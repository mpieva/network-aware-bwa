/* Estimation of insert sizes.  We want to do this incrementally and
 * globally, and separately for the read groups.
 *
 * - infer_isize starts with some filtering, then determines quartiles
 *   to select a hard cutoff.
 *   - This can be done on a small sample (~1000 reads), because
 *     precision is not really important.
 * - it then estimates average, variance, kurtosis and skewness, but
 *   only from the stuff between the cutoffs
 *   - Therefore, incremental computation is out the window and an
 *     estimate from the same subsample must be sufficient.
 * - then a cutoff is computed so that a fraction ap_prior is above the
 *   cutoff [contains Voodoo]
 * - then the actual fraction above is counted and used instead of
 *   ap_prior, it it is bigger
 *   - This will need a sample of the larger ISIZE fraction at higher
 *     resolution.  What constitutes the higher fraction can be defined
 *     based on an early eastimate of avg and std.
 *
 * Another option is to simply collect insert sizes until our sample is
 * big enough, do the estimate and be done with it.  Limit ISIZE to
 * 65536 and it's actually pretty compact.  And the old code can be
 * recycled (mwahaha!).  This way, we need a few hundred thousand
 * samples at most.
 *
 * Another option: keep the histogram, maybe limited to 64k entries per
 * bin.  Put the upper limit to 64k, we're looking at only 128k of
 * storage per read group.  Fine for thousands of read groups.
 *
 *
 * Therefore, the detailed plan:  Our knowledge of insert sizes is a
 * hash map with one entry per read group.  We start out with a
 * histogram, which is an array 0..65535 of shorts, counting the
 * respective insert size.  At the end of the input, or when one of the
 * bins overflows, we estimate the insert size.  Now we throw the
 * histogram away and store an isize_info_t.  Those get broadcast and
 * distributed with the config to the workers.
 *
 * We do not store insert sizes < 4; this makes sure the first 8 bytes
 * of our histogram are 0.  That way a histogram can be distinguished 
 */

#include "bwape.h"
#include <math.h>

#define MAX_ISIZE 100000

// ap_prior?  L?
static int infer_isize_hist( const char *rg, isize_info_t *ii, double ap_prior, int64_t L)
{
	uint64_t x, n_ap = 0;
    // max len becomes maximum read length -- not ISIZE!
    // p25, p50, p75 becomes quartiles of ISIZE distribution
	int n, i, tot=0, cum=0, p25=0, p50=0, p75=0, tmp;
	double skewness = 0.0, kurtosis = 0.0, y;

    // avg becomes average of stuff between low and high.  (is this
    // still unbiased?  i think so.)
	ii->avg = ii->std = -1.0;
	ii->low = ii->high = ii->high_bayesian = 0;

    for( i=0; i != MAX_ISIZE; ++i ) tot += ii->hist[i] ;

	if (tot < 20) {
        fprintf(stderr, "[infer_isize] %s: too few good pairs\n", rg?rg:"(null)");
		free(ii->hist);
        ii->hist=0;
		return -1;
	}

    for( i=0; i != MAX_ISIZE; ++i ) {
        int cum2 = cum + ii->hist[i] ;
        if( cum <= tot*0.25+0.5 && cum2 > tot*0.25+0.5 ) p25 = i ;
        if( cum <= tot*0.50+0.5 && cum2 > tot*0.50+0.5 ) p50 = i ;
        if( cum <= tot*0.75+0.5 && cum2 > tot*0.75+0.5 ) p75 = i ;
        cum = cum2 ;
    }

	tmp  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
    // Clamp to positive number.  Used to be the maximum read length,
    // but that doesn't make sense if inserts can actually be shorter
    // than a read (and they usually can).
	ii->low = tmp > 1 ? tmp : 1 ;
	ii->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);

    // first pass: number in range and average
	for (i = 0, x = n = 0; i < MAX_ISIZE; ++i) {
        if (i >= ii->low && i <= ii->high) {
            n += ii->hist[i] ;
            x += ii->hist[i] * i ;
        }
    }
	ii->avg = (double)x / n;        // ...gives average of everything between low and high

    // second pass; higher moments
	for (i = 0; i < MAX_ISIZE; ++i) {
		if (i >= ii->low && i <= ii->high) {
			double tmp = (i - ii->avg) * (i - ii->avg);
			ii->std += tmp * ii->hist[i];
			skewness += tmp * (i - ii->avg) * ii->hist[i];
			kurtosis += tmp * tmp * ii->hist[i];
		}
	}
    // everything that follows is stolen directly from infer_isize
	kurtosis = kurtosis/n / (ii->std / n * ii->std / n) - 3;
	ii->std = sqrt(ii->std / n);
	skewness = skewness / n / (ii->std * ii->std * ii->std);

    // Chose high_bayesian such that P( ISIZE > high_bayesian ) ==
    // ap_prior, assuming a normal distribution.
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (bwtint_t)(y * ii->std + ii->avg + .499);
    
	for (i = 0; i < MAX_ISIZE; ++i)
		if (i > ii->high_bayesian)
            n_ap += ii->hist[i];  // count how many are too high

    free(ii->hist); ii->hist=0;
    
    // estimate rate of abnormal pairs, use the estimate iff it is
    // higher than the prior.  (Is this expressed in percent?  And we
    // add a hundredth of a percent to get the rounding right or
    // something?)
	ii->ap_prior = .01 * (n_ap + .01) / tot;
	if (ii->ap_prior < ap_prior) ii->ap_prior = ap_prior;

	fprintf(stderr, "[infer_isize] %s: qu(%d, %d, %d)", rg?rg:"(null)", p25, p50, p75);
	if (isnan(ii->std) || p75 > MAX_ISIZE) {
		ii->low = ii->high = ii->high_bayesian = 0; ii->avg = ii->std = -1.0;
		fprintf(stderr, " -- not useable\n");
		return -1;
	}

	fprintf(stderr, " bound(%d,%d), num/avg/std/kur/skw %d/%.3lf/%.3lf/%.3lf/%.3lf, ap %.2e, max %d, %.2lf sigma\n",
            ii->low, ii->high, n, ii->avg, ii->std, skewness, kurtosis, ii->ap_prior, ii->high_bayesian, y);
	return 0;
}

void improve_isize_est(khash_t(isize_infos) *iinfos, bam_pair_t *p, double ap_prior, int64_t L)
{
    bwa_seq_t *s = p->bwa_seq ;
    if( p->kind<1 || s[0].mapQ < 20 ) return ;
    if( p->kind>1 && s[1].mapQ < 20 ) return ;

    int len = p->kind==1 ? s[0].len :       // single read: isize is length (trimming/merging is assumed)
        s[0].pos < s[1].pos ? s[1].pos + s[1].len - s[0].pos :
        s[0].pos + s[0].len - s[1].pos;
    if( len < 0 || len >= MAX_ISIZE ) return ;

    int ret = 0 ;
    const char *rg = bam_get_rg(p->bam_rec);
    khiter_t it = kh_get(isize_infos, iinfos, rg) ;
    isize_info_t *ii = &kh_value(iinfos, it) ;
    if( it==kh_end(iinfos) ) { 
        it = kh_put(isize_infos, iinfos, strdup(rg), &ret) ;
        ii = &kh_value(iinfos, it) ;
        memset( ii, 0, sizeof(isize_info_t) ) ;
        ii->hist = calloc( MAX_ISIZE, sizeof(unsigned short) ) ;
    }
    if( !ii->hist ) return ;
    // if we hit the ceiling when counting, infer isize and be done
    if( ++ii->hist[len] == (-1) ) infer_isize_hist(rg, ii, ap_prior, L) ;
}

void infer_all_isizes( khash_t(isize_infos) *iinfos, double ap_prior, int64_t L)
{
    khiter_t k;
    for(k = kh_begin(iinfos); k != kh_end(iinfos); ++k)
        if(kh_exist(iinfos,k) && kh_value(iinfos,k).hist )
            infer_isize_hist(kh_key(iinfos,k), &kh_value(iinfos,k), ap_prior, L) ;
}
