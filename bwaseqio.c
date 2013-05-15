#include <zlib.h>
#include <ctype.h>
#include "bwtaln.h"
#include "utils.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];
static char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

struct __bwa_seqio_t {
	// for BAM input
	int is_bam, which; // 1st bit: read1, 2nd bit: read2, 3rd: SE
	bamFile fp;
	// for fastq input
	kseq_t *ks;
    // in case sai failes available
    FILE *sai[3];
};

bwa_seqio_t *bwa_bam_open(const char *fn, int which, char **saif,
                          gap_opt_t *o0, bam_header_t **hh)
{
    int c, b=0;
	bwa_seqio_t *bs;
	bam_header_t *h;
	bs = (bwa_seqio_t*)calloc(1, sizeof(bwa_seqio_t));
	bs->is_bam = 1;
	bs->which = which;
    bs->fp = (fn[0]!='-' || fn[1]) ? bam_open(fn, "r") : bam_dopen(0, "r") ;
	h = bam_header_read(bs->fp);
	if(hh) *hh=h; else bam_header_destroy(h);

    if( saif ) for(c=0;c!=3;++c)
    {
        gap_opt_t opt;
        if( saif[c] ) {
            bs->sai[c] = xopen(saif[c], "r");
            if( 1 > fread(&opt, sizeof(gap_opt_t), 1, bs->sai[c]) )
            {
                fclose(bs->sai[c]); 
                bs->sai[c] = 0;
            }
            opt.n_threads=o0->n_threads;
            if(o0) {
                if(b) {
                    opt.mode=o0->mode;
                    if( memcmp(o0, &opt, sizeof(gap_opt_t)) ) {
                        fprintf( stderr, "[bwa_bam_open] options from sai file \"%s\" conflict with others.\n", saif[c] ) ;
                        exit(1);
                    }
                    fprintf( stderr, "[bwa_bam_open] options from sai file \"%s\" match.\n", saif[c] ) ;
                }
                else {
                    fprintf( stderr, "[bwa_bam_open] recovered options from sai file \"%s\".\n", saif[c] ) ;
                    memcpy(o0, &opt, sizeof(gap_opt_t));
                    b=1;
                }
            }
        }
    }
	return bs;
}

bwa_seqio_t *bwa_seq_open(const char *fn)
{
	gzFile fp;
	bwa_seqio_t *bs;
	bs = (bwa_seqio_t*)calloc(1, sizeof(bwa_seqio_t));
	fp = xzopen(fn, "r");
	bs->ks = kseq_init(fp);
	return bs;
}

void bwa_seq_close(bwa_seqio_t *bs)
{
    int i;
	if (bs == 0) return;
	if (bs->is_bam) bam_close(bs->fp);
	else {
		gzclose(bs->ks->f->f);
		kseq_destroy(bs->ks);
	}
    for(i=0;i!=3;++i)
        if(bs->sai[i])
            fclose(bs->sai[i]);
	free(bs);
}

void seq_reverse(int len, ubyte_t *seq, int is_comp)
{
	int i;
	if (is_comp) {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			if (tmp < 4) tmp = 3 - tmp;
			seq[len-1-i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
			seq[i] = tmp;
		}
		if (len&1) seq[i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
	} else {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			seq[len-1-i] = seq[i]; seq[i] = tmp;
		}
	}
}

int bwa_trim_read(int trim_qual, bwa_seq_t *p)
{
	int s = 0, l, max = 0, max_l = p->len - 1;
	if (trim_qual < 1 || p->qual == 0) return 0;
	for (l = p->len - 1; l >= BWA_MIN_RDLEN - 1; --l) {
		s += trim_qual - (p->qual[l] - 33);
		if (s < 0) break;
		if (s > max) {
			max = s; max_l = l;
		}
	}
	p->clip_len = p->len = max_l + 1;
	return p->full_len - p->len;
}

bwa_seq_t *bwa_read_bam(bwa_seqio_t *bs, int n_needed, int *n, int is_comp, int trim_qual)
{
	bwa_seq_t *seqs, *p;
	int n_seqs, l, i;
	long n_trimmed = 0, n_tot = 0;
	bam1_t *b;

	b = bam_init1();
	n_seqs = 0;
	seqs = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
	while (bam_read1(bs->fp, b) >= 0) {
		uint8_t *s, *q;
		int go = 0;
		if ((bs->which & 1) &&  (b->core.flag & BAM_FPAIRED) && (b->core.flag & BAM_FREAD1)) go = 1;
		if ((bs->which & 2) &&  (b->core.flag & BAM_FPAIRED) && (b->core.flag & BAM_FREAD2)) go = 1;
		if ((bs->which & 4) && !(b->core.flag & BAM_FPAIRED)) go = 1;
		if (go == 0) continue;
		l = b->core.l_qseq;
		p = &seqs[n_seqs++];
		p->tid = -1; // no assigned to a thread
		p->qual = 0;
		p->full_len = p->clip_len = p->len = l;
		n_tot += p->full_len;
		s = bam1_seq(b); q = bam1_qual(b);
		p->seq = (ubyte_t*)calloc(p->len + 1, 1);
		p->qual = (ubyte_t*)calloc(p->len + 1, 1);
		for (i = 0; i != p->full_len; ++i) {
			p->seq[i] = bam_nt16_nt4_table[(int)bam1_seqi(s, i)];
			p->qual[i] = q[i] + 33 < 126? q[i] + 33 : 126;
		}
		if (bam1_strand(b)) { // then reverse 
			seq_reverse(p->len, p->seq, 1);
			seq_reverse(p->len, p->qual, 0);
		}
		if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
		p->rseq = (ubyte_t*)calloc(p->full_len, 1);
		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
		seq_reverse(p->len, p->rseq, is_comp);
		p->name = strdup((const char*)bam1_qname(b));
		if (n_seqs == n_needed) break;
	}
	*n = n_seqs;
	if (n_seqs && trim_qual >= 1)
		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed/n_tot);
	if (n_seqs == 0) {
		free(seqs);
		bam_destroy1(b);
		return 0;
	}
	bam_destroy1(b);
	return seqs;
}

#define BARCODE_LOW_QUAL 13

bwa_seq_t *bwa_read_seq(bwa_seqio_t *bs, int n_needed, int *n, int mode, int trim_qual)
{
	bwa_seq_t *seqs, *p;
	kseq_t *seq = bs->ks;
	int n_seqs, l, i, is_comp = mode&BWA_MODE_COMPREAD, is_64 = mode&BWA_MODE_IL13, l_bc = mode>>24;
	long n_trimmed = 0, n_tot = 0;

	if (l_bc > BWA_MAX_BCLEN) {
		fprintf(stderr, "[%s] the maximum barcode length is %d.\n", __func__, BWA_MAX_BCLEN);
		return 0;
	}
	if (bs->is_bam) return bwa_read_bam(bs, n_needed, n, is_comp, trim_qual); // l_bc has no effect for BAM input
	n_seqs = 0;
	seqs = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
	while ((l = kseq_read(seq)) >= 0) {
		if ((mode & BWA_MODE_CFY) && (seq->comment.l != 0)) {
			// skip reads that are marked to be filtered by Casava
			char *s = index(seq->comment.s, ':');
			if (s && *(++s) == 'Y') {
				continue;
			}
		}
		if (is_64 && seq->qual.l)
			for (i = 0; i < seq->qual.l; ++i) seq->qual.s[i] -= 31;
		if (seq->seq.l <= l_bc) continue; // sequence length equals or smaller than the barcode length
		p = &seqs[n_seqs++];
		if (l_bc) { // then trim barcode
			for (i = 0; i < l_bc; ++i)
				p->bc[i] = (seq->qual.l && seq->qual.s[i]-33 < BARCODE_LOW_QUAL)? tolower(seq->seq.s[i]) : toupper(seq->seq.s[i]);
			p->bc[i] = 0;
			for (; i < seq->seq.l; ++i)
				seq->seq.s[i - l_bc] = seq->seq.s[i];
			seq->seq.l -= l_bc; seq->seq.s[seq->seq.l] = 0;
			if (seq->qual.l) {
				for (i = l_bc; i < seq->qual.l; ++i)
					seq->qual.s[i - l_bc] = seq->qual.s[i];
				seq->qual.l -= l_bc; seq->qual.s[seq->qual.l] = 0;
			}
			l = seq->seq.l;
		} else p->bc[0] = 0;
		p->tid = -1; // no assigned to a thread
		p->qual = 0;
		p->full_len = p->clip_len = p->len = l;
		n_tot += p->full_len;
		p->seq = (ubyte_t*)calloc(p->len, 1);
		for (i = 0; i != p->full_len; ++i)
			p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
		if (seq->qual.l) { // copy quality
			p->qual = (ubyte_t*)strdup((char*)seq->qual.s);
			if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
		}
		p->rseq = (ubyte_t*)calloc(p->full_len, 1);
		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
		seq_reverse(p->len, p->rseq, is_comp);
		p->name = strdup((const char*)seq->name.s);
		{ // trim /[12]$
			int t = strlen(p->name);
			if (t > 2 && p->name[t-2] == '/' && (p->name[t-1] == '1' || p->name[t-1] == '2')) p->name[t-2] = '\0';
		}
		if (n_seqs == n_needed) break;
	}
	*n = n_seqs;
	if (n_seqs && trim_qual >= 1)
		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed/n_tot);
	if (n_seqs == 0) {
		free(seqs);
		return 0;
	}
	return seqs;
}

void bwa_free_read_seq1(bwa_seq_t *p)
{
	int j;
    for (j = 0; j < p->n_multi; ++j)
        if (p->multi[j].cigar) free(p->multi[j].cigar);
    free(p->name);
    free(p->seq); free(p->rseq); free(p->qual); free(p->aln); free(p->md); free(p->multi);
    free(p->cigar);
}

void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs)
{
	int i;
	for (i = 0; i != n_seqs; ++i)
        bwa_free_read_seq1( seqs+i ) ;
	free(seqs);
}

// Mostly stolen from bwa_read_bam.
void bam1_to_seq(bam1_t *raw, bwa_seq_t *p, int is_comp, int trim_qual)
{
    // long n_trimmed = 0;

    uint8_t *s, *q;
    int i, l = raw->core.l_qseq;
    p->tid = -1; // no assigned to a thread
    p->qual = 0;
    p->full_len = p->clip_len = p->len = l;
    // n_tot += p->full_len;
    s = bam1_seq(raw); q = bam1_qual(raw);
    p->seq = (ubyte_t*)calloc(p->len + 1, 1);
    p->qual = (ubyte_t*)calloc(p->len + 1, 1);
    for (i = 0; i != p->full_len; ++i) {
        p->seq[i] = bam_nt16_nt4_table[(int)bam1_seqi(s, i)];
        p->qual[i] = q[i] + 33 < 126? q[i] + 33 : 126;
    }
    if (bam1_strand(raw)) { // then reverse 
        seq_reverse(p->len, p->seq, 1);
        seq_reverse(p->len, p->qual, 0);
    }
    if (trim_qual >= 1) /* n_trimmed += */ bwa_trim_read(trim_qual, p);
    p->rseq = (ubyte_t*)calloc(p->full_len, 1);
    memcpy(p->rseq, p->seq, p->len);
    seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
    seq_reverse(p->len, p->rseq, is_comp);
    p->max_entries = 0 ;

    // We don't set a name, it's contained in the original record
    // anyway.
    // p->name = strdup((const char*)bam1_qname(raw));

    // No place to put the tally right now.
    // if (n_seqs && trim_qual >= 1)
    // fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed/n_tot);
}

static void memswap(void *pp, void *qq, size_t s)
{
    char *p = (char*)pp, *q = (char*)qq;
    while(s)
    {
        char x = *p ;
        *p = *q ;
        *q = x ;
        --s ;
        ++p ; 
        ++q ;
    }
}

int try_get_sai( FILE **f, int c, int *naln, bwt_aln1_t **aln )
{
    if( f && f[c] ) {
        if( 1 == fread(naln, 4, 1, f[c]) ) {
            *aln = (bwt_aln1_t*)calloc(*naln, sizeof(bwt_aln1_t));
            if( *naln == fread(*aln, sizeof(bwt_aln1_t), *naln, f[c]) ) return 1 ;
            free(*aln) ;
        }
        fprintf( stderr, "[read_bam_pair] note: sai file %d has ended.\n", c ) ;
        fclose(f[c]);
        f[c]=0;
    }
    *aln = 0 ;
    *naln = 0 ;
    return 0 ;
}

/* Read one pair from a bam file.  Returns 1 if we got a singleton, 2 if
 * we got a pair, 0 if we reached EOF, -1 if something outside our
 * control went wrong, -2 if we got something unexpected (missing mate,
 * fragment with unexpected PE flags).
 */
static int read_bam_pair_core(bwa_seqio_t *bs, bam_pair_t *pair, int allow_broken)
{
    static int num_wrong_pair = 128 ;

    memset(pair, 0, sizeof(bam_pair_t)) ;
	if (bam_read1(bs->fp, &pair->bam_rec[0]) < 0) return 0 ;
    while(1) {
        if (pair->bam_rec[0].core.flag & BAM_FPAIRED) { // paired read, get another
            if (bam_read1(bs->fp, &pair->bam_rec[1]) >= 0) {
                uint32_t flag1 = pair->bam_rec[0].core.flag & (BAM_FPAIRED|BAM_FREAD1|BAM_FREAD2);
                uint32_t flag2 = pair->bam_rec[1].core.flag & (BAM_FPAIRED|BAM_FREAD1|BAM_FREAD2);
                if (!strcmp(bam1_qname(&pair->bam_rec[0]), bam1_qname(&pair->bam_rec[1]))) { // actual mates
                    if( flag1 == (BAM_FPAIRED|BAM_FREAD1) && flag2 == (BAM_FPAIRED|BAM_FREAD2) ) { // correct order
                        pair->kind = proper_pair ;
                        return 2 ;
                    } else if (flag2 == (BAM_FPAIRED|BAM_FREAD1) && flag1 == (BAM_FPAIRED|BAM_FREAD2) ) { // reverse order
                        memswap(&pair->bam_rec[0], &pair->bam_rec[1], sizeof(bam1_t));
                        pair->kind = proper_pair ;
                        return 2 ;
                    } else {
                        fprintf( stderr, "[read_bam_pair] got a pair, but the flags are wrong (%s).\n", bam1_qname(&pair->bam_rec[0]) ) ;
                        if( allow_broken ) {
                            pair->bam_rec[0].core.flag &= ~BAM_FREAD2;
                            pair->bam_rec[0].core.flag |= BAM_FPAIRED|BAM_FREAD1;
                            pair->bam_rec[1].core.flag &= ~BAM_FREAD1;
                            pair->bam_rec[1].core.flag |= BAM_FPAIRED|BAM_FREAD2;
                            pair->kind = proper_pair ;
                            return 2 ;
                        }
                        else return -2 ;
                    }
                } else {
                    // This is arguably wrong, we discard a lone mate.  But what else could we do?  Buffering it
                    // somewhere to wait is too hard for the time being, returning it as a single means we need to buffer the
                    // next one.  Not very appealing.  So only two options remain: discard it or bail out.
                    if( num_wrong_pair ) {
                        fprintf( stderr, "[read_bam_pair] got two reads, but the names don't match (%s,%s).\n",
                                bam1_qname(&pair->bam_rec[0]), bam1_qname(&pair->bam_rec[1]) ) ;
                        --num_wrong_pair ;
                        if( !num_wrong_pair ) 
                            fprintf( stderr, "[read_bam_pair] too many mismatched names, not reporting anymore.\n" ) ;
                    }
                    try_get_sai( bs->sai, flag1 & BAM_FREAD1 ? 1 : 2, &pair->bwa_seq[0].n_aln, &pair->bwa_seq[0].aln ) ;
                    free(pair->bam_rec[0].data);
                    if(pair->bwa_seq[0].n_aln) free(pair->bwa_seq[0].aln);
                    if( !allow_broken ) {
                        free(pair->bam_rec[1].data);
                        if(pair->bwa_seq[0].n_aln) free(pair->bwa_seq[0].aln);
                        return -2 ;
                    }
                    memmove(&pair->bam_rec[0], &pair->bam_rec[1], sizeof(bam1_t));
                    memset(&pair->bam_rec[1], 0, sizeof(bam1_t));
                }
            } else {
                fprintf( stderr, "[read_bam_pair] got a paired read and hit EOF.\n" ) ;
                free(pair->bam_rec[0].data);
                if(pair->bwa_seq[0].n_aln) free(pair->bwa_seq[0].aln);
                return allow_broken ? 0 : -2 ;
            }
        } else { // singleton read
            pair->kind = singleton ;
            return 1 ;
        }
    }
}

// Erase unwanted tagged fields:
// AM NM CM SM MD X0 X1 XA XC XG XM XN XO XT YQ
static void erase_unwanted_tags(bam1_t *out)
{
    int total = out->core.n_cigar*4 + out->core.l_qname +
        out->core.l_qseq + (out->core.l_qseq + 1)/2;
    uint8_t *p = out->data + total;
    uint8_t *q = p ;

    while (p < out->data + out->data_len) {
        int keep = 1 ;
        switch(p[0]) {
            case 'A':
            case 'S': 
            case 'C':
            case 'N': keep = p[1] != 'M' ; break ;
            case 'M': keep = p[1] != 'D' ; break ;
            case 'X': keep = !strchr("01ACGMNOT", p[1]) ; break ;
            case 'Y': keep = p[1] != 'Q' ; break ;
        }

        int len = 3 ;
        int count ;
        switch( p[2] & ~32 ) {
            case 'C':
            case 'A': len++; break;
            case 'S': len+=2; break;
            case 'I':
            case 'F': len+=4; break;
            case 'D': len+=8; break;
            case 'Z':
            case 'H': while( p[len] ) len++ ; len++; break;
            case 'B': count=(int)p[4] << 0 | (int)p[5] << 8 | (int)p[6] << 16 | (int)p[7] << 24 ;
                      len += 5;
                      switch( p[3] & ~32 ) {
                          case 'C':
                          case 'A': len += count ; break;
                          case 'S': len += 2*count; break;
                          case 'I':
                          case 'F': len += 4*count ; break;
                          case 'D': len += 8*count ; break;
                      }
                      break;
        }

        if(keep) {
            memmove(q, p, len);
            q+=len;
            total+=len;
        }
        p+=len;
	}
    out->data_len=total;
}

int read_bam_pair(bwa_seqio_t *bs, bam_pair_t *pair, int allow_broken)
{
    int i, r = read_bam_pair_core( bs, pair, allow_broken ) ;

    if( pair->kind == singleton ) {
        if( try_get_sai( bs->sai, 0, &pair->bwa_seq[0].n_aln, &pair->bwa_seq[0].aln ) )
            pair->phase = aligned ;
    }
    else if( pair->kind == proper_pair ) {
        if( try_get_sai( bs->sai, 1, &pair->bwa_seq[0].n_aln, &pair->bwa_seq[0].aln ) 
                + try_get_sai( bs->sai, 2, &pair->bwa_seq[1].n_aln, &pair->bwa_seq[1].aln ) == 2 )
            pair->phase = aligned ;
    }

    // make sure that either none or both pass QC
    if( pair->kind == proper_pair ) {
        pair->bam_rec[0].core.flag |= pair->bam_rec[1].core.flag & SAM_FQC ;
        pair->bam_rec[1].core.flag |= pair->bam_rec[0].core.flag & SAM_FQC ;
    }

    for( i = 0 ; i != pair->kind ; ++i )
        erase_unwanted_tags( &pair->bam_rec[i] ) ;
    return r ;
}
