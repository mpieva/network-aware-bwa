#ifndef INCLUDED_BWAPE_H
#define INCLUDED_BWAPE_H

#include "bwtaln.h"
#include "kvec.h"
#include "bntseq.h"

typedef struct {
	int n;
	bwtint_t *a;
} poslist_t;

// When running in two passes, 'hist' is first filled in, but then
// cleared once there is enough data to compute the other fields.
// Therefore, there is no need to pass an array to pass two.
typedef struct {
    unsigned short *hist ;
	double avg, std, ap_prior;
	bwtint_t low, high, high_bayesian;
} isize_info_t;

#include "khash.h"
#include "ksort.h"

typedef struct {
	kvec_t(uint64_t) arr;
	kvec_t(uint64_t) pos[2];
	kvec_t(bwt_aln1_t) aln[2];
} pe_data_t;

#define MIN_HASH_WIDTH 1000

// for normal distribution, this is about 3std
#define OUTLIER_BOUND 2.0

#define SW_MIN_MATCH_LEN 20
#define SW_MIN_MAPQ 17


extern int g_log_n[256]; // in bwase.c

typedef struct {
	kvec_t(bwt_aln1_t) aln;
} aln_buf_t;

int pairing(bwa_seq_t *p[2], pe_data_t *d, const pe_opt_t *opt, int s_mm, const isize_info_t *ii) ;
ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *pac, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii);
pe_opt_t *bwa_init_pe_opt() ;
void bwa_paired_sw1(const bntseq_t *bns, const ubyte_t *pacseq, bwa_seq_t *p[2], const pe_opt_t *popt, const isize_info_t *ii,uint64_t n_tot[2], uint64_t n_mapped[2]);

KHASH_MAP_INIT_STR(isize_infos, isize_info_t) ;

void improve_isize_est(khash_t(isize_infos) *iinfos, bam_pair_t *p, double ap_prior, int64_t);
void infer_all_isizes( khash_t(isize_infos) *iinfos, double ap_prior, int64_t);
int iinfo_encoded_size( khash_t(isize_infos) *iinfos ) ;
char *encode_iinfo( char* p, khash_t(isize_infos) *iinfos ) ;
khash_t(isize_infos) *decode_iinfo( char *p, char *q ) ;

#endif

