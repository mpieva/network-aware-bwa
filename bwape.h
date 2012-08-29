#ifndef INCLUDED_BWAPE_H
#define INCLUDED_BWAPE_H

#include "bwtaln.h"
#include "kvec.h"
#include "bntseq.h"

typedef struct {
	int n;
	bwtint_t *a;
} poslist_t;

typedef struct {
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

extern int g_log_n[256]; // in bwase.c

typedef struct {
	kvec_t(bwt_aln1_t) aln;
} aln_buf_t;

int pairing(bwa_seq_t *p[2], pe_data_t *d, const pe_opt_t *opt, int s_mm, const isize_info_t *ii) ;
ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *pac, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii);
pe_opt_t *bwa_init_pe_opt() ;

#endif

