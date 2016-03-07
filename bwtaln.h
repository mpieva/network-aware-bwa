#ifndef BWTALN_H
#define BWTALN_H

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "bwt.h"
#include "bamlite.h"

#define BWA_TYPE_NO_MATCH 0
#define BWA_TYPE_UNIQUE 1
#define BWA_TYPE_REPEAT 2
#define BWA_TYPE_MATESW 3

#define SAM_FPD   1 // paired
#define SAM_FPP   2 // properly paired
#define SAM_FSU   4 // self-unmapped
#define SAM_FMU   8 // mate-unmapped
#define SAM_FSR  16 // self on the reverse strand
#define SAM_FMR  32 // mate on the reverse strand
#define SAM_FR1  64 // this is read one
#define SAM_FR2 128 // this is read two
#define SAM_FSC 256 // secondary alignment
#define SAM_FQC 512 // fails quality control
#define SAM_FDP 1024 // is a duplicate

#define BWA_AVG_ERR 0.02
#define BWA_MIN_RDLEN 35 // for read trimming

#define BWA_MAX_BCLEN 63 // maximum barcode length; 127 is the maximum

#ifndef bns_pac
#define bns_pac(pac, k) ((pac)[(k)>>2] >> ((~(k)&3)<<1) & 3)
#endif

typedef struct {
	bwtint_t w;
	int bid;
} bwt_width_t;

typedef struct {
	uint32_t n_mm:8, n_gapo:8, n_gape:8, a:1;
	bwtint_t k, l;
	int score;
} bwt_aln1_t;

typedef uint16_t bwa_cigar_t;
/* rgoya: If changing order of bytes, beware of operations like:
 *     s->cigar[0] += s->full_len - s->len;
 */
#define CIGAR_OP_SHIFT 14
#define CIGAR_LN_MASK 0x3fff

#define __cigar_op(__cigar) ((__cigar)>>CIGAR_OP_SHIFT)
#define __cigar_len(__cigar) ((__cigar)&CIGAR_LN_MASK)
#define __cigar_create(__op, __len) ((__op)<<CIGAR_OP_SHIFT | (__len))

typedef struct {
	uint32_t pos;
	uint32_t n_cigar:15, gap:8, mm:8, strand:1;
	bwa_cigar_t *cigar;
} bwt_multi1_t;

typedef struct {
	char *name;
	ubyte_t *seq, *rseq, *qual;
	uint32_t len:20, strand:1, type:2, dummy:1, extra_flag:8;
	uint32_t n_mm:8, n_gapo:8, n_gape:8, mapQ:8;
	int score;
	int clip_len;
	// alignments in SA coordinates
	int n_aln;
	bwt_aln1_t *aln;
	// multiple hits
	int n_multi;
	bwt_multi1_t *multi;
	// alignment information
	bwtint_t sa, pos;
	uint64_t c1:28, c2:28, seQ:8; // number of top1 and top2 hits; single-end mapQ
	int n_cigar;
	bwa_cigar_t *cigar;
	// for multi-threading only
	int tid;
	// barcode
	char bc[BWA_MAX_BCLEN+1]; // null terminated; up to BWA_MAX_BCLEN bases
	// NM and MD tags
	uint32_t full_len:20, nm:12;
	char *md;
    int max_entries;
} bwa_seq_t;

/* The structure of the existing code is sufficiently twisted and warped
 * that I don't see how to untangle it.  So, for bam2bam, reboot the
 * structure.  We start be reading a logical record from the input file:
 * read a sequence, and if it's single ended, keep it as such.  If it's
 * a first or second read, get another one, it should be the other half.
 * Exchange them if necessary, bomb out if one is missing.
 *
 * We further want to keep everything from the original bam file, so we
 * can reproduce it on output.  Therefore, we simply keep the original
 * record.  Here's the new data structure.  If kind is set to singleton,
 * the second read is invalid.
 */

/** Describes the kind of data contained in bam_pair_t, which makes up a
 * logical record.
 */
enum pair_kind {
    eof_marker=0,   // end-of-file marker, no paylod
    singleton=1,    // singleton, only bam_rec[0] and bwa_seq[0] are valid
    proper_pair=2   // a pair, {bam_rec,bwa_seq}[0] is forward, the other reverse
} ;

/** Describes the processing phase of a logcal bam record (bam_pair_t),
 * and thereby determines which fields are valid.
 */
enum pair_phase { 
    pristine=0,     // fresh input, only bam_rec[] is valid
    aligned=1,      // aligned, bam_rec[] and bwa_seq[].aln are valid
    positioned=2,   // coordinates computed, bam_rec[] and part of bwa_seq[] are valid
    finished=3      // everything done, only bam_rec[] is valid
} ;

typedef struct {
    uint64_t recno ;             // to put them back into the correct order
    enum pair_kind kind ;
    enum pair_phase phase ;
    bam1_t bam_rec[2] ;
    bwa_seq_t bwa_seq[2] ;
} bam_pair_t ;

#define BWA_MODE_GAPE       0x01
#define BWA_MODE_COMPREAD   0x02
#define BWA_MODE_LOGGAP     0x04
#define BWA_MODE_CFY        0x08
#define BWA_MODE_NONSTOP    0x10
#define BWA_MODE_BAM        0x20
#define BWA_MODE_BAM_SE     0x40
#define BWA_MODE_BAM_READ1  0x80
#define BWA_MODE_BAM_READ2  0x100
#define BWA_MODE_IL13       0x200

typedef struct {
	int s_mm, s_gapo, s_gape;
	int mode; // bit 24-31 are the barcode length
	int indel_end_skip, max_del_occ, max_entries;
	float fnr;
	int max_diff, max_gapo, max_gape;
	int max_seed_diff, seed_len;
	int n_threads;
	int max_top2;
	int trim_qual;
} gap_opt_t;

#define BWA_PET_STD   1
#define BWA_PET_SOLID 2

typedef struct {
	int max_isize, force_isize;
	int max_occ, max_occ_se;
	int n_multi, N_multi;
	int type, is_sw, is_preload;
	double ap_prior;
} pe_opt_t;

struct __bwa_seqio_t;
typedef struct __bwa_seqio_t bwa_seqio_t;

#ifdef __cplusplus
extern "C" {
#endif

	gap_opt_t *gap_init_opt();
	void bwa_aln_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt, int nskip);

	bwa_seqio_t *bwa_seq_open(const char *fn);
	bwa_seqio_t *bwa_bam_open(const char *fn, int which, char **saif, gap_opt_t *o0, bam_header_t **hh);
	void bwa_seq_close(bwa_seqio_t *bs);
	void seq_reverse(int len, ubyte_t *seq, int is_comp);
    int read_bam_pair(bwa_seqio_t *bs, bam_pair_t *pair, int allow_broken, int drop_aligned);
    void bam1_to_seq(bam1_t *raw, bwa_seq_t *p, int is_comp, int trim_qual);
	bwa_seq_t *bwa_read_seq(bwa_seqio_t *seq, int n_needed, int *n, int mode, int trim_qual);
    void bwa_free_read_seq1(bwa_seq_t *p);
	void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs);

	int bwa_cal_maxdiff(int l, double err, double thres);
	void bwa_cal_sa_reg_gap(bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);

	void bwa_cs2nt_core(bwa_seq_t *p, bwtint_t l_pac, ubyte_t *pac);


	/* rgoya: Temporary clone of aln_path2cigar to accomodate for bwa_cigar_t,
	__cigar_op and __cigar_len while keeping stdaln stand alone */
#include "stdaln.h"

	bwa_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar);

#ifdef __cplusplus
}
#endif

static inline void bam_init_pair( bam_pair_t *p )
{
    memset(p, 0, sizeof(bam_pair_t));
}

static inline void bam_destroy_pair( bam_pair_t *p )
{
    int i ;
    if(p) {
        for( i = 0 ; i != p->kind ; ++i ) {
            free(p->bam_rec[i].data);
            switch( p->phase ) {
                case aligned:
                case positioned:
                    bwa_free_read_seq1( &p->bwa_seq[i] ) ;
                case pristine:
                case finished:
                    break ;
            }
        }
    }
}


#endif
