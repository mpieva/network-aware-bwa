/* Combine aln, samse, sampe into one workflow.  We steal the workhorse
 * functions from the other subprograms, then wrap them in new top-level
 * functions.  Comes with a command line that combines all the other
 * options, maybe even some more.
 */

static const int chunksize = 0x100000;
static const int loudness  = 0;
static const int ring_size = 0x280000;
static const int timeout   = 30;

#include "bamlite.h"
#include "bwtaln.h"
#include "bwase.h"
#include "bwape.h"
#include "khash.h"
#include "main.h"
#include "bgzf.h"
#include "utils.h"

#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/utsname.h>

#include <zmq.h>
#include <pthread.h>
#include <signal.h>
#include <stdarg.h>
#include <zlib.h>

static const int the_hwm = 64;
static const int the_linger = 2000;

KHASH_MAP_INIT_INT64(64, poslist_t)

struct option longopts[] = {
    { "num-diff",               1, 0, 'n' },
    { "max-gap-open",           1, 0, 'o' },
    { "max-gap-extensions",     1, 0, 'e' },
    { "indel-near-end",         1, 0, 'i' },
    { "deletion-occurences",    1, 0, 'd' },
    { "seed-length",            1, 0, 'l' },
    { "seed-mismatches",        1, 0, 'k' },
    { "queue-size",             1, 0, 'm' },
    { "num-threads",            1, 0, 't' },
    { "mismatch-penalty",       1, 0, 'M' },
    { "gap-open-penalty",       1, 0, 'O' },
    { "gap-extension-penalty",  1, 0, 'E' },
    { "max-best-hits",          1, 0, 'R' },
    { "trim-quality",           1, 0, 'q' },
    { "log-gap-penalty",        0, 0, 'L' },
    { "non-iterative",          0, 0, 'N' },
    { "barcode-length",         1, 0, 'B' },
    { "output",                 1, 0, 'f' },
    { "only-aligned",           0, 0, 128 },
    { "debug-bam",              0, 0, 129 },
    { "broken-input",           0, 0, 130 },
    { "skip-duplicates",        0, 0, 131 },
    { "temp-dir",               1, 0, 132 },
    { "max-insert-size",        1, 0, 'a' },
    { "max-occurences",         1, 0, 'C' },
    { "max-occurences-se",      1, 0, 'D' },
    { "max-hits",               1, 0, 'h' },
    { "max-discordant-hits",    1, 0, 'H' },
    { "chimeric-rate",          1, 0, 'c' },
    { "disable-sw",             1, 0, 's' },
    { "disable-isize-estimate", 1, 0, 'A' },
    { "listen-port",            1, 0, 'p' },
    { 0,0,0,0 }
} ;

struct option workeropts[] = {
    { "num-threads",            1, 0, 't' },
    { "host",                   1, 0, 'h' },
    { "port",                   1, 0, 'p' },
    { "timeout",                1, 0, 'T' },
    { 0,0,0,0 }
} ;

// Yes, global variables.  Not very nice, but we'll only load one genome
// and passing it everywhere gets old after a while.

static bwt_t     *bwt[2]      = {0,0};
static bntseq_t  *bns         = 0;
static bntseq_t  *ntbns       = 0;
static ubyte_t   *pac         = 0;
static bwtint_t genome_length = 0;

static gap_opt_t *gap_opt     = 0;
static pe_opt_t  *pe_opt      = 0;
static int       only_aligned = 0;
static int       debug_bam    = 0;
static int       broken_input = 0;
static int       max_run_time = 90;
static int    skip_duplicates = 0;

static void      *zmq_context = 0;
static int       listen_port  = 0;

static isize_info_t null_ii = {0} ;
static khash_t(isize_infos) * volatile g_iinfos = 0;

static uint8_t revcom1[256] ;
static void init_revcom1()
{
    int i,j;
    for(i=0 ; i!=256; ++i)
    {
        j=0 ;
        if(i&0x1) j|=0x80 ;
        if(i&0x2) j|=0x40 ;
        if(i&0x4) j|=0x20 ;
        if(i&0x8) j|=0x10 ;
        if(i&0x10) j|=0x8 ;
        if(i&0x20) j|=0x4 ;
        if(i&0x40) j|=0x2 ;
        if(i&0x80) j|=0x1 ;
        revcom1[i] = j ;
    }
}

static volatile int s_interrupted = 0;
static void s_signal_handler (int signal_value)
{
    s_interrupted = 1;
}

static void s_catch_signals (void)
{
    struct sigaction action;
    action.sa_handler = s_signal_handler;
    action.sa_flags = 0;
    sigemptyset (&action.sa_mask);
    sigaction (SIGINT, &action, NULL);
    sigaction (SIGTERM, &action, NULL);
}

int zcheck_at( int rc, int testterm, int status, const char *label, int line, char *format, ... )
{
    va_list ap;
    if (rc != -1) return 0;
    if (testterm && (zmq_errno() == ETERM || zmq_errno() == EINTR)) return 1;
    fprintf( stderr, "[%s:%d] ", label, line );
    va_start(ap,format);
    vfprintf( stderr, format, ap );
    va_end(ap);
    fprintf( stderr, "%s%s\n", zmq_errno()?": ":"", zmq_strerror(zmq_errno()) );
    if (status) exit(status) ;
    return 0;
}

#define zcheck(rc, status, ...) (void)zcheck_at(rc,0,status,__FUNCTION__,__LINE__,__VA_ARGS__)
#define zterm(rc, status, ...)     if(zcheck_at(rc,1,status,__FUNCTION__,__LINE__,__VA_ARGS__))

// We take the old header, remove the lines we generate ourselves, and
// incorporate the rest verbatim.  More specifically, HD is added, PG is added
// and linked to the previous one using PG-PP, SQ and HD lines are removed.
int32_t bwa_print_header_text( char* buf, int32_t len, bntseq_t *bns,
                               const char* oldhdr, const char *pptag, const char* myid, int argc, char *argv[] )
{
    int i;
	int32_t c = snprintf( buf, len, "@HD\tVN:1.4\n@PG\tID:%s%s%s\tPN:bwa\tVN:%s%s",
                          myid, pptag?"\tPP:":"", pptag?pptag:"", PACKAGE_VERSION,
                          argc?"\tCL:":"" );
    for (i = 0; i < argc ; ++i) 
        c += snprintf( buf+c, len>c?len-c:0, "%s%c", argv[i], i==argc-1?'\n':' ' );
        
	for (i = 0; i < bns->n_seqs; ++i)
        c += snprintf( buf+c, len>c?len-c:0, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);

    while (*oldhdr) {
        // check for an interesting header: HD and SQ are uninteresting
        int interesting = 1 ;
        if( oldhdr[0] == '@' && oldhdr[1] && oldhdr[2] ) {
            if( oldhdr[1] == 'S' && oldhdr[2] == 'Q' ) interesting = 0 ;
            if( oldhdr[1] == 'H' && oldhdr[2] == 'D' ) interesting = 0 ;
        }

        // scan old header, copy only interesting lines
        while( *oldhdr && *oldhdr != '\n' ) {
            if( interesting ) {
                if( c < len ) buf[c] = *oldhdr ;
                ++c;
            }
            ++oldhdr;
        }
        if( interesting ) {
            if( c < len ) buf[c] = '\n' ;
            ++c;
        }
        if( *oldhdr ) ++oldhdr;
    }
    return c;
}

KHASH_SET_INIT_STR(words)

// Generate useful program id and pp-tag.  We scan the PG lines of the
// header looking for one that is present and not linked to (the first
// one is fine, there should be only one).  Then we try our program id
// ("bwa"), and if that's already present, we append numbers until it
// isn't anymore.
//
// Returns two strings in (*pp) and (*id), both of which must be
// free()d.
void find_pp_tag( const char* h, char** pp, char** id )
{
    const char *he, *k ;
    int r, n;
    khint_t iter;
    kh_words_t *present = kh_init_words() ;
    kh_words_t *linked = kh_init_words() ;

    *pp=0;
    while( *h ) {
        if( h[0] == '@' && h[1] == 'P' && h[2] == 'G' ) {
            while( *h && *h != '\n' ) {
                if( ( (h[0] == 'I' && h[1] == 'D') 
                            || (h[0] == 'P' && h[1] == 'P') )
                        && h[2] == ':' )
                {
                    h += 3 ;
                    he = h ;
                    while( *he && *he != '\n' && *he != '\t' ) ++he ;
                    k = strndup(h, he-h) ;
                    kh_put_words(h[-3]=='I' ? present : linked, k, &r) ;
                    h=he ;
                }
                while( *h && *h != '\n' && *h != '\t' ) ++h ;
                if( *h == '\t' ) ++h ;
            }
        }
        while( *h && *h != '\n' ) ++h ;
        if( *h ) ++h ;
    }

	for (iter = kh_begin(present); iter != kh_end(present); ++iter) {
		if (kh_exist(present, iter)) {
            khint_t j = kh_get(words, linked, kh_key(present, iter)) ;
            if( j == kh_end(linked) || !kh_exist(linked,j) ) {
                *pp = strdup( kh_key(present, iter) ) ;
                break ;
            }
        }
    }

    char myid[50] = "bwa" ;
    for( n = 1; 1 ; ++n ) {
        khint_t j = kh_get(words, present, myid) ;
        if( j == kh_end(present) || !kh_exist(present,j) ) {
            *id = strdup( myid ) ;
            break ;
        }
        snprintf( myid, 50, "bwa-%d", n ) ;
    }

	for (iter = kh_begin(linked); iter != kh_end(linked); ++iter) 
		if (kh_exist(linked, iter)) 
            free((void*)kh_key(linked,iter)) ;
	for (iter = kh_begin(present); iter != kh_end(present); ++iter) 
		if (kh_exist(present, iter)) 
            free((void*)kh_key(present,iter)) ;
    kh_destroy_words( linked ) ;
    kh_destroy_words( present ) ;
}

int bwa_print_bam_header( BGZF* output, bntseq_t *bns,
                          char* oldhdr, int argc, char *argv[] )
{
#define DO(x) if(r>=0) { int r1 = x; if(r1>=0) r+=r1; else r=r1; } else {}
	int i, r=0;
    char magic[4] = { 'B', 'A', 'M', 1 } ;
    char *pptag, *myid ;
    find_pp_tag( oldhdr, &pptag, &myid ) ;

    int32_t h_len = bwa_print_header_text( 0, 0, bns, oldhdr, pptag, myid, argc, argv );
    char *buf = malloc( h_len+1 ) ;
    bwa_print_header_text( buf, h_len+1, bns, oldhdr, pptag, myid, argc, argv );
    DO( bgzf_write( output, magic, 4 ) ) ;
    DO( bgzf_write( output, &h_len, 4 ) ) ;
    DO( bgzf_write( output, buf, h_len ) ) ;
    DO( bgzf_write( output, &bns->n_seqs, 4 ) ) ;

	for (i = 0; i < bns->n_seqs; ++i)
    {
        int32_t nlen = 1+strlen( bns->anns[i].name ) ;
        DO( bgzf_write( output, &nlen, 4 ) ) ;
        DO( bgzf_write( output, bns->anns[i].name, nlen ) ) ;
        DO( bgzf_write( output, &bns->anns[i].len, 4 ) ) ;
    }
    free(buf) ;
    free(pptag) ;
    free(myid) ;
    return r ;
}

int bwa_print_bam1( BGZF* output, bam1_t *b )
{
    uint32_t blocksize = sizeof(bam1_core_t) + b->data_len ;
    uint32_t y = (int)b->core.bin << 16 | (int)b->core.qual << 8 | (int)b->core.l_qname;
    uint32_t z = (int)b->core.flag << 16 | (int)b->core.n_cigar;

    int r = bgzf_write( output, &blocksize, 4);
	DO( bgzf_write( output, &b->core.tid, 4));
	DO( bgzf_write( output, &b->core.pos, 4));
	DO( bgzf_write( output, &y, 4 ));
	DO( bgzf_write( output, &z, 4 ));
	DO( bgzf_write( output, &b->core.l_qseq, 4));
	DO( bgzf_write( output, &b->core.mtid, 4));
	DO( bgzf_write( output, &b->core.mpos, 4));
	DO( bgzf_write( output, &b->core.isize, 4));
    DO( bgzf_write( output, b->data, b->data_len));
    return r;
}
#undef DO

static inline int bam_reg2bin(uint32_t beg, uint32_t end)
{
	--end;
	if (beg>>14 == end>>14) return 4681 + (beg>>14);
	if (beg>>17 == end>>17) return  585 + (beg>>17);
	if (beg>>20 == end>>20) return   73 + (beg>>20);
	if (beg>>23 == end>>23) return    9 + (beg>>23);
	if (beg>>26 == end>>26) return    1 + (beg>>26);
	return 0;
}

static void revcom_bam1( bam1_t *b )
{
    uint8_t *p,*q;
    
    // flip flag
    b->core.flag ^= SAM_FSR;

    // revcom sequence (incl. up to one nybble of padding)
    for(p=bam1_seq(b), q=bam1_seq(b) + (b->core.l_qseq+1) / 2 - 1; p<q ; ++p,--q)
    {
        uint8_t x = *q ; *q = revcom1[*p] ; *p = revcom1[x] ;
    }
    if(p==q) *p = revcom1[*q] ;

    // shift one nybble if necessary
    if(b->core.l_qseq & 1) {
        for(p=bam1_seq(b), q=bam1_seq(b) + (b->core.l_qseq+1) / 2 - 1; p<q ; ++p)
            *p = (p[0] & 0xf) << 4 | (p[1] & 0xf0) >> 4 ;
        *p = (p[0] & 0xf) << 4 ;
    }


    // reverse quality
    for(p=bam1_qual(b), q=bam1_qual(b) + b->core.l_qseq -1; p<q ; ++p,--q)
    {
        uint8_t x = *q ; *q = *p ; *p = x ;
    }
}

// Make room for at least n more bytes.
static void bam_realloc( int n, bam1_t *b )
{
    if( b->m_data - b->data_len < n ) {
        b->m_data = b->data_len+n ;
        kroundup32(b->m_data) ;
        b->data = realloc( b->data, b->m_data ) ;
    }
}

// add a tag with an int content.  Right now, encode as unsigned int, need to
// add the other tags if time permits.
static void bam_push_int( bam1_t *b, char u, char v, int x )
{
    bam_realloc( 7, b );
    b->data[b->data_len + 0] = u;
    b->data[b->data_len + 1] = v;
    b->data[b->data_len + 2] = 'i';
    *(uint32_t*)(b->data + b->data_len + 3) = x;
    b->data_len += 7;
}

static void bam_push_char( bam1_t *b, char u, char v, char c )
{
    bam_realloc( 4, b );
    b->data[b->data_len + 0] = u;
    b->data[b->data_len + 1] = v;
    b->data[b->data_len + 2] = 'A';
    b->data[b->data_len + 3] = c;
    b->data_len += 4;
}

static void bam_push_string( bam1_t *b, char u, char v, char *s )
{
    int len = strlen(s) ;
    bam_realloc( 4+len, b );
    b->data[b->data_len + 0] = u;
    b->data[b->data_len + 1] = v;
    b->data[b->data_len + 2] = 'Z';
    strcpy((char*)b->data + b->data_len + 3, s) ;
    b->data_len += 4+len;
}

static void bam_resize_cigar( bam1_t *b, int n_cigar )
{
    bam_realloc( 4*(n_cigar - b->core.n_cigar), b ) ;
    memmove( bam1_cigar(b) + n_cigar, bam1_cigar(b) + b->core.n_cigar,
            (b->data + b->data_len) - (uint8_t*)(bam1_cigar(b) + b->core.n_cigar) ) ;
    b->data_len += 4*(n_cigar - b->core.n_cigar) ;
    b->core.n_cigar = n_cigar ;
}


// The old record must be cleaned (remove the alignment), then the new
// alignment gets spliced in.  Alignment related tagged fields must also be
// overwritten or removed.  Specifically:
//
// - refID, pos, bin, mapq must be set
// - next_refID, next_pos, tlen must be set
// - cigar must be generated or cleared
// - seq/qual may need to be reverse-complemented
// - some flags must be cleared: proper pair, unmapped, mate unmapped,
//   reversed, mate reversed, secondary alignment
// - tags must be copied, but those we generate must be removed and
//   regenerated: XC, NM, XN, SM, AM, X0, X1, XM, XO, XG, XA, YQ, XT, MD

void bwa_update_bam1(bam1_t *out, const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{
    if (p->clip_len < p->full_len) bam_push_int( out, 'X', 'C', p->clip_len );
    if (p->max_entries && debug_bam) bam_push_int( out, 'Y', 'Q', p->max_entries );

	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
		int seqid, am = 0, nn, j;
		if (p->type == BWA_TYPE_NO_MATCH) {
			p->pos = mate->pos;
			p->strand = mate->strand;
			p->extra_flag |= SAM_FSU;
			j = 1;
		} else j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

        // reversed now but originally wasn't or not reversed now but was originally?
        // Note that this also inverts SAM_FSR if necessary,
		if (p->strand != ((out->core.flag & SAM_FSR) != 0) ) revcom_bam1( out ) ;
            
        out->core.flag &= ~( SAM_FPP | SAM_FSU | SAM_FMU | SAM_FSC | SAM_FMR ) ;
        out->core.flag |= p->extra_flag;

		// get seqid
		nn = bns_coor_pac2real(bns, p->pos, j, &seqid);
		if (p->type != BWA_TYPE_NO_MATCH && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
        {
			out->core.flag |= SAM_FSU;  // flag UNMAP and not PROPERLY PAIRED as this alignment
            out->core.flag &= ~SAM_FPP; // bridges two adjacent reference sequences
            p->mapQ = 0;
        }

		// update other flags in output record
        out->core.tid = seqid ; 
        out->core.pos = p->pos - bns->anns[seqid].offset ;
        out->core.bin = bam_reg2bin( p->pos - bns->anns[seqid].offset, pos_end(p) - bns->anns[seqid].offset ) ;
        out->core.qual = p->mapQ ;

		if (p->cigar) {
            bam_resize_cigar( out, p->n_cigar ) ;
			for (j = 0; j != p->n_cigar; ++j)
                bam1_cigar(out)[j] = __cigar_len(p->cigar[j]) << 4 | "\000\001\002\004"[__cigar_op(p->cigar[j])] ;

		} else if (p->type == BWA_TYPE_NO_MATCH) {
            bam_resize_cigar( out, 0 ) ;
        } else {
            bam_resize_cigar( out, 1 ) ;
            *bam1_cigar(out) = p->len << 4 ;
        }

		if (mate && mate->type != BWA_TYPE_NO_MATCH) {
			int m_seqid, m_j;
			am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
            // also, add up the nn values so self and mate get the same
            // XN field
			nn += bns_coor_pac2real(bns, mate->pos, mate->len, &m_seqid);

            m_j = pos_end(mate) - mate->pos; // m_j is the length of the reference in the alignment
            if( mate->pos + m_j - bns->anns[m_seqid].offset > bns->anns[m_seqid].len ) {
                out->core.flag |= SAM_FMU;  // flag MUNMAP and not PROPERLY PAIRED as the mate's alignment
                out->core.flag &= ~SAM_FPP; // bridges two adjacent reference sequences
            }
            if (mate->strand) out->core.flag |= SAM_FMR;
            out->core.mtid = m_seqid ;
            out->core.mpos = mate->pos - bns->anns[m_seqid].offset ;
			if (p->type == BWA_TYPE_NO_MATCH) out->core.isize = 0;
            else out->core.isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
		} else if (mate) {
            out->core.flag |= SAM_FMU;
            out->core.flag &= ~SAM_FPP;
            out->core.mtid = seqid ;
            out->core.mpos = p->pos - bns->anns[seqid].offset ;
            out->core.isize = 0;
        } else {
            out->core.mtid = -1 ;
            out->core.mpos = -1 ;
            out->core.isize = 0;
        }

		if (p->type != BWA_TYPE_NO_MATCH) {
			int i;
			// calculate XT tag
			char XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';

			// print tags
            bam_push_char( out, 'X', 'T', XT ) ;
            if( mode & BWA_MODE_COMPREAD ) bam_push_int( out, 'N', 'M', p->nm ) ;
            else bam_push_int( out, 'C', 'M', p->nm ) ;

            if (nn) bam_push_int( out, 'X', 'N', nn ) ;

            if (mate) {
                bam_push_int( out, 'S', 'M', p->seQ ) ;
                bam_push_int( out, 'A', 'M', am ) ;
            }

			if (p->type != BWA_TYPE_MATESW) { // X0 and X1 are not available for this type of alignment
                bam_push_int( out, 'X', '0', p->c1 ) ;
				if (p->c1 <= max_top2) 
                    bam_push_int( out, 'X', '1', p->c2 ) ;
			}
            bam_push_int( out, 'X', 'M', p->n_mm ) ;
            bam_push_int( out, 'X', 'O', p->n_gapo ) ;
            bam_push_int( out, 'X', 'G', p->n_gapo+p->n_gape ) ;

            if (p->md) bam_push_string( out, 'M', 'D', p->md ) ;
            
			// print multiple hits
			if (p->n_multi) for(;;) {
                char* outp = (char*)out->data + out->data_len ;
                char* pe   = (char*)out->data + out->m_data ;

#define mysnprintf(...) outp+=snprintf(outp,pe>outp?pe-outp:0,__VA_ARGS__)

                mysnprintf("XAZ");
				for (i = 0; i < p->n_multi; ++i) {
					bwt_multi1_t *q = p->multi + i;
					int k;
					j = pos_end_multi(q, p->len) - q->pos;
					bns_coor_pac2real(bns, q->pos, j, &seqid);
                    mysnprintf("%s,%c%d,", bns->anns[seqid].name, q->strand? '-' : '+',
						   (int)(q->pos - bns->anns[seqid].offset + 1));
					if (q->cigar) {
						for (k = 0; k < q->n_cigar; ++k)
							mysnprintf("%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
					} else mysnprintf("%dM", p->len);
					mysnprintf(",%d;", q->gap + q->mm);
				}
                if(outp<pe) *outp=0;
                outp++;

#undef mysnprintf
                if(outp<=pe) {
                    out->data_len = outp - (char*)out->data;
                    break;
                }
                out->m_data = outp - (char*)out->data;
                kroundup32(out->m_data) ;
                out->data = realloc( out->data, out->m_data ) ;
			}
		}
	} else { // this read has no match
        out->core.tid = -1 ; 
        out->core.pos = -1 ; 
        out->core.bin = 0 ;
        out->core.qual = 0 ;
        out->core.mtid = -1 ; 
        out->core.mpos = -1 ; 
        out->core.isize = 0 ;

        out->core.flag &= ~( SAM_FPP | SAM_FMU | SAM_FSC ) ;
        out->core.flag |= SAM_FSU ;
		if (mate && mate->type == BWA_TYPE_NO_MATCH) out->core.flag |= SAM_FMU;

        bam_resize_cigar( out, 0 ) ;

        // if the mate has an XN tag, we need to reproduce it here
		if (mate && mate->type != BWA_TYPE_NO_MATCH) {
			int m_seqid, nn ;
			nn = bns_coor_pac2real(bns, mate->pos, mate->len, &m_seqid);
            if (nn) bam_push_int( out, 'X', 'N', nn ) ;
        }
	}
}

static int unique(bam_pair_t *p)
{
    if( skip_duplicates ) {
        switch (p->kind) {
            case eof_marker:  return 0;
            case singleton:   return !(p->bam_rec[0].core.flag & SAM_FDP) ;
            case proper_pair: return !(p->bam_rec[0].core.flag & SAM_FDP)
                                  && !(p->bam_rec[1].core.flag & SAM_FDP) ;
        }
    }
    return 1;
}

static void aln_singleton( bam_pair_t *raw )
{
    // If we already have alignments (could come from sai file or
    // because the firs phase has already run), we skip the alignment.
    // (This point shouldn't actually be reached in that case.)
    if(raw->phase == pristine) {
        if(unique(raw)) {
            bam1_to_seq(&raw->bam_rec[0], &raw->bwa_seq[0], 1, gap_opt->trim_qual);
            bwa_cal_sa_reg_gap(bwt, 1, &raw->bwa_seq[0], gap_opt);
        }
        raw->phase = aligned ;
    }
}

static void posn_singleton( bam_pair_t *raw )
{
    int j ;
    if( raw->phase == aligned ) {
        if(unique(raw)) {
            // from bwa_sai2sam_se_core
            bwa_seq_t *p = &raw->bwa_seq[0] ;
            bwa_aln2seq_core(p->n_aln, p->aln, p, 1, pe_opt->max_occ_se );

            // from bwa_cal_pac_pos
            bwa_cal_pac_pos_core(bwt[0], bwt[1], p, gap_opt->max_diff, gap_opt->fnr);
            for (j = 0; j < p->n_multi; ++j) {
                bwt_multi1_t *q = raw->bwa_seq[0].multi + j;
                if (q->strand) q->pos = bwt_sa(bwt[0], q->pos);
                else q->pos = bwt[1]->seq_len - (bwt_sa(bwt[1], q->pos) + p->len);
            }
        }
        raw->phase = positioned ;
    }
}

static void finish_singleton( bam_pair_t *raw )
{
    if( raw->phase == positioned ) {
        if(unique(raw)) {
            bwa_seq_t *p = &raw->bwa_seq[0] ;
            if( !p->seq ) bam1_to_seq(&raw->bam_rec[0], p, 1, gap_opt->trim_qual);
            bwa_refine_gapped(bns, 1, p, pac, ntbns);
            bwa_update_bam1(&raw->bam_rec[0], bns, p, 0, gap_opt->mode, gap_opt->max_top2);
            bwa_free_read_seq1(p);
        }
        raw->phase = finished ;
    }
}

// First stage of paired end alignment: this computes the arrays 'aln1'
// and 'aln2', the two positions, strands and single map qualities.
// pos, strand, mapQ are stored directly into the bam record.
static void aln_pair( bam_pair_t *raw )
{
    int j ;
    if( raw->phase == pristine ) {
        if(unique(raw)) {
            for( j = 0 ; j != 2 ; ++j ) {
                bam1_to_seq(&raw->bam_rec[j], &raw->bwa_seq[j], 1, gap_opt->trim_qual);

                // BEGIN from bwa_sai2sam_pe_core
                // BEGIN from bwa_cal_pac_pos_pe
                // note that seqs is not an array of two arrays anymore!!
                // also note that we don't track this
                // cnt_chg... whatever it was good
                // for.
                // cnt_chg = bwa_cal_pac_pos_pe(seqs, &ii, &last_ii);

                bwa_cal_sa_reg_gap(bwt, 1, &raw->bwa_seq[j], gap_opt);
            }
        }
        raw->phase = aligned ;
    }
}

static void posn_pair( bam_pair_t *raw )
{
    int j ;
    if( raw->phase == aligned ) {
        if(unique(raw)) {
            for( j = 0 ; j != 2 ; ++j ) {
                // SE part.  Almost, but not quite the
                // same as for singletons.  (*sigh*)
                raw->bwa_seq[j].n_multi = 0 ;

                // generate SE alignment and mapping quality
                bwa_aln2seq(raw->bwa_seq[j].n_aln, raw->bwa_seq[j].aln, &raw->bwa_seq[j]);

                // Computes pos, seQ, mapQ.  need to store only those to avoid
                // repeated computation!
                bwa_cal_pac_pos_core(bwt[0], bwt[1], &raw->bwa_seq[j], gap_opt->max_diff, gap_opt->fnr);
            }
        }
        raw->phase = positioned ;
    }
}

static void finish_pair(
        bam_pair_t *raw, khash_t(isize_infos) *iinfos,
        uint64_t n_tot[2], uint64_t n_mapped[2], kh_64_t *my_hash )
{
    if( raw->phase != positioned ) return ;
    if( unique(raw) ) {
        bwa_seq_t *p[2] = { raw->bwa_seq, raw->bwa_seq+1 } ;
        pe_data_t d;
        int j;

        khiter_t it = kh_get(isize_infos, iinfos, bam_get_rg(raw->bam_rec)) ;
        isize_info_t *ii = it == kh_end(iinfos) ? &null_ii : &kh_val(iinfos,it) ;

        memset(&d, 0, sizeof(pe_data_t));
        for (j = 0; j < 2; ++j) 
        {
            if( !raw->bwa_seq[j].seq ) bam1_to_seq(&raw->bam_rec[j], &raw->bwa_seq[j], 1, gap_opt->trim_qual);
            d.aln[j].a = p[j]->aln;     // cheapting, but we do not really need to copy here
            d.aln[j].n = p[j]->n_aln; 
        }

        if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
                && (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT))
        { // only when both ends mapped
            uint64_t x;
            int j, k;
            long long n_occ[2];
            for (j = 0; j < 2; ++j) {
                n_occ[j] = 0;
                for (k = 0; k < d.aln[j].n; ++k)
                    n_occ[j] += d.aln[j].a[k].l - d.aln[j].a[k].k + 1;
            }
            if (n_occ[0] <= pe_opt->max_occ && n_occ[1] <= pe_opt->max_occ) {
                d.arr.n = 0;
                for (j = 0; j < 2; ++j) {
                    for (k = 0; k < d.aln[j].n; ++k) {
                        bwt_aln1_t *r = d.aln[j].a + k;
                        bwtint_t l;
                        if (r->l - r->k + 1 >= MIN_HASH_WIDTH) { // then check hash table
                            uint64_t key = (uint64_t)r->k<<32 | r->l;
                            int ret;
                            khint_t iter = kh_put(64, my_hash, key, &ret);
                            if (ret) { // not in the hash table; ret must equal 1 as we never remove elements
                                poslist_t *z = &kh_val(my_hash, iter);
                                z->n = r->l - r->k + 1;
                                z->a = (bwtint_t*)malloc(sizeof(bwtint_t) * z->n);
                                for (l = r->k; l <= r->l; ++l)
                                    z->a[l - r->k] = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);
                            }
                            for (l = 0; l < kh_val(my_hash, iter).n; ++l) {
                                x = kh_val(my_hash, iter).a[l];
                                x = x<<32 | k<<1 | j;
                                kv_push(uint64_t, d.arr, x);
                            }
                        } else { // then calculate on the fly
                            for (l = r->k; l <= r->l; ++l) {
                                x = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);
                                x = x<<32 | k<<1 | j;
                                kv_push(uint64_t, d.arr, x);
                            }
                        }
                    }
                }
                pairing(p, &d, pe_opt, gap_opt->s_mm, ii);
            }
        }

        if (pe_opt->N_multi || pe_opt->n_multi) {
            for (j = 0; j < 2; ++j) {
                if (p[j]->type != BWA_TYPE_NO_MATCH) {
                    int k;
                    if (!(p[j]->extra_flag&SAM_FPP) && p[1-j]->type != BWA_TYPE_NO_MATCH) {
                        bwa_aln2seq_core(
                                d.aln[j].n, d.aln[j].a, p[j], 0,
                                p[j]->c1+p[j]->c2-1 > pe_opt->N_multi ?
                                pe_opt->n_multi : pe_opt->N_multi );
                    }
                    else
                        bwa_aln2seq_core(d.aln[j].n, d.aln[j].a, p[j], 0, pe_opt->n_multi);
                    for (k = 0; k < p[j]->n_multi; ++k) {
                        bwt_multi1_t *q = p[j]->multi + k;
                        q->pos = q->strand? bwt_sa(bwt[0], q->pos) : bwt[1]->seq_len - (bwt_sa(bwt[1], q->pos) + p[j]->len);
                    }
                }
            }
        }

        kv_destroy(d.arr);
        kv_destroy(d.pos[0]); kv_destroy(d.pos[1]);

        // END from bwa_cal_pac_pos_pe
        bwa_paired_sw1(bns, pac, p, pe_opt, ii, n_tot, n_mapped);

        // END from bwa_sai2sam_pe_core
        bwa_refine_gapped(bns, 1, p[0], pac, ntbns);
        bwa_refine_gapped(bns, 1, p[1], pac, ntbns);

        // For PE reads, stock BWA would have concatenated their
        // barcodes.  We don't, for once because we don't identify
        // barcodes, but also because the idea feels wrong.
        bwa_update_bam1( &raw->bam_rec[0], bns, &raw->bwa_seq[0], &raw->bwa_seq[1], gap_opt->mode, gap_opt->max_top2);
        bwa_update_bam1( &raw->bam_rec[1], bns, &raw->bwa_seq[1], &raw->bwa_seq[0], gap_opt->mode, gap_opt->max_top2);
        bwa_free_read_seq1(&raw->bwa_seq[1]);
        bwa_free_read_seq1(&raw->bwa_seq[0]);
    }
    raw->phase = finished ;
}

static inline double tdiff( struct timeval* tv1, struct timeval *tv2 )
{
    return (double)(tv2->tv_sec-tv1->tv_sec) + 0.000001 * (double)(tv2->tv_usec-tv1->tv_usec) ;
}

/* XXX  This is creating difficulties in a cluster environment.
 *
 * We want the genome index to live in shared memory, in case the
 * cluster scheduler starts up more than one instance on a single node.
 * So we need mmap().  We also want to have the genome index on a
 * network file system, to avoid silly duplication.
 *
 * The genome index on NFS already causes substantial traffic, combined
 * with mmap() it seems to get a lot worse.  The way out could be a
 * "shared memory object", which is filled explicitly and then mmapped,
 * maybe combined with a multicast protocol for distribution.  But then
 * we need to dispose of the shared memory object after some timeout.
 *
 * So, here's the clean idea:  Operate on a named shared memory object
 * ("index") and two named semaphores ("ready", "alive").  To initialize
 * the genome index, grab "ready", mmap "index" and put back "ready".
 * At regular intervals, put "alive", when done, unmap "index" (but do
 * not unlink anything).
 *
 * If the three don't exist, create them, and fork a process to
 * initialize them, then proceed as above.  The child loads the genome
 * index, then puts "ready".  It then does a timed wait on "alive".  If
 * it runs into a timeout, it unlinks everything and exits.
 *
 * To be implemented when I feel like it.
 */
void init_genome_index( const char* prefix, int touch )
{
    struct timeval tv, tv1 ;
    gettimeofday( &tv, 0 ) ;

    char *str = (char*)calloc(strlen(prefix) + 10, 1);
    bwase_initialize();
    init_revcom1();

    fprintf(stderr, "[init_genome_index] loading index... ");
    strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str,touch);
    strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str,touch);
    strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0],touch);
    strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1],touch);
    free(str);

    // if (!(gap_opt.mode & BWA_MODE_COMPREAD)) {  // in color space; initialize ntpac
    //	pe_opt->type = BWA_PET_SOLID;
    //	ntbns = bwa_open_nt(prefix);
    // }	
    pac = bwt_restore_pac( bns,touch ) ;
    gettimeofday( &tv1, 0 ) ;
    fprintf(stderr, "%.2f sec\n", tdiff(&tv, &tv1));
    genome_length = bwt[0]->seq_len ;
}

void init_genome_params( const char* prefix ) 
{
    bwtint_t foo[5] ;
    char *fn = alloca(strlen(prefix) + 10);
    strcpy(fn, prefix); strcat(fn, ".bwt");
    FILE *fp = xopen(fn, "rb");
	err_fread(foo, sizeof(bwtint_t), 5, fp);
    genome_length = foo[4] ;
    fclose(fp) ;
}


void pair_aln(bam_pair_t *p)
{
    switch (p->kind) {
        case eof_marker:  break ;
        case singleton:   aln_singleton(p) ; break ;
        case proper_pair: aln_pair(p) ; break ;
    }
}

void pair_posn(bam_pair_t *p)
{
    switch (p->kind) {
        case eof_marker: break ;
        case singleton:   posn_singleton(p) ; break ;
        case proper_pair: posn_pair(p) ; break ;
    }
}

void pair_finish(
        bam_pair_t *p, khash_t(isize_infos) *iinfos,
        uint64_t n_tot[2], uint64_t n_mapped[2], kh_64_t *my_hash )
{
    switch (p->kind) {
        case eof_marker: break ;
        case singleton:   finish_singleton(p) ; break ;
        case proper_pair: finish_pair(p,iinfos,n_tot,n_mapped,my_hash) ; break ;
    }
}

void pair_print_bam(BGZF *output, bam_pair_t *p)
{
    int i ;
    // If "only aligned" is requested and at least one mate is unmapped,
    // we skip the output.  (Conserves time and space where we don't
    // need to pass through the unaligned sequences.)
    if( only_aligned )
        for( i = 0 ; i != p->kind ; ++i )
            if( p->bam_rec[i].core.flag & SAM_FSU )
                return ;

    for( i = 0 ; i != p->kind ; ++i )
        bwa_print_bam1(output, &p->bam_rec[i]);
}

inline void put_int( unsigned char **p, int x )
{
    (*p)[0] = x >>  0 & 0xff ;
    (*p)[1] = x >>  8 & 0xff ;
    (*p)[2] = x >> 16 & 0xff ;
    (*p)[3] = x >> 24 & 0xff ;
    *p += 4 ;
}

inline void put_block( unsigned char **p, void *q, size_t s )
{
    memcpy( *p, q, s ) ;
    *p += s ;
}

/** 
 * Serializes a bam_pair_t, no matter what phase it is in.
 * We only serialize the valid fields, depending on phase.
 */
void msg_init_from_pair(zmq_msg_t *m, bam_pair_t *p)
{
    // recno, kind, phase, 0-2 bam_recs w/ data_len
    int rc, i, len = 6 ;
    for( i = 0 ; i != p->kind ; ++i ) {
        len += sizeof( bam1_core_t ) + 4 + p->bam_rec[i].data_len ;
        switch( p->phase ) {
            case pristine: break ;
            case finished: break ;
            case positioned: len += 38 + p->bwa_seq[i].n_multi * sizeof(bwt_multi1_t) ;
                             // fallthrough!
            case aligned: len += 8 + p->bwa_seq[i].n_aln * sizeof(bwt_aln1_t) ;
                          break ;
        }
    }

    rc = zmq_msg_init_size(m, len) ;
    xassert (rc == 0, "fail to init message");

    unsigned char *q = zmq_msg_data(m);
    put_int( &q, p->recno ) ;
    *q++ = p->kind ;
    *q++ = p->phase ;

    for( i = 0 ; i != p->kind ; ++i ) {
        put_block( &q, &p->bam_rec[i].core, sizeof( bam1_core_t ) ) ;
        put_int(   &q,  p->bam_rec[i].data_len ) ;
        put_block( &q,  p->bam_rec[i].data, p->bam_rec[i].data_len ) ;

        switch( p->phase ) {
            case pristine:   break ;
            case finished:   break ;
            case positioned: *q++ = p->bwa_seq[i].strand << 4 | p->bwa_seq[i].type ;
                             *q++ = p->bwa_seq[i].n_mm ;
                             *q++ = p->bwa_seq[i].n_gapo ;
                             *q++ = p->bwa_seq[i].n_gape ;
                             *q++ = p->bwa_seq[i].seQ ;
                             *q++ = p->bwa_seq[i].mapQ ;
                             put_int( &q, p->bwa_seq[i].len ) ;
                             put_int( &q, p->bwa_seq[i].clip_len ) ;
                             put_int( &q, p->bwa_seq[i].score ) ;
                             put_int( &q, p->bwa_seq[i].sa ) ;
                             put_int( &q, p->bwa_seq[i].c1 ) ;
                             put_int( &q, p->bwa_seq[i].c2 ) ;
                             put_int( &q, p->bwa_seq[i].pos ) ;
                             put_int( &q, p->bwa_seq[i].n_multi ) ;
                             put_block( &q, p->bwa_seq[i].multi, p->bwa_seq[i].n_multi * sizeof(bwt_multi1_t) ) ;
                             // fallthrough!
            case aligned:    put_int(   &q, p->bwa_seq[i].max_entries ) ;
                             put_int(   &q, p->bwa_seq[i].n_aln ) ;
                             put_block( &q, p->bwa_seq[i].aln, p->bwa_seq[i].n_aln * sizeof(bwt_aln1_t) ) ;
                             break ;
        }
    }
    xassert( q == zmq_msg_data(m) + len, "logic error in msg_init_from_pair" ) ;
}

inline int get_int( unsigned char **p )
{
    *p += 4 ;
    return (int)(*p)[-1] << 24 |
           (int)(*p)[-2] << 16 |
           (int)(*p)[-3] <<  8 |
           (int)(*p)[-4] <<  0 ;
}

inline void get_block( unsigned char **p, void *q, size_t s )
{
    memcpy( q, *p, s ) ;
    *p += s ;
}


void pair_init_from_msg(bam_pair_t *p, zmq_msg_t *m)
{
    int i ;
    unsigned char *q = zmq_msg_data(m);

    memset( p, 0, sizeof( bam_pair_t ) ) ;
    p->recno = get_int( &q ) ;
    p->kind  = *q++ ;
    p->phase = *q++ ;

    for( i = 0 ; i != p->kind ; ++i ) {
        get_block( &q, &p->bam_rec[i].core, sizeof( bam1_core_t ) ) ;
        p->bam_rec[i].data_len = p->bam_rec[i].m_data = get_int( &q ) ;
        p->bam_rec[i].data = malloc( p->bam_rec[i].data_len ) ;
        get_block( &q,  p->bam_rec[i].data, p->bam_rec[i].data_len ) ;

        switch( p->phase ) {
            case pristine:   break ;
            case finished:   break ;
            case positioned: p->bwa_seq[i].strand = *q >> 4 ;
                             p->bwa_seq[i].type   = *q++ & 0xff ;
                             p->bwa_seq[i].n_mm   = *q++ ;
                             p->bwa_seq[i].n_gapo = *q++ ;
                             p->bwa_seq[i].n_gape = *q++ ;
                             p->bwa_seq[i].seQ    = *q++ ;
                             p->bwa_seq[i].mapQ   = *q++ ;
                             p->bwa_seq[i].len    = get_int( &q ) ;
                             p->bwa_seq[i].clip_len=get_int( &q ) ;
                             p->bwa_seq[i].score  = get_int( &q ) ;
                             p->bwa_seq[i].sa     = get_int( &q ) ;
                             p->bwa_seq[i].c1     = get_int( &q ) ;
                             p->bwa_seq[i].c2     = get_int( &q ) ;
                             p->bwa_seq[i].pos    = get_int( &q ) ;
                             p->bwa_seq[i].n_multi= get_int( &q ) ;
                             p->bwa_seq[i].multi  = malloc( p->bwa_seq[i].n_multi * sizeof(bwt_multi1_t) ) ;
                             get_block( &q, p->bwa_seq[i].multi, p->bwa_seq[i].n_multi * sizeof(bwt_multi1_t) ) ;
                             // fallthrough!
            case aligned:    p->bwa_seq[i].max_entries = get_int( &q ) ;
                             p->bwa_seq[i].n_aln       = get_int( &q ) ;
                             p->bwa_seq[i].aln         = malloc( p->bwa_seq[i].n_aln * sizeof(bwt_aln1_t) ) ;
                             get_block( &q,  p->bwa_seq[i].aln, p->bwa_seq[i].n_aln * sizeof(bwt_aln1_t) ) ;
                             break ;
        }
    }
    if( q != (unsigned char*)zmq_msg_data(m) + zmq_msg_size(m) ) {
        int i ;
        q = (unsigned char*)zmq_msg_data(m) ;
        fprintf( stderr, "error decoding message: %d bytes, phase %d, kind %d, recno %d, hexdump follows",
                (int)zmq_msg_size(m), p->phase, p->kind, p->recno ) ;
        for( i = 0 ; i != zmq_msg_size(m) ; ++i )
        {
            if( !(i & 31) ) fputc( '\n', stderr ) ;
            fprintf( stderr, " %02x", (int)(q[i]) ) ;
        }
        fputc( '\n', stderr ) ;
        exit(1) ;
    }
}

void pair_print_custom( gzFile f, bam_pair_t *p )
{
    zmq_msg_t m ;
    msg_init_from_pair( &m, p ) ;
    uint32_t len = zmq_msg_size( &m ) ;
    int r, l1 = gzwrite( f, &len, sizeof(len) ) ;
    int l2 = gzwrite( f, zmq_msg_data( &m ), len ) ;
    if( l1 != sizeof(len) || l2 != len ) {
        fprintf(stderr, "[%s] error writing temporary file (%s,%s)\n",
                __FUNCTION__, gzerror(f,&r), strerror(errno) ) ;
        exit(1);
    }
    zmq_msg_close(&m) ;
}

int read_pair_custom( gzFile f, bam_pair_t *p )
{
    uint32_t len ;
    int r ;
    if( sizeof(len) == gzread( f, &len, sizeof(len) ) ) 
    {
        zmq_msg_t m ;
        zmq_msg_init_size( &m, len ) ;
        int l1 = gzread( f, zmq_msg_data(&m), len ) ;
        if( l1 != len ) goto hell ;
        pair_init_from_msg(p, &m);
        zmq_msg_close(&m);
        return p->kind;
    }
    else if( gzeof(f) )
    {
        memset(p, 0, sizeof(bam_pair_t)) ;
        return 0;
    }
hell:
    fprintf(stderr, "[%s] error reading temporary file (%s,%s)\n",
            __FUNCTION__, gzerror(f,&r), strerror(errno) ) ;
    exit(1);
}

// This is the simple sequential loop:  We run it if exactly one thread
// and no listening port is requested.  Doesn't require 0MQ.  As opposed
// to "classic" BWA, we stop that blockwise nonsense and get exactly one
// singleton or pair per iteration.
void sequential_loop_pass1( bwa_seqio_t *ks, gzFile temporary, khash_t(isize_infos) *iinfos )
{
    struct timeval tv0, tv1 ;
    gettimeofday( &tv0, 0 ) ;

    int tot_seqs = 0, rc=0 ;
    bam_pair_t raw ;
    for(;;) {
        rc=read_bam_pair(ks, &raw, broken_input);
        if(rc<0) {
            fprintf(stderr, "[%s] error reading input BAM%s\n",
                    __FUNCTION__, rc==-2 ? " (lone mate)" : "" ) ;
            exit(1);
        }
        if(rc==0) break ;

        tot_seqs += raw.kind ;
        if(tot_seqs%0x20000 < raw.kind) {
            gettimeofday( &tv1, 0 ) ;
            fprintf(stderr, "[%s] %d sequences processed in %.2f sec\n",
                    __FUNCTION__, tot_seqs, tdiff( &tv0, &tv1 ));
        }
        pair_aln(&raw);
        pair_posn(&raw);
        if(unique(&raw)) improve_isize_est(iinfos, &raw,pe_opt->ap_prior,genome_length);

        pair_print_custom(temporary,&raw);
        bam_destroy_pair(&raw);
    }
    gettimeofday( &tv1, 0 ) ;
    fprintf(stderr, "[%s] %d sequences processed in %.2f sec\n",
            __FUNCTION__, tot_seqs, tdiff( &tv0, &tv1 ));
    fprintf(stderr, "[%s] finished cleanly.\n", __FUNCTION__);
}

void sequential_loop_pass2( gzFile temporary, BGZF* output, khash_t(isize_infos) *iinfos )
{
    uint64_t n_tot[2]={0,0}, n_mapped[2]={0,0};
    struct timeval tv0, tv1 ;
    gettimeofday( &tv0, 0 ) ;

    int tot_seqs = 0, rc=0 ;
    khiter_t iter ;
    kh_64_t *my_hash ;
	my_hash = kh_init(64);
    bam_pair_t raw ;
    for(;;) {
        rc=read_pair_custom(temporary, &raw);
        if(rc<0) {
            fprintf(stderr, "[%s] error reading intermediate file\n", __FUNCTION__);
            exit(1);
        }
        if(rc==0) break ;

        tot_seqs += raw.kind ;
        if(tot_seqs%0x100000 < raw.kind) {
            gettimeofday( &tv1, 0 ) ;
            fprintf(stderr, "[%s] %d sequences processed in %.2f sec\n",
                    __FUNCTION__, tot_seqs, tdiff( &tv0, &tv1 ));
        }
        pair_finish(&raw, iinfos, n_tot, n_mapped, my_hash);
        pair_print_bam(output,&raw);
        bam_destroy_pair(&raw);
    }
    gettimeofday( &tv1, 0 ) ;
    fprintf(stderr, "[%s] %d sequences processed in %.2f sec\n"
                    "[%s] finished cleanly, shutting down.\n"
	                "[bwa_paired_sw] %lld out of %lld Q%d singletons are mated.\n"
	                "[bwa_paired_sw] %lld out of %lld Q%d discordant pairs are fixed.\n",
                    __FUNCTION__, tot_seqs, tdiff( &tv0, &tv1 ), __FUNCTION__,
                    (long long)n_mapped[1], (long long)n_tot[1], SW_MIN_MAPQ,
                    (long long)n_mapped[0], (long long)n_tot[0], SW_MIN_MAPQ );
	for (iter = kh_begin(my_hash); iter != kh_end(my_hash); ++iter)
		if (kh_exist(my_hash, iter))
            free(kh_val(my_hash, iter).a);
	kh_destroy(64, my_hash);
}

void set_sockopts( void *socket )
{
    zcheck( zmq_setsockopt( socket, ZMQ_SNDHWM, &the_hwm, sizeof(int) ),
            1, "zmq_setsockopt failed" ) ;
    zcheck( zmq_setsockopt( socket, ZMQ_RCVHWM, &the_hwm, sizeof(int) ),
            1, "zmq_setsockopt failed" ) ;
}
void zmq_xclose( void *socket )
{
    if( socket ) {
        int the_linger = 0;
        zmq_setsockopt( socket, ZMQ_LINGER, &the_linger, sizeof(int) ) ;
        zmq_close(socket);
    }
}


void *run_config_service( void *prefix )
{
    char buf[100] ;
    void *socket = zmq_socket( zmq_context, ZMQ_REP ) ;
    xassert(socket, "fail to create socket" );
    set_sockopts( socket );
    snprintf( buf, sizeof(buf), "tcp://*:%d", listen_port ) ;
    zcheck( zmq_bind( socket, buf ), 1, "zmq_bind to %s failed", buf );

    fprintf( stderr, "[run_config_service] Listening on port %d.\n", listen_port ) ;
    while (!s_interrupted) {
        zmq_msg_t m ;
        zmq_msg_init(&m) ;

        zterm( zmq_msg_recv(&m, socket, 0), 1, "zmq_recv failed" ) break;
        if( zmq_msg_size(&m) > 0 ) {
            char key = *(char*)zmq_msg_data(&m) ;
            if( key == 0 ) {
                zmq_msg_t m2 ;
                fprintf( stderr, "[run_config_service] Received hello: %.*s.\n", 
                        (int)zmq_msg_size(&m)-1, (char*)zmq_msg_data(&m)+1 ) ;

                zmq_msg_init_size(&m2, sizeof(gap_opt_t) + sizeof(pe_opt_t) + strlen(prefix)) ;
                memcpy(zmq_msg_data(&m2),                                       gap_opt, sizeof(gap_opt_t));
                memcpy(zmq_msg_data(&m2) + sizeof(gap_opt_t),                    pe_opt, sizeof(pe_opt_t));
                memcpy(zmq_msg_data(&m2) + sizeof(gap_opt_t) + sizeof(pe_opt_t), prefix, strlen(prefix) );
                zterm( zmq_msg_send(&m2, socket, 0), 1, "zmq_send failed" ) break;
                zmq_msg_close(&m2) ;
            }
            else if( key == 1 ) {
                zmq_msg_t m2 ;
                fprintf( stderr, "[run_config_service] Received update request: %.*s.\n", 
                        (int)zmq_msg_size(&m)-1, (char*)zmq_msg_data(&m)+1 ) ;

                zmq_msg_init_size(&m2, g_iinfos ? iinfo_encoded_size(g_iinfos) : 0 ) ;
                if( g_iinfos ) {
                    char *q = encode_iinfo( zmq_msg_data(&m2), g_iinfos ) ;
                    xassert( q == zmq_msg_data(&m2) + zmq_msg_size(&m2), "encoding of config failed" ) ;
                }
                zterm( zmq_msg_send(&m2, socket, 0), 1, "zmq_send failed" ) break;
                zmq_msg_close(&m2) ;
            }

        }
        zmq_msg_close(&m) ;
    }
    zmq_xclose( socket ) ;
    return 0;
}

    
struct reader_thread_args {
    void *socket ;
    bwa_seqio_t *ks;
    gzFile gzf ;
} ;


void *run_reader_thread( void* vargs )
{
    struct reader_thread_args *args = vargs ;

    int recno = 0, rc=0;
    bam_pair_t raw;
    zmq_msg_t msg;
    bam_init_pair(&raw);
    while (!s_interrupted) {
        if( args->ks ) rc=read_bam_pair(args->ks, &raw, broken_input);
        else rc=read_pair_custom(args->gzf, &raw);

        if(!rc) break ;
        if(rc<0) {
            fprintf( stderr, "[run_reader_thread] error reading input BAM %s\n",
                    rc==-2 ? "(lone mate)" : "" ) ;
            exit(1);
        }

        raw.recno = recno++; 
        if(!(recno%chunksize)) fprintf( stderr, "[run_reader_thread] %d records read so far.\n", recno );
        msg_init_from_pair( &msg, &raw );
        bam_destroy_pair(&raw);
        zterm( zmq_msg_send(&msg, args->socket, 0), 1, "zmq_send failed" ) break;
        zmq_msg_close(&msg);
    }

    if(!s_interrupted) {
        bam_init_pair(&raw);
        raw.kind = eof_marker;
        raw.recno = recno++;
        msg_init_from_pair( &msg, &raw );
        bam_destroy_pair(&raw);
        zterm( zmq_msg_send(&msg, args->socket, 0), 0, "zmq_send failed" );
        fprintf( stderr, "[run_reader_thread] finished, %d records in total.\n", recno );
        zmq_msg_close(&msg);
    }
    zmq_close(args->socket);
    return 0;
}
 
void *run_output_thread( BGZF *ks, gzFile output, void* socket, khash_t(isize_infos) *iinfos )
{
    struct timeval tv0, tv1 ;
    int lastrn = 0 ;
    double rate = -1 ;

    gettimeofday( &tv0, 0 ) ;
    while (!s_interrupted) 
    {
        zmq_msg_t msg;
        bam_pair_t pair;
        zmq_msg_init(&msg);
        zterm( zmq_msg_recv(&msg, socket, 0), 1, "zmq_recv failed" ) break;
        pair_init_from_msg(&pair, &msg);
        zmq_msg_close(&msg);
        if(!pair.kind) break ;
        if(!pair.kind) fprintf( stderr, "[run_output_thread] %d records received in total.\n", pair.recno );
        else if(!(pair.recno%0x100)) {
            gettimeofday( &tv1, 0 );
            double sec = tdiff( &tv0, &tv1 ) ;
            if( sec >= 10 ) {
                if( rate < 0 ) rate = (pair.recno-lastrn) / (1000*sec) ;
                else rate = ((pair.recno-lastrn) / (1000*sec) + 15*rate) * 0.0625 ;
                fprintf( stderr, "[run_output_thread] %d records received in %0.2fs, rate = %.1f kHz.\n",
                        pair.recno-lastrn, sec, rate ) ;
                lastrn=pair.recno;
                memmove( &tv0, &tv1, sizeof(struct timeval) ) ;
            }
        }
        if(unique(&pair)) improve_isize_est(iinfos, &pair, pe_opt->ap_prior, genome_length);
        if( ks ) pair_print_bam(ks,&pair);
        else pair_print_custom(output,&pair);
        bam_destroy_pair(&pair);
    }
    zmq_close( socket );
    fprintf(stderr, "[run_output_thread] finished cleanly, shutting down\n");
    return 0;
}

// Hm.  Should be rather easy:  read message, align, write message.
// Both sockets are bound, and always to the inproc endpoints.
//
// We bind to inproc://work_{in,out}.  When running standalone, the
// other end is the I/O multiplexor, when running as worker process, the
// other end is a streaming device.
//
// Could be split into a lean-memory and a full version, but that's
// probably not worth the hassle when we're moving on to 0.7 anyway. 
void *run_worker_thread( void *arg )
{
    zmq_msg_t msg;
    bam_pair_t pair;
    int failure_count = 0 ;
    int rc;
    uint64_t n_tot[2], n_mapped[2] ;
    kh_64_t *my_hash ;
    khiter_t iter ;
	my_hash = kh_init(64);

    void *incoming = zmq_socket( zmq_context, ZMQ_PULL ) ;
    void *outgoing = zmq_socket( zmq_context, ZMQ_PUSH ) ;
    xassert(incoming, "error creating socket");
    xassert(outgoing, "error creating socket");
    set_sockopts( incoming ) ;
    set_sockopts( outgoing ) ;

    zcheck( zmq_connect( incoming, "inproc://work_out" ), 1, "zmq_connect failed" );
    zcheck( zmq_connect( outgoing, "inproc://work_in"  ), 1, "zmq_connect failed" );
    do {
        zmq_msg_init(&msg);
        if( 0>zmq_msg_recv(&msg, incoming, 0) ) {
            // clean exit if we got terminated here
            if( zmq_errno() == ETERM || zmq_errno() == EINTR ) break ;
            fprintf( stderr, "zmq_msg_recv failed: %s\n", zmq_strerror( zmq_errno() ) ) ;
            exit(1) ;
        }
        pair_init_from_msg(&pair, &msg);
        zmq_msg_close(&msg);

        switch( pair.phase ) {
            case pristine: pair_aln(&pair);     // first two steps can run together
            case aligned:  pair_posn(&pair);
                           break ;      // but then we must send the result back
                                        // until we receive it again, then we finish it
            case positioned: if( g_iinfos ) pair_finish(&pair, g_iinfos, n_tot, n_mapped, my_hash) ;
                             else ++failure_count ;    
            case finished: break ;              // boring...
        }

        msg_init_from_pair(&msg,&pair);
        bam_destroy_pair(&pair);
        rc = zmq_msg_send(&msg, outgoing, 0);
        zmq_msg_close(&msg);
        zterm( rc, 1, "zmq_msg_send failed" ) break;
        if( failure_count >= 1024 ) {
            fprintf( stderr, "[run_worker_thread] Lots of failures due to missing insert size information.\n" ) ;
            fprintf( stderr, "[run_worker_thread] Terminating due to suspected communication problem.\n" ) ;
            s_interrupted = 1 ;
        }
    } while(pair.kind) ;
    fprintf( stderr, "[run_worker_thread] exiting.\n" ) ;
	for (iter = kh_begin(my_hash); iter != kh_end(my_hash); ++iter)
		if (kh_exist(my_hash, iter))
            free(kh_val(my_hash, iter).a);
	kh_destroy(64, my_hash);
    zmq_close(outgoing);
    zmq_close(incoming);
    return 0;
}

// Control loop around the worker threads.  Here we manage a sliding
// window of outstanding records, multiplex between input, output, the
// network (or worker threads), and sort the output again.
//
// The ring buffer contains records of type bam_pair_t.  A record with
// kind==0 and isdone==0 is invalid.  A record with kind==0 and
// isdone==done is the end marker.  Other records are either unacknowledged
// or acknowledged work packets.
//
// FIXME: For the time being, the sliding window has a fixed size.
//        Don't know yet how to make this more dynamic.
struct muxer_args {
    void *workio[2] ;
    enum pair_phase end_phase ;
    const char *addr_in, *addr_out ;
    int cleanup ;
} ;

void *run_io_multiplexor( void* vargs )
{
    struct muxer_args *args = vargs ;
    // Let's see... we're multiplexing between four sockets:
    // - input from BAM         (incoming_bam)
    // - output to BAM          (outgoing_bam)
    // - fan-out to workers     (work_out)
    // - fan-in from workers    (work_in)
    void *incoming_bam = zmq_socket( zmq_context, ZMQ_PULL ) ;
    void *outgoing_bam = zmq_socket( zmq_context, ZMQ_PUSH ) ;

    set_sockopts( incoming_bam ) ;
    set_sockopts( outgoing_bam ) ;

    zcheck( zmq_connect( incoming_bam, args->addr_in ), 1, "zmq_connect failed" );
    zcheck( zmq_connect( outgoing_bam, args->addr_out ), 1, "zmq_connect failed" );

    // We use an array as a ring buffer.  Three pointers point at
    // - the first free slot (index),
    // - the next record to be sent on if it is done (recno!)
    // - the next record to be resent if necessary (index).

    // The ring buffer is empty if (next_free == next_expected), full if
    // (next_free == next_expected-1).

    bam_pair_t *ringbuf = calloc( ring_size, sizeof(bam_pair_t) );
    int next_free = 0, next_expected = 0, next_resend = 0 ;
    int total_out = 0, total_resends = 0, total_in = 0, total_dups = 0 ;
    int current_undone = 0, current_done = 0, current_invalid = ring_size ;
    int rn ;
    struct timeval tv0, tv1;
    gettimeofday( &tv0, 0 ); tv0.tv_sec -= 5;

    // When do we poll what?
    // - Output is not polled, we send when we're ready and block if
    //   necessary.
    // - We poll input if the network is ready, our queue is not full
    //   and input has not ended yet.
    // - We always poll for fan-out.
    // - We always poll for fan-in.

    zmq_pollitem_t polls[] = {
        { socket :    incoming_bam, fd : 0, events :  ZMQ_POLLIN, revents : 0 },
        { socket : args->workio[1], fd : 0, events :           0, revents : 0 },
        { socket : args->workio[0], fd : 0, events :  ZMQ_POLLIN, revents : 0 }
    } ;

    while (!s_interrupted) {
        xassert( current_done + current_undone + current_invalid == ring_size, "logic error in queueing thread (1)" ) ;
        gettimeofday( &tv1, 0 ) ;
        if( loudness >= 2 || tdiff(&tv0,&tv1) >= 10 ) {
            fprintf( stderr, "[run_io_multiplexor] polling -- %d records received, "
                             "%d written, %d resends; %d done in queue, %d undone in queue.\n",
                     total_in, total_out, total_resends, current_done, current_undone ) ;
            memmove( &tv0, &tv1, sizeof(struct timeval) );
        }
        if( loudness >= 3 )
            fprintf( stderr, "next expected: %d, next free: %d\n", next_expected % ring_size, next_free % ring_size ) ; 

        // OMFG, this is becoming unmaintainable quickly...
        zterm( zmq_poll( incoming_bam ? polls : polls+1, incoming_bam ? 3 : 2, -1 ), 1, "zmq_poll failed" ) break;
        if( s_interrupted ) break ;

        // If something is ready, what do we do?
        // - If fresh input arrives, send it to the network and put it in
        //   the ring buffer.  (If we did poll, the net was already
        //   guaranteed to be ready.)
        if( polls[0].revents & ZMQ_POLLIN) { // input available
            // next time, do not poll for input (we change that later
            // once we know we can handle more input)
            polls[0].events = 0 ;

            // the buffer should not be full (otherwise it's a logic error)
            xassert( (next_free+1) % ring_size != next_expected % ring_size, "logic error in queueing thread (2)" ) ;
            if( loudness >= 3 ) fprintf( stderr, "got input record.\n" ) ;
            zmq_msg_t msg ;
            zmq_msg_init( &msg ) ;
            zterm( zmq_msg_recv( &msg, incoming_bam, 0 ), 1, "zmq_recv failed" ) break;

            // input comes ordered, so the correct slot should be
            // (next_free)
            rn = *(int*)zmq_msg_data(&msg) ;
            if( next_free != rn )
            {
                bam_pair_t p;
                pair_init_from_msg( &p, &msg ) ;
                fprintf( stderr, "%d != %d\n", next_free, *(int*)zmq_msg_data(&msg) ) ;
                fprintf( stderr, "kind == %d, isdone == %d, recno == %d\n", p.kind, p.phase, p.recno ) ;
                exit(2);
            }

            xassert( !ringbuf[next_free % ring_size].kind && ringbuf[next_free % ring_size].phase < args->end_phase,
                     "logic error in queueing thread (3)" ) ;
            pair_init_from_msg( &ringbuf[next_free % ring_size], &msg ) ;
            current_invalid-- ;
            if( ringbuf[next_free % ring_size].kind ) {
                total_in += ringbuf[next_free % ring_size].kind ;
                if( !current_undone ) next_resend = next_free % ring_size;
                current_undone++ ;
                if( loudness >= 3 )
                    fprintf( stderr, "received recno %d (%d), will send.\n",
                             ringbuf[next_free % ring_size].recno, ringbuf[next_free % ring_size].kind ) ;
                zterm( zmq_msg_send( &msg, args->workio[1], 0 ), 1, "zmq_send failed" ) break;

                polls[1].events = ZMQ_POLLOUT ;
                next_free++;
                zmq_msg_close(&msg) ;
            } else {
                // How to shut down?
                // - If EOF arrives from input, we shut down input and enqueue the
                //   EOF marker (as already acknowledged).
                ringbuf[next_free % ring_size].phase = finished ;
                xassert( incoming_bam, "repeated closing of incoming_bam" ) ;
                zmq_close( incoming_bam ) ;
                incoming_bam = 0 ;
                polls[0].events = 0 ;
                polls[0].revents = 0 ;
                current_done++ ;
                next_free++;
                zmq_msg_close(&msg) ;
                goto received_one ;
            }
        }

        // - If the fan-out is ready and the queue is full, resend an
        //   unacknowledged packet.
        // - If we're out of input, resend an unacknowledged packet.  We
        //   will be sending out packets rapidly, duplicating some work,
        //   but we never get stuck this way.  We must make sure not to
        //   resend if the queue is full of finished stuff.
        if(polls[1].revents & ZMQ_POLLOUT) { // network ready
            if( (next_free+1) % ring_size == next_expected % ring_size || !incoming_bam ) {
                if( current_undone ) {
                    if( loudness >= 3 ) fprintf( stderr, "ready to send, will send recno %d (%d).\n",
                            ringbuf[next_resend].recno, ringbuf[next_resend].kind ) ;
                    zmq_msg_t msg ;
                    msg_init_from_pair( &msg, &ringbuf[next_resend] ) ;
                    zterm( zmq_msg_send( &msg, args->workio[1], 0 ), 1, "zmq_send failed" ) break;

                    do next_resend = (next_resend+1) % ring_size ;
                    while( ringbuf[next_resend].phase >= args->end_phase || !ringbuf[next_resend].kind ) ;
                    total_resends++ ;
                }
            }
            // next time, poll for input
            else polls[0].events = ZMQ_POLLIN ;
        }

        // - If fan-in arrives, put it in the queue.  If the next packet to
        //   be sent out arrived, start sending what is ready.
        if(polls[2].revents & ZMQ_POLLIN) { // work done
            if( loudness >= 3 ) fprintf( stderr, "work done.\n" ) ;
            zmq_msg_t msg ;
            zmq_msg_init( &msg ) ;
            zterm( zmq_msg_recv( &msg, args->workio[0], 0 ), 1, "zmq_recv failed" ) break;
            
            rn = *(int*)zmq_msg_data(&msg) ;
            if( loudness >= 3 ) fprintf( stderr, "[run_io_multiplexor] received record %d.\n", rn ) ;
            if( rn < next_expected ) {
                if( loudness >= 1 ) fprintf( stderr, "[run_io_multiplexor] this is old shit: %d.\n", rn ) ;
                total_dups++ ;
                zmq_msg_close(&msg ) ;
            }
            else if( rn >= next_free ) {
                if( loudness >= 1 ) fprintf( stderr, "[run_io_multiplexor] this comes from the future: %d.\n", rn ) ;
                zmq_msg_close(&msg ) ;
            }
            else if( ringbuf[rn % ring_size].phase >= args->end_phase ) {
                if( loudness >= 1 ) fprintf( stderr, "[run_io_multiplexor] this is a duplicate: %d.\n", rn ) ;
                total_dups++ ;
                zmq_msg_close(&msg ) ;
            }
            else {
                if( loudness >= 3 ) fprintf( stderr, "this is not a duplicate, dealing with it\n" ) ;
                bam_destroy_pair( &ringbuf[ rn % ring_size ] ) ;
                pair_init_from_msg( &ringbuf[ rn % ring_size ], &msg ) ;
                // I think this is implicit now:
                xassert( ringbuf[ rn % ring_size ].kind, "wrong record type in queue" ) ;
                zmq_msg_close(&msg ) ;
                if( ringbuf[ rn % ring_size ].phase >= args->end_phase ) {
                    current_undone-- ;
                    current_done++ ;
                }

received_one:
                if( current_undone && next_resend == rn % ring_size ) 
                    do next_resend = (next_resend+1) % ring_size ;
                    while( ringbuf[next_resend].phase >= args->end_phase || !ringbuf[next_resend].kind ) ;

                int foo = 0 ;
                while( next_expected % ring_size != next_free % ring_size &&
                        ringbuf[next_expected % ring_size].phase >= args->end_phase )
                {
                    int mykind = ringbuf[next_expected % ring_size].kind ;
                    if( loudness >= 3 ) fprintf( stderr, "next expected: %d, next free: %d\n",
                            next_expected % ring_size, next_free % ring_size ) ;

                    // make message from pair, then invalidate the slot
                    // so it doesn't get sent again.
                    msg_init_from_pair( &msg, &ringbuf[next_expected % ring_size] ) ;
                    bam_destroy_pair( &ringbuf[next_expected % ring_size] ) ;
                    ringbuf[next_expected % ring_size].kind = 0 ; 
                    ringbuf[next_expected % ring_size].phase = pristine ; 
                    // shouldn't this be implicit now:
                    // ringbuf[next_expected % ring_size].phase = pristine ;

                    // this should be dispatched based on message phase
                    zterm( zmq_msg_send( &msg, outgoing_bam, 0 ), 1, "zmq_send failed" ) goto forced_exit ;
                    zmq_msg_close(&msg) ;
                    total_out+=mykind;
                    current_done--;
                    current_invalid++;

                    // How to shut down?
                    // - If EOF is ready to be sent out, shut everything
                    //   down and send the EOF.
                    if( !mykind ) { if( args->cleanup ) goto forced_exit ; else goto clean_exit ; }
                    next_expected++ ;
                    if( loudness >= 4 ) fprintf( stderr, "spinning(1)\n" ) ;
                    xassert( foo++ < ring_size, "logic error in queueing thread (4)" ) ;
                }

                if( current_undone ) {
                    int foo = 0 ;
                    while( ringbuf[next_resend].phase >= args->end_phase ) {
                        next_resend = (next_resend+1) % ring_size ;
                        xassert( foo++ < ring_size, "logic error in queueing thread (5)" ) ;
                    }
                }
            }
        }
    }
    
forced_exit:
    zmq_xclose( args->workio[0] ) ;
    zmq_xclose( args->workio[1] ) ;

clean_exit:
    fprintf( stderr, "[run_io_multiplexor] finished: %d records in total, %d resends, %d dups.\n",
            total_out, total_resends, total_dups ) ;
    if( incoming_bam ) zmq_close( incoming_bam ) ;
    zmq_close( outgoing_bam ) ;
    free( ringbuf ) ;
    return 0 ;
}

// This is the driver loop for bam2bam.  Depending on setting, we run
// different strategies:
//
// - One thread and no listening port: Plain loop.
// - No threads, but a listening port: Don't load index, run I/O threads.
// - Threads and no listening port: Load index, run both I/O and worker
//   threads.
// - Threads and listening port: Load index, run both I/O and worker
//   threads, also bind to listening port.
// - Worker mode: connect, run worker threads and two streamer devices.

void bwa_bam2bam_core( const char *prefix, char* tmpdir, bwa_seqio_t *ks, BGZF *output )
{
    struct muxer_args muxargs ;
    void *broadcast ;
    khint_t iter;
    char *tmpname = alloca( strlen(tmpdir) + 11 ) ;
    strcpy( tmpname, tmpdir ) ;
    strcat( tmpname, "/bwaXXXXXX" ) ;

    // Initialize genome index.  Needed if and only if we're going to
    // have a worker thread.  Else we still need a handful of numbers
    // from the genome index.
    if( gap_opt->n_threads ) init_genome_index( prefix,1 ) ;
    else init_genome_params( prefix ) ;
    fprintf( stderr, "[bwa_bam2bam_core] genome length is %ld\n", (long)genome_length ) ;

    khash_t(isize_infos) *iinfos = kh_init( isize_infos ) ;
    srand48(bns->seed);

    if( !gap_opt->n_threads && !listen_port ) {
        fprintf( stderr, "[bwa_bam2bam_core] No threads and no listening port specified, nothing to do.\n" ) ;
        exit(1);
    } 
    
    int tmpfd1 = mkstemp( tmpname ) ;
    xassert( tmpfd1, "could not create temporary file" ) ;
    fprintf( stderr, "[bwa_bam2bam_core] buffering data in %s\n", tmpname ) ;
    int tmpfd2 = dup(tmpfd1) ;
    xassert( tmpfd1, "could not duplicate file descriptor" ) ;

    gzFile tmpout = gzdopen( tmpfd1, "wb2" ) ;
    
    // Can we get away with a simple loop?
    if( gap_opt->n_threads == 1 && !listen_port ) {
        sequential_loop_pass1( ks, tmpout, iinfos ) ;
        if( gzclose( tmpout ) != Z_OK ) {
            fprintf( stderr, "[%s] error closing temporary file (%s)", __FUNCTION__, strerror(errno) ) ;
            exit(1) ;
        }

        infer_all_isizes( iinfos, pe_opt->ap_prior, genome_length ) ;
        g_iinfos = iinfos ;
        off_t new_off = lseek(tmpfd2, 0, SEEK_SET) ;
        xassert( new_off == 0, "failed to seek in temporary file" ) ;

        gzFile tmpin = gzdopen( tmpfd2, "rb" ) ;
        sequential_loop_pass2( tmpin, output, iinfos ) ;
        if( gzclose( tmpin ) != Z_OK ) {
            fprintf( stderr, "[%s] error closing temporary file (%s)", __FUNCTION__, strerror(errno) ) ;
            exit(1) ;
        }
    } else {
        int n;
        pthread_t config_service_tid ;
        pthread_t reader_tid, reader_tid2 ;
        pthread_t mux_tid, mux_tid2 ;

        // Okay, we need 0MQ.
        s_catch_signals();
        zmq_context = zmq_init(1) ;
        xassert( zmq_context, "0MQ init failed" ) ;

        broadcast = zmq_socket( zmq_context, ZMQ_PUB ) ;
        xassert(broadcast, "couldn't create socket");

        // Are we going to be network aware?  Then we need a
        // configuration service.  And a broadcast channel.
        if( listen_port ) {
            char addr[100] ;
            pthread_create( &config_service_tid, 0, run_config_service, (void*)prefix ) ;
            pthread_detach( config_service_tid ) ;
            snprintf( addr, 100, "tcp://*:%d", listen_port+3 ) ;
            zcheck( zmq_bind( broadcast, addr ), 1, "zmq_bind failed" );
            set_sockopts( broadcast ) ;
        }

        muxargs.workio[0] = zmq_socket( zmq_context, ZMQ_PULL ) ;
        muxargs.workio[1] = zmq_socket( zmq_context, ZMQ_PUSH ) ;
        
        set_sockopts(muxargs.workio[1]);
        set_sockopts(muxargs.workio[0]);

        zcheck( zmq_bind( muxargs.workio[1], "inproc://work_out" ), 1, "zmq_bind failed" );
        zcheck( zmq_bind( muxargs.workio[0], "inproc://work_in" ), 1, "zmq_bind" );
        
        if( listen_port ) {
            char addr[100] ;
            snprintf( addr, 100, "tcp://*:%d", listen_port+1 ) ;
            zcheck( zmq_bind( muxargs.workio[1], addr ), 1, "zmq_bind failed" );
            snprintf( addr, 100, "tcp://*:%d", listen_port+2 ) ;
            zcheck( zmq_bind( muxargs.workio[0], addr ), 1, "zmq_bind failed" );
        }

        void *bam_in = zmq_socket( zmq_context, ZMQ_PUSH ) ;
        xassert(bam_in, "couldn't create socket");
        set_sockopts(bam_in);
        zcheck( zmq_bind( bam_in, "inproc://bam_in" ), 1, "zmq_bind failed" );

        // Start reading in the background...
        struct reader_thread_args rargs = { bam_in, ks, 0 } ;
        pthread_create( &reader_tid, 0, run_reader_thread, &rargs );
        pthread_detach( reader_tid ) ;

        // Are threads requested?  Then fire them up.
        for(n=0;n!=gap_opt->n_threads;++n) {
            pthread_t worker_tid ;
            pthread_create(&worker_tid, 0, run_worker_thread, 0 ) ;
            pthread_detach(worker_tid);
        }

        void *tmp_out = zmq_socket( zmq_context, ZMQ_PULL ) ;
        xassert(tmp_out, "couldn't create socket");
        set_sockopts(tmp_out);
        zcheck( zmq_bind( tmp_out, "inproc://tmp_out" ), 1, "zmq_bind failed" );

        // ...multiplex records between I/O and workers...
        muxargs.end_phase = positioned ;
        muxargs.addr_in = "inproc://bam_in" ;
        muxargs.addr_out = "inproc://tmp_out" ;
        muxargs.cleanup = 0 ;
        pthread_create( &mux_tid, 0, run_io_multiplexor, &muxargs );
        pthread_detach( mux_tid ) ;

        // ...and shift to output in the foreground.
        run_output_thread( 0, tmpout, tmp_out, iinfos );

        // when we get here, we got the last record and an EOF marker,
        // so round 1 is finished
        if( gzclose( tmpout ) != Z_OK ) {
            fprintf( stderr, "[%s] error closing temporary file (%s)", __FUNCTION__, strerror(errno) ) ;
            exit(1) ;
        }
        close( tmpfd1 ) ;

        if( !s_interrupted ) {
            infer_all_isizes( iinfos, pe_opt->ap_prior, genome_length ) ;
            g_iinfos = iinfos ;

            // broadcast here iff we're networked.  doing it without the
            // network would be no harm, but why do unnecessary work?
            if( listen_port ) {
                zmq_msg_t iinfo_msg ;
                zmq_msg_init_size( &iinfo_msg, 1+iinfo_encoded_size(iinfos) ) ;
                char *p = zmq_msg_data(&iinfo_msg) ; 
                char *q = p + zmq_msg_size(&iinfo_msg) ;
                *p = 2 ;
                xassert( q == encode_iinfo( p+1, iinfos ), "failed encoding isize info" ) ;
                fprintf( stderr, "[bwa_bam2bam_core] broadcasting %d bytes of insert size info\n", (int)zmq_msg_size(&iinfo_msg)-1 ) ;
                zcheck( zmq_msg_send( &iinfo_msg, broadcast, 0 ), 1, "zmq_send failed" ) ;
            }

            off_t new_off = lseek(tmpfd2, 0, SEEK_SET) ;
            xassert( new_off == 0, "failed to seek in temporary file" ) ;
            gzFile tmpin = gzdopen( tmpfd2, "rb" ) ;

            // At this point, all the workers and muxers are still connected
            // and have been supplied with isize information.  The timing
            // may screw this up, but if workers receive the isize info
            // late, their work will simply be repeated.  All we need is
            // another loop like the one above... but it should be
            // parameterized, not cloned.

            void *tmp_in = zmq_socket( zmq_context, ZMQ_PUSH ) ;
            xassert(tmp_in, "couldn't create socket");
            set_sockopts(tmp_in);
            zcheck( zmq_bind( tmp_in, "inproc://tmp_in" ), 1, "zmq_bind failed" );

            // Start reading in the background...
            struct reader_thread_args rargs2 = { tmp_in, 0, tmpin } ;
            pthread_create( &reader_tid2, 0, run_reader_thread, &rargs2 );
            pthread_detach( reader_tid2 ) ;

            void *bam_out = zmq_socket( zmq_context, ZMQ_PULL ) ;
            xassert(bam_out, "couldn't create socket");
            set_sockopts(bam_out);
            zcheck( zmq_bind( bam_out, "inproc://bam_out" ), 1, "zmq_bind failed" );

            // ...multiplex records between I/O and workers...
            muxargs.end_phase = finished ;
            muxargs.addr_in = "inproc://tmp_in" ;
            muxargs.addr_out = "inproc://bam_out" ;
            muxargs.cleanup = 1 ;
            pthread_create( &mux_tid2, 0, run_io_multiplexor, &muxargs );
            pthread_detach( mux_tid2 ) ;

            run_output_thread( output, 0, bam_out, iinfos ) ;
            if( gzclose( tmpin ) != Z_OK ) {
                fprintf( stderr, "[%s] error closing temporary file (%s)", __FUNCTION__, strerror(errno) ) ;
                exit(1) ;
            }
            close( tmpfd2 ) ;
        }

        // Done, tell everyone and get rid of 0MQ.
        zmq_msg_t finmsg ;
        zmq_msg_init_data( &finmsg, "\1go away", 8, 0, 0 ) ;
        fprintf( stderr, "[bwa_bam2bam_core] sending termination signal\n" ) ;
        zcheck( zmq_msg_send( &finmsg, broadcast, 0 ), 1, "zmq_send failed" ) ;
        fprintf( stderr, "[bwa_bam2bam_core] done, waiting for clean shutdown\n" ) ;

        zmq_xclose( broadcast ) ;
        zmq_term( zmq_context ) ;
    }

    unlink( tmpname ) ;
	if (pac)    bwt_destroy_pac(pac,bns);
	if (ntbns)  bns_destroy(ntbns);
	if (bns)    bns_destroy(bns);
	if (bwt[0]) bwt_destroy(bwt[0]);
    if (bwt[1]) bwt_destroy(bwt[1]);

	for (iter = kh_begin(iinfos); iter != kh_end(iinfos); ++iter)
    {
        if (kh_exist(iinfos, iter)) {
            free((char*)kh_key(iinfos, iter));
            free(kh_val(iinfos, iter).hist);
        }
    }
    kh_destroy(isize_infos, iinfos);
}

int bwa_bam_to_bam( int argc, char *argv[] )
{
	int c, opte = -1;
    char *saif[3] = {0,0,0} ;
    char *ofile = 0;
    char *tmpname = "/var/tmp" ;
    char *prefix = 0;

	gap_opt = gap_init_opt();
	pe_opt = bwa_init_pe_opt();

	while ((c = getopt_long(argc, argv, "g:n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:C:D:a:sc:h:H:Ap:0:1:2:", longopts, 0)) >= 0) {
		switch (c) {
        case 'g': prefix = optarg; break;
		case 'n':
			if (strstr(optarg, ".")) gap_opt->fnr = atof(optarg), gap_opt->max_diff = -1;
			else gap_opt->max_diff = atoi(optarg), gap_opt->fnr = -1.0;
			break;
		case 'o': gap_opt->max_gapo = atoi(optarg); break;
		case 'e': opte = atoi(optarg); break;
		case 'M': gap_opt->s_mm = atoi(optarg); break;
		case 'O': gap_opt->s_gapo = atoi(optarg); break;
		case 'E': gap_opt->s_gape = atoi(optarg); break;
		case 'd': gap_opt->max_del_occ = atoi(optarg); break;
		case 'i': gap_opt->indel_end_skip = atoi(optarg); break;
		case 'l': gap_opt->seed_len = atoi(optarg); break;
		case 'k': gap_opt->max_seed_diff = atoi(optarg); break;
		case 'm': gap_opt->max_entries = atoi(optarg); break;
		case 't': gap_opt->n_threads = atoi(optarg); break;
		case 'L': gap_opt->mode |= BWA_MODE_LOGGAP; break;
		case 'R': gap_opt->max_top2 = atoi(optarg); break;
		case 'q': gap_opt->trim_qual = atoi(optarg); break;
		case 'N': gap_opt->mode |= BWA_MODE_NONSTOP; gap_opt->max_top2 = 0x7fffffff; break;
		case 'f': ofile = optarg; break;
		case 'C': pe_opt->max_occ = atoi(optarg); break;
		case 'D': pe_opt->max_occ_se = atoi(optarg); break;
		case 'a': pe_opt->max_isize = atoi(optarg); break;
		case 's': pe_opt->is_sw = 0; break;
		case 'c': pe_opt->ap_prior = atof(optarg); break;
		case 'A': pe_opt->force_isize = 1; break;
		case 'h': pe_opt->n_multi = atoi(optarg); break;
		case 'H': pe_opt->N_multi = atoi(optarg); break;
		case 'p': listen_port = atoi(optarg); break;
        case '0':
        case '1':
        case '2': saif[c-'0'] = optarg; break;
        case 128: only_aligned = 1; break;
        case 129: debug_bam = 1; break;
        case 130: broken_input = 1; break;
        case 131: skip_duplicates = 1; break;
        case 132: tmpname = *optarg ? optarg : "." ; break ;
		default: return 1;
		}
	}
	if (opte > 0) {
		gap_opt->max_gape = opte;
		gap_opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 1 > argc || !prefix) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa bam2bam [options] <in.bam>\n\n");
		fprintf(stderr, "Options: -g, --genome PREFIX               prefix of genome index files [%s]\n", prefix);
        fprintf(stderr, "         -f, --output FILE                 file to write output to instead of stdout\n");
        fprintf(stderr, "\n");
		fprintf(stderr, "         -n, --num-diff NUM                max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n", BWA_AVG_ERR, gap_opt->fnr);
		fprintf(stderr, "         -o, --max-gap-open INT            maximum number or fraction of gap opens [%d]\n", gap_opt->max_gapo);
		fprintf(stderr, "         -e, --max-gap-extensions INT      maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr, "         -i, --indel-near-end INT          do not put an indel within INT bp towards the ends [%d]\n", gap_opt->indel_end_skip);
		fprintf(stderr, "         -d, --deletion-occurences INT     maximum occurrences for extending a long deletion [%d]\n", gap_opt->max_del_occ);
		fprintf(stderr, "         -l, --seed-length INT             seed length [%d]\n", gap_opt->seed_len);
		fprintf(stderr, "         -k, --seed-mismatches INT         maximum differences in the seed [%d]\n", gap_opt->max_seed_diff);
		fprintf(stderr, "         -M, --mismatch-penalty INT        mismatch penalty [%d]\n", gap_opt->s_mm);
		fprintf(stderr, "         -O, --gap-open-penalty INT        gap open penalty [%d]\n", gap_opt->s_gapo);
		fprintf(stderr, "         -E, --gap-extension-penalty INT   gap extension penalty [%d]\n", gap_opt->s_gape);
        fprintf(stderr, "\n");
		fprintf(stderr, "         -m, --queue-size INT              maximum entries in the queue [%d]\n", gap_opt->max_entries);
		fprintf(stderr, "         -R, --max-best-hits INT           stop searching when there are >INT equally best hits [%d]\n", gap_opt->max_top2);
		fprintf(stderr, "         -q, --trim-quality INT            quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, gap_opt->trim_qual);
		fprintf(stderr, "         -B, --barcode-length INT          length of barcode\n");
		fprintf(stderr, "         -L, --log-gap-penalty             log-scaled gap penalty for long deletions\n");
		fprintf(stderr, "         -N, --non-iterative               non-iterative mode: search for all n-difference hits (slooow)\n");
		fprintf(stderr, "         -a, --max-insert-size INT         maximum insert size [%d]\n", pe_opt->max_isize);
		fprintf(stderr, "         -C, --max-occurences INT          maximum occurrences for one end of a pair [%d]\n", pe_opt->max_occ);
		fprintf(stderr, "         -D, --max-occurences-se INT       maximum occurrences for a single ended read [%d]\n", pe_opt->max_occ_se);
		fprintf(stderr, "         -h, --max-hits INT                maximum hits to output [%d]\n", pe_opt->n_multi);
		fprintf(stderr, "         -H, --max-discordant-hits INT     maximum hits to output for discordant pairs [%d]\n", pe_opt->N_multi);
        fprintf(stderr, "         -c, --chimeric-rate FLOAT         prior of chimeric rate (lower bound) [%.1le]\n", pe_opt->ap_prior);
		fprintf(stderr, "         -s, --disable-sw                  disable Smith-Waterman for the unmapped mate\n");
		fprintf(stderr, "         -A, --disable-isize-estimate      disable insert size estimate (force -s)\n");
        fprintf(stderr, "             --only-aligned                output only aligned reads or pairs\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "         -p, --listen-port PORT            listen for workers on PORT [%d]\n", listen_port);
		fprintf(stderr, "         -t, --num-threads INT             number of threads [%d]\n", gap_opt->n_threads);
        fprintf(stderr, "             --debug-bam                   add additional fields to BAM output to aid debugging\n");
        fprintf(stderr, "             --broken-input                ignore lone mates in input BAM (not recommended)\n");
        fprintf(stderr, "             --skip-duplicates             do not bother mapping reads marked as duplicates\n");
        fprintf(stderr, "             --temp-dir                    location of intermediate file [%s]\n", tmpname);
		fprintf(stderr, "         -0, -1, -2                        provide up to three sai files to resume from (not recommended)\n");
		fprintf(stderr, "\n");
        if( !prefix ) 
            fprintf(stderr, "No genome prefix specified.\n\n");
		return 1;
	}
	if (gap_opt->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, gap_opt->fnr);
			if (l != k) fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l);
			k = l;
		}
	}

    bam_header_t *hdr;
	bwa_seqio_t *ks = bwa_bam_open(argv[optind], 7, saif, gap_opt, &hdr);
    bns = bns_restore(prefix);

    BGZF* output = ofile ? bgzf_open(ofile, "w2") : bgzf_fdopen(1, "w2") ;
    if( 0>bwa_print_bam_header(output, bns, hdr->text, argc, argv) ) {
        fprintf( stderr, "[bwa_bam2bam_core] Error writing BAM header.\n" ) ;
        exit(1);
    }
    bam_header_destroy(hdr);

    bwa_bam2bam_core(prefix, tmpname, ks, output);
	bwa_seq_close(ks);

    bgzf_close( output );
    if( !s_interrupted ) final_rename( "bam2bam", ofile ) ;

	free(pe_opt);
	free(gap_opt);
	return 0;
}

void *run_device( void *args )
{
    void **socks = args ;
    zmq_device( ZMQ_STREAMER, socks[0], socks[1] ) ;
    return 0 ;
}

int handle_broadcast( zmq_msg_t *m, void* sock )
{
    zmq_msg_recv( m, sock, 0 );
    fprintf( stderr, "[bwa_worker_core] Received broadcast of %d bytes.\n", (int)zmq_msg_size(m) ) ;
    if( zmq_msg_size(m) == 0 ) return 1 ;
    if( *(char*)zmq_msg_data(m) == 1 ) return 0 ; // the termination code
    if( *(char*)zmq_msg_data(m) == 2 ) {
        // received isize info!!
        fprintf( stderr, "[bwa_worker_core] Received broadcast of %d bytes of isize info.\n", (int)zmq_msg_size(m)-1 ) ;
        char *p = zmq_msg_data(m)+1 ;
        char *q = zmq_msg_data(m)+zmq_msg_size(m) ;
        if( p != q ) {
            if( g_iinfos ) kh_destroy(isize_infos, g_iinfos);
            g_iinfos = decode_iinfo( p, q ) ;
        }
        return 1 ;
    }
    return -1 ; // Huh?
}

void bwa_worker_core( int nthreads, char* host, int port ) 
{
    time_t start_time = time(0) ;
    s_catch_signals();

    // Need a streamer that writes into work out and one that reads from
    // work_in.  (Those would mimic the role of the multiplexor.)
    void *socks[5] ;
    socks[0] = zmq_socket( zmq_context, ZMQ_PULL ) ;
    socks[1] = zmq_socket( zmq_context, ZMQ_PUSH ) ;
    socks[2] = zmq_socket( zmq_context, ZMQ_PULL ) ;
    socks[3] = zmq_socket( zmq_context, ZMQ_PUSH ) ;
    socks[4] = zmq_socket( zmq_context, ZMQ_SUB ) ;

    set_sockopts(socks[0]);
    set_sockopts(socks[1]);
    set_sockopts(socks[2]);
    set_sockopts(socks[3]);
    set_sockopts(socks[4]);

    char addr[100] ;
    snprintf( addr, 100, "tcp://%s:%d", host, port+1 ) ;
    zcheck( zmq_connect( socks[0], addr ), 1, "zmq_connect failed" );
    zcheck( zmq_bind( socks[1], "inproc://work_out" ), 1, "zmq_bind failed" );
    zcheck( zmq_bind( socks[2], "inproc://work_in" ), 1, "zmq_bind failed" );
    snprintf( addr, 100, "tcp://%s:%d", host, port+2 ) ;
    zcheck( zmq_connect( socks[3], addr ), 1, "zmq_connect failed" );
    snprintf( addr, 100, "tcp://%s:%d", host, port+3 ) ;
    zcheck( zmq_connect( socks[4], addr ), 1, "zmq_connect failed" );

    zcheck( zmq_setsockopt( socks[4], ZMQ_SUBSCRIBE, "", 0 ), 1, "subscription failed" ) ;

    int n;
    for(n=0;n!=nthreads;++n) {
        pthread_t worker_tid ;
        pthread_create(&worker_tid, 0, run_worker_thread, 0 ) ;
        pthread_detach(worker_tid);
    }

    // This device relays stuff from work_out back to the network.
    pthread_t fan_out ;
    pthread_create( &fan_out, 0, run_device, socks ) ;

    // Now run the streamer for the fan-in.  We just ping-pong between
    // receive and send, but if either times out, we simply terminate.
    zmq_pollitem_t pitems[] = {
        { socks[2], 0, ZMQ_POLLIN, 0 },
        { socks[4], 0, ZMQ_POLLIN, 0 },
        { socks[3], 0, ZMQ_POLLOUT, 0 } } ;

    int npoll=0;
    zmq_msg_t m, m2 ;
    zmq_msg_init(&m) ;
    zmq_msg_init(&m2) ;
    while (!s_interrupted && time(0) - start_time < max_run_time*60) {
        npoll = zmq_poll( pitems, 2, 1E3 * timeout ) ;
        zterm( npoll, 1, "zmq_poll failed" ) break;

        if( npoll == 0 ) break ;
        if( pitems[1].revents & ZMQ_POLLIN ) 
            if( !handle_broadcast( &m2, socks[4] ) ) break ;
        zterm( zmq_msg_recv( &m, socks[2], 0 ), 1, "zmq_recv failed" ) break;

        npoll = zmq_poll( pitems+1, 2, 1E3 * timeout ) ;
        zterm( npoll, 1, "zmq_poll failed" ) break;

        if( npoll == 0 ) break ;
        if( pitems[1].revents & ZMQ_POLLIN ) 
            if( !handle_broadcast( &m2, socks[4] ) ) break ;

        zterm( zmq_msg_send( &m, socks[3], 0 ), 1, "zmq_send failed" ) break;
    }

    if( npoll==0 )
        fprintf( stderr, "[bwa_worker_core] No work delivered in %ds.  Terminating.\n", timeout ) ;
    else if( pitems[1].revents & ZMQ_POLLIN ) {
        fprintf( stderr, "[bwa_worker_core] Received termination signal: \"%.*s\".\n",
                (int)zmq_msg_size(&m2)-1, (char*)zmq_msg_data(&m2)+1 ) ;
    } else if( s_interrupted ) 
        fprintf( stderr, "[bwa_worker_core] Received interrupt signal.\n" ) ;
    else 
        fprintf( stderr, "[bwa_worker_core] I've been going for %d minutes, I'm tired now.\n",
                (int)(time(0) - start_time) / 60 ) ;

    zmq_msg_close(&m);
    zmq_msg_close(&m2);

    pthread_cancel( fan_out ) ;
    pthread_join( fan_out, 0 ) ;
    zmq_xclose(socks[4]) ;
    zmq_xclose(socks[3]) ;
    zmq_close(socks[2]) ;
    zmq_close(socks[1]) ;
    zmq_xclose(socks[0]) ;

    // Termination here with initialization outside is a bit unclean,
    // but we need to make sure nobody accesses shared data anymore when
    // we free it.
    zmq_term(zmq_context);

    if (pac) bwt_destroy_pac(pac,bns);
    if (ntbns) bns_destroy(ntbns);
    if (bns) bns_destroy(bns);
    if (bwt[0]) bwt_destroy(bwt[0]);
    if (bwt[1]) bwt_destroy(bwt[1]);
}

int bwa_worker( int argc, char *argv[] )
{
    int c, nthreads = 1, port = 0 ;
    char *host = "localhost" ;
    char *prefix = 0;

	while ((c = getopt_long(argc, argv, "t:h:p:T:", workeropts, 0)) >= 0) {
		switch (c) {
            case 't': nthreads = atoi(optarg) ; break ;
            case 'h': host = optarg ; break ;
            case 'p': port = atoi(optarg) ; break ;
            case 'T': max_run_time = atoi(optarg) ; break ;
            default: return 1;
        }
    }

    if( optind != argc ) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa worker [options]\n\n") ;
		fprintf(stderr, "Options: -t, --num-threads NUM             number of worker threads [%d]\n", nthreads);
        fprintf(stderr, "         -h, --host HOST                   host to connect to [%s]\n", host);
        fprintf(stderr, "         -p, --port NUM                    port to connect to [%d]\n", port);
        fprintf(stderr, "         -T, --timeout NUM                 terminate after NUM minutes [%d]\n", max_run_time);
        return 1;
    }

    zmq_context = zmq_init(1) ;
    xassert( zmq_context, "0MQ init failed" ) ;

    void *conf_sock = zmq_socket( zmq_context, ZMQ_REQ );
    xassert( conf_sock, "couldn't create socket" );

    char addr[100] ;
    snprintf( addr, 100, "tcp://%s:%d", host, port ) ;
    zcheck( zmq_connect( conf_sock, addr ), 1, "zmq_connect failed" );

    struct utsname un ;
    uname( &un ) ;

    zmq_msg_t m;
    zmq_msg_init_size(&m,strlen(un.nodename)+1);
    *(char*)zmq_msg_data(&m) = 0 ;
    memcpy( (char*)zmq_msg_data(&m)+1, un.nodename, strlen(un.nodename) ) ;
    zcheck( zmq_msg_send(&m, conf_sock, 0), 1, "zmq_send failed" );
    // zmq_msg_close(&m);
    // zmq_msg_init(&m);
    zcheck( zmq_msg_recv(&m, conf_sock, 0), 1, "zmq_recv failed" );

	gap_opt = gap_init_opt();
	pe_opt = bwa_init_pe_opt();

    memcpy( gap_opt, zmq_msg_data(&m), sizeof(gap_opt_t) );
    memcpy( pe_opt, zmq_msg_data(&m) + sizeof(gap_opt_t), sizeof(pe_opt_t) );
    prefix = strndup( zmq_msg_data(&m) + sizeof(gap_opt_t) + sizeof(pe_opt_t),
                      zmq_msg_size(&m) - sizeof(gap_opt_t) - sizeof(pe_opt_t) ) ;
    zmq_msg_close(&m);

    // Initialize genome index.  This is basically an mmap, but since
    // the data is read lazily, startup is slow afterwards.  Therefore,
    // we touch all of the index (by looping over the data structures
    // and reading a word every so often), and *then* open the streaming
    // sockets.  Otherwise we risk running into a timeout before the
    // index is was even loaded.
    bns = bns_restore(prefix);
    init_genome_index( prefix,1 ) ;
    srand48(bns->seed);

    zmq_msg_init_size(&m,strlen(un.nodename)+1);
    *(char*)zmq_msg_data(&m) = 1 ;
    memcpy( (char*)zmq_msg_data(&m)+1, un.nodename, strlen(un.nodename) ) ;
    zcheck( zmq_msg_send(&m, conf_sock, 0), 1, "zmq_send failed" );
    // zmq_msg_close(&m);
    // zmq_msg_init(&m);
    zcheck( zmq_msg_recv(&m, conf_sock, 0), 1, "zmq_recv failed" );
    zmq_xclose(conf_sock);
    char *iinfo_start = zmq_msg_data(&m) ;
    char *iinfo_end   = zmq_msg_data(&m) + zmq_msg_size(&m) ;

    if( iinfo_start != iinfo_end ) 
    {
        xassert( !g_iinfos, "got mysterious iinfo from nowhere" ) ;
        g_iinfos = decode_iinfo( iinfo_start, iinfo_end ) ;
    }

    zmq_msg_close(&m);
    bwa_worker_core( nthreads, host, port ) ;

    free(prefix);
	free(pe_opt);
	free(gap_opt);
    return 0;
}

