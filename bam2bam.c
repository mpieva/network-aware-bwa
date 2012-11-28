#define CHUNKSIZE 0x100000

static const int loudness = 0;
static const int ring_size = 0x100000;
static const int timeout = 30;


/* Combine aln, samse, sampe into one workflow.  We steal the workhorse
 * functions from the other subprograms, then wrap them in new top-level
 * functions.  Comes with a command line that combines all the other
 * options, maybe even some more.
 */

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

#include <zmq.h>
#include <pthread.h>
#include <signal.h>
#include <stdarg.h>

static const uint64_t the_hwm = 64;
static const int the_linger = 5000;

KHASH_MAP_INIT_INT64(64, poslist_t)
static kh_64_t *g_hash;

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

// Yes, global variables.  Not very nice, but we'll only load one genome
// and passing it everywhere gets old after a while.

static bwt_t     *bwt[2]      = {0,0};
static bntseq_t  *bns         = 0;
static bntseq_t  *ntbns       = 0;
static ubyte_t   *pac         = 0;

static gap_opt_t *gap_opt     = 0;
static pe_opt_t  *pe_opt      = 0;
static char      *prefix      = 0;
static int       only_aligned = 0;
static int       debug_bam    = 0;
static int       broken_input = 0;

static void      *zmq_context = 0;
static void      *work_out    = 0;
static void      *work_in     = 0;
static int       listen_port  = 0;

static bwa_seqio_t *ks;
static khint_t iter;
static isize_info_t ii;

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

static int s_interrupted = 0;
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

int zcheck_at( int rc, int testterm, int status, const char *label, char *format, ... )
{
    va_list ap;
    if (rc>=0) return 0;
    if (testterm && (zmq_errno() == ETERM || zmq_errno() == EINTR)) return 1;
    fprintf( stderr, "[%s] ", label );
    va_start(ap,format);
    vfprintf( stderr, format, ap );
    va_end(ap);
    fprintf( stderr, "%s%s\n", zmq_errno()?": ":"", zmq_strerror(zmq_errno()) );
    if (status) exit(status) ;
    return 0;
}

#define zcheck(rc, status, ...) (void)zcheck_at(rc,0,status,__FUNCTION__,__VA_ARGS__)
#define zterm(rc, status, ...)     if(zcheck_at(rc,1,status,__FUNCTION__,__VA_ARGS__))

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

    // revcom sequence (incl. one nybble of padding)
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

// Erase unwanted tagged fields:
// AM NM CM SM MD X0 X1 XA XC XG XM XN XO XT YQ
void erase_unwanted_tags(bam1_t *out)
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
            case 'Y': keep = p[1] != 'Q' || !debug_bam; break ;
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
	int j;
    if (p->clip_len < p->full_len) bam_push_int( out, 'X', 'C', p->clip_len );
    if (p->max_entries && debug_bam) bam_push_int( out, 'Y', 'Q', p->max_entries );

	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
		int seqid, nn, am = 0;
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
			out->core.flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences
            p->mapQ = 0;
        }

		// update other flags in output record
		if (mate) {
			if (mate->type != BWA_TYPE_NO_MATCH) {
				if (mate->strand) out->core.flag |= SAM_FMR;
			} else out->core.flag |= SAM_FMU;
		}
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
			int m_seqid;
			am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
			bns_coor_pac2real(bns, mate->pos, mate->len, &m_seqid);

            out->core.mtid = m_seqid ;
            out->core.mpos = mate->pos - bns->anns[m_seqid].offset ;
			if (p->type == BWA_TYPE_NO_MATCH) out->core.isize = 0;
            else out->core.isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
		} else if (mate) {
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
					nn = bns_coor_pac2real(bns, q->pos, j, &seqid);
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
	}
}

static void aln_singleton( bam_pair_t *raw )
{
    int j;
    bwa_seq_t seq_first ;
    erase_unwanted_tags(&raw->first);
    memset( &seq_first, 0, sizeof(bwa_seq_t) ) ;
    bam1_to_seq(&raw->first, &seq_first, 1);

    // Horrible hack: if we have alignments from an external source
    // (recovered from sai file), we use them instead of calculating
    // them.
    if(raw->n_aln1) {
        seq_first.n_aln = raw->n_aln1-1 ;
        seq_first.aln = raw->aln1 ;
        raw->n_aln1 = 0 ;
        raw->aln1 = 0 ;
    }
    else bwa_cal_sa_reg_gap(bwt, 1, &seq_first, gap_opt);

    // from bwa_sai2sam_se_core
    bwa_seq_t *p = &seq_first ;
    bwa_aln2seq_core(p->n_aln, p->aln, p, 1, pe_opt->max_occ_se );

    // from bwa_cal_pac_pos
    bwa_cal_pac_pos_core(bwt[0], bwt[1], p, gap_opt->max_diff, gap_opt->fnr);
    for (j = 0; j < p->n_multi; ++j) {
        bwt_multi1_t *q = seq_first.multi + j;
        if (q->strand) q->pos = bwt_sa(bwt[0], q->pos);
        else q->pos = bwt[1]->seq_len - (bwt_sa(bwt[1], q->pos) + p->len);
    }
    bwa_refine_gapped(bns, 1, p, pac, ntbns);
    bwa_update_bam1(&raw->first, bns, &seq_first, 0, gap_opt->mode, gap_opt->max_top2);
    bwa_free_read_seq1(&seq_first);
}

static void aln_pair( bam_pair_t *raw )
{
    int j;
    pe_data_t d;
    bwa_seq_t seq_first, seq_second ;
    memset( &seq_first, 0, sizeof(bwa_seq_t) ) ;
    memset( &seq_second, 0, sizeof(bwa_seq_t) ) ;
    erase_unwanted_tags(&raw->first);
    erase_unwanted_tags(&raw->second);
    bam1_to_seq(&raw->first, &seq_first, 1);
    bam1_to_seq(&raw->second, &seq_second, 1);

    // Strategy for paired end data: Originally, bwa operated on a block
    // at a time, did alignments, inferred insert size.  This is no
    // longer possible.  Instead, we will collect sequences per read
    // group until we get the insert size estimate, then complete what
    // we collected and use the insert size for everything else.  But
    // right now, we don't do the estimate at all.

    // BEGIN from bwa_sai2sam_pe_core
    // BEGIN from bwa_cal_pac_pos_pe
    // note that seqs is not an array of two arrays anymore!!
    // also note that we don't track this
    // cnt_chg... whatever it was good
    // for.
    // cnt_chg = bwa_cal_pac_pos_pe(seqs, &ii, &last_ii);

    // Horrible hack: if we have alignments from an external source
    // (recovered from sai file or delivered by small worker), we use
    // them instead of calculating them.
    if(raw->n_aln1) {
        seq_first.n_aln = raw->n_aln1-1 ;
        seq_first.aln = raw->aln1 ;
        raw->n_aln1 = 0 ;
        raw->aln1 = 0 ;
    }
    else bwa_cal_sa_reg_gap(bwt, 1, &seq_first, gap_opt);

    if(raw->n_aln2) {
        seq_second.n_aln = raw->n_aln2-1 ;
        seq_second.aln = raw->aln2 ;
        raw->n_aln2 = 0 ;
        raw->aln2 = 0 ;
    }
    else bwa_cal_sa_reg_gap(bwt, 1, &seq_second, gap_opt);

    memset(&d, 0, sizeof(pe_data_t));

    // SE part.  Almost, but not quite the
    // same as for singletons.  (*sigh*)
    seq_first.n_multi = 0 ;
    seq_second.n_multi = 0 ;

    bwa_seq_t *p[2] = { &seq_first, &seq_second } ;
    // fill d.aln from p.aln
    for (j = 0; j < 2; ++j) 
    {
        if(p[j]->n_aln > kv_max(d.aln[j]))
            kv_resize(bwt_aln1_t, d.aln[j], p[j]->n_aln);
        memcpy(d.aln[j].a, p[j]->aln, p[j]->n_aln*sizeof(bwt_aln1_t));
        d.aln[j].n = p[j]->n_aln; 
    }

    // generate SE alignment and mapping quality
    bwa_aln2seq(seq_first.n_aln,  d.aln[0].a, &seq_first);
    bwa_aln2seq(seq_second.n_aln, d.aln[1].a, &seq_second);

    bwa_cal_pac_pos_core(bwt[0], bwt[1], &seq_first,  gap_opt->max_diff, gap_opt->fnr);
    bwa_cal_pac_pos_core(bwt[0], bwt[1], &seq_second, gap_opt->max_diff, gap_opt->fnr);

    // XXX inferring isize is thoroughly broken!
    // infer_isize(n_seqs, seqs, &ii, pe_opt->ap_prior, bwt[0]->seq_len);
    // if (ii.avg < 0.0 && last_ii.avg > 0.0) ii = last_ii;
    // if (pe_opt->force_isize) {
    // fprintf(stderr, "[%s] discard insert size estimate as user's request.\n", __func__);
    // ii.low = ii.high = 0; ii.avg = ii.std = -1.0;
    // }

    // PE

    // WTF?  Restore the old crap?  Who broke it in the first place?
    for (j = 0; j < 2; ++j) 
    {
        if(p[j]->n_aln > kv_max(d.aln[j]))
            kv_resize(bwt_aln1_t, d.aln[j], p[j]->n_aln);
        memcpy(d.aln[j].a, p[j]->aln, p[j]->n_aln*sizeof(bwt_aln1_t));
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
                        khint_t iter = kh_put(64, g_hash, key, &ret);
                        if (ret) { // not in the hash table; ret must equal 1 as we never remove elements
                            poslist_t *z = &kh_val(g_hash, iter);
                            z->n = r->l - r->k + 1;
                            z->a = (bwtint_t*)malloc(sizeof(bwtint_t) * z->n);
                            for (l = r->k; l <= r->l; ++l)
                                z->a[l - r->k] = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);
                        }
                        for (l = 0; l < kh_val(g_hash, iter).n; ++l) {
                            x = kh_val(g_hash, iter).a[l];
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
            /* cnt_chg += */
            pairing(p, &d, pe_opt, gap_opt->s_mm, &ii);
        }
    }

    if (pe_opt->N_multi || pe_opt->n_multi) {
        for (j = 0; j < 2; ++j) {
            if (p[j]->type != BWA_TYPE_NO_MATCH) {
                int k;
                if (!(p[j]->extra_flag&SAM_FPP) && p[1-j]->type != BWA_TYPE_NO_MATCH) {
                    bwa_aln2seq_core(d.aln[j].n, d.aln[j].a, p[j], 0, p[j]->c1+p[j]->c2-1 > pe_opt->N_multi? pe_opt->n_multi : pe_opt->N_multi);
                } else bwa_aln2seq_core(d.aln[j].n, d.aln[j].a, p[j], 0, pe_opt->n_multi);
                for (k = 0; k < p[j]->n_multi; ++k) {
                    bwt_multi1_t *q = p[j]->multi + k;
                    q->pos = q->strand? bwt_sa(bwt[0], q->pos) : bwt[1]->seq_len - (bwt_sa(bwt[1], q->pos) + p[j]->len);
                }
            }
        }
    }

    kv_destroy(d.arr);
    kv_destroy(d.pos[0]); kv_destroy(d.pos[1]);
    kv_destroy(d.aln[0]); kv_destroy(d.aln[1]);

    // END from bwa_cal_pac_pos_pe
    bwa_paired_sw(bns, pac, 1, p, pe_opt, &ii);
    // last_ii = ii;

    // END from bwa_sai2sam_pe_core
    bwa_refine_gapped(bns, 1, p[0], pac, ntbns);
    bwa_refine_gapped(bns, 1, p[1], pac, ntbns);

    // For PE reads, BWA would have concatenated their barcodes.  We
    // don't, for once because we don't identify barcodes, but also
    // because the idea feels wrong.
    bwa_update_bam1( &raw->first,  bns, &seq_first, &seq_second, gap_opt->mode, gap_opt->max_top2);
    bwa_update_bam1( &raw->second, bns, &seq_second, &seq_first, gap_opt->mode, gap_opt->max_top2);
    bwa_free_read_seq1(&seq_second);
    bwa_free_read_seq1(&seq_first);
}

static inline double tdiff( struct timeval* tv1, struct timeval *tv2 )
{
    return (double)(tv2->tv_sec-tv1->tv_sec) + 0.000001 * (double)(tv2->tv_usec-tv1->tv_usec) ;
}

void init_genome_index( const char* prefix )
{
    struct timeval tv, tv1 ;
    gettimeofday( &tv, 0 ) ;

    char *str = (char*)calloc(strlen(prefix) + 10, 1);
    bwase_initialize();
    init_revcom1();

    fprintf(stderr, "[init_genome_index] loading index... ");
    strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
    strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
    strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);
    strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1]);
    free(str);

    // if (!(gap_opt.mode & BWA_MODE_COMPREAD)) {  // in color space; initialize ntpac
    //	pe_opt->type = BWA_PET_SOLID;
    //	ntbns = bwa_open_nt(prefix);
    // }	
    pac = bwt_restore_pac( bns ) ;
    gettimeofday( &tv1, 0 ) ;
    fprintf(stderr, "%.2f sec\n", tdiff(&tv, &tv1));
}

void pair_aln(bam_pair_t *p)
{
    switch (p->kind) {
        case eof_marker: break ;
        case singleton: aln_singleton(p) ; break ;
        case proper_pair: aln_pair(p) ; break ;
        default: fprintf(stderr, "[%s] got impossible type of input records, bailing out.\n", __func__);
                 exit(1);
    }
}

void pair_print_bam(BGZF *output, bam_pair_t *p)
{
    // If "only aligned" is requested and at least one mate is unmapped,
    // we skip the output.  (Conserves time and space where we don't
    // need to pass through the unaligned sequences.)
    if( only_aligned && p->kind > 0 && (p->first.core.flag  & SAM_FSU) ) return ;
    if( only_aligned && p->kind > 1 && (p->second.core.flag & SAM_FSU) ) return ;

    if(p->kind > 0) bwa_print_bam1(output, &p->first);
    if(p->kind > 1) bwa_print_bam1(output, &p->second);
}

// This is the simple sequential loop:  We run it if exactly one thread
// and no listening port is requested.  Doesn't require 0MQ.  As opposed
// to "classic" BWA, we stop that blockwise nonsense and get exactly one
// singleton or pair per iteration.
void sequential_loop( BGZF* output )
{
    struct timeval tv0, tv1 ;
    gettimeofday( &tv0, 0 ) ;

    int tot_seqs = 0, rc=0 ;
    bam_pair_t raw ;
    for(;;) {
        rc=read_bam_pair(ks, &raw, broken_input);
        if(rc<0) {
            fprintf( stderr, "[run_reader_thread] error reading input BAM %s\n",
                    rc==-2 ? "(lone mate)" : "" ) ;
            exit(1);
        }
        if(rc==0) break ;

        tot_seqs += raw.kind ;
        if(tot_seqs%0x20000 < raw.kind) {
            gettimeofday( &tv1, 0 ) ;
            fprintf(stderr, "[sequential_loop] %d sequences processed in %.2f sec\n",
                    tot_seqs, tdiff( &tv0, &tv1 ));
        }
        pair_aln(&raw);
        pair_print_bam(output,&raw);
        bam_destroy_pair(&raw);
    }
    gettimeofday( &tv1, 0 ) ;
    fprintf(stderr, "[sequential_loop] %d sequences processed in %.2f sec\n",
            tot_seqs, tdiff( &tv0, &tv1 ));
    fprintf(stderr, "[sequential_loop] finished cleanly, shutting down\n");
}

void set_sockopts( void *socket )
{
    zcheck( zmq_setsockopt( socket, ZMQ_HWM, &the_hwm, sizeof(uint64_t) ),
            1, "zmq_setsockopt failed" ) ;
}
void zmq_xclose( void *socket )
{
    int the_linger = 0;
    zmq_setsockopt( socket, ZMQ_LINGER, &the_linger, sizeof(int) ) ;
    zmq_close(socket);
}


void *run_config_service( void *arg )
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

        zterm( zmq_recv(socket, &m, 0), 1, "zmq_recv failed" ) break;
        
        fprintf( stderr, "[run_config_service] Received hello: %.*s.\n", 
                (int)zmq_msg_size(&m), (char*)zmq_msg_data(&m) ) ;
        zmq_msg_close(&m) ;

        zmq_msg_init_size(&m, sizeof(gap_opt_t) + sizeof(pe_opt_t) + strlen(prefix));
        memcpy(zmq_msg_data(&m),                                       gap_opt, sizeof(gap_opt_t));
        memcpy(zmq_msg_data(&m) + sizeof(gap_opt_t),                    pe_opt, sizeof(pe_opt_t));
        memcpy(zmq_msg_data(&m) + sizeof(gap_opt_t) + sizeof(pe_opt_t), prefix, strlen(prefix));
        zterm( zmq_send(socket, &m, 0), 1, "zmq_send failed" ) break;
        zmq_msg_close(&m) ;
    }

    zmq_xclose( socket ) ;
    return 0;
}

void msg_init_from_pair(zmq_msg_t *m, bam_pair_t *p)
{
    int rc, len = (char*)&p->first - (char*)p ;
    if( p->kind > 0 ) 
        len += sizeof(bam1_core_t) + p->first.data_len + sizeof(int)
            + (p->n_aln1 ? sizeof(bwt_aln1_t) * (p->n_aln1-1) : 0) ;

    if( p->kind > 1 )
        len += sizeof(bam1_core_t) + p->second.data_len + sizeof(int)
            + (p->n_aln2 ? sizeof(bwt_aln1_t) * (p->n_aln2-1) : 0) ;

    rc = zmq_msg_init_size(m, len) ;
    xassert (rc == 0, "fail to init message");

    char *q = zmq_msg_data(m);
    memcpy(q, p, (char*)&p->first - (char*)p);
    q += (char*)&p->first - (char*)p;
    if( p->kind > 0 ) {
        memcpy(q, &p->first.core, sizeof(bam1_core_t)) ;
        q += sizeof(bam1_core_t) ;
        memcpy(q, &p->first.data_len, sizeof(int)) ;
        q += sizeof(int) ;
        memcpy(q, p->first.data, p->first.data_len);
        q += p->first.data_len;
        if(p->n_aln1) {
            memcpy(q, p->aln1, sizeof(bwt_aln1_t) * (p->n_aln1-1)) ;
            q += sizeof(bwt_aln1_t) * (p->n_aln1-1) ;
        }
    }
    if( p->kind > 1 ) {
        memcpy(q, &p->second.core, sizeof(bam1_core_t)) ;
        q += sizeof(bam1_core_t) ;
        memcpy(q, &p->second.data_len, sizeof(int)) ;
        q += sizeof(int) ;
        memcpy(q, p->second.data, p->second.data_len);
        q += p->second.data_len;
        if(p->n_aln2) {
            memcpy(q, p->aln2, sizeof(bwt_aln1_t) * (p->n_aln2-1)) ;
            q += sizeof(bwt_aln1_t) * (p->n_aln2-1) ;
        }
    }
    xassert( q == zmq_msg_data(m) + zmq_msg_size(m), "error encoding message" ) ;
}

void pair_init_from_msg(bam_pair_t *p, zmq_msg_t *m)
{
    char *q = zmq_msg_data(m);
    memcpy(p, q, (char*)&p->first - (char*)p);
    q += (char*)&p->first - (char*)p;
    if( p->kind > 0 ) {
        memcpy(&p->first.core, q, sizeof(bam1_core_t)) ;
        q += sizeof(bam1_core_t) ;
        memcpy(&p->first.data_len, q, sizeof(int)) ;
        q += sizeof(int) ;
        p->first.data = malloc( p->first.data_len );
        p->first.m_data = p->first.data_len ;
        memcpy(p->first.data, q, p->first.data_len);
        q += p->first.data_len;
        if(p->n_aln1) {
            p->aln1 = malloc(sizeof(bwt_aln1_t) * (p->n_aln1-1)) ;
            memcpy(p->aln1, q, sizeof(bwt_aln1_t) * (p->n_aln1-1)) ;
            q += sizeof(bwt_aln1_t) * (p->n_aln1-1) ;
        }
    }
    if( p->kind > 1 ) {
        memcpy(&p->second.core, q, sizeof(bam1_core_t)) ;
        q += sizeof(bam1_core_t) ;
        memcpy(&p->second.data_len, q, sizeof(int)) ;
        q += sizeof(int) ;
        p->second.data = malloc( p->second.data_len );
        p->second.m_data = p->second.data_len ;
        memcpy(p->second.data, q, p->second.data_len);
        q += p->second.data_len;
        if(p->n_aln2) {
            p->aln2 = malloc(sizeof(bwt_aln1_t) * (p->n_aln2-1)) ;
            memcpy(p->aln2, q, sizeof(bwt_aln1_t) * (p->n_aln2-1)) ;
            q += sizeof(bwt_aln1_t) * (p->n_aln2-1) ;
        }
    }
    xassert( q == zmq_msg_data(m) + zmq_msg_size(m), "error decoding message" ) ;
}

void *run_reader_thread( void *socket )
{
    int recno = 0, rc=0;
    bam_pair_t raw;
    zmq_msg_t msg;
    bam_init_pair(&raw);
    while (!s_interrupted) {
        rc=read_bam_pair(ks, &raw, broken_input);
        if(rc<0) {
            fprintf( stderr, "[run_reader_thread] error reading input BAM %s\n",
                    rc==-2 ? "(lone mate)" : "" ) ;
            exit(1);
        }
        if(rc==0) break ;

        raw.recno = recno++; 
        if(!(recno%CHUNKSIZE)) fprintf( stderr, "[run_reader_thread] %d records read so far.\n", recno );
        msg_init_from_pair( &msg, &raw );
        bam_destroy_pair(&raw);
        zterm( zmq_send(socket, &msg, 0), 1, "zmq_send failed" ) break;
        zmq_msg_close(&msg);
    }

    bam_init_pair(&raw);
    raw.kind = eof_marker;
    raw.recno = recno++;
    msg_init_from_pair( &msg, &raw );
    bam_destroy_pair(&raw);
    zterm( zmq_send(socket, &msg, 0), 0, "zmq_send failed" );
    fprintf( stderr, "[run_reader_thread] finished, %d records in total.\n", recno );
    zmq_msg_close(&msg);
    zmq_close(socket);
    return 0;
}
 
void *run_output_thread( BGZF *output, void* socket )
{
    struct timeval tv0, tv1 ;
    gettimeofday( &tv0, 0 ) ;
    int lastrn = 0 ;
    double rate = -1 ;
    while (!s_interrupted) 
    {
        zmq_msg_t msg;
        bam_pair_t pair;
        zmq_msg_init(&msg);
        zterm( zmq_recv(socket, &msg, 0), 1, "zmq_recv failed" ) break;
        pair_init_from_msg(&pair, &msg);
        zmq_msg_close(&msg);
        if(!pair.kind) break ;
        if(!pair.kind) fprintf( stderr, "[run_output_thread] %d records received in total.\n", pair.recno );
        else if(!(pair.recno%0x100)) {
            gettimeofday( &tv1, 0 );
            double sec = tdiff( &tv0, &tv1 ) ;
            if( sec >= 10 ) {
                if( rate < 0 ) rate = (pair.recno-lastrn) / (1000*sec) ;
                else rate = ((pair.recno-lastrn) / (1000*sec) + 3*rate) * 0.25 ;
                fprintf( stderr, "[run_output_thread] %d records received in %0.2fs, rate = %.1f kHz.\n",
                        pair.recno-lastrn, sec, rate ) ;
                lastrn=pair.recno;
                memmove( &tv0, &tv1, sizeof(struct timeval) ) ;
            }
        }
        pair_print_bam(output, &pair);
        bam_destroy_pair(&pair);
    }
    fprintf(stderr, "[run_output_thread] finished cleanly, shutting down\n");
    zmq_close(socket);
    return 0;
}

// Hm.  Should be rather easy:  read message, align, write message.
// Both sockets are bound, and always to the inproc endpoints.
//
// We bind to inproc://work_{in,out}.  When running standalone, the
// other end is the I/O multiplexor, when running as worker process, the
// other end is a streaming device.
void *run_worker_thread( void *arg )
{
    zmq_msg_t msg;
    bam_pair_t pair;
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
        if( 0>zmq_recv(incoming, &msg, 0) ) {
            // clean exit if we got terminated here
            if( zmq_errno() == ETERM || zmq_errno() == EINTR ) break ;
            fprintf( stderr, "zmq_recv failed: %s\n", zmq_strerror( zmq_errno() ) ) ;
            exit(1) ;
        }
        pair_init_from_msg(&pair, &msg);
        // fprintf( stderr, "[run_worker_thread] got record %d.\n", pair.recno ) ;
        zmq_msg_close(&msg);
        pair_aln(&pair);
        msg_init_from_pair(&msg,&pair);
        bam_destroy_pair(&pair);
        zterm( zmq_send(outgoing, &msg, 0), 1, "zmq_send failed" ) break;
        zmq_msg_close(&msg);
    } while(pair.kind) ;
    fprintf( stderr, "[run_worker_thread] exiting.\n" ) ;
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
// isdone==1 is the end marker.  Other records are either unacknowledged
// or acknowledged work packets.
//
// FIXME: For the time being, the sliding window has a fixed size.
//        Don't know yet how to make this more dynamic.
void *run_io_multiplexor( void* arg )
{
    // Let's see... we're multiplexing between four sockets:
    // - input from BAM         (incoming_bam)
    // - output to BAM          (outgoing_bam)
    // - fan-out to workers     (work_out)
    // - fan-in from workers    (work_in)
    void *incoming_bam = zmq_socket( zmq_context, ZMQ_PULL ) ;
    void *outgoing_bam = zmq_socket( zmq_context, ZMQ_PUSH ) ;

    set_sockopts( incoming_bam ) ;
    set_sockopts( outgoing_bam ) ;

    zcheck( zmq_connect( incoming_bam, "inproc://bam_in"  ), 1, "zmq_connect failed" );
    zcheck( zmq_connect( outgoing_bam, "inproc://bam_out" ), 1, "zmq_connect failed" );

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
        { socket : incoming_bam, fd : 0, events :  ZMQ_POLLIN, revents : 0 },
        { socket :     work_out, fd : 0, events :           0, revents : 0 },
        { socket :      work_in, fd : 0, events :  ZMQ_POLLIN, revents : 0 }
    } ;

    while (!s_interrupted) {
        xassert( current_done + current_undone + current_invalid == ring_size, "logic error in queueing thread" ) ;
        gettimeofday( &tv1, 0 ) ;
        if( loudness >= 2 || tdiff(&tv0,&tv1) >= 10 ) {
            fprintf( stderr, "[run_io_multiplexor] polling -- %d records received, "
                             "%d written, %d resends; %d done in queue, %d undone in queue.\n",
                     total_in, total_out, total_resends, current_done, current_undone ) ;
            memmove( &tv0, &tv1, sizeof(struct timeval) );
        }
        if( loudness >= 3 )
            fprintf( stderr, "next expected: %d, next free: %d\n", next_expected % ring_size, next_free % ring_size ) ; 

        zterm( zmq_poll( polls, 3, -1 ), 1, "zmq_poll failed" ) break;

        // next time, do not poll for input
        polls[0].events = 0 ;

        // If something is ready, what do we do?
        // - If fresh input arrives, send it to the network and put it in
        //   the ring buffer.  (If we did poll, the net was already
        //   guaranteed to be ready.)
        if( polls[0].revents & ZMQ_POLLIN) { // input available
            // the buffer should not be full (otherwise it's a logic error)
            xassert( (next_free+1) % ring_size != next_expected % ring_size, "logic error in queueing thread" ) ;
            if( loudness >= 3 ) fprintf( stderr, "got input record.\n" ) ;
            zmq_msg_t msg ;
            zmq_msg_init( &msg ) ;
            zterm( zmq_recv( incoming_bam, &msg, 0 ), 1, "zmq_recv failed" ) break;

            // input comes ordered, so the correct slot should be
            // (next_free)
            rn = *(int*)zmq_msg_data(&msg) ;
            if( next_free != rn )
            {
                bam_pair_t p;
                pair_init_from_msg( &p, &msg ) ;
                fprintf( stderr, "%d != %d\n", next_free, *(int*)zmq_msg_data(&msg) ) ;
                fprintf( stderr, "kind == %d, isdone == %d, recno == %d\n", p.kind, p.isdone, p.recno ) ;
                exit(2);
            }

            xassert( !ringbuf[next_free % ring_size].kind && !ringbuf[next_free % ring_size].isdone, "logic error in queueing thread" ) ;
            pair_init_from_msg( &ringbuf[next_free % ring_size], &msg ) ;
            current_invalid-- ;
            if( ringbuf[next_free % ring_size].kind ) {
                total_in++ ;
                if( !current_undone ) next_resend = next_free % ring_size;
                current_undone++ ;
                if( loudness >= 3 )
                    fprintf( stderr, "received recno %d (%d), will send.\n",
                             ringbuf[next_free % ring_size].recno, ringbuf[next_free % ring_size].kind ) ;
                zterm( zmq_send( work_out, &msg, 0 ), 1, "zmq_send failed" ) break;

                polls[1].events = ZMQ_POLLOUT ;
                next_free++;
                zmq_msg_close(&msg) ;
            } else {
                // How to shut down?
                // - If EOF arrives from input, we shut down input and enqueue the
                //   EOF marker (as already acknowledged).
                ringbuf[next_free % ring_size].isdone = done ;
                zmq_close( incoming_bam ) ;
                incoming_bam = 0 ;
                polls[0].socket = 0 ;
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
                    zterm( zmq_send( work_out, &msg, 0 ), 1, "zmq_send failed" ) break;

                    do next_resend = (next_resend+1) % ring_size ;
                    while( ringbuf[next_resend].isdone == done || !ringbuf[next_resend].kind ) ;
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
            zterm( zmq_recv( work_in, &msg, 0 ), 1, "zmq_recv failed" ) break;
            
            rn = *(int*)zmq_msg_data(&msg) ;
            if( loudness >= 3 ) fprintf( stderr, "[run_io_multiplexor] received record %d.\n", rn ) ;
            if( rn < next_expected ) {
                if( loudness >= 1 ) fprintf( stderr, "[run_io_multiplexor] this is old shit: %d.\n", rn ) ;
                total_dups++ ;
            }
            else if( rn >= next_free ) {
                if( loudness >= 1 ) fprintf( stderr, "[run_io_multiplexor] this comes from the future: %d.\n", rn ) ;
            }
            else if( ringbuf[rn % ring_size].isdone ) {
                if( loudness >= 1 ) fprintf( stderr, "[run_io_multiplexor] this is a duplicate: %d.\n", rn ) ;
                total_dups++ ;
            }
            else {
                if( loudness >= 3 ) fprintf( stderr, "this is not a duplicate, dealing with it\n" ) ;
                bam_destroy_pair( &ringbuf[ rn % ring_size ] ) ;
                pair_init_from_msg( &ringbuf[ rn % ring_size ], &msg ) ;
                ringbuf[ rn % ring_size ].isdone = done ;
                xassert( ringbuf[ rn % ring_size ].kind, "wrong record type in queue" ) ;
                zmq_msg_close(&msg ) ;
                current_undone-- ;
                current_done++ ;

received_one:
                if( current_undone && next_resend == rn % ring_size ) 
                    do next_resend = (next_resend+1) % ring_size ;
                    while( ringbuf[next_resend].isdone == done || !ringbuf[next_resend].kind ) ;

                int foo = 0 ;
                while( next_expected % ring_size != next_free % ring_size && ringbuf[next_expected % ring_size].isdone == done )
                {
                    int mykind = ringbuf[next_expected % ring_size].kind ;
                    if( loudness >= 3 ) fprintf( stderr, "next expected: %d, next free: %d\n",
                            next_expected % ring_size, next_free % ring_size ) ;

                    // make message from pair, then invalidate the slot
                    // so it doesn't get sent again.
                    msg_init_from_pair( &msg, &ringbuf[next_expected % ring_size] ) ;
                    bam_destroy_pair( &ringbuf[next_expected % ring_size] ) ;
                    ringbuf[next_expected % ring_size].kind = 0 ; 
                    ringbuf[next_expected % ring_size].isdone = 0 ;

                    zterm( zmq_send( outgoing_bam, &msg, 0 ), 1, "zmq_send failed" ) goto break2;
                    zmq_msg_close(&msg) ;
                    total_out++;
                    current_done--;
                    current_invalid++;

                    // How to shut down?
                    // - If EOF is ready to be sent out, shut everything
                    //   down and send the EOF.
                    if( !mykind ) goto break2 ;
                    next_expected++ ;
                    if( loudness >= 4 ) fprintf( stderr, "spinning(1)\n" ) ;
                    xassert( foo++ < ring_size, "logic error in queueing thread" ) ;
                }

                if( current_undone ) {
                    int foo = 0 ;
                    while( ringbuf[next_resend].isdone == done ) {
                        next_resend = (next_resend+1) % ring_size ;
                        xassert( foo++ < ring_size, "logic error in queueing thread" ) ;
                    }
                }
            }
        }
    }
break2:
    fprintf( stderr, "[run_io_multiplexor] finished: %d records in total, %d resends, %d dups.\n",
            total_out, total_resends, total_dups ) ;
    if( listen_port ) {
        zmq_xclose( work_out ) ;
        zmq_xclose( work_in ) ;
    }
    else {
        zmq_close( work_out ) ;
        zmq_close( work_in ) ;
    }
    zmq_xclose( incoming_bam ) ;
    zmq_xclose( outgoing_bam ) ;
    return 0 ;
}

// FIXME:  Don't know what to do about that weird trimming stuff.
//
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

void bwa_bam2bam_core( const char *prefix, BGZF *output )
{
	// cnt_chg=0;
	// isize_info_t last_ii = {0} ;

    // Initialize genome index.  Needed if and only if we're going to
    // have a worker thread.
    if( gap_opt->n_threads ) init_genome_index( prefix ) ;

	g_hash = kh_init(64);
	// last_ii.avg = -1.0;
    ii.low = ii.high = 0; ii.avg = ii.std = -1.0;
    srand48(bns->seed);

    // Can we get away with a simple loop?
    if( !gap_opt->n_threads && !listen_port ) {
        fprintf( stderr, "[bwa_bam2bam_core] No threads and no listening port specified, nothing to do.\n" ) ;
        exit(1);
    } else if( gap_opt->n_threads == 1 && !listen_port ) {
        sequential_loop( output ) ;
    } else {
        pthread_t config_service_tid ;
        pthread_t reader_tid ;
        pthread_t mux_tid ;
        pthread_t *worker_tid ;

        // Okay, we need 0MQ.
        s_catch_signals();
        zmq_context = zmq_init(1) ;
        xassert( zmq_context, "0MQ init failed" ) ;

        void *bam_out = zmq_socket( zmq_context, ZMQ_PULL ) ;
        xassert(bam_out, "couldn't create socket");
        set_sockopts(bam_out);
        zcheck( zmq_bind( bam_out, "inproc://bam_out" ), 1, "zmq_bind failed" );

        void *broadcast = zmq_socket( zmq_context, ZMQ_PUB ) ;
        xassert(broadcast, "couldn't create socket");

        // Are we going to be network aware?  Then we need a
        // configuration service.  And a broadcast channel.
        if( listen_port ) {
            char addr[100] ;
            pthread_create( &config_service_tid, 0, run_config_service, 0 ) ;
            snprintf( addr, 100, "tcp://*:%d", listen_port+3 ) ;
            zcheck( zmq_bind( broadcast, addr ), 1, "zmq_bind failed" );
            set_sockopts( broadcast ) ;
        }

        // Start reading in the background...
        void *bam_in = zmq_socket( zmq_context, ZMQ_PUSH ) ;
        xassert(bam_in, "couldn't create socket");
        set_sockopts(bam_in);

        zcheck( zmq_bind( bam_in, "inproc://bam_in" ), 1, "zmq_bind failed" );
        pthread_create( &reader_tid, 0, run_reader_thread, bam_in );

        work_out = zmq_socket( zmq_context, ZMQ_PUSH ) ;
        work_in  = zmq_socket( zmq_context, ZMQ_PULL ) ;
        
        set_sockopts(work_out);
        set_sockopts(work_in);

        zcheck( zmq_bind( work_out, "inproc://work_out" ), 1, "zmq_bind failed" );
        zcheck( zmq_bind( work_in, "inproc://work_in" ), 1, "zmq_bind" );
        
        if( listen_port ) {
            char addr[100] ;
            snprintf( addr, 100, "tcp://*:%d", listen_port+1 ) ;
            zcheck( zmq_bind( work_out, addr ), 1, "zmq_bind failed" );
            snprintf( addr, 100, "tcp://*:%d", listen_port+2 ) ;
            zcheck( zmq_bind( work_in, addr ), 1, "zmq_bind failed" );
        }

        // Are threads requested?  Then fire them up.
        if( gap_opt->n_threads ) {
            int n;
            worker_tid = calloc( gap_opt->n_threads, sizeof(pthread_t) ) ;
            xassert( worker_tid, "out of momory" ) ;

            for(n=0;n!=gap_opt->n_threads;++n)
                // Hmm, arguments?
                pthread_create(worker_tid+n, 0, run_worker_thread, 0 ) ;
        }

        // ...multiplex records between I/O and workers...
        pthread_create( &mux_tid, 0, run_io_multiplexor, 0 );

        // ...and shift to output in the foreground.
        run_output_thread( output, bam_out );

        // when we get here, we got the last record and an EOF marker,
        // so we can safely tear down everything

        // Done, tell everyone and get rid of 0MQ.
        zmq_msg_t finmsg ;
        zmq_msg_init_data( &finmsg, "go away", 7, 0, 0 ) ;
        fprintf( stderr, "[bwa_bam2bam_core] sending termination signal\n" ) ;
        zcheck( zmq_send( broadcast, &finmsg, 0 ), 1, "zmq_send failed" ) ;
        fprintf( stderr, "[bwa_bam2bam_core] done, waiting for clean shutdown\n" ) ;
        zmq_xclose( broadcast ) ;
        zmq_term( zmq_context ) ;
    }

	if (pac) bwt_destroy_pac(pac,bns);
	if (ntbns) bns_destroy(ntbns);
	if (bns) bns_destroy(bns);
	if (bwt[0]) bwt_destroy(bwt[0]);
    if (bwt[1]) bwt_destroy(bwt[1]);

	for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
		if (kh_exist(g_hash, iter))
            free(kh_val(g_hash, iter).a);
	kh_destroy(64, g_hash);
}

int bwa_bam_to_bam( int argc, char *argv[] )
{
	int c, opte = -1;
    char *saif[3] = {0,0,0} ;
    char *ofile = 0;

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
	ks = bwa_bam_open(argv[optind], 7, saif, gap_opt, &hdr);
    bns = bns_restore(prefix);

    BGZF* output = ofile ? bgzf_open(ofile, "w2") : bgzf_fdopen(1, "w2") ;
    if( 0>bwa_print_bam_header(output, bns, hdr->text, argc, argv) ) {
        fprintf( stderr, "[bwa_bam2bam_core] Error writing BAM header.\n" ) ;
        exit(1);
    }
    bam_header_destroy(hdr);

    bwa_bam2bam_core(prefix, output);
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

void bwa_worker_core( int nthreads, char* host, int port ) 
{
    time_t start_time = time(0) ;
    s_catch_signals();
	// int tot_seqs=0; // cnt_chg=0;
	// isize_info_t last_ii = {0} ;

    // Initialize genome index.
    bns = bns_restore(prefix);
    init_genome_index( prefix ) ;

	g_hash = kh_init(64);
	// last_ii.avg = -1.0;
    ii.low = ii.high = 0; ii.avg = ii.std = -1.0;

    srand48(bns->seed);

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

    pthread_t *worker_tid ;
    int n;
    worker_tid = calloc( nthreads, sizeof(pthread_t) ) ;
    xassert( worker_tid, "out of memory" ) ;

    for(n=0;n!=nthreads;++n)
        pthread_create(worker_tid+n, 0, run_worker_thread, 0 ) ;

    pthread_t fan_out ;
    pthread_create( &fan_out, 0, run_device, socks ) ;

    // Now run the streamer for the fan-in.  We just ping-pong between
    // receive and send, but if either times out, we simply terminate.
    zmq_pollitem_t pitems[] = {
        { socks[2], 0, ZMQ_POLLIN, 0 },
        { socks[4], 0, ZMQ_POLLIN, 0 },
        { socks[3], 0, ZMQ_POLLOUT, 0 } } ;

    int npoll=0;
    zmq_msg_t m ;
    zmq_msg_init(&m) ;
    while (!s_interrupted && time(0) - start_time < 90*60) {
        npoll = zmq_poll( pitems, 2, 1E6 * timeout ) ;
        zterm( npoll, 1, "zmq_poll failed" ) break;

        if( npoll == 0 || (pitems[1].revents & ZMQ_POLLIN) ) break ;
        zterm( zmq_recv( socks[2], &m, 0 ), 1, "zmq_recv failed" ) break;

        npoll = zmq_poll( pitems+1, 2, 1E6 * timeout ) ;
        zterm( npoll, 1, "zmq_poll failed" ) break;

        if( npoll == 0 || (pitems[1].revents & ZMQ_POLLIN) ) break ;
        zterm( zmq_send( socks[3], &m, 0 ), 1, "zmq_send failed" ) break;
    }

    if( npoll==0 )
        fprintf( stderr, "[bwa_worker_core] No work delivered in %ds.  Terminating.\n", timeout ) ;
    else if( pitems[1].revents & ZMQ_POLLIN ) {
        zmq_recv( socks[4], &m, 0 );
        fprintf( stderr, "[bwa_worker_core] Received termination signal: \"%.*s\".\n",
                (int)zmq_msg_size(&m), (char*)zmq_msg_data(&m) ) ;
    } else if( s_interrupted ) 
        fprintf( stderr, "[bwa_worker_core] Received interrupt signal.\n" ) ;
    else 
        fprintf( stderr, "[bwa_worker_core] I've been going for %d minutes, I'm tired now.\n",
                (time(0) - start_time) / 60 ) ;

    zmq_msg_close(&m);
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

    for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
        if (kh_exist(g_hash, iter))
            free(kh_val(g_hash, iter).a);
    kh_destroy(64, g_hash);
}

int bwa_worker( int argc, char *argv[] )
{
    if( argc != 4 ) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa worker <nthreads> <host> <port>\n\n");
        return 1;
    }

    int nthreads = atoi( argv[1] ) ;
    char *host = argv[2] ;
    int port = atoi( argv[3] ) ;

    zmq_context = zmq_init(1) ;
    xassert( zmq_context, "0MQ init failed" ) ;

    void *conf_sock = zmq_socket( zmq_context, ZMQ_REQ );
    xassert( conf_sock, "couldn't create socket" );

    char addr[100] ;
    snprintf( addr, 100, "tcp://%s:%d", host, port ) ;
    zcheck( zmq_connect( conf_sock, addr ), 1, "zmq_connect failed" );

    zmq_msg_t m;
    zmq_msg_init_size(&m,20);
    gethostname( zmq_msg_data(&m), 20 ) ;
    zcheck( zmq_send(conf_sock, &m, 0), 1, "zmq_send failed" );
    zmq_msg_close(&m);
    zmq_msg_init(&m);
    zcheck( zmq_recv(conf_sock, &m, 0), 1, "zmq_recv failed" );
    zmq_xclose(conf_sock);

    int plen = zmq_msg_size(&m) - sizeof(gap_opt_t) - sizeof(pe_opt_t) ;
    zcheck( plen, 1, "received truncated configuration information" ) ;

	gap_opt = gap_init_opt();
	pe_opt = bwa_init_pe_opt();
    prefix = malloc( 1+plen ) ;

    memcpy( gap_opt, zmq_msg_data(&m), sizeof(gap_opt_t) );
    memcpy( pe_opt, zmq_msg_data(&m) + sizeof(gap_opt_t), sizeof(pe_opt_t) );
    memcpy( prefix, zmq_msg_data(&m) + sizeof(gap_opt_t) + sizeof(pe_opt_t), plen ) ;
    prefix[plen] = 0 ;

    zmq_msg_close(&m);
    bwa_worker_core( nthreads, host, port ) ;

    free(prefix);
	free(pe_opt);
	free(gap_opt);
    return 0;
}

