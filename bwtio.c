#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MMAP
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#endif

#include "bwt.h"
#include "utils.h"

void bwt_dump_bwt(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	fwrite(bwt->bwt, sizeof(bwtint_t), bwt->bwt_size, fp);
	fclose(fp);
}

void bwt_dump_sa(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	fwrite(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	fwrite(&bwt->seq_len, sizeof(bwtint_t), 1, fp);
	fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	fclose(fp);
}

#ifdef USE_MMAP
void *mmap_t(int touch, void *addr, size_t length, int prot, int flags,
        int fd, off_t offset)
{
#ifdef MAP_POPULATE
    if( touch ) flags |= MAP_POPULATE ;
#endif
    void *p = mmap(addr, length, prot, flags, fd, offset);
    if( p != (void*)(-1) && touch ) {
        char *q = p, *qe = p + length ;
        int x, ps = sysconf(_SC_PAGE_SIZE) ;
        for( ; q < qe ; q += ps ) x += *q ;
    }
    return p ;
}

ubyte_t *bwt_restore_pac( const bntseq_t *bns, int touch)
{
    ubyte_t *pac = mmap_t( touch, 0, bns->l_pac/4+1, PROT_READ, MAP_PRIVATE, fileno(bns->fp_pac), 0 ) ;
    xassert(pac != (void*)-1, "failed to mmap pac file");
    madvise(pac, bns->l_pac/4+1, MADV_WILLNEED);
    return pac;
}

void bwt_destroy_pac( ubyte_t *pac, const bntseq_t *bns )
{
    if (pac && 0>munmap( pac, bns->l_pac/4+1))
        fprintf(stderr, "Warning: munmap() of pax failed (%s).  Memory leak?\n", strerror(errno));
}

void bwt_restore_sa(const char *fn, bwt_t *bwt, int touch)
{
    int fd = open(fn, O_RDONLY);
    xassert(fd>=0, "failed to open sa file");
    off_t len = lseek(fd, 0, SEEK_END);
    xassert(len>=0, "failed to seek in sa file");
    bwtint_t *raw = mmap_t( touch, 0, len, PROT_READ, MAP_PRIVATE, fd, 0);
    xassert(raw != (void*)-1, "failed to mmap sa file");
    madvise(raw, len, MADV_WILLNEED);

	xassert(raw[0] == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	xassert(raw[6] == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");
	bwt->sa_intv = raw[5];
	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = raw+6 ;
    bwt->mmap_sa = 1;
	/* bwt->sa[0] = -1; */
	close(fd);
    xassert((bwt->n_sa+6) * sizeof(bwtint_t)==len, "SA-BWT inconsistency: length does not match.");
}

bwt_t *bwt_restore_bwt(const char *fn, int touch )
{
    int i, fd = open(fn, O_RDONLY);
    xassert(fd>=0, "failed to open bwt file");
    off_t len = lseek(fd, 0, SEEK_END);
    xassert(len>=0, "failed to seek in bwt file");
	bwt_t *bwt = (bwt_t*)calloc(1, sizeof(bwt_t));

	bwt->bwt_size = (len - sizeof(bwtint_t) * 5) >> 2;

    bwtint_t *raw = mmap_t( touch, 0, len, PROT_READ, MAP_PRIVATE, fd, 0);
    xassert(raw != (void*)-1, "failed to mmap bwt file");
    madvise(raw, len, MADV_WILLNEED);
    bwt->primary = raw[0] ;
    for(i=1;i!=5;++i) bwt->L2[i] = raw[i] ;

	bwt->bwt = (uint32_t*)(raw+5) ;
    bwt->mmap_bwt = 1;
	bwt->seq_len = bwt->L2[4];
	close(fd);
	bwt_gen_cnt_table(bwt);
	return bwt;
}

void bwt_destroy_bwt(bwt_t *bwt, bwtint_t bwt_size)
{
    if (bwt->bwt) {
        if(bwt->mmap_bwt) {
            if( 0>munmap(((bwtint_t*)bwt->bwt)-5, 5*sizeof(bwtint_t) + (bwt_size << 2)) )
                fprintf(stderr, "Warning: munmap() of bwt->bwt failed (%s).  Memory leak?\n", strerror(errno));
        }
        else free(bwt->bwt);
    }
}

void bwt_destroy_sa(bwt_t *bwt)
{
    if (bwt->sa) {
        if( bwt->mmap_sa ) {
            if (0>munmap(bwt->sa-6, (bwt->n_sa+6) * sizeof(bwtint_t)))
                fprintf(stderr, "Warning: munmap() of bwt->sa failed (%s).  Memory leak?\n", strerror(errno));
        }
        else free(bwt->sa);
    }
}

void bwt_destroy(bwt_t *bwt)
{
    if (bwt) {
        bwt_destroy_sa(bwt);
        bwt_destroy_bwt(bwt, bwt->bwt_size);
        free(bwt);
    }
}

#else

ubyte_t *bwt_restore_pac( const bntseq_t *bns )
{
    ubyte_t *pac = (ubyte_t*)malloc(bns->l_pac/4+1);
    rewind(bns->fp_pac);
    err_fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
    return pac;
}

void bwt_destroy_pac( ubyte_t *pac, const bntseq_t *bns )
{
    free(pac) ;
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = xopen(fn, "rb");
	err_fread(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	err_fread(skipped, sizeof(bwtint_t), 4, fp); // skip
	err_fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fread(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	/* bwt->sa[0] = -1; */

	err_fread(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	fclose(fp);
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = xopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	fseek(fp, 0, SEEK_SET);
	err_fread(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fread(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	err_fread(bwt->bwt, 4, bwt->bwt_size, fp);
	bwt->seq_len = bwt->L2[4];
	fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}

void bwt_destroy_bwt(bwt_t *bwt, bwtint_t bwt_size) { free(bwt->bwt); }
void bwt_destroy_sa(bwt_t *bwt) { free(bwt->sa); }


#endif
