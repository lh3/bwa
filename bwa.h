#ifndef BWA_H_
#define BWA_H_

#include <stdint.h>
#include "bntseq.h"
#include "bwt.h"

#define BWA_IDX_BWT 0x1
#define BWA_IDX_BNS 0x2
#define BWA_IDX_PAC 0x4
#define BWA_IDX_ALL 0x7

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

typedef struct {
	int l_seq;
	char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

extern int bwa_verbose;

#ifdef __cplusplus
extern "C" {
#endif

	bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);

	uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar);

	char *bwa_idx_infer_prefix(const char *hint);
	bwt_t *bwa_idx_load_bwt(const char *hint);
	bwaidx_t *bwa_idx_load(const char *hint, int which);
	void bwa_idx_destroy(bwaidx_t *idx);

#ifdef __cplusplus
}
#endif

#endif
