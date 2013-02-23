#ifndef BWA_H_
#define BWA_H_

#include <stdint.h>

typedef struct {
	int l_seq;
	char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

#ifdef __cplusplus
extern "C" {
#endif

	bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);

	uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar);

#ifdef __cplusplus
}
#endif

#endif
