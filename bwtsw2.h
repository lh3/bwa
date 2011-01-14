#ifndef LH3_BWTSW2_H
#define LH3_BWTSW2_H

#include <stdint.h>
#include "bntseq.h"
#include "bwt_lite.h"
#include "bwt.h"

typedef struct {
	int a, b, q, r, t, qr, bw;
	int z, is, t_seeds, hard_clip;
	float yita, mask_level, coef;
	int n_threads, chunk_size;
} bsw2opt_t;

typedef struct {
	uint32_t k, l, flag:18, n_seeds:14;
	int len, G, G2;
	int beg, end;
} bsw2hit_t;

typedef struct {
	int n, max;
	bsw2hit_t *hits;
	int *n_cigar;
	uint32_t **cigar;
} bwtsw2_t;

typedef struct {
	void *stack;
	int max_l;
	uint8_t *aln_mem;
} bsw2global_t;

#ifdef __cplusplus
extern "C" {
#endif

	bsw2opt_t *bsw2_init_opt();
	bwtsw2_t **bsw2_core(const bsw2opt_t *opt, const bwtl_t *target, const bwt_t *query, bsw2global_t *pool);
	void bsw2_aln(const bsw2opt_t *opt, const bntseq_t *bns, bwt_t * const target[2], const char *fn);
	void bsw2_destroy(bwtsw2_t *b);

	bsw2global_t *bsw2_global_init();
	void bsw2_global_destroy(bsw2global_t *_pool);

#ifdef __cplusplus
}
#endif

#endif
