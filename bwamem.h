#ifndef BWAMEM_H_
#define BWAMEM_H_

#include "bwt.h"

typedef struct {
	const bwt_t *bwt;
	const uint8_t *query;
	int start, len, min_intv;
	bwtintv_v *tmpvec[2], *matches;
} smem_i;

typedef struct {
	int a, b, q, r, w;
	int min_seed_len, max_occ;
} memopt_t;

#ifdef __cplusplus
extern "C" {
#endif

smem_i *smem_itr_init(const bwt_t *bwt);
void smem_itr_destroy(smem_i *itr);
void smem_set_query(smem_i *itr, int min_intv, int len, const uint8_t *query);
int smem_next(smem_i *itr);

memopt_t *mem_opt_init(void);

#ifdef __cplusplus
}
#endif

#endif
