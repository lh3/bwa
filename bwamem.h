#ifndef BWAMEM_H_
#define BWAMEM_H_

#include "bwt.h"

struct __smem_i;
typedef struct __smem_i smem_i;

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
} memseed_t;

typedef struct {
	int a, b, q, r, w;
	int min_seed_len, max_occ, max_chain_gap;
} memopt_t;

typedef struct {
	int n, m;
	int64_t pos;
	memseed_t *seeds;
} memchain1_t;

typedef struct {
	int n, m;
	memchain1_t *chains;
} memchain_t;

#ifdef __cplusplus
extern "C" {
#endif

smem_i *smem_itr_init(const bwt_t *bwt);
void smem_itr_destroy(smem_i *itr);
void smem_set_query(smem_i *itr, int len, const uint8_t *query);
const bwtintv_v *smem_next(smem_i *itr, int split_len);

memopt_t *mem_opt_init(void);

memchain_t mem_chain(const memopt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq);

#ifdef __cplusplus
}
#endif

#endif
