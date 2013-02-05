#ifndef BWAMEM_H_
#define BWAMEM_H_

#include "bwt.h"

struct __smem_i;
typedef struct __smem_i smem_i;

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
} mem_seed_t;

typedef struct {
	int a, b, q, r, w;
	int min_seed_len, max_occ, max_chain_gap;
	int8_t mat[25]; // scoring matrix; mat[0] == 0 if unset
	float mask_level, chain_drop_ratio;
} mem_opt_t;

typedef struct {
	int n, m;
	int64_t pos;
	mem_seed_t *seeds;
} mem_chain1_t;

typedef struct {
	int n, m;
	mem_chain1_t *chains;
} mem_chain_t;

typedef struct {
	int64_t pos, rb, re;
	int n_cigar, len, score, qb, qe, is_all;
	uint32_t *cigar;
} mem_aln_t;

#ifdef __cplusplus
extern "C" {
#endif

smem_i *smem_itr_init(const bwt_t *bwt);
void smem_itr_destroy(smem_i *itr);
void smem_set_query(smem_i *itr, int len, const uint8_t *query);
const bwtintv_v *smem_next(smem_i *itr, int split_len);

mem_opt_t *mem_opt_init(void);
void mem_fill_scmat(int a, int b, int8_t mat[25]);

mem_chain_t mem_chain(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq);
void mem_chain_flt(const mem_opt_t *opt, mem_chain_t *chn);
void mem_chain2aln(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain1_t *c, mem_aln_t *a);

#ifdef __cplusplus
}
#endif

#endif
