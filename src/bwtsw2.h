#ifndef LH3_BWTSW2_H
#define LH3_BWTSW2_H

#include <stdint.h>
#include "bntseq.h"
#include "bwt_lite.h"
#include "bwt.h"

#define BSW2_FLAG_MATESW  0x100
#define BSW2_FLAG_TANDEM  0x200
#define BSW2_FLAG_MOVED   0x400
#define BSW2_FLAG_RESCUED 0x800

typedef struct {
	int skip_sw:8, cpy_cmt:8, hard_clip:16;
	int a, b, q, r, t, qr, bw, max_ins, max_chain_gap;
	int z, is, t_seeds, multi_2nd;
	float mask_level, coef;
	int n_threads, chunk_size;
} bsw2opt_t;

typedef struct {
	bwtint_t k, l;
	uint32_t flag:18, n_seeds:13, is_rev:1;
	int len, G, G2;
	int beg, end;
} bsw2hit_t;

typedef struct {
	int flag, nn, n_cigar, chr, pos, qual, mchr, mpos, pqual, isize, nm;
	uint32_t *cigar;
} bsw2aux_t;

typedef struct {
	int n, max;
	bsw2hit_t *hits;
	bsw2aux_t *aux;
} bwtsw2_t;

typedef struct {
	void *stack;
	int max_l;
	uint8_t *aln_mem;
} bsw2global_t;

typedef struct {
	int l, tid;
	char *name, *seq, *qual, *sam, *comment;
} bsw2seq1_t;

#ifdef __cplusplus
extern "C" {
#endif

	bsw2opt_t *bsw2_init_opt();
	bwtsw2_t **bsw2_core(const bntseq_t *bns, const bsw2opt_t *opt, const bwtl_t *target, const bwt_t *query, bsw2global_t *pool);
	void bsw2_aln(const bsw2opt_t *opt, const bntseq_t *bns, bwt_t * const target, const char *fn, const char *fn2);
	void bsw2_destroy(bwtsw2_t *b);

	bsw2global_t *bsw2_global_init();
	void bsw2_global_destroy(bsw2global_t *_pool);

	void bsw2_pair(const bsw2opt_t *opt, int64_t l_pac, const uint8_t *pac, int n, bsw2seq1_t *seq, bwtsw2_t **hit);

#ifdef __cplusplus
}
#endif

#endif
