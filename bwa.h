#ifndef BWA_H_
#define BWA_H_

#include <stdint.h>

#define BWA_DEF_MAX_SCORE 2048
#define BWA_MAX_QUERY_LEN 1024

struct bwa_idx_t;
typedef struct bwa_idx_t bwa_idx_t;

struct bwa_aux_t;
typedef struct bwa_aux_t bwa_aux_t;

typedef struct {
	int s_gapo, s_gape; // the mismatch penalty is fixed at 3
	int max_diff, max_gapo, max_gape;
	int seed_len, max_seed_diff;
	float fnr;
} bwa_opt_t;

typedef struct {
	uint32_t n_mm:16, n_gapo:8, n_gape:8;
	int score;
	uint64_t k, l;
} bwa_alnpre_t;

typedef struct {
	uint32_t n_cigar:15, gap:8, mm:8, strand:1;
	uint32_t ref_id;
	uint64_t offset;
	uint32_t *cigar;
} bwa_aln_t;

extern bwa_opt_t bwa_def_opt;

#ifdef __cplusplus
extern "C" {
#endif

	bwa_idx_t *bwa_idx_load(const char *prefix);
	void bwa_idx_destroy(bwa_idx_t *p);
	bwa_aux_t *bwa_aux_init(const bwa_opt_t *opt, int max_score);
	void bwa_aux_destroy(bwa_aux_t *p);
	bwa_alnpre_t *bwa_aln_pre(const bwa_idx_t *idx, bwa_aux_t *aux, const char *seq, int *n_aln);

#ifdef __cplusplus
}
#endif

#endif
