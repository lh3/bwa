#ifndef BWA_H_
#define BWA_H_

#include <stdint.h>

#define BWA_DEF_MAX_SCORE 2048
#define BWA_MAX_QUERY_LEN 1024

// BWA index
struct bwa_idx_t;
typedef struct bwa_idx_t bwa_idx_t;

// Buffer for BWA alignment
struct bwa_buf_t;
typedef struct bwa_buf_t bwa_buf_t;

// BWA alignment options
typedef struct {
	int s_gapo, s_gape;               // gap open and extension penalties; the mismatch penalty is fixed at 3
	int max_diff, max_gapo, max_gape; // max differences (-1 to use fnr for length-adjusted max diff), gap opens and gap extensions
	int seed_len, max_seed_diff;      // seed length and max differences allowed in the seed
	float fnr;                        // parameter for automatic length-adjusted max differences
} bwa_opt_t;

// default BWA alignment options
extern bwa_opt_t bwa_def_opt; // = { 11, 4, -1, 1, 6, 32, 2, 0.04 }

// an interval hit in the SA coordinate; basic unit in .sai files
typedef struct {
	uint32_t n_mm:16, n_gapo:8, n_gape:8;
	int score;
	uint64_t k, l; // [k,l] is the SA interval; each interval has l-k+1 hits
} bwa_sai1_t;

// all interval hits in the SA coordinate
typedef struct {
	int n; // number of interval hits
	bwa_sai1_t *sai;
} bwa_sai_t;

// an alignment
typedef struct {
	uint32_t n_n:8, n_gap:12, n_mm:12; // number of ambiguous bases, gaps and mismatches in the alignment
	int32_t ref_id;                    // referece sequence index (the first seq is indexed by 0)
	uint32_t offset;                   // coordinate on the reference; zero-based
	uint32_t n_cigar:16, flag:16;      // number of CIGAR operations; SAM flag
	uint32_t *cigar;                   // CIGAR in the BAM 28+4 encoding; having n_cigar operations
} bwa_aln_t;

typedef struct {
	int mapQs, mapQ, c1, c2;
	uint64_t sa;
	bwa_sai1_t *which;
	bwa_sai_t sai;
	bwa_aln_t one;
} bwa_one_t;

typedef struct {
	double avg, std, ap_prior;
	uint64_t low, high, high_bayesian;
} bwa_pestat_t;

#ifdef __cplusplus
extern "C" {
#endif

	// load a BWA index
	bwa_idx_t *bwa_idx_load(const char *prefix);
	void bwa_idx_destroy(bwa_idx_t *p);

	// allocate a BWA alignment buffer; if unsure, set opt to &bwa_def_opt and max_score to BWA_DEF_MAX_SCORE
	bwa_buf_t *bwa_buf_init(const bwa_opt_t *opt, int max_score);
	void bwa_buf_destroy(bwa_buf_t *p);

	/**
	 * Find all the SA intervals
	 *
	 * @param idx    BWA index; multiple threads can share the same index
	 * @param buf    BWA alignment buffer; each thread should have its own buffer
	 * @param seq    NULL terminated C string, consisting of A/C/G/T/N only
	 *
	 * @return       SA intervals seq is matched to
	 */
	bwa_sai_t bwa_sai(const bwa_idx_t *idx, bwa_buf_t *buf, const char *seq);

	/**
	 * Construct an alignment in the base-pair coordinate
	 *
	 * @param idx     BWA index
	 * @param buf     BWA alignment buffer
	 * @param seq     NULL terinated C string
	 * @param sa      Suffix array value
	 * @param n_gaps  Number of gaps (typically equal to bwa_sai1_t::n_gapo + bwa_sai1_t::n_gape
	 *
	 * @return        An alignment
	 */
	void bwa_sa2aln(const bwa_idx_t *idx, bwa_buf_t *buf, const char *seq, uint64_t sa, int n_gaps, bwa_aln_t *aln);

	bwa_one_t *bwa_se(const bwa_idx_t *idx, bwa_buf_t *buf, const char *seq, int gen_cigar);

	void bwa_one_destroy(bwa_one_t *one);

#ifdef __cplusplus
}
#endif

#endif
