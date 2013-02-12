#ifndef BWAMEM_H_
#define BWAMEM_H_

#include "bwt.h"
#include "bntseq.h"
#include "utils.h"

struct __smem_i;
typedef struct __smem_i smem_i;

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
} mem_seed_t;

typedef struct {
	int a, b, q, r, w;
	int split_width;
	int min_seed_len, max_occ, max_chain_gap;
	int n_threads, chunk_size;
	int pe_dir, is_pe;
	int is_hard; // if to use hard clip
	float mask_level, chain_drop_ratio;
	int max_ins; // maximum insert size
	int8_t mat[25]; // scoring matrix; mat[0] == 0 if unset
} mem_opt_t;

typedef struct {
	int n, m;
	int64_t pos;
	mem_seed_t *seeds;
} mem_chain_t;

typedef struct {
	int64_t rb, re;
	int score, qb, qe, seedcov, sub, csub; // sub: suboptimal score; csub: suboptimal inside the chain
	int sub_n, is_primary;
} mem_alnreg_t;

typedef struct {
	int low, high, failed;
	double avg, std;
} mem_pestat_t;

typedef struct {
	int64_t rb, re;
	int qb, qe, flag, qual;
	// optional info
	int score, sub;
	int64_t mb, me; // mb: mate start; -1 if single-end; -2 if mate unmapped
} bwahit_t;

typedef struct { size_t n, m; mem_chain_t *a;  } mem_chain_v;
typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;

extern int mem_verbose;

#ifdef __cplusplus
extern "C" {
#endif

smem_i *smem_itr_init(const bwt_t *bwt);
void smem_itr_destroy(smem_i *itr);
void smem_set_query(smem_i *itr, int len, const uint8_t *query);
const bwtintv_v *smem_next(smem_i *itr, int split_len, int split_width);

mem_opt_t *mem_opt_init(void);
void mem_fill_scmat(int a, int b, int8_t mat[25]);

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq);
int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *chains);
void mem_chain2aln(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_t *a);
uint32_t *mem_gen_cigar(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar);

int mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int n, bseq1_t *seqs);

void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4]);

#ifdef __cplusplus
}
#endif

#endif
