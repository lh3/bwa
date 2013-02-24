#ifndef BWAMEM_H_
#define BWAMEM_H_

#include "bwt.h"
#include "bntseq.h"
#include "bwa.h"

#define MEM_MAPQ_COEF 40.0
#define MEM_MAPQ_MAX  60

struct __smem_i;
typedef struct __smem_i smem_i;

#define MEM_F_HARDCLIP  0x1
#define MEM_F_PE        0x2
#define MEM_F_NOPAIRING 0x4

typedef struct {
	int a, b, q, r, w;
	int flag;
	int split_width;
	int min_seed_len, max_occ, max_chain_gap;
	int n_threads, chunk_size;
	int pe_dir;
	float mask_level;
	float chain_drop_ratio;
	float split_factor; // split into a seed if MEM is longer than min_seed_len*split_factor
	int pen_unpaired; // phred-scaled penalty for unpaired reads
	int max_ins; // maximum insert size
	int8_t mat[25]; // scoring matrix; mat[0] == 0 if unset
} mem_opt_t;

typedef struct {
	int64_t rb, re;
	int score, qb, qe, seedcov, sub, csub; // sub: suboptimal score; csub: suboptimal inside the chain
	int sub_n; // approximate number of suboptimal hits
	int secondary; // non-negative if the hit is secondary
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
} bwahit_t;

typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;

#ifdef __cplusplus
extern "C" {
#endif

	smem_i *smem_itr_init(const bwt_t *bwt);
	void smem_itr_destroy(smem_i *itr);
	void smem_set_query(smem_i *itr, int len, const uint8_t *query);
	const bwtintv_v *smem_next(smem_i *itr, int split_len, int split_width);

	mem_opt_t *mem_opt_init(void);
	void mem_fill_scmat(int a, int b, int8_t mat[25]);

	int mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int n, bseq1_t *seqs);

	void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4]);

#ifdef __cplusplus
}
#endif

#endif
