#ifndef BWAMEM_H_
#define BWAMEM_H_

#include "bwt.h"
#include "bntseq.h"
#include "bwa.h"

#define MEM_MAPQ_COEF 30.0
#define MEM_MAPQ_MAX  60

struct __smem_i;
typedef struct __smem_i smem_i;

#define MEM_F_HARDCLIP  0x1
#define MEM_F_PE        0x2
#define MEM_F_NOPAIRING 0x4
#define MEM_F_ALL       0x8
#define MEM_F_NO_MULTI  0x10
#define MEM_F_NO_RESCUE 0x20

typedef struct {
	int a, b, q, r;         // match score, mismatch penalty and gap open/extension penalty. A gap of size k costs q+k*r
	int pen_unpaired;       // phred-scaled penalty for unpaired reads
	int pen_clip;           // clipping penalty. This score is not deducted from the DP score.
	int w;                  // band width
	int zdrop;              // Z-dropoff

	int T;                  // output score threshold; only affecting output
	int flag;               // see MEM_F_* macros
	int min_seed_len;       // minimum seed length
	float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
	int split_width;        // split into a seed if its occurence is smaller than this value
	int max_occ;            // skip a seed if its occurence is larger than this value
	int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
	int n_threads;          // number of threads
	int chunk_size;         // process chunk_size-bp sequences in a batch
	float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
	float chain_drop_ratio; // drop a chain if its seed coverage is below chain_drop_ratio times the seed coverage of a better chain overlapping with the small chain
	int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
	int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
	int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
} mem_opt_t;

typedef struct {
	int64_t rb, re; // [rb,re): reference sequence in the alignment
	int qb, qe;     // [qb,qe): query sequence in the alignment
	int score;      // best local SW score
	int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
	int sub;        // 2nd best SW score
	int csub;       // SW score of a tandem hit
	int sub_n;      // approximate number of suboptimal hits
	int w;          // actual band width used in extension
	int seedcov;    // length of regions coverged by seeds
	int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
} mem_alnreg_t;

typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;

typedef struct {
	int low, high, failed;
	double avg, std;
} mem_pestat_t;

typedef struct { // This struct is only used for the convenience of API.
	int64_t pos;     // forward strand 5'-end mapping position
	int rid;         // reference sequence index in bntseq_t; <0 for unmapped
	int flag;        // extra flag
	uint32_t is_rev:1, mapq:8, NM:23; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
	int n_cigar;     // number of CIGAR operations
	uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234

	int score, sub;
} mem_aln_t;

#ifdef __cplusplus
extern "C" {
#endif

	smem_i *smem_itr_init(const bwt_t *bwt);
	void smem_itr_destroy(smem_i *itr);
	void smem_set_query(smem_i *itr, int len, const uint8_t *query);
	const bwtintv_v *smem_next(smem_i *itr, int split_len, int split_width);

	mem_opt_t *mem_opt_init(void);
	void mem_fill_scmat(int a, int b, int8_t mat[25]);

	/**
	 * Align a batch of sequences and generate the alignments in the SAM format
	 *
	 * This routine requires $seqs[i].{l_seq,seq,name} and write $seqs[i].sam.
	 * Note that $seqs[i].sam may consist of several SAM lines if the
	 * corresponding sequence has multiple primary hits.
	 *
	 * In the paired-end mode (i.e. MEM_F_PE is set in $opt->flag), query
	 * sequences must be interleaved: $n must be an even number and the 2i-th
	 * sequence and the (2i+1)-th sequence constitute a read pair. In this
	 * mode, there should be enough (typically >50) unique pairs for the
	 * routine to infer the orientation and insert size.
	 *
	 * @param opt    alignment parameters
	 * @param bwt    FM-index of the reference sequence
	 * @param bns    Information of the reference
	 * @param pac    2-bit encoded reference
	 * @param n      number of query sequences
	 * @param seqs   query sequences; $seqs[i].seq/sam to be modified after the call
	 */
	void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int n, bseq1_t *seqs);

	/**
	 * Find the aligned regions for one query sequence
	 *
	 * Note that this routine does not generate CIGAR. CIGAR should be
	 * generated later by mem_reg2aln() below.
	 *
	 * @param opt    alignment parameters
	 * @param bwt    FM-index of the reference sequence
	 * @param bns    Information of the reference
	 * @param pac    2-bit encoded reference
	 * @param l_seq  length of query sequence
	 * @param seq    query sequence
	 *
	 * @return       list of aligned regions.
	 */
	mem_alnreg_v mem_align1(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq);

	/**
	 * Generate CIGAR and forward-strand position from alignment region
	 *
	 * @param opt    alignment parameters
	 * @param bns    Information of the reference
	 * @param pac    2-bit encoded reference
	 * @param l_seq  length of query sequence
	 * @param seq    query sequence
	 * @param ar     one alignment region
	 *
	 * @return       CIGAR, strand, mapping quality and forward-strand position
	 */
	mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, const mem_alnreg_t *ar);

	/**
	 * Infer the insert size distribution from interleaved alignment regions
	 *
	 * This function can be called after mem_align1(), as long as paired-end
	 * reads are properly interleaved.
	 *
	 * @param opt    alignment parameters
	 * @param l_pac  length of concatenated reference sequence
	 * @param n      number of query sequences; must be an even number
	 * @param regs   region array of size $n; 2i-th and (2i+1)-th elements constitute a pair
	 * @param pes    inferred insert size distribution (output)
	 */
	void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4]);

#ifdef __cplusplus
}
#endif

#endif
