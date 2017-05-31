#ifndef BWTALN_H
#define BWTALN_H

#include <stdint.h>
#include "bwt.h"

#define BWA_TYPE_NO_MATCH 0
#define BWA_TYPE_UNIQUE 1
#define BWA_TYPE_REPEAT 2
#define BWA_TYPE_MATESW 3

#define SAM_FPD   1 // paired
#define SAM_FPP   2 // properly paired
#define SAM_FSU   4 // self-unmapped
#define SAM_FMU   8 // mate-unmapped
#define SAM_FSR  16 // self on the reverse strand
#define SAM_FMR  32 // mate on the reverse strand
#define SAM_FR1  64 // this is read one
#define SAM_FR2 128 // this is read two
#define SAM_FSC 256 // secondary alignment

#define BWA_AVG_ERR 0.02
#define BWA_MIN_RDLEN 35 // for read trimming

#define BWA_MAX_BCLEN 63 // maximum barcode length; 127 is the maximum

#ifndef bns_pac
#define bns_pac(pac, k) ((pac)[(k)>>2] >> ((~(k)&3)<<1) & 3)
#endif

#define FROM_M 0
#define FROM_I 1
#define FROM_D 2
#define FROM_S 3

#define SAI_MAGIC "SAI\1"

typedef struct {
	bwtint_t w;
	int bid;
} bwt_width_t;

typedef struct {
	uint64_t n_mm:8, n_gapo:8, n_gape:8, score:20, n_ins:10, n_del:10;
	bwtint_t k, l;
} bwt_aln1_t;

typedef uint16_t bwa_cigar_t;
/* rgoya: If changing order of bytes, beware of operations like:
 *     s->cigar[0] += s->full_len - s->len;
 */
#define CIGAR_OP_SHIFT 14
#define CIGAR_LN_MASK 0x3fff

#define __cigar_op(__cigar) ((__cigar)>>CIGAR_OP_SHIFT)
#define __cigar_len(__cigar) ((__cigar)&CIGAR_LN_MASK)
#define __cigar_create(__op, __len) ((__op)<<CIGAR_OP_SHIFT | (__len))

typedef struct {
	uint32_t n_cigar:15, gap:8, mm:8, strand:1;
	int ref_shift;
	bwtint_t pos;
	bwa_cigar_t *cigar;
} bwt_multi1_t;

typedef struct {
	char *name;
	ubyte_t *seq, *rseq, *qual;
	uint32_t len:20, strand:1, type:2, dummy:1, extra_flag:8;
	uint32_t n_mm:8, n_gapo:8, n_gape:8, mapQ:8;
	int score;
	int clip_len;
	// alignments in SA coordinates
	int n_aln;
	bwt_aln1_t *aln;
	// multiple hits
	int n_multi;
	bwt_multi1_t *multi;
	// alignment information
	bwtint_t sa, pos;
	uint64_t c1:28, c2:28, seQ:8; // number of top1 and top2 hits; single-end mapQ
	int ref_shift;
	int n_cigar;
	bwa_cigar_t *cigar;
	// for multi-threading only
	int tid;
	// barcode
	char bc[BWA_MAX_BCLEN+1]; // null terminated; up to BWA_MAX_BCLEN bases
	// NM and MD tags
	uint32_t full_len:20, nm:12;
	char *md;
} bwa_seq_t;

#define BWA_MODE_GAPE       0x01
#define BWA_MODE_COMPREAD   0x02
#define BWA_MODE_LOGGAP     0x04
#define BWA_MODE_CFY        0x08
#define BWA_MODE_NONSTOP    0x10
#define BWA_MODE_BAM        0x20
#define BWA_MODE_BAM_SE     0x40
#define BWA_MODE_BAM_READ1  0x80
#define BWA_MODE_BAM_READ2  0x100
#define BWA_MODE_IL13       0x200

typedef struct {
	int s_mm, s_gapo, s_gape;
	int mode; // bit 24-31 are the barcode length
	int indel_end_skip, max_del_occ, max_entries;
	float fnr;
	int max_diff, max_gapo, max_gape;
	int max_seed_diff, seed_len;
	int n_threads;
	int max_top2;
	int trim_qual;
} gap_opt_t;

#define BWA_PET_STD   1

typedef struct {
	int max_isize, force_isize;
	int max_occ;
	int n_multi, N_multi;
	int type, is_sw, is_preload;
	double ap_prior;
} pe_opt_t;

struct __bwa_seqio_t;
typedef struct __bwa_seqio_t bwa_seqio_t;

#ifdef __cplusplus
extern "C" {
#endif

	gap_opt_t *gap_init_opt();
	void bwa_aln_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt);

	bwa_seqio_t *bwa_seq_open(const char *fn);
	bwa_seqio_t *bwa_bam_open(const char *fn, int which);
	void bwa_seq_close(bwa_seqio_t *bs);
	void seq_reverse(int len, ubyte_t *seq, int is_comp);
	bwa_seq_t *bwa_read_seq(bwa_seqio_t *seq, int n_needed, int *n, int mode, int trim_qual);
	void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs);

	int bwa_cal_maxdiff(int l, double err, double thres);
	void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);

	void bwa_cs2nt_core(bwa_seq_t *p, bwtint_t l_pac, ubyte_t *pac);

#ifdef __cplusplus
}
#endif

#endif
