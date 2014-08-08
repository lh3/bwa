#ifndef BWA_H_
#define BWA_H_

#include <stdint.h>
#include "bntseq.h"
#include "bwt.h"
#ifdef USE_HTSLIB
#include <htslib/sam.h>
#endif

#define BWA_IDX_BWT 0x1
#define BWA_IDX_BNS 0x2
#define BWA_IDX_PAC 0x4
#define BWA_IDX_ALL 0x7

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

#ifdef USE_HTSLIB
typedef struct {
	int l, m;
	bam1_t **bams;
} bams_t;
#endif

typedef struct {
	int l_seq;
#ifdef USE_HTSLIB
	char *name, *comment, *seq, *qual;
	bams_t *bams;
#else
	char *name, *comment, *seq, *qual, *sam;
#endif
} bseq1_t;

// This is here to faciliate passing around HTSLIB's bam_hdr_t structure when we are not compiling with HTSLIB
#ifndef USE_HTSLIB
typedef struct {
	void *ptr;
} bam_hdr_t; // DO NOT USE
#endif


extern int bwa_verbose;
extern char bwa_rg_id[256];

#ifdef __cplusplus
extern "C" {
#endif

	bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);

	void bwa_fill_scmat(int a, int b, int8_t mat[25]);
	uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM);
	uint32_t *bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM);

	char *bwa_idx_infer_prefix(const char *hint);
	bwt_t *bwa_idx_load_bwt(const char *hint);

	bwaidx_t *bwa_idx_load(const char *hint, int which);
	void bwa_idx_destroy(bwaidx_t *idx);

	void bwa_print_sam_hdr(const bntseq_t *bns, const char *rg_line);
#ifdef USE_HTSLIB
	void bwa_format_sam_hdr(const bntseq_t *bns, const char *rg_line, kstring_t *str);
#endif
	char *bwa_set_rg(const char *s);

#ifdef USE_HTSLIB
	bams_t *bams_init();
	void bams_add(bams_t *bams, bam1_t *b);
	void bams_destroy(bams_t *bams);
#endif

#ifdef __cplusplus
}
#endif

#endif
