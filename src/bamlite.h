#ifndef BAMLITE_H_
#define BAMLITE_H_

#include <stdint.h>
#include <zlib.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define USE_VERBOSE_ZLIB_WRAPPERS

typedef gzFile bamFile;
#ifdef USE_VERBOSE_ZLIB_WRAPPERS
/* These print error messages on failure */
#  define bam_open(fn, mode)      bamlite_gzopen(fn, mode)
#  define bam_dopen(fd, mode)     gzdopen(fd, mode)
#  define bam_close(fp)           bamlite_gzclose(fp)
#  define bam_read(fp, buf, size) bamlite_gzread(fp, buf, size)
#else
#  define bam_open(fn, mode)      gzopen(fn, mode)
#  define bam_dopen(fd, mode)     gzdopen(fd, mode)
#  define bam_close(fp)           gzclose(fp)
#  define bam_read(fp, buf, size) gzread(fp, buf, size)
#endif /* USE_VERBOSE_ZLIB_WRAPPERS */

typedef struct {
	int32_t n_targets;
	char **target_name;
	uint32_t *target_len;
	size_t l_text, n_text;
	char *text;
} bam_header_t;

#define BAM_FPAIRED        1
#define BAM_FPROPER_PAIR   2
#define BAM_FUNMAP         4
#define BAM_FMUNMAP        8
#define BAM_FREVERSE      16
#define BAM_FMREVERSE     32
#define BAM_FREAD1        64
#define BAM_FREAD2       128
#define BAM_FSECONDARY   256
#define BAM_FQCFAIL      512
#define BAM_FDUP        1024

#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6

typedef struct {
	int32_t tid;
	int32_t pos;
	uint32_t bin:16, qual:8, l_qname:8;
	uint32_t flag:16, n_cigar:16;
	int32_t l_qseq;
	int32_t mtid;
	int32_t mpos;
	int32_t isize;
} bam1_core_t;

typedef struct {
	bam1_core_t core;
	int l_aux, data_len, m_data;
	uint8_t *data;
} bam1_t;

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define bam1_strand(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define bam1_mstrand(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
#define bam1_qname(b) ((char*)((b)->data))
#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_qual(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define bam1_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define bam1_aux(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2)

#define bam_init1() ((bam1_t*)calloc(1, sizeof(bam1_t)))
#define bam_destroy1(b) do {					\
		if (b) { free((b)->data); free(b); }	\
	} while (0)

extern int bam_is_be;

#ifdef __cplusplus
extern "C" {
#endif

	bam_header_t *bam_header_init(void);
	void bam_header_destroy(bam_header_t *header);
	bam_header_t *bam_header_read(bamFile fp);
	int bam_read1(bamFile fp, bam1_t *b);

#ifdef USE_VERBOSE_ZLIB_WRAPPERS
	gzFile bamlite_gzopen(const char *fn, const char *mode);
	int bamlite_gzread(gzFile file, void *ptr, unsigned int len);
	int bamlite_gzclose(gzFile file);
#endif /* USE_VERBOSE_ZLIB_WRAPPERS */

#ifdef __cplusplus
}
#endif

#endif
