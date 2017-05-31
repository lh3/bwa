#ifndef BWT_LITE_H_
#define BWT_LITE_H_

#include <stdint.h>

typedef struct {
	uint32_t seq_len, bwt_size, n_occ;
	uint32_t primary;
	uint32_t *bwt, *occ, *sa, L2[5];
	uint32_t cnt_table[256];
} bwtl_t;

#define bwtl_B0(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

#ifdef __cplusplus
extern "C" {
#endif

	bwtl_t *bwtl_seq2bwtl(int len, const uint8_t *seq);
	uint32_t bwtl_occ(const bwtl_t *bwt, uint32_t k, uint8_t c);
	void bwtl_occ4(const bwtl_t *bwt, uint32_t k, uint32_t cnt[4]);
	void bwtl_2occ4(const bwtl_t *bwt, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4]);
	void bwtl_destroy(bwtl_t *bwt);

#ifdef __cplusplus
}
#endif

#endif
