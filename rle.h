#ifndef RLE6_H_
#define RLE6_H_

#include <stdint.h>

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#else
#define LIKELY(x) (x)
#endif
#ifdef __cplusplus

extern "C" {
#endif

	int rle_insert_cached(uint8_t *block, int64_t x, int a, int64_t rl, int64_t cnt[6], const int64_t ec[6], int *beg, int64_t bc[6]);
	int rle_insert(uint8_t *block, int64_t x, int a, int64_t rl, int64_t cnt[6], const int64_t end_cnt[6]);
	void rle_split(uint8_t *block, uint8_t *new_block);
	void rle_count(const uint8_t *block, int64_t cnt[6]);
	void rle_rank2a(const uint8_t *block, int64_t x, int64_t y, int64_t *cx, int64_t *cy, const int64_t ec[6]);
	#define rle_rank1a(block, x, cx, ec) rle_rank2a(block, x, -1, cx, 0, ec)

	void rle_print(const uint8_t *block, int expand);

#ifdef __cplusplus
}
#endif

/******************
 *** 43+3 codec ***
 ******************/

const uint8_t rle_auxtab[8];

#define RLE_MIN_SPACE 18
#define rle_nptr(block) ((uint16_t*)(block))

// decode one run (c,l) and move the pointer p
#define rle_dec1(p, c, l) do { \
		(c) = *(p) & 7; \
		if (LIKELY((*(p)&0x80) == 0)) { \
			(l) = *(p)++ >> 3; \
		} else if (LIKELY(*(p)>>5 == 6)) { \
			(l) = (*(p)&0x18L)<<3L | ((p)[1]&0x3fL); \
			(p) += 2; \
		} else { \
			int n = ((*(p)&0x10) >> 2) + 4; \
			(l) = *(p)++ >> 3 & 1; \
			while (--n) (l) = ((l)<<6) | (*(p)++&0x3fL); \
		} \
	} while (0)

static inline int rle_enc1(uint8_t *p, int c, int64_t l)
{
	if (l < 1LL<<4) {
		*p = l << 3 | c;
		return 1;
	} else if (l < 1LL<<8) {
		*p = 0xC0 | l >> 6 << 3 | c;
		p[1] = 0x80 | (l & 0x3f);
		return 2;
	} else if (l < 1LL<<19) {
		*p = 0xE0 | l >> 18 << 3 | c;
		p[1] = 0x80 | (l >> 12 & 0x3f);
		p[2] = 0x80 | (l >> 6 & 0x3f);
		p[3] = 0x80 | (l & 0x3f);
		return 4;
	} else {
		int i, shift = 36;
		*p = 0xF0 | l >> 42 << 3 | c;
		for (i = 1; i < 8; ++i, shift -= 6)
			p[i] = 0x80 | (l>>shift & 0x3f);
		return 8;
	}
}

#endif
