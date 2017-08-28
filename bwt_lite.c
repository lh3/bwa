#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bwt_lite.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int is_sa(const uint8_t *T, int *SA, int n);
int is_bwt(uint8_t *T, int n);

bwtl_t *bwtl_seq2bwtl(int len, const uint8_t *seq)
{
	bwtl_t *b;
	int i;
	b = (bwtl_t*)calloc(1, sizeof(bwtl_t));
	b->seq_len = len;

	{ // calculate b->bwt
		uint8_t *s;
		b->sa = (uint32_t*)calloc(len + 1, 4);
		is_sa(seq, (int*)b->sa, len);
		s = (uint8_t*)calloc(len + 1, 1);
		for (i = 0; i <= len; ++i) {
			if (b->sa[i] == 0) b->primary = i;
			else s[i] = seq[b->sa[i] - 1];
		}
		for (i = b->primary; i < len; ++i) s[i] = s[i + 1];
		b->bwt_size = (len + 15) / 16;
		b->bwt = (uint32_t*)calloc(b->bwt_size, 4);
		for (i = 0; i < len; ++i)
			b->bwt[i>>4] |= s[i] << ((15 - (i&15)) << 1);
		free(s);
	}
	{ // calculate b->occ
		uint32_t c[4];
		b->n_occ = (len + 15) / 16 * 4;
		b->occ = (uint32_t*)calloc(b->n_occ, 4);
		memset(c, 0, 16);
		for (i = 0; i < len; ++i) {
			if (i % 16 == 0)
				memcpy(b->occ + (i/16) * 4, c, 16);
			++c[bwtl_B0(b, i)];
		}
		memcpy(b->L2+1, c, 16);
		for (i = 2; i < 5; ++i) b->L2[i] += b->L2[i-1];
	}
	{ // generate cnt_table
		for (i = 0; i != 256; ++i) {
			uint32_t j, x = 0;
			for (j = 0; j != 4; ++j)
				x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
			b->cnt_table[i] = x;
		}
	}
	return b;
}
uint32_t bwtl_occ(const bwtl_t *bwt, uint32_t k, uint8_t c)
{
	uint32_t n, b;
	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (uint32_t)(-1)) return 0;
	if (k >= bwt->primary) --k; // because $ is not in bwt
	n = bwt->occ[k/16<<2|c];
	b = bwt->bwt[k/16] & ~((1U<<((15-(k&15))<<1)) - 1);
	n += (bwt->cnt_table[b&0xff] + bwt->cnt_table[b>>8&0xff]
		  + bwt->cnt_table[b>>16&0xff] + bwt->cnt_table[b>>24]) >> (c<<3) & 0xff;
	if (c == 0) n -= 15 - (k&15); // corrected for the masked bits
	return n;
}
void bwtl_occ4(const bwtl_t *bwt, uint32_t k, uint32_t cnt[4])
{
	uint32_t x, b;
	if (k == (uint32_t)(-1)) {
		memset(cnt, 0, 16);
		return;
	}
	if (k >= bwt->primary) --k; // because $ is not in bwt
	memcpy(cnt, bwt->occ + (k>>4<<2), 16);
	b = bwt->bwt[k>>4] & ~((1U<<((~k&15)<<1)) - 1);
	x = bwt->cnt_table[b&0xff] + bwt->cnt_table[b>>8&0xff]
		+ bwt->cnt_table[b>>16&0xff] + bwt->cnt_table[b>>24];
	x -= 15 - (k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}
void bwtl_2occ4(const bwtl_t *bwt, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4])
{
	bwtl_occ4(bwt, k, cntk);
	bwtl_occ4(bwt, l, cntl);
}
void bwtl_destroy(bwtl_t *bwt)
{
	if (bwt) {
		free(bwt->occ); free(bwt->bwt); free(bwt->sa);
		free(bwt);
	}
}
