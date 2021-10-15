#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "rle.h"

const uint8_t rle_auxtab[8] = { 0x01, 0x11, 0x21, 0x31, 0x03, 0x13, 0x07, 0x17 };

// insert symbol $a after $x symbols in $str; marginal counts added to $cnt; returns the size increase
int rle_insert_cached(uint8_t *block, int64_t x, int a, int64_t rl, int64_t cnt[6], const int64_t ec[6], int *beg, int64_t bc[6])
{
	uint16_t *nptr = (uint16_t*)block;
	int diff;

	block += 2; // skip the first 2 counting bytes
	if (*nptr == 0) {
		memset(cnt, 0, 48);
		diff = rle_enc1(block, a, rl);
	} else {
		uint8_t *p, *end = block + *nptr, *q;
		int64_t pre, z, l = 0, tot, beg_l;
		int c = -1, n_bytes = 0, n_bytes2, t = 0;
		uint8_t tmp[24];
		beg_l = bc[0] + bc[1] + bc[2] + bc[3] + bc[4] + bc[5];
		tot   = ec[0] + ec[1] + ec[2] + ec[3] + ec[4] + ec[5];
		if (x < beg_l) {
			beg_l = 0, *beg = 0;
			memset(bc, 0, 48);
		}
		if (x == beg_l) {
			p = q = block + (*beg); z = beg_l;
			memcpy(cnt, bc, 48);
		} else if (x - beg_l <= ((tot-beg_l)>>1) + ((tot-beg_l)>>3)) { // forward
			z = beg_l; p = block + (*beg);
			memcpy(cnt, bc, 48);
			while (z < x) {
				rle_dec1(p, c, l);
				z += l; cnt[c] += l;
			}
			for (q = p - 1; *q>>6 == 2; --q);
		} else { // backward
			memcpy(cnt, ec, 48);
			z = tot; p = end;
			while (z >= x) {
				--p;
				if (*p>>6 != 2) {
					l |= *p>>7? (int64_t)rle_auxtab[*p>>3&7]>>4 << t : *p>>3;
					z -= l; cnt[*p&7] -= l;
					l = 0; t = 0;
				} else {
					l |= (*p&0x3fL) << t;
					t += 6;
				}
			}
			q = p;
			rle_dec1(p, c, l);
			z += l; cnt[c] += l;
		}
		*beg = q - block;
		memcpy(bc, cnt, 48);
		bc[c] -= l;
		n_bytes = p - q;
		if (x == z && a != c && p < end) { // then try the next run
			int tc;
			int64_t tl;
			q = p;
			rle_dec1(q, tc, tl);
			if (a == tc)
				c = tc, n_bytes = q - p, l = tl, z += l, p = q, cnt[tc] += tl;
		}
		if (z != x) cnt[c] -= z - x;
		pre = x - (z - l); p -= n_bytes;
		if (a == c) { // insert to the same run
			n_bytes2 = rle_enc1(tmp, c, l + rl);
		} else if (x == z) { // at the end; append to the existing run
			p += n_bytes; n_bytes = 0;
			n_bytes2 = rle_enc1(tmp, a, rl);
		} else { // break the current run
			n_bytes2 = rle_enc1(tmp, c, pre);
			n_bytes2 += rle_enc1(tmp + n_bytes2, a, rl);
			n_bytes2 += rle_enc1(tmp + n_bytes2, c, l - pre);
		}
		if (n_bytes != n_bytes2 && end != p + n_bytes) // size changed
			memmove(p + n_bytes2, p + n_bytes, end - p - n_bytes);
		memcpy(p, tmp, n_bytes2);
		diff = n_bytes2 - n_bytes;
	}
	return (*nptr += diff);
}

int rle_insert(uint8_t *block, int64_t x, int a, int64_t rl, int64_t cnt[6], const int64_t ec[6])
{
	int beg = 0;
	int64_t bc[6];
	memset(bc, 0, 48);
	return rle_insert_cached(block, x, a, rl, cnt, ec, &beg, bc);
}

void rle_split(uint8_t *block, uint8_t *new_block)
{
	int n = *(uint16_t*)block;
	uint8_t *end = block + 2 + n, *q = block + 2 + (n>>1);
	while (*q>>6 == 2) --q;
	memcpy(new_block + 2, q, end - q);
	*(uint16_t*)new_block = end - q;
	*(uint16_t*)block = q - block - 2;
}

void rle_count(const uint8_t *block, int64_t cnt[6])
{
	const uint8_t *q = block + 2, *end = q + *(uint16_t*)block;
	while (q < end) {
		int c;
		int64_t l;
		rle_dec1(q, c, l);
		cnt[c] += l;
	}
}

void rle_print(const uint8_t *block, int expand)
{
	const uint16_t *p = (const uint16_t*)block;
	const uint8_t *q = block + 2, *end = block + 2 + *p;
	while (q < end) {
		int c;
		int64_t l, x;
		rle_dec1(q, c, l);
		if (expand) for (x = 0; x < l; ++x) putchar("$ACGTN"[c]);
		else printf("%c%ld", "$ACGTN"[c], (long)l);
	}
	putchar('\n');
}

void rle_rank2a(const uint8_t *block, int64_t x, int64_t y, int64_t *cx, int64_t *cy, const int64_t ec[6])
{
	int a;
	int64_t tot, cnt[6];
	const uint8_t *p;

	y = y >= x? y : x;
	tot = ec[0] + ec[1] + ec[2] + ec[3] + ec[4] + ec[5];
	if (tot == 0) return;
	if (x <= (tot - y) + (tot>>3)) {
		int c = 0;
		int64_t l, z = 0;
		memset(cnt, 0, 48);
		p = block + 2;
		while (z < x) {
			rle_dec1(p, c, l);
			z += l; cnt[c] += l;
		}
		for (a = 0; a != 6; ++a) cx[a] += cnt[a];
		cx[c] -= z - x;
		if (cy) {
			while (z < y) {
				rle_dec1(p, c, l);
				z += l; cnt[c] += l;
			}
			for (a = 0; a != 6; ++a) cy[a] += cnt[a];
			cy[c] -= z - y;
		}
	} else {
#define move_backward(_x) \
		while (z >= (_x)) { \
			--p; \
			if (*p>>6 != 2) { \
				l |= *p>>7? (int64_t)rle_auxtab[*p>>3&7]>>4 << t : *p>>3; \
				z -= l; cnt[*p&7] -= l; \
				l = 0; t = 0; \
			} else { \
				l |= (*p&0x3fL) << t; \
				t += 6; \
			} \
		} \

		int t = 0;
		int64_t l = 0, z = tot;
		memcpy(cnt, ec, 48);
		p = block + 2 + *(const uint16_t*)block;
		if (cy) {
			move_backward(y)
			for (a = 0; a != 6; ++a) cy[a] += cnt[a];
			cy[*p&7] += y - z;
		}
		move_backward(x)
		for (a = 0; a != 6; ++a) cx[a] += cnt[a];
		cx[*p&7] += x - z;

#undef move_backward
	}
}
