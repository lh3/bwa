#ifndef SCALAR_SSE_H
#define SCALAR_SSE_H

#include <assert.h>
#include <stdint.h>
#include <string.h>

typedef union m128i {
	uint8_t u8[16];
	int16_t i16[8];
} __m128i;

static inline __m128i _mm_set1_epi32(int32_t n) {
	assert(n >= 0 && n <= 255);
	__m128i r; memset(&r, n, sizeof r); return r;
}

static inline __m128i _mm_load_si128(const __m128i *ptr) { __m128i r; memcpy(&r, ptr, sizeof r); return r; }
static inline void _mm_store_si128(__m128i *ptr, __m128i a) { memcpy(ptr, &a, sizeof a); }

static inline int m128i_allzero(__m128i a) {
	static const char zero[] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
	return memcmp(&a, zero, sizeof a) == 0;
}

static inline __m128i _mm_slli_si128(__m128i a, int n) {
	int i;
	memmove(&a.u8[n], &a.u8[0], 16 - n);
	for (i = 0; i < n; i++) a.u8[i] = 0;
	return a;
}

static inline __m128i _mm_adds_epu8(__m128i a, __m128i b) {
	int i;
	for (i = 0; i < 16; i++) {
		uint16_t aa = a.u8[i];
		aa += b.u8[i];
		a.u8[i] = (aa < 256)? aa : 255;
	}
	return a;
}

static inline __m128i _mm_max_epu8(__m128i a, __m128i b) {
	int i;
	for (i = 0; i < 16; i++)
		if (a.u8[i] < b.u8[i]) a.u8[i] = b.u8[i];
	return a;
}

static inline uint8_t m128i_max_u8(__m128i a) {
	uint8_t max = 0;
	int i;
	for (i = 0; i < 16; i++)
		if (max < a.u8[i]) max = a.u8[i];
	return max;
}

static inline __m128i _mm_set1_epi8(int8_t n) { __m128i r; memset(&r, n, sizeof r); return r; }

static inline __m128i _mm_subs_epu8(__m128i a, __m128i b) {
	int i;
	for (i = 0; i < 16; i++) {
		int16_t aa = a.u8[i];
		aa -= b.u8[i];
		a.u8[i] = (aa >= 0)? aa : 0;
	}
	return a;
}

static inline __m128i _mm_adds_epi16(__m128i a, __m128i b) {
	int i;
	for (i = 0; i < 8; i++) {
		int32_t aa = a.i16[i];
		aa += b.i16[i];
		a.i16[i] = (aa < 32768)? aa : 32767;
	}
	return a;
}

static inline __m128i _mm_cmpgt_epi16(__m128i a, __m128i b) {
	int i;
	for (i = 0; i < 8; i++)
		a.i16[i] = (a.i16[i] > b.i16[i])? 0xffff : 0x0000;
	return a;
}

static inline __m128i _mm_max_epi16(__m128i a, __m128i b) {
	int i;
	for (i = 0; i < 8; i++)
		if (a.i16[i] < b.i16[i]) a.i16[i] = b.i16[i];
	return a;
}

static inline __m128i _mm_set1_epi16(int16_t n) {
	__m128i r;
	r.i16[0] = r.i16[1] = r.i16[2] = r.i16[3] =
	r.i16[4] = r.i16[5] = r.i16[6] = r.i16[7] = n;
	return r;
}

static inline int16_t m128i_max_s16(__m128i a) {
	int16_t max = -32768;
	int i;
	for (i = 0; i < 8; i++)
		if (max < a.i16[i]) max = a.i16[i];
	return max;
}

static inline __m128i _mm_subs_epu16(__m128i a, __m128i b) {
	int i;
	for (i = 0; i < 8; i++) {
		int32_t aa = a.i16[i];
		aa -= b.i16[i];
		a.i16[i] = (aa >= 0)? aa : 0;
	}
	return a;
}

#endif
