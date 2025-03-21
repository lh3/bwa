#ifndef NEON_SSE_H
#define NEON_SSE_H

#include <arm_neon.h>

typedef uint8x16_t __m128i;

static inline __m128i _mm_load_si128(const __m128i *ptr) { return vld1q_u8((const uint8_t *) ptr); }
static inline __m128i _mm_set1_epi32(int n) { return vreinterpretq_u8_s32(vdupq_n_s32(n)); }
static inline void _mm_store_si128(__m128i *ptr, __m128i a) { vst1q_u8((uint8_t *) ptr, a); }

static inline __m128i _mm_adds_epu8(__m128i a, __m128i b) { return vqaddq_u8(a, b); }
static inline __m128i _mm_max_epu8(__m128i a, __m128i b) { return vmaxq_u8(a, b); }
static inline __m128i _mm_set1_epi8(int8_t n) { return vreinterpretq_u8_s8(vdupq_n_s8(n)); }
static inline __m128i _mm_subs_epu8(__m128i a, __m128i b) { return vqsubq_u8(a, b); }

#define M128I(a)  vreinterpretq_u8_s16((a))
#define UM128I(a) vreinterpretq_u8_u16((a))
#define S16(a)    vreinterpretq_s16_u8((a))
#define U16(a)    vreinterpretq_u16_u8((a))

static inline __m128i _mm_adds_epi16(__m128i a, __m128i b) { return M128I(vqaddq_s16(S16(a), S16(b))); }
static inline __m128i _mm_cmpgt_epi16(__m128i a, __m128i b) { return UM128I(vcgtq_s16(S16(a), S16(b))); }
static inline __m128i _mm_max_epi16(__m128i a, __m128i b) { return M128I(vmaxq_s16(S16(a), S16(b))); }
static inline __m128i _mm_set1_epi16(int16_t n) { return vreinterpretq_u8_s16(vdupq_n_s16(n)); }
static inline __m128i _mm_subs_epu16(__m128i a, __m128i b) { return UM128I(vqsubq_u16(U16(a), U16(b))); }

#undef M128I
#undef UM128I
#undef S16
#undef U16

#endif
