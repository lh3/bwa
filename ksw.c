/* The MIT License

   Copyright (c) 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdlib.h>
#include <stdint.h>
#include "ksw.h"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef _NO_SSE2
#include <emmintrin.h>

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/***************
 *** SSE2 SW ***
 ***************/

struct _ksw_query_t {
	int qlen, slen;
	uint8_t shift, mdiff, max, size;
	__m128i *qp, *H0, *H1, *E, *Hmax;
};

ksw_query_t *ksw_qinit(int size, int qlen, const uint8_t *query, int m, const int8_t *mat)
{
	ksw_query_t *q;
	int slen, a, tmp, p;

	size = size > 1? 2 : 1;
	p = 8 * (3 - size); // # values per __m128i
	slen = (qlen + p - 1) / p; // segmented length
	q = malloc(sizeof(ksw_query_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
	q->qp = (__m128i*)(((size_t)q + sizeof(ksw_query_t) + 15) >> 4 << 4); // align memory
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
	// compute shift
	tmp = m * m;
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
	q->max = q->mdiff;
	q->shift = 256 - q->shift; // NB: q->shift is uint8_t
	q->mdiff += q->shift; // this is the difference between the min and max scores
	// An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
		int8_t *t = (int8_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
		}
	} else {
		int16_t *t = (int16_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]);
		}
	}
	return q;
}

int ksw_sse2_16(ksw_query_t *q, int tlen, const uint8_t *target, ksw_aux_t *a) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0;
	uint64_t *b;
	__m128i zero, gapoe, gape, shift, *H0, *H1, *E, *Hmax;

#define __max_16(ret, xx) do { \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
    	(ret) = _mm_extract_epi16((xx), 0) & 0x00ff; \
	} while (0)

	// initialization
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	gapoe = _mm_set1_epi8(a->gapo + a->gape);
	gape = _mm_set1_epi8(a->gape);
	shift = _mm_set1_epi8(q->shift);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, cmp, imax;
		__m128i e, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
		for (j = 0; LIKELY(j < slen); ++j) {
			/* SW cells are computed in the following order:
			 *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */
			// compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));
			h = _mm_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
			e = _mm_load_si128(E + j); // e=E'(i,j)
			h = _mm_max_epu8(h, e);
			h = _mm_max_epu8(h, f); // h=H'(i,j)
			max = _mm_max_epu8(max, h); // set max
			_mm_store_si128(H1 + j, h); // save to H'(i,j)
			// now compute E'(i+1,j)
			h = _mm_subs_epu8(h, gapoe); // h=H'(i,j)-gapo
			e = _mm_subs_epu8(e, gape); // e=E'(i,j)-gape
			e = _mm_max_epu8(e, h); // e=E'(i+1,j)
			_mm_store_si128(E + j, e); // save to E'(i+1,j)
			// now compute F'(i,j+1)
			f = _mm_subs_epu8(f, gape);
			f = _mm_max_epu8(f, h);
			// get H'(i-1,j) and prepare for the next j
			h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
		}
		// NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
		for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epu8(h, f); // h=H'(i,j)
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, gapoe);
				f = _mm_subs_epu8(f, gape);
				cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h), zero));
				if (UNLIKELY(cmp == 0xffff)) goto end_loop16;
			}
		}
end_loop16:
		//int k;for (k=0;k<16;++k)printf("%d ", ((uint8_t*)&max)[k]);printf("\n");
		__max_16(imax, max); // imax is the maximum number in max
		if (imax >= a->T) { // write the b array; this condition adds branching unfornately
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i; // te is the end position on the target
			for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax + q->shift >= 255) break;
		}
		S = H1; H1 = H0; H0 = S; // swap H0 and H1
	}
	a->score = gmax; a->te = te;
	{ // get a->qe, the end of query match; find the 2nd best score
		int max = -1, low, high, qlen = slen * 16;
		uint8_t *t = (uint8_t*)Hmax;
		for (i = 0, a->qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, a->qe = i / 16 + i % 16 * slen;
		//printf("%d,%d\n", max, gmax);
		i = (a->score + q->max - 1) / q->max;
		low = te - i; high = te + i;
		for (i = 0, a->score2 = 0; i < n_b; ++i) {
			int e = (int32_t)b[i];
			if ((e < low || e > high) && b[i]>>32 > (uint32_t)a->score2)
				a->score2 = b[i]>>32, a->te2 = e;
		}
	}
	free(b);
	return a->score + q->shift >= 255? 255 : a->score;
}

int ksw_sse2_8(ksw_query_t *q, int tlen, const uint8_t *target, ksw_aux_t *a) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0;
	uint64_t *b;
	__m128i zero, gapoe, gape, *H0, *H1, *E, *Hmax;

#define __max_8(ret, xx) do { \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    	(ret) = _mm_extract_epi16((xx), 0); \
	} while (0)

	// initialization
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	gapoe = _mm_set1_epi16(a->gapo + a->gape);
	gape = _mm_set1_epi16(a->gape);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 2);
		for (j = 0; LIKELY(j < slen); ++j) {
			h = _mm_adds_epi16(h, *S++);
			e = _mm_load_si128(E + j);
			h = _mm_max_epi16(h, e);
			h = _mm_max_epi16(h, f);
			max = _mm_max_epi16(max, h);
			_mm_store_si128(H1 + j, h);
			h = _mm_subs_epu16(h, gapoe);
			e = _mm_subs_epu16(e, gape);
			e = _mm_max_epi16(e, h);
			_mm_store_si128(E + j, e);
			f = _mm_subs_epu16(f, gape);
			f = _mm_max_epi16(f, h);
			h = _mm_load_si128(H0 + j);
		}
		for (k = 0; LIKELY(k < 16); ++k) {
			f = _mm_slli_si128(f, 2);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epi16(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu16(h, gapoe);
				f = _mm_subs_epu16(f, gape);
				if(UNLIKELY(!_mm_movemask_epi8(_mm_cmpgt_epi16(f, h)))) goto end_loop8;
			}
		}
end_loop8:
		__max_8(imax, max);
		if (imax >= a->T) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
		}
		S = H1; H1 = H0; H0 = S;
	}
	a->score = gmax; a->te = te;
	{
		int max = -1, low, high, qlen = slen * 8;
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, a->qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, a->qe = i / 8 + i % 8 * slen;
		i = (a->score + q->max - 1) / q->max;
		low = te - i; high = te + i;
		for (i = 0, a->score2 = 0; i < n_b; ++i) {
			int e = (int32_t)b[i];
			if ((e < low || e > high) && b[i]>>32 > (uint32_t)a->score2)
				a->score2 = b[i]>>32, a->te2 = e;
		}
	}
	free(b);
	return a->score;
}

int ksw_sse2(ksw_query_t *q, int tlen, const uint8_t *target, ksw_aux_t *a)
{
	if (q->size == 1) return ksw_sse2_16(q, tlen, target, a);
	else return ksw_sse2_8(q, tlen, target, a);
}

#endif // _NO_SSE2

/********************
 *** SW extension ***
 ********************/

typedef struct {
	int32_t h, e;
} eh_t;

int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int h0, int *_qle, int *_tle)
{
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, gapoe = gapo + gape, beg, end, max, max_i, max_j, max_gap;
	if (h0 < 0) h0 = 0;
	// allocate memory
	qp = malloc(qlen * m);
	eh = calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = h0; eh[1].h = h0 > gapoe? h0 - gapoe : 0;
	for (j = 2; j <= qlen && eh[j-1].h > gape; ++j)
		eh[j].h = eh[j-1].h - gape;
	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_gap = (int)((double)(qlen * max - gapo) / gape + 1.);
	max_gap = max_gap > 1? max_gap : 1;
	w = w < max_gap? w : max_gap;
	// DP loop
	max = h0, max_i = max_j = -1;
	beg = 0, end = qlen;
	for (i = 0; LIKELY(i < tlen); ++i) {
		int f = 0, h1, m = 0, mj = -1;
		int8_t *q = &qp[target[i] * qlen];
		// compute the first column
		h1 = h0 - (gapo + gape * (i + 1));
		if (h1 < 0) h1 = 0;
		// apply the band and the constraint (if provided)
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen;
		for (j = beg; LIKELY(j < end); ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int h = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			h += q[j];
			h = h > e? h : e;
			h = h > f? h : f;
			h1 = h;             // save H(i,j) to h1 for the next column
			mj = m > h? mj : j;
			m = m > h? m : h;   // m is stored at eh[mj+1]
			h -= gapoe;
			h = h > 0? h : 0;
			e -= gape;
			e = e > h? e : h;   // computed E(i+1,j)
			p->e = e;           // save E(i+1,j) for the next row
			f -= gape;
			f = f > h? f : h;   // computed F(i,j+1)
		}
		eh[end].h = h1; eh[end].e = 0;
		if (m == 0) break;
		if (m > max) max = m, max_i = i, max_j = mj;
		// update beg and end for the next round
		for (j = mj; j >= beg && eh[j].h; --j);
		beg = j + 1;
		for (j = mj + 2; j <= end && eh[j].h; ++j);
		end = j;
		//beg = 0; end = qlen; // uncomment this line for debugging
	}
	free(eh); free(qp);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	return max;
}

/********************
 * Global alignment *
 ********************/

#define MINUS_INF -0x40000000

static inline uint32_t *push_cigar(int *n_cigar, int *m_cigar, uint32_t *cigar, int op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = realloc(cigar, (*m_cigar) << 4);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

int ksw_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *n_cigar_, uint32_t **cigar_)
{
	eh_t *eh;
	int8_t *qp; // query profile
	int i, j, k, gapoe = gapo + gape, score, n_col;
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	if (n_cigar_) *n_cigar_ = 0;
	// allocate memory
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	z = malloc(n_col * tlen);
	qp = malloc(qlen * m);
	eh = calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = 0; eh[0].e = MINUS_INF;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(gapo + gape * j), eh[j].e = MINUS_INF;
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	// DP loop
	for (i = 0; LIKELY(i < tlen); ++i) { // target sequence is in the outer loop
		int32_t f = MINUS_INF, h1, beg, end;
		int8_t *q = &qp[target[i] * qlen];
		uint8_t *zi = &z[i * n_col];
		beg = i > w? i - w : 0;
		end = i + w + 1 < qlen? i + w + 1 : qlen; // only loop through [beg,end) of the query sequence
		h1 = beg == 0? -(gapo + gape * (i + 1)) : MINUS_INF;
		for (j = beg; LIKELY(j < end); ++j) {
			// This loop is organized in a similar way to ksw_extend() and ksw_sse2(), except:
			// 1) not checking h>0; 2) recording direction for backtracking
			eh_t *p = &eh[j];
			int32_t h = p->h, e = p->e;
			uint8_t d; // direction
			p->h = h1;
			h += q[j];
			d = h > e? 0 : 1;
			h = h > e? h : e;
			d = h > f? d : 2;
			h = h > f? h : f;
			h1 = h;
			h -= gapoe;
			e -= gape;
			d |= e > h? 1<<2 : 0;
			e = e > h? e : h;
			p->e = e;
			f -= gape;
			d |= f > h? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
			f = f > h? f : h;
			zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
		}
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}
	score = eh[qlen].h;
	if (n_cigar_ && cigar_) { // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1; k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1; // (i,k) points to the last cell
		while (i >= 0 && k >= 0) {
			which = z[i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			if (which == 0)      cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --k;
			else if (which == 1) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i;
			else                 cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --k;
		}
		if (i >= 0) push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		if (k >= 0) push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;
	}
	free(eh); free(qp); free(z);
	return score;
}

/*******************************************
 * Main function (not compiled by default) *
 *******************************************/

#if defined(_KSW_MAIN) && !defined(_NO_SSE2)

#include <unistd.h>
#include <stdio.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int main(int argc, char *argv[])
{
	int c, sa = 1, sb = 3, i, j, k, forward_only = 0, size = 2;
	int8_t mat[25];
	ksw_aux_t a;
	gzFile fpt, fpq;
	kseq_t *kst, *ksq;
	// parse command line
	a.gapo = 5; a.gape = 2; a.T = 10;
	while ((c = getopt(argc, argv, "a:b:q:r:ft:s:")) >= 0) {
		switch (c) {
			case 'a': sa = atoi(optarg); break;
			case 'b': sb = atoi(optarg); break;
			case 'q': a.gapo = atoi(optarg); break;
			case 'r': a.gape = atoi(optarg); break;
			case 't': a.T = atoi(optarg); break;
			case 'f': forward_only = 1; break;
			case 's': size = atoi(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: ksw [-s%d] [-a%d] [-b%d] [-q%d] [-r%d] <target.fa> <query.fa>\n", size, sa, sb, a.gapo, a.gape);
		return 1;
	}
	// initialize scoring matrix
	for (i = k = 0; i < 5; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? sa : -sb;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
	// open file
	fpt = gzopen(argv[optind],   "r"); kst = kseq_init(fpt);
	fpq = gzopen(argv[optind+1], "r"); ksq = kseq_init(fpq);
	// all-pair alignment
	while (kseq_read(ksq) > 0) {
		ksw_query_t *q[2];
		for (i = 0; i < ksq->seq.l; ++i) ksq->seq.s[i] = seq_nt4_table[(int)ksq->seq.s[i]];
		q[0] = ksw_qinit(size, ksq->seq.l, (uint8_t*)ksq->seq.s, 5, mat);
		if (!forward_only) { // reverse
			for (i = 0; i < ksq->seq.l/2; ++i) {
				int t = ksq->seq.s[i];
				ksq->seq.s[i] = ksq->seq.s[ksq->seq.l-1-i];
				ksq->seq.s[ksq->seq.l-1-i] = t;
			}
			for (i = 0; i < ksq->seq.l; ++i)
				ksq->seq.s[i] = ksq->seq.s[i] == 4? 4 : 3 - ksq->seq.s[i];
			q[1] = ksw_qinit(size, ksq->seq.l, (uint8_t*)ksq->seq.s, 5, mat);
		} else q[1] = 0;
		gzrewind(fpt); kseq_rewind(kst);
		while (kseq_read(kst) > 0) {
			int s;
			for (i = 0; i < kst->seq.l; ++i) kst->seq.s[i] = seq_nt4_table[(int)kst->seq.s[i]];
			s = ksw_sse2(q[0], kst->seq.l, (uint8_t*)kst->seq.s, &a);
			printf("%s\t%s\t+\t%d\t%d\t%d\n", ksq->name.s, kst->name.s, s, a.te+1, a.qe+1);
			if (q[1]) {
				s = ksw_sse2(q[1], kst->seq.l, (uint8_t*)kst->seq.s, &a);
				printf("%s\t%s\t-\t%d\t%d\t%d\n", ksq->name.s, kst->name.s, s, a.te+1, a.qe+1);
			}
		}
		free(q[0]); free(q[1]);
	}
	kseq_destroy(kst); gzclose(fpt);
	kseq_destroy(ksq); gzclose(fpq);
	return 0;
}
#endif // _KSW_MAIN
