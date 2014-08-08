/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_VEC_H
#define _FAST_EXTEND_VEC_H

#include <immintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#ifdef USE_SSE4
#include <tmmintrin.h>
#include <smmintrin.h>

#define BLENDV(__v1, __v2, __pred)		\
  _mm_blendv_epi8(__v1, __v2, __pred)

#else // #ifdef USE_SSE4

#define BLENDV(__v1, __v2, __pred)					\
  _mm_or_si128(_mm_andnot_si128(__pred, __v1), _mm_and_si128(__v2, __pred))

#endif // #ifdef USE_SSE4

using std::cout ;
using std::endl ;

inline void printVec16(const __m128i& v) {
  
  cout << "{" ;
  cout << " " << _mm_extract_epi16(v, 0) ;
  cout << " " << _mm_extract_epi16(v, 1) ;
  cout << " " << _mm_extract_epi16(v, 2) ;
  cout << " " << _mm_extract_epi16(v, 3) ;
  cout << " " <<_mm_extract_epi16(v, 4) ;
  cout << " " <<_mm_extract_epi16(v, 5) ;
  cout << " " <<_mm_extract_epi16(v, 6) ;
  cout << " " <<_mm_extract_epi16(v, 7) ;

  cout << "}" << endl ;
}



#endif // #ifdef _FAST_EXTEND_VEC_H
