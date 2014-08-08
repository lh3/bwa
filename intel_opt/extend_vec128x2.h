/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <immintrin.h>

#include <iostream>
#include <iomanip>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

using std::cout ;
using std::endl ;

class EC16BitSSE_16 {
  
public:
  
  typedef int16_t Word ;

  struct Vec {
    __m128i v[2]; 
    
    Vec(){}
    explicit Vec(__m128i v1, __m128i v0) {
      v[1] = v1 ;
      v[0] = v0 ;
    }
  } ;


  static const int MAX_VEC_LEN = 16 ;
  static const int BAD_SCORE = -1000 ;
  
  static inline Vec vec_set(Word* arr) {
    return Vec(_mm_set_epi16(arr[15], arr[14], arr[13], arr[12], arr[11], arr[10], arr[9], arr[8]),
	       _mm_set_epi16(arr[7], arr[6], arr[5], arr[4], arr[3], arr[2], arr[1], arr[0])) ;
  }

  static inline Vec vec_set1(Word val) {
    return Vec(_mm_set1_epi16(val), _mm_set1_epi16(val)) ;
  }

  static inline Vec vec_add(const Vec& v1, const Vec& v2) {
    return Vec(_mm_add_epi16(v1.v[1], v2.v[1]), 
	       _mm_add_epi16(v1.v[0], v2.v[0])) ;
  }

  static inline Vec vec_or(const Vec& v1, const Vec& v2) {
    return Vec(_mm_or_si128(v1.v[1], v2.v[1]),
	       _mm_or_si128(v1.v[0], v2.v[0])) ;
  }

  static inline Vec vec_max(const Vec& v1, const Vec& v2) {
    return Vec(_mm_max_epi16(v1.v[1], v2.v[1]), 
	       _mm_max_epi16(v1.v[0], v2.v[0])) ;
  }

  static inline Vec vec_compare_eq(const Vec& v1, const Vec& v2) {
    return Vec(_mm_cmpeq_epi16(v1.v[1], v2.v[1]), 
	       _mm_cmpeq_epi16(v1.v[0], v2.v[0])) ;
  }

  static inline Vec vec_compare_gt(const Vec& v1, const Vec& v2) {
    return Vec(_mm_cmpgt_epi16(v1.v[1], v2.v[1]), 
	       _mm_cmpgt_epi16(v1.v[0], v2.v[0])) ;
  }
  
  static inline Vec vec_blend(const Vec& vFalse, const Vec& vTrue, const Vec& mask) {
    return Vec(_mm_blendv_epi8(vFalse.v[1], vTrue.v[1], mask.v[1]), 
	       _mm_blendv_epi8(vFalse.v[0], vTrue.v[0], mask.v[0])) ;
  }

  static inline bool vec_test_all_ones(const Vec& vec, const Vec& allOnesV) {
    return _mm_test_all_ones(vec.v[1]) && _mm_test_all_ones(vec.v[0]) ;
  }

  // 0 is inserted to LSW
  static inline Vec vec_shift_left (const Vec& vec) {
    return Vec(_mm_alignr_epi8(vec.v[1], vec.v[0], 14),
	       _mm_slli_si128(vec.v[0], 2)) ;
  }

  // The most significant word (MSW) in words is inserted to LSW of vec
  static inline Vec vec_shift_left_and_insert (const Vec& vec, const Vec& words) {
    return Vec(_mm_alignr_epi8(vec.v[1], vec.v[0], 14),
	       _mm_alignr_epi8(vec.v[0], words.v[1], 14)) ;
  }

  // 0 is inserted to MSW
  static inline Vec vec_shift_right (const Vec& vec) {
    return Vec(_mm_srli_si128(vec.v[1], 2),
	       _mm_alignr_epi8(vec.v[1], vec.v[0], 2)) ;
  }
  
  // The LSW of words is inserted to the MSW of vec
  static inline Vec vec_shift_right_and_insert (const Vec& vec, const Vec& words) {
    return Vec(_mm_alignr_epi8(words.v[0], vec.v[1], 2),
	       _mm_alignr_epi8(vec.v[1], vec.v[0], 2)) ;
  }

  static void vec_print(const Vec& vec) ;

} ;

void EC16BitSSE_16::vec_print(const Vec& vec) {
  for(int i = 7; i >= 0; i--)
    cout << std::setw(5) << ((const int16_t*) (&vec.v[1]))[i] << " " ;

  for(int i = 7; i >= 0; i--)
    cout << std::setw(5) << ((const int16_t*) (&vec.v[0]))[i] << " " ;

  cout << endl ;
}

