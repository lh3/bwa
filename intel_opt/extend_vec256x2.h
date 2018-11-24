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

class EC16BitAvx2_32 {

public:
  
  typedef int16_t Word ;

  struct Vec {
    __m256i v[2]; 

    Vec(){}
    explicit Vec(__m256i v1, __m256i v0) {
      v[1] = v1 ;
      v[0] = v0 ;
    }

  } ;


  static const int MAX_VEC_LEN = 32 ;
  static const int BAD_SCORE = -1000 ;
  
  static inline Vec vec_set(Word* arr) {

    return Vec(_mm256_set_epi16(arr[31], arr[30], arr[29], arr[28], arr[27],
				arr[26], arr[25], arr[24], arr[23], arr[22], arr[21],
				arr[20], arr[19], arr[18], arr[17], arr[16]),
	       _mm256_set_epi16(arr[15], arr[14], arr[13], arr[12], arr[11],
				arr[10], arr[9], arr[8], arr[7], arr[6], arr[5],
				arr[4], arr[3], arr[2], arr[1], arr[0])) ;
  }

  static inline Vec vec_set1(Word val) {
    return Vec(_mm256_set1_epi16(val), _mm256_set1_epi16(val)) ;
  }

  static inline Vec vec_add(const Vec& v1, const Vec& v2) {
    return Vec(_mm256_add_epi16(v1.v[1], v2.v[1]), _mm256_add_epi16(v1.v[0], v2.v[0])) ;
  }

  static inline Vec vec_or(const Vec& v1, const Vec& v2) {
    return Vec(_mm256_or_si256(v1.v[1], v2.v[1]), _mm256_or_si256(v1.v[0], v2.v[0])) ;
  }

  static inline Vec vec_max(const Vec& v1, const Vec& v2) {
    return Vec(_mm256_max_epi16(v1.v[1], v2.v[1]), _mm256_max_epi16(v1.v[0], v2.v[0])) ;
  }

  static inline Vec vec_compare_eq(const Vec& v1, const Vec& v2) {
    return Vec(_mm256_cmpeq_epi16(v1.v[1], v2.v[1]), _mm256_cmpeq_epi16(v1.v[0], v2.v[0])) ;
  }

  static inline Vec vec_compare_gt(const Vec& v1, const Vec& v2) {
    return Vec(_mm256_cmpgt_epi16(v1.v[1], v2.v[1]), _mm256_cmpgt_epi16(v1.v[0], v2.v[0])) ;
  }
  
  static inline Vec vec_blend(const Vec& vFalse, const Vec& vTrue, const Vec& mask) {
    return Vec(_mm256_blendv_epi8(vFalse.v[1], vTrue.v[1], mask.v[1]), 
	       _mm256_blendv_epi8(vFalse.v[0], vTrue.v[0], mask.v[0])) ;
  }

  static inline bool vec_test_all_ones(const Vec& vec, const Vec& allOnesV) {
    return _mm256_testc_si256(vec.v[0], allOnesV.v[0]) && 
      _mm256_testc_si256(vec.v[1], allOnesV.v[0]) ;
  }


  // 0 is inserted to LSW
  static inline Vec vec_shift_left (const Vec& vec) {
    return Vec(_mm256_alignr_epi8(vec.v[1], _mm256_permute2x128_si256(vec.v[1], vec.v[0], 0x03), 14),
	       _mm256_alignr_epi8(vec.v[0], _mm256_permute2x128_si256(vec.v[0], vec.v[0], 0x04), 14)) ;
  }

  // The most significant word (MSW) in words is inserted to LSW of vec
  static inline Vec vec_shift_left_and_insert (const Vec& vec, const Vec& words) {
    return Vec(_mm256_alignr_epi8(vec.v[1], _mm256_permute2x128_si256(vec.v[1], vec.v[0], 0x03), 14),
	       _mm256_alignr_epi8(vec.v[0], _mm256_permute2x128_si256(vec.v[0], words.v[1], 0x03), 14)) ;
  }
  
  // 0 is inserted to LSW
  static inline Vec vec_shift_right (const Vec& vec) {
    return Vec(_mm256_alignr_epi8(_mm256_permute2x128_si256(vec.v[1], vec.v[1], 0x41), vec.v[1], 2),
	       _mm256_alignr_epi8(_mm256_permute2x128_si256(vec.v[0], vec.v[1], 0x21), vec.v[0], 2)) ;
  }

  // The LSW of words is inserted to the MSW of vec
  static inline Vec vec_shift_right_and_insert (const Vec& vec, const Vec& words) {
    return Vec(_mm256_alignr_epi8(_mm256_permute2x128_si256(vec.v[1], words.v[0], 0x21), vec.v[1], 2),
	       _mm256_alignr_epi8(_mm256_permute2x128_si256(vec.v[0], vec.v[1], 0x21), vec.v[0], 2));
  }

  static void vec_print(const Vec& vec) ;

} ;

void EC16BitAvx2_32::vec_print(const Vec& vec) {
  for(int i = 15; i >= 0; i--)
    cout << std::setw(5) << ((int16_t *)&vec.v[1])[i] << " " ;
  for(int i = 15; i >= 0; i--)
    cout << std::setw(5) << ((int16_t *)&vec.v[0])[i] << " " ;
  cout << endl ;
}

