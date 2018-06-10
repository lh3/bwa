/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_VEC256_H
#define _FAST_EXTEND_VEC256_H

#include <limits>

#include "fast_extend_vec.h"
#include "fast_extend_bitv128.h"
#include "fast_extend_bitv128x2.h"
#include "fast_extend_bitv256.h"

//
// Edit distance vector for a less-than-255-base query.
// The probes are at every 16 bases. Each probe is 16 bits.
// There are 16 probes in total.
//
class EDVec256Every16 {
  __m256i vec_ ;
  static const int WORD_CNT_ = 16 ;
  static const int PERIOD_ = 16 ;

 public:
  EDVec256Every16() {}

  explicit EDVec256Every16(int16_t val) {
    vec_ = _mm256_set1_epi16(val) ;
  }

  // gcc complains about undefined _mm256_set_m128i function. 
  // Commented out this adaptor because it's not really needed.
  //explicit EDVec256Every16(const BitVec128x2& bv) {
  //  vec_ = _mm256_set_m128i(bv.bitV[1], bv.bitV[0]) ;
  //}

#ifdef USE_AVX2
  explicit EDVec256Every16(__m256i v) {
    vec_ = v ;
  }

  explicit EDVec256Every16(const BitVec256& bv) {
    vec_ = bv.bitV ;
  }
#endif

  inline void setAll(int16_t val) {
    vec_ = _mm256_set1_epi16(val) ;
  }

  static int16_t getMaxWordVal() { return std::numeric_limits<int16_t>::max() ;}

  inline int getWord(int index) const {
    switch(index) {
    case 0: return _mm_extract_epi16(_mm256_castsi256_si128(vec_), 0) ;
    case 1: return _mm_extract_epi16(_mm256_castsi256_si128(vec_), 1) ;
    case 2: return _mm_extract_epi16(_mm256_castsi256_si128(vec_), 2) ;
    case 3: return _mm_extract_epi16(_mm256_castsi256_si128(vec_), 3) ;
    case 4: return _mm_extract_epi16(_mm256_castsi256_si128(vec_), 4) ;
    case 5: return _mm_extract_epi16(_mm256_castsi256_si128(vec_), 5) ;
    case 6: return _mm_extract_epi16(_mm256_castsi256_si128(vec_), 6) ;
    case 7: return _mm_extract_epi16(_mm256_castsi256_si128(vec_), 7) ;

    case  8: return _mm_extract_epi16(_mm256_extractf128_si256(vec_,1), 0) ;
    case  9: return _mm_extract_epi16(_mm256_extractf128_si256(vec_,1), 1) ;
    case 10: return _mm_extract_epi16(_mm256_extractf128_si256(vec_,1), 2) ;
    case 11: return _mm_extract_epi16(_mm256_extractf128_si256(vec_,1), 3) ;
    case 12: return _mm_extract_epi16(_mm256_extractf128_si256(vec_,1), 4) ;
    case 13: return _mm_extract_epi16(_mm256_extractf128_si256(vec_,1), 5) ;
    case 14: return _mm_extract_epi16(_mm256_extractf128_si256(vec_,1), 6) ;
    case 15: return _mm_extract_epi16(_mm256_extractf128_si256(vec_,1), 7) ;

    default: assert(false); return 0 ;
    }
  }

  bool allLessThanOrEqualTo(const EDVec256Every16& other) const {
    // Each elt in pred will be set to 0xFFFF if the element in this vector 
    // is greater than the corresponding element in other
    __m256i pred = _mm256_cmpgt_epi16(vec_, other.vec_) ; 

    int flags = _mm256_movemask_epi8(pred) ; // creates a mask from MSB of each 8-bit elt in pred

    return flags == 0 ;
    // if there's at least one elt greater than "other", flags != 0
  }

  inline void setWords (int* v) {
    vec_ = _mm256_set_epi16(v[15], v[14], v[13], v[12], v[11], v[10], v[9], v[8],
			    v[7],  v[6],  v[5],  v[4],  v[3],  v[2],  v[1], v[0]) ;
  }

  inline void setFirstWord(int val) {
    vec_ = _mm256_insertf128_si256(vec_,
				   _mm_insert_epi16(_mm256_castsi256_si128(vec_), val, 0),
				   0) ;
  }

  inline void setAllWords (uint16_t val) {
    vec_ = _mm256_set1_epi16(val) ;
  }

  inline void setWordsAsMask () {
    vec_ = _mm256_set1_epi16(0x1) ;
  }

  static int getLastWordIndexFor (int queryLen) {
    int baseIndex = queryLen - 1 ;
    return baseIndex/PERIOD_ ;
  }

  static int getProbeOffsetFor (int queryLen) {
    return (queryLen-1) % PERIOD_ ;
  }

  void setWordsAsEndBonus(int queryLen, int endBonus) {
    int v[WORD_CNT_] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    v[lastWordIndex] = endBonus ;

    setWords(v) ;
  }

  void setWordsAsBadScore(int queryLen, int thr, int infScore) {
    int v[WORD_CNT_] = {thr, thr, thr, thr, thr, thr, thr, thr, thr, thr, thr, thr, thr, thr, thr, thr} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    for (int i=lastWordIndex+1; i < WORD_CNT_; ++i) {
      v[i] = infScore ;
    }

    setWords(v) ;
  }

  static int getDistAt(int wordIndex, int queryLen) {
    return 1 + wordIndex * PERIOD_ + getProbeOffsetFor(queryLen) ;
  }

  void setWordsAsDist (int queryLen, int distWt) {
    int v[WORD_CNT_] ;
    int bitOffset = getProbeOffsetFor(queryLen) ;

    for (int wi=0; wi < WORD_CNT_; ++wi) {
      v[wi] = distWt * (1 + wi * PERIOD_ + bitOffset) ;
    }
    setWords(v) ;

  }

  void setMin(const EDVec256Every16& other) {
    vec_ = _mm256_min_epi16(vec_, other.vec_) ;
  }


  void setMax(const EDVec256Every16& other, EDVec256Every16& bestIndices, const EDVec256Every16& otherIndices) {
    __m256i pred = _mm256_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this

    vec_ = _mm256_blendv_epi8(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other

    bestIndices.vec_ = _mm256_blendv_epi8(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec256Every16& other, EDVec256Every16& bestIndices, const EDVec256Every16& otherIndices, EDVec256Every16& bestAccumDist, const EDVec256Every16& otherAccumDist) {
 
   __m256i pred = _mm256_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this

    vec_ = _mm256_blendv_epi8(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other

    bestIndices.vec_ = _mm256_blendv_epi8(bestIndices.vec_, otherIndices.vec_, pred) ;

    bestAccumDist.vec_ = _mm256_blendv_epi8(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;
  }

  // The return value will be nonzero iff at least one of the entries is updated.
  // 2 bits per word is returned because of SSE limitations.
  uint32_t setMaxAndReturnFlag(const EDVec256Every16& other, EDVec256Every16& bestIndices, const EDVec256Every16& otherIndices) {
    __m256i pred = _mm256_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    
    vec_ = _mm256_blendv_epi8(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    
    bestIndices.vec_ = _mm256_blendv_epi8(bestIndices.vec_, otherIndices.vec_, pred) ;
    
    return (uint32_t) _mm256_movemask_epi8(pred) ;
  }
  
  // The return value will be nonzero iff at least one of the entries is updated.
  // 2 bits per word is returned because of SSE limitations.
  uint32_t setMaxAndReturnFlag(const EDVec256Every16& other, EDVec256Every16& bestIndices, const EDVec256Every16& otherIndices, EDVec256Every16& bestAccumDist, const EDVec256Every16& otherAccumDist) {
    __m256i pred = _mm256_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this

    vec_ = _mm256_blendv_epi8(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other

    bestIndices.vec_ = _mm256_blendv_epi8(bestIndices.vec_, otherIndices.vec_, pred) ;

    bestAccumDist.vec_ = _mm256_blendv_epi8(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;
    
    return ((uint32_t) _mm256_movemask_epi8(pred)) ;
  }

  void addThirdIfFirstGTSecond(const EDVec256Every16& first, const EDVec256Every16& second, 
			       const EDVec256Every16& third, const EDVec256Every16& zero) {
    
    __m256i pred = _mm256_cmpgt_epi16(first.vec_, second.vec_) ; // pred is 0xFFFF if first > second

    // If pred is 0, then blendv will add zero; o/w it will add third
    vec_ = _mm256_add_epi16(vec_, _mm256_blendv_epi8(zero.vec_, third.vec_, pred)) ;
  }
 
  EDVec256Every16 subSat(const EDVec256Every16& other) const {
    return EDVec256Every16(_mm256_subs_epu16(this->vec_, other.vec_)) ;
  }

  inline void shiftWordsLeftByOne() {
    // Shift the vector by 2 bytes (16-bit word) across 128-bit boundary
    vec_ = _mm256_or_si256(_mm256_slli_si256(vec_, 2),
			   _mm256_permute4x64_epi64(_mm256_srli_si256(vec_,14), 0x45)) ;
  }
  
  inline bool operator == (const EDVec256Every16& other) {
    return (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi16(vec_, other.vec_)) == 0xFFFFFFFF ;
  }
  
  inline EDVec256Every16& operator += (const EDVec256Every16& other) {
    vec_ = _mm256_add_epi16(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec256Every16& operator -= (const EDVec256Every16& other) {
    vec_ = _mm256_sub_epi16(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec256Every16 operator+ (const EDVec256Every16& other) const {
    return EDVec256Every16(_mm256_add_epi16(vec_, other.vec_)) ;
  }

  inline EDVec256Every16 operator- (const EDVec256Every16& other) const {
    return EDVec256Every16(_mm256_sub_epi16(vec_, other.vec_)) ;
  }

  inline EDVec256Every16 operator* (const EDVec256Every16& other) const {
    return EDVec256Every16(_mm256_mullo_epi16(vec_, other.vec_)) ;
  }

  inline EDVec256Every16 operator& (const EDVec256Every16& other) const {
    return EDVec256Every16(_mm256_and_si256(vec_, other.vec_)) ;
  }

  static int WORD_CNT() {
    return WORD_CNT_ ;
  }

  static int PERIOD() {
    return PERIOD_ ;
  }

  friend std::ostream& operator<< (std::ostream& os, const EDVec256Every16& v) {
    os << "{" ;
    for (int i=0; i < WORD_CNT_; ++i) {
      os << v.getWord(i) << " " ;
    }
    os << "}" ;
    return os ;
  }

} ;



#endif // #ifndef _FAST_EXTEND_VEC256_H
