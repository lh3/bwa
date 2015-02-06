/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_VEC128x2_H
#define _FAST_EXTEND_VEC128x2_H

#include <limits>

#include "fast_extend_vec.h"
#include "fast_extend_bitv128.h"
#include "fast_extend_bitv128x2.h"

#ifdef USE_AVX2
#include "fast_extend_bitv256.h"
#endif

//
// Edit distance vector for a less-than-255-base query.
// The probes are at every 16 bases. Each probe is 16 bits.
// There are 16 probes in total.
//
class EDVec128x2Every16 {
  __m128i vec_[2] ;
  static const int WORD_CNT_ = 16 ;
  static const int PERIOD_ = 16 ;

 public:
  EDVec128x2Every16() {}

  explicit EDVec128x2Every16(__m128i v0, __m128i v1) {
    vec_[0] = v0 ;
    vec_[1] = v1 ;
  }

  explicit EDVec128x2Every16(int16_t val) {
    vec_[0] = _mm_set1_epi16(val) ;
    vec_[1] = _mm_set1_epi16(val) ; 
  }

  explicit EDVec128x2Every16(const BitVec128x2& bv) {
    vec_[0] = bv.bitV[0] ;
    vec_[1] = bv.bitV[1] ;
  }

#ifdef USE_AVX2
  explicit EDVec128x2Every16(const BitVec256& bv) {
    vec_[0] = _mm256_castsi256_si128(bv.bitV) ;
    vec_[1] = _mm256_extractf128_si256(bv.bitV, 1) ;
  }
#endif

  inline void setAll(int16_t val) {
    vec_[0] = _mm_set1_epi16(val) ;    
    vec_[1] = _mm_set1_epi16(val) ;    
  }

  static int16_t getMaxWordVal() { return std::numeric_limits<int16_t>::max() ;}

  inline int getWord(int index) const {
    switch(index) {
    case 0: return _mm_extract_epi16(vec_[0], 0) ;
    case 1: return  _mm_extract_epi16(vec_[0], 1) ;
    case 2: return  _mm_extract_epi16(vec_[0], 2) ;
    case 3: return  _mm_extract_epi16(vec_[0], 3) ;
    case 4: return  _mm_extract_epi16(vec_[0], 4) ;
    case 5: return  _mm_extract_epi16(vec_[0], 5) ;
    case 6: return  _mm_extract_epi16(vec_[0], 6) ;
    case 7: return  _mm_extract_epi16(vec_[0], 7) ;
    case 8: return  _mm_extract_epi16(vec_[1], 0) ;
    case 9: return  _mm_extract_epi16(vec_[1], 1) ;
    case 10: return  _mm_extract_epi16(vec_[1], 2) ;
    case 11: return  _mm_extract_epi16(vec_[1], 3) ;
    case 12: return  _mm_extract_epi16(vec_[1], 4) ;
    case 13: return  _mm_extract_epi16(vec_[1], 5) ;
    case 14: return  _mm_extract_epi16(vec_[1], 6) ;
    case 15: return  _mm_extract_epi16(vec_[1], 7) ;
    default: assert(false); return 0 ;
    }
  }

  bool allLessThanOrEqualTo(const EDVec128x2Every16& other) const {
    // Each elt in pred will be set to 0xFFFF if the element in this vector 
    // is greater than the corresponding element in other
    __m128i pred0 = _mm_cmpgt_epi16(vec_[0], other.vec_[0]) ; 
    __m128i pred1 = _mm_cmpgt_epi16(vec_[1], other.vec_[1]) ; 
    int flags0 = _mm_movemask_epi8(pred0) ; // creates a mask from MSB of each 8-bit elt in pred
    int flags1 = _mm_movemask_epi8(pred1) ; // creates a mask from MSB of each 8-bit elt in pred

    return (flags0 | flags1) == 0 ; 
    // if there's at least one elt greater than "other", flags != 0
  }

  inline void setWords (int* v) {
    vec_[0] = _mm_set_epi16(v[7],  v[6],  v[5],  v[4],  v[3],  v[2],  v[1], v[0]) ;
    vec_[1] = _mm_set_epi16(v[15], v[14], v[13], v[12], v[11], v[10], v[9], v[8]) ;
  }

  inline void setFirstWord(int val) {
    vec_[0] = _mm_insert_epi16(vec_[0], val, 0) ;
  }

  inline void setAllWords (uint16_t val) {
    vec_[0] = _mm_set1_epi16(val) ;
    vec_[1] = _mm_set1_epi16(val) ;
  }

  inline void setWordsAsMask () {
    vec_[0] = _mm_set1_epi16(0x1) ;
    vec_[1] = _mm_set1_epi16(0x1) ;
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

  void setMin(const EDVec128x2Every16& other) {
    vec_[0] = _mm_min_epi16(vec_[0], other.vec_[0]) ;
    vec_[1] = _mm_min_epi16(vec_[1], other.vec_[1]) ;
  }

  void setMin(const EDVec128x2Every16& other, EDVec128x2Every16& bestIndices, const EDVec128x2Every16& otherIndices) {
    __m128i pred0 = _mm_cmplt_epi16(other.vec_[0], vec_[0]) ; // pred is 0xFFFF if other < this
    __m128i pred1 = _mm_cmplt_epi16(other.vec_[1], vec_[1]) ; // pred is 0xFFFF if other < this

    vec_[0] = BLENDV(vec_[0], other.vec_[0], pred0) ; // if pred is 0, choose this; o/w choose other
    vec_[1] = BLENDV(vec_[1], other.vec_[1], pred1) ; // if pred is 0, choose this; o/w choose other

    bestIndices.vec_[0] = BLENDV(bestIndices.vec_[0], otherIndices.vec_[0], pred0) ;
    bestIndices.vec_[1] = BLENDV(bestIndices.vec_[1], otherIndices.vec_[1], pred1) ;
  }

  void setMax(const EDVec128x2Every16& other, EDVec128x2Every16& bestIndices, const EDVec128x2Every16& otherIndices) {
    __m128i pred0 = _mm_cmpgt_epi16(other.vec_[0], vec_[0]) ; // pred is 0xFFFF if other > this
    __m128i pred1 = _mm_cmpgt_epi16(other.vec_[1], vec_[1]) ; // pred is 0xFFFF if other > this

    vec_[0] = BLENDV(vec_[0], other.vec_[0], pred0) ; // if pred is 0, choose this; o/w choose other
    vec_[1] = BLENDV(vec_[1], other.vec_[1], pred1) ; // if pred is 0, choose this; o/w choose other

    bestIndices.vec_[0] = BLENDV(bestIndices.vec_[0], otherIndices.vec_[0], pred0) ;
    bestIndices.vec_[1] = BLENDV(bestIndices.vec_[1], otherIndices.vec_[1], pred1) ;
  }

  void setMax(const EDVec128x2Every16& other, EDVec128x2Every16& bestIndices, const EDVec128x2Every16& otherIndices, EDVec128x2Every16& bestAccumDist, const EDVec128x2Every16& otherAccumDist) {
 
   __m128i pred0 = _mm_cmpgt_epi16(other.vec_[0], vec_[0]) ; // pred is 0xFFFF if other > this
   __m128i pred1 = _mm_cmpgt_epi16(other.vec_[1], vec_[1]) ; // pred is 0xFFFF if other > this

    vec_[0] = BLENDV(vec_[0], other.vec_[0], pred0) ; // if pred is 0, choose this; o/w choose other
    vec_[1] = BLENDV(vec_[1], other.vec_[1], pred1) ; // if pred is 0, choose this; o/w choose other

    bestIndices.vec_[0] = BLENDV(bestIndices.vec_[0], otherIndices.vec_[0], pred0) ;
    bestIndices.vec_[1] = BLENDV(bestIndices.vec_[1], otherIndices.vec_[1], pred1) ;

    bestAccumDist.vec_[0] = BLENDV(bestAccumDist.vec_[0], otherAccumDist.vec_[0], pred0) ;
    bestAccumDist.vec_[1] = BLENDV(bestAccumDist.vec_[1], otherAccumDist.vec_[1], pred1) ;
  }

  // The return value will be nonzero iff at least one of the entries is updated.
  // 2 bits per word is returned because of SSE limitations.
  uint32_t setMaxAndReturnFlag(const EDVec128x2Every16& other, EDVec128x2Every16& bestIndices, const EDVec128x2Every16& otherIndices) {
    __m128i pred0 = _mm_cmpgt_epi16(other.vec_[0], vec_[0]) ; // pred is 0xFFFF if other > this
    __m128i pred1 = _mm_cmpgt_epi16(other.vec_[1], vec_[1]) ; // pred is 0xFFFF if other > this
    
    vec_[0] = BLENDV(vec_[0], other.vec_[0], pred0) ; // if pred is 0, choose this; o/w choose other
    vec_[1] = BLENDV(vec_[1], other.vec_[1], pred1) ; // if pred is 0, choose this; o/w choose other
    
    bestIndices.vec_[0] = BLENDV(bestIndices.vec_[0], otherIndices.vec_[0], pred0) ;
    bestIndices.vec_[1] = BLENDV(bestIndices.vec_[1], otherIndices.vec_[1], pred1) ;
    
    return (((uint32_t) _mm_movemask_epi8(pred1)) << 16) | _mm_movemask_epi8(pred0) ;
  }
  
  // The return value will be nonzero iff at least one of the entries is updated.
  // 2 bits per word is returned because of SSE limitations.
  uint32_t setMaxAndReturnFlag(const EDVec128x2Every16& other, EDVec128x2Every16& bestIndices, const EDVec128x2Every16& otherIndices, EDVec128x2Every16& bestAccumDist, const EDVec128x2Every16& otherAccumDist) {
    __m128i pred0 = _mm_cmpgt_epi16(other.vec_[0], vec_[0]) ; // pred is 0xFFFF if other > this
    __m128i pred1 = _mm_cmpgt_epi16(other.vec_[1], vec_[1]) ; // pred is 0xFFFF if other > this

    vec_[0] = BLENDV(vec_[0], other.vec_[0], pred0) ; // if pred is 0, choose this; o/w choose other
    vec_[1] = BLENDV(vec_[1], other.vec_[1], pred1) ; // if pred is 0, choose this; o/w choose other

    bestIndices.vec_[0] = BLENDV(bestIndices.vec_[0], otherIndices.vec_[0], pred0) ;
    bestIndices.vec_[1] = BLENDV(bestIndices.vec_[1], otherIndices.vec_[1], pred1) ;

    bestAccumDist.vec_[0] = BLENDV(bestAccumDist.vec_[0], otherAccumDist.vec_[0], pred0) ;
    bestAccumDist.vec_[1] = BLENDV(bestAccumDist.vec_[1], otherAccumDist.vec_[1], pred1) ;
    
    return (((uint32_t) _mm_movemask_epi8(pred1)) << 16) | _mm_movemask_epi8(pred0) ;
  }

  void addThirdIfFirstGTSecond(const EDVec128x2Every16& first, const EDVec128x2Every16& second, 
			       const EDVec128x2Every16& third, const EDVec128x2Every16& zero) {
    
    __m128i pred0 = _mm_cmpgt_epi16(first.vec_[0], second.vec_[1]) ; // pred is 0xFFFF if first > second
    __m128i pred1 = _mm_cmpgt_epi16(first.vec_[0], second.vec_[1]) ; // pred is 0xFFFF if first > second

    // If pred is 0, then blendv will add zero; o/w it will add third
    vec_[0] = _mm_add_epi16(vec_[0], BLENDV(zero.vec_[0], third.vec_[0], pred0)) ;
    vec_[1] = _mm_add_epi16(vec_[1], BLENDV(zero.vec_[1], third.vec_[1], pred1)) ;
  }
 
  EDVec128x2Every16 subSat(const EDVec128x2Every16& other) const {
    return EDVec128x2Every16(_mm_subs_epu16(this->vec_[0], other.vec_[0]),
			     _mm_subs_epu16(this->vec_[1], other.vec_[1])) ;
  }

  inline void shiftWordsLeftByOne() {
    vec_[1] = _mm_alignr_epi8(vec_[1], vec_[0], 14) ;
    vec_[0] = _mm_slli_si128(vec_[0], 2) ; // 16-bit shift = 2 * 8-bit shift
  }
  
  inline bool operator == (const EDVec128x2Every16& other) {
    return (_mm_movemask_epi8(_mm_cmpeq_epi16(vec_[0], other.vec_[0])) &
	    _mm_movemask_epi8(_mm_cmpeq_epi16(vec_[1], other.vec_[1]))) == 0xFFFF ;
  }
  
  inline EDVec128x2Every16& operator += (const EDVec128x2Every16& other) {
    vec_[0] = _mm_add_epi16(vec_[0], other.vec_[0]) ;
    vec_[1] = _mm_add_epi16(vec_[1], other.vec_[1]) ;
    return *this ;
  }

  inline EDVec128x2Every16& operator -= (const EDVec128x2Every16& other) {
    vec_[0] = _mm_sub_epi16(vec_[0], other.vec_[0]) ;
    vec_[1] = _mm_sub_epi16(vec_[1], other.vec_[1]) ;
    return *this ;
  }

  inline EDVec128x2Every16 operator+ (const EDVec128x2Every16& other) const {
    return EDVec128x2Every16(_mm_add_epi16(vec_[0], other.vec_[0]),
			     _mm_add_epi16(vec_[1], other.vec_[1])) ;
  }

  inline EDVec128x2Every16 operator- (const EDVec128x2Every16& other) const {
    return EDVec128x2Every16(_mm_sub_epi16(vec_[0], other.vec_[0]),
			     _mm_sub_epi16(vec_[1], other.vec_[1])) ;
  }

  inline EDVec128x2Every16 operator* (const EDVec128x2Every16& other) const {
    return EDVec128x2Every16(_mm_mullo_epi16(vec_[0], other.vec_[0]),
			     _mm_mullo_epi16(vec_[1], other.vec_[1])) ;
  }

  inline EDVec128x2Every16 operator& (const EDVec128x2Every16& other) const {
    return EDVec128x2Every16(_mm_and_si128(vec_[0], other.vec_[0]),
			     _mm_and_si128(vec_[1], other.vec_[1])) ;
  }

  static int WORD_CNT() {
    return WORD_CNT_ ;
  }

  static int PERIOD() {
    return PERIOD_ ;
  }

  friend std::ostream& operator<< (std::ostream& os, const EDVec128x2Every16& v) {
    os << "{" ;
    for (int i=0; i < WORD_CNT_; ++i) {
      os << v.getWord(i) << " " ;
    }
    os << "}" ;
    return os ;
  }

} ;



#endif // #ifndef _FAST_EXTEND_VEC128x2_H
