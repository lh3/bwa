/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_BITV128_H
#define _FAST_EXTEND_BITV128_H

#include <stdint.h>
#include <stdlib.h>
#include <iostream>

using namespace std ;

class BitVec128 {

  __m128i bitV ;

public:
  
  BitVec128():  bitV(_mm_setzero_si128()) {}

  explicit BitVec128(const __m128i& v): bitV(v) {}

  inline uint16_t getLow16Bits() const {
    return (uint16_t) _mm_extract_epi16(bitV, 0) ;
  }

  inline void setAll64BitWords(uint64_t val) {
    bitV = _mm_set1_epi64x(val) ;
  }

  inline void setAllOnes () {
    bitV = _mm_set1_epi8(0xFF) ;
  }

  inline void setAllZeroes () {
    bitV = _mm_setzero_si128() ;
  }

  inline void setAllBits (bool bit) {
    if (bit) setAllOnes() ;
    else setAllZeroes() ;
  }

  inline void setLSBClearRest (bool bit) {
    bitV = _mm_set_epi64x(0x0, bit) ;
  }

  inline void shiftLeftAndInsert (bool newBit) {

    static const __m128i lsbSetVec = _mm_set_epi64x(0x0, 0x1) ; // only the 0th bit is 1, rest are 0
    
    bitV = _mm_or_si128(_mm_slli_epi64(bitV, 1),
			_mm_srli_epi64(_mm_slli_si128(bitV, 8), 63)) ;
    if (newBit)
      bitV = _mm_or_si128(bitV, lsbSetVec) ;

  }

  inline void shiftLeft () {
    bitV = _mm_or_si128(_mm_slli_epi64(bitV, 1),
			_mm_srli_epi64(_mm_slli_si128(bitV, 8), 63)) ;

  }

 
  inline BitVec128 shiftProbesRight(int probeOffset, const BitVec128& shiftCnt) const {
    
    BitVec128 r ;
    if (probeOffset == 15)
      r.bitV = _mm_srli_si128(bitV, 2) ; // shift right by 16 bits (across 64-bit boundary)
    else
      r.bitV = _mm_srl_epi64(bitV, shiftCnt.bitV) ;

    return r ;
  }

  void setWord (int index, uint64_t w) {
    if (index == 0)
      bitV = _mm_insert_epi64(bitV, w, 0) ;
    else
      bitV = _mm_insert_epi64(bitV, w, 1) ;
  }

  inline BitVec128 operator~ () const {
    static const __m128i onesVec = _mm_set1_epi8(0xFF) ;

    return BitVec128(_mm_andnot_si128(bitV, onesVec)) ;
  }

  inline BitVec128 operator& (const BitVec128& other) const {

    return BitVec128(_mm_and_si128(bitV, other.bitV)) ;
  }

  inline BitVec128 andnot (const BitVec128& other) const {

    return BitVec128(_mm_andnot_si128(other.bitV, bitV)) ;
  }

  inline BitVec128 operator| (const BitVec128& other) const {

    return BitVec128(_mm_or_si128(bitV, other.bitV)) ;
  }

  inline BitVec128 operator^ (const BitVec128& other) const {

    return BitVec128(_mm_xor_si128(bitV, other.bitV)) ;
  }

  inline BitVec128 operator+ (const BitVec128& other) const {


    BitVec128 sum ;
    sum.bitV = _mm_add_epi64(this->bitV, other.bitV) ;
    
    // We should add the carryout from the lower word (if any) to upper word

    /*
    static __m128i lsbSetVec = _mm_set1_epi64x(0x1) ;
    static __m128i msbSetVec = _mm_set1_epi64x(0x8000000000000000) ;

    // set pred to all 1s (0xFFFFFF...) in case of overflow
    __m128i pred = _mm_cmpgt_epi64(_mm_xor_si128(this->bitV, msbSetVec),
				   _mm_xor_si128(sum.bitV, msbSetVec)) ;
    // Note: Since SSE does not have unsigned comparison opn, we first xor the MSBs,
    // and then do the signed comparison

    // If lower 64-bit of pred is set, add 0x1 to the upper word
    sum.bitV = _mm_add_epi64(sum.bitV, _mm_and_si128(lsbSetVec, _mm_slli_si128(pred, 8))) ;
    */
    /* Alternative implementation: */
    
    uint64_t sumLower = _mm_extract_epi64(sum.bitV, 0) ;
    uint64_t thisLower = _mm_extract_epi64(this->bitV, 0) ;
    if (sumLower < thisLower) {
      static __m128i carryVec = _mm_set_epi64x(0x1, 0x0) ;
      sum.bitV = _mm_add_epi64(sum.bitV, carryVec) ;
    }
    

    return sum ;
  }

  static bool isQueryLengthOk(int qlen) {
    return qlen <= 127 ;
  }

  friend std::ostream& operator<< (std::ostream& os, const BitVec128& b) {
    os << "{" << _mm_extract_epi64(b.bitV, 0) << " " << _mm_extract_epi64(b.bitV, 1) << "}" ;
    return os ;
  }

  friend class EDVec128Every16 ;

} ;


#endif // #ifndef _FAST_EXTEND_BITV128_H
