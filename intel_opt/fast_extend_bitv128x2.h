/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_BITV128x2_H
#define _FAST_EXTEND_BITV128x2_H

#include <stdint.h>
#include <stdlib.h>
#include <iostream>

using namespace std ;

class BitVec128x2 {

  __m128i bitV[2] ;

public:
  
  BitVec128x2() {
    bitV[0] = _mm_setzero_si128() ;
    bitV[1] = _mm_setzero_si128() ;
  }

  explicit BitVec128x2(const __m128i& v0, const __m128i& v1) {
    bitV[0] = v0 ;
    bitV[1] = v1 ;
  }

  inline uint16_t getLow16Bits() const {
    return (uint16_t) _mm_extract_epi16(bitV[0], 0) ;
  }

  inline void setAll64BitWords(uint64_t val) {
    bitV[0] = _mm_set1_epi64x(val) ;
    bitV[1] = _mm_set1_epi64x(val) ;
  }

  inline void setAllOnes () {
    bitV[0] = _mm_set1_epi8(0xFF) ;
    bitV[1] = _mm_set1_epi8(0xFF) ;
  }

  inline void setAllZeroes () {
    bitV[0] = _mm_setzero_si128() ;
    bitV[1] = _mm_setzero_si128() ;
  }

  inline void setAllBits (bool bit) {
    if (bit) setAllOnes() ;
    else setAllZeroes() ;
  }

  inline void setLSBClearRest (bool bit) {
    bitV[0] = _mm_set_epi64x(0x0, bit) ;
    bitV[1] = _mm_set_epi64x(0x0, 0x0) ;
  }

  inline void shiftLeftAndInsert (bool newBit) {

    static const __m128i lsbSetVec = _mm_set_epi64x(0x0, 0x1) ; // only the 0th bit is 1, rest are 0
    
    // Shift the MSBs to LSB positions
    __m128i msb2lsb0 = _mm_srli_epi64(bitV[0], 63) ;
    __m128i msb2lsb1 = _mm_srli_epi64(bitV[1], 63) ;
    
    bitV[0] = _mm_or_si128(_mm_slli_epi64(bitV[0], 1),
			   _mm_slli_si128(msb2lsb0, 8)) ;
    
    bitV[1] = _mm_or_si128(_mm_slli_epi64(bitV[1], 1),
			   _mm_alignr_epi8(msb2lsb1, msb2lsb0, 8)) ;

    if (newBit)
      bitV[0] = _mm_or_si128(bitV[0], lsbSetVec) ;
  }

  inline void shiftLeft () {
    // Shift the MSBs to LSB positions
    __m128i msb2lsb0 = _mm_srli_epi64(bitV[0], 63) ;
    __m128i msb2lsb1 = _mm_srli_epi64(bitV[1], 63) ;
    
    bitV[0] = _mm_or_si128(_mm_slli_epi64(bitV[0], 1),
			   _mm_slli_si128(msb2lsb0, 8)) ;
    
    bitV[1] = _mm_or_si128(_mm_slli_epi64(bitV[1], 1),
			   _mm_alignr_epi8(msb2lsb1, msb2lsb0, 8)) ;
  }

  inline BitVec128x2 shiftProbesRight(int probeOffset, const BitVec128x2& shiftCnt) const {
    
    BitVec128x2 r ;
    if (probeOffset == 15) {
      // shift right by 16 bits (across 128-bit boundary)
      r.bitV[0] = _mm_or_si128(_mm_srli_si128(bitV[0], 2),
			       _mm_slli_si128(bitV[1], 14)) ;
      r.bitV[1] = _mm_srli_si128(bitV[1], 2) ;
    }
    else {
      // shift right by shiftCnt within each 16-bit word
      r.bitV[0] = _mm_srl_epi64(bitV[0], shiftCnt.bitV[0]) ;
      r.bitV[1] = _mm_srl_epi64(bitV[1], shiftCnt.bitV[1]) ;
    }

    return r ;
  }

  void setWord (int index, uint64_t w) {
    switch(index) {
    case 0: bitV[0] = _mm_insert_epi64(bitV[0], w, 0) ; break ;
    case 1: bitV[0] = _mm_insert_epi64(bitV[0], w, 1) ; break ;
    case 2: bitV[1] = _mm_insert_epi64(bitV[1], w, 0) ; break ;
    case 3: bitV[1] = _mm_insert_epi64(bitV[1], w, 1) ; break ;
    default: assert(false) ;
    }
  }

  inline BitVec128x2 operator~ () const {
    static const __m128i onesVec = _mm_set1_epi8(0xFF) ;

    return BitVec128x2(_mm_andnot_si128(bitV[0], onesVec),
		       _mm_andnot_si128(bitV[1], onesVec)) ;
  }

  inline BitVec128x2 operator& (const BitVec128x2& other) const {

    return BitVec128x2(_mm_and_si128(bitV[0], other.bitV[0]),
		       _mm_and_si128(bitV[1], other.bitV[1])) ;
  }

  inline BitVec128x2 andnot (const BitVec128x2& other) const {

    return BitVec128x2(_mm_andnot_si128(other.bitV[0], bitV[0]),
		       _mm_andnot_si128(other.bitV[1], bitV[1])) ;
  }

  inline BitVec128x2 operator| (const BitVec128x2& other) const {

    return BitVec128x2(_mm_or_si128(bitV[0], other.bitV[0]),
		       _mm_or_si128(bitV[1], other.bitV[1])) ;
  }

  inline BitVec128x2 operator^ (const BitVec128x2& other) const {

    return BitVec128x2(_mm_xor_si128(bitV[0], other.bitV[0]),
		       _mm_xor_si128(bitV[1], other.bitV[1])) ;
  }

  inline BitVec128x2 operator+ (const BitVec128x2& other) const {

    BitVec128x2 sum ;
    sum.bitV[0] = _mm_add_epi64(this->bitV[0], other.bitV[0]) ;
    sum.bitV[1] = _mm_add_epi64(this->bitV[1], other.bitV[1]) ;
    
    // We should add the carryout from the lower word (if any) to upper word

    uint64_t sum0 = _mm_extract_epi64(sum.bitV[0], 0) ;
    uint64_t this0 = _mm_extract_epi64(this->bitV[0], 0) ;
    uint carry0 = (sum0 < this0) ? 1 : 0 ;

    uint64_t sum1 = _mm_extract_epi64(sum.bitV[0], 1) ;
    uint64_t this1 = _mm_extract_epi64(this->bitV[0], 1) ;
    uint carry1 = (sum1 < this1) ? 1 : 0 ;

    uint64_t sum2 = _mm_extract_epi64(sum.bitV[1], 0) ;
    uint64_t this2 = _mm_extract_epi64(this->bitV[1], 0) ;
    uint carry2 = (sum2 < this2) ? 1 : 0 ;

    uint64_t sum3 = _mm_extract_epi64(sum.bitV[1], 1) ;

    sum1 += carry0 ;
    sum2 += carry1 ;
    sum3 += carry2 ;

    carry2 = carry1 && (sum2 == 0) ;
    carry1 = carry0 && (sum1 == 0) ;

    sum2 += carry1 ;
    sum3 += carry2 ;
    
    carry2 = carry1 && (sum2 == 0) ;
    sum3 += carry2 ;

    sum.bitV[0] = _mm_set_epi64x(sum1, sum0) ;
    sum.bitV[1] = _mm_set_epi64x(sum3, sum2) ;

    return sum ;
  }

  static bool isQueryLengthOk(int qlen) {
    return qlen <= 255 ;
  }

  friend std::ostream& operator<< (std::ostream& os, const BitVec128x2& b) {
    os << "{" << _mm_extract_epi64(b.bitV[0], 0) << " " << _mm_extract_epi64(b.bitV[0], 1) << " "
       << _mm_extract_epi64(b.bitV[1], 0) << " " << _mm_extract_epi64(b.bitV[1], 1) << "}" ;
    return os ;
  }

  friend class EDVec128x2Every16 ;
  friend class EDVec256Every16 ;
  
} ;



#endif // #ifndef _FAST_EXTEND_BITV128_H
