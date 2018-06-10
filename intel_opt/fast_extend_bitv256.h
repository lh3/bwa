/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_BITV256_H
#define _FAST_EXTEND_BITV256_H

#include <immintrin.h>

#include <stdint.h>
#include <stdlib.h>
#include <iostream>

using namespace std ;

class BitVec256 {

  __m256i bitV ;

public:
  
  BitVec256() {
    bitV = _mm256_setzero_si256() ;
  }

  explicit BitVec256(const __m256i& v) {
    bitV = v ;
  }

  inline uint16_t getLow16Bits() const {
    return (uint16_t) _mm_extract_epi16(_mm256_castsi256_si128(bitV), 0) ;
  }

  inline void setAll64BitWords(uint64_t val) {
    bitV = _mm256_set1_epi64x(val) ;
  }

  inline void setAllOnes () {
    bitV = _mm256_set1_epi8(0xFF) ;
  }

  inline void setAllZeroes () {
    bitV = _mm256_setzero_si256() ;
  }

  inline void setAllBits (bool bit) {
    if (bit) setAllOnes() ;
    else setAllZeroes() ;
  }

  inline void setLSBClearRest (bool bit) {
    bitV = _mm256_set_epi64x(0x0, 0x0, 0x0, bit) ;
  }

  inline void shiftLeftAndInsert (bool newBit) {

    static const __m256i lsbSetVec = _mm256_set_epi64x(0x0, 0x0, 0x0, 0x1) ; // only the 0th bit is 1, rest are 0
    static const __m256i lsbClearVec = _mm256_set_epi64x(0xFFFFFFFFFFFFFFFF, 
							 0xFFFFFFFFFFFFFFFF, 
							 0xFFFFFFFFFFFFFFFF, 
							 0xFFFFFFFFFFFFFFFE) ; // only the 0th bit is 0, rest are 1
    
    // Shift the MSBs to LSB positions
    __m256i msb2lsb = _mm256_srli_epi64(bitV, 63) ;

    // Permute 64-bit words as follows to shift words left
    // v[0] = v[0] ;
    // v[1] = v[0] ;
    // v[2] = v[1] ;
    // v[3] = v[2] ;
    // Selector bits: 10 01 00 00 = 0x90
    bitV = _mm256_or_si256(_mm256_slli_epi64(bitV, 1),
			   _mm256_permute4x64_epi64(msb2lsb, 0x90)) ;

    // The LSB needs to be set or cleared based on input
    if (newBit)
      bitV = _mm256_or_si256(bitV, lsbSetVec) ;
    else
      bitV = _mm256_and_si256(bitV, lsbClearVec) ;
  }

  inline void shiftLeft () {

    static const __m256i lsbClearVec = _mm256_set_epi64x(0xFFFFFFFFFFFFFFFF, 
							 0xFFFFFFFFFFFFFFFF, 
							 0xFFFFFFFFFFFFFFFF, 
							 0xFFFFFFFFFFFFFFFE) ; // only the 0th bit is 0, rest are 1
    
    // Shift the MSBs to LSB positions
    __m256i msb2lsb = _mm256_srli_epi64(bitV, 63) ;

    // Permute 64-bit words as follows to shift words left
    // v[0] = v[0] ;
    // v[1] = v[0] ;
    // v[2] = v[1] ;
    // v[3] = v[2] ;
    // Selector bits: 10 01 00 00 = 0x90
    bitV = _mm256_or_si256(_mm256_slli_epi64(bitV, 1),
			   _mm256_permute4x64_epi64(msb2lsb, 0x90)) ;

    // The LSB needs to be cleared
    bitV = _mm256_and_si256(bitV, lsbClearVec) ;
  }
  
  inline BitVec256 shiftProbesRight(int probeOffset, const BitVec256& shiftCnt) const {
    
    BitVec256 r ;
    if (probeOffset == 15) {
      // shift right by 16 bits (across 128-bit boundary)
      r.bitV = _mm256_or_si256(_mm256_srli_si256(bitV, 2),
			       _mm256_permute4x64_epi64(_mm256_slli_si256(bitV,14), 0x0C)) ;
    }
    else {
      // shift right by shiftCnt within each 16-bit word
      r.bitV = _mm256_srl_epi64(bitV, _mm256_castsi256_si128(shiftCnt.bitV)) ;
    }

    return r ;
  }

  void setWord (int index, uint64_t w) {
    switch(index) {
    case 0: bitV = _mm256_insertf128_si256(bitV, 
					   _mm_insert_epi64(_mm256_castsi256_si128(bitV), w, 0),
					   0) ; 
      break ;

    case 1: bitV = _mm256_insertf128_si256(bitV, 
					   _mm_insert_epi64(_mm256_castsi256_si128(bitV), w, 1),
					   0) ; 
      break ;
      
    case 2: bitV = _mm256_insertf128_si256(bitV, 
					   _mm_insert_epi64(_mm256_extractf128_si256(bitV, 1), 
							    w, 0),
					   1) ; 
      break ;

    case 3: bitV = _mm256_insertf128_si256(bitV, 
					   _mm_insert_epi64(_mm256_extractf128_si256(bitV, 1), 
							    w, 1),
					   1) ; 
      break ;

    default: assert(false) ;
    }
  }

  uint64_t getWord (int index) const {

    switch(index) {
    case 0: 
      return _mm_extract_epi64(_mm256_castsi256_si128(bitV), 0) ;
    case 1:
      return _mm_extract_epi64(_mm256_castsi256_si128(bitV), 1) ;
    case 2:
      return _mm_extract_epi64(_mm256_extractf128_si256(bitV,1), 0) ;
    case 3:
      return _mm_extract_epi64(_mm256_extractf128_si256(bitV,1), 1) ;
    default:
      assert(false) ;
    }

    return 0 ;
  }


  
  inline BitVec256 operator~ () const {
    static const __m256i onesVec = _mm256_set1_epi8(0xFF) ;

    return BitVec256(_mm256_andnot_si256(bitV, onesVec)) ;
  }

  inline BitVec256 operator& (const BitVec256& other) const {

    return BitVec256(_mm256_and_si256(bitV, other.bitV)) ;
  }

  inline BitVec256 andnot (const BitVec256& other) const {

    return BitVec256(_mm256_andnot_si256(other.bitV, bitV)) ;
  }

  inline BitVec256 operator| (const BitVec256& other) const {

    return BitVec256(_mm256_or_si256(bitV, other.bitV)) ;
  }

  inline BitVec256 operator^ (const BitVec256& other) const {

    return BitVec256(_mm256_xor_si256(bitV, other.bitV)) ;
  }

  inline BitVec256 operator+ (const BitVec256& other) const {

    BitVec256 sum ;
    sum.bitV = _mm256_add_epi64(this->bitV, other.bitV) ;
    
    // We should add the carryout from the lower word (if any) to upper word
    static __m256i zeroVec = _mm256_set1_epi64x(0x0) ;
    static __m256i lsbSetVec1 = _mm256_set_epi64x(0x1, 0x1, 0x1, 0x0) ;
    static __m256i lsbSetVec2 = _mm256_set_epi64x(0x1, 0x1, 0x0, 0x0) ;
    static __m256i lsbSetVec3 = _mm256_set_epi64x(0x1, 0x0, 0x0, 0x0) ;
    static __m256i msbSetVec = _mm256_set1_epi64x(0x8000000000000000) ;

    // set pred to all 1s (0xFFFFFF...) in case of overflow
    __m256i pred = _mm256_cmpgt_epi64(_mm256_xor_si256(this->bitV, msbSetVec),
				      _mm256_xor_si256(sum.bitV, msbSetVec)) ;
    // Note: Since SSE/avx do not have unsigned comparison opn, we first xor the MSBs,
    // and then do the signed comparison

    __m256i carry1 = _mm256_permute4x64_epi64(pred, 0x90) ;

    // If lower 64-bit of pred is set, add 0x1 to the upper word
    sum.bitV = _mm256_add_epi64(sum.bitV, _mm256_and_si256(lsbSetVec1, carry1)) ;

    // It seems rare that the second and third round carries are needed.
    // We disable the code below for better performance.
    return sum ; 

    pred = _mm256_cmpeq_epi64(sum.bitV, zeroVec) ;
    __m256i carry2 = _mm256_permute4x64_epi64(_mm256_and_si256(pred, carry1), 0x90) ;
    sum.bitV = _mm256_add_epi64(sum.bitV, _mm256_and_si256(lsbSetVec2, carry2)) ;

    pred = _mm256_cmpeq_epi64(sum.bitV, zeroVec) ;
    __m256i carry3 = _mm256_permute4x64_epi64(_mm256_and_si256(pred, carry2), 0x90) ;
    sum.bitV = _mm256_add_epi64(sum.bitV, _mm256_and_si256(lsbSetVec3, carry3)) ;
					 
    
    return sum ; // zzz
    
  }

  static bool isQueryLengthOk(int qlen) {
    return qlen <= 255 ;
  }

  friend std::ostream& operator<< (std::ostream& os, const BitVec256& b) {
    os << "{" <<  b.getWord(0) << " " << b.getWord(1) << " " << b.getWord(2) << " "
       << b.getWord(3) << "}" ;

    return os ;
  }

  friend class EDVec128x2Every16 ;
  friend class EDVec256Every16 ;
  
} ;



#endif // #ifndef _FAST_EXTEND_BITV128_H
