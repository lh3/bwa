/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#include <stdint.h>
#include <iostream>


class BitCount8 {
  static uint8_t lut_[256] ;

 public:

  static void init_lut() ;
  
  static int getCnt(uint8_t bitVec) {
    return lut_[bitVec] ;
  }

} ;


class BitCount16 {
  static uint8_t lut_[256] ;

 public:

  static void init_lut() ;
  
  static int getCnt(uint16_t bitVec) {
    return BitCount8::getCnt((uint8_t) (bitVec >> 8)) + 
      BitCount8::getCnt((uint8_t) bitVec) ;
  }
  
} ;



// Fine-grain alignment for up to the first 16 bases.
// Considers only alignments where qlen=rlen.
class FineAlignment16 {

  int queryLen_, endBonus_ ;
  
  int bestAlignedLenLow_, bestAlignedLenHigh_ ;
  int maxScoreLow_, maxScoreHigh_ ;

  int mismatchWtP1Low_, mismatchWtP1High_ ;
  int ambigWtP1_ ;
  

 public:

 FineAlignment16(int queryLen, int endBonus, int mismatchWtLow, int mismatchWtHigh, int ambigWt)
   : queryLen_(queryLen), endBonus_(endBonus), 
    mismatchWtP1Low_(mismatchWtLow+1), mismatchWtP1High_(mismatchWtHigh+1), 
    ambigWtP1_(ambigWt+1),
    bestAlignedLenLow_(0), bestAlignedLenHigh_(0), maxScoreLow_(0), maxScoreHigh_(0) {}
  
  void update(uint16_t UP0, uint16_t UP1, int partialLen, int ambigCnt) {

    uint16_t mask = ((uint32_t)(1) << partialLen) - 1 ; // clear higher bits beyond partialLen
    UP0 &= mask ;
    UP1 &= mask ;

    int editDist = (partialLen <= 8) ? 
      (BitCount8::getCnt((uint8_t)UP0) - BitCount8::getCnt((uint8_t)UP1)) :
      (BitCount16::getCnt(UP0) - BitCount16::getCnt(UP1)) ;

    int bonus = (partialLen == queryLen_) ? endBonus_ : 0 ;

    int scoreLow = partialLen - editDist * mismatchWtP1High_ - ambigCnt * ambigWtP1_ + bonus ;
    int scoreHigh = partialLen - editDist * mismatchWtP1Low_ - ambigCnt * ambigWtP1_ + bonus ;

    //cout << std::hex << mask << " " << UP0 << " " << UP1 << std::dec << endl ;
    //cout << BitCount16::getCnt(UP0) << " " << BitCount16::getCnt(UP1) << endl ;
    //cout << std::hex << (int) ((uint8_t) (UP0 >> 8)) << " " << (int) ((uint8_t) UP0) << endl ;
    //cout << "###################len = " << partialLen << ": " << ambigCnt << " " << editDist
    //	 << " " << scoreLow << " " << scoreHigh << endl ;

    if (scoreLow > maxScoreLow_) {
      maxScoreLow_ = scoreLow ;
      bestAlignedLenLow_ = partialLen ;
    }

    if (scoreHigh > maxScoreHigh_) {
      maxScoreHigh_ = scoreHigh ;
      bestAlignedLenHigh_ = partialLen ;
    }

  }

  void getBestLow(int& alignedLen, int& score) {
    alignedLen = bestAlignedLenLow_ ;
    score = maxScoreLow_ ;
  }

  void getBestHigh(int& alignedLen, int& score) {
    alignedLen = bestAlignedLenHigh_ ;
    score = maxScoreHigh_ ;
  }

  friend ostream& operator<< (ostream& os, const FineAlignment16& a) {
    os << "{maxScores = (" << a.maxScoreLow_ << ", " << a.maxScoreHigh_ << "); best cols: "
       << a.bestAlignedLenLow_ << ", " << a.bestAlignedLenHigh_ << "}" ;

    return os ;
  }

} ;
