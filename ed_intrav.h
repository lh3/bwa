/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef ED_INTRAV_H
#define ED_INTRAV_H

#include "ed_intravED.h"

const int costMatrixRowCnt = 5 ;
const int8_t costMatrix[] = {
  1, -4, -4, -4, -1,
  -4, 1, -4, -4, -1,
  -4, -4, 1, -4, -1,
  -4, -4, -4, 1, -1,
  -1, -1, -1, -1, -1 
} ;

const int gapo = 6 ;
const int gape = 1 ; 

int run_ksw_extend(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
		   int bandW, int initScore, int endBonus, int zdrop, int costMatrixRowCnt,
		   const int8_t* costMatrix, int gapo, int gape,
		   int& alignedQLen, int& alignedRLen) ;

struct SWFeedback {
  int maxQLen ;
  int maxRLen ;
  int maxBand ;

  SWFeedback():maxQLen(0), maxRLen(0), maxBand(0) {} ;

  SWFeedback(int qlen, int rlen, int band): maxQLen(qlen), maxRLen(rlen), maxBand(band) {}
  void updateMax(int qlen, int rlen, int band) {
    maxQLen = std::max(maxQLen, qlen) ;
    maxRLen = std::max(maxRLen, rlen) ;
    maxBand = std::max(maxBand, band) ;
  }
} ;

inline ostream& operator<<(ostream& os, const SWFeedback& swfb) {
  os << "{" << swfb.maxQLen << " " << swfb.maxRLen << " " << swfb.maxBand << "}" ;
  return os ;
}

template<class BitVec, class EDVec>
class DistVec {

  EDVec dist_, mask_ ;
  int queryLen_, msWordIndex_, probeOffset_ ;
  
 public: 
  
 DistVec(int queryLen):
  queryLen_(queryLen), msWordIndex_(EDVec::getLastWordIndexFor(queryLen)),
    probeOffset_(EDVec::getProbeOffsetFor(queryLen)) {

    mask_.setWordsAsMask() ;
  }

  const EDVec& getVec() const {
    return dist_ ;
  }

  void initDist(int initU0, int initU1) {
    dist_.setWordsAsDist(queryLen_, initU0 - initU1) ;
  }

  void addDist (const BitVec& LP0, const BitVec& LP1) {

    /*
    cout << std::hex << "+++++" << LP0 << " " << LP1 << " " << std::dec << probeOffset_ << endl
	 << std::hex << (EDVec(LP0.shiftProbesRight(probeOffset_))) << " "
	 << ((EDVec(LP0.shiftProbesRight(probeOffset_))) & mask_) << endl 
	 << std::hex << (EDVec(LP1.shiftProbesRight(probeOffset_))) << " "
	 << ((EDVec(LP1.shiftProbesRight(probeOffset_))) & mask_)
	 << std::dec << endl ;
    cout << "Dist before: " << dist_ << endl ;
    */

    dist_ += (EDVec(LP0.shiftProbesRight(probeOffset_)) & mask_) ;
    dist_ -= (EDVec(LP1.shiftProbesRight(probeOffset_)) & mask_) ;

    //cout << "Dist after: " << dist_ << endl ;

  }

  EDVec getDeltaDistVec(const BitVec& LP0, const BitVec& LP1) {
    EDVec r = EDVec(LP0.shiftProbesRight(probeOffset_)) & mask_ ;
    r -= EDVec(LP1.shiftProbesRight(probeOffset_)) & mask_ ;

    return r ;
  }

  void addDist (const EDVec& deltaDist) {
    dist_ += deltaDist ;
  }

  int getValidDistCnt() const {
    return msWordIndex_ + 1 ;
  }

  int getProbeColumn(int index) const {
    return EDVec::getDistAt(index, queryLen_) - 1 ;
  }

  int getDist(int index) const {
    return (int) dist_.getWord(index) ;
  }

  int getTotalDist() const {
    return (int) dist_.getWord(msWordIndex_) ;
  }

  void setMin (const DistVec& other) {
    this->dist_.setMin(other.dist_) ;
  }

  friend std::ostream& operator<< (std::ostream& os, const DistVec& d) {
    os << d.dist_ ;
    return os ;
  }

} ;

#endif // #ifndef ED_INTRAV_H
