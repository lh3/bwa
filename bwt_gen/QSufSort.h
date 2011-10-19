/* QSufSort.h

   Header file for QSufSort.c

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

   Modified by Wong Chi-Kwong, 2004

   Changes summary:	- Used long variable and function names
					- Removed global variables
					- Replace pointer references with array references
					- Used insertion sort in place of selection sort and increased insertion sort threshold
					- Reconstructing suffix array from inverse becomes an option
					- Add handling where end-of-text symbol is not necessary < all characters
					- Removed codes for supporting alphabet size > number of characters
  
  No warrenty is given regarding the quality of the modifications.

*/

#ifndef __QSUFSORT_H__
#define __QSUFSORT_H__

#include <stdint.h>

#define KEY(V, I, p, h)					( V[ I[p] + h ] )
#define INSERT_SORT_NUM_ITEM	16

typedef int64_t qsint_t;
#define QSINT_MAX INT64_MAX

void QSufSortSuffixSort(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t numChar, const qsint_t largestInputSymbol, 
						const qsint_t smallestInputSymbol, const int skipTransform);
void QSufSortGenerateSaFromInverse(const qsint_t *V, qsint_t* __restrict I, const qsint_t numChar);


#endif
