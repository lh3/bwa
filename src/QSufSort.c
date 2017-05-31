/* QSufSort.c

   Original source from qsufsort.c

   Copyright 1999, N. Jesper Larsson, all rights reserved.

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


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "QSufSort.h"

#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define swap(a, b, t);							t = a; a = b; b = t;

// Static functions
static void QSufSortSortSplit(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t lowestPos, 
							  const qsint_t highestPos, const qsint_t numSortedChar);
static qsint_t QSufSortChoosePivot(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t lowestPos, 
							   const qsint_t highestPos, const qsint_t numSortedChar);
static void QSufSortInsertSortSplit(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t lowestPos, 
									const qsint_t highestPos, const qsint_t numSortedChar);
static void QSufSortBucketSort(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t numChar, const qsint_t alphabetSize);
static qsint_t QSufSortTransform(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t numChar, const qsint_t largestInputSymbol, 
							 const qsint_t smallestInputSymbol, const qsint_t maxNewAlphabetSize, qsint_t *numSymbolAggregated);

/* Makes suffix array p of x. x becomes inverse of p. p and x are both of size
   n+1. Contents of x[0...n-1] are integers in the range l...k-1. Original
   contents of x[n] is disregarded, the n-th symbol being regarded as
   end-of-string smaller than all other symbols.*/
void QSufSortSuffixSort(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t numChar, const qsint_t largestInputSymbol, 
						const qsint_t smallestInputSymbol, const int skipTransform)
{
	qsint_t i, j;
	qsint_t s, negatedSortedGroupLength;
	qsint_t numSymbolAggregated;
	qsint_t numSortedPos = 1;
	qsint_t newAlphabetSize;
   
	if (!skipTransform) {
		/* bucketing possible*/
		newAlphabetSize = QSufSortTransform(V, I, numChar, largestInputSymbol, smallestInputSymbol, 
											numChar, &numSymbolAggregated);
		QSufSortBucketSort(V, I, numChar, newAlphabetSize);
		I[0] = -1;
		V[numChar] = 0;
		numSortedPos = numSymbolAggregated;
	}

	while ((qsint_t)(I[0]) >= -(qsint_t)numChar) {
		i = 0;
		negatedSortedGroupLength = 0;
		do {
			s = I[i];
			if (s < 0) {
				i -= s;						/* skip over sorted group.*/
				negatedSortedGroupLength += s;
			} else {
				if (negatedSortedGroupLength) {
					I[i+negatedSortedGroupLength] = negatedSortedGroupLength;	/* combine preceding sorted groups */
					negatedSortedGroupLength = 0;
				}
				j = V[s] + 1;
				QSufSortSortSplit(V, I, i, j - 1, numSortedPos);
				i = j;
			}
		} while (i <= numChar);
		if (negatedSortedGroupLength) {
			/* array ends with a sorted group.*/
			I[i+negatedSortedGroupLength] = negatedSortedGroupLength;	/* combine sorted groups at end of I.*/
		}
		numSortedPos *= 2;	/* double sorted-depth.*/
	}
}

void QSufSortGenerateSaFromInverse(const qsint_t* V, qsint_t* __restrict I, const qsint_t numChar)
{
	qsint_t i;
	for (i=0; i<=numChar; i++)
		I[V[i]] = i + 1;
}

/* Sorting routine called for each unsorted group. Sorts the array of integers
   (suffix numbers) of length n starting at p. The algorithm is a ternary-split
   quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
   Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
   function is based on Program 7.*/
static void QSufSortSortSplit(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t lowestPos, 
							  const qsint_t highestPos, const qsint_t numSortedChar) {

	qsint_t a, b, c, d;
	qsint_t l, m;
	qsint_t f, v, s, t;
	qsint_t tmp;
	qsint_t numItem;

	numItem = highestPos - lowestPos + 1;

	if (numItem <= INSERT_SORT_NUM_ITEM) {
		QSufSortInsertSortSplit(V, I, lowestPos, highestPos, numSortedChar);
		return;
	}

	v = QSufSortChoosePivot(V, I, lowestPos, highestPos, numSortedChar);

	a = b = lowestPos;
	c = d = highestPos;

	while (1) {
		while (c >= b && (f = KEY(V, I, b, numSortedChar)) <= v) {
			if (f == v) {
				swap(I[a], I[b], tmp);
				a++;
			}
			b++;
		}
		while (c >= b && (f = KEY(V, I, c, numSortedChar)) >= v) {
			if (f == v) {
				swap(I[c], I[d], tmp);
				d--;
			}
			c--;
		}
		if (b > c)
			break;
		swap(I[b], I[c], tmp);
		b++;
		c--;
	}

	s = a - lowestPos;
	t = b - a;
	s = min(s, t);
	for (l = lowestPos, m = b - s; m < b; l++, m++) {
		swap(I[l], I[m], tmp);
	}

	s = d - c;
	t = highestPos - d;
	s = min(s, t);
	for (l = b, m = highestPos - s + 1; m <= highestPos; l++, m++) {
		swap(I[l], I[m], tmp);
	}

	s = b - a;
	t = d - c;
	if (s > 0)
		QSufSortSortSplit(V, I, lowestPos, lowestPos + s - 1, numSortedChar);

	// Update group number for equal portion
	a = lowestPos + s;
	b = highestPos - t;
	if (a == b) {
		// Sorted group
		V[I[a]] = a;
		I[a] = -1;
	} else {
		// Unsorted group
		for (c=a; c<=b; c++)
			V[I[c]] = b;
	}

	if (t > 0)
		QSufSortSortSplit(V, I, highestPos - t + 1, highestPos, numSortedChar);

}

/* Algorithm by Bentley & McIlroy.*/
static qsint_t QSufSortChoosePivot(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t lowestPos, 
							   const qsint_t highestPos, const qsint_t numSortedChar) {

	qsint_t m;
	qsint_t keyl, keym, keyn;
	qsint_t key1, key2, key3;
	qsint_t s;
	qsint_t numItem;

	numItem = highestPos - lowestPos + 1;

	m = lowestPos + numItem / 2;

	s = numItem / 8;
	key1 = KEY(V, I, lowestPos, numSortedChar);
	key2 = KEY(V, I, lowestPos+s, numSortedChar);
	key3 = KEY(V, I, lowestPos+2*s, numSortedChar);
	keyl = med3(key1, key2, key3);
	key1 = KEY(V, I, m-s, numSortedChar);
	key2 = KEY(V, I, m, numSortedChar);
	key3 = KEY(V, I, m+s, numSortedChar);
	keym = med3(key1, key2, key3);
	key1 = KEY(V, I, highestPos-2*s, numSortedChar);
	key2 = KEY(V, I, highestPos-s, numSortedChar);
	key3 = KEY(V, I, highestPos, numSortedChar);
	keyn = med3(key1, key2, key3);

	return med3(keyl, keym, keyn);


}

/* Quadratic sorting method to use for small subarrays. */
static void QSufSortInsertSortSplit(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t lowestPos, 
									const qsint_t highestPos, const qsint_t numSortedChar)
{
	qsint_t i, j;
	qsint_t tmpKey, tmpPos;
	qsint_t numItem;
	qsint_t key[INSERT_SORT_NUM_ITEM], pos[INSERT_SORT_NUM_ITEM];
	qsint_t negativeSortedLength;
	qsint_t groupNum;

	numItem = highestPos - lowestPos + 1;

	for (i=0; i<numItem; i++) {
		pos[i] = I[lowestPos + i];
		key[i] = V[pos[i] + numSortedChar];
	}

	for (i=1; i<numItem; i++) {
		tmpKey = key[i];
		tmpPos = pos[i];
		for (j=i; j>0 && key[j-1] > tmpKey; j--) {
			key[j] = key[j-1];
			pos[j] = pos[j-1];
		}
		key[j] = tmpKey;
		pos[j] = tmpPos;
	}

	negativeSortedLength = -1;

	i = numItem - 1;
	groupNum = highestPos;
	while (i > 0) {
		I[i+lowestPos] = pos[i];
		V[I[i+lowestPos]] = groupNum;
		if (key[i-1] == key[i]) {
			negativeSortedLength = 0;
		} else {
			if (negativeSortedLength < 0)
				I[i+lowestPos] = negativeSortedLength;
			groupNum = i + lowestPos - 1;
			negativeSortedLength--;
		}
		i--;
	}

	I[lowestPos] = pos[0];
	V[I[lowestPos]] = groupNum;
	if (negativeSortedLength < 0)
		I[lowestPos] = negativeSortedLength;
}

/* Bucketsort for first iteration.

   Input: x[0...n-1] holds integers in the range 1...k-1, all of which appear
   at least once. x[n] is 0. (This is the corresponding output of transform.) k
   must be at most n+1. p is array of size n+1 whose contents are disregarded.

   Output: x is V and p is I after the initial sorting stage of the refined
   suffix sorting algorithm.*/
      
static void QSufSortBucketSort(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t numChar, const qsint_t alphabetSize)
{
	qsint_t i, c;
	qsint_t d;
	qsint_t groupNum;
	qsint_t currentIndex;

	// mark linked list empty
	for (i=0; i<alphabetSize; i++)
		I[i] = -1;

	// insert to linked list
	for (i=0; i<=numChar; i++) {
		c = V[i];
		V[i] = (qsint_t)(I[c]);
		I[c] = i;
	}

	currentIndex = numChar;
	for (i=alphabetSize; i>0; i--) {
		c = I[i-1];
		d = (qsint_t)(V[c]);
		groupNum = currentIndex;
		V[c] = groupNum;
		if (d >= 0) {
			I[currentIndex] = c;
			while (d >= 0) {
				c = d;
				d = V[c];
				V[c] = groupNum;
				currentIndex--;
				I[currentIndex] = c;
			}
		} else {
			// sorted group
			I[currentIndex] = -1;
		}
		currentIndex--;
	}
}

/* Transforms the alphabet of x by attempting to aggregate several symbols into
   one, while preserving the suffix order of x. The alphabet may also be
   compacted, so that x on output comprises all integers of the new alphabet
   with no skipped numbers.

   Input: x is an array of size n+1 whose first n elements are positive
   integers in the range l...k-1. p is array of size n+1, used for temporary
   storage. q controls aggregation and compaction by defining the maximum intue
   for any symbol during transformation: q must be at least k-l; if q<=n,
   compaction is guaranteed; if k-l>n, compaction is never done; if q is
   INT_MAX, the maximum number of symbols are aggregated into one.
   
   Output: Returns an integer j in the range 1...q representing the size of the
   new alphabet. If j<=n+1, the alphabet is compacted. The global variable r is
   set to the number of old symbols grouped into one. Only x[n] is 0.*/
static qsint_t QSufSortTransform(qsint_t* __restrict V, qsint_t* __restrict I, const qsint_t numChar, const qsint_t largestInputSymbol, 
							 const qsint_t smallestInputSymbol, const qsint_t maxNewAlphabetSize, qsint_t *numSymbolAggregated)
{
	qsint_t c, i, j;
	qsint_t a;	// numSymbolAggregated
	qsint_t mask;
	qsint_t minSymbolInChunk = 0, maxSymbolInChunk = 0;
	qsint_t newAlphabetSize;
	qsint_t maxNumInputSymbol, maxNumBit, maxSymbol;

	maxNumInputSymbol = largestInputSymbol - smallestInputSymbol + 1;

	for (maxNumBit = 0, i = maxNumInputSymbol; i; i >>= 1) ++maxNumBit;
	maxSymbol = QSINT_MAX >> maxNumBit;

	c = maxNumInputSymbol;
	for (a = 0; a < numChar && maxSymbolInChunk <= maxSymbol && c <= maxNewAlphabetSize; a++) {
		minSymbolInChunk = (minSymbolInChunk << maxNumBit) | (V[a] - smallestInputSymbol + 1);
		maxSymbolInChunk = c;
		c = (maxSymbolInChunk << maxNumBit) | maxNumInputSymbol;
	}

	mask = (1 << (a-1) * maxNumBit) - 1;	/* mask masks off top old symbol from chunk.*/
	V[numChar] = smallestInputSymbol - 1;	/* emulate zero terminator.*/

	/* bucketing possible, compact alphabet.*/
	for (i=0; i<=maxSymbolInChunk; i++)
		I[i] = 0;	/* zero transformation table.*/
	c = minSymbolInChunk;
	for (i=a; i<=numChar; i++) {
		I[c] = 1;			/* mark used chunk symbol.*/
		c = ((c & mask) << maxNumBit) | (V[i] - smallestInputSymbol + 1);	/* shift in next old symbol in chunk.*/
	}
	for (i=1; i<a; i++) {	/* handle last r-1 positions.*/
		I[c] = 1;			/* mark used chunk symbol.*/
		c = (c & mask) << maxNumBit;	/* shift in next old symbol in chunk.*/
	}
	newAlphabetSize = 1;
	for (i=0; i<=maxSymbolInChunk; i++) {
		if (I[i]) {
			I[i] = newAlphabetSize;
			newAlphabetSize++;
		}
	}
	c = minSymbolInChunk;
	for (i=0, j=a; j<=numChar; i++, j++) {
		V[i] = I[c];						/* transform to new alphabet.*/
		c = ((c & mask) << maxNumBit) | (V[j] - smallestInputSymbol + 1);	/* shift in next old symbol in chunk.*/
	}
	for (; i<numChar; i++) {	/* handle last a-1 positions.*/
		V[i] = I[c];			/* transform to new alphabet.*/
		c = (c & mask) << maxNumBit;	/* shift right-end zero in chunk.*/
	}

	V[numChar] = 0;		/* end-of-string symbol is zero.*/

    *numSymbolAggregated = a;
	return newAlphabetSize;
}
