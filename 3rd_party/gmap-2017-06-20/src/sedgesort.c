static char rcsid[] = "$Id: sedgesort.c 196273 2016-08-12 15:15:06Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "sedgesort.h"


#define CUTOFF 50
#define SWAP(x,y) {temp = (x); (x) = (y); (y) = temp;}

#if 0
/* insertion sort */
static void
insort_uint4 (register unsigned int array[], register int len) {
  register int i, j;
  register unsigned int temp;

  for (i = 1; i < len; i++) {
    j = i;
    temp = array[j];
    while (j > 0 && array[j-1] > temp) {
      array[j] = array[j-1];
      j--;
    }
    array[j] = temp;
  }

  return;
}
#endif


#if 0
static void
insertion_sort_uint4 (register unsigned int array[], register int len) {
  register int i, j;
  register unsigned int *hi, *lo;
  register unsigned int temp;
  int thresh;

  thresh = (len < CUTOFF + 1) ? len : CUTOFF + 1;

  /* Find smallest element in first threshold and place it at array
     beginning.  This is the smallest array element, and the operation
     speeds up the inner loop of insertion sort */

  j = 0;
  for (i = 1; i < thresh; i++) {
    if (array[i] < array[j]) {
      j = i;
    }
  }
  SWAP(array[j],array[0]);


#if 0
  /* Insertion sort, running from left to right */
  i = 1;
  while (++i < len) {
    temp = array[i];
    j = i - 1;
    while (temp < array[j]) {
      j--;
    }
    j++;
    
    if (j != i) {
      hi = lo = &(array[i]);
      while (--lo >= &(array[j])) {
	*hi = *lo;
	hi = lo;
      }
      *hi = temp;
    }
  }
#else
  /* Faster under -O3 */
  for (i = 1; i < len; i++) {
    j = i;
    temp = array[j];
    while (array[j-1] > temp) {
      array[j] = array[j-1];
      j--;
    }
    array[j] = temp;
  }
#endif


  return;
}
#endif


static void
insertion_sort_ptr_uint4 (register unsigned int array[], register int len) {
  register unsigned int *ptri, *ptrj, *ptrk;
  register unsigned int temp;
  int thresh;

  thresh = (len < CUTOFF + 1) ? len : CUTOFF + 1;

  /* Find smallest element in first threshold and place it at array
     beginning.  This is the smallest array element, and the operation
     speeds up the inner loop of insertion sort */

  ptrj = &(array[0]);
  for (ptri = &(array[1]); ptri < &(array[thresh]); ptri++) {
    if (*ptri < *ptrj) {
      ptrj = ptri;
    }
  }
  SWAP(array[0],*ptrj);

  for (ptri = &(array[1]); ptri < &(array[len]); ptri++) {
    temp = *ptri;
    ptrk = ptri;
    ptrj = ptri - 1;
    while (temp < *ptrj) {
      *ptrk = *ptrj;
      ptrk = ptrj--;
    }
    *ptrk = temp;
  }

  return;
}


#if 0
static void
partial_quickersort_uint4 (register unsigned int array[], register int lower,
			   register int upper) {
  register int i, j;
  register unsigned int temp, pivot;

  if (upper - lower > CUTOFF) {
    SWAP(array[lower], array[(upper+lower)/2]);
    i = lower;
    j = upper + 1;
    pivot = array[lower];

    while (1) {
      do i++; while (array[i] < pivot);
      do j--; while (array[j] > pivot);
      if (j < i) break;
      SWAP(array[i], array[j]);
    }

    SWAP(array[lower], array[j]);
    partial_quickersort_uint4(array, lower, j-1);
    partial_quickersort_uint4(array, i, upper);
  }

  return;
}
#endif


static void
partial_quickersort_ptr_uint4 (register unsigned int array[], register int lower,
			       register int upper) {
  register unsigned int *ptri, *ptrj;
  register unsigned int temp, pivot;

  if (upper - lower > CUTOFF) {
    SWAP(array[lower], array[(upper+lower)/2]);
    ptri = &(array[lower]);
    ptrj = &(array[upper + 1]);
    pivot = *ptri;

    while (1) {
      do ptri++; while (*ptri < pivot);
      do ptrj--; while (*ptrj > pivot);
      if (ptrj < ptri) break;
      SWAP(*ptri, *ptrj);
    }

    SWAP(array[lower], *ptrj);
    partial_quickersort_ptr_uint4(array, lower, (ptrj-1) - array);
    partial_quickersort_ptr_uint4(array, ptri - array, upper);
  }

  return;
}


/* Requires array[len] to be available */
void
Sedgesort_uint4 (register unsigned int array[], register int len) {

  array[len] = -1U;

  partial_quickersort_ptr_uint4(array, 0, len-1);
  insertion_sort_ptr_uint4(array, len);
  return;
}


  

