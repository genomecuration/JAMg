static char rcsid[] = "$Id: sedgesort.c 216943 2018-10-10 07:13:39Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "sedgesort.h"
#include "mem.h"
#if 0
#include <stdio.h>
#include <stdlib.h>
#endif
#include <limits.h>		/* For INT_MAX */


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


#ifdef LARGE_GENOMES
static void
insertion_sort_ptr_uint8 (register UINT8 array[], register int len) {
  register UINT8 *ptri, *ptrj, *ptrk;
  register UINT8 temp;
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
#endif



static void
insertion_sort_order (register int order[], register unsigned int array[], register int len) {
  register int i, j;
  register int temp;
  /* register unsigned int *hi, *lo; */
  register unsigned int temporary;
  int thresh;

  thresh = (len < CUTOFF + 1) ? len : CUTOFF + 1;

  /* Find smallest element in first threshold and place it at array
     beginning.  This is the smallest array element, and the operation
     speeds up the inner loop of insertion sort */

  j = 0;
  for (i = 1; i < thresh; i++) {
    if (array[order[i]] < array[order[j]]) {
      j = i;
    }
  }
  SWAP(order[j],order[0]);


#if 0
  /* Insertion sort, running from left to right */
  i = 1;
  while (++i < len) {
    temporary = array[order[i]];
    j = i - 1;
    while (temporary < array[order[j]]) {
      j--;
    }
    j++;
    
    if (j != i) {
      hi = lo = &(array[order[i]]);
      while (--lo >= &(array[order[j]])) {
	*hi = *lo;
	hi = lo;
      }
      *hi = temporary;
    }
  }
#else
  /* Faster under -O3 */
  for (i = 1; i < len; i++) {
    j = i;
    temp = order[j];
    temporary = array[temp];
    while (array[order[j-1]] > temporary) {
      order[j] = order[j-1];
      j--;
    }
    order[j] = temp;
  }
#endif


  return;
}

static void
insertion_sort_order_int (register int order[], register int array[], register int len) {
  register int i, j;
  register int temp;
  /* register int *hi, *lo; */
  register int temporary;
  int thresh;

  thresh = (len < CUTOFF + 1) ? len : CUTOFF + 1;

  /* Find smallest element in first threshold and place it at array
     beginning.  This is the smallest array element, and the operation
     speeds up the inner loop of insertion sort */

  j = 0;
  for (i = 1; i < thresh; i++) {
    if (array[order[i]] < array[order[j]]) {
      j = i;
    }
  }
  SWAP(order[j],order[0]);


#if 0
  /* Insertion sort, running from left to right */
  i = 1;
  while (++i < len) {
    temporary = array[order[i]];
    j = i - 1;
    while (temporary < array[order[j]]) {
      j--;
    }
    j++;
    
    if (j != i) {
      hi = lo = &(array[order[i]]);
      while (--lo >= &(array[order[j]])) {
	*hi = *lo;
	hi = lo;
      }
      *hi = temporary;
    }
  }
#else
  /* Faster under -O3 */
  for (i = 1; i < len; i++) {
    j = i;
    temp = order[j];
    temporary = array[temp];
    while (array[order[j-1]] > temporary) {
      order[j] = order[j-1];
      j--;
    }
    order[j] = temp;
  }
#endif


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
partial_quickersort_order (register int order[], register unsigned int array[],
			   register int lower, register int upper) {
  register int i, j;
  register int temp;
  register unsigned int pivot;

  if (upper - lower > CUTOFF) {
    SWAP(order[lower], order[(upper+lower)/2]);
    i = lower;
    j = upper + 1;
    pivot = array[order[lower]];

    while (1) {
      do i++; while (array[order[i]] < pivot);
      do j--; while (array[order[j]] > pivot);
      if (j < i) break;
      SWAP(order[i], order[j]);
    }

    SWAP(order[lower], order[j]);
    partial_quickersort_order(order, array, lower, j-1);
    partial_quickersort_order(order, array, i, upper);
  }

  return;
}


static void
partial_quickersort_order_int (register int order[], register int array[],
			   register int lower, register int upper) {
  register int i, j;
  register int temp;
  register int pivot;

  if (upper - lower > CUTOFF) {
    SWAP(order[lower], order[(upper+lower)/2]);
    i = lower;
    j = upper + 1;
    pivot = array[order[lower]];

    while (1) {
      do i++; while (array[order[i]] < pivot);
      do j--; while (array[order[j]] > pivot);
      if (j < i) break;
      SWAP(order[i], order[j]);
    }

    SWAP(order[lower], order[j]);
    partial_quickersort_order_int(order, array, lower, j-1);
    partial_quickersort_order_int(order, array, i, upper);
  }

  return;
}



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

#ifdef LARGE_GENOMES
static void
partial_quickersort_ptr_uint8 (register UINT8 array[], register int lower,
			       register int upper) {
  register UINT8 *ptri, *ptrj;
  register UINT8 temp, pivot;

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
    partial_quickersort_ptr_uint8(array, lower, (ptrj-1) - array);
    partial_quickersort_ptr_uint8(array, ptri - array, upper);
  }

  return;
}
#endif


#if 0
static void
partial_quickersort_ptr_order (register int order[], register unsigned int array[],
			       register int lower, register int upper) {
  register int *ptri, *ptrj, temp;
  register unsigned int pivot;

  if (upper - lower > CUTOFF) {
    SWAP(order[lower], order[(upper+lower)/2]);
    ptri = &(order[lower]);
    ptrj = &(order[upper + 1]);
    pivot = array[*ptri];

    while (1) {
      do ptri++; while (array[*ptri] < pivot);
      do ptrj--; while (array[*ptrj] > pivot);
      if (ptrj < ptri) break;
      SWAP(*ptri, *ptrj);
    }

    SWAP(order[lower], *ptrj);
    partial_quickersort_ptr_order(order, array, lower, (ptrj-1) - order);
    partial_quickersort_ptr_order(order, array, ptri - order, upper);
  }

  return;
}
#endif


/* Requires array[len] to be available */
void
Sedgesort_uint4 (register unsigned int array[], register int len) {

  array[len] = (unsigned int) -1;

  partial_quickersort_ptr_uint4(array, 0, len-1);
  insertion_sort_ptr_uint4(array, len);
  return;
}

#ifdef LARGE_GENOMES
/* Requires array[len] to be available */
void
Sedgesort_uint8 (register UINT8 array[], register int len) {

  array[len] = (UINT8) -1;

  partial_quickersort_ptr_uint8(array, 0, len-1);
  insertion_sort_ptr_uint8(array, len);
  return;
}
#endif


int *
Sedgesort_order (register unsigned int array[], register int len) {
  int *order;
  int i;

  order = (int *) MALLOC((len + 1)*sizeof(int));
  for (i = 0; i <= len; i++) {
    order[i] = i;
  }
  array[len] = (unsigned int) -1;

  partial_quickersort_order(order, array, 0, len-1);
  insertion_sort_order(order, array, len);
  return order;
}

int *
Sedgesort_order_int (register int array[], register int len) {
  int *order;
  int i;

  order = (int *) MALLOC((len + 1)*sizeof(int));
  for (i = 0; i <= len; i++) {
    order[i] = i;
  }
  array[len] = INT_MAX;

  partial_quickersort_order_int(order, array, 0, len-1);
  insertion_sort_order_int(order, array, len);
  return order;
}


#if 0
int
main (int argc, char *argv[]) {
  int *order;
  int array[10+1];
  int i;

  array[0] = 3;
  array[1] = 1;
  array[2] = 4;
  array[3] = 1;
  array[4] = 5;
  array[5] = 9;
  array[6] = 2;
  array[7] = 6;
  array[8] = 5;
  array[9] = 3;

  order = Sedgesort_order(array,10);

  for (i = 0; i < 10; i++) {
    printf("%d ",order[i]);
  }
  printf("\n");

  return 0;
}
#endif


  

