static char rcsid[] = "$Id: merge-diagonals-heap.c 216918 2018-10-09 02:34:49Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "merge-diagonals-heap.h"
#include "assert.h"
#include "mem.h"
#include "list.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif


#define PARENT(i) (i >> 1)
#define LEFT(i) (i << 1)
#define RIGHT(i) ((i << 1) | 1)

#define GETPOS(high,low) (((Univcoord_T) high << 32) + low)


#if !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
static void
min_heap_insert (Univcoord_T **heap, int heapi, Univcoord_T *diagonals) {
  int i;
  Univcoord_T diagonal;

  i = heapi;
  diagonal = diagonals[0];
  while (i > 1 && (*(heap[PARENT(i)]) > diagonal)) {
    heap[i] = heap[PARENT(i)];
    i = PARENT(i);
  }
  heap[i] = diagonals;

  return;
}
#endif


#if !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
/* Provide ancestori as inserti */
static void
heapify (Univcoord_T **heap, Univcoord_T *diagonals
#ifdef DEBUG6
	 , int heapsize
#endif
	 ) {
  Univcoord_T diagonal;
  int inserti, smallesti, righti;
#ifdef DEBUG6
  int i;
  Univcoord_T *ptr;
#endif

  diagonal = *diagonals;

  debug6(printf("Starting heapify with %llu\n",(unsigned long long) diagonal));
#ifdef DEBUG6
  for (i = 1; i <= heapsize; i++) {
    printf("%d: ",i);
    ptr = heap[i];
    while (*ptr < (Univcoord_T) -1) {
      printf(" %u",*ptr);
      ptr++;
    }
    printf("\n");
  }
  printf("\n");
#endif

  inserti = 1;
  smallesti = (*(heap[3]) < *(heap[2])) ? 3 : 2;
  debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		2,3,(unsigned long long) *(heap[2]),(unsigned long long) *(heap[3])));
  while (diagonal > *(heap[smallesti])) {
    heap[inserti] = heap[smallesti];
    inserti = smallesti;
    smallesti = LEFT(inserti);
    righti = smallesti+1;
    debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		  smallesti,righti,(unsigned long long) *(heap[smallesti]),
		  (unsigned long long) *(heap[righti])));
    if (*(heap[righti]) < *(heap[smallesti])) {
      smallesti = righti;
    }
  }
  heap[inserti] = diagonals;
  debug6(printf("Inserting at %d\n\n",inserti));
  return;
}
#endif


#if !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
/* No sorting.  Does not use SIMD merge, so non-destructive.  However,
   we need to add a sentinel to the end of each stream. */
static Univcoord_T **
make_diagonals_heap (List_T *free_list, unsigned char **stream_high_array, UINT4 **stream_low_array,
		     int *streamsize_array, int *diagterm_array, int nstreams) {
  unsigned char *stream_high;
  UINT4 *stream_low;
  Univcoord_T **heap, *copy, *storage, diagterm;
  int heapsize, heapi;
  int nelts, l;
  int streami;


  heapsize = 2*nstreams + 1;
  heap = (Univcoord_T **) CALLOC((heapsize + 1),sizeof(Univcoord_T *));

  debug6(printf("nstreams %d, heapsize %d, heap defined from 0..%d\n",
		nstreams,heapsize,heapsize));
  
  *free_list = (List_T) NULL;
  heapi = 1;
  for (streami = 0; streami < nstreams; streami++) {
    stream_high = stream_high_array[streami];
    stream_low = stream_low_array[streami];
    nelts = streamsize_array[streami];
    diagterm = diagterm_array[streami];

    copy = (Univcoord_T *) MALLOC((nelts+1)*sizeof(Univcoord_T));
    *free_list = List_push(*free_list,(void *) copy);

    for (l = 0; l < nelts; l++) {
      assert(GETPOS(stream_high[l],stream_low[l]) >= (Univcoord_T) -diagterm);
      copy[l] = GETPOS(stream_high[l],stream_low[l]) + diagterm;
    }

    copy[nelts] = (Univcoord_T) -1;	/* sentinel */
    
    min_heap_insert(heap,heapi,copy);
    heapi++;
  }
		    
  /* Set up rest of heap */
  storage = (Univcoord_T *) MALLOC(sizeof(Univcoord_T));
  storage[0] = (Univcoord_T) -1; 	/* sentinel */
  *free_list = List_push(*free_list,(void *) storage);
  while (heapi <= heapsize) {
    heap[heapi] = storage;
    heapi++;
  }

  return heap;
}
#endif


#if !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
/* Output: Not aligned, unlike SIMD version */
Univcoord_T *
Merge_diagonals_large (int *nelts, unsigned char **stream_high_array, UINT4 **stream_low_array,
		       int *streamsize_array, int *diagterm_array, int nstreams) {
  unsigned char *stream_high;
  UINT4 *stream_low;
  Univcoord_T *result, *out, **heap, *diagonals, diagonal, diagterm;
  int streami;
  List_T free_list, p;
  int l;

  if (nstreams == 0) {
    *nelts = 0;
    return (Univcoord_T *) NULL;

  } else if (nstreams == 1) {
    *nelts = streamsize_array[0];
    stream_high = stream_high_array[0];
    stream_low = stream_low_array[0];
    diagterm = diagterm_array[0];

    result = MALLOC((*nelts)*sizeof(Univcoord_T));
    for (l = 0; l < *nelts; l++) {
      assert(GETPOS(stream_high[l],stream_low[l]) >= (Univcoord_T) -diagterm);
      result[l] = GETPOS(stream_high[l],stream_low[l]) + diagterm;
    }

    return result;

  } else {
    *nelts = 0;
    for (streami = 0; streami < nstreams; streami++) {
      *nelts += streamsize_array[streami];
    }
    out = result = MALLOC((*nelts)*sizeof(Univcoord_T));

    heap = make_diagonals_heap(&free_list,stream_high_array,stream_low_array,
			       streamsize_array,diagterm_array,nstreams);

    while ((diagonal = *(heap[1])) < (Univcoord_T) -1) {
      *out++ = diagonal;
      diagonals = ++(heap[1]);	/* Advance pointer */
#ifdef DEBUG6
      heapify(heap,diagonals,/*heapsize*/2*nstreams+1);
#else
      heapify(heap,diagonals);
#endif
    }

    for (p = free_list; p != NULL; p = List_next(p)) {
      diagonals = (Univcoord_T *) List_head(p);
      FREE(diagonals);
    }
    List_free(&free_list);
    FREE(heap);

    return result;
  }

}
#endif


#if !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
/* No sorting.  Does not use SIMD merge, so non-destructive.  However,
   we need to add a sentinel to the end of each stream. */
static Univcoord_T **
make_univdiagonals_heap (List_T *free_list, Univcoord_T **stream_array,
			 int *streamsize_array, int nstreams) {
  Univcoord_T *stream;
  Univcoord_T **heap, *copy, *storage;
  int heapsize, heapi;
  int nelts;
  int streami;


  heapsize = 2*nstreams + 1;
  heap = (Univcoord_T **) CALLOC((heapsize + 1),sizeof(Univcoord_T *));

  debug6(printf("nstreams %d, heapsize %d, heap defined from 0..%d\n",
		nstreams,heapsize,heapsize));
  
  *free_list = (List_T) NULL;
  heapi = 1;
  for (streami = 0; streami < nstreams; streami++) {
    stream = stream_array[streami];
    nelts = streamsize_array[streami];

    copy = (Univcoord_T *) MALLOC((nelts+1)*sizeof(Univcoord_T));
    *free_list = List_push(*free_list,(void *) copy);

    memcpy(copy,stream,nelts*sizeof(Univcoord_T));
    copy[nelts] = (Univcoord_T) -1;	/* sentinel */
    
    min_heap_insert(heap,heapi,copy);
    heapi++;
  }
		    
  /* Set up rest of heap */
  storage = (Univcoord_T *) MALLOC(sizeof(Univcoord_T));
  storage[0] = (Univcoord_T) -1; 	/* sentinel */
  *free_list = List_push(*free_list,(void *) storage);
  while (heapi <= heapsize) {
    heap[heapi] = storage;
    heapi++;
  }

  return heap;
}
#endif


/* For AVX2 and AVX512 version, see merge-diagonals-simd-uint8.c */
#if !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
/* Output: Not aligned, unlike SIMD version */
Univcoord_T *
Merge_diagonals_uint8 (int *nelts, Univcoord_T **stream_array, int *streamsize_array,
		       int nstreams) {
  Univcoord_T *stream;
  Univcoord_T *result, *out, **heap, *diagonals, diagonal;
  List_T free_list, p;
  int streami;

  if (nstreams == 0) {
    *nelts = 0;
    return (Univcoord_T *) NULL;

  } else if (nstreams == 1) {
    *nelts = streamsize_array[0];
    stream = stream_array[0];

    result = MALLOC((*nelts)*sizeof(Univcoord_T));
    memcpy(result,stream,(*nelts)*sizeof(Univcoord_T));
    return result;

  } else {
    *nelts = 0;
    for (streami = 0; streami < nstreams; streami++) {
      *nelts += streamsize_array[streami];
    }
    out = result = MALLOC((*nelts)*sizeof(Univcoord_T));

    heap = make_univdiagonals_heap(&free_list,stream_array,streamsize_array,nstreams);

    while ((diagonal = *(heap[1])) < (Univcoord_T) -1) {
      *out++ = diagonal;
      diagonals = ++(heap[1]);	/* Advance pointer */
#ifdef DEBUG6
      heapify(heap,diagonals,/*heapsize*/2*nstreams+1);
#else
      heapify(heap,diagonals);
#endif
    }

    for (p = free_list; p != NULL; p = List_next(p)) {
      diagonals = (Univcoord_T *) List_head(p);
      FREE(diagonals);
    }
    List_free(&free_list);
    FREE(heap);

    return result;
  }
}
#endif

