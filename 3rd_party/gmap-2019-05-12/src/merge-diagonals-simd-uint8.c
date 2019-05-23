static char rcsid[] = "$Id: merge-diagonals-simd-uint8.c 216918 2018-10-09 02:34:49Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "merge-diagonals-simd-uint8.h"
#include "assert.h"
#include "mem.h"
#include "popcount.h"		/* For clz_table */
#include "sedgesort.h"
#include "merge-uint8.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */


#if defined(HAVE_SSE4_1)
#include <smmintrin.h>
#endif
#if defined(HAVE_AVX2)
#include <immintrin.h>
#endif
#if defined(HAVE_AVX512)
#include <immintrin.h>
#endif


/* #define PYRAMID_SIZE 4 */

#define CUTOFF 1000
#define PYRAMID_SIZE 32

#define GETPOS(high,low) (((Univcoord_T) high << 32) + low)


#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif



#define PARENT(i) (i >> 1)
#define LEFT(i) (i << 1)
#define RIGHT(i) ((i << 1) | 1)

#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
static int
pyramid_merge (UINT8 **heap, int nstreams, int heapsize, int *nelts,
	       int pyramid_start, int pyramid_end) {
  int nodei;
#ifdef DEBUG
  int i;
#endif

  while (pyramid_end > pyramid_start) {
    debug(printf("Merging level: %d..%d for heapsize %d\n",pyramid_start,pyramid_end,heapsize));

    if (pyramid_end > heapsize) {
      nodei = heapsize;
    } else {
      nodei = pyramid_end;
    }

    while (nodei >= pyramid_start) {
      debug2(printf("Merging nodes %d (%d elts) and %d (%d elts) => %d\n",
		    nodei-1,nelts[nodei-1],nodei,nelts[nodei],PARENT(nodei)));
      heap[PARENT(nodei)] = Merge_uint8(/*dest*/NULL,heap[nodei-1],heap[nodei],nelts[nodei-1],nelts[nodei]);
      CHECK_ALIGN(heap[PARENT(nodei)]);
      nelts[PARENT(nodei)] = nelts[nodei-1] + nelts[nodei];
      debug2(printf("Created list %p of length %d at node %d\n",
		    heap[PARENT(nodei)],nelts[PARENT(nodei)],PARENT(nodei)));

#ifdef DEBUG
      for (i = 0; i < nelts[PARENT(nodei)]; i++) {
	printf("%u\n",heap[PARENT(nodei)][i]);
      }
#endif

      /* Don't free original lists (when nodei >= nstreams) */
      debug(printf("Freeing nodes %d and %d\n",nodei-1,nodei));
      if (nodei < nstreams) {
	FREE_ALIGN(heap[nodei]);
      }
      if (nodei-1 < nstreams) {
	FREE_ALIGN(heap[nodei-1]);
      }
      nodei -= 2;
    }

    pyramid_end = PARENT(pyramid_end);
    pyramid_start = PARENT(pyramid_start);
  }

  debug(printf("Returning ancestor %d\n\n",pyramid_start));
  return pyramid_start;
}
#endif


#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
static UINT8 **
make_diagonals_heap (int **nelts, int *ncopied,
		     unsigned char **stream_high_array, UINT4 **stream_low_array,
		     int *streamsize_array, int *diagterm_array, int nstreams) {
  UINT8 **heap, **combined, *out, diagterm;
  int heapsize, heapi, ncombined;
  int *order, *totals, total, n, i, j, k, l;

  /* Put smallest streams at the end of the heap, to increase speed */
  /* Note: Extra space needed for Sedgesort, but querylength > (query_lastpos+1) + 1 */
  order = Sedgesort_order_int(streamsize_array,nstreams);

  combined = (UINT8 **) MALLOC(nstreams*sizeof(UINT8 *));
  totals = (int *) MALLOC(nstreams*sizeof(int));
  ncombined = 0;

  /* Combine the smallest streams, up to CUTOFF, and just use Sedgesort */
  i = 0;
  total = 1;			/* To start the loop */
  while (i < nstreams && total > 0) {
    total = 0;
    j = i;
    while (j < nstreams && total + streamsize_array[order[j]] < CUTOFF) {
      total += streamsize_array[order[j++]];
    }

    if (total > 0) {
      debug(printf("Merging from %d to %d with %d total elements\n",i,j-1,total));
      totals[ncombined] = total;
      /* Need an extra value for Sedgesort_uint8 */
      out = combined[ncombined] = (UINT8 *) MALLOC_ALIGN((total+1)*sizeof(UINT8));
      for (k = i; k < j; k++) {
	n = streamsize_array[order[k]];

	/* Combine stream_high and stream_low and add diagterm (Could use SIMD here) */
	diagterm = (UINT8) diagterm_array[order[k]];
	for (l = 0; l < n; l++) {
	  out[l] = GETPOS(stream_high_array[order[k]][l],stream_low_array[order[k]][l]) + diagterm;
	}

	out += n;
      }
      Sedgesort_uint8(combined[ncombined++],total);
      i = j;
    }
  }
  
  *ncopied = (nstreams - i) + ncombined;
  heapsize = 2*(*ncopied) - 1;
  
  heap = (UINT8 **) CALLOC((heapsize + 1),sizeof(UINT8 *));
  *nelts = (int *) CALLOC((heapsize + 1),sizeof(int));
  
  heapi = heapsize;
  /* Handle individual contents: Start with i value from before */
  while (i < nstreams) {
    n = (*nelts)[heapi] = streamsize_array[order[i]];
    
    /* Copy to make the merging process non-destructive */
    out = heap[heapi] = MALLOC_ALIGN(n*sizeof(UINT8));
    CHECK_ALIGN(heap[heapi]);

    /* Combine stream_high and stream_low and add diagterm (Could use SIMD here) */
    diagterm = (UINT8) diagterm_array[order[i]];
    for (l = 0; l < n; l++) {
      out[l] = GETPOS(stream_high_array[order[i]][l],stream_low_array[order[i]][l]) + diagterm;
    }

#ifdef DEBUG
    printf("Assigning node %d with %d elts:",heapi,(*nelts)[heapi]);
    for (k = 0; k < (*nelts)[heapi]; k++) {
      printf(" %u",heap[heapi][k]);
    }
    printf("\n");
#endif
    heapi--;
    i++;
  }

  /* Handle combined contents */
  for (i = 0; i < ncombined; i++) {
    heap[heapi] = combined[i];
    (*nelts)[heapi] = totals[i];
#ifdef DEBUG
    printf("Assigning node %d with %d elts:",heapi,(*nelts)[heapi]);
    for (k = 0; k < (*nelts)[heapi]; k++) {
      printf(" %u",heap[heapi][k]);
    }
    printf("\n");
#endif
    heapi--;
  }

  FREE(totals);
  FREE(combined);
  FREE(order);

  return heap;
}
#endif


/* For non-AVX2, non-AVX512 code, see merge-diagonals-heap.c */
#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
/* Input: streams might be aligned or not.  Caller will have to free appropriately */
/* Output: aligned */
Univcoord_T *
Merge_diagonals_large (int *nelts1, unsigned char **stream_high_array, UINT4 **stream_low_array,
		       int *streamsize_array, int *diagterm_array, int nstreams) {
  Univcoord_T *result, **heap, diagterm;
  unsigned char *stream_high;
  UINT4 *stream_low;
  int *nelts, l;
  int ncopied, heapi, heapsize, base, ancestori, pyramid_start, pyramid_end;
  int bits;
#ifdef DEBUG
  int i;
#endif


  if (nstreams == 0) {
    *nelts1 = 0;
    return (Univcoord_T *) NULL;

  } else if (nstreams == 1) {
    *nelts1 = streamsize_array[0];
    stream_high = stream_high_array[0];
    stream_low = stream_low_array[0];
    diagterm = diagterm_array[0];

    result = MALLOC_ALIGN((*nelts1)*sizeof(UINT8)); /* Output must be aligned */
    /* Combine stream_high and stream_low and add diagterm (Could use SIMD here) */
    for (l = 0; l < *nelts1; l++) {
      assert(GETPOS(stream_high[l],stream_low[l]) >= (Univcoord_T) -diagterm);
      result[l] = GETPOS(stream_high[l],stream_low[l]) + diagterm;
    }

    return result;

  } else {
    heap = make_diagonals_heap(&nelts,&ncopied,stream_high_array,stream_low_array,
			       streamsize_array,diagterm_array,nstreams);
    if (ncopied == 1) {
      *nelts1 = nelts[1];
      result = heap[1];
    } else {
      heapsize = 2*ncopied - 1;	/* also index of last node */
#ifdef HAVE_BUILTIN_CLZ
      bits = 31 - __builtin_clz((unsigned int) heapsize);
#elif defined(HAVE_ASM_BSR)
      asm("bsr %1,%0" : "=r"(bits) : "r"(heapsize));
#else
      bits = 31 - ((heapsize >> 16) ? clz_table[heapsize >> 16] : 16 + clz_table[heapsize]); 
#endif

      base = (1 << bits);
      debug(printf("nstreams %d, ncopied %d, heapsize %d, clz %d, bits %d, base %d\n",
		   nstreams,ncopied,heapsize,__builtin_clz(heapsize),bits,base));
      
      /* Middle pyramids */
      while (base > PYRAMID_SIZE) {
	for (pyramid_start = 2*base - PYRAMID_SIZE, pyramid_end = 2*base - 1; pyramid_start >= base;
	     pyramid_start -= PYRAMID_SIZE, pyramid_end -= PYRAMID_SIZE) {
	  debug(printf("diagonals: pyramid_start %d, pyramid_end %d, ncopied %d\n",pyramid_start,pyramid_end,ncopied));
	  ancestori = pyramid_merge(heap,ncopied,heapsize,nelts,pyramid_start,pyramid_end);
	}
	base = ancestori;
      }

      /* Last pyramid */
      pyramid_start = base;
      pyramid_end = 2*base - 1;
      debug(printf("diagonals: pyramid_start %d, pyramid_end %d, ncopied %d\n",pyramid_start,pyramid_end,ncopied));
      /* base = */ pyramid_merge(heap,ncopied,heapsize,nelts,pyramid_start,pyramid_end);

      *nelts1 = nelts[1];
      result = heap[1];

      for (heapi = heapsize; heapi > heapsize - ncopied; heapi--) {
	FREE_ALIGN(heap[heapi]);
      }
    }

    FREE(heap);
    FREE(nelts);
  }


#ifdef DEBUG
  printf("Merge_diagonals_large returning result of length %d\n",*nelts1);
  for (i = 0; i < *nelts1; i++) {
    printf("%u\n",result[i]);
  }
#endif

  return result;
}
#endif


#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
static Univcoord_T **
make_univdiagonals_heap (int **nelts, int *ncopied, Univcoord_T **stream_array,
			 int *streamsize_array, int nstreams) {
  Univcoord_T **heap, **combined, *out;
  int heapsize, heapi, ncombined;
  int *order, *totals, total, n, i, j, k;

  /* Put smallest streams at the end of the heap, to increase speed */
  /* Note: Extra space needed for Sedgesort, but querylength > (query_lastpos+1) + 1 */
  order = Sedgesort_order_int(streamsize_array,nstreams);

  combined = (Univcoord_T **) MALLOC(nstreams*sizeof(Univcoord_T *));
  totals = (int *) MALLOC(nstreams*sizeof(int));
  ncombined = 0;

  /* Combine the smallest streams, up to CUTOFF, and just use Sedgesort */
  i = 0;
  total = 1;			/* To start the loop */
  while (i < nstreams && total > 0) {
    total = 0;
    j = i;
    while (j < nstreams && total + streamsize_array[order[j]] < CUTOFF) {
      total += streamsize_array[order[j++]];
    }

    if (total > 0) {
      debug(printf("Merging from %d to %d with %d total elements\n",i,j-1,total));
      totals[ncombined] = total;
      /* Need an extra value for Sedgesort_uint8 */
      out = combined[ncombined] = (Univcoord_T *) MALLOC_ALIGN((total+1)*sizeof(Univcoord_T));
      for (k = i; k < j; k++) {
	n = streamsize_array[order[k]];
	memcpy(out,stream_array[order[k]],n*sizeof(Univcoord_T));
	out += n;
      }
      Sedgesort_uint8(combined[ncombined++],total);
      i = j;
    }
  }
  
  *ncopied = (nstreams - i) + ncombined;
  heapsize = 2*(*ncopied) - 1;
  
  heap = (Univcoord_T **) CALLOC((heapsize + 1),sizeof(Univcoord_T *));
  *nelts = (int *) CALLOC((heapsize + 1),sizeof(int));
  
  heapi = heapsize;
  /* Handle individual contents: Start with i value from before */
  while (i < nstreams) {
    n = (*nelts)[heapi] = streamsize_array[order[i]];
    
    /* Copy to make the merging process non-destructive */
    out = heap[heapi] = MALLOC_ALIGN(n*sizeof(Univcoord_T));
    CHECK_ALIGN(heap[heapi]);
    memcpy(heap[heapi],stream_array[order[i]],n*sizeof(Univcoord_T));

#ifdef DEBUG
    printf("Assigning node %d with %d elts:",heapi,(*nelts)[heapi]);
    for (k = 0; k < (*nelts)[heapi]; k++) {
      printf(" %u",heap[heapi][k]);
    }
    printf("\n");
#endif
    heapi--;
    i++;
  }

  /* Handle combined contents */
  for (i = 0; i < ncombined; i++) {
    heap[heapi] = combined[i];
    (*nelts)[heapi] = totals[i];
#ifdef DEBUG
    printf("Assigning node %d with %d elts:",heapi,(*nelts)[heapi]);
    for (k = 0; k < (*nelts)[heapi]; k++) {
      printf(" %u",heap[heapi][k]);
    }
    printf("\n");
#endif
    heapi--;
  }

  FREE(totals);
  FREE(combined);
  FREE(order);

  return heap;
}
#endif


#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
/* kmer-search.c calls Merge_diagonals_large with stream_high_list and
   stream_low_list.  localdb.c calls Merge_diagonals_uint4.  
   No procedure calls Merge_diagonals_uint8 with stream_list. */
Univcoord_T *
Merge_diagonals_uint8 (int *nelts1, Univcoord_T **stream_array, int *streamsize_array,
		       int nstreams) {
  Univcoord_T *result, **heap, *stream;
  int *nelts;
  int ncopied, heapi, heapsize, base, ancestori, pyramid_start, pyramid_end;
  int bits;
#ifdef DEBUG
  int i;
#endif


  if (nstreams == 0) {
    *nelts1 = 0;
    return (Univcoord_T *) NULL;

  } else if (nstreams == 1) {
    *nelts1 = streamsize_array[0];
    stream = stream_array[0];

    result = MALLOC_ALIGN((*nelts1)*sizeof(Univcoord_T)); /* Output must be aligned */
    memcpy(result,stream,(*nelts1)*sizeof(Univcoord_T));
    return result;

  } else {
    heap = make_univdiagonals_heap(&nelts,&ncopied,stream_array,streamsize_array,nstreams);
    if (ncopied == 1) {
      *nelts1 = nelts[1];
      result = heap[1];
    } else {
      heapsize = 2*ncopied - 1;	/* also index of last node */
#ifdef HAVE_BUILTIN_CLZ
      bits = 31 - __builtin_clz((unsigned int) heapsize);
#elif defined(HAVE_ASM_BSR)
      asm("bsr %1,%0" : "=r"(bits) : "r"(heapsize));
#else
      bits = 31 - ((heapsize >> 16) ? clz_table[heapsize >> 16] : 16 + clz_table[heapsize]); 
#endif

      base = (1 << bits);
      debug(printf("nstreams %d, ncopied %d, heapsize %d, clz %d, bits %d, base %d\n",
		   nstreams,ncopied,heapsize,__builtin_clz(heapsize),bits,base));
      
      /* Middle pyramids */
      while (base > PYRAMID_SIZE) {
	for (pyramid_start = 2*base - PYRAMID_SIZE, pyramid_end = 2*base - 1; pyramid_start >= base;
	     pyramid_start -= PYRAMID_SIZE, pyramid_end -= PYRAMID_SIZE) {
	  debug(printf("diagonals: pyramid_start %d, pyramid_end %d, ncopied %d\n",pyramid_start,pyramid_end,ncopied));
	  ancestori = pyramid_merge(heap,ncopied,heapsize,nelts,pyramid_start,pyramid_end);
	}
	base = ancestori;
      }

      /* Last pyramid */
      pyramid_start = base;
      pyramid_end = 2*base - 1;
      debug(printf("diagonals: pyramid_start %d, pyramid_end %d, ncopied %d\n",pyramid_start,pyramid_end,ncopied));
      /* base = */ pyramid_merge(heap,ncopied,heapsize,nelts,pyramid_start,pyramid_end);

      *nelts1 = nelts[1];
      result = heap[1];

      for (heapi = heapsize; heapi > heapsize - ncopied; heapi--) {
	FREE_ALIGN(heap[heapi]);
      }
    }

    FREE(heap);
    FREE(nelts);
  }


#ifdef DEBUG
  printf("Merge_diagonals_uint8 returning result of length %d\n",*nelts1);
  for (i = 0; i < *nelts1; i++) {
    printf("%llu\n",result[i]);
  }
#endif

  return result;
}
#endif


