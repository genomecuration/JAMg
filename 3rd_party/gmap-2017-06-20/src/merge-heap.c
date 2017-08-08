static char rcsid[] = "$Id: merge-heap.c 201745 2016-12-16 16:51:24Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "merge-heap.h"
#include "assert.h"
#include "mem.h"
#include "popcount.h"

#include <stdio.h>
#include <stdlib.h>

#define PYRAMID_SIZE 32
#define KEY_MASK (~0U << 5)


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

static void
min_heap_insert (UINT4 *heap, int *heapsize, UINT4 diagonal) {
  int i;

  i = ++(*heapsize);
  while (i > 1 && (heap[PARENT(i)] > diagonal)) {
    heap[i] = heap[PARENT(i)];
    i = PARENT(i);
  }
  heap[i] = diagonal;

  return;
}


/* Provide ancestori as inserti */
static void
heapify (unsigned int *heap, unsigned int diagonal, int merge_heap_size) {
  int inserti, smallesti, righti;
  int i;

  debug6(printf("Starting heapify with %llu\n",(unsigned long long) diagonal));
#ifdef DEBUG6
  for (i = 1; i <= 2*merge_heap_size + 1; i++) {
    printf("%d %u\n",i,heap[i]);
  }
  printf("\n");
#endif

  inserti = 1;
  smallesti = (heap[3] < heap[2]) ? 3 : 2;
  debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		2,3,(unsigned long long) heap[2],(unsigned long long)heap[3]));
  while (diagonal > heap[smallesti]) {
    heap[inserti] = heap[smallesti];
    inserti = smallesti;
    smallesti = LEFT(inserti);
    righti = smallesti+1;
    debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		  smallesti,righti,(unsigned long long) heap[smallesti],
		  (unsigned long long) heap[righti]));
    if (heap[righti] < heap[smallesti]) {
      smallesti = righti;
    }
  }
  heap[inserti] = diagonal;
  debug6(printf("Inserting at %d\n\n",inserti));
  return;
}


static int
pyramid_merge_full (Record_T **record_heap, unsigned int **key_streams, unsigned int *merge_heap,
		    int node_start, int ancestori, int merge_heap_size) {
  int nelts = 0;
  unsigned int diagonal;
  int streami, k;
  int ptrs[PYRAMID_SIZE];

  k = 0;
  memset(ptrs,0,PYRAMID_SIZE*sizeof(int));
  while ((diagonal = merge_heap[1]) < -1U) {
    /* Convert integer to structure */
    streami = diagonal & ~KEY_MASK;
    record_heap[ancestori][k++] = record_heap[node_start + streami][ptrs[streami]];
    debug(printf("Writing %u (stream %d): %u\n",diagonal,streami,record_heap[ancestori][k-1]->diagonal));

    /* Advance pointer and get next value */
    diagonal = key_streams[streami][++ptrs[streami]];
    heapify(merge_heap,diagonal,merge_heap_size);
    nelts += 1;
  }

  return nelts;
}


static Record_T **
make_record_heap (int **nelts, List_T stream_list, Intlist_T streamsize_list, 
		  Intlist_T querypos_list, Intlist_T diagterm_list, int nstreams,
		  int base, struct Record_T *all_records) {
  Record_T **record_heap;
  UINT4 *diagonals;
  int heapsize, null_pyramid_start, heapi, basei;
  int querypos, diagterm;
  int i, k;

  heapsize = 2*nstreams - 1;
  null_pyramid_start = (heapsize + PYRAMID_SIZE - 1)/PYRAMID_SIZE * PYRAMID_SIZE; /* full or partial pyramid for entries below this */

  /* Add 4 to handle partial pyramid */
  record_heap = (Record_T **) CALLOC(heapsize + PYRAMID_SIZE,sizeof(Record_T *));
  *nelts = (int *) CALLOC(heapsize + PYRAMID_SIZE,sizeof(int));

  /* Process as (base - 1) downto nstreams, then heapsize downto base,
     because stream_list is in reverse order of elts */
  k = 0;
  for (heapi = base - 1; heapi >= PARENT(null_pyramid_start); heapi--) {
    /* Put all information into penultimate row */
    stream_list = List_pop(stream_list,(void *) &diagonals); /* already padded */
    streamsize_list = Intlist_pop(streamsize_list,&((*nelts)[heapi]));
    querypos_list = Intlist_pop(querypos_list,&querypos);
    diagterm_list = Intlist_pop(diagterm_list,&diagterm);
    record_heap[heapi] = (Record_T *) MALLOC(((*nelts)[heapi]) * sizeof(Record_T));
    debug(printf("NULL: Assigning node %d with %d elts (%p)",heapi,(*nelts)[heapi],record_heap[heapi]));

    for (i = 0; i < (*nelts)[heapi]; i++) {
      /* Process in forward order to keep records in order */
      all_records[k].diagonal = diagonals[i] + diagterm;
      all_records[k].querypos = querypos;
      record_heap[heapi][i] = &(all_records[k]);
      debug(printf(" %u+%d",diagonals[i],querypos));
      k++;
    }
    debug(printf("\n"));
  }
    
  for ( ; heapi >= nstreams; heapi--) {
    /* Move all information down to left child */
    basei = LEFT(heapi);
    stream_list = List_pop(stream_list,(void *) &diagonals); /* already padded */
    streamsize_list = Intlist_pop(streamsize_list,&((*nelts)[basei]));
    querypos_list = Intlist_pop(querypos_list,&querypos);
    diagterm_list = Intlist_pop(diagterm_list,&diagterm);
    record_heap[basei] = (Record_T *) MALLOC(((*nelts)[basei]) * sizeof(Record_T));
    debug(printf("PART: Assigning node %d => %d with %d elts (%p)",heapi,basei,(*nelts)[basei],record_heap[basei]));

    for (i = 0; i < (*nelts)[basei]; i++) {
      /* Process in forward order to keep records in order */
      all_records[k].diagonal = diagonals[i] + diagterm;
      all_records[k].querypos = querypos;
      record_heap[basei][i] = &(all_records[k]);
      debug(printf(" %u+%d",diagonals[i],querypos));
      k++;
    }
    debug(printf("\n"));
  }

  for (heapi = heapsize; heapi >= base; heapi--) {
    /* Put all information into base row */
    stream_list = List_pop(stream_list,(void *) &diagonals); /* already padded */
    streamsize_list = Intlist_pop(streamsize_list,&((*nelts)[heapi]));
    querypos_list = Intlist_pop(querypos_list,&querypos);
    diagterm_list = Intlist_pop(diagterm_list,&diagterm);
    record_heap[heapi] = (Record_T *) MALLOC(((*nelts)[heapi]) * sizeof(Record_T));
    debug(printf("FULL: Assigning node %d with %d elts (%p)",heapi,(*nelts)[heapi],record_heap[heapi]));

    for (i = 0; i < (*nelts)[heapi]; i++) {
      /* Process in forward order to keep records in order */
      all_records[k].diagonal = diagonals[i] + diagterm;
      all_records[k].querypos = querypos;
      record_heap[heapi][i] = &(all_records[k]);
      debug(printf(" %u+%d",diagonals[i],querypos));
      k++;
    }
    debug(printf("\n"));
  }

  return record_heap;
}


/* For initializing heap, there are three categories:
   base..(heapsize % PYRAMID_SIZE) + PYRAMID_SIZE: Fill bottom row
   straddling heapsize: Pull down some nodes to bottom row
   heapsize..(2*base - 1): Fill penultimate row */
Record_T *
Merge_records_heap (int *nelts1, List_T stream_list, Intlist_T streamsize_list,
		    Intlist_T querypos_list, Intlist_T diagterm_list, 
		    struct Record_T *all_records) {
  Record_T *result, **record_heap, curr;
  UINT4 *key_streams[PYRAMID_SIZE];
  UINT4 merge_heap[2*PYRAMID_SIZE+1+1]; /* Add second 1 because top node is at 1 */
  UINT4 *storage;
  int *nelts, nalloc;
  int nstreams, heapsize, base, ancestori, pyramid_start, pyramid_end,
    node_start, node_end, start, end;
  int merge_heap_size;
  int bits;
  int heapi, streami, i, j;

  debug(printf("Entered Merge_records\n"));

  if ((nstreams = List_length(stream_list)) == 0) {
    *nelts1 = 0;
    return (Record_T *) NULL;

  } else {
    heapsize = 2*nstreams - 1;	/* also index of last node */
#ifdef HAVE_BUILTIN_CLZ
    bits = 31 - __builtin_clz(heapsize);
#elif defined(HAVE_ASM_BSR)
    asm("bsr %1,%0" : "=r"(bits) : "r"(heapsize));
#else
    bits = 31 - ((heapsize >> 16) ? clz_table[heapsize >> 16] : 16 + clz_table[heapsize]); 
#endif
    base = (1 << bits);
    debug(printf("nstreams %d, heapsize %d, base %d\n",nstreams,heapsize,base));
    record_heap = make_record_heap(&nelts,stream_list,streamsize_list,querypos_list,diagterm_list,
				   nstreams,base,all_records);
  }


  while (base > 1) {
    if (base < PYRAMID_SIZE) {
      pyramid_start = base;
      pyramid_end = 2*base - 1;

      ancestori = 1;
      debug(printf("records: pyramid_start %d, pyramid_end %d, nstreams %d\n",pyramid_start,pyramid_end,nstreams));

      /* Allocate memory for the pyramid key_streams */
      nalloc = 0;
      for (heapi = pyramid_start; heapi <= pyramid_end; heapi++) {
	nalloc += (nelts[heapi] + 1);
      }
      storage = (UINT4 *) MALLOC(nalloc * sizeof(UINT4));

      /* Convert structures to integers (key_streams) */
      nalloc = 0;
      merge_heap_size = 0;
      for (heapi = pyramid_start, streami = 0; heapi <= pyramid_end; heapi++, streami++) {
	key_streams[streami] = &(storage[nalloc]);
	for (i = 0; i < nelts[heapi]; i++) {
	  key_streams[streami][i] = (record_heap[heapi][i]->diagonal & KEY_MASK) + streami;
	}
	key_streams[streami][i] = -1U;
	nalloc += (i + 1);	/* nelts[heapi] + 1 */

	min_heap_insert(merge_heap,&merge_heap_size,key_streams[streami][0]);
      }

      /* Set up bounds of heap (sentinels) */
      assert(merge_heap_size <= PYRAMID_SIZE);
      debug(printf("merge_heap_size is %d\n",merge_heap_size));
      for (i = merge_heap_size+1; i <= 2*merge_heap_size+1; i++) {
	merge_heap[i] = -1U;
      }

      /* Merge and convert integers to structures */
      record_heap[1] = (Record_T *) MALLOC(nalloc * sizeof(Record_T));
      nelts[1] = pyramid_merge_full(record_heap,key_streams,merge_heap,pyramid_start,ancestori,merge_heap_size);

      /* Free base heaps */
      for (heapi = pyramid_start, streami = 0; heapi <= pyramid_end; heapi++, streami++) {
	FREE(record_heap[pyramid_start + streami]);
      }

      /* Free key_streams storage */
      FREE(storage);

    } else {
      for (pyramid_start = 2*base - PYRAMID_SIZE, pyramid_end = 2*base - 1; pyramid_start >= base;
	   pyramid_start -= PYRAMID_SIZE, pyramid_end -= PYRAMID_SIZE) {
	debug(printf("records: pyramid_start %d, pyramid_end %d, nstreams %d",pyramid_start,pyramid_end,nstreams));

	if (pyramid_start > heapsize) {
	  node_start = PARENT(pyramid_start);
	  node_end = PARENT(pyramid_end);
	  debug(printf(" => node_start %d, node_end %d\n",node_start,node_end));
	} else {
	  node_start = pyramid_start;
	  node_end = pyramid_end;
	}
	debug(printf("\n"));

	/* Determine ancestori */
	start = node_start;
	end = node_end;
	while ((start = PARENT(start)) < (end = PARENT(end))) ;
	ancestori = start;

	/* Allocate memory for the pyramid key_streams */
	nalloc = 0;
	for (heapi = node_start; heapi <= node_end; heapi++) {
	  nalloc += (nelts[heapi] + 1);
	}
	storage = (UINT4 *) MALLOC(nalloc * sizeof(UINT4));

	/* Convert structures to integers (key_streams) */
	nalloc = 0;
	merge_heap_size = 0;
	for (heapi = node_start, streami = 0; heapi <= node_end; heapi++, streami++) {
	  key_streams[streami] = &(storage[nalloc]);
	  for (i = 0; i < nelts[heapi]; i++) {
	    key_streams[streami][i] = (record_heap[heapi][i]->diagonal & KEY_MASK) + streami;
	  }
	  key_streams[streami][i] = -1U;
	  nalloc += (i + 1);	/* nelts[heapi] + 1 */

	  min_heap_insert(merge_heap,&merge_heap_size,key_streams[streami][0]);
	}

#ifdef DEBUG
	for (heapi = node_start, streami = 0; heapi <= node_end; heapi++, streami++) {
	  printf("key_stream %d:",streami);
	  for (i = 0; i <= nelts[heapi]; i++) {
	    printf(" %u",key_streams[streami][i]);
	  }
	  printf("\n");
	}
#endif

	/* Set up bounds of heap (sentinels) */
	assert(merge_heap_size <= PYRAMID_SIZE);
	debug(printf("merge_heap_size is %d\n",merge_heap_size));
	for (i = merge_heap_size+1; i <= 2*merge_heap_size+1; i++) {
	  merge_heap[i] = -1U;
	}

	/* Merge and convert integers to structures */
	record_heap[ancestori] = (Record_T *) MALLOC(nalloc * sizeof(Record_T));
	nelts[ancestori] = pyramid_merge_full(record_heap,key_streams,merge_heap,node_start,ancestori,merge_heap_size);

	/* Free base heaps */
	for (heapi = node_start; heapi <= node_end; heapi++) {
	  FREE(record_heap[heapi]);
	}

	/* Free key_streams storage */
	FREE(storage);
      }
    }

    base = ancestori;
  }

  *nelts1 = nelts[1];
  result = record_heap[1];

  /* Final insertion sort to correct for truncation of keys */
  for (j = 1; j < *nelts1; j++) {
    curr = result[j];
    i = j - 1;
    /* For a stable merge sort, is the second condition possible? */
    while (i >= 0 && (result[i]->diagonal > curr->diagonal ||
		     (result[i]->diagonal == curr->diagonal &&
		      result[i]->querypos > curr->querypos))) {
      assert(result[i]->diagonal > curr->diagonal);
      result[i+1] = result[i];
      i--;
    }
    result[i+1] = curr;
  }


  FREE(nelts);
  FREE(record_heap);

#ifdef DEBUG0
  printf("Merge_records returning result of length %d\n",*nelts1);
  for (i = 0; i < *nelts1; i++) {
    printf("%u %d\n",result[i]->diagonal,result[i]->querypos);
  }
#endif

  return result;
}


