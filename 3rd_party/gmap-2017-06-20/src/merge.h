#ifndef MERGE_INCLUDED
#define MERGE_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "types.h"
#include "list.h"
#include "intlist.h"


/* Pad lengths at end for row-based storage */
#ifdef HAVE_AVX512
#define PAD_UINT4(x) (((x + 15)/16) * 16)
#elif defined(HAVE_AVX2)
#define PAD_UINT4(x) (((x + 7)/8) * 8)
#elif defined(HAVE_SSE4_1)
#define PAD_UINT4(x) (((x + 3)/4) * 4)
#else
#define PAD_UINT4(x) (x)
#endif


typedef struct Record_T *Record_T;
struct Record_T {
  Univcoord_T diagonal;		/* Primary sort */
  int querypos;			/* Secondary sort */
};


extern unsigned int *
Merge_uint4 (unsigned int *__restrict__ dest, unsigned int *__restrict__ A,
	     unsigned int *__restrict__ B, int nA, int nB);

extern UINT4 *
Merge_diagonals (int *nelts1, List_T stream_list, Intlist_T streamsize_list);

extern Record_T *
Merge_records (int *nelts1, List_T stream_list, Intlist_T streamsize_list,
	       Intlist_T querypos_list, Intlist_T diagterm_list,
	       struct Record_T *all_records);

#endif


