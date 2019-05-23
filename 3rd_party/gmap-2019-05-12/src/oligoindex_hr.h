/* $Id: oligoindex_hr.h 218167 2019-01-17 06:23:27Z twu $ */
#ifndef OLIGOINDEX_HR_INCLUDED
#define OLIGOINDEX_HR_INCLUDED

#include "bool.h"
#include "types.h"
#include "univcoord.h"

#include "mode.h"
#include "genomicpos.h"
#include "list.h"
#include "diagpool.h"


#if 0
/* Old code, no longer used */
#define OVERABUNDANCE_CHECK 50
#define OVERABUNDANCE_PCT 0.97
#define OVERABUNDANCE_MIN 200
#endif

/* Attempted to use int, so we could use i32gather_epi32.  However, SIMD is much faster on bytes than on ints */
typedef unsigned char Count_T;
typedef unsigned char Inquery_T;
#define INQUERY_FALSE 0x00
#define INQUERY_TRUE  0xFF
#define INCR_COUNT(counts) counts += 1;

#if defined(HAVE_AVX512)
#define SIMD_NELTS 64		/* 64 bytes in 256 bits */
#elif defined(HAVE_AVX2)
#define SIMD_NELTS 32		/* 32 bytes in 256 bits */
#elif defined(HAVE_SSE2)
#define SIMD_NELTS 16		/* 16 bytes in 128 bits */
#endif


/* These need to be here, because the struct is here */
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSSE3
#include <tmmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#ifdef HAVE_AVX512
#include <immintrin.h>
#endif


#define T Oligoindex_T
typedef struct T *T;
struct T {

  int indexsize;
  Shortoligomer_T mask;

  int diag_lookback;
  int suffnconsecutive;

  /* bool query_evaluated_p; */

  Oligospace_T oligospace;
#if defined(HAVE_AVX512)
  __m512i *inquery_allocated;
  __m512i *counts_allocated;
#elif defined(HAVE_AVX2)
  __m256i *inquery_allocated;
  __m256i *counts_allocated;
#elif defined(HAVE_SSE2)
  __m128i *inquery_allocated;
  __m128i *counts_allocated;
#endif
  Inquery_T *inquery;
  Count_T *counts;

  Chrpos_T *table;
  UINT4 *positions;
  /* UINT4 *pointers; */
  /* UINT4 *pointers_allocated; */
};



typedef struct Oligoindex_array_T *Oligoindex_array_T;

extern void
Oligoindex_hr_setup (Genomecomp_T *ref_blocks_in, Mode_T mode_in);

extern int
Oligoindex_indexsize (T this);

extern int
Oligoindex_array_length (Oligoindex_array_T oligoindices);
extern T
Oligoindex_array_elt (Oligoindex_array_T oligoindices, int source);

extern Oligoindex_array_T
Oligoindex_array_new_major (int max_querylength, int max_genomiclength);

extern Oligoindex_array_T
Oligoindex_array_new_minor (int max_querylength, int max_genomiclength);

extern Chrpos_T *
Oligoindex_allocate_positions (UINT4 *__restrict__ positions,
			       Inquery_T *__restrict__ inquery, Count_T *counts, int oligospace);


extern double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, char *queryuc_ptr, int querystart, int queryend, bool trimp);
extern void
Oligoindex_hr_tally (T this, Univcoord_T mappingstart, Univcoord_T mappingend, bool plusp,
		     char *queryuc_ptr, int querystart, int queryend, Chrpos_T chrpos, int genestrand);
extern void
Oligoindex_untally (T this);
extern void
Oligoindex_clear_inquery (T this, char *queryuc_ptr, int querystart, int queryend);
extern void
Oligoindex_array_free(Oligoindex_array_T *old);

extern List_T
Oligoindex_get_mappings (List_T diagonals, bool *coveredp, Chrpos_T **mappings, int *npositions,
			 int *totalpositions, bool *oned_matrix_p, int *maxnconsecutive, 
			 Oligoindex_array_T array, T this, char *queryuc_ptr,
			 int querystart, int queryend, int querylength,
			 Chrpos_T chrstart, Chrpos_T chrend,
			 Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
			 Diagpool_T diagpool);

#undef T
#endif

