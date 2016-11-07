static char rcsid[] = "$Id: sarray-read.c 197551 2016-09-08 01:16:14Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "sarray-read.h"

#ifdef WORDS_BIGENDIAN
#define CONVERT(x) Bigendian_convert_uint(x)
#include "bigendian.h"
#else
#define CONVERT(x) (x)
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */

#include "mem.h"
#include "bool.h"
#include "assert.h"
#include "access.h"
#include "types.h"
#include "listdef.h"
#include "list.h"
#include "genome128_hr.h"
#include "splice.h"
#include "indel.h"
#include "stage3hr.h"
#include "bytecoding.h"
#include "bitpack64-read.h"
#include "bitpack64-readtwo.h"
#include "bitpack64-access.h"

#include "comp.h"
#include "diagdef.h"
#include "diag.h"
#include "univdiagdef.h"
#include "univdiag.h"
#include "substring.h"
#include "junction.h"
#include "stage3hr.h"
#include "sedgesort.h"


#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#else
#include <emmintrin.h>
#endif
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSSE3)
#else
#include <tmmintrin.h>
#endif
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_AVX2)
#else
#include <immintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip popcnt */
#elif defined(HAVE_POPCNT)
#include <immintrin.h>
#elif defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif


#define MIN_ENDLENGTH 12
#define MIN_INTRONLEN 9

#define MAX_HITS_FOR_BEST_ELT 1000

/* A value of 10000 misses various splices, although they are caught by GSNAP algorithm */
#define EXCESS_SARRAY_HITS 100000
#define LOCALSPLICING_NMATCHES_SLOP 1
#define LOCALSPLICING_PROB_SLOP 0.05

#define USE_SHUFFLE_MASK 1	/* Alternative requires AVX, and that part of the code isn't called much */
/* #define REQUIRE_ALIGNMENT 1 */

#define GUESS_ALLOCATION 10

/* #define USE_SEPARATE_BUCKETS 1 */

/* Results of each suffix array search */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define MAX_DEBUG1_HITS 100

/* Details of suffix array search */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Search through saindex */
#ifdef DEBUG1A
#define debug1a(x) x
#else
#define debug1a(x)
#endif

/* get_child */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Compressed suffix array */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Compressed suffix array: comparison with sarray */
#ifdef DEBUG3A
#define debug3a(x) x
#else
#define debug3a(x)
#endif

/* Compressed suffix array: comparison with csa phi */
#ifdef DEBUG3B
#define debug3b(x) x
#else
#define debug3b(x)
#endif

/* known splicing */
#ifdef DEBUG4S
#define debug4s(x) x
#else
#define debug4s(x)
#endif

/* find_multimiss_iter */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* find_multimiss_iter details */
#ifdef DEBUG7A
#define debug7a(x) x
#else
#define debug7a(x)
#endif

/* SIMD new filtering */
#ifdef DEBUG7B
#define debug7b(x) x
#else
#define debug7b(x)
#endif


/* Comparing SIMD with non-SIMD */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Sorting of diagonals */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* GMAP */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif

/* Oligoindex fillin */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

/* Compare separate buckets with a single one */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif


#ifdef DEBUG7B
static void
print_vector_hex (__m128i x) {
  UINT4 *s = (UINT4 *) &x;

  /* printf("%08X %08X %08X %08X\n",s[0],s[1],s[2],s[3]); */
  printf("%08X %08X %08X %08X\n",s[3],s[2],s[1],s[0]);
  return;
}

static void
print_vector_uint (__m128i x) {
  UINT4 *s = (UINT4 *) &x;

  /* printf("%d %d %d %d\n",s[0],s[1],s[2],s[3]); */
  printf("%u %u %u %u\n",s[3],s[2],s[1],s[0]);
  return;
}

#ifdef HAVE_AVX2
static void
print_vector_hex_256 (__m256i x) {
  UINT4 *s = (UINT4 *) &x;

  /* printf("%d %d %d %d\n",s[0],s[1],s[2],s[3]); */
  printf("%08X %08X %08X %08X %08X %08X %08X %08X\n",s[7],s[6],s[5],s[4],s[3],s[2],s[1],s[0]);
  return;
}

static void
print_vector_uint_256 (__m256i x) {
  UINT4 *s = (UINT4 *) &x;

  /* printf("%d %d %d %d\n",s[0],s[1],s[2],s[3]); */
  printf("%u %u %u %u %u %u %u %u\n",s[7],s[6],s[5],s[4],s[3],s[2],s[1],s[0]);
  return;
}
#endif
#endif



#define T Sarray_T
struct T {
  Univcoord_T n;
  Univcoord_T n_plus_one;

  /* Old format */
  int array_shmid;
  key_t array_key;
  Univcoord_T *array;

  int lcpchilddc_shmid;
  key_t lcpchilddc_key;
  unsigned char *lcpchilddc;

  int lcp_guide_shmid;
  key_t lcp_guide_key;
  int lcp_exceptions_shmid;
  key_t lcp_exceptions_key;
  UINT4 *lcp_guide;
  UINT4 *lcp_exceptions;
  int n_lcp_exceptions;		/* Won't be necessary if we change lcpchilddc to use guide array */
  /* int lcp_guide_interval; -- Always use 1024 */
  
  int child_guide_shmid;
  key_t child_guide_key;
  int child_exceptions_shmid;
  key_t child_exceptions_key;
  UINT4 *child_guide;
  UINT4 *child_exceptions;
  /* int n_child_exceptions; */
  int child_guide_interval; /* Always use 1024 */

#if 0
  Sarrayptr_T initindexi[4];	/* For A, C, G, T */
  Sarrayptr_T initindexj[4];	/* For A, C, G, T */
#endif

  int indexsize;
  UINT4 indexspace;		/* 4^indexsize.  Used by sarray_search to detect when we have a poly-T oligo shorter than indexsize */
#ifdef DEBUG15
  UINT4 *indexi_ptrs, *indexi_comp, *indexj_ptrs, *indexj_comp; /* bucket array: oligomer lookup into suffix array */
  UINT4 *indexij_ptrs, *indexij_comp;
#elif defined(USE_SEPARATE_BUCKETS)
  UINT4 *indexi_ptrs, *indexi_comp, *indexj_ptrs, *indexj_comp; /* bucket array: oligomer lookup into suffix array */
#else
  int indexij_ptrs_shmid;
  key_t indexij_ptrs_key;
  int indexij_comp_shmid;
  key_t indexij_comp_key;
  UINT4 *indexij_ptrs, *indexij_comp;
#endif

  Access_T array_access; int array_fd; size_t array_len;

#ifdef DEBUG15
  int indexi_ptrs_fd; size_t indexi_ptrs_len; int indexi_comp_fd; size_t indexi_comp_len;
  int indexj_ptrs_fd; size_t indexj_ptrs_len; int indexj_comp_fd; size_t indexj_comp_len;
  int indexij_ptrs_fd; size_t indexij_ptrs_len; int indexij_comp_fd; size_t indexij_comp_len;
#elif defined(USE_SEPARATE_BUCKETS)
  int indexi_ptrs_fd; size_t indexi_ptrs_len; int indexi_comp_fd; size_t indexi_comp_len;
  int indexj_ptrs_fd; size_t indexj_ptrs_len; int indexj_comp_fd; size_t indexj_comp_len;
#else
  Access_T indexij_ptrs_access; int indexij_ptrs_fd; size_t indexij_ptrs_len;
  Access_T indexij_comp_access; int indexij_comp_fd; size_t indexij_comp_len;
#endif

  Access_T lcpchilddc_access; int lcpchilddc_fd; size_t lcpchilddc_len;

  Access_T lcp_guide_access; int lcp_guide_fd; size_t lcp_guide_len;
  Access_T lcp_exceptions_access; int lcp_exceptions_fd; size_t lcp_exceptions_len;

  Access_T child_guide_access; int child_guide_fd; size_t child_guide_len;
  Access_T child_exceptions_access; int child_exceptions_fd; size_t child_exceptions_len;

};


/* For benchmarking */
Univcoord_T
Sarray_size (Sarray_T this) {
  return this->n_plus_one;
}


static Sarray_T sarray_fwd;
static Sarray_T sarray_rev;
static Genome_T genome;
static bool *circularp;

static char conversion_fwd[128];
static char conversion_rev[128];

static Univ_IIT_T chromosome_iit;
static int circular_typeint;
static int splicing_penalty;

static Chrpos_T overall_max_distance;
static Chrpos_T shortsplicedist;
static Chrpos_T max_deletionlen;
static Chrpos_T max_insertionlen_default;
static int max_end_deletions;
static int max_middle_insertions_default;

/* Splicing */
static Univcoord_T *splicesites;
static Splicetype_T *splicetypes;
static Chrpos_T *splicedists;
static int nsplicesites;


#if defined(HAVE_AVX2) && !defined(WORDS_BIGENDIAN)
static __m256i epi32_convert_256;	/* For converting unsigned ints to signed ints */
#endif

#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
static __m128i epi32_convert;	/* For converting unsigned ints to signed ints */
#endif

#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN) && defined(USE_SHUFFLE_MASK)
static __m128i shuffle_mask16[16];
#endif


#if 0
/* Simplified from sarray_search_simple in sarray-write.c */
static void
sarray_search_char (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, char desired_char,
		    UINT4 *SA, UINT4 n, char *chartable) {
  Sarrayptr_T low, high, mid;
  Univcoord_T pos;
  char c;

  low = 1;
  high = n + 1;

  while (low < high) {
#if 0
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif
#ifdef WORDS_BIGENDIAN
    pos = Bigendian_convert_uint(SA[mid]);
#else
    pos = SA[mid];
#endif
    c = Genome_get_char_lex(genome,pos,n,chartable);
    if (desired_char > c) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  *initptr = low;

  low--;
  high = n;
  while (low < high) {
#if 1
    /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 || high % 2 == 1) {
      mid += 1;
    }
#else
    /* This does not work for ceiling */
    mid = low + ((high - low) / 2);
#endif
#ifdef WORDS_BIGENDIAN
    pos = Bigendian_convert_uint(SA[mid]);
#else
    pos = SA[mid];
#endif
    c = Genome_get_char_lex(genome,pos,n,chartable);
    if (desired_char >= c) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  *finalptr = high;
  return;
}
#endif


void
Sarray_setup (T sarray_fwd_in, T sarray_rev_in, Genome_T genome_in, Mode_T mode,
	      Univ_IIT_T chromosome_iit_in, int circular_typeint_in, bool *circularp_in,
	      Chrpos_T shortsplicedist_in, int splicing_penalty_in,
	      int max_deletionlength, int max_end_deletions_in,
	      int max_middle_insertions_in, int max_end_insertions,
	      Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
	      Chrpos_T *splicedists_in, int nsplicesites_in) {
  int i;

  sarray_fwd = sarray_fwd_in;
  sarray_rev = sarray_rev_in;
  genome = genome_in;
  circularp = circularp_in;

  for (i = 0; i < 128; i++) {
    conversion_fwd[i] = i;
    conversion_rev[i] = i;
  }
  if (mode == STANDARD) {
    /* Don't change conversion */
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    conversion_fwd['C'] = 'T';	/* CT */
    conversion_rev['G'] = 'A';	/* GA */
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    conversion_fwd['A'] = 'G';	/* AG */
    conversion_rev['T'] = 'C';	/* TC */
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    conversion_fwd['T'] = 'C';	/* TC */
    conversion_rev['A'] = 'G';	/* AG */
  }

  chromosome_iit = chromosome_iit_in;
  circular_typeint = circular_typeint_in;
  shortsplicedist = shortsplicedist_in;
  splicing_penalty = splicing_penalty_in;

  max_deletionlen = max_deletionlength;
  max_end_deletions = max_end_deletions_in;
  max_middle_insertions_default = max_middle_insertions_in;
  if (max_middle_insertions_in > max_end_insertions) {
    max_insertionlen_default = max_middle_insertions_in;
  } else {
    max_insertionlen_default = max_end_insertions;
  }

  if (shortsplicedist > max_deletionlen) {
    overall_max_distance = shortsplicedist;
  } else {
    overall_max_distance = max_deletionlen;
  }

  splicesites = splicesites_in;
  splicetypes = splicetypes_in;
  splicedists = splicedists_in;
  nsplicesites = nsplicesites_in;

#if 0
  sarray_search_char(&(sarray->initindexi[0]),&(sarray->initindexj[0]),/*desired_char*/'A',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[1]),&(sarray->initindexj[1]),/*desired_char*/'C',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[2]),&(sarray->initindexj[2]),/*desired_char*/'G',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[3]),&(sarray->initindexj[3]),/*desired_char*/'T',sarray->array,sarray->n);
#endif

#if 0
  printf("A => %u %u\n",sarray->initindexi[0],sarray->initindexj[0]);
  printf("C => %u %u\n",sarray->initindexi[1],sarray->initindexj[1]);
  printf("G => %u %u\n",sarray->initindexi[2],sarray->initindexj[2]);
  printf("T => %u %u\n",sarray->initindexi[3],sarray->initindexj[3]);
#endif

#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
  epi32_convert = _mm_set1_epi32(2147483648); /* 2^31 */
#endif

#if defined(HAVE_AVX2) && !defined(WORDS_BIGENDIAN)
  epi32_convert_256 = _mm256_set1_epi32(2147483648); /* 2^31 */
#endif

#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN) && defined(USE_SHUFFLE_MASK)
  /* Used by fill_positions_filtered_first */
  shuffle_mask16[0] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1);
  shuffle_mask16[1] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,  3, 2, 1, 0);
  shuffle_mask16[2] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,  7, 6, 5, 4);
  shuffle_mask16[3] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1,  7, 6, 5, 4,  3, 2, 1, 0);
  shuffle_mask16[4] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, 11,10, 9, 8);
  shuffle_mask16[5] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 11,10, 9, 8,  3, 2, 1, 0);
  shuffle_mask16[6] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 11,10, 9, 8,  7, 6, 5, 4);
  shuffle_mask16[7] =  _mm_set_epi8(-1,-1,-1,-1, 11,10, 9, 8,  7, 6, 5, 4,  3, 2, 1, 0);
  shuffle_mask16[8] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, 15,14,13,12);
  shuffle_mask16[9] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 15,14,13,12,  3, 2, 1, 0);
  shuffle_mask16[10] = _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 15,14,13,12,  7, 6, 5, 4);
  shuffle_mask16[11] = _mm_set_epi8(-1,-1,-1,-1, 15,14,13,12,  7, 6, 5, 4,  3, 2, 1, 0);
  shuffle_mask16[12] = _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 15,14,13,12, 11,10, 9, 8);
  shuffle_mask16[13] = _mm_set_epi8(-1,-1,-1,-1, 15,14,13,12, 11,10, 9, 8,  3, 2, 1, 0);
  shuffle_mask16[14] = _mm_set_epi8(-1,-1,-1,-1, 15,14,13,12, 11,10, 9, 8,  7, 6, 5, 4);
  shuffle_mask16[15] = _mm_set_epi8(15,14,13,12, 11,10, 9, 8,  7, 6, 5, 4,  3, 2, 1, 0);
#endif
  
  return;
}


static int
log4 (int result) {
  int exponent = 0;

  while (result > 1) {
    result /= 4;
    exponent++;
  }

  return exponent;
}

static UINT4
power (int base, int exponent) {
  UINT4 result = 1;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }

  return result;
}


#if 0
void
Sarray_shmem_remove (char *dir, char *fileroot, char *snps_root, Mode_T mode, bool fwdp) {
  char *mode_prefix;
  char *sarrayfile;
  char *lcpchilddcfile;
  char *lcp_guidefile, *lcp_exceptionsfile;
  char *child_guidefile, *child_exceptionsfile;
  char *indexij_ptrsfile, *indexij_compfile;

  if (mode == STANDARD) {
    mode_prefix = ".";
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".metct.";
    } else {
      mode_prefix = ".metga.";
    }
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2iag.";
    } else {
      mode_prefix = ".a2itc.";
    }
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2itc.";
    } else {
      mode_prefix = ".a2iag.";
    }
  }

  sarrayfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sarray")+1,sizeof(char));
  sprintf(sarrayfile,"%s/%s%ssarray",dir,fileroot,mode_prefix);

  lcpchilddcfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpchilddc")+1,sizeof(char));
  sprintf(lcpchilddcfile,"%s/%s%ssalcpchilddc",dir,fileroot,mode_prefix);

  lcp_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpguide1024")+1,sizeof(char));
  sprintf(lcp_guidefile,"%s/%s%ssalcpguide1024",dir,fileroot,mode_prefix);
  lcp_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpexc")+1,sizeof(char));
  sprintf(lcp_exceptionsfile,"%s/%s%ssalcpexc",dir,fileroot,mode_prefix);

  child_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildguide1024")+1,sizeof(char));
  sprintf(child_guidefile,"%s/%s%ssachildguide1024",dir,fileroot,mode_prefix);
  child_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildexc")+1,sizeof(char));
  sprintf(child_exceptionsfile,"%s/%s%ssachildexc",dir,fileroot,mode_prefix);

  indexij_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64meta")+1,sizeof(char));
  sprintf(indexij_ptrsfile,"%s/%s%ssaindex64meta",dir,fileroot,mode_prefix);
  indexij_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64strm")+1,sizeof(char));
  sprintf(indexij_compfile,"%s/%s%ssaindex64strm",dir,fileroot,mode_prefix);

  Access_shmem_remove(indexij_ptrsfile);
  Access_shmem_remove(indexij_compfile);

  Access_shmem_remove(sarrayfile);
  Access_shmem_remove(lcpchilddcfile);
  Access_shmem_remove(lcp_guidefile);
  Access_shmem_remove(lcp_exceptionsfile);

  Access_shmem_remove(child_guidefile);
  Access_shmem_remove(child_exceptionsfile);

  FREE(child_exceptionsfile);
  FREE(child_guidefile);

  FREE(lcp_exceptionsfile);
  FREE(lcp_guidefile);

  FREE(lcpchilddcfile);

  FREE(sarrayfile);

  return;
}
#endif


#ifdef USE_CSA

static Univcoord_T
csa_lookup (T sarray, Sarrayptr_T i) {
  Univcoord_T nhops = 0, expected_sa_i;
  Sarrayptr_T expected_i;
  __m128i converted, cmp;
  int matchbits;

  debug3(printf("Entered csa_lookup for %u:",i));
#ifdef DEBUG3A
  expected_sa_i = sarray->array[i];
#endif

  if (
#ifdef DEBUG3A
      0 && 
#endif
      sarray->array != NULL) {
    debug3(printf("Returning %u\n",sarray->array[i]));
    return sarray->array[i];
  } else {
    while ((i % sarray->sa_sampling) != 0) {
      debug3(printf(",%u",i));
#ifdef DEBUG3B
      expected_i = sarray->csa[i];
#endif

#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
      converted = _mm_sub_epi32(_mm_set1_epi32(i),epi32_convert);
      cmp = _mm_cmpgt_epi32(converted,sarray->indices0); /* To use cmpgt, sarray->indices0 is shifted down by 1 */
      matchbits = _mm_movemask_ps(_mm_castsi128_ps(cmp));
      /* assert(matchbits == 0 || matchbits == 1 || matchbits == 3 || matchbits == 7 || matchbits == 15); */
      debug3(printf("(%d)",matchbits));
      i = Bitpack64_read_one(i - sarray->index0[matchbits],sarray->csa0ptrs[matchbits],sarray->csa0comp[matchbits]);
#else
      if (i >= sarray->indexX) {
	assert(matchbits == 15);
	printf("X");
	i = Bitpack64_read_one(i-sarray->indexX,sarray->csaXptrs,sarray->csaXcomp);
      } else if (i >= sarray->indexT) {
	assert(matchbits == 7);
	printf("T");
	i = Bitpack64_read_one(i-sarray->indexT,sarray->csaTptrs,sarray->csaTcomp);
      } else if (i >= sarray->indexG) {
	assert(matchbits == 3);
	printf("G");
	i = Bitpack64_read_one(i-sarray->indexG,sarray->csaGptrs,sarray->csaGcomp);
      } else if (i >= sarray->indexC) {
	assert(matchbits == 1);
	printf("C");
	i = Bitpack64_read_one(i-sarray->indexC,sarray->csaCptrs,sarray->csaCcomp);
      } else {
	assert(matchbits == 0);
	printf("A");
	i = Bitpack64_read_one(i-sarray->indexA,sarray->csaAptrs,sarray->csaAcomp);
      }
#endif

      debug3b(assert(i == expected_i));
      nhops += 1;
    }

    debug3(printf("\n"));
    debug3(printf("Returning %u = %u - nhops %u\n",
		   sarray->array_samples[i/sarray->sa_sampling] - nhops,
		   sarray->array_samples[i/sarray->sa_sampling],nhops));
    
    debug3a(assert(sarray->array_samples[i/sarray->sa_sampling] - nhops == expected_sa_i));

    return sarray->array_samples[i/sarray->sa_sampling] - nhops;
  }
}

#elif defined(WORDS_BIGENDIAN)

#define csa_lookup(sarray,i) Bigendian_convert_uint(sarray->array[i])

#else

#define csa_lookup(sarray,i) sarray->array[i]

#endif


T
Sarray_new (char *dir, char *fileroot, Access_mode_T sarray_access, Access_mode_T lcp_access,
	    Access_mode_T guideexc_access, Access_mode_T indexij_access, bool sharedp, Mode_T mode, bool fwdp) {
  T new;
  char *comma1;
  double seconds;
  int npages;

  bool old_format_p;
  char *sarrayfile;		/* Old format */

  char *lcpchilddcfile;
  char *lcp_guidefile, *lcp_exceptionsfile;
  char *child_guidefile, *child_exceptionsfile;
#ifdef DEBUG15
  char *indexi_ptrsfile, *indexi_compfile;
  char *indexj_ptrsfile, *indexj_compfile;
  char *indexij_ptrsfile, *indexij_compfile;
#elif defined(USE_SEPARATE_BUCKETS)
  char *indexi_ptrsfile, *indexi_compfile;
  char *indexj_ptrsfile, *indexj_compfile;
#else
  char *indexij_ptrsfile, *indexij_compfile;
#endif

  char *mode_prefix;

  if (mode == STANDARD) {
    mode_prefix = ".";
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".metct.";
    } else {
      mode_prefix = ".metga.";
    }
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2iag.";
    } else {
      mode_prefix = ".a2itc.";
    }
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2itc.";
    } else {
      mode_prefix = ".a2iag.";
    }
  }

  /* Old format */
  sarrayfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sarray")+1,sizeof(char));
  sprintf(sarrayfile,"%s/%s%ssarray",dir,fileroot,mode_prefix);

  lcpchilddcfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpchilddc")+1,sizeof(char));
  sprintf(lcpchilddcfile,"%s/%s%ssalcpchilddc",dir,fileroot,mode_prefix);

  lcp_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpguide1024")+1,sizeof(char));
  sprintf(lcp_guidefile,"%s/%s%ssalcpguide1024",dir,fileroot,mode_prefix);
  lcp_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpexc")+1,sizeof(char));
  sprintf(lcp_exceptionsfile,"%s/%s%ssalcpexc",dir,fileroot,mode_prefix);

  child_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildguide1024")+1,sizeof(char));
  sprintf(child_guidefile,"%s/%s%ssachildguide1024",dir,fileroot,mode_prefix);
  child_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildexc")+1,sizeof(char));
  sprintf(child_exceptionsfile,"%s/%s%ssachildexc",dir,fileroot,mode_prefix);

#ifdef DEBUG15
  indexi_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64meta")+1,sizeof(char));
  sprintf(indexi_ptrsfile,"%s/%s%ssaindexi64meta",dir,fileroot,mode_prefix);
  indexi_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64strm")+1,sizeof(char));
  sprintf(indexi_compfile,"%s/%s%ssaindexi64strm",dir,fileroot,mode_prefix);
  indexj_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64meta")+1,sizeof(char));
  sprintf(indexj_ptrsfile,"%s/%s%ssaindexj64meta",dir,fileroot,mode_prefix);
  indexj_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64strm")+1,sizeof(char));
  sprintf(indexj_compfile,"%s/%s%ssaindexj64strm",dir,fileroot,mode_prefix);
  indexij_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindex64meta")+1,sizeof(char));
  sprintf(indexij_ptrsfile,"%s/%s%ssaindex64meta",dir,fileroot,mode_prefix);
  indexij_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindex64strm")+1,sizeof(char));
  sprintf(indexij_compfile,"%s/%s%ssaindex64strm",dir,fileroot,mode_prefix);
#elif defined(USE_SEPARATE_BUCKETS)
  indexi_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64meta")+1,sizeof(char));
  sprintf(indexi_ptrsfile,"%s/%s%ssaindexi64meta",dir,fileroot,mode_prefix);
  indexi_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64strm")+1,sizeof(char));
  sprintf(indexi_compfile,"%s/%s%ssaindexi64strm",dir,fileroot,mode_prefix);
  indexj_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64meta")+1,sizeof(char));
  sprintf(indexj_ptrsfile,"%s/%s%ssaindexj64meta",dir,fileroot,mode_prefix);
  indexj_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64strm")+1,sizeof(char));
  sprintf(indexj_compfile,"%s/%s%ssaindexj64strm",dir,fileroot,mode_prefix);
#else
  indexij_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64meta")+1,sizeof(char));
  sprintf(indexij_ptrsfile,"%s/%s%ssaindex64meta",dir,fileroot,mode_prefix);
  indexij_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64strm")+1,sizeof(char));
  sprintf(indexij_compfile,"%s/%s%ssaindex64strm",dir,fileroot,mode_prefix);
#endif

  if (Access_file_exists_p(sarrayfile) == false) {
    fprintf(stderr,"No suffix array for genome\n");
    new = (T) NULL;

  } else if (Access_file_exists_p(lcpchilddcfile) == false) {
    fprintf(stderr,"Enhanced suffix array file %s does not exist.  The genome was built using an obsolete version\n",
	    lcpchilddcfile);
    new = (T) NULL;
    exit(9);

  } else {
    new = (T) MALLOC_KEEP(sizeof(*new));
    old_format_p = true;

    if (sarray_access == USE_MMAP_PRELOAD) {
      if (old_format_p == true) {
	fprintf(stderr,"Pre-loading suffix array...");
	new->array = (UINT4 *) Access_mmap_and_preload(&new->array_fd,&new->array_len,&npages,&seconds,sarrayfile,
						       sizeof(UINT4));
	new->n_plus_one = new->array_len/sizeof(UINT4); /* Should be genomiclength + 1*/
	new->n = new->n_plus_one - 1;

	comma1 = Genomicpos_commafmt(new->array_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      }
      new->array_access = MMAPPED;

    } else if (sarray_access == USE_MMAP_ONLY) {
      if (old_format_p == true) {
	new->array = (UINT4 *) Access_mmap(&new->array_fd,&new->array_len,sarrayfile,/*randomp*/true);
	new->n_plus_one = new->array_len/sizeof(UINT4); /* Should be genomiclength + 1*/
	new->n = new->n_plus_one - 1;
      }
      new->array_access = MMAPPED;

    } else if (sarray_access == USE_ALLOCATE) {
      if (old_format_p == true) {
	fprintf(stderr,"Allocating memory for suffix array...");
	if (sharedp == true) {
	  new->array = (UINT4 *) Access_allocate_shared(&new->array_access,&new->array_shmid,&new->array_key,
							&new->array_fd,&new->array_len,&seconds,sarrayfile,sizeof(UINT4));
	} else {
	  new->array = (UINT4 *) Access_allocate_private(&new->array_access,&new->array_len,&seconds,sarrayfile,sizeof(UINT4));
	}
	new->n_plus_one = new->array_len/sizeof(UINT4); /* Should be genomiclength + 1*/
	new->n = new->n_plus_one - 1;
	comma1 = Genomicpos_commafmt(new->array_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      }
    }

#ifdef DEBUG15
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    new->indexi_ptrs = (UINT4 *) Access_allocate_private(&new->indexi_ptrs_len,&seconds,indexi_ptrsfile,sizeof(UINT4));
    new->indexi_comp = (UINT4 *) Access_allocate_private(&new->indexi_comp_len,&seconds,indexi_compfile,sizeof(UINT4));
    new->indexj_ptrs = (UINT4 *) Access_allocate_private(&new->indexj_ptrs_len,&seconds,indexj_ptrsfile,sizeof(UINT4));
    new->indexj_comp = (UINT4 *) Access_allocate_private(&new->indexj_comp_len,&seconds,indexj_compfile,sizeof(UINT4));
    new->indexij_ptrs = (UINT4 *) Access_allocate_private(&new->indexij_ptrs_len,&seconds,indexij_ptrsfile,sizeof(UINT4));
    new->indexij_comp = (UINT4 *) Access_allocate_private(&new->indexij_comp_len,&seconds,indexij_compfile,sizeof(UINT4));
    new->indexsize = 3 + log4(((new->indexij_ptrs_len - 8)/sizeof(UINT4)/2)/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#elif defined(USE_SEPARATE_BUCKETS)
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    new->indexi_ptrs = (UINT4 *) Access_allocate_private(&new->indexi_ptrs_access,&new->indexi_ptrs_len,&seconds,indexi_ptrsfile,sizeof(UINT4));
    new->indexi_comp = (UINT4 *) Access_allocate_private(&new->indexi_comp_access,&new->indexi_comp_len,&seconds,indexi_compfile,sizeof(UINT4));
    new->indexj_ptrs = (UINT4 *) Access_allocate_private(&new->indexj_ptrs_access,&new->indexj_ptrs_len,&seconds,indexj_ptrsfile,sizeof(UINT4));
    new->indexj_comp = (UINT4 *) Access_allocate_private(&new->indexj_comp_access,&new->indexj_comp_len,&seconds,indexj_compfile,sizeof(UINT4));
    new->indexsize = 3 + log4(((new->indexi_ptrs_len - 8)/sizeof(UINT4))/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#else
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    if (indexij_access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading indexij ptrs...");
      new->indexij_ptrs = (UINT4 *) Access_mmap_and_preload(&new->indexij_ptrs_fd,&new->indexij_ptrs_len,&npages,&seconds,indexij_ptrsfile,
							    sizeof(UINT4));
      comma1 = Genomicpos_commafmt(new->indexij_ptrs_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
      FREE(comma1);

      fprintf(stderr,"Pre-loading indexij comp...");
      new->indexij_comp = (UINT4 *) Access_mmap_and_preload(&new->indexij_comp_fd,&new->indexij_comp_len,&npages,&seconds,indexij_compfile,
							    sizeof(UINT4));
      comma1 = Genomicpos_commafmt(new->indexij_comp_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
      FREE(comma1);

      new->indexij_ptrs_access = MMAPPED;
      new->indexij_comp_access = MMAPPED;

    } else if (indexij_access == USE_MMAP_ONLY) {
      new->indexij_ptrs = (UINT4 *) Access_mmap(&new->indexij_ptrs_fd,&new->indexij_ptrs_len,indexij_ptrsfile,/*randomp*/true);
      new->indexij_comp = (UINT4 *) Access_mmap(&new->indexij_comp_fd,&new->indexij_comp_len,indexij_compfile,/*randomp*/true);

      new->indexij_ptrs_access = MMAPPED;
      new->indexij_comp_access = MMAPPED;

    } else if (indexij_access == USE_ALLOCATE) {
      if (sharedp == true) {
	fprintf(stderr,"Allocating memory for indexij ptrs...");
	new->indexij_ptrs = (UINT4 *) Access_allocate_shared(&new->indexij_ptrs_access,&new->indexij_ptrs_shmid,&new->indexij_ptrs_key,
							     &new->indexij_ptrs_fd,&new->indexij_ptrs_len,&seconds,indexij_ptrsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->indexij_ptrs_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);

	fprintf(stderr,"Allocating memory for indexij comp...");
	new->indexij_comp = (UINT4 *) Access_allocate_shared(&new->indexij_comp_access,&new->indexij_comp_shmid,&new->indexij_comp_key,
							     &new->indexij_comp_fd,&new->indexij_comp_len,&seconds,indexij_compfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->indexij_comp_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      } else {
	fprintf(stderr,"Allocating memory for indexij ptrs...");
	new->indexij_ptrs = (UINT4 *) Access_allocate_private(&new->indexij_ptrs_access,&new->indexij_ptrs_len,&seconds,indexij_ptrsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->indexij_ptrs_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);

	fprintf(stderr,"Allocating memory for indexij comp...");
	new->indexij_comp = (UINT4 *) Access_allocate_private(&new->indexij_comp_access,&new->indexij_comp_len,&seconds,indexij_compfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->indexij_comp_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      }

    }

    new->indexsize = 3 + log4(((new->indexij_ptrs_len - 8)/sizeof(UINT4)/2)/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#endif
    new->indexspace = power(4,new->indexsize);

    if (lcp_access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading LCP/child/DC arrays...");
      new->lcpchilddc = (unsigned char *) Access_mmap_and_preload(&new->lcpchilddc_fd,&new->lcpchilddc_len,&npages,&seconds,
								  lcpchilddcfile,sizeof(unsigned char));
      new->lcpchilddc_access = MMAPPED;
      comma1 = Genomicpos_commafmt(new->lcpchilddc_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
      FREE(comma1);
    } else if (lcp_access == USE_MMAP_ONLY) {
      new->lcpchilddc = (unsigned char *) Access_mmap(&new->lcpchilddc_fd,&new->lcpchilddc_len,lcpchilddcfile,/*randomp*/true);
      new->lcpchilddc_access = MMAPPED;
    } else if (lcp_access == USE_ALLOCATE) {
      fprintf(stderr,"Allocating memory for lcpchildc...");
      if (sharedp == true) {
	new->lcpchilddc = (unsigned char *) Access_allocate_shared(&new->lcpchilddc_access,&new->lcpchilddc_shmid,&new->lcpchilddc_key,
								   &new->lcpchilddc_fd,&new->lcpchilddc_len,&seconds,lcpchilddcfile,sizeof(unsigned char));
      } else {
	new->lcpchilddc = (unsigned char *) Access_allocate_private(&new->lcpchilddc_access,&new->lcpchilddc_len,&seconds,lcpchilddcfile,sizeof(unsigned char));
      }
      comma1 = Genomicpos_commafmt(new->lcpchilddc_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
      FREE(comma1);
    }

    if (guideexc_access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading guide/exceptions...");
      new->lcp_guide = (UINT4 *) Access_mmap_and_preload(&new->lcp_guide_fd,&new->lcp_guide_len,&npages,&seconds,
							 lcp_guidefile,sizeof(UINT4));
      new->lcp_exceptions = (UINT4 *) Access_mmap_and_preload(&new->lcp_exceptions_fd,&new->lcp_exceptions_len,&npages,&seconds,
							 lcp_exceptionsfile,sizeof(UINT4));
      new->child_guide = (UINT4 *) Access_mmap_and_preload(&new->child_guide_fd,&new->child_guide_len,&npages,&seconds,
							 child_guidefile,sizeof(UINT4));
      new->child_exceptions = (UINT4 *) Access_mmap_and_preload(&new->child_exceptions_fd,&new->child_exceptions_len,&npages,&seconds,
							 child_exceptionsfile,sizeof(UINT4));
      new->lcp_guide_access = MMAPPED;
      new->lcp_exceptions_access = MMAPPED;
      new->child_guide_access = MMAPPED;
      new->child_exceptions_access = MMAPPED;
      fprintf(stderr,"done\n");

    } else if (guideexc_access == USE_MMAP_ONLY) {
      new->lcp_guide = (UINT4 *) Access_mmap(&new->lcp_guide_fd,&new->lcp_guide_len,
					     lcp_guidefile,/*randomp*/true);
      new->lcp_exceptions = (UINT4 *) Access_mmap(&new->lcp_exceptions_fd,&new->lcp_exceptions_len,
						  lcp_exceptionsfile,/*randomp*/true);
      new->child_guide = (UINT4 *) Access_mmap(&new->child_guide_fd,&new->child_guide_len,
					       child_guidefile,/*randomp*/true);
      new->child_exceptions = (UINT4 *) Access_mmap(&new->child_exceptions_fd,&new->child_exceptions_len,
							 child_exceptionsfile,/*randomp*/true);
      new->lcp_guide_access = MMAPPED;
      new->lcp_exceptions_access = MMAPPED;
      new->child_guide_access = MMAPPED;
      new->child_exceptions_access = MMAPPED;

    } else if (guideexc_access == USE_ALLOCATE) {
      fprintf(stderr,"Allocating memory for lcp guide...");
      if (sharedp == true) {
	new->lcp_guide = (UINT4 *) Access_allocate_shared(&new->lcp_guide_access,&new->lcp_guide_shmid,&new->lcp_guide_key,
							  &new->lcp_guide_fd,&new->lcp_guide_len,&seconds,lcp_guidefile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->lcp_guide_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      
	fprintf(stderr,"Allocating memory for lcp exceptions...");
	new->lcp_exceptions = (UINT4 *) Access_allocate_shared(&new->lcp_exceptions_access,&new->lcp_exceptions_shmid,&new->lcp_exceptions_key,
							       &new->lcp_exceptions_fd,&new->lcp_exceptions_len,&seconds,lcp_exceptionsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->lcp_exceptions_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
	
	fprintf(stderr,"Allocating memory for child guide...");
	new->child_guide = (UINT4 *) Access_allocate_shared(&new->child_guide_access,&new->child_guide_shmid,&new->child_guide_key,
							    &new->child_guide_fd,&new->child_guide_len,&seconds,child_guidefile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->child_guide_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);

	fprintf(stderr,"Allocating memory for child exceptions...");
	new->child_exceptions = (UINT4 *) Access_allocate_shared(&new->child_exceptions_access,&new->child_exceptions_shmid,&new->child_exceptions_key,
								 &new->child_exceptions_fd,&new->child_exceptions_len,&seconds,child_exceptionsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->child_exceptions_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      
      } else {
	new->lcp_guide = (UINT4 *) Access_allocate_private(&new->lcp_guide_access,&new->lcp_guide_len,&seconds,lcp_guidefile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->lcp_guide_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      
	fprintf(stderr,"Allocating memory for lcp exceptions...");
	new->lcp_exceptions = (UINT4 *) Access_allocate_private(&new->lcp_exceptions_access,&new->lcp_exceptions_len,&seconds,
								lcp_exceptionsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->lcp_exceptions_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
	
	fprintf(stderr,"Allocating memory for child guide...");
	new->child_guide = (UINT4 *) Access_allocate_private(&new->child_guide_access,&new->child_guide_len,&seconds,child_guidefile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->child_guide_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);

	fprintf(stderr,"Allocating memory for child exceptions...");
	new->child_exceptions = (UINT4 *) Access_allocate_private(&new->child_exceptions_access,&new->child_exceptions_len,&seconds,
								  child_exceptionsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->child_exceptions_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      }
    }

    new->n_lcp_exceptions = new->lcp_exceptions_len/(sizeof(UINT4) + sizeof(UINT4));
    new->child_guide_interval = 1024;
  }


  FREE(child_exceptionsfile);
  FREE(child_guidefile);

  FREE(lcp_exceptionsfile);
  FREE(lcp_guidefile);

  FREE(lcpchilddcfile);

#ifdef DEBUG15
  FREE(indexi_compfile);
  FREE(indexi_ptrsfile);
  FREE(indexj_compfile);
  FREE(indexj_ptrsfile);
  FREE(indexij_compfile);
  FREE(indexij_ptrsfile);
#elif defined(USE_SEPARATE_BUCKETS)
  FREE(indexi_compfile);
  FREE(indexi_ptrsfile);
  FREE(indexj_compfile);
  FREE(indexj_ptrsfile);
#else
  FREE(indexij_compfile);
  FREE(indexij_ptrsfile);
#endif

  FREE(sarrayfile);

  return new;
}


void
Sarray_free (T *old) {
  if (*old) {
#ifdef DEBUG15
    FREE((*old)->indexi_ptrs);
    FREE((*old)->indexi_comp);
    FREE((*old)->indexj_ptrs);
    FREE((*old)->indexj_comp);
    FREE((*old)->indexij_ptrs);
    FREE((*old)->indexij_comp);
#elif defined(USE_SEPARATE_BUCKETS)
    FREE((*old)->indexi_ptrs);
    FREE((*old)->indexi_comp);
    FREE((*old)->indexj_ptrs);
    FREE((*old)->indexj_comp);
#else
    if ((*old)->indexij_ptrs_access == MMAPPED) {
      munmap((void *) (*old)->indexij_ptrs,(*old)->indexij_ptrs_len);
      close((*old)->indexij_ptrs_fd);
    } else if ((*old)->indexij_ptrs_access == ALLOCATED_PRIVATE) {
      FREE((*old)->indexij_ptrs);
    } else if ((*old)->indexij_ptrs_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->indexij_ptrs,(*old)->indexij_ptrs_shmid,(*old)->indexij_ptrs_key);
    }
    if ((*old)->indexij_comp_access == MMAPPED) {
      munmap((void *) (*old)->indexij_comp,(*old)->indexij_comp_len);
      close((*old)->indexij_comp_fd);
    } else if ((*old)->indexij_comp_access == ALLOCATED_PRIVATE) {
      FREE((*old)->indexij_comp);
    } else if ((*old)->indexij_comp_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->indexij_comp,(*old)->indexij_comp_shmid,(*old)->indexij_comp_key);
    }
#endif

    if ((*old)->lcp_guide_access == MMAPPED) {
      munmap((void *) (*old)->lcp_guide,(*old)->lcp_guide_len);
      close((*old)->lcp_guide_fd);
    } else if ((*old)->lcp_guide_access == ALLOCATED_PRIVATE) {
      FREE((*old)->lcp_guide);
    } else if ((*old)->lcp_guide_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->lcp_guide,(*old)->lcp_guide_shmid,(*old)->lcp_guide_key);
    }

    if ((*old)->lcp_exceptions_access == MMAPPED) {
      munmap((void *) (*old)->lcp_exceptions,(*old)->lcp_exceptions_len);
      close((*old)->lcp_exceptions_fd);
    } else if ((*old)->lcp_exceptions_access == ALLOCATED_PRIVATE) {
      FREE((*old)->lcp_exceptions);
    } else if ((*old)->lcp_exceptions_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->lcp_exceptions,(*old)->lcp_exceptions_shmid,(*old)->lcp_exceptions_key);
    }

    if ((*old)->child_guide_access == MMAPPED) {
      munmap((void *) (*old)->child_guide,(*old)->child_guide_len);
      close((*old)->child_guide_fd);
    } else if ((*old)->child_guide_access == ALLOCATED_PRIVATE) {
      FREE((*old)->child_guide);
    } else if ((*old)->child_guide_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->child_guide,(*old)->child_guide_shmid,(*old)->child_guide_key);
    }

    if ((*old)->child_exceptions_access == MMAPPED) {
      munmap((void *) (*old)->child_exceptions,(*old)->child_exceptions_len);
      close((*old)->child_exceptions_fd);
    } else if ((*old)->child_exceptions_access == ALLOCATED_PRIVATE) {
      FREE((*old)->child_exceptions);
    } else if ((*old)->child_exceptions_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->child_exceptions,(*old)->child_exceptions_shmid,(*old)->child_exceptions_key);
    }

    if ((*old)->lcpchilddc_access == MMAPPED) {
      munmap((void *) (*old)->lcpchilddc,(*old)->lcpchilddc_len);
      close((*old)->lcpchilddc_fd);
    } else if ((*old)->lcpchilddc_access == ALLOCATED_PRIVATE) {
      FREE((*old)->lcpchilddc);
    } else if ((*old)->lcpchilddc_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->lcpchilddc,(*old)->lcpchilddc_shmid,(*old)->lcpchilddc_key);
    }

    if ((*old)->array_access == MMAPPED) {
      munmap((void *) (*old)->array,(*old)->array_len);
      close((*old)->array_fd);
    } else if ((*old)->array_access == ALLOCATED_PRIVATE) {
      FREE((*old)->array);
    } else if ((*old)->array_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->array,(*old)->array_shmid,(*old)->array_key);
    }

    FREE_KEEP(*old);
  }

  return;
}



#if 0
/* Old search method.  O(m*(log n)), where m is the querylength and n
   is the size of the suffix array searched */
static Sarrayptr_T
sarray_search_init (char *query, int querylength, int queryoffset, Compress_T query_compress, bool plusp,
		    Sarrayptr_T low, Sarrayptr_T high, Univcoord_T nmatches_low, Univcoord_T nmatches_high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  Univcoord_T nmatches_mid, fasti;
  char c;
  UINT4 sa_low, sa_mid;
  UINT4 lcp_low, lcp_mid;

  assert(querylength > 0);

  debug1(printf("sarray_search_init on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
#if 0
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif

    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

#ifdef WORDS_BIGENDIAN
    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/Bigendian_convert_uint(sarray->array[mid])-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand);
    pos = Bigendian_convert_uint(sarray->array[mid]) + fasti;
#else
    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand);
    pos = sarray->array[mid] + fasti;
#endif
    c = Genome_get_char_lex(genome,pos,sarray->n,chartable);

    if (fasti == (Univcoord_T) querylength || c > query[fasti]) {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
#ifdef WORDS_BIGENDIAN
      sa_mid = Bigendian_convert_uint(sarray->array[mid]);
#else
      sa_mid = sarray->array[mid];
#endif
      lcp_mid = Bitpack64_read_one(sa_mid,sarray->plcp_ptrs,sarray->plcp_comp) - sa_mid;
#ifdef USE_LCP
      if (lcp_mid != sarray->lcp[mid]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
    } else {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
#ifdef WORDS_BIGENDIAN
      sa_low = Bigendian_convert_uint(sarray->array[low]);
#else
      sa_low = sarray->array[low];
#endif
      lcp_low = Bitpack64_read_one(sa_low,sarray->plcp_ptrs,sarray->plcp_comp) - sa_low;
#ifdef USE_LCP
      if (lcp_low != sarray->lcp[low]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
    }

    debug1(printf("sarray_search_init with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_init ended.  Returning low %u+1\n\n",low));
  return low + 1;
}
#endif


#if 0
/* Old search method.  O(m*(log n)), where m is the querylength and n
   is the size of the suffix array searched */
static Sarrayptr_T
sarray_search_final (char *query, int querylength, int queryoffset, Compress_T query_compress, bool plusp,
		     Sarrayptr_T low, Sarrayptr_T high, Univcoord_T nmatches_low, Univcoord_T nmatches_high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  Univcoord_T nmatches_mid, fasti;
  UINT4 sa_low, sa_mid;
  UINT4 lcp_low, lcp_mid;
  char c;

  assert(querylength > 0);

  debug1(printf("sarray_search_final on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
#if 0
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif
    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

#ifdef WORDS_BIGENDIAN
    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/Bigendian_convert_uint(sarray->array[mid])-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand);
    pos = Bigendian_convert_uint(sarray->array[mid]) + fasti;
#else
    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand);
    pos = sarray->array[mid] + fasti;
#endif
    c = Genome_get_char_lex(genome,pos,sarray->n,chartable);

    if (fasti == (Univcoord_T) querylength || c < query[fasti]) {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
#ifdef WORDS_BIGENDIAN
      sa_low = Bigendian_convert_uint(sarray->array[low]);
#else
      sa_low = sarray->array[low];
#endif
      lcp_low = Bitpack64_read_one(sa_low,sarray->plcp_ptrs,sarray->plcp_comp) - sa_low;
#ifdef USE_LCP
      if (lcp_low != sarray->lcp[low]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
    } else {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
#ifdef WORDS_BIGENDIAN
      sa_mid = Bigendian_convert_uint(sarray->array[mid]);
#else
      sa_mid = sarray->array[mid];
#endif
      lcp_mid = Bitpack64_read_one(sa_mid,sarray->plcp_ptrs,sarray->plcp_comp) - sa_mid;
#ifdef USE_LCP
      if (lcp_mid != sarray->lcp[mid]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
    }

    debug1(printf("sarray_search_final with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_final ended.  Returning high %u-1\n\n",high-1));
  return high - 1;
}
#endif


int
nt_querylength (char *query, int querylength) {
  int i;
  char c;

  i = 0;
  while (i < querylength && ((c = query[i]) == 'A' || c == 'C' || c == 'G' || c == 'T')) {
    i++;
  }

  return i;
}


Oligospace_T
nt_oligo (char *query, int indexsize) {
  Oligospace_T oligo = 0U;
  int i;

  for (i = 0; i < indexsize; i++) {
    oligo *= 4;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default:
      fprintf(stderr,"Saw N in nt_oligo\n");
      abort();
    }
  }

  return oligo;
}

Oligospace_T
nt_oligo_truncate (char *query, int truncsize, int indexsize, int subst_value) {
  Oligospace_T oligo = 0U;
  int i;

  for (i = 0; i < truncsize; i++) {
    oligo *= 4;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default:
      fprintf(stderr,"Saw N in nt_oligo\n");
      abort();
    }
  }

  for ( ; i < indexsize; i++) {
    oligo *= 4;
    oligo += subst_value;
  }

  return oligo;
}



/* For child[index+1].up, just calling child[index] */
#define decode_up(index,child_bytes,child_guide,child_exceptions,child_guide_interval) index - Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval)
#define decode_down(index,child_bytes,child_guide,child_exceptions,child_guide_interval) Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval) + index + 1
#define decode_next(index,child_bytes,child_guide,child_exceptions,child_guide_interval) Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval) + index + 1

#if 0
/*                                      0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F */
static char discrim_char_before[16] = {'?','$','$','$','$','$','A','A','A','A','C','C','C','G','G','T'};
static char discrim_char_after[16]  = {'?','A','C','G','T','X','C','G','T','X','G','T','X','T','X','X'};
#endif

static bool
get_child_given_first (Sarrayptr_T *l, Sarrayptr_T *r, Sarrayptr_T i, Sarrayptr_T j, char desired_char,
		       T sarray, unsigned char *lcpchilddc, UINT4 lcp_whole, UINT4 nextl) {
  char c1, c2;
  UINT4 child_next;

  debug2(printf("Getting children for l-interval from %u to %u, char %c\n",i,j,desired_char));

#if 0
  /* First child already given */
  debug1(printf("lcp-interval %u..%u\n",i,j));
  up = decode_up(j,sarray->child_bytes,sarray->child_guide,sarray->child_exceptions,sarray->child_guide_interval);
  if (i < up && up <= j) {
    nextl = up;
    debug2(printf("nextl is up: %u\n",nextl));
  } else {
    nextl = decode_down(i,sarray->child_bytes,sarray->child_guide,sarray->child_exceptions,sarray->child_guide_interval); /* down */
    debug2(printf("nextl is down: %u\n",nextl));
  }
#endif

  /* Test first child: Use discrim_chars, rather than looking up S[SA[i] + lcp_whole] */
  c2 = Bytecoding_lcpchilddc_dc(&c1,nextl,lcpchilddc);
  debug2(printf("First child: %u to %u, discrim chars %c and %c\n",i,nextl-1,c1,c2));

  if (desired_char < c1) {
    debug2(printf("1.  Returning false, because desired %c < c1 %c\n",desired_char,c1));
    return false;
  } else if (desired_char == c1) {
    *l = i;
    *r = nextl - 1;
    debug2(printf("Returning true\n\n"));
    return true;
  } else if (desired_char < c2) {
    debug2(printf("1.  Returning false, because desired %c < c2 %c\n",desired_char,c2));
    return false;
  } else {
    /* Advance to middle children or final child */
    debug2(printf("1.  Advancing\n"));
  }

  /* Test for child[i] being down: lcp[child[i]] > lcp[i] */
  /* Test for child[i] being next_lindex: lcp[child[i]] == lcp[i] */
  /* Test middle children */
  while (nextl < j && Bytecoding_lcpchilddc_lcp_next(&child_next,nextl,/*bytes*/lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						     sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions) == lcp_whole) {
    /* Already tested for desired_char < c2 */
    if (desired_char == c2) {
      *l = nextl;
#if 0
      *r = Bytecoding_lcpchilddc_child_next(nextl,lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					    sarray->child_guide_interval) - 1; /* child[nextl] - 1 */
#else
      *r = child_next - 1;
#endif
      debug2(printf("Child: %u to %u, c2 %c\n",nextl,*r,c2));
      debug2(printf("Returning true\n\n"));
      return true;
    } else {
      debug2(printf("Child: %u",nextl));
#if 0
      nextl = Bytecoding_lcpchilddc_child_next(nextl,lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					       sarray->child_guide_interval); /* child[nextl] */
#else
      nextl = child_next;
#endif
      c2 = Bytecoding_lcpchilddc_dc(&c1,nextl,lcpchilddc);
      debug2(printf(" to %u, discrim chars %c and %c\n",nextl-1,c1,c2));

      if (desired_char < c2) {
	debug2(printf("M.  Returning false, because desired %c < c2 %c\n",desired_char,c2));
	return false;
      } else {
	debug2(printf("M.  Advancing\n"));
      }
    }
  }

  /* Test last child */
  /* Already tested for desired_char < c2 */
  debug2(printf("Final child: %u to %u, c2 %c\n",nextl,j,c2));
  if (desired_char == c2) {
    *l = nextl;
    *r = j;
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    debug2(printf("3.  Returning false, because desired %c != c2 %c\n",desired_char,c2));
    return false;
  }
}


static UINT4
find_longest_match (UINT4 nmatches, Sarrayptr_T *initptr, Sarrayptr_T *finalptr,
		    Sarrayptr_T i, Sarrayptr_T j, char *query, UINT4 querylength,
		    int queryoffset, Compress_T query_compress, T sarray, bool plusp,
		    int genestrand, char conversion[]) {
  UINT4 lcp_whole, nextl, up;
  UINT4 minlength;
  UINT4 l, r;
  Univcoord_T SA_i;

  while (nmatches < querylength) {
    if (i == j) {
      /* Singleton interval */
      debug1(printf("Singleton interval %u..%u\n",i,j));
      SA_i = csa_lookup(sarray,i);
      nmatches +=
	Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
					     /*pos5*/queryoffset+nmatches,/*pos3*/queryoffset+querylength,
					     plusp,genestrand);
      *initptr = i;
      *finalptr = j;
      return nmatches;

    } else {
      /* First child */
      debug1(printf("lcp-interval %u..%u\n",i,j));
      up = Bytecoding_lcpchilddc_child_up(j,sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					  sarray->child_guide_interval);
      if (i < up && up <= j) {
	nextl = up;
	debug2(printf("nextl is up: %u\n",nextl));
      } else {
	nextl = Bytecoding_lcpchilddc_child_next(i,sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval); /* really down */
	debug2(printf("nextl is down: %u\n",nextl));
      }

      lcp_whole = Bytecoding_lcpchilddc_lcp(nextl,sarray->lcpchilddc,sarray->lcp_exceptions,
					    sarray->n_lcp_exceptions); /* lcp(i,j) */
      debug1(printf("lcp_whole for %u..%u is %d, compared with nmatches %d\n",i,j,lcp_whole,nmatches));

      if (lcp_whole > nmatches) {
	/* Check only up to minlength, so we validate the entire interval */
	minlength = (lcp_whole < querylength) ? lcp_whole : querylength;
	debug1(printf("Looking up genome for query from %d .. %d - 1\n",nmatches,minlength));
	SA_i = csa_lookup(sarray,i);
	nmatches +=
	  Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
					       /*pos5*/queryoffset+nmatches,/*pos3*/queryoffset+minlength,
					       plusp,genestrand);
	if (nmatches < minlength) {
	  *initptr = i;
	  *finalptr = j;
	  return nmatches;

	} else if (nmatches >= querylength) {
	  debug1(printf("nmatches is now %d >= querylength %d => success\n",nmatches,querylength));
	  *initptr = i;
	  *finalptr = j;
	  return nmatches;
	}
      }
	
      debug1(printf("nmatches is now %d => desired_char is %c => %c\n",
		    nmatches,query[nmatches],conversion[query[nmatches]]));
      if (get_child_given_first(&l,&r,i,j,/*desired_char*/conversion[(int) query[nmatches]],
				sarray,sarray->lcpchilddc,lcp_whole,nextl) == false) {
	*initptr = i;
	*finalptr = j;
	return nmatches;
      } else {
	nmatches += 1;
	i = l;
	j = r;
      }
    }
  }

  *initptr = i;
  *finalptr = j;
  return nmatches;
}



/* Searches using LCP and child arrays.  Should be O(m * |Sigma|),
   where m wis the querylength and |Sigma| is the size of the alphabet
   (4 for DNA) */
/* query is a substring of the original, starting with queryoffset */
static void
sarray_search (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, bool *successp,
	       UINT4 *nmatches, char *query, UINT4 querylength, int queryoffset,
	       Compress_T query_compress, T sarray, bool plusp, int genestrand,
	       char conversion[]) {
  int effective_querylength;	/* length to first N */
  Oligospace_T oligo;
  UINT4 l, r;

#ifdef DEBUG1
  Univcoord_T SA_i, hit, child_next;
  int k = 0;
  UINT4 recount, lcp_prev, lcp_next, lcp_i, max_lcp;
  char Buffer[1000+1], c1, c2;
  bool failp;
#endif

  debug1(printf("sarray_search on %.*s, querylength %d, plusp %d\n",querylength,query,querylength,plusp));

  /* Find initial lcp-interval */
  effective_querylength = nt_querylength(query,querylength);

  *nmatches = 0;
  if (effective_querylength == 0) {
    *initptr = *finalptr = 0;
    *successp = false;
    return;

  } else if (effective_querylength < sarray->indexsize) {
    debug1(printf("string %.*s with effective querylength %d is shorter than indexsize",
		  querylength,query,effective_querylength));
    l = 1;
    r = sarray->n;

  } else {
    oligo = nt_oligo(query,sarray->indexsize);
#ifdef DEBUG15
    if ((l = Bitpack64_read_two(&r,oligo*2,sarray->indexij_ptrs,sarray->indexij_comp)) !=
	Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp)) {
      abort();
    } else if (r - 1 != Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp)) {
      printf("For oligo %u, separate buckets give %u and %u, while single bucket gives %u and %u\n",
	     oligo,
	     Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp),
	     Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp),
	     l,r);
      abort();
    }
    r--;			/* Because interleaved writes r+1 to maintain monotonicity */
#elif defined(USE_SEPARATE_BUCKETS)
    l = Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp);
    r = Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp);
#else
    l = Bitpack64_read_two(&r,oligo*2,sarray->indexij_ptrs,sarray->indexij_comp);
    r--;			/* Because interleaved writes r+1 to maintain monotonicity */
#endif
    debug1(printf("string %.*s is equal/longer than indexsize %d => oligo %u => interval %u..%u",
		  querylength,query,sarray->indexsize,oligo,l,r));
    if (l <= r) {
      debug1(printf(" (good)\n"));
      *nmatches = sarray->indexsize;
      /* i = l; */
      /* j = r; */
    } else {
      /* The entire lcp-interval [1,sarray->n] should also work without initindex */
      l = 1;
      r = sarray->n;
      debug1(printf(" (bad) => entire lcp-interval: %u..%u\n",l,r));
    }
  }

  if (l > r) {
    /* Did not find a match using saindex or one letter */
    *initptr = l;
    *finalptr = r;
  } else {
    *nmatches = find_longest_match(*nmatches,&(*initptr),&(*finalptr),/*i*/l,/*j*/r,
				   query,querylength,queryoffset,query_compress,sarray,
				   plusp,genestrand,conversion);
  }

  /* Search through suffix tree */
  debug1(printf("initptr gets %u, finalptr gets %u\n",*initptr,*finalptr));

  if (*nmatches < querylength) {
    *successp = false;
    debug1(printf("%s fail at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  } else {
    *successp = true;
    debug1(printf("%s success at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  }

#ifdef DEBUG1
  failp = false;

  /* Before */
  if (*nmatches > 0 && *initptr > 0U) {
    SA_i = csa_lookup(sarray,(*initptr)-1);
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand);
    printf("%d\t%u\t%u\t",recount,(*initptr)-1,SA_i/*+ 1U*/);
    c2 = Bytecoding_lcpchilddc_dc(&c1,(*initptr)-1,sarray->lcpchilddc);
    printf("%c%c\t",c1,c2);
    lcp_i = Bytecoding_lcpchilddc_lcp((*initptr)-1,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    printf("%u\t",lcp_i);
    lcp_next = Bytecoding_lcpchilddc_lcp((*initptr),/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    printf("%u\t",Bytecoding_lcpchilddc_lcp_next(&child_next,(*initptr)-1,/*bytes*/sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(SA_i,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(SA_i,recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(SA_i,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(SA_i,recount+1,Buffer);
      }
    }
    printf("%s\n",Buffer);
    if (recount >= *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false negative: recount %d at %u before init does equal expected nmatches %d\n",
	     recount,SA_i,*nmatches);
      failp = true;
    }
  }
  printf("\n");


  /* Hits */
  lcp_prev = lcp_i;
  for (k = 0; k < (int) (*finalptr - *initptr + 1) && k < MAX_DEBUG1_HITS; k++) {
    SA_i = csa_lookup(sarray,(*initptr)+k);
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand);
    printf("%d\t%u\t%u\t",recount,(*initptr)+k,SA_i/*+ 1U*/);
    c2 = Bytecoding_lcpchilddc_dc(&c1,(*initptr)+k,sarray->lcpchilddc);
    printf("%c%c\t",c1,c2);
    lcp_i = Bytecoding_lcpchilddc_lcp((*initptr)+k,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    lcp_next = Bytecoding_lcpchilddc_lcp((*initptr)+k+1,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    printf("%u\t",lcp_i);
    printf("%u\t",Bytecoding_lcpchilddc_lcp_next(&child_next,(*initptr)+k,/*bytes*/sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    max_lcp = lcp_i;
    if (lcp_prev > max_lcp) {
      max_lcp = lcp_prev;
    }
    if (lcp_next > max_lcp) {
      max_lcp = lcp_next;
    }
    if (max_lcp > 1000) {
      max_lcp = 1000;
    }

    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(SA_i,max_lcp+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(SA_i,max_lcp+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(SA_i,max_lcp+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(SA_i,max_lcp+1,Buffer);
      }
    }
    printf("%s\n",Buffer);
    if (recount != *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false positive: recount %d at %u does not equal expected nmatches %d\n",
	     recount,csa_lookup(sarray,(*initptr)),*nmatches);
      failp = true;
    }

    lcp_prev = lcp_i;
  }

  if (k < (int) (*finalptr - *initptr + 1)) {
    /* Overflow */
    printf("...\n");
    k = (int) (*finalptr - *initptr);
    hit = csa_lookup(sarray,(*initptr)+k);
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/hit-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand);
    printf("%d\t%u\t%u\t",recount,(*initptr)+k,hit /*+ 1U*/);
    c2 = Bytecoding_lcpchilddc_dc(&c1,(*initptr)+k,sarray->lcpchilddc);
    printf("%c%c\t",c1,c2);
    lcp_i = Bytecoding_lcpchilddc_lcp((*initptr)+k,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    lcp_next = Bytecoding_lcpchilddc_lcp((*initptr)+k+1,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    printf("%u\t",lcp_i);
    printf("%u\t",Bytecoding_lcpchilddc_lcp_next(&child_next,(*initptr)+k,/*bytes*/sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(hit,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(hit,recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(hit,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(hit,recount+1,Buffer);
      }
    }
    printf("%s\n",Buffer);
    if (recount != *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false positive: recount %d at %u does not equal expected nmatches %d\n",
	     recount,csa_lookup(sarray,*initptr),*nmatches);
      failp = true;
    }
    /* hits[k] = sarray->array[(*initptr)++]; */
  }


  /* After */
  if (*nmatches > 0 && (SA_i = csa_lookup(sarray,(*finalptr)+1)) > 0U) {
    printf("\n");
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand);
    printf("%d\t%u\t%u\t",recount,(*finalptr)+1,SA_i/*+ 1U*/);
    c2 = Bytecoding_lcpchilddc_dc(&c1,(*finalptr)+1,sarray->lcpchilddc);
    printf("%c%c\t",c1,c2);
    printf("%u\t",Bytecoding_lcpchilddc_lcp((*finalptr)+1,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    printf("%u\t",Bytecoding_lcpchilddc_lcp_next(&child_next,(*finalptr)+1,/*bytes*/sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(SA_i,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(SA_i,recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(SA_i,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(SA_i,recount+1,Buffer);
      }
    }
    printf("%s\n",Buffer);
    if (recount >= *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false negative: recount %d at %u after (*finalptr) does equal expected nmatches %d\n",
	     recount,SA_i,*nmatches);
      failp = true;
    }
  }

  if (failp == true) {
    /* Can happen because $ ranks below 0 */
    /* Can also happen with CMET or ATOI, since genome128_hr procedures find genome-to-query mismatches */
    /* abort(); */
  }
#endif

  return;
}


/* For fill_positions_all: ELT_VIRGIN -> ELT_FILLED */
/* For fill_positions_filtered: ELT_VIRGIN -(1st call)-> ELT_UNSORTED -(2nd call)-> ELT_SORTED */
typedef enum {ELT_VIRGIN, ELT_FILLED, ELT_UNSORTED, ELT_SORTED} Elt_status_T;


/* Simplified version of Spanningelt_T */
typedef struct Elt_T *Elt_T;
struct Elt_T {
  int querystart;
  int queryend;

  int querystart_leftward; /* Modified when we extend matches leftward */
  int queryend_leftward; /* Modified when we extend matches leftward */

  int nmatches;

  Sarrayptr_T initptr;			/* in sarray */
  Sarrayptr_T finalptr;
  Sarrayptr_T nptr;

  Univcoord_T *positions_allocated; /* all or filtered positions needed */
  Univcoord_T *positions;
  int npositions_allocated;
  int npositions;		/* from goal to high */

  bool temporaryp;
  bool fillin_p;		/* Created by oligoindex algorithm */

  /* filled/sorted by Elt_fill_positions_filtered to speed up on multiple calls */
  Univcoord_T *all_positions;
  int n_all_positions;

  Elt_status_T status;
};


static void
Elt_reset (Elt_T this) {
  this->querystart_leftward = this->querystart;
  this->queryend_leftward = this->queryend;
  return;
}


static Elt_T
Elt_new (int querypos, int nmatches, Sarrayptr_T initptr, Sarrayptr_T finalptr, bool temporaryp) {
  Elt_T new = (Elt_T) MALLOC(sizeof(*new));

  new->querystart = new->querystart_leftward = querypos;
  new->queryend = new->queryend_leftward = querypos + nmatches - 1;
  new->nmatches = nmatches;

  new->initptr = initptr;
  new->finalptr = finalptr;
  new->nptr = new->finalptr - new->initptr + 1;

  /* new->positions is a pointer that advances to goal */
  new->positions_allocated = new->positions = (Univcoord_T *) NULL;
  new->npositions_allocated = new->npositions = 0;

  new->temporaryp = temporaryp;
  new->fillin_p = false;

  new->all_positions = (Univcoord_T *) NULL;
  new->n_all_positions = 0;

  new->status = ELT_VIRGIN;

  return new;
}

#if 0
static Elt_T
Elt_new_fillin (int querystart, int queryend, int indexsize, Univcoord_T left) {
  Elt_T new = (Elt_T) MALLOC(sizeof(*new));

  new->querystart = new->querystart_leftward = querystart;
  new->queryend = new->queryend_leftward = queryend + indexsize - 1;
  new->nmatches = new->queryend - querystart + 1;

  new->initptr = 0;
  new->finalptr = 0;
  new->nptr = 0;

  new->npositions = 1;
  new->positions_allocated = new->positions = (Univcoord_T *) MALLOC(sizeof(Univcoord_T));
  new->positions[0] = left;

  new->temporaryp = true;
  new->fillin_p = true;

  new->all_positions = (Univcoord_T *) NULL;
  new->n_all_positions = 0;

  new->status = ELT_VIRGIN;

  return new;
}
#endif

#if 0
static void
Elt_replace (Elt_T this, int querypos, int nmatches, Sarrayptr_T initptr, Sarrayptr_T finalptr) {
  this->querystart = querypos;
  this->queryend = querypos + nmatches - 1;
  this->nmatches = nmatches;

  this->initptr = initptr;
  this->finalptr = finalptr;

  if (this->positions_allocated != NULL) {
    FREE(this->positions_allocated);
  }
  this->positions_allocated = this->positions = (Univcoord_T *) NULL;
  this->npositions_allocated = this->npositions = 0;


  if (this->all_positions != NULL) {
    FREE(this->all_positions);
  }
  this->all_positions = (Univcoord_T *) NULL;
  this->n_all_positions = 0;

  this->status = ELT_VIRGIN;

  return;
}
#endif


static void
Elt_free (Elt_T *old) {

  if ((*old)->positions_allocated != NULL) {
    FREE((*old)->positions_allocated);
  }
  if ((*old)->all_positions != NULL) {
    FREE((*old)->all_positions);
  }
  FREE(*old);
  return;
}


#if 0
static int
Elt_nmatches_cmp (const void *a, const void *b) {
  Elt_T x = * (Elt_T *) a;
  Elt_T y = * (Elt_T *) b;

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else {
    return 0;
  }
}
#endif

#if 0
static int
Elt_querypos_ascending_cmp (const void *a, const void *b) {
  Elt_T x = * (Elt_T *) a;
  Elt_T y = * (Elt_T *) b;

  if (x->querystart < y->querystart) {
    return -1;
  } else if (y->querystart < x->querystart) {
    return +1;
  } else {
    return 0;
  }
}
#endif

#if 0
static int
Elt_querypos_descending_cmp (const void *a, const void *b) {
  Elt_T x = * (Elt_T *) a;
  Elt_T y = * (Elt_T *) b;

  if (x->querystart > y->querystart) {
    return -1;
  } else if (y->querystart > x->querystart) {
    return +1;
  } else {
    return 0;
  }
}
#endif


static int
Elt_extend_leftward (int *min_leftward, Elt_T elt, Compress_T query_compress,
		     bool plusp, int genestrand, int skip_left) {
  int max_leftward, nmatches;
  int i;

  if (elt->npositions == 0) {
    *min_leftward = 0;
    return 0;
  } else {
    max_leftward = *min_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/elt->positions[0],
								       /*pos5*/0,/*pos3*/elt->querystart - skip_left,
								       plusp,genestrand);
    for (i = 1; i < elt->npositions; i++) {
      if ((nmatches = Genome_consecutive_matches_leftward(query_compress,/*left*/elt->positions[i],
							  /*pos5*/0,/*pos3*/elt->querystart,
							  plusp,genestrand)) > max_leftward) {
	max_leftward = nmatches;
      } else if (nmatches < *min_leftward) {
	*min_leftward = nmatches;
      }
    }
    return max_leftward;
  }
}


static void
Elt_fill_positions_all (Elt_T this, T sarray) {
  Sarrayptr_T ptr;
  Univcoord_T pos;
  int i;

  debug7(printf("Entering Elt_fill_positions_all on %p\n",this));
  if (this->positions_allocated != NULL) {
    debug7(printf("  positions_allocated is already non-NULL, so skipping\n"));
    /* Don't free positions_allocated.  Use it. */

  } else {
    this->npositions_allocated = this->npositions = this->finalptr - this->initptr + 1;
    debug7(printf("  filling %d positions\n",this->npositions));

    if (this->nmatches == 0 || this->npositions > EXCESS_SARRAY_HITS) {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
      this->npositions_allocated = this->npositions = 0;
    } else {
#ifdef USE_QSORT
      this->positions_allocated = this->positions = (Univcoord_T *) MALLOC(this->npositions * sizeof(Univcoord_T));
#else
      this->positions_allocated = this->positions = (Univcoord_T *) MALLOC((this->npositions + 1) * sizeof(Univcoord_T));
#endif
      i = 0;
      ptr = this->initptr;
      while (ptr <= this->finalptr) {
	if ((pos = csa_lookup(sarray,ptr++)) >= (Univcoord_T) this->querystart) {
	  this->positions[i++] = pos - this->querystart;
	}
      }
      this->npositions = i;
#ifdef USE_QSORT
      qsort(this->positions,this->npositions,sizeof(Univcoord_T),Univcoord_compare);
#else
      Sedgesort_uint4(this->positions,this->npositions);
#endif
    }
  }

  this->status = ELT_FILLED;
  return;
}


#ifdef DEBUG7
static void
print_vector (__m128i x, char *label) {
  __m128i a[1];
  unsigned int *s = a;

  _mm_store_si128(a,x);
  printf("%s: %u %u %u %u\n",label,s[0],s[1],s[2],s[3]);
  return;
}

static void
print_vector_looking (__m128i x, Univcoord_T low, Univcoord_T high) {
  __m128i a[1];
  unsigned int *s = a;

  _mm_store_si128(a,x);
  printf("Looking at value %u, relative to low %u and high %u\n",s[0],low,high);
  printf("Looking at value %u, relative to low %u and high %u\n",s[1],low,high);
  printf("Looking at value %u, relative to low %u and high %u\n",s[2],low,high);
  printf("Looking at value %u, relative to low %u and high %u\n",s[3],low,high);
  return;
}
#endif


#ifdef DEBUG8
/* Non-SIMD methods for comparison */
static void
positions_compare (Univcoord_T *positions, int npositions,
		   Univcoord_T *positions_std, int npositions_std) {
  int i;
  bool problemp = false;

  if (npositions != npositions_std) {
    fprintf(stderr,"npositions %d != npositions_std %d\n",npositions,npositions_std);
    for (i = 0; i < npositions; i++) {
      printf("%u\n",positions[i]);
    }
    printf("\n");

    for (i = 0; i < npositions_std; i++) {
      printf("%u\n",positions_std[i]);
    }
    printf("\n");
    abort();

  } else {
    qsort(positions,npositions,sizeof(Univcoord_T),Univcoord_compare);
    qsort(positions_std,npositions,sizeof(Univcoord_T),Univcoord_compare);
    for (i = 0; i < npositions; i++) {
      if (positions[i] != positions_std[i]) {
	fprintf(stderr,"At %d, positions %u != positions_std %u\n",i,positions[i],positions_std[i]);
	problemp = true;
      }
    }
    if (problemp == true) {
      abort();
    }
  }

  return;
}
#endif


#ifdef DEBUG8
static Univcoord_T *
fill_positions_std (int *npositions, Univcoord_T low_adj, Univcoord_T high_adj,
		    Sarrayptr_T initptr, Sarrayptr_T finalptr,
		    int querystart, Univcoord_T *array) {
  Univcoord_T *more_positions;
  Univcoord_T *positions, value;
  Sarrayptr_T ptr, lastptr;
  int i;

  positions = (Univcoord_T *) MALLOC(GUESS_ALLOCATION * sizeof(Univcoord_T)); /* Return value, so cannot use alloca */

  *npositions = 0;
  ptr = initptr;      

  while (ptr <= finalptr) {
    debug7a(printf("Std: Looking at value %u, relative to low %u and high %u\n",array[ptr],low_adj,high_adj));
    if ((value = CONVERT(array[ptr++])) < low_adj) {
      /* Skip */
    } else if (value > high_adj) {
      /* Skip */
    } else if (*npositions < GUESS_ALLOCATION) {
      debug7(printf("Std: Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
      positions[(*npositions)++] = value - querystart;
    } else {
      debug7(printf("Std: Found position %u between low %u and high %u, but exceeds allocation\n",value,low_adj,high_adj));
      (*npositions)++;
      lastptr = ptr;		/* saves us from going through the entire sarray below */
    }
  }

  debug7(printf("Std method found %d positions\n",*npositions));
  if (*npositions > GUESS_ALLOCATION) {
    /* Copy the positions we have stored so far */
    more_positions = (Univcoord_T *) MALLOC((*npositions) * sizeof(Univcoord_T));
    memcpy(more_positions,positions,GUESS_ALLOCATION*sizeof(Univcoord_T));
    FREE(positions);
    positions = more_positions;
    
    i = GUESS_ALLOCATION;	/* Start count with the number stored */
    ptr = lastptr;	/* One past the last ptr with a result */

    while (i < *npositions) {
      if ((value = CONVERT(array[--ptr])) < low_adj) {
	/* Skip */
      } else if (value > high_adj) {
	/* Skip */
      } else {
	positions[i++] = value - querystart;
      }
    }
  }

  return positions;
}
#endif



/* Call fill_positions_filtered_first for first time, which is
   linear in number of entries or O(n), then on second call, do sort with O(n*log n),
   plus O(log n) for each additional call */

#ifdef HAVE_ALLOCA

#if defined(HAVE_AVX2) && !defined(WORDS_BIGENDIAN)

/* Using pext method, because _mm256_shuffle_epi32 doesn't work well
because it works only within lanes, and MASTER_CONTROL does not extend
well to 256 bits */

static void
fill_positions_filtered_first (Elt_T this, T sarray, Univcoord_T low, Univcoord_T high) {
  Univcoord_T low_adj, high_adj;
  Univcoord_T *array = sarray->array, value0;
  Sarrayptr_T *array_stop, *array_end, *array_ptr;
  Univcoord_T *positions_temp;
  Univcoord_T *out;
  __m256i converted, adjusted, match;
  __m256i floor, ceiling, values, adj, p;
  unsigned int mask, selector;
  unsigned int univ_selector = 0x76543210;
  __m256i _selector, _varshifts;
#if defined(REQUIRE_ALIGNMENT)
  int n_prealign, k;
#endif
#if defined(DEBUG) || defined(DEBUG7)
  int nmatches;
#endif
#ifdef DEBUG7
  UINT8 pointer;
  int i;
#endif
#ifdef DEBUG8
  Univcoord_T *positions_std;
  int npositions_std;
#endif


  debug(printf("Entered fill_positions_filtered_first with low %u and high %u, initptr %u and finalptr %u (n = %d), nmatches %d\n",
	       low,high,this->initptr,this->finalptr,this->finalptr - this->initptr + 1,this->nmatches));
  debug7(printf("Entered fill_positions_filtered_first with low %u and high %u, initptr %u and finalptr %u (n = %d), nmatches %d\n",
		low,high,this->initptr,this->finalptr,this->finalptr - this->initptr + 1,this->nmatches));
  
  if (this->positions_allocated != NULL) {
    /* Filled from a previous call */
    FREE(this->positions_allocated);
  }

  if ((this->n_all_positions = this->finalptr - this->initptr + 1) == 0 /*|| this->n_all_positions > EXCESS_SARRAY_HITS*/) {
    this->all_positions = (Univcoord_T *) NULL;

  } else {
    /* Function surrounded by HAVE_ALLOCA */
#ifdef USE_QSORT
    positions_temp = out = (Univcoord_T *) MALLOCA((this->finalptr - this->initptr + 1) * sizeof(Univcoord_T));
#else
    positions_temp = out = (Univcoord_T *) MALLOCA((this->finalptr - this->initptr + 1 + 1) * sizeof(Univcoord_T));
#endif

    low_adj = low + this->querystart;
    high_adj = high + this->querystart;
  
    floor = _mm256_set1_epi32(low_adj - 1 - 2147483648);
    ceiling = _mm256_set1_epi32(high_adj + 1 - 2147483648);
    adj = _mm256_set1_epi32(this->querystart);
    _varshifts = _mm256_set_epi32(28,24,20,16,12,8,4,0);

    this->npositions_allocated = this->npositions = 0;
#if defined(REQUIRE_ALIGNMENT)
    array_ptr = &(array[this->initptr]);
    
    /* Initial part */
    n_prealign = ((32 - ((UINT8) array_ptr & 0x1F))/4) & 0x7;
    debug7(printf("Initial ptr is at location %p.  Need %d to get to 256-bit boundary\n",pointer,n_prealign));

    debug7(printf("Initial part:\n"));
    if (n_prealign > this->finalptr - this->initptr + 1) {
      n_prealign = this->finalptr - this->initptr + 1;
    }
    for (k = 0; k < n_prealign; k++) {
      debug7a(printf("Looking at value %u, relative to low %u and high %u\n",CONVERT(array[ptr]),low_adj,high_adj));
      if ((value0 = *array_ptr++) >= low_adj && value0 <= high_adj) {
	*out++ = value0 - this->querystart;
      }
    }
#else
    array_ptr = &(array[this->initptr]);
#endif	/* REQUIRE_ALIGNMENT */


    /* Aligned part */
    if (this->finalptr < 8) {
      array_stop = &(array[0]);
    } else {
      array_stop = &(array[this->finalptr - 8]);
    }
    array_end = &(array[this->finalptr]);

    while (array_ptr < array_stop) {

#if defined(REQUIRE_ALIGNMENT)
      /* Use stream_load to avoid polluting the cache with suffix array entries */
      values = _mm256_stream_load_si256((__m256i *) array_ptr);
#else
      /* It looks like loadu is just as fast as load */
      values = _mm256_loadu_si256((__m256i *) array_ptr);
#endif
      debug7b(print_vector_uint_256(values));

      converted = _mm256_sub_epi32(values,epi32_convert_256);
      /* match = _mm256_andnot_si256(_mm256_cmpgt_epi32(floor,converted),_mm256_cmpgt_epi32(ceiling,converted)); -- This is off by 1 at floor */
      match = _mm256_and_si256(_mm256_cmpgt_epi32(converted,floor),_mm256_cmpgt_epi32(ceiling,converted));
      /* Example:
	 {0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0x00000000} (8 x 32-bit ints) */
      debug7b(print_vector_hex_256(match));

      /* Get most significant from each byte from each byte to obtain a mask in terms of 4-bit nibbles */
      mask = (unsigned int) _mm256_movemask_epi8(match);
      /* 11110000 11110000 11110000 11110000  (32-bit int) */

      if (mask) {
	adjusted = _mm256_sub_epi32(values,adj);

	selector = _pext_u32(univ_selector,mask); /* Requires compilation with -mbmi2 */
	/* 01110110 01010100 00110010 00010000  (32-bit int, univ_selector) */
	/* 00000000 00000000 01110101 00110001  (32-bit int, selector), equiv to 0x00007531 */
	debug7b(printf("selector: %08X\n",selector));

	_selector = _mm256_set1_epi32(selector);
	/* {0x00007531, 0x00007531, 0x00007531, 0x00007531, 0x00007531, 0x00007531, 0x00007531, 0x00007531} */

	_selector = _mm256_srlv_epi32(_selector,_varshifts);
	/* {0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000007, 0x00000075, 0x00000753, 0x00007531} */
	debug7b(print_vector_hex_256(_selector));

	p = _mm256_permutevar8x32_epi32(adjusted,_selector);
	_mm256_storeu_si256((__m256i *) out, p);

	/* Divide popcount(mask) by 4 to get number of matches */
#ifdef HAVE_POPCNT
	out += _popcnt32(mask)/4;
	debug7b(printf("mask: %08X (%d ones)\n",mask,_popcnt32(mask)/4));
#elif defined HAVE_MM_POPCNT
	out += _mm_popcnt_u32(mask)/4;
	debug7b(printf("mask: %08X (%d ones)\n",mask,_mm_popcnt_u32(mask)/4));
#else
	out += __builtin_popcount(mask)/4;
	debug7b(printf("mask: %08X (%d ones)\n",mask,__builtin_popcount(mask)/4));
#endif
	debug7b(print_vector_uint_256(p));
      }

      array_ptr += 8;
    }

    /* Partial block at end; do scalar */
    debug7(printf("\nFinal part:\n"));
    while (array_ptr <= array_end) {
      if ((value0 = *array_ptr++) >= low_adj && value0 <= high_adj) {
	*out++ = value0 - this->querystart;
      }
    }

    this->npositions_allocated = this->npositions = out - positions_temp;
    debug7(printf("SIMD method found %d positions\n",this->npositions));

    /* Copy the positions into heap from temp in stack */
    if (this->npositions == 0) {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
    } else {
      debug7(printf("Sorting %d positions\n",this->npositions));
#ifdef USE_QSORT
      qsort(positions_temp,this->npositions,sizeof(Univcoord_T),Univcoord_compare);
#else
      Sedgesort_uint4(positions_temp,this->npositions);
#endif

      /* Need to copy positions before the goal */
#ifdef USE_QSORT
      this->positions_allocated = this->positions = MALLOC(this->npositions * sizeof(Univcoord_T));
#else
      this->positions_allocated = this->positions = MALLOC((this->npositions + 1) * sizeof(Univcoord_T));
#endif
      memcpy(this->positions,positions_temp,this->npositions * sizeof(Univcoord_T));
#ifdef DEBUG7
      for (i = 0; i < this->npositions; i++) {
	printf("%u\n",this->positions[i]);
      }
#endif

#if 0
      /* Not sure why we were doing this.  We will find collinear set of diagonals later. */
      /* Advance pointer to goal (note: do not want goal_adj, since we have already subtracted this->querystart) */
      /* Have tested positions[i] <= goal, but want positions[-1] to be < goal, or positions[0] >= goal */
      /* ? Replace with a binary search */
      i = 0;
      while (i < this->npositions && positions_temp[i] < goal) {
	debug7(printf("1 Skipping position %u (%u) < goal %u (%u)\n",
		      positions_temp[i],positions_temp[i] - chroffset,goal,goal - chroffset));
	i++;
      }
      this->positions += i;
      this->npositions -= i;
      debug7(printf("Remaining: %d positions\n",this->npositions));
#endif
    }
    
    /* Function surrounded by HAVE_ALLOCA */
    FREEA(positions_temp);
  }

  return;
}


#elif defined(HAVE_SSSE3) && !defined(WORDS_BIGENDIAN)
/* SSSE3 needed for _mm_shuffle_epi8 */

/* Prefer shuffle_mask, because MASTER_CONTROL requires AVX command
   _mm_permutevar_ps.  Cannot use pext approach, because that requires
   BMI2. */

/* Nevertheless, here is an explanation of MASTER_CONTROL:

   For blocks of 128 bits:

 MASTER_CONTROL =
    _mm_set_epi8(0x10, 0x12, 0x13, 0x12, 0x40, 0x68, 0x7C, 0x6B,
                 0x00, 0x80, 0xC0, 0xBC, 0x00, 0x00, 0x00, 0xC0);

  matchbits = _mm_movemask_ps(_mm_castsi128_ps(match));
  p = _mm_permutevar_ps(input,_mm_srli_epi32(MASTER_CONTROL,matchbits*2));

  These values come from rotating the control values horizontally, like this:

  15 14 13 12  11 10 09 08  07 06 05 04  03 02 01 00  (matchbits)

  00 01 00 10  00 01 00 11  00 01 00 10  00 01 00 XX  (control 15 downto 0, bits 1 and 0)
  01 10 10 11  01 11 11 XX  01 10 10 XX  01 XX XX XX  (control 15 downto 0, bits 3 and 2)
  10 11 11 XX  11 XX XX XX  10 XX XX XX  XX XX XX XX  (control 15 downto 0, bits 5 and 4)
  11 XX XX XX  XX XX XX XX  XX XX XX XX  XX XX XX XX  (control 15 downto 0, bits 7 and 6)

  and then reading in bytes from right to left for each row.

  Alternatively, create the following table, and read in bits upward
  starting from the rightmost column:

  matchbits 00: XX XX XX XX ^
  matchbits 01: XX XX XX 00 |
  matchbits 02: XX XX XX 01 |
  matchbits 03: XX XX 01 00 _
  matchbits 04: XX XX XX 10 ^
  matchbits 05: XX XX 10 00 |
  matchbits 06: XX XX 10 01 |
  matchbits 07: XX 10 01 00 _
  matchbits 08: XX XX XX 11 .
  matchbits 09: XX XX 11 00 .
  matchbits 10: XX XX 11 01 .
  matchbits 11: XX 11 01 00
  matchbits 12: XX XX 11 10
  matchbits 13: XX 11 10 00
  matchbits 14: XX 11 10 01
  matchbits 15: 11 10 01 00
*/

static void
fill_positions_filtered_first (Elt_T this, T sarray, Univcoord_T low, Univcoord_T high) {
  Univcoord_T low_adj, high_adj;
  Univcoord_T *array = sarray->array, value0;
  Sarrayptr_T *array_stop, *array_end, *array_ptr;
  Univcoord_T *positions_temp;
  Univcoord_T *out;
  __m128i converted, adjusted, match;
  __m128i floor, ceiling, values, adj, p;
  int matchbits;
#if defined(REQUIRE_ALIGNMENT)
  int n_prealign, k;
#endif
#ifndef USE_SHUFFLE_MASK
  __m128i MASTER_CONTROL;
#endif
#ifdef DEBUG7
  int i;
#endif
#ifdef DEBUG8
  Univcoord_T *positions_std;
  int npositions_std;
#endif


  debug(printf("Entered fill_positions_filtered_first with low %u and high %u, initptr %u and finalptr %u (n = %d), nmatches %d\n",
	       low,high,this->initptr,this->finalptr,this->finalptr - this->initptr + 1,this->nmatches));
  debug7(printf("Entered fill_positions_filtered_first with low %u and high %u, initptr %u and finalptr %u (n = %d), nmatches %d\n",
		low,high,this->initptr,this->finalptr,this->finalptr - this->initptr + 1,this->nmatches));
  
  if (this->positions_allocated != NULL) {
    /* Filled from a previous call */
    FREE(this->positions_allocated);
  }

  if ((this->n_all_positions = this->finalptr - this->initptr + 1) == 0 /*|| this->n_all_positions > EXCESS_SARRAY_HITS*/) {
    this->all_positions = (Univcoord_T *) NULL;

  } else {
    /* Function surrounded by HAVE_ALLOCA */
    positions_temp = out = (Univcoord_T *) MALLOCA((this->finalptr - this->initptr + 1) * sizeof(Univcoord_T));

    low_adj = low + this->querystart;
    high_adj = high + this->querystart;
  
    floor = _mm_set1_epi32(low_adj - 1 - 2147483648);
    ceiling = _mm_set1_epi32(high_adj + 1 - 2147483648);
    adj = _mm_set1_epi32(this->querystart);

    this->npositions_allocated = this->npositions = 0;
#if defined(REQUIRE_ALIGNMENT)
    array_ptr = &(array[this->initptr]);
    
    /* Initial part */
#ifdef HAVE_64_BIT
    n_prealign = ((16 - ((UINT8) array_ptr & 0xF))/4) & 0x3;
#else
    n_prealign = ((16 - ((UINT4) array_ptr & 0xF))/4) & 0x3;
#endif
    debug7(printf("Initial ptr is at location %p.  Need %d to get to 128-bit boundary\n",pointer,n_prealign));

    debug7(printf("Initial part:\n"));
    if (n_prealign > this->finalptr - this->initptr + 1) {
      n_prealign = this->finalptr - this->initptr + 1;
    }
    for (k = 0; k < n_prealign; k++) {
      debug7a(printf("Looking at value %u, relative to low %u and high %u\n",CONVERT(array[ptr]),low_adj,high_adj));
      if ((value0 = *array_ptr++) >= low_adj && value0 <= high_adj) {
	*out++ = value0 - this->querystart;
      }
    }
#else
    array_ptr = &(array[this->initptr]);
#endif	/* REQUIRE_ALIGNMENT */


    /* Aligned part */
    if (this->finalptr < 4) {
      array_stop = &(array[0]);
    } else {
      array_stop = &(array[this->finalptr - 4]);
    }
    array_end = &(array[this->finalptr]);

#ifndef USE_SHUFFLE_MASK
    MASTER_CONTROL = _mm_setr_epi8(0x10, 0x12, 0x13, 0x12, 0x40, 0x68, 0x7C, 0x6B,
				   0x00, 0x80, 0xC0, 0xBC, 0x00, 0x00, 0x00, 0xC0);
#endif

    while (array_ptr < array_stop) {
#if defined(REQUIRE_ALIGNMENT)

#ifdef HAVE_SSE4_1      
      /* Use stream_load to avoid polluting the cache with suffix array entries */
      values = _mm_stream_load_si128((__m128i *) array_ptr);
#else
      values = _mm_load_si128((__m128i *) array_ptr);
#endif

#else
      /* It looks like loadu is just as fast as load */
      values = _mm_loadu_si128((__m128i *) array_ptr);
#endif
      debug7b(print_vector_uint(values));

      converted = _mm_sub_epi32(values,epi32_convert);
      /* match = _mm_andnot_si128(_mm_cmpgt_epi32(floor,converted),_mm_cmpgt_epi32(ceiling,converted)); -- This is off by 1 at floor */
      match = _mm_and_si128(_mm_cmpgt_epi32(converted,floor),_mm_cmplt_epi32(converted,ceiling));
      debug7b(print_vector_hex(match));

      matchbits = _mm_movemask_ps(_mm_castsi128_ps(match));
      if (matchbits) {
	adjusted = _mm_sub_epi32(values,adj);
#ifdef USE_SHUFFLE_MASK
	p = _mm_shuffle_epi8(adjusted, shuffle_mask16[matchbits]);
#else
	p = _mm_castps_si128(_mm_permutevar_ps(_mm_castsi128_ps(adjusted),_mm_srli_epi32(MASTER_CONTROL,matchbits*2)));
#endif
	_mm_storeu_si128((__m128i *) out, p);

#if !defined(HAVE_SSE4_2)
	out += __builtin_popcount(matchbits);
	debug7b(printf("matchbits: %08X (%d ones)\n",matchbits,__builtin_popcount(matchbits)));
#elif defined(HAVE_POPCNT)
	out += _popcnt32(matchbits);
	debug7b(printf("matchbits: %08X (%d ones)\n",matchbits,_popcnt32(matchbits)));
#elif defined HAVE_MM_POPCNT
	out += _mm_popcnt_u32(matchbits);
	debug7b(printf("matchbits: %08X (%d ones)\n",matchbits,_mm_popcnt_u32(matchbits)));
#else
	out += __builtin_popcount(matchbits);
	debug7b(printf("matchbits: %08X (%d ones)\n",matchbits,__builtin_popcount(matchbits)));
#endif
	debug7b(print_vector_hex(shuffle_mask16[matchbits]));
	debug7b(print_vector_uint(p));
      }

      array_ptr += 4;
    }

    /* Partial block at end; do scalar */
    debug7(printf("\nFinal part:\n"));
    while (array_ptr <= array_end) {
      if ((value0 = *array_ptr++) >= low_adj && value0 <= high_adj) {
	*out++ = value0 - this->querystart;
      }
    }

    this->npositions_allocated = this->npositions = out - positions_temp;
    debug7(printf("SIMD method found %d positions\n",this->npositions));

    /* Copy the positions into heap from temp in stack */
    if (this->npositions == 0) {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
    } else {
      debug7(printf("Sorting %d positions\n",this->npositions));
#ifdef USE_QSORT
      qsort(positions_temp,this->npositions,sizeof(Univcoord_T),Univcoord_compare);
#else
      Sedgesort_uint4(positions_temp,this->npositions);
#endif

      /* Need to copy positions before the goal */
#ifdef USE_QSORT
      this->positions_allocated = this->positions = MALLOC(this->npositions * sizeof(Univcoord_T));
#else
      this->positions_allocated = this->positions = MALLOC((this->npositions + 1) * sizeof(Univcoord_T));
#endif
      memcpy(this->positions,positions_temp,this->npositions * sizeof(Univcoord_T));
#ifdef DEBUG7
      for (i = 0; i < this->npositions; i++) {
	printf("%u\n",this->positions[i]);
      }
#endif

#if 0
      /* Not sure why we were doing this.  We will find collinear set of diagonals later. */
      /* Advance pointer to goal (note: do not want goal_adj, since we have already subtracted this->querystart) */
      /* Have tested positions[i] <= goal, but want positions[-1] to be < goal, or positions[0] >= goal */
      /* ? Replace with a binary search */
      i = 0;
      while (i < this->npositions && positions_temp[i] < goal) {
	debug7(printf("1 Skipping position %u (%u) < goal %u (%u)\n",
		      positions_temp[i],positions_temp[i] - chroffset,goal,goal - chroffset));
	i++;
      }
      this->positions += i;
      this->npositions -= i;
      debug7(printf("Remaining: %d positions\n",this->npositions));
#endif
    }
    
    /* Function surrounded by HAVE_ALLOCA */
    FREEA(positions_temp);
  }

  return;
}


#else
/* Bigendian or missing SSSE3 */

static void
fill_positions_filtered_first (Elt_T this, T sarray, Univcoord_T low, Univcoord_T high) {
  Sarrayptr_T ptr;
  Univcoord_T low_adj, high_adj;
  Univcoord_T *array = sarray->array;
  Univcoord_T value3, value2, value1, value0;
  Univcoord_T *positions_temp;
#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
#ifdef HAVE_64_BIT
  UINT8 pointer;
#else
  UINT4 pointer;
#endif
  __m128i floor, ceiling, values, compare;
  int n_prealign, k;
#endif


  debug7(printf("Entered fill_positions_filtered_first with low %u and high %u, initptr %u and finalptr %u (n = %d), nmatches %d\n",
		low,high,this->initptr,this->finalptr,this->finalptr - this->initptr + 1,this->nmatches));
  
  if (this->positions_allocated != NULL) {
    /* Filled from a previous call */
    FREE(this->positions_allocated);
  }

  if ((this->n_all_positions = this->finalptr - this->initptr + 1) == 0 /*|| this->n_all_positions > EXCESS_SARRAY_HITS*/) {
    this->all_positions = (Univcoord_T *) NULL;

  } else {
    /* Function surrounded by HAVE_ALLOCA */
#ifdef USE_QSORT
    positions_temp = (Univcoord_T *) MALLOCA((this->finalptr - this->initptr + 1) * sizeof(Univcoord_T));
#else
    positions_temp = (Univcoord_T *) MALLOCA((this->finalptr - this->initptr + 1 + 1) * sizeof(Univcoord_T));
#endif

    low_adj = low + this->querystart;
    high_adj = high + this->querystart;

    this->npositions_allocated = this->npositions = 0;
    ptr = this->initptr;
#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
    if (ptr + 3 > this->finalptr) { /* ptr + 4 > (this->finalptr + 1) */
      /* Handle in normal manner */
      debug7(printf("Small batch, because %u + 3 <= %u\n",ptr,this->finalptr));
      while (ptr <= this->finalptr) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",csa_lookup(sarray,ptr),low_adj,high_adj));
	if ((value0 = csa_lookup(sarray,ptr++)) < low_adj) {
	  /* Skip */
	} else if (value0 > high_adj) {
	  /* Skip */
	} else {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	  positions_temp[this->npositions++] = value0 - this->querystart;
	}
      }

    } else {
#ifdef HAVE_64_BIT
      pointer = (UINT8) &(array[ptr]);
#else
      pointer = (UINT4) &(array[ptr]);
#endif
      n_prealign = ((16 - (pointer & 0xF))/4) & 0x3;
      debug7(printf("Initial ptr is at location %p.  Need %d to get to 128-bit boundary\n",
		    &(array[ptr]),n_prealign));

      /* Initial part */
      debug7(printf("Initial part:\n"));
      for (k = 0; k < n_prealign; k++) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",CONVERT(array[ptr]),low_adj,high_adj));
	if ((value0 = CONVERT(array[ptr++])) < low_adj) {
	  /* Skip */
	} else if (value0 > high_adj) {
	  /* Skip */
	} else {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	  positions_temp[this->npositions++] = value0 - this->querystart;
	}
      }

      /* Aligned part */
      debug7(printf("\nAligned part:\n"));
      /* Since compare operations not available for unsigned ints, using the fact that
	 unsigned_gt(a,b) is equivalent to signed_gt(a - 2^31, b - 2^31) */
      floor = _mm_set1_epi32(low_adj - 1 - 2147483648);
      ceiling = _mm_set1_epi32(high_adj + 1 - 2147483648);
      while (ptr + 3 <= this->finalptr) { /* ptr + 4 < this->finalptr + 1 */
	/* Missing SSSE3 so cannot do stream load */
	values = _mm_load_si128((__m128i *) &(array[ptr]));
	debug7a(print_vector_looking(values,low_adj,high_adj));
	values = _mm_sub_epi32(values,epi32_convert);
	compare = _mm_and_si128(_mm_cmpgt_epi32(values,floor),_mm_cmplt_epi32(values,ceiling));
	if (/*cmp*/_mm_movemask_epi8(compare) == 0x0000) {
	  /* All results are false, indicating no values between low_adj and high_adj (most common case) */
	  ptr += 4;
	} else {
	  value3 = CONVERT(array[ptr++]);
	  if (value3 < low_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u < low %u\n",value3,low_adj));
	  } else if (value3 > high_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u > high %u\n",value3,high_adj));
	  } else {
	    debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value3,low_adj,high_adj));
	    positions_temp[this->npositions++] = value3 - this->querystart;
	  }

	  value2 = CONVERT(array[ptr++]);
	  if (value2 < low_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u < low %u\n",value2,low_adj));
	  } else if (value2 > high_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u > high %u\n",value2,high_adj));
	  } else {
	    debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value2,low_adj,high_adj));
	    positions_temp[this->npositions++] = value2 - this->querystart;
	  }

	  value1 = CONVERT(array[ptr++]);
	  if (value1 < low_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u < low %u\n",value1,low_adj));
	  } else if (value1 > high_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u > high %u\n",value1,high_adj));
	  } else {
	    debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value1,low_adj,high_adj));
	    positions_temp[this->npositions++] = value1 - this->querystart;
	  }

	  value0 = CONVERT(array[ptr++]);
	  if (value0 < low_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u < low %u\n",value0,low_adj));
	  } else if (value0 > high_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u > high %u\n",value0,high_adj));
	  } else {
	    debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	    positions_temp[this->npositions++] = value0 - this->querystart;
	  }
	}
      }

      /* Final part */
      debug7(printf("\nFinal part:\n"));
      while (ptr <= this->finalptr) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",csa_lookup(sarray,ptr),low_adj,high_adj));
	if ((value0 = csa_lookup(sarray,ptr++)) < low_adj) {
	  /* Skip */
	} else if (value0 > high_adj) {
	  /* Skip */
	} else {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	  positions_temp[this->npositions++] = value0 - this->querystart;
	}
      }
    }

#else

    while (ptr <= this->finalptr) {
      debug7a(printf("Looking at value %u, relative to low %u and high %u\n",csa_lookup(sarray,ptr),low_adj,high_adj));
      if ((value0 = csa_lookup(sarray,ptr++)) < low_adj) {
	/* Skip */
      } else if (value0 > high_adj) {
	/* Skip */
      } else {
	debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	positions_temp[this->npositions++] = value0 - this->querystart;
      }
    }
#endif

    debug7(printf("SIMD method found %d positions\n",this->npositions));

    /* Copy the positions from temp */
    if (this->npositions == 0) {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
    } else {
      debug7(printf("Sorting %d positions\n",this->npositions));
#ifdef USE_QSORT
      qsort(positions_temp,this->npositions,sizeof(Univcoord_T),Univcoord_compare);
#else
      Sedgesort_uint4(positions_temp,this->npositions);
#endif

      /* Need to copy positions before the goal */
      this->positions_allocated = this->positions = MALLOC(this->npositions * sizeof(Univcoord_T));
      memcpy(this->positions,positions_temp,this->npositions * sizeof(Univcoord_T));

#if 0
      /* Not sure why we were doing this.  We will find collinear set of diagonals later. */
      /* Advance pointer to goal (note: do not want goal_adj, since we have already subtracted this->querystart) */
      /* Have tested positions[i] <= goal, but want positions[-1] to be < goal, or positions[0] >= goal */
      /* ? Replace with a binary search */
      i = 0;
      while (i < this->npositions && positions_temp[i] < goal) {
	debug7(printf("2 Skipping position %u < goal %u\n",positions_temp[i] - chroffset,goal - chroffset));
	i++;
      }
      this->positions += i;
      this->npositions -= i;
      debug7(printf("Remaining: %d positions\n",this->npositions));
#endif
    }
    
    /* Function surrounded by HAVE_ALLOCA */
    FREEA(positions_temp);
  }

  return;
}
#endif


#else
/* Non-ALLOCA version */

static void
fill_positions_filtered_first (Elt_T this, T sarray, Univcoord_T low, Univcoord_T high) {
  Sarrayptr_T ptr, lastptr;
  int nmatches;
  int i;
  Univcoord_T low_adj, high_adj;
  Univcoord_T value3, value2, value1, value0;
  Univcoord_T *array = sarray->array;
  Univcoord_T *more_positions;
#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
#ifdef HAVE_64_BIT
  UINT8 pointer;
#else
  UINT4 pointer;
#endif
  __m128i floor, ceiling, values, compare;
  int n_prealign, k;
#endif


  debug7(printf("Entered fill_positions_filtered_first with low %u and high %u, initptr %u and finalptr %u (n = %d), nmatches %d\n",
		low,high,this->initptr,this->finalptr,this->finalptr - this->initptr + 1,this->nmatches));
  
  if (this->positions_allocated != NULL) {
    /* Filled from a previous call */
    FREE(this->positions_allocated);
  }

  if ((this->n_all_positions = this->finalptr - this->initptr + 1) == 0 /*|| this->n_all_positions > EXCESS_SARRAY_HITS*/) {
    this->all_positions = (Univcoord_T *) NULL;

  } else {
    /* Guess at allocation size */
#ifdef USE_QSORT
    this->positions_allocated = this->positions = (Univcoord_T *) MALLOC(GUESS_ALLOCATION * sizeof(Univcoord_T));
#else
    this->positions_allocated = this->positions = (Univcoord_T *) MALLOC((GUESS_ALLOCATION + 1) * sizeof(Univcoord_T));
#endif

    low_adj = low + this->querystart;
    high_adj = high + this->querystart;

    this->npositions_allocated = this->npositions = 0;
    ptr = this->initptr;
#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
    if (ptr + 3 > this->finalptr) { /* ptr + 4 > (this->finalptr + 1) */
      /* Handle in normal manner */
      debug7(printf("Small batch, because %u + 3 <= %u\n",ptr,this->finalptr));
      while (ptr <= this->finalptr) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",csa_lookup(sarray,ptr),low_adj,high_adj));
	if ((value0 = csa_lookup(sarray,ptr++)) < low_adj) {
	  /* Skip */
	} else if (value0 > high_adj) {
	  /* Skip */
	} else if (this->npositions < GUESS_ALLOCATION) {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	  this->positions[this->npositions++] = value0 - this->querystart;
	} else {
	  debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value0,low_adj,high_adj));
	  this->npositions++;
	  lastptr = ptr;		/* saves us from going through the entire sarray below */
	}
      }

    } else {
#ifdef HAVE_64_BIT
      pointer = (UINT8) &(array[ptr]);
#else
      pointer = (UINT4) &(array[ptr]);
#endif
      n_prealign = ((16 - (pointer & 0xF))/4) & 0x3;
      debug7(printf("Initial ptr is at location %p.  Need %d to get to 128-bit boundary\n",
		    &(array[ptr]),n_prealign));

      /* Initial part */
      debug7(printf("Initial part:\n"));
      for (k = 0; k < n_prealign; k++) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",CONVERT(array[ptr]),low_adj,high_adj));
	if ((value0 = CONVERT(array[ptr++])) < low_adj) {
	  /* Skip */
	} else if (value0 > high_adj) {
	  /* Skip */
	} else if (this->npositions < GUESS_ALLOCATION) {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	  this->positions[this->npositions++] = value0 - this->querystart;
	} else {
	  debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value0,low_adj,high_adj));
	  this->npositions++;
	  lastptr = ptr;		/* saves us from going through the entire sarray below */
	}
      }

      /* Aligned part */
      debug7(printf("\nAligned part:\n"));
      /* Since compare operations not available for unsigned ints, using the fact that
	 unsigned_gt(a,b) is equivalent to signed_gt(a - 2^31, b - 2^31) */
      floor = _mm_set1_epi32(low_adj - 1 - 2147483648);
      ceiling = _mm_set1_epi32(high_adj + 1 - 2147483648);
      while (ptr + 3 <= this->finalptr) { /* ptr + 4 < this->finalptr + 1 */
	values = _mm_load_si128((__m128i *) &(array[ptr]));
	debug7a(print_vector_looking(values,low_adj,high_adj));
	values = _mm_sub_epi32(values,epi32_convert);
	compare = _mm_and_si128(_mm_cmpgt_epi32(values,floor),_mm_cmplt_epi32(values,ceiling));
	if (/*cmp*/_mm_movemask_epi8(compare) == 0x0000) {
	  /* All results are false, indicating no values between low_adj and high_adj (most common case) */
	  ptr += 4;
	} else {
	  value3 = CONVERT(array[ptr++]);
	  if (value3 < low_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u < low %u\n",value3,low_adj));
	  } else if (value3 > high_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u > high %u\n",value3,high_adj));
	  } else if (this->npositions < GUESS_ALLOCATION) {
	    debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value3,low_adj,high_adj));
	    this->positions[this->npositions++] = value3 - this->querystart;
	  } else {
	    debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value3,low_adj,high_adj));
	    this->npositions++;
	    lastptr = ptr;		/* saves us from going through the entire sarray below */
	  }

	  value2 = CONVERT(array[ptr++]);
	  if (value2 < low_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u < low %u\n",value2,low_adj));
	  } else if (value2 > high_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u > high %u\n",value2,high_adj));
	  } else if (this->npositions < GUESS_ALLOCATION) {
	    debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value2,low_adj,high_adj));
	    this->positions[this->npositions++] = value2 - this->querystart;
	  } else {
	    debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value2,low_adj,high_adj));
	    this->npositions++;
	    lastptr = ptr;		/* saves us from going through the entire sarray below */
	  }

	  value1 = CONVERT(array[ptr++]);
	  if (value1 < low_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u < low %u\n",value1,low_adj));
	  } else if (value1 > high_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u > high %u\n",value1,high_adj));
	  } else if (this->npositions < GUESS_ALLOCATION) {
	    debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value1,low_adj,high_adj));
	    this->positions[this->npositions++] = value1 - this->querystart;
	  } else {
	    debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value1,low_adj,high_adj));
	    this->npositions++;
	    lastptr = ptr;		/* saves us from going through the entire sarray below */
	  }

	  value0 = CONVERT(array[ptr++]);
	  if (value0 < low_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u < low %u\n",value0,low_adj));
	  } else if (value0 > high_adj) {
	    /* Skip */
	    debug7(printf("Skipping position %u > high %u\n",value0,high_adj));
	  } else if (this->npositions < GUESS_ALLOCATION) {
	    debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	    this->positions[this->npositions++] = value0 - this->querystart;
	  } else {
	    debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value0,low_adj,high_adj));
	    this->npositions++;
	    lastptr = ptr;		/* saves us from going through the entire sarray below */
	  }
	}
      }

      /* Final part */
      debug7(printf("\nFinal part:\n"));
      while (ptr <= this->finalptr) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",csa_lookup(sarray,ptr),low_adj,high_adj));
	if ((value0 = csa_lookup(sarray,ptr++)) < low_adj) {
	  /* Skip */
	} else if (value0 > high_adj) {
	  /* Skip */
	} else if (this->npositions < GUESS_ALLOCATION) {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	  this->positions[this->npositions++] = value0 - this->querystart;
	} else {
	  debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value0,low_adj,high_adj));
	  this->npositions++;
	  lastptr = ptr;		/* saves us from going through the entire sarray below */
	}
      }
    }

#else

    while (ptr <= this->finalptr) {
      debug7a(printf("Looking at value %u, relative to low %u and high %u\n",csa_lookup(sarray,ptr),low_adj,high_adj));
      if ((value0 = csa_lookup(sarray,ptr++)) < low_adj) {
	/* Skip */
      } else if (value0 > high_adj) {
	/* Skip */
      } else if (this->npositions < GUESS_ALLOCATION) {
	debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value0,low_adj,high_adj));
	this->positions[this->npositions++] = value0 - this->querystart;
      } else {
	debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value0,low_adj,high_adj));
	this->npositions++;
	lastptr = ptr;		/* saves us from going through the entire sarray below */
      }
    }
#endif

    debug7(printf("SIMD method found %d positions\n",this->npositions));
    if (this->npositions > GUESS_ALLOCATION) {
      /* Handle the case if we exceeded GUESS_ALLOCATION */

      /* Copy the positions we have stored so far */
#ifdef USE_QSORT
      more_positions = (Univcoord_T *) MALLOC(this->npositions * sizeof(Univcoord_T));
#else
      more_positions = (Univcoord_T *) MALLOC((this->npositions + 1) * sizeof(Univcoord_T));
#endif
      memcpy(more_positions,this->positions,GUESS_ALLOCATION*sizeof(Univcoord_T));
      FREE(this->positions_allocated);
      this->positions_allocated = this->positions = more_positions;

      i = GUESS_ALLOCATION;	/* Start count with the number stored */
      ptr = lastptr;		/* One past the last ptr with a result */
#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
      if (this->initptr + 4 < ptr) {
	while (i < this->npositions) {
	  if ((value0 = csa_lookup(sarray,--ptr)) < low_adj) {
	    /* Skip */
	  } else if (value0 > high_adj) {
	    /* Skip */
	  } else {
	    this->positions[i++] = value0 - this->querystart;
	  }
	}

      } else {
#ifdef HAVE_64_BIT
	pointer = (UINT8) &(array[ptr]);
#else
	pointer = (UINT4) &(array[ptr]);
#endif
	n_prealign = ((pointer & 0xF)/4) & 0x3;
	debug7(printf("Initial ptr is at location %p.  Need %d to get to 128-bit boundary\n",
		      &(array[ptr]),n_prealign));

	/* Initial part */
	while (i < this->npositions) {
	  if ((value0 = CONVERT(array[--ptr])) < low_adj) {
	    /* Skip */
	  } else if (value0 > high_adj) {
	    /* Skip */
	  } else {
	    this->positions[i++] = value0 - this->querystart;
	  }
	}

	/* Aligned part */
	while (i < this->npositions && this->initptr + 4 < ptr) {
	  values = _mm_load_si128((__m128i *) &(array[ptr-4]));
	  values = _mm_sub_epi32(values,epi32_convert);
	  compare = _mm_and_si128(_mm_cmpgt_epi32(values,floor),_mm_cmplt_epi32(values,ceiling));
	  if (/*cmp*/_mm_movemask_epi8(compare) == 0x0000) {
	    /* All results are false, indicating no values between low_adj and high_adj (most common case) */
	    ptr -= 4;
	  } else {
	    value0 = CONVERT(array[--ptr]);
	    if (value0 < low_adj) {
	      /* Skip */
	    } else if (value0 > high_adj) {
	      /* Skip */
	    } else {
	      this->positions[i++] = value0 - this->querystart;
	    }

	    value1 = CONVERT(array[--ptr]);
	    if (value1 < low_adj) {
	      /* Skip */
	    } else if (value1 > high_adj) {
	      /* Skip */
	    } else {
	      this->positions[i++] = value1 - this->querystart;
	    }

	    value2 = CONVERT(array[--ptr]);
	    if (value2 < low_adj) {
	      /* Skip */
	    } else if (value2 > high_adj) {
	      /* Skip */
	    } else {
	      this->positions[i++] = value2 - this->querystart;
	    }

	    value3 = CONVERT(array[--ptr]);
	    if (value3 < low_adj) {
	      /* Skip */
	    } else if (value3 > high_adj) {
	      /* Skip */
	    } else {
	      this->positions[i++] = value3 - this->querystart;
	    }
	  }
  	}
	  
	/* Last part */
	while (i < this->npositions) {
	  if ((value0 = csa_lookup(sarray,--ptr)) < low_adj) {
	    /* Skip */
	  } else if (value0 > high_adj) {
	    /* Skip */
	  } else {
	    this->positions[i++] = value0 - this->querystart;
	  }
	}
      }

#else

      while (i < this->npositions) {
	if ((value0 = csa_lookup(sarray,--ptr)) < low_adj) {
	  /* Skip */
	} else if (value0 > high_adj) {
	  /* Skip */
	} else {
	  this->positions[i++] = value0 - this->querystart;
	}
      }
#endif
    }

#ifdef USE_QSORT
    qsort(this->positions,this->npositions,sizeof(Univcoord_T),Univcoord_compare);
#else
    Sedgesort_uint4(this->positions,this->npositions);
#endif
    debug7(printf("Sorting %d positions\n",this->npositions));

#if 0
    /* Not sure why we were doing this.  We will find collinear set of diagonals later. */
    /* Advance pointer to goal (note: do not want goal_adj, since we have already subtracted this->querystart) */
    /* Have tested positions[i] <= goal, but want positions[-1] to be < goal, or positions[0] >= goal */
    i = 0;
    while (i < this->npositions && this->positions[i] < goal) {
      debug7(printf("3 Skipping position %u < goal %u\n",this->positions[i] - chroffset,goal - chroffset));
      i++;
    }
    this->positions += i;
    this->npositions -= i;
    debug7(printf("Remaining: %d positions\n",this->npositions));
#endif
  }

  return;
}
  
#endif


/* ? Returns first entry that is >= goal */
static int
binary_search (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi-1,positions[highi-1],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}



/* Used upon second call to Elt_fill_positions_filtered */
static void
fill_positions_filtered_again (Elt_T this, Univcoord_T low, Univcoord_T high) {
  int lowi, highi, i;


  debug(printf("Entered fill_positions_filtered_again with low %u and high %u, initptr %u and finalptr %u (n = %d), nmatches %d\n",
	       low,high,this->initptr,this->finalptr,this->finalptr - this->initptr + 1,this->nmatches));

  if (this->positions_allocated != NULL) {
    /* Filled from a previous call */
    FREE(this->positions_allocated);
  }

  if (this->n_all_positions == 0) {
    this->positions_allocated = this->positions = (Univcoord_T *) NULL;
    this->npositions_allocated = this->npositions = 0;

  } else {
    /* low_adj and high_adj are inclusive */
    lowi = binary_search(/*lowi*/0,/*highi*/this->n_all_positions,this->all_positions,/*goal*/low + this->querystart);
    highi = binary_search(lowi,/*highi*/this->n_all_positions,this->all_positions,/*goal*/high + this->querystart + 1) - 1;
    if ((this->npositions_allocated = this->npositions = highi - lowi + 1) == 0) {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;

    } else {
      this->positions_allocated = this->positions = (Univcoord_T *) MALLOC(this->npositions * sizeof(Univcoord_T));
      memcpy(this->positions,&(this->all_positions[lowi]),this->npositions*sizeof(Univcoord_T));
      for (i = 0; i < this->npositions; i++) {
	this->positions[i] -= this->querystart;
      }
    }
  }

  return;
}



static void
Elt_fill_positions_filtered (Elt_T this, T sarray, Univcoord_T goal, Univcoord_T low, Univcoord_T high,
			     Compress_T query_compress, bool plusp, int genestrand, bool multiplep) {
  int nmatches;
#ifdef DEBUG8
  Univcoord_T *positions_std;
  int npositions_std;
#endif
#ifdef WORDS_BIGENDIAN
  int i;
#endif


  if (this->nmatches == 0 || this->finalptr - this->initptr + 1 > EXCESS_SARRAY_HITS) {
    /* Check for an extension */
    nmatches = Genome_consecutive_matches_rightward(query_compress,/*left*/goal,/*pos5*/this->querystart,
						    /*pos3*/this->queryend + 1,plusp,genestrand);
    debug7(printf("rightward at goal %u from %d to %d shows %d matches (want %d)\n",goal,this->querystart,this->queryend,
		  nmatches,this->queryend - this->querystart + 1));

    if (this->positions_allocated != NULL) {
      /* Filled from a previous call */
      FREE(this->positions_allocated);
    }

    if (nmatches == this->queryend - this->querystart + 1) {
      /* Create a position that works */
      this->positions_allocated = this->positions = (Univcoord_T *) CALLOC(1,sizeof(Univcoord_T));
      this->positions[0] = goal;
      this->npositions_allocated = this->npositions = 1;
    } else {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
      this->npositions_allocated = this->npositions = 0;
    }
    return;			/* Don't even try other methods */

  } else if (multiplep == true) {
    if (this->status == ELT_VIRGIN) {
      /* Just go directly to sorting method, and skip SIMD filtering method */
      this->status = ELT_UNSORTED;
    }
  }

  if (this->status == ELT_VIRGIN) {
    fill_positions_filtered_first(this,sarray,low,high);
    if (this->finalptr - this->initptr + 1 > EXCESS_SARRAY_HITS) {
      /* Just keep filtering using SIMD method */
      this->all_positions = (Univcoord_T *) NULL;
      this->n_all_positions = 0;
    } else {
      this->status = ELT_UNSORTED;
    }

  } else if (this->status == ELT_UNSORTED) {
    if ((this->n_all_positions = this->finalptr - this->initptr + 1) == 0 /*|| this->npositions > EXCESS_SARRAY_HITS*/) {
      this->all_positions = (Univcoord_T *) NULL;
      this->n_all_positions = 0;
    } else {
#ifdef USE_QSORT
      this->all_positions = (Univcoord_T *) MALLOC(this->n_all_positions*sizeof(Univcoord_T));
#else
      this->all_positions = (Univcoord_T *) MALLOC((this->n_all_positions + 1) * sizeof(Univcoord_T));
#endif
#ifdef WORDS_BIGENDIAN
      for (i = 0; i < this->n_all_positions; i++) {
	this->all_positions[i] = Bigendian_convert_uint(sarray->array[this->initptr+i]);
      }
#else
      memcpy(this->all_positions,&(sarray->array[this->initptr]),this->n_all_positions*sizeof(Univcoord_T));
#endif
#ifdef USE_QSORT
      qsort(this->all_positions,this->n_all_positions,sizeof(Univcoord_T),Univcoord_compare);
#else
      Sedgesort_uint4(this->all_positions,this->n_all_positions);
#endif
    }
#ifdef DEBUG10
    for (i = 0; i < this->n_all_positions; i++) {
      printf("%d: %u\n",i,this->all_positions[i]);
    }
    printf("\n");
#endif

    fill_positions_filtered_again(this,low,high);
    this->status = ELT_SORTED;

  } else {
    /* ELT_SORTED */
    fill_positions_filtered_again(this,low,high);
  }

#ifdef DEBUG8
  positions_std = fill_positions_std(&npositions_std,/*low_adj*/low + this->querystart,
				     /*high_adj*/high + this->querystart,
				     this->initptr,this->finalptr,this->querystart,sarray->array);
  positions_compare(this->positions_allocated,this->npositions_allocated,positions_std,npositions_std);
  FREE(positions_std);
#endif

  return;
}


#if 0
static void
Elt_dump_list (List_T list) {
  List_T p;
  Elt_T elt;
  int maxn = 0, k;

  for (p = list; p != NULL; p = p->rest) {
    elt = (Elt_T) p->first;
    if (elt->npositions > maxn) {
      maxn = elt->npositions;
    }
  }

  for (k = 0; k < maxn /* && k < 100 */; k++) {
    for (p = list; p != NULL; p = p->rest) {
      elt = (Elt_T) p->first;
      if (k >= elt->npositions) {
	printf("\t");
      } else {
	printf("%d..%d:%u\t",elt->querystart,elt->queryend,elt->positions[k]);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}
#endif

#ifdef DEBUG
static void
Elt_dump (Elt_T elt) {
  int k;

  printf("Elt %d..%d (SA %u+%d) with %d positions:\n",
	 elt->querystart,elt->queryend,elt->initptr,elt->finalptr - elt->initptr,elt->npositions);
  for (k = 0; k < elt->npositions; k++) {
    printf("  %u\n",elt->positions[k]);
  }
  printf("\n");

  return;
}
#endif



#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


#if 0
/* Copied to stage1hr.c */
static int
donor_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length_orig(Stage3end_substring_donor(x));
  int y_length = Substring_match_length_orig(Stage3end_substring_donor(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    return 0;
  }
}
#endif


#if 0
/* Copied to stage1hr.c */
static int
acceptor_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length_orig(Stage3end_substring_acceptor(x));
  int y_length = Substring_match_length_orig(Stage3end_substring_acceptor(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    return 0;
  }
}
#endif


/* Also defined in stage1hr.c */
#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))

#if 0
/* Previously called collect_elt_matches */
static bool
solve_twopart (int *found_score, List_T *subs, List_T *indels, List_T *ambiguous, List_T *singlesplicing,
	       int querystart_same, int queryend_same,
	       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
	       Chrpos_T chrlength, Univcoord_T goal, List_T rightward_set, List_T leftward_set,
	       int querylength, Compress_T query_compress,
	       bool plusp, int genestrand, bool first_read_p, int nmisses_allowed) {
  bool twopartp = false;
  List_T set, p;
  Stage3end_T hit, *hitarray;
  Elt_T elt;
  Univcoord_T left, left1, left2, *array;
  Uintlist_T difflist = NULL;	/* Won't work with LARGE_GENOMES */
  int nmismatches, nindels;
  int nsame, ndiff;
  int querystart_diff, queryend_diff, indel_pos;
#if 0
  int nmismatches1, nmismatches2;
#endif

  List_T accepted_hits, rejected_hits;
  List_T spliceends_sense, spliceends_antisense, lowprob;
  List_T donor_hits, acceptor_hits;
  int donor_length, acceptor_length;
  int nhits, nspliceends_sense, nspliceends_antisense, n_good_spliceends;
  int best_nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob;
  Substring_T donor, acceptor;

  Uintlist_T ambcoords;
  Intlist_T amb_knowni, amb_nmismatches;
  Doublelist_T amb_probs;

  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;
  int k, j, i, n;
  bool segmenti_usedp, segmentj_usedp;
  bool foundp;
  int *segmenti_donor_knownpos, *segmentj_acceptor_knownpos,  *segmentj_antidonor_knownpos, *segmenti_antiacceptor_knownpos,
    *segmenti_donor_knowni, *segmentj_acceptor_knowni, *segmentj_antidonor_knowni, *segmenti_antiacceptor_knowni;


  /* Potential success */
  debug7(printf("  successful candidate found\n"));
  if (goal < (Univcoord_T) querylength) {
    debug7(printf("  Goes over beginning of chromosome\n"));
    return false;
  } else if (goal + querylength > chrhigh) {
    debug7(printf("  Goes over end of chromosome\n"));
    return false;
  } else {
    left = goal /* - querylength */;

#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      segmenti_donor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
      segmentj_acceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
      segmentj_antidonor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
      segmenti_antiacceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
      segmenti_donor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
      segmentj_acceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
      segmentj_antidonor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
      segmenti_antiacceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    } else {
      segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
      segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
      segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
      segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
      segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
      segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
      segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
      segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    }
#else
    segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
#endif
  }


  nsame = ndiff = 0;
  querystart_diff = querylength;
  queryend_diff = 0;
  for (set = rightward_set; set /* != NULL */; set = set->rest) {
    elt = (Elt_T) set->first;
    debug7(printf("%d..%d:%u vs %u: ",elt->querystart,elt->queryend,elt->positions[-1],goal));
    /* assert(elt->status != ELT_VIRGIN); */
    if (elt->positions[-1] == goal) {
      debug7(printf("same\n"));
      if (elt->querystart < querystart_same) {
	querystart_same = elt->querystart;
      }
      if (elt->queryend > queryend_same) {
	queryend_same = elt->queryend;
      }
      nsame++;

    } else {
#if 0
      /* Assertion holds because of values for low and high given to Elt_fill_positions_filtered */
      assert(elt->positions[-1] + max_insertionlen + overall_max_distance > goal &&
	     elt->positions[-1] < goal + max_insertionlen + overall_max_distance);
#endif

      debug7(printf("diff (npositions %d)\n",elt->npositions));
      debug7(printf("Pushing position %u\n",elt->positions[-1]));
      difflist = Uintlist_push(difflist,elt->positions[-1]);
      for (i = 0; i < elt->npositions; i++) {
	debug7(printf("Pushing position %u\n",elt->positions[i]));
	difflist = Uintlist_push(difflist,elt->positions[i]);
      }
      if (elt->querystart < querystart_diff) {
	querystart_diff = elt->querystart;
      }
      if (elt->queryend > queryend_diff) {
	queryend_diff = elt->queryend;
      }
      ndiff++;
    }
  }

  for (set = leftward_set; set /* != NULL */; set = set->rest) {
    elt = (Elt_T) set->first;
    debug7(printf("%d..%d:%u vs %u: ",elt->querystart,elt->queryend,elt->positions[-1],goal));
    /* assert(elt->status != ELT_VIRGIN); */
    if (elt->positions[-1] == goal) {
      debug7(printf("same\n"));
      if (elt->querystart < querystart_same) {
	querystart_same = elt->querystart;
      }
      if (elt->queryend > queryend_same) {
	queryend_same = elt->queryend;
      }
      nsame++;

    } else {
#if 0
      /* Assertion holds because of values for low and high given to Elt_fill_positions_filtered */
      assert(elt->positions[-1] + max_insertionlen + overall_max_distance > goal &&
	     elt->positions[-1] < goal + max_insertionlen + overall_max_distance);
#endif

      debug7(printf("diff (npositions %d)\n",elt->npositions));
      debug7(printf("Pushing position %u\n",elt->positions[-1]));
      difflist = Uintlist_push(difflist,elt->positions[-1]);
      for (i = 0; i < elt->npositions; i++) {
	debug7(printf("Pushing position %u\n",elt->positions[i]));
	difflist = Uintlist_push(difflist,elt->positions[i]);
      }
      if (elt->querystart < querystart_diff) {
	querystart_diff = elt->querystart;
      }
      if (elt->queryend > queryend_diff) {
	queryend_diff = elt->queryend;
      }
      ndiff++;
    }
  }

  debug7(printf("Got %d same, %d diff\n",nsame,ndiff));

  if (ndiff == 0) {
    /* sub */
    debug7(printf("  Testing in entire query\n"));
    nmismatches = Genome_count_mismatches_substring(query_compress,left,/*pos5*/0,/*pos3*/querylength,
						    plusp,genestrand);
    debug7(printf("nmismatches = %d (vs %d misses allowed)\n",nmismatches,nmisses_allowed));

    if (nmismatches > nmisses_allowed) {
      debug7(printf("Result: too many mismatches\n"));

    } else {
      debug7(printf("Result: successful hit saved\n"));
      if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					    left,/*genomiclength*/querylength,
					    query_compress,plusp,genestrand,first_read_p,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*sarrayp*/true)) != NULL) {
	debug1(printf("1. Reporting hit with %d mismatches vs %d allowed\n",nmismatches,nmisses_allowed));
	*subs = List_push(*subs,(void *) hit);
	twopartp = true;
      }
    }
    assert(difflist == NULL);

  } else if (querystart_same == 0 && queryend_diff == querylength - 1) {
    left1 = left;
    indel_pos = queryend_same + 1;
    debug7(printf("same is at %u from %d to %d\n",left,querystart_same,queryend_same));

    n = Uintlist_length(difflist);
#ifdef HAVE_ALLOCA
#ifdef USE_QSORT
    array = (UINT4 *) MALLOCA(n * sizeof(UINT4));
#else
    array = (UINT4 *) MALLOCA((n + 1) * sizeof(UINT4));
#endif
#else
#ifdef USE_QSORT
    array = (UINT4 *) MALLOC(n * sizeof(UINT4));
#else
    array = (UINT4 *) MALLOC((n + 1) * sizeof(UINT4));
#endif
#endif

    Uintlist_fill_array_and_free(array,&difflist);
#ifdef USE_QSORT
    qsort(array,n,sizeof(Univcoord_T),Univcoord_compare);
#else
    Sedgesort(array,n);
#endif
    debug7(printf("Have %d matching diffs\n",n));

    spliceends_sense = spliceends_antisense = (List_T) NULL;
    lowprob = (List_T) NULL;
    for (i = 0; i < n; i++) {
      left2 = array[i];
      debug7(printf("diff %d/%d is at %u, from %d to %d\n",i,n,left2,querystart_diff - 1,queryend_diff));

      if (i > 0 && left2 == array[i-1]) {
	/* Already processed */

      } else if (left2 + querylength >= chrhigh) {
	/* Splice or deletion would extend to next chromosome */

      } else if (left2 > left1 + max_deletionlen) {
	debug7(printf("A splice..."));

	segmenti_donor_nknown = segmenti_antiacceptor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left1,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left1);
	  while (j < nsplicesites && splicesites[j] < left1 + querylength) {
	    if (splicetypes[j] == DONOR) {
	      debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - left1;
	      segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - left1;
	      segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmenti_donor_knownpos[segmenti_donor_nknown] = querylength + 100;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength + 100;
	  
	segmentj_acceptor_nknown = segmentj_antidonor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left2,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left2);
	  while (j < nsplicesites && splicesites[j] < left2 + querylength) {
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - left2;
	      segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - left2;
	      segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength + 100;
	segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength + 100;

	/* nspliceends = 0; */
	assert(segmenti_donor_knownpos[0] == querylength);
	assert(segmentj_acceptor_knownpos[0] == querylength);
	assert(segmentj_antidonor_knownpos[0] == querylength);
	assert(segmenti_antiacceptor_knownpos[0] == querylength);

	spliceends_sense =
	  Splice_solve_single_sense(&(*found_score),&nspliceends_sense,spliceends_sense,&lowprob,
				    &segmenti_usedp,&segmentj_usedp,
				    /*segmenti_left*/left1,/*segmentj_left*/left2,
				    chrnum,chroffset,chrhigh,chrlength,
				    chrnum,chroffset,chrhigh,chrlength,
				    querylength,query_compress,
				    segmenti_donor_knownpos,segmentj_acceptor_knownpos,
				    segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
				    segmenti_donor_knowni,segmentj_acceptor_knowni,
				    segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
				    segmenti_donor_nknown,segmentj_acceptor_nknown,
				    segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
				    splicing_penalty,/*max_mismatches_allowed*/1000,
				    plusp,genestrand,first_read_p,/*subs_or_indels_p*/false,
				    /*sarrayp*/true);

	assert(segmenti_donor_knownpos[0] == querylength);
	assert(segmentj_acceptor_knownpos[0] == querylength);
	assert(segmentj_antidonor_knownpos[0] == querylength);
	assert(segmenti_antiacceptor_knownpos[0] == querylength);

	spliceends_antisense =
	  Splice_solve_single_antisense(&(*found_score),&nspliceends_antisense,spliceends_antisense,&lowprob,
					&segmenti_usedp,&segmentj_usedp,
					/*segmenti_left*/left1,/*segmentj_left*/left2,
					chrnum,chroffset,chrhigh,chrlength,
					chrnum,chroffset,chrhigh,chrlength,
					querylength,query_compress,
					segmenti_donor_knownpos,segmentj_acceptor_knownpos,
					segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					segmenti_donor_knowni,segmentj_acceptor_knowni,
					segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
					segmenti_donor_nknown,segmentj_acceptor_nknown,
					segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
					splicing_penalty,/*max_mismatches_allowed*/1000,
					plusp,genestrand,first_read_p,/*subs_or_indels_p*/false,
					/*sarrayp*/true);

      } else if (left2 > left1) {
	nindels = left2 - left1;
	debug7(printf("B deletion of %d bp relative to max_deletionlen %d (nmisses allowed %d)...",
		      nindels,max_deletionlen,nmisses_allowed));
	if ((indel_pos < 17 || querylength - indel_pos < 17) && nindels > max_end_deletions) {
	  /* Allow regular GSNAP algorithm to find this */
	  debug7(printf("too long for end deletion"));
	} else {
#if 0
	  nmismatches1 = Genome_count_mismatches_substring(query_compress,left1,/*pos5*/0,/*pos3*/indel_pos,
							   plusp,genestrand);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos,
							   /*pos3*/querylength,plusp,genestrand);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos;
	  }
	  /* genomiclength = querylength+nindels; */
	  if ((hit = Stage3end_new_deletion(&(*found_score),nindels,query_indel_pos,
					    nmismatches1,nmismatches2,left1,
					    query_compress,querylength,plusp,genestrand,first_read_p,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	    twopartp = true;
	  }
#else
	  *indels = Indel_solve_middle_deletion(&foundp,&(*found_score),&nhits,*indels,
						/*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						/*indels*/-nindels,query_compress,querylength,
						nmisses_allowed,plusp,genestrand,/*sarray*/true);
	  debug7(
		 if (foundp == true) {
		   printf("successful");
		 }
		 );
#endif
	}
	debug7(printf("\n"));
      
      } else if (left2 < left1) {
	nindels = left1 - left2;
	if (nindels >= indel_pos || indel_pos + nindels >= querylength) {
	  debug7(printf("X insertion of %d bp too long\n",nindels));
	} else {
	  debug7(printf("C insertion of %d bp (nmisses allowed %d)...",nindels,nmisses_allowed));
#if 0
	  nmismatches1 = Genome_count_mismatches_substring(query_compress,left1,/*pos5*/0,/*pos3*/indel_pos-nindels,
							   plusp,genestrand);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos+nindels,
							   /*pos3*/querylength,plusp,genestrand);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos - nindels;
	  }
	  /* genomiclength = querylength-nindels; */
	  if ((hit = Stage3end_new_insertion(&(*found_score),nindels,query_indel_pos,
					     nmismatches1,nmismatches2,left1,
					     query_compress,querylength,plusp,genestrand,first_read_p,
					     chrnum,chroffset,chrhigh,chrlength,
					     /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	    twopartp = true;
	  }
#else
	  *indels = Indel_solve_middle_insertion(&foundp,&(*found_score),&nhits,*indels,
						 /*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						 /*indels*/+nindels,query_compress,querylength,nmisses_allowed,
						 plusp,genestrand,/*sarrayp*/true);
	  debug7(
		 if (foundp == true) {
		   printf("successful");
		 }
		 );
#endif
	  debug7(printf("\n"));
	}
      }
    }

    if (spliceends_sense != NULL) {
      /* nmismatches here may be different for spliceends from Splice_solve, so pick based on prob and nmismatches */
      best_nmismatches = querylength;
      best_prob = 0.0;
      for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
		      Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
		      Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
		      Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
		      Stage3end_nmismatches_whole(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
		      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	}
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      accepted_hits = rejected_hits = (List_T) NULL;
      for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
	    Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	  accepted_hits = List_push(accepted_hits,(void *) hit);
	} else {
	  rejected_hits = List_push(rejected_hits,(void *) hit);
	}
      }

      if (n_good_spliceends == 0) {
	/* Conjunction is too strict.  Allow for disjunction instead. */
	List_free(&rejected_hits);
	for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
	      Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	    debug7(printf("accepting distance %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			  Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	    n_good_spliceends += 1;
	    accepted_hits = List_push(accepted_hits,(void *) hit);
	  } else {
	    rejected_hits = List_push(rejected_hits,(void *) hit);
	  }
	}
      }

      for (p = rejected_hits; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	Stage3end_free(&hit);
      }
      List_free(&rejected_hits);
      List_free(&spliceends_sense);

      if (n_good_spliceends == 1) {
	*singlesplicing = List_push(*singlesplicing,List_head(accepted_hits));
	nhits += 1;
	List_free(&accepted_hits);

      } else {
	/* 1.  Multiple hits, sense, left1 */
	debug7(printf("multiple hits with best prob, sense\n"));
	donor_hits = acceptor_hits = (List_T) NULL;
	if (plusp == true) {
	  for (p = accepted_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
	    acceptor = Stage3end_substring_acceptor(hit);
	    if (Substring_genomicstart(donor) == left1) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicstart(acceptor) == left1) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	} else {
	  for (p = accepted_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
	    acceptor = Stage3end_substring_acceptor(hit);
	    if (Substring_genomicend(donor) == left1) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicend(acceptor) == left1) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	}

	if (donor_hits != NULL) {
	  hitarray = (Stage3end_T *) List_to_array_n(&n,donor_hits);
	  qsort(hitarray,n,sizeof(Stage3end_T),donor_match_length_cmp);
	  i = 0;
	  while (i < n) {
	    hit = hitarray[i];
	    donor = Stage3end_substring_donor(hit);
	    donor_length = Substring_match_length_orig(donor);
	    j = i + 1;
	    while (j < n && Substring_match_length_orig(Stage3end_substring_donor(hitarray[j])) == donor_length) {
	      j++;
	    }
	    if (j == i + 1) {
	      *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    } else {
	      ambcoords = (Uintlist_T) NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (k = i; k < j; k++) {
		acceptor = Stage3end_substring_acceptor(hitarray[k]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_A(acceptor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_A(acceptor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
	      }

	      nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	      prob = best_prob - Substring_siteD_prob(donor);
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								   donor,/*acceptor*/NULL,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*ambcoords_donor*/NULL,ambcoords,
								   /*amb_knowni_donor*/NULL,amb_knowni,
								   /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								   /*amb_probs_donor*/NULL,amb_probs,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								   Stage3end_sensedir(hit),/*sarrayp*/true));
	      twopartp = true;
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_nmismatches);
	      Intlist_free(&amb_knowni);
	      Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	      for (k = i; k < j; k++) {
		hit = hitarray[k];
		Stage3end_free(&hit);
	      }
	    }

	    i = j;
	  }
	  FREE(hitarray);
	  List_free(&donor_hits);
	}

	if (acceptor_hits != NULL) {
	  hitarray = (Stage3end_T *) List_to_array_n(&n,acceptor_hits);
	  qsort(hitarray,n,sizeof(Stage3end_T),acceptor_match_length_cmp);
	  i = 0;
	  while (i < n) {
	    hit = hitarray[i];
	    acceptor = Stage3end_substring_acceptor(hit);
	    acceptor_length = Substring_match_length_orig(acceptor);
	    j = i + 1;
	    while (j < n && Substring_match_length_orig(Stage3end_substring_acceptor(hitarray[j])) == acceptor_length) {
	      j++;
	    }
	    if (j == i + 1) {
	      *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    } else {
	      ambcoords = (Uintlist_T) NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (k = i; k < j; k++) {
		donor = Stage3end_substring_donor(hitarray[k]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_D(donor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_D(donor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
	      }
	    
	      nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	      prob = best_prob - Substring_siteA_prob(acceptor);
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								   /*donor*/NULL,acceptor,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   ambcoords,/*ambcoords_acceptor*/NULL,
								   amb_knowni,/*amb_knowni_acceptor*/NULL,
								   amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								   amb_probs,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								   Stage3end_sensedir(hit),/*sarrayp*/true));
	      twopartp = true;
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_nmismatches);
	      Intlist_free(&amb_knowni);
	      Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	      for (k = i; k < j; k++) {
		hit = hitarray[k];
		Stage3end_free(&hit);
	      }
	    }

	    i = j;
	  }
	  FREE(hitarray);
	  List_free(&acceptor_hits);
	}

	List_free(&accepted_hits);
      }
    }

    if (spliceends_antisense != NULL) {
      /* nmismatches here may be different for spliceends from Splice_solve, so pick based on prob and nmismatches */
      best_nmismatches = querylength;
      best_prob = 0.0;
      for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
		      Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
		      Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
		      Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
		      Stage3end_nmismatches_whole(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
		      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	}
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      accepted_hits = rejected_hits = (List_T) NULL;
      for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
	    Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	  debug7(printf("accepting distance %d, donor length %d and acceptor length %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			Substring_siteD_prob(Stage3end_substring_donor(hit)),
			Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	  accepted_hits = List_push(accepted_hits,(void *) hit);
	} else {
	  rejected_hits = List_push(rejected_hits,(void *) hit);
	}
      }

      if (n_good_spliceends == 0) {
	/* Conjunction is too strict.  Allow for disjunction instead. */
	List_free(&rejected_hits);
	for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
	      Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	    debug7(printf("accepting distance %d, donor length %d and acceptor length %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			  Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			  Substring_siteD_prob(Stage3end_substring_donor(hit)),
			  Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	    n_good_spliceends += 1;
	    accepted_hits = List_push(accepted_hits,(void *) hit);
	  } else {
	    rejected_hits = List_push(rejected_hits,(void *) hit);
	  }
	}
      }

      for (p = rejected_hits; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	Stage3end_free(&hit);
      }
      List_free(&rejected_hits);
      List_free(&spliceends_antisense);

      if (n_good_spliceends == 1) {
	*singlesplicing = List_push(*singlesplicing,List_head(accepted_hits));
	nhits += 1;
	List_free(&accepted_hits);

      } else {
	/* 2.  Multiple hits, antisense, left1 */
	debug7(printf("multiple hits with best prob, antisense\n"));
	donor_hits = acceptor_hits = (List_T) NULL;
	if (plusp == true) {
	  for (p = accepted_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
	    acceptor = Stage3end_substring_acceptor(hit);
	    if (Substring_genomicstart(donor) == left1) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicstart(acceptor) == left1) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	} else {
	  for (p = accepted_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
	    acceptor = Stage3end_substring_acceptor(hit);
	    if (Substring_genomicend(donor) == left1) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicend(acceptor) == left1) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	}

	if (donor_hits != NULL) {
	  hitarray = (Stage3end_T *) List_to_array_n(&n,donor_hits);
	  qsort(hitarray,n,sizeof(Stage3end_T),donor_match_length_cmp);
	  i = 0;
	  while (i < n) {
	    hit = hitarray[i];
	    donor = Stage3end_substring_donor(hit);
	    donor_length = Substring_match_length_orig(donor);
	    j = i + 1;
	    while (j < n && Substring_match_length_orig(Stage3end_substring_donor(hitarray[j])) == donor_length) {
	      j++;
	    }
	    if (j == i + 1) {
	      *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    } else {
	      ambcoords = (Uintlist_T) NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;
	      
	      for (k = i; k < j; k++) {
		acceptor = Stage3end_substring_acceptor(hitarray[k]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_A(acceptor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_A(acceptor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
	      }
	      
	      nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	      prob = best_prob - Substring_siteD_prob(donor);
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								   donor,/*acceptor*/NULL,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*ambcoords_donor*/NULL,ambcoords,
								   /*amb_knowni_donor*/NULL,amb_knowni,
								   /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								   /*amb_probs_donor*/NULL,amb_probs,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								   Stage3end_sensedir(hit),/*sarrayp*/true));
	      twopartp = true;
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_nmismatches);
	      Intlist_free(&amb_knowni);
	      Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	      for (k = i; k < j; k++) {
		hit = hitarray[k];
		Stage3end_free(&hit);
	      }
	    }

	    i = j;
	  }
	  FREE(hitarray);
	  List_free(&donor_hits);
	}

	if (acceptor_hits != NULL) {
	  hitarray = (Stage3end_T *) List_to_array_n(&n,acceptor_hits);
	  qsort(hitarray,n,sizeof(Stage3end_T),acceptor_match_length_cmp);
	  i = 0;
	  while (i < n) {
	    hit = hitarray[i];
	    acceptor = Stage3end_substring_acceptor(hit);
	    acceptor_length = Substring_match_length_orig(acceptor);
	    j = i + 1;
	    while (j < n && Substring_match_length_orig(Stage3end_substring_acceptor(hitarray[j])) == acceptor_length) {
	      j++;
	    }
	    if (j == i + 1) {
	      *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    } else {
	      ambcoords = (Uintlist_T) NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (k = i; k < j; k++) {
		donor = Stage3end_substring_donor(hitarray[k]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_D(donor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_D(donor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
	      }
	    
	      nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	      prob = best_prob - Substring_siteA_prob(acceptor);
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								   /*donor*/NULL,acceptor,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   ambcoords,/*ambcoords_acceptor*/NULL,
								   amb_knowni,/*amb_knowni_acceptor*/NULL,
								   amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								   amb_probs,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								   Stage3end_sensedir(hit),/*sarrayp*/true));
	      twopartp = true;
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_nmismatches);
	      Intlist_free(&amb_knowni);
	      Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	      for (k = i; k < j; k++) {
		hit = hitarray[k];
		Stage3end_free(&hit);
	      }
	    }

	    i = j;
	  }
	  FREE(hitarray);
	  List_free(&acceptor_hits);
	}

	List_free(&accepted_hits);
      }
    }

    /* Don't use lowprob in suffix array stage */
    debug7(printf("freeing lowprobs\n"));
    for (p = lowprob; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      Stage3end_free(&hit);
    }
    List_free(&lowprob);

#ifdef HAVE_ALLOCA
    FREEA(array);
#else
    FREE(array);
#endif

  } else if (querystart_diff == 0 && queryend_same == querylength - 1) {
    left2 = left;
    indel_pos = querystart_same;
    debug7(printf("same is at %u from %d to %d\n",left,querystart_same,queryend_same));
    
    n = Uintlist_length(difflist);
#ifdef HAVE_ALLOCA
#ifdef USE_QSORT
    array = (UINT4 *) MALLOCA(n * sizeof(UINT4));
#else
    array = (UINT4 *) MALLOCA((n + 1) * sizeof(UINT4));
#endif
#else
#ifdef USE_QSORT
    array = (UINT4 *) MALLOC(n * sizeof(UINT4));
#else
    array = (UINT4 *) MALLOC((n + 1) * sizeof(UINT4));
#endif
#endif

    Uintlist_fill_array_and_free(array,&difflist);
#ifdef USE_QSORT
    qsort(array,n,sizeof(Univcoord_T),Univcoord_compare);
#else
    Sedgesort_uint4(array,n);
#endif
    debug7(printf("Have %d matching diffs\n",n));

    spliceends_sense = spliceends_antisense = (List_T) NULL;
    lowprob = (List_T) NULL;
    for (i = 0; i < n; i++) {
      left1 = array[i];
      debug7(printf("diff %d/%d is at %u, from %d to %d\n",i,n,left1,querystart_diff,queryend_diff));

      if (i > 0 && left1 == array[i-1]) {
	/* Already processed */

      } else if (left2 + querylength >= chrhigh) {
	/* Splice or deletion would extend to next chromosome */

      } else if (left2 > left1 + max_deletionlen) {
	debug7(printf("A splice..."));

	segmenti_donor_nknown = segmenti_antiacceptor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left1,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left1);
	  while (j < nsplicesites && splicesites[j] < left1 + querylength) {
	    if (splicetypes[j] == DONOR) {
	      debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - left1;
	      segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - left1;
	      segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmenti_donor_knownpos[segmenti_donor_nknown] = querylength + 100;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength + 100;
	  
	segmentj_acceptor_nknown = segmentj_antidonor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left2,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left2);
	  while (j < nsplicesites && splicesites[j] < left2 + querylength) {
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - left2;
	      segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - left2;
	      segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength + 100;
	segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength + 100;

	/* nspliceends = 0; */
	spliceends_sense =
	  Splice_solve_single_sense(&(*found_score),&nspliceends_sense,spliceends_sense,&lowprob,
				    &segmenti_usedp,&segmentj_usedp,
				    /*segmenti_left*/left1,/*segmentj_left*/left2,
				    chrnum,chroffset,chrhigh,chrlength,
				    chrnum,chroffset,chrhigh,chrlength,
				    querylength,query_compress,
				    segmenti_donor_knownpos,segmentj_acceptor_knownpos,
				    segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
				    segmenti_donor_knowni,segmentj_acceptor_knowni,
				    segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
				    segmenti_donor_nknown,segmentj_acceptor_nknown,
				    segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
				    splicing_penalty,/*max_mismatches_allowed*/1000,
				    plusp,genestrand,first_read_p,/*subs_or_indels_p*/false,
				    /*sarrayp*/true);
	spliceends_antisense =
	  Splice_solve_single_antisense(&(*found_score),&nspliceends_antisense,spliceends_antisense,&lowprob,
				    &segmenti_usedp,&segmentj_usedp,
				    /*segmenti_left*/left1,/*segmentj_left*/left2,
				    chrnum,chroffset,chrhigh,chrlength,
				    chrnum,chroffset,chrhigh,chrlength,
				    querylength,query_compress,
				    segmenti_donor_knownpos,segmentj_acceptor_knownpos,
				    segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
				    segmenti_donor_knowni,segmentj_acceptor_knowni,
				    segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
				    segmenti_donor_nknown,segmentj_acceptor_nknown,
				    segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
				    splicing_penalty,/*max_mismatches_allowed*/1000,
				    plusp,genestrand,first_read_p,/*subs_or_indels_p*/false,
				    /*sarrayp*/true);

      } else if (left2 > left1) {
	nindels = left2 - left1;
	debug7(printf("B deletion of %d bp relative to max_deletionlen %d (nmisses allowed %d)...",
		      nindels,max_deletionlen,nmisses_allowed));
	if ((indel_pos < 17 || querylength - indel_pos < 17) && nindels > max_end_deletions) {
	  /* Allow regular GSNAP algorithm to find this */
	  debug7(printf("too long for end deletion"));
	} else {
#if 0
	  nmismatches1 = Genome_count_mismatches_substring(query_compress,left1,/*pos5*/0,/*pos3*/indel_pos,
							   plusp,genestrand);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos,
							   /*pos3*/querylength,plusp,genestrand);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos;
	  }
	  /* genomiclength = querylength+nindels; */
	  if ((hit = Stage3end_new_deletion(&(*found_score),nindels,query_indel_pos,
					    nmismatches1,nmismatches2,left1,
					    query_compress,querylength,plusp,genestrand,first_read_p,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	    twopartp = true;
	  }
#else
	  *indels = Indel_solve_middle_deletion(&foundp,&(*found_score),&nhits,*indels,
						/*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						/*indels*/-nindels,query_compress,querylength,
						nmisses_allowed,plusp,genestrand,/*sarray*/true);
	  debug7(
		 if (foundp == true) {
		   printf("successful");
		 }
		 );
#endif
	}
	debug7(printf("\n"));
      
      } else if (left2 < left1) {
	nindels = left1 - left2;
	if (nindels >= indel_pos || indel_pos + nindels >= querylength) {
	  debug7(printf("X insertion of %d bp too long\n",nindels));
	} else {
	  debug7(printf("C insertion of %d bp (nmisses allowed %d)...",nindels,nmisses_allowed));
#if 0      
	  nmismatches1 = Genome_count_mismatches_substring(query_compress,left1,/*pos5*/0,/*pos3*/indel_pos-nindels,
							   plusp,genestrand);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos+nindels,
							   /*pos3*/querylength,plusp,genestrand);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos - nindels;
	  }
	  /*genomiclength = querylength-nindels; */
	  if ((hit = Stage3end_new_insertion(&(*found_score),nindels,query_indel_pos,
					     nmismatches1,nmismatches2,left1,
					     query_compress,querylength,plusp,genestrand,first_read_p,
					     chrnum,chroffset,chrhigh,chrlength,
					     /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	    twopartp = true;
	  }
#else
	  *indels = Indel_solve_middle_insertion(&foundp,&(*found_score),&nhits,*indels,
						 /*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						 /*indels*/+nindels,query_compress,querylength,nmisses_allowed,
						 plusp,genestrand,/*sarrayp*/true);
	  debug7(
		 if (foundp == true) {
		   printf("successful");
		 }
		 );
#endif
	  debug7(printf("\n"));
	}
      }
    }

    if (spliceends_sense != NULL) {
      /* nmismatches here may be different for spliceends from Splice_solve, so pick based on prob and nmismatches */
      best_nmismatches = querylength;
      best_prob = 0.0;
      for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
		      Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
		      Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
		      Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
		      Stage3end_nmismatches_whole(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
		      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	}
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      accepted_hits = rejected_hits = (List_T) NULL;
      for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
	    Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	  accepted_hits = List_push(accepted_hits,(void *) hit);
	} else {
	  rejected_hits = List_push(rejected_hits,(void *) hit);
	}
      }
      
      if (n_good_spliceends == 0) {
	/* Conjunction is too strict.  Allow for disjunction instead. */
	List_free(&rejected_hits);
	for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
	      Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	    debug7(printf("accepting distance %d, donor length %d and acceptor length %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			  Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			  Substring_siteD_prob(Stage3end_substring_donor(hit)),
			  Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	    n_good_spliceends += 1;
	    accepted_hits = List_push(accepted_hits,(void *) hit);
	  } else {
	    rejected_hits = List_push(rejected_hits,(void *) hit);
	  }
	}
      }

      for (p = rejected_hits; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	Stage3end_free(&hit);
      }
      List_free(&rejected_hits);
      List_free(&spliceends_sense);

      if (n_good_spliceends == 1) {
	*singlesplicing = List_push(*singlesplicing,List_head(accepted_hits));
	nhits += 1;
	List_free(&accepted_hits);

      } else {
	/* 3.  Multiple hits, sense, left2 */
	debug7(printf("multiple hits with best prob, sense\n"));
	donor_hits = acceptor_hits = (List_T) NULL;
	if (plusp == true) {
	  for (p = accepted_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
	    acceptor = Stage3end_substring_acceptor(hit);
	    if (Substring_genomicstart(donor) == left2) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicstart(acceptor) == left1) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	} else {
	  for (p = accepted_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
	    acceptor = Stage3end_substring_acceptor(hit);
	    if (Substring_genomicend(donor) == left2) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicend(acceptor) == left1) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	}
	  
	if (donor_hits != NULL) {
	  hitarray = (Stage3end_T *) List_to_array_n(&n,donor_hits);
	  qsort(hitarray,n,sizeof(Stage3end_T),donor_match_length_cmp);
	  i = 0;
	  while (i < n) {
	    hit = hitarray[i];
	    donor = Stage3end_substring_donor(hit);
	    donor_length = Substring_match_length_orig(donor);
	    j = i + 1;
	    while (j < n && Substring_match_length_orig(Stage3end_substring_donor(hitarray[j])) == donor_length) {
	      j++;
	    }
	    if (j == i + 1) {
	      *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    } else {
	      ambcoords = (Uintlist_T) NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (k = i; k < j; k++) {
		acceptor = Stage3end_substring_acceptor(hitarray[k]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_A(acceptor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_A(acceptor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
	      }

	      nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	      prob = best_prob - Substring_siteD_prob(donor);
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								   donor,/*acceptor*/NULL,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*ambcoords_donor*/NULL,ambcoords,
								   /*amb_knowni_donor*/NULL,amb_knowni,
								   /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								   /*amb_probs_donor*/NULL,amb_probs,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								   Stage3end_sensedir(hit),/*sarrayp*/true));
	      twopartp = true;
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_nmismatches);
	      Intlist_free(&amb_knowni);
	      Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */


	      for (k = i; k < j; k++) {
		hit = hitarray[k];
		Stage3end_free(&hit);
	      }
	    }

	    i = j;
	  }
	  FREE(hitarray);
	  List_free(&donor_hits);
	}

	if (acceptor_hits != NULL) {
	  hitarray = (Stage3end_T *) List_to_array_n(&n,acceptor_hits);
	  qsort(hitarray,n,sizeof(Stage3end_T),acceptor_match_length_cmp);
	  i = 0;
	  while (i < n) {
	    hit = hitarray[i];
	    acceptor = Stage3end_substring_acceptor(hit);
	    acceptor_length = Substring_match_length_orig(acceptor);
	    j = i + 1;
	    while (j < n && Substring_match_length_orig(Stage3end_substring_acceptor(hitarray[j])) == acceptor_length) {
	      j++;
	    }
	    if (j == i + 1) {
	      *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    } else {
	      ambcoords = (Uintlist_T) NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (k = i; k < j; k++) {
		donor = Stage3end_substring_donor(hitarray[k]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_D(donor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_D(donor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
	      }

	      nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	      prob = best_prob - Substring_siteA_prob(acceptor);
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								   /*donor*/NULL,acceptor,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   ambcoords,/*ambcoords_acceptor*/NULL,
								   amb_knowni,/*amb_knowni_acceptor*/NULL,
								   amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								   amb_probs,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								   Stage3end_sensedir(hit),/*sarrayp*/true));
	      twopartp = true;
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_nmismatches);
	      Intlist_free(&amb_knowni);
	      Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	      for (k = i; k < j; k++) {
		hit = hitarray[k];
		Stage3end_free(&hit);
	      }
	    }

	    i = j;
	  }
	  FREE(hitarray);
	  List_free(&acceptor_hits);
	}

	List_free(&accepted_hits);
      }
    }

    if (spliceends_antisense != NULL) {
      /* nmismatches here may be different for spliceends from Splice_solve, so pick based on prob and nmismatches */
      best_nmismatches = querylength;
      best_prob = 0.0;
      for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
		      Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
		      Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
		      Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
		      Stage3end_nmismatches_whole(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
		      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	}
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      accepted_hits = rejected_hits = (List_T) NULL;
      for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
	    Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	  accepted_hits = List_push(accepted_hits,(void *) hit);
	} else {
	  rejected_hits = List_push(rejected_hits,(void *) hit);
	}
      }
      
      if (n_good_spliceends == 0) {
	/* Conjunction is too strict.  Allow for disjunction instead. */
	List_free(&rejected_hits);
	for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
	      Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	    debug7(printf("accepting distance %d, donor length %d and acceptor length %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			  Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			  Substring_siteD_prob(Stage3end_substring_donor(hit)),
			  Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	    n_good_spliceends += 1;
	    accepted_hits = List_push(accepted_hits,(void *) hit);
	  } else {
	    rejected_hits = List_push(rejected_hits,(void *) hit);
	  }
	}
      }

      for (p = rejected_hits; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	Stage3end_free(&hit);
      }
      List_free(&rejected_hits);
      List_free(&spliceends_antisense);

      if (n_good_spliceends == 1) {
	*singlesplicing = List_push(*singlesplicing,List_head(accepted_hits));
	nhits += 1;
	List_free(&accepted_hits);

      } else {
	/* 4.  Multiple hits, antisense, left2 */
	debug7(printf("multiple hits with best prob, antisense\n"));
	donor_hits = acceptor_hits = (List_T) NULL;
	if (plusp == true) {
	  for (p = accepted_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
	    acceptor = Stage3end_substring_acceptor(hit);
	    if (Substring_genomicstart(donor) == left2) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicstart(acceptor) == left2) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	} else {
	  for (p = accepted_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
	    acceptor = Stage3end_substring_acceptor(hit);
	    if (Substring_genomicend(donor) == left2) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicend(acceptor) == left2) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	}

	if (donor_hits != NULL) {
	  hitarray = (Stage3end_T *) List_to_array_n(&n,donor_hits);
	  qsort(hitarray,n,sizeof(Stage3end_T),donor_match_length_cmp);
	  i = 0;
	  while (i < n) {
	    hit = hitarray[i];
	    donor = Stage3end_substring_donor(hit);
	    donor_length = Substring_match_length_orig(donor);
	    j = i + 1;
	    while (j < n && Substring_match_length_orig(Stage3end_substring_donor(hitarray[j])) == donor_length) {
	      j++;
	    }
	    if (j == i + 1) {
	      *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    } else {
	      ambcoords = (Uintlist_T) NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (k = i; k < j; k++) {
		acceptor = Stage3end_substring_acceptor(hitarray[k]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_A(acceptor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_A(acceptor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
	      }

	      nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	      prob = best_prob - Substring_siteD_prob(donor);
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								   donor,/*acceptor*/NULL,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*ambcoords_donor*/NULL,ambcoords,
								   /*amb_knowni_donor*/NULL,amb_knowni,
								   /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								   /*amb_probs_donor*/NULL,amb_probs,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								   Stage3end_sensedir(hit),/*sarrayp*/true));
	      twopartp = true;
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_nmismatches);
	      Intlist_free(&amb_knowni);
	      Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	      for (k = i; k < j; k++) {
		hit = hitarray[k];
		Stage3end_free(&hit);
	      }
	    }

	    i = j;
	  }
	  FREE(hitarray);
	  List_free(&donor_hits);
	}

	if (acceptor_hits != NULL) {
	  hitarray = (Stage3end_T *) List_to_array_n(&n,acceptor_hits);
	  qsort(hitarray,n,sizeof(Stage3end_T),acceptor_match_length_cmp);
	  i = 0;
	  while (i < n) {
	    hit = hitarray[i];
	    acceptor = Stage3end_substring_acceptor(hit);
	    acceptor_length = Substring_match_length_orig(acceptor);
	    j = i + 1;
	    while (j < n && Substring_match_length_orig(Stage3end_substring_acceptor(hitarray[j])) == acceptor_length) {
	      j++;
	    }
	    if (j == i + 1) {
	      *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    } else {
	      ambcoords = (Uintlist_T) NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (k = i; k < j; k++) {
		donor = Stage3end_substring_donor(hitarray[k]);
#ifdef LARGE_GENOMES
		ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_D(donor));
#else
		ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_D(donor));
#endif
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
	      }

	      nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	      prob = best_prob - Substring_siteA_prob(acceptor);
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								   /*donor*/NULL,acceptor,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   ambcoords,/*ambcoords_acceptor*/NULL,
								   amb_knowni,/*amb_knowni_acceptor*/NULL,
								   amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								   amb_probs,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								   Stage3end_sensedir(hit),/*sarrayp*/true));
	      twopartp = true;
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_nmismatches);
	      Intlist_free(&amb_knowni);
	      Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	      for (k = i; k < j; k++) {
		hit = hitarray[k];
		Stage3end_free(&hit);
	      }
	    }

	    i = j;
	  }
	  FREE(hitarray);
	  List_free(&acceptor_hits);
	}

	List_free(&accepted_hits);
      }
    }


    /* Don't use lowprob in suffix array stage */
    debug7(printf("freeing lowprobs\n"));
    for (p = lowprob; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      Stage3end_free(&hit);
    }
    List_free(&lowprob);

#ifdef HAVE_ALLOCA
    FREEA(array);
#else
    FREE(array);
#endif

  } else {
    Uintlist_free(&difflist);
  }


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(segmenti_donor_knownpos);
    FREEA(segmentj_acceptor_knownpos);
    FREEA(segmentj_antidonor_knownpos);
    FREEA(segmenti_antiacceptor_knownpos);
    FREEA(segmenti_donor_knowni);
    FREEA(segmentj_acceptor_knowni);
    FREEA(segmentj_antidonor_knowni);
    FREEA(segmenti_antiacceptor_knowni);
  } else {
    FREE(segmenti_donor_knownpos);
    FREE(segmentj_acceptor_knownpos);
    FREE(segmentj_antidonor_knownpos);
    FREE(segmenti_antiacceptor_knownpos);
    FREE(segmenti_donor_knowni);
    FREE(segmentj_acceptor_knowni);
    FREE(segmentj_antidonor_knowni);
    FREE(segmenti_antiacceptor_knowni);
  }
#else
  FREE(segmenti_donor_knownpos);
  FREE(segmentj_acceptor_knownpos);
  FREE(segmentj_antidonor_knownpos);
  FREE(segmenti_antiacceptor_knownpos);
  FREE(segmenti_donor_knowni);
  FREE(segmentj_acceptor_knowni);
  FREE(segmentj_antidonor_knowni);
  FREE(segmenti_antiacceptor_knowni);
#endif

  return twopartp;
}
#endif


static int
get_diagonals (Univdiag_T *middle_diagonal, List_T *best_right_diagonals, List_T *best_left_diagonals, 
	       List_T *all_right_diagonals, List_T *all_left_diagonals,
	       T sarray, char *queryptr, int querylength, Compress_T query_compress,
	       Univcoord_T chroffset, Univcoord_T chrhigh,
	       Univcoord_T goal, Elt_T *original_elt_array, int best_i, int nelts,
	       bool plusp, int genestrand, char conversion[]) {
  int best_score_right, best_score_left, best_score, score;
  Elt_T elt, right_elt;
  List_T *elt_tree;
  Univcoord_T low, high;
  int max_leftward, min_leftward, skip_left;
  int querystart, queryend;

  Sarrayptr_T initptr, finalptr;
  bool successp;
  UINT4 nmatches;

  int i, j, k;
  List_T p;

  Univdiag_T *diagonal_array, diagonal, prev_diagonal;
  int querypos;
  int ndiagonals;
  List_T left_diagonals, right_diagonals;

#ifdef SUBDIVIDE_NOMATCHES
  Chrpos_T low_chrpos, high_chrpos;
  List_T sub_diagonals;
  Diag_T sub_diagonal;
  int nfound;
  Univcoord_T mappingstart, mappingend;
  int maxnconsecutive = 0;

  Oligoindex_T oligoindex;
  Chrpos_T **mappings, chrstart, chrend;
  bool *coveredp;
  int *npositions, totalpositions = 0;
#endif

  int max_insertionlen;

  if (max_middle_insertions_default >= 0) {
    max_insertionlen = max_insertionlen_default;
  } else {
    max_insertionlen = querylength;
  }

  debug13(printf("\n***Entered get_diagonals, plusp %d, with goal %u\n",plusp,goal));

  /* Make elt tree, which allows for subdivisions of an elt */
  elt_tree = (List_T *) MALLOC(nelts*sizeof(List_T));
  for (i = 0; i < nelts; i++) {
    elt_tree[i] = List_push(NULL,(void *) original_elt_array[i]);
  }


  /* Compute leftward extensions for right side */
  debug13(printf("Performing leftward extensions for right side.  Subtracting %d and adding %d\n",
		 max_insertionlen,overall_max_distance));
  low = subtract_bounded(goal,/*minusterm*/max_insertionlen,chroffset);
  high = add_bounded(goal,/*plusterm*/overall_max_distance,chrhigh);
  for (i = best_i + 1; i < nelts; i++) {
    elt = (Elt_T) elt_tree[i]->first;
    Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
    if (elt->npositions > 0) {
      /* Success: Update low and high for next search */
      low = subtract_bounded(elt->positions[0],/*minusterm*/max_insertionlen,chroffset);
      high = add_bounded(elt->positions[elt->npositions-1],/*plusterm*/overall_max_distance,chrhigh);
    } else {
      debug13(printf("Elt %d..%d (leftward %d..%d) has no positions, so trying to reduce elt->queryend\n",
		     elt->querystart,elt->queryend,elt->querystart_leftward,elt->queryend));
      if (i + 1 < nelts) {
	/* A.  Try moving boundary to the left */
	right_elt = (Elt_T) elt_tree[i+1]->first;
	Elt_fill_positions_filtered(right_elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	if ((max_leftward = Elt_extend_leftward(&min_leftward,right_elt,query_compress,
						plusp,genestrand,/*skip_left*/0)) > 0) {
	  debug13(printf("Can extend %d..%d leftward by max %d, min %d\n",
			 right_elt->querystart,right_elt->queryend,max_leftward,min_leftward));
	  right_elt->querystart_leftward -= min_leftward; /* Using min_leftward is conservative */
	  queryend = right_elt->querystart_leftward - 2;

	  j = i;
	  while (j >= best_i && ((Elt_T) elt_tree[j]->first)->querystart_leftward >= queryend) {
	    debug13(printf("Left-extension of elt %d..%d => %d..%d obliterates elt %d..%d => %d..%d\n",
			   right_elt->querystart,right_elt->queryend,right_elt->querystart_leftward,right_elt->queryend_leftward,
			   ((Elt_T) elt_tree[j]->first)->querystart,((Elt_T) elt_tree[j]->first)->queryend,((Elt_T) elt_tree[j]->first)->querystart_leftward,queryend));
	    --j;
	  }

	  if (j >= best_i) {
	    /* Create a new elt with new positions */
	    querystart = ((Elt_T) elt_tree[j]->first)->querystart_leftward;
	    /* queryend was computed above */
	    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querystart]),
			  /*querylength*/(queryend + 1) - querystart,/*queryoffset*/querystart,
			  query_compress,sarray,plusp,genestrand,conversion);
	    elt_tree[j] = List_pop(elt_tree[j],(void **) &elt);
	    if (elt->temporaryp == true) {
	      Elt_free(&elt);
	    }
	    elt = Elt_new(querystart,nmatches,initptr,finalptr,/*temporaryp*/true);
	    elt_tree[j] = List_push(NULL,(void *) elt);
	    Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	  }
	}
      }

      if (elt->npositions > 0) {
	/* Success: Update low and high for next search */
	low = subtract_bounded(elt->positions[0],/*minusterm*/max_insertionlen,chroffset);
	high = add_bounded(elt->positions[elt->npositions-1],/*plusterm*/overall_max_distance,chrhigh);
      }
    }
  }


  /* Compute leftward extensions for left side */
  debug13(printf("Performing leftward extensions for left side.  Subtracting %d and adding %d\n",
		 overall_max_distance,max_insertionlen));
  low = subtract_bounded(goal,/*minusterm*/overall_max_distance,chroffset);
  high = add_bounded(goal,/*plusterm*/max_insertionlen,chrhigh);
  for (i = best_i - 1; i >= 0; --i) {
    elt = (Elt_T) elt_tree[i]->first;
    Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
    if (elt->npositions > 0) {
      /* Success: Update low and high for next search */
      low = subtract_bounded(elt->positions[0],/*minusterm*/overall_max_distance,chroffset);
      high = add_bounded(elt->positions[elt->npositions-1],/*plusterm*/max_insertionlen,chrhigh);
    } else {
      /* A.  Try moving boundary to the left */
      debug13(printf("Elt %d..%d has no positions, so trying to reduce elt->queryend\n",
		     elt->querystart,elt->queryend));
      if (i + 1 < nelts) {
	right_elt = (Elt_T) elt_tree[i+1]->first;
	skip_left = 0;
	if ((max_leftward = Elt_extend_leftward(&min_leftward,right_elt,query_compress,
						plusp,genestrand,/*skip_left*/0)) == 0) {
	  skip_left = 1;
	  max_leftward = Elt_extend_leftward(&min_leftward,right_elt,query_compress,
					     plusp,genestrand,skip_left);
	  debug13(printf("On second try, min_leftward is %d, max_leftward is %d\n",min_leftward,max_leftward));
	}

	if (max_leftward > 0) {
	  debug13(printf("Can extend %d..%d leftward by max %d, min %d\n",
			 right_elt->querystart,right_elt->queryend,max_leftward,min_leftward));
	  right_elt->querystart_leftward -= min_leftward + skip_left; /* Using min_leftward is conservative */
	  queryend = right_elt->querystart_leftward - 2;
	  
	  j = i;
	  while (j >= best_i && ((Elt_T) elt_tree[j]->first)->querystart_leftward >= queryend) {
	    debug13(printf("Left-extension of elt %d..%d => %d..%d obliterates elt %d..%d => %d..%d\n",
			   right_elt->querystart,right_elt->queryend,right_elt->querystart_leftward,right_elt->querystart_leftward,
			   ((Elt_T) elt_tree[j]->first)->querystart,((Elt_T) elt_tree[j]->first)->queryend,((Elt_T) elt_tree[j]->first)->querystart_leftward,queryend));
	    --j;
	  }
	  
	  if (j >= 0) {
	    /* Create a new elt with new positions */
	    querystart = ((Elt_T) elt_tree[j]->first)->querystart_leftward;
	    /* queryend was computed above */
	    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querystart]),
			  /*querylength*/(queryend + 1) - querystart,/*queryoffset*/querystart,
			  query_compress,sarray,plusp,genestrand,conversion);
	    elt_tree[j] = List_pop(elt_tree[j],(void **) &elt);
	    if (elt->temporaryp == true) {
	      Elt_free(&elt);
	    }
	    elt = Elt_new(querystart,nmatches,initptr,finalptr,/*temporaryp*/true);
	    elt_tree[j] = List_push(NULL,(void *) elt);
	    Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	  }
	}
      }

      if (elt->npositions > 0) {
	/* Success: Update low and high for next search */
	low = subtract_bounded(elt->positions[0],/*minusterm*/overall_max_distance,chroffset);
	high = add_bounded(elt->positions[elt->npositions-1],/*plusterm*/max_insertionlen,chrhigh);
      }
    }
  }
    
#ifdef SUBDIVIDE_NOMATCHES
  /* Try to subdivide elts that have no matches */
#ifdef HAVE_ALLOCA
  coveredp = (bool *) CALLOCA(querylength,sizeof(bool));
  mappings = (Chrpos_T **) ALLOCA(querylength * sizeof(Chrpos_T *));
  npositions = (int *) CALLOCA(querylength,sizeof(int));
#else
  coveredp = (bool *) CALLOC(querylength,sizeof(bool));
  mappings = (Chrpos_T **) MALLOC(querylength * sizeof(Chrpos_T *));
  npositions = (int *) CALLOC(querylength,sizeof(int));
#endif
  oligoindex = Oligoindex_array_elt(oligoindices_minor,/*source*/0);
  indexsize = Oligoindex_indexsize(oligoindex);


  debug13(printf("Starting subdivisions on right side\n"));
  low = subtract_bounded(goal,/*minusterm*/max_insertionlen,chroffset);
  high = add_bounded(goal,/*plusterm*/overall_max_distance,chrhigh);
  i = best_i + 1;
  while (i < nelts) {
    elt = (Elt_T) elt_tree[i]->first;
    debug13(printf("Elt #%d at %d..%d has %d matching positions\n",i,elt->querystart,elt->queryend,elt->npositions));

    if (elt->npositions > 0) {
      low = subtract_bounded(elt->positions[0],/*minusterm*/max_insertionlen,chroffset);
      high = add_bounded(elt->positions[elt->npositions-1],/*plusterm*/overall_max_distance,chrhigh);
      i++;
    } else {
      j = i;
      querystart = elt->querystart_leftward;
      while (j + 1 < nelts && ((Elt_T) elt_tree[j+1]->first)->npositions <= 0) {
	j = j + 1;
      }
      elt = (Elt_T) elt_tree[j]->first;
      queryend = elt->queryend_leftward;
      debug13(printf("Elts from %d through %d have no matching positions\n",i,j));

#if 0
      nfound = 0;
      /* B.  Try subdividing elt using 16-mers every 8 */
      debug13(printf("B.  Try to subdivide elt region at %d..%d\n",querystart,queryend));
      for (querypos = queryend - 16; querypos >= querystart; querypos -= 8) {
	sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querypos]),
		      /*querylength*/16,/*queryoffset*/querypos,
		      query_compress,sarray,plusp,genestrand,conversion);
	elt = Elt_new(querypos,nmatches,initptr,finalptr,/*temporaryp*/true);
	elt_tree[i] = List_push(elt_tree[i],(void *) elt);
	Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	nfound += elt->npositions;
	debug13(printf("Subelt at %d..%d has %d matching positions\n",elt->querystart,elt->queryend,elt->npositions));
      }

      if (nfound == 0) {
	/* C.  Try subdividing elt using 16-mers every 1 */
	debug13(printf("C.  Try to subdivide elt region at %d..%d\n",querystart,queryend));
	for (querypos = queryend - 16; querypos >= querystart; querypos -= 1) {
	  sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querypos]),
			/*querylength*/16,/*queryoffset*/querypos,
			query_compress,sarray,plusp,genestrand,conversion);
	  elt = Elt_new(querypos,nmatches,initptr,finalptr,/*temporaryp*/true);
	  elt_tree[i] = List_push(elt_tree[i],(void *) elt);
	  Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	  nfound += elt->npositions;
	  debug13(printf("Subelt at %d..%d has %d matching positions\n",elt->querystart,elt->queryend,elt->npositions));
	}
      }

      if (nfound == 0) {
	/* D.  Try subdividing elt using 8-mers every 1 */
	debug13(printf("D.  Try to subdivide elt region at %d..%d\n",querystart,queryend));
	for (querypos = queryend - 8; querypos >= querystart; querypos -= 1) {
	  sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querypos]),
			/*querylength*/8,/*queryoffset*/querypos,
			query_compress,sarray,plusp,genestrand,conversion);
	  elt = Elt_new(querypos,nmatches,initptr,finalptr,/*temporaryp*/true);
	  elt_tree[i] = List_push(elt_tree[i],(void *) elt);
	  Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	  nfound += elt->npositions;
	  debug13(printf("Subelt at %d..%d has %d matching positions\n",elt->querystart,elt->queryend,elt->npositions));
	}
      }

#else

      mappingstart = low + querystart;
      mappingend = high + queryend;
      chrstart = mappingstart - chroffset;
      chrend = mappingend - chroffset;

      Oligoindex_hr_tally(oligoindex,mappingstart,mappingend,/*plusp:true*/true,
			  queryptr,querystart,queryend,/*chrpos*/chrstart,genestrand);
      sub_diagonals = Oligoindex_get_mappings(NULL,coveredp,mappings,npositions,&totalpositions,
					      &oned_matrix_p,&maxnconsecutive,oligoindices_minor,oligoindex,
					      queryptr,querystart,queryend,querylength,
					      chrstart,chrend,chroffset,chrhigh,/*plusp:true*/true,diagpool);
      Oligoindex_untally(oligoindex,queryptr,querylength);

      debug14(printf("Got %d sub diagonals\n",List_length(sub_diagonals)));
      for (p = sub_diagonals; p != NULL; p = List_next(p)) {
	sub_diagonal = (Diag_T) List_head(p);
	debug14(printf("%d..%d %u\n",sub_diagonal->querystart,sub_diagonal->queryend + indexsize - 1,chrstart + sub_diagonal->diagonal));
	elt = Elt_new_fillin(sub_diagonal->querystart,sub_diagonal->queryend,indexsize,chroffset + chrstart + sub_diagonal->diagonal);
	elt_tree[i] = List_push(elt_tree[i],(void *) elt);
      }

#endif

      i = j + 1;
    }
  }


  debug13(printf("Starting subdivisions on left side\n"));
  low = subtract_bounded(goal,/*minusterm*/overall_max_distance,chroffset);
  high = add_bounded(goal,/*plusterm*/max_insertionlen,chrhigh);
  i = best_i - 1;
  while (i >= 0) {
    elt = (Elt_T) elt_tree[i]->first;
    debug13(printf("Elt #%d at %d..%d has %d matching positions\n",i,elt->querystart,elt->queryend,elt->npositions));

    if (elt->npositions > 0) {
      low = subtract_bounded(elt->positions[0],/*minusterm*/overall_max_distance,chroffset);
      high = add_bounded(elt->positions[elt->npositions-1],/*plusterm*/max_insertionlen,chrhigh);
      --i;

    } else {
      j = i;
      queryend = elt->queryend_leftward;
      while (j - 1 >= 0 && ((Elt_T) elt_tree[j-1]->first)->npositions <= 0) {
	j = j - 1;
      }
      elt = (Elt_T) elt_tree[j]->first;
      querystart = elt->querystart_leftward;
      debug13(printf("Elts from %d through %d have no matching positions\n",i,j));

#if 0
      nfound = 0;
      /* B.  Try subdividing elt using 16-mers every 8 */
      debug13(printf("B.  Try to subdivide elt region at %d..%d\n",querystart,queryend));
      for (querypos = queryend - 16; querypos >= querystart; querypos -= 8) {
	sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querystart]),
		      /*querylength*/16,/*queryoffset*/querystart,
		      query_compress,sarray,plusp,genestrand,conversion);
	elt = Elt_new(querystart,nmatches,initptr,finalptr,/*temporaryp*/true);
	elt_tree[i] = List_push(elt_tree[i],(void *) elt);
	Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	nfound += elt->npositions;
	debug13(printf("Subelt at %d..%d has %d matching positions\n",elt->querystart,elt->queryend,elt->npositions));
      }

      if (nfound == 0) {
	/* C.  Try subdividing elt using 16-mers every 1 */
	debug13(printf("C.  Try to subdivide elt region at %d..%d\n",querystart,queryend));
	for (querypos = queryend - 16; querypos >= querystart; querypos -= 1) {
	  sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querypos]),
			/*querylength*/16,/*queryoffset*/querypos,
			query_compress,sarray,plusp,genestrand,conversion);
	  elt = Elt_new(querypos,nmatches,initptr,finalptr,/*temporaryp*/true);
	  elt_tree[i] = List_push(elt_tree[i],(void *) elt);
	  Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	  nfound += elt->npositions;
	  debug13(printf("Subelt at %d..%d has %d matching positions\n",elt->querystart,elt->queryend,elt->npositions));
	}
      }

      if (nfound == 0) {
	/* D.  Try subdividing elt using 8-mers every 1 */
	debug13(printf("D.  Try to subdivide elt region at %d..%d\n",querystart,queryend));
	for (querypos = queryend - 8; querypos >= querystart; querypos -= 1) {
	  sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querypos]),
			/*querylength*/8,/*queryoffset*/querypos,
			query_compress,sarray,plusp,genestrand,conversion);
	  elt = Elt_new(querypos,nmatches,initptr,finalptr,/*temporaryp*/true);
	  elt_tree[i] = List_push(elt_tree[i],(void *) elt);
	  Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,/*multiplep*/false);
	  nfound += elt->npositions;
	  debug13(printf("Subelt at %d..%d has %d matching positions\n",elt->querystart,elt->queryend,elt->npositions));
	}
      }

#else

      mappingstart = low + querystart;
      mappingend = high + queryend;
      chrstart = mappingstart - chroffset;
      chrend = mappingend - chroffset;

      Oligoindex_hr_tally(oligoindex,mappingstart,mappingend,/*plusp:true*/true,
			  queryptr,querystart,queryend,/*chrpos*/chrstart,genestrand);
      sub_diagonals = Oligoindex_get_mappings(NULL,coveredp,mappings,npositions,&totalpositions,
					      &oned_matrix_p,&maxnconsecutive,oligoindices_minor,oligoindex,
					      queryptr,querystart,queryend,querylength,
					      chrstart,chrend,chroffset,chrhigh,/*plusp:true*/true,diagpool);
      Oligoindex_untally(oligoindex,queryptr,querylength);

      debug14(printf("Got %d sub diagonals\n",List_length(sub_diagonals)));
      for (p = sub_diagonals; p != NULL; p = List_next(p)) {
	sub_diagonal = (Diag_T) List_head(p);
	debug14(printf("%d..%d %u\n",sub_diagonal->querystart,sub_diagonal->queryend + indexsize - 1,chrstart + sub_diagonal->diagonal));
	elt = Elt_new_fillin(sub_diagonal->querystart,sub_diagonal->queryend,indexsize,chroffset + chrstart + sub_diagonal->diagonal);
	elt_tree[i] = List_push(elt_tree[i],(void *) elt);
      }
#endif

      i = j - 1;
    }
  }
#endif


  /* Create diagonals.  We give a bonus of +1 for being on the same
     diagonal.  This means that we should count consecutive regions
     within each diagonal as 2 points.  Then an indel or gap will
     give only 1 point, or a relative penalty. */
  assert(List_length(elt_tree[best_i]) == 1);
  elt = (Elt_T) elt_tree[best_i]->first;
  /* Don't use leftward values */
  *middle_diagonal = Univdiag_new(elt->querystart,elt->queryend,/*univdiagonal*/goal);
  (*middle_diagonal)->intscore = 2*(elt->queryend - elt->querystart + 1);
  debug13(printf("Creating middle diagonal: query %d..%d, diagonal %u = goal %u - chroffset %u\n",
		 elt->querystart,elt->queryend,goal - chroffset,goal,chroffset));
  if (elt->temporaryp == true) {
    Elt_free(&elt);
  } else {
    Elt_reset(elt);
  }
  List_free(&(elt_tree[best_i]));


  right_diagonals = (List_T) NULL;
  for (i = nelts - 1; i > best_i; --i) { /* Go in this order to avoid reversing list at the end */
    for (p = elt_tree[i]; p != NULL; p = List_next(p)) {
      elt = (Elt_T) p->first;
      if (elt->fillin_p == true) {
	/* Created by oligoindex */
	diagonal = Univdiag_new(elt->querystart_leftward,elt->queryend_leftward,/*univdiagonal*/elt->positions[0]);
	diagonal->nmismatches_known_p = false;
	right_diagonals = List_push(right_diagonals,(void *) diagonal);
      } else if (elt->querystart_leftward < elt->queryend_leftward) {
	for (j = elt->npositions - 1; j >= 0; --j) {  /* Go in this order to avoid reversing list at the end */
	  debug13(printf("Creating right diagonal: query %d..%d (leftward %d..%d), diagonal %u\n",
			 elt->querystart,elt->queryend,elt->querystart_leftward,elt->queryend_leftward,elt->positions[j] - chroffset));
	  right_diagonals = List_push(right_diagonals,Univdiag_new(elt->querystart_leftward,elt->queryend_leftward,
								   /*univdiagonal*/elt->positions[j]));
	}
      }
      if (elt->temporaryp == true) {
	Elt_free(&elt);
      } else {
	Elt_reset(elt);
      }
    }
    List_free(&(elt_tree[i]));
  }


  left_diagonals = (List_T) NULL;
  for (i = 0; i < best_i; i++) { /* Go in this order to avoid reversing list at the end */
    for (p = elt_tree[i]; p != NULL; p = List_next(p)) {
      elt = (Elt_T) p->first;
      if (elt->fillin_p == true) {
	/* Created by oligoindex */
	diagonal = Univdiag_new(elt->querystart_leftward,elt->queryend_leftward,/*univdiagonal*/elt->positions[0]);
	diagonal->nmismatches_known_p = false; /* Signifies that we don't know the number of mismatches */
	left_diagonals = List_push(left_diagonals,(void *) diagonal);
      } else if (elt->querystart_leftward < elt->queryend_leftward) {
	for (j = 0; j < elt->npositions; j++) {	/* Go in this order to avoid reversing list at the end */
	  debug13(printf("Creating left diagonal: query %d..%d (leftward %d..%d), diagonal %u\n",
			 elt->querystart,elt->queryend,elt->querystart_leftward,elt->queryend_leftward,elt->positions[j] - chroffset));
	  left_diagonals = List_push(left_diagonals,Univdiag_new(elt->querystart_leftward,elt->queryend_leftward,
								 /*univdiagonal*/elt->positions[j]));
	}
      }
      if (elt->temporaryp == true) {
	Elt_free(&elt);
      } else {
	Elt_reset(elt);
      }
    }
    List_free(&(elt_tree[i]));
  }

  FREE(elt_tree);



  /* A.  Compute right diagonals */
  /* A1.  Scoring for dynamic programming */
  diagonal_array = (Univdiag_T *) List_to_array_n(&ndiagonals,right_diagonals);
  List_free(&right_diagonals);
#ifdef DEBUG12
  printf("Right side before consolidating\n");
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    printf("%d..%d at %u\n",diagonal->querystart,diagonal->queryend,diagonal->univdiagonal);
  }
#endif

  *all_right_diagonals = (List_T) NULL;
  qsort(diagonal_array,ndiagonals,sizeof(Univdiag_T),Univdiag_diagonal_cmp);
  i = 0;
  while (i < ndiagonals) {
    j = i;
    while (j < ndiagonals && diagonal_array[j]->univdiagonal == diagonal_array[i]->univdiagonal) {
      j++;
    }
    if (j == i) {
      *all_right_diagonals = List_push(*all_right_diagonals,(void *) diagonal_array[i]);
    } else {
      *all_right_diagonals = List_push(*all_right_diagonals,
				      (void *) Univdiag_new(diagonal_array[i]->querystart,
							    diagonal_array[j-1]->queryend,
							    diagonal_array[i]->univdiagonal));
      for (k = i; k < j; k++) {
	Univdiag_free(&(diagonal_array[k]));
      }
    }
    i = j;
  }
  FREE(diagonal_array);

  /* TODO: May be able to skip this sorting step */
  diagonal_array = (Univdiag_T *) List_to_array_n(&ndiagonals,*all_right_diagonals);
  qsort(diagonal_array,ndiagonals,sizeof(Univdiag_T),Univdiag_ascending_cmp);
#ifdef DEBUG12
  printf("Right side after consolidating and sorting\n");
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    printf("%d..%d at %u\n",diagonal->querystart,diagonal->queryend,diagonal->univdiagonal);
  }
#endif


  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    debug13(printf("%d: %d..%d at %u\n",i,diagonal->querystart,diagonal->queryend,diagonal->univdiagonal));

    low = subtract_bounded(diagonal->univdiagonal,overall_max_distance,chroffset);
    high = add_bounded(diagonal->univdiagonal,max_insertionlen,chrhigh);
    querypos = diagonal->querystart;
    best_score = 0;

    for (j = i - 1; j >= 0; --j) {
      prev_diagonal = diagonal_array[j];
      debug13(printf("  %d: %d..%d at %u  ",j,prev_diagonal->querystart,prev_diagonal->queryend,prev_diagonal->univdiagonal));

      if (prev_diagonal->queryend >= querypos) {
	debug13(printf("Skipping because queryend %d >= querypos %d\n",prev_diagonal->queryend,querypos));
      } else if (prev_diagonal->univdiagonal < low) {
	debug13(printf("Skipping because diagonal %u < low_chrpos %u\n",prev_diagonal->univdiagonal,low));
      } else if (prev_diagonal->univdiagonal > high) {
	debug13(printf("Skipping because diagonal %u > high_chrpos %u\n",prev_diagonal->univdiagonal,high));
      } else {
	score = prev_diagonal->intscore;
	if (prev_diagonal->univdiagonal == diagonal->univdiagonal) {
	  score += 1;
	}
	if (score <= best_score) {
	  debug13(printf("Skipping because score %d <= best_score %d\n",score,best_score));
	} else {
	  best_score = score;
	  diagonal->prev = prev_diagonal;
	  debug13(printf("Updating best score to be %d.  Prev diagonal is %d..%d at %u\n",
			 best_score,prev_diagonal->querystart,prev_diagonal->queryend,prev_diagonal->univdiagonal));
	}
      }
    }

    /* Handle links to middle diagonal */
    prev_diagonal = *middle_diagonal;
    debug13(printf("  Middle: %d..%d at %u  ",prev_diagonal->querystart,prev_diagonal->queryend,prev_diagonal->univdiagonal));
    if (prev_diagonal->queryend >= querypos) {
      debug13(printf("Skipping because queryend %d >= querypos %d\n",prev_diagonal->queryend,querypos));
    } else if (prev_diagonal->univdiagonal < low) {
      debug13(printf("Skipping because diagonal %u < low_chrpos %u\n",prev_diagonal->univdiagonal,low));
    } else if (prev_diagonal->univdiagonal > high) {
      debug13(printf("Skipping because diagonal %u > high_chrpos %u\n",prev_diagonal->univdiagonal,high));
    } else {
      score = prev_diagonal->intscore;
      if (prev_diagonal->univdiagonal == diagonal->univdiagonal) {
	score += 1;		/* This bonus means we should double count contiguous region within each segment */
      }
      if (score <= best_score) {
	debug13(printf("Skipping because score %d <= best_score %d\n",score,best_score));
      } else {
	best_score = score;
	/* diagonal->prev = (Univdiag_T) NULL; */
	debug13(printf("Updating best score (for link to middle diagonal) to be %d\n",best_score));
      }
    }

    diagonal->intscore = best_score + 2*diagonal->nconsecutive;
    debug13(printf("Right diagonal %d..%d at %u gets score %d\n",
		   diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->intscore));
  }
  FREE(diagonal_array);


  /* A2.  Optimizing for dynamic programming */
  best_score_right = 0;
  *best_right_diagonals = (List_T) NULL;
  for (p = *all_right_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    if (diagonal->intscore > best_score_right) {
      best_score_right = diagonal->intscore;
      List_free(&(*best_right_diagonals));
      *best_right_diagonals = List_push(NULL,(void *) diagonal);
    } else if (diagonal->intscore == best_score_right) {
      *best_right_diagonals = List_push(*best_right_diagonals,(void *) diagonal);
    }
  }


  /* C.  Compute left diagonals */
  /* C1.  Scoring for dynamic programming */
  diagonal_array = (Univdiag_T *) List_to_array_n(&ndiagonals,left_diagonals);
  List_free(&left_diagonals);
#ifdef DEBUG12
  printf("Left side before consolidating\n");
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    printf("%d..%d at %u\n",diagonal->querystart,diagonal->queryend,diagonal->univdiagonal);
  }
#endif

  *all_left_diagonals = (List_T) NULL;
  qsort(diagonal_array,ndiagonals,sizeof(Univdiag_T),Univdiag_diagonal_cmp);
  i = 0;
  while (i < ndiagonals) {
    j = i;
    while (j < ndiagonals && diagonal_array[j]->univdiagonal == diagonal_array[i]->univdiagonal) {
      j++;
    }
    if (j == i) {
      *all_left_diagonals = List_push(*all_left_diagonals,(void *) diagonal_array[i]);
    } else {
      *all_left_diagonals = List_push(*all_left_diagonals,
				      (void *) Univdiag_new(diagonal_array[i]->querystart,
							    diagonal_array[j-1]->queryend,
							    diagonal_array[i]->univdiagonal));
      for (k = i; k < j; k++) {
	Univdiag_free(&(diagonal_array[k]));
      }
    }
    i = j;
  }
  FREE(diagonal_array);

  /* TODO: May be able to skip this sorting step */
  diagonal_array = (Univdiag_T *) List_to_array_n(&ndiagonals,*all_left_diagonals);
  qsort(diagonal_array,ndiagonals,sizeof(Univdiag_T),Univdiag_descending_cmp);
#ifdef DEBUG12
  printf("Left side after consolidating and sorting\n");
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    printf("%d..%d at %u\n",diagonal->querystart,diagonal->queryend,diagonal->univdiagonal);
  }
#endif


  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    debug13(printf("%d: %d..%d at %u\n",i,diagonal->querystart,diagonal->queryend,diagonal->univdiagonal));

    low = subtract_bounded(diagonal->univdiagonal,max_insertionlen,chroffset);
    high = add_bounded(diagonal->univdiagonal,overall_max_distance,chrhigh);
    querypos = diagonal->queryend;
    best_score = 0;

    for (j = i - 1; j >= 0; --j) {
      prev_diagonal = diagonal_array[j];
      debug13(printf("  %d: %d..%d at %u  ",j,prev_diagonal->querystart,prev_diagonal->queryend,prev_diagonal->univdiagonal));

      if (prev_diagonal->querystart <= querypos) {
	debug13(printf("Skipping because querystart %d <= querypos %d\n",prev_diagonal->querystart,querypos));
      } else if (prev_diagonal->univdiagonal < low) {
	debug13(printf("Skipping because diagonal %u < low %u\n",prev_diagonal->univdiagonal,low));
      } else if (prev_diagonal->univdiagonal > high) {
	debug13(printf("Skipping because diagonal %u > high %u\n",prev_diagonal->univdiagonal,high));
      } else {
	score = prev_diagonal->intscore;
	if (prev_diagonal->univdiagonal == diagonal->univdiagonal) {
	  score += 1;
	}
	if (score <= best_score) {
	  debug13(printf("Skipping because score %d <= best_score %d\n",score,best_score));
	} else {
	  best_score = score;
	  diagonal->prev = prev_diagonal;
	  debug13(printf("Updating best score to be %d.  Prev diagonal is %d..%d at %u\n",
			 best_score,prev_diagonal->querystart,prev_diagonal->queryend,prev_diagonal->univdiagonal));
	}
      }
    }

    /* Handle links to middle diagonal */
    prev_diagonal = *middle_diagonal;
    debug13(printf("  Middle: %d..%d at %u  ",prev_diagonal->querystart,prev_diagonal->queryend,prev_diagonal->univdiagonal));
    if (prev_diagonal->querystart <= querypos) {
      debug13(printf("Skipping because querystart %d <= querypos %d\n",prev_diagonal->querystart,querypos));
    } else if (prev_diagonal->univdiagonal < low) {
      debug13(printf("Skipping because diagonal %u < low_chrpos %u\n",prev_diagonal->univdiagonal,low));
    } else if (prev_diagonal->univdiagonal > high) {
      debug13(printf("Skipping because diagonal %u > high_chrpos %u\n",prev_diagonal->univdiagonal,high));
    } else {
      score = prev_diagonal->intscore;
      if (prev_diagonal->univdiagonal == diagonal->univdiagonal) {
	score += 1;		/* This bonus means we should double count contiguous region within each segment */
      }
      if (score <= best_score) {
	debug13(printf("Skipping because score %d <= best_score %d\n",prev_diagonal->intscore,best_score));
      } else {
	best_score = score;
	/* diagonal->prev = (Univdiag_T) NULL; */
	debug13(printf("Updating best score (for link to middle diagonal) to be %d\n",best_score));
      }
    }

    diagonal->intscore = best_score + 2*diagonal->nconsecutive;
    debug13(printf("Left diagonal %d..%d at %u gets score %d\n",
		   diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->intscore));
  }
  FREE(diagonal_array);


  /* C2.  Optimizing for dynamic programming */
  best_score_left = 0;
  *best_left_diagonals = (List_T) NULL;
  for (p = *all_left_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    if (diagonal->intscore > best_score_left) {
      best_score_left = diagonal->intscore;
      List_free(&(*best_left_diagonals));
      *best_left_diagonals = List_push(NULL,(void *) diagonal);
    } else if (diagonal->intscore == best_score_left) {
      *best_left_diagonals = List_push(*best_left_diagonals,(void *) diagonal);
    }
  }

#if 0
  printf("Best on the left\n");
  for (p = *best_left_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    printf("Score %d: %d..%d at %u\n",diagonal->intscore,diagonal->querystart,diagonal->queryend,diagonal->diagonal);
  }
#endif


  if (best_score_left == 0 && best_score_right == 0) {
    return (*middle_diagonal)->intscore;
  } else if (best_score_left == 0) {
    return best_score_right;
  } else if (best_score_right == 0) {
    return best_score_left;
  } else {
    /* middle_diagonal score is double counted */
    return best_score_left + best_score_right - (*middle_diagonal)->intscore;
  }
}


static List_T
find_best_path (List_T *right_paths, Intlist_T *right_endpoints_sense, Intlist_T *right_endpoints_antisense,
		Intlist_T *right_queryends_sense, Intlist_T *right_queryends_antisense,
		Uintlist_T *right_ambcoords_sense, Uintlist_T *right_ambcoords_antisense,
		Intlist_T *right_amb_knowni_sense, Intlist_T *right_amb_knowni_antisense,
		Intlist_T *right_amb_nmismatchesi_sense, Intlist_T *right_amb_nmismatchesi_antisense,
		Intlist_T *right_amb_nmismatchesj_sense, Intlist_T *right_amb_nmismatchesj_antisense,
		Doublelist_T *right_amb_probsi_sense, Doublelist_T *right_amb_probsi_antisense,
		Doublelist_T *right_amb_probsj_sense, Doublelist_T *right_amb_probsj_antisense,

		List_T *left_paths, Intlist_T *left_endpoints_sense, Intlist_T *left_endpoints_antisense,
		Intlist_T *left_querystarts_sense, Intlist_T *left_querystarts_antisense,
		Uintlist_T *left_ambcoords_sense, Uintlist_T *left_ambcoords_antisense,
		Intlist_T *left_amb_knowni_sense, Intlist_T *left_amb_knowni_antisense,
		Intlist_T *left_amb_nmismatchesi_sense, Intlist_T *left_amb_nmismatchesi_antisense,
		Intlist_T *left_amb_nmismatchesj_sense, Intlist_T *left_amb_nmismatchesj_antisense,
		Doublelist_T *left_amb_probsi_sense, Doublelist_T *left_amb_probsi_antisense,
		Doublelist_T *left_amb_probsj_sense, Doublelist_T *left_amb_probsj_antisense,

		List_T *fillin_diagonals,

		Univdiag_T middle_diagonal, List_T best_right_diagonals, List_T best_left_diagonals,

		int querylength, Compress_T query_compress, Univcoord_T chroffset,
		bool plusp, int genestrand, int max_mismatches_allowed) {
  List_T middle_path;
  List_T p;

  List_T diagonal_path, ambig_path;
  Univdiag_T diagonal, common_diagonal, prev_diagonal, right_indel_diagonal = NULL, left_indel_diagonal = NULL;
  int nbest;

#ifdef SUBDIVIDE_ENDS
  int indexsize;
  bool oned_matrix_p;
  Chrpos_T **mappings, chrstart, chrend;
  int maxnconsecutive = 0;
  int *npositions, totalpositions = 0;
  bool *coveredp;
  int querystart, queryend;
  Univcoord_T mappingstart, mappingend;
  List_T sub_diagonals;
  Diag_T sub_diagonal;
  Oligoindex_T oligoindex;
#endif
  Univcoord_T left, prev_left;
  
  /* Chrpos_T splice_distance; */
  int splice_pos;
  int best_knowni_i, best_knowni_j, best_nmismatches_i, best_nmismatches_j;
  double best_prob_i, best_prob_j;

  Chrpos_T first_dist_sense, second_dist_sense, first_dist_antisense, second_dist_antisense;
  double first_prob_sense, second_prob_sense, first_prob_antisense, second_prob_antisense;
  int firsti_sense, secondi_sense, firsti_antisense, secondi_antisense;
  int sensei, antisensei;

  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;
  int *segmenti_donor_knownpos, *segmentj_acceptor_knownpos, *segmentj_antidonor_knownpos, *segmenti_antiacceptor_knownpos,
    *segmenti_donor_knowni, *segmentj_acceptor_knowni, *segmentj_antidonor_knowni, *segmenti_antiacceptor_knowni;
  int j;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    segmenti_donor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


  debug13(printf("***Entered find_best_path\n"));

#ifdef SUBDIVIDE_ENDS
  mappings = (Chrpos_T **) MALLOCA(querylength * sizeof(Chrpos_T *));
  indexsize = Oligoindex_indexsize(oligoindex);
  coveredp = (bool *) CALLOCA(querylength,sizeof(bool));
  npositions = (int *) CALLOCA(querylength,sizeof(int));
  oligoindex = Oligoindex_array_elt(oligoindices_minor,/*source*/0);
#endif

  /* A3.  Traceback for dynamic programming */
  *right_endpoints_sense = *right_endpoints_antisense = (Intlist_T) NULL;
  *right_queryends_sense = *right_queryends_antisense = (Intlist_T) NULL;
  *right_ambcoords_sense = *right_ambcoords_antisense = (Uintlist_T) NULL;
  *right_amb_knowni_sense = *right_amb_knowni_antisense = (Intlist_T) NULL;
  *right_amb_nmismatchesi_sense = *right_amb_nmismatchesi_antisense = (Intlist_T) NULL;
  *right_amb_nmismatchesj_sense = *right_amb_nmismatchesj_antisense = (Intlist_T) NULL;
  *right_amb_probsi_sense = *right_amb_probsi_antisense = (Doublelist_T) NULL;
  *right_amb_probsj_sense = *right_amb_probsj_antisense = (Doublelist_T) NULL;

  *right_paths = (List_T) NULL;
  if ((nbest = List_length(best_right_diagonals)) == 0) {
    common_diagonal = (Univdiag_T) NULL;

    /* querystart = middle_diagonal->queryend + 1; */
    left = middle_diagonal->univdiagonal;

  } else if (nbest == 1) {
    common_diagonal = (Univdiag_T) List_head(best_right_diagonals);

    /* querystart = common_diagonal->queryend + 1; */
    left = common_diagonal->univdiagonal;

  } else {
    debug13(printf("Multiple (%d) best right diagonals\n",nbest));

    /* Distinguish between common and divergent diagonals */
    for (p = best_right_diagonals; p != NULL; p = List_next(p)) {
      diagonal = (Univdiag_T) List_head(p);
      while (diagonal != NULL) {
	diagonal->nlinked += 1;
	diagonal = diagonal->prev;
      }
    }

    /* Handle divergent diagonals */
    /* Now that we are running oligoindex, we may need to obtain only the last common_diagonal */
    for (p = best_right_diagonals; p != NULL; p = List_next(p)) {
      ambig_path = (List_T) NULL;
      diagonal = (Univdiag_T) List_head(p);
      while (diagonal != NULL && diagonal->nlinked < nbest) {
	ambig_path = List_push(ambig_path,(void *) diagonal);
	diagonal = diagonal->prev;
      }
      *right_paths = List_push(*right_paths,(void *) ambig_path);

      common_diagonal = diagonal; /* Last elt on prev path.  Save for later */
    }

    if (common_diagonal == NULL) {
      /* All paths connect directly to the middle diagonal, so there is no common diagonal */
      prev_diagonal = middle_diagonal;
      /* querystart = middle_diagonal->queryend + 1; */
      prev_left = middle_diagonal->univdiagonal;
    } else {
      prev_diagonal = common_diagonal;
      /* querystart = common_diagonal->queryend + 1; */
      prev_left = common_diagonal->univdiagonal;
    }

    /* Distinguish right paths by looking for indel (which wins) or splicing */
    debug13(printf("Have %d right_paths.  Distinguish by looking for indels\n",List_length(*right_paths)));
    for (p = *right_paths; p != NULL; p = List_next(p)) {
      ambig_path = (List_T) List_head(p);
      diagonal = (Univdiag_T) List_head(ambig_path);
      left = diagonal->univdiagonal;
      debug13(printf("left %u, prev_left %u, difference %d\n",left,prev_left,(int) left - prev_left));
      if (left < prev_left) {
	/* Insertion */
	debug13(printf("Found insertion\n"));
	right_indel_diagonal = diagonal;
      } else if (left - prev_left < MIN_INTRONLEN) {
	/* Deletion */
	debug13(printf("Found deletion\n"));
	right_indel_diagonal = diagonal;
      }
    }

    if (right_indel_diagonal != NULL) {
      /* Push onto middle path later */
      /* querystart = right_indel_diagonal->queryend + 1; */
      left = right_indel_diagonal->univdiagonal;

    } else {
      debug13(printf("Still have %d right_paths.  Distinguish by looking for best splice\n",List_length(*right_paths)));
      first_dist_sense = second_dist_sense = 0;
      first_prob_sense = second_prob_sense = 0.0;
      firsti_sense = secondi_sense = -1;
      first_dist_antisense = second_dist_antisense = 0;
      first_prob_antisense = second_prob_antisense = 0.0;
      firsti_antisense = secondi_antisense = -1;
      sensei = antisensei = 0;

      for (p = *right_paths; p != NULL; p = List_next(p)) {
	ambig_path = (List_T) List_head(p);
	diagonal = (Univdiag_T) List_head(ambig_path);
	left = diagonal->univdiagonal;

	segmenti_donor_nknown = segmenti_antiacceptor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(prev_left,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,prev_left);
	  while (j < nsplicesites && splicesites[j] < prev_left + querylength) {
	    if (splicetypes[j] == DONOR) {
	      debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - prev_left;
	      segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - prev_left;
	      segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmenti_donor_knownpos[segmenti_donor_nknown] = querylength + 100;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength + 100;
	  
	segmentj_acceptor_nknown = segmentj_antidonor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left);
	  while (j < nsplicesites && splicesites[j] < left + querylength) {
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - left;
	      segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - left;
	      segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength + 100;
	segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength + 100;
      
	/* splice_distance = left - prev_left; */
#if 0
	max_mismatches_allowed = (diagonal->querystart - prev_diagonal->queryend - 1);
	debug13(printf("max_mismatches %d = %d - %d - 1\n",max_mismatches_allowed,diagonal->querystart,prev_diagonal->queryend));
	if (prev_diagonal->intscore > 0) {
	  max_mismatches_allowed += 1;
	}
	if (diagonal->intscore > 0) {
	  max_mismatches_allowed += 1;
	}
#endif
      
	if ((splice_pos = Splice_resolve_sense(&best_knowni_i,&best_knowni_j,&best_nmismatches_i,&best_nmismatches_j,
					       &best_prob_i,&best_prob_j,
					       /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
					       prev_diagonal->querystart,diagonal->queryend+1,querylength,query_compress,
					       segmenti_donor_knownpos,segmentj_acceptor_knownpos,
					       segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					       segmenti_donor_knowni,segmentj_acceptor_knowni,
					       segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
					       segmenti_donor_nknown,segmentj_acceptor_nknown,
					       segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
					       max_mismatches_allowed,plusp,genestrand)) >= 0) {
	  debug13(printf("Found sense splice_pos %d with probs %f and %f\n",splice_pos,best_prob_i,best_prob_j));
	  if (best_prob_i + best_prob_j > first_prob_sense) {
	    second_dist_sense = first_dist_sense;
	    second_prob_sense = first_prob_sense;
	    secondi_sense = firsti_sense;
	    first_dist_sense = left - prev_left;
	    first_prob_sense = best_prob_i + best_prob_j;
	    firsti_sense = sensei;
	  } else if (best_prob_i + best_prob_j > second_prob_sense) {
	    second_dist_sense = left - prev_left;
	    second_prob_sense = best_prob_i + best_prob_j;
	    secondi_sense = sensei;
	  }

	  *right_endpoints_sense = Intlist_push(*right_endpoints_sense,splice_pos);
	  *right_queryends_sense = Intlist_push(*right_queryends_sense,diagonal->queryend + 1);
	  *right_ambcoords_sense = Uintlist_push(*right_ambcoords_sense,left + splice_pos);
	  *right_amb_knowni_sense = Intlist_push(*right_amb_knowni_sense,best_knowni_j);
	  *right_amb_nmismatchesi_sense = Intlist_push(*right_amb_nmismatchesi_sense,best_nmismatches_i);
	  *right_amb_nmismatchesj_sense = Intlist_push(*right_amb_nmismatchesj_sense,best_nmismatches_j);
	  *right_amb_probsi_sense = Doublelist_push(*right_amb_probsi_sense,best_prob_i);
	  *right_amb_probsj_sense = Doublelist_push(*right_amb_probsj_sense,best_prob_j);
	  sensei++;
	}

	if ((splice_pos = Splice_resolve_antisense(&best_knowni_i,&best_knowni_j,&best_nmismatches_i,&best_nmismatches_j,
						   &best_prob_i,&best_prob_j,
						   /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
						   prev_diagonal->querystart,diagonal->queryend+1,querylength,query_compress,
						   segmenti_donor_knownpos,segmentj_acceptor_knownpos,
						   segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
						   segmenti_donor_knowni,segmentj_acceptor_knowni,
						   segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
						   segmenti_donor_nknown,segmentj_acceptor_nknown,
						   segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
						   max_mismatches_allowed,plusp,genestrand)) >= 0) {
	  debug13(printf("Found antisense splice_pos %d with probs %f and %f\n",splice_pos,best_prob_i,best_prob_j));
	  if (best_prob_i + best_prob_j > first_prob_antisense) {
	    second_dist_antisense = first_dist_antisense;
	    second_prob_antisense = first_prob_antisense;
	    secondi_antisense = firsti_antisense;
	    first_dist_antisense = left - prev_left;
	    first_prob_antisense = best_prob_i + best_prob_j;
	    firsti_antisense = antisensei;
	  } else if (best_prob_i + best_prob_j > second_prob_antisense) {
	    second_dist_antisense = left - prev_left;
	    second_prob_antisense = best_prob_i + best_prob_j;
	    secondi_antisense = antisensei;
	  }

	  *right_endpoints_antisense = Intlist_push(*right_endpoints_antisense,splice_pos);
	  *right_queryends_antisense = Intlist_push(*right_queryends_antisense,diagonal->queryend + 1);
	  *right_ambcoords_antisense = Uintlist_push(*right_ambcoords_antisense,left + splice_pos);
	  *right_amb_knowni_antisense = Intlist_push(*right_amb_knowni_antisense,best_knowni_j);
	  *right_amb_nmismatchesi_antisense = Intlist_push(*right_amb_nmismatchesi_antisense,best_nmismatches_i);
	  *right_amb_nmismatchesj_antisense = Intlist_push(*right_amb_nmismatchesj_antisense,best_nmismatches_j);
	  *right_amb_probsi_antisense = Doublelist_push(*right_amb_probsi_antisense,best_prob_i);
	  *right_amb_probsj_antisense = Doublelist_push(*right_amb_probsj_antisense,best_prob_j);
	  antisensei++;
	}
      }

      if (Intlist_length(*right_endpoints_sense) > 1) {
	if (first_dist_sense < second_dist_sense/2) {
	  debug13(printf("first dist sense %u is significantly shorter than second dist sense %u.  Keeping %d from end\n",
			 first_dist_sense,second_dist_sense,firsti_sense));
	  firsti_sense = Intlist_length(*right_endpoints_sense) - 1 - firsti_sense; /* Because we don't reverse lists */
	  *right_endpoints_sense = Intlist_keep_one(*right_endpoints_sense,firsti_sense);
	  *right_queryends_sense = Intlist_keep_one(*right_queryends_sense,firsti_sense);
	  *right_ambcoords_sense = Uintlist_keep_one(*right_ambcoords_sense,firsti_sense);
	  *right_amb_knowni_sense = Intlist_keep_one(*right_amb_knowni_sense,firsti_sense);
	  *right_amb_nmismatchesi_sense = Intlist_keep_one(*right_amb_nmismatchesi_sense,firsti_sense);
	  *right_amb_nmismatchesj_sense = Intlist_keep_one(*right_amb_nmismatchesj_sense,firsti_sense);
	  *right_amb_probsi_sense = Doublelist_keep_one(*right_amb_probsi_sense,firsti_sense);
	  *right_amb_probsj_sense = Doublelist_keep_one(*right_amb_probsj_sense,firsti_sense);
	}
      }

      if (Intlist_length(*right_endpoints_antisense) > 1) {
	if (first_dist_antisense < second_dist_antisense/2) {
	  debug13(printf("first dist antisense %u is significantly shorter than second dist antisense %u.  Keeping %d from end\n",
			 first_dist_antisense,second_dist_antisense,firsti_antisense));
	  firsti_antisense = Intlist_length(*right_endpoints_antisense) - 1 - firsti_antisense; /* Because we don't reverse lists */
	  *right_endpoints_antisense = Intlist_keep_one(*right_endpoints_antisense,firsti_antisense);
	  *right_queryends_antisense = Intlist_keep_one(*right_queryends_antisense,firsti_antisense);
	  *right_ambcoords_antisense = Uintlist_keep_one(*right_ambcoords_antisense,firsti_antisense);
	  *right_amb_knowni_antisense = Intlist_keep_one(*right_amb_knowni_antisense,firsti_antisense);
	  *right_amb_nmismatchesi_antisense = Intlist_keep_one(*right_amb_nmismatchesi_antisense,firsti_antisense);
	  *right_amb_nmismatchesj_antisense = Intlist_keep_one(*right_amb_nmismatchesj_antisense,firsti_antisense);
	  *right_amb_probsi_antisense = Doublelist_keep_one(*right_amb_probsi_antisense,firsti_antisense);
	  *right_amb_probsj_antisense = Doublelist_keep_one(*right_amb_probsj_antisense,firsti_antisense);
	}
      }
    }
  }


#ifdef SUBDIVIDE_ENDS
  sub_diagonals = (List_T) NULL;

  if (querystart + MIN_ENDLENGTH >= querylength) {
  } else {
    /* Run oligoindex here to right of common_diagonal */
    mappingstart = subtract_bounded(left + querystart,/*minusterm*/max_insertionlen,chroffset);
    mappingend = add_bounded(left + querylength,/*plusterm*/overall_max_distance,chrhigh);
    chrstart = mappingstart - chroffset;
    chrend = mappingend - chroffset;

    Oligoindex_hr_tally(oligoindex,mappingstart,mappingend,/*plusp:true*/true,
			queryptr,querystart,/*queryend*/querylength,/*chrpos*/chrstart,genestrand);
    sub_diagonals = Oligoindex_get_mappings(NULL,coveredp,mappings,npositions,&totalpositions,
					    &oned_matrix_p,&maxnconsecutive,oligoindices_minor,oligoindex,
					    queryptr,querystart,/*queryend*/querylength,querylength,
					    chrstart,chrend,chroffset,chrhigh,/*plusp:true*/true,diagpool);
    Oligoindex_untally(oligoindex,queryptr,querylength);

    debug14(printf("Got %d sub diagonals\n",List_length(sub_diagonals)));
#ifdef DEBUG14
    for (p = sub_diagonals; p != NULL; p = List_next(p)) {
      sub_diagonal = (Diag_T) List_head(p);
      /* Need to alter oligoindex diagonal for our needs */
      printf("%d..%d %u\n",sub_diagonal->querystart,sub_diagonal->queryend + indexsize - 1,chrstart + sub_diagonal->diagonal);
    }
#endif

#if 0
    /* Perform dynamic programming on these diagonals */
    for (p = sub_diagonals; p != NULL; p = List_next(p)) {
      diagonal = List_head(p);
      querypos = diagonal->querystart;
      best_score = 0;

      for (q = sub_diagonals; q != p; q = List_next(q)) {
	prev_diagonal = List_head(q);
	if (prev_diagonal->queryend >= querypos) {
	  debug13(printf("Skipping because queryend %d >= querypos %d\n",prev_diagonal->queryend,querypos));
	} else if (prev_diagonal->univdiagonal < low) {
	  debug13(printf("Skipping because diagonal %u < low_chrpos %u\n",prev_diagonal->diagonal,low_chrpos));
	} else if (prev_diagonal->diagonal > high_chrpos) {
	  debug13(printf("Skipping because diagonal %u > high_chrpos %u\n",prev_diagonal->diagonal,high_chrpos));
	} else {
	  score = prev_diagonal->intscore;
	  if (prev_diagonal->diagonal == diagonal->diagonal) {
	    score += 1;
	  }
	  if (score <= best_score) {
	    debug13(printf("Skipping because score %d <= best_score %d\n",score,best_score));
	  } else {
	    best_score = score;
	    diagonal->prev = prev_diagonal;
	    debug13(printf("Updating best score to be %d.  Prev diagonal is %d..%d at %u\n",
			   best_score,prev_diagonal->querystart,prev_diagonal->queryend,prev_diagonal->diagonal));
	  }
	}
      }
    }
#endif

  }
#endif	/* SUBDIVIDE_ENDS */


  *fillin_diagonals = (List_T) NULL;
  middle_path = (List_T) NULL;

#ifdef SUBDIVIDE_ENDS
  /* Without SUBDIVIDE_ENDS, sub_diagonals is guaranteed to be NULL */
  /* A4.  Process oligoindex diagonals from right */
  if (List_length(sub_diagonals) == 0) {
    /* Skip */
  } else if (List_length(sub_diagonals) == 1) {
    sub_diagonal = List_head(sub_diagonals);
    diagonal = Univdiag_new_fillin(sub_diagonal->querystart,sub_diagonal->queryend,indexsize,
				   /*univdiagonal*/chroffset + chrstart + sub_diagonal->diagonal);
    *fillin_diagonals = List_push(*fillin_diagonals,(void *) diagonal);
    middle_path = List_push(middle_path,(void *) diagonal);
  } else {
#ifdef DEBUG13
    printf("Have %d sub_diagonals\n",List_length(sub_diagonals));
    for (p = sub_diagonals; p != NULL; p = List_next(p)) {
      sub_diagonal = List_head(p);
      printf("%d..%d %u\n",sub_diagonal->querystart,sub_diagonal->queryend,chrstart + sub_diagonal->diagonal);
    }
#endif
  }
#endif

  if (right_indel_diagonal != NULL) {
    debug13(printf("Pushing right indel diagonal onto middle: query %d..%d, diagonal %u\n",
		   right_indel_diagonal->querystart,right_indel_diagonal->queryend,right_indel_diagonal->univdiagonal - chroffset));
    middle_path = List_push(middle_path,(void *) right_indel_diagonal);
  }

  /* A5. Process common diagonal from right */
  while (common_diagonal != NULL) {
    middle_path = List_push(middle_path,(void *) common_diagonal);
    debug13(printf("Pushing common diagonal onto middle: query %d..%d, diagonal %u\n",
		   common_diagonal->querystart,common_diagonal->queryend,common_diagonal->univdiagonal - chroffset));
    common_diagonal = common_diagonal->prev;
  }

  /* B. Process original middle diagonal */
  middle_path = List_push(middle_path,(void *) middle_diagonal);
  debug13(printf("Pushing middle diagonal onto middle: query %d..%d, diagonal %u\n",
		 middle_diagonal->querystart,middle_diagonal->queryend,middle_diagonal->univdiagonal - chroffset));


  /* C3.  Traceback for dynamic programming */
  *left_endpoints_sense = *left_endpoints_antisense = (Intlist_T) NULL;
  *left_querystarts_sense = *left_querystarts_antisense = (Intlist_T) NULL;
  *left_ambcoords_sense = *left_ambcoords_antisense = (Uintlist_T) NULL;
  *left_amb_knowni_sense = *left_amb_knowni_antisense = (Intlist_T) NULL;
  *left_amb_nmismatchesi_sense = *left_amb_nmismatchesi_antisense = (Intlist_T) NULL;
  *left_amb_nmismatchesj_sense = *left_amb_nmismatchesj_antisense = (Intlist_T) NULL;
  *left_amb_probsi_sense = *left_amb_probsi_antisense = (Doublelist_T) NULL;
  *left_amb_probsj_sense = *left_amb_probsj_antisense = (Doublelist_T) NULL;

  *left_paths = (List_T) NULL;
  debug13(printf("On left, have %d best_left_diagonals\n",List_length(best_left_diagonals)));
  if ((nbest = List_length(best_left_diagonals)) == 0) {
    common_diagonal = (Univdiag_T) NULL;

    /* queryend = middle_diagonal->querystart; */
    left = middle_diagonal->univdiagonal;

  } else if (nbest == 1) {
    common_diagonal = (Univdiag_T) List_head(best_left_diagonals);

    /* queryend = common_diagonal->querystart; */
    left = common_diagonal->univdiagonal;

  } else {
    debug13(printf("Multiple (%d) best left diagonals\n",nbest));

    /* Distinguish between common and divergent diagonals */
    for (p = best_left_diagonals; p != NULL; p = List_next(p)) {
      diagonal = (Univdiag_T) List_head(p);
      while (diagonal != NULL) {
	diagonal->nlinked += 1;
	diagonal = diagonal->prev;
      }
    }

    /* Handle divergent diagonals */
    /* Now that we are running oligoindex, we may need to obtain only the last common_diagonal */
    for (p = best_left_diagonals; p != NULL; p = List_next(p)) {
      ambig_path = (List_T) NULL;
      diagonal = (Univdiag_T) List_head(p);
      while (diagonal != NULL && diagonal->nlinked < nbest) {
	ambig_path = List_push(ambig_path,(void *) diagonal);
	diagonal = diagonal->prev;
      }
      *left_paths = List_push(*left_paths,(void *) ambig_path);

      common_diagonal = diagonal; /* Last elt on prev path.  Save for later */
    }

    if (common_diagonal == NULL) {
      /* All paths connect directly to the middle diagonal, so there is no common diagonal */
      diagonal = middle_diagonal;
      /* queryend = middle_diagonal->querystart; */
      left = middle_diagonal->univdiagonal;
    } else {
      diagonal = common_diagonal;
      /* queryend = common_diagonal->querystart; */
      left = common_diagonal->univdiagonal;
    }

    /* Distinguish left paths by looking for indel (which wins) or splicing */
    debug13(printf("Have %d left_paths.  Distinguish by looking for indel\n",List_length(*left_paths)));
    for (p = *left_paths; p != NULL; p = List_next(p)) {
      ambig_path = (List_T) List_head(p);
      prev_diagonal = (Univdiag_T) List_head(ambig_path);
      prev_left = prev_diagonal->univdiagonal;
      debug13(printf("left %u, prev_left %u, difference %d\n",left,prev_left,(int) left - prev_left));
      if (left < prev_left) {
	/* Insertion */
	debug13(printf("Found insertion\n"));
	left_indel_diagonal = prev_diagonal;
      } else if (left - prev_left < MIN_INTRONLEN) {
	/* Deletion */
	debug13(printf("Found deletion\n"));
	left_indel_diagonal = prev_diagonal;
      }
    }

    if (left_indel_diagonal != NULL) {
      /* Push onto middle path later */
      left = left_indel_diagonal->univdiagonal;
      /* queryend = left_indel_diagonal->querystart; */

    } else {
      debug13(printf("Still have %d left_paths.  Distinguish by looking for best splice\n",List_length(*left_paths)));
      first_dist_sense = second_dist_sense = 0;
      first_prob_sense = second_prob_sense = 0.0;
      firsti_sense = secondi_sense = -1;
      first_dist_antisense = second_dist_antisense = 0;
      first_prob_antisense = second_prob_antisense = 0.0;
      firsti_antisense = secondi_antisense = -1;
      sensei = antisensei = 0;

      for (p = *left_paths; p != NULL; p = List_next(p)) {
	ambig_path = (List_T) List_head(p);
	prev_diagonal = (Univdiag_T) List_head(ambig_path);
	prev_left = prev_diagonal->univdiagonal;

	segmenti_donor_nknown = segmenti_antiacceptor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(prev_left,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,prev_left);
	  while (j < nsplicesites && splicesites[j] < prev_left + querylength) {
	    if (splicetypes[j] == DONOR) {
	      debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - prev_left;
	      segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - prev_left;
	      segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmenti_donor_knownpos[segmenti_donor_nknown] = querylength + 100;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength + 100;
	  
	segmentj_acceptor_nknown = segmentj_antidonor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left);
	  while (j < nsplicesites && splicesites[j] < left + querylength) {
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - left;
	      segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - left;
	      segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength + 100;
	segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength + 100;
      
	/* splice_distance = left - prev_left; */
#if 0
	max_mismatches_allowed = (diagonal->querystart - prev_diagonal->queryend - 1);
	debug13(printf("max_mismatches %d = %d - %d - 1\n",max_mismatches_allowed,diagonal->querystart,prev_diagonal->queryend));
	if (prev_diagonal->intscore > 0) {
	  max_mismatches_allowed += 1;
	}
	if (diagonal->intscore > 0) {
	  max_mismatches_allowed += 1;
	}
#endif
      
	if ((splice_pos = Splice_resolve_sense(&best_knowni_i,&best_knowni_j,&best_nmismatches_i,&best_nmismatches_j,
					       &best_prob_i,&best_prob_j,
					       /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
					       prev_diagonal->querystart,diagonal->queryend+1,querylength,query_compress,
					       segmenti_donor_knownpos,segmentj_acceptor_knownpos,
					       segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					       segmenti_donor_knowni,segmentj_acceptor_knowni,
					       segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
					       segmenti_donor_nknown,segmentj_acceptor_nknown,
					       segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
					       max_mismatches_allowed,plusp,genestrand)) >= 0) {
	  debug13(printf("Found sense splice_pos %d with probs %f and %f\n",splice_pos,best_prob_i,best_prob_j));
	  if (best_prob_i + best_prob_j > first_prob_sense) {
	    second_dist_sense = first_dist_sense;
	    second_prob_sense = first_prob_sense;
	    secondi_sense = firsti_sense;
	    first_dist_sense = left - prev_left;
	    first_prob_sense = best_prob_i + best_prob_j;
	    firsti_sense = sensei;
	  } else if (best_prob_i + best_prob_j > second_prob_sense) {
	    second_dist_sense = left - prev_left;
	    second_prob_sense = best_prob_i + best_prob_j;
	    secondi_sense = sensei;
	  }

	  *left_endpoints_sense = Intlist_push(*left_endpoints_sense,splice_pos);
	  *left_querystarts_sense = Intlist_push(*left_querystarts_sense,prev_diagonal->querystart);
	  *left_ambcoords_sense = Uintlist_push(*left_ambcoords_sense,prev_left + splice_pos);
	  *left_amb_knowni_sense = Intlist_push(*left_amb_knowni_sense,best_knowni_i);
	  *left_amb_nmismatchesi_sense = Intlist_push(*left_amb_nmismatchesi_sense,best_nmismatches_i);
	  *left_amb_nmismatchesj_sense = Intlist_push(*left_amb_nmismatchesj_sense,best_nmismatches_j);
	  *left_amb_probsi_sense = Doublelist_push(*left_amb_probsi_sense,best_prob_i);
	  *left_amb_probsj_sense = Doublelist_push(*left_amb_probsj_sense,best_prob_j);
	  sensei++;
	}

	if ((splice_pos = Splice_resolve_antisense(&best_knowni_i,&best_knowni_j,&best_nmismatches_i,&best_nmismatches_j,
						   &best_prob_i,&best_prob_j,
						   /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
						   prev_diagonal->querystart,diagonal->queryend+1,querylength,query_compress,
						   segmenti_donor_knownpos,segmentj_acceptor_knownpos,
						   segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
						   segmenti_donor_knowni,segmentj_acceptor_knowni,
						   segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
						   segmenti_donor_nknown,segmentj_acceptor_nknown,
						   segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
						   max_mismatches_allowed,plusp,genestrand)) >= 0) {
	  debug13(printf("Found antisense splice_pos %d with probs %f and %f\n",splice_pos,best_prob_i,best_prob_j));
	  if (best_prob_i + best_prob_j > first_prob_antisense) {
	    second_dist_antisense = first_dist_antisense;
	    second_prob_antisense = first_prob_antisense;
	    secondi_antisense = firsti_antisense;
	    first_dist_antisense = left - prev_left;
	    first_prob_antisense = best_prob_i + best_prob_j;
	    firsti_antisense = antisensei;
	  } else if (best_prob_i + best_prob_j > second_prob_antisense) {
	    second_dist_antisense = left - prev_left;
	    second_prob_antisense = best_prob_i + best_prob_j;
	    secondi_antisense = antisensei;
	  }

	  *left_endpoints_antisense = Intlist_push(*left_endpoints_antisense,splice_pos);
	  *left_querystarts_antisense = Intlist_push(*left_querystarts_antisense,prev_diagonal->querystart);
	  *left_ambcoords_antisense = Uintlist_push(*left_ambcoords_antisense,prev_left + splice_pos);
	  *left_amb_knowni_antisense = Intlist_push(*left_amb_knowni_antisense,best_knowni_i);
	  *left_amb_nmismatchesi_antisense = Intlist_push(*left_amb_nmismatchesi_antisense,best_nmismatches_i);
	  *left_amb_nmismatchesj_antisense = Intlist_push(*left_amb_nmismatchesj_antisense,best_nmismatches_j);
	  *left_amb_probsi_antisense = Doublelist_push(*left_amb_probsi_antisense,best_prob_i);
	  *left_amb_probsj_antisense = Doublelist_push(*left_amb_probsj_antisense,best_prob_j);
	  antisensei++;
	}
      }

      if (Intlist_length(*left_endpoints_sense) > 1) {
	if (first_dist_sense < second_dist_sense/2) {
	  debug13(printf("first dist sense %u is significantly shorter than second dist sense %u.  Keeping %d from end\n",
			 first_dist_sense,second_dist_sense,firsti_sense));
	  firsti_sense = Intlist_length(*left_endpoints_sense) - 1 - firsti_sense; /* Because we don't reverse lists */
	  *left_endpoints_sense = Intlist_keep_one(*left_endpoints_sense,firsti_sense);
	  *left_querystarts_sense = Intlist_keep_one(*left_querystarts_sense,firsti_sense);
	  *left_ambcoords_sense = Uintlist_keep_one(*left_ambcoords_sense,firsti_sense);
	  *left_amb_knowni_sense = Intlist_keep_one(*left_amb_knowni_sense,firsti_sense);
	  *left_amb_nmismatchesi_sense = Intlist_keep_one(*left_amb_nmismatchesi_sense,firsti_sense);
	  *left_amb_nmismatchesj_sense = Intlist_keep_one(*left_amb_nmismatchesj_sense,firsti_sense);
	  *left_amb_probsi_sense = Doublelist_keep_one(*left_amb_probsi_sense,firsti_sense);
	  *left_amb_probsj_sense = Doublelist_keep_one(*left_amb_probsj_sense,firsti_sense);
	}
      }

      if (Intlist_length(*left_endpoints_antisense) > 1) {
	if (first_dist_antisense < second_dist_antisense/2) {
	  debug13(printf("first dist antisense %u is significantly shorter than second dist antisense %u.  Keeping %d from end\n",
			 first_dist_antisense,second_dist_antisense,firsti_antisense));
	  firsti_antisense = Intlist_length(*left_endpoints_antisense) - 1 - firsti_antisense; /* Because we don't reverse lists */
	  *left_endpoints_antisense = Intlist_keep_one(*left_endpoints_antisense,firsti_antisense);
	  *left_querystarts_antisense = Intlist_keep_one(*left_querystarts_antisense,firsti_antisense);
	  *left_ambcoords_antisense = Uintlist_keep_one(*left_ambcoords_antisense,firsti_antisense);
	  *left_amb_knowni_antisense = Intlist_keep_one(*left_amb_knowni_antisense,firsti_antisense);
	  *left_amb_nmismatchesi_antisense = Intlist_keep_one(*left_amb_nmismatchesi_antisense,firsti_antisense);
	  *left_amb_nmismatchesj_antisense = Intlist_keep_one(*left_amb_nmismatchesj_antisense,firsti_antisense);
	  *left_amb_probsi_antisense = Doublelist_keep_one(*left_amb_probsi_antisense,firsti_antisense);
	  *left_amb_probsj_antisense = Doublelist_keep_one(*left_amb_probsj_antisense,firsti_antisense);
	}
      }
    }
  }


#ifdef SUBDIVIDE_ENDS
  sub_diagonals = (List_T) NULL;

  /* Run oligoindex here to left of common_diagonal */
  if (queryend < MIN_ENDLENGTH) {
  } else {
    mappingstart = subtract_bounded(left + 0,/*minusterm*/overall_max_distance,chroffset);
    mappingend = add_bounded(left + queryend,/*plusterm*/max_insertionlen,chrhigh);
    chrstart = mappingstart - chroffset;
    chrend = mappingend - chroffset;

    Oligoindex_hr_tally(oligoindex,mappingstart,mappingend,/*plusp:true*/true,
			queryptr,/*querystart*/0,queryend,/*chrpos*/chrstart,genestrand);
    sub_diagonals = Oligoindex_get_mappings(NULL,coveredp,mappings,npositions,&totalpositions,
					    &oned_matrix_p,&maxnconsecutive,oligoindices_minor,oligoindex,
					    queryptr,/*querystart*/0,queryend,querylength,
					    chrstart,chrend,chroffset,chrhigh,/*plusp:true*/true,diagpool);
    Oligoindex_untally(oligoindex,queryptr,querylength);

    debug14(printf("Got %d sub diagonals\n",List_length(sub_diagonals)));
#ifdef DEBUG14
    for (p = sub_diagonals; p != NULL; p = List_next(p)) {
      sub_diagonal = (Diag_T) List_head(p);
      /* Need to alter oligoindex diagonal for our needs */
      printf("%d..%d %u\n",sub_diagonal->querystart,sub_diagonal->queryend + indexsize - 1,chrstart + sub_diagonal->diagonal);
    }
#endif
    /* Need to perform dynamic programming on these diagonals, or select one */
  }
#endif	/* SUBDIVIDE_ENDS */


  diagonal_path = (List_T) NULL;

  /* C5. Process left diagonals in reverse */
  while (common_diagonal != NULL) {
    diagonal_path = List_push(diagonal_path,(void *) common_diagonal);
    common_diagonal = common_diagonal->prev;
  }
  /* Pops off in reverse */
  for (p = diagonal_path; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    debug13(printf("Pushing common diagonal onto middle: query %d..%d, diagonal %u\n",
		   diagonal->querystart,diagonal->queryend,diagonal->univdiagonal - chroffset));
    middle_path = List_push(middle_path,(void *) diagonal);
  }
  List_free(&diagonal_path);


  if (left_indel_diagonal != NULL) {
    debug13(printf("Pushing left indel diagonal onto middle: query %d..%d, diagonal %u\n",
		   left_indel_diagonal->querystart,left_indel_diagonal->queryend,left_indel_diagonal->univdiagonal - chroffset));
    middle_path = List_push(middle_path,(void *) left_indel_diagonal);
  }


#ifdef SUBDIVIDE_ENDS
  /* Without SUBDIVIDE_ENDS, sub_diagonals is guaranteed to be NULL */
  /* C4. Process oligoindex diagonals from left */
  if (List_length(sub_diagonals) == 0) {
    /* Skip */
  } else if (List_length(sub_diagonals) == 1) {
    sub_diagonal = List_head(sub_diagonals);
    diagonal = Univdiag_new_fillin(sub_diagonal->querystart,sub_diagonal->queryend,indexsize,
				   /*univdiagonal*/chroffset + chrstart + sub_diagonal->diagonal);
    *fillin_diagonals = List_push(*fillin_diagonals,(void *) diagonal);
    middle_path = List_push(middle_path,(void *) diagonal);
  } else {
#ifdef DEBUG13
    printf("Have %d sub_diagonals\n",List_length(sub_diagonals));
    for (p = sub_diagonals; p != NULL; p = List_next(p)) {
      sub_diagonal = (Diag_T) List_head(p);
      printf("%d..%d %u\n",sub_diagonal->querystart,sub_diagonal->queryend,chrstart + sub_diagonal->diagonal);
    }
#endif
  }
#endif

  debug13(printf("***Exiting find_best_path\n"));

#ifdef SUBDIVIDE_ENDS
#ifdef HAVE_ALLOCA
  FREEA(npositions);
  FREEA(coveredp);
  FREEA(mappings);
#else
  FREE(npositions);
  FREE(coveredp);
  FREE(mappings);
#endif
#endif


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(segmenti_donor_knownpos);
    FREEA(segmentj_acceptor_knownpos);
    FREEA(segmentj_antidonor_knownpos);
    FREEA(segmenti_antiacceptor_knownpos);
    FREEA(segmenti_donor_knowni);
    FREEA(segmentj_acceptor_knowni);
    FREEA(segmentj_antidonor_knowni);
    FREEA(segmenti_antiacceptor_knowni);
  } else {
    FREE(segmenti_donor_knownpos);
    FREE(segmentj_acceptor_knownpos);
    FREE(segmentj_antidonor_knownpos);
    FREE(segmenti_antiacceptor_knownpos);
    FREE(segmenti_donor_knowni);
    FREE(segmentj_acceptor_knowni);
    FREE(segmentj_antidonor_knowni);
    FREE(segmenti_antiacceptor_knowni);
  }
#else
  FREE(segmenti_donor_knownpos);
  FREE(segmentj_acceptor_knownpos);
  FREE(segmentj_antidonor_knownpos);
  FREE(segmenti_antiacceptor_knownpos);
  FREE(segmenti_donor_knowni);
  FREE(segmentj_acceptor_knowni);
  FREE(segmentj_antidonor_knowni);
  FREE(segmenti_antiacceptor_knowni);
#endif

  return middle_path;
}



/* Note: This GMAP from sarray suffers from relying on middle_path and
end paths to get stage2.  Would be better to run oligoindex_hr to get
a better stage2, or to run GMAP from GSNAP or pairsearch */

#if 0
static List_T
run_gmap_plus (List_T gmap, List_T middle_path, List_T start_paths, List_T end_paths,
	       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
	       Chrpos_T chrlength, char *queryuc_ptr, int querylength,
	       int genestrand, bool first_read_p,
	       int maxpeelback, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	       Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  Stage3end_T hit;
  List_T stage2pairs, all_stage2_starts, all_stage2_ends;
  List_T p, q;

  int sensedir;

  struct Pair_T *pairarray;
  List_T pairs;
  List_T diagonal_path;
  Univdiag_T diagonal, prev_diagonal;
  int querypos;
  Chrpos_T genomepos;
  int c;

  int npairs, goodness, cdna_direction, matches, nmatches_posttrim,
    max_match_length, ambig_end_length_5, ambig_end_length_3,
    unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double ambig_prob_5, ambig_prob_3, min_splice_prob;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T start, end;
  int nsegments, nmismatches_whole, nindels, nintrons, nindelbreaks;


  /* D.  Make all_stage2_starts (paths) */
  all_stage2_starts = (List_T) NULL;
  diagonal = (Univdiag_T) List_head(middle_path);
  for (q = start_paths; q != NULL; q = List_next(q)) {
    q->first = diagonal_path = List_reverse((List_T) List_head(q));
    prev_diagonal = (Univdiag_T) List_head(diagonal_path);
    if (diagonal->univdiagonal > prev_diagonal->univdiagonal) {
      debug13(printf("START, PLUS\n"));
      stage2pairs = (List_T) NULL;
      for (p = diagonal_path; p != NULL; p = List_next(p)) {
	diagonal = (Univdiag_T) List_head(p);
	debug13(printf("Diagonal %d..%d at %u [%u]\n",
		       diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->univdiagonal - chroffset));
	querypos = diagonal->querystart;
	genomepos = diagonal->univdiagonal + diagonal->querystart - chroffset;
	while (querypos <= diagonal->queryend) {
	  c = queryuc_ptr[querypos];
	  stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				      /*cdna*/c,MATCH_COMP,/*genome*/c,/*genomealt*/c,
				      /*dynprogindex*/0);
	  debug13(printf("Pushing %c | %c at %d,%d\n",queryuc_ptr[querypos],queryuc_ptr[querypos],querypos,genomepos));
	  querypos++;
	  genomepos++;
	}
	debug13(printf("\n"));
      }
      all_stage2_starts = List_push(all_stage2_starts,(void *) stage2pairs);
    }
  }


  /* E.  Make all_stage2_ends (pairs) */
  all_stage2_ends = (List_T) NULL;
  prev_diagonal = (Univdiag_T) List_last_value(middle_path);
  for (q = end_paths; q != NULL; q = List_next(q)) {
    diagonal_path = (List_T) List_head(q);
    diagonal = (Univdiag_T) List_head(diagonal_path);
    if (diagonal->univdiagonal > prev_diagonal->univdiagonal) {
      debug13(printf("END, PLUS\n"));
      stage2pairs = (List_T) NULL;
      for (p = diagonal_path; p != NULL; p = List_next(p)) {
	diagonal = (Univdiag_T) List_head(p);
	debug13(printf("Diagonal %d..%d at %u [%u]\n",
		       diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->univdiagonal - chroffset));
	querypos = diagonal->querystart;
	genomepos = diagonal->univdiagonal + diagonal->querystart - chroffset;
	while (querypos <= diagonal->queryend) {
	  c = queryuc_ptr[querypos];
	  stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				      /*cdna*/c,MATCH_COMP,/*genome*/c,/*genomealt*/c,
				      /*dynprogindex*/0);
	  debug13(printf("Pushing %c | %c at %d,%d\n",queryuc_ptr[querypos],queryuc_ptr[querypos],querypos,genomepos));
	  querypos++;
	  genomepos++;
	}
	debug13(printf("\n"));
      }
      all_stage2_ends = List_push(all_stage2_ends,(void *) List_reverse(stage2pairs));
    }
  }


#ifdef DEBUG13
  printf("MIDDLE DIAGONALS, PLUS\n");
  for (p = middle_path; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    printf("Diagonal %d..%d at %u [%u]\n",
	   diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->univdiagonal - chroffset);
  }
#endif

  /* F.  Make stage2pairs */
  stage2pairs = (List_T) NULL;
  for (p = middle_path; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    querypos = diagonal->querystart;
    genomepos = diagonal->univdiagonal + diagonal->querystart - chroffset;
    while (querypos <= diagonal->queryend) {
      c = queryuc_ptr[querypos];
      stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				  /*cdna*/c,MATCH_COMP,/*genome*/c,/*genomealt*/c,
				  /*dynprogindex*/0);
      debug13(printf("Pushing %c | %c at %d,%d\n",queryuc_ptr[querypos],queryuc_ptr[querypos],querypos,genomepos));
      querypos++;
      genomepos++;
    }
    debug13(printf("\n"));
  }


  knownsplice_limit_high = ((Pair_T) stage2pairs->first)->genomepos + chroffset;
  stage2pairs = List_reverse(stage2pairs);
  knownsplice_limit_low = ((Pair_T) stage2pairs->first)->genomepos + chroffset;

  if ((pairarray = Stage3_compute(&pairs,&npairs,&goodness,&cdna_direction,&sensedir,
				  &matches,&nmatches_posttrim,&max_match_length,
				  &ambig_end_length_5,&ambig_end_length_3,
				  &ambig_splicetype_5,&ambig_splicetype_3,
				  &ambig_prob_5,&ambig_prob_3,
				  &unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
				  &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
				  stage2pairs,all_stage2_starts,all_stage2_ends,
#ifdef END_KNOWNSPLICING_SHORTCUT
				  cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
				  watsonp ? query_compress_fwd : query_compress_rev,
#endif
				  /*queryseq_ptr*/queryuc_ptr,queryuc_ptr,querylength,/*skiplength*/0,
#ifdef EXTRACT_GENOMICSEG
				  /*query_subseq_offset*/0,
#else
				  /*query_subseq_offset*/0,
#endif
				  chrnum,chroffset,chrhigh,
				  knownsplice_limit_low,knownsplice_limit_high,/*plusp*/true,genestrand,
				  /*jump_late_p*/false,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				  /*sense_try*/0,/*sense_filter*/0,
				  oligoindices_minor,diagpool,cellpool)) == NULL) {

  } else {
    nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
				     pairarray,npairs);
    start = subtract_bounded(chroffset + Pair_genomepos(&(pairarray[0])),
			     /*minusterm*/Pair_querypos(&(pairarray[0])),chroffset);
    end = add_bounded(chroffset + Pair_genomepos(&(pairarray[npairs-1])),
		      /*plusterm*/querylength - 1 - Pair_querypos(&(pairarray[npairs-1])),chrhigh);
    if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim,max_match_length,
				  ambig_end_length_5,ambig_end_length_3,
				  ambig_splicetype_5,ambig_splicetype_3,
				  min_splice_prob,
				  pairarray,npairs,nsegments,nintrons,nindelbreaks,
				  /*left*/start,/*genomiclength*/end - start + 1,
				  /*plusp*/true,genestrand,first_read_p,
				  /*accession*/NULL,querylength,chrnum,chroffset,chrhigh,chrlength,
				  cdna_direction,sensedir,/*sarrayp*/true)) == NULL) {
      FREE_OUT(pairarray);
    } else {
      gmap = List_push(gmap,(void *) hit);
    }
  }

  List_free(&all_stage2_ends);
  List_free(&all_stage2_starts);

  return gmap;
}
#endif


#if 0
static List_T
run_gmap_minus (List_T gmap, List_T middle_path, List_T start_paths, List_T end_paths,
		Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		Chrpos_T chrlength, char *queryuc_ptr, int querylength,
		int genestrand, bool first_read_p,
		int maxpeelback, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  Stage3end_T hit;
  List_T stage2pairs, all_stage2_starts, all_stage2_ends;

  List_T p, q;

  int sensedir;

  struct Pair_T *pairarray;
  List_T pairs;
  List_T diagonal_path;
  Univdiag_T diagonal, prev_diagonal;
  int querypos;
  Chrpos_T genomepos;
  int c;

  int npairs, goodness, cdna_direction, matches, nmatches_posttrim,
    max_match_length, ambig_end_length_5, ambig_end_length_3,
    unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double ambig_prob_5, ambig_prob_3, min_splice_prob;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T start, end;
  int nsegments, nmismatches_whole, nindels, nintrons, nindelbreaks;


  /* D.  Make all_stage2_starts (paths) */
  all_stage2_starts = (List_T) NULL;
  diagonal = (Univdiag_T) List_head(middle_path);
  for (q = start_paths; q != NULL; q = List_next(q)) {
    q->first = diagonal_path = List_reverse((List_T) List_head(q));
    prev_diagonal = (Univdiag_T) List_head(diagonal_path);
    if (diagonal->univdiagonal < prev_diagonal->univdiagonal) {
      debug13(printf("START, MINUS\n"));
      stage2pairs = (List_T) NULL;
      for (p = diagonal_path; p != NULL; p = List_next(p)) {
	diagonal = (Univdiag_T) List_head(p);
	debug13(printf("Diagonal %d..%d at %u [%u]\n",
		       diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->univdiagonal - chroffset));
	querypos = querylength - 1 - diagonal->queryend;
	genomepos = chrhigh - (diagonal->univdiagonal + diagonal->queryend);
	while (querypos <= querylength - 1 - diagonal->querystart) {
	  c = queryuc_ptr[querypos];
	  stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				      /*cdna*/c,MATCH_COMP,/*genome*/c,/*genomealt*/c,
				      /*dynprogindex*/0);
	  debug13(printf("Pushing %c | %c at %d,%d\n",queryuc_ptr[querypos],queryuc_ptr[querypos],querypos,genomepos));
	  querypos++;
	  genomepos++;
	}
	debug13(printf("\n"));
      }
      all_stage2_starts = List_push(all_stage2_starts,(void *) stage2pairs);
    }
  }


  /* E.  Make all_stage2_ends (pairs) */
  all_stage2_ends = (List_T) NULL;
  prev_diagonal = (Univdiag_T) List_last_value(middle_path);
  for (q = end_paths; q != NULL; q = List_next(q)) {
    diagonal_path = (List_T) List_head(q);
    diagonal = (Univdiag_T) List_head(diagonal_path);
    if (diagonal->univdiagonal < prev_diagonal->univdiagonal) {
      debug13(printf("END, MINUS\n"));
      stage2pairs = (List_T) NULL;
      for (p = diagonal_path; p != NULL; p = List_next(p)) {
	diagonal = (Univdiag_T) List_head(p);
	debug13(printf("Diagonal %d..%d at %u [%u]\n",
		       diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->univdiagonal - chroffset));
	querypos = querylength - 1 - diagonal->queryend;
	genomepos = chrhigh - (diagonal->univdiagonal + diagonal->queryend);
	while (querypos <= querylength - 1 - diagonal->querystart) {
	  c = queryuc_ptr[querypos];
	  stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				      /*cdna*/c,MATCH_COMP,/*genome*/c,/*genomealt*/c,
				      /*dynprogindex*/0);
	  debug13(printf("Pushing %c | %c at %d,%d\n",queryuc_ptr[querypos],queryuc_ptr[querypos],querypos,genomepos));
	  querypos++;
	  genomepos++;
	}
	debug13(printf("\n"));
      }
      all_stage2_ends = List_push(all_stage2_ends,(void *) List_reverse(stage2pairs));
    }
  }


#ifdef DEBUG13
  printf("MIDDLE DIAGONALS, MINUS\n");
  for (p = middle_path; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    printf("Diagonal %d..%d at %u [%u]\n",
	   diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->univdiagonal - chroffset);
  }
#endif

  /* F.  Make stage2pairs */
  stage2pairs = (List_T) NULL;
  middle_path = List_reverse(middle_path); /* For minus */
  for (p = middle_path; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    querypos = querylength - 1 - diagonal->queryend;
    assert(chrhigh > diagonal->univdiagonal + diagonal->queryend);
    genomepos = chrhigh - (diagonal->univdiagonal + diagonal->queryend);
    while (querypos <= querylength - 1 - diagonal->querystart) {
      c = queryuc_ptr[querypos];
      stage2pairs = Pairpool_push(stage2pairs,pairpool,querypos,genomepos,
				  /*cdna*/c,MATCH_COMP,/*genome*/c,/*genomealt*/c,
				  /*dynprogindex*/0);
      debug13(printf("Pushing %c | %c at %d,%d\n",queryuc_ptr[querypos],queryuc_ptr[querypos],querypos,genomepos));
      querypos++;
      genomepos++;
    }
    debug13(printf("\n"));
  }


  knownsplice_limit_low = ((Pair_T) stage2pairs->first)->genomepos + chroffset;
  stage2pairs = List_reverse(stage2pairs);
  knownsplice_limit_high = ((Pair_T) stage2pairs->first)->genomepos + chroffset;


  if ((pairarray = Stage3_compute(&pairs,&npairs,&goodness,&cdna_direction,&sensedir,
				  &matches,&nmatches_posttrim,&max_match_length,
				  &ambig_end_length_5,&ambig_end_length_3,
				  &ambig_splicetype_5,&ambig_splicetype_3,
				  &ambig_prob_5,&ambig_prob_3,
				  &unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
				  &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
				  stage2pairs,all_stage2_starts,all_stage2_ends,
#ifdef END_KNOWNSPLICING_SHORTCUT
				  cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
				  watsonp ? query_compress_fwd : query_compress_rev,
#endif
				  /*queryseq_ptr*/queryuc_ptr,queryuc_ptr,querylength,/*skiplength*/0,
#ifdef EXTRACT_GENOMICSEG
				  /*query_subseq_offset*/0,
#else
				  /*query_subseq_offset*/0,
#endif
				  chrnum,chroffset,chrhigh,
				  knownsplice_limit_low,knownsplice_limit_high,/*plusp*/false,genestrand,
				  /*jump_late_p*/true,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				  /*sense_try*/0,/*sense_filter*/0,
				  oligoindices_minor,diagpool,cellpool)) == NULL) {

  } else {
    nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
				     pairarray,npairs);
    start = add_bounded(chroffset + Pair_genomepos(&(pairarray[0])),
			/*plusterm*/Pair_querypos(&(pairarray[0])),chrhigh);
    end = subtract_bounded(chroffset + Pair_genomepos(&(pairarray[npairs-1])),
			   /*minusterm*/querylength - 1 - Pair_querypos(&(pairarray[npairs-1])),chroffset);
    if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim,max_match_length,
				  ambig_end_length_5,ambig_end_length_3,
				  ambig_splicetype_5,ambig_splicetype_3,
				  min_splice_prob,
				  pairarray,npairs,nsegments,nintrons,nindelbreaks,
				  /*left*/end,/*genomiclength*/start - end + 1,
				  /*plusp*/false,genestrand,first_read_p,
				  /*accession*/NULL,querylength,chrnum,chroffset,chrhigh,chrlength,
				  cdna_direction,sensedir,/*sarrayp*/true)) == NULL) {
      FREE_OUT(pairarray);
    } else {
      gmap = List_push(gmap,(void *) hit);
    }
  }

  List_free(&all_stage2_ends);
  List_free(&all_stage2_starts);

  return gmap;
}
#endif


#if 0
static bool
find_sense (int *sensedir, List_T sense_junctions, List_T antisense_junctions,
	    Intlist_T sense_endpoints, Intlist_T antisense_endpoints) {
  bool sense_acceptable_p = true, antisense_acceptable_p = true;
  double sense_prob = 0.0, antisense_prob = 0.0;
  Junction_T sense_junction, antisense_junction;
  List_T p;
  Intlist_T a;
  int last_endpoint;

  last_endpoint = -1;
  for (a = sense_endpoints; a != NULL; a = Intlist_next(a)) {
    if (Intlist_head(a) <= last_endpoint) {
      sense_acceptable_p = false;
    }
    last_endpoint = Intlist_head(a);
  }

  last_endpoint = -1;
  for (a = antisense_endpoints; a != NULL; a = Intlist_next(a)) {
    if (Intlist_head(a) <= last_endpoint) {
      antisense_acceptable_p = false;
    }
    last_endpoint = Intlist_head(a);
  }

  for (p = sense_junctions; p != NULL; p = List_next(p)) {
    sense_junction = (Junction_T) List_head(p);
    if (sense_junction == NULL) {
      sense_acceptable_p = false;
    } else if (Junction_type(sense_junction) == AMB_JUNCTION) {
      /* Ignore */
    } else {
      sense_prob += Junction_prob(sense_junction);
    }
  }

  for (p = antisense_junctions; p != NULL; p = List_next(p)) {
    antisense_junction = (Junction_T) List_head(p);
    if (antisense_junction == NULL) {
      antisense_acceptable_p = false;
    } else if (Junction_type(antisense_junction) == AMB_JUNCTION) {
      /* Ignore */
    } else {
      antisense_prob += Junction_prob(antisense_junction);
    }
  }

  if (sense_acceptable_p == false && antisense_acceptable_p == false) {
    return false;
  } else if (sense_acceptable_p == false) {
    *sensedir = SENSE_ANTI;
    return true;
  } else if (antisense_acceptable_p == false) {
    *sensedir = SENSE_FORWARD;
    return true;
  } else if (sense_prob > antisense_prob) {
    *sensedir = SENSE_FORWARD;
    return true;
  } else if (antisense_prob > sense_prob) {
    *sensedir = SENSE_ANTI;
    return true;
  } else {
    *sensedir = SENSE_NULL;
    return true;
  }
}
#endif


static bool
endpoints_acceptable_p (bool *intronp, List_T junctions, Intlist_T endpoints) {
  bool acceptable_p = true;
  Junction_T junction;
  List_T p;
  Intlist_T a;
  int last_endpoint;

  last_endpoint = -1;
  for (a = endpoints; a != NULL; a = Intlist_next(a)) {
    if (Intlist_head(a) <= last_endpoint) {
      acceptable_p = false;
    }
    last_endpoint = Intlist_head(a);
  }

  *intronp = false;
  for (p = junctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (junction == NULL) {
      acceptable_p = false;
    } else if (Junction_type(junction) == SPLICE_JUNCTION) {
      *intronp = true;
    }
  }

  return acceptable_p;
}



#if 0
static bool
incomplete_result_p (List_T middle_path, int querylength) {
  Univdiag_T diagonal;
  int querystart, queryend;

  diagonal = (Univdiag_T) List_head(middle_path);
  querystart = diagonal->querystart;

  diagonal = (Univdiag_T) List_last_value(middle_path);
  queryend = diagonal->queryend;

  if (querystart > 8 || queryend < querylength - 8) {
    return true;
  } else {
    return false;
  }
}
#endif


/* Always solves against plus strand of genome.  Just provide either
   queryuc/query_compress_fwd (coords measured from beginning of
   sequence) or queryrc/query_compress_rev (coords measured from end
   of sequence).  All coordinates measured from low end.
   Sense/antisense is with respect to the plus strand.  But to
   interface with Stage3end_new_substring command, need to flip
   coordinates for case where queryrc aligns to plus strand. */

static List_T
solve_via_segments (int *found_score, bool *completep, List_T hits, List_T middle_path,

		    Intlist_T right_endpoints_sense, Intlist_T right_endpoints_antisense,
		    Intlist_T right_queryends_sense, Intlist_T right_queryends_antisense,
		    Uintlist_T right_ambcoords_sense, Uintlist_T right_ambcoords_antisense,
		    Intlist_T right_amb_knowni_sense, Intlist_T right_amb_knowni_antisense,
		    Intlist_T right_amb_nmismatchesi_sense, Intlist_T right_amb_nmismatchesi_antisense,
		    Intlist_T right_amb_nmismatchesj_sense, Intlist_T right_amb_nmismatchesj_antisense,
		    Doublelist_T right_amb_probsi_sense, Doublelist_T right_amb_probsi_antisense,
		    Doublelist_T right_amb_probsj_sense, Doublelist_T right_amb_probsj_antisense,

		    Intlist_T left_endpoints_sense, Intlist_T left_endpoints_antisense,
		    Intlist_T left_querystarts_sense, Intlist_T left_querystarts_antisense,
		    Uintlist_T left_ambcoords_sense, Uintlist_T left_ambcoords_antisense,
		    Intlist_T left_amb_knowni_sense, Intlist_T left_amb_knowni_antisense,
		    Intlist_T left_amb_nmismatchesi_sense, Intlist_T left_amb_nmismatchesi_antisense,
		    Intlist_T left_amb_nmismatchesj_sense, Intlist_T left_amb_nmismatchesj_antisense,
		    Doublelist_T left_amb_probsi_sense, Doublelist_T left_amb_probsi_antisense,
		    Doublelist_T left_amb_probsj_sense, Doublelist_T left_amb_probsj_antisense,

		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    Chrpos_T chrlength, int querylength, Compress_T query_compress,
		    bool plusp, int genestrand, int max_mismatches_allowed) {
  List_T super_path;
  Stage3end_T hit;
  int sensedir, sense_sensedir, antisense_sensedir;

  List_T p;
  Univdiag_T diagonal, prev_diagonal, new_diagonal;
  Chrpos_T splice_distance;
  int querystart_for_merge, querystart, queryend;
  int max_leftward;
  int nmismatches, prev_nmismatches;
  bool fillin_p;

  int indel_pos;
  int nindels;
  Univcoord_T deletionpos;

  int splice_pos;
  double donor_prob, acceptor_prob;

  bool sense_acceptable_p, antisense_acceptable_p, sense_intronp, antisense_intronp;
  Univcoord_T left, prev_left;
  Uintlist_T sense_lefts = NULL, antisense_lefts = NULL;
  Intlist_T sense_nmismatches = NULL, antisense_nmismatches = NULL;
  Intlist_T sense_endpoints = NULL, antisense_endpoints = NULL;
  List_T sense_junctions = NULL, antisense_junctions = NULL;

  int best_knowni_i, best_knowni_j, best_nmismatches_i, best_nmismatches_j;
  double best_prob_i, best_prob_j;

  Substring_T right_ambig_sense, right_ambig_antisense,
    left_ambig_sense, left_ambig_antisense;
  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;
  int *segmenti_donor_knownpos, *segmentj_acceptor_knownpos, *segmentj_antidonor_knownpos, *segmenti_antiacceptor_knownpos,
    *segmenti_donor_knowni, *segmentj_acceptor_knowni, *segmentj_antidonor_knowni, *segmenti_antiacceptor_knowni;
  int j;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    segmenti_donor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


#ifdef DEBUG13
  printf("\n");
  printf("Original diagonals:\n");
  for (p = middle_path; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    printf("%d..%d at %u [%u]\n",diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,diagonal->univdiagonal - chroffset);
  }
  printf("\n");
#endif

  /*  Step 1:  Handle mismatches */
  *completep = false;
  super_path = (List_T) NULL;

  p = middle_path;
  prev_diagonal = (Univdiag_T) List_head(p);
  querystart_for_merge = prev_diagonal->querystart;
  prev_left = prev_diagonal->univdiagonal;
  nmismatches = 0;
  fillin_p = false;

  for (p = List_next(p); p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    if ((left = diagonal->univdiagonal) == prev_left) {
      /* Mismatch */
      nmismatches += (diagonal->querystart - prev_diagonal->queryend - 1); /* This could be an overestimate */
      debug13(printf("We have mismatch or mismatches between %d..%d and %d..%d.  Incrementing mismatches by %d => %d\n",
		     prev_diagonal->querystart,prev_diagonal->queryend,diagonal->querystart,diagonal->queryend,
		     (diagonal->querystart - prev_diagonal->queryend - 1),nmismatches));
      if (diagonal->nmismatches_known_p == false) {
	fillin_p = true;
      }
	     
    } else {
      /* Indel or splice */

      /* Handle previous segment (for prev_left) */
      new_diagonal = Univdiag_new(querystart_for_merge,prev_diagonal->queryend,prev_diagonal->univdiagonal);
      if (fillin_p == true || prev_diagonal->nmismatches_known_p == false) {
	new_diagonal->intscore = 100; /* Positive score allows for many mismatches in indel/splice routines */
      } else {
	new_diagonal->intscore = nmismatches;
      }
      super_path = List_push(super_path,(void *) new_diagonal);

      prev_left = left;
      querystart_for_merge = diagonal->querystart;
      nmismatches = 0;
      fillin_p = false;
    }

    prev_diagonal = diagonal;
  }

  new_diagonal = Univdiag_new(querystart_for_merge,prev_diagonal->queryend,prev_diagonal->univdiagonal);
  if (fillin_p == true || prev_diagonal->nmismatches_known_p == false) {
    new_diagonal->intscore = 100; /* Positive score allows for many mismatches in indel/splice routines */
  } else {
    new_diagonal->intscore = nmismatches;
  }
  super_path = List_push(super_path,(void *) new_diagonal);

  super_path = List_reverse(super_path);

#ifdef DEBUG13
  printf("\n");
  printf("Super diagonals on chrnum %d:\n",chrnum);
  for (p = super_path; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    printf("%d..%d at %u [%u] with %d mismatches\n",
	   diagonal->querystart,diagonal->queryend,diagonal->univdiagonal,
	   diagonal->univdiagonal - chroffset,diagonal->intscore);
  }
  printf("\n");
#endif


  /*  Step 2:  Handle indels and splices */

  p = super_path;
  prev_diagonal = (Univdiag_T) List_head(p);
  prev_left = prev_diagonal->univdiagonal;

  debug13(printf("Diagonal %d..%d for left %u [%u]\n",
		 prev_diagonal->querystart,prev_diagonal->queryend,prev_left,prev_left - chroffset));

  sense_endpoints = Intlist_push(NULL,prev_diagonal->querystart);
  antisense_endpoints = Intlist_push(NULL,prev_diagonal->querystart);

  /* Previously pushed prev_diagonal->intscore, but that is not
     correct (unless there are no mismatches).  Pushing -1 indicates
     that we need to compute the value */
  if (prev_diagonal->intscore == 0) {
    sense_nmismatches = Intlist_push(NULL,0);
    antisense_nmismatches = Intlist_push(NULL,0);
  } else {
    sense_nmismatches = Intlist_push(NULL,-1);
    antisense_nmismatches = Intlist_push(NULL,-1);
  }

  for (p = List_next(p); p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    left = diagonal->univdiagonal;
    assert(left != prev_left); /* Because we already handled mismatches above */

    debug13(printf("Diagonal %d..%d at leftpos %u [%u], diff %d\n",
		   diagonal->querystart,diagonal->queryend,left,left - chroffset,left - prev_left));

    if (left < prev_left) {
      /* Insertion */
      nindels = prev_left - left;
#if 0
      max_mismatches_allowed = (diagonal->querystart - prev_diagonal->queryend - 1);
      debug13(printf("max_mismatches %d = %d - %d - 1\n",max_mismatches_allowed,diagonal->querystart,prev_diagonal->queryend));
      if (prev_diagonal->intscore > 0) {
	max_mismatches_allowed += 1;
      }
      if (diagonal->intscore > 0) {
	max_mismatches_allowed += 1;
      }
#endif
      if ((indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
						      /*left*/prev_left,/*indels*/+nindels,query_compress,
						      prev_diagonal->querystart,diagonal->queryend,querylength,
						      max_mismatches_allowed,/*plusp:true*/true,genestrand)) < 0) {
	sense_junctions = List_push(sense_junctions,NULL);
	antisense_junctions = List_push(antisense_junctions,NULL);
      } else {
	sense_junctions = List_push(sense_junctions,Junction_new_insertion(nindels));
	antisense_junctions = List_push(antisense_junctions,Junction_new_insertion(nindels));
	diagonal->querystart += nindels; /* Needed for subsequent indel computation */
      }

      if ((prev_nmismatches = Intlist_head(sense_nmismatches)) < 0) {
	/* Still need to compute */
	sense_nmismatches = Intlist_push(sense_nmismatches,best_nmismatches_j);
      } else if (Intlist_head(sense_endpoints) != prev_diagonal->querystart) {
	/* Endpoints not quite the same, so need to recompute */
	Intlist_head_set(sense_nmismatches,-1);
	sense_nmismatches = Intlist_push(sense_nmismatches,-1);
      } else {
	Intlist_head_set(sense_nmismatches,best_nmismatches_i /*+ prev_nmismatches*/);
	sense_nmismatches = Intlist_push(sense_nmismatches,best_nmismatches_j);
      }

      if ((prev_nmismatches = Intlist_head(antisense_nmismatches)) < 0) {
	/* Still need to compute */
	antisense_nmismatches = Intlist_push(antisense_nmismatches,best_nmismatches_j);
      } else if (Intlist_head(antisense_endpoints) != prev_diagonal->querystart) {
	/* Endpoints not quite the same, so need to recompute */
	Intlist_head_set(antisense_nmismatches,-1);
	antisense_nmismatches = Intlist_push(antisense_nmismatches,-1);
      } else {
	Intlist_head_set(antisense_nmismatches,best_nmismatches_i /*+ prev_nmismatches*/);
	antisense_nmismatches = Intlist_push(antisense_nmismatches,best_nmismatches_j);
      }

      sense_lefts = Uintlist_push(sense_lefts,prev_left);
      antisense_lefts = Uintlist_push(antisense_lefts,prev_left);

      sense_endpoints = Intlist_push(sense_endpoints,indel_pos);
      antisense_endpoints = Intlist_push(antisense_endpoints,indel_pos);
      debug13(printf("insertion pos in range %d..%d is %d with nmismatches %d+%d\n",
		     prev_diagonal->querystart,diagonal->queryend,indel_pos,best_nmismatches_i,best_nmismatches_j));
      
    } else if (left <= prev_left + max_deletionlen) {
      /* Deletion */
      nindels = left - prev_left;
#if 0
      max_mismatches_allowed = (diagonal->querystart - prev_diagonal->queryend - 1);
      debug13(printf("max_mismatches %d = %d - %d - 1\n",max_mismatches_allowed,diagonal->querystart,prev_diagonal->queryend));
      if (prev_diagonal->intscore > 0) {
	max_mismatches_allowed += 1;
      }
      if (diagonal->intscore > 0) {
	max_mismatches_allowed += 1;
      }
#endif
      if ((indel_pos = Indel_resolve_middle_deletion(&best_nmismatches_i,&best_nmismatches_j,
						     /*left*/prev_left,/*indels*/-nindels,query_compress,
						     prev_diagonal->querystart,diagonal->queryend,querylength,
						     max_mismatches_allowed,/*plusp:true*/true,genestrand)) < 0) {
	sense_junctions = List_push(sense_junctions,NULL);
	antisense_junctions = List_push(antisense_junctions,NULL);
      } else {
	deletionpos = prev_left + indel_pos;
	sense_junctions = List_push(sense_junctions,Junction_new_deletion(nindels,deletionpos));
	antisense_junctions = List_push(antisense_junctions,Junction_new_deletion(nindels,deletionpos));
      }

      if ((prev_nmismatches = Intlist_head(sense_nmismatches)) < 0) {
	/* Still need to compute */
	sense_nmismatches = Intlist_push(sense_nmismatches,best_nmismatches_j);
      } else if (Intlist_head(sense_endpoints) != prev_diagonal->querystart) {
	/* Endpoints not quite the same, so need to recompute */
	Intlist_head_set(sense_nmismatches,-1);
	sense_nmismatches = Intlist_push(sense_nmismatches,-1);
      } else {
	Intlist_head_set(sense_nmismatches,best_nmismatches_i /*+ prev_nmismatches*/);
	sense_nmismatches = Intlist_push(sense_nmismatches,best_nmismatches_j);
      }

      if ((prev_nmismatches = Intlist_head(antisense_nmismatches)) < 0) {
	/* Still need to compute */
	antisense_nmismatches = Intlist_push(antisense_nmismatches,best_nmismatches_j);
      } else if (Intlist_head(antisense_endpoints) != prev_diagonal->querystart) {
	/* Endpoints not quite the same, so need to recompute */
	Intlist_head_set(antisense_nmismatches,-1);
	antisense_nmismatches = Intlist_push(antisense_nmismatches,-1);
      } else {
	Intlist_head_set(antisense_nmismatches,best_nmismatches_i /*+ prev_nmismatches*/);
	antisense_nmismatches = Intlist_push(antisense_nmismatches,best_nmismatches_j);
      }

      sense_lefts = Uintlist_push(sense_lefts,prev_left);
      antisense_lefts = Uintlist_push(antisense_lefts,prev_left);

      sense_endpoints = Intlist_push(sense_endpoints,indel_pos);
      antisense_endpoints = Intlist_push(antisense_endpoints,indel_pos);
      debug13(printf("deletion pos in range %d..%d is %d with nmismatches %d+%d\n",
		     prev_diagonal->querystart,diagonal->queryend,indel_pos,best_nmismatches_i,best_nmismatches_j));
      
    } else {
      /* Splice */
      segmenti_donor_nknown = segmenti_antiacceptor_nknown = 0;
      if (nsplicesites > 0 &&
	  Splicetrie_splicesite_p(prev_left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search(0,nsplicesites,splicesites,prev_left);
	while (j < nsplicesites && splicesites[j] < prev_left + querylength) {
	  if (splicetypes[j] == DONOR) {
	    debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	    segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - prev_left;
	    segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIACCEPTOR) {
	    debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	    segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - prev_left;
	    segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	  }
	  j++;
	}
      }
      segmenti_donor_knownpos[segmenti_donor_nknown] = querylength + 100;
      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength + 100;
	  
      segmentj_acceptor_nknown = segmentj_antidonor_nknown = 0;
      if (nsplicesites > 0 &&
	  Splicetrie_splicesite_p(left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search(0,nsplicesites,splicesites,left);
	while (j < nsplicesites && splicesites[j] < left + querylength) {
	  if (splicetypes[j] == ACCEPTOR) {
	    debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	    segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - left;
	    segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIDONOR) {
	    debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	    segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - left;
	    segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
	  }
	  j++;
	}
      }
      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength + 100;
      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength + 100;
      
      splice_distance = left - prev_left;
#if 0
      max_mismatches_allowed = (diagonal->querystart - prev_diagonal->queryend - 1);
      debug13(printf("max_mismatches %d = %d - %d - 1\n",max_mismatches_allowed,diagonal->querystart,prev_diagonal->queryend));
      if (prev_diagonal->intscore > 0) {
	max_mismatches_allowed += 1;
      }
      if (diagonal->intscore > 0) {
	max_mismatches_allowed += 1;
      }
#endif

      if ((splice_pos = Splice_resolve_sense(&best_knowni_i,&best_knowni_j,&best_nmismatches_i,&best_nmismatches_j,
					     &best_prob_i,&best_prob_j,
					     /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
					     prev_diagonal->querystart,diagonal->queryend+1,querylength,query_compress,
					     segmenti_donor_knownpos,segmentj_acceptor_knownpos,
					     segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					     segmenti_donor_knowni,segmentj_acceptor_knowni,
					     segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
					     segmenti_donor_nknown,segmentj_acceptor_nknown,
					     segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
					     max_mismatches_allowed,plusp,genestrand)) < 0) {
	sense_endpoints = Intlist_push(sense_endpoints,-1); /* Mark as invalid */
	sense_junctions = List_push(sense_junctions,NULL);
      } else if (plusp == true) {
	sense_endpoints = Intlist_push(sense_endpoints,splice_pos);
	sense_junctions = List_push(sense_junctions,Junction_new_splice(splice_distance,SENSE_FORWARD,
									/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
      } else {
	sense_endpoints = Intlist_push(sense_endpoints,splice_pos);
	sense_junctions = List_push(sense_junctions,Junction_new_splice(splice_distance,SENSE_FORWARD,
									/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
      }
      debug13(printf("sense splice_pos in range %d..%d is %d with mismatches %d+%d\n",
		     prev_diagonal->querystart,diagonal->queryend,splice_pos,best_nmismatches_i,best_nmismatches_j));
      if ((prev_nmismatches = Intlist_head(sense_nmismatches)) < 0) {
	/* Still need to compute */
	sense_nmismatches = Intlist_push(sense_nmismatches,best_nmismatches_j);
      } else if (Intlist_head(sense_endpoints) != prev_diagonal->querystart) {
	/* Endpoints not quite the same, so need to recompute */
	Intlist_head_set(sense_nmismatches,-1);
	sense_nmismatches = Intlist_push(sense_nmismatches,-1);
      } else {
	Intlist_head_set(sense_nmismatches,best_nmismatches_i /*+ prev_nmismatches*/);
	sense_nmismatches = Intlist_push(sense_nmismatches,best_nmismatches_j);
      }
      sense_lefts = Uintlist_push(sense_lefts,prev_left);

      if ((splice_pos = Splice_resolve_antisense(&best_knowni_i,&best_knowni_j,&best_nmismatches_i,&best_nmismatches_j,
						 &best_prob_i,&best_prob_j,
						 /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
						 prev_diagonal->querystart,diagonal->queryend+1,querylength,query_compress,
						 segmenti_donor_knownpos,segmentj_acceptor_knownpos,
						 segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
						 segmenti_donor_knowni,segmentj_acceptor_knowni,
						 segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
						 segmenti_donor_nknown,segmentj_acceptor_nknown,
						 segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
						 max_mismatches_allowed,plusp,genestrand)) < 0) {
	antisense_endpoints = Intlist_push(antisense_endpoints,-1); /* Mark as invalid */
	antisense_junctions = List_push(antisense_junctions,NULL);
      } else if (plusp == true) {
	antisense_endpoints = Intlist_push(antisense_endpoints,splice_pos);
	antisense_junctions = List_push(antisense_junctions,Junction_new_splice(splice_distance,SENSE_ANTI,
										/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
      } else {
	antisense_endpoints = Intlist_push(antisense_endpoints,splice_pos);
	antisense_junctions = List_push(antisense_junctions,Junction_new_splice(splice_distance,SENSE_ANTI,
										/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
      }
      debug13(printf("antisense splice_pos in range %d..%d is %d with nmismatches %d+%d\n",
		     prev_diagonal->querystart,diagonal->queryend,splice_pos,best_nmismatches_i,best_nmismatches_j));
      if ((prev_nmismatches = Intlist_head(antisense_nmismatches)) < 0) {
	/* Still need to compute */
	antisense_nmismatches = Intlist_push(antisense_nmismatches,best_nmismatches_j);
      } else if (Intlist_head(antisense_endpoints) != prev_diagonal->querystart) {
	/* Endpoints not quite the same, so need to recompute */
	Intlist_head_set(antisense_nmismatches,-1);
	antisense_nmismatches = Intlist_push(antisense_nmismatches,-1);
      } else {
	Intlist_head_set(antisense_nmismatches,best_nmismatches_i /*+ prev_nmismatches*/);
	antisense_nmismatches = Intlist_push(antisense_nmismatches,best_nmismatches_j);
      }
      antisense_lefts = Uintlist_push(antisense_lefts,prev_left);
    }

    /* Handle previous segment (for prev_left) */
    prev_left = left;
    prev_diagonal = diagonal;
  }

  /* Finish up lists */
  sense_lefts = Uintlist_push(sense_lefts,prev_left);
  antisense_lefts = Uintlist_push(antisense_lefts,prev_left);
  sense_endpoints = Intlist_push(sense_endpoints,prev_diagonal->queryend + 1);
  antisense_endpoints = Intlist_push(antisense_endpoints,prev_diagonal->queryend + 1);


  debug13(printf("After step 2 (indels and splices)\n"));
  debug13(printf("sense (wrt plus): %s\n",Intlist_to_string(sense_endpoints)));
  debug13(printf("antisense (wrt plus): %s\n",Intlist_to_string(antisense_endpoints)));
  debug13(printf("sense nmismatches: %s\n",Intlist_to_string(sense_nmismatches)));
  debug13(printf("antisense nmismatches: %s\n",Intlist_to_string(antisense_nmismatches)));


  /*  Step 3:  Handle ambiguous ends on right */
  right_ambig_sense = (Substring_T) NULL;
  if (circularp[chrnum] == true) {
    /* Skip */

  } else if (right_endpoints_sense == NULL) {
    /* Skip */

  } else if (Intlist_length(right_endpoints_sense) == 1) {
    /* Only one splice on right */
    debug13(printf("Only one sense splice on right, which should have %d and %d mismatches\n",
		   Intlist_head(right_amb_nmismatchesi_sense),Intlist_head(right_amb_nmismatchesj_sense)));
    splice_pos = Intlist_head(right_endpoints_sense);
    queryend = Intlist_head(right_queryends_sense);
    left = Uintlist_head(right_ambcoords_sense) - splice_pos;
    splice_distance = left - prev_left;
    if (plusp == true) {
      donor_prob = Doublelist_head(right_amb_probsi_sense);
      acceptor_prob = Doublelist_head(right_amb_probsj_sense);
    } else {
      acceptor_prob = Doublelist_head(right_amb_probsi_sense);
      donor_prob = Doublelist_head(right_amb_probsj_sense);
    }

    if (Intlist_head(sense_endpoints) != splice_pos) {
      Intlist_head_set(sense_nmismatches,-1);
      sense_nmismatches = Intlist_push(sense_nmismatches,-1);
      Intlist_head_set(sense_endpoints,splice_pos);
    } else {
      /* Only distal nmismatches is reliable */
      /* Intlist_head_set(sense_nmismatches,Intlist_head(right_amb_nmismatchesi_sense)); */
      sense_nmismatches = Intlist_push(sense_nmismatches,Intlist_head(right_amb_nmismatchesj_sense));
    }
    sense_endpoints = Intlist_push(sense_endpoints,queryend);

    sense_lefts = Uintlist_push(sense_lefts,left);
    sense_junctions = List_push(sense_junctions,Junction_new_splice(splice_distance,SENSE_FORWARD,
								    donor_prob,acceptor_prob));

  } else if (Intlist_vary(right_endpoints_sense) == true) {
    /* Skip */
  } else {
    /* Ambiguous substring on right */
    debug13(printf("Ambiguous substring on right\n"));
    splice_pos = Intlist_head(right_endpoints_sense);
    queryend = Intlist_head(right_queryends_sense); /* Should all be the same */

    if (Intlist_head(sense_endpoints) != splice_pos) {
      Intlist_head_set(sense_nmismatches,-1);
      Intlist_head_set(sense_endpoints,splice_pos);
    }
    /* sense_endpoints = Intlist_push(sense_endpoints,queryend); */

    if (plusp == true) {
      right_ambig_sense = Substring_new_ambig_A(/*querystart*/splice_pos,queryend,
						/*splice_pos*/splice_pos,querylength,
						chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
						right_ambcoords_sense,right_amb_knowni_sense,
						right_amb_nmismatchesj_sense,right_amb_probsj_sense,
						/*amb_common_prob*/Doublelist_head(right_amb_probsi_sense),
						/*substring1p*/false);
    } else {
      right_ambig_sense = Substring_new_ambig_D(/*querystart*/querylength - queryend,querylength - splice_pos,
						/*splice_pos*/querylength - splice_pos,querylength,
						chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
						right_ambcoords_sense,right_amb_knowni_sense,
						right_amb_nmismatchesj_sense,right_amb_probsj_sense,
						/*amb_common_prob*/Doublelist_head(right_amb_probsi_sense),
						/*substring1p*/true);
    }
  }

  if (right_ambig_sense != NULL) {
    /* Endpoints end before ambiguous substring */
  } else if (Intlist_head(sense_endpoints) == querylength) {
    /* Last substring already goes to the end */
  } else {
    Intlist_head_set(sense_endpoints,querylength);
    Intlist_head_set(sense_nmismatches,-1); /* Recalculate */
  }


  right_ambig_antisense = (Substring_T) NULL;
  if (circularp[chrnum] == true) {
    /* Skip */

  } else if (right_endpoints_antisense == NULL) {
    /* Skip */

  } else if (Intlist_length(right_endpoints_antisense) == 1) {
    /* Only one splice on right */
    debug13(printf("Only one antisense splice on right, which should have %d and %d mismatches\n",
		   Intlist_head(right_amb_nmismatchesi_antisense),Intlist_head(right_amb_nmismatchesj_antisense)));
    splice_pos = Intlist_head(right_endpoints_antisense);
    queryend = Intlist_head(right_queryends_antisense);
    left = Uintlist_head(right_ambcoords_antisense) - splice_pos;
    splice_distance = left - prev_left;
    if (plusp == true) {
      acceptor_prob = Doublelist_head(right_amb_probsi_antisense);
      donor_prob = Doublelist_head(right_amb_probsj_antisense);
    } else {
      donor_prob = Doublelist_head(right_amb_probsi_antisense);
      acceptor_prob = Doublelist_head(right_amb_probsj_antisense);
    }

    if (Intlist_head(antisense_endpoints) != splice_pos) {
      Intlist_head_set(antisense_nmismatches,-1);
      antisense_nmismatches = Intlist_push(antisense_nmismatches,-1);
      Intlist_head_set(antisense_endpoints,splice_pos);
    } else {
      /* Only distal nmismatches is reliable */
      /* Intlist_head_set(antisense_nmismatches,Intlist_head(right_amb_nmismatchesi_antisense)); */
      antisense_nmismatches = Intlist_push(antisense_nmismatches,Intlist_head(right_amb_nmismatchesj_antisense));
    }
    antisense_endpoints = Intlist_push(antisense_endpoints,queryend);

    antisense_lefts = Uintlist_push(antisense_lefts,left);
    antisense_junctions = List_push(antisense_junctions,Junction_new_splice(splice_distance,SENSE_ANTI,
									    donor_prob,acceptor_prob));

  } else if (Intlist_vary(right_endpoints_antisense) == true) {
    /* Skip */
  } else {
    /* Ambiguous substring on right */
    debug13(printf("Ambiguous substring on right\n"));
    splice_pos = Intlist_head(right_endpoints_antisense);
    queryend = Intlist_head(right_queryends_antisense); /* Should all be the same */

    if (Intlist_head(antisense_endpoints) != splice_pos) {
      Intlist_head_set(antisense_nmismatches,-1);
      Intlist_head_set(antisense_endpoints,splice_pos);
    }
    /* antisense_endpoints = Intlist_push(antisense_endpoints,queryend); */

    if (plusp == true) {
      right_ambig_antisense = Substring_new_ambig_D(/*querystart*/splice_pos,queryend,
						    /*splice_pos*/splice_pos,querylength,
						    chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
						    right_ambcoords_antisense,right_amb_knowni_antisense,
						    right_amb_nmismatchesj_antisense,right_amb_probsj_antisense,
						    /*amb_common_prob*/Doublelist_head(right_amb_probsi_antisense),
						    /*substring1p*/false);
    } else {
      right_ambig_antisense = Substring_new_ambig_A(/*querystart*/querylength - queryend,querylength - splice_pos,
						    /*splice_pos*/querylength - splice_pos,querylength,
						    chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
						    right_ambcoords_antisense,right_amb_knowni_antisense,
						    right_amb_nmismatchesj_antisense,right_amb_probsj_antisense,
						    /*amb_common_prob*/Doublelist_head(right_amb_probsi_antisense),
						    /*substring1p*/true);
    }
  }

  if (right_ambig_antisense != NULL) {
    /* Endpoints end before ambiguous substring */
  } else if (Intlist_head(antisense_endpoints) == querylength) {
    /* Last substring already goes to the end */
  } else {
    Intlist_head_set(antisense_endpoints,querylength);
    Intlist_head_set(antisense_nmismatches,-1); /* Recalculate */
  }


  debug13(printf("After step 3 (ambiguous ends on right)\n"));
  debug13(printf("sense (wrt plus): %s\n",Intlist_to_string(sense_endpoints)));
  debug13(printf("antisense (wrt plus): %s\n",Intlist_to_string(antisense_endpoints)));
  debug13(printf("sense nmismatches: %s\n",Intlist_to_string(sense_nmismatches)));
  debug13(printf("antisense nmismatches: %s\n",Intlist_to_string(antisense_nmismatches)));

  /*  Step 4:  Reverse sense and antisense alignments */
  sense_nmismatches = Intlist_reverse(sense_nmismatches);
  antisense_nmismatches = Intlist_reverse(antisense_nmismatches);
  sense_lefts = Uintlist_reverse(sense_lefts);
  antisense_lefts = Uintlist_reverse(antisense_lefts);
  sense_endpoints = Intlist_reverse(sense_endpoints);
  antisense_endpoints = Intlist_reverse(antisense_endpoints);
  sense_junctions = List_reverse(sense_junctions);
  antisense_junctions = List_reverse(antisense_junctions);

  debug13(printf("After step 4 (reverse alignments)\n"));
  debug13(printf("sense (wrt plus): %s\n",Intlist_to_string(sense_endpoints)));
  debug13(printf("antisense (wrt plus): %s\n",Intlist_to_string(antisense_endpoints)));
  debug13(printf("sense nmismatches: %s\n",Intlist_to_string(sense_nmismatches)));
  debug13(printf("antisense nmismatches: %s\n",Intlist_to_string(antisense_nmismatches)));


  /*  Step 5:  Handle ambiguous ends on left */
  left_ambig_sense = (Substring_T) NULL;
  if (circularp[chrnum] == true) {
    /* Skip */
    debug13(printf("Sense: Chrnum %d is circular, so not computing mismatches\n",chrnum));

  } else if (left_endpoints_sense == NULL) {
    /* Skip, but extend leftward */
    debug13(printf("Sense: Skip, but extend leftward\n"));
    if (Intlist_head(sense_endpoints) > 0) {
      querystart = Intlist_head(sense_endpoints);
      if ((max_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/Uintlist_head(sense_lefts),
							      /*pos5*/0,/*pos3*/querystart,plusp,genestrand)) > 0) {
	Intlist_head_set(sense_endpoints,querystart - max_leftward);
	Intlist_head_set(sense_nmismatches,-1);
      } else if ((max_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/Uintlist_head(sense_lefts),
								     /*pos5*/0,/*pos3*/querystart-1,plusp,genestrand)) > 0) {
	Intlist_head_set(sense_endpoints,querystart - max_leftward - 1);
	Intlist_head_set(sense_nmismatches,-1);
      } else {
	/* Keep value as querystart */
      }
    }

  } else if (Intlist_length(left_endpoints_sense) == 1) {
    /* Only one splice on left */
    debug13(printf("Only one sense splice on left, which should have %d and %d mismatches, plusp %d\n",
		   Intlist_head(left_amb_nmismatchesi_sense),Intlist_head(left_amb_nmismatchesj_sense),plusp));
    prev_left = Uintlist_head(sense_lefts);
    splice_pos = Intlist_head(left_endpoints_sense);
    querystart = Intlist_head(left_querystarts_sense);
    left = Uintlist_head(left_ambcoords_sense) - splice_pos;
    splice_distance = prev_left - left;
    if (plusp == true) {
      donor_prob = Doublelist_head(left_amb_probsi_sense);
      acceptor_prob = Doublelist_head(left_amb_probsj_sense);
    } else {
      acceptor_prob = Doublelist_head(left_amb_probsi_sense);
      donor_prob = Doublelist_head(left_amb_probsj_sense);
    }

    if (Intlist_head(sense_endpoints) != splice_pos) {
      Intlist_head_set(sense_nmismatches,-1);
      sense_nmismatches = Intlist_push(sense_nmismatches,-1);
      Intlist_head_set(sense_endpoints,splice_pos);
    } else {
      /* Only distal nmismatches is reliable */
      /* Intlist_head_set(sense_nmismatches,Intlist_head(left_amb_nmismatchesj_sense)); */
      sense_nmismatches = Intlist_push(sense_nmismatches,Intlist_head(left_amb_nmismatchesi_sense));
    }
    sense_endpoints = Intlist_push(sense_endpoints,querystart);

    sense_lefts = Uintlist_push(sense_lefts,left);
    sense_junctions = List_push(sense_junctions,Junction_new_splice(splice_distance,SENSE_FORWARD,
								    donor_prob,acceptor_prob));

  } else if (Intlist_vary(left_endpoints_sense) == true) {
    /* Skip, but extend leftward */
    querystart = Intlist_head(sense_endpoints);
    if ((max_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/Uintlist_head(sense_lefts),
							    /*pos5*/0,/*pos3*/querystart,plusp,genestrand)) > 0) {
      Intlist_head_set(sense_endpoints,querystart - max_leftward);
      Intlist_head_set(sense_nmismatches,-1);
    } else if ((max_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/Uintlist_head(sense_lefts),
								   /*pos5*/0,/*pos3*/querystart-1,plusp,genestrand)) > 0) {
      Intlist_head_set(sense_endpoints,querystart - max_leftward - 1);
      Intlist_head_set(sense_nmismatches,-1);
    } else {
      /* Keep value as querystart */
    }

  } else {
    /* Ambiguous substring on left */
    debug13(printf("Ambiguous substring on left\n"));
    splice_pos = Intlist_head(left_endpoints_sense);
    querystart = Intlist_head(left_querystarts_sense); /* Should all be the same */

    if (Intlist_head(sense_endpoints) != splice_pos) {
      Intlist_head_set(sense_nmismatches,-1);
      Intlist_head_set(sense_endpoints,splice_pos);
    }
    /* sense_endpoints = Intlist_push(sense_endpoints,querystart); */

    if (plusp == true) {
      left_ambig_sense = Substring_new_ambig_D(querystart,/*queryend*/splice_pos,
					       /*splice_pos*/splice_pos,querylength,
					       chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
					       left_ambcoords_sense,left_amb_knowni_sense,
					       left_amb_nmismatchesi_sense,left_amb_probsi_sense,
					       /*amb_common_prob*/Doublelist_head(left_amb_probsj_sense),
					       /*substring1p*/true);
    } else {
      left_ambig_sense = Substring_new_ambig_A(querylength - splice_pos,/*queryend*/querylength - querystart,
					       /*splice_pos*/querylength - splice_pos,querylength,
					       chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
					       left_ambcoords_sense,left_amb_knowni_sense,
					       left_amb_nmismatchesi_sense,left_amb_probsi_sense,
					       /*amb_common_prob*/Doublelist_head(left_amb_probsj_sense),
					       /*substring1p*/false);
    }
  }

  if (left_ambig_sense != NULL) {
    /* Endpoints begin after ambiguous substring */
    debug13(printf("Sense: Endpoints begin after ambiguous substring\n"));
  } else if (Intlist_head(sense_endpoints) == 0) {
    /* First substring already goes to the beginning */
    debug13(printf("Sense: First substring already goes to the beginning\n"));
  } else {
    Intlist_head_set(sense_endpoints,0);
    Intlist_head_set(sense_nmismatches,-1); /* Recalculate */
  }


  left_ambig_antisense = (Substring_T) NULL;
  if (circularp[chrnum] == true) {
    /* Skip */
    debug13(printf("Antisense: Chrnum %d is circular, so not computing mismatches\n",chrnum));

  } else if (left_endpoints_antisense == NULL) {
    /* Skip, but extend leftward */
    debug13(printf("Antisense: Skip, but extend leftward\n"));
    querystart = Intlist_head(antisense_endpoints);
    if (Intlist_head(antisense_endpoints) > 0) {
      if ((max_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/Uintlist_head(antisense_lefts),
							      /*pos5*/0,/*pos3*/querystart,plusp,genestrand)) > 0) {
	Intlist_head_set(antisense_endpoints,querystart - max_leftward);
	Intlist_head_set(antisense_nmismatches,-1);
      } else if ((max_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/Uintlist_head(antisense_lefts),
								     /*pos5*/0,/*pos3*/querystart-1,plusp,genestrand)) > 0) {
	Intlist_head_set(antisense_endpoints,querystart - max_leftward - 1);
	Intlist_head_set(antisense_nmismatches,-1);
      } else {
	/* Keep value as querystart */
      }
    }

  } else if (Intlist_length(left_endpoints_antisense) == 1) {
    /* Only one splice on left */
    debug13(printf("Only one antisense splice on left, which should have %d and %d mismatches\n",
		   Intlist_head(left_amb_nmismatchesi_antisense),Intlist_head(left_amb_nmismatchesj_antisense)));
    prev_left = Uintlist_head(antisense_lefts);
    splice_pos = Intlist_head(left_endpoints_antisense);
    querystart = Intlist_head(left_querystarts_antisense);
    left = Uintlist_head(left_ambcoords_antisense) - splice_pos;
    splice_distance = prev_left - left;
    if (plusp == true) {
      acceptor_prob = Doublelist_head(left_amb_probsi_antisense);
      donor_prob = Doublelist_head(left_amb_probsj_antisense);
    } else {
      donor_prob = Doublelist_head(left_amb_probsi_antisense);
      acceptor_prob = Doublelist_head(left_amb_probsj_antisense);
    }

    if (Intlist_head(antisense_endpoints) != splice_pos) {
      Intlist_head_set(antisense_nmismatches,-1);
      antisense_nmismatches = Intlist_push(antisense_nmismatches,-1);
      Intlist_head_set(antisense_endpoints,splice_pos);
    } else {
      /* Only distal nmismatches is reliable */
      /* Intlist_head_set(antisense_nmismatches,Intlist_head(left_amb_nmismatchesj_antisense)); */
      antisense_nmismatches = Intlist_push(antisense_nmismatches,Intlist_head(left_amb_nmismatchesi_antisense));
    }
    antisense_endpoints = Intlist_push(antisense_endpoints,querystart);

    antisense_lefts = Uintlist_push(antisense_lefts,left);
    antisense_junctions = List_push(antisense_junctions,Junction_new_splice(splice_distance,SENSE_ANTI,
									    donor_prob,acceptor_prob));

  } else if (Intlist_vary(left_endpoints_antisense) == true) {
    /* Skip, but extend leftward */
    querystart = Intlist_head(antisense_endpoints);
    if ((max_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/Uintlist_head(antisense_lefts),
							    /*pos5*/0,/*pos3*/querystart,plusp,genestrand)) > 0) {
      Intlist_head_set(antisense_endpoints,querystart - max_leftward);
      Intlist_head_set(antisense_endpoints,-1);
    } else if ((max_leftward = Genome_consecutive_matches_leftward(query_compress,/*left*/Uintlist_head(antisense_lefts),
								   /*pos5*/0,/*pos3*/querystart-1,plusp,genestrand)) > 0) {
      Intlist_head_set(antisense_endpoints,querystart - max_leftward - 1);
      Intlist_head_set(antisense_endpoints,-1);
    } else {
      /* Keep value as querystart */
    }

  } else {
    /* Ambiguous substring on left */
    debug13(printf("Ambiguous substring on left\n"));
    splice_pos = Intlist_head(left_endpoints_antisense);
    querystart = Intlist_head(left_querystarts_antisense); /* Should all be the same */

    if (Intlist_head(antisense_endpoints) != splice_pos) {
      Intlist_head_set(antisense_nmismatches,-1);
      Intlist_head_set(antisense_endpoints,splice_pos);
    }
    /* antisense_endpoints = Intlist_push(antisense_endpoints,querystart); */

    if (plusp == true) {
      left_ambig_antisense = Substring_new_ambig_A(querystart,/*queryend*/splice_pos,
						   /*splice_pos*/splice_pos,querylength,
						   chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
						   left_ambcoords_antisense,left_amb_knowni_antisense,
						   left_amb_nmismatchesi_antisense,left_amb_probsi_antisense,
						   /*amb_common_prob*/Doublelist_head(left_amb_probsj_antisense),
						   /*substring1p*/true);
    } else {
      left_ambig_antisense = Substring_new_ambig_D(querylength - splice_pos,/*queryend*/querylength - querystart,
						   /*splice_pos*/querylength - splice_pos,querylength,
						   chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
						   left_ambcoords_antisense,left_amb_knowni_antisense,
						   left_amb_nmismatchesi_antisense,left_amb_probsi_antisense,
						   /*amb_common_prob*/Doublelist_head(left_amb_probsj_antisense),
						   /*substring1p*/false);
    }
  }

  if (left_ambig_antisense != NULL) {
    /* Endpoints begin after ambiguous substring */
    debug13(printf("Antisense: Endpoints begin after ambiguous substring\n"));
  } else if (Intlist_head(antisense_endpoints) == 0) {
    /* First substring already goes to the beginning */
    debug13(printf("Antisense: First substring already goes to the beginning\n"));
  } else {
    Intlist_head_set(antisense_endpoints,0);
    Intlist_head_set(antisense_nmismatches,-1); /* Recalculate */
  }


  debug13(printf("After step 5 (ambiguous ends on left)\n"));
  debug13(printf("sense (wrt plus): %s\n",Intlist_to_string(sense_endpoints)));
  debug13(printf("antisense (wrt plus): %s\n",Intlist_to_string(antisense_endpoints)));
  debug13(printf("sense nmismatches: %s\n",Intlist_to_string(sense_nmismatches)));
  debug13(printf("antisense nmismatches: %s\n",Intlist_to_string(antisense_nmismatches)));

#ifdef DEBUG13
  printf("Sense junctions\n");
  for (p = sense_junctions; p != NULL; p = List_next(p)) {
    Junction_print(List_head(p));
  }
  printf("\n");
  printf("Antisense junctions\n");
  for (p = antisense_junctions; p != NULL; p = List_next(p)) {
    Junction_print(List_head(p));
  }
  printf("\n");
#endif


  /* Need to rely on probability filtering in splice.c to get correct
     results for sense and antisense */
  sense_acceptable_p = endpoints_acceptable_p(&sense_intronp,sense_junctions,sense_endpoints);
  antisense_acceptable_p = endpoints_acceptable_p(&antisense_intronp,antisense_junctions,
						  antisense_endpoints);
  if (sense_acceptable_p == true && antisense_acceptable_p == true) {
    if (sense_intronp == true || right_ambig_sense != NULL || left_ambig_sense != NULL) {
      sense_sensedir = SENSE_FORWARD;
    } else {
      sense_sensedir = SENSE_NULL;
    }
    if (antisense_intronp == true || right_ambig_antisense != NULL || left_ambig_antisense != NULL) {
      antisense_sensedir = SENSE_ANTI;
    } else {
      antisense_sensedir = SENSE_NULL;
    }

    if (sense_sensedir == SENSE_NULL && antisense_sensedir == SENSE_NULL) {
      /* Create just one hit */
      if ((hit = Stage3end_new_substrings(&(*found_score),sense_endpoints,sense_lefts,
					  sense_nmismatches,sense_junctions,querylength,query_compress,
					  /*right_ambig*/NULL,/*left_ambig*/NULL,plusp,genestrand,/*sensedir*/SENSE_NULL,
					  chrnum,chroffset,chrhigh,chrlength,/*sarrayp*/true)) == NULL) {
	Substring_free(&right_ambig_sense);
	Substring_free(&left_ambig_sense);
	/* Junction_gc(&sense_junctions); -- Done by Stage3end_new_substrings */
	Substring_free(&right_ambig_antisense);
	Substring_free(&left_ambig_antisense);
      } else {
	if (Stage3end_substrings_querystart(hit) < 8 &&
	    Stage3end_substrings_queryend(hit) >= querylength - 8) {
	  *completep = true;
	}
	hits = List_push(hits,(void *) hit);
      }
      Junction_gc(&antisense_junctions);

    } else {
      /* Create just both sense and antisense hits */
      if ((hit = Stage3end_new_substrings(&(*found_score),sense_endpoints,sense_lefts,
					  sense_nmismatches,sense_junctions,querylength,query_compress,
					  right_ambig_sense,left_ambig_sense,plusp,genestrand,sense_sensedir,
					  chrnum,chroffset,chrhigh,chrlength,/*sarrayp*/true)) == NULL) {
	Substring_free(&right_ambig_sense);
	Substring_free(&left_ambig_sense);
	/* Junction_gc(&sense_junctions); -- Done by Stage3end_new_substrings */
      } else {
	if (Stage3end_substrings_querystart(hit) < 8 &&
	    Stage3end_substrings_queryend(hit) >= querylength - 8) {
	  *completep = true;
	}
	hits = List_push(hits,(void *) hit);
      }

      if ((hit = Stage3end_new_substrings(&(*found_score),antisense_endpoints,antisense_lefts,
					  antisense_nmismatches,antisense_junctions,querylength,query_compress,
					  right_ambig_antisense,left_ambig_antisense,plusp,genestrand,antisense_sensedir,
					  chrnum,chroffset,chrhigh,chrlength,/*sarrayp*/true)) == NULL) {
	Substring_free(&right_ambig_antisense);
	Substring_free(&left_ambig_antisense);
	/* Junction_gc(&antisense_junctions); -- Done by Stage3end_new_substrings */
      } else {
	if (Stage3end_substrings_querystart(hit) < 8 &&
	    Stage3end_substrings_queryend(hit) >= querylength - 8) {
	  *completep = true;
	}
	hits = List_push(hits,(void *) hit);
      }
    }
    
  } else if (sense_acceptable_p == true) {
    if (sense_intronp == true || right_ambig_sense != NULL || left_ambig_sense != NULL) {
      sensedir = SENSE_FORWARD;
    } else {
      sensedir = SENSE_NULL;
    }
    if ((hit = Stage3end_new_substrings(&(*found_score),sense_endpoints,sense_lefts,
					sense_nmismatches,sense_junctions,querylength,query_compress,
					right_ambig_sense,left_ambig_sense,plusp,genestrand,sensedir,
					chrnum,chroffset,chrhigh,chrlength,/*sarrayp*/true)) == NULL) {
      Substring_free(&right_ambig_sense);
      Substring_free(&left_ambig_sense);
      /* Junction_gc(&sense_junctions); -- Done by Stage3end_new_substrings */
    } else {
      if (Stage3end_substrings_querystart(hit) < 8 &&
	  Stage3end_substrings_queryend(hit) >= querylength - 8) {
	*completep = true;
      }
      hits = List_push(hits,(void *) hit);
    }

    Substring_free(&right_ambig_antisense);
    Substring_free(&left_ambig_antisense);
    Junction_gc(&antisense_junctions);

  } else if (antisense_acceptable_p == true) {
    if (antisense_intronp == true || right_ambig_antisense != NULL || left_ambig_antisense != NULL) {
      sensedir = SENSE_ANTI;
    } else {
      sensedir = SENSE_NULL;
    }
    if ((hit = Stage3end_new_substrings(&(*found_score),antisense_endpoints,antisense_lefts,
					antisense_nmismatches,antisense_junctions,querylength,query_compress,
					right_ambig_antisense,left_ambig_antisense,plusp,genestrand,sensedir,
					chrnum,chroffset,chrhigh,chrlength,/*sarrayp*/true)) == NULL) {
      Substring_free(&right_ambig_antisense);
      Substring_free(&left_ambig_antisense);
      /* Junction_gc(&antisense_junctions); -- Done by Stage3end_new_substrings */
    } else {
      if (Stage3end_substrings_querystart(hit) < 8 &&
	  Stage3end_substrings_queryend(hit) >= querylength - 8) {
	*completep = true;
      }
      hits = List_push(hits,(void *) hit);
    }

    Substring_free(&right_ambig_sense);
    Substring_free(&left_ambig_sense);
    Junction_gc(&sense_junctions);

  } else {
    /* Neither set of junctions/endpoints works */
    Substring_free(&right_ambig_sense);
    Substring_free(&left_ambig_sense);
    Substring_free(&right_ambig_antisense);
    Substring_free(&left_ambig_antisense);

    Junction_gc(&sense_junctions);
    Junction_gc(&antisense_junctions);
  }


  Intlist_free(&sense_nmismatches);
  Intlist_free(&antisense_nmismatches);
  Uintlist_free(&sense_lefts);
  Uintlist_free(&antisense_lefts);
  Intlist_free(&sense_endpoints);
  Intlist_free(&antisense_endpoints);
  
  for (p = super_path; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    Univdiag_free(&diagonal);
  }
  List_free(&super_path);

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(segmenti_donor_knownpos);
    FREEA(segmentj_acceptor_knownpos);
    FREEA(segmentj_antidonor_knownpos);
    FREEA(segmenti_antiacceptor_knownpos);
    FREEA(segmenti_donor_knowni);
    FREEA(segmentj_acceptor_knowni);
    FREEA(segmentj_antidonor_knowni);
    FREEA(segmenti_antiacceptor_knowni);
  } else {
    FREE(segmenti_donor_knownpos);
    FREE(segmentj_acceptor_knownpos);
    FREE(segmentj_antidonor_knownpos);
    FREE(segmenti_antiacceptor_knownpos);
    FREE(segmenti_donor_knowni);
    FREE(segmentj_acceptor_knowni);
    FREE(segmentj_antidonor_knowni);
    FREE(segmenti_antiacceptor_knowni);
  }
#else
  FREE(segmenti_donor_knownpos);
  FREE(segmentj_acceptor_knownpos);
  FREE(segmentj_antidonor_knownpos);
  FREE(segmenti_antiacceptor_knownpos);
  FREE(segmenti_donor_knowni);
  FREE(segmentj_acceptor_knowni);
  FREE(segmentj_antidonor_knowni);
  FREE(segmenti_antiacceptor_knowni);
#endif

  return hits;
}




List_T
Sarray_search_greedy (int *found_score, char *queryuc_ptr, char *queryrc, int querylength,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev, 
		      int nmisses_allowed, int genestrand) {
  List_T hits = NULL;
  List_T plus_set, minus_set, p;
  List_T rightward_set_plus = NULL, leftward_set_plus = NULL, rightward_set_minus = NULL, leftward_set_minus = NULL;
  Elt_T best_plus_elt, best_minus_elt, elt, *plus_elt_array, *minus_elt_array;
  UINT4 best_plus_nmatches, best_minus_nmatches, nmatches;
  Sarrayptr_T initptr, finalptr;
  bool successp, completep;
  int plus_querypos, minus_querypos;
  int i;
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh, left;
  Chrpos_T chrlength;
  T plus_sarray, minus_sarray;
  char *plus_conversion, *minus_conversion;

  int nseeds_plus, nseeds_minus;
  int *scores_plus = NULL, *scores_minus = NULL;
  int niter, best_plus_i, best_minus_i, nplus, nminus;
  List_T *middle_path_plus = NULL, *right_paths_plus = NULL, *left_paths_plus = NULL,
    *middle_path_minus = NULL, *right_paths_minus = NULL, *left_paths_minus = NULL;
  Univdiag_T *middle_diagonals_plus = NULL, *middle_diagonals_minus = NULL;
  List_T *best_right_diagonals_plus = NULL, *best_left_diagonals_plus = NULL,
    *all_right_diagonals_plus = NULL, *all_left_diagonals_plus = NULL,
    *fillin_diagonals_plus = NULL, *fillin_diagonals_minus = NULL,
    *best_right_diagonals_minus = NULL, *best_left_diagonals_minus = NULL,
    *all_right_diagonals_minus = NULL, *all_left_diagonals_minus = NULL;

  Intlist_T right_endpoints_sense, right_endpoints_antisense,
    left_endpoints_sense, left_endpoints_antisense;
  Intlist_T right_queryends_sense, right_queryends_antisense,
    left_querystarts_sense, left_querystarts_antisense;
  Uintlist_T right_ambcoords_sense, right_ambcoords_antisense,
    left_ambcoords_sense, left_ambcoords_antisense;
  Intlist_T right_amb_knowni_sense, right_amb_knowni_antisense,
    left_amb_knowni_sense, left_amb_knowni_antisense;
  Intlist_T right_amb_nmismatchesi_sense, right_amb_nmismatchesi_antisense,
    right_amb_nmismatchesj_sense, right_amb_nmismatchesj_antisense,
    left_amb_nmismatchesi_sense, left_amb_nmismatchesi_antisense,
    left_amb_nmismatchesj_sense, left_amb_nmismatchesj_antisense;
  Doublelist_T right_amb_probsi_sense, right_amb_probsi_antisense,
    right_amb_probsj_sense, right_amb_probsj_antisense,
    left_amb_probsi_sense, left_amb_probsi_antisense,
    left_amb_probsj_sense, left_amb_probsj_antisense;

  List_T diagonal_path;
  Univdiag_T diagonal;
#if 0
  bool *coveredp;
  Chrpos_T **mappings;
  int *npositions;
  Oligoindex_T oligoindex;
  int indexsize;
  int best_score;
#endif


  debug(printf("\nStarting Sarray_search_greedy with querylength %d and indexsize %d and nmisses_allowed %d, genestrand %d\n",
	       querylength,sarray_fwd->indexsize,nmisses_allowed,genestrand));
  if (nmisses_allowed < 0) {
    nmisses_allowed = 0;
#if 0
  } else {
    /* It is possible that this makes GSNAP too slow */
    nmisses_allowed = querylength;
#endif
  }
  *found_score = querylength;

  if (genestrand == +2) {
    plus_conversion = conversion_rev;
    minus_conversion = conversion_fwd;
    plus_sarray = sarray_rev;
    minus_sarray = sarray_fwd;
  } else {
    plus_conversion = conversion_fwd;
    minus_conversion = conversion_rev;
    plus_sarray = sarray_fwd;
    minus_sarray = sarray_rev;
  }


  /* I.  Race from plus and minus start to end */
  plus_set = minus_set = (List_T) NULL;
  best_plus_nmatches = best_minus_nmatches = 0;
  best_plus_elt = best_minus_elt = (Elt_T) NULL;
  plus_querypos = 0;
  minus_querypos = 0;
  niter = 0;
  while (niter < nmisses_allowed && plus_querypos < querylength && minus_querypos < querylength) {
    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryuc_ptr[plus_querypos]),
		  querylength - plus_querypos,/*queryoffset*/plus_querypos,
		  query_compress_fwd,plus_sarray,/*plusp*/true,genestrand,plus_conversion);
    elt = Elt_new(plus_querypos,nmatches,initptr,finalptr,/*temporaryp*/false);
    if (nmatches > best_plus_nmatches && elt->nptr <= MAX_HITS_FOR_BEST_ELT) {
      best_plus_elt = elt;
      best_plus_nmatches = nmatches;
      best_plus_i = niter;
    }
    plus_set = List_push(plus_set,elt);
    plus_querypos += nmatches;
    plus_querypos += 1;		/* To skip the presumed mismatch */

    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryrc[minus_querypos]),
		  querylength - minus_querypos,/*queryoffset*/minus_querypos,
		  query_compress_rev,minus_sarray,/*plusp*/false,genestrand,minus_conversion);
    elt = Elt_new(minus_querypos,nmatches,initptr,finalptr,/*temporaryp*/false);
    if (nmatches > best_minus_nmatches && elt->nptr < MAX_HITS_FOR_BEST_ELT) {
      best_minus_elt = elt;
      best_minus_nmatches = nmatches;
      best_minus_i = niter;
    }
    minus_set = List_push(minus_set,elt);
    minus_querypos += nmatches;
    minus_querypos += 1;		/* To skip the presumed mismatch */

    niter++;
  }

#ifdef DEBUG
  printf("niter %d vs %d allowed, plus 0..%d, minus 0..%d\n",niter,nmisses_allowed,plus_querypos,minus_querypos);
  if (best_plus_elt != NULL) {
    printf("best plus %d..%d (SA %u+%d)\n",
	   best_plus_elt->querystart,best_plus_elt->queryend,best_plus_elt->initptr,best_plus_elt->finalptr - best_plus_elt->initptr);
  }
  if (best_minus_elt != NULL) {
    printf("best minus %d..%d (SA %u+%d)\n",
	 best_minus_elt->querystart,best_minus_elt->queryend,best_minus_elt->initptr,best_minus_elt->finalptr - best_minus_elt->initptr);
  }
  printf("plus set (positions not yet filled):\n");
  for (p = plus_set; p != NULL; p = List_next(p)) {
    Elt_dump((Elt_T) List_head(p));
  }
  printf("\n");
  printf("minus set (positions not yet filled):\n");
  for (p = minus_set; p != NULL; p = List_next(p)) {
    Elt_dump((Elt_T) List_head(p));
  }
#endif

  if (plus_querypos < querylength) {
    debug(printf("Plus: could not find large pieces\n"));
    nseeds_plus = 0;

  } else if (best_plus_elt == NULL) {
    debug(printf("Plus: No best elt\n"));
    nseeds_plus = 0;

  } else {
    Elt_fill_positions_all(best_plus_elt,plus_sarray);
    if (best_plus_elt->npositions == 0) {
      /* Could happen if there are too many positions */
      debug(printf("Plus: Best elt has no positions\n"));
      nseeds_plus = 0;

    } else {
      plus_set = List_reverse(plus_set);
      plus_elt_array = (Elt_T *) List_to_array_n(&nplus,plus_set);

#ifdef DEBUG
      printf("LEFT\n");
      for (i = 0; i < best_plus_i; i++) {
	Elt_dump(plus_elt_array[i]);
      }
      printf("MIDDLE\n");
      Elt_dump(plus_elt_array[best_plus_i]);
      printf("RIGHT\n");
      for (i = best_plus_i + 1; i < nplus; i++) {
	Elt_dump(plus_elt_array[i]);
      }
#endif

      nseeds_plus = best_plus_elt->npositions;
      scores_plus = (int *) MALLOC(nseeds_plus*sizeof(int));
      /* Assigned only if score is high */
      middle_path_plus = (List_T *) CALLOC(nseeds_plus,sizeof(List_T));
      right_paths_plus = (List_T *) CALLOC(nseeds_plus,sizeof(List_T));
      left_paths_plus = (List_T *) CALLOC(nseeds_plus,sizeof(List_T));

      middle_diagonals_plus = (Univdiag_T *) MALLOC(nseeds_plus*sizeof(Univdiag_T));
      best_right_diagonals_plus = (List_T *) MALLOC(nseeds_plus*sizeof(List_T));
      best_left_diagonals_plus = (List_T *) MALLOC(nseeds_plus*sizeof(List_T));
      all_right_diagonals_plus = (List_T *) MALLOC(nseeds_plus*sizeof(List_T));
      all_left_diagonals_plus = (List_T *) MALLOC(nseeds_plus*sizeof(List_T));
      fillin_diagonals_plus = (List_T *) CALLOC(nseeds_plus,sizeof(List_T));

      chrnum = 1;
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,/*chrnum*/1,circular_typeint);
      for (i = 0; i < nseeds_plus; i++) {
	left = best_plus_elt->positions[i];
	if (left > chrhigh) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  /* *chrhigh += 1U; */
	}
	/* May not want to solve for best_right_diagonals and best_left_diagonals.  Use oligoindex instead. */
	scores_plus[i] = get_diagonals(&(middle_diagonals_plus[i]),
				       &(best_right_diagonals_plus[i]),&(best_left_diagonals_plus[i]),
				       &(all_right_diagonals_plus[i]),&(all_left_diagonals_plus[i]),
				       plus_sarray,/*queryptr*/queryuc_ptr,querylength,query_compress_fwd,
				       chroffset,chrhigh,/*goal*/left,plus_elt_array,
				       best_plus_i,nplus,/*plusp*/true,genestrand,plus_conversion);
	debug(printf("Got plus score %d\n",scores_plus[i]));
      }

      FREE(plus_elt_array);
    }
  }

  if (minus_querypos < querylength) {
    debug(printf("Minus: Could not find large pieces\n"));
    nseeds_minus = 0;
    
  } else if (best_minus_elt == NULL) {
    debug(printf("Minus: No best elt\n"));
    nseeds_minus = 0;

  } else {
    Elt_fill_positions_all(best_minus_elt,minus_sarray);
    if (best_minus_elt->npositions == 0) {
      /* Could happen if there are too many positions */
      debug(printf("Minus: Best elt has no positions\n"));
      nseeds_minus = 0;

    } else {
      minus_set = List_reverse(minus_set);
      minus_elt_array = (Elt_T *) List_to_array_n(&nminus,minus_set);

#ifdef DEBUG
      printf("LEFT\n");
      for (i = 0; i < best_minus_i; i++) {
	Elt_dump(minus_elt_array[i]);
      }
      printf("MIDDLE\n");
      Elt_dump(minus_elt_array[best_minus_i]);
      printf("RIGHT\n");
      for (i = best_minus_i + 1; i < nminus; i++) {
	Elt_dump(minus_elt_array[i]);
      }
#endif

      nseeds_minus = best_minus_elt->npositions;
      scores_minus = (int *) MALLOC(nseeds_minus*sizeof(int));
      /* Assigned only if score is high */
      middle_path_minus = (List_T *) CALLOC(nseeds_minus,sizeof(List_T));
      right_paths_minus = (List_T *) CALLOC(nseeds_minus,sizeof(List_T));
      left_paths_minus = (List_T *) CALLOC(nseeds_minus,sizeof(List_T));

      middle_diagonals_minus = (Univdiag_T *) MALLOC(nseeds_minus*sizeof(Univdiag_T));
      best_right_diagonals_minus = (List_T *) MALLOC(nseeds_minus*sizeof(List_T));
      best_left_diagonals_minus = (List_T *) MALLOC(nseeds_minus*sizeof(List_T));
      all_right_diagonals_minus = (List_T *) MALLOC(nseeds_minus*sizeof(List_T));
      all_left_diagonals_minus = (List_T *) MALLOC(nseeds_minus*sizeof(List_T));
      fillin_diagonals_minus = (List_T *) CALLOC(nseeds_minus,sizeof(List_T));

      chrnum = 1;
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,/*chrnum*/1,circular_typeint);
      for (i = 0; i < nseeds_minus; i++) {
	left = best_minus_elt->positions[i];
	if (left > chrhigh) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  /* *chrhigh += 1U; */
	}
	/* May not want to solve for best_right_diagonals and best_left_diagonals.  Use oligoindex instead. */
	scores_minus[i] = get_diagonals(&(middle_diagonals_minus[i]),
					&(best_right_diagonals_minus[i]),&(best_left_diagonals_minus[i]),
					&(all_right_diagonals_minus[i]),&(all_left_diagonals_minus[i]),
					minus_sarray,/*queryptr*/queryrc,querylength,query_compress_rev,
					chroffset,chrhigh,/*goal*/left,minus_elt_array,
					best_minus_i,nminus,/*plusp*/false,genestrand,minus_conversion);
	debug(printf("Got minus score %d\n",scores_minus[i]));
      }

      FREE(minus_elt_array);
    }
  }

#if 0
  /* Because we don't always left-extend, we cannot trust best_score */
  best_score = 0;
  for (i = 0; i < nseeds_plus; i++) {
    if (scores_plus[i] > best_score) {
      best_score = scores_plus[i];
    }
  }
  for (i = 0; i < nseeds_minus; i++) {
    if (scores_minus[i] > best_score) {
      best_score = scores_minus[i];
    }
  }
#endif

  debug(printf("Have %d nseeds_plus and %d nseeds_minus\n",nseeds_plus,nseeds_minus));

#if 0
  coveredp = (bool *) CALLOCA(querylength,sizeof(bool));
  mappings = (Chrpos_T **) MALLOCA(querylength * sizeof(Chrpos_T *));
  npositions = (int *) CALLOCA(querylength,sizeof(int));
  oligoindex = Oligoindex_array_elt(oligoindices_minor,/*source*/0);
  indexsize = Oligoindex_indexsize(oligoindex);
#endif

  /* *sarray_gmap = (List_T) NULL; */

  chrnum = 1;
  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,/*chrnum*/1,circular_typeint);
  for (i = 0; i < nseeds_plus; i++) {
    if (1 /*|| scores_plus[i] > best_score - 20*/) {
      diagonal = middle_diagonals_plus[i];
      left = diagonal->univdiagonal;
      if (left > chrhigh) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* *chrhigh += 1U; */
      }
      middle_path_plus[i] = find_best_path(&(right_paths_plus[i]),&right_endpoints_sense,&right_endpoints_antisense,
					   &right_queryends_sense,&right_queryends_antisense,
					   &right_ambcoords_sense,&right_ambcoords_antisense,
					   &right_amb_knowni_sense,&right_amb_knowni_antisense,
					   &right_amb_nmismatchesi_sense,&right_amb_nmismatchesi_antisense,
					   &right_amb_nmismatchesj_sense,&right_amb_nmismatchesj_antisense,
					   &right_amb_probsi_sense,&right_amb_probsi_antisense,
					   &right_amb_probsj_sense,&right_amb_probsj_antisense,
					   &(left_paths_plus[i]),&left_endpoints_sense,&left_endpoints_antisense,
					   &left_querystarts_sense,&left_querystarts_antisense,
					   &left_ambcoords_sense,&left_ambcoords_antisense,
					   &left_amb_knowni_sense,&left_amb_knowni_antisense,
					   &left_amb_nmismatchesi_sense,&left_amb_nmismatchesi_antisense,
					   &left_amb_nmismatchesj_sense,&left_amb_nmismatchesj_antisense,
					   &left_amb_probsi_sense,&left_amb_probsi_antisense,
					   &left_amb_probsj_sense,&left_amb_probsj_antisense,
					   &(fillin_diagonals_plus[i]),diagonal,best_right_diagonals_plus[i],best_left_diagonals_plus[i],
					   querylength,query_compress_fwd,chroffset,
					   /*plusp*/true,genestrand,/*nmismatches_allowed*/nmisses_allowed);

      hits = solve_via_segments(&(*found_score),&completep,hits,middle_path_plus[i],
				right_endpoints_sense,right_endpoints_antisense,
				right_queryends_sense,right_queryends_antisense,
				right_ambcoords_sense,right_ambcoords_antisense,
				right_amb_knowni_sense,right_amb_knowni_antisense,
				right_amb_nmismatchesi_sense,right_amb_nmismatchesi_antisense,
				right_amb_nmismatchesj_sense,right_amb_nmismatchesj_antisense,
				right_amb_probsi_sense,right_amb_probsi_antisense,
				right_amb_probsj_sense,right_amb_probsj_antisense,

				left_endpoints_sense,left_endpoints_antisense,
				left_querystarts_sense,left_querystarts_antisense,
				left_ambcoords_sense,left_ambcoords_antisense,
				left_amb_knowni_sense,left_amb_knowni_antisense,
				left_amb_nmismatchesi_sense,left_amb_nmismatchesi_antisense,
				left_amb_nmismatchesj_sense,left_amb_nmismatchesj_antisense,
				left_amb_probsi_sense,left_amb_probsi_antisense,
				left_amb_probsj_sense,left_amb_probsj_antisense,

				chrnum,chroffset,chrhigh,chrlength,
				querylength,query_compress_fwd,/*plusp*/true,genestrand,
				/*max_mismatches_allowed*/nmisses_allowed);

#if 0
      if (0 && completep == false) {
	*sarray_gmap = run_gmap_plus(*sarray_gmap,middle_path_plus[i],/*start_paths*/left_paths_plus[i],/*end_paths*/right_paths_plus[i],
				     chrnum,chroffset,chrhigh,chrlength,queryuc_ptr,querylength,
				     genestrand,first_read_p,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				     oligoindices_minor,diagpool,cellpool);
      }
#endif

      Intlist_free(&right_endpoints_sense); Intlist_free(&right_endpoints_antisense);
      Intlist_free(&right_queryends_sense); Intlist_free(&right_queryends_antisense);
      Uintlist_free(&right_ambcoords_sense); Uintlist_free(&right_ambcoords_antisense);
      Intlist_free(&right_amb_knowni_sense); Intlist_free(&right_amb_knowni_antisense);
      Intlist_free(&right_amb_nmismatchesi_sense); Intlist_free(&right_amb_nmismatchesi_antisense);
      Intlist_free(&right_amb_nmismatchesj_sense); Intlist_free(&right_amb_nmismatchesj_antisense);
      Doublelist_free(&right_amb_probsi_sense); Doublelist_free(&right_amb_probsi_antisense);
      Doublelist_free(&right_amb_probsj_sense); Doublelist_free(&right_amb_probsj_antisense);

      Intlist_free(&left_endpoints_sense); Intlist_free(&left_endpoints_antisense);
      Intlist_free(&left_querystarts_sense); Intlist_free(&left_querystarts_antisense);
      Uintlist_free(&left_ambcoords_sense); Uintlist_free(&left_ambcoords_antisense);
      Intlist_free(&left_amb_knowni_sense); Intlist_free(&left_amb_knowni_antisense);
      Intlist_free(&left_amb_nmismatchesi_sense); Intlist_free(&left_amb_nmismatchesi_antisense);
      Intlist_free(&left_amb_nmismatchesj_sense); Intlist_free(&left_amb_nmismatchesj_antisense);
      Doublelist_free(&left_amb_probsi_sense); Doublelist_free(&left_amb_probsi_antisense);
      Doublelist_free(&left_amb_probsj_sense); Doublelist_free(&left_amb_probsj_antisense);
    }
  }

  chrnum = 1;
  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,/*chrnum*/1,circular_typeint);
  for (i = 0; i < nseeds_minus; i++) {
    if (1 /*|| scores_minus[i] > best_score - 20*/) {
      diagonal = middle_diagonals_minus[i];
      left = diagonal->univdiagonal;
      if (left > chrhigh) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* *chrhigh += 1U; */
      }
      middle_path_minus[i] = find_best_path(&(right_paths_minus[i]),&right_endpoints_sense,&right_endpoints_antisense,
					    &right_queryends_sense,&right_queryends_antisense,
					    &right_ambcoords_sense,&right_ambcoords_antisense,
					    &right_amb_knowni_sense,&right_amb_knowni_antisense,
					    &right_amb_nmismatchesi_sense,&right_amb_nmismatchesi_antisense,
					    &right_amb_nmismatchesj_sense,&right_amb_nmismatchesj_antisense,
					    &right_amb_probsi_sense,&right_amb_probsi_antisense,
					    &right_amb_probsj_sense,&right_amb_probsj_antisense,
					    &(left_paths_minus[i]),&left_endpoints_sense,&left_endpoints_antisense,
					    &left_querystarts_sense,&left_querystarts_antisense,
					    &left_ambcoords_sense,&left_ambcoords_antisense,
					    &left_amb_knowni_sense,&left_amb_knowni_antisense,
					    &left_amb_nmismatchesi_sense,&left_amb_nmismatchesi_antisense,
					    &left_amb_nmismatchesj_sense,&left_amb_nmismatchesj_antisense,
					    &left_amb_probsi_sense,&left_amb_probsi_antisense,
					    &left_amb_probsj_sense,&left_amb_probsj_antisense,
					    &(fillin_diagonals_minus[i]),diagonal,best_right_diagonals_minus[i],best_left_diagonals_minus[i],
					    querylength,query_compress_rev,chroffset,
					    /*plusp*/false,genestrand,/*nmismatches_allowed*/nmisses_allowed);
      
      hits = solve_via_segments(&(*found_score),&completep,hits,middle_path_minus[i],
				right_endpoints_sense,right_endpoints_antisense,
				right_queryends_sense,right_queryends_antisense,
				right_ambcoords_sense,right_ambcoords_antisense,
				right_amb_knowni_sense,right_amb_knowni_antisense,
				right_amb_nmismatchesi_sense,right_amb_nmismatchesi_antisense,
				right_amb_nmismatchesj_sense,right_amb_nmismatchesj_antisense,
				right_amb_probsi_sense,right_amb_probsi_antisense,
				right_amb_probsj_sense,right_amb_probsj_antisense,

				left_endpoints_sense,left_endpoints_antisense,
				left_querystarts_sense,left_querystarts_antisense,
				left_ambcoords_sense,left_ambcoords_antisense,
				left_amb_knowni_sense,left_amb_knowni_antisense,
				left_amb_nmismatchesi_sense,left_amb_nmismatchesi_antisense,
				left_amb_nmismatchesj_sense,left_amb_nmismatchesj_antisense,
				left_amb_probsi_sense,left_amb_probsi_antisense,
				left_amb_probsj_sense,left_amb_probsj_antisense,
				
				chrnum,chroffset,chrhigh,chrlength,
				querylength,query_compress_rev,/*plusp*/false,genestrand,
				/*max_mismatches_allowed*/nmisses_allowed);

#if 0
      if (0 && completep == false) {
	*sarray_gmap = run_gmap_minus(*sarray_gmap,middle_path_minus[i],/*start_paths*/right_paths_minus[i],/*end_paths*/left_paths_minus[i],
				      chrnum,chroffset,chrhigh,chrlength,queryuc_ptr,querylength,
				      genestrand,first_read_p,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				      oligoindices_minor,diagpool,cellpool);
      }
#endif

      Intlist_free(&right_endpoints_sense); Intlist_free(&right_endpoints_antisense);
      Intlist_free(&right_queryends_sense); Intlist_free(&right_queryends_antisense);
      Uintlist_free(&right_ambcoords_sense); Uintlist_free(&right_ambcoords_antisense);
      Intlist_free(&right_amb_knowni_sense); Intlist_free(&right_amb_knowni_antisense);
      Intlist_free(&right_amb_nmismatchesi_sense); Intlist_free(&right_amb_nmismatchesi_antisense);
      Intlist_free(&right_amb_nmismatchesj_sense); Intlist_free(&right_amb_nmismatchesj_antisense);
      Doublelist_free(&right_amb_probsi_sense); Doublelist_free(&right_amb_probsi_antisense);
      Doublelist_free(&right_amb_probsj_sense); Doublelist_free(&right_amb_probsj_antisense);

      Intlist_free(&left_endpoints_sense); Intlist_free(&left_endpoints_antisense);
      Intlist_free(&left_querystarts_sense); Intlist_free(&left_querystarts_antisense);
      Uintlist_free(&left_ambcoords_sense); Uintlist_free(&left_ambcoords_antisense);
      Intlist_free(&left_amb_knowni_sense); Intlist_free(&left_amb_knowni_antisense);
      Intlist_free(&left_amb_nmismatchesi_sense); Intlist_free(&left_amb_nmismatchesi_antisense);
      Intlist_free(&left_amb_nmismatchesj_sense); Intlist_free(&left_amb_nmismatchesj_antisense);
      Doublelist_free(&left_amb_probsi_sense); Doublelist_free(&left_amb_probsi_antisense);
      Doublelist_free(&left_amb_probsj_sense); Doublelist_free(&left_amb_probsj_antisense);

    }
  }


#if 0
  /* Salvage using gmap */
  chrnum = 1;
  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,/*chrnum*/1,circular_typeint);
  for (i = 0; i < nseeds_plus; i++) {
    if (incomplete_result_p(middle_path_plus[i],querylength) == true) {
      left = best_plus_elt->positions[i];
      if (left > chrhigh) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* chrhigh += 1U; */
      }
      *sarray_gmap = run_gmap_plus(*sarray_gmap,middle_path_plus[i],/*start_paths*/left_paths_plus[i],/*end_paths*/right_paths_plus[i],
				   chrnum,chroffset,chrhigh,chrlength,queryuc_ptr,querylength,
				   genestrand,first_read_p,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				   oligoindices_minor,diagpool,cellpool);
    }
  }

  chrnum = 1;
  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,/*chrnum*/1,circular_typeint);
  for (i = 0; i < nseeds_minus; i++) {
    if (incomplete_result_p(middle_path_minus[i],querylength) == true) {
      left = best_minus_elt->positions[i];
      if (left > chrhigh) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* chrhigh += 1U; */
      }
      *sarray_gmap = run_gmap_minus(*sarray_gmap,middle_path_minus[i],/*start_paths*/right_paths_minus[i],/*end_paths*/left_paths_minus[i],
				    chrnum,chroffset,chrhigh,chrlength,queryuc_ptr,querylength,
				    genestrand,first_read_p,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				    oligoindices_minor,diagpool,cellpool);
    }
  }
#endif


  if (nseeds_minus > 0) {
    FREE(scores_minus);
    for (i = 0; i < nseeds_minus; i++) {
      for (p = right_paths_minus[i]; p != NULL; p = List_next(p)) {
	diagonal_path = (List_T) List_head(p);
	List_free(&diagonal_path);
      }
      for (p = left_paths_minus[i]; p != NULL; p = List_next(p)) {
	diagonal_path = (List_T) List_head(p);
	List_free(&diagonal_path);
      }
      List_free(&(middle_path_minus[i]));
      List_free(&(left_paths_minus[i]));
      List_free(&(right_paths_minus[i]));
	     

      Univdiag_free(&(middle_diagonals_minus[i]));
      List_free(&(best_right_diagonals_minus[i]));
      List_free(&(best_left_diagonals_minus[i]));
      Univdiag_gc(&(all_right_diagonals_minus[i]));
      Univdiag_gc(&(all_left_diagonals_minus[i]));
      Univdiag_gc(&(fillin_diagonals_minus[i]));
    }
    FREE(middle_diagonals_minus);
    FREE(best_right_diagonals_minus);
    FREE(best_left_diagonals_minus);
    FREE(all_right_diagonals_minus);
    FREE(all_left_diagonals_minus);
    FREE(fillin_diagonals_minus);

    FREE(middle_path_minus);
    FREE(right_paths_minus);
    FREE(left_paths_minus);
  }

  if (nseeds_plus > 0) {
    FREE(scores_plus);
    for (i = 0; i < nseeds_plus; i++) {
      for (p = right_paths_plus[i]; p != NULL; p = List_next(p)) {
	diagonal_path = (List_T) List_head(p);
	List_free(&diagonal_path);
      }
      for (p = left_paths_plus[i]; p != NULL; p = List_next(p)) {
	diagonal_path = (List_T) List_head(p);
	List_free(&diagonal_path);
      }
      List_free(&(middle_path_plus[i]));
      List_free(&(left_paths_plus[i]));
      List_free(&(right_paths_plus[i]));

      Univdiag_free(&(middle_diagonals_plus[i]));
      List_free(&(best_right_diagonals_plus[i]));
      List_free(&(best_left_diagonals_plus[i]));
      Univdiag_gc(&(all_right_diagonals_plus[i]));
      Univdiag_gc(&(all_left_diagonals_plus[i]));
      Univdiag_gc(&(fillin_diagonals_plus[i]));
    }
    FREE(middle_diagonals_plus);
    FREE(best_right_diagonals_plus);
    FREE(best_left_diagonals_plus);
    FREE(all_right_diagonals_plus);
    FREE(all_left_diagonals_plus);
    FREE(fillin_diagonals_plus);

    FREE(middle_path_plus);
    FREE(right_paths_plus);
    FREE(left_paths_plus);
  }

  List_free(&leftward_set_minus);
  List_free(&rightward_set_minus);
  List_free(&leftward_set_plus);
  List_free(&rightward_set_plus);

  for (p = plus_set; p != NULL; p = p->rest) {
    elt = (Elt_T) p->first;
    Elt_free(&elt);
  }
  List_free(&plus_set);

  for (p = minus_set; p != NULL; p = p->rest) {
    elt = (Elt_T) p->first;
    Elt_free(&elt);
  }
  List_free(&minus_set);

  debug(printf("Found %d hits\n",List_length(hits)));

  return hits;
}

