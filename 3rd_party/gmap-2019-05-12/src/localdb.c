static char rcsid[] = "$Id: localdb.c 218530 2019-03-04 23:52:01Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "localdb.h"
#include "localdbdef.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */

#include <sys/mman.h>		/* For munmap */

#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For lseek and close */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#include "assert.h"
#include "mem.h"
#include "fopen.h"
#include "types.h"		/* For Oligospace_T */
#include "filesuffix.h"

#include "epu16-bitpack64-read.h"
#include "epu16-bitpack64-readtwo.h"
#include "genome128_hr.h"

#ifndef LARGE_GENOMES
#include "merge-diagonals-simd-uint4.h"
#elif !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
#include "merge-diagonals-heap.h" /* For Merge_diagonals_uint8 */
#else
#include "merge-diagonals-simd-uint8.h" /* For Merge_diagonals_uint8 */
#endif

#ifdef GSNAP
#include "univdiag.h"
#include "univdiagdef.h"
#endif


#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#if defined(HAVE_SSE2)
#include <emmintrin.h>
#endif
#if defined(HAVE_SSE4_1)
#include <smmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#ifdef HAVE_AVX512
#include <immintrin.h>
#endif



#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* fill procedures */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif


#define T Localdb_T


/************************************************************************
 *   Debugging procedures
 ************************************************************************/

#if 0
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
#endif



#ifndef PMAP

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

#if (defined(DEBUG0) || defined(DEBUG1) || defined(DEBUG2))
static char *
shortoligo_nt (Oligospace_T oligo, Width_T oligosize) {
  char *nt;
  Width_T i, j;
  Oligospace_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return nt;
}
#endif

#endif



#define LOCTABLE_SIZE 2
#define DIFFERENTIAL_METAINFO_SIZE 2


#if 0
static int
Localdb_count (T this, Oligospace_T oligo, int regioni) {
  UINT4 region_strm_start;
  UINT2 *offsetsmeta, *offsetsstrm;
  UINT4 ptr0, end0;

  debug0(printf("Localdb_count: oligo = %06X, regioni %d\n",oligo,regioni));

  offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);
  debug0(printf("Offsetsmeta for regionsize %d showing %u for block start and %u for prefix sum\n",
		this->regionsize,offsetsmeta[0],offsetsmeta[1]));

  region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
  offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

#ifdef DEBUG0
  printf("Loctable showing strm starts at %u registers => %u shorts\n",
	 region_strm_start,region_strm_start * 8);
  printf("Pointing bitpack to %04X %04X %04X %04X %04X %04X %04X %04X\n",
	 offsetsstrm[0],offsetsstrm[1],offsetsstrm[2],offsetsstrm[3],
	 offsetsstrm[4],offsetsstrm[5],offsetsstrm[5],offsetsstrm[7]);
#endif


  ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
  debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

  return (int) (end0 - ptr0);
}
#endif


#ifndef UTILITYP
/* Have to use UINT4, and not UINT2, because highi could be 65536 */
static int
binary_search_with_term (UINT4 lowi, UINT4 highi, UINT2 *positions, UINT4 term, UINT4 goal) {
  UINT4 middlei;

  debug10(printf("entered binary search_with_term with lowi=%d, highi=%d, term=%u, goal=%u\n",
		 lowi,highi,term,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi] + term,middlei,positions[middlei] + term,
		   highi-1,positions[highi-1] + term,goal));
    if (goal < positions[middlei] + term) {
      highi = middlei;
    } else if (goal > positions[middlei] + term) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}
#endif


#if 0
static int
Localdb_count_with_bounds (T this, Oligospace_T oligo, Univcoord_T low, Univcoord_T high) {
  int total_n, n;
  int low_regioni, high_regioni, regioni;
  UINT4 region_strm_start, region_positions_start;
  Chrpos_T region_term;
  UINT2 *offsetsmeta, *offsetsstrm, *region_positions;
  UINT4 ptr0, end0, ptr1, end1;

  debug0(printf("Localdb_read: oligo = %06X\n",oligo));

  low_regioni = low/65536;
  high_regioni = high/65536;

  if (high_regioni == low_regioni) {
    regioni = low_regioni;

    offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);

    region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
    offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

    ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
    debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

    if (end0 == ptr0) {
      return 0;

    } else {
      region_positions_start = this->loctable[regioni * LOCTABLE_SIZE + 1];
      region_positions = &(this->locpositions[region_positions_start]);
      region_term = regioni * 65536;

      ptr1 = binary_search_with_term(ptr0,end0,region_positions,region_term,low);
      end1 = ptr1;
      while (end1 < end0 && region_positions[end1] + region_term <= high) {
	end1++;
      }

      return (int) (end1 - ptr1);
    }

  } else {
    total_n = 0;

    /* low region */
    regioni = low_regioni;
    offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);

    region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
    offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

    ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
    debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

    if ((n = end0 - ptr0) == 0) {
      /* Skip low region */

    } else {
      region_positions_start = this->loctable[regioni * LOCTABLE_SIZE + 1];
      region_positions = &(this->locpositions[region_positions_start]);
      region_term = regioni * 65536;

      ptr1 = binary_search_with_term(ptr0,end0,region_positions,region_term,low);
      end1 = end0;		/* Because we are below high_regioni */

      total_n += (int) (end1 - ptr1);
    }

    /* middle regions */
    for (regioni = low_regioni + 1; regioni < high_regioni; regioni++) {
      offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);

      region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
      offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

      ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
      debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

      /* Count whole region, because we are above low_regioni and high_regioni */
      total_n += (int) (end0 - ptr0);
    }

    /* high region */
    regioni = high_regioni;
    offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);

    region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
    offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

    ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
    debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

    if (end0 == ptr0) {
      /* Skip high region */

    } else {
      region_positions_start = this->loctable[regioni * LOCTABLE_SIZE + 1];
      region_positions = &(this->locpositions[region_positions_start]);
      region_term = regioni * 65536;

      ptr1 = ptr0;		/* Because we are above low_regioni */
      end1 = ptr1;
      while (end1 < end0 && region_term + region_positions[end1] <= high) {
	end1++;
      }

      total_n += (int) (end1 - ptr1);
    }

    return total_n;
  }
}
#endif


#if 0
/* Called by Localdb_test, for debugging only */
static UINT4 *
Localdb_read_uint4 (int *nentries, T this, Oligospace_T oligo, int regioni) {
  UINT4 region_strm_start, region_positions_start;
  UINT2 *offsetsmeta, *offsetsstrm, *region_positions, *in;
  UINT4 *diagonals, region_term;
  UINT4 ptr0, end0;
#ifdef LARGE_GENOMES
  UINT4 ptr;
#elif defined(HAVE_AVX2)
  __m256i _region_term, _in_epu32, _zeroes;
  UINT4 *out;
#elif defined(HAVE_SSE2)
  __m128i _region_term, _in_epu32, _in_epu16, _zeroes;
  UINT4 *out;
#else
  UINT4 ptr;
#endif
  int i;


  debug0(printf("Localdb_read: oligo = %06X, regioni %d\n",oligo,regioni));

  offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);
  debug0(printf("Offsetsmeta for regionsize %d showing %u for block start and %u for prefix sum\n",
		this->regionsize,offsetsmeta[0],offsetsmeta[1]));

  region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
  offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

#ifdef DEBUG0
  printf("Loctable showing strm starts at %u registers => %u shorts\n",
	 region_strm_start,region_strm_start * 8);
  printf("Pointing bitpack to %04X %04X %04X %04X %04X %04X %04X %04X\n",
	 offsetsstrm[0],offsetsstrm[1],offsetsstrm[2],offsetsstrm[3],
	 offsetsstrm[4],offsetsstrm[5],offsetsstrm[5],offsetsstrm[7]);
#endif


  ptr0 = Epu16_bitpack64_read_one(oligo,offsetsmeta,offsetsstrm);
  end0 = Epu16_bitpack64_read_one(oligo+1,offsetsmeta,offsetsstrm);
  debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (UINT4 *) NULL;
  } else {
    region_positions_start = this->loctable[regioni * LOCTABLE_SIZE + 1];
    region_positions = &(this->locpositions[region_positions_start]);
    region_term = regioni * 65536;

    diagonals = (UINT4 *) MALLOC((*nentries)*sizeof(UINT4));

#ifdef LARGE_GENOMES
    for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
      diagonals[i++] = (UINT4) (region_positions[ptr]) + region_term;
    }

#elif defined(HAVE_AVX2)
    _region_term = _mm256_set1_epi32(region_term);
    _zeroes = _mm256_setzero_si256();

    in = &(region_positions[ptr0]);
    out = &(diagonals[0]);
    while (in + 8 < &(region_positions[end0])) {
      _in_epu32 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *) in)); /* Read and convert 8 shorts to ints */
      _mm256_storeu_si256((__m256i *) out, _mm256_add_epi32(_in_epu32, _region_term));
      in += 8;
      out += 8;
    }
    while (in < &(region_positions[end0])) {
      *out++ = (UINT4) (*in++) + region_term;
    }

#elif defined(HAVE_SSE2)
    _region_term = _mm_set1_epi32(region_term);
    _zeroes = _mm_setzero_si128();

    in = &(region_positions[ptr0]);
    out = &(diagonals[0]);
    while (in + 4 < &(region_positions[end0])) {
      _in_epu16 = _mm_loadu_si128((__m128i *) in); /* Read 8 shorts, but we unpack only the first 4 */
      _in_epu32 = _mm_unpacklo_epi16(_in_epu16,_zeroes);
      _mm_storeu_si128((__m128i *) out, _mm_add_epi32(_in_epu32, _region_term));
      in += 4;
      out += 4;
    }
    while (in < &(region_positions[end0])) {
      *out++ = (UINT4) (*in++) + region_term;
    }

#else
    for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
      diagonals[i++] = (UINT4) (region_positions[ptr]) + region_term;
    }
#endif

    debug0(
	   printf("%d entries:",*nentries);
	   for (i = 0; i < *nentries; i++) {
	     printf(" %u",diagonals[i]);
	   }
	   printf("\n");
	   );

    return diagonals;
  }

}
#endif

#if 0
void
Localdb_test (T this, Oligospace_T oligo, int regioni) {
  static UINT4 *positions;
  int nentries, k;

  positions = Localdb_read_uint4(&nentries,this,oligo,regioni);
  printf("%d entries:\n",nentries);
  for (k = 0; k < nentries; k++) {
    printf("%u\n",positions[k]);
  }
  printf("\n");

  return;
}
#endif



#ifndef UTILITYP
/* Output for GSNAP: Returns (UINT4 *) aligned memory because of call to
   Merge_diagonals_uint4.  Need to free with FREE_ALIGN */
/* Output for GSNAPL: Returns (UINT8 *) non-aligned memory.  Free with FREE */
static Univcoord_T *
read_with_bounds (int *nentries, T this, Localspace_T oligo, Univcoord_T low, Univcoord_T high,
		  int diagterm, Univcoord_T **region_stream_alloc, int *region_streamsize_alloc,
		  Univcoord_T *region_alloc) {
  int low_regioni, high_regioni, regioni;
  UINT4 region_strm_start, region_positions_start;
  UINT2 *offsetsmeta, *offsetsstrm, *region_positions, *in;
  Univcoord_T *diagonals, *out, region_term;
  UINT4 ptr0, end0, ptr1, end1;
  int total_n, n;
#ifdef LARGE_GENOMES
  UINT4 ptr;
  int i, k;

#elif defined(HAVE_AVX2)
  __m256i _region_term, _in_epu32, _zeroes;
#ifdef DEBUG0
  int i;
#endif

#elif defined(HAVE_SSE2)
  __m128i _region_term, _in_epu32, _in_epu16, _zeroes;
#ifdef DEBUG0
  int i;
#endif

#else
  UINT4 ptr;
  int i, k;
#endif
  Univcoord_T **stream_array, *stream_ptr, *stream;
  int *streamsize_array;
  int streami = 0;


  debug0(printf("Localdb_read: oligo = %06X, diagterm %d\n",oligo,diagterm));

  low_regioni = low/65536;
  high_regioni = high/65536;

  if (high_regioni == low_regioni) {
    regioni = low_regioni;

    offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);

    region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
    offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

    ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
    debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

    if ((*nentries = end0 - ptr0) == 0) {
      return (Univcoord_T *) NULL;

    } else {
      region_positions_start = this->loctable[regioni * LOCTABLE_SIZE + 1];
      region_positions = &(this->locpositions[region_positions_start]);
      region_term = regioni * 65536;

      ptr1 = binary_search_with_term(ptr0,end0,region_positions,region_term,low);
      end1 = ptr1;
      while (end1 < end0 && region_positions[end1] + region_term <= high) {
	end1++;
      }

      if ((*nentries = end1 - ptr1) == 0) {
	return (Univcoord_T *) NULL;
      } else {
	assert(diagterm <= 0);
	if ((Univcoord_T) -diagterm > region_term) {
	  /* Check for possible negative coordinates */
	  while (ptr1 < end1 && region_positions[ptr1] + region_term < (Univcoord_T) -diagterm) {
	    ptr1++;
	  }
	  if ((*nentries = end1 - ptr1) == 0) {
	    return (Univcoord_T *) NULL;
	  }
	}

#ifndef LARGE_GENOMES
	diagonals = (Univcoord_T *) MALLOC_ALIGN((*nentries)*sizeof(Univcoord_T));
	debug_align(printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
			   diagonals,diagonals,__FILE__,__LINE__));
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
	diagonals = (Univcoord_T *) MALLOC_ALIGN((*nentries)*sizeof(Univcoord_T));
#else
	diagonals = (Univcoord_T *) MALLOC((*nentries)*sizeof(Univcoord_T));
#endif


#ifdef LARGE_GENOMES
	for (ptr = ptr1, i = 0; ptr < end1; ptr++) {
	  diagonals[i++] = (Univcoord_T) (region_positions[ptr]) + region_term + diagterm;
	}

#elif defined(HAVE_AVX2)
	_region_term = _mm256_set1_epi32(region_term + diagterm);
	_zeroes = _mm256_setzero_si256();
	
	in = &(region_positions[ptr1]);
	out = &(diagonals[0]);
	while (in + 8 < &(region_positions[end1])) {
	  _in_epu32 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *) in)); /* Read and convert 8 shorts to ints */
	  _mm256_storeu_si256((__m256i *) out, _mm256_add_epi32(_in_epu32, _region_term));
	  in += 8;
	  out += 8;
	}
	while (in < &(region_positions[end1])) {
	  *out++ = (Univcoord_T) (*in++) + region_term + diagterm;
	}

#elif defined(HAVE_SSE2)
	_region_term = _mm_set1_epi32(region_term + diagterm);
	_zeroes = _mm_setzero_si128();

	in = &(region_positions[ptr1]);
	out = &(diagonals[0]);
	while (in + 4 < &(region_positions[end1])) {
	  _in_epu16 = _mm_loadu_si128((__m128i *) in); /* Read 8 shorts, but we unpack only the first 4 */
	  _in_epu32 = _mm_unpacklo_epi16(_in_epu16,_zeroes);
	  _mm_storeu_si128((__m128i *) out, _mm_add_epi32(_in_epu32, _region_term));
	  in += 4;
	  out += 4;
	}
	while (in < &(region_positions[end1])) {
	  *out++ = (Univcoord_T) (*in++) + region_term + diagterm;
	}

#else
	for (ptr = ptr1, i = 0; ptr < end1; ptr++) {
	  diagonals[i++] = (Univcoord_T) (region_positions[ptr]) + region_term + diagterm;
	}
#endif

	debug0(
	       printf("%d entries:",*nentries);
	       for (i = 0; i < *nentries; i++) {
		 printf(" %u",diagonals[i]);
	       }
	       printf("\n");
	       );

	return diagonals;
      }
    }

  } else {
    /* stream_array = (Univcoord_T **) MALLOC((high_regioni - low_regioni + 1)*sizeof(Univcoord_T *)); */
    /* streamsize_array = (int *) MALLOC((high_regioni - low_regioni + 1 + 1)*sizeof(int)); -- Add extra space for Sedgesort */
    stream_array = region_stream_alloc;
    streamsize_array = region_streamsize_alloc;
    stream_ptr = region_alloc;

    streami = 0;
    total_n = 0;

    /* low region */
    regioni = low_regioni;
    offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);

    region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
    offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

    ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
    debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

    if ((n = end0 - ptr0) == 0) {
      /* Skip stream */

    } else {
      region_positions_start = this->loctable[regioni * LOCTABLE_SIZE + 1];
      region_positions = &(this->locpositions[region_positions_start]);
      region_term = regioni * 65536;

      ptr1 = binary_search_with_term(ptr0,end0,region_positions,region_term,low);
      end1 = end0;		/* Because we are below high_regioni */

      if ((n = end1 - ptr1) == 0) {
	/* Skip stream */
      } else {
	assert(diagterm <= 0);
	if ((Univcoord_T) -diagterm > region_term) {
	  /* Check for possible negative coordinates */
	  while (ptr1 < end1 && region_positions[ptr1] + region_term < (Univcoord_T) -diagterm) {
	    ptr1++;
	  }
	  n = end1 - ptr1;
	}
	if (n == 0) {
	  /* Skip stream */
	} else {
	  debug0(printf("Pushing a stream of size %d\n",n));
	  /* stream = (Univcoord_T *) MALLOC(n*sizeof(Univcoord_T)); */
	  stream = stream_ptr;
	  stream_ptr += n;

#ifdef LARGE_GENOMES
	  k = 0;
	  for (ptr = ptr1; ptr < end1; ptr++) {
	    stream[k++] = region_term + diagterm + (Univcoord_T) (region_positions[ptr]);
	  }

#elif defined(HAVE_AVX2)
	  _region_term = _mm256_set1_epi32(region_term + diagterm);
	  _zeroes = _mm256_setzero_si256();
	
	  in = &(region_positions[ptr1]);
	  out = &(stream[0]);
	  while (in + 8 < &(region_positions[end1])) {
	    _in_epu32 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *) in)); /* Read and convert 8 shorts to ints */
	    _mm256_storeu_si256((__m256i *) out, _mm256_add_epi32(_in_epu32, _region_term));
	    in += 8;
	    out += 8;
	  }
	  while (in < &(region_positions[end1])) {
	    *out++ = (Univcoord_T) (*in++) + region_term + diagterm;
	  }

#elif defined(HAVE_SSE2)
	  _region_term = _mm_set1_epi32(region_term + diagterm);
	  _zeroes = _mm_setzero_si128();

	  in = &(region_positions[ptr1]);
	  out = &(stream[0]);
	  while (in + 4 < &(region_positions[end1])) {
	    _in_epu16 = _mm_loadu_si128((__m128i *) in); /* Read 8 shorts, but we unpack only the first 4 */
	    _in_epu32 = _mm_unpacklo_epi16(_in_epu16,_zeroes);
	    _mm_storeu_si128((__m128i *) out, _mm_add_epi32(_in_epu32, _region_term));
	    in += 4;
	    out += 4;
	  }
	  while (in < &(region_positions[end1])) {
	    *out++ = (Univcoord_T) (*in++) + region_term + diagterm;
	  }

#else	  
	  k = 0;
	  for (ptr = ptr1; ptr < end1; ptr++) {
	    stream[k++] = region_term + diagterm + (Univcoord_T) (region_positions[ptr]);
	  }
#endif

	  stream_array[streami] = stream;
	  streamsize_array[streami] = n;
	  total_n += n;
	  streami++;
	}
      }
    }

    /* middle regions */
    for (regioni = low_regioni + 1; regioni < high_regioni; regioni++) {
      offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);

      region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
      offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

      ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
      debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

      if ((n = end0 - ptr0) == 0) {
	/* Skip stream */

      } else {
	region_positions_start = this->loctable[regioni * LOCTABLE_SIZE + 1];
	region_positions = &(this->locpositions[region_positions_start]);
	region_term = regioni * 65536;

	ptr1 = ptr0;		/* Because we are above low_regioni */
	end1 = end0;		/* Because we are below high_regioni */

	if ((n = end1 - ptr1) == 0) {
	  /* Skip stream */
	} else {
	  assert(diagterm <= 0);
	  if ((Univcoord_T) -diagterm > region_term) {
	    /* Check for possible negative coordinates */
	    while (ptr1 < end1 && region_positions[ptr1] + region_term < (Univcoord_T) -diagterm) {
	      ptr1++;
	    }
	    n = end1 - ptr1;
	  }
	  if (n == 0) {
	    /* Skip stream */
	  } else {
	    debug0(printf("Pushing a stream of size %d\n",n));
	    /* stream = (Univcoord_T *) MALLOC(n*sizeof(Univcoord_T)); */
	    stream = stream_ptr;
	    stream_ptr += n;

#ifdef LARGE_GENOMES
	  k = 0;
	  for (ptr = ptr1; ptr < end1; ptr++) {
	    stream[k++] = region_term + diagterm + (Univcoord_T) (region_positions[ptr]);
	  }

#elif defined(HAVE_AVX2)
	  _region_term = _mm256_set1_epi32(region_term + diagterm);
	  _zeroes = _mm256_setzero_si256();
	
	  in = &(region_positions[ptr1]);
	  out = &(stream[0]);
	  while (in + 8 < &(region_positions[end1])) {
	    _in_epu32 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *) in)); /* Read and convert 8 shorts to ints */
	    _mm256_storeu_si256((__m256i *) out, _mm256_add_epi32(_in_epu32, _region_term));
	    in += 8;
	    out += 8;
	  }
	  while (in < &(region_positions[end1])) {
	    *out++ = (Univcoord_T) (*in++) + region_term + diagterm;
	  }

#elif defined(HAVE_SSE2)
	  _region_term = _mm_set1_epi32(region_term + diagterm);
	  _zeroes = _mm_setzero_si128();

	  in = &(region_positions[ptr1]);
	  out = &(stream[0]);
	  while (in + 4 < &(region_positions[end1])) {
	    _in_epu16 = _mm_loadu_si128((__m128i *) in); /* Read 8 shorts, but we unpack only the first 4 */
	    _in_epu32 = _mm_unpacklo_epi16(_in_epu16,_zeroes);
	    _mm_storeu_si128((__m128i *) out, _mm_add_epi32(_in_epu32, _region_term));
	    in += 4;
	    out += 4;
	  }
	  while (in < &(region_positions[end1])) {
	    *out++ = (Univcoord_T) (*in++) + region_term + diagterm;
	  }

#else	  
	    k = 0;
	    for (ptr = ptr1; ptr < end1; ptr++) {
	      stream[k++] = region_term + diagterm + (Univcoord_T) (region_positions[ptr]);
	    }
#endif

	    stream_array[streami] = stream;
	    streamsize_array[streami] = n;
	    total_n += n;
	    streami++;
	  }
	}
      }
    }

    /* high region */
    regioni = high_regioni;
    offsetsmeta = &(this->locoffsetsmeta[regioni * this->regionsize * DIFFERENTIAL_METAINFO_SIZE]);

    region_strm_start = this->loctable[regioni * LOCTABLE_SIZE];
    offsetsstrm = &(this->locoffsetsstrm[region_strm_start * 8]); /* 8 shorts per 128-bit register */

    ptr0 = Epu16_bitpack64_read_two(&end0,oligo,offsetsmeta,offsetsstrm);
    debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

    if ((n = end0 - ptr0) == 0) {
      /* Skip stream */

    } else {
      region_positions_start = this->loctable[regioni * LOCTABLE_SIZE + 1];
      region_positions = &(this->locpositions[region_positions_start]);
      region_term = regioni * 65536;

      ptr1 = ptr0;		/* Because we are above low_regioni */
      end1 = ptr1;
      while (end1 < end0 && region_term + region_positions[end1] <= high) {
	end1++;
      }

      if ((n = end1 - ptr1) == 0) {
	/* Skip stream */
      } else {
	assert(diagterm <= 0);
	if ((Univcoord_T) -diagterm > region_term) {
	  /* Check for possible negative coordinates */
	  while (ptr1 < end1 && region_positions[ptr1] + region_term < (Univcoord_T) -diagterm) {
	    ptr1++;
	  }
	  n = end1 - ptr1;
	}
	if (n == 0) {
	  /* Skip stream */
	} else {
	  debug0(printf("Pushing a stream of size %d\n",n));
	  /* stream = (Univcoord_T *) MALLOC(n*sizeof(Univcoord_T)); */
	  stream = stream_ptr;
	  stream_ptr += n;

#ifdef LARGE_GENOMES
	  k = 0;
	  for (ptr = ptr1; ptr < end1; ptr++) {
	    stream[k++] = region_term + diagterm + (Univcoord_T) (region_positions[ptr]);
	  }

#elif defined(HAVE_AVX2)
	  _region_term = _mm256_set1_epi32(region_term + diagterm);
	  _zeroes = _mm256_setzero_si256();
	
	  in = &(region_positions[ptr1]);
	  out = &(stream[0]);
	  while (in + 8 < &(region_positions[end1])) {
	    _in_epu32 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *) in)); /* Read and convert 8 shorts to ints */
	    _mm256_storeu_si256((__m256i *) out, _mm256_add_epi32(_in_epu32, _region_term));
	    in += 8;
	    out += 8;
	  }
	  while (in < &(region_positions[end1])) {
	    *out++ = (Univcoord_T) (*in++) + region_term + diagterm;
	  }

#elif defined(HAVE_SSE2)
	  _region_term = _mm_set1_epi32(region_term + diagterm);
	  _zeroes = _mm_setzero_si128();

	  in = &(region_positions[ptr1]);
	  out = &(stream[0]);
	  while (in + 4 < &(region_positions[end1])) {
	    _in_epu16 = _mm_loadu_si128((__m128i *) in); /* Read 8 shorts, but we unpack only the first 4 */
	    _in_epu32 = _mm_unpacklo_epi16(_in_epu16,_zeroes);
	    _mm_storeu_si128((__m128i *) out, _mm_add_epi32(_in_epu32, _region_term));
	    in += 4;
	    out += 4;
	  }
	  while (in < &(region_positions[end1])) {
	    *out++ = (Univcoord_T) (*in++) + region_term + diagterm;
	  }

#else	  
	  k = 0;
	  for (ptr = ptr1; ptr < end1; ptr++) {
	    stream[k++] = region_term + diagterm + (Univcoord_T) (region_positions[ptr]);
	  }
#endif

	  stream_array[streami] = stream;
	  streamsize_array[streami] = n;
	  total_n += n;
	  streami++;
	}
      }
    }

    /* Merge regions */
    if (streami == 0) {
      /* FREE(streamsize_array); -- Using region_streamsize_alloc from caller */
      /* FREE(stream_array); -- Using region_stream_alloc from caller */
      *nentries = 0;
      return (Univcoord_T *) NULL;

    } else {
#ifdef LARGE_GENOMES
      diagonals = Merge_diagonals_uint8(&(*nentries),stream_array,streamsize_array,/*nstreams*/streami);
#else
      diagonals = Merge_diagonals_uint4(&(*nentries),stream_array,streamsize_array,/*nstreams*/streami);
#endif

#if 0
      /* No longer necessary with region_alloc, allocated by caller */
      while (--streami >= 0) {
	stream = stream_array[streami];
	FREE(stream);		/* Not aligned */
      }
#endif
      /* FREE(streamsize_array); -- Using region_streamsize_alloc from caller */
      /* FREE(stream_array); -- Using region_stream_alloc from caller */

#ifdef DEBUG0
      printf("%d entries:",*nentries);
      for (i = 0; i < *nentries; i++) {
	printf(" %u",diagonals[i]);
      }
      printf("\n");
#endif

      return diagonals;
    }
  }
}
#endif


#define VALUE_A 0
#define VALUE_C 1
#define VALUE_G 2
#define VALUE_T 3

static int oligosize;

/* Number of possible genestrands: 3.  Number of possible nucleotides: 4 */
static Localspace_T forward_conv[3][4];
static Localspace_T revcomp_conv[3][4];

#ifdef GSNAP
static int char_to_int[128] =
  { -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 1..10 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 11..20 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 21..30 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 31..40 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 41..50 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 51..60 */
    -1, -1, -1, -1,

  /* A,  B,  C,  D,  E,  F,  G,  H,  I,  J, */
     0, -1,  1, -1, -1, -1,  2, -1, -1, -1, /* 65..74 */

  /* K,  L,  M,  N,  O,  P,  Q,  R,  S,  T  */
    -1, -1, -1, -1, -1, -1, -1, -1, -1,  3, /* 75..84 */

  /* U,  V,  W,  X,  Y,  Z  */
    -1, -1, -1, -1, -1, -1,

    -1, -1, -1, -1, -1, -1,

  /* a,  b,  c,  d,  e,  f,  g,  h,  i,  j, */
     0, -1,  1, -1, -1, -1,  2, -1, -1, -1, /* 97..106 */

  /* k,  l,  m,  n,  o,  p,  q,  r,  s,  t  */
    -1, -1, -1, -1, -1, -1, -1, -1, -1,  3, /* 107..116 */

  /* u,  v,  w,  x,  y,  z  */
    -1, -1, -1, -1, -1, -1,

    -1, -1, -1, -1, -1};
#endif


#ifdef GSNAP
static Localspace_T
nt_oligo_plus (bool *validp, char *queryptr, int genestrand) {
  Localspace_T oligo = 0U;
  int i, d;
  char c;

  *validp = true;
  for (i = 0; i < oligosize; i++) {
    oligo <<= 2;

    c = queryptr[i];
    if ((d = char_to_int[(int) c]) < 0) {
      *validp = false;
      return 0;
    } else {
      oligo += forward_conv[genestrand][d];
    }
  }

  return oligo;
}
#endif


#ifdef GSNAP
static Localspace_T
nt_oligo_minus (bool *validp, char *queryptr, int genestrand) {
  Localspace_T oligo = 0U;
  int i, d;
  char c;

  *validp = true;
  for (i = 0; i < oligosize; i++) {
    oligo <<= 2;

    c = queryptr[i];
    if ((d = char_to_int[(int) c]) < 0) {
      *validp = false;
      return 0;
    } else {
      oligo += revcomp_conv[genestrand][d];
    }
  }

  return oligo;
}
#endif



#if 0

#define MAX_NEIGHBORS 0		/* was 3 */

/* Modified from kmer-search.c */
static Univdiag_T *
find_local_sets (int *nloci, List_T **left_diagonals, List_T **right_diagonals,
		 int subopt_count, Univcoord_T *values, int nvalues,
		 Compress_T query_compress, int pos5, int pos3, Genome_T genomebits,
		 int local1part, bool plusp) {
  Univdiag_T *middle_diagonals, *out;
  Univcoord_T *ptr, *end, *first, nearby_value, left, diagonal;
#if 0
  int trim5, trim3;
#else
  int querystart, queryend;
#endif
  int nmismatches5, nmismatches3;
  int *indices, *outi;
  int nneighbors;
  int i, k;

#ifdef TRIM_CHECK
  int trim5_old, trim3_old;
  int nmismatches5_old, nmismatches3_old;
#endif


  out = middle_diagonals = (Univdiag_T *) MALLOC(nvalues*sizeof(Univdiag_T));
  outi = indices = (int *) MALLOC(nvalues*sizeof(int));

  /* Get all values that have subopt_count or more */
  ptr = values;
  end = &(values[nvalues]);
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u",*first));
    if (ptr + subopt_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[subopt_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[subopt_count - 1] != ptr[subopt_count - 2]) {
	/* Can jump immediately */
	ptr += subopt_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      /* left = *first - querylength; */
      left = *first;
#if 0
      trim_ends(&trim5,&trim3,&nmismatches5,&nmismatches3,poly_p,query_compress,left,querylength,genomebits,plusp,15);
#else
      querystart = Genome_first_kmer_left(&nmismatches5,genomebits,query_compress,left,pos5,pos3,
					  plusp,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/local1part);
      queryend = Genome_first_kmer_right(&nmismatches3,genomebits,query_compress,left,pos5,pos3,
					 plusp,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/local1part);
#endif

      debug0(printf(" => MIDDLE DIAGONAL at left %u has query %d..%d\n",left,querystart,queryend));
      *out++ = Univdiag_new(querystart,queryend,/*univdiagonal*/left);
      *outi++ = ptr - values;	/* Save index so we can find local diagonals later */

      ptr += subopt_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
    }
  }

  *nloci = out - &(middle_diagonals[0]);
  *left_diagonals = (List_T *) CALLOC(*nloci,sizeof(List_T));
  *right_diagonals = (List_T *) CALLOC(*nloci,sizeof(List_T));

  /* Add nearby diagonals */
  /* TODO: Check for too many nearby diagonals, indicating a repetitive oligo */
  for (i = 0; i < *nloci; i++) {
    /* diagonal = middle_diagonals[i]->univdiagonal + querylength; */
    diagonal = middle_diagonals[i]->univdiagonal;

    k = indices[i];
    nneighbors = 0;
    while (nneighbors <= MAX_NEIGHBORS && k - 1 >= 0 && diagonal - values[k - 1] < 200000) {
      nearby_value = values[--k];
      debug0(printf("Potential before %u is %u\n",diagonal,nearby_value));

      /* left = nearby_value - querylength; */
      left = nearby_value;
#if 0
      trim_ends(&trim5,&trim3,&nmismatches5,&nmismatches3,poly_p,query_compress,left,querylength,genomebits,plusp,15);
#else
      querystart = Genome_first_kmer_left(&nmismatches5,genomebits,query_compress,left,pos5,pos3,
					  plusp,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/local1part);
      queryend = Genome_first_kmer_right(&nmismatches3,genomebits,query_compress,left,pos5,pos3,
					 plusp,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/local1part);
#endif

      debug0(printf("LEFT DIAGONAL at left %u has query %d..%d",left,querystart,queryend));
      if (querystart >= middle_diagonals[i]->querystart) {
	debug0(printf(" => skipping because querystart %d >= middle querystart %d",
		      querystart,middle_diagonals[i]->querystart));
      } else {
	(*left_diagonals)[i] =
	  List_push((*left_diagonals)[i],(void *) Univdiag_new(querystart,queryend,/*univdiagonal*/left));
	nneighbors += 1;
      }
      debug0(printf("\n"));
      while (k - 1 >= 0 && values[k - 1] == nearby_value) {
	k--;
      }
    }

    k = indices[i] + subopt_count - 1;
    while (k + 1 < nvalues && values[k + 1] == diagonal) {
      k++;
    }

    nneighbors = 0;
    while (nneighbors <= MAX_NEIGHBORS && k + 1 < nvalues && values[k + 1] - diagonal < 200000) {
      nearby_value = values[++k];
      debug0(printf("Potential after %u is %u\n",diagonal,nearby_value));
      
      /* left = nearby_value - querylength; */
      left = nearby_value;
#if 0
      trim_ends(&trim5,&trim3,&nmismatches5,&nmismatches3,poly_p,query_compress,left,querylength,genomebits,plusp,15);
#else
      querystart = Genome_first_kmer_left(&nmismatches5,genomebits,query_compress,left,pos5,pos3,
					  plusp,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/local1part);
      queryend = Genome_first_kmer_right(&nmismatches3,genomebits,query_compress,left,pos5,pos3,
					 plusp,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/local1part);
#endif

      debug0(printf("RIGHT DIAGONAL at left %u has query %d..%d",left,querystart,queryend));
      if (queryend <= middle_diagonals[i]->queryend) {
	debug0(printf(" => skipping because queryend %d <= middle queryend %d",
		      queryend,middle_diagonals[i]->queryend));
      } else {
	(*right_diagonals)[i] =
	  List_push((*right_diagonals)[i],(void *) Univdiag_new(querystart,queryend,/*univdiagonal*/left));
	nneighbors += 1;
      }
      debug0(printf("\n"));
      while (k + 1 < nvalues && values[k + 1] == nearby_value) {
	k++;
      }
    }
  }

  FREE(indices);
  return middle_diagonals;
}

#endif


#if 0

/* Copied from kmer-search.c */
static int
most_prevalent_uint (int *nloci, Univcoord_T *values, int nvalues) {
  int max_count, count;
  Univcoord_T *out, *ptr, *end, *first;

  assert(nvalues > 0);

  ptr = out = &(values[0]);	/* Reset */
  end = &(values[nvalues]);

  max_count = 1;
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u",*first));
    if (ptr + max_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[max_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[max_count - 1] != ptr[max_count - 2]) {
	/* Can jump immediately */
	ptr += max_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      ptr += max_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
      debug0(printf(" => Count %ld, index %d => printing\n",ptr - first,first - values));
      if ((count = ptr - first) > max_count) {
	out = &(values[0]);	/* Reset */
	max_count = count;
      }
      *out++ = *first;
    }
  }

  *nloci = out - &(values[0]);
  return max_count;
}

#endif


#ifdef GSNAP
/* Copied from kmer-search.c */
static int
most_prevalent_count (Univcoord_T *values, int nvalues) {
  int max_count, count;
  Univcoord_T *ptr, *end, *first;

  assert(nvalues > 0);

  ptr = values;
  end = &(values[nvalues]);

  max_count = 1;
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u",*first));
    if (ptr + max_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[max_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[max_count - 1] != ptr[max_count - 2]) {
	/* Can jump immediately */
	ptr += max_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      ptr += max_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
      debug0(printf(" => Count %ld\n",ptr - first));
      if ((count = ptr - first) > max_count) {
	max_count = count;
      }
    }
  }

  return max_count;
}


/* Modified from kmer-search.c */
/* Overwrites values */
static void
find_prevalent (int *nloci, Univcoord_T *values, int nvalues, int subopt_count) {
  Univcoord_T *out, *ptr, *end, *first;
  /* int querystart, queryend; */
  /* int nmismatches5, nmismatches3; */

  out = &(values[0]);		/* overwrite */

  /* Get all values that have subopt_count or more */
  ptr = values;
  end = &(values[nvalues]);
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u",*first));
    if (ptr + subopt_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[subopt_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[subopt_count - 1] != ptr[subopt_count - 2]) {
	/* Can jump immediately */
	ptr += subopt_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      *out++ = *first;
      
      ptr += subopt_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
    }
  }

  *nloci = out - &(values[0]);
  return;
}


#define SUBOPT 2


/* Output for GSNAP: Returns (UINT4 *) aligned memory because of call
   to Merge_diagonals_uint4.  Internal streams are also aligned from
   call to read_with_bounds */
/* Output for GSNAP: Returns (UINT8 *) non-aligned memory.  Internal
   streams are also non-aligned from call to read_with_bounds */

/* queryptr is either queryuc_ptr (plus) or queryrc (minus) */
Univcoord_T *
Localdb_get_diagonals (int *nentries, T this, char *queryptr, int pos5, int pos3,
		       Univcoord_T low, Univcoord_T high, bool plusp, int genestrand,
		       Univcoord_T **stream_alloc, int *streamsize_alloc,
		       bool remove_repetitive_p) {
  Univcoord_T *diagonals, *stream, *p, *q, last_diagonal;
  Univcoord_T **region_stream_alloc, *region_alloc;
  int *region_streamsize_alloc;
  int low_regioni, high_regioni, nregions;
  Localspace_T oligo;
  bool validp;
  int querypos;
  int ndiagonals, total_n = 0, n;
  int max_count, subopt_count;
  int streami = 0;
#ifdef DEBUG0
  int i;
#endif
  

  debug0(printf("Entered Localdb_get_diagonals with queryptr %s\n",queryptr));
  debug0(printf("  querypos %d..%d, low %u, and high %u\n",pos5,pos3,low,high));

  if (low >= high) {
    /* Have to check for this because a substring that is at the
       beginning or end of a chromosome may attempt to find a local
       piece that is beyond the chromosome, indicated by low == high */
    debug0(printf("  exiting because low %u >= high %u\n",low,high));
    *nentries = 0;
    return (Univcoord_T *) NULL;

  } else {
    low_regioni = low/65536;
    high_regioni = high/65536;
    nregions = high_regioni - low_regioni + 1;

    region_stream_alloc = (Univcoord_T **) MALLOC(nregions*sizeof(Univcoord_T *));
    region_streamsize_alloc = (int *) MALLOC((nregions + 1)*sizeof(int)); /* Add extra space for Sedgesort */

    /* At most, each stream can contain 65536 positions */
    region_alloc = (Univcoord_T *) MALLOC(nregions*65536*sizeof(Univcoord_T));
  }

  if (plusp) {
    for (querypos = pos5; querypos <= pos3 - this->local1part; querypos++) {
      oligo = nt_oligo_plus(&validp,&(queryptr[querypos]),genestrand);
      debug0(printf("Oligo: %s\n",shortoligo_nt(oligo,this->local1part)));

      if (validp == true) {
	stream = read_with_bounds(&n,this,oligo,low,high,/*diagterm*/-querypos,
				  region_stream_alloc,region_streamsize_alloc,region_alloc);
#ifdef DEBUG0
	printf("At querypos %d, %d entries:",querypos,n);
	for (i = 0; i < n; i++) {
	  printf(" %u",stream[i]);
	}
	printf("\n");
#endif
      
	if (n == 0) {
	  /* Skip stream */
	} else {
	  stream_alloc[streami] = stream;
	  streamsize_alloc[streami] = n;
	  total_n += n;
	  streami++;
	}
      }
    }

  } else {
    for (querypos = pos5; querypos <= pos3 - this->local1part; querypos++) {
      oligo = nt_oligo_minus(&validp,&(queryptr[querypos]),genestrand);
      debug0(printf("Oligo: %s\n",shortoligo_nt(oligo,this->local1part)));

      if (validp == true) {
	stream = read_with_bounds(&n,this,oligo,low,high,/*diagterm*/-querypos,
				  region_stream_alloc,region_streamsize_alloc,region_alloc);
#ifdef DEBUG0
	printf("At querypos %d, %d entries:",querypos,n);
	for (i = 0; i < n; i++) {
	  printf(" %u",stream[i]);
	}
	printf("\n");
#endif
      
	if (n == 0) {
	  /* Skip stream */
	} else {
	  stream_alloc[streami] = stream;
	  streamsize_alloc[streami] = n;
	  total_n += n;
	  streami++;
	}
      }
    }
  }

  FREE(region_alloc);
  FREE(region_streamsize_alloc);
  FREE(region_stream_alloc);



  /* Merge streams */
  if (streami == 0) {
    *nentries = 0;
    return (Univcoord_T *) NULL;
  } else {
#ifdef LARGE_GENOMES
    diagonals = Merge_diagonals_uint8(&ndiagonals,stream_alloc,streamsize_alloc,/*nstreams*/streami);
#else
    diagonals = Merge_diagonals_uint4(&ndiagonals,stream_alloc,streamsize_alloc,/*nstreams*/streami);
#endif
    while (--streami >= 0) {
      stream = stream_alloc[streami];
#ifndef LARGE_GENOMES
      FREE_ALIGN(stream);	/* Aligned */
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
      FREE_ALIGN(stream);	/* Aligned */
#else
      FREE(stream);
#endif
    }

#if 0
    /* max_count = */ most_prevalent_uint(&(*nentries),diagonals,ndiagonals);
#else
    max_count = most_prevalent_count(diagonals,ndiagonals);
    if ((subopt_count = max_count - SUBOPT) < 1) {
      subopt_count = 1;
    }
    find_prevalent(&(*nentries),diagonals,ndiagonals,subopt_count);
#endif

    if (remove_repetitive_p == true) {
#ifdef DEBUG0
      printf("\n");
      printf("%d entries before removing repetitive:",*nentries);
      for (i = 0; i < *nentries; i++) {
	printf(" %u",diagonals[i]);
      }
      printf("\n");
#endif

      /* Skip repetitive diagonals */
      q = diagonals;
      last_diagonal = *q;
      for (p = &(diagonals[1]); p < &(diagonals[*nentries]); p++) {
	assert(*p > last_diagonal);
	if (*p == last_diagonal + 1) {
	  /* Skip: repetitive */
	} else {
	  *++q = *p;
	}
	last_diagonal = *p;
      }
      *nentries = (q - diagonals) + 1;
    }
    
#ifdef DEBUG0
    printf("\n");
    printf("%d entries:",*nentries);
    for (i = 0; i < *nentries; i++) {
      printf(" %u",diagonals[i]);
    }
    printf("\n");
#endif

    return diagonals;
  }
}
#endif


#ifdef PMAP
#define BASE_KMER_SAMPLING 3   /* e.g., 677 */
#define KMER_SAMPLING 2   /* e.g., 77 */
#else
#define BASE_KMER_SAMPLING 5   /* e.g., 12153 */
#define KMER_SAMPLING 3   /* e.g., 153 */
#endif


static Localdb_filenames_T
Localdb_filenames_new (char *loctable_filename, char *locpointers_filename,
		       char *locoffsets_filename, char *locpositions_filename,
		       char *loctable_basename_ptr, char *locpointers_basename_ptr,
		       char *locoffsets_basename_ptr, char *locpositions_basename_ptr,
		       char *loctable_local1info_ptr, char *locpointers_local1info_ptr,
		       char *locoffsets_local1info_ptr, char *locpositions_local1info_ptr) {
  Localdb_filenames_T new = (Localdb_filenames_T) MALLOC(sizeof(*new));

  new->loctable_filename = loctable_filename;
  new->locpointers_filename = locpointers_filename;
  new->locoffsets_filename = locoffsets_filename;
  new->locpositions_filename = locpositions_filename;

  new->loctable_basename_ptr = loctable_basename_ptr;
  new->locpointers_basename_ptr = locpointers_basename_ptr;
  new->locoffsets_basename_ptr = locoffsets_basename_ptr;
  new->locpositions_basename_ptr = locpositions_basename_ptr;

  new->loctable_local1info_ptr = loctable_local1info_ptr;
  new->locpointers_local1info_ptr = locpointers_local1info_ptr;
  new->locoffsets_local1info_ptr = locoffsets_local1info_ptr;
  new->locpositions_local1info_ptr = locpositions_local1info_ptr;

  return new;
}

void
Localdb_filenames_free (Localdb_filenames_T *old) {

  FREE((*old)->loctable_filename);
  FREE((*old)->locpointers_filename);
  FREE((*old)->locoffsets_filename);
  FREE((*old)->locpositions_filename);

  FREE(*old);

  return;
}


Localdb_filenames_T
Localdb_get_filenames (
#ifdef PMAP
		       Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
		       Width_T *local1part, Width_T *local1interval,
		       char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
		       Width_T required_local1part, Width_T required_interval, bool offsets_only_p) {
  Blocksize_T blocksize = 64;

  char *loctable_filename, *locoffsetsmeta_filename,
    *locoffsetsstrm_filename, *locpositions_filename,
    *loctable_basename_ptr, *locoffsetsmeta_basename_ptr,
    *locoffsetsstrm_basename_ptr, *locpositions_basename_ptr,
    *loctable_local1info_ptr, *locoffsetsmeta_local1info_ptr,
    *locoffsetsstrm_local1info_ptr, *locpositions_local1info_ptr;

  char *base_filename, *filename;
#ifdef PMAP
  char *pattern1, *pattern2, *a;
  int patternlength1, patternlength2, alphabet_strlen;
  Alphabet_T found_alphabet;
#else
  char *pattern;
  char tens;
#endif
  char interval_char, digit_string[2], *p, *q;
  Width_T found_local1part = 0, found_interval = 0;
  int rootlength, patternlength;

  char ones;
  char *loctable_suffix, *locoffsetsstrm_suffix, *locoffsetsmeta_suffix, *locpositions_suffix;
  struct dirent *entry;
  DIR *dp;


  if (snps_root == NULL) {
    if (blocksize == 32) {
      loctable_suffix = "loctable";
      locoffsetsmeta_suffix = "locoffsets32meta";
      locoffsetsstrm_suffix = "locoffsets32strm";
      locpositions_suffix = "locpositions";
    } else if (blocksize == 64) {
      loctable_suffix = "loctable";
      locoffsetsmeta_suffix = "locoffsets64meta";
      locoffsetsstrm_suffix = "locoffsets64strm";
      locpositions_suffix = "locpositions";
    } else {
      fprintf(stderr,"Unexpected blocksize %d\n",blocksize);
      abort();
    }
  } else {
    if (blocksize == 32) {
      loctable_suffix = (char *) CALLOC(strlen("loctable.")+strlen(snps_root)+1,sizeof(char));
      locoffsetsmeta_suffix = (char *) CALLOC(strlen("locoffsets32meta.")+strlen(snps_root)+1,sizeof(char));
      locoffsetsstrm_suffix = (char *) CALLOC(strlen("locoffsets32strm.")+strlen(snps_root)+1,sizeof(char));
      locpositions_suffix = (char *) CALLOC(strlen("locpositions.")+strlen(snps_root)+1,sizeof(char));

      sprintf(loctable_suffix,"loctable.%s",snps_root);
      sprintf(locoffsetsmeta_suffix,"locoffsets32meta.%s",snps_root);
      sprintf(locoffsetsstrm_suffix,"locoffsets32strm.%s",snps_root);
      sprintf(locpositions_suffix,"locpositions.%s",snps_root);

    } else if (blocksize == 64) {
      loctable_suffix = (char *) CALLOC(strlen("loctable.")+strlen(snps_root)+1,sizeof(char));
      locoffsetsmeta_suffix = (char *) CALLOC(strlen("locoffsets64meta.")+strlen(snps_root)+1,sizeof(char));
      locoffsetsstrm_suffix = (char *) CALLOC(strlen("locoffsets64strm.")+strlen(snps_root)+1,sizeof(char));
      locpositions_suffix = (char *) CALLOC(strlen("locpositions.")+strlen(snps_root)+1,sizeof(char));

      sprintf(loctable_suffix,"loctable.%s",snps_root);
      sprintf(locoffsetsmeta_suffix,"locoffsets64meta.%s",snps_root);
      sprintf(locoffsetsstrm_suffix,"locoffsets64strm.%s",snps_root);
      sprintf(locpositions_suffix,"locpositions.%s",snps_root);

    } else {
      fprintf(stderr,"Unexpected blocksize %d\n",blocksize);
      abort();
    }
  }


#ifdef PMAP
  *alphabet = NALPHABETS + 1;
#endif
  *local1part = 0;
  *local1interval = 1000;
  base_filename = (char *) NULL;

  if ((dp = opendir(genomesubdir)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",genomesubdir);
    exit(9);
  }

#ifdef PMAP
  pattern1 = (char *) CALLOC(strlen(fileroot)+strlen(".")+1,sizeof(char)); /* e.g., "hg19." */
  sprintf(pattern1,"%s.",fileroot);
  patternlength1 = strlen(pattern1);

  pattern2 = (char *) CALLOC(strlen(".")+strlen(idx_filesuffix)+1,sizeof(char)); /* e.g., ".pr" */
  sprintf(pattern2,".%s",idx_filesuffix);
  patternlength2 = strlen(pattern2);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern1,patternlength1)) {
      a = &(filename[strlen(pattern1)]); /* Points after fileroot, e.g., "hg19." */
      if ((p = strstr(a,pattern2)) != NULL && (q = strstr(p,locoffsetsstrm_suffix)) != NULL && !strcmp(q,locoffsetsstrm_suffix)) {
	if ((found_alphabet = Alphabet_find(a)) != AA0) {
	  alphabet_strlen = p - a;
	  p += patternlength2;

	  if (q - p == KMER_SAMPLING) {
	    /* Latest style, e.g., pf77 */
	    if (sscanf(p,"%c%c",&ones,&interval_char) == 2) {
	      digit_string[0] = ones;
	      found_local1part = atoi(digit_string);

	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    } else {
	      abort();
	    }

#if 0
	  } else if (q - p == BASE_KMER_SAMPLING) {
	    /* Previous style, e.g., pf677.  No longer supported. */
	    if (sscanf(p,"%c%c%c",&ones0,&ones,&interval_char) == 3) {
	      digit_string[0] = ones0;
	      found_basesize = atoi(digit_string);

	      digit_string[0] = ones;
	      found_local1part = atoi(digit_string);

	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    } else {
	      abort();
	    }
#endif

	  } else {
	    /* fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename); */
	    if (snps_root != NULL) {
	      FREE(loctable_suffix);
	      FREE(locoffsetsstrm_suffix);
	      FREE(locoffsetsmeta_suffix);
	      FREE(locpositions_suffix);
	    }
	    return (Filenames_T) NULL;
	  }

	  if ((required_alphabet == AA0 || found_alphabet == required_alphabet) &&
	      (required_local1part == 0 || found_local1part == required_local1part) &&
	      (required_interval == 0 || found_interval == required_interval)) {
	    if (required_alphabet == AA0 && found_alphabet > *alphabet) {
	      /* Skip, since we have already found an earlier alphabet */
	    } else if (required_local1part == 0 && found_local1part < *local1part) {
	      /* Skip, since we have already found a larger local1part */
	    } else if (required_interval == 0 && found_interval > *local1interval) {
	      /* Skip, since we have already found a smaller interval */
	    } else {
	      patternlength = patternlength1 + alphabet_strlen + patternlength2;
	      /* *basesize = found_basesize; */
	      *local1part = found_local1part;
	      *local1interval = found_interval;
	      *alphabet = found_alphabet;
	      FREE(base_filename);
	      base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	      strcpy(base_filename,filename);
	    }
	  }
	}
      }
    }
  }

  FREE(pattern2);
  FREE(pattern1);

#else

  pattern = (char *) CALLOC(strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+1,sizeof(char));
  sprintf(pattern,"%s.%s",fileroot,idx_filesuffix);
  patternlength = strlen(pattern);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern,patternlength)) {
      p = &(filename[strlen(pattern)]); /* Points after idx_filesuffix, e.g., "ref" */
      if ((q = strstr(p,locoffsetsstrm_suffix)) != NULL && !strcmp(q,locoffsetsstrm_suffix)) {

	if (q - p == KMER_SAMPLING) {
	  /* New style, e.g., ref153 */
	  if (sscanf(p,"%c%c%c",&tens,&ones,&interval_char) == 3) {
	    digit_string[0] = tens;
	    found_local1part = 10*atoi(digit_string);
	    digit_string[0] = ones;
	    found_local1part += atoi(digit_string);

	    digit_string[0] = interval_char;
	    found_interval = atoi(digit_string);
	  } else {
	    abort();
	  }

	} else {
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s: found %ld characters, expecting %d\n",
		  idx_filesuffix,filename,q-p,BASE_KMER_SAMPLING);
	  abort();
	  if (snps_root != NULL) {
	    FREE(loctable_suffix);
	    FREE(locoffsetsstrm_suffix);
	    FREE(locoffsetsmeta_suffix);
	    FREE(locpositions_suffix);
	  }
	  return (Localdb_filenames_T) NULL;
	}

	if ((required_local1part == 0 || found_local1part == required_local1part) &&
	    (required_interval == 0 || found_interval == required_interval)) {
	  if (required_local1part == 0 && found_local1part < *local1part) {
	    /* Skip, since we have already found a larger local1part */
	  } else if (required_interval == 0 && found_interval > *local1interval) {
	    /* Skip, since we have already found a smaller interval */
	  } else {
	    *local1part = found_local1part;
	    *local1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	}
      }
    }
  }

  FREE(pattern);
#endif


  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }

  /* Construct full filenames */
  if (base_filename == NULL) {
    /* offsetspages_filename = (char *) NULL; */
    /* offsetsmeta_filename = (char *) NULL; */
    /* offsetsstrm_filename = (char *) NULL; */
    /* positions_high_filename = (char *) NULL; */
    /* positions_low_filename = (char *) NULL; */
    if (snps_root != NULL) {
      FREE(loctable_suffix);
      FREE(locoffsetsstrm_suffix);
      FREE(locoffsetsmeta_suffix);
      FREE(locpositions_suffix);
    }
    return (Localdb_filenames_T) NULL;

  } else {
    locoffsetsstrm_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(base_filename)+1,sizeof(char));
    locoffsetsstrm_basename_ptr = &(locoffsetsstrm_filename[strlen(genomesubdir)+strlen("/")]);
    locoffsetsstrm_local1info_ptr = &(locoffsetsstrm_basename_ptr[patternlength]);

    sprintf(locoffsetsstrm_filename,"%s/%s",genomesubdir,base_filename);
    if (Access_file_exists_p(locoffsetsstrm_filename) == false) {
      fprintf(stderr,"Offsets filename %s does not exist\n",locoffsetsstrm_filename);
      FREE(locoffsetsstrm_filename);
      /* offsetsstrm_filename = (char *) NULL; */
      /* positions_high_filename = (char *) NULL; */
      /* positions_low_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(loctable_suffix);
	FREE(locoffsetsstrm_suffix);
	FREE(locoffsetsmeta_suffix);
	FREE(locpositions_suffix);
      }
      return (Localdb_filenames_T) NULL;
    }


    if ((q = strstr(base_filename,locoffsetsstrm_suffix)) == NULL) {
      abort();
    } else {
      rootlength = q - base_filename;
    }

    loctable_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(loctable_suffix)+1,sizeof(char));
    loctable_basename_ptr = &(loctable_filename[strlen(genomesubdir)+strlen("/")]);
    loctable_local1info_ptr = &(loctable_basename_ptr[patternlength]);

    sprintf(loctable_filename,"%s/",genomesubdir);
    strncpy(loctable_basename_ptr,base_filename,rootlength);
    strcpy(&(loctable_basename_ptr[rootlength]),loctable_suffix);


    locoffsetsmeta_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(locoffsetsmeta_suffix)+1,sizeof(char));
    locoffsetsmeta_basename_ptr = &(locoffsetsmeta_filename[strlen(genomesubdir)+strlen("/")]);
    locoffsetsmeta_local1info_ptr = &(locoffsetsmeta_basename_ptr[patternlength]);

    sprintf(locoffsetsmeta_filename,"%s/",genomesubdir);
    strncpy(locoffsetsmeta_basename_ptr,base_filename,rootlength);
    strcpy(&(locoffsetsmeta_basename_ptr[rootlength]),locoffsetsmeta_suffix);

    if (Access_file_exists_p(locoffsetsmeta_filename) == false) {
      fprintf(stderr,"Offsetsmeta filename %s does not exist\n",locoffsetsmeta_filename);
      FREE(locoffsetsstrm_filename);
      /* offsetsstrm_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(loctable_suffix);
	FREE(locoffsetsstrm_suffix);
	FREE(locoffsetsmeta_suffix);
	FREE(locpositions_suffix);
      }
      return (Localdb_filenames_T) NULL;
    }


    locpositions_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(locpositions_suffix)+1,sizeof(char));
    locpositions_basename_ptr = &(locpositions_filename[strlen(genomesubdir)+strlen("/")]);
    locpositions_local1info_ptr = &(locpositions_basename_ptr[patternlength]);

    sprintf(locpositions_filename,"%s/",genomesubdir);
    strncpy(locpositions_basename_ptr,base_filename,rootlength);
    strcpy(&(locpositions_basename_ptr[rootlength]),locpositions_suffix);

    if (offsets_only_p == true) {
      /* Do not look for a positions file */
    } else if (Access_file_exists_p(locpositions_filename) == false) {
      /* Try newer naming scheme: ref153positions instead of ref12153positions */
      sprintf(locpositions_filename,"%s/",genomesubdir);
      strncpy(locpositions_basename_ptr,base_filename,rootlength-BASE_KMER_SAMPLING); /* e.g., skip "12153" */
      strncpy(&(locpositions_basename_ptr[rootlength-BASE_KMER_SAMPLING]),&(base_filename[rootlength-KMER_SAMPLING]),KMER_SAMPLING);
      strcpy(&(locpositions_basename_ptr[rootlength+KMER_SAMPLING-BASE_KMER_SAMPLING]),locpositions_suffix);

      if (Access_file_exists_p(locpositions_filename) == false) {
	fprintf(stderr,"Positions filename %s does not exist\n",locpositions_filename);
	FREE(loctable_filename);
	FREE(locoffsetsmeta_filename);
	FREE(locoffsetsstrm_filename);
	FREE(locpositions_filename);
	FREE(base_filename);
	if (snps_root != NULL) {
	  FREE(loctable_suffix);
	  FREE(locoffsetsstrm_suffix);
	  FREE(locoffsetsmeta_suffix);
	  FREE(locpositions_suffix);
	}
	return (Localdb_filenames_T) NULL;
      }
    }

    if (snps_root != NULL) {
      FREE(loctable_suffix);
      FREE(locoffsetsstrm_suffix);
      FREE(locoffsetsmeta_suffix);
      FREE(locpositions_suffix);
    }

    FREE(base_filename);

    fprintf(stderr,"Looking for local files in directory %s\n",genomesubdir);
    fprintf(stderr,"  Table file is %s\n",loctable_basename_ptr);
    fprintf(stderr,"  Pointers file is %s\n",locoffsetsmeta_basename_ptr);
    fprintf(stderr,"  Offsets file is %s\n",locoffsetsstrm_basename_ptr);
    fprintf(stderr,"  Positions file is %s\n",locpositions_basename_ptr);
    return Localdb_filenames_new(loctable_filename,locoffsetsmeta_filename,
				 locoffsetsstrm_filename,locpositions_filename,
				 loctable_basename_ptr,locoffsetsmeta_basename_ptr,
				 locoffsetsstrm_basename_ptr,locpositions_basename_ptr,
				 loctable_local1info_ptr,locoffsetsmeta_local1info_ptr,
				 locoffsetsstrm_local1info_ptr,locpositions_local1info_ptr);

  }
}



static Oligospace_T
power (int base, Width_T exponent) {
#ifdef OLIGOSPACE_NOT_LONG
  Oligospace_T result = 1U;
#else
  Oligospace_T result = 1UL;
#endif
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


T
Localdb_new_genome (Width_T *local1part, Width_T *local1interval,
		    char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		    Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		    Width_T required_local1part, Width_T required_interval,
		    Access_mode_T locoffsetsstrm_access, Access_mode_T locpositions_access, bool sharedp,
		    bool multiple_sequences_p, bool unload_shared_memory_p) {
  T new = (T) MALLOC(sizeof(*new));

  Localdb_filenames_T filenames;

  /* Oligospace_T basespace, base; */
  /* Oligospace_T poly_T; */
  /* Positionsptr_T ptr0; -- UINT8 or UINT4 */
  /* Positionsptr_T end0; -- UINT8 or UINT4 */

  char *comma;
  double seconds;
#ifdef HAVE_MMAP
  int npages;
#endif


  if ((filenames = Localdb_get_filenames(
#ifdef PMAP
					 &(*alphabet),required_alphabet,
#endif
					 &new->local1part,&new->local1interval,
					 genomesubdir,fileroot,idx_filesuffix,snps_root,
					 required_local1part,required_interval,
					 /*offsets_only_p*/false)) == NULL) {
    /* For backward compatibility and large genomes, allow alignment without localdb */
    /* fprintf(stderr,"Cannot find genomic index files in either current or old format.  Looking for files containing %s\n",idx_filesuffix); */
    FREE(new);
    return (T) NULL;

  } else {

    /* Try bitpack compression  */
    *local1part = new->local1part;
    *local1interval = new->local1interval;

    new->compression_type = BITPACK64_COMPRESSION;

#ifdef PMAP
    *alphabet_size = Alphabet_get_size(*alphabet);
#endif
    new->regionsize = power(4,(*local1part) - 3); /* Subtracting 3 for blocksize of 64 */
    new->blocksize = 64;  /* Used to be determined by 4^(kmer - basesize), but now fixed at 64 */

    /* offsetsmeta and offsetsmeta always ALLOCATED */
    if (snps_root) {
      fprintf(stderr,"Allocating memory for %s (%s) offset pointers, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->local1part,new->local1interval);
    } else {
      fprintf(stderr,"Allocating memory for %s offset pointers, kmer %d, interval %d...",
	      idx_filesuffix,new->local1part,new->local1interval);
    }
#ifdef HAVE_MMAP
    if (multiple_sequences_p == false && unload_shared_memory_p == false) {
      new->locoffsetsmeta = (UINT2 *) Access_mmap(&new->locoffsetsmeta_fd,&new->locoffsetsmeta_len,
						  filenames->locpointers_filename,/*randomp*/false);
      new->locoffsetsmeta_access = MMAPPED;
    } else
#endif
      if (sharedp == true) {
	new->locoffsetsmeta = (UINT2 *) Access_allocate_shared(&new->locoffsetsmeta_access,&new->locoffsetsmeta_shmid,&new->locoffsetsmeta_key,
							       &new->locoffsetsmeta_fd,&new->locoffsetsmeta_len,&seconds,
							       filenames->locpointers_filename,sizeof(UINT2));
      } else {
	new->locoffsetsmeta = (UINT2 *) Access_allocate_private(&new->locoffsetsmeta_access,&new->locoffsetsmeta_len,&seconds,
								filenames->locpointers_filename,sizeof(UINT2));
      }

    comma = Genomicpos_commafmt(new->locoffsetsmeta_len);
    fprintf(stderr,"done (%s bytes",comma);
    FREE(comma);
    if (multiple_sequences_p == true) {
      fprintf(stderr,", %.2f sec",seconds);
    }
    fprintf(stderr,")\n");

    /* offsetsstrm could be ALLOCATED or MMAPPED +/- PRELOAD */
    if (locoffsetsstrm_access == USE_ALLOCATE) {
      if (snps_root) {
	fprintf(stderr,"Allocating memory for %s (%s) offsets, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->local1part,new->local1interval);
      } else {
	fprintf(stderr,"Allocating memory for %s offsets, kmer %d, interval %d...",
		idx_filesuffix,new->local1part,new->local1interval);
      }
      if (sharedp == true) {
	new->locoffsetsstrm = (UINT2 *) Access_allocate_shared(&new->locoffsetsstrm_access,&new->locoffsetsstrm_shmid,&new->locoffsetsstrm_key,
							       &new->locoffsetsstrm_fd,&new->locoffsetsstrm_len,&seconds,
							       filenames->locoffsets_filename,sizeof(UINT2));
      } else {
	new->locoffsetsstrm = (UINT2 *) Access_allocate_private(&new->locoffsetsstrm_access,&new->locoffsetsstrm_len,&seconds,
								filenames->locoffsets_filename,sizeof(UINT2));
      }
      if (new->locoffsetsstrm == NULL) {
	fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	exit(9);
      } else {
	comma = Genomicpos_commafmt(new->locoffsetsstrm_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);
      }

#ifdef HAVE_MMAP
    } else if (locoffsetsstrm_access == USE_MMAP_PRELOAD) {
      if (snps_root) {
	fprintf(stderr,"Pre-loading %s (%s) offsets, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->local1part,new->local1interval);
      } else {
	fprintf(stderr,"Pre-loading %s offsets, kmer %d, interval %d...",
		idx_filesuffix,new->local1part,new->local1interval);
      }
      new->locoffsetsstrm = (UINT2 *) Access_mmap_and_preload(&new->locoffsetsstrm_fd,&new->locoffsetsstrm_len,&npages,&seconds,
							      filenames->locoffsets_filename,sizeof(UINT2));
      if (new->locoffsetsstrm == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead, but program may not run)\n");
#ifdef PMAP
	new->locoffsetsstrm_access = FILEIO;
#else
	exit(9);
#endif
      } else {
	comma = Genomicpos_commafmt(new->locoffsetsstrm_len);
	fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	FREE(comma);
	new->locoffsetsstrm_access = MMAPPED;
      }

    } else if (locoffsetsstrm_access == USE_MMAP_ONLY) {
      new->locoffsetsstrm = (UINT2 *) Access_mmap(&new->locoffsetsstrm_fd,&new->locoffsetsstrm_len,
						  filenames->locoffsets_filename,/*randomp*/false);
      if (new->locoffsetsstrm == NULL) {
	fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		filenames->locoffsets_filename);
#ifdef PMAP
	new->locoffsetsstrm_access = FILEIO;
#else
	exit(9);
#endif
      } else {
	new->locoffsetsstrm_access = MMAPPED;
      }
#endif

    } else if (locoffsetsstrm_access == USE_FILEIO) {
      fprintf(stderr,"Offsetsstrm file I/O access of %s not allowed\n",filenames->locoffsets_filename);
      exit(9);

    } else {
      fprintf(stderr,"Don't recognize offsetsstrm_access type %d\n",locoffsetsstrm_access);
      abort();
    }

#if 0
    printf("Localdb strm is %04X %04X %04X %04X %04X %04X %04X %04X\n",
	   new->locoffsetsstrm[0],new->locoffsetsstrm[1],new->locoffsetsstrm[2],new->locoffsetsstrm[3],
	   new->locoffsetsstrm[4],new->locoffsetsstrm[5],new->locoffsetsstrm[5],new->locoffsetsstrm[7]);
#endif

#ifdef PMAP
#else
#ifdef HAVE_MMAP
    if (multiple_sequences_p == false && unload_shared_memory_p == false) {
      new->loctable = (UINT4 *) Access_mmap(&new->loctable_fd,&new->loctable_len,
					    filenames->loctable_filename,/*randomp*/false);
      new->loctable_access = MMAPPED;
    } else
#endif
      if (sharedp == true) {
	new->loctable = (UINT4 *) Access_allocate_shared(&new->loctable_access,&new->loctable_shmid,&new->loctable_key,
							 &new->loctable_fd,&new->loctable_len,&seconds,
							 filenames->loctable_filename,sizeof(UINT4));
      } else {
	new->loctable = (UINT4 *) Access_allocate_private(&new->loctable_access,&new->loctable_len,&seconds,
							  filenames->loctable_filename,sizeof(UINT4));
      }

#endif	/* PMAP */

    if (new->loctable == NULL) {
      fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->loctable_len);
      fprintf(stderr,"done (%s bytes",comma);
      FREE(comma);
      if (multiple_sequences_p == true) {
	fprintf(stderr,", %.2f sec",seconds);
      }
      fprintf(stderr,")\n");
    }
  }


  /* Read or memory map positions file */
  if (locpositions_access == USE_ALLOCATE) {
    if (snps_root) {
      fprintf(stderr,"Allocating memory for %s (%s) positions, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->local1part,new->local1interval);
    } else {
      fprintf(stderr,"Allocating memory for %s positions, kmer %d, interval %d...",
	      idx_filesuffix,new->local1part,new->local1interval);
    }

#ifdef HAVE_MMAP
    if (multiple_sequences_p == false && unload_shared_memory_p == false) {
      new->locpositions = (UINT2 *) Access_mmap(&new->locpositions_fd,&new->locpositions_len,
						filenames->locpositions_filename,/*randomp*/false);
      new->locpositions_access = MMAPPED;
    } else
#endif

      if (sharedp == true) {
	new->locpositions = (UINT2 *) Access_allocate_shared(&new->locpositions_access,&new->locpositions_shmid,&new->locpositions_key,
							     &new->locpositions_fd,&new->locpositions_len,&seconds,
							     filenames->locpositions_filename,sizeof(UINT2));
      } else {
	new->locpositions = (UINT2 *) Access_allocate_private(&new->locpositions_access,&new->locpositions_len,&seconds,
							      filenames->locpositions_filename,sizeof(UINT2));
      }

    if (new->locpositions == NULL) {
      fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->locpositions_len);
      fprintf(stderr,"done (%s bytes",comma);
      FREE(comma);
      if (multiple_sequences_p == true) {
	fprintf(stderr,", %.2f sec",seconds);
      }
      fprintf(stderr,")\n");
    }
    

#ifdef HAVE_MMAP
  } else if (locpositions_access == USE_MMAP_PRELOAD) {

    if (snps_root) {
      fprintf(stderr,"Pre-loading %s (%s) positions, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->local1part,new->local1interval);
    } else {
      fprintf(stderr,"Pre-loading %s positions, kmer %d, interval %d...",
	      idx_filesuffix,new->local1part,new->local1interval);
    }

    new->locpositions = (UINT2 *) Access_mmap_and_preload(&new->locpositions_fd,&new->locpositions_len,&npages,&seconds,
							  filenames->locpositions_filename,sizeof(UINT2));
    if (new->locpositions == NULL) {
      fprintf(stderr,"insufficient memory (will use disk file instead, but program will be slow)\n");
      new->locpositions_access = FILEIO;
    } else {
      comma = Genomicpos_commafmt(new->locpositions_len);
      fprintf(stderr,"done (%s bytes",comma);
      FREE(comma);
      if (multiple_sequences_p == true) {
	fprintf(stderr,", %.2f sec",seconds);
      }
      fprintf(stderr,")\n");
      new->locpositions_access = MMAPPED;
    }

  } else if (locpositions_access == USE_MMAP_ONLY) {
    new->locpositions = (UINT2 *) Access_mmap(&new->locpositions_fd,&new->locpositions_len,
					      filenames->locpositions_filename,/*randomp*/true);

    if (new->locpositions == NULL) {
      fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program will be slow)\n",
	      filenames->locpositions_filename);
      new->locpositions_access = FILEIO;
    } else {
      new->locpositions_access = MMAPPED;
    }

#endif

  } else {
    fprintf(stderr,"Don't recognize positions_access %d\n",locpositions_access);
    abort();
  }

  Localdb_filenames_free(&filenames);

  return new;
}


void
Localdb_free (T *old) {
  if (*old) {

    if ((*old)->locpositions_access == ALLOCATED_PRIVATE) {
      FREE((*old)->locpositions);

    } else if ((*old)->locpositions_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->locpositions,(*old)->locpositions_shmid,(*old)->locpositions_key);

#ifdef HAVE_MMAP
    } else if ((*old)->locpositions_access == MMAPPED) {
      munmap((void *) (*old)->locpositions,(*old)->locpositions_len);
      close((*old)->locpositions_fd);
#endif
    }

    if ((*old)->locoffsetsstrm_access == ALLOCATED_PRIVATE) {
      FREE((*old)->locoffsetsstrm);

    } else if ((*old)->locoffsetsstrm_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->locoffsetsstrm,(*old)->locoffsetsstrm_shmid,(*old)->locoffsetsstrm_key);

#ifdef HAVE_MMAP
    } else if ((*old)->locoffsetsstrm_access == MMAPPED) {
      munmap((void *) (*old)->locoffsetsstrm,(*old)->locoffsetsstrm_len);
      close((*old)->locoffsetsstrm_fd);
#endif
    }
      
    if ((*old)->locoffsetsmeta_access == ALLOCATED_PRIVATE) {
      FREE((*old)->locoffsetsmeta);
    } else if ((*old)->locoffsetsmeta_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->locoffsetsmeta,(*old)->locoffsetsmeta_shmid,(*old)->locoffsetsmeta_key);
#ifdef HAVE_MMAP
    } else if ((*old)->locoffsetsmeta_access == MMAPPED) {
      munmap((void *) (*old)->locoffsetsmeta,(*old)->locoffsetsmeta_len);
      close((*old)->locoffsetsmeta_fd);
#endif
    }

    if ((*old)->loctable_access == ALLOCATED_PRIVATE) {
      FREE((*old)->loctable);
    } else if ((*old)->loctable_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->loctable,(*old)->loctable_shmid,(*old)->loctable_key);
#ifdef HAVE_MMAP
    } else if ((*old)->loctable_access == MMAPPED) {
      munmap((void *) (*old)->loctable,(*old)->loctable_len);
      close((*old)->loctable_fd);
#endif
    }

    FREE(*old);
  }
  return;
}


/* Modified from Oligo_setup */
void
Localdb_setup (int local1part, int mode) {
  oligosize = local1part;

  if (mode == STANDARD) {
    memcpy(forward_conv[+0],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));
    memcpy(revcomp_conv[+0],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));

  } else if (mode == CMET_STRANDED) {
    memcpy(forward_conv[+0],((Localspace_T[]){/*   */VALUE_A, /*C>T*/VALUE_T, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));
    memcpy(revcomp_conv[+0],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*G>A*/VALUE_A, /*   */VALUE_T}),4*sizeof(Localspace_T));

  } else if (mode == CMET_NONSTRANDED) {
    /* genestrand: 0 or +1 */
    /* Queryuc_ptr from gplus strand: (C==T).  Matches genomicfwd kmer */
    /* Queryrc from gminus strand: (G==A),  Then revcomp to match genomic genomicfwd kmer: (C==T) */

    /* genestrand: +2 (pcr step) */
    /* Queryuc_ptr from gminus strand: (C==T). Then pcr: (G==A).  Matches genomicfwd kmer */
    /* Queryrc from gplus strand: (G==A).  Then pcr: (C==T).  Then revcomp to match genomicfwd kmer: (G==A) */

    memcpy(forward_conv[+1],((Localspace_T[]){/*   */VALUE_A, /*C>T*/VALUE_T, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));
    memcpy(forward_conv[+2],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*G>A*/VALUE_A, /*   */VALUE_T}),4*sizeof(Localspace_T));

    memcpy(revcomp_conv[+1],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*G>A*/VALUE_A, /*   */VALUE_T}),4*sizeof(Localspace_T));
    memcpy(revcomp_conv[+2],((Localspace_T[]){/*   */VALUE_A, /*C>T*/VALUE_T, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));

  } else if (mode == ATOI_STRANDED) {
    memcpy(forward_conv[+0],((Localspace_T[]){/*A>G*/VALUE_G, /*   */VALUE_C, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));
    memcpy(revcomp_conv[+0],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*   */VALUE_G, /*T>C*/VALUE_C}),4*sizeof(Localspace_T));

  } else if (mode == ATOI_NONSTRANDED) {
    memcpy(forward_conv[+1],((Localspace_T[]){/*A>G*/VALUE_G, /*   */VALUE_C, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));
    memcpy(forward_conv[+2],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*   */VALUE_G, /*T>C*/VALUE_C}),4*sizeof(Localspace_T));

    memcpy(revcomp_conv[+1],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*   */VALUE_G, /*T>C*/VALUE_C}),4*sizeof(Localspace_T));
    memcpy(revcomp_conv[+2],((Localspace_T[]){/*A>G*/VALUE_G, /*   */VALUE_C, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));

  } else if (mode == TTOC_STRANDED) {
    memcpy(forward_conv[+0],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*   */VALUE_G, /*T>C*/VALUE_C}),4*sizeof(Localspace_T));
    memcpy(revcomp_conv[+0],((Localspace_T[]){/*A>G*/VALUE_G, /*   */VALUE_C, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));

  } else if (mode == TTOC_NONSTRANDED) {
    memcpy(forward_conv[+1],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*   */VALUE_G, /*T>C*/VALUE_C}),4*sizeof(Localspace_T));
    memcpy(forward_conv[+2],((Localspace_T[]){/*A>G*/VALUE_G, /*   */VALUE_C, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));

    memcpy(revcomp_conv[+1],((Localspace_T[]){/*A>G*/VALUE_G, /*   */VALUE_C, /*   */VALUE_G, /*   */VALUE_T}),4*sizeof(Localspace_T));
    memcpy(revcomp_conv[+2],((Localspace_T[]){/*   */VALUE_A, /*   */VALUE_C, /*   */VALUE_G, /*T>C*/VALUE_C}),4*sizeof(Localspace_T));

  } else {
    fprintf(stderr,"Unexpected mode %d\n",mode);
    exit(9);
  }
}

