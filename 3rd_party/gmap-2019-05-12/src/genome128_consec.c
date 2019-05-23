static char rcsid[] = "$Id: genome128_consec.c 218286 2019-01-23 16:46:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "genome128_consec.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For tolower() */

#include "assert.h"
#include "except.h"
#include "compress.h"
#include "popcount.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


/* Consecutive_matches_rightward and leftward */
/* Slower with shift and wrap, perhaps because we need to extract integers from the SIMD object */
/* #define USE_SHIFT_FIRST_MISMATCH 1 */
/* #define USE_WRAP_FIRST_MISMATCH 1 */

/* Genome_mismatches_right and left */
/* Slower with shift and wrap, probably because we need to loop over the SIMD object */
/* #define USE_SHIFT_MISMATCH_POSITIONS 1 */
/* #define USE_WRAP_MISMATCH_POSITIONS 1 */

/* Genome_count_mismatches_substring */
/* Faster with shift and wrap.  Does not involve any loops. */
#define USE_SHIFT_POPCOUNT 1
#define USE_WRAP_POPCOUNT 1

/* Genome_mismatches_right_trim and left_trim */
/* Slower with shift and wrap */
/* #define USE_SHIFT_TRIM 1 */
/* #define USE_WRAP_TRIM 1 */


/* Faster to use a straight shift, and _mm_bsrli_si128 is not defined in gcc 4.7 */
/* #define USE_SHIFT_HILO 1 */


#ifdef HAVE_SSE2
#define QUERY_NEXTCOL 1		/* high0, high1, high2, high3 */
#define QUERY_NEXTROW 8
#else
#define QUERY_NEXTCOL 3		/* high, low, flags */
/* #define QUERY_NEXTROW 0 */
#endif

#define GENOME_NEXTCOL 1
#define GENOME_NEXTROW 8


#ifdef WORDS_BIGENDIAN
/* Do not use SIMD */
#elif defined(HAVE_SSE2)
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

#if !defined(HAVE_SSE4_2)
/* Skip popcnt, which comes after SSE4.2 */
#elif defined(HAVE_POPCNT)
#include <immintrin.h>
#elif defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip lzcnt and tzcnt, which come after SSE4.2 */
#elif defined(HAVE_LZCNT) || defined(HAVE_TZCNT)
#include <immintrin.h>
#endif

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Fragments */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Genome_consecutive_matches_pair */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* clear_highbit and clear_lowbit */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* count_leading_zeroes and count_trailing_zeroes */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* mark mismatches */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* 32-bit shortcuts */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif


/* Genome_consec code starts here */

#if defined(DEBUG) || defined(DEBUG5)
#ifdef HAVE_SSE4_1
static void
print_vector_hex (__m128i x) {
  printf("%08X %08X %08X %08X\n",
	 _mm_extract_epi32(x,0),_mm_extract_epi32(x,1),_mm_extract_epi32(x,2),_mm_extract_epi32(x,3));
  return;
}

static void
print_vector_dec (__m128i x) {
  printf("%u %u %u %u\n",
	 _mm_extract_epi32(x,0),_mm_extract_epi32(x,1),_mm_extract_epi32(x,2),_mm_extract_epi32(x,3));
  return;
}

#elif defined(HAVE_SSE2)
static void
print_vector_hex (__m128i x) {
  printf("%08X %08X %08X %08X\n",
	 (_mm_extract_epi16(x,1) << 16) | (_mm_extract_epi16(x,0) & 0x0000FFFF),
	 (_mm_extract_epi16(x,3) << 16) | (_mm_extract_epi16(x,2) & 0x0000FFFF),
	 (_mm_extract_epi16(x,5) << 16) | (_mm_extract_epi16(x,4) & 0x0000FFFF),
	 (_mm_extract_epi16(x,7) << 16) | (_mm_extract_epi16(x,6) & 0x0000FFFF));
  return;
}

static void
print_vector_dec (__m128i x) {
  printf("%u %u %u %u\n",
	 (_mm_extract_epi16(x,1) << 16) | (_mm_extract_epi16(x,0) & 0x0000FFFF),
	 (_mm_extract_epi16(x,3) << 16) | (_mm_extract_epi16(x,2) & 0x0000FFFF),
	 (_mm_extract_epi16(x,5) << 16) | (_mm_extract_epi16(x,4) & 0x0000FFFF),
	 (_mm_extract_epi16(x,7) << 16) | (_mm_extract_epi16(x,6) & 0x0000FFFF));
  return;
}
#endif

#ifdef HAVE_AVX2
static void
print_vector_256_hex (__m256i x) {
  printf("%08X %08X %08X %08X %08X %08X %08X %08X\n",
	 _mm256_extract_epi32(x,0),_mm256_extract_epi32(x,1),_mm256_extract_epi32(x,2),_mm256_extract_epi32(x,3),
	 _mm256_extract_epi32(x,4),_mm256_extract_epi32(x,5),_mm256_extract_epi32(x,6),_mm256_extract_epi32(x,7));
  return;
}

static void
print_vector_256_dec (__m256i x) {
  printf("%u %u %u %u %u %u %u %u\n",
	 _mm256_extract_epi32(x,0),_mm256_extract_epi32(x,1),_mm256_extract_epi32(x,2),_mm256_extract_epi32(x,3),
	 _mm256_extract_epi32(x,4),_mm256_extract_epi32(x,5),_mm256_extract_epi32(x,6),_mm256_extract_epi32(x,7));
  return;
}
#endif

#ifdef HAVE_AVX512
static void
print_vector_512_hex (__m512i x) {
  unsigned int array[16];

  _mm512_store_si512((__m512i *) array,x);
  printf("%08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X\n",
	 array[0],array[1],array[2],array[3],array[4],array[5],array[6],array[7],
	 array[8],array[9],array[10],array[11],array[12],array[13],array[14],array[15]);
  return;
}

static void
print_vector_512_dec (__m512i x) {
  unsigned int array[16];

  _mm512_store_si512((__m512i *) array,x);
  printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u\n",
	 array[0],array[1],array[2],array[3],array[4],array[5],array[6],array[7],
	 array[8],array[9],array[10],array[11],array[12],array[13],array[14],array[15]);
  return;
}
#endif

#endif


#if 0
/* Note: outer unshuffle operation, as implemented below has twice the
   instruction count as lookup of reduce_nt array */
static inline UINT4
reduce_nt_unshuffle (UINT4 xhigh, UINT4 xlow) {
  UINT8 x, t;

  x = (UINT8) xhigh;
  x <<= 32;
  x |= xlow;

  t = (x ^ (x >> 1))  & 0x2222222222222222;  x = x ^ t ^ (t << 1);
  t = (x ^ (x >> 2))  & 0x0C0C0C0C0C0C0C0C;  x = x ^ t ^ (t << 2);
  t = (x ^ (x >> 4))  & 0x00F000F000F000F0;  x = x ^ t ^ (t << 4);
  t = (x ^ (x >> 8))  & 0x0000FF000000FF00;  x = x ^ t ^ (t << 8);
  t = (x ^ (x >> 16)) & 0x00000000FFFF0000;  x = x ^ t ^ (t << 16);

  return (UINT4) ((x >> 32) | x);
}
#endif


#if defined(DEBUG) || defined(DEBUG2) || defined(DEBUG5)
static void
write_chars (Genomecomp_T high, Genomecomp_T low, Genomecomp_T flags) {
  char Buffer[33];
  int i;

  Buffer[32] = '\0';
  /* printf("%08X %08X %08X => ",high,low,flags); */

  for (i = 0; i < 32; i++) {
    switch (((high & 0x01) << 1) | (low & 0x01)) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    high >>= 1;
    low >>= 1;
  }

  if (flags != 0U) {
    for (i = 0; i < 32; i++) {
      if (flags & 0x01) {
	Buffer[i] = 'N';
      }
      flags >>= 1;
    }
  }

  printf("%s",Buffer);
  return;
}
#endif



#if defined(DEBUG) || defined(DEBUG2) || defined(DEBUG5)
static void
print_blocks (Genomecomp_T *blocks, Univcoord_T startpos, Univcoord_T endpos) {
  /* Chrpos_T length = endpos - startpos; */
  Genomecomp_T *ptr;
  Univcoord_T startblocki, endblocki;
  int startcolumni, endcolumni;
  int startdiscard32, enddiscard32;
  Genomecomp_T high, low, flags;
  int i;

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  startcolumni = (startpos % 128) / 32;
  startblocki = startpos/128U*12 + startcolumni;

  endcolumni = (endpos % 128) / 32;
  endblocki = endpos/128U*12 + endcolumni;
#else
  startcolumni = (startpos % 128) / 32;
  endcolumni = (endpos % 128) / 32;

  startblocki = startpos/128U*12;
  endblocki = endpos/128U*12;
#endif

  startdiscard32 = startpos % 32;
  enddiscard32 = endpos % 32;


  ptr = &(blocks[startblocki]);
  while (ptr <= &(blocks[endblocki])) {
#if defined(WORDS_BIGENDIAN)
    high = Bigendian_convert_uint(ptr[0]);
    low = Bigendian_convert_uint(ptr[4]);
    flags = Bigendian_convert_uint(ptr[8]);
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    write_chars(high,low,flags);
    printf("\n");

    ptr += GENOME_NEXTCOL; if (++startcolumni == 4) {ptr += GENOME_NEXTROW; startcolumni = 0;}
#elif !defined(HAVE_SSE2)
    high = ptr[0]; low = ptr[4]; flags = ptr[8];
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    write_chars(high,low,flags);
    printf("\n");

    ptr += GENOME_NEXTCOL; if (++startcolumni == 4) {ptr += GENOME_NEXTROW; startcolumni = 0;}

#else
    if (startcolumni == 0) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < startdiscard32; i++) {
	printf("*");
      }
      for ( ; i < 32; i++) {
	printf(" ");
      }
      printf("\n");
      startcolumni = -1;
    }
    high = ptr[0]; low = ptr[4]; flags = ptr[8];
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    write_chars(high,low,flags);
    printf("\n");
    if (ptr == &(blocks[endblocki]) && endcolumni == 0) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < enddiscard32; i++) {
	printf(" ");
      }
      for ( ; i < 32; i++) {
	printf("*");
      }
      printf("\n");
    }


    if (startcolumni == 1) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < startdiscard32; i++) {
	printf("*");
      }
      for ( ; i < 32; i++) {
	printf(" ");
      }
      printf("\n");
      startcolumni = -1;
    }
    high = ptr[1]; low = ptr[5]; flags = ptr[9];
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    write_chars(high,low,flags);
    printf("\n");
    if (ptr == &(blocks[endblocki]) && endcolumni == 1) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < enddiscard32; i++) {
	printf(" ");
      }
      for ( ; i < 32; i++) {
	printf("*");
      }
      printf("\n");
    }

    if (startcolumni == 2) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < startdiscard32; i++) {
	printf("*");
      }
      for ( ; i < 32; i++) {
	printf(" ");
      }
      printf("\n");
      startcolumni = -1;
    }
    high = ptr[2]; low = ptr[6]; flags = ptr[10];
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    write_chars(high,low,flags);
    printf("\n");
    if (ptr == &(blocks[endblocki]) && endcolumni == 2) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < enddiscard32; i++) {
	printf(" ");
      }
      for ( ; i < 32; i++) {
	printf("*");
      }
      printf("\n");
    }

    if (startcolumni == 3) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < startdiscard32; i++) {
	printf("*");
      }
      for ( ; i < 32; i++) {
	printf(" ");
      }
      printf("\n");
      startcolumni = -1;
    }
    high = ptr[3]; low = ptr[7]; flags = ptr[11];
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    write_chars(high,low,flags);
    printf("\n");
    if (ptr == &(blocks[endblocki]) && endcolumni == 3) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < enddiscard32; i++) {
	printf(" ");
      }
      for ( ; i < 32; i++) {
	printf("*");
      }
      printf("\n");
    }

    printf("\n");
    ptr += 12;
#endif
  }

  return;
}
#endif


#if defined(USE_SHIFT_HILO) && defined(HAVE_SSE2)
static inline void
read_128_shift_lo (__m128i *__restrict__ high, __m128i *__restrict__ low, __m128i *__restrict__ flags, UINT4 *__restrict__ ptr,
		    int startcolumni) {
  __m128i a, b, c;
  
  ptr -= startcolumni;
  a = _mm_load_si128((__m128i *) ptr); ptr += 4;
  b = _mm_load_si128((__m128i *) ptr); ptr += 4;
  c = _mm_load_si128((__m128i *) ptr); ptr += 4;

  switch (startcolumni) {
  case 0:
    *high = _mm_bsrli_si128(a, 0);
    *low = _mm_bsrli_si128(b, 0);
    *flags = _mm_bsrli_si128(c, 0);
    break;
  case 1:
    *high = _mm_bsrli_si128(a, 4);
    *low = _mm_bsrli_si128(b, 4);
    *flags = _mm_bsrli_si128(c, 4);
    break;
  case 2:
    *high = _mm_bsrli_si128(a, 8);
    *low = _mm_bsrli_si128(b, 8);
    *flags = _mm_bsrli_si128(c, 8);
    break;
  default:
    *high = _mm_bsrli_si128(a, 12);
    *low = _mm_bsrli_si128(b, 12);
    *flags = _mm_bsrli_si128(c, 12);
    break;
  }

  return;
}

static inline void
read_128_shift_hi (__m128i *__restrict__ high, __m128i *__restrict__ low, __m128i *__restrict__ flags, UINT4 *__restrict__ ptr,
		   int endcolumni) {
  __m128i a, b, c;
  
  ptr -= endcolumni;
  a = _mm_load_si128((__m128i *) ptr); ptr += 4;
  b = _mm_load_si128((__m128i *) ptr); ptr += 4;
  c = _mm_load_si128((__m128i *) ptr); ptr += 4;

  switch (endcolumni) {
  case 0:
    *high = _mm_bslli_si128(a, 12);
    *low = _mm_bslli_si128(b, 12);
    *flags = _mm_bslli_si128(c, 12);
    break;
  case 1:
    *high = _mm_bslli_si128(a, 8);
    *low = _mm_bslli_si128(b, 8);
    *flags = _mm_bslli_si128(c, 8);
    break;
  case 2:
    *high = _mm_bslli_si128(a, 4);
    *low = _mm_bslli_si128(b, 4);
    *flags = _mm_bslli_si128(c, 4);
    break;
  default:
    *high = _mm_bslli_si128(a, 0);
    *low = _mm_bslli_si128(b, 0);
    *flags = _mm_bslli_si128(c, 0);
    break;
  }

  return;
}
#endif


#ifdef HAVE_SSSE3
static inline void
read_128_wrap_lo (__m128i *__restrict__ high, __m128i *__restrict__ low, __m128i *__restrict__ flags, UINT4 *__restrict__ ptr,
		  int startcolumni) {
  __m128i a, b, c, d, e, f;
  
  ptr -= startcolumni;
  a = _mm_load_si128((__m128i *) ptr); ptr += 4;
  b = _mm_load_si128((__m128i *) ptr); ptr += 4;
  c = _mm_load_si128((__m128i *) ptr); ptr += 4;
  d = _mm_load_si128((__m128i *) ptr); ptr += 4;
  e = _mm_load_si128((__m128i *) ptr); ptr += 4;
  f = _mm_load_si128((__m128i *) ptr);

  switch (startcolumni) {
  case 0:
    *high = _mm_alignr_epi8(d, a, 0);
    *low = _mm_alignr_epi8(e, b, 0);
    *flags = _mm_alignr_epi8(f, c, 0);
    break;
  case 1:
    *high = _mm_alignr_epi8(d, a, 4);
    *low = _mm_alignr_epi8(e, b, 4);
    *flags = _mm_alignr_epi8(f, c, 4);
    break;
  case 2:
    *high = _mm_alignr_epi8(d, a, 8);
    *low = _mm_alignr_epi8(e, b, 8);
    *flags = _mm_alignr_epi8(f, c, 8);
    break;
  default:
    *high = _mm_alignr_epi8(d, a, 12);
    *low = _mm_alignr_epi8(e, b, 12);
    *flags = _mm_alignr_epi8(f, c, 12);
    break;
  }

  return;
}

static inline void
read_128_wrap_hi (__m128i *__restrict__ high, __m128i *__restrict__ low, __m128i *__restrict__ flags, UINT4 *__restrict__ ptr,
		  int endcolumni) {
  __m128i a, b, c, d, e, f;
  
  ptr -= endcolumni;
  ptr -= 12;
  a = _mm_load_si128((__m128i *) ptr); ptr += 4;
  b = _mm_load_si128((__m128i *) ptr); ptr += 4;
  c = _mm_load_si128((__m128i *) ptr); ptr += 4;
  d = _mm_load_si128((__m128i *) ptr); ptr += 4;
  e = _mm_load_si128((__m128i *) ptr); ptr += 4;
  f = _mm_load_si128((__m128i *) ptr);
    
  switch (endcolumni) {
  case 0:
    *high = _mm_alignr_epi8(d, a, 4);
    *low = _mm_alignr_epi8(e, b, 4);
    *flags = _mm_alignr_epi8(f, c, 4);
    break;
  case 1:
    *high = _mm_alignr_epi8(d, a, 8);
    *low = _mm_alignr_epi8(e, b, 8);
    *flags = _mm_alignr_epi8(f, c, 8);
    break;
  case 2:
    *high = _mm_alignr_epi8(d, a, 12);
    *low = _mm_alignr_epi8(e, b, 12);
    *flags = _mm_alignr_epi8(f, c, 12);
    break;
  default:
    *high = _mm_alignr_epi8(d, a, 16);
    *low = _mm_alignr_epi8(e, b, 16);
    *flags = _mm_alignr_epi8(f, c, 16);
    break;
  }

  return;
}
#endif


#ifdef HAVE_AVX2
static inline void
read_256 (__m256i *__restrict__ high, __m256i *__restrict__ low, __m256i *__restrict__ flags, UINT4 *__restrict__ ptr) {
  __m256i a, b, c;
  a = _mm256_loadu_si256((__m256i *) ptr); /* query0_high, query0_low */
  b = _mm256_loadu_si256((__m256i *) &(ptr[8])); /* query0_flags, query1_high */
  c = _mm256_loadu_si256((__m256i *) &(ptr[16])); /* query1_low, query1_flags */

  *high = _mm256_permute2x128_si256(a, b, 0x30);
  *low = _mm256_permute2x128_si256(a, c, 0x21);
  *flags = _mm256_permute2x128_si256(b, c, 0x30);

  return;
}
#endif

#ifdef HAVE_AVX512
static inline void
read_512 (__m512i *__restrict__ high, __m512i *__restrict__ low, __m512i *__restrict__ flags, UINT4 *__restrict__ ptr) {
  __m512i a, b, c, d, e, f;
  a = _mm512_loadu_si512((__m512i *) ptr); /* query0_high, query0_low, query0_flags, query1_high */
  b = _mm512_loadu_si512((__m512i *) &(ptr[16])); /* query1_low, query1_flags, query2_high, query2_low */
  c = _mm512_loadu_si512((__m512i *) &(ptr[32])); /* query2_flags, query3_high, query3_low, query3_flags */

  d = _mm512_permutex2var_epi32(a, _mm512_setr_epi32(0, 1, 2, 3, 12, 13, 14, 15,
						     4, 5, 6, 7, 16+0, 16+1, 16+2, 16+3), b);
  e = _mm512_permutex2var_epi32(b, _mm512_setr_epi32(8, 9, 10, 11, 16+4, 16+5, 16+6, 16+7,
						     12, 13, 14, 15, 16+8, 16+9, 16+10, 16+11), c);
  f = _mm512_permutex2var_epi32(a, _mm512_setr_epi32(8, 9, 10, 11, 16+4, 16+5, 16+6, 16+7,
						     12, 13, 14, 15, 16+8, 16+9, 16+10, 16+11), b);

  *high = _mm512_permutex2var_epi64(d, _mm512_setr_epi64(0, 1, 2, 3, 8+0, 8+1, 8+2, 8+3), e);
  *low = _mm512_permutex2var_epi64(d, _mm512_setr_epi64(4, 5, 6, 7, 8+4, 8+5, 8+6, 8+7), e);
  *flags = _mm512_permutex2var_epi64(f, _mm512_setr_epi64(0, 1, 2, 3, 8+0, 8+1, 8+6, 8+7), c);

  return;
}
#endif


/* These are global values, used for alignment.  Previously for
   trimming, treated query N's as mismatches, but this is not correct
   for query N's in the middle of the read.  Also, trimming query N's
   can affect the --clip-overlap feature. */
static bool query_unk_mismatch_p = false;
static bool genome_unk_mismatch_p = true;

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#define STEP_SIZE 32
#else
/* Holds for SSE2, AVX2, and AVX512 */
#define STEP_SIZE 128
#endif


static UINT4
block_diff_standard_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  UINT4 diff;

  debug(printf("Comparing high: query %08X with genome %08X ",query_shifted[0],ref_ptr[0]));
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#else
  debug(printf("Comparing low: query %08X with genome %08X ",query_shifted[4],ref_ptr[4]));
#endif

#ifdef WORDS_BIGENDIAN
  diff = (query_shifted[0] ^ Bigendian_convert_uint(ref_ptr[0])) | (query_shifted[1] ^ Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff = (query_shifted[0] ^ ref_ptr[0]) | (query_shifted[1] ^ ref_ptr[4]);
#else
  diff = (query_shifted[0] ^ ref_ptr[0]) | (query_shifted[4] ^ ref_ptr[4]);
#endif

  /* Query Ns */
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  if (query_unk_mismatch_local_p) {
    /* Query: Considering N as a mismatch */
    diff |= query_shifted[2];
  } else {
    /* Query: Considering N as a wildcard */
    diff &= ~(query_shifted[2]);
  }
#else
  if (query_unk_mismatch_local_p) {
    /* Query: Considering N as a mismatch */
    diff |= query_shifted[8];
  } else {
    /* Query: Considering N as a wildcard */
    diff &= ~(query_shifted[8]);
  }
#endif

  /* Genome Ns */
  if (genome_unk_mismatch_p) {
    /* Genome: Considering N as a mismatch */
#ifdef WORDS_BIGENDIAN
    diff |= Bigendian_convert_uint(ref_ptr[8]);
#else
    diff |= ref_ptr[8];
#endif
  } else {
    /* Genome: Considering N as a wildcard */
#ifdef WORDS_BIGENDIAN
    diff &= ~(Bigendian_convert_uint(ref_ptr[8]));
#else
    diff &= ~(ref_ptr[8]);
#endif
  }

  debug(printf(" => diff %08X\n",diff));

  return diff;
}


#ifdef HAVE_SSE2
static __m128i
block_diff_standard_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			 bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_shifted);
  _query_low = _mm_load_si128((__m128i *) &(query_shifted[4]));
  _ref_high = _mm_load_si128((__m128i *) ref_ptr);
  _ref_low = _mm_load_si128((__m128i *) &(ref_ptr[4]));

  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) &(query_shifted[8]));
  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) &(ref_ptr[8]));
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_standard_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				  bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				  int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_shift_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_standard_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				  bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				  int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_shift_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif
#endif

#ifdef HAVE_SSSE3
static __m128i
block_diff_standard_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				 bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				 int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_wrap_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_standard_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				 bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				 int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_wrap_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX2
static __m256i
block_diff_standard_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			 bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  __m256i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_256(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_256(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  _diff = _mm256_or_si256(_mm256_xor_si256(_query_high, _ref_high), _mm256_xor_si256(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm256_or_si256(_query_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm256_or_si256(_ref_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX512
static __m512i
block_diff_standard_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			 bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  __m512i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_512(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_512(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  _diff = _mm512_or_si512(_mm512_xor_si512(_query_high, _ref_high), _mm512_xor_si512(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm512_or_si512(_query_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm512_or_si512(_ref_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_ref_flags, _diff);
  }

  return _diff;
}
#endif


/************************************************************************
 *   CMET
 ************************************************************************/

static UINT4
block_diff_metct_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		     bool query_unk_mismatch_local_p) {
  UINT4 diff;

  /* sarrayp == true */
  /* Convert everything to 3-nucleotide space */
  diff = 0U;

  /* Compare reduced C->T nts  */
#ifdef WORDS_BIGENDIAN
  diff |= ((query_shifted[0] | query_shifted[1]) ^ (Bigendian_convert_uint(ref_ptr[0]) | Bigendian_convert_uint(ref_ptr[4]))) |
    (query_shifted[1] ^ Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff |= ((query_shifted[0] | query_shifted[1]) ^ (ref_ptr[0] | ref_ptr[4])) | (query_shifted[1] ^ ref_ptr[4]);
#else
  diff |= ((query_shifted[0] | query_shifted[4]) ^ (ref_ptr[0] | ref_ptr[4])) | (query_shifted[4] ^ ref_ptr[4]);
#endif
  debug(printf(" => diff %08X\n",diff));


  /* Flags: Considering N as a mismatch */
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  if (query_unk_mismatch_local_p) {
    debug(printf("Marking query flags: query %08X ",query_shifted[8]));
    diff |= query_shifted[8];
  } else {
    debug(printf("Clearing query flags: query %08X ",query_shifted[8]));
    diff &= ~(query_shifted[8]);
  }
#else
  if (query_unk_mismatch_local_p) {
    debug(printf("Marking query flags: query %08X ",query_shifted[8]));
    diff |= query_shifted[8];
  } else {
    debug(printf("Clearing query flags: query %08X ",query_shifted[8]));
    diff &= ~(query_shifted[8]);
  }
#endif

  if (genome_unk_mismatch_p) {
    debug(printf("Marking genome flags: genome %08X ",ref_ptr[8]));
#ifdef WORDS_BIGENDIAN
    diff |= Bigendian_convert_uint(ref_ptr[8]);
#else
    diff |= (ref_ptr[8]);
#endif
  } else {
    debug(printf("Clearing genome flags: genome %08X ",ref_ptr[8]));
#ifdef WORDS_BIGENDIAN
    diff &= ~(Bigendian_convert_uint(ref_ptr[8]));
#else
    diff &= ~(ref_ptr[8]);
#endif
  }
  debug(printf(" => diff %08X\n",diff));

  return diff;
}


#ifdef HAVE_SSE2
/* Convert C to T: high/low (A) 0 0 => new high 0; (C) 0 1 => 1; (G) 1 0 => 1; (T) 1 0 => 1 */
/* new high = high | low */
static __m128i
block_diff_metct_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_shifted);
  _query_low = _mm_load_si128((__m128i *) &(query_shifted[4]));
  _ref_high = _mm_load_si128((__m128i *) ref_ptr);
  _ref_low = _mm_load_si128((__m128i *) &(ref_ptr[4]));

  /* sarrayp == true */
  /* Ignore genome-T to query-C mismatches.  Convert everything to 3-nucleotide space */
  _diff = _mm_setzero_si128();

  /* Compare reduced C->T nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_or_si128(_query_high, _query_low), _mm_or_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) &(query_shifted[8]));
  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) &(ref_ptr[8]));
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_metct_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			       bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_shift_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  /* sarrayp == true */
  /* Ignore genome-T to query-C mismatches.  Convert everything to 3-nucleotide space */
  _diff = _mm_setzero_si128();

  /* Compare reduced C->T nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_or_si128(_query_high, _query_low), _mm_or_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_metct_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			       bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_shift_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  /* sarrayp == true */
  /* Ignore genome-T to query-C mismatches.  Convert everything to 3-nucleotide space */
  _diff = _mm_setzero_si128();

  /* Compare reduced C->T nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_or_si128(_query_high, _query_low), _mm_or_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif
#endif

#ifdef HAVE_SSSE3
/* Convert C to T: high/low (A) 0 0 => new high 0; (C) 0 1 => 1; (G) 1 0 => 1; (T) 1 0 => 1 */
/* new high = high | low */
static __m128i
block_diff_metct_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_wrap_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  /* sarrayp == true */
  /* Ignore genome-T to query-C mismatches.  Convert everything to 3-nucleotide space */
  _diff = _mm_setzero_si128();

  /* Compare reduced C->T nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_or_si128(_query_high, _query_low), _mm_or_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_metct_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_wrap_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  /* sarrayp == true */
  /* Ignore genome-T to query-C mismatches.  Convert everything to 3-nucleotide space */
  _diff = _mm_setzero_si128();

  /* Compare reduced C->T nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_or_si128(_query_high, _query_low), _mm_or_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif


#ifdef HAVE_AVX2
/* Convert C to T: high/low (A) 0 0 => new high 0; (C) 0 1 => 1; (G) 1 0 => 1; (T) 1 0 => 1 */
/* new high = high | low */
static __m256i
block_diff_metct_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m256i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_256(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_256(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  /* sarrayp == true */
  /* Ignore genome-T to query-C mismatches.  Convert everything to 3-nucleotide space */
  _diff = _mm256_setzero_si256();

  /* Compare reduced C->T nts  */
  _diff = _mm256_or_si256(_diff, _mm256_xor_si256(_mm256_or_si256(_query_high, _query_low), _mm256_or_si256(_ref_high, _ref_low)));
  _diff = _mm256_or_si256(_diff, _mm256_xor_si256(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm256_or_si256(_query_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm256_or_si256(_ref_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX512
/* Convert C to T: high/low (A) 0 0 => new high 0; (C) 0 1 => 1; (G) 1 0 => 1; (T) 1 0 => 1 */
/* new high = high | low */
static __m512i
block_diff_metct_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m512i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_512(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_512(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  /* sarrayp == true */
  /* Ignore genome-T to query-C mismatches.  Convert everything to 3-nucleotide space */
  _diff = _mm512_setzero_si512();

  /* Compare reduced C->T nts  */
  _diff = _mm512_or_si512(_diff, _mm512_xor_si512(_mm512_or_si512(_query_high, _query_low), _mm512_or_si512(_ref_high, _ref_low)));
  _diff = _mm512_or_si512(_diff, _mm512_xor_si512(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm512_or_si512(_query_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm512_or_si512(_ref_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_ref_flags, _diff);
  }

  return _diff;
}
#endif


static UINT4
block_diff_metga_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		     bool query_unk_mismatch_local_p) {
  UINT4 diff;

  /* sarrayp == true */
  /* Ignore genome-A to query-G mismatches.  Convert everything to 3-nucleotide space. */
  diff = 0U;

  /* Compare reduced G->A nts  */
#ifdef WORDS_BIGENDIAN
  diff |= ((query_shifted[0] & query_shifted[1]) ^ (Bigendian_convert_uint(ref_ptr[0]) & Bigendian_convert_uint(ref_ptr[4]))) |
    (query_shifted[1] ^ Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff |= ((query_shifted[0] & query_shifted[1]) ^ (ref_ptr[0] & ref_ptr[4])) | (query_shifted[1] ^ ref_ptr[4]);
#else
  diff |= ((query_shifted[0] & query_shifted[4]) ^ (ref_ptr[0] & ref_ptr[4])) | (query_shifted[4] ^ ref_ptr[4]);
#endif
  debug(printf(" => diff %08X\n",diff));


  /* Flags: Considering N as a mismatch */
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  if (query_unk_mismatch_local_p) {
    debug(printf("Marking query flags: query %08X ",query_shifted[2]));
    diff |= query_shifted[2];
  } else {
    debug(printf("Clearing query flags: query %08X ",query_shifted[2]));
    diff &= ~(query_shifted[2]);
  }
#else
  if (query_unk_mismatch_local_p) {
    debug(printf("Marking query flags: query %08X ",query_shifted[8]));
    diff |= query_shifted[8];
  } else {
    debug(printf("Clearing query flags: query %08X ",query_shifted[8]));
    diff &= ~(query_shifted[8]);
  }
#endif

  if (genome_unk_mismatch_p) {
    debug(printf("Marking genome flags: genome %08X ",ref_ptr[8]));
#ifdef WORDS_BIGENDIAN
    diff |= Bigendian_convert_uint(ref_ptr[8]);
#else
    diff |= (ref_ptr[8]);
#endif
  } else {
    debug(printf("Clearing genome flags: genome %08X ",ref_ptr[8]));
#ifdef WORDS_BIGENDIAN
    diff &= ~(Bigendian_convert_uint(ref_ptr[8]));
#else
    diff &= ~(ref_ptr[8]);
#endif
  }
  debug(printf(" => diff %08X\n",diff));

  return diff;
}


#ifdef HAVE_SSE2
/* Convert G to A: high/low (A) 0 0 => new high 0; (C) 0 1 => 0; (G) 1 0 => 0; (T) 1 0 => 1 */
/* new high = high & low */
static __m128i
block_diff_metga_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_shifted);
  _query_low = _mm_load_si128((__m128i *) &(query_shifted[4]));
  _ref_high = _mm_load_si128((__m128i *) ref_ptr);
  _ref_low = _mm_load_si128((__m128i *) &(ref_ptr[4]));

  /* sarrayp == true */
  /* Ignore genome-A to query-G mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced G->A nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_and_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) &(query_shifted[8]));
  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) &(ref_ptr[8]));
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_metga_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			       bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_shift_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  /* sarrayp == true */
  /* Ignore genome-A to query-G mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced G->A nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_and_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_metga_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			       bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_shift_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  /* sarrayp == true */
  /* Ignore genome-A to query-G mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced G->A nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_and_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif
#endif

#ifdef HAVE_SSSE3
/* Convert G to A: high/low (A) 0 0 => new high 0; (C) 0 1 => 0; (G) 1 0 => 0; (T) 1 0 => 1 */
/* new high = high & low */
static __m128i
block_diff_metga_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_wrap_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  /* sarrayp == true */
  /* Ignore genome-A to query-G mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced G->A nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_and_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_metga_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_wrap_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  /* sarrayp == true */
  /* Ignore genome-A to query-G mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced G->A nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_and_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX2
/* Convert G to A: high/low (A) 0 0 => new high 0; (C) 0 1 => 0; (G) 1 0 => 0; (T) 1 0 => 1 */
/* new high = high & low */
static __m256i
block_diff_metga_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m256i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_256(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_256(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  /* sarrayp == true */
  /* Ignore genome-A to query-G mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm256_setzero_si256();

  /* Compare reduced G->A nts  */
  _diff = _mm256_or_si256(_diff, _mm256_xor_si256(_mm256_and_si256(_query_high, _query_low), _mm256_and_si256(_ref_high, _ref_low)));
  _diff = _mm256_or_si256(_diff, _mm256_xor_si256(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm256_or_si256(_query_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm256_or_si256(_ref_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX512
/* Convert G to A: high/low (A) 0 0 => new high 0; (C) 0 1 => 0; (G) 1 0 => 0; (T) 1 0 => 1 */
/* new high = high & low */
static __m512i
block_diff_metga_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m512i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_512(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_512(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  /* sarrayp == true */
  /* Ignore genome-A to query-G mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm512_setzero_si512();

  /* Compare reduced G->A nts  */
  _diff = _mm512_or_si512(_diff, _mm512_xor_si512(_mm512_and_si512(_query_high, _query_low), _mm512_and_si512(_ref_high, _ref_low)));
  _diff = _mm512_or_si512(_diff, _mm512_xor_si512(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm512_or_si512(_query_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm512_or_si512(_ref_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_ref_flags, _diff);
  }

  return _diff;
}
#endif


static UINT4
block_diff_cmet_sarray_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			   bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_metct_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_metga_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}

#ifdef HAVE_SSE2
static __m128i
block_diff_cmet_sarray_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_metct_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_metga_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_cmet_sarray_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				     bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				     int startcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    } else {
      return block_diff_metct_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_metct_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    } else {
      return block_diff_metga_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    }
  }
}

static __m128i
block_diff_cmet_sarray_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				     bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				     int endcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    } else {
      return block_diff_metct_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_metct_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    } else {
      return block_diff_metga_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    }
  }
}
#endif
#endif

#ifdef HAVE_SSSE3
static __m128i
block_diff_cmet_sarray_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				    bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				    int startcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    } else {
      return block_diff_metct_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_metct_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    } else {
      return block_diff_metga_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    }
  }
}

static __m128i
block_diff_cmet_sarray_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				    bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				    int endcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    } else {
      return block_diff_metct_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_metct_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    } else {
      return block_diff_metga_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    }
  }
}
#endif

#ifdef HAVE_AVX2
static __m256i
block_diff_cmet_sarray_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_metct_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_metga_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}
#endif

#ifdef HAVE_AVX512
static __m512i
block_diff_cmet_sarray_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_metct_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_metga_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}
#endif



/************************************************************************
 *   ATOI
 ************************************************************************/

static UINT4
block_diff_a2iag_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		     bool query_unk_mismatch_local_p) {
  UINT4 diff;

  /* sarrayp == true */
  /* Ignore genome-G to query-A mismatches.  Convert everything to 3-nucleotide space. */
  diff = 0U;

  /* Compare reduced A->G nts  */
#ifdef WORDS_BIGENDIAN
  diff |= ((query_shifted[0] | ~(query_shifted[1])) ^ (Bigendian_convert_uint(ref_ptr[0]) | ~(Bigendian_convert_uint(ref_ptr[4])))) |
    (query_shifted[1] ^ Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff |= ((query_shifted[0] | ~(query_shifted[1])) ^ (ref_ptr[0] | ~(ref_ptr[4]))) | (query_shifted[1] ^ ref_ptr[4]);
  /* Because (a ^ b) = (~a ^ ~b), this is equivalent to 
  diff |= ((~query_shifted[0] & query_shifted[1]) ^ (~ref_ptr[0] & ref_ptr[4])) | (query_shifted[1] ^ ref_ptr[4]);
  */
#else
  diff |= ((query_shifted[0] | ~(query_shifted[4])) ^ (ref_ptr[0] | ~(ref_ptr[4]))) | (query_shifted[4] ^ ref_ptr[4]);
#endif
  debug(printf(" => diff %08X\n",diff));

  /* Flags: Considering N as a mismatch */
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  if (query_unk_mismatch_local_p) {
    debug(printf("Marking query flags: query %08X ",query_shifted[2]));
    diff |= query_shifted[2];
  } else {
    debug(printf("Clearing query flags: query %08X ",query_shifted[2]));
    diff &= ~(query_shifted[2]);
  }
#else
  if (query_unk_mismatch_local_p) {
    debug(printf("Marking query flags: query %08X ",query_shifted[8]));
    diff |= query_shifted[8];
  } else {
    debug(printf("Clearing query flags: query %08X ",query_shifted[8]));
    diff &= ~(query_shifted[8]);
  }
#endif

  if (genome_unk_mismatch_p) {
    debug(printf("Marking genome flags: genome %08X ",ref_ptr[8]));
#ifdef WORDS_BIGENDIAN
    diff |= Bigendian_convert_uint(ref_ptr[8]);
#else
    diff |= (ref_ptr[8]);
#endif
  } else {
    debug(printf("Clearing genome flags: genome %08X ",ref_ptr[8]));
#ifdef WORDS_BIGENDIAN
    diff &= ~(Bigendian_convert_uint(ref_ptr[8]));
#else
    diff &= ~(ref_ptr[8]);
#endif
  }
  debug(printf(" => diff %08X\n",diff));

  return diff;
}


#ifdef HAVE_SSE2
/* Convert A->G: new high = high | ~low */
static __m128i
block_diff_a2iag_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_shifted);
  _query_low = _mm_load_si128((__m128i *) &(query_shifted[4]));
  _ref_high = _mm_load_si128((__m128i *) ref_ptr);
  _ref_low = _mm_load_si128((__m128i *) &(ref_ptr[4]));

  /* sarrayp == true */
  /* Ignore genome-G to query-A mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced A->G nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) &(query_shifted[8]));
  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) &(ref_ptr[8]));
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_a2iag_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			       bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_shift_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  /* sarrayp == true */
  /* Ignore genome-G to query-A mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced A->G nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_a2iag_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			       bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_shift_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  /* sarrayp == true */
  /* Ignore genome-G to query-A mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced A->G nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif
#endif

#ifdef HAVE_SSSE3
/* Convert A->G: new high = high | ~low */
static __m128i
block_diff_a2iag_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_wrap_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  /* sarrayp == true */
  /* Ignore genome-G to query-A mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced A->G nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_a2iag_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_wrap_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  /* sarrayp == true */
  /* Ignore genome-G to query-A mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm_setzero_si128();

  /* Compare reduced A->G nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX2
/* Convert A->G: new high = high | ~low */
static __m256i
block_diff_a2iag_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m256i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_256(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_256(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  /* sarrayp == true */
  /* Ignore genome-G to query-A mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm256_setzero_si256();

  /* Compare reduced A->G nts  */
  _diff = _mm256_or_si256(_diff, _mm256_xor_si256(_mm256_andnot_si256(_query_high, _query_low), _mm256_andnot_si256(_ref_high, _ref_low)));
  _diff = _mm256_or_si256(_diff, _mm256_xor_si256(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm256_or_si256(_query_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm256_or_si256(_ref_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX512
/* Convert A->G: new high = high | ~low */
static __m512i
block_diff_a2iag_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m512i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_512(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_512(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  /* sarrayp == true */
  /* Ignore genome-G to query-A mismatches.  Convert everything to 3-nucleotide space. */
  _diff = _mm512_setzero_si512();

  /* Compare reduced A->G nts  */
  _diff = _mm512_or_si512(_diff, _mm512_xor_si512(_mm512_andnot_si512(_query_high, _query_low), _mm512_andnot_si512(_ref_high, _ref_low)));
  _diff = _mm512_or_si512(_diff, _mm512_xor_si512(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm512_or_si512(_query_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm512_or_si512(_ref_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_ref_flags, _diff);
  }

  return _diff;
}
#endif


static UINT4
block_diff_a2itc_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		     bool query_unk_mismatch_local_p) {
  UINT4 diff;

  /* sarrayp == true */
  /* Ignore genome-C to query-T mismatches */
  diff = 0U;

  /* Compare reduced T->C nts  */
#ifdef WORDS_BIGENDIAN
  diff |= ((query_shifted[0] & ~(query_shifted[1])) ^ (Bigendian_convert_uint(ref_ptr[0]) & ~(Bigendian_convert_uint(ref_ptr[4])))) |
    (query_shifted[1] ^ Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff |= ((query_shifted[0] & ~(query_shifted[1])) ^ (ref_ptr[0] & ~(ref_ptr[4]))) | (query_shifted[1] ^ ref_ptr[4]);
#else
  diff |= ((query_shifted[0] & ~(query_shifted[4])) ^ (ref_ptr[0] & ~(ref_ptr[4]))) | (query_shifted[4] ^ ref_ptr[4]);
#endif
  debug(printf(" => diff %08X\n",diff));

  /* Flags: Considering N as a mismatch */
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  if (query_unk_mismatch_local_p) {
    debug(printf("Marking query flags: query %08X ",query_shifted[8]));
    diff |= query_shifted[8];
  } else {
    debug(printf("Clearing query flags: query %08X ",query_shifted[8]));
    diff &= ~(query_shifted[8]);
  }
#else
  if (query_unk_mismatch_local_p) {
    debug(printf("Marking query flags: query %08X ",query_shifted[2]));
    diff |= query_shifted[2];
  } else {
    debug(printf("Clearing query flags: query %08X ",query_shifted[2]));
    diff &= ~(query_shifted[2]);
  }
#endif

  if (genome_unk_mismatch_p) {
    debug(printf("Marking genome flags: genome %08X ",ref_ptr[8]));
#ifdef WORDS_BIGENDIAN
    diff |= Bigendian_convert_uint(ref_ptr[8]);
#else
    diff |= (ref_ptr[8]);
#endif
  } else {
    debug(printf("Clearing genome flags: genome %08X ",ref_ptr[8]));
#ifdef WORDS_BIGENDIAN
    diff &= ~(Bigendian_convert_uint(ref_ptr[8]));
#else
    diff &= ~(ref_ptr[8]);
#endif
  }
  debug(printf(" => diff %08X\n",diff));

  return diff;
}


#ifdef HAVE_SSE2
/* Convert T->C: new high = high & ~low */
static __m128i
block_diff_a2itc_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_shifted);
  _query_low = _mm_load_si128((__m128i *) &(query_shifted[4]));
  _ref_high = _mm_load_si128((__m128i *) ref_ptr);
  _ref_low = _mm_load_si128((__m128i *) &(ref_ptr[4]));

  /* sarrayp == true */
  /* Ignore genome-C to query-T mismatches */
  _diff = _mm_setzero_si128();

  /* Compare reduced T->C nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_low, _query_high), _mm_andnot_si128(_ref_low, _ref_high)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) &(query_shifted[8]));
  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) &(ref_ptr[8]));
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_a2itc_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			       bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_shift_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  /* sarrayp == true */
  /* Ignore genome-C to query-T mismatches */
  _diff = _mm_setzero_si128();

  /* Compare reduced T->C nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_low, _query_high), _mm_andnot_si128(_ref_low, _ref_high)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_a2itc_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			       bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_shift_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_shift_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  /* sarrayp == true */
  /* Ignore genome-C to query-T mismatches */
  _diff = _mm_setzero_si128();

  /* Compare reduced T->C nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_low, _query_high), _mm_andnot_si128(_ref_low, _ref_high)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif
#endif

#ifdef HAVE_SSSE3
/* Convert T->C: new high = high & ~low */
static __m128i
block_diff_a2itc_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_wrap_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);

  /* sarrayp == true */
  /* Ignore genome-C to query-T mismatches */
  _diff = _mm_setzero_si128();

  /* Compare reduced T->C nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_low, _query_high), _mm_andnot_si128(_ref_low, _ref_high)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}

static __m128i
block_diff_a2itc_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_128_wrap_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_wrap_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);

  /* sarrayp == true */
  /* Ignore genome-C to query-T mismatches */
  _diff = _mm_setzero_si128();

  /* Compare reduced T->C nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_low, _query_high), _mm_andnot_si128(_ref_low, _ref_high)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX2
/* Convert T->C: new high = high & ~low */
static __m256i
block_diff_a2itc_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m256i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_256(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_256(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  /* sarrayp == true */
  /* Ignore genome-C to query-T mismatches */
  _diff = _mm256_setzero_si256();

  /* Compare reduced T->C nts  */
  _diff = _mm256_or_si256(_diff, _mm256_xor_si256(_mm256_andnot_si256(_query_low, _query_high), _mm256_andnot_si256(_ref_low, _ref_high)));
  _diff = _mm256_or_si256(_diff, _mm256_xor_si256(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm256_or_si256(_query_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm256_or_si256(_ref_flags, _diff);
  } else {
    _diff = _mm256_andnot_si256(_ref_flags, _diff);
  }

  return _diff;
}
#endif

#ifdef HAVE_AVX512
/* Convert T->C: new high = high & ~low */
static __m512i
block_diff_a2itc_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
		      bool query_unk_mismatch_local_p) {
  __m512i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  read_512(&_query_high,&_query_low,&_query_flags,query_shifted);
  read_512(&_ref_high,&_ref_low,&_ref_flags,ref_ptr);

  /* sarrayp == true */
  /* Ignore genome-C to query-T mismatches */
  _diff = _mm512_setzero_si512();

  /* Compare reduced T->C nts  */
  _diff = _mm512_or_si512(_diff, _mm512_xor_si512(_mm512_andnot_si512(_query_low, _query_high), _mm512_andnot_si512(_ref_low, _ref_high)));
  _diff = _mm512_or_si512(_diff, _mm512_xor_si512(_query_low, _ref_low));

  if (query_unk_mismatch_local_p) {
    _diff = _mm512_or_si512(_query_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_query_flags, _diff);
  }

  if (genome_unk_mismatch_p) {
    _diff = _mm512_or_si512(_ref_flags, _diff);
  } else {
    _diff = _mm512_andnot_si512(_ref_flags, _diff);
  }

  return _diff;
}
#endif


static UINT4
block_diff_atoi_sarray_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			   bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2iag_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2itc_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}

#ifdef HAVE_SSE2
static __m128i
block_diff_atoi_sarray_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2iag_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2itc_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_atoi_sarray_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				     bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				     int startcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    } else {
      return block_diff_a2iag_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    } else {
      return block_diff_a2itc_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    }
  }
}

static __m128i
block_diff_atoi_sarray_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				     bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				     int endcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    } else {
      return block_diff_a2iag_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    } else {
      return block_diff_a2itc_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    }
  }
}
#endif
#endif

#ifdef HAVE_SSSE3
static __m128i
block_diff_atoi_sarray_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				    bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				    int startcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    } else {
      return block_diff_a2iag_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    } else {
      return block_diff_a2itc_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    }
  }
}

static __m128i
block_diff_atoi_sarray_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				    bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				    int endcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    } else {
      return block_diff_a2iag_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    } else {
      return block_diff_a2itc_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    }
  }
}
#endif

#ifdef HAVE_AVX2
static __m256i
block_diff_atoi_sarray_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2iag_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2itc_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}
#endif

#ifdef HAVE_AVX512
static __m512i
block_diff_atoi_sarray_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2iag_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2itc_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}
#endif



/************************************************************************
 *  TTOC
 ************************************************************************/

static UINT4
block_diff_ttoc_sarray_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			   bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2itc_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2iag_32(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}

#ifdef HAVE_SSE2
static __m128i
block_diff_ttoc_sarray_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2itc_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2iag_128(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_ttoc_sarray_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				     bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				     int startcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    } else {
      return block_diff_a2itc_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    } else {
      return block_diff_a2iag_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   startcolumni);
    }
  }
}

static __m128i
block_diff_ttoc_sarray_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				     bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				     int endcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    } else {
      return block_diff_a2itc_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    } else {
      return block_diff_a2iag_128_shift_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    }
  }
}
#endif
#endif

#ifdef HAVE_SSSE3
static __m128i
block_diff_ttoc_sarray_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				    bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				    int startcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    } else {
      return block_diff_a2itc_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    } else {
      return block_diff_a2iag_128_wrap_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  startcolumni);
    }
  }
}

static __m128i
block_diff_ttoc_sarray_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
				    bool plusp, int genestrand, bool query_unk_mismatch_local_p,
				    int endcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    } else {
      return block_diff_a2itc_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    } else {
      return block_diff_a2iag_128_wrap_hi(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					  endcolumni);
    }
  }
}
#endif

#ifdef HAVE_AVX2
static __m256i
block_diff_ttoc_sarray_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2itc_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2iag_256(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}
#endif

#ifdef HAVE_AVX512
static __m512i
block_diff_ttoc_sarray_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			    bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2itc_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    } else {
      return block_diff_a2iag_512(query_shifted,ref_ptr,query_unk_mismatch_local_p);
    }
  }
}
#endif


/* query_shifted, ref_ptr, plusp, genestrand, query_unk_mismatch_local_p */
#ifdef HAVE_AVX512
typedef __m512i (*Diffproc_512_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool);
#endif

#ifdef HAVE_AVX2
typedef __m256i (*Diffproc_256_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool);
#endif

#ifdef HAVE_SSSE3
typedef __m128i (*Diffproc_128_wrap_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool, int);
#endif

#ifdef HAVE_SSE2
#ifdef USE_SHIFT_HILO
typedef __m128i (*Diffproc_128_shift_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool, int);
#endif

typedef __m128i (*Diffproc_128_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool);
#endif

typedef UINT4 (*Diffproc_32_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool);


/* For CMET and ATOI, ignores genome-to-query mismatches.  Used by
   Genome_consecutive procedures, called only by sarray-read.c */
#ifdef HAVE_AVX512
static Diffproc_512_T block_diff_sarray_512;
#endif
#ifdef HAVE_AVX2
static Diffproc_256_T block_diff_sarray_256; 
#endif
#ifdef HAVE_SSSE3
static Diffproc_128_wrap_T block_diff_sarray_128_wrap_lo;
static Diffproc_128_wrap_T block_diff_sarray_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
static Diffproc_128_T block_diff_sarray_128; 
#ifdef USE_SHIFT_HILO
static Diffproc_128_shift_T block_diff_sarray_128_shift_lo;
static Diffproc_128_shift_T block_diff_sarray_128_shift_hi;
#endif
#endif
static Diffproc_32_T block_diff_sarray_32; 

#ifdef HAVE_AVX512
static __m512i _BOUND_HIGH_512;
static __m512i _BOUND_LOW_512;
#endif
#ifdef HAVE_AVX2
static __m256i _BOUND_HIGH_256;
static __m256i _BOUND_LOW_256;
#endif
#ifdef HAVE_SSE2
static __m128i _BOUND_HIGH;
static __m128i _BOUND_LOW;
#endif

void
Genome_consec_setup (bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		     Mode_T mode) {

#ifdef HAVE_AVX512
  _BOUND_HIGH_512 = _mm512_set_epi32(512,480,448,416, 384,352,320,288, 256,224,192,160, 128,96,64,32);
  _BOUND_LOW_512 = _mm512_set_epi32(480,448,416,384, 352,320,288,256, 224,192,160,128, 96,64,32,0);
#endif
#ifdef HAVE_AVX2
  _BOUND_HIGH_256 = _mm256_set_epi32(256,224,192,160,128,96,64,32);
  _BOUND_LOW_256 = _mm256_set_epi32(224,192,160,128,96,64,32,0);
#endif
#ifdef HAVE_SSE2
  _BOUND_HIGH = _mm_set_epi32(128,96,64,32);
  _BOUND_LOW = _mm_set_epi32(96,64,32,0);
#endif

  query_unk_mismatch_p = query_unk_mismatch_p_in;
  genome_unk_mismatch_p = genome_unk_mismatch_p_in;

  switch (mode) {
  case STANDARD:
#ifdef HAVE_AVX512
    block_diff_sarray_512 = block_diff_standard_512;
#endif
#ifdef HAVE_AVX2
    block_diff_sarray_256 = block_diff_standard_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_sarray_128_wrap_lo = block_diff_standard_128_wrap_lo;
    block_diff_sarray_128_wrap_hi = block_diff_standard_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_sarray_128 = block_diff_standard_128;
#ifdef USE_SHIFT_HILO
    block_diff_sarray_128_shift_lo = block_diff_standard_128_shift_lo;
    block_diff_sarray_128_shift_hi = block_diff_standard_128_shift_hi;
#endif
#endif
    block_diff_sarray_32 = block_diff_standard_32;
    break;

  case CMET_STRANDED: case CMET_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_sarray_512 = block_diff_cmet_sarray_512;
#endif
#ifdef HAVE_AVX2
    block_diff_sarray_256 = block_diff_cmet_sarray_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_sarray_128_wrap_lo = block_diff_cmet_sarray_128_wrap_lo;
    block_diff_sarray_128_wrap_hi = block_diff_cmet_sarray_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_sarray_128 = block_diff_cmet_sarray_128;
#ifdef USE_SHIFT_HILO
    block_diff_sarray_128_shift_lo = block_diff_cmet_128_shift_lo;
    block_diff_sarray_128_shift_hi = block_diff_cmet_128_shift_hi;
#endif
#endif
    block_diff_sarray_32 = block_diff_cmet_sarray_32;
    break;

  case ATOI_STRANDED: case ATOI_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_sarray_512 = block_diff_atoi_sarray_512;
#endif
#ifdef HAVE_AVX2
    block_diff_sarray_256 = block_diff_atoi_sarray_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_sarray_128_wrap_lo = block_diff_atoi_sarray_128_wrap_lo;
    block_diff_sarray_128_wrap_hi = block_diff_atoi_sarray_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_sarray_128 = block_diff_atoi_sarray_128;
#ifdef USE_SHIFT_HILO
    block_diff_sarray_128_shift_lo = block_diff_atoi_128_shift_lo;
    block_diff_sarray_128_shift_hi = block_diff_atoi_128_shift_hi;
#endif
#endif
    block_diff_sarray_32 = block_diff_atoi_sarray_32;
    break;

  case TTOC_STRANDED: case TTOC_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_sarray_512 = block_diff_ttoc_sarray_512;
#endif
#ifdef HAVE_AVX2
    block_diff_sarray_256 = block_diff_ttoc_sarray_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_sarray_128_wrap_lo = block_diff_ttoc_sarray_128_wrap_lo;
    block_diff_sarray_128_wrap_hi = block_diff_ttoc_sarray_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_sarray_128 = block_diff_ttoc_sarray_128;
#ifdef USE_SHIFT_HILO
    block_diff_sarray_128_shift_lo = block_diff_ttoc_128_shift_lo;
    block_diff_sarray_128_shift_hi = block_diff_ttoc_128_shift_hi;
#endif
#endif
    block_diff_sarray_32 = block_diff_ttoc_sarray_32;
    break;

  default: fprintf(stderr,"Mode %d not recognized\n",mode); abort();
  }

  return;
}


/************************************************************************/

#ifdef HAVE_AVX512
/* Need to implement.  Extract procedures not available. */
#endif


#ifdef HAVE_AVX2
#define nonzero_p_256(diff) !_mm256_testz_si256(diff,diff)


#if defined(HAVE_POPCNT)
#define popcount_ones_256(_diff) (_popcnt64(_mm256_extract_epi64(_diff,0)) + _popcnt64(_mm256_extract_epi64(_diff,1)) + popcnt64(_mm256_extract_epi64(_diff,2)) + _popcnt64(_mm256_extract_epi64(_diff,3)))
#elif defined(HAVE_MM_POPCNT)
#define popcount_ones_256(_diff) (_mm_popcnt_u64(_mm256_extract_epi64(_diff,0)) + _mm_popcnt_u64(_mm256_extract_epi64(_diff,1)) + _mm_popcnt_u64(_mm256_extract_epi64(_diff,2)) + _mm_popcnt_u64(_mm256_extract_epi64(_diff,3)))
#else
#define popcount_ones_256(_diff) (__builtin_popcountll(_mm256_extract_epi64(_diff,0)) + __builtin_popcountll(_mm256_extract_epi64(_diff,1)) + __builtin_popcountll(_mm256_extract_epi64(_diff,2)) + __builtin_popcountll(_mm256_extract_epi64(_diff,3)))
#endif

static int
count_leading_zeroes_256 (__m256i _diff) {
  debug4(printf("Entered count_leading_zeroes with "));
  debug4(print_vector_256_hex(_diff));
  UINT8 x;

#ifdef HAVE_LZCNT
  if ((x = _mm256_extract_epi64(_diff,3)) != 0) {
    return (int) _lzcnt_u64(x);
  } else if ((x = _mm256_extract_epi64(_diff,2)) != 0) {
    return 64 + (int) _lzcnt_u64(x);
  } else if ((x = _mm256_extract_epi64(_diff,1)) != 0) {
    return 128 + (int) _lzcnt_u64(x);
  } else {
    return 192 + (int) _lzcnt_u64(_mm256_extract_epi64(_diff,0));
  }

#elif defined(HAVE_BUILTIN_CLZ)
  if ((x = _mm256_extract_epi64(_diff,3)) != 0) {
    return (int) __builtin_clzll(x);
  } else if ((x = _mm256_extract_epi64(_diff,2)) != 0) {
    return 64 + (int) __builtin_clzll(x);
  } else if ((x = _mm256_extract_epi64(_diff,1)) != 0) {
    return 128 + (int) __builtin_clzll(x);
  } else {
    return 192 + (int) __builtin_clzll(_mm256_extract_epi64(_diff,0));
  }

#else
  abort();
#endif
}

static int
count_trailing_zeroes_256 (__m256i _diff) {
  debug4(printf("Entered count_trailing_zeroes with "));
  debug4(print_vector_256_hex(_diff));
  UINT8 x;

#ifdef HAVE_TZCNT
  if ((x = _mm256_extract_epi64(_diff,0)) != 0) {
    return (int) _tzcnt_u64(x);
  } else if ((x = _mm256_extract_epi64(_diff,1)) != 0) {
    return 64 + (int) _tzcnt_u64(x);
  } else if ((x = _mm256_extract_epi64(_diff,2)) != 0) {
    return 128 + (int) _tzcnt_u64(x);
  } else {
    return 192 + (int) _tzcnt_u64(_mm256_extract_epi64(_diff,3));
  }

#elif defined(HAVE_BUILTIN_CTZ)
  if ((x = _mm256_extract_epi64(_diff,0)) != 0) {
    return (int) __builtin_ctzll(x);
  } else if ((x = _mm256_extract_epi64(_diff,1)) != 0) {
    return 64 + (int) __builtin_ctzll(x);
  } else if ((x = _mm256_extract_epi64(_diff,2)) != 0) {
    return 128 + (int) __builtin_ctzll(x);
  } else {
    return 192 + (int) __builtin_ctzll(_mm256_extract_epi64(_diff,3));
  }

#else
  abort();
#endif

}

static __m256i
clear_highbit_256 (__m256i _diff, int leading_zeroes) {
  __m256i _subtract, _relpos;
  int relpos;

  relpos = 255 - leading_zeroes;
  debug3(printf("Clearing high bit at relpos %d\n",relpos));

  _subtract = _mm256_slli_epi32(_mm256_set1_epi32(1), relpos % 32);
  _relpos = _mm256_set1_epi32(relpos);
  _subtract = _mm256_and_si256(_mm256_cmpgt_epi32(_BOUND_HIGH_256, _relpos), _subtract);
  _subtract = _mm256_andnot_si256(_mm256_cmpgt_epi32(_BOUND_LOW_256, _relpos), _subtract);

  debug3(printf("Subtract: "));
  debug3(print_vector_256_hex(_subtract));
#if 0
  /* latency 1, throughput: 0.5 */
  return _mm256_sub_epi32(_diff, _subtract);
#else
  /* _mm256_xor_si128 also works if all other bits are 0.  latency 1, throughput: 0.33 */
  return _mm256_xor_si256(_diff, _subtract);
#endif
}

/* relpos is equal to trailing_zeroes */
static __m256i
clear_lowbit_256 (__m256i _diff, int relpos) {
  __m256i _subtract, _relpos;

  debug3(printf("Clearing low bit at relpos %d\n",relpos));

  _subtract = _mm256_slli_epi32(_mm256_set1_epi32(1), relpos % 32);
  _relpos = _mm256_set1_epi32(relpos);
  _subtract = _mm256_and_si256(_mm256_cmpgt_epi32(_BOUND_HIGH_256, _relpos), _subtract);
  _subtract = _mm256_andnot_si256(_mm256_cmpgt_epi32(_BOUND_LOW_256, _relpos), _subtract);

  debug3(printf("Subtract: "));
  debug3(print_vector_256_hex(_subtract));
#if 0
  /* latency 1, throughput: 0.5 */
  return _mm256_sub_epi32(_diff, _subtract);
#else
  /* _mm256_xor_si128 also works if all other bits are 0.  latency 1, throughput: 0.33 */
  return _mm256_xor_si256(_diff, _subtract);
#endif
}

#endif


#ifdef HAVE_SSE2

#ifdef HAVE_SSE4_1
#define nonzero_p_128(diff) !_mm_testz_si128(diff,diff)
#else
#define nonzero_p_128(diff) _mm_movemask_epi8(_mm_cmpeq_epi8(diff,_mm_setzero_si128())) != 0xFFFF
#endif

#if (defined(USE_SHIFT_FIRST_MISMATCH) && defined(HAVE_SSE2)) || (defined(USE_WRAP_FIRST_MISMATCH) && defined(HAVE_SSSE3))
static __m128i
clear_start_128 (__m128i _diff, int startdiscard) {
  __m128i _mask, _startdiscard;
#ifdef DEBUG
  __m128i _result;
#endif

  debug(printf("Clearing start to startdiscard %d\n",startdiscard));
  debug(printf("Before: "));
  debug(print_vector_hex(_diff));

#ifdef DEFECTIVE_SSE2_COMPILER
  _mask = _mm_sll_epi32(_mm_set1_epi32(~0U), _mm_setr_epi32(startdiscard % 32,0,0,0));
#else
  _mask = _mm_slli_epi32(_mm_set1_epi32(~0U), startdiscard % 32);
#endif
  _startdiscard = _mm_set1_epi32(startdiscard);
  _mask = _mm_or_si128(_mask, _mm_cmplt_epi32(_startdiscard, _BOUND_LOW));
  _mask = _mm_and_si128(_mask, _mm_cmplt_epi32(_startdiscard, _BOUND_HIGH));

#ifdef DEBUG
  _result = _mm_and_si128(_mask, _diff);
  debug(printf("After:  "));
  print_vector_hex(_result);
#endif

  return _mm_and_si128(_mask, _diff);
}
#endif

#if (defined(USE_SHIFT_FIRST_MISMATCH) && defined(HAVE_SSE2)) || (defined(USE_WRAP_FIRST_MISMATCH) && defined(HAVE_SSSE3))
static __m128i
clear_end_128 (__m128i _diff, int enddiscard) {
  __m128i _mask, _enddiscard;
#ifdef DEBUG
  __m128i _result;
#endif

  debug(printf("Clearing end from enddiscard %d\n",enddiscard));
  debug(printf("Before: "));
  debug(print_vector_hex(_diff));

#ifdef DEFECTIVE_SSE2_COMPILER
  _mask = _mm_sll_epi32(_mm_set1_epi32(~0U), _mm_setr_epi32(enddiscard % 32,0,0,0));
#else
  _mask = _mm_slli_epi32(_mm_set1_epi32(~0U), enddiscard % 32);
#endif
  _enddiscard = _mm_set1_epi32(enddiscard);
  _mask = _mm_or_si128(_mm_cmplt_epi32(_enddiscard, _BOUND_LOW), _mask);
  _mask = _mm_and_si128(_mm_cmplt_epi32(_enddiscard, _BOUND_HIGH), _mask);

#ifdef DEBUG
  _result = _mm_andnot_si128(_mask, _diff);
  debug(printf("After:  "));
  print_vector_hex(_result);
#endif

  return _mm_andnot_si128(_mask, _diff);
}
#endif



#if !defined(HAVE_SSE4_2) || !defined(HAVE_MM_EXTRACT_EPI64)

#if 0
/* Naive method for pre-SSE4.2.  Requires four popcount operations. */
static int
popcount_ones_128 (__m128i _diff) {
  UINT4 diff[4];

  _mm_store_si128((__m128i *) diff,_diff);

  return __builtin_popcount(diff[0]) + __builtin_popcount(diff[1]) + __builtin_popcount(diff[2]) + __builtin_popcount(diff[3]);
}
#endif


/************************************************************************
 *  Method for pre-SSE4.2: Using Harley's method to reduce number of
 *  popcount operations when we need to compute four 32-bit popcounts
 *  in a 128-bit register.
 *
 *  The naive method is to do _popcnt32(diff[0]) + _popcnt32(diff[1]) + _popcnt32(diff[2]) + _popcnt32(diff[3]);
 *
 *  Harley's method uses a carry-save adder (CSA) to reduce the number
 *  of popcount operations for 3 elements from 3 to 2.
 ************************************************************************/

#define CSA(h,l, a,b,c, u,v) u = a ^ b; v = c; h = (a & b) | (u & v); l = u ^ v;

static int
popcount_ones_128 (__m128i _diff) {
  UINT4 ones, twos, u, v;
  UINT4 diff[4];

  _mm_store_si128((__m128i *) diff,_diff);

  CSA(twos, ones, diff[0], diff[1], diff[2], u, v);

  return 2*__builtin_popcount(twos) + __builtin_popcount(ones) + __builtin_popcount(diff[3]);
}


#elif defined(HAVE_POPCNT)
#define popcount_ones_128(_diff) (_popcnt64(_mm_extract_epi64(_diff,0)) + _popcnt64(_mm_extract_epi64(_diff,1)))
#elif defined(HAVE_MM_POPCNT) && defined(HAVE_MM_POPCNT_U64)
#define popcount_ones_128(_diff) (_mm_popcnt_u64(_mm_extract_epi64(_diff,0)) + _mm_popcnt_u64(_mm_extract_epi64(_diff,1)))
#else
#define popcount_ones_128(_diff) (__builtin_popcountll(_mm_extract_epi64(_diff,0)) + __builtin_popcountll(_mm_extract_epi64(_diff,1)))
#endif


static int
count_leading_zeroes_128 (__m128i _diff) {
  debug4(printf("Entered count_leading_zeroes_128 with "));
  debug4(print_vector_hex(_diff));

#if defined(HAVE_SSE4_2) && defined(HAVE_LZCNT) && defined(HAVE_MM_EXTRACT_EPI64)
  UINT8 x;

  if ((x = _mm_extract_epi64(_diff,1)) != 0) {
    return (int) _lzcnt_u64(x);
  } else {
    return 64 + (int) _lzcnt_u64(_mm_extract_epi64(_diff,0));
  }

#elif defined(HAVE_SSE4_1) && defined(HAVE_MM_EXTRACT_EPI64)
  UINT8 x;

  if ((x = _mm_extract_epi64(_diff,1)) != 0) {
    return __builtin_clzll(x);
  } else {
    return 64 + __builtin_clzll(_mm_extract_epi64(_diff,0));
  }

#else
  UINT4 x;

  if ((x = (_mm_extract_epi16(_diff,7) << 16) | (_mm_extract_epi16(_diff,6) & 0x0000FFFF)) != 0) {
    debug4(printf("word 3 is non-empty, so returning %d vs %d\n",__builtin_clz(x),__builtin_clz(((UINT4 *) &_diff)[3])));
    return __builtin_clz(x);
  } else if ((x = (_mm_extract_epi16(_diff,5) << 16) | (_mm_extract_epi16(_diff,4) & 0x0000FFFF)) != 0) {
    debug4(printf("word 2 is non-empty, so returning 32 + %d vs %d\n",__builtin_clz(x),__builtin_clz(((UINT4 *) &_diff)[2])));
    return 32 + __builtin_clz(x);
  } else if ((x = (_mm_extract_epi16(_diff,3) << 16) | (_mm_extract_epi16(_diff,2) & 0x0000FFFF)) != 0) {
    debug4(printf("word 1 is non-empty, so returning 64 + %d vs %d\n",__builtin_clz(x),__builtin_clz(((UINT4 *) &_diff)[1])));
    return 64 + __builtin_clz(x);
  } else {
    x = (_mm_extract_epi16(_diff,1) << 16) | (_mm_extract_epi16(_diff,0) & 0x0000FFFF);
    debug4(printf("word 0 is non-empty, so returning 96 + %d vs %d\n",__builtin_clz(x),__builtin_clz(((UINT4 *) &_diff)[0])));
    return 96 + __builtin_clz(x);
  }
#endif
}

static int
count_trailing_zeroes_128 (__m128i _diff) {
  debug4(printf("Entered count_trailing_zeroes_128 with "));
  debug4(print_vector_hex(_diff));

#if defined(HAVE_SSE4_2) && defined(HAVE_TZCNT) && defined(HAVE_MM_EXTRACT_EPI64)
  UINT8 x;

  if ((x = _mm_extract_epi64(_diff,0)) != 0) {
    return (int) _tzcnt_u64(x);
  } else {
    return 64 + (int) _tzcnt_u64(_mm_extract_epi64(_diff,1));
  }

#elif defined(HAVE_SSE4_1) && defined(HAVE_MM_EXTRACT_EPI64)
  UINT8 x;

  if ((x = _mm_extract_epi64(_diff,0)) != 0) {
    return __builtin_ctzll(x);
  } else {
    return 64 + __builtin_ctzll(_mm_extract_epi64(_diff,1));
  }

#else
  UINT4 x;

  if ((x = (_mm_extract_epi16(_diff,1) << 16) | (_mm_extract_epi16(_diff,0) & 0x0000FFFF)) != 0) {
    debug4(printf("word 0 is non-empty, so returning %d vs %d\n",__builtin_ctz(x),__builtin_ctz(((UINT4 *) &_diff)[0])));
    return __builtin_ctz(x);
  } else if ((x = (_mm_extract_epi16(_diff,3) << 16) | (_mm_extract_epi16(_diff,2) & 0x0000FFFF)) != 0) {
    debug4(printf("word 1 is non-empty, so returning 32 + %d vs %d\n",__builtin_ctz(x),__builtin_ctz(((UINT4 *) &_diff)[1])));
    return 32 + __builtin_ctz(x);
  } else if ((x = (_mm_extract_epi16(_diff,5) << 16) | (_mm_extract_epi16(_diff,4) & 0x0000FFFF)) != 0) {
    debug4(printf("word 2 is non-empty, so returning 64 + %d vs %d\n",__builtin_ctz(x),__builtin_ctz(((UINT4 *) &_diff)[2])));
    return 64 + __builtin_ctz(x);
  } else {
    x = (_mm_extract_epi16(_diff,7) << 16) | (_mm_extract_epi16(_diff,6) & 0x0000FFFF);
    debug4(printf("word 3 is non-empty, so returning 96 + %d vs %d\n",__builtin_ctz(x),__builtin_ctz(((UINT4 *) &_diff)[3])));
    return 96 + __builtin_ctz(x);
  }
#endif
}

#endif


/*                 76543210 */
#define HIGH_BIT 0x80000000

#define nonzero_p_32(diff) diff

#define clear_start_32(diff,startdiscard) (diff & (~0U << (startdiscard)))
#define clear_end_32(diff,enddiscard) (diff & ~(~0U << (enddiscard)))

/* For trimming */
#define set_start_32(diff,startdiscard) (diff | ~(~0U << startdiscard))
#define set_end_32(diff,enddiscard) (diff | (~0U << enddiscard))

/* For fragment functions that evaluate only the end 16-mer */
#define clear_start_mask(startdiscard) (~0U << (startdiscard))
#define clear_end_mask(enddiscard) (~(~0U << (enddiscard)))


/* Same speed: clear_highbit(diff,relpos) (diff - (HIGH_BIT >> relpos)) */
/* Note: xor assumes that bit at relpos was on */
#define clear_highbit_32(diff,relpos) (diff ^ (HIGH_BIT >> relpos))

/* Slower: clear_lowbit(diff,relpos) diff -= (1 << relpos) */
#define clear_lowbit_32(diff,relpos) (diff & (diff - 1));


#if !defined(HAVE_SSE4_2)
#define popcount_ones_32(diff) (count_bits[diff & 0x0000FFFF] + count_bits[diff >> 16])
#elif defined(HAVE_POPCNT)
#define popcount_ones_32(diff) (_popcnt32(diff))
#elif defined(HAVE_MM_POPCNT)
#define popcount_ones_32(diff) (_mm_popcnt_u32(diff))
#elif defined(HAVE_BUILTIN_POPCOUNT)
#define popcount_ones_32(diff) (__builtin_popcount(diff))
#else
#define popcount_ones_32(diff) (count_bits[diff & 0x0000FFFF] + count_bits[diff >> 16])
#endif

#if !defined(HAVE_SSE4_2)
#define count_leading_zeroes_32(diff) ((diff >> 16) ? clz_table[diff >> 16] : 16 + clz_table[diff])
#elif defined(HAVE_LZCNT)
#define count_leading_zeroes_32(diff) _lzcnt_u32(diff)
#elif defined(HAVE_BUILTIN_CLZ)
#define count_leading_zeroes_32(diff) __builtin_clz(diff)
#else
#define count_leading_zeroes_32(diff) ((diff >> 16) ? clz_table[diff >> 16] : 16 + clz_table[diff])
#endif

#if !defined(HAVE_SSE4_2)
#define count_trailing_zeroes_32(diff) mod_37_bit_position[(-diff & diff) % 37]
#elif defined(HAVE_TZCNT)
#define count_trailing_zeroes_32(diff) _tzcnt_u32(diff)
#elif defined(HAVE_BUILTIN_CTZ)
#define count_trailing_zeroes_32(diff) __builtin_ctz(diff)
#else
/* lowbit = -diff & diff */
#define count_trailing_zeroes_32(diff) mod_37_bit_position[(-diff & diff) % 37]
#endif


/* Counts matches from pos5 to pos3 up to first mismatch.  Modified from mismatches_left */
int
Genome_consecutive_matches_rightward (Genome_T ome, Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				      bool plusp, int genestrand) {
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  int mismatch_position, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *query_shifted, *ptr, *endptr;
  int relpos;
  int startcolumni, endcolumni;
  UINT4 diff_32;
#ifdef HAVE_SSE2
  __m128i diff_128;
#endif
#ifdef HAVE_AVX2
  __m256i diff_256;
#endif

  debug(
	printf("\n\n");
	printf("Genome (in consecutive_matches_rightward):\n");
	print_blocks(ref_blocks,left+pos5,left+pos3);
	printf("\n");
	);


  startblocki = (left+pos5)/128U*12;
  startcolumni = ((left+pos5) % 128) / 32;
  startblocki_32 = startblocki + startcolumni;

  endblocki = (left+pos3)/128U*12;
  endcolumni = ((left+pos3) % 128) / 32;
  endblocki_32 = endblocki + endcolumni;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % STEP_SIZE;
  query_shifted = Compress_shift(query_compress,nshift);
  debug(printf("Query shifted %d:\n",nshift));
  debug(Compress_print_blocks(query_shifted,nshift,pos5,pos3));
  query_shifted += (nshift+pos5)/STEP_SIZE*COMPRESS_BLOCKSIZE;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#else
  query_shifted += startcolumni;
#endif

  startdiscard = (left+pos5) % 32;
  enddiscard = (left+pos3) % 32;
  offset = -startdiscard + pos5;
  ptr = &(ref_blocks[startblocki_32]);
  endptr = &(ref_blocks[endblocki_32]);

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
      return (mismatch_position - pos5);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_FIRST_MISMATCH) && defined(HAVE_SSE2)
    /* Shift */
    enddiscard += (endcolumni - startcolumni)*32;
    assert(startdiscard == ((left+pos5) % 128) - startcolumni*32);
    assert(enddiscard == ((left+pos3) % 128) - startcolumni*32);

    diff_128 = (block_diff_sarray_128_shift_lo)(query_shifted,ptr,plusp,genestrand,/*query_unk_mismatch_local_p*/true,
						startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    if (nonzero_p_128(diff_128)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
      debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
      return (mismatch_position - pos5);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

#else
    /* Start block */
    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
                                     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_start_32(diff_32,startdiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
      return (mismatch_position - pos5);
    }
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				       plusp,genestrand,/*query_unk_mismatch_local_p*/true);

      if (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
	return (mismatch_position - pos5);
      }
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
                                     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_end_32(diff_32,enddiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
      return (mismatch_position - pos5);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }
#endif

#if defined(USE_WRAP_FIRST_MISMATCH) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    assert(startdiscard == ((left+pos5) % 128) - startcolumni*32);
    assert(enddiscard == ((left+pos3) % 128) + (4 - startcolumni)*32);

    diff_128 = (block_diff_sarray_128_wrap_lo)(query_shifted,ptr,plusp,genestrand,/*query_unk_mismatch_local_p*/true,
					       startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    if (nonzero_p_128(diff_128)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
      debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
      return (mismatch_position - pos5);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
                                     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_start_32(diff_32,startdiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
      return (mismatch_position - pos5);
    }
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				       plusp,genestrand,/*query_unk_mismatch_local_p*/true);
      if (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
	return (mismatch_position - pos5);
      }
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ptr + 24 <= endptr) {
      diff_256 = (block_diff_sarray_256)(query_shifted,ptr,plusp,genestrand,/*query_unk_mismatch_local_p*/true);

      if (nonzero_p_256(diff_256)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_256(diff_256));
	debug(printf("returning %d - %d consecutive matches\n",mismatch_position,pos5));
	return (mismatch_position - pos5);
      }
      query_shifted += 24; ptr += 24;
      offset += 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr + 12 <= endptr) {
      diff_128 = (block_diff_sarray_128)(query_shifted,ptr,plusp,genestrand,/*query_unk_mismatch_local_p*/true);

      if (nonzero_p_128(diff_128)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
	debug(printf("returning %d - %d consecutive matches\n",mismatch_position,pos5));
	return (mismatch_position - pos5);
      }
      query_shifted += 12; ptr += 12;
      offset += 128;
    }
#else
    while (ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
					 plusp,genestrand,/*query_unk_mismatch_local_p*/true);
	if (nonzero_p_32(diff_32)) {
	  mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	  debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
	  return (mismatch_position - pos5);
	}
	query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
	offset += 32;
      }
      /* query_shifted += QUERY_NEXTROW; */ ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ptr < endptr) {
      diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				       plusp,genestrand,/*query_unk_mismatch_local_p*/true);

      if (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
	return (mismatch_position - pos5);
      }
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_end_32(diff_32,enddiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      debug(printf("Would return %d - %d consecutive matches\n",mismatch_position,pos5));
      return (mismatch_position - pos5);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }
  }
}


/* Counts matches from pos3 to pos5 up to first mismatch.  Modified from mismatches_right */
int
Genome_consecutive_matches_leftward (Genome_T ome, Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				     bool plusp, int genestrand) {
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  int mismatch_position, offset, relpos, nshift;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *query_shifted, *ptr, *startptr;
#ifndef HAVE_BUILTIN_CLZ
  Genomecomp_T top;
#endif
  int startcolumni, endcolumni;
  UINT4 diff_32;
#ifdef HAVE_SSE2
  __m128i diff_128;
#endif
#ifdef HAVE_AVX2
  __m256i diff_256;
#endif

  /* static int ncalls = 0; */
  /* printf("Number of calls to leftward: %d\n",++ncalls); */

  debug(
	printf("\n\n");
	printf("Genome (in consecutive_matches_leftward):\n");
	print_blocks(ref_blocks,left+pos5,left+pos3);
	printf("\n");
	);


  startblocki = (left+pos5)/128U*12;
  startcolumni = ((left+pos5) % 128) / 32;
  startblocki_32 = startblocki + startcolumni;

  endblocki = (left+pos3)/128U*12;
  endcolumni = ((left+pos3) % 128) / 32;
  endblocki_32 = endblocki + endcolumni;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % STEP_SIZE;
  query_shifted = Compress_shift(query_compress,nshift);
  debug(printf("Query shifted %d:\n",nshift));
  debug(Compress_print_blocks(query_shifted,nshift,pos5,pos3));
  query_shifted += (nshift+pos3)/STEP_SIZE*COMPRESS_BLOCKSIZE;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#else
  query_shifted += endcolumni;
#endif

  startdiscard = (left+pos5) % 32;
  enddiscard = (left+pos3) % 32;
  offset = (pos3 - 1) - enddiscard + 32;
  ptr = &(ref_blocks[endblocki_32]);
  startptr = &(ref_blocks[startblocki_32]);

  if (startblocki_32 == endblocki_32) {
    /* Single block */
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
      return (pos3 - mismatch_position - 1);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

  } else if (startblocki == endblocki) {
#if defined(USE_SHIFT_FIRST_MISMATCH) && defined(HAVE_SSE2)
    /* Shift */
    startdiscard += 96 - (endcolumni - startcolumni)*32;
    enddiscard += 96;
    assert(startdiscard == ((left+pos5) % 128) + (3 - endcolumni)*32);
    assert(enddiscard == ((left+pos3) % 128) + (3 - endcolumni)*32);

    diff_128 = (block_diff_sarray_128_shift_hi)(query_shifted,ptr,plusp,genestrand,/*query_unk_mismatch_local_p*/true,
						endcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    if (nonzero_p_128(diff_128)) {
      mismatch_position = offset - (relpos = count_leading_zeroes_128(diff_128));
      debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
      return (pos3 - mismatch_position - 1);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

#else
    /* End block */
    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
                                     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_end_32(diff_32,enddiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
      return (pos3 - mismatch_position - 1);
    }
    query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    offset -= 32;

    /* Single row */
    while (--endcolumni > startcolumni) {
      diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				       plusp,genestrand,/*query_unk_mismatch_local_p*/true);

      if (nonzero_p_32(diff_32)) {
	mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
	return (pos3 - mismatch_position - 1);
      }
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }

    /* Start block */
    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
                                     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_start_32(diff_32,startdiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
      return (pos3 - mismatch_position - 1);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }
#endif

#if defined(USE_WRAP_FIRST_MISMATCH) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    startdiscard += (startcolumni - endcolumni - 1)*32;
    enddiscard += 96;
    assert(startdiscard == ((left+pos5) % 128) - (endcolumni + 1)*32);
    assert(enddiscard == ((left+pos3) % 128) + (3 - endcolumni)*32);

    diff_128 = (block_diff_sarray_128_wrap_hi)(query_shifted,ptr,plusp,genestrand,/*query_unk_mismatch_local_p*/true,
					       endcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    if (nonzero_p_128(diff_128)) {
      mismatch_position = offset - (relpos = count_leading_zeroes_128(diff_128));
      debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
      return (pos3 - mismatch_position - 1);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }
#endif

  } else {
    /* End block */
    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_end_32(diff_32,enddiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
      return (pos3 - mismatch_position - 1);
    }
    query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    offset -= 32;

    /* End row */
    while (--endcolumni >= 0) {
      diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				       plusp,genestrand,/*query_unk_mismatch_local_p*/true);

      if (nonzero_p_32(diff_32)) {
	mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
	return (pos3 - mismatch_position - 1);
      }
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }
#ifdef HAVE_SSE2
    query_shifted -= QUERY_NEXTROW;
#endif
    ptr -= GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ptr >= startptr + 24) {
      diff_256 = (block_diff_sarray_256)(&(query_shifted[-15]),&(ptr[-15]),plusp,genestrand,/*query_unk_mismatch_local_p*/true);

      if (nonzero_p_256(diff_256)) {
	mismatch_position = offset - (relpos = count_leading_zeroes_256(diff_256));
	debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
	return (pos3 - mismatch_position - 1);
      }
      query_shifted -= 24; ptr -= 24;
      offset -= 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr >= startptr + 12) {
      diff_128 = (block_diff_sarray_128)(&(query_shifted[-3]),&(ptr[-3]),plusp,genestrand,/*query_unk_mismatch_local_p*/true);

      if (nonzero_p_128(diff_128)) {
	mismatch_position = offset - (relpos = count_leading_zeroes_128(diff_128));
	debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
	return (pos3 - mismatch_position - 1);
      }
      query_shifted -= 12; ptr -= 12;
      offset -= 128;
    }
#else
    while (ptr >= startptr + 12) {
     for (endcolumni = 3; endcolumni >= 0; --endcolumni) {
       diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
					plusp,genestrand,/*query_unk_mismatch_local_p*/true);
       
       if (nonzero_p_32(diff_32)) {
	 mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	 debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
	 return (pos3 - mismatch_position - 1);
       }
       query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
       offset -= 32;
     }
     /* query_shifted -= QUERY_NEXTROW; */ ptr -= GENOME_NEXTROW;
    }
#endif

    /* Start row */
    while (ptr > startptr) {
      diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
				       plusp,genestrand,/*query_unk_mismatch_local_p*/true);
      if (nonzero_p_32(diff_32)) {
	mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
	return (pos3 - mismatch_position - 1);
      }
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }

    /* Start block */
    diff_32 = (block_diff_sarray_32)(query_shifted,ptr,
                                     plusp,genestrand,/*query_unk_mismatch_local_p*/true);
    diff_32 = clear_start_32(diff_32,startdiscard);

    if (nonzero_p_32(diff_32)) {
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      debug(printf("returning %d - %d - 1 consecutive matches\n",pos3,mismatch_position));
      return (pos3 - mismatch_position - 1);
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }
  }
}


/* Modified from Genome_consecutive_matches_rightward.  Intended for
   writing LCP only, so coords are UINT4, and does not need SIMD
   instructions. */

int
Genome_consecutive_matches_pair (Genome_T ome, UINT4 lefta, UINT4 leftb, UINT4 genomelength) {
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  int mismatch_position, offset, nshift, rightshift, nshift1, nshift2, nshifta, nshiftb;
  int startdiscard, enddiscard, nblocks;
  UINT4 left1, left2;
  UINT4 startblocki_1, startblocki_2, endblocki;
  int startcolumni_1, startcolumni_2, endcolumni;
  Genomecomp_T *ptr1, *ptr2, *end, *ptr1_prev;
  Genomecomp_T diff;
  Genomecomp_T shifted1[3];
  int relpos;

  nshifta = lefta % 32;		/* Not STEP_SIZE */
  nshiftb = leftb % 32;		/* Not STEP_SIZE */
  if (nshifta < nshiftb) {
    left1 = lefta;
    left2 = leftb;
    nshift1 = nshifta;
    nshift2 = nshiftb;

  } else {
    left1 = leftb;
    left2 = lefta;
    nshift1 = nshiftb;
    nshift2 = nshifta;
  }
  nshift = nshift2 - nshift1;
  rightshift = 32 - nshift;	/* Not STEP_SIZE */
  startdiscard = nshift2;

  /* Non-SSE2 code */
  startcolumni_1 = (left1 % 128) / 32;
  startblocki_1 = left1/128U*12 + startcolumni_1;

  startcolumni_2 = (left2 % 128) / 32;
  startblocki_2 = left2/128U*12 + startcolumni_2;

  endcolumni = (genomelength % 128) / 32;
  endblocki = genomelength/128U*12 + endcolumni;


  ptr1 = &(ref_blocks[startblocki_1]);
  ptr2 = &(ref_blocks[startblocki_2]);
  end = &(ref_blocks[endblocki]);

  offset = -startdiscard;
  
  debug2(printf("\n\n"));
  debug2(printf("left1 = %u, left2 = %u, startblocki_1 = %u, startblocki_2 = %u, endblocki = %u\n",
	       left1,left2,startblocki_1,startblocki_2,endblocki));
  debug2(printf("nshift1 = %d, nshift2 = %d, nshift = %d, startdiscard = %u\n",
	       nshift1,nshift2,nshift,startdiscard));

  if (ptr1 == end) {
    /* Single block */
    enddiscard = genomelength % 32; /* Not STEP_SIZE */
    if (nshift + enddiscard < 32) { /* Not STEP_SIZE */
      enddiscard = nshift + enddiscard;

      ptr1 = &(ref_blocks[startblocki_1]);
      ptr2 = &(ref_blocks[startblocki_2]);
#ifdef WORDS_BIGENDIAN
      shifted1[0] = Bigendian_convert_uint(ptr1[0]) << nshift;
      shifted1[1] = Bigendian_convert_uint(ptr1[4]) << nshift;
      shifted1[2] = Bigendian_convert_uint(ptr1[8]) << nshift;
#else
      shifted1[0] = ptr1[0] << nshift;
      shifted1[1] = ptr1[4] << nshift;
      shifted1[2] = ptr1[8] << nshift;
#endif
      debug2(Compress_print_one_block(ptr1));
      debug2(Compress_print_one_block(ptr2));
      debug2(Compress_print_one_block(shifted1));

#ifdef WORDS_BIGENDIAN
      diff = (shifted1[0] ^ Bigendian_convert_uint(ptr2[0])) | (shifted1[1] ^ Bigendian_convert_uint(ptr2[4])) | (shifted1[2] ^ Bigendian_convert_uint(ptr2[8]));
#else
      diff = (shifted1[0] ^ ptr2[0]) | (shifted1[1] ^ ptr2[4]) | (shifted1[2] ^ ptr2[8]);
#endif
      diff = clear_start_32(diff,startdiscard);
      diff = clear_end_32(diff,enddiscard);

      if (diff /* != 0U */) {
#ifdef HAVE_BUILTIN_CTZ
	mismatch_position = offset + (relpos = __builtin_ctz(diff));
#else
	mismatch_position = offset + mod_37_bit_position[(-diff & diff) % 37];
#endif
	debug2(printf("case 1: returning %d matches\n",mismatch_position));
	return mismatch_position;
      } else {
	debug2(printf("case 2: returning %d - %d matches\n",enddiscard,startdiscard));
	return (enddiscard - startdiscard);
      }

    } else {
      /* Two blocks */
      if (nshift > 0) {
	enddiscard -= (32 - nshift); /* Not STEP_SIZE */
      }

      /* Block 1 */
      ptr1 = &(ref_blocks[startblocki_1]);
      ptr2 = &(ref_blocks[startblocki_2]);
#ifdef WORDS_BIGENDIAN
      shifted1[0] = Bigendian_convert_uint(ptr1[0]) << nshift;
      shifted1[1] = Bigendian_convert_uint(ptr1[4]) << nshift;
      shifted1[2] = Bigendian_convert_uint(ptr1[8]) << nshift;
#else
      shifted1[0] = ptr1[0] << nshift;
      shifted1[1] = ptr1[4] << nshift;
      shifted1[2] = ptr1[8] << nshift;
#endif
      debug2(Compress_print_one_block(ptr1));
      debug2(Compress_print_one_block(ptr2));
      debug2(Compress_print_one_block(shifted1));

#ifdef WORDS_BIGENDIAN
      diff = (shifted1[0] ^ Bigendian_convert_uint(ptr2[0])) | (shifted1[1] ^ Bigendian_convert_uint(ptr2[4])) | (shifted1[2] ^ Bigendian_convert_uint(ptr2[8]));
#else
      diff = (shifted1[0] ^ ptr2[0]) | (shifted1[1] ^ ptr2[4]) | (shifted1[2] ^ ptr2[8]);
#endif
      diff = clear_start_32(diff,startdiscard);

      if (diff /* != 0U */) {
#ifdef HAVE_BUILTIN_CTZ
	mismatch_position = offset + (relpos = __builtin_ctz(diff));
#else
	mismatch_position = offset + mod_37_bit_position[(-diff & diff) % 37];
#endif
	debug2(printf("case 3: returning %d matches\n",mismatch_position));
	return mismatch_position;
      } else {
	ptr1_prev = ptr1;
	ptr1 += 1; if (++startcolumni_1 == 4) {ptr1 += 8; startcolumni_1 = 0;}
	ptr2 += 1; if (++startcolumni_2 == 4) {ptr2 += 8; startcolumni_2 = 0;}
	offset += 32; /* Not STEP_SIZE */
      }

      /* Block 2 */
      if (nshift == 0) {
	/* rightshift of 32 is a no-op */
#ifdef WORDS_BIGENDIAN
	shifted1[0] = Bigendian_convert_uint(ptr1[0]); shifted1[1] = Bigendian_convert_uint(ptr1[4]); shifted1[2] = Bigendian_convert_uint(ptr1[8]);
#else
	shifted1[0] = ptr1[0]; shifted1[1] = ptr1[4]; shifted1[2] = ptr1[8];
#endif
      } else {
#ifdef WORDS_BIGENDIAN
	shifted1[0] = (Bigendian_convert_uint(ptr1[0]) << nshift) | (Bigendian_convert_uint(ptr1_prev[0]) >> rightshift);
	shifted1[1] = (Bigendian_convert_uint(ptr1[4]) << nshift) | (Bigendian_convert_uint(ptr1_prev[4]) >> rightshift);
	shifted1[2] = (Bigendian_convert_uint(ptr1[8]) << nshift) | (Bigendian_convert_uint(ptr1_prev[8]) >> rightshift);
#else
	shifted1[0] = (ptr1[0] << nshift) | (ptr1_prev[0] >> rightshift);
	shifted1[1] = (ptr1[4] << nshift) | (ptr1_prev[4] >> rightshift);
	shifted1[2] = (ptr1[8] << nshift) | (ptr1_prev[8] >> rightshift);
#endif
      }
      debug2(Compress_print_one_block(ptr1));
      debug2(Compress_print_one_block(ptr2));
      debug2(Compress_print_one_block(shifted1));

#ifdef WORDS_BIGENDIAN
      diff = (shifted1[0] ^ Bigendian_convert_uint(ptr2[0])) | (shifted1[1] ^ Bigendian_convert_uint(ptr2[4])) | (shifted1[2] ^ Bigendian_convert_uint(ptr2[8]));
#else
      diff = (shifted1[0] ^ ptr2[0]) | (shifted1[1] ^ ptr2[4]) | (shifted1[2] ^ ptr2[8]);
#endif
      diff = clear_end_32(diff,enddiscard);

      if (diff /* != 0U */) {
#ifdef HAVE_BUILTIN_CTZ
	mismatch_position = offset + (relpos = __builtin_ctz(diff));
#else
	mismatch_position = offset + mod_37_bit_position[(-diff & diff) % 37];
#endif
	debug2(printf("case 4: returning %d matches\n",mismatch_position));
	return mismatch_position;
      } else {
	debug2(printf("case 5: returning offset %d + enddiscard %d matches\n",offset,enddiscard));
	return offset + enddiscard;
      }
    }

  } else if (ptr2 == end) {
    /* Single block */
    enddiscard = genomelength % 32; /* Not STEP_SIZE */

    ptr1 = &(ref_blocks[startblocki_1]);
    ptr2 = &(ref_blocks[startblocki_2]);
#ifdef WORDS_BIGENDIAN
    shifted1[0] = Bigendian_convert_uint(ptr1[0]) << nshift;
    shifted1[1] = Bigendian_convert_uint(ptr1[4]) << nshift;
    shifted1[2] = Bigendian_convert_uint(ptr1[8]) << nshift;
#else
    shifted1[0] = ptr1[0] << nshift;
    shifted1[1] = ptr1[4] << nshift;
    shifted1[2] = ptr1[8] << nshift;
#endif
    debug2(Compress_print_one_block(ptr1));
    debug2(Compress_print_one_block(ptr2));
    debug2(Compress_print_one_block(shifted1));

#ifdef WORDS_BIGENDIAN
    diff = (shifted1[0] ^ Bigendian_convert_uint(ptr2[0])) | (shifted1[1] ^ Bigendian_convert_uint(ptr2[4])) | (shifted1[2] ^ Bigendian_convert_uint(ptr2[8]));
#else
    diff = (shifted1[0] ^ ptr2[0]) | (shifted1[1] ^ ptr2[4]) | (shifted1[2] ^ ptr2[8]);
#endif
    diff = clear_start_32(diff,startdiscard);
    diff = clear_end_32(diff,enddiscard);

    if (diff /* != 0U */) {
#ifdef HAVE_BUILTIN_CTZ
      mismatch_position = offset + (relpos = __builtin_ctz(diff));
#else
      mismatch_position = offset + mod_37_bit_position[(-diff & diff) % 37];
#endif
      debug2(printf("case 6: returning %d matches\n",mismatch_position));
      return mismatch_position;
    } else {
      debug2(printf("case 7: returning %d - %d matches\n",enddiscard,startdiscard));
      return (enddiscard - startdiscard);
    }

  } else {

    /* Startblock */
    ptr1 = &(ref_blocks[startblocki_1]);
    ptr2 = &(ref_blocks[startblocki_2]);
#ifdef WORDS_BIGENDIAN
    shifted1[0] = (Bigendian_convert_uint(ptr1[0]) << nshift);
    shifted1[1] = (Bigendian_convert_uint(ptr1[4]) << nshift);
    shifted1[2] = (Bigendian_convert_uint(ptr1[8]) << nshift);
#else
    shifted1[0] = (ptr1[0] << nshift);
    shifted1[1] = (ptr1[4] << nshift);
    shifted1[2] = (ptr1[8] << nshift);
#endif
    debug2(Compress_print_one_block(ptr1));
    debug2(Compress_print_one_block(ptr2));
    debug2(Compress_print_one_block(shifted1));

#ifdef WORDS_BIGENDIAN
    diff = (shifted1[0] ^ Bigendian_convert_uint(ptr2[0])) | (shifted1[1] ^ Bigendian_convert_uint(ptr2[4])) | (shifted1[2] ^ Bigendian_convert_uint(ptr2[8]));
#else
    diff = (shifted1[0] ^ ptr2[0]) | (shifted1[1] ^ ptr2[4]) | (shifted1[2] ^ ptr2[8]);
#endif
    diff = clear_start_32(diff,startdiscard);

    if (diff /* != 0U */) {
#ifdef HAVE_BUILTIN_CTZ
      mismatch_position = offset + (relpos = __builtin_ctz(diff));
#else
      mismatch_position = offset + mod_37_bit_position[(-diff & diff) % 37];
#endif
      debug2(printf("case 8: returning %d matches\n",mismatch_position));
      return mismatch_position;
    } else {
      ptr1_prev = ptr1;
      ptr1 += 1; if (++startcolumni_1 == 4) {ptr1 += 8; startcolumni_1 = 0;}
      ptr2 += 1; if (++startcolumni_2 == 4) {ptr2 += 8; startcolumni_2 = 0;}
      offset += 32;		/* Not STEP_SIZE */
    }

    while (ptr1 < end && ptr2 < end) {
      if (nshift == 0) {
	/* rightshift of 32 is a no-op */
#ifdef WORDS_BIGENDIAN
	shifted1[0] = Bigendian_convert_uint(ptr1[0]); shifted1[1] = Bigendian_convert_uint(ptr1[4]); shifted1[2] = Bigendian_convert_uint(ptr1[8]);
#else
	shifted1[0] = ptr1[0]; shifted1[1] = ptr1[4]; shifted1[2] = ptr1[8];
#endif
      } else {
#ifdef WORDS_BIGENDIAN
	shifted1[0] = (Bigendian_convert_uint(ptr1[0]) << nshift) | (Bigendian_convert_uint(ptr1_prev[0]) >> rightshift);
	shifted1[1] = (Bigendian_convert_uint(ptr1[4]) << nshift) | (Bigendian_convert_uint(ptr1_prev[4]) >> rightshift);
	shifted1[2] = (Bigendian_convert_uint(ptr1[8]) << nshift) | (Bigendian_convert_uint(ptr1_prev[8]) >> rightshift);
#else
	shifted1[0] = (ptr1[0] << nshift) | (ptr1_prev[0] >> rightshift);
	shifted1[1] = (ptr1[4] << nshift) | (ptr1_prev[4] >> rightshift);
	shifted1[2] = (ptr1[8] << nshift) | (ptr1_prev[8] >> rightshift);
#endif
      }
      debug2(Compress_print_one_block(ptr1));
      debug2(Compress_print_one_block(ptr2));
      debug2(Compress_print_one_block(shifted1));

#ifdef WORDS_BIGENDIAN
      diff = (shifted1[0] ^ Bigendian_convert_uint(ptr2[0])) | (shifted1[1] ^ Bigendian_convert_uint(ptr2[4])) | (shifted1[2] ^ Bigendian_convert_uint(ptr2[8]));
#else
      diff = (shifted1[0] ^ ptr2[0]) | (shifted1[1] ^ ptr2[4]) | (shifted1[2] ^ ptr2[8]);
#endif
      if (diff /* != 0U */) {
#ifdef HAVE_BUILTIN_CTZ
	mismatch_position = offset + (relpos = __builtin_ctz(diff));
#else
	mismatch_position = offset + mod_37_bit_position[(-diff & diff) % 37];
#endif
	debug2(printf("case 9: returning %d matches\n",mismatch_position));
	return mismatch_position;
      } else {
	ptr1_prev = ptr1;
	ptr1 += 1; if (++startcolumni_1 == 4) {ptr1 += 8; startcolumni_1 = 0;}
	ptr2 += 1; if (++startcolumni_2 == 4) {ptr2 += 8; startcolumni_2 = 0;}
	offset += 32;		/* Not STEP_SIZE */
      }
    }

    /* Last block of entire genome */
    enddiscard = genomelength % 32; /* Not STEP_SIZE */
    if (ptr2 == end) {
      debug2(printf("ptr2 == end\n"));
      /* Keep enddiscard */
      nblocks = 1;
    } else if (nshift + enddiscard < 32) {
      debug2(printf("ptr1 == end and nshift %d + enddiscard %d < 32\n",nshift,enddiscard));
      enddiscard = nshift + enddiscard;
      nblocks = 1;
    } else if (nshift > 0) {
      debug2(printf("ptr1 == end and nshift %d + enddiscard %d >= 32\n",nshift,enddiscard));
      enddiscard -= (32 - nshift);
      nblocks = 2;
    } else {
      debug2(printf("ptr1 == end and nshift %d + enddiscard %d >= 32\n",nshift,enddiscard));
      /* Keep enddiscard */
      nblocks = 2;
    }

    /* Block 1 */
    if (nshift == 0) {
      /* rightshift of 32 is a no-op */
#ifdef WORDS_BIGENDIAN
      shifted1[0] = Bigendian_convert_uint(ptr1[0]); shifted1[1] = Bigendian_convert_uint(ptr1[4]); shifted1[2] = Bigendian_convert_uint(ptr1[8]);
#else
      shifted1[0] = ptr1[0]; shifted1[1] = ptr1[4]; shifted1[2] = ptr1[8];
#endif
    } else {
#ifdef WORDS_BIGENDIAN
      shifted1[0] = (Bigendian_convert_uint(ptr1[0]) << nshift) | (Bigendian_convert_uint(ptr1_prev[0]) >> rightshift);
      shifted1[1] = (Bigendian_convert_uint(ptr1[4]) << nshift) | (Bigendian_convert_uint(ptr1_prev[4]) >> rightshift);
      shifted1[2] = (Bigendian_convert_uint(ptr1[8]) << nshift) | (Bigendian_convert_uint(ptr1_prev[8]) >> rightshift);
#else
      shifted1[0] = (ptr1[0] << nshift) | (ptr1_prev[0] >> rightshift);
      shifted1[1] = (ptr1[4] << nshift) | (ptr1_prev[4] >> rightshift);
      shifted1[2] = (ptr1[8] << nshift) | (ptr1_prev[8] >> rightshift);
#endif
    }
    debug2(Compress_print_one_block(ptr1));
    debug2(Compress_print_one_block(ptr2));
    debug2(Compress_print_one_block(shifted1));

#ifdef WORDS_BIGENDIAN
    diff = (shifted1[0] ^ Bigendian_convert_uint(ptr2[0])) | (shifted1[1] ^ Bigendian_convert_uint(ptr2[4])) | (shifted1[2] ^ Bigendian_convert_uint(ptr2[8]));
#else
    diff = (shifted1[0] ^ ptr2[0]) | (shifted1[1] ^ ptr2[4]) | (shifted1[2] ^ ptr2[8]);
#endif
    if (nblocks == 1) {
      diff = clear_end_32(diff,enddiscard);
    }

    if (diff /* != 0U */) {
#ifdef HAVE_BUILTIN_CTZ
      mismatch_position = offset + (relpos = __builtin_ctz(diff));
#else
      mismatch_position = offset + mod_37_bit_position[(-diff & diff) % 37];
#endif
      debug2(printf("case 10: returning %d matches\n",mismatch_position));
      return mismatch_position;
    } else if (nblocks == 1) {
      debug2(printf("case 11: returning offset %d + enddiscard %d matches\n",offset,enddiscard));
      return offset + enddiscard;
    } else {
      ptr1_prev = ptr1;
      ptr1 += 1; if (++startcolumni_1 == 4) {ptr1 += 8; startcolumni_1 = 0;}
      ptr2 += 1; if (++startcolumni_2 == 4) {ptr2 += 8; startcolumni_2 = 0;}
      offset += 32;		/* Not STEP_SIZE */
    }

    /* Block 2 */
#ifdef WORDS_BIGENDIAN
    shifted1[0] = (Bigendian_convert_uint(ptr1_prev[0]) >> rightshift);
    shifted1[1] = (Bigendian_convert_uint(ptr1_prev[4]) >> rightshift);
    shifted1[2] = (Bigendian_convert_uint(ptr1_prev[8]) >> rightshift);
#else
    shifted1[0] = (ptr1_prev[0] >> rightshift);
    shifted1[1] = (ptr1_prev[4] >> rightshift);
    shifted1[2] = (ptr1_prev[8] >> rightshift);
#endif
    debug2(Compress_print_one_block(ptr1));
    debug2(Compress_print_one_block(ptr2));
    debug2(Compress_print_one_block(shifted1));

#ifdef WORDS_BIGENDIAN
    diff = (shifted1[0] ^ Bigendian_convert_uint(ptr2[0])) | (shifted1[1] ^ Bigendian_convert_uint(ptr2[4])) | (shifted1[2] ^ Bigendian_convert_uint(ptr2[8]));
#else
    diff = (shifted1[0] ^ ptr2[0]) | (shifted1[1] ^ ptr2[4]) | (shifted1[2] ^ ptr2[8]);
#endif
    diff = clear_end_32(diff,enddiscard);

    if (diff /* != 0U */) {
#ifdef HAVE_BUILTIN_CTZ
      mismatch_position = offset + (relpos = __builtin_ctz(diff));
#else
      mismatch_position = offset + mod_37_bit_position[(-diff & diff) % 37];
#endif
      debug2(printf("case 12: returning %d matches\n",mismatch_position));
      return mismatch_position;
    } else {
      debug2(printf("case 13: returning offset %d + enddiscard %d matches\n",offset,enddiscard));
      return offset + enddiscard;
    }
  }
}

