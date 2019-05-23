static char rcsid[] = "$Id: genome128_hr.c 218673 2019-03-16 01:23:08Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "genome128_hr.h"

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
      printf("  startdiscard:%d\n",startdiscard32);
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
      printf("  enddiscard:%d\n",enddiscard32);
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
      printf("  startdiscard:%d\n",startdiscard32);
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
      printf("  enddiscard:%d\n",enddiscard32);
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
      printf("  startdiscard:%d\n",startdiscard32);
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
      printf("  enddiscard:%d\n",enddiscard32);
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
      printf("  startdiscard:%d\n",startdiscard32);
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
      printf("  enddiscard:%d\n",enddiscard32);
    }

    printf("\n");
    ptr += 12;
#endif
  }

  return;
}
#endif

#if defined(DEBUG) || defined(DEBUG5)
static void
print_blocks_snp (Genomecomp_T *blocks, Genomecomp_T *snp_blocks, Univcoord_T startpos, Univcoord_T endpos) {
  /* Chrpos_T length = endpos - startpos; */
  Genomecomp_T *ref_ptr, *snp_ptr;
  Univcoord_T startblocki, endblocki;
  int startcolumni, endcolumni;
  int startdiscard32, enddiscard32;
  Genomecomp_T high, low, flags, snpmask;

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  startcolumni = (startpos % 128) / 32;
  startblocki = startpos/128U*12 + startcolumni;

  endcolumni = (endpos % 128) / 32;
  endblocki = endpos/128U*12 + endcolumni;
#else
  startblocki = startpos/128U*12;
  endblocki = endpos/128U*12;
#endif

  startdiscard32 = startpos % 32;
  enddiscard32 = endpos % 32;
  
  ref_ptr = &(blocks[startblocki]);
  snp_ptr = &(snp_blocks[startblocki]);
  while (ref_ptr <= &(blocks[endblocki])) {
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
    high = ref_ptr[0]; low = ref_ptr[4]; flags = ref_ptr[8]; snpmask = snp_ptr[8];
    printf("high: %08X  low: %08X  flags: %08X  snpmask: %08X\n",high,low,flags,snpmask);
    ref_ptr += GENOME_NEXTCOL; snp_ptr += GENOME_NEXTCOL; if (++startcolumni == 4) {ref_ptr += GENOME_NEXTROW; snp_ptr += GENOME_NEXTROW; startcolumni = 0;}
#else
    high = ref_ptr[0]; low = ref_ptr[4]; flags = ref_ptr[8]; snpmask = snp_ptr[8];
    printf("high: %08X  low: %08X  flags: %08X  snpmask: %08X\n",high,low,flags,snpmask);

    high = ref_ptr[1]; low = ref_ptr[5]; flags = ref_ptr[9]; snpmask = snp_ptr[9];
    printf("high: %08X  low: %08X  flags: %08X  snpmask: %08X\n",high,low,flags,snpmask);

    high = ref_ptr[2]; low = ref_ptr[6]; flags = ref_ptr[10]; snpmask = snp_ptr[10];
    printf("high: %08X  low: %08X  flags: %08X  snpmask: %08X\n",high,low,flags,snpmask);

    high = ref_ptr[3]; low = ref_ptr[7]; flags = ref_ptr[11]; snpmask = snp_ptr[11];
    printf("high: %08X  low: %08X  flags: %08X  snpmask: %08X\n",high,low,flags,snpmask);

    ref_ptr += 12; snp_ptr += 12;
#endif
  }

  return;
}
#endif


static Genomecomp_T *ref_blocks;
static Genomecomp_T *snp_blocks;

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

#if 0
  printf("query high: ");
  print_vector_hex(_query_high);
  printf("query low:  ");
  print_vector_hex(_query_low);

  printf("  ref high: ");
  print_vector_hex(_ref_high);
  printf("  ref low:  ");
  print_vector_hex(_ref_low);
#endif

  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

#if 0
  printf("      diff: ");
  print_vector_hex(_diff);
  printf("\n");
#endif

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


static UINT4
block_diff_standard_wildcard_32 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
				   bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  UINT4 diff, non_wildcard;

  /* Taken from block_diff_standard */
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

  /* Add difference relative to SNP */
#ifdef WORDS_BIGENDIAN
  diff &= (query_shifted[0] ^ Bigendian_convert_uint(snp_ptr[0])) | (query_shifted[1] ^ Bigendian_convert_uint(snp_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff &= (query_shifted[0] ^ snp_ptr[0]) | (query_shifted[1] ^ snp_ptr[4]);
#else
  diff &= (query_shifted[0] ^ snp_ptr[0]) | (query_shifted[4] ^ snp_ptr[4]);
#endif

  /* Test for equality of ref and alt */
  debug(printf("Equality high: ref genome %08X with alt genome %08X ",ref_ptr[0],snp_ptr[0]));
#ifdef WORDS_BIGENDIAN
  non_wildcard = (Bigendian_convert_uint(ref_ptr[0]) ^ Bigendian_convert_uint(snp_ptr[0])) |
    (Bigendian_convert_uint(ref_ptr[4]) ^ Bigendian_convert_uint(snp_ptr[4]));
#else
  non_wildcard = (ref_ptr[0] ^ snp_ptr[0]) | (ref_ptr[4] ^ snp_ptr[4]);
#endif
  debug(printf(" => diff %08X\n",non_wildcard));
  
  /* Ref flags */
  debug(printf("Wildcard add ref flags: ref genome %08X and alt genome %08X ",ref_ptr[8],snp_ptr[8]));
#ifdef WORDS_BIGENDIAN
  non_wildcard |= Bigendian_convert_uint(ref_ptr[8]);
#else
  non_wildcard |= ref_ptr[8];
#endif

  /* Alt flags */
  debug(printf("Wildcard add alt flags: ref genome %08X and alt genome %08X ",ref_ptr[8],snp_ptr[8]));
#ifdef WORDS_BIGENDIAN
  non_wildcard |= ~(Bigendian_convert_uint(snp_ptr[8]));
#else
  non_wildcard |= ~(snp_ptr[8]);
#endif
  debug(printf(" => non_wildcard %08X\n",non_wildcard));

  return diff & non_wildcard;
}


/* wildcard if ref == alt && ref_flag == 0 && alt_flag == 1 */
/* not wildcard if ref != alt || ref_flag == 1 || alt_flag == 0 */
/* diffs are (query ^ ref) & (query ^ alt) & ~wildcard */
/* snp_ptr here is alt_ptr */
#ifdef HAVE_SSE2
static __m128i
block_diff_standard_wildcard_128 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
				  bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  __m128i _diff, _wildcard, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _snp_high, _snp_low, _snp_flags;

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
  /* End of (query ^ ref) */


  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
  _snp_high = _mm_load_si128((__m128i *) snp_ptr);
  _snp_low = _mm_load_si128((__m128i *) &(snp_ptr[4]));

  _diff = _mm_and_si128(_diff, _mm_or_si128(_mm_xor_si128(_query_high, _snp_high), _mm_xor_si128(_query_low, _snp_low)));


  /* Test for equality of ref and alt */
  _snp_flags = _mm_load_si128((__m128i *) &(snp_ptr[8]));
  _wildcard = _mm_andnot_si128(_ref_flags, _snp_flags);
  _wildcard = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(_ref_high, _snp_high), _mm_xor_si128(_ref_low, _snp_low)), _wildcard);

  _diff = _mm_andnot_si128(_wildcard, _diff);

  return _diff;
}

#ifdef USE_SHIFT_HILO
static __m128i
block_diff_standard_wildcard_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
					   bool plusp, int genestrand, bool query_unk_mismatch_local_p,
					   int startcolumni) {
  __m128i _diff, _wildcard, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _snp_high, _snp_low, _snp_flags;

  read_128_shift_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_shift_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);
  read_128_shift_lo(&_snp_high,&_snp_low,&_snp_flags,snp_ptr,startcolumni);

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
  /* End of (query ^ ref) */

  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
  _diff = _mm_and_si128(_diff, _mm_or_si128(_mm_xor_si128(_query_high, _snp_high), _mm_xor_si128(_query_low, _snp_low)));

  /* Test for equality of ref and alt */
  _wildcard = _mm_andnot_si128(_ref_flags, _snp_flags);
  _wildcard = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(_ref_high, _snp_high), _mm_xor_si128(_ref_low, _snp_low)), _wildcard);

  _diff = _mm_andnot_si128(_wildcard, _diff);

  return _diff;
}

static __m128i
block_diff_standard_wildcard_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
					   bool plusp, int genestrand, bool query_unk_mismatch_local_p,
					   int endcolumni) {
  __m128i _diff, _wildcard, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _snp_high, _snp_low, _snp_flags;

  read_128_shift_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_shift_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);
  read_128_shift_hi(&_snp_high,&_snp_low,&_snp_flags,snp_ptr,endcolumni);

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
  /* End of (query ^ ref) */

  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
  _diff = _mm_and_si128(_diff, _mm_or_si128(_mm_xor_si128(_query_high, _snp_high), _mm_xor_si128(_query_low, _snp_low)));

  /* Test for equality of ref and alt */
  _wildcard = _mm_andnot_si128(_ref_flags, _snp_flags);
  _wildcard = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(_ref_high, _snp_high), _mm_xor_si128(_ref_low, _snp_low)), _wildcard);

  _diff = _mm_andnot_si128(_wildcard, _diff);

  return _diff;
}
#endif
#endif

#if defined(GSNAP) && defined(HAVE_SSSE3)
static __m128i
block_diff_standard_wildcard_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
					  bool plusp, int genestrand, bool query_unk_mismatch_local_p, int startcolumni) {
  __m128i _diff, _wildcard, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _snp_high, _snp_low, _snp_flags;

  read_128_wrap_lo(&_query_high,&_query_low,&_query_flags,query_shifted,startcolumni);
  read_128_wrap_lo(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,startcolumni);
  read_128_wrap_lo(&_snp_high,&_snp_low,&_snp_flags,snp_ptr,startcolumni);

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
  /* End of (query ^ ref) */


  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
  _diff = _mm_and_si128(_diff, _mm_or_si128(_mm_xor_si128(_query_high, _snp_high), _mm_xor_si128(_query_low, _snp_low)));


  /* Test for equality of ref and alt */
  _wildcard = _mm_andnot_si128(_ref_flags, _snp_flags);
  _wildcard = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(_ref_high, _snp_high), _mm_xor_si128(_ref_low, _snp_low)), _wildcard);

  _diff = _mm_andnot_si128(_wildcard, _diff);

  return _diff;
}
#endif

#if defined(GSNAP) && defined(HAVE_SSSE3)
static __m128i
block_diff_standard_wildcard_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
					  bool plusp, int genestrand, bool query_unk_mismatch_local_p, int endcolumni) {
  __m128i _diff, _wildcard, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _snp_high, _snp_low, _snp_flags;

  read_128_wrap_hi(&_query_high,&_query_low,&_query_flags,query_shifted,endcolumni);
  read_128_wrap_hi(&_ref_high,&_ref_low,&_ref_flags,ref_ptr,endcolumni);
  read_128_wrap_hi(&_snp_high,&_snp_low,&_snp_flags,snp_ptr,endcolumni);

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
  /* End of (query ^ ref) */


  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
  _diff = _mm_and_si128(_diff, _mm_or_si128(_mm_xor_si128(_query_high, _snp_high), _mm_xor_si128(_query_low, _snp_low)));


  /* Test for equality of ref and alt */
  _wildcard = _mm_andnot_si128(_ref_flags, _snp_flags);
  _wildcard = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(_ref_high, _snp_high), _mm_xor_si128(_ref_low, _snp_low)), _wildcard);

  _diff = _mm_andnot_si128(_wildcard, _diff);

  return _diff;
}
#endif

/* wildcard if ref == alt && ref_flag == 0 && alt_flag == 1 */
/* not wildcard if ref != alt || ref_flag == 1 || alt_flag == 0 */
/* diffs are (query ^ ref) & (query ^ alt) & ~wildcard */
/* snp_ptr here is alt_ptr */
#ifdef HAVE_AVX2
static __m256i
block_diff_standard_wildcard_256 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
				  bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  __m256i _diff, _wildcard, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _snp_high, _snp_low, _snp_flags;

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
  /* End of (query ^ ref) */

  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
  read_256(&_snp_high,&_snp_low,&_snp_flags,snp_ptr);

  _diff = _mm256_and_si256(_diff, _mm256_or_si256(_mm256_xor_si256(_query_high, _snp_high), _mm256_xor_si256(_query_low, _snp_low)));

  /* Test for equality of ref and alt */
  _wildcard = _mm256_andnot_si256(_ref_flags, _snp_flags);
  _wildcard = _mm256_andnot_si256(_mm256_or_si256(_mm256_xor_si256(_ref_high, _snp_high), _mm256_xor_si256(_ref_low, _snp_low)), _wildcard);

  _diff = _mm256_andnot_si256(_wildcard, _diff);

  return _diff;
}
#endif

/* wildcard if ref == alt && ref_flag == 0 && alt_flag == 1 */
/* not wildcard if ref != alt || ref_flag == 1 || alt_flag == 0 */
/* diffs are (query ^ ref) & (query ^ alt) & ~wildcard */
/* snp_ptr here is alt_ptr */
#ifdef HAVE_AVX512
static __m512i
block_diff_standard_wildcard_512 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
				  bool plusp, int genestrand, bool query_unk_mismatch_local_p) {
  __m512i _diff, _wildcard, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _snp_high, _snp_low, _snp_flags;

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
  /* End of (query ^ ref) */

  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
  read_512(&_snp_high,&_snp_low,&_snp_flags,snp_ptr);

  _diff = _mm512_and_si512(_diff, _mm512_or_si512(_mm512_xor_si512(_query_high, _snp_high), _mm512_xor_si512(_query_low, _snp_low)));

  /* Test for equality of ref and alt */
  _wildcard = _mm512_andnot_si512(_ref_flags, _snp_flags);
  _wildcard = _mm512_andnot_si512(_mm512_or_si512(_mm512_xor_si512(_ref_high, _snp_high), _mm512_xor_si512(_ref_low, _snp_low)), _wildcard);

  _diff = _mm512_andnot_si512(_wildcard, _diff);

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

  /* sarrayp == false */
  /* Mark genome-T to query-C mismatches */
#ifdef WORDS_BIGENDIAN
  diff = (~(query_shifted[0]) & query_shifted[1]) &
    (Bigendian_convert_uint(ref_ptr[0]) & Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff = (~(query_shifted[0]) & query_shifted[1]) & (ref_ptr[0] & ref_ptr[4]);
#else
  diff = (~(query_shifted[0]) & query_shifted[4]) & (ref_ptr[0] & ref_ptr[4]);
#endif
  debug(printf(" => diff %08X\n",diff));

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

  /* sarrayp == false */
  /* Mark genome-T to query-C mismatches */
  _diff = _mm_and_si128(_mm_andnot_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-T to query-C mismatches */
  _diff = _mm_and_si128(_mm_andnot_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-T to query-C mismatches */
  _diff = _mm_and_si128(_mm_andnot_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-T to query-C mismatches */
  _diff = _mm_and_si128(_mm_andnot_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-T to query-C mismatches */
  _diff = _mm_and_si128(_mm_andnot_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-T to query-C mismatches */
  _diff = _mm256_and_si256(_mm256_andnot_si256(_query_high, _query_low), _mm256_and_si256(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-T to query-C mismatches */
  _diff = _mm512_and_si512(_mm512_andnot_si512(_query_high, _query_low), _mm512_and_si512(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-A to query-G mismatches */
#ifdef WORDS_BIGENDIAN
  diff = (query_shifted[0] & ~(query_shifted[1])) &
    ~(Bigendian_convert_uint(ref_ptr[0]) | Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff = (query_shifted[0] & ~(query_shifted[1])) & ~(ref_ptr[0] | ref_ptr[4]);
#else
  diff = (query_shifted[0] & ~(query_shifted[4])) & ~(ref_ptr[0] | ref_ptr[4]);
#endif
  debug(printf(" => diff %08X\n",diff));

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

  /* sarrayp == false */
  /* Mark genome-A to query-G mismatches */
  _diff = _mm_andnot_si128(_query_low, _query_high);
  _diff = _mm_andnot_si128(_ref_high, _diff);
  _diff = _mm_andnot_si128(_ref_low, _diff);

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

  /* sarrayp == false */
  /* Mark genome-A to query-G mismatches */
  _diff = _mm_andnot_si128(_query_low, _query_high);
  _diff = _mm_andnot_si128(_ref_high, _diff);
  _diff = _mm_andnot_si128(_ref_low, _diff);

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

  /* sarrayp == false */
  /* Mark genome-A to query-G mismatches */
  _diff = _mm_andnot_si128(_query_low, _query_high);
  _diff = _mm_andnot_si128(_ref_high, _diff);
  _diff = _mm_andnot_si128(_ref_low, _diff);

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

  /* sarrayp == false */
  /* Mark genome-A to query-G mismatches */
  _diff = _mm_andnot_si128(_query_low, _query_high);
  _diff = _mm_andnot_si128(_ref_high, _diff);
  _diff = _mm_andnot_si128(_ref_low, _diff);

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

  /* sarrayp == false */
  /* Mark genome-A to query-G mismatches */
  _diff = _mm_andnot_si128(_query_low, _query_high);
  _diff = _mm_andnot_si128(_ref_high, _diff);
  _diff = _mm_andnot_si128(_ref_low, _diff);

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

  /* sarrayp == false */
  /* Mark genome-A to query-G mismatches */
  _diff = _mm256_andnot_si256(_query_low, _query_high);
  _diff = _mm256_andnot_si256(_ref_high, _diff);
  _diff = _mm256_andnot_si256(_ref_low, _diff);

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

  /* sarrayp == false */
  /* Mark genome-A to query-G mismatches */
  _diff = _mm512_andnot_si512(_query_low, _query_high);
  _diff = _mm512_andnot_si512(_ref_high, _diff);
  _diff = _mm512_andnot_si512(_ref_low, _diff);

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
block_diff_cmet_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_cmet_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_cmet_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_cmet_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_cmet_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_cmet_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_cmet_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_cmet_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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


#ifdef GSNAP
/* Ignores snp_ptr */
static UINT4
block_diff_cmet_snp_32 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
#endif


#if defined(GSNAP) && defined(HAVE_SSE2)
/* Ignores snp_ptr */
static __m128i
block_diff_cmet_snp_128 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_cmet_snp_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_cmet_snp_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

#if defined(GSNAP) && defined(HAVE_SSSE3)
/* Ignores snp_ptr */
static __m128i
block_diff_cmet_snp_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_cmet_snp_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

#if defined(GSNAP) && defined(HAVE_AVX2)
/* Ignores snp_ptr */
static __m256i
block_diff_cmet_snp_256 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

#if defined(GSNAP) && defined(HAVE_AVX512)
/* Ignores snp_ptr */
static __m512i
block_diff_cmet_snp_512 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

  /* sarrayp == false */
  /* Mark genome-G to query-A mismatches */
#ifdef WORDS_BIGENDIAN
  diff = ~(query_shifted[0] | query_shifted[1]) &
    (Bigendian_convert_uint(ref_ptr[0]) & ~Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff = ~(query_shifted[0] | query_shifted[1]) & (ref_ptr[0] & ~(ref_ptr[4]));
#else
  diff = ~(query_shifted[0] | query_shifted[4]) & (ref_ptr[0] & ~(ref_ptr[4]));
#endif
  debug(printf(" => diff %08X\n",diff));

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

  /* sarrayp == false */
  /* Mark genome-G to query-A mismatches */
  _diff = _mm_andnot_si128(_mm_or_si128(_query_high, _query_low), _mm_andnot_si128(_ref_low, _ref_high));

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

  /* sarrayp == false */
  /* Mark genome-G to query-A mismatches */
  _diff = _mm_andnot_si128(_mm_or_si128(_query_high, _query_low), _mm_andnot_si128(_ref_low, _ref_high));

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

  /* sarrayp == false */
  /* Mark genome-G to query-A mismatches */
  _diff = _mm_andnot_si128(_mm_or_si128(_query_high, _query_low), _mm_andnot_si128(_ref_low, _ref_high));

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

  /* sarrayp == false */
  /* Mark genome-G to query-A mismatches */
  _diff = _mm_andnot_si128(_mm_or_si128(_query_high, _query_low), _mm_andnot_si128(_ref_low, _ref_high));

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

  /* sarrayp == false */
  /* Mark genome-G to query-A mismatches */
  _diff = _mm_andnot_si128(_mm_or_si128(_query_high, _query_low), _mm_andnot_si128(_ref_low, _ref_high));

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

  /* sarrayp == false */
  /* Mark genome-G to query-A mismatches */
  _diff = _mm256_andnot_si256(_mm256_or_si256(_query_high, _query_low), _mm256_andnot_si256(_ref_low, _ref_high));

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

  /* sarrayp == false */
  /* Mark genome-G to query-A mismatches */
  _diff = _mm512_andnot_si512(_mm512_or_si512(_query_high, _query_low), _mm512_andnot_si512(_ref_low, _ref_high));

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

  /* sarrayp == false */
  /* Mark genome-C to query-T mismatches */
#ifdef WORDS_BIGENDIAN
  diff = (query_shifted[0] & query_shifted[1]) &
    (~(Bigendian_convert_uint(ref_ptr[0])) & Bigendian_convert_uint(ref_ptr[4]));
#elif !defined(HAVE_SSE2)
  diff = (query_shifted[0] & query_shifted[1]) & (~(ref_ptr[0]) & ref_ptr[4]);
#else
  diff = (query_shifted[0] & query_shifted[4]) & (~(ref_ptr[0]) & ref_ptr[4]);
#endif
  debug(printf(" => diff %08X\n",diff));

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

  /* sarrayp == false */
  /* Mark genome-C to query-T mismatches */
  _diff = _mm_and_si128(_mm_and_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-C to query-T mismatches */
  _diff = _mm_and_si128(_mm_and_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-C to query-T mismatches */
  _diff = _mm_and_si128(_mm_and_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-C to query-T mismatches */
  _diff = _mm_and_si128(_mm_and_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-C to query-T mismatches */
  _diff = _mm_and_si128(_mm_and_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-C to query-T mismatches */
  _diff = _mm256_and_si256(_mm256_and_si256(_query_high, _query_low), _mm256_andnot_si256(_ref_high, _ref_low));

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

  /* sarrayp == false */
  /* Mark genome-C to query-T mismatches */
  _diff = _mm512_and_si512(_mm512_and_si512(_query_high, _query_low), _mm512_andnot_si512(_ref_high, _ref_low));

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
block_diff_atoi_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_atoi_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_atoi_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_atoi_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
			      bool plusp, int genestrand, bool query_unk_mismatch_local_p,
			      int endcolumni) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    } else {
      return block_diff_a2iag_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    } else {
      return block_diff_a2itc_128_shift_lo(query_shifted,ref_ptr,query_unk_mismatch_local_p,
					   endcolumni);
    }
  }
}
#endif
#endif

#ifdef HAVE_SSSE3
static __m128i
block_diff_atoi_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_atoi_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_atoi_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_atoi_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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


#ifdef GSNAP
/* Ignores snp_ptr */
static UINT4
block_diff_atoi_snp_32 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
#endif

#if defined(GSNAP) && defined(HAVE_SSE2)
/* Ignores snp_ptr */
static __m128i
block_diff_atoi_snp_128 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_atoi_snp_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
				  bool plusp, int genestrand, bool query_unk_mismatch_local_p, int startcolumni) {
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
block_diff_atoi_snp_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
				  bool plusp, int genestrand, bool query_unk_mismatch_local_p, int endcolumni) {
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

#if defined(GSNAP) && defined(HAVE_SSSE3)
/* Ignores snp_ptr */
static __m128i
block_diff_atoi_snp_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_atoi_snp_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

#if defined(GSNAP) && defined(HAVE_AVX2)
/* Ignores snp_ptr */
static __m256i
block_diff_atoi_snp_256 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

#if defined(GSNAP) && defined(HAVE_AVX512)
/* Ignores snp_ptr */
static __m512i
block_diff_atoi_snp_512 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_32 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_128 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_256 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_512 (Genomecomp_T *query_shifted, Genomecomp_T *ref_ptr,
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


#ifdef GSNAP
/* Ignores snp_ptr */
static UINT4
block_diff_ttoc_snp_32 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
#endif

#if defined(GSNAP) && defined(HAVE_SSE2)
/* Ignores snp_ptr */
static __m128i
block_diff_ttoc_snp_128 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_snp_128_shift_lo (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_snp_128_shift_hi (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

#if defined(GSNAP) && defined(HAVE_SSSE3)
/* Ignores snp_ptr */
static __m128i
block_diff_ttoc_snp_128_wrap_lo (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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
block_diff_ttoc_snp_128_wrap_hi (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

#if defined(GSNAP) && defined(HAVE_AVX2)
/* Ignores snp_ptr */
static __m256i
block_diff_ttoc_snp_256 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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

#if defined(GSNAP) && defined(HAVE_AVX512)
/* Ignores snp_ptr */
static __m512i
block_diff_ttoc_snp_512 (Genomecomp_T *query_shifted, Genomecomp_T *snp_ptr, Genomecomp_T *ref_ptr,
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


/* query_shifted, (snp_ptr,) ref_ptr, plusp, genestrand, query_unk_mismatch_local_p */
#ifdef HAVE_AVX512
typedef __m512i (*Diffproc_512_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool);
typedef __m512i (*Diffproc_snp_512_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *, bool, int, bool);
static Diffproc_512_T block_diff_512;
static Diffproc_snp_512_T block_diff_snp_512;
#endif

#ifdef HAVE_AVX2
typedef __m256i (*Diffproc_256_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool);
typedef __m256i (*Diffproc_snp_256_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *, bool, int, bool);
static Diffproc_256_T block_diff_256;
static Diffproc_snp_256_T block_diff_snp_256;
#endif

#ifdef HAVE_SSSE3
typedef __m128i (*Diffproc_128_wrap_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool, int);
typedef __m128i (*Diffproc_snp_128_wrap_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *, bool, int, bool, int);
static Diffproc_128_wrap_T block_diff_128_wrap_lo;
static Diffproc_snp_128_wrap_T block_diff_snp_128_wrap_lo;
static Diffproc_128_wrap_T block_diff_128_wrap_hi;
static Diffproc_snp_128_wrap_T block_diff_snp_128_wrap_hi;
#endif

#ifdef HAVE_SSE2
#ifdef USE_SHIFT_HILO
typedef __m128i (*Diffproc_128_shift_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool, int);
typedef __m128i (*Diffproc_snp_128_shift_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *, bool, int, bool, int);
static Diffproc_128_shift_T block_diff_128_shift_lo;
static Diffproc_snp_128_shift_T block_diff_snp_128_shift_lo;
static Diffproc_128_shift_T block_diff_128_shift_hi;
static Diffproc_snp_128_shift_T block_diff_snp_128_shift_hi;
#endif

typedef __m128i (*Diffproc_128_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool);
typedef __m128i (*Diffproc_snp_128_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *, bool, int, bool);
static Diffproc_128_T block_diff_128;
static Diffproc_snp_128_T block_diff_snp_128;
#endif

typedef UINT4 (*Diffproc_32_T) (Genomecomp_T *, Genomecomp_T *, bool, int, bool);
typedef UINT4 (*Diffproc_snp_32_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *, bool, int, bool);

static Diffproc_32_T block_diff_32;
static Diffproc_snp_32_T block_diff_snp_32;

/* For CMET and ATOI, ignores genome-to-query mismatches.  Used by
   Genome_consecutive procedures, called only by sarray-read.c */

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
Genome_hr_setup (Genomecomp_T *ref_blocks_in, Genomecomp_T *snp_blocks_in,
		 bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
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

  ref_blocks = ref_blocks_in;
  snp_blocks = snp_blocks_in;
  query_unk_mismatch_p = query_unk_mismatch_p_in;
  genome_unk_mismatch_p = genome_unk_mismatch_p_in;

  switch (mode) {
  case STANDARD:
#ifdef HAVE_AVX512
    block_diff_512 = block_diff_standard_512;
#endif
#ifdef HAVE_AVX2
    block_diff_256 = block_diff_standard_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_128_wrap_lo = block_diff_standard_128_wrap_lo;
    block_diff_128_wrap_hi = block_diff_standard_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_128 = block_diff_standard_128;
#ifdef USE_SHIFT_HILO
    block_diff_128_shift_lo = block_diff_standard_128_shift_lo;
    block_diff_128_shift_hi = block_diff_standard_128_shift_hi;
#endif
#endif
    block_diff_32 = block_diff_standard_32;
    break;

  case CMET_STRANDED: case CMET_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_512 = block_diff_cmet_512;
#endif
#ifdef HAVE_AVX2
    block_diff_256 = block_diff_cmet_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_128_wrap_lo = block_diff_cmet_128_wrap_lo;
    block_diff_128_wrap_hi = block_diff_cmet_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_128 = block_diff_cmet_128;
#ifdef USE_SHIFT_HILO
    block_diff_128_shift_lo = block_diff_cmet_128_shift_lo;
    block_diff_128_shift_hi = block_diff_cmet_128_shift_hi;
#endif
#endif
    block_diff_32 = block_diff_cmet_32;
    break;

  case ATOI_STRANDED: case ATOI_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_512 = block_diff_atoi_512;
#endif
#ifdef HAVE_AVX2
    block_diff_256 = block_diff_atoi_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_128_wrap_lo = block_diff_atoi_128_wrap_lo;
    block_diff_128_wrap_hi = block_diff_atoi_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_128 = block_diff_atoi_128;
#ifdef USE_SHIFT_HILO
    block_diff_128_shift_lo = block_diff_atoi_128_shift_lo;
    block_diff_128_shift_hi = block_diff_atoi_128_shift_hi;
#endif
#endif
    block_diff_32 = block_diff_atoi_32;
    break;

  case TTOC_STRANDED: case TTOC_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_512 = block_diff_ttoc_512;
#endif
#ifdef HAVE_AVX2
    block_diff_256 = block_diff_ttoc_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_128_wrap_lo = block_diff_ttoc_128_wrap_lo;
    block_diff_128_wrap_hi = block_diff_ttoc_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_128 = block_diff_ttoc_128;
#ifdef USE_SHIFT_HILO
    block_diff_128_shift_lo = block_diff_ttoc_128_shift_lo;
    block_diff_128_shift_hi = block_diff_ttoc_128_shift_hi;
#endif
#endif
    block_diff_32 = block_diff_ttoc_32;
    break;
  default: fprintf(stderr,"Mode %d not recognized\n",mode); abort();
  }

#ifndef GSNAP
#ifdef HAVE_AVX512
  block_diff_snp_512 = block_diff_standard_wildcard_512;
#endif
#ifdef HAVE_AVX2
  block_diff_snp_256 = block_diff_standard_wildcard_256;
#endif
#ifdef HAVE_SSE2
  block_diff_snp_128 = block_diff_standard_wildcard_128;
#endif
  block_diff_snp_32 = block_diff_standard_wildcard_32;

#else
  switch (mode) {
  case STANDARD:
#ifdef HAVE_AVX512
    block_diff_snp_512 = block_diff_standard_wildcard_512;
#endif
#ifdef HAVE_AVX2
    block_diff_snp_256 = block_diff_standard_wildcard_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_snp_128_wrap_lo = block_diff_standard_wildcard_128_wrap_lo;
    block_diff_snp_128_wrap_hi = block_diff_standard_wildcard_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_snp_128 = block_diff_standard_wildcard_128;
#ifdef USE_SHIFT_HILO
    block_diff_snp_128_shift_lo = block_diff_standard_wildcard_128_shift_lo;
    block_diff_snp_128_shift_hi = block_diff_standard_wildcard_128_shift_hi;
#endif
#endif
    block_diff_snp_32 = block_diff_standard_wildcard_32;
    break;

  case CMET_STRANDED: case CMET_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_snp_512 = block_diff_cmet_snp_512;
#endif
#ifdef HAVE_AVX2
    block_diff_snp_256 = block_diff_cmet_snp_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_snp_128_wrap_lo = block_diff_cmet_snp_128_wrap_lo;
    block_diff_snp_128_wrap_hi = block_diff_cmet_snp_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_snp_128 = block_diff_cmet_snp_128;
#ifdef USE_SHIFT_HILO
    block_diff_snp_128_shift_lo = block_diff_cmet_snp_128_shift_lo;
    block_diff_snp_128_shift_hi = block_diff_cmet_snp_128_shift_hi;
#endif
#endif
    block_diff_snp_32 = block_diff_cmet_snp_32;
    break;

  case ATOI_STRANDED: case ATOI_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_snp_512 = block_diff_atoi_snp_512;
#endif
#ifdef HAVE_AVX2
    block_diff_snp_256 = block_diff_atoi_snp_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_snp_128_wrap_lo = block_diff_atoi_snp_128_wrap_lo;
    block_diff_snp_128_wrap_hi = block_diff_atoi_snp_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_snp_128 = block_diff_atoi_snp_128;
#ifdef USE_SHIFT_HILO
    block_diff_snp_128_shift_lo = block_diff_atoi_snp_128_shift_lo;
    block_diff_snp_128_shift_hi = block_diff_atoi_snp_128_shift_hi;
#endif
#endif
    block_diff_snp_32 = block_diff_atoi_snp_32;
    break;

  case TTOC_STRANDED: case TTOC_NONSTRANDED:
#ifdef HAVE_AVX512
    block_diff_snp_512 = block_diff_ttoc_snp_512;
#endif
#ifdef HAVE_AVX2
    block_diff_snp_256 = block_diff_ttoc_snp_256;
#endif
#ifdef HAVE_SSSE3
    block_diff_snp_128_wrap_lo = block_diff_ttoc_snp_128_wrap_lo;
    block_diff_snp_128_wrap_hi = block_diff_ttoc_snp_128_wrap_hi;
#endif
#ifdef HAVE_SSE2
    block_diff_snp_128 = block_diff_ttoc_snp_128;
#ifdef USE_SHIFT_HILO
    block_diff_snp_128_shift_lo = block_diff_ttoc_snp_128_shift_lo;
    block_diff_snp_128_shift_hi = block_diff_ttoc_snp_128_shift_hi;
#endif
#endif
    block_diff_snp_32 = block_diff_ttoc_snp_32;
    break;
  default: fprintf(stderr,"Mode %d not recognized\n",mode); abort();
  }
#endif

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
#elif defined(HAVE_MM_POPCNT) && defined(HAVE_MM_POPCNT_U64)
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
/* ? Alternative: _mm_movemask_ps(_mm_castsi128_ps(diff)) != 0 */
#define nonzero_p_128(diff) _mm_movemask_epi8(_mm_cmpeq_epi8(diff,_mm_setzero_si128())) != 0xFFFF
#endif

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

#if (defined(USE_SHIFT_TRIM) && defined(HAVE_SSE2)) || (defined(USE_WRAP_TRIM) && defined(HAVE_SSSE3))
/* Based on clear_end */
static __m128i
set_start_128 (__m128i _diff, int startdiscard) {
  __m128i _mask, _startdiscard;

  debug(printf("Setting start at startdiscard %d\n",startdiscard));

#ifdef DEFECTIVE_SSE2_COMPILER
  _mask = _mm_sll_epi32(_mm_set1_epi32(~0U), _mm_setr_epi32(startdiscard % 32,0,0,0));
#else
  _mask = _mm_slli_epi32(_mm_set1_epi32(~0U), startdiscard % 32);
#endif
  _startdiscard = _mm_set1_epi32(startdiscard);
  _mask = _mm_or_si128(_mm_cmplt_epi32(_startdiscard, _BOUND_LOW), _mask);
  _mask = _mm_and_si128(_mm_cmplt_epi32(_startdiscard, _BOUND_HIGH), _mask);

  _mask = _mm_xor_si128(_mask, _mm_set1_epi32(~0U)); /* Take complement of _mask */

  return _mm_or_si128(_mask, _diff);
}
#endif


#if (defined(USE_SHIFT_TRIM) && defined(HAVE_SSE2)) || (defined(USE_WRAP_TRIM) && defined(HAVE_SSSE3))/* Based on clear_start */
static __m128i
set_end_128 (__m128i _diff, int enddiscard) {
  __m128i _mask, _enddiscard;

  debug(printf("Setting end at enddiscard %d\n",enddiscard));

#ifdef DEFECTIVE_SSE2_COMPILER
  _mask = _mm_sll_epi32(_mm_set1_epi32(~0U), _mm_setr_epi32(enddiscard % 32,0,0,0));
#else
  _mask = _mm_slli_epi32(_mm_set1_epi32(~0U), enddiscard % 32);
#endif
  _enddiscard = _mm_set1_epi32(enddiscard);
  _mask = _mm_or_si128(_mask, _mm_cmplt_epi32(_enddiscard, _BOUND_LOW));
  _mask = _mm_and_si128(_mask, _mm_cmplt_epi32(_enddiscard, _BOUND_HIGH));

  return _mm_or_si128(_mask, _diff);
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
#elif defined(HAVE_MM_POPCNT)
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

static __m128i
clear_highbit_128 (__m128i _diff, int leading_zeroes) {
  __m128i _subtract, _relpos;
  int relpos;

  relpos = 127 - leading_zeroes;
  debug3(printf("Clearing high bit at relpos %d\n",relpos));

#ifdef DEFECTIVE_SSE2_COMPILER
  _subtract = _mm_sll_epi32(_mm_set1_epi32(1), _mm_setr_epi32(relpos % 32,0,0,0));
#else
  _subtract = _mm_slli_epi32(_mm_set1_epi32(1), relpos % 32);
#endif
  _relpos = _mm_set1_epi32(relpos);
  _subtract = _mm_and_si128(_mm_cmplt_epi32(_relpos, _BOUND_HIGH), _subtract);
  _subtract = _mm_andnot_si128(_mm_cmplt_epi32(_relpos, _BOUND_LOW), _subtract);

  debug3(printf("Subtract: "));
  debug3(print_vector_hex(_subtract));
#if 0
  /* latency 1, throughput: 0.5 */
  return _mm_sub_epi32(_diff, _subtract);
#else
  /* _mm_xor_si128 also works if all other bits are 0.  latency 1, throughput: 0.33 */
  return _mm_xor_si128(_diff, _subtract);
#endif
}

/* relpos is equal to trailing_zeroes */
static __m128i
clear_lowbit_128 (__m128i _diff, int relpos) {
  __m128i _subtract, _relpos;

  debug3(printf("Clearing low bit at relpos %d\n",relpos));

#ifdef DEFECTIVE_SSE2_COMPILER
  _subtract = _mm_sll_epi32(_mm_set1_epi32(1), _mm_setr_epi32(relpos % 32,0,0,0));
#else
  _subtract = _mm_slli_epi32(_mm_set1_epi32(1), relpos % 32);
#endif
  _relpos = _mm_set1_epi32(relpos);
  _subtract = _mm_and_si128(_mm_cmplt_epi32(_relpos, _BOUND_HIGH), _subtract);
  _subtract = _mm_andnot_si128(_mm_cmplt_epi32(_relpos, _BOUND_LOW), _subtract);

  debug3(printf("Subtract: "));
  debug3(print_vector_hex(_subtract));
#if 0
  /* latency 1, throughput: 0.5 */
  return _mm_sub_epi32(_diff, _subtract);
#else
  /* _mm_xor_si128 also works if all other bits are 0.  latency 1, throughput: 0.33 */
  return _mm_xor_si128(_diff, _subtract);
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



static int
count_mismatches_limit (Compress_T query_compress, Univcoord_T left, 
			int pos5, int pos3, int max_mismatches, bool plusp, int genestrand) {
  int nmismatches = 0;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *ptr, *endptr;
  Genomecomp_T *query_shifted;
  int nshift;
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
	printf("Genome (in count_mismatches_limit) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
  ptr = &(ref_blocks[startblocki_32]);
  endptr = &(ref_blocks[endblocki_32]);

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    debug(printf("** Single block **\n"));
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    return popcount_ones_32(diff_32);

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_POPCOUNT) && defined(HAVE_SSE2)
    /* Shift */
#ifdef USE_SHIFT_HILO
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_shift_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					 startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#else
    /* Faster */
    startdiscard += startcolumni*32;
    enddiscard += endcolumni*32;

    diff_128 = (block_diff_128)(query_shifted - startcolumni,ptr - startcolumni,plusp,genestrand,query_unk_mismatch_p);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#endif

#else
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    nmismatches = popcount_ones_32(diff_32);
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    return (nmismatches + popcount_ones_32(diff_32));
#endif

#if defined(USE_WRAP_POPCOUNT) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_wrap_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    nmismatches = popcount_ones_32(diff_32);
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ptr + 24 <= endptr) {
      diff_256 = (block_diff_256)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_256(diff_256);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += 24; ptr += 24;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr + 12 <= endptr) {
      diff_128 = (block_diff_128)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_128(diff_128);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += 12; ptr += 12;
    }
#else
    while (ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
	nmismatches += popcount_ones_32(diff_32);
	if (nmismatches > max_mismatches) {
	  return nmismatches;
	}
	query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      }
      /* query_shifted += QUERY_NEXTROW; */ ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ptr < endptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    return (nmismatches + popcount_ones_32(diff_32));
  }
}


static int
count_mismatches_limit_snps (Compress_T query_compress, Univcoord_T left, 
			     int pos5, int pos3, int max_mismatches, bool plusp, int genestrand) {
  int nmismatches = 0;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *ref_ptr, *alt_ptr, *endptr;
  Genomecomp_T *query_shifted;
  int nshift;
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
	printf("Genome (in count_mismatches_limit) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
  ref_ptr = &(ref_blocks[startblocki_32]);
  alt_ptr = &(snp_blocks[startblocki_32]);
  endptr = &(ref_blocks[endblocki_32]);

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    debug(printf("** Single block **\n"));
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    return popcount_ones_32(diff_32);

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_POPCOUNT) && defined(HAVE_SSE2)
    /* Shift */
#ifdef USE_SHIFT_HILO
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_shift_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					     startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#else
    /* Faster */
    startdiscard += startcolumni*32;
    enddiscard += endcolumni*32;

    diff_128 = (block_diff_snp_128)(query_shifted - startcolumni,alt_ptr - startcolumni,ref_ptr - startcolumni,
				    plusp,genestrand,query_unk_mismatch_p);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#endif

#else
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    nmismatches = popcount_ones_32(diff_32);
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    return (nmismatches + popcount_ones_32(diff_32));
#endif

#if defined(USE_WRAP_POPCOUNT) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_wrap_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					    startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    nmismatches = popcount_ones_32(diff_32);
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ref_ptr + 24 <= endptr) {
      diff_256 = (block_diff_snp_256)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_256(diff_256);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += 24; ref_ptr += 24; alt_ptr += 24;
    }
#endif

#ifdef HAVE_SSE2
    while (ref_ptr + 12 <= endptr) {
      diff_128 = (block_diff_snp_128)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_128(diff_128);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += 12; ref_ptr += 12; alt_ptr += 12;
    }
#else
    while (ref_ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
	nmismatches += popcount_ones_32(diff_32);
	if (nmismatches > max_mismatches) {
	  return nmismatches;
	}
	query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      }
      /* query_shifted += QUERY_NEXTROW; */ ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ref_ptr < endptr) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    return (nmismatches + popcount_ones_32(diff_32));
  }
}



int
Genome_count_mismatches_limit (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
			       int max_mismatches, bool plusp, int genestrand) {

#if 0
  if (dibasep) {
    debug(printf("Dibase_count_mismatches_limit from %u+%d to %u+%d with max_mismatches %d:\n",
		 left,pos5,left,pos3,max_mismatches));

    return Dibase_count_mismatches_limit(&(*ncolordiffs),query,pos5,pos3,
					 /*startpos*/left+pos5,/*endpos*/left+pos3,max_mismatches);
  }
#endif

  if (snp_blocks == NULL) {
    return count_mismatches_limit(query_compress,left,pos5,pos3,max_mismatches,plusp,genestrand);
  } else {
    return count_mismatches_limit_snps(query_compress,left,pos5,pos3,max_mismatches,plusp,genestrand);
  }
}



int
Genome_count_mismatches_substring_ref (Genome_T ome, Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				       bool plusp, int genestrand) {
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  int nmismatches = 0;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *ptr, *endptr;
  Genomecomp_T *query_shifted;
  int nshift;
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
	printf("Genome (in count_mismatches_substring) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
  ptr = &(ref_blocks[startblocki_32]);
  endptr = &(ref_blocks[endblocki_32]);

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    debug(printf("(1) Returning %d mismatches\n",popcount_ones_32(diff_32)));
    return popcount_ones_32(diff_32);

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_POPCOUNT) && defined(HAVE_SSE2)
    /* Shift */
#ifdef USE_SHIFT_HILO
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_shift_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					 startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    debug(printf("(2) Returning %d mismatches\n",popcounts_ones_128(diff_128)));
    return popcount_ones_128(diff_128);
#else
    /* Faster */
    startdiscard += startcolumni*32;
    enddiscard += endcolumni*32;

    diff_128 = (block_diff_128)(query_shifted - startcolumni,ptr - startcolumni,plusp,genestrand,query_unk_mismatch_p);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    debug(printf("(3) Returning %d mismatches\n",popcount_ones_128(diff_128)));
    return popcount_ones_128(diff_128);
#endif

#else
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    nmismatches = popcount_ones_32(diff_32);
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("(4) Returning %d mismatches\n",nmismatches + popcount_ones_32(diff_32)));
    return (nmismatches + popcount_ones_32(diff_32));
#endif

#if defined(USE_WRAP_POPCOUNT) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_wrap_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    debug(printf("(5) Returning %d mismatches\n",popcount_ones_128(diff_128)));
    return popcount_ones_128(diff_128);
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    nmismatches = popcount_ones_32(diff_32);
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ptr + 24 <= endptr) {
      diff_256 = (block_diff_256)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_256(diff_256);
      query_shifted += 24; ptr += 24;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr + 12 <= endptr) {
      diff_128 = (block_diff_128)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_128(diff_128);
      query_shifted += 12; ptr += 12;
    }
#else
    while (ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
	nmismatches += popcount_ones_32(diff_32);
	query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      }
      /* query_shifted += QUERY_NEXTROW; */ ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ptr < endptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("(6) Returning %d mismatches\n",nmismatches + popcount_ones_32(diff_32)));
    return (nmismatches + popcount_ones_32(diff_32));
  }
}

static int
count_mismatches_substring_snps (Genome_T ome, Genome_T ome_alt, Compress_T query_compress, Univcoord_T left,
				 int pos5, int pos3, bool plusp, int genestrand) {
#ifdef DEBUG14
  int answer;
#endif
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  Genomecomp_T *snp_blocks = Genome_blocks(ome_alt);
  int nmismatches;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *ref_ptr, *alt_ptr, *endptr;
  Genomecomp_T *query_shifted;
  int nshift;
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
	printf("Genome (in count_mismatches_substring_snps) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	print_blocks_snp(ref_blocks,snp_blocks,left+pos5,left+pos3);
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
  ref_ptr = &(ref_blocks[startblocki_32]);
  alt_ptr = &(snp_blocks[startblocki_32]);
  endptr = &(ref_blocks[endblocki_32]);

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    debug(printf("** Single block **\n"));
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    return popcount_ones_32(diff_32);

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_POPCOUNT) && defined(HAVE_SSE2)
    /* Shift */
#ifdef USE_SHIFT_HILO
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_shift_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					     startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#else
    /* Faster */
    startdiscard += startcolumni*32;
    enddiscard += endcolumni*32;

    diff_128 = (block_diff_snp_128)(query_shifted - startcolumni,alt_ptr - startcolumni,ref_ptr - startcolumni,
                                    plusp,genestrand,query_unk_mismatch_p);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#endif

#else
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    nmismatches = popcount_ones_32(diff_32);
    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    return (nmismatches + popcount_ones_32(diff_32));
#endif

#if defined(USE_WRAP_POPCOUNT) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_wrap_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					    startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    return popcount_ones_128(diff_128);
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    nmismatches = popcount_ones_32(diff_32);
    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ref_ptr + 24 <= endptr) {
      diff_256 = (block_diff_snp_256)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_256(diff_256);
      query_shifted += 24; ref_ptr += 24; alt_ptr += 24;
    }
#endif

#ifdef HAVE_SSE2
    while (ref_ptr + 12 <= endptr) {
      diff_128 = (block_diff_snp_128)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_128(diff_128);
      query_shifted += 12; ref_ptr += 12; alt_ptr += 12;
    }
#else
    while (ref_ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
	nmismatches += popcount_ones_32(diff_32);
	query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      }
      /* query_shifted += QUERY_NEXTROW; */  ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ref_ptr < endptr) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      nmismatches += popcount_ones_32(diff_32);
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    return (nmismatches + popcount_ones_32(diff_32));
  }
}


/* left is where the start of the query matches.  pos5 is where we
   want to start comparing in the query.  pos3 is just after where we
   want to stop comparing in the query, i.e., stop at (pos3-1)
   inclusive */
int
Genome_count_mismatches_substring (Genome_T ome, Genome_T ome_alt, Compress_T query_compress, Univcoord_T left,
				   int pos5, int pos3, bool plusp, int genestrand) {

#if 0
  if (dibasep) {
    Dibase_count_mismatches_substring(&ncolordiffs,query,pos5,pos3,
				      /*startpos*/left+pos5,/*endpos*/left+pos3);
  }
#endif

  if (snp_blocks == NULL) {
    return Genome_count_mismatches_substring_ref(ome,query_compress,left,pos5,pos3,plusp,genestrand);
  } else {
    return count_mismatches_substring_snps(ome,ome_alt,query_compress,left,pos5,pos3,plusp,genestrand);
  }
}


/* pos5 is where we want to start comparing in the query.  pos3 is
   just after where we want to stop comparing in the query, i.e., stop
   at (pos3-1) inclusive */
int
Genome_count_mismatches_fragment_left (Compress_T query_compress, int pos5, int pos3,
				       Genomecomp_T ref_fragment, Genomecomp_T alt_fragment) {
  Genomecomp_T diff, alt_diff, mask;
  int startdiscard;
  Genomecomp_T query_high, query_low, query_flags;
  Genomecomp_T ref_high, ref_low, alt_high, alt_low;

  Compress_get_16mer_left(&query_high,&query_low,&query_flags,query_compress,pos3);
  startdiscard = 16 - (pos3 - pos5);

  mask = clear_start_mask(startdiscard);
  mask &= 0x0000FFFF;		/* Therefore, result of Compress does not need masking */
  debug1(printf("Mask for startdiscard %d: %08X\n",startdiscard,mask));


  /* Unpack genomic fragments */
  ref_high = ref_fragment >> 16;
  ref_low = ref_fragment /* & 0x0000FFFF */;

  alt_high = alt_fragment >> 16;
  alt_low = alt_fragment /* & 0x0000FFFF */;


  debug1(printf("Comparing: query high %08X, low %08X with ref fragment high %08X, %08X\n",query_high & 0xFFFF,query_low & 0xFFFF,ref_high & 0xFFFF,ref_low & 0xFFFF));

  /* Taken from block_diff */
  diff = (query_high ^ ref_high) | (query_low ^ ref_low);
  debug1(printf(" => ref_diff %04X",(unsigned short) diff));

  alt_diff = (query_high ^ alt_high) | (query_low ^ alt_low);
  debug1(printf(" and alt_diff %04X\n",(unsigned short) alt_diff));

  diff &= alt_diff;

  diff |= query_flags;

  diff &= mask;

  assert(diff <= 0x0000FFFF);

#if !defined(HAVE_SSE4_2)
  debug1(printf("nmismatches %08X => %d\n",diff,count_bits[diff]));
  return count_bits[diff];
#elif defined(HAVE_POPCNT)
  debug1(printf("nmismatches %08X => %d\n",diff,_popcnt32(diff)));
  return _popcnt32(diff);
#elif defined(HAVE_MM_POPCNT)
  debug1(printf("nmismatches %08X => %d\n",diff,_popcnt32(diff)));
  return _mm_popcnt_u32(diff);
#elif defined(HAVE_BUILTIN_POPCOUNT)
  debug1(printf("nmismatches %08X => %d\n",diff,__builtin_popcount(diff)));
  return __builtin_popcount(diff);
#else
  debug1(printf("nmismatches %08X => %d\n",diff,count_bits[diff]));
  return count_bits[diff];
#endif
}


/* pos5 is where we want to start comparing in the query.  pos3 is
   just after where we want to stop comparing in the query, i.e., stop
   at (pos3-1) inclusive */
int
Genome_count_mismatches_fragment_right (Compress_T query_compress, int pos5, int pos3,
					Genomecomp_T ref_fragment, Genomecomp_T alt_fragment) {
  Genomecomp_T diff, alt_diff, mask;
  int enddiscard;
  Genomecomp_T query_high, query_low, query_flags;
  Genomecomp_T ref_high, ref_low, alt_high, alt_low;

  Compress_get_16mer_right(&query_high,&query_low,&query_flags,query_compress,pos5);
  enddiscard = pos3 - pos5;

  mask = clear_end_mask(enddiscard);
  mask &= 0x0000FFFF;		/* Therefore, result of Compress does not need masking */
  debug1(printf("Mask for enddiscard %d: %08X\n",enddiscard,mask));


  /* Unpack genomic fragments */
  ref_high = ref_fragment >> 16;
  ref_low = ref_fragment /* & 0x0000FFFF */;

  alt_high = alt_fragment >> 16;
  alt_low = alt_fragment /* & 0x0000FFFF */;


  debug1(printf("Comparing: query high %08X, low %08X with ref fragment high %08X, %08X\n",query_high & 0xFFFF,query_low & 0xFFFF,ref_high & 0xFFFF,ref_low & 0xFFFF));

  /* Taken from block_diff */
  diff = (query_high ^ ref_high) | (query_low ^ ref_low);
  debug1(printf(" => ref_diff %08X",diff));

  alt_diff = (query_high ^ alt_high) | (query_low ^ alt_low);
  debug1(printf(" and alt_diff %08X\n",alt_diff));

  diff &= alt_diff;

  diff |= query_flags;

  diff &= mask;

  assert(diff <= 0x0000FFFF);

#if !defined(HAVE_SSE4_2)
  debug1(printf("nmismatches %08X => %d\n",diff,count_bits[diff]));
  return count_bits[diff];
#elif defined(HAVE_POPCNT)
  debug1(printf("nmismatches %08X => %d\n",diff,_popcnt32(diff)));
  return _popcnt32(diff);
#elif defined(HAVE_MM_POPCNT)
  debug1(printf("nmismatches %08X => %d\n",diff,_popcnt32(diff)));
  return _mm_popcnt_u32(diff);
#elif defined(HAVE_BUILTIN_POPCOUNT)
  debug1(printf("nmismatches %08X => %d\n",diff,__builtin_popcount(diff)));
  return __builtin_popcount(diff);
#else
  debug1(printf("nmismatches %08X => %d\n",diff,count_bits[diff]));
  return count_bits[diff];
#endif
}



/* ome here needs to be genomebits */
static int
mismatches_left (int *mismatch_positions, int max_mismatches, Genome_T ome, Compress_T query_compress,
		 Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand,
		 bool query_unk_mismatch_local_p) {
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  int nmismatches = 0, offset, nshift;
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
	printf("Entered mismatches_left with %d max_mismatches\n",max_mismatches);
	printf("Genome (in mismatches_left) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    return nmismatches;

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT) && defined(HAVE_SSE2)
    /* Shift */
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_shift_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p,
					 startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_128(diff_128));
      diff_128 = clear_lowbit_128(diff_128,relpos);
    }
    return nmismatches;

#else
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }

    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    return nmismatches;
#endif

#if defined(USE_WRAP) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_wrap_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p,
					startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_128(diff_128));
      diff_128 = clear_lowbit_128(diff_128,relpos);
    }
    return nmismatches;
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
    
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
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
      diff_256 = (block_diff_256)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_256(diff_256) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_256(diff_256));
	diff_256 = clear_lowbit_256(diff_256,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += 24; ptr += 24;
      offset += 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr + 12 <= endptr) {
      diff_128 = (block_diff_128)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_128(diff_128));
	diff_128 = clear_lowbit_128(diff_128,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += 12; ptr += 12;
      offset += 128;
    }
#else
    while (ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

	while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	  mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	  diff_32 = clear_lowbit_32(diff_32,relpos);
	}
	if (nmismatches > max_mismatches) {
	  return nmismatches;
	}
      
	query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
	offset += 32;
      }
      /* query_shifted += QUERY_NEXTROW; */ ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ptr < endptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    return nmismatches;
  }
}

/* Returns nmismatches_both */
static int
mismatches_left_snps (int *mismatch_positions, int max_mismatches, Genome_T ome, Genome_T ome_alt,
		      Compress_T query_compress,
		      Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand,
		      bool query_unk_mismatch_local_p) {
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  Genomecomp_T *snp_blocks = Genome_blocks(ome_alt);
  int nmismatches = 0, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *query_shifted, *ref_ptr, *alt_ptr, *endptr;
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
	printf("Entered mismatches_left_snps with %d max_mismatches\n",max_mismatches);
	printf("Genome (in mismatches_left_snps) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
  ref_ptr = &(ref_blocks[startblocki_32]);
  alt_ptr = &(snp_blocks[startblocki_32]);
  endptr = &(ref_blocks[endblocki_32]);

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    return nmismatches;

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT) && defined(HAVE_SSE2)
    /* Shift */
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_shift_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p,
					     startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_128(diff_128));
      diff_128 = clear_lowbit_128(diff_128,relpos);
    }
    return nmismatches;

#else
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }

    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    return nmismatches;
#endif

#if defined(USE_WRAP) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_wrap_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p,
					    startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_128(diff_128));
      diff_128 = clear_lowbit_128(diff_128,relpos);
    }
    return nmismatches;
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
    
    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      offset += 32;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ref_ptr + 24 <= endptr) {
      diff_256 = (block_diff_snp_256)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_256(diff_256) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_256(diff_256));
	diff_256 = clear_lowbit_256(diff_256,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += 24; ref_ptr += 24; alt_ptr += 24;
      offset += 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ref_ptr + 12 <= endptr) {
      diff_128 = (block_diff_snp_128)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_128(diff_128));
	diff_128 = clear_lowbit_128(diff_128,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += 12; ref_ptr += 12; alt_ptr += 12;
      offset += 128;
    }
#else
    while (ref_ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
	
	while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	  mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	  diff_32 = clear_lowbit_32(diff_32,relpos);
	}
	if (nmismatches > max_mismatches) {
	  return nmismatches;
	}
	
	query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
	offset += 32;
      }
      /* query_shifted += QUERY_NEXTROW; */ ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ref_ptr < endptr) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    return nmismatches;
  }
}


/* Returns mismatch_positions[0..nmismatches], where nmismatches <= max_mismatches + 1 */
/* If request max_mismatches 3, could return 0, 1, 2, 3, or 4.  Array
   contains one more element than the return value.  Therefore, need
   to supply array with (max_mismatches + 2) spaces */
int
Genome_mismatches_left (int *mismatch_positions, int max_mismatches, Genome_T ome, Genome_T ome_alt,
			Compress_T query_compress,
			Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand) {
  int nmismatches;
#ifdef DEBUG
  int i;
#endif

#if 0
  if (dibasep) {
    debug(printf("Dibase_mismatches_left from %u+%d to %u+%d:\n",left,pos5,left,pos3));

    nmismatches = Dibase_mismatches_left(&(*mismatch_positions),&(*colordiffs),max_mismatches,query,
					 pos5,pos3,/*startpos*/left+pos5,/*endpos*/left+pos3);
    mismatch_positions[nmismatches] = pos3 + 1;	/* Need +1 because of starting assumed nt */

  }
#endif

  if (ome_alt == NULL) {
    nmismatches = mismatches_left(&(*mismatch_positions),max_mismatches,ome,query_compress,
				  left,pos5,pos3,plusp,genestrand,query_unk_mismatch_p);
    mismatch_positions[nmismatches] = pos3;
  } else {
    nmismatches = mismatches_left_snps(&(*mismatch_positions),max_mismatches,ome,ome_alt,query_compress,
				       left,pos5,pos3,plusp,genestrand,query_unk_mismatch_p);
    mismatch_positions[nmismatches] = pos3;
  }
  debug(
	printf("%d mismatches on left: ",nmismatches);
	for (i = 0; i <= nmismatches; i++) {
	  printf("%d ",mismatch_positions[i]);
	}
	printf("\n");
	);
  
  return nmismatches;
}


/* Returns mismatch_positions[0..nmismatches], where nmismatches <= max_mismatches + 1 */
/* If request max_mismatches 3, could return m0, m1, m2, m3, m4 */
/* See note above about why we set query_unk_mismatch_p to false */
int
Genome_mismatches_left_trim (int *mismatch_positions, int max_mismatches, Genome_T ome, Genome_T ome_alt,
			     Compress_T query_compress,
			     Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand) {
  int nmismatches;
#ifdef DEBUG
  int i;
#endif

#if 0
  if (dibasep) {
    debug(printf("Dibase_mismatches_left from %u+%d to %u+%d:\n",left,pos5,left,pos3));

    nmismatches = Dibase_mismatches_left(&(*mismatch_positions),&(*colordiffs),max_mismatches,query,
					 pos5,pos3,/*startpos*/left+pos5,/*endpos*/left+pos3);
    mismatch_positions[nmismatches] = pos3 + 1;	/* Need +1 because of starting assumed nt */

  }
#endif

  if (ome_alt == NULL) {
    nmismatches = mismatches_left(&(*mismatch_positions),max_mismatches,ome,query_compress,
				  left,pos5,pos3,plusp,genestrand,/*query_unk_mismatch_p*/false);
    mismatch_positions[nmismatches] = pos3;
  } else {
    nmismatches = mismatches_left_snps(&(*mismatch_positions),max_mismatches,ome,ome_alt,query_compress,
				       left,pos5,pos3,plusp,genestrand,/*query_unk_mismatch_p*/false);
    mismatch_positions[nmismatches] = pos3;
  }
  debug(
	printf("%d mismatches on left: ",nmismatches);
	for (i = 0; i <= nmismatches; i++) {
	  printf("%d ",mismatch_positions[i]);
	}
	printf("\n");
	);
  
  return nmismatches;
}


static int
mismatches_right (int *mismatch_positions, int max_mismatches, Genome_T ome, Compress_T query_compress,
		  Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand,
		  bool query_unk_mismatch_local_p) {
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  int nmismatches = 0, offset, relpos, nshift;
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


  debug(
	printf("\n\n");
	printf("Entered mismatches_right with %d max_mismatches\n",max_mismatches);
	printf("Genome (in mismatches_right) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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

    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    return nmismatches;

  } else if (startblocki == endblocki) {
#if defined(USE_SHIFT_MISMATCH_POSITIONS) && defined(HAVE_SSE2)
    /* Shift */
    startdiscard += 96 - (endcolumni - startcolumni)*32;
    enddiscard += 96;
    diff_128 = (block_diff_128_shift_hi)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p,
					 endcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_128(diff_128));
      diff_128 = clear_highbit_128(diff_128,relpos);
    }
    return nmismatches;

#else
    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }

    query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    offset -= 32;

    /* Single row */
    while (--endcolumni > startcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }

      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }

    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    return nmismatches;
#endif

#if defined(USE_WRAP_MISMATCH_POSITIONS) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    startdiscard += (startcolumni - endcolumni - 1)*32;
    enddiscard += 96;
    diff_128 = (block_diff_128_wrap_hi)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p,
					endcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_128(diff_128));
      diff_128 = clear_highbit_128(diff_128,relpos);
    }
    return nmismatches;
#endif

  } else {
    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
    query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    offset -= 32;

    /* End row */
    while (--endcolumni >= 0) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
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
      diff_256 = (block_diff_256)(&(query_shifted[-15]),&(ptr[-15]),plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_256(diff_256) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_256(diff_256));
	diff_256 = clear_highbit_256(diff_256,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted -= 24; ptr -= 24;
      offset -= 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr >= startptr + 12) {
      diff_128 = (block_diff_128)(&(query_shifted[-3]),&(ptr[-3]),plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_128(diff_128));
	diff_128 = clear_highbit_128(diff_128,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted -= 12; ptr -= 12;
      offset -= 128;
    }
#else
    while (ptr >= startptr + 12) {
      for (endcolumni = 3; endcolumni >= 0; --endcolumni) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

	while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	  mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	  diff_32 = clear_highbit_32(diff_32,relpos);
	}
	if (nmismatches > max_mismatches) {
	  return nmismatches;
	}
	query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
	offset -= 32;
      }
      /* query_shifted -= QUERY_NEXTROW; */ ptr -= GENOME_NEXTROW;
    }
#endif

    /* Start row */
    while (ptr > startptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }

    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    return nmismatches;
  }
}

/* Returns nmismatches_both */
static int
mismatches_right_snps (int *mismatch_positions, int max_mismatches, Genome_T ome, Genome_T ome_alt,
		       Compress_T query_compress,
		       Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand,
		       bool query_unk_mismatch_local_p) {
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  Genomecomp_T *snp_blocks = Genome_blocks(ome_alt);
  int nmismatches = 0, offset, relpos, nshift;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *query_shifted, *ref_ptr, *alt_ptr, *startptr;
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


  debug(
	printf("\n\n");
	printf("Entered mismatches_right_snps with %d max_mismatches\n",max_mismatches);
	printf("Genome (in mismatches_right_snps) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
  ref_ptr = &(ref_blocks[endblocki_32]);
  alt_ptr = &(snp_blocks[endblocki_32]);
  startptr = &(ref_blocks[startblocki_32]);

  if (startblocki_32 == endblocki_32) {
    /* Single block */
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    return nmismatches;

  } else if (startblocki == endblocki) {
#if defined(USE_SHIFT_MISMATCH_POSITIONS) && defined(HAVE_SSE2)
    /* Shift */
    startdiscard += 96 - (endcolumni - startcolumni)*32;
    enddiscard += 96;
    diff_128 = (block_diff_snp_128_shift_hi)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p,
					     endcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_128(diff_128));
      diff_128 = clear_highbit_128(diff_128,relpos);
    }
    return nmismatches;

#else
    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }

    query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
    offset -= 32;

    /* Single row */
    while (--endcolumni > startcolumni) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }

      query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }

    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    return nmismatches;
#endif

#if defined(USE_WRAP_MISMATCH_POSITIONS) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    startdiscard += (startcolumni - endcolumni - 1)*32;
    enddiscard += 96;
    diff_128 = (block_diff_snp_128_wrap_hi)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p,
					    endcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_128(diff_128));
      diff_128 = clear_highbit_128(diff_128,relpos);
    }
    return nmismatches;
#endif

  } else {
    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
    query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
    offset -= 32;

    /* End row */
    while (--endcolumni >= 0) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }
#ifdef HAVE_SSE2
    query_shifted -= QUERY_NEXTROW;
#endif
    ref_ptr -= GENOME_NEXTROW; alt_ptr -= GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ref_ptr >= startptr + 24) {
      diff_256 = (block_diff_snp_256)(&(query_shifted[-15]),&(alt_ptr[-15]),&(ref_ptr[-15]),plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_256(diff_256) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_256(diff_256));
	diff_256 = clear_highbit_256(diff_256,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted -= 24; ref_ptr -= 24; alt_ptr -= 24;
      offset -= 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ref_ptr >= startptr + 12) {
      diff_128 = (block_diff_snp_128)(&(query_shifted[-3]),&(alt_ptr[-3]),&(ref_ptr[-3]),plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_128(diff_128) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_128(diff_128));
	diff_128 = clear_highbit_128(diff_128,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted -= 12; ref_ptr -= 12; alt_ptr -= 12;
      offset -= 128;
    }
#else
    while (ref_ptr >= startptr + 12) {
      for (endcolumni = 3; endcolumni >= 0; --endcolumni) {
	diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

	while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	  mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	  diff_32 = clear_highbit_32(diff_32,relpos);
	}
	if (nmismatches > max_mismatches) {
	  return nmismatches;
	}
	query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
	offset -= 32;
      }
      /* query_shifted -= QUERY_NEXTROW; */ ref_ptr -= GENOME_NEXTROW; alt_ptr -= GENOME_NEXTROW;
    }
#endif

    /* Start row */
    while (ref_ptr > startptr) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
      query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }

    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    return nmismatches;
  }
}



/* Returns mismatch_positions[0..nmismatches], where nmismatches <= max_mismatches */
int
Genome_mismatches_right (int *mismatch_positions, int max_mismatches, Genome_T ome, Genome_T ome_alt,
			 Compress_T query_compress,
			 Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand) {
  int nmismatches;
#ifdef DEBUG
  int i;
#endif

#if 0
  if (dibasep) {
    debug(printf("Dibase_mismatches_right from %u+%d to %u+%d:\n",left,pos5,left,pos3));

    nmismatches = Dibase_mismatches_right(&(*mismatch_positions),&(*colordiffs),max_mismatches,query,
					  pos5,pos3,/*startpos*/left+pos5,/*endpos*/left+pos3);
  }
#endif

  if (ome_alt == NULL) {
    nmismatches = mismatches_right(&(*mismatch_positions),max_mismatches,ome,query_compress,
				   left,pos5,pos3,plusp,genestrand,query_unk_mismatch_p);
  } else {
    nmismatches = mismatches_right_snps(&(*mismatch_positions),max_mismatches,ome,ome_alt,query_compress,
					left,pos5,pos3,plusp,genestrand,query_unk_mismatch_p);
  }
  mismatch_positions[nmismatches] = -1;
  debug(
	printf("%d mismatches on right: ",nmismatches);
	for (i = 0; i <= nmismatches; i++) {
	  printf("%d ",mismatch_positions[i]);
	}
	printf("\n");
	);
  return nmismatches;
}


/* Returns mismatch_positions[0..nmismatches], where nmismatches <= max_mismatches */
/* See note above about why we set query_unk_mismatch_p to false */
int
Genome_mismatches_right_trim (int *mismatch_positions, int max_mismatches, Genome_T ome, Genome_T ome_alt,
			      Compress_T query_compress,
			      Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand) {
  int nmismatches;
#ifdef DEBUG
  int i;
#endif

#if 0
  if (dibasep) {
    debug(printf("Dibase_mismatches_right from %u+%d to %u+%d:\n",left,pos5,left,pos3));

    nmismatches = Dibase_mismatches_right(&(*mismatch_positions),&(*colordiffs),max_mismatches,query,
					  pos5,pos3,/*startpos*/left+pos5,/*endpos*/left+pos3);
  }
#endif

  if (ome_alt == NULL) {
    nmismatches = mismatches_right(&(*mismatch_positions),max_mismatches,ome,query_compress,
				   left,pos5,pos3,plusp,genestrand,/*query_unk_mismatch_p*/false);
  } else {
    nmismatches = mismatches_right_snps(&(*mismatch_positions),max_mismatches,ome,ome_alt,query_compress,
					left,pos5,pos3,plusp,genestrand,/*query_unk_mismatch_p*/false);
  }
  mismatch_positions[nmismatches] = -1;
  debug(
	printf("%d mismatches on right: ",nmismatches);
	for (i = 0; i <= nmismatches; i++) {
	  printf("%d ",mismatch_positions[i]);
	}
	printf("\n");
	);
  return nmismatches;
}


int
Genome_first_kmer_left (int *nmismatches5, Genome_T ome, Compress_T query_compress,
			Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand,
			bool query_unk_mismatch_local_p, int kmer) {
  int mismatch_position, last_mismatch;
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  int offset, nshift;
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
	printf("Entered Genome_first_kmer_left with pos5 %d and pos3 %d\n",pos5,pos3);
	printf("Genome (in first_kmer_left) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	print_blocks(ref_blocks,left+pos5,left+pos3);
	printf("\n");
	);


  *nmismatches5 = -1;

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

  last_mismatch = mismatch_position = pos5 - 1;

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      *nmismatches5 += 1;
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    debug(printf("1 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    } else {
      *nmismatches5 += 1;
      return (mismatch_position + 1);
    }

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT) && defined(HAVE_SSE2)
    /* Shift */
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_shift_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p,
					 startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
      *nmismatches5 += 1;
      diff_128 = clear_lowbit_128(diff_128,relpos);
    }
    debug(printf("2 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    } else {
      *nmismatches5 += 1;
      return (mismatch_position + 1);
    }

#else
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      *nmismatches5 += 1;
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    debug(printf("3 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    }

    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	*nmismatches5 += 1;
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      debug(printf("4 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (mismatch_position - last_mismatch > kmer) {
	return (last_mismatch + 1);
      }
      
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      *nmismatches5 += 1;
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    debug(printf("5 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    } else {
      *nmismatches5 += 1;
      return (mismatch_position + 1);
    }
#endif

#if defined(USE_WRAP) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_wrap_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p,
					startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
      *nmismatches5 += 1;
      diff_128 = clear_lowbit_128(diff_128,relpos);
    }
    debug(printf("6 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    } else {
      *nmismatches5 += 1;
      return (mismatch_position + 1);
    }
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      *nmismatches5 += 1;
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    debug(printf("7 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    }
    
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	*nmismatches5 += 1;
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      debug(printf("8 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (mismatch_position - last_mismatch > kmer) {
	return (last_mismatch + 1);
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
      diff_256 = (block_diff_256)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_256(diff_256) && mismatch_position - last_mismatch <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset + (relpos = count_trailing_zeroes_256(diff_256));
	*nmismatches5 += 1;
	diff_256 = clear_lowbit_256(diff_256,relpos);
      }
      debug(printf("9 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (mismatch_position - last_mismatch > kmer) {
	return (last_mismatch + 1);
      }
      
      query_shifted += 24; ptr += 24;
      offset += 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr + 12 <= endptr) {
      diff_128 = (block_diff_128)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_128(diff_128) && mismatch_position - last_mismatch <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
	*nmismatches5 += 1;
	diff_128 = clear_lowbit_128(diff_128,relpos);
      }
      debug(printf("10 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (mismatch_position - last_mismatch > kmer) {
	return (last_mismatch + 1);
      }
      
      query_shifted += 12; ptr += 12;
      offset += 128;
    }
#else
    while (ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

	while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
	  last_mismatch = mismatch_position;
	  mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	  *nmismatches5 += 1;
	  diff_32 = clear_lowbit_32(diff_32,relpos);
	}
	debug(printf("11 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
	if (mismatch_position - last_mismatch > kmer) {
	  return (last_mismatch + 1);
	}
      
	query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
	offset += 32;
      }
      /* query_shifted += QUERY_NEXTROW; */ ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ptr < endptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	*nmismatches5 += 1;
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      debug(printf("12 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (mismatch_position - last_mismatch > kmer) {
	return (last_mismatch + 1);
      }
      
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      *nmismatches5 += 1;
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    debug(printf("13 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    } else {
      *nmismatches5 += 1;
      return (mismatch_position + 1);
    }
  }
}




int
Genome_first_kmer_right (int *nmismatches3, Genome_T ome, Compress_T query_compress,
			 Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand,
			 bool query_unk_mismatch_local_p, int kmer) {
  int mismatch_position, last_mismatch;
  Genomecomp_T *ref_blocks = Genome_blocks(ome);
  int offset, relpos, nshift;
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


  debug(
	printf("\n\n");
	printf("Entered Genome_first_kmer_right with pos5 %d and pos3 %d\n",pos5,pos3);
	printf("Genome (in first_kmer_right) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	print_blocks(ref_blocks,left+pos5,left+pos3);
	printf("\n");
	);

  *nmismatches3 = -1;

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

  last_mismatch = mismatch_position = pos3;

  if (startblocki_32 == endblocki_32) {
    /* Single block */
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      *nmismatches3 += 1;
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    debug(printf("1 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    } else {
      *nmismatches3 += 1;
      return mismatch_position;
    }

  } else if (startblocki == endblocki) {
#if defined(USE_SHIFT_MISMATCH_POSITIONS) && defined(HAVE_SSE2)
    /* Shift */
    startdiscard += 96 - (endcolumni - startcolumni)*32;
    enddiscard += 96;
    diff_128 = (block_diff_128_shift_hi)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p,
					 endcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_128(diff_128));
      *nmismatches3 += 1;
      diff_128 = clear_highbit_128(diff_128,relpos);
    }
    debug(printf("2 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    } else {
      *nmismatches3 += 1;
      return mismatch_position;
    }

#else
    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      *nmismatches3 += 1;
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    debug(printf("3 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    }
    query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    offset -= 32;

    /* Single row */
    while (--endcolumni > startcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	*nmismatches3 += 1;
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      debug(printf("4 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (last_mismatch - mismatch_position > kmer) {
	return last_mismatch;
      }
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }

    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      *nmismatches3 += 1;
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    debug(printf("5 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    } else {
      *nmismatches3 += 1;
      return mismatch_position;
    }
#endif

#if defined(USE_WRAP_MISMATCH_POSITIONS) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    startdiscard += (startcolumni - endcolumni - 1)*32;
    enddiscard += 96;
    diff_128 = (block_diff_128_wrap_hi)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p,
					endcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_128(diff_128));
      *nmismatches3 += 1;
      diff_128 = clear_highbit_128(diff_128,relpos);
    }
    debug(printf("6 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    } else {
      *nmismatches3 += 1;
      return mismatch_position;
    }
#endif

  } else {
    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      *nmismatches3 += 1;
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    debug(printf("7 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    }
    query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    offset -= 32;

    /* End row */
    while (--endcolumni >= 0) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	*nmismatches3 += 1;
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      debug(printf("8 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (last_mismatch - mismatch_position > kmer) {
	return last_mismatch;
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
      diff_256 = (block_diff_256)(&(query_shifted[-15]),&(ptr[-15]),plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_256(diff_256) && last_mismatch - mismatch_position <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset - (relpos = count_leading_zeroes_256(diff_256));
	*nmismatches3 += 1;
	diff_256 = clear_highbit_256(diff_256,relpos);
      }
      debug(printf("9 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (last_mismatch - mismatch_position > kmer) {
	return last_mismatch;
      }
      query_shifted -= 24; ptr -= 24;
      offset -= 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr >= startptr + 12) {
      diff_128 = (block_diff_128)(&(query_shifted[-3]),&(ptr[-3]),plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_128(diff_128) && last_mismatch - mismatch_position <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset - (relpos = count_leading_zeroes_128(diff_128));
	*nmismatches3 += 1;
	diff_128 = clear_highbit_128(diff_128,relpos);
      }
      debug(printf("10 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (last_mismatch - mismatch_position > kmer) {
	return last_mismatch;
      }
      query_shifted -= 12; ptr -= 12;
      offset -= 128;
    }
#else
    while (ptr >= startptr + 12) {
      for (endcolumni = 3; endcolumni >= 0; --endcolumni) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

	while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
	  last_mismatch = mismatch_position;
	  mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	  *nmismatches3 += 1;
	  diff_32 = clear_highbit_32(diff_32,relpos);
	}
	debug(printf("11 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
	if (last_mismatch - mismatch_position > kmer) {
	  return last_mismatch;
	}
	query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
	offset -= 32;
      }
      /* query_shifted -= QUERY_NEXTROW; */ ptr -= GENOME_NEXTROW;
    }
#endif

    /* Start row */
    while (ptr > startptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);

      while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	*nmismatches3 += 1;
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      debug(printf("12 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (last_mismatch - mismatch_position > kmer) {
	return last_mismatch;
      }
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
      offset -= 32;
    }

    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_local_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      *nmismatches3 += 1;
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    debug(printf("13 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    } else {
      *nmismatches3 += 1;
      return mismatch_position;
    }
  }
}



/************************************************************************
 *  Marking
 ************************************************************************/

/* Derived from mismatches_left() */
int
Genome_mark_mismatches_ref (char *genomic, int querylength, Compress_T query_compress,
			    Univcoord_T left, int pos5, int pos3,
			    bool plusp, int genestrand) {
  int mismatch_position;
  int nmismatches = 0, offset, nshift;
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


  debug5(
	printf("\n\n");
	printf("genomic = %s\n",genomic);
	printf("Genome (in mark_mismatches_ref) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	print_blocks(ref_blocks,left+pos5,left+pos3);
	printf("\n");
	);


  startblocki = (left+pos5)/128U*12;
  startcolumni = ((left+pos5) % 128) / 32;
  startblocki_32 = startblocki + startcolumni;

  endblocki = (left+pos3)/128U*12;
  endcolumni = ((left+pos3) % 128) / 32;
  endblocki_32 = endblocki + endcolumni;

  debug5(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u, plusp %d, step_size %d\n",
		left,pos5,pos3,startblocki,endblocki,plusp,STEP_SIZE));

  nshift = left % STEP_SIZE;
  query_shifted = Compress_shift(query_compress,nshift);
  debug5(printf("Query shifted %d:\n",nshift));
  debug5(Compress_print_blocks(query_shifted,nshift,pos5,pos3));
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
    debug5(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    debug5(printf("genomic = %s\n",genomic));
    return nmismatches;

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_MISMATCH_POSITIONS) && defined(HAVE_SSE2)
    /* Shift */
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_shift_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					 startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
      diff_128 = clear_lowbit_128(diff_128,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    return nmismatches;

#else
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }

    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    
    return nmismatches;
#endif

#if defined(USE_WRAP_MISMATCH_POSITIONS) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_wrap_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
      diff_128 = clear_lowbit_128(diff_128,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    return nmismatches;
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }

    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
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
      diff_256 = (block_diff_256)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_256(diff_256)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_256(diff_256));
	diff_256 = clear_lowbit_256(diff_256,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += 24; ptr += 24;
      offset += 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr + 12 <= endptr) {
      diff_128 = (block_diff_128)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_128(diff_128)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
	diff_128 = clear_lowbit_128(diff_128,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += 12; ptr += 12;
      offset += 128;
    }
#else
    while (ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

	while (nonzero_p_32(diff_32)) {
	  mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	  diff_32 = clear_lowbit_32(diff_32,relpos);
	  if (plusp == false) {
	    mismatch_position = (querylength - 1) - mismatch_position;
	  }
	  genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	  nmismatches++;
	}
      
	query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
	offset += 32;
      }
      /* query_shifted += QUERY_NEXTROW; */ ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ptr < endptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    return nmismatches;
  }
}


/* Derived from mismatches_left_snps() */
/* Returns nmismatches_both */
static int
mark_mismatches_snps (char *genomic, int querylength, Compress_T query_compress,
		      Univcoord_T left, int pos5, int pos3,
		      bool plusp, int genestrand) {
  int mismatch_position;
  int nmismatches = 0, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *query_shifted, *ref_ptr, *alt_ptr, *endptr;
  int relpos;
  int startcolumni, endcolumni;
  UINT4 diff_32;
#ifdef HAVE_SSE2
  __m128i diff_128;
#endif
#ifdef HAVE_AVX2
  __m256i diff_256;
#endif


  debug5(
	printf("\n\n");
	printf("genomic = %s\n",genomic);
	printf("Genome (in mark_mismatches_ref) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	print_blocks(ref_blocks,left+pos5,left+pos3);
	printf("\n");
	);


  startblocki = (left+pos5)/128U*12;
  startcolumni = ((left+pos5) % 128) / 32;
  startblocki_32 = startblocki + startcolumni;

  endblocki = (left+pos3)/128U*12;
  endcolumni = ((left+pos3) % 128) / 32;
  endblocki_32 = endblocki + endcolumni;

  debug5(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u, plusp %d, step_size %d\n",
		left,pos5,pos3,startblocki,endblocki,plusp,STEP_SIZE));

  nshift = left % STEP_SIZE;
  query_shifted = Compress_shift(query_compress,nshift);
  debug5(printf("Query shifted %d:\n",nshift));
  debug5(Compress_print_blocks(query_shifted,nshift,pos5,pos3));
  query_shifted += (nshift+pos5)/STEP_SIZE*COMPRESS_BLOCKSIZE;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
#else
  query_shifted += startcolumni;
#endif

  startdiscard = (left+pos5) % 32;
  enddiscard = (left+pos3) % 32;
  offset = -startdiscard + pos5;
  ref_ptr = &(ref_blocks[startblocki_32]);
  alt_ptr = &(snp_blocks[startblocki_32]);
  endptr = &(ref_blocks[endblocki_32]);

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    debug5(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    debug5(printf("genomic = %s\n",genomic));
    return nmismatches;

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_MISMATCH_POSITIONS) && defined(HAVE_SSE2)
    /* Shift */
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_shift_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					     startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
      diff_128 = clear_lowbit_128(diff_128,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    debug5(printf("genomic = %s\n",genomic));
    return nmismatches;

#else
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }

    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    
    return nmismatches;
#endif

#if defined(USE_WRAP_MISMATCH_POSITIONS) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_wrap_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					    startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard);
    diff_128 = clear_end_128(diff_128,enddiscard);

    while (nonzero_p_128(diff_128)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
      diff_128 = clear_lowbit_128(diff_128,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    debug5(printf("genomic = %s\n",genomic));
    return nmismatches;
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }

    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    offset += 32;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      offset += 32;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ref_ptr + 24 <= endptr) {
      diff_256 = (block_diff_snp_256)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_256(diff_256)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_256(diff_256));
	diff_256 = clear_lowbit_256(diff_256,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += 24; ref_ptr += 24; alt_ptr += 24;
      offset += 256;
    }
#endif

#ifdef HAVE_SSE2
    while (ref_ptr + 12 <= endptr) {
      diff_128 = (block_diff_snp_128)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_128(diff_128)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_128(diff_128));
	diff_128 = clear_lowbit_128(diff_128,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += 12; ref_ptr += 12; alt_ptr += 12;
      offset += 128;
    }
#else
    while (ref_ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

	while (nonzero_p_32(diff_32)) {
	  mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	  diff_32 = clear_lowbit_32(diff_32,relpos);
	  if (plusp == false) {
	    mismatch_position = (querylength - 1) - mismatch_position;
	  }
	  genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	  nmismatches++;
	}
      
	query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
	offset += 32;
      }
      /* query_shifted += QUERY_NEXTROW; */ ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ref_ptr < endptr) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      while (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
	if (plusp == false) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	}
	genomic[mismatch_position] = tolower(genomic[mismatch_position]);
	nmismatches++;
      }
      
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      offset += 32;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      if (plusp == false) {
	mismatch_position = (querylength - 1) - mismatch_position;
      }
      genomic[mismatch_position] = tolower(genomic[mismatch_position]);
      nmismatches++;
    }
    return nmismatches;
  }
}


int
Genome_mark_mismatches (char *genomic, int querylength, Compress_T query_compress,
			Univcoord_T left, int pos5, int pos3,
			bool plusp, int genestrand) {

#if 0
  if (dibasep) {
    fprintf(stderr,"Not implemented\n");
#if 0
    debug5(printf("Dibase_mismatches_left from %u+%d to %u+%d:\n",left,pos5,left,pos3));

    nmismatches = Dibase_mismatches_left(&(*mismatch_positions),&(*colordiffs),max_mismatches,query,
					 pos5,pos3,/*startpos*/left+pos5,/*endpos*/left+pos3);
    mismatch_positions[nmismatches] = pos3 + 1;	/* Need +1 because of starting assumed nt */
#endif
    return 0;
  }
#endif

  if (snp_blocks == NULL) {
    return Genome_mark_mismatches_ref(&(*genomic),querylength,query_compress,
				      left,pos5,pos3,plusp,genestrand);
  } else {
    return mark_mismatches_snps(&(*genomic),querylength,query_compress,
				left,pos5,pos3,plusp,genestrand);
  }
}


/************************************************************************
 *  Trimming
 ************************************************************************/

#if 0
static int
trim_left_substring (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
		     bool plusp, int genestrand) {
  int startdiscard, enddiscard, offset;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *ptr, *startptr;
  Genomecomp_T *query_shifted;
  int nshift;
  int startcolumni, endcolumni;
  UINT4 diff_32;
#ifdef HAVE_SSE2
  __m128i diff_128;
  int i;
#endif
#ifdef HAVE_AVX2
  __m256i diff_256;
#endif
#ifdef HAVE_AVX2
  unsigned short array[16];
#elif defined(HAVE_SSE2)
  unsigned short array[8];
#endif



  int totalscore, bestscore, score;
  int trimpos;
  Genomecomp_T p;

  debug(
	printf("\n\n");
	printf("Genome (in trim_left_substring) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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

    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    diff_32 = set_start_32(diff_32,startdiscard);  /* puts 1 (mismatches) at start */

    p = 3*(diff_32 >> 16);
    bestscore = score_high[p];
    trimpos = offset - score_high[p+1];
    totalscore = score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    /* totalscore += score_high[p+2]; */
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset -= 16; */

    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */

  } else if (startblocki == endblocki) {
#if defined(USE_SHIFT_TRIM) && defined(HAVE_SSE2)
    /* Shift */
    startdiscard += 96 - (endcolumni - startcolumni)*32;
    enddiscard += 96;
    diff_128 = (block_diff_128_shift_hi)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					 endcolumni);
    diff_128 = clear_end_128(diff_128,enddiscard);
    diff_128 = set_start_128(diff_128,startdiscard);  /* puts 1 (mismatches) at start */
    _mm_store_si128((__m128i *) array,diff_128);

    bestscore = -100;
    for (i = 7; i >= 0; --i) {
      p = 3*array[i];
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
    }

    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */

#else
    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    p = 3*(diff_32 >> 16);
    bestscore = score_high[p];
    trimpos = offset - score_high[p+1];
    totalscore = score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    totalscore += score_high[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;
    query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;

    /* Single row */
    while (--endcolumni > startcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 >> 16);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
      
      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    }

    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = set_start_32(diff_32,startdiscard);  /* puts 1 (mismatches) at start */

    p = 3*(diff_32 >> 16);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    totalscore += score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;
    
    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    /* totalscore += score_high[p+2]; */
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset -= 16; */
    
    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */
#endif

#if defined(USE_WRAP_TRIM) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    startdiscard += (startcolumni - endcolumni - 1)*32;
    enddiscard += 96;
    diff_128 = (block_diff_128_wrap_hi)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					endcolumni);
    diff_128 = clear_end_128(diff_128,enddiscard);
    diff_128 = set_start_128(diff_128,startdiscard);  /* puts 1 (mismatches) at start */
    _mm_store_si128((__m128i *) array,diff_128);

    bestscore = -100;
    for (i = 7; i >= 0; --i) {
      p = 3*array[i];
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
    }

    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */
#endif

  } else {
    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    p = 3*(diff_32 >> 16);
    bestscore = score_high[p];
    trimpos = offset - score_high[p+1];
    totalscore = score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    /* totalscore += score_high[p+2]; */
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;

    /* End row */
    while (--endcolumni >= 0) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 >> 16);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore = score_high[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    }
#ifdef HAVE_SSE2
    query_shifted -= QUERY_NEXTROW;
#endif
    ptr -= GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ptr >= startptr + 24) {
      diff_256 = (block_diff_256)(&(query_shifted[-15]),&(ptr[-15]),plusp,genestrand,query_unk_mismatch_p);
      _mm256_store_si256((__m256i *) array,diff_256);

      for (i = 15; i >= 0; --i) {
	p = 3*array[i];
	if ((score = score_high[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset - score_high[p+1];
	}
	totalscore += score_high[p+2];
	debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset -= 16;
      }
      query_shifted -= 24; ptr -= 24;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr >= startptr + 12) {
      diff_128 = (block_diff_128)(&(query_shifted[-3]),&(ptr[-3]),plusp,genestrand,query_unk_mismatch_p);
      _mm_store_si128((__m128i *) array,diff_128);

      for (i = 7; i >= 0; --i) {
	p = 3*array[i];
	if ((score = score_high[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset - score_high[p+1];
	}
	totalscore += score_high[p+2];
	debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset -= 16;
      }
      query_shifted -= 12; ptr -= 12;
    }
#else
    while (ptr >= startptr + 12) {
      for (endcolumni = 3; endcolumni >= 0; --endcolumni) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

	p = 3*(diff_32 >> 16);
	if ((score = score_high[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset - score_high[p+1];
	}
	totalscore = score_high[p+2];
	debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset -= 16;
	
	p = 3*(diff_32 & 0x0000FFFF);
	if ((score = score_high[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset - score_high[p+1];
	}
	totalscore += score_high[p+2];
	debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset -= 16;
	query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
      }
      /* query_shifted -= QUERY_NEXTROW; */ ptr -= GENOME_NEXTROW;
    }
#endif

    /* Start row */
    while (ptr > startptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 >> 16);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore = score_high[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
      query_shifted -= QUERY_NEXTCOL; ptr -= GENOME_NEXTCOL;
    }

    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = set_start_32(diff_32,startdiscard);  /* puts 1 (mismatches) at start */

    p = 3*(diff_32 >> 16);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    totalscore += score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;
    
    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    /* totalscore += score_high[p+2]; */
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset -= 16; */
    
    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */
  }
}
#endif


#if 0
static int
trim_left_substring_snps (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
			  bool plusp, int genestrand) {
  int startdiscard, enddiscard, offset;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *ref_ptr, *alt_ptr, *startptr;
  Genomecomp_T *query_shifted;
  int nshift;
  int startcolumni, endcolumni;
  UINT4 diff_32;
#ifdef HAVE_SSE2
  __m128i diff_128;
  int i;
#endif
#ifdef HAVE_AVX2
  __m256i diff_256;
#endif
#ifdef HAVE_AVX2
  unsigned short array[16];
#elif defined(HAVE_SSE2)
  unsigned short array[8];
#endif

  int totalscore, bestscore, score;
  int trimpos;
  Genomecomp_T p;

  debug(
	printf("\n\n");
	printf("Genome (in trim_left_substring_snps) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
  ref_ptr = &(ref_blocks[endblocki_32]);
  alt_ptr = &(snp_blocks[endblocki_32]);
  startptr = &(ref_blocks[startblocki_32]);

  if (startblocki_32 == endblocki_32) {
    /* Single block */
    debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    diff_32 = set_start_32(diff_32,startdiscard);  /* puts 1 (mismatches) at start */

    p = 3*(diff_32 >> 16);
    bestscore = score_high[p];
    trimpos = offset - score_high[p+1];
    totalscore = score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    /* totalscore += score_high[p+2]; */
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset -= 16; */

    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */

  } else if (startblocki == endblocki) {
#if defined(USE_SHIFT_TRIM) && defined(HAVE_SSE2)
    /* Shift */
    startdiscard += 96 - (endcolumni - startcolumni)*32;
    enddiscard += 96;
    diff_128 = (block_diff_snp_128_shift_hi)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					     endcolumni);
    diff_128 = clear_end_128(diff_128,enddiscard);
    diff_128 = set_start_128(diff_128,startdiscard);  /* puts 1 (mismatches) at start */
    _mm_store_si128((__m128i *) array,diff_128);

    bestscore = -100;
    for (i = 7; i >= 0; --i) {
      p = 3*array[i];
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
    }
    
    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */

#else
    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    p = 3*(diff_32 >> 16);
    bestscore = score_high[p];
    trimpos = offset - score_high[p+1];
    totalscore = score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    totalscore += score_high[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;
    query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;

    /* Single row */
    while (--endcolumni > startcolumni) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 >> 16);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
      
      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
      query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
    }

    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = set_start_32(diff_32,startdiscard);  /* puts 1 (mismatches) at start */

    p = 3*(diff_32 >> 16);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    totalscore += score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;
    
    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    /* totalscore += score_high[p+2]; */
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset -= 16; */
    
    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */
#endif

#if defined(USE_WRAP_TRIM) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    startdiscard += (startcolumni - endcolumni - 1)*32;
    enddiscard += 96;
    diff_128 = (block_diff_snp_128_wrap_hi)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					    endcolumni);
    diff_128 = clear_end_128(diff_128,enddiscard);
    diff_128 = set_start_128(diff_128,startdiscard);  /* puts 1 (mismatches) at start */
    _mm_store_si128((__m128i *) array,diff_128);

    bestscore = -100;
    for (i = 7; i >= 0; --i) {
      p = 3*array[i];
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
    }
    
    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */
#endif

  } else {
    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);

    p = 3*(diff_32 >> 16);
    bestscore = score_high[p];
    trimpos = offset - score_high[p+1];
    totalscore = score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    /* totalscore += score_high[p+2]; */
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;
    query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;

    /* End row */
    while (--endcolumni >= 0) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 >> 16);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore = score_high[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
      query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
    }
#ifdef HAVE_SSE2
    query_shifted -= QUERY_NEXTROW;
#endif
    ref_ptr -= GENOME_NEXTROW; alt_ptr -= GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ref_ptr >= startptr + 24) {
      diff_256 = (block_diff_snp_256)(&(query_shifted[-15]),&(ref_ptr[-15]),alt_ptr,plusp,genestrand,query_unk_mismatch_p);
      _mm256_store_si256((__m256i *) array,diff_256);

      for (i = 15; i >= 0; --i) {
	p = 3*array[i];
	if ((score = score_high[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset - score_high[p+1];
	}
	totalscore += score_high[p+2];
	debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset -= 16;
      }
      query_shifted -= 24; ref_ptr -= 24; alt_ptr -= 24;
    }
#endif

#ifdef HAVE_SSE2
    while (ref_ptr >= startptr + 12) {
      diff_128 = (block_diff_snp_128)(&(query_shifted[-3]),&(ref_ptr[-3]),alt_ptr,plusp,genestrand,query_unk_mismatch_p);
      _mm_store_si128((__m128i *) array,diff_128);

      for (i = 7; i >= 0; --i) {
	p = 3*array[i];
	if ((score = score_high[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset - score_high[p+1];
	}
	totalscore += score_high[p+2];
	debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset -= 16;
      }
      query_shifted -= 12; ref_ptr -= 12; alt_ptr -= 12;
    }
#else
    while (ref_ptr >= startptr + 12) {
      for (endcolumni = 3; endcolumni >= 0; --endcolumni) {
	diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

	p = 3*(diff_32 >> 16);
	if ((score = score_high[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset - score_high[p+1];
	}
	totalscore = score_high[p+2];
	debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset -= 16;
	
	p = 3*(diff_32 & 0x0000FFFF);
	if ((score = score_high[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset - score_high[p+1];
	}
	totalscore += score_high[p+2];
	debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset -= 16;
	query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
      }
      /* query_shifted -= QUERY_NEXTROW; */ ref_ptr -= GENOME_NEXTROW; alt_ptr -= GENOME_NEXTROW;
    }
#endif

    /* Start row */
    while (ref_ptr > startptr) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 >> 16);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore = score_high[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_high[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset - score_high[p+1];
      }
      totalscore += score_high[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset -= 16;
      query_shifted -= QUERY_NEXTCOL; ref_ptr -= GENOME_NEXTCOL; alt_ptr -= GENOME_NEXTCOL;
    }

    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = set_start_32(diff_32,startdiscard);  /* puts 1 (mismatches) at start */

    p = 3*(diff_32 >> 16);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    totalscore += score_high[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset -= 16;
    
    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_high[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset - score_high[p+1];
    }
    /* totalscore += score_high[p+2]; */
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset -= 16; */
    
    return (trimpos - 1);	/* trimpos-1 is on side of mismatch */
  }
}
#endif


#if 0
static int
trim_right_substring (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
		      bool plusp, int genestrand) {
  int startdiscard, enddiscard, offset;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *ptr, *endptr;
  Genomecomp_T *query_shifted;
  int nshift;
  int startcolumni, endcolumni;
  UINT4 diff_32;
#ifdef HAVE_SSE2
  __m128i diff_128;
  int i;
#endif
#ifdef HAVE_AVX2
  __m256i diff_256;
#endif
#ifdef HAVE_AVX2
  unsigned short array[16];
#elif defined(HAVE_SSE2)
  unsigned short array[8];
#endif

  int totalscore, bestscore, score;
  int trimpos;
  Genomecomp_T p;

  debug(
	printf("\n\n");
	printf("Genome (in trim_right_substring) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
    debug5(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard); /* puts 0 (matches) at start */
    diff_32 = set_end_32(diff_32,enddiscard);  /* puts 1 (mismatches) at end */

    p = 3*(diff_32 & 0x0000FFFF);
    bestscore = score_low[p];
    trimpos = offset + score_low[p+1];
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    /* totalscore += score_low[p+2]; */
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset += 16; */
    
    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_TRIM) && defined(HAVE_SSE2)
    /* Shift */
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_shift_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					 startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard); /* puts 0 (matches) at start */
    diff_128 = set_end_128(diff_128,enddiscard);  /* puts 1 (mismatches) at end */
    _mm_store_si128((__m128i *) array,diff_128);

    bestscore = -100;
    for (i = 0; i < 8; i++) {
      p = 3*array[i];
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
    }

    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */
#else
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    p = 3*(diff_32 & 0x0000FFFF);
    bestscore = score_low[p];
    trimpos = offset + score_low[p+1];
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    totalscore += score_low[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore = score_low[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      
      p = 3*(diff_32 >> 16);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = set_end_32(diff_32,enddiscard);  /* puts 1 (mismatches) at end */

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    /* totalscore += score_low[p+2]; */
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset += 16; */
    
    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */
#endif

#if defined(USE_WRAP_TRIM) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_128_wrap_lo)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p,
					startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard); /* puts 0 (matches) at start */
    diff_128 = set_end_128(diff_128,enddiscard);  /* puts 1 (mismatches) at end */
    _mm_store_si128((__m128i *) array,diff_128);

    bestscore = -100;
    for (i = 0; i < 8; i++) {
      p = 3*array[i];
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
    }

    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */
#endif

  } else {
    /* Start block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    p = 3*(diff_32 & 0x0000FFFF);
    bestscore = score_low[p];
    trimpos = offset + score_low[p+1];
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    totalscore += score_low[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;
    query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore = score_low[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      
      p = 3*(diff_32 >> 16);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ptr + 24 <= endptr) {
      diff_256 = (block_diff_256)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      _mm256_store_si256((__m256i *) array,diff_256);

      for (i = 0; i < 16; i++) {
	p = 3*array[i];
	if ((score = score_low[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset + score_low[p+1];
	}
	totalscore += score_low[p+2];
	debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset += 16;
      }
      query_shifted += 24; ptr += 24;
    }
#endif

#ifdef HAVE_SSE2
    while (ptr + 12 <= endptr) {
      diff_128 = (block_diff_128)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
      _mm_store_si128((__m128i *) array,diff_128);

      for (i = 0; i < 8; i++) {
	p = 3*array[i];
	if ((score = score_low[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset + score_low[p+1];
	}
	totalscore += score_low[p+2];
	debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset += 16;
      }
      query_shifted += 12; ptr += 12;
    }
#else
    while (ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

	p = 3*(diff_32 & 0x0000FFFF);
	if ((score = score_low[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset + score_low[p+1];
	}
	totalscore = score_low[p+2];
	debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset += 16;
	
	p = 3*(diff_32 >> 16);
	if ((score = score_low[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset + score_low[p+1];
	}
	totalscore += score_low[p+2];
	debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset += 16;
	query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
      }
      /* query_shifted += QUERY_NEXTROW; */ ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ptr < endptr) {
      diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore = score_low[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      
      p = 3*(diff_32 >> 16);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      query_shifted += QUERY_NEXTCOL; ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_32)(query_shifted,ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = set_end_32(diff_32,enddiscard);  /* puts 1 (mismatches) at end */

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    /* totalscore += score_low[p+2]; */
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset += 16; */
    
    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */
  }
}
#endif


#if 0
static int
trim_right_substring_snps (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
			   bool plusp, int genestrand) {
  int startdiscard, enddiscard, offset;
  Univcoord_T startblocki, endblocki, startblocki_32, endblocki_32;
  Genomecomp_T *ref_ptr, *alt_ptr, *endptr;
  Genomecomp_T *query_shifted;
  int nshift;
  int startcolumni, endcolumni;
  UINT4 diff_32;
#ifdef HAVE_SSE2
  __m128i diff_128;
  int i;
#endif
#ifdef HAVE_AVX2
  __m256i diff_256;
#endif
#ifdef HAVE_AVX2
  unsigned short array[16];
#elif defined(HAVE_SSE2)
  unsigned short array[8];
#endif

  int totalscore, bestscore, score;
  int trimpos;
  Genomecomp_T p;

  debug(
	printf("\n\n");
	printf("Genome (in trim_right_substring) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
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
  ref_ptr = &(ref_blocks[startblocki_32]);
  alt_ptr = &(snp_blocks[startblocki_32]);
  endptr = &(ref_blocks[endblocki_32]);

  if (endblocki_32 == startblocki_32) {
    /* Single block */
    debug5(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard); /* puts 0 (matches) at start */
    diff_32 = set_end_32(diff_32,enddiscard);  /* puts 1 (mismatches) at end */

    p = 3*(diff_32 & 0x0000FFFF);
    bestscore = score_low[p];
    trimpos = offset + score_low[p+1];
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    /* totalscore += score_low[p+2]; */
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset += 16; */
    
    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */

  } else if (endblocki == startblocki) {
#if defined(USE_SHIFT_TRIM) && defined(HAVE_SSE2)
    /* Shift */
    enddiscard += (endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_shift_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					     startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard); /* puts 0 (matches) at start */
    diff_128 = set_end_128(diff_128,enddiscard);  /* puts 1 (mismatches) at end */
    _mm_store_si128((__m128i *) array,diff_128);

    bestscore = -100;
    for (i = 0; i < 8; i++) {
      p = 3*array[i];
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
    }

    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */

#else
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    p = 3*(diff_32 & 0x0000FFFF);
    bestscore = score_low[p];
    trimpos = offset + score_low[p+1];
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    totalscore += score_low[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;
    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;

    /* Single row */
    while (++startcolumni < endcolumni) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore = score_low[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      
      p = 3*(diff_32 >> 16);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = set_end_32(diff_32,enddiscard);  /* puts 1 (mismatches) at end */

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    /* totalscore += score_low[p+2]; */
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset += 16; */
    
    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */
#endif

#if defined(USE_WRAP_TRIM) && defined(HAVE_SSSE3)
  } else if (endblocki == startblocki + 12 && endcolumni < startcolumni) {
    /* Wrap */
    enddiscard += (4 + endcolumni - startcolumni)*32;
    diff_128 = (block_diff_snp_128_wrap_lo)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p,
					    startcolumni);
    diff_128 = clear_start_128(diff_128,startdiscard); /* puts 0 (matches) at start */
    diff_128 = set_end_128(diff_128,enddiscard);  /* puts 1 (mismatches) at end */
    _mm_store_si128((__m128i *) array,diff_128);

    bestscore = -100;
    for (i = 0; i < 8; i++) {
      p = 3*array[i];
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
    }
#endif

    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */

  } else {
    /* Start block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);

    p = 3*(diff_32 & 0x0000FFFF);
    bestscore = score_low[p];
    trimpos = offset + score_low[p+1];
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    totalscore += score_low[p+2];
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;
    query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;

    /* Start row */
    while (++startcolumni < 4) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore = score_low[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      
      p = 3*(diff_32 >> 16);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }
#ifdef HAVE_SSE2
    query_shifted += QUERY_NEXTROW;
#endif
    ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;


    /* Middle rows */
#ifdef HAVE_AVX2
    while (ref_ptr + 24 <= endptr) {
      diff_256 = (block_diff_snp_256)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      _mm256_store_si256((__m256i *) array,diff_256);

      for (i = 0; i < 16; i++) {
	p = 3*array[i];
	if ((score = score_low[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset + score_low[p+1];
	}
	totalscore += score_low[p+2];
	debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset += 16;
      }
      query_shifted += 24; ref_ptr += 24; alt_ptr += 24;
    }
#endif

#ifdef HAVE_SSE2
    while (ref_ptr + 12 <= endptr) {
      diff_128 = (block_diff_snp_128)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
      _mm_store_si128((__m128i *) array,diff_128);

      for (i = 0; i < 8; i++) {
	p = 3*array[i];
	if ((score = score_low[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset + score_low[p+1];
	}
	totalscore += score_low[p+2];
	debug(printf("diff piece %d %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     i,array[i],score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset += 16;
      }
      query_shifted += 12; ref_ptr += 12; alt_ptr += 12;
    }
#else
    while (ref_ptr + 12 <= endptr) {
      for (startcolumni = 0; startcolumni < 4; startcolumni++) {
	diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

	p = 3*(diff_32 & 0x0000FFFF);
	if ((score = score_low[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset + score_low[p+1];
	}
	totalscore = score_low[p+2];
	debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset += 16;
	
	p = 3*(diff_32 >> 16);
	if ((score = score_low[p] + totalscore) > bestscore) {
	  bestscore = score;
	  trimpos = offset + score_low[p+1];
	}
	totalscore += score_low[p+2];
	debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		     diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
	offset += 16;
	query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
      }
      /* query_shifted += QUERY_NEXTROW; */ ref_ptr += GENOME_NEXTROW; alt_ptr += GENOME_NEXTROW;
    }
#endif

    /* End row */
    while (ref_ptr < endptr) {
      diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);

      p = 3*(diff_32 & 0x0000FFFF);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore = score_low[p+2];
      debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      
      p = 3*(diff_32 >> 16);
      if ((score = score_low[p] + totalscore) > bestscore) {
	bestscore = score;
	trimpos = offset + score_low[p+1];
      }
      totalscore += score_low[p+2];
      debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		   diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
      offset += 16;
      query_shifted += QUERY_NEXTCOL; ref_ptr += GENOME_NEXTCOL; alt_ptr += GENOME_NEXTCOL;
    }

    /* End block */
    diff_32 = (block_diff_snp_32)(query_shifted,alt_ptr,ref_ptr,plusp,genestrand,query_unk_mismatch_p);
    diff_32 = set_end_32(diff_32,enddiscard);  /* puts 1 (mismatches) at end */

    p = 3*(diff_32 & 0x0000FFFF);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    totalscore = score_low[p+2];
    debug(printf("diff low %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 & 0x0000FFFF,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    offset += 16;

    p = 3*(diff_32 >> 16);
    if ((score = score_low[p] + totalscore) > bestscore) {
      bestscore = score;
      trimpos = offset + score_low[p+1];
    }
    /* totalscore += score_low[p+2]; */
    debug(printf("diff high %04X => bestscore %d at pos %d, offset %d, trimpos %d, totalscore %d\n",
		 diff_32 >> 16,score_high[p],score_high[p+1],offset,trimpos,totalscore));
    /* offset += 16; */
    
    return (trimpos + 1);	/* trimpos+1 is on side of mismatch */
  }
}
#endif


#if 0
/* Not currently called by any procedure */
int
Genome_trim_left (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
		  bool plusp, int genestrand) {
#if 0
  if (dibasep) {
    /* Not implemented */
    return 0;
  }
#endif

  if (snp_blocks == NULL) {
    return trim_left_substring(query_compress,left,pos5,pos3,plusp,genestrand);
  } else {
    return trim_left_substring_snps(query_compress,left,pos5,pos3,plusp,genestrand);
  }
}
#endif

#if 0
/* Not currently called by any procedure */
int
Genome_trim_right (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
		   bool plusp, int genestrand) {
#if 0
  if (dibasep) {
    /* Not implemented */
    return 0;
  }
#endif

  if (snp_blocks == NULL) {
    return trim_right_substring(query_compress,left,pos5,pos3,plusp,genestrand);
  } else {
    return trim_right_substring_snps(query_compress,left,pos5,pos3,plusp,genestrand);
  }
}
#endif




