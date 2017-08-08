static char rcsid[] = "$Id: merge.c 205967 2017-05-04 00:49:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "merge.h"
#include "assert.h"
#include "mem.h"
#include "popcount.h"

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
/* #define KEY_MASK (~0U << 2) */

#define PYRAMID_SIZE 32
#define KEY_MASK (~0U << 5)

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


#ifdef DEBUG
#ifdef HAVE_SSE4_1
static void
print_vector (__m128i x, char *label) {
  unsigned int *s = (unsigned int *) &x;

  printf("%s: %u %u %u %u\n",label,s[0],s[1],s[2],s[3]);
  return;
}
#endif

#ifdef HAVE_AVX2
static void
print_vector_256 (__m256i x, char *label) {
  unsigned int *s = (unsigned int *) &x;

  printf("%s: %u %u %u %u %u %u %u %u\n",label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7]);
  return;
}
#endif

#ifdef HAVE_AVX512
static void
print_vector_512 (__m512i x, char *label) {
  unsigned int *s = (unsigned int *) &x;

  printf("%s: %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u\n",
	 label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],
	 s[8],s[9],s[10],s[11],s[12],s[13],s[14],s[15]);
  return;
}
#endif
#endif



#ifdef HAVE_SSE4_1
/* The min and max procedures require SSE4.1, which makes SSE4.1 the minimum requirement for SIMD-based merge */
static void
merge_4x4 (__m128i *__restrict__ vMergedA, __m128i *__restrict__ vMergedB, __m128i vA, __m128i vB) {
  __m128i vTmp, vMin, vMax;

  vMin = _mm_min_epu32(vA, vB);
  vMax = _mm_max_epu32(vA, vB);
  /* print_vector(vMin,"Min 1"); */
  /* print_vector(vMax,"Max 1"); */

  vTmp = _mm_alignr_epi8(vMin, vMin, 4); /* Rotate Min by 4 */
  vMin = _mm_min_epu32(vTmp, vMax);
  vMax = _mm_max_epu32(vTmp, vMax);
  /* print_vector(vTmp,"Tmp 2"); */
  /* print_vector(vMin,"Min 2"); */
  /* print_vector(vMax,"Max 2"); */

  vTmp = _mm_alignr_epi8(vMin, vMin, 4);
  vMin = _mm_min_epu32(vTmp, vMax);
  vMax = _mm_max_epu32(vTmp, vMax);
  /* print_vector(vTmp,"Tmp 3"); */
  /* print_vector(vMin,"Min 3"); */
  /* print_vector(vMax,"Max 3"); */

  vTmp = _mm_alignr_epi8(vMin, vMin, 4);
  vMin = _mm_min_epu32(vTmp, vMax);
  /* print_vector(vTmp,"Tmp 4"); */
  /* print_vector(vMin,"Min 4"); */

  *vMergedB = _mm_max_epu32(vTmp, vMax);
  *vMergedA = _mm_alignr_epi8(vMin, vMin, 4);

  return;
}


#ifndef HAVE_AVX2
static void
merge_8x8_network (__m128i *__restrict__ vMergedA, __m128i *__restrict__ vMergedB,
		   __m128i *__restrict__ vMergedC, __m128i *__restrict__ vMergedD,
		   __m128i vA0, __m128i vA1, __m128i vB0, __m128i vB1) {
  merge_4x4(&(*vMergedA),&(*vMergedB),vA0,vB0);
  merge_4x4(&(*vMergedC),&(*vMergedD),vA1,vB1);

  merge_4x4(&(*vMergedB),&(*vMergedC),*vMergedC,*vMergedB);
  return;
}
#endif
#endif


#ifdef HAVE_AVX2
/* The problem is that _mm256_alignr_epi8 rotates within 128-bit lanes */
/* So use _mm256_permutevar8x32_epi32, which shuffles across lanes */
static void
merge_8x8 (__m256i *__restrict__ vMergedA, __m256i *__restrict__ vMergedB, __m256i vA, __m256i vB) {
  __m256i vTmp, vMin, vMax;
  __m256i vRot;

  vRot = _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 0);

  /* print_vector_256(vA,"vA"); */
  /* print_vector_256(*vB,"vB"); */

  /* 1 */
  vMin = _mm256_min_epu32(vA, vB);
  vMax = _mm256_max_epu32(vA, vB);
  /* print_vector_256(vMin,"Min 1"); */
  /* print_vector_256(vMax,"Max 1"); */

  /* 2 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = _mm256_min_epu32(vTmp, vMax);
  vMax = _mm256_max_epu32(vTmp, vMax);

  /* 3 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = _mm256_min_epu32(vTmp, vMax);
  vMax = _mm256_max_epu32(vTmp, vMax);

  /* 4 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = _mm256_min_epu32(vTmp, vMax);
  vMax = _mm256_max_epu32(vTmp, vMax);

  /* 5 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = _mm256_min_epu32(vTmp, vMax);
  vMax = _mm256_max_epu32(vTmp, vMax);

  /* 6 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = _mm256_min_epu32(vTmp, vMax);
  vMax = _mm256_max_epu32(vTmp, vMax);

  /* 7 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = _mm256_min_epu32(vTmp, vMax);
  vMax = _mm256_max_epu32(vTmp, vMax);

  /* 8 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = _mm256_min_epu32(vTmp, vMax);
  /* print_vector_256(vTmp,"Tmp 8"); */
  /* print_vector_256(vMin,"Min 8"); */

  *vMergedB = _mm256_max_epu32(vTmp, vMax);
  *vMergedA = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  /* print_vector_256(*vMergedA,"vMergedA"); */
  /* print_vector_256(*vMergedB,"vMergedB"); */
  /* printf("\n"); */

  return;
}

#ifndef HAVE_AVX512
static void
merge_16x16_network (__m256i *__restrict__ vMergedA, __m256i *__restrict__ vMergedB,
		     __m256i *__restrict__ vMergedC, __m256i *__restrict__ vMergedD,
		     __m256i vA0, __m256i vA1, __m256i vB0, __m256i vB1) {
  merge_8x8(&(*vMergedA),&(*vMergedB),vA0,vB0);
  merge_8x8(&(*vMergedC),&(*vMergedD),vA1,vB1);

  merge_8x8(&(*vMergedB),&(*vMergedC),*vMergedC,*vMergedB);
  return;
}
#endif
#endif


#ifdef HAVE_AVX512
static void
merge_16x16 (__m512i *__restrict__ vMergedA, __m512i *__restrict__ vMergedB, __m512i vA, __m512i vB) {
  __m512i vTmp, vMin, vMax;
  __m512i vRot;
  int i;

  vRot = _mm512_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0);

  /* print_vector_512(vA,"vA"); */
  /* print_vector_512(vB,"vB"); */

  /* 1 */
  vMin = _mm512_min_epu32(vA, vB);
  vMax = _mm512_max_epu32(vA, vB);
  /* print_vector_512(vMin,"Min 1"); */
  /* print_vector_512(vMax,"Max 1"); */

  /* 2..15 */
  for (i = 0; i < 14; i++) {
    vTmp = _mm512_permutexvar_epi32(vRot, vMin); /* Rotate Min by ints */
    vMin = _mm512_min_epu32(vTmp, vMax);
    vMax = _mm512_max_epu32(vTmp, vMax);
    /* print_vector_512(vTmp,"Tmp 2"); */
    /* print_vector_512(vMin,"Min 2"); */
    /* print_vector_512(vMax,"Max 2"); */
  }

  /* 16 */
  vTmp = _mm512_permutexvar_epi32(vRot, vMin); /* Rotate Min by ints */
  vMin = _mm512_min_epu32(vTmp, vMax);
  /* print_vector_512(vTmp,"Tmp 16"); */
  /* print_vector_512(vMin,"Min 16"); */

  *vMergedB = _mm512_max_epu32(vTmp, vMax);
  *vMergedA = _mm512_permutexvar_epi32(vRot, vMin); /* Rotate Min by ints */
  /* print_vector_512(*vMergedA,"vMergedA"); */
  /* print_vector_512(*vMergedB,"vMergedB"); */
  /* printf("\n"); */

  return;
}

static void
merge_32x32_network (__m512i *__restrict__ vMergedA, __m512i *__restrict__ vMergedB,
		     __m512i *__restrict__ vMergedC, __m512i *__restrict__ vMergedD,
		     __m512i vA0, __m512i vA1, __m512i vB0, __m512i vB1) {
  merge_16x16(&(*vMergedA),&(*vMergedB),vA0,vB0);
  merge_16x16(&(*vMergedC),&(*vMergedD),vA1,vB1);

  merge_16x16(&(*vMergedB),&(*vMergedC),*vMergedC,*vMergedB);
  return;
}
#endif



/* Assumes padding to nearest 4 uints, and alignment to nearest 16 bytes */
/* If dest is NULL, then allocates and returns memory.  Otherwise, fills in at dest */
unsigned int *
Merge_uint4 (unsigned int *__restrict__ dest, unsigned int *__restrict__ A,
	     unsigned int *__restrict__ B, int nA, int nB) {
  unsigned int *C0, *C, *Aend, *Bend;
  unsigned int nextA, nextB;
  int nC;
#ifdef HAVE_AVX512
  __m512i vMerged512, vMerged512_0, vMerged512_1,
    vOld512, vNew512, vOld512_0, vOld512_1, vNew512_0, vNew512_1;
#endif
#ifdef HAVE_AVX2
  __m256i vMerged256, vMerged256_0, vMerged256_1,
    vOld256, vNew256, vOld256_0, vOld256_1, vNew256_0, vNew256_1;
#endif
#ifdef HAVE_SSE4_1
  __m128i vMerged128, vMerged128_0, vMerged128_1,
    vOld128, vNew128, vOld128_0, vOld128_1, vNew128_0, vNew128_1;
#endif


  if ((nC = nA + nB) == 0) {
    return (unsigned int *) NULL;
  } else if (dest) {
    C0 = C = dest;
  } else {
#if defined(HAVE_SSE4_1)
    C0 = C = (unsigned int *) MALLOC_ALIGN(nC * sizeof(unsigned int));
#else
    C0 = C = (unsigned int *) MALLOC(nC * sizeof(unsigned int));
#endif
  }

  Aend = &(A[nA]);
  Bend = &(B[nB]);

#ifdef HAVE_AVX512
  if (A < Aend - 32 && B < Bend - 32) {
    /* 32 ints = 1024 bits */
    if ((nextA = A[32]) < (nextB = B[32])) {
      vOld512_0 = _mm512_load_si512((__m512i *) B); B += 16;
      vOld512_1 = _mm512_load_si512((__m512i *) B); B += 16;
      vNew512_0 = _mm512_load_si512((__m512i *) A); A += 16;
      vNew512_1 = _mm512_load_si512((__m512i *) A); A += 16;
    } else {
      vOld512_0 = _mm512_load_si512((__m512i *) A); A += 16;
      vOld512_1 = _mm512_load_si512((__m512i *) A); A += 16;
      vNew512_0 = _mm512_load_si512((__m512i *) B); B += 16;
      vNew512_1 = _mm512_load_si512((__m512i *) B); B += 16;
    }
    merge_32x32_network(&vMerged512_0,&vMerged512_1,&vOld512_0,&vOld512_1,
			vOld512_0,vOld512_1,vNew512_0,vNew512_1);
    _mm512_stream_si512((__m512i *) C,vMerged512_0); C += 16;
    _mm512_stream_si512((__m512i *) C,vMerged512_1); C += 16;

    while (A < Aend - 32 && B < Bend - 32) {
      if (nextA < nextB) {
	vNew512_0 = _mm512_load_si512((__m512i *) A); A += 16;
	vNew512_1 = _mm512_load_si512((__m512i *) A); A += 16;
	nextA = *A;
      } else {
	vNew512_0 = _mm512_load_si512((__m512i *) B); B += 16;
	vNew512_1 = _mm512_load_si512((__m512i *) B); B += 16;
	nextB = *B;
      }
      merge_32x32_network(&vMerged512_0,&vMerged512_1,&vOld512_0,&vOld512_1,
			  vOld512_0,vOld512_1,vNew512_0,vNew512_1);
      _mm512_stream_si512((__m512i *) C,vMerged512_0); C += 16;
      _mm512_stream_si512((__m512i *) C,vMerged512_1); C += 16;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 16; _mm512_store_si512((__m512i *) B,vOld512_1);
      B -= 16; _mm512_store_si512((__m512i *) B,vOld512_0);
    } else {
      A -= 16; _mm512_store_si512((__m512i *) A,vOld512_1);
      A -= 16; _mm512_store_si512((__m512i *) A,vOld512_0);
    }
  }
#endif


#ifdef HAVE_AVX512
  if (A < Aend - 16 && B < Bend - 16) {
    /* 512 bits */
    if ((nextA = A[16]) < (nextB = B[16])) {
      vOld512 = _mm512_load_si512((__m512i *) B); B += 16;
      vNew512 = _mm512_load_si512((__m512i *) A); A += 16;
    } else {
      vOld512 = _mm512_load_si512((__m512i *) A); A += 16;
      vNew512 = _mm512_load_si512((__m512i *) B); B += 16;
    }
    merge_16x16(&vMerged512,&vOld512,vOld512,vNew512);
    _mm512_stream_si512((__m512i *) C,vMerged512); C += 16;

    while (A < Aend - 16 && B < Bend - 16) {
      if (nextA < nextB) {
	vNew512 = _mm512_load_si512((__m512i *) A); A += 16; nextA = *A;
      } else {
	vNew512 = _mm512_load_si512((__m512i *) B); B += 16; nextB = *B;
      }
      merge_16x16(&vMerged512,&vOld512,vOld512,vNew512);
      _mm512_stream_si512((__m512i *) C,vMerged512); C += 16;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 16; _mm512_store_si512((__m512i *) B,vOld512);
    } else {
      A -= 16; _mm512_store_si512((__m512i *) A,vOld512);
    }
  }

#elif defined(HAVE_AVX2)
  if (A < Aend - 16 && B < Bend - 16) {
    if ((nextA = A[16]) < (nextB = B[16])) {
      vOld256_0 = _mm256_load_si256((__m256i *) B); B += 8;
      vOld256_1 = _mm256_load_si256((__m256i *) B); B += 8;
      vNew256_0 = _mm256_load_si256((__m256i *) A); A += 8;
      vNew256_1 = _mm256_load_si256((__m256i *) A); A += 8;
    } else {
      vOld256_0 = _mm256_load_si256((__m256i *) A); A += 8;
      vOld256_1 = _mm256_load_si256((__m256i *) A); A += 8;
      vNew256_0 = _mm256_load_si256((__m256i *) B); B += 8;
      vNew256_1 = _mm256_load_si256((__m256i *) B); B += 8;
    }
    merge_16x16_network(&vMerged256_0,&vMerged256_1,&vOld256_0,&vOld256_1,
			vOld256_0,vOld256_1,vNew256_0,vNew256_1);
    _mm256_stream_si256((__m256i *) C,vMerged256_0); C += 8;
    _mm256_stream_si256((__m256i *) C,vMerged256_1); C += 8;

    while (A < Aend - 16 && B < Bend - 16) {
      if (nextA < nextB) {
	vNew256_0 = _mm256_load_si256((__m256i *) A); A += 8;
	vNew256_1 = _mm256_load_si256((__m256i *) A); A += 8;
	nextA = *A;
      } else {
	vNew256_0 = _mm256_load_si256((__m256i *) B); B += 8;
	vNew256_1 = _mm256_load_si256((__m256i *) B); B += 8;
	nextB = *B;
      }
      merge_16x16_network(&vMerged256_0,&vMerged256_1,&vOld256_0,&vOld256_1,
			  vOld256_0,vOld256_1,vNew256_0,vNew256_1);
      _mm256_stream_si256((__m256i *) C,vMerged256_0); C += 8;
      _mm256_stream_si256((__m256i *) C,vMerged256_1); C += 8;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 8; _mm256_store_si256((__m256i *) B,vOld256_1);
      B -= 8; _mm256_store_si256((__m256i *) B,vOld256_0);
    } else {
      A -= 8; _mm256_store_si256((__m256i *) A,vOld256_1);
      A -= 8; _mm256_store_si256((__m256i *) A,vOld256_0);
    }
  }
#endif


#ifdef HAVE_AVX2
  if (A < Aend - 8 && B < Bend - 8) {
    /* 256 bits */
    if ((nextA = A[8]) < (nextB = B[8])) {
      vOld256 = _mm256_load_si256((__m256i *) B); B += 8;
      vNew256 = _mm256_load_si256((__m256i *) A); A += 8;
    } else {
      vOld256 = _mm256_load_si256((__m256i *) A); A += 8;
      vNew256 = _mm256_load_si256((__m256i *) B); B += 8;
    }
    merge_8x8(&vMerged256,&vOld256,vOld256,vNew256);
    _mm256_stream_si256((__m256i *) C,vMerged256); C += 8;

    while (A < Aend - 8 && B < Bend - 8) {
      if (nextA < nextB) {
	vNew256 = _mm256_load_si256((__m256i *) A); A += 8; nextA = *A;
      } else {
	vNew256 = _mm256_load_si256((__m256i *) B); B += 8; nextB = *B;
      }
      merge_8x8(&vMerged256,&vOld256,vOld256,vNew256);
      _mm256_stream_si256((__m256i *) C,vMerged256); C += 8;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 8; _mm256_store_si256((__m256i *) B,vOld256);
    } else {
      A -= 8; _mm256_store_si256((__m256i *) A,vOld256);
    }
  }

#elif defined(HAVE_SSE4_1)
  if (A < Aend - 8 && B < Bend - 8) {
    if ((nextA = A[8]) < (nextB = B[8])) {
      vOld128_0 = _mm_load_si128((__m128i *) B); B += 4;
      vOld128_1 = _mm_load_si128((__m128i *) B); B += 4;
      vNew128_0 = _mm_load_si128((__m128i *) A); A += 4;
      vNew128_1 = _mm_load_si128((__m128i *) A); A += 4;
    } else {
      vOld128_0 = _mm_load_si128((__m128i *) A); A += 4;
      vOld128_1 = _mm_load_si128((__m128i *) A); A += 4;
      vNew128_0 = _mm_load_si128((__m128i *) B); B += 4;
      vNew128_1 = _mm_load_si128((__m128i *) B); B += 4;
    }
    merge_8x8_network(&vMerged128_0,&vMerged128_1,&vOld128_0,&vOld128_1,
		      vOld128_0,vOld128_1,vNew128_0,vNew128_1);
    _mm_stream_si128((__m128i *) C,vMerged128_0); C += 4;
    _mm_stream_si128((__m128i *) C,vMerged128_1); C += 4;

    while (A < Aend - 8 && B < Bend - 8) {
      if (nextA < nextB) {
	vNew128_0 = _mm_load_si128((__m128i *) A); A += 4;
	vNew128_1 = _mm_load_si128((__m128i *) A); A += 4;
	nextA = *A;
      } else {
	vNew128_0 = _mm_load_si128((__m128i *) B); B += 4;
	vNew128_1 = _mm_load_si128((__m128i *) B); B += 4;
	nextB = *B;
      }
      merge_8x8_network(&vMerged128_0,&vMerged128_1,&vOld128_0,&vOld128_1,
			vOld128_0,vOld128_1,vNew128_0,vNew128_1);
      _mm_stream_si128((__m128i *) C,vMerged128_0); C += 4;
      _mm_stream_si128((__m128i *) C,vMerged128_1); C += 4;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 4; _mm_store_si128((__m128i *) B,vOld128_1);
      B -= 4; _mm_store_si128((__m128i *) B,vOld128_0);
    } else {
      A -= 4; _mm_store_si128((__m128i *) A,vOld128_1);
      A -= 4; _mm_store_si128((__m128i *) A,vOld128_0);
    }
  }
#endif


#ifdef HAVE_SSE4_1
  if (A < Aend - 4 && B < Bend - 4) {
    /* 128 bits */
    if ((nextA = A[4]) < (nextB = B[4])) {
      vOld128 = _mm_load_si128((__m128i *) B); B += 4;
      vNew128 = _mm_load_si128((__m128i *) A); A += 4;
    } else {
      vOld128 = _mm_load_si128((__m128i *) A); A += 4;
      vNew128 = _mm_load_si128((__m128i *) B); B += 4;
    }
    merge_4x4(&vMerged128,&vOld128,vOld128,vNew128);
    _mm_stream_si128((__m128i *) C,vMerged128); C += 4;

    while (A < Aend - 4 && B < Bend - 4) {
      if (nextA < nextB) {
	vNew128 = _mm_load_si128((__m128i *) A); A += 4; nextA = *A;
      } else {
	vNew128 = _mm_load_si128((__m128i *) B); B += 4; nextB = *B;
      }
      merge_4x4(&vMerged128,&vOld128,vOld128,vNew128);
      _mm_stream_si128((__m128i *) C,vMerged128); C += 4;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 4; _mm_store_si128((__m128i *) B,vOld128);
    } else {
      A -= 4; _mm_store_si128((__m128i *) A,vOld128);
    }
  }
#endif

  /* Serial */
  while (A < Aend && B < Bend) {
    if (*A < *B) {
      *C++ = *A++;
    } else {
      *C++ = *B++;
    }
  }

  memcpy(C,A,(Aend - A) * sizeof(unsigned int));
  memcpy(C,B,(Bend - B) * sizeof(unsigned int));

  return C0;
}



#define PARENT(i) (i >> 1)
#define LEFT(i) (i << 1)
#define RIGHT(i) ((i << 1) | 1)

static int
pyramid_merge (unsigned int **heap, int nstreams, int heapsize, int *nelts,
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
      heap[PARENT(nodei)] = Merge_uint4(/*dest*/NULL,heap[nodei-1],heap[nodei],nelts[nodei-1],nelts[nodei]);
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


/* Assumes heapi < base put into LEFT(heapi) */
static int
pyramid_merge_prealloc (unsigned int **heap, unsigned int *curr_storage, unsigned int *prev_storage,
			int *nelts, int pyramid_start, int pyramid_end) {
  unsigned int *temp;
  int nodei;
  int nalloc;
#ifdef HAVE_SSE4_1
  int n;
#endif

  while (pyramid_end > pyramid_start) {
    debug2(printf("Merging level: %d..%d\n",pyramid_start,pyramid_end));
    nalloc = 0;

    nodei = pyramid_end;
    while (nodei >= pyramid_start) {
      debug2(printf("Merging nodes %d (%d elts) and %d (%d elts) => %d\n",
		   nodei-1,nelts[nodei-1],nodei,nelts[nodei],PARENT(nodei)));
      heap[PARENT(nodei)] = Merge_uint4(/*dest*/&(curr_storage[nalloc]),heap[nodei-1],heap[nodei],nelts[nodei-1],nelts[nodei]);
      CHECK_ALIGN(heap[PARENT(nodei)]);
      /* Have to align start of each entry curr_storage[nalloc], regardless of end padding */
#ifdef HAVE_SSE4_1
      n = nelts[PARENT(nodei)] = nelts[nodei-1] + nelts[nodei];
      nalloc += PAD_UINT4(n);
#else
      nalloc += (nelts[PARENT(nodei)] = nelts[nodei-1] + nelts[nodei]);
#endif
      debug2(printf("Created list %p of length %d at node %d\n",
		    heap[PARENT(nodei)],nelts[PARENT(nodei)],PARENT(nodei)));

#ifdef DEBUG2
      for (i = 0; i < nelts[PARENT(nodei)]; i++) {
	printf("%u\n",heap[PARENT(nodei)][i]);
      }
#endif

      /* Freeing memory one row at a time, so don't do it here */
      nodei -= 2;
    }

    /* Swap memory spaces */
    debug(printf("Swapping storage spaces\n"));
    temp = prev_storage;
    prev_storage = curr_storage;
    curr_storage = temp;

    /* Go up a level */
    pyramid_end = PARENT(pyramid_end);
    pyramid_start = PARENT(pyramid_start);
  }

  debug(printf("Returning ancestor %d\n\n",pyramid_start));
  return pyramid_start;
}


static UINT4 **
make_diagonals_heap (int *ncopied, int **nelts, List_T stream_list, Intlist_T streamsize_list, int nstreams) {
  UINT4 **heap, *stream;
  int heapsize, heapi, n;

#ifdef DEBUG
  int i;
#endif

  *ncopied = 0;
  heapsize = 2*nstreams - 1;

  heap = (UINT4 **) CALLOC((heapsize + 1),sizeof(UINT4 *));
  *nelts = (int *) CALLOC((heapsize + 1),sizeof(int));

  /* Process in reverse order, because stream_list is in reverse order of elts */
  heapi = heapsize;
  while (stream_list != NULL) {
    streamsize_list = Intlist_pop(streamsize_list,&n);
    (*nelts)[heapi] = n;

#if 0
    stream_list = List_pop(stream_list,(void *) &(heap[heapi])); /* already padded */
#else
    /* Copy to make the merging process non-destructive */
    heap[heapi] = MALLOC_ALIGN(n*sizeof(UINT4));
    stream_list = List_pop(stream_list,(void *) &stream);
    memcpy(heap[heapi],stream,n*sizeof(UINT4));
    *ncopied += 1;
#endif

    CHECK_ALIGN(heap[heapi]);
    debug(printf("Assigning node %d with %d elts",heapi,(*nelts)[heapi]));
#ifdef DEBUG
    for (i = 0; i < (*nelts)[heapi]; i++) {
      printf(" %u",heap[heapi][i]);
    }
#endif
    debug(printf("\n"));
    heapi--;
  }

  return heap;
}


UINT4 *
Merge_diagonals (int *nelts1, List_T stream_list, Intlist_T streamsize_list) {
  UINT4 *result, **heap, *stream;
  int *nelts;
  int nstreams, ncopied, heapi, heapsize, base, ancestori, pyramid_start, pyramid_end;
  int bits;
#ifdef DEBUG
  int i;
#endif


  if ((nstreams = List_length(stream_list)) == 0) {
    *nelts1 = 0;
    return (UINT4 *) NULL;

  } else if (nstreams == 1) {
    streamsize_list = Intlist_pop(streamsize_list,&(*nelts1));
    stream_list = List_pop(stream_list,(void *) &stream);
    result = MALLOC_ALIGN((*nelts1)*sizeof(UINT4));
    memcpy(result,stream,(*nelts1)*sizeof(UINT4));
    return result;

  } else {
    heapsize = 2*nstreams - 1;	/* also index of last node */
#ifdef HAVE_BUILTIN_CLZ
    bits = 31 - __builtin_clz((unsigned int) heapsize);
#elif defined(HAVE_ASM_BSR)
    asm("bsr %1,%0" : "=r"(bits) : "r"(heapsize));
#else
    bits = 31 - ((heapsize >> 16) ? clz_table[heapsize >> 16] : 16 + clz_table[heapsize]); 
#endif

    base = (1 << bits);
    heap = make_diagonals_heap(&ncopied,&nelts,stream_list,streamsize_list,nstreams);
    debug(printf("nstreams %d, heapsize %d, clz %d, bits %d, base %d\n",nstreams,heapsize,__builtin_clz(heapsize),bits,base));
  }

  /* Middle pyramids */
  while (base > PYRAMID_SIZE) {
    for (pyramid_start = 2*base - PYRAMID_SIZE, pyramid_end = 2*base - 1; pyramid_start >= base;
	 pyramid_start -= PYRAMID_SIZE, pyramid_end -= PYRAMID_SIZE) {
      debug(printf("diagonals: pyramid_start %d, pyramid_end %d, nstreams %d\n",pyramid_start,pyramid_end,nstreams));
      ancestori = pyramid_merge(heap,nstreams,heapsize,nelts,pyramid_start,pyramid_end);
    }
    base = ancestori;
  }

  /* Last pyramid */
  pyramid_start = base;
  pyramid_end = 2*base - 1;
  debug(printf("diagonals: pyramid_start %d, pyramid_end %d, nstreams %d\n",pyramid_start,pyramid_end,nstreams));
  /* base = */ pyramid_merge(heap,nstreams,heapsize,nelts,pyramid_start,pyramid_end);

  *nelts1 = nelts[1];
  result = heap[1];

  for (heapi = heapsize; heapi > heapsize - ncopied; heapi--) {
    FREE_ALIGN(heap[heapi]);
  }

  FREE(heap);
  FREE(nelts);

#ifdef DEBUG
  printf("Merge_diagonals returning result of length %d\n",*nelts1);
  for (i = 0; i < *nelts1; i++) {
    printf("%u\n",result[i]);
  }
#endif

  return result;
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

  /* Add PYRAMID_SIZE to handle partial pyramid */
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
    debug2(printf("NULL: Assigning node %d with %d elts (%p)",heapi,(*nelts)[heapi],record_heap[heapi]));

    for (i = 0; i < (*nelts)[heapi]; i++) {
      /* Process in forward order to keep records in order */
      all_records[k].diagonal = diagonals[i] + diagterm;
      all_records[k].querypos = querypos;
      record_heap[heapi][i] = &(all_records[k]);
      debug2(printf(" %u+%d",diagonals[i],querypos));
      k++;
    }
    debug2(printf("\n"));
  }
    
  for ( ; heapi >= nstreams; heapi--) {
    /* Move all information down to left child */
    basei = LEFT(heapi);
    stream_list = List_pop(stream_list,(void *) &diagonals); /* already padded */
    streamsize_list = Intlist_pop(streamsize_list,&((*nelts)[basei]));
    querypos_list = Intlist_pop(querypos_list,&querypos);
    diagterm_list = Intlist_pop(diagterm_list,&diagterm);
    record_heap[basei] = (Record_T *) MALLOC(((*nelts)[basei]) * sizeof(Record_T));
    debug2(printf("PART: Assigning node %d => %d with %d elts (%p)",heapi,basei,(*nelts)[basei],record_heap[basei]));

    for (i = 0; i < (*nelts)[basei]; i++) {
      /* Process in forward order to keep records in order */
      all_records[k].diagonal = diagonals[i] + diagterm;
      all_records[k].querypos = querypos;
      record_heap[basei][i] = &(all_records[k]);
      debug2(printf(" %u+%d",diagonals[i],querypos));
      k++;
    }
    debug2(printf("\n"));
  }

  for (heapi = heapsize; heapi >= base; heapi--) {
    /* Put all information into base row */
    stream_list = List_pop(stream_list,(void *) &diagonals); /* already padded */
    streamsize_list = Intlist_pop(streamsize_list,&((*nelts)[heapi]));
    querypos_list = Intlist_pop(querypos_list,&querypos);
    diagterm_list = Intlist_pop(diagterm_list,&diagterm);
    record_heap[heapi] = (Record_T *) MALLOC(((*nelts)[heapi]) * sizeof(Record_T));
    debug2(printf("FULL: Assigning node %d with %d elts (%p)",heapi,(*nelts)[heapi],record_heap[heapi]));

    for (i = 0; i < (*nelts)[heapi]; i++) {
      /* Process in forward order to keep records in order */
      all_records[k].diagonal = diagonals[i] + diagterm;
      all_records[k].querypos = querypos;
      record_heap[heapi][i] = &(all_records[k]);
      debug2(printf(" %u+%d",diagonals[i],querypos));
      k++;
    }
    debug2(printf("\n"));
  }

  return record_heap;
}



/* For initializing heap, there are three categories:
   base..(heapsize % PYRAMID_SIZE) + PYRAMID_SIZE: Fill bottom row
   straddling heapsize: Pull down some nodes to bottom row
   heapsize..(2*base - 1): Fill penultimate row */
Record_T *
Merge_records (int *nelts1, List_T stream_list, Intlist_T streamsize_list,
	       Intlist_T querypos_list, Intlist_T diagterm_list, 
	       struct Record_T *all_records) {
  Record_T *result, **record_heap, curr;
  UINT4 **key_heap, *prev_storage, *curr_storage;
  int ptrs[PYRAMID_SIZE];
  int *nelts, nalloc;
  int nstreams, heapsize, base, ancestori, pyramid_start, pyramid_end,
    node_start, node_end;
  int bits;
  int heapi, streami, i, j, k;

  debug2(printf("Entered Merge_records\n"));

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
    debug2(printf("nstreams %d, heapsize %d, base %d\n",nstreams,heapsize,base));
    record_heap = make_record_heap(&nelts,stream_list,streamsize_list,querypos_list,diagterm_list,
				   nstreams,base,all_records);
  }

  if (nstreams == 1) {
    *nelts1 = nelts[1];
    result = record_heap[1];

    FREE(nelts);
    FREE(record_heap);

#ifdef DEBUG2
    printf("Merge_records returning result of length %d\n",*nelts1);
    for (i = 0; i < *nelts1; i++) {
      printf("%u %d\n",result[i]->diagonal,result[i]->querypos);
    }
#endif

    return result;
  }


  key_heap = (UINT4 **) CALLOC(heapsize + PYRAMID_SIZE,sizeof(UINT4 *));

  /* Middle pyramids */
  while (base > PYRAMID_SIZE) {
    for (pyramid_start = 2*base - PYRAMID_SIZE, pyramid_end = 2*base - 1; pyramid_start >= base;
	 pyramid_start -= PYRAMID_SIZE, pyramid_end -= PYRAMID_SIZE) {
      debug2(printf("records: pyramid_start %d, pyramid_end %d, nstreams %d",pyramid_start,pyramid_end,nstreams));

      if (pyramid_start > heapsize) {
	node_start = PARENT(pyramid_start);
	node_end = PARENT(pyramid_end);
	debug2(printf(" => node_start %d, node_end %d\n",node_start,node_end));
      } else {
	node_start = pyramid_start;
	node_end = pyramid_end;
      }
      debug2(printf("\n"));

      /* Allocate memory for the pyramid */
      nalloc = 0;
      /* Have to align start of each entry prev_storage[nalloc] and curr_storage[nalloc], regardless of end padding */
#ifdef HAVE_SSE4_1
      for (heapi = node_start; heapi <= node_end; heapi++) {
	nalloc += PAD_UINT4(nelts[heapi]);
      }
      prev_storage = (UINT4 *) MALLOC_ALIGN(nalloc * sizeof(UINT4));
      curr_storage = (UINT4 *) MALLOC_ALIGN(nalloc * sizeof(UINT4));
#else
      for (heapi = node_start; heapi <= node_end; heapi++) {
        nalloc += nelts[heapi];
      }
      prev_storage = (UINT4 *) MALLOC(nalloc * sizeof(UINT4));
      curr_storage = (UINT4 *) MALLOC(nalloc * sizeof(UINT4));
#endif
      
      /* Convert structures to integers (key_heap) */
      nalloc = 0;
      for (heapi = node_start, streami = 0; heapi <= node_end; heapi++, streami++) {
	debug2(printf("Creating key node %d from %p\n",heapi,record_heap[heapi]));
	/* key_heap[heapi] = (UINT4 *) MALLOC((npadded + 1) * sizeof(UINT4)); */
	key_heap[heapi] = &(prev_storage[nalloc]);
	for (i = 0; i < nelts[heapi]; i++) {
	  key_heap[heapi][i] = (record_heap[heapi][i]->diagonal & KEY_MASK) + streami;
	}
        /* Had to align start of each entry prev_storage[nalloc], regardless of end padding */
#ifdef HAVE_SSE4_1
	nalloc += PAD_UINT4(nelts[heapi]);
#else
	nalloc += nelts[heapi];
#endif
      }

      ancestori = pyramid_merge_prealloc(key_heap,curr_storage,prev_storage,nelts,
                                	 node_start,node_end);

      /* Convert integers to structures */
      record_heap[ancestori] = (Record_T *) MALLOC(nelts[ancestori] * sizeof(Record_T));
      memset(ptrs,0,PYRAMID_SIZE*sizeof(int));
      k = 0;
      for (i = 0; i < nelts[ancestori]; i++) {
	streami = key_heap[ancestori][i] & ~KEY_MASK;
	record_heap[ancestori][k++] = record_heap[node_start + streami][ptrs[streami]++];
      }
      
      /* Free base heaps */
      for (heapi = node_start; heapi <= node_end; heapi++) {
	FREE(record_heap[heapi]);
      }
      
      /* Free key_heap storage */
      FREE_ALIGN(prev_storage);
      FREE_ALIGN(curr_storage);

    }
    base = ancestori;
  }

  /* Last pyramid */
  pyramid_start = base;
  pyramid_end = 2*base - 1;
  debug2(printf("records: pyramid_start %d, pyramid_end %d, nstreams %d\n",pyramid_start,pyramid_end,nstreams));

  /* Allocate memory for the pyramid */
  nalloc = 0;

  /* Have to align start of each entry prev_storage[nalloc] and curr_storage[nalloc], regardless of end padding */
#ifdef HAVE_SSE4_1
  for (heapi = pyramid_start; heapi <= pyramid_end; heapi++) {
    nalloc += PAD_UINT4(nelts[heapi]);
  }
  prev_storage = (UINT4 *) MALLOC_ALIGN(nalloc * sizeof(UINT4));
  curr_storage = (UINT4 *) MALLOC_ALIGN(nalloc * sizeof(UINT4));
#else
  for (heapi = pyramid_start; heapi <= pyramid_end; heapi++) {
    nalloc += nelts[heapi];
  }
  prev_storage = (UINT4 *) MALLOC(nalloc * sizeof(UINT4));
  curr_storage = (UINT4 *) MALLOC(nalloc * sizeof(UINT4));
#endif

  /* Convert structures to integers (key_heap) */
  nalloc = 0;
  for (heapi = pyramid_start, streami = 0; heapi <= pyramid_end; heapi++, streami++) {
    /* key_heap[heapi] = (UINT4 *) MALLOC((npadded + 1) * sizeof(UINT4)); */
    key_heap[heapi] = &(prev_storage[nalloc]);
    for (i = 0; i < nelts[heapi]; i++) {
      key_heap[heapi][i] = (record_heap[heapi][i]->diagonal & KEY_MASK) + streami;
    }
    /* Had to align start each entry prev_storage[nalloc], regardless of end padding */
#ifdef HAVE_SSE4_1
    nalloc += PAD_UINT4(nelts[heapi]);
#else
    nalloc += nelts[heapi];
#endif
  }

  ancestori = pyramid_merge_prealloc(key_heap,curr_storage,prev_storage,
				     nelts,pyramid_start,pyramid_end);
  /* ancestori should be 1 */

  /* Convert integers to structures */
  record_heap[ancestori] = (Record_T *) MALLOC(nelts[ancestori] * sizeof(Record_T));
  memset(ptrs,0,PYRAMID_SIZE*sizeof(int));
  k = 0;
  for (i = 0; i < nelts[ancestori]; i++) {
    streami = key_heap[ancestori][i] & ~KEY_MASK;
    record_heap[ancestori][k++] = record_heap[pyramid_start + streami][ptrs[streami]++];
  }

  /* Free base heaps (unless pyramid_start == 1, implying that base == 1) */
  for (heapi = pyramid_start, streami = 0; heapi <= pyramid_end; heapi++, streami++) {
    FREE(record_heap[pyramid_start + streami]);
  }

  /* Free key_heap storage */
  FREE_ALIGN(prev_storage);
  FREE_ALIGN(curr_storage);


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


  FREE(key_heap);
  FREE(nelts);
  FREE(record_heap);

#ifdef DEBUG2
  printf("Merge_records returning result of length %d\n",*nelts1);
  for (i = 0; i < *nelts1; i++) {
    printf("%u %d\n",result[i]->diagonal,result[i]->querypos);
  }
#endif

  return result;
}


