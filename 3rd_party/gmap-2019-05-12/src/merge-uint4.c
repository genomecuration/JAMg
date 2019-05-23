static char rcsid[] = "$Id: merge-uint4.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "merge-uint4.h"
#include "mem.h"

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


/* Turning off, because it appears to be slower */
/* #define USE_AVX2 1 */


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


#if !defined(HAVE_AVX2) || !defined(USE_AVX2)
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


#if defined(HAVE_AVX2) && defined(USE_AVX2)
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
#if defined(HAVE_AVX2) && defined(USE_AVX2)
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
    debug_align(printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
		       C,C,__FILE__,__LINE__));
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

#elif defined(HAVE_AVX2) && defined(USE_AVX2)
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


#if defined(HAVE_AVX2) && defined(USE_AVX2)
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


