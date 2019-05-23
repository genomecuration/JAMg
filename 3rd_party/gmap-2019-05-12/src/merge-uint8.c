static char rcsid[] = "$Id: merge-uint8.c 219217 2019-05-12 22:24:49Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "merge-uint8.h"
#include "mem.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */


#if defined(HAVE_AVX2)
#include <immintrin.h>
#endif
#if defined(HAVE_AVX512)
#include <immintrin.h>
#endif


#ifdef DEBUG
static void
print_vector_256 (__m256i x, char *label) {
  UINT8 *s = (UINT8 *) &x;

  printf("%s: %llu %llu %llu %llu\n",label,s[0],s[1],s[2],s[3]);
  return;
}

static void
print_vector_512 (__m512i x, char *label) {
  UINT8 *s = (UINT8 *) &x;

  printf("%s: %llu %llu %llu %llu %llu %llu %llu %llu\n",
	   label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7]);
  return;
}
#endif


#ifdef HAVE_AVX512
/* The min and max procedures for epu64 require AVX512VL + AVX512F */
static void
merge_4x4 (__m256i *__restrict__ vMergedA, __m256i *__restrict__ vMergedB, __m256i vA, __m256i vB) {
  __m256i vTmp, vMin, vMax;

  vMin = _mm256_min_epu64(vA, vB);
  vMax = _mm256_max_epu64(vA, vB);
  /* print_vector(vMin,"Min 1"); */
  /* print_vector(vMax,"Max 1"); */

  vTmp = _mm256_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64-bits */
  vMin = _mm256_min_epu64(vTmp, vMax);
  vMax = _mm256_max_epu64(vTmp, vMax);
  /* print_vector(vTmp,"Tmp 2"); */
  /* print_vector(vMin,"Min 2"); */
  /* print_vector(vMax,"Max 2"); */

  vTmp = _mm256_alignr_epi64(vMin, vMin, 1);
  vMin = _mm256_min_epu64(vTmp, vMax);
  vMax = _mm256_max_epu64(vTmp, vMax);
  /* print_vector(vTmp,"Tmp 3"); */
  /* print_vector(vMin,"Min 3"); */
  /* print_vector(vMax,"Max 3"); */

  vTmp = _mm256_alignr_epi64(vMin, vMin, 1);
  vMin = _mm256_min_epu64(vTmp, vMax);
  /* print_vector(vTmp,"Tmp 4"); */
  /* print_vector(vMin,"Min 4"); */

  *vMergedB = _mm256_max_epu64(vTmp, vMax);
  *vMergedA = _mm256_alignr_epi64(vMin, vMin, 1);

  return;
}
#elif defined(HAVE_AVX2)
/* AVX has _mm256_min_pd and _mm256_max_pd.  If coordinates < 2^52,
   then floating-point arithmetic will see 0 for the sign and 0 for
   the exponent, and min/max results should be the same */
/* Another problem is that _mm256_alignr_epi8 rotates within 128-bit
   lanes.  So use _mm256_permutevar8x32_epi32, which shuffles across
   lanes */

static void
merge_4x4 (__m256i *__restrict__ vMergedA, __m256i *__restrict__ vMergedB, __m256i vA, __m256i vB) {
  __m256i vTmp, vMin, vMax;
  __m256i vRot;

  vRot = _mm256_setr_epi32(2, 3, 4, 5, 6, 7, 0, 1);

  /* print_vector_256(vA,"vA"); */
  /* print_vector_256(*vB,"vB"); */

  /* 1 */
  vMin = (__m256i) _mm256_min_pd((__m256d) vA, (__m256d) vB);
  vMax = (__m256i) _mm256_max_pd((__m256d) vA, (__m256d) vB);
  /* print_vector_256(vMin,"Min 1"); */
  /* print_vector_256(vMax,"Max 1"); */

  /* 2 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = (__m256i) _mm256_min_pd((__m256d) vTmp, (__m256d) vMax);
  vMax = (__m256i) _mm256_max_pd((__m256d) vTmp, (__m256d) vMax);

  /* 3 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = (__m256i) _mm256_min_pd((__m256d) vTmp, (__m256d) vMax);
  vMax = (__m256i) _mm256_max_pd((__m256d) vTmp, (__m256d) vMax);

  /* 4 */
  vTmp = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  vMin = (__m256i) _mm256_min_pd((__m256d) vTmp, (__m256d) vMax);

  *vMergedB = (__m256i) _mm256_max_pd((__m256d) vTmp, (__m256d) vMax);
  *vMergedA = _mm256_permutevar8x32_epi32(vMin, vRot); /* Rotate Min by ints */
  /* print_vector_256(*vMergedA,"vMergedA"); */
  /* print_vector_256(*vMergedB,"vMergedB"); */
  /* printf("\n"); */

  return;
}
#endif


#ifdef HAVE_AVX512
/* If we have AVX512, don't need the network method; just use merge_8x8 */
static void
merge_8x8 (__m512i *__restrict__ vMergedA, __m512i *__restrict__ vMergedB, __m512i vA, __m512i vB) {
  __m512i vTmp, vMin, vMax;


  /* print_vector_256(vA,"vA"); */
  /* print_vector_256(*vB,"vB"); */

  /* 1 */
  vMin = _mm512_min_epu64(vA, vB);
  vMax = _mm512_max_epu64(vA, vB);
  /* print_vector_256(vMin,"Min 1"); */
  /* print_vector_256(vMax,"Max 1"); */

  /* 2 */
  vTmp = _mm512_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64 bits */
  vMin = _mm512_min_epu64(vTmp, vMax);
  vMax = _mm512_max_epu64(vTmp, vMax);

  /* 3 */
  vTmp = _mm512_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64 bits */
  vMin = _mm512_min_epu64(vTmp, vMax);
  vMax = _mm512_max_epu64(vTmp, vMax);

  /* 4 */
  vTmp = _mm512_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64 bits */
  vMin = _mm512_min_epu64(vTmp, vMax);
  vMax = _mm512_max_epu64(vTmp, vMax);

  /* 5 */
  vTmp = _mm512_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64 bits */
  vMin = _mm512_min_epu64(vTmp, vMax);
  vMax = _mm512_max_epu64(vTmp, vMax);

  /* 6 */
  vTmp = _mm512_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64 bits */
  vMin = _mm512_min_epu64(vTmp, vMax);
  vMax = _mm512_max_epu64(vTmp, vMax);

  /* 7 */
  vTmp = _mm512_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64 bits */
  vMin = _mm512_min_epu64(vTmp, vMax);
  vMax = _mm512_max_epu64(vTmp, vMax);

  /* 8 */
  vTmp = _mm512_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64 bits */
  vMin = _mm512_min_epu64(vTmp, vMax);
  /* print_vector_256(vTmp,"Tmp 8"); */
  /* print_vector_256(vMin,"Min 8"); */

  *vMergedB = _mm512_max_epu64(vTmp, vMax);
  *vMergedA = _mm512_alignr_epi64(vMin, vMin, 1); /* Rotate Min by 64 bits */

  /* print_vector_256(*vMergedA,"vMergedA"); */
  /* print_vector_256(*vMergedB,"vMergedB"); */
  /* printf("\n"); */

  return;
}

#elif defined(HAVE_AVX2)
static void
merge_8x8_network (__m256i *__restrict__ vMergedA, __m256i *__restrict__ vMergedB,
		   __m256i *__restrict__ vMergedC, __m256i *__restrict__ vMergedD,
		   __m256i vA0, __m256i vA1, __m256i vB0, __m256i vB1) {
  merge_4x4(&(*vMergedA),&(*vMergedB),vA0,vB0);
  merge_4x4(&(*vMergedC),&(*vMergedD),vA1,vB1);

  merge_4x4(&(*vMergedB),&(*vMergedC),*vMergedC,*vMergedB);
  return;
}
#endif



#ifdef HAVE_AVX512
static void
merge_16x16_network (__m512i *__restrict__ vMergedA, __m512i *__restrict__ vMergedB,
		     __m512i *__restrict__ vMergedC, __m512i *__restrict__ vMergedD,
		     __m512i vA0, __m512i vA1, __m512i vB0, __m512i vB1) {
  merge_8x8(&(*vMergedA),&(*vMergedB),vA0,vB0);
  merge_8x8(&(*vMergedC),&(*vMergedD),vA1,vB1);

  merge_8x8(&(*vMergedB),&(*vMergedC),*vMergedC,*vMergedB);
  return;
}
#endif



#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
/* Assumes padding to nearest 4 uints, and alignment to nearest 16 bytes */
/* If dest is NULL, then allocates and returns memory.  Otherwise, fills in at dest */
UINT8 *
Merge_uint8 (UINT8 *__restrict__ dest, UINT8 *__restrict__ A,
	     UINT8 *__restrict__ B, int nA, int nB) {
  UINT8 *C0, *C, *Aend, *Bend;
  UINT8 nextA, nextB;
  int nC;
#ifdef HAVE_AVX512
  __m512i vMerged512, vMerged512_0, vMerged512_1,
    vOld512, vNew512, vOld512_0, vOld512_1, vNew512_0, vNew512_1;
#endif
  __m256i vMerged256, vMerged256_0, vMerged256_1,
    vOld256, vNew256, vOld256_0, vOld256_1, vNew256_0, vNew256_1;


  if ((nC = nA + nB) == 0) {
    return (UINT8 *) NULL;
  } else if (dest) {
    C0 = C = dest;
  } else {
    C0 = C = (UINT8 *) MALLOC_ALIGN(nC * sizeof(UINT8));
  }

  Aend = &(A[nA]);
  Bend = &(B[nB]);


#ifdef HAVE_AVX512
  if (A < Aend - 16 && B < Bend - 16) {
    /* 16 UINT8s = 1024 bits */
    if ((nextA = A[16]) < (nextB = B[16])) {
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
    merge_16x16_network(&vMerged512_0,&vMerged512_1,&vOld512_0,&vOld512_1,
			vOld512_0,vOld512_1,vNew512_0,vNew512_1);
    _mm512_stream_si512((__m512i *) C,vMerged512_0); C += 16;
    _mm512_stream_si512((__m512i *) C,vMerged512_1); C += 16;

    while (A < Aend - 16 && B < Bend - 16) {
      if (nextA < nextB) {
	vNew512_0 = _mm512_load_si512((__m512i *) A); A += 16;
	vNew512_1 = _mm512_load_si512((__m512i *) A); A += 16;
	nextA = *A;
      } else {
	vNew512_0 = _mm512_load_si512((__m512i *) B); B += 16;
	vNew512_1 = _mm512_load_si512((__m512i *) B); B += 16;
	nextB = *B;
      }
      merge_16x16_network(&vMerged512_0,&vMerged512_1,&vOld512_0,&vOld512_1,
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
  if (A < Aend - 8 && B < Bend - 8) {
    /* 512 bits */
    if ((nextA = A[8]) < (nextB = B[8])) {
      vOld512 = _mm512_load_si512((__m512i *) B); B += 8;
      vNew512 = _mm512_load_si512((__m512i *) A); A += 8;
    } else {
      vOld512 = _mm512_load_si512((__m512i *) A); A += 8;
      vNew512 = _mm512_load_si512((__m512i *) B); B += 8;
    }
    merge_8x8(&vMerged512,&vOld512,vOld512,vNew512);
    _mm512_stream_si512((__m512i *) C,vMerged512); C += 8;

    while (A < Aend - 8 && B < Bend - 8) {
      if (nextA < nextB) {
	vNew512 = _mm512_load_si512((__m512i *) A); A += 8; nextA = *A;
      } else {
	vNew512 = _mm512_load_si512((__m512i *) B); B += 8; nextB = *B;
      }
      merge_8x8(&vMerged512,&vOld512,vOld512,vNew512);
      _mm512_stream_si512((__m512i *) C,vMerged512); C += 8;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 8; _mm512_store_si512((__m512i *) B,vOld512);
    } else {
      A -= 8; _mm512_store_si512((__m512i *) A,vOld512);
    }
  }

#elif defined(HAVE_AVX2)
  if (A < Aend - 8 && B < Bend - 8) {
    if ((nextA = A[8]) < (nextB = B[8])) {
      vOld256_0 = _mm256_load_si256((__m256i *) B); B += 4;
      vOld256_1 = _mm256_load_si256((__m256i *) B); B += 4;
      vNew256_0 = _mm256_load_si256((__m256i *) A); A += 4;
      vNew256_1 = _mm256_load_si256((__m256i *) A); A += 4;
    } else {
      vOld256_0 = _mm256_load_si256((__m256i *) A); A += 4;
      vOld256_1 = _mm256_load_si256((__m256i *) A); A += 4;
      vNew256_0 = _mm256_load_si256((__m256i *) B); B += 4;
      vNew256_1 = _mm256_load_si256((__m256i *) B); B += 4;
    }
    merge_8x8_network(&vMerged256_0,&vMerged256_1,&vOld256_0,&vOld256_1,
		      vOld256_0,vOld256_1,vNew256_0,vNew256_1);
    _mm256_stream_si256((__m256i *) C,vMerged256_0); C += 4;
    _mm256_stream_si256((__m256i *) C,vMerged256_1); C += 4;

    while (A < Aend - 8 && B < Bend - 8) {
      if (nextA < nextB) {
	vNew256_0 = _mm256_load_si256((__m256i *) A); A += 4;
	vNew256_1 = _mm256_load_si256((__m256i *) A); A += 4;
	nextA = *A;
      } else {
	vNew256_0 = _mm256_load_si256((__m256i *) B); B += 4;
	vNew256_1 = _mm256_load_si256((__m256i *) B); B += 4;
	nextB = *B;
      }
      merge_8x8_network(&vMerged256_0,&vMerged256_1,&vOld256_0,&vOld256_1,
			vOld256_0,vOld256_1,vNew256_0,vNew256_1);
      _mm256_stream_si256((__m256i *) C,vMerged256_0); C += 4;
      _mm256_stream_si256((__m256i *) C,vMerged256_1); C += 4;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 4; _mm256_store_si256((__m256i *) B,vOld256_1);
      B -= 4; _mm256_store_si256((__m256i *) B,vOld256_0);
    } else {
      A -= 4; _mm256_store_si256((__m256i *) A,vOld256_1);
      A -= 4; _mm256_store_si256((__m256i *) A,vOld256_0);
    }
  }
#endif

#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
  if (A < Aend - 4 && B < Bend - 4) {
    /* 256 bits */
    if ((nextA = A[4]) < (nextB = B[4])) {
      vOld256 = _mm256_load_si256((__m256i *) B); B += 4;
      vNew256 = _mm256_load_si256((__m256i *) A); A += 4;
    } else {
      vOld256 = _mm256_load_si256((__m256i *) A); A += 4;
      vNew256 = _mm256_load_si256((__m256i *) B); B += 4;
    }
    merge_4x4(&vMerged256,&vOld256,vOld256,vNew256);
    _mm256_stream_si256((__m256i *) C,vMerged256); C += 4;

    while (A < Aend - 4 && B < Bend - 4) {
      if (nextA < nextB) {
	vNew256 = _mm256_load_si256((__m256i *) A); A += 4; nextA = *A;
      } else {
	vNew256 = _mm256_load_si256((__m256i *) B); B += 4; nextB = *B;
      }
      merge_4x4(&vMerged256,&vOld256,vOld256,vNew256);
      _mm256_stream_si256((__m256i *) C,vMerged256); C += 4;
    }

    /* Re-insert before largest element */
    if (nextA < nextB) {
      B -= 4; _mm256_store_si256((__m256i *) B,vOld256);
    } else {
      A -= 4; _mm256_store_si256((__m256i *) A,vOld256);
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

  memcpy(C,A,(Aend - A) * sizeof(UINT8));
  memcpy(C,B,(Bend - B) * sizeof(UINT8));

  return C0;
}
#endif



