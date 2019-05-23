static char rcsid[] = "$Id: epu16-bitpack64-readtwo.c 218164 2019-01-17 06:09:06Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "epu16-bitpack64-readtwo.h"

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#elif defined(HAVE_SSE2)
#include <emmintrin.h>
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#ifdef DEBUG
#include "epu16-bitpack64-read.h"
#endif


#if 0
#ifdef HAVE_SSE2
/* For debugging */
static void
print_vector_hex (__m128i x) {
#if 1
  UINT2 *s = (UINT2 *) &x;
#else
  UINT2 s[8];
  _mm_store_si128((__m128i *) s,x);
#endif

  printf("%04X %04X %04X %04X %04X %04X %04X %04X\n",
	 s[7],s[6],s[5],s[4],s[3],s[2],s[1],s[0]);
  return;
}

static void
print_vector (__m128i x) {
  UINT2 *s = (UINT2 *) &x;

  printf("%u %u %u %u %u %u %u %u\n",
	 s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7]);
  return;
}
#endif
#endif


static void
unpack_00_fwd_0_0 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i zero = _mm_set1_epi16(0U);

  debug(printf("Entering unpack_00_fwd_0_0\n"));

  _mm_store_si128(out++, zero);
  _mm_store_si128(out++, zero);

  return;
}



static void
unpack_02_fwd_1_2  (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  debug(printf("Entering unpack_02_fwd_1_2\n"));

    /* 1 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask2);
  _mm_store_si128(out++, OutReg);


  /* 2 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_02_fwd_2_3  (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  debug(printf("Entering unpack_02_fwd_2_3\n"));

  /* 2 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask2);
  _mm_store_si128(out++, OutReg);


  /* 3 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_02_fwd_3_4  (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  debug(printf("Entering unpack_02_fwd_3_4\n"));

  /* 3 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask2);
  _mm_store_si128(out++, OutReg);


  /* 4 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,6) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_02_fwd_4_5  (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  debug(printf("Entering unpack_02_fwd_4_5\n"));

  /* 4 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,6) , mask2);
  _mm_store_si128(out++, OutReg);


  /* 5 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_02_fwd_5_6  (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  debug(printf("Entering unpack_02_fwd_5_6\n"));

  /* 5 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask2);
  _mm_store_si128(out++, OutReg);


  /* 6 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,10) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_02_fwd_6_7  (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  debug(printf("Entering unpack_02_fwd_6_7\n"));

  /* 6 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,10) , mask2);
  _mm_store_si128(out++, OutReg);


  /* 7 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,12) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_02_fwd_7_8  (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  debug(printf("Entering unpack_02_fwd_7_8\n"));

  /* 7 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,12) , mask2);
  _mm_store_si128(out++, OutReg);


  /* 8 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,14) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_02_fwd_8_1  (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  debug(printf("Entering unpack_02_fwd_8_1\n"));

  /* 8 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,14) , mask2);
  _mm_store_si128(out++, OutReg);


  /* 1 */
  OutReg =  _mm_and_si128( InReg , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}




static void
unpack_04_fwd_1_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  debug(printf("Entering unpack_04_fwd_1_2\n"));

  /* 1 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask4);
  _mm_store_si128(out++, OutReg);


  /* 2 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_2_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  debug(printf("Entering unpack_04_fwd_2_3\n"));

  /* 2 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask4);
  _mm_store_si128(out++, OutReg);


  /* 3 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_3_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  debug(printf("Entering unpack_04_fwd_3_4\n"));

  /* 3 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask4);
  _mm_store_si128(out++, OutReg);


  /* 4 */
  OutReg =  _mm_srli_epi16(InReg,12);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_4_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  debug(printf("Entering unpack_04_fwd_4_5\n"));

  /* 4 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_srli_epi16(InReg,12);
  _mm_store_si128(out++, OutReg);


  /* 5 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128( InReg , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_5_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  debug(printf("Entering unpack_04_fwd_5_6\n"));

  /* 5 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128( InReg , mask4);
  _mm_store_si128(out++, OutReg);


  /* 6 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_04_fwd_6_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  debug(printf("Entering unpack_04_fwd_6_7\n"));

  /* 6 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask4);
  _mm_store_si128(out++, OutReg);


  /* 7 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_7_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  debug(printf("Entering unpack_04_fwd_7_8\n"));

  /* 7 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask4);
  _mm_store_si128(out++, OutReg);


  /* 8 */
  OutReg =  _mm_srli_epi16(InReg,12);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_04_fwd_8_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  debug(printf("Entering unpack_04_fwd_8_1\n"));

  /* 1 first */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask4);
  _mm_store_si128(&(out[1]), OutReg);


  /* 8 second */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_srli_epi16(InReg,12);
  _mm_store_si128(out++, OutReg);

  return;
}



static void
unpack_06_fwd_1_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  debug(printf("Entering unpack_06_fwd_1_2\n"));

  /* 1 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask6);
  _mm_store_si128(out++, OutReg);


  /* 2 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,6) , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_2_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  debug(printf("Entering unpack_06_fwd_2_3\n"));

  /* 2 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,6) , mask6);
  _mm_store_si128(out++, OutReg);


  /* 3 */
  OutReg =  _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 6 - 2), mask6));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_3_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  debug(printf("Entering unpack_06_fwd_3_4\n"));

  /* 3 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 6 - 2), mask6));
  _mm_store_si128(out++, OutReg);


  /* 4 */
  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_06_fwd_4_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  debug(printf("Entering unpack_06_fwd_4_5\n"));

  /* 4 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask6);
  _mm_store_si128(out++, OutReg);


  /* 5 */
  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_06_fwd_5_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  debug(printf("Entering unpack_06_fwd_5_6\n"));

  /* 5 */
  InReg = _mm_load_si128(++in);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask6);
  _mm_store_si128(out++, OutReg);


  /* 6 */
  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 6 - 4), mask6));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_6_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  debug(printf("Entering unpack_06_fwd_6_7\n"));

  /* 6 */
  InReg = _mm_load_si128(++in);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 6 - 4), mask6));
  _mm_store_si128(out++, OutReg);


  /* 7 */
  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_7_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  debug(printf("Entering unpack_06_fwd_7_8\n"));

  /* 7 */
  in += 2;
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask6);
  _mm_store_si128(out++, OutReg);


  /* 8 */
  OutReg =   _mm_srli_epi16(InReg,10) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_06_fwd_8_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  debug(printf("Entering unpack_06_fwd_8_1\n"));

  /* 1 first */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask6);
  _mm_store_si128(&(out[1]), OutReg);


  /* 8 second */
  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  _mm_store_si128(out++, OutReg);

  return;
}



static void
unpack_08_fwd_1_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  debug(printf("Entering unpack_08_fwd_1_2\n"));

  /* 1 */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);


  /* 2 */
  OutReg =   _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_08_fwd_2_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  debug(printf("Entering unpack_08_fwd_2_3\n"));

  /* 2 */
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);


  /* 3 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_3_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  debug(printf("Entering unpack_08_fwd_3_4\n"));

  /* 3 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);


  /* 4 */
  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_4_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  debug(printf("Entering unpack_08_fwd_4_5\n"));

  /* 4 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);


  /* 5 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_5_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  debug(printf("Entering unpack_08_fwd_5_6\n"));

  /* 5 */
  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  /* 6 */
  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_6_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  debug(printf("Entering unpack_08_fwd_6_7\n"));

  /* 6 */
  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);


  /* 7 */
  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_7_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  debug(printf("Entering unpack_08_fwd_7_8\n"));

  /* 7 */
  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);


  /* 8 */
  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_08_fwd_8_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  debug(printf("Entering unpack_08_fwd_8_1\n"));

  /* 1 first */
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(&(out[1]), OutReg);


  /* 8 second */
  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_1_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  debug(printf("Entering unpack_10_fwd_1_2\n"));

  /* 1 */
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask10);
  _mm_store_si128(out++, OutReg);


  /* 2 */
  OutReg =   _mm_srli_epi16(InReg,10) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 4), mask10));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_2_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  debug(printf("Entering unpack_10_fwd_2_3\n"));

  /* 2 */
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 4), mask10));
  _mm_store_si128(out++, OutReg);


  /* 3 */
  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask10);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_3_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  debug(printf("Entering unpack_10_fwd_3_4\n"));

  /* 3 */
  InReg = _mm_load_si128(++in);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask10);
  _mm_store_si128(out++, OutReg);


  /* 4 */
  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 8), mask10));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_4_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  debug(printf("Entering unpack_10_fwd_4_5\n"));

  /* 4 */
  InReg = _mm_load_si128(++in);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 8), mask10));
  _mm_store_si128(out++, OutReg);


  /* 5 */
  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 2), mask10));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_5_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  debug(printf("Entering unpack_10_fwd_5_6\n"));

  /* 5 */
  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 2), mask10));
  _mm_store_si128(out++, OutReg);


  /* 6 */
  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask10);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_6_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  debug(printf("Entering unpack_10_fwd_6_7\n"));

  /* 6 */
  in += 3;
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask10);
  _mm_store_si128(out++, OutReg);


  /* 7 */
  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 6), mask10));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_7_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  debug(printf("Entering unpack_10_fwd_7_8\n"));

  /* 7 */
  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 6), mask10));
  _mm_store_si128(out++, OutReg);


  /* 8 */
  OutReg =   _mm_srli_epi16(InReg,6) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_8_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  debug(printf("Entering unpack_10_fwd_8_1\n"));

  /* 1 first */
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask10);
  _mm_store_si128(&(out[1]), OutReg);


  /* 8 second */
  in += 4;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,6) ;
  _mm_store_si128(out++, OutReg);

  return;
}



static void
unpack_12_fwd_1_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  debug(printf("Entering unpack_12_fwd_1_2\n"));

  /* 1 */
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask12);
  _mm_store_si128(out++, OutReg);


  /* 2 */
  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 8), mask12));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_12_fwd_2_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  debug(printf("Entering unpack_12_fwd_2_3\n"));

  /* 2 */
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 8), mask12));
  _mm_store_si128(out++, OutReg);


  /* 3 */
  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 4), mask12));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_12_fwd_3_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  debug(printf("Entering unpack_12_fwd_3_4\n"));

  /* 3 */
  InReg = _mm_load_si128(++in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 4), mask12));
  _mm_store_si128(out++, OutReg);


  /* 4 */
  OutReg =   _mm_srli_epi16(InReg,4) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_12_fwd_4_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  debug(printf("Entering unpack_12_fwd_4_5\n"));

  /* 4 */
  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  _mm_store_si128(out++, OutReg);


  /* 5 */
  InReg = _mm_load_si128(++in);

  OutReg = _mm_and_si128( InReg , mask12);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_12_fwd_5_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  debug(printf("Entering unpack_12_fwd_5_6\n"));

  /* 5 */
  in += 3;
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask12);
  _mm_store_si128(out++, OutReg);


  /* 6 */
  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 8), mask12));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_12_fwd_6_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  debug(printf("Entering unpack_12_fwd_6_7\n"));

  /* 6 */
  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 8), mask12));
  _mm_store_si128(out++, OutReg);


  /* 7 */
  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 4), mask12));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_12_fwd_7_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  debug(printf("Entering unpack_12_fwd_7_8\n"));

  /* 7 */
  in += 4;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 4), mask12));
  _mm_store_si128(out++, OutReg);


  /* 8 */
  OutReg =   _mm_srli_epi16(InReg,4) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_12_fwd_8_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  debug(printf("Entering unpack_12_fwd_8_1\n"));

  /* 1 first */
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask12);
  _mm_store_si128(&(out[1]), OutReg);


  /* 8 second */
  in += 5;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  _mm_store_si128(out++, OutReg);

  return;
}



static void
unpack_14_fwd_1_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  debug(printf("Entering unpack_14_fwd_1_2\n"));

  /* 1 */
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask14);
  _mm_store_si128(out++, OutReg);


  /* 2 */
  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 12), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_14_fwd_2_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  debug(printf("Entering unpack_14_fwd_2_3\n"));

  /* 2 */
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 12), mask14));
  _mm_store_si128(out++, OutReg);


  /* 3 */
  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 10), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_14_fwd_3_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  debug(printf("Entering unpack_14_fwd_3_4\n"));

  /* 3 */
  InReg = _mm_load_si128(++in);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 10), mask14));
  _mm_store_si128(out++, OutReg);


  /* 4 */
  OutReg =   _mm_srli_epi16(InReg,10) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 8), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_14_fwd_4_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  debug(printf("Entering unpack_14_fwd_4_5\n"));

  /* 4 */
  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 8), mask14));
  _mm_store_si128(out++, OutReg);


  /* 5 */
  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 6), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_14_fwd_5_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  debug(printf("Entering unpack_14_fwd_5_6\n"));

  /* 5 */
  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 6), mask14));
  _mm_store_si128(out++, OutReg);


  /* 6 */
  OutReg =   _mm_srli_epi16(InReg,6) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 4), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_14_fwd_6_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  debug(printf("Entering unpack_14_fwd_6_7\n"));

  /* 6 */
  in += 4;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,6) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 4), mask14));
  _mm_store_si128(out++, OutReg);


  /* 7 */
  OutReg =   _mm_srli_epi16(InReg,4) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 2), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_14_fwd_7_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  debug(printf("Entering unpack_14_fwd_7_8\n"));

  /* 7 */
  in += 5;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 2), mask14));
  _mm_store_si128(out++, OutReg);


  /* 8 */
  OutReg =   _mm_srli_epi16(InReg,2) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_14_fwd_8_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  debug(printf("Entering unpack_14_fwd_8_1\n"));

  /* 1 first */
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask14);
  _mm_store_si128(&(out[1]), OutReg);


  /* 8 second */
  in += 6;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,2) ;
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_16_fwd_1_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  debug(printf("Entering unpack_16_fwd_1_2\n"));

  /* 1 */
  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  /* 2 */
  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_2_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  debug(printf("Entering unpack_16_fwd_2_3\n"));

  /* 2 */
  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  /* 3 */
  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_3_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  debug(printf("Entering unpack_16_fwd_3_4\n"));

  /* 3 */
  in += 2;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);


  /* 4 */
  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  return;
}



static void
unpack_16_fwd_4_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  debug(printf("Entering unpack_16_fwd_4_5\n"));

  /* 4 */
  in += 3;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);


  /* 5 */
  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_5_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  debug(printf("Entering unpack_16_fwd_5_6\n"));

  /* 5 */
  in += 4;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);


  /* 6 */
  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_6_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  debug(printf("Entering unpack_16_fwd_6_7\n"));

  /* 6 */
  in += 5;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);


  /* 7 */
  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_7_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  debug(printf("Entering unpack_16_fwd_7_8\n"));

  /* 7 */
  in += 6;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);


  /* 8 */
  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_16_fwd_8_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  debug(printf("Entering unpack_16_fwd_8_1\n"));

  /* 1 first */
  OutReg = _mm_load_si128(in);
  _mm_store_si128(&(out[1]), OutReg);

  /* 8 second */
  in += 7;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  return;
}


#if !defined(HAVE_SSE2)
typedef void (*Unpacker_T) (UINT2* __restrict__, const UINT2* __restrict__);
#else
typedef void (*Unpacker_T) (__m128i* __restrict__, const __m128i* __restrict__);
#endif


#if !defined(HAVE_SSE2)
static Unpacker_T unpacker_all_table[17] =
  {unpack_00,
   unpack_00, unpack_02, unpack_00, unpack_04,
   unpack_00, unpack_06, unpack_00, unpack_08,
   unpack_00, unpack_10, unpack_00, unpack_12,
   unpack_00, unpack_14, unpack_00, unpack_16};

#else

/* Entry 16 in each packsize handles remainder == 64 => quarter_block == 4, column 3, row -1 */
static Unpacker_T unpacker_table[9][9] = 
  {{unpack_00_fwd_0_0, unpack_00_fwd_0_0, unpack_00_fwd_0_0, unpack_00_fwd_0_0,
    unpack_00_fwd_0_0, unpack_00_fwd_0_0, unpack_00_fwd_0_0, unpack_00_fwd_0_0,
    unpack_00_fwd_0_0},

   {unpack_02_fwd_1_2, unpack_02_fwd_2_3, unpack_02_fwd_3_4, unpack_02_fwd_4_5,
    unpack_02_fwd_5_6, unpack_02_fwd_6_7, unpack_02_fwd_7_8, unpack_02_fwd_8_1,
    unpack_00_fwd_0_0},

   {unpack_04_fwd_1_2, unpack_04_fwd_2_3, unpack_04_fwd_3_4, unpack_04_fwd_4_5,
    unpack_04_fwd_5_6, unpack_04_fwd_6_7, unpack_04_fwd_7_8, unpack_04_fwd_8_1,
    unpack_00_fwd_0_0},

   {unpack_06_fwd_1_2, unpack_06_fwd_2_3, unpack_06_fwd_3_4, unpack_06_fwd_4_5,
    unpack_06_fwd_5_6, unpack_06_fwd_6_7, unpack_06_fwd_7_8, unpack_06_fwd_8_1,
    unpack_00_fwd_0_0},

   {unpack_08_fwd_1_2, unpack_08_fwd_2_3, unpack_08_fwd_3_4, unpack_08_fwd_4_5,
    unpack_08_fwd_5_6, unpack_08_fwd_6_7, unpack_08_fwd_7_8, unpack_08_fwd_8_1,
    unpack_00_fwd_0_0},

   {unpack_10_fwd_1_2, unpack_10_fwd_2_3, unpack_10_fwd_3_4, unpack_10_fwd_4_5,
    unpack_10_fwd_5_6, unpack_10_fwd_6_7, unpack_10_fwd_7_8, unpack_10_fwd_8_1,
    unpack_00_fwd_0_0},

   {unpack_12_fwd_1_2, unpack_12_fwd_2_3, unpack_12_fwd_3_4, unpack_12_fwd_4_5,
    unpack_12_fwd_5_6, unpack_12_fwd_6_7, unpack_12_fwd_7_8, unpack_12_fwd_8_1,
    unpack_00_fwd_0_0},

   {unpack_14_fwd_1_2, unpack_14_fwd_2_3, unpack_14_fwd_3_4, unpack_14_fwd_4_5,
    unpack_14_fwd_5_6, unpack_14_fwd_6_7, unpack_14_fwd_7_8, unpack_14_fwd_8_1,
    unpack_00_fwd_0_0},

   {unpack_16_fwd_1_2, unpack_16_fwd_2_3, unpack_16_fwd_3_4, unpack_16_fwd_4_5,
    unpack_16_fwd_5_6, unpack_16_fwd_6_7, unpack_16_fwd_7_8, unpack_16_fwd_8_1,
    unpack_00_fwd_0_0},

};
#endif



#define BLOCKSIZE 64
#define METAINFO_SIZE 2

#define get_column(s) (s) & 7 /* Not s % 8, which fails on negative values */
#define get_row(s) (s) >> 3 /* Not s / 8, which fails on negative values */


/* bitpackpages: A list of b-mers (12-mers by default), ending with -1U */
/* Returning UINT4 instead of UINT2, to avoid the problem where end0 could be 65536, which is represented by 0 */
UINT4
Epu16_bitpack64_read_two (UINT4 *end0, Localspace_T oligo, UINT2 *bitpackptrs, UINT2 *bitpackcomp) {
  UINT4 ptr;
  Localspace_T bmer;
  UINT2 *info, nwritten4, packsize_div2;
  UINT2 nwritten8;
  UINT2 offset0;
  int column;
#if !defined(HAVE_SSE2)
  int remainder, row, k;
  UINT2 diffs[BLOCKSIZE+1], *bitpack, offset1;
#else
  __m128i diffs[2];  /* Need to provide space for ptr and for end0 */
  int remainder0;
  int delta, row0, row1;
  __m128i *bitpack;
  UINT2 *_diffs;
  int i;
#endif
#ifdef DEBUG
  UINT2 offsets[BLOCKSIZE+1];
#endif


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * METAINFO_SIZE]);

  debug(printf("Entered Epu16_bitpack64_read_two with oligo %llu => bmer %u\n",oligo,bmer));
#if 0
  printf("At entry, bitpack points to %04X %04X %04X %04X %04X %04X %04X %04X\n",
	 bitpackcomp[0],bitpackcomp[1],bitpackcomp[2],bitpackcomp[3],
	 bitpackcomp[4],bitpackcomp[5],bitpackcomp[5],bitpackcomp[7]);
#endif

  /* Because for localdb, we not writing 0 0 at the start */
  if (bmer == 0) {
    nwritten4 = 0;
    offset0 = 0;
  } else {
    nwritten4 = info[-METAINFO_SIZE];
    offset0 = info[-METAINFO_SIZE+1];
  }
  nwritten8 = 8 * nwritten4;	       /* In 16-bit shorts */
  packsize_div2 = info[0] - nwritten4;

#if !defined(HAVE_SSE2)
  bitpack = (UINT2 *) &(bitpackcomp[nwritten8]);
  offset1 = info[1];
#else
  bitpack = (__m128i *) &(bitpackcomp[nwritten8]);
#endif


  debug(printf("nwritten %u registers or %u shorts, packsize %d\n",
	       nwritten4,nwritten8,packsize_div2 * 2));
#ifdef DEBUG
  printf("bitpack now points to %04X %04X %04X %04X %04X %04X %04X %04X\n",
	 bitpackcomp[nwritten8+0],bitpackcomp[nwritten8+1],
	 bitpackcomp[nwritten8+2],bitpackcomp[nwritten8+3],
	 bitpackcomp[nwritten8+4],bitpackcomp[nwritten8+5],
	 bitpackcomp[nwritten8+5],bitpackcomp[nwritten8+7]);
#endif

  debug(Epu16_bitpack64_unpack(&(bitpackcomp)[nwritten8],packsize_div2 * 2));

#if !defined(HAVE_SSE2)
  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);

#ifdef DEBUG
  printf("oligo: %08X, remainder %d, offset0 %u, offset1 %u\n",
	 oligo,oligo % BLOCKSIZE,info[1],info[METAINFO_SIZE+1]);
  printf("bitpack:\n");

  for (i = 1; i <= BLOCKSIZE; i++) {
    printf("%d ",diffs[i]);
    if (i % (BLOCKSIZE/8) == 0) {
      printf("\n");
    }
  }
  printf("\n");
  printf("end of diffs\n");
#endif  

  ptr = (UINT4) offset0;
  if ((remainder = oligo % BLOCKSIZE) == 0) {
    /* Just return offset0 */

  } else {
    column = (remainder - 1) % 8; /* Goes from 0 to 7 */
    row = (remainder - 1) / 8;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column + 1, i = 0; i <= row; k += BLOCKSIZE/8, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }
  }

  remainder++;
  if (remainder == 64) {
    *end0 = offset1;

  } else {
    /* Compute necessary cumulative sums */
    *end0 = offset0;

    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    diffs[0] = 0;
    column = (remainder - 1) % 8; /* Goes from 0 to 7 */
    row = (remainder - 1) / 8;
    debug(printf("column %d, row %d\n",column,row));
    
    for (k = column + 1, i = 0; i <= row; k += BLOCKSIZE/8, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      *end0 += diffs[k];
    }
  }

#else			    /* SSE2 */
  remainder0 = oligo % BLOCKSIZE;
  /* remainder1 = remainder0 + 1; */

  _diffs = (UINT2 *) diffs;	/* Assumes a dummy register in diffs[0] */

  delta = remainder0 - 1;
  column = get_column(delta);
  row0 = get_row(delta);
  row1 = get_row(delta + 1);
  debug(printf("Column is %d, Rows are %d and %d\n",column,row0,row1));
  debug(printf("Calling unpacker_table[%d][%d]\n",packsize_div2,column));
  assert(packsize_div2 <= 9);
  (unpacker_table[packsize_div2][column])(diffs,bitpack);

  debug(printf("offset0 is %u\n",offset0));

  _diffs = (UINT2 *) &(diffs[1]);
  debug(printf("diffs for end: %u %u %u %u %u %u %u %u\n",
	       _diffs[0],_diffs[1],_diffs[2],_diffs[3],_diffs[4],_diffs[5],_diffs[6],_diffs[7]));
  *end0 = (UINT4) offset0;
  for (i = 0; i <= row1; i++) {
    *end0 += _diffs[i];
    debug(printf("For end, adding diffs[%d] = %u\n",i,_diffs[i]));
  }
  debug(printf("end0 is %u\n",*end0));

  _diffs = (UINT2 *) &(diffs[0]);
  debug(printf("diffs for ptr: %u %u %u %u %u %u %u %u\n",
	       _diffs[0],_diffs[1],_diffs[2],_diffs[3],_diffs[4],_diffs[5],_diffs[6],_diffs[7]));
  ptr = (UINT4) offset0;
  for (i = 0; i <= row0; i++) {
    ptr += _diffs[i];
    debug(printf("For ptr, adding diffs[%d] = %u\n",i,_diffs[i]));
  }
  debug(printf("ptr0 is %u\n",ptr));

#endif	/* HAVE_SSE2 */

  return ptr;

}
