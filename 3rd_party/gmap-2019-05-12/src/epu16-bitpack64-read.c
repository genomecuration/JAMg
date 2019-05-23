static char rcsid[] = "$Id: epu16-bitpack64-read.c 218164 2019-01-17 06:09:06Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "epu16-bitpack64-read.h"

#include <stdio.h>
#include <stdlib.h>

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

/* block offsets */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
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


#if 0
#ifdef HAVE_SSE2
#ifdef ALLOW_ODD_PACKSIZES
static __m128i mask1, mask2, mask3, mask4, mask5, mask6, mask7, mask8,
  mask9, mask10, mask11, mask12, mask13, mask14, mask15;
#else
static __m128i mask2, mask4, mask6, mask8, mask10, mask12, mask14;
#endif
#endif
#endif


#define BLOCKSIZE 64

#if 0
void
Epu16_bitpack64_read_setup () {

#ifdef HAVE_SSE2
#ifdef ALLOW_ODD_PACKSIZES
  mask1 = _mm_set1_epi16(1U);
  mask3 =  _mm_set1_epi16(7U);
  mask5 =  _mm_set1_epi16(31U);
  mask7 =  _mm_set1_epi16(127U);
  mask9 =  _mm_set1_epi16(511U);
  mask11 =  _mm_set1_epi16(2047U);
  mask13 =  _mm_set1_epi16(8191U);
  mask15 =  _mm_set1_epi16(32767U);
#endif
  mask2 = _mm_set1_epi16(3U);
  mask4 =  _mm_set1_epi16(15U);
  mask6 =  _mm_set1_epi16(63U);
  mask8 =  _mm_set1_epi16(255U);
  mask10 =  _mm_set1_epi16(1023U);
  mask12 =  _mm_set1_epi16(4095U);
  mask14 =  _mm_set1_epi16(16383U);
#endif

  return;
}
#endif



#if !defined(HAVE_SSE2)
static void
unpack_00 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  int i;

  for (i = 0; i < BLOCKSIZE; i++) {
    *out++ = 0;
  }

  return;
}

#else
static void
unpack_00 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i total = _mm_set1_epi16(0U);

  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);

  return;
}

static void
unpack_00_0 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i zero = _mm_set1_epi16(0U);

  _mm_store_si128(out++, zero);

  return;
}
#endif




#if !defined(HAVE_SSE2)
static void
unpack_02 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  unsigned int column;
  const UINT2 *bitpack = in;

  for (column = 0; column < 8; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  2  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  4  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  6  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 2 ) ;
    out++;
  }

  return;
}

#else

/* Note: the commented-out lines for total are vestiges from the
   original code, which was designed for serial decoding and assumed a
   vertical block structure, rather than a columnar one */

static void
unpack_02_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg = _mm_load_si128(in);
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  OutReg = _mm_and_si128( InReg , mask2);
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask2);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask2);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,6) , mask2);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask2);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,10) , mask2);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,12) , mask2);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,14) , mask2);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_02_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_02_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_02_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_02_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,6) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_02_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}
  
static void
unpack_02_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,10) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_02_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,12) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_02_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask2 = _mm_set1_epi16(3U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,14) , mask2);
  _mm_store_si128(out++, OutReg);

  return;
}
#endif



#if !defined(HAVE_SSE2)
static void
unpack_04 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  unsigned int column;
  const UINT2 *bitpack = in;

  for (column = 0; column < 8; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 4 ) ;
    out++;
    *out = ( (*in) >>  4  )   % (1U << 4 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 4 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 4 ) ;
    out++;
    in += 8;
    *out = ( (*in) >>  0  )   % (1U << 4 ) ;
    out++;
    *out = ( (*in) >>  4  )   % (1U << 4 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 4 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 4 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_04_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg = _mm_load_si128(in);
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  OutReg = _mm_and_si128( InReg , mask4);
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask4);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask4);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_srli_epi16(InReg,12);
  InReg = _mm_load_si128(++in);

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128( InReg , mask4);
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask4);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask4);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_srli_epi16(InReg,12);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  InReg = _mm_load_si128(in);

  OutReg =  _mm_srli_epi16(InReg,12);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128( InReg , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask4 = _mm_set1_epi16(15U);

  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask4);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_04_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  InReg = _mm_load_si128(++in);

  OutReg =  _mm_srli_epi16(InReg,12);
  _mm_store_si128(out++, OutReg);

  return;
}
#endif


#if !defined(HAVE_SSE2)
static void
unpack_06 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  unsigned int column;
  const UINT2 *bitpack = in;

  for (column = 0; column < 8; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  6  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 6 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 2 ))<<( 6 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 6 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 4 ))<<( 6 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 6 ) ;
    in += 8;
    out++;
  }

  return;
}

#else

static void
unpack_06_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg = _mm_load_si128(in);
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  OutReg = _mm_and_si128( InReg , mask6);
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,6) , mask6);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 6 - 2), mask6));
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask6);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask6);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 6 - 4), mask6));
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask6);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,6) , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 6 - 2), mask6));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  InReg = _mm_load_si128(++in);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,8) , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  InReg = _mm_load_si128(++in);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 6 - 4), mask6));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask6 =  _mm_set1_epi16(63U);

  in += 2;
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask6);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_06_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  _mm_store_si128(out++, OutReg);

  return;
}
#endif



#if !defined(HAVE_SSE2)
static void
unpack_08 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  unsigned int column;
  const UINT2 *bitpack = in;

  for (column = 0; column < 8; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 8 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 8 ) ;
    out++;
    in += 8;
    *out = ( (*in) >>  0  )   % (1U << 8 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 8 ) ;
    out++;
    in += 8;
    *out = ( (*in) >>  0  )   % (1U << 8 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 8 ) ;
    out++;
    in += 8;
    *out = ( (*in) >>  0  )   % (1U << 8 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 8 ) ;
    out++;
  }

  return;
}

#else

static void
unpack_08_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg = _mm_load_si128(in);
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  OutReg = _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128( InReg , mask8);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128( InReg , mask8);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128( InReg , mask8);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);
  
  return;
}


static void
unpack_08_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  InReg = _mm_load_si128(++in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  InReg = _mm_load_si128(++in);

  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask8 =  _mm_set1_epi16(255U);

  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =  _mm_and_si128( InReg , mask8);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_08_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =  _mm_srli_epi16(InReg,8) ;
  _mm_store_si128(out++, OutReg);

  return;
}
#endif


#if !defined(HAVE_SSE2)
static void
unpack_10 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  unsigned int column;
  const UINT2 *bitpack = in;

  for (column = 0; column < 8; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 10 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 4 ))<<( 10 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 10 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 8 ))<<( 10 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 10 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 2 ))<<( 10 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 10 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 6 ))<<( 10 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 10 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_10_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg = _mm_load_si128(in);
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  OutReg = _mm_and_si128( InReg , mask10);
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 4), mask10));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask10);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 8), mask10));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 2), mask10));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask10);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 6), mask10));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,6) ;
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_10_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask10);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_10_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 4), mask10));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_10_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  InReg = _mm_load_si128(++in);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,4) , mask10);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_10_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  InReg = _mm_load_si128(++in);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 8), mask10));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_10_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 2), mask10));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_10_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  in += 3;
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128(  _mm_srli_epi16(InReg,2) , mask10);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_10_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask10 =  _mm_set1_epi16(1023U);

  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 10 - 6), mask10));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_10_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  in += 4;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,6) ;
  _mm_store_si128(out++, OutReg);

  return;
}
#endif


#if !defined(HAVE_SSE2)

static void
unpack_12 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  unsigned int column;
  const UINT2 *bitpack = in;

  for (column = 0; column < 8; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 12 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 8 ))<<( 12 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 12 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 4 ))<<( 12 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 12 ) ;
    in += 8;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 12 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 8 ))<<( 12 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 12 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 4 ))<<( 12 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 12 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_12_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg = _mm_load_si128(in);
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  OutReg = _mm_and_si128( InReg , mask12);
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 8), mask12));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 4), mask12));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  InReg = _mm_load_si128(++in);

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_and_si128( InReg , mask12);
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 8), mask12));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 4), mask12));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_12_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask12);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_12_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 8), mask12));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_12_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  InReg = _mm_load_si128(++in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 4), mask12));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_12_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_12_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  in += 3;
  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask12);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_12_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 8), mask12));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_12_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask12 =  _mm_set1_epi16(4095U);

  in += 4;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 12 - 4), mask12));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_12_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  in += 5;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  _mm_store_si128(out++, OutReg);

  return;
}
#endif


#if !defined(HAVE_SSE2)

static void
unpack_14 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  unsigned int column;
  const UINT2 *bitpack = in;

  for (column = 0; column < 8; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 14 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 12 ))<<( 14 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 14 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 10 ))<<( 14 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 14 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 8 ))<<( 14 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 14 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 6 ))<<( 14 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 14 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 4 ))<<( 14 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 14 ) ;
    in += 8;
    *out |= ((*in) % (1U<< 2 ))<<( 14 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 14 ) ;
    out++;
  }

  return;
}

#else
static void
unpack_14_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg = _mm_load_si128(in);
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  OutReg = _mm_and_si128( InReg , mask14);
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 12), mask14));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 10), mask14));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 8), mask14));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 6), mask14));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,6) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 4), mask14));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 2), mask14));

  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  OutReg =   _mm_srli_epi16(InReg,2) ;
  /* total = _mm_add_epi16(total, OutReg); */
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_14_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  InReg = _mm_load_si128(in);

  OutReg = _mm_and_si128( InReg , mask14);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_14_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,14) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 12), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_14_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  InReg = _mm_load_si128(++in);

  OutReg =   _mm_srli_epi16(InReg,12) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 10), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_14_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  in += 2;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,10) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 8), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_14_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  in += 3;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,8) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 6), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_14_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  in += 4;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,6) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 4), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}


static void
unpack_14_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;
  const __m128i mask14 =  _mm_set1_epi16(16383U);

  in += 5;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,4) ;
  InReg = _mm_load_si128(++in);

  OutReg = _mm_or_si128(OutReg, _mm_and_si128(_mm_slli_epi16(InReg, 14 - 2), mask14));
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_14_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i InReg;
  __m128i OutReg;

  in += 6;
  InReg = _mm_load_si128(in);

  OutReg =   _mm_srli_epi16(InReg,2) ;
  _mm_store_si128(out++, OutReg);

  return;
}
#endif


#if !defined(HAVE_SSE2)
static void
unpack_16 (UINT2* __restrict__ out, const UINT2* __restrict__ in) {
  unsigned int column;
  const UINT2 *bitpack = in;

  for (column = 0; column < 8; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   ;
    in += 8;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 8;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 8;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 8;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 8;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 8;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 8;
    out++;
    *out = ( (*in) >>  0  )   ;
    out++;
  }

  return;
}

#else
static void
unpack_16_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  OutReg = _mm_load_si128(in++);
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_load_si128(in++);
  /* total = _mm_add_epi16(total, _mm_load_si128(in++)); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_load_si128(in++);
  /* total = _mm_add_epi16(total, _mm_load_si128(in++)); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_load_si128(in++);
  /* total = _mm_add_epi16(total, _mm_load_si128(in++)); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_load_si128(in++);
  /* total = _mm_add_epi16(total, _mm_load_si128(in++)); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_load_si128(in++);
  /* total = _mm_add_epi16(total, _mm_load_si128(in++)); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_load_si128(in++);
  /* total = _mm_add_epi16(total, _mm_load_si128(in++)); */
  _mm_store_si128(out++, OutReg);

  OutReg = _mm_load_si128(in++);
  /* total = _mm_add_epi16(total, _mm_load_si128(in++)); */
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_1 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_2 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  OutReg = _mm_load_si128(++in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_3 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  in += 2;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_4 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  in += 3;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_5 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  in += 4;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_6 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  in += 5;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_7 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  in += 6;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  return;
}

static void
unpack_16_fwd_8 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i OutReg;

  in += 7;

  OutReg = _mm_load_si128(in);
  _mm_store_si128(out++, OutReg);

  return;
}
#endif


/* Note: output is UINT4, but input is UINT2 */
static void
vertical_order (UINT4 *vertical, UINT2 *columnar) {

  vertical[0] = columnar[0];
  vertical[8] = columnar[1];
  vertical[16] = columnar[2];
  vertical[24] = columnar[3];
  vertical[32] = columnar[4];
  vertical[40] = columnar[5];
  vertical[48] = columnar[6];
  vertical[56] = columnar[7];

  vertical[1] = columnar[8];
  vertical[9] = columnar[9];
  vertical[17] = columnar[10];
  vertical[25] = columnar[11];
  vertical[33] = columnar[12];
  vertical[41] = columnar[13];
  vertical[49] = columnar[14];
  vertical[57] = columnar[15];

  vertical[2] = columnar[16];
  vertical[10] = columnar[17];
  vertical[18] = columnar[18];
  vertical[26] = columnar[19];
  vertical[34] = columnar[20];
  vertical[42] = columnar[21];
  vertical[50] = columnar[22];
  vertical[58] = columnar[23];

  vertical[3] = columnar[24];
  vertical[11] = columnar[25];
  vertical[19] = columnar[26];
  vertical[27] = columnar[27];
  vertical[35] = columnar[28];
  vertical[43] = columnar[29];
  vertical[51] = columnar[30];
  vertical[59] = columnar[31];

  vertical[4] = columnar[32];
  vertical[12] = columnar[33];
  vertical[20] = columnar[34];
  vertical[28] = columnar[35];
  vertical[36] = columnar[36];
  vertical[44] = columnar[37];
  vertical[52] = columnar[38];
  vertical[60] = columnar[39];

  vertical[5] = columnar[40];
  vertical[13] = columnar[41];
  vertical[21] = columnar[42];
  vertical[29] = columnar[43];
  vertical[37] = columnar[44];
  vertical[45] = columnar[45];
  vertical[53] = columnar[46];
  vertical[61] = columnar[47];

  vertical[6] = columnar[48];
  vertical[14] = columnar[49];
  vertical[22] = columnar[50];
  vertical[30] = columnar[51];
  vertical[38] = columnar[52];
  vertical[46] = columnar[53];
  vertical[54] = columnar[54];
  vertical[62] = columnar[55];

  vertical[7] = columnar[56];
  vertical[15] = columnar[57];
  vertical[23] = columnar[58];
  vertical[31] = columnar[59];
  vertical[39] = columnar[60];
  vertical[47] = columnar[61];
  vertical[55] = columnar[62];
  vertical[63] = columnar[63];

  return;
}



#if !defined(HAVE_SSE2)
typedef void (*Unpacker_T) (UINT2* __restrict__, const UINT2* __restrict__);
#else
typedef void (*Unpacker_T) (__m128i* __restrict__, const __m128i* __restrict__);
#endif


#ifdef ALLOW_ODD_PACKSIZES
static Unpacker_T unpacker_table[17] =
  {unpack_00,
   unpack_01, unpack_02, unpack_03, unpack_04,
   unpack_05, unpack_06, unpack_07, unpack_08,
   unpack_09, unpack_10, unpack_11, unpack_12,
   unpack_13, unpack_14, unpack_15, unpack_16};

#elif !defined(HAVE_SSE2)
static Unpacker_T unpacker_all_table[17] =
  {unpack_00,
   unpack_00, unpack_02, unpack_00, unpack_04,
   unpack_00, unpack_06, unpack_00, unpack_08,
   unpack_00, unpack_10, unpack_00, unpack_12,
   unpack_00, unpack_14, unpack_00, unpack_16};

#else
static Unpacker_T unpacker_all_table[17] =
  {unpack_00,
   unpack_00, unpack_02_fwd, unpack_00, unpack_04_fwd,
   unpack_00, unpack_06_fwd, unpack_00, unpack_08_fwd,
   unpack_00, unpack_10_fwd, unpack_00, unpack_12_fwd,
   unpack_00, unpack_14_fwd, unpack_00, unpack_16_fwd};

/* Entry 9 in each packsize handles remainder == 64 => quarter_block == 4, column 3, row -1 */
static Unpacker_T unpacker_table[9][9] = 
  {{unpack_00_0, unpack_00_0, unpack_00_0, unpack_00_0,
    unpack_00_0, unpack_00_0, unpack_00_0, unpack_00_0,
    unpack_00_0},

   {unpack_02_fwd_1, unpack_02_fwd_2, unpack_02_fwd_3, unpack_02_fwd_4, 
    unpack_02_fwd_5, unpack_02_fwd_6, unpack_02_fwd_7, unpack_02_fwd_8, 
    unpack_00_0},

   {unpack_04_fwd_1, unpack_04_fwd_2, unpack_04_fwd_3, unpack_04_fwd_4, 
    unpack_04_fwd_5, unpack_04_fwd_6, unpack_04_fwd_7, unpack_04_fwd_8, 
    unpack_00_0},

   {unpack_06_fwd_1, unpack_06_fwd_2, unpack_06_fwd_3, unpack_06_fwd_4, 
    unpack_06_fwd_5, unpack_06_fwd_6, unpack_06_fwd_7, unpack_06_fwd_8, 
    unpack_00_0},

   {unpack_08_fwd_1, unpack_08_fwd_2, unpack_08_fwd_3, unpack_08_fwd_4, 
    unpack_08_fwd_5, unpack_08_fwd_6, unpack_08_fwd_7, unpack_08_fwd_8, 
    unpack_00_0},

   {unpack_10_fwd_1, unpack_10_fwd_2, unpack_10_fwd_3, unpack_10_fwd_4, 
    unpack_10_fwd_5, unpack_10_fwd_6, unpack_10_fwd_7, unpack_10_fwd_8, 
    unpack_00_0},

   {unpack_12_fwd_1, unpack_12_fwd_2, unpack_12_fwd_3, unpack_12_fwd_4, 
    unpack_12_fwd_5, unpack_12_fwd_6, unpack_12_fwd_7, unpack_12_fwd_8, 
    unpack_00_0},

   {unpack_14_fwd_1, unpack_14_fwd_2, unpack_14_fwd_3, unpack_14_fwd_4, 
    unpack_14_fwd_5, unpack_14_fwd_6, unpack_14_fwd_7, unpack_14_fwd_8, 
    unpack_00_0},

   {unpack_16_fwd_1, unpack_16_fwd_2, unpack_16_fwd_3, unpack_16_fwd_4, 
    unpack_16_fwd_5, unpack_16_fwd_6, unpack_16_fwd_7, unpack_16_fwd_8, 
    unpack_00_0},

};
   
#endif


#define return_sum_fwd(offset0,diffs,row) switch (row) {	\
  case -1: return offset0; break;					\
  case 0: return offset0 + diffs[0]; break;				\
  case 1: return offset0 + diffs[0] + diffs[1]; break;			\
  case 2: return offset0 + diffs[0] + diffs[1] + diffs[2]; break;	\
  case 3: return offset0 + diffs[0] + diffs[1] + diffs[2] + diffs[3]; break; \
  case 4: return offset0 + diffs[0] + diffs[1] + diffs[2] + diffs[3] + diffs[4]; break; \
  case 5: return offset0 + diffs[0] + diffs[1] + diffs[2] + diffs[3] + diffs[4] + diffs[5]; break; \
  case 6: return offset0 + diffs[3] + diffs[4] + diffs[5] + diffs[6]; break; \
  case 7: return offset0 + diffs[4] + diffs[5] + diffs[6] + diffs[7]; break; \
  default: abort();							\
  }


#define DIFFERENTIAL_METAINFO_SIZE 2
#define METAINFO_SIZE 2

#define get_column(s) (s) & 7 /* Not s % 8, which fails on negative values */
#define get_row(s) (s) >> 3 /* Not s / 8, which fails on negative values */


UINT4
Epu16_bitpack64_read_one (Localspace_T oligo, UINT2 *bitpackptrs, UINT2 *bitpackcomp) {
  UINT4 ptr;
  Localspace_T bmer;
  UINT2 *info, nwritten4, packsize_div2;
  UINT8 nwritten8;
  int delta, remainder, column, row, i;
  UINT2 offset0;
#if !defined(HAVE_SSE2)
  UINT2 diffs[BLOCKSIZE+1], *bitpack;
  int k;
#else
  __m128i diffs[2];
  __m128i *bitpack;
  UINT2 *_diffs;
#endif


  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * METAINFO_SIZE]);

  debug(printf("Entered Bitpack64_read_one with oligo %u => bmer %u\n",oligo,bmer));

  /* Because for localdb, we not writing 0 0 at the start */
  if (bmer == 0) {
    nwritten4 = 0;
    offset0 = 0;
  } else {
    nwritten4 = info[-METAINFO_SIZE];
    offset0 = info[-METAINFO_SIZE+1];
  }
  debug(printf("Starting with offset0 %u\n",offset0));
  nwritten8 = 8 * nwritten4;	       /* In 16-bit shorts */
  packsize_div2 = info[0] - nwritten4;

#if !defined(HAVE_SSE2)
  bitpack = (UINT2 *) &(bitpackcomp[nwritten8]);
#else
  bitpack = (__m128i *) &(bitpackcomp[nwritten8]);
#endif

  remainder = oligo % BLOCKSIZE;

#if !defined(HAVE_SSE2)

  ptr = (UINT4) offset0;
  if (remainder == 0) {
    /* Just use offset0 */

  } else {
    delta = remainder - 1;
    column = get_column(delta);
    row = get_row(delta);
    debug(printf("column %d, row %d\n",column,row));
    
    /* Unpack all 64 diffs for non-SIMD */
    (unpacker_all_table[packsize_div2*2])(&(diffs[1]),bitpack);
    for (k = column + 1, i = 0; i <= row; k += BLOCKSIZE/8, i++) {
      debug(printf("Adding diffs[%d] = %u\n",k,diffs[k]));
      ptr += diffs[k];
    }
  }

#else  /* SSE2 */
  _diffs = (UINT2 *) diffs;	/* Assumes a dummy register in diffs[0] */

  delta = remainder - 1;
  column = get_column(delta);
  row = get_row(delta);
  debug(printf("delta %d, column %d, row %d\n",delta,column,row));

  (unpacker_table[packsize_div2][column*8])(diffs,bitpack);

  ptr = (UINT4) offset0;
  for (i = 0; i <= row; i++) {
    ptr += (UINT4) _diffs[i];
  }

#endif	/* littleendian and SSE2 */

  debug(printf("Returning %u\n",ptr));
  return ptr;

}


/* Unpack all offsets.  Can treat offset0 as a special case */
/* DONE: Corrected for shift of meta up by 1 */
/* offsets needs to be UINT4 to hold 65536 as the last value */
void
Epu16_bitpack64_block_offsets (UINT4 *offsets, Localspace_T oligo,
			       UINT2 *bitpackptrs, UINT2 *bitpackcomp) {
  Localspace_T bmer;
  UINT2 *info, nwritten4;
  UINT8 nwritten8;
  UINT4 offset0;
  int packsize, k;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  int column, row;
  UINT2 diffs[BLOCKSIZE], columnar[BLOCKSIZE], *bitpack, *vertical;
#else
  __m128i diffs[8], *bitpack;
  UINT2 *_diffs;
#endif
#ifdef DEBUG1
  UINT4 offset1;
  int i;
#endif

  bmer = oligo/BLOCKSIZE;
  info = &(bitpackptrs[bmer * DIFFERENTIAL_METAINFO_SIZE]);
  if (bmer == 0) {
    nwritten4 = 0;
    offset0 = 0;
  } else {
    /* For localdb, meta is shifted by 1 (we don't write the first offset of 0 0) */
    nwritten4 = info[-DIFFERENTIAL_METAINFO_SIZE]; /* In 128-bit registers */
    offset0 = (UINT4) info[-DIFFERENTIAL_METAINFO_SIZE+1];
  }
  nwritten8 = 8 * (UINT8) nwritten4; /* In 16-bit shorts */
  packsize = (info[0] - nwritten4)*2;

#if !defined(HAVE_SSE2)
  bitpack = (UINT2 *) &(bitpackcomp[nwritten8]);
#else
  bitpack = (__m128i *) &(bitpackcomp[nwritten8]);
#endif

#ifdef DEBUG1
  offset1 = (UINT4) info[1];
  printf("oligo: %08X, nwritten %u, offset0 %u, offset1 %u, packsize %d\n",
	 oligo,nwritten4,offset0,offset1,packsize);
#endif


#if !defined(HAVE_SSE2)
  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize])(&(diffs[0]),bitpack);

#ifdef DEBUG1
  for (i = 0; i < BLOCKSIZE; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   diffs[i],diffs[i+1],diffs[i+2],diffs[i+3],
	   diffs[i+4],diffs[i+5],diffs[i+6],diffs[i+7]);
  }
  printf("end of diffs horizontal (because non-SIMD unpackers are horizontal)\n");
#endif

  /* Convert to columnar */
  vertical = &(diffs[0]);
  for (column = 0; column < 8; column++) {
    k = column;
    for (row = 0; row < BLOCKSIZE/8; row++) {
      columnar[k] = *vertical++;
      k += 8;
    }
  }

  /* Convert to vertical, shifting by 1 */
  vertical_order(&(offsets[1]),columnar);

#ifdef DEBUG1
  printf("%u\n",offset0);
  for (i = 1; i <= BLOCKSIZE; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   offsets[i],offsets[i+1],offsets[i+2],offsets[i+3],
	   offsets[i+4],offsets[i+5],offsets[i+6],offsets[i+7]);
  }
  printf("end of diffs vertical\n");
#endif

#else  /* SSE2 */

#if 0
  printf("bitpack:\n");
  for (i = 0; i < packsize/2; i++) {
    print_vector_hex(bitpack[i]);
  }
  printf("\n");
#endif  

  _diffs = (UINT2 *) &(diffs[0]);

  /* Unpack fwd 64 cumulative sums under SIMD */
  (unpacker_all_table[packsize])(&(diffs[0]),bitpack);
  vertical_order(&(offsets[1]),_diffs);

#ifdef DEBUG1
  printf("%u\n",offsets[i]);
  for (i = 1; i <= BLOCKSIZE; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   offsets[i],offsets[i+1],offsets[i+2],offsets[i+3],
	   offsets[i+4],offsets[i+5],offsets[i+6],offsets[i+7]);
  }
  printf("end of diffs vertical\n");
#endif

#endif	/* SSE2 */

  /* Perform cumulative sum */
  offsets[0] = offset0;
  offsets[1] += offset0;
  offsets[2] += offset0;
  offsets[3] += offset0;
  offsets[4] += offset0;
  offsets[5] += offset0;
  offsets[6] += offset0;
  offsets[7] += offset0;
  offsets[8] += offset0;

  for (k = 9; k <= 64; k++) {
    offsets[k] += offsets[k-8];
  }

#ifdef DEBUG1
  printf("%u\n",offsets[0]);
  for (i = 1; i <= 64; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   offsets[i],offsets[i+1],offsets[i+2],offsets[i+3],
	   offsets[i+4],offsets[i+5],offsets[i+6],offsets[i+7]);
  }
  printf("end of offsets\n");
#endif

  return;
}



#if 0
/* Used only for debugging */
void
Epu16_bitpack64_unpack (UINT2 *bitpackptr, int packsize) {
  UINT4 offsets[BLOCKSIZE+1], offset0 = 0;
  int k;
#if defined(WORDS_BIGENDIAN) || !defined(HAVE_SSE2)
  int column, row;
  UINT2 diffs[BLOCKSIZE], columnar[BLOCKSIZE], *bitpack, *vertical;
#else
  __m128i diffs[8], *bitpack;
  UINT2 *_diffs;
#endif
  int i;


#if !defined(HAVE_SSE2)
  bitpack = (UINT2 *) bitpackptr;

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_all_table[packsize])(&(diffs[0]),bitpack);

  for (i = 0; i < BLOCKSIZE; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   diffs[i],diffs[i+1],diffs[i+2],diffs[i+3],
	   diffs[i+4],diffs[i+5],diffs[i+6],diffs[i+7]);
  }
  printf("end of diffs horizontal (because non-SIMD unpackers are horizontal)\n");

  /* Convert to columnar */
  vertical = &(diffs[0]);
  for (column = 0; column < 8; column++) {
    k = column;
    for (row = 0; row < BLOCKSIZE/8; row++) {
      columnar[k] = *vertical++;
      k += 8;
    }
  }

  /* Convert to vertical, shifting by 1 */
  vertical_order(&(offsets[1]),columnar);

  for (i = 1; i <= BLOCKSIZE; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   offsets[i],offsets[i+1],offsets[i+2],offsets[i+3],
	   offsets[i+4],offsets[i+5],offsets[i+6],offsets[i+7]);
  }
  printf("end of diffs vertical\n");

#else  /* SSE2 */
  bitpack = (__m128i *) bitpackptr;

  printf("bitpack:\n");
  for (i = 0; i < packsize/2; i++) {
    print_vector_hex(bitpack[i]);
  }
  printf("\n");

  _diffs = (UINT2 *) &(diffs[0]);

  /* Unpack fwd 64 cumulative sums under SIMD */
  (unpacker_all_table[packsize])(&(diffs[0]),bitpack);
  vertical_order(&(offsets[1]),_diffs);

  for (i = 1; i <= BLOCKSIZE; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   offsets[i],offsets[i+1],offsets[i+2],offsets[i+3],
	   offsets[i+4],offsets[i+5],offsets[i+6],offsets[i+7]);
  }
  printf("end of diffs vertical\n");

#endif	/* SSE2 */

  /* Perform cumulative sum */
  offset0 = 0;
  offsets[0] = offset0;
  offsets[1] += offset0;
  offsets[2] += offset0;
  offsets[3] += offset0;
  offsets[4] += offset0;
  offsets[5] += offset0;
  offsets[6] += offset0;
  offsets[7] += offset0;
  offsets[8] += offset0;

  for (k = 9; k <= 64; k++) {
    offsets[k] += offsets[k-8];
  }

#ifdef DEBUG
  printf("%u\n",offsets[0]);
  for (i = 1; i <= 64; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   offsets[i],offsets[i+1],offsets[i+2],offsets[i+3],
	   offsets[i+4],offsets[i+5],offsets[i+6],offsets[i+7]);
  }
  printf("end of offsets\n");
#endif

  return;
}
#endif

