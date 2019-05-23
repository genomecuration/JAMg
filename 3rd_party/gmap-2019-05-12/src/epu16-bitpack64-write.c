static char rcsid[] = "$Id: epu16-bitpack64-write.c 216773 2018-10-03 21:27:20Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "epu16-bitpack64-write.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"		/* For FWRITE_USHORTS */
#else
#include "littleendian.h"	/* For FWRITE_USHORTS */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include "mem.h"
#include "assert.h"
#include "fopen.h"
#include "popcount.h"
#include "epu16-bitpack64-access.h"	/* For Epu16_bitpack64_extract */

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif


#define CHECK 1

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* #define ALLOW_ODD_PACKSIZES 1 */

/* #define USE_ONE_FILE_FOR_FIXED 1 */

#define DIFFERENTIAL_METAINFO_SIZE 2
#define PAIRED_METAINFO_SIZE 3
#define DIRECT_METAINFO_SIZE 1
#define BLOCKSIZE 64

#define BUFFER_SIZE 1000000



#ifdef HAVE_SSE2
static int
write_reg_buffered_vert (FILE *strm_fp, UINT2 *strm_buffer,
			 int strm_buffer_size, int strm_buffer_i, __m128i OutReg) {

  debug(printf("Entered write_reg_buffered_vert, SSE2\n"));

#if 0
  /* Type casting method (when we passed in pointer to OutReg).  Needs a memory fence. */
  UINT4 *buffer = (UINT4 *) OutReg;
  _mm_lfence();  /* Needed to avoid storing incorrect values into strm_buffer */
#else
  /* Storing method.  Safer.  */
  UINT2 buffer[8];
  _mm_store_si128((__m128i *) buffer,OutReg);
#endif

  /* printf("Writing %08X %08X %08X %08X\n",buffer[0],buffer[1],buffer[2],buffer[3]); */

  strm_buffer[strm_buffer_i++] = buffer[0];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[1];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[2];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[3];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[4];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[5];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[6];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  strm_buffer[strm_buffer_i++] = buffer[7];
  if (strm_buffer_i == strm_buffer_size) {
    FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
    strm_buffer_i = 0;
  }

  return strm_buffer_i;
}
#else
static int
write_reg_buffered_vert (FILE *strm_fp, UINT2 *strm_buffer,
			 int strm_buffer_size, int strm_buffer_i,
			 UINT2 *horizontal, int nwritten) {
  UINT2 vertical[64];
  int nrows = nwritten/8, row, column, k;

  debug(printf("Entered write_reg_buffered_vert, non-SIMD\n"));

  /* Convert to vertical */
  for (column = 0; column < 8; column++) {
    k = column;
    for (row = 0; row < nrows; row++) {
      vertical[k] = *horizontal++;
      k += 8;
    }
  }
    
#ifdef DEBUG
  for (k = 0; k < nwritten; k += 8) {
    printf("%04X %04X %04X %04X %04X %04X %04X %04X\n",
	   vertical[k+7],vertical[k+6],vertical[k+5],vertical[k+4],
	   vertical[k+3],vertical[k+2],vertical[k+1],vertical[k+0]);
  }
#endif

  /* Send to output buffer */
  for (k = 0; k < nwritten; k++) {
    strm_buffer[strm_buffer_i++] = vertical[k];
    if (strm_buffer_i == strm_buffer_size) {
      FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
      strm_buffer_i = 0;
    }
  }

  return strm_buffer_i;
}
#endif



#if 0
static int
write_reg_buffered_horiz (FILE *strm_fp, UINT2 *strm_buffer,
			  int strm_buffer_size, int strm_buffer_i,
			  UINT2 *values, int nwritten) {
  int k;

  /* Send to output buffer */
  for (k = 0; k < nwritten; k++) {
    /* printf("Writing %08X\n",values[k]); */
    strm_buffer[strm_buffer_i++] = values[k];
    if (strm_buffer_i == strm_buffer_size) {
      FWRITE_USHORTS(strm_buffer,strm_buffer_size,strm_fp);
      strm_buffer_i = 0;
    }
  }

  return strm_buffer_i;
}
#endif




#ifdef HAVE_SSE2
static __m128i mask1, mask2, mask3, mask4, mask5, mask6, mask7, mask8,
  mask9, mask10, mask11, mask12, mask13, mask14, mask15;
#endif


static void
write_setup () {

#ifdef HAVE_SSE2
  mask1 = _mm_set1_epi16(1U);
  mask2 = _mm_set1_epi16(3U);
  mask3 =  _mm_set1_epi16(7U);
  mask4 =  _mm_set1_epi16(15U);
  mask5 =  _mm_set1_epi16(31U);
  mask6 =  _mm_set1_epi16(63U);
  mask7 =  _mm_set1_epi16(127U);
  mask8 =  _mm_set1_epi16(255U);
  mask9 =  _mm_set1_epi16(511U);
  mask10 =  _mm_set1_epi16(1023U);
  mask11 =  _mm_set1_epi16(2047U);
  mask12 =  _mm_set1_epi16(4095U);
  mask13 =  _mm_set1_epi16(8191U);
  mask14 =  _mm_set1_epi16(16383U);
  mask15 =  _mm_set1_epi16(32767U);
#endif

  return;
}



#ifdef HAVE_SSE2
/* nwritten = 1 * 8 = 8 unsigned ints */
static int
write_02_vert (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT2 *_in) {
  const __m128i *in = (const __m128i *) _in;
  __m128i OutReg;

  debug(printf("Entered write_02_vert\n"));

  __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask2);
  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 2));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 4));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 6));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 10));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 12));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 14));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  return strm_buffer_i;
}
#endif

static int
pack_02_horiz (UINT2 *out, const UINT2 *in) {
  int column;

  debug(printf("Entered pack_02_horiz\n"));

  for (column = 0; column < 8; column++) {
    *out |= (*in)   % (1U << 2 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  14 ;
    ++out;
    ++in;
  }

  return 8;
}


#ifdef HAVE_SSE2
/* nwritten = 2 * 8 = 16 unsigned shorts */
static int
write_04_vert (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT2 *_in) {
  const __m128i *in = (const __m128i *) _in;
  __m128i OutReg;

  debug(printf("Entered write_04_vert\n"));

  __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask4);
  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 4));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 12));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 4));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 12));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  return strm_buffer_i;
}
#endif


static int
pack_04_horiz (UINT2 *out, const UINT2 *in) {
  int column;

  debug(printf("Entered pack_04_horiz\n"));

  for (column = 0; column < 8; column++) {
    *out |= (*in)   % (1U << 4 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  12 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 4 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  12 ;
    ++out;
    ++in;
  }

  return 16;
}



#ifdef HAVE_SSE2
/* nwritten = 3 * 8 = 24 unsigned shorts */
static int
write_06_vert (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT2 *_in) {
  const __m128i *in = (const __m128i *) _in;
  __m128i OutReg;

  debug(printf("Entered write_06_vert\n"));

  __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask6);
  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 6));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 12));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 6 - 2);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 2));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 14));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 6 - 4);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 4));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 10));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  return strm_buffer_i;
}
#endif

static int
pack_06_horiz (UINT2 *out, const UINT2 *in) {
  int column;

  debug(printf("Entered pack_06_horiz\n"));

  for (column = 0; column < 8; column++) {

    *out |= (*in)   % (1U << 6 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 6 ) ) >> ( 6  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 6 ) ) >> ( 6  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  10 ;
    ++out;
    ++in;
  }

  return 24;
}




#ifdef HAVE_SSE2
/* nwritten = 4 * 8 = 32 unsigned shorts */
static int
write_08_vert (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT2 *_in) {
  const __m128i *in = (const __m128i *) _in;
  __m128i OutReg;

  __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask8);

  debug(printf("Entered write_08_vert\n"));

  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  return strm_buffer_i;
}
#endif

static int
pack_08_horiz (UINT2 *out, const UINT2 *in) {
  int column;

  debug(printf("Entered pack_08_horiz\n"));

  for (column = 0; column < 8; column++) {
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++out;
    ++in;
  }

  return 32;
}




#ifdef HAVE_SSE2
/* nwritten = 4 * 8 = 40 unsigned shorts */
static int
write_10_vert (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT2 *_in) {
  const __m128i *in = (const __m128i *) _in;
  __m128i OutReg;

  debug(printf("Entered write_10_vert\n"));

  __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask10);
  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 10));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 10 - 4);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 4));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 14));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 10 - 8);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 10 - 2);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 2));
  InReg = _mm_and_si128(_mm_load_si128(++in), mask10);
  
  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 12));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 10 - 6);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 6));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);
  return strm_buffer_i;
}
#endif


static int
pack_10_horiz (UINT2 *out, const UINT2 *in) {
  int column;

  debug(printf("Entered pack_10_horiz\n"));

  for (column = 0; column < 8; column++) {
    *out |= (*in)   % (1U << 10 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  10 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  6 ;
    ++out;
    ++in;
  }

  return 40;
}




#ifdef HAVE_SSE2
/* nwritten = 6 * 8 = 48 unsigned shorts */
static int
write_12_vert (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT2 *_in) {
  const __m128i *in = (const __m128i *) _in;
  __m128i OutReg;

  debug(printf("Entered write_12_vert\n"));

  __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask12);
  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 12));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 12 - 8);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 12 - 4);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 4));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 12));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 12 - 8);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 12 - 4);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 4));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  return strm_buffer_i;
}
#endif

static int
pack_12_horiz (UINT2 *out, const UINT2 *in) {
  int column;

  debug(printf("Entered pack_12_horiz\n"));

  for (column = 0; column < 8; column++) {
    *out |= (*in)   % (1U << 12 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  4 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 12 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  4 ;
    ++out;
    ++in;
  }

  return 48;
}



#ifdef HAVE_SSE2
/* nwritten = 7 * 8 = 56 unsigned shorts */
static int
write_14_vert (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT2 *_in) {
  const __m128i *in = (const __m128i *) _in;
  __m128i OutReg;

  debug(printf("Entered write_14_vert\n"));

  __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask14);
  OutReg = InReg;
  InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 14));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 14 - 12);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 12));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 14 - 10);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 10));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 14 - 8);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 8));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 14 - 6);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 6));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 14 - 4);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 4));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  OutReg = _mm_srli_epi16(InReg, 14 - 2);
  InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

  OutReg =  _mm_or_si128(OutReg,_mm_slli_epi16(InReg, 2));
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  return strm_buffer_i;
}
#endif


static int
pack_14_horiz (UINT2 *out, const UINT2 *in) {
  int column;

  debug(printf("Entered pack_14_horiz\n"));

  for (column = 0; column < 8; column++) {
    *out |= (*in)   % (1U << 14 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  10 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  6 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  4 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  2 ;
    ++out;
    ++in;
  }

  return 56;
}


#ifdef HAVE_SSE2
/* nwritten = 8 * 8 = 64 unsigned shorts */
static int
write_16_vert (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i, const UINT2 *_in) {
  const __m128i *in = (const __m128i *) _in;
  __m128i OutReg;

  debug(printf("Entered write_16_vert\n"));

  __m128i InReg = _mm_load_si128(in);
  OutReg = InReg;
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_load_si128(++in);

  OutReg = InReg;
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_load_si128(++in);

  OutReg = InReg;
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_load_si128(++in);

  OutReg = InReg;
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_load_si128(++in);

  OutReg = InReg;
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_load_si128(++in);

  OutReg = InReg;
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_load_si128(++in);

  OutReg = InReg;
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  InReg = _mm_load_si128(++in);

  OutReg = InReg;
  strm_buffer_i = write_reg_buffered_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,
					  OutReg);

  return strm_buffer_i;
}
#endif


static int
pack_16_horiz (UINT2 *out, const UINT2 *in) {
  int column;

  debug(printf("Entered pack_16_horiz\n"));

  for (column = 0; column < 8; column++) {
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
  }

  return 64;
}



#ifdef HAVE_SSE2

/* Columnar order allows just the necessary values in a block to be decoded */
static void
columnar_order (UINT2 *columnar, const UINT2 *vertical) {

  columnar[0] = vertical[0];
  columnar[1] = vertical[8];
  columnar[2] = vertical[16];
  columnar[3] = vertical[24];
  columnar[4] = vertical[32];
  columnar[5] = vertical[40];
  columnar[6] = vertical[48];
  columnar[7] = vertical[56];

  columnar[8] = vertical[1];
  columnar[9] = vertical[9];
  columnar[10] = vertical[17];
  columnar[11] = vertical[25];
  columnar[12] = vertical[33];
  columnar[13] = vertical[41];
  columnar[14] = vertical[49];
  columnar[15] = vertical[57];

  columnar[16] = vertical[2];
  columnar[17] = vertical[10];
  columnar[18] = vertical[18];
  columnar[19] = vertical[26];
  columnar[20] = vertical[34];
  columnar[21] = vertical[42];
  columnar[22] = vertical[50];
  columnar[23] = vertical[58];

  columnar[24] = vertical[3];
  columnar[25] = vertical[11];
  columnar[26] = vertical[19];
  columnar[27] = vertical[27];
  columnar[28] = vertical[35];
  columnar[29] = vertical[43];
  columnar[30] = vertical[51];
  columnar[31] = vertical[59];

  columnar[32] = vertical[4];
  columnar[33] = vertical[12];
  columnar[34] = vertical[20];
  columnar[35] = vertical[28];
  columnar[36] = vertical[36];
  columnar[37] = vertical[44];
  columnar[38] = vertical[52];
  columnar[39] = vertical[60];

  columnar[40] = vertical[5];
  columnar[41] = vertical[13];
  columnar[42] = vertical[21];
  columnar[43] = vertical[29];
  columnar[44] = vertical[37];
  columnar[45] = vertical[45];
  columnar[46] = vertical[53];
  columnar[47] = vertical[61];

  columnar[48] = vertical[6];
  columnar[49] = vertical[14];
  columnar[50] = vertical[22];
  columnar[51] = vertical[30];
  columnar[52] = vertical[38];
  columnar[53] = vertical[46];
  columnar[54] = vertical[54];
  columnar[55] = vertical[62];

  columnar[56] = vertical[7];
  columnar[57] = vertical[15];
  columnar[58] = vertical[23];
  columnar[59] = vertical[31];
  columnar[60] = vertical[39];
  columnar[61] = vertical[47];
  columnar[62] = vertical[55];
  columnar[63] = vertical[63];

  return;
}

static int
write_columnar (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i,
		const UINT2 *_in, int packsize) {
  UINT2 columnar[BLOCKSIZE];

#ifdef DEBUG
  int i;

  printf("Entering write_columnar, SSE2, with packsize %d\n",packsize);
  for (i = 0; i < BLOCKSIZE; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  columnar_order(columnar,_in);

  switch (packsize) {
  case 0: return strm_buffer_i;
  case 2: return write_02_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 4: return write_04_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 6: return write_06_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 8: return write_08_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 10: return write_10_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 12: return write_12_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 14: return write_14_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);
  case 16: return write_16_vert(strm_fp,strm_buffer,strm_buffer_size,strm_buffer_i,columnar);

  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }
}

#else

/* Non-SIMD code cannot write vertical format easily, so using
   horizontal code and conversions */

static int
write_columnar (FILE *strm_fp, UINT2 *strm_buffer, int strm_buffer_size, int strm_buffer_i,
		const UINT2 *horizontal, int packsize) {
  int nwritten;
  UINT2 buffer[BLOCKSIZE];
#ifdef DEBUG
  int i;
#endif
  /* UINT2 vertical[BLOCKSIZE]; */
  /* UINT2 columnar[BLOCKSIZE]; */

#ifdef DEBUG
  printf("Entering write_columnar, non-SIMD, with packsize %d\n",packsize);
  for (i = 0; i < BLOCKSIZE; i++) {
    printf("%d ",horizontal[i]);
  }
  printf("\n");
#endif

#if 0
  /* No need since we there are no issues with bidirectional formats */
  columnar_order(columnar,horizontal);
  reorder_values_vertically(vertical,columnar);
#endif
  memset((void *) buffer,0,BLOCKSIZE*sizeof(UINT2));

#if 0
  printf("diffs in horizontal order\n");
  for (i = 0; i < 64; i += 8) {
    printf("%u %u %u %u %u %u %u %u\n",
	   horizontal[i],horizontal[i+1],horizontal[i+2],horizontal[i+3],
	   horizontal[i+4],horizontal[i+5],horizontal[i+6],horizontal[i+7]);
  }
#endif

  switch (packsize) {
  case 0: return strm_buffer_i;
  case 2: nwritten = pack_02_horiz(buffer,&(horizontal[0])); break;
  case 4: nwritten = pack_04_horiz(buffer,&(horizontal[0])); break;
  case 6: nwritten = pack_06_horiz(buffer,&(horizontal[0])); break;
  case 8: nwritten = pack_08_horiz(buffer,&(horizontal[0])); break;
  case 10: nwritten = pack_10_horiz(buffer,&(horizontal[0])); break;
  case 12: nwritten = pack_12_horiz(buffer,&(horizontal[0])); break;
  case 14: nwritten = pack_14_horiz(buffer,&(horizontal[0])); break;
  case 16: nwritten = pack_16_horiz(buffer,&(horizontal[0])); break;
  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }

  return write_reg_buffered_vert(strm_fp,strm_buffer,
				 strm_buffer_size,strm_buffer_i,
				 buffer,nwritten);
}

#endif


/* Processes 64 values at a time.  Returns packsize. */
int
compute_q8_diffs (UINT2 *diffs, UINT2 *values) {
  UINT2 packsize;
  int i;
  UINT2 maxdiff = 0;
  int firstbit;
#ifdef HAVE_BUILTIN_CLZ
#elif defined(HAVE_ASM_BSR)
  int msb;
#endif

#ifdef CHECK
  for (i = 0; i < 64; i++) {
    assert(values[i+1] >= values[i]);
  }
#endif

  for (i = 63; i >= 8; i--) {
    maxdiff |= (diffs[i] = values[i+1] - values[i+1-8]);
  }
  maxdiff |= (diffs[7] = values[8] - values[0]);
  maxdiff |= (diffs[6] = values[7] - values[0]);
  maxdiff |= (diffs[5] = values[6] - values[0]);
  maxdiff |= (diffs[4] = values[5] - values[0]);
  maxdiff |= (diffs[3] = values[4] - values[0]);
  maxdiff |= (diffs[2] = values[3] - values[0]);
  maxdiff |= (diffs[1] = values[2] - values[0]);
  maxdiff |= (diffs[0] = values[1] - values[0]);

  if (maxdiff == 0) {
    /* __builtin_clz() behaves oddly on zero */
    return 0;

  } else {
#ifdef HAVE_BUILTIN_CLZ
    firstbit = __builtin_clz(maxdiff);
    packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
    asm("bsr %1,%0" : "=r"(msb) : "r"(maxdiff));
    packsize = msb + 1;
#else
    firstbit = ((maxdiff >> 16) ? clz_table[maxdiff >> 16] : 16 + clz_table[maxdiff]);
    packsize = 32 - firstbit;
#endif

#ifdef ALLOW_ODD_PACKSIZES
    return packsize;
#else
    return (packsize + 1) & ~1;	/* Converts packsizes to the next multiple of 2 */
#endif
  }
}


static UINT2
compute_ascending (UINT2 *ascending, UINT2 *counts) {
  int i;

  ascending[0] = 0;
  for (i = 1; i <= 64; i++) {
    ascending[i] = ascending[i-1] + counts[i-1];
  }

  return ascending[64];
}

/* We want to store values 0..n, with final value at ascending[n]
   possibly stored as the final metainfo value */
/* Stored in columnar order */
UINT4
Epu16_bitpack64_append_differential (UINT4 *totalcount, FILE *ptrs_fp, FILE *comp_fp,
				     char *packsizes, UINT2 **bitpacks, Localspace_T n) {
  UINT4 nregisters;
  UINT2 *ptrs, *p;
  size_t nptrs;
  Localspace_T positioni, bmer;

  /* Buffer is used to avoid frequent writes to the file */
  UINT2 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT2 diffs[BLOCKSIZE], ascending[BLOCKSIZE+1], counts[BLOCKSIZE];
  /* UINT2 last_block[BLOCKSIZE]; */
  int packsize;


  write_setup();

  debug(printf("Entered Epu16_bitpack64_append_differential with n %llu\n",n));

  /* 2 metavalues: nwritten (pointer) and cumulative sum for block.
     Packsize can be computed from difference between successive
     pointers, if only even packsizes are allowed */
  p = ptrs = (UINT2 *) CALLOC(((n + BLOCKSIZE)/BLOCKSIZE + 1) * DIFFERENTIAL_METAINFO_SIZE,sizeof(UINT2));

  buffer = (UINT2 *) CALLOC(buffer_size,sizeof(UINT2));
  buffer_i = 0;

  nregisters = 0;

  /* Last value of ascending is at ascending[n] */
  /* Use <= n instead of < n, because we do want to write the final pointers */
  *totalcount = 0;
  for (positioni = 0, bmer = 0; positioni + BLOCKSIZE <= n; positioni += BLOCKSIZE, bmer++) {
    debug(printf("Looping at positioni %d out of %d\n",positioni,n));

    /* Pack block of 64 diffs */
    Epu16_bitpack64_extract(counts,packsizes[bmer],bitpacks[bmer]);
    *totalcount += (UINT4) compute_ascending(ascending,counts);
    packsize = compute_q8_diffs(diffs,ascending); /* Note: This packsize may differ from packsizes[bmer] */
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);
    nregisters += packsize / 2;

    /* Do not write first "0 0".  Different from Bitpack64_write_differential */
    /* Pointer */
    *p++ = (UINT2) nregisters;	/* In 128-bit registers */

    /* Value for start of block */
    *p++ = (UINT2) (*totalcount);
  }

#if 0
  /* Old: Check for positioni < n, because if positioni == n, ascending[n] will be taken care of as metainfo */
  /* Use <= n instead of < n, because we do want to write final pointers */
  /* For nucleotides, expect a single final block where positioni == n */
  debug(printf("Deciding whether to finish last block.  positioni is %d\n",positioni));
  if (positioni <= n) {
    debug(printf("Finishing last block of 64, because positioni is %d < n %d\n",positioni,n));

#if 0
    /* Finish last block of 64 */
    *p++ = (UINT2) nregisters;	/* In 128-bit registers */

    /* Value for start of block */
    *p++ = (UINT2) (*totalcount);
#endif

    if (positioni == n) {
      /* Don't have a bitpack at [bmerspace].  Just fills counts with zeroes. */
      Epu16_bitpack64_extract(counts,/*packsize*/0,/*bitpack*/NULL);
    } else {
      Epu16_bitpack64_extract(counts,packsizes[bmer],bitpacks[bmer]);
    }

    /* For differential, want <=.  For direct, want < */
    for (i = 0; i <= (int) (n - positioni); i++) {
      last_block[i] = counts[i];
    }
    for ( ; i < BLOCKSIZE; i++) {
      /* Copy last value for rest of block */
      last_block[i] = 0;
    }

    /* Pack block of < 64 diffs */
    *totalcount += (UINT4) compute_ascending(ascending,last_block);
    packsize = compute_q8_diffs(diffs,ascending);
    buffer_i = write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

    nregisters += packsize / 2;

    /* Write final pointers */
    *p++ = (UINT2) nregisters;	/* In 128-bit registers */

    /* Value for end of block */
    *p++ = (UINT2) *totalcount;
  }
#endif

  nptrs = p - ptrs;
  FWRITE_USHORTS(ptrs,nptrs,ptrs_fp);
  FREE(ptrs);
    
  /* Empty buffer */
  if (buffer_i > 0) {
    FWRITE_USHORTS(buffer,buffer_i,comp_fp);	
    buffer_i = 0;
  }
  FREE(buffer);

  return nregisters;
}


