static char rcsid[] = "$Id: bitpack64-access.c 180341 2015-12-07 18:29:40Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bitpack64-access.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#define CONVERT(x) Bigendian_convert_uint(x)
#else
#define CONVERT(x) x
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define BLOCKSIZE 64

/* Vertical access is slightly more efficient than horizontal */

#ifdef HORIZONTAL
#define WORD_INCR 1		/* 1 for horizontal; 4 for vertical */
#else
#define WORD_INCR 4
#endif


static UINT4
access_00 (const UINT4 *in) {
  return 0U;
}



static UINT4
access_02_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 2 ) ;
}

static UINT4
access_02_01 (const UINT4 *in) {
  return ( CONVERT(*in) >>  2  )   % (1U << 2 ) ;
}

static UINT4
access_02_02 (const UINT4 *in) {
  return ( CONVERT(*in) >>  4  )   % (1U << 2 ) ;
}

static UINT4
access_02_03 (const UINT4 *in) {
  return ( CONVERT(*in) >>  6  )   % (1U << 2 ) ;
}

static UINT4
access_02_04 (const UINT4 *in) {
  return ( CONVERT(*in) >>  8  )   % (1U << 2 ) ;
}

static UINT4
access_02_05 (const UINT4 *in) {
  return ( CONVERT(*in) >>  10  )   % (1U << 2 ) ;
}

static UINT4
access_02_06 (const UINT4 *in) {
  return ( CONVERT(*in) >>  12  )   % (1U << 2 ) ;
}

static UINT4
access_02_07 (const UINT4 *in) {
  return ( CONVERT(*in) >>  14  )   % (1U << 2 ) ;
}

static UINT4
access_02_08 (const UINT4 *in) {
  return ( CONVERT(*in) >>  16  )   % (1U << 2 ) ;
}

static UINT4
access_02_09 (const UINT4 *in) {
  return ( CONVERT(*in) >>  18  )   % (1U << 2 ) ;
}

static UINT4
access_02_10 (const UINT4 *in) {
  return ( CONVERT(*in) >>  20  )   % (1U << 2 ) ;
}

static UINT4
access_02_11 (const UINT4 *in) {
  return ( CONVERT(*in) >>  22  )   % (1U << 2 ) ;
}

static UINT4
access_02_12 (const UINT4 *in) {
  return ( CONVERT(*in) >>  24  )   % (1U << 2 ) ;
}

static UINT4
access_02_13 (const UINT4 *in) {
  return ( CONVERT(*in) >>  26  )   % (1U << 2 ) ;
}

static UINT4
access_02_14 (const UINT4 *in) {
  return ( CONVERT(*in) >>  28  )   % (1U << 2 ) ;
}

static UINT4
access_02_15 (const UINT4 *in) {
  return ( CONVERT(*in) >>  30  )   % (1U << 2 ) ;
}



static UINT4
access_04_00 (const UINT4 *in) {
  return ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;
}

static UINT4
access_04_01 (const UINT4 *in) {
  return ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;
}

static UINT4
access_04_02 (const UINT4 *in) {
  return ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;
}

static UINT4
access_04_03 (const UINT4 *in) {
  return ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;
}

static UINT4
access_04_04 (const UINT4 *in) {
  return ( CONVERT(*in) >> 16 )   % (1U << 4 ) ;
}

static UINT4
access_04_05 (const UINT4 *in) {
  return ( CONVERT(*in) >> 20 )   % (1U << 4 ) ;
}

static UINT4
access_04_06 (const UINT4 *in) {
  return ( CONVERT(*in) >> 24 )   % (1U << 4 ) ;
}

static UINT4
access_04_07 (const UINT4 *in) {
  return ( CONVERT(*in) >> 28 )   % (1U << 4 ) ;
}

static UINT4
access_04_08 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;
}

static UINT4
access_04_09 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;
}

static UINT4
access_04_10 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;
}

static UINT4
access_04_11 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;
}

static UINT4
access_04_12 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 4 ) ;
}

static UINT4
access_04_13 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 20 )   % (1U << 4 ) ;
}

static UINT4
access_04_14 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 24 )   % (1U << 4 ) ;
}

static UINT4
access_04_15 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 28 )   % (1U << 4 ) ;
}


static UINT4
access_06_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 6 ) ;
}

static UINT4
access_06_01 (const UINT4 *in) {
  return ( CONVERT(*in) >>  6  )   % (1U << 6 ) ;
}

static UINT4
access_06_02 (const UINT4 *in) {
  return ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;
}

static UINT4
access_06_03 (const UINT4 *in) {
  return ( CONVERT(*in) >>  18  )   % (1U << 6 ) ;
}

static UINT4
access_06_04 (const UINT4 *in) {
  return ( CONVERT(*in) >>  24  )   % (1U << 6 ) ;
}

static UINT4
access_06_05 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  30  )   % (1U << 6 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 );
  return out;
}

static UINT4
access_06_06 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 6 ) ;
}

static UINT4
access_06_07 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  10  )   % (1U << 6 ) ;
}

static UINT4
access_06_08 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  16  )   % (1U << 6 ) ;
}

static UINT4
access_06_09 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  22  )   % (1U << 6 ) ;
}

static UINT4
access_06_10 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 6 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 );
  return out;
}

static UINT4
access_06_11 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 6 ) ;
}

static UINT4
access_06_12 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 6 ) ;
}

static UINT4
access_06_13 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;
}

static UINT4
access_06_14 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  20  )   % (1U << 6 ) ;
}

static UINT4
access_06_15 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  26  )   % (1U << 6 ) ;
}


static UINT4
access_08_00 (const UINT4 *in) {
  return ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
}

static UINT4
access_08_01 (const UINT4 *in) {
  return ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
}

static UINT4
access_08_02 (const UINT4 *in) {
  return ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
}

static UINT4
access_08_03 (const UINT4 *in) {
  return ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
}

static UINT4
access_08_04 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
}

static UINT4
access_08_05 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
}

static UINT4
access_08_06 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
}

static UINT4
access_08_07 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
}

static UINT4
access_08_08 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
}

static UINT4
access_08_09 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
}

static UINT4
access_08_10 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
}

static UINT4
access_08_11 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
}

static UINT4
access_08_12 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
}

static UINT4
access_08_13 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
}

static UINT4
access_08_14 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
}

static UINT4
access_08_15 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
}


static UINT4
access_10_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 10 ) ;
}

static UINT4
access_10_01 (const UINT4 *in) {
  return ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;
}

static UINT4
access_10_02 (const UINT4 *in) {
  return ( CONVERT(*in) >>  20  )   % (1U << 10 ) ;
}

static UINT4
access_10_03 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  30  )   % (1U << 10 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 );
  return out;
}

static UINT4
access_10_04 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;
}

static UINT4
access_10_05 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  18  )   % (1U << 10 ) ;
}

static UINT4
access_10_06 (const UINT4 *in) {
  UINT4 out;
  
  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 10 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 );
  return out;
}

static UINT4
access_10_07 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  6  )   % (1U << 10 ) ;
}

static UINT4
access_10_08 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  16  )   % (1U << 10 ) ;
}

static UINT4
access_10_09 (const UINT4 *in) {
  UINT4 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  26  )   % (1U << 10 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 );
  return out;
}

static UINT4
access_10_10 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 10 ) ;
}

static UINT4
access_10_11 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;
}

static UINT4
access_10_12 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 10 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 );
  return out;
}

static UINT4
access_10_13 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 10 ) ;
}

static UINT4
access_10_14 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;
}

static UINT4
access_10_15 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  22  )   % (1U << 10 ) ;
}


static UINT4
access_12_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
}

static UINT4
access_12_01 (const UINT4 *in) {
  return ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
}

static UINT4
access_12_02 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  return out;
}

static UINT4
access_12_03 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
}

static UINT4
access_12_04 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  16  )   % (1U << 12 ) ;
}

static UINT4
access_12_05 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  return out;
}

static UINT4
access_12_06 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
}

static UINT4
access_12_07 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  20  )   % (1U << 12 ) ;
}

static UINT4
access_12_08 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
}

static UINT4
access_12_09 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
}

static UINT4
access_12_10 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  return out;
}

static UINT4
access_12_11 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
}

static UINT4
access_12_12 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  16  )   % (1U << 12 ) ;
}

static UINT4
access_12_13 (const UINT4 *in) {
  UINT4 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  return out;
}

static UINT4
access_12_14 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
}

static UINT4
access_12_15 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >>  20  )   % (1U << 12 ) ;
}


static UINT4
access_14_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 14 ) ;
}

static UINT4
access_14_01 (const UINT4 *in) {
  return ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;
}

static UINT4
access_14_02 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  28  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 );
  return out;
}

static UINT4
access_14_03 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;
}

static UINT4
access_14_04 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 );
  return out;
}

static UINT4
access_14_05 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;
}

static UINT4
access_14_06 (const UINT4 *in) {
  UINT4 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  20  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 );
  return out;
}

static UINT4
access_14_07 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 14 ) ;
}

static UINT4
access_14_08 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  16  )   % (1U << 14 ) ;
}

static UINT4
access_14_09 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  30  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 );
  return out;
}

static UINT4
access_14_10 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;
}

static UINT4
access_14_11 (const UINT4 *in) {
  UINT4 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  26  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 );
  return out;
}

static UINT4
access_14_12 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;
}

static UINT4
access_14_13 (const UINT4 *in) {
  UINT4 out;

  in += 5 * WORD_INCR;
  out = ( CONVERT(*in) >>  22  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 );
  return out;
}

static UINT4
access_14_14 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;
}

static UINT4
access_14_15 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >>  18  )   % (1U << 14 ) ;
}


static UINT4
access_16_00 (const UINT4 *in) {
  return ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
}

static UINT4
access_16_01 (const UINT4 *in) {
  return ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
}

static UINT4
access_16_02 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
}

static UINT4
access_16_03 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
}

static UINT4
access_16_04 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
}

static UINT4
access_16_05 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
}

static UINT4
access_16_06 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
}

static UINT4
access_16_07 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
}

static UINT4
access_16_08 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
}

static UINT4
access_16_09 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
}

static UINT4
access_16_10 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
}

static UINT4
access_16_11 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
}

static UINT4
access_16_12 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
}

static UINT4
access_16_13 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
}

static UINT4
access_16_14 (const UINT4 *in) {
  in += 7 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
}

static UINT4
access_16_15 (const UINT4 *in) {
  in += 7 * WORD_INCR;
  return ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
}


static UINT4
access_18_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 18 ) ;
}

static UINT4
access_18_01 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  18  )   % (1U << 18 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 18 - 4 );
  return out;
}

static UINT4
access_18_02 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 18 ) ;
}

static UINT4
access_18_03 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  22  )   % (1U << 18 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 18 - 8 );
  return out;
}

static UINT4
access_18_04 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 18 ) ;
}

static UINT4
access_18_05 (const UINT4 *in) {
  UINT4 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  26  )   % (1U << 18 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 18 - 12 );
  return out;
}

static UINT4
access_18_06 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  12  )   % (1U << 18 ) ;
}

static UINT4
access_18_07 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  30  )   % (1U << 18 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 18 - 16 );
  return out;
}

static UINT4
access_18_08 (const UINT4 *in) {
  UINT4 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 18 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 18 - 2 );
  return out;
}

static UINT4
access_18_09 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 18 ) ;
}

static UINT4
access_18_10 (const UINT4 *in) {
  UINT4 out;

  in += 5 * WORD_INCR;
  out = ( CONVERT(*in) >>  20  )   % (1U << 18 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 6 ))<<( 18 - 6 );
  return out;
}

static UINT4
access_18_11 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >>  6  )   % (1U << 18 ) ;
}

static UINT4
access_18_12 (const UINT4 *in) {
  UINT4 out;

  in += 6 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 18 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 10 ))<<( 18 - 10 );
  return out;
}

static UINT4
access_18_13 (const UINT4 *in) {
  in += 7 * WORD_INCR;
  return ( CONVERT(*in) >>  10  )   % (1U << 18 ) ;
}

static UINT4
access_18_14 (const UINT4 *in) {
  UINT4 out;

  in += 7 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 18 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 14 ))<<( 18 - 14 );
  return out;
}

static UINT4
access_18_15 (const UINT4 *in) {
  in += 8 * WORD_INCR;
  return ( CONVERT(*in) >>  14  )   % (1U << 18 ) ;
}


static UINT4
access_20_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 20 ) ;
}

static UINT4
access_20_01 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 );
  return out;
}

static UINT4
access_20_02 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 20 ) ;
}

static UINT4
access_20_03 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 );
  return out;
}

static UINT4
access_20_04 (const UINT4 *in) {
  UINT4 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 );
  return out;
}

static UINT4
access_20_05 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 20 ) ;
}

static UINT4
access_20_06 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 );
  return out;
}

static UINT4
access_20_07 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  12  )   % (1U << 20 ) ;
}

static UINT4
access_20_08 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >>  0  )   % (1U << 20 ) ;
}

static UINT4
access_20_09 (const UINT4 *in) {
  UINT4 out;

  in += 5 * WORD_INCR;
  out = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 );
  return out;
}

static UINT4
access_20_10 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 20 ) ;
}

static UINT4
access_20_11 (const UINT4 *in) {
  UINT4 out;

  in += 6 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 );
  return out;
}

static UINT4
access_20_12 (const UINT4 *in) {
  UINT4 out;

  in += 7 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 );
  return out;
}

static UINT4
access_20_13 (const UINT4 *in) {
  in += 8 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 20 ) ;
}

static UINT4
access_20_14 (const UINT4 *in) {
  UINT4 out;

  in += 8 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 );
  return out;
}

static UINT4
access_20_15 (const UINT4 *in) {
  in += 9 * WORD_INCR;
  return ( CONVERT(*in) >>  12  )   % (1U << 20 ) ;
}


static UINT4
access_22_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 22 ) ;
}

static UINT4
access_22_01 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  22  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 22 - 12 );
  return out;
}

static UINT4
access_22_02 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  12  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 22 - 2 );
  return out;
}

static UINT4
access_22_03 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 22 ) ;
}

static UINT4
access_22_04 (const UINT4 *in) {
  UINT4 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 14 ))<<( 22 - 14 );
  return out;
}

static UINT4
access_22_05 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  14  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 22 - 4 );
  return out;
}

static UINT4
access_22_06 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 22 ) ;
}

static UINT4
access_22_07 (const UINT4 *in) {
  UINT4 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  26  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 22 - 16 );
  return out;
}

static UINT4
access_22_08 (const UINT4 *in) {
  UINT4 out;

  in += 5 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 6 ))<<( 22 - 6 );
  return out;
}

static UINT4
access_22_09 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >>  6  )   % (1U << 22 ) ;
}

static UINT4
access_22_10 (const UINT4 *in) {
  UINT4 out;

  in += 6 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 18 ))<<( 22 - 18 );
  return out;
}

static UINT4
access_22_11 (const UINT4 *in) {
  UINT4 out;

  in += 7 * WORD_INCR;
  out = ( CONVERT(*in) >>  18  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 22 - 8 );
  return out;
}

static UINT4
access_22_12 (const UINT4 *in) {
  in += 8 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 22 ) ;
}

static UINT4
access_22_13 (const UINT4 *in) {
  UINT4 out;

  in += 8 * WORD_INCR;
  out = ( CONVERT(*in) >>  30  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 20 ))<<( 22 - 20 );
  return out;
}

static UINT4
access_22_14 (const UINT4 *in) {
  UINT4 out;

  in += 9 * WORD_INCR;
  out = ( CONVERT(*in) >>  20  )   % (1U << 22 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 10 ))<<( 22 - 10 );
  return out;
}

static UINT4
access_22_15 (const UINT4 *in) {
  in += 10 * WORD_INCR;
  return ( CONVERT(*in) >>  10  )   % (1U << 22 ) ;
}



static UINT4
access_24_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
}

static UINT4
access_24_01 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  return out;
}

static UINT4
access_24_02 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  return out;
}

static UINT4
access_24_03 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
}

static UINT4
access_24_04 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
}

static UINT4
access_24_05 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  return out;
}

static UINT4
access_24_06 (const UINT4 *in) {
  UINT4 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  return out;
}

static UINT4
access_24_07 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
}

static UINT4
access_24_08 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
}

static UINT4
access_24_09 (const UINT4 *in) {
  UINT4 out;

  in += 6 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  return out;
}

static UINT4
access_24_10 (const UINT4 *in) {
  UINT4 out;

  in += 7 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  return out;
}

static UINT4
access_24_11 (const UINT4 *in) {
  in += 8 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
}

static UINT4
access_24_12 (const UINT4 *in) {
  in += 9 * WORD_INCR;
  return ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
}

static UINT4
access_24_13 (const UINT4 *in) {
  UINT4 out;

  in += 9 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  return out;
}

static UINT4
access_24_14 (const UINT4 *in) {
  UINT4 out;

  in += 10 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  return out;
}

static UINT4
access_24_15 (const UINT4 *in) {
  in += 11 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
}



static UINT4
access_26_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 26 ) ;
}

static UINT4
access_26_01 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  26  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 20 ))<<( 26 - 20 );
  return out;
}

static UINT4
access_26_02 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  20  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 14 ))<<( 26 - 14 );
  return out;
}

static UINT4
access_26_03 (const UINT4 *in) {
  UINT4 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  14  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 26 - 8 );
  return out;
}

static UINT4
access_26_04 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  8  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 26 - 2 );
  return out;
}

static UINT4
access_26_05 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 26 ) ;
}

static UINT4
access_26_06 (const UINT4 *in) {
  UINT4 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 22 ))<<( 26 - 22 );
  return out;
}

static UINT4
access_26_07 (const UINT4 *in) {
  UINT4 out;

  in += 5 * WORD_INCR;
  out = ( CONVERT(*in) >>  22  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 26 - 16 );
  return out;
}

static UINT4
access_26_08 (const UINT4 *in) {
  UINT4 out;

  in += 6 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 10 ))<<( 26 - 10 );
  return out;
}

static UINT4
access_26_09 (const UINT4 *in) {
  UINT4 out;

  in += 7 * WORD_INCR;
  out = ( CONVERT(*in) >>  10  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 26 - 4 );
  return out;
}

static UINT4
access_26_10 (const UINT4 *in) {
  in += 8 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 26 ) ;
}

static UINT4
access_26_11 (const UINT4 *in) {
  UINT4 out;

  in += 8 * WORD_INCR;
  out = ( CONVERT(*in) >>  30  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 24 ))<<( 26 - 24 );
  return out;
}

static UINT4
access_26_12 (const UINT4 *in) {
  UINT4 out;

  in += 9 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 18 ))<<( 26 - 18 );
  return out;
}

static UINT4
access_26_13 (const UINT4 *in) {
  UINT4 out;

  in += 10 * WORD_INCR;
  out = ( CONVERT(*in) >>  18  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 26 - 12 );
  return out;
}

static UINT4
access_26_14 (const UINT4 *in) {
  UINT4 out;

  in += 11 * WORD_INCR;
  out = ( CONVERT(*in) >>  12  )   % (1U << 26 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 6 ))<<( 26 - 6 );
  return out;
}

static UINT4
access_26_15 (const UINT4 *in) {
  in += 12 * WORD_INCR;
  return ( CONVERT(*in) >>  6  )   % (1U << 26 ) ;
}


static UINT4
access_28_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 28 ) ;
}

static UINT4
access_28_01 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 );
  return out;
}

static UINT4
access_28_02 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 );
  return out;
}

static UINT4
access_28_03 (const UINT4 *in) {
  UINT4 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 );
  return out;
}

static UINT4
access_28_04 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 );
  return out;
}

static UINT4
access_28_05 (const UINT4 *in) {
  UINT4 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 );
  return out;
}

static UINT4
access_28_06 (const UINT4 *in) {
  UINT4 out;

  in += 5 * WORD_INCR;
  out = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 );
  return out;
}

static UINT4
access_28_07 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 28 ) ;
}

static UINT4
access_28_08 (const UINT4 *in) {
  in += 7 * WORD_INCR;
  return ( CONVERT(*in) >>  0  )   % (1U << 28 ) ;
}

static UINT4
access_28_09 (const UINT4 *in) {
  UINT4 out;

  in += 7 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 );
  return out;
}

static UINT4
access_28_10 (const UINT4 *in) {
  UINT4 out;

  in += 8 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 );
  return out;
}

static UINT4
access_28_11 (const UINT4 *in) {
  UINT4 out;

  in += 9 * WORD_INCR;
  out = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 );
  return out;
}

static UINT4
access_28_12 (const UINT4 *in) {
  UINT4 out;

  in += 10 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 );
  return out;
}

static UINT4
access_28_13 (const UINT4 *in) {
  UINT4 out;

  in += 11 * WORD_INCR;
  out = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 );
  return out;
}

static UINT4
access_28_14 (const UINT4 *in) {
  UINT4 out;

  in += 12 * WORD_INCR;
  out = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 );
  return out;
}

static UINT4
access_28_15 (const UINT4 *in) {
  in += 13 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 28 ) ;
}


static UINT4
access_30_00 (const UINT4 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 30 ) ;
}

static UINT4
access_30_01 (const UINT4 *in) {
  UINT4 out;

  out = ( CONVERT(*in) >>  30  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 28 ))<<( 30 - 28 );
  return out;
}

static UINT4
access_30_02 (const UINT4 *in) {
  UINT4 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  28  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 26 ))<<( 30 - 26 );
  return out;
}

static UINT4
access_30_03 (const UINT4 *in) {
  UINT4 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  26  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 24 ))<<( 30 - 24 );
  return out;
}

static UINT4
access_30_04 (const UINT4 *in) {
  UINT4 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  24  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 22 ))<<( 30 - 22 );
  return out;
}

static UINT4
access_30_05 (const UINT4 *in) {
  UINT4 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  22  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 20 ))<<( 30 - 20 );
  return out;
}

static UINT4
access_30_06 (const UINT4 *in) {
  UINT4 out;

  in += 5 * WORD_INCR;
  out = ( CONVERT(*in) >>  20  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 18 ))<<( 30 - 18 );
  return out;
}

static UINT4
access_30_07 (const UINT4 *in) {
  UINT4 out;

  in += 6 * WORD_INCR;
  out = ( CONVERT(*in) >>  18  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 16 ))<<( 30 - 16 );
  return out;
}

static UINT4
access_30_08 (const UINT4 *in) {
  UINT4 out;

  in += 7 * WORD_INCR;
  out = ( CONVERT(*in) >>  16  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 14 ))<<( 30 - 14 );
  return out;
}

static UINT4
access_30_09 (const UINT4 *in) {
  UINT4 out;

  in += 8 * WORD_INCR;
  out = ( CONVERT(*in) >>  14  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 30 - 12 );
  return out;
}

static UINT4
access_30_10 (const UINT4 *in) {
  UINT4 out;

  in += 9 * WORD_INCR;
  out = ( CONVERT(*in) >>  12  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 10 ))<<( 30 - 10 );
  return out;
}

static UINT4
access_30_11 (const UINT4 *in) {
  UINT4 out;

  in += 10 * WORD_INCR;
  out = ( CONVERT(*in) >>  10  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 30 - 8 );
  return out;
}

static UINT4
access_30_12 (const UINT4 *in) {
  UINT4 out;

  in += 11 * WORD_INCR;
  out = ( CONVERT(*in) >>  8  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 6 ))<<( 30 - 6 );
  return out;
}

static UINT4
access_30_13 (const UINT4 *in) {
  UINT4 out;

  in += 12 * WORD_INCR;
  out = ( CONVERT(*in) >>  6  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 30 - 4 );
  return out;
}

static UINT4
access_30_14 (const UINT4 *in) {
  UINT4 out;

  in += 13 * WORD_INCR;
  out = ( CONVERT(*in) >>  4  )   % (1U << 30 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 30 - 2 );
  return out;
}

static UINT4
access_30_15 (const UINT4 *in) {
  in += 14 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 30 ) ;
}


static UINT4
access_32_00 (const UINT4 *in) {
  return CONVERT(*in);
}

static UINT4
access_32_01 (const UINT4 *in) {
  in += 1 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_02 (const UINT4 *in) {
  in += 2 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_03 (const UINT4 *in) {
  in += 3 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_04 (const UINT4 *in) {
  in += 4 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_05 (const UINT4 *in) {
  in += 5 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_06 (const UINT4 *in) {
  in += 6 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_07 (const UINT4 *in) {
  in += 7 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_08 (const UINT4 *in) {
  in += 8 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_09 (const UINT4 *in) {
  in += 9 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_10 (const UINT4 *in) {
  in += 10 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_11 (const UINT4 *in) {
  in += 11 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_12 (const UINT4 *in) {
  in += 12 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_13 (const UINT4 *in) {
  in += 13 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_14 (const UINT4 *in) {
  in += 14 * WORD_INCR;
  return CONVERT(*in);
}

static UINT4
access_32_15 (const UINT4 *in) {
  in += 15 * WORD_INCR;
  return CONVERT(*in);
}



typedef UINT4 (*Accessor_T) (const UINT4 *);

static Accessor_T accessor_table[272] =
  {access_00, access_00, access_00, access_00,
   access_00, access_00, access_00, access_00,
   access_00, access_00, access_00, access_00,
   access_00, access_00, access_00, access_00,

   access_02_00, access_02_01, access_02_02, access_02_03,
   access_02_04, access_02_05, access_02_06, access_02_07,
   access_02_08, access_02_09, access_02_10, access_02_11,
   access_02_12, access_02_13, access_02_14, access_02_15,

   access_04_00, access_04_01, access_04_02, access_04_03,
   access_04_04, access_04_05, access_04_06, access_04_07,
   access_04_08, access_04_09, access_04_10, access_04_11,
   access_04_12, access_04_13, access_04_14, access_04_15,

   access_06_00, access_06_01, access_06_02, access_06_03,
   access_06_04, access_06_05, access_06_06, access_06_07,
   access_06_08, access_06_09, access_06_10, access_06_11,
   access_06_12, access_06_13, access_06_14, access_06_15,

   access_08_00, access_08_01, access_08_02, access_08_03,
   access_08_04, access_08_05, access_08_06, access_08_07,
   access_08_08, access_08_09, access_08_10, access_08_11,
   access_08_12, access_08_13, access_08_14, access_08_15,

   access_10_00, access_10_01, access_10_02, access_10_03,
   access_10_04, access_10_05, access_10_06, access_10_07,
   access_10_08, access_10_09, access_10_10, access_10_11,
   access_10_12, access_10_13, access_10_14, access_10_15,

   access_12_00, access_12_01, access_12_02, access_12_03,
   access_12_04, access_12_05, access_12_06, access_12_07,
   access_12_08, access_12_09, access_12_10, access_12_11,
   access_12_12, access_12_13, access_12_14, access_12_15,

   access_14_00, access_14_01, access_14_02, access_14_03,
   access_14_04, access_14_05, access_14_06, access_14_07,
   access_14_08, access_14_09, access_14_10, access_14_11,
   access_14_12, access_14_13, access_14_14, access_14_15,

   access_16_00, access_16_01, access_16_02, access_16_03,
   access_16_04, access_16_05, access_16_06, access_16_07,
   access_16_08, access_16_09, access_16_10, access_16_11,
   access_16_12, access_16_13, access_16_14, access_16_15,

   access_18_00, access_18_01, access_18_02, access_18_03,
   access_18_04, access_18_05, access_18_06, access_18_07,
   access_18_08, access_18_09, access_18_10, access_18_11,
   access_18_12, access_18_13, access_18_14, access_18_15,

   access_20_00, access_20_01, access_20_02, access_20_03,
   access_20_04, access_20_05, access_20_06, access_20_07,
   access_20_08, access_20_09, access_20_10, access_20_11,
   access_20_12, access_20_13, access_20_14, access_20_15,

   access_22_00, access_22_01, access_22_02, access_22_03,
   access_22_04, access_22_05, access_22_06, access_22_07,
   access_22_08, access_22_09, access_22_10, access_22_11,
   access_22_12, access_22_13, access_22_14, access_22_15,

   access_24_00, access_24_01, access_24_02, access_24_03,
   access_24_04, access_24_05, access_24_06, access_24_07,
   access_24_08, access_24_09, access_24_10, access_24_11,
   access_24_12, access_24_13, access_24_14, access_24_15,

   access_26_00, access_26_01, access_26_02, access_26_03,
   access_26_04, access_26_05, access_26_06, access_26_07,
   access_26_08, access_26_09, access_26_10, access_26_11,
   access_26_12, access_26_13, access_26_14, access_26_15,

   access_28_00, access_28_01, access_28_02, access_28_03,
   access_28_04, access_28_05, access_28_06, access_28_07,
   access_28_08, access_28_09, access_28_10, access_28_11,
   access_28_12, access_28_13, access_28_14, access_28_15,

   access_30_00, access_30_01, access_30_02, access_30_03,
   access_30_04, access_30_05, access_30_06, access_30_07,
   access_30_08, access_30_09, access_30_10, access_30_11,
   access_30_12, access_30_13, access_30_14, access_30_15,

   access_32_00, access_32_01, access_32_02, access_32_03,
   access_32_04, access_32_05, access_32_06, access_32_07,
   access_32_08, access_32_09, access_32_10, access_32_11,
   access_32_12, access_32_13, access_32_14, access_32_15,
  };
  


#define DIRECT_METAINFO_SIZE 1

#ifdef HORIZONTAL

UINT4
Bitpack64_access (Oligospace_T oligo, UINT4 *ptrs, UINT4 *comp) {
  UINT4 *info, start;
  int nwritten, remainder;
  UINT4 *bitpack;
  int index, row;
#ifdef DEBUG
  int packsize, i;
#endif

  info = &(ptrs[oligo/BLOCKSIZE * DIRECT_METAINFO_SIZE]);

#ifdef WORDS_BIGENDIAN
  start = Bigendian_convert_uint(info[0]);
  bitpack = (UINT4 *) &(comp[start*4]);
  nwritten = Bigendian_convert_uint(info[1]) - start;	/* In 128-bit registers */
#else
  start = info[0];
  bitpack = (UINT4 *) &(comp[start*4]);
  nwritten = info[1] - start;	/* In 128-bit registers */
#endif

  remainder = oligo % BLOCKSIZE;
  index = nwritten*16 + remainder % 16;
  row = (remainder / 16) * (packsize / 2);   /* Complexity of this calculation makes horizontal format slower */

#ifdef DEBUG
  packsize = nwritten*2;
  printf("Entered Bitpack64_access with oligo %u, packsize %d, remainder %d, row %d, index %d\n",
	 oligo,packsize,remainder,row,index);
  printf("bitpack:\n");
  for (i = 0; i < nwritten*4; i += 4) {
    printf("%08X %08X %08X %08X\n",bitpack[i],bitpack[i+1],bitpack[i+2],bitpack[i+3]);
  }
  printf("\n");
#endif
  
  return (accessor_table[index])(&(bitpack[row]));
}

#else

UINT4
Bitpack64_access (Oligospace_T oligo, UINT4 *ptrs, UINT4 *comp) {
  UINT4 *info, start;
  int nwritten, remainder;
  UINT4 *bitpack;
  int index, column;
#ifdef DEBUG
  int packsize, i;
#endif

  info = &(ptrs[oligo/BLOCKSIZE * DIRECT_METAINFO_SIZE]);

#ifdef WORDS_BIGENDIAN
  start = Bigendian_convert_uint(info[0]);
  bitpack = (UINT4 *) &(comp[start*4]);
  nwritten = Bigendian_convert_uint(info[1]) - start;	/* In 128-bit registers */
#else
  start = info[0];
  bitpack = (UINT4 *) &(comp[start*4]);
  nwritten = info[1] - start;	/* In 128-bit registers */
#endif

  remainder = oligo % BLOCKSIZE;
  index = nwritten*16 + remainder/4;
  column = remainder % 4;

#ifdef DEBUG
  packsize = nwritten*2;
  printf("Entered Bitpack64_access with oligo %u, packsize %d, remainder %d, column %d, index %d\n",
	 oligo,packsize,remainder,column,index);
  printf("bitpack:\n");
  for (i = 0; i < nwritten*4; i += 4) {
    printf("%08X %08X %08X %08X\n",bitpack[i],bitpack[i+1],bitpack[i+2],bitpack[i+3]);
  }
  printf("\n");
#endif

  return (accessor_table[index])(&(bitpack[column]));
}

UINT4
Bitpack64_access_bitpack (Oligospace_T oligo, int packsize, UINT4 *bitpack) {
  int nregisters, remainder;	/* nregisters is same as nwritten */
  int index, column;

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*16 + remainder/4;
  column = remainder % 4;

  return (accessor_table[index])(&(bitpack[column]));
}

static UINT4 threshold[17] = 
  {/*0*/0U, /*2*/3U, /*4*/15U, /*6*/63U,
   /*8*/255U, /*10*/1023U, /*12*/4095U, /*14*/16383U,
   /*16*/65535U, /*18*/262143U, /*20*/1048575U, /*22*/4194303U,
   /*24*/16777215U, /*26*/67108863U, /*28*/268435455U, /*30*/1073741823U,
   /*32*/4294967295U};


bool
Bitpack64_access_filledp (Oligospace_T oligo, int packsize, UINT4 *bitpack) {
  int nregisters, remainder;	/* nregisters is same as nwritten */
  int index, column;
  UINT4 value;

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*16 + remainder/4;
  column = remainder % 4;

  value = (accessor_table[index])(&(bitpack[column]));
  if (value == threshold[nregisters]) {
    return true;
  } else {
    return false;
  }
}

int
Bitpack64_access_new_packsize (Oligospace_T oligo, int old_packsize, UINT4 *bitpack, int increment) {
  int new_packsize;
  int nregisters, remainder;	/* nregisters is same as nwritten */
  int index, column;
  UINT4 value;
  Blocksize_T firstbit;
#ifdef HAVE_BUILTIN_CLZ
#elif defined(HAVE_ASM_BSR)
  int msb;
#endif


  nregisters = old_packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*16 + remainder/4;
  column = remainder % 4;

  value = (accessor_table[index])(&(bitpack[column]));
  value += increment;

  if (value == 0) {
    new_packsize = 0;
  } else {
#ifdef HAVE_BUILTIN_CLZ
    firstbit = __builtin_clz(value);
    new_packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
    asm("bsr %1,%0" : "=r"(msb) : "r"(value));
    new_packsize = msb + 1;
#else
    firstbit = ((value >> 16) ? clz_table[value >> 16] : 16 + clz_table[value]);
    new_packsize = 32 - firstbit;
#endif
  }

  new_packsize = (new_packsize + 1) & ~1; /* Converts packsizes to the next multiple of 2 */
  if (new_packsize > old_packsize) {
    return new_packsize;
  } else {
    return old_packsize;
  }
}



#endif



static void
extract_00 (UINT4 *out, const UINT4 *in) {
  /* 00 */
  out[0] = 0;

  /* 01 */
  out[4] = 0;

  /* 02 */
  out[8] = 0;

  /* 03 */
  out[12] = 0;

  /* 04 */
  out[16] = 0;

  /* 05 */
  out[20] = 0;

  /* 06 */
  out[24] = 0;

  /* 07 */
  out[28] = 0;

  /* 08 */
  out[32] = 0;

  /* 09 */
  out[36] = 0;

  /* 10 */
  out[40] = 0;

  /* 11 */
  out[44] = 0;

  /* 12 */
  out[48] = 0;

  /* 13 */
  out[52] = 0;

  /* 14 */
  out[56] = 0;

  /* 15 */
  out[60] = 0;

  return;
}

static void
extract_02 (UINT4 *out, const UINT4 *in) {
  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 2 ) ;

  /* 01 */
  out[4] = ( CONVERT(*in) >>  2  )   % (1U << 2 ) ;

  /* 02 */
  out[8] = ( CONVERT(*in) >>  4  )   % (1U << 2 ) ;

  /* 03 */
  out[12] = ( CONVERT(*in) >>  6  )   % (1U << 2 ) ;

  /* 04 */
  out[16] = ( CONVERT(*in) >>  8  )   % (1U << 2 ) ;

  /* 05 */
  out[20] = ( CONVERT(*in) >>  10  )   % (1U << 2 ) ;

  /* 06 */
  out[24] = ( CONVERT(*in) >>  12  )   % (1U << 2 ) ;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  14  )   % (1U << 2 ) ;

  /* 08 */
  out[32] = ( CONVERT(*in) >>  16  )   % (1U << 2 ) ;

  /* 09 */
  out[36] = ( CONVERT(*in) >>  18  )   % (1U << 2 ) ;

  /* 10 */
  out[40] = ( CONVERT(*in) >>  20  )   % (1U << 2 ) ;

  /* 11 */
  out[44] = ( CONVERT(*in) >>  22  )   % (1U << 2 ) ;

  /* 12 */
  out[48] = ( CONVERT(*in) >>  24  )   % (1U << 2 ) ;

  /* 13 */
  out[52] = ( CONVERT(*in) >>  26  )   % (1U << 2 ) ;

  /* 14 */
  out[56] = ( CONVERT(*in) >>  28  )   % (1U << 2 ) ;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  30  )   % (1U << 2 ) ;

  return;
}

static void
extract_04 (UINT4 *out, const UINT4 *in) {

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 4 ) ;

  /* 01 */
  out[4] = ( CONVERT(*in) >>  4  )   % (1U << 4 ) ;

  /* 02 */
  out[8] = ( CONVERT(*in) >>  8  )   % (1U << 4 ) ;

  /* 03 */
  out[12] = ( CONVERT(*in) >>  12  )   % (1U << 4 ) ;

  /* 04 */
  out[16] = ( CONVERT(*in) >>  16 )   % (1U << 4 ) ;

  /* 05 */
  out[20] = ( CONVERT(*in) >>  20  )   % (1U << 4 ) ;

  /* 06 */
  out[24] = ( CONVERT(*in) >>  24  )   % (1U << 4 ) ;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  28  )   % (1U << 4 ) ;

  /* 08 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >>  0  )   % (1U << 4 ) ;

  /* 09 */
  out[36] = ( CONVERT(*in) >>  4  )   % (1U << 4 ) ;

  /* 10 */
  out[40] = ( CONVERT(*in) >>  8  )   % (1U << 4 ) ;

  /* 11 */
  out[44] = ( CONVERT(*in) >>  12  )   % (1U << 4 ) ;

  /* 12 */
  out[48] = ( CONVERT(*in) >>  16  )   % (1U << 4 ) ;

  /* 13 */
  out[52] = ( CONVERT(*in) >>  20  )   % (1U << 4 ) ;

  /* 14 */
  out[56] = ( CONVERT(*in) >>  24  )   % (1U << 4 ) ;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  28  )   % (1U << 4 ) ;

  return;
}


static void
extract_06 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 6 ) ;

  /* 01 */
  out[4] = ( CONVERT(*in) >>  6  )   % (1U << 6 ) ;

  /* 02 */
  out[8] = ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;

  /* 03 */
  out[12] = ( CONVERT(*in) >>  18  )   % (1U << 6 ) ;

  /* 04 */
  out[16] = ( CONVERT(*in) >>  24  )   % (1U << 6 ) ;

  /* 05 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 6 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 );
  out[20] = value;

  /* 06 */
  out[24] = ( CONVERT(*in) >>  4  )   % (1U << 6 ) ;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  10  )   % (1U << 6 ) ;

  /* 08 */
  out[32] = ( CONVERT(*in) >>  16  )   % (1U << 6 ) ;

  /* 09 */
  out[36] = ( CONVERT(*in) >>  22  )   % (1U << 6 ) ;

  /* 10 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 6 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 );
  out[40] = value;

  /* 11 */
  out[44] = ( CONVERT(*in) >>  2  )   % (1U << 6 ) ;

  /* 12 */
  out[48] = ( CONVERT(*in) >>  8  )   % (1U << 6 ) ;

  /* 13 */
  out[52] = ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;

  /* 14 */
  out[56] = ( CONVERT(*in) >>  20  )   % (1U << 6 ) ;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  26  )   % (1U << 6 ) ;

  return;
}


static void
extract_08 (UINT4 *out, const UINT4 *in) {

  /* 00 */
  out[0] = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;

  /* 01 */
  out[4] = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;

  /* 02 */
  out[8] = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;

  /* 03 */
  out[12] = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;

  /* 04 */
  in += WORD_INCR;
  out[16] = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;

  /* 05 */
  out[20] = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;

  /* 06 */
  out[24] = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;

  /* 07 */
  out[28] = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;

  /* 08 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;

  /* 09 */
  out[36] = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;

  /* 10 */
  out[40] = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;

  /* 11 */
  out[44] = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;

  /* 12 */
  in += WORD_INCR;
  out[48] = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;

  /* 13 */
  out[52] = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;

  /* 14 */
  out[56] = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;

  /* 15 */
  out[60] = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;

  return;
}

static void
extract_10 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 10 ) ;

  /* 01 */
  out[4] = ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;

  /* 02 */
  out[8] = ( CONVERT(*in) >>  20  )   % (1U << 10 ) ;

  /* 03 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 );
  out[12] = value;

  /* 04 */
  out[16] = ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;

  /* 05 */
  out[20] = ( CONVERT(*in) >>  18  )   % (1U << 10 ) ;

  /* 06 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 );
  out[24] = value;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  6  )   % (1U << 10 ) ;

  /* 08 */
  out[32] = ( CONVERT(*in) >>  16  )   % (1U << 10 ) ;

  /* 09 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 );
  out[36] = value;

  /* 10 */
  out[40] = ( CONVERT(*in) >>  4  )   % (1U << 10 ) ;

  /* 11 */
  out[44] = ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;

  /* 12 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 );
  out[48] = value;

  /* 13 */
  out[52] = ( CONVERT(*in) >>  2  )   % (1U << 10 ) ;

  /* 14 */
  out[56] = ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  22  )   % (1U << 10 ) ;

  return;
}

static void
extract_12 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;

  /* 01 */
  out[4] = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;

  /* 02 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  out[8] = value;

  /* 03 */
  out[12] = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;

  /* 04 */
  out[16] = ( CONVERT(*in) >>  16  )   % (1U << 12 ) ;

  /* 05 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  out[20] = value;

  /* 06 */
  out[24] = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  20  )   % (1U << 12 ) ;

  /* 08 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;

  /* 09 */
  out[36] = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;

  /* 10 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  out[40] = value;

  /* 11 */
  out[44] = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;

  /* 12 */
  out[48] = ( CONVERT(*in) >>  16  )   % (1U << 12 ) ;

  /* 13 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  out[52] = value;

  /* 14 */
  out[56] = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  20  )   % (1U << 12 ) ;

  return;
}


static void
extract_14 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 14 ) ;

  /* 01 */
  out[4] = ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;

  /* 02 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 );
  out[8] = value;

  /* 03 */
  out[12] = ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;

  /* 04 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 );
  out[16] = value;

  /* 05 */
  out[20] = ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;

  /* 06 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 );
  out[24] = value;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  2  )   % (1U << 14 ) ;

  /* 08 */
  out[32] = ( CONVERT(*in) >>  16  )   % (1U << 14 ) ;

  /* 09 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 );
  out[36] = value;

  /* 10 */
  out[40] = ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;

  /* 11 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 );
  out[44] = value;

  /* 12 */
  out[48] = ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;

  /* 13 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 );
  out[52] = value;

  /* 14 */
  out[56] = ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  18  )   % (1U << 14 ) ;

  return;
}

static void
extract_16 (UINT4 *out, const UINT4 *in) {

  /* 00 */
  out[0] = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;

  /* 01 */
  out[4] = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;

  /* 02 */
  in += WORD_INCR;
  out[8] = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;

  /* 03 */
  out[12] = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;

  /* 04 */
  in += WORD_INCR;
  out[16] = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;

  /* 05 */
  out[20] = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;

  /* 06 */
  in += WORD_INCR;
  out[24] = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;

  /* 07 */
  out[28] = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;

  /* 08 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;

  /* 09 */
  out[36] = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;

  /* 10 */
  in += WORD_INCR;
  out[40] = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;

  /* 11 */
  out[44] = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;

  /* 12 */
  in += WORD_INCR;
  out[48] = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;

  /* 13 */
  out[52] = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;

  /* 14 */
  in += WORD_INCR;
  out[56] = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;

  /* 15 */
  out[60] = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;

  return;
}


static void
extract_18 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 18 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 18 - 4 );
  out[4] = value;

  /* 02 */
  out[8] = ( CONVERT(*in) >>  4  )   % (1U << 18 ) ;

  /* 03 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 18 - 8 );
  out[12] = value;

  /* 04 */
  out[16] = ( CONVERT(*in) >>  8  )   % (1U << 18 ) ;

  /* 05 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 18 - 12 );
  out[20] = value;

  /* 06 */
  out[24] = ( CONVERT(*in) >>  12  )   % (1U << 18 ) ;

  /* 07 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 18 - 16 );
  out[28] = value;

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 18 - 2 );
  out[32] = value;

  /* 09 */
  out[36] = ( CONVERT(*in) >>  2  )   % (1U << 18 ) ;

  /* 10 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 18 - 6 );
  out[40] = value;

  /* 11 */
  out[44] = ( CONVERT(*in) >>  6  )   % (1U << 18 ) ;

  /* 12 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 18 - 10 );
  out[48] = value;

  /* 13 */
  out[52] = ( CONVERT(*in) >>  10  )   % (1U << 18 ) ;

  /* 14 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 14 ))<<( 18 - 14 );
  out[56] = value;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  14  )   % (1U << 18 ) ;

  return;
}

static void
extract_20 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 20 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 );
  out[4] = value;

  /* 02 */
  out[8] = ( CONVERT(*in) >>  8  )   % (1U << 20 ) ;

  /* 03 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 );
  out[12] = value;

  /* 04 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 );
  out[16] = value;

  /* 05 */
  out[20] = ( CONVERT(*in) >>  4  )   % (1U << 20 ) ;

  /* 06 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 );
  out[24] = value;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  12  )   % (1U << 20 ) ;

  /* 08 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >>  0  )   % (1U << 20 ) ;

  /* 09 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 );
  out[36] = value;

  /* 10 */
  out[40] = ( CONVERT(*in) >>  8  )   % (1U << 20 ) ;

  /* 11 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 );
  out[44] = value;

  /* 12 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 );
  out[48] = value;

  /* 13 */
  out[52] = ( CONVERT(*in) >>  4  )   % (1U << 20 ) ;

  /* 14 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 );
  out[56] = value;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  12  )   % (1U << 20 ) ;

  return;
}

static void
extract_22 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 22 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 22 - 12 );
  out[4] = value;

  /* 02 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 22 - 2 );
  out[8] = value;

  /* 03 */
  out[12] = ( CONVERT(*in) >>  2  )   % (1U << 22 ) ;

  /* 04 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 14 ))<<( 22 - 14 );
  out[16] = value;

  /* 05 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 22 - 4 );
  out[20] = value;

  /* 06 */
  out[24] = value = ( CONVERT(*in) >>  4  )   % (1U << 22 ) ;

  /* 07 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 22 - 16 );
  out[28] = value;

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 22 - 6 );
  out[32] = value;

  /* 09 */
  out[36] = ( CONVERT(*in) >>  6  )   % (1U << 22 ) ;

  /* 10 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 18 ))<<( 22 - 18 );
  out[40] = value;

  /* 11 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 22 - 8 );
  out[44] = value;

  /* 12 */
  out[48] = ( CONVERT(*in) >>  8  )   % (1U << 22 ) ;

  /* 13 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 22 - 20 );
  out[52] = value;

  /* 14 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 22 - 10 );
  out[56] = value;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  10  )   % (1U << 22 ) ;

  return;
}

static void
extract_24 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  out[4] = value;

  /* 02 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  out[8] = value;

  /* 03 */
  out[12] = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;

  /* 04 */
  in += WORD_INCR;
  out[16] = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;

  /* 05 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  out[20] = value;

  /* 06 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  out[24] = value;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;

  /* 08 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;

  /* 09 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  out[36] = value;

  /* 10 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  out[40] = value;

  /* 11 */
  out[44] = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;

  /* 12 */
  in += WORD_INCR;
  out[48] = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;

  /* 13 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  out[52] = value;

  /* 14 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  out[56] = value;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;

  return;
}

static void
extract_26 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 26 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 26 - 20 );
  out[4] = value;

  /* 02 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 14 ))<<( 26 - 14 );
  out[8] = value;

  /* 03 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 26 - 8 );
  out[12] = value;

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 26 - 2 );
  out[16] = value;

  /* 05 */
  out[20] = ( CONVERT(*in) >>  2  )   % (1U << 26 ) ;

  /* 06 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 22 ))<<( 26 - 22 );
  out[24] = value;

  /* 07 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 26 - 16 );
  out[28] = value;

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 26 - 10 );
  out[32] = value;

  /* 09 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 26 - 4 );
  out[36] = value;

  /* 10 */
  out[40] = ( CONVERT(*in) >>  4  )   % (1U << 26 ) ;

  /* 11 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 24 ))<<( 26 - 24 );
  out[44] = value;

  /* 12 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 18 ))<<( 26 - 18 );
  out[48] = value;

  /* 13 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 26 - 12 );
  out[52] = value;

  /* 14 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 26 - 6 );
  out[56] = value;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  6  )   % (1U << 26 ) ;

  return;
}

static void
extract_28 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 28 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 );
  out[4] = value;

  /* 02 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 );
  out[8] = value;

  /* 03 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 );
  out[12] = value;

  /* 04 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 );
  out[16] = value;

  /* 05 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 );
  out[20] = value;

  /* 06 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 );
  out[24] = value;

  /* 07 */
  out[28] = ( CONVERT(*in) >>  4  )   % (1U << 28 ) ;

  /* 08 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >>  0  )   % (1U << 28 ) ;

  /* 09 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 );
  out[36] = value;

  /* 10 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 );
  out[40] = value;

  /* 11 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 );
  out[44] = value;

  /* 12 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 );
  out[48] = value;

  /* 13 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 );
  out[52] = value;

  /* 14 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 );
  out[56] = value;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  4  )   % (1U << 28 ) ;

  return;
}


static void
extract_30 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 30 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 28 ))<<( 30 - 28 );
  out[4] = value;

  /* 02 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 26 ))<<( 30 - 26 );
  out[8] = value;

  /* 03 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 24 ))<<( 30 - 24 );
  out[12] = value;

  /* 04 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 22 ))<<( 30 - 22 );
  out[16] = value;

  /* 05 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 30 - 20 );
  out[20] = value;

  /* 06 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 18 ))<<( 30 - 18 );
  out[24] = value;

  /* 07 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 30 - 16 );
  out[28] = value;

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 14 ))<<( 30 - 14 );
  out[32] = value;

  /* 09 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 30 - 12 );
  out[36] = value;

  /* 10 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 30 - 10 );
  out[40] = value;

  /* 11 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 30 - 8 );
  out[44] = value;

  /* 12 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 30 - 6 );
  out[48] = value;

  /* 13 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 30 - 4 );
  out[52] = value;

  /* 14 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 30 - 2 );
  out += WORD_INCR;
  out[56] = value;

  /* 15 */
  out[60] = ( CONVERT(*in) >>  2  )   % (1U << 30 ) ;

  return;
}


static void
extract_32 (UINT4 *out, const UINT4 *in) {

  /* 00 */
  out[0] = CONVERT(*in);

  /* 01 */
  in += WORD_INCR;
  out[4] = CONVERT(*in);

  /* 02 */
  in += WORD_INCR;
  out[8] = CONVERT(*in);

  /* 03 */
  in += WORD_INCR;
  out[12] = CONVERT(*in);

  /* 04 */
  in += WORD_INCR;
  out[16] = CONVERT(*in);

  /* 05 */
  in += WORD_INCR;
  out[20] = CONVERT(*in);

  /* 06 */
  in += WORD_INCR;
  out[24] = CONVERT(*in);

  /* 07 */
  in += WORD_INCR;
  out[28] = CONVERT(*in);

  /* 08 */
  in += WORD_INCR;
  out[32] = CONVERT(*in);

  /* 09 */
  in += WORD_INCR;
  out[36] = CONVERT(*in);

  /* 10 */
  in += WORD_INCR;
  out[40] = CONVERT(*in);

  /* 11 */
  in += WORD_INCR;
  out[44] = CONVERT(*in);

  /* 12 */
  in += WORD_INCR;
  out[48] = CONVERT(*in);

  /* 13 */
  in += WORD_INCR;
  out[52] = CONVERT(*in);

  /* 14 */
  in += WORD_INCR;
  out[56] = CONVERT(*in);

  /* 15 */
  in += WORD_INCR;
  out[60] = CONVERT(*in);

  return;
}


typedef void (*Extractor_T) (UINT4 *, const UINT4 *);

static Extractor_T extractor_table[17] =
  {extract_00, extract_02, extract_04, extract_06,
   extract_08, extract_10, extract_12, extract_14,
   extract_16, extract_18, extract_20, extract_22,
   extract_24, extract_26, extract_28, extract_30,
   extract_32};

void
Bitpack64_extract_bitpack (UINT4 *out, int packsize, UINT4 *bitpack) {
  int nregisters;

  nregisters = packsize / 2;
  (extractor_table[nregisters])(&(out[0]),&(bitpack[/*column*/0]));
  (extractor_table[nregisters])(&(out[1]),&(bitpack[/*column*/1]));
  (extractor_table[nregisters])(&(out[2]),&(bitpack[/*column*/2]));
  (extractor_table[nregisters])(&(out[3]),&(bitpack[/*column*/3]));

  return;
}
