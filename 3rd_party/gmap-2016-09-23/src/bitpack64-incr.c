static char rcsid[] = "$Id: bitpack64-incr.c 197549 2016-09-08 01:14:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bitpack64-incr.h"
#include "bitpack64-access.h"

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#define CONVERT(x) Bigendian_convert_uint(x)
#else
#define CONVERT(x) x
#endif


#define CHECK 1


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define BLOCKSIZE 64

/* Vertical access is slightly more efficient than horizontal */
/* Choice needs to match that in bitpack64-access */

#ifdef HORIZONTAL
#define WORD_INCR 1		/* 1 for horizontal; 4 for vertical */
#else
#define WORD_INCR 4
#endif


#define BIT0  0x1;
#define BIT2  0x4;
#define BIT4  0x10;
#define BIT6  0x40;
#define BIT8  0x100;
#define BIT10 0x400;
#define BIT12 0x1000;
#define BIT14 0x4000;
#define BIT16 0x10000;
#define BIT18 0x40000;
#define BIT20 0x100000;
#define BIT22 0x400000;
#define BIT24 0x1000000;
#define BIT26 0x4000000;
#define BIT28 0x10000000;
#define BIT30 0x40000000;



static UINT4
incr_00 (UINT4 *in) {
  abort();
  return 0U;
}



static UINT4
incr_02_00 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  0  )   % (1U << 2 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_02_01 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  2  )   % (1U << 2 ) ;
  *in += BIT2;
  return prev;
}

static UINT4
incr_02_02 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  4  )   % (1U << 2 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_02_03 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  6  )   % (1U << 2 ) ;
  *in += BIT6;
  return prev;
}

static UINT4
incr_02_04 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  8  )   % (1U << 2 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_02_05 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  10  )   % (1U << 2 ) ;
  *in += BIT10;
  return prev;
}

static UINT4
incr_02_06 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  12  )   % (1U << 2 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_02_07 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  14  )   % (1U << 2 ) ;
  *in += BIT14;
  return prev;
}

static UINT4
incr_02_08 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  16  )   % (1U << 2 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_02_09 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  18  )   % (1U << 2 ) ;
  *in += BIT18;
  return prev;
}

static UINT4
incr_02_10 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  20  )   % (1U << 2 ) ;
  *in += BIT20;
  return prev;
}

static UINT4
incr_02_11 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  22  )   % (1U << 2 ) ;
  *in += BIT22;
  return prev;
}

static UINT4
incr_02_12 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  24  )   % (1U << 2 ) ;
  *in += BIT24;
  return prev;
}

static UINT4
incr_02_13 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  26  )   % (1U << 2 ) ;
  *in += BIT26;
  return prev;
}

static UINT4
incr_02_14 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  28  )   % (1U << 2 ) ;
  *in += BIT28;
  return prev;
}

static UINT4
incr_02_15 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >>  30  )   % (1U << 2 ) ;
  *in += BIT30;
  return prev;
}



static UINT4
incr_04_00 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_04_01 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_04_02 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_04_03 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_04_04 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >> 16 )   % (1U << 4 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_04_05 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >> 20 )   % (1U << 4 ) ;
  *in += BIT20;
  return prev;
}

static UINT4
incr_04_06 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >> 24 )   % (1U << 4 ) ;
  *in += BIT24;
  return prev;
}

static UINT4
incr_04_07 (UINT4 *in) {
  UINT4 prev = ( CONVERT(*in) >> 28 )   % (1U << 4 ) ;
  *in += BIT28;
  return prev;
}

static UINT4
incr_04_08 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_04_09 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_04_10 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_04_11 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_04_12 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 4 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_04_13 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 20 )   % (1U << 4 ) ;
  *in += BIT20;
  return prev;
}

static UINT4
incr_04_14 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 24 )   % (1U << 4 ) ;
  *in += BIT24;
  return prev;
}

static UINT4
incr_04_15 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 28 )   % (1U << 4 ) ;
  *in += BIT28;
  return prev;
}


static UINT4
incr_06_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 6 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_06_01 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  6  )   % (1U << 6 ) ;
  *in += BIT6;
  return prev;
}

static UINT4
incr_06_02 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_06_03 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  18  )   % (1U << 6 ) ;
  *in += BIT18;
  return prev;
}

static UINT4
incr_06_04 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  24  )   % (1U << 6 ) ;
  *in += BIT24;
  return prev;
}

static UINT4
incr_06_05 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 6 ) ;
  carry = ( (prev + 1) >> (6 - 4) ); 
  *in += BIT30;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_06_06 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 6 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_06_07 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 6 ) ;
  *in += BIT10;
  return prev;
}

static UINT4
incr_06_08 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 6 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_06_09 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 6 ) ;
  *in += BIT22;
  return prev;
}

static UINT4
incr_06_10 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 6 ) ;
  carry = ( (prev + 1) >> (6 - 2) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 );
  *in += carry;

  return prev;
}

static UINT4
incr_06_11 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 6 ) ;
  *in += BIT2;
  return prev;
}

static UINT4
incr_06_12 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 6 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_06_13 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;
  *in += BIT14;
  return prev;
}

static UINT4
incr_06_14 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 6 ) ;
  *in += BIT20;
  return prev;
}

static UINT4
incr_06_15 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 6 ) ;
  *in += BIT26;
  return prev;
}


static UINT4
incr_08_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_08_01 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_08_02 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_08_03 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
  *in += BIT24;
  return prev;
}

static UINT4
incr_08_04 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_08_05 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_08_06 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_08_07 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
  *in += BIT24;
  return prev;
}

static UINT4
incr_08_08 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_08_09 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_08_10 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_08_11 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
  *in += BIT24;
  return prev;
}

static UINT4
incr_08_12 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_08_13 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_08_14 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_08_15 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
  *in += BIT24;
  return prev;
}



static UINT4
incr_10_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 10 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_10_01 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;
  *in += BIT10;
  return prev;
}

static UINT4
incr_10_02 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  20  )   % (1U << 10 ) ;
  *in += BIT20;
  return prev;
}

static UINT4
incr_10_03 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 10 ) ;
  carry = ( (prev + 1) >> (10 - 8) );
  *in += BIT30;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_10_04 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_10_05 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 10 ) ;
  *in += BIT18;
  return prev;
}

static UINT4
incr_10_06 (UINT4 *in) {
  UINT4 prev, carry;
  
  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 10 ) ;
  carry = ( (prev + 1) >> (10 - 6) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 );
  *in += carry;

  return prev;
}

static UINT4
incr_10_07 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 10 ) ;
  *in += BIT6;
  return prev;
}

static UINT4
incr_10_08 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 10 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_10_09 (UINT4 *in) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 10 ) ;
  carry = ( (prev + 1) >> (10 - 4) );
  *in += BIT26;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_10_10 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 10 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_10_11 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;
  *in += BIT14;
  return prev;
}

static UINT4
incr_10_12 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 10 ) ;
  carry = ( (prev + 1) >> (10 - 2) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 );
  *in += carry;

  return prev;
}

static UINT4
incr_10_13 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 10 ) ;
  *in += BIT2;
  return prev;
}

static UINT4
incr_10_14 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_10_15 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 10 ) ;
  *in += BIT22;
  return prev;
}



static UINT4
incr_12_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_12_01 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_12_02 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  carry = ( (prev + 1) >> (12 - 4) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_12_03 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_12_04 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 12 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_12_05 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  carry = ( (prev + 1) >> (12 - 8));
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_12_06 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_12_07 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 12 ) ;
  *in += BIT20;
  return prev;
}

static UINT4
incr_12_08 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_12_09 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_12_10 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  carry = ( (prev + 1) >> (12 - 4) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_12_11 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_12_12 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 12 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_12_13 (UINT4 *in) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  carry = ( (prev + 1) >> (12 - 8));
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_12_14 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_12_15 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 12 ) ;
  *in += BIT20;
  return prev;
}


static UINT4
incr_14_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 14 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_14_01 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;
  *in += BIT14;
  return prev;
}

static UINT4
incr_14_02 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  28  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 10) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 );
  *in += carry;

  return prev;
}

static UINT4
incr_14_03 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;
  *in += BIT10;
  return prev;
}

static UINT4
incr_14_04 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 6) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 );
  *in += carry;

  return prev;
}

static UINT4
incr_14_05 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;
  *in += BIT6;
  return prev;
}

static UINT4
incr_14_06 (UINT4 *in) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 2) );
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 );
  *in += carry;

  return prev;
}

static UINT4
incr_14_07 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 14 ) ;
  *in += BIT2;
  return prev;
}

static UINT4
incr_14_08 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 14 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_14_09 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 12) );
  *in += BIT30;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_14_10 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_14_11 (UINT4 *in) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 8) );
  *in += BIT26;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_14_12 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_14_13 (UINT4 *in) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 4) );
  *in += BIT22;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_14_14 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_14_15 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 14 ) ;
  *in += BIT18;
  return prev;
}


static UINT4
incr_16_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_16_01 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_16_02 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_16_03 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_16_04 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_16_05 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_16_06 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_16_07 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_16_08 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_16_09 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_16_10 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_16_11 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_16_12 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_16_13 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *in += BIT16;
  return prev;
}

static UINT4
incr_16_14 (UINT4 *in) {
  UINT4 prev;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_16_15 (UINT4 *in) {
  UINT4 prev;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *in += BIT16;
  return prev;
}


static UINT4
incr_18_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 18 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_18_01 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  18  )   % (1U << 18 ) ;
  carry = ( (prev + 1) >> (18 - 4) );
  *in += BIT18;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 18 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_18_02 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 18 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_18_03 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 18 ) ;
  carry = ( (prev + 1) >> (18 - 8) );
  *in += BIT22;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 18 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_18_04 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 18 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_18_05 (UINT4 *in) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 18 ) ;
  carry = ( (prev + 1) >> (18 - 12) );
  *in += BIT26;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 18 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_18_06 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 18 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_18_07 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 18 ) ;
  carry = ( (prev + 1) >> (18 - 16) );
  *in += BIT30;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 18 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_18_08 (UINT4 *in) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 18 ) ;
  carry = ( (prev + 1) >> (18 - 2) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 18 - 2 );
  *in += carry;

  return prev;
}

static UINT4
incr_18_09 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 18 ) ;
  *in += BIT2;
  return prev;
}

static UINT4
incr_18_10 (UINT4 *in) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 18 ) ;
  carry = ( (prev + 1) >> (18 - 6) );
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 6 ))<<( 18 - 6 );
  *in += carry;

  return prev;
}

static UINT4
incr_18_11 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 18 ) ;
  *in += BIT6;
  return prev;
}

static UINT4
incr_18_12 (UINT4 *in) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 18 ) ;
  carry = ( (prev + 1) >> (18 - 10) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 10 ))<<( 18 - 10 );
  *in += carry;

  return prev;
}

static UINT4
incr_18_13 (UINT4 *in) {
  UINT4 prev;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 18 ) ;
  *in += BIT10;
  return prev;
}

static UINT4
incr_18_14 (UINT4 *in) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 18 ) ;
  carry = ( (prev + 1) >> (18 - 14) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 14 ))<<( 18 - 14 );
  *in += carry;

  return prev;
}

static UINT4
incr_18_15 (UINT4 *in) {
  UINT4 prev;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 18 ) ;
  *in += BIT14;
  return prev;
}


static UINT4
incr_20_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 20 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_20_01 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  carry = ( (prev + 1) >> (20 - 8) );
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_20_02 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 20 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_20_03 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  carry = ( (prev + 1) >> (20 - 16) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_20_04 (UINT4 *in) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  carry = ( (prev + 1) >> (20 - 4) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_20_05 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 20 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_20_06 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  carry = ( (prev + 1) >> (20 - 12) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_20_07 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 20 ) ;
  *in += BIT12;
  return prev;
}

static UINT4
incr_20_08 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  0  )   % (1U << 20 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_20_09 (UINT4 *in) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  carry = ( (prev + 1) >> (20 - 8) );
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_20_10 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 20 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_20_11 (UINT4 *in) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  carry = ( (prev + 1) >> (20 - 16) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_20_12 (UINT4 *in) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  carry = ( (prev + 1) >> (20 - 4) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_20_13 (UINT4 *in) {
  UINT4 prev;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 20 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_20_14 (UINT4 *in) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  carry = ( (prev + 1) >> (20 - 12) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_20_15 (UINT4 *in) {
  UINT4 prev;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 20 ) ;
  *in += BIT12;
  return prev;
}


static UINT4
incr_22_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 22 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_22_01 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  22  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 12) );
  *in += BIT22;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 22 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_02 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 2) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 22 - 2 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_03 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 22 ) ;
  *in += BIT2;
  return prev;
}

static UINT4
incr_22_04 (UINT4 *in) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 14) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 14 ))<<( 22 - 14 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_05 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 4) );
  *in += BIT14;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 22 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_06 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 22 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_22_07 (UINT4 *in) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 16) );
  *in += BIT26;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 22 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_08 (UINT4 *in) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 6) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 6 ))<<( 22 - 6 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_09 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 22 ) ;
  *in += BIT6;
  return prev;
}

static UINT4
incr_22_10 (UINT4 *in) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 18) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 18 ))<<( 22 - 18 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_11 (UINT4 *in) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 8) );
  *in += BIT18;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 22 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_12 (UINT4 *in) {
  UINT4 prev;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 22 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_22_13 (UINT4 *in) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 20) );
  *in += BIT30;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 20 ))<<( 22 - 20 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_14 (UINT4 *in) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 22 ) ;
  carry = ( (prev + 1) >> (22 - 10) );
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 10 ))<<( 22 - 10 );
  *in += carry;

  return prev;
}

static UINT4
incr_22_15 (UINT4 *in) {
  UINT4 prev;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 22 ) ;
  *in += BIT10;
  return prev;
}



static UINT4
incr_24_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_24_01 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev + 1) >> (24 - 16) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_24_02 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev + 1) >> (24 - 8) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_24_03 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_24_04 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_24_05 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev + 1) >> (24 - 16) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_24_06 (UINT4 *in) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev + 1) >> (24 - 8) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_24_07 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_24_08 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_24_09 (UINT4 *in) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev + 1) >> (24 - 16) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_24_10 (UINT4 *in) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev + 1) >> (24 - 8) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_24_11 (UINT4 *in) {
  UINT4 prev;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
  *in += BIT8;
  return prev;
}

static UINT4
incr_24_12 (UINT4 *in) {
  UINT4 prev;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_24_13 (UINT4 *in) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev + 1) >> (24 - 16) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_24_14 (UINT4 *in) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev + 1) >> (24 - 8) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_24_15 (UINT4 *in) {
  UINT4 prev;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
  *in += BIT8;
  return prev;
}



static UINT4
incr_26_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 26 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_26_01 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  26  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 20) );
  *in += BIT26;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 20 ))<<( 26 - 20 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_02 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 14) );
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 14 ))<<( 26 - 14 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_03 (UINT4 *in) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 8) );
  *in += BIT14;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 26 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_04 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 2) );
  *in += BIT8;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 26 - 2 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_05 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 26 ) ;
  *in += BIT2;
  return prev;
}

static UINT4
incr_26_06 (UINT4 *in) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 22) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 22 ))<<( 26 - 22 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_07 (UINT4 *in) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 16) );
  *in += BIT22;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 26 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_08 (UINT4 *in) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 10) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 10 ))<<( 26 - 10 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_09 (UINT4 *in) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 4) );
  *in += BIT10;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 26 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_10 (UINT4 *in) {
  UINT4 prev;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 26 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_26_11 (UINT4 *in) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 24) );
  *in += BIT30;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 24 ))<<( 26 - 24 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_12 (UINT4 *in) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 18) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 18 ))<<( 26 - 18 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_13 (UINT4 *in) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 12) );
  *in += BIT18;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 26 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_14 (UINT4 *in) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 26 ) ;
  carry = ( (prev + 1) >> (26 - 6) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 6 ))<<( 26 - 6 );
  *in += carry;

  return prev;
}

static UINT4
incr_26_15 (UINT4 *in) {
  UINT4 prev;

  in += 12 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 26 ) ;
  *in += BIT6;
  return prev;
}


static UINT4
incr_28_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 28 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_28_01 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 24) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_02 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 20) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_03 (UINT4 *in) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 16) );
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_04 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 12) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_05 (UINT4 *in) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 8) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_06 (UINT4 *in) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 4) );
  *in += BIT8;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_07 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 28 ) ;
  *in += BIT4;
  return prev;
}

static UINT4
incr_28_08 (UINT4 *in) {
  UINT4 prev;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  0  )   % (1U << 28 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_28_09 (UINT4 *in) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 24) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_10 (UINT4 *in) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 20) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_11 (UINT4 *in) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 16) );
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_12 (UINT4 *in) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 12) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_13 (UINT4 *in) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 8) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_14 (UINT4 *in) {
  UINT4 prev, carry;

  in += 12 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  carry = ( (prev + 1) >> (28 - 4) );
  *in += BIT8;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_28_15 (UINT4 *in) {
  UINT4 prev;

  in += 13 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 28 ) ;
  *in += BIT4;
  return prev;
}


static UINT4
incr_30_00 (UINT4 *in) {
  UINT4 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 30 ) ;
  *in += BIT0;
  return prev;
}

static UINT4
incr_30_01 (UINT4 *in) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 28) );
  *in += BIT30;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 28 ))<<( 30 - 28 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_02 (UINT4 *in) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 26) );
  *in += BIT28;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 26 ))<<( 30 - 26 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_03 (UINT4 *in) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 24) );
  *in += BIT26;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 24 ))<<( 30 - 24 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_04 (UINT4 *in) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 22) );
  *in += BIT24;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 22 ))<<( 30 - 22 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_05 (UINT4 *in) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 20) );
  *in += BIT22;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 20 ))<<( 30 - 20 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_06 (UINT4 *in) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 18));
  *in += BIT20;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 18 ))<<( 30 - 18 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_07 (UINT4 *in) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 16) );
  *in += BIT18;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 16 ))<<( 30 - 16 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_08 (UINT4 *in) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 14) );
  *in += BIT16;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 14 ))<<( 30 - 14 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_09 (UINT4 *in) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 12) );
  *in += BIT14;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 30 - 12 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_10 (UINT4 *in) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 10) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 10 ))<<( 30 - 10 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_11 (UINT4 *in) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 8) );
  *in += BIT10;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 30 - 8 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_12 (UINT4 *in) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 6) );
  *in += BIT8;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 6 ))<<( 30 - 6 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_13 (UINT4 *in) {
  UINT4 prev, carry;

  in += 12 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 4) );
  *in += BIT6;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 30 - 4 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_14 (UINT4 *in) {
  UINT4 prev, carry;

  in += 13 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 30 ) ;
  carry = ( (prev + 1) >> (30 - 2) );
  *in += BIT4;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 30 - 2 );
  *in += carry;

  return prev;
}

static UINT4
incr_30_15 (UINT4 *in) {
  UINT4 prev;

  in += 14 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 30 ) ;
  *in += BIT2;
  return prev;
}


static UINT4
incr_32_00 (UINT4 *in) {
  UINT4 prev;

  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_01 (UINT4 *in) {
  UINT4 prev;

  in += 1 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_02 (UINT4 *in) {
  UINT4 prev;

  in += 2 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_03 (UINT4 *in) {
  UINT4 prev;

  in += 3 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_04 (UINT4 *in) {
  UINT4 prev;

  in += 4 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_05 (UINT4 *in) {
  UINT4 prev;

  in += 5 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_06 (UINT4 *in) {
  UINT4 prev;

  in += 6 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}


static UINT4
incr_32_07 (UINT4 *in) {
  UINT4 prev;

  in += 7 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_08 (UINT4 *in) {
  UINT4 prev;

  in += 8 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_09 (UINT4 *in) {
  UINT4 prev;

  in += 9 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_10 (UINT4 *in) {
  UINT4 prev;

  in += 10 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_11 (UINT4 *in) {
  UINT4 prev;

  in += 11 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_12 (UINT4 *in) {
  UINT4 prev;

  in += 12 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_13 (UINT4 *in) {
  UINT4 prev;

  in += 13 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_14 (UINT4 *in) {
  UINT4 prev;

  in += 14 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT4
incr_32_15 (UINT4 *in) {
  UINT4 prev;

  in += 15 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}




typedef UINT4 (*Incrementor_T) (UINT4 *);

static Incrementor_T incrementor_table[272] =
  {incr_00, incr_00, incr_00, incr_00,
   incr_00, incr_00, incr_00, incr_00,
   incr_00, incr_00, incr_00, incr_00,
   incr_00, incr_00, incr_00, incr_00,

   incr_02_00, incr_02_01, incr_02_02, incr_02_03,
   incr_02_04, incr_02_05, incr_02_06, incr_02_07,
   incr_02_08, incr_02_09, incr_02_10, incr_02_11,
   incr_02_12, incr_02_13, incr_02_14, incr_02_15,

   incr_04_00, incr_04_01, incr_04_02, incr_04_03,
   incr_04_04, incr_04_05, incr_04_06, incr_04_07,
   incr_04_08, incr_04_09, incr_04_10, incr_04_11,
   incr_04_12, incr_04_13, incr_04_14, incr_04_15,

   incr_06_00, incr_06_01, incr_06_02, incr_06_03,
   incr_06_04, incr_06_05, incr_06_06, incr_06_07,
   incr_06_08, incr_06_09, incr_06_10, incr_06_11,
   incr_06_12, incr_06_13, incr_06_14, incr_06_15,

   incr_08_00, incr_08_01, incr_08_02, incr_08_03,
   incr_08_04, incr_08_05, incr_08_06, incr_08_07,
   incr_08_08, incr_08_09, incr_08_10, incr_08_11,
   incr_08_12, incr_08_13, incr_08_14, incr_08_15,

   incr_10_00, incr_10_01, incr_10_02, incr_10_03,
   incr_10_04, incr_10_05, incr_10_06, incr_10_07,
   incr_10_08, incr_10_09, incr_10_10, incr_10_11,
   incr_10_12, incr_10_13, incr_10_14, incr_10_15,

   incr_12_00, incr_12_01, incr_12_02, incr_12_03,
   incr_12_04, incr_12_05, incr_12_06, incr_12_07,
   incr_12_08, incr_12_09, incr_12_10, incr_12_11,
   incr_12_12, incr_12_13, incr_12_14, incr_12_15,

   incr_14_00, incr_14_01, incr_14_02, incr_14_03,
   incr_14_04, incr_14_05, incr_14_06, incr_14_07,
   incr_14_08, incr_14_09, incr_14_10, incr_14_11,
   incr_14_12, incr_14_13, incr_14_14, incr_14_15,

   incr_16_00, incr_16_01, incr_16_02, incr_16_03,
   incr_16_04, incr_16_05, incr_16_06, incr_16_07,
   incr_16_08, incr_16_09, incr_16_10, incr_16_11,
   incr_16_12, incr_16_13, incr_16_14, incr_16_15,

   incr_18_00, incr_18_01, incr_18_02, incr_18_03,
   incr_18_04, incr_18_05, incr_18_06, incr_18_07,
   incr_18_08, incr_18_09, incr_18_10, incr_18_11,
   incr_18_12, incr_18_13, incr_18_14, incr_18_15,

   incr_20_00, incr_20_01, incr_20_02, incr_20_03,
   incr_20_04, incr_20_05, incr_20_06, incr_20_07,
   incr_20_08, incr_20_09, incr_20_10, incr_20_11,
   incr_20_12, incr_20_13, incr_20_14, incr_20_15,

   incr_22_00, incr_22_01, incr_22_02, incr_22_03,
   incr_22_04, incr_22_05, incr_22_06, incr_22_07,
   incr_22_08, incr_22_09, incr_22_10, incr_22_11,
   incr_22_12, incr_22_13, incr_22_14, incr_22_15,

   incr_24_00, incr_24_01, incr_24_02, incr_24_03,
   incr_24_04, incr_24_05, incr_24_06, incr_24_07,
   incr_24_08, incr_24_09, incr_24_10, incr_24_11,
   incr_24_12, incr_24_13, incr_24_14, incr_24_15,

   incr_26_00, incr_26_01, incr_26_02, incr_26_03,
   incr_26_04, incr_26_05, incr_26_06, incr_26_07,
   incr_26_08, incr_26_09, incr_26_10, incr_26_11,
   incr_26_12, incr_26_13, incr_26_14, incr_26_15,

   incr_28_00, incr_28_01, incr_28_02, incr_28_03,
   incr_28_04, incr_28_05, incr_28_06, incr_28_07,
   incr_28_08, incr_28_09, incr_28_10, incr_28_11,
   incr_28_12, incr_28_13, incr_28_14, incr_28_15,

   incr_30_00, incr_30_01, incr_30_02, incr_30_03,
   incr_30_04, incr_30_05, incr_30_06, incr_30_07,
   incr_30_08, incr_30_09, incr_30_10, incr_30_11,
   incr_30_12, incr_30_13, incr_30_14, incr_30_15,

   incr_32_00, incr_32_01, incr_32_02, incr_32_03,
   incr_32_04, incr_32_05, incr_32_06, incr_32_07,
   incr_32_08, incr_32_09, incr_32_10, incr_32_11,
   incr_32_12, incr_32_13, incr_32_14, incr_32_15,
  };
  


static void
add_00 (UINT4 *in, UINT4 increment) {
  abort();
  return;
}



static void
add_02_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_02_01 (UINT4 *in, UINT4 increment) {
  *in += (increment << 2);
  return;
}

static void
add_02_02 (UINT4 *in, UINT4 increment) {
  *in += (increment << 4);
  return;
}

static void
add_02_03 (UINT4 *in, UINT4 increment) {
  *in += (increment << 6);
  return;
}

static void
add_02_04 (UINT4 *in, UINT4 increment) {
  *in += (increment << 8);
  return;
}

static void
add_02_05 (UINT4 *in, UINT4 increment) {
  *in += (increment << 10);
  return;
}

static void
add_02_06 (UINT4 *in, UINT4 increment) {
  *in += (increment << 12);
  return;
}

static void
add_02_07 (UINT4 *in, UINT4 increment) {
  *in += (increment << 14);
  return;
}

static void
add_02_08 (UINT4 *in, UINT4 increment) {
  *in += (increment << 16);
  return;
}

static void
add_02_09 (UINT4 *in, UINT4 increment) {
  *in += (increment << 18);
  return;
}

static void
add_02_10 (UINT4 *in, UINT4 increment) {
  *in += (increment << 20);
  return;
}

static void
add_02_11 (UINT4 *in, UINT4 increment) {
  *in += (increment << 22);
  return;
}

static void
add_02_12 (UINT4 *in, UINT4 increment) {
  *in += (increment << 24);
  return;
}

static void
add_02_13 (UINT4 *in, UINT4 increment) {
  *in += (increment << 26);
  return;
}

static void
add_02_14 (UINT4 *in, UINT4 increment) {
  *in += (increment << 28);
  return;
}

static void
add_02_15 (UINT4 *in, UINT4 increment) {
  *in += (increment << 30);
  return;
}


static void
add_04_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_04_01 (UINT4 *in, UINT4 increment) {
  *in += (increment << 4);
  return;
}

static void
add_04_02 (UINT4 *in, UINT4 increment) {
  *in += (increment << 8);
  return;
}

static void
add_04_03 (UINT4 *in, UINT4 increment) {
  *in += (increment << 12);
  return;
}

static void
add_04_04 (UINT4 *in, UINT4 increment) {
  *in += (increment << 16);
  return;
}

static void
add_04_05 (UINT4 *in, UINT4 increment) {
  *in += (increment << 20);
  return;
}

static void
add_04_06 (UINT4 *in, UINT4 increment) {
  *in += (increment << 24);
  return;
}

static void
add_04_07 (UINT4 *in, UINT4 increment) {
  *in += (increment << 28);
  return;
}

static void
add_04_08 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_04_09 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_04_10 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_04_11 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 12);
  return;
}

static void
add_04_12 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_04_13 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 20);
  return;
}

static void
add_04_14 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 24);
  return;
}

static void
add_04_15 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 28);
  return;
}


static void
add_06_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_06_01 (UINT4 *in, UINT4 increment) {
  *in += (increment << 6);
  return;
}

static void
add_06_02 (UINT4 *in, UINT4 increment) {
  *in += (increment << 12);
  return;
}

static void
add_06_03 (UINT4 *in, UINT4 increment) {
  *in += (increment << 18);
  return;
}

static void
add_06_04 (UINT4 *in, UINT4 increment) {
  *in += (increment << 24);
  return;
}

static void
add_06_05 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 6 ) ;
  carry = ( (prev + increment) >> (6 - 4) ); 
  *in += (increment << 30);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 ); */
  *in += carry;

  return;
}

static void
add_06_06 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_06_07 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 10);
  return;
}

static void
add_06_08 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_06_09 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 22);
  return;
}

static void
add_06_10 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 6 ) ;
  carry = ( (prev + increment) >> (6 - 2) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 ); */
  *in += carry;

  return;
}

static void
add_06_11 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 2);
  return;
}

static void
add_06_12 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_06_13 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 14);
  return;
}

static void
add_06_14 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 20);
  return;
}

static void
add_06_15 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 26);
  return;
}


static void
add_08_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_08_01 (UINT4 *in, UINT4 increment) {
  *in += (increment << 8);
  return;
}

static void
add_08_02 (UINT4 *in, UINT4 increment) {
  *in += (increment << 16);
  return;
}

static void
add_08_03 (UINT4 *in, UINT4 increment) {
  *in += (increment << 24);
  return;
}

static void
add_08_04 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_08_05 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_08_06 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_08_07 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 24);
  return;
}

static void
add_08_08 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_08_09 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_08_10 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_08_11 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 24);
  return;
}

static void
add_08_12 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_08_13 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_08_14 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_08_15 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 24);
  return;
}



static void
add_10_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_10_01 (UINT4 *in, UINT4 increment) {
  *in += (increment << 10);
  return;
}

static void
add_10_02 (UINT4 *in, UINT4 increment) {
  *in += (increment << 20);
  return;
}

static void
add_10_03 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 10 ) ;
  carry = ( (prev + increment) >> (10 - 8) );
  *in += (increment << 30);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 ); */
  *in += carry;

  return;
}

static void
add_10_04 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_10_05 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 18);
  return;
}

static void
add_10_06 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;
  
  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 10 ) ;
  carry = ( (prev + increment) >> (10 - 6) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 ); */
  *in += carry;

  return;
}

static void
add_10_07 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 6);
  return;
}

static void
add_10_08 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_10_09 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 10 ) ;
  carry = ( (prev + increment) >> (10 - 4) );
  *in += (increment << 26);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 ); */
  *in += carry;

  return;
}

static void
add_10_10 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_10_11 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 14);
  return;
}

static void
add_10_12 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 10 ) ;
  carry = ( (prev + increment) >> (10 - 2) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 ); */
  *in += carry;

  return;
}

static void
add_10_13 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 2);
  return;
}

static void
add_10_14 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 12);
  return;
}

static void
add_10_15 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 22);
  return;
}



static void
add_12_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_12_01 (UINT4 *in, UINT4 increment) {
  *in += (increment << 12);
  return;
}

static void
add_12_02 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  carry = ( (prev + increment) >> (12 - 4) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 ); */
  *in += carry;

  return;
}

static void
add_12_03 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_12_04 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_12_05 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  carry = ( (prev + increment) >> (12 - 8));
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 ); */
  *in += carry;

  return;
}

static void
add_12_06 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_12_07 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 20);
  return;
}

static void
add_12_08 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_12_09 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 12);
  return;
}

static void
add_12_10 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  carry = ( (prev + increment) >> (12 - 4) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 ); */
  *in += carry;

  return;
}

static void
add_12_11 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_12_12 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_12_13 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  carry = ( (prev + increment) >> (12 - 8));
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 ); */
  *in += carry;

  return;
}

static void
add_12_14 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_12_15 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 20);
  return;
}


static void
add_14_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_14_01 (UINT4 *in, UINT4 increment) {
  *in += (increment << 14);
  return;
}

static void
add_14_02 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  28  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 10) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 ); */
  *in += carry;

  return;
}

static void
add_14_03 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 10);
  return;
}

static void
add_14_04 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 6) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 ); */
  *in += carry;

  return;
}

static void
add_14_05 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 6);
  return;
}

static void
add_14_06 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 2) );
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 ); */
  *in += carry;

  return;
}

static void
add_14_07 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 2);
  return;
}

static void
add_14_08 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_14_09 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 12) );
  *in += (increment << 30);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 ); */
  *in += carry;

  return;
}

static void
add_14_10 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 12);
  return;
}

static void
add_14_11 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 8) );
  *in += (increment << 26);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 ); */
  *in += carry;

  return;
}

static void
add_14_12 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_14_13 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 4) );
  *in += (increment << 22);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 ); */
  *in += carry;

  return;
}

static void
add_14_14 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_14_15 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 18);
  return;
}


static void
add_16_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_16_01 (UINT4 *in, UINT4 increment) {
  *in += (increment << 16);
  return;
}

static void
add_16_02 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_16_03 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_16_04 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_16_05 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_16_06 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_16_07 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_16_08 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_16_09 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_16_10 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_16_11 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_16_12 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_16_13 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 16);
  return;
}

static void
add_16_14 (UINT4 *in, UINT4 increment) {
  in += 7 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_16_15 (UINT4 *in, UINT4 increment) {
  in += 7 * WORD_INCR;
  *in += (increment << 16);
  return;
}


static void
add_18_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_18_01 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  18  )   % (1U << 18 ) ;
  carry = ( (prev + increment) >> (18 - 4) );
  *in += (increment << 18);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 18 - 4 ); */
  *in += carry;

  return;
}

static void
add_18_02 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_18_03 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 18 ) ;
  carry = ( (prev + increment) >> (18 - 8) );
  *in += (increment << 22);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 18 - 8 ); */
  *in += carry;

  return;
}

static void
add_18_04 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_18_05 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 18 ) ;
  carry = ( (prev + increment) >> (18 - 12) );
  *in += (increment << 26);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 18 - 12 ); */
  *in += carry;

  return;
}

static void
add_18_06 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 12);
  return;
}

static void
add_18_07 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 18 ) ;
  carry = ( (prev + increment) >> (18 - 16) );
  *in += (increment << 30);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 18 - 16 ); */
  *in += carry;

  return;
}

static void
add_18_08 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 18 ) ;
  carry = ( (prev + increment) >> (18 - 2) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 18 - 2 ); */
  *in += carry;

  return;
}

static void
add_18_09 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 2);
  return;
}

static void
add_18_10 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 18 ) ;
  carry = ( (prev + increment) >> (18 - 6) );
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 6 ))<<( 18 - 6 ); */
  *in += carry;

  return;
}

static void
add_18_11 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 6);
  return;
}

static void
add_18_12 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 18 ) ;
  carry = ( (prev + increment) >> (18 - 10) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 10 ))<<( 18 - 10 ); */
  *in += carry;

  return;
}

static void
add_18_13 (UINT4 *in, UINT4 increment) {
  in += 7 * WORD_INCR;
  *in += (increment << 10);
  return;
}

static void
add_18_14 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 18 ) ;
  carry = ( (prev + increment) >> (18 - 14) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 14 ))<<( 18 - 14 ); */
  *in += carry;

  return;
}

static void
add_18_15 (UINT4 *in, UINT4 increment) {
  in += 8 * WORD_INCR;
  *in += (increment << 14);
  return;
}


static void
add_20_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_20_01 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  carry = ( (prev + increment) >> (20 - 8) );
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 ); */
  *in += carry;

  return;
}

static void
add_20_02 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_20_03 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  carry = ( (prev + increment) >> (20 - 16) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 ); */
  *in += carry;

  return;
}

static void
add_20_04 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  carry = ( (prev + increment) >> (20 - 4) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 ); */
  *in += carry;

  return;
}

static void
add_20_05 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_20_06 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  carry = ( (prev + increment) >> (20 - 12) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 ); */
  *in += carry;

  return;
}

static void
add_20_07 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 12);
  return;
}

static void
add_20_08 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_20_09 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  carry = ( (prev + increment) >> (20 - 8) );
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 ); */
  *in += carry;

  return;
}

static void
add_20_10 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_20_11 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  carry = ( (prev + increment) >> (20 - 16) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 ); */
  *in += carry;

  return;
}

static void
add_20_12 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  carry = ( (prev + increment) >> (20 - 4) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 ); */
  *in += carry;

  return;
}

static void
add_20_13 (UINT4 *in, UINT4 increment) {
  in += 8 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_20_14 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  carry = ( (prev + increment) >> (20 - 12) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 ); */
  *in += carry;

  return;
}

static void
add_20_15 (UINT4 *in, UINT4 increment) {
  in += 9 * WORD_INCR;
  *in += (increment << 12);
  return;
}


static void
add_22_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_22_01 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  22  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 12) );
  *in += (increment << 22);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 22 - 12 ); */
  *in += carry;

  return;
}

static void
add_22_02 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 2) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 22 - 2 ); */
  *in += carry;

  return;
}

static void
add_22_03 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 2);
  return;
}

static void
add_22_04 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 14) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 14 ))<<( 22 - 14 ); */
  *in += carry;

  return;
}

static void
add_22_05 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 4) );
  *in += (increment << 14);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 22 - 4 ); */
  *in += carry;

  return;
}

static void
add_22_06 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_22_07 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 16) );
  *in += (increment << 26);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 22 - 16 ); */
  *in += carry;

  return;
}

static void
add_22_08 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 6) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 6 ))<<( 22 - 6 ); */
  *in += carry;

  return;
}

static void
add_22_09 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 6);
  return;
}

static void
add_22_10 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 18) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 18 ))<<( 22 - 18 ); */
  *in += carry;

  return;
}

static void
add_22_11 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 8) );
  *in += (increment << 18);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 22 - 8 ); */
  *in += carry;

  return;
}

static void
add_22_12 (UINT4 *in, UINT4 increment) {
  in += 8 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_22_13 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 20) );
  *in += (increment << 30);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 20 ))<<( 22 - 20 ); */
  *in += carry;

  return;
}

static void
add_22_14 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 22 ) ;
  carry = ( (prev + increment) >> (22 - 10) );
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 10 ))<<( 22 - 10 ); */
  *in += carry;

  return;
}

static void
add_22_15 (UINT4 *in, UINT4 increment) {
  in += 10 * WORD_INCR;
  *in += (increment << 10);
  return;
}



static void
add_24_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_24_01 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev + increment) >> (24 - 16) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 ); */
  *in += carry;

  return;
}

static void
add_24_02 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev + increment) >> (24 - 8) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 ); */
  *in += carry;

  return;
}

static void
add_24_03 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_24_04 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_24_05 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev + increment) >> (24 - 16) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 ); */
  *in += carry;

  return;
}

static void
add_24_06 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev + increment) >> (24 - 8) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 ); */
  *in += carry;

  return;
}

static void
add_24_07 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_24_08 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_24_09 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev + increment) >> (24 - 16) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 ); */
  *in += carry;

  return;
}

static void
add_24_10 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev + increment) >> (24 - 8) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 ); */
  *in += carry;

  return;
}

static void
add_24_11 (UINT4 *in, UINT4 increment) {
  in += 8 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_24_12 (UINT4 *in, UINT4 increment) {
  in += 9 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_24_13 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev + increment) >> (24 - 16) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 ); */
  *in += carry;

  return;
}

static void
add_24_14 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev + increment) >> (24 - 8) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 ); */
  *in += carry;

  return;
}

static void
add_24_15 (UINT4 *in, UINT4 increment) {
  in += 11 * WORD_INCR;
  *in += (increment << 8);
  return;
}



static void
add_26_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_26_01 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  26  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 20) );
  *in += (increment << 26);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 20 ))<<( 26 - 20 ); */
  *in += carry;

  return;
}

static void
add_26_02 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 14) );
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 14 ))<<( 26 - 14 ); */
  *in += carry;

  return;
}

static void
add_26_03 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 8) );
  *in += (increment << 14);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 26 - 8 ); */
  *in += carry;

  return;
}

static void
add_26_04 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 2) );
  *in += (increment << 8);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 26 - 2 ); */
  *in += carry;

  return;
}

static void
add_26_05 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 2);
  return;
}

static void
add_26_06 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 22) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 22 ))<<( 26 - 22 ); */
  *in += carry;

  return;
}

static void
add_26_07 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 16) );
  *in += (increment << 22);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 26 - 16 ); */
  *in += carry;

  return;
}

static void
add_26_08 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 10) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 10 ))<<( 26 - 10 ); */
  *in += carry;

  return;
}

static void
add_26_09 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 4) );
  *in += (increment << 10);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 26 - 4 ); */
  *in += carry;

  return;
}

static void
add_26_10 (UINT4 *in, UINT4 increment) {
  in += 8 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_26_11 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 24) );
  *in += (increment << 30);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 24 ))<<( 26 - 24 ); */
  *in += carry;

  return;
}

static void
add_26_12 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 18) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 18 ))<<( 26 - 18 ); */
  *in += carry;

  return;
}

static void
add_26_13 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 12) );
  *in += (increment << 18);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 26 - 12 ); */
  *in += carry;

  return;
}

static void
add_26_14 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 26 ) ;
  carry = ( (prev + increment) >> (26 - 6) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 6 ))<<( 26 - 6 ); */
  *in += carry;

  return;
}

static void
add_26_15 (UINT4 *in, UINT4 increment) {
  in += 12 * WORD_INCR;
  *in += (increment << 6);
  return;
}


static void
add_28_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_28_01 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 24) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 ); */
  *in += carry;

  return;
}

static void
add_28_02 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 20) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 ); */
  *in += carry;

  return;
}

static void
add_28_03 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 16) );
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 ); */
  *in += carry;

  return;
}

static void
add_28_04 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 12) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 ); */
  *in += carry;

  return;
}

static void
add_28_05 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 8) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 ); */
  *in += carry;

  return;
}

static void
add_28_06 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 4) );
  *in += (increment << 8);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 ); */
  *in += carry;

  return;
}

static void
add_28_07 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_28_08 (UINT4 *in, UINT4 increment) {
  in += 7 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_28_09 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 24) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 ); */
  *in += carry;

  return;
}

static void
add_28_10 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 20) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 ); */
  *in += carry;

  return;
}

static void
add_28_11 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 16) );
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 ); */
  *in += carry;

  return;
}

static void
add_28_12 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 12) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 ); */
  *in += carry;

  return;
}

static void
add_28_13 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 8) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 ); */
  *in += carry;

  return;
}

static void
add_28_14 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 12 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  carry = ( (prev + increment) >> (28 - 4) );
  *in += (increment << 8);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 ); */
  *in += carry;

  return;
}

static void
add_28_15 (UINT4 *in, UINT4 increment) {
  in += 13 * WORD_INCR;
  *in += (increment << 4);
  return;
}


static void
add_30_00 (UINT4 *in, UINT4 increment) {
  *in += (increment << 0);
  return;
}

static void
add_30_01 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 28) );
  *in += (increment << 30);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 28 ))<<( 30 - 28 ); */
  *in += carry;

  return;
}

static void
add_30_02 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 26) );
  *in += (increment << 28);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 26 ))<<( 30 - 26 ); */
  *in += carry;

  return;
}

static void
add_30_03 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 24) );
  *in += (increment << 26);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 24 ))<<( 30 - 24 ); */
  *in += carry;

  return;
}

static void
add_30_04 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 22) );
  *in += (increment << 24);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 22 ))<<( 30 - 22 ); */
  *in += carry;

  return;
}

static void
add_30_05 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 20) );
  *in += (increment << 22);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 20 ))<<( 30 - 20 ); */
  *in += carry;

  return;
}

static void
add_30_06 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 18));
  *in += (increment << 20);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 18 ))<<( 30 - 18 ); */
  *in += carry;

  return;
}

static void
add_30_07 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 16) );
  *in += (increment << 18);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 16 ))<<( 30 - 16 ); */
  *in += carry;

  return;
}

static void
add_30_08 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 14) );
  *in += (increment << 16);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 14 ))<<( 30 - 14 ); */
  *in += carry;

  return;
}

static void
add_30_09 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 12) );
  *in += (increment << 14);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 30 - 12 ); */
  *in += carry;

  return;
}

static void
add_30_10 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 10) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 10 ))<<( 30 - 10 ); */
  *in += carry;

  return;
}

static void
add_30_11 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 8) );
  *in += (increment << 10);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 30 - 8 ); */
  *in += carry;

  return;
}

static void
add_30_12 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 6) );
  *in += (increment << 8);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 6 ))<<( 30 - 6 ); */
  *in += carry;

  return;
}

static void
add_30_13 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 12 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 4) );
  *in += (increment << 6);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 30 - 4 ); */
  *in += carry;

  return;
}

static void
add_30_14 (UINT4 *in, UINT4 increment) {
  UINT4 prev, carry;

  in += 13 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 30 ) ;
  carry = ( (prev + increment) >> (30 - 2) );
  *in += (increment << 4);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 30 - 2 ); */
  *in += carry;

  return;
}

static void
add_30_15 (UINT4 *in, UINT4 increment) {
  in += 14 * WORD_INCR;
  *in += (increment << 2);
  return;
}


static void
add_32_00 (UINT4 *in, UINT4 increment) {
  *in += increment;
  return;
}

static void
add_32_01 (UINT4 *in, UINT4 increment) {
  in += 1 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_02 (UINT4 *in, UINT4 increment) {
  in += 2 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_03 (UINT4 *in, UINT4 increment) {
  in += 3 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_04 (UINT4 *in, UINT4 increment) {
  in += 4 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_05 (UINT4 *in, UINT4 increment) {
  in += 5 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_06 (UINT4 *in, UINT4 increment) {
  in += 6 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_07 (UINT4 *in, UINT4 increment) {
  in += 7 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_08 (UINT4 *in, UINT4 increment) {
  in += 8 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_09 (UINT4 *in, UINT4 increment) {
  in += 9 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_10 (UINT4 *in, UINT4 increment) {
  in += 10 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_11 (UINT4 *in, UINT4 increment) {
  in += 11 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_12 (UINT4 *in, UINT4 increment) {
  in += 12 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_13 (UINT4 *in, UINT4 increment) {
  in += 13 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_14 (UINT4 *in, UINT4 increment) {
  in += 14 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_32_15 (UINT4 *in, UINT4 increment) {
  in += 15 * WORD_INCR;
  *in += increment;
  return;
}


typedef void (*Adder_T) (UINT4 *, UINT4);

static Adder_T adder_table[272] =
  {add_00, add_00, add_00, add_00,
   add_00, add_00, add_00, add_00,
   add_00, add_00, add_00, add_00,
   add_00, add_00, add_00, add_00,

   add_02_00, add_02_01, add_02_02, add_02_03,
   add_02_04, add_02_05, add_02_06, add_02_07,
   add_02_08, add_02_09, add_02_10, add_02_11,
   add_02_12, add_02_13, add_02_14, add_02_15,

   add_04_00, add_04_01, add_04_02, add_04_03,
   add_04_04, add_04_05, add_04_06, add_04_07,
   add_04_08, add_04_09, add_04_10, add_04_11,
   add_04_12, add_04_13, add_04_14, add_04_15,

   add_06_00, add_06_01, add_06_02, add_06_03,
   add_06_04, add_06_05, add_06_06, add_06_07,
   add_06_08, add_06_09, add_06_10, add_06_11,
   add_06_12, add_06_13, add_06_14, add_06_15,

   add_08_00, add_08_01, add_08_02, add_08_03,
   add_08_04, add_08_05, add_08_06, add_08_07,
   add_08_08, add_08_09, add_08_10, add_08_11,
   add_08_12, add_08_13, add_08_14, add_08_15,

   add_10_00, add_10_01, add_10_02, add_10_03,
   add_10_04, add_10_05, add_10_06, add_10_07,
   add_10_08, add_10_09, add_10_10, add_10_11,
   add_10_12, add_10_13, add_10_14, add_10_15,

   add_12_00, add_12_01, add_12_02, add_12_03,
   add_12_04, add_12_05, add_12_06, add_12_07,
   add_12_08, add_12_09, add_12_10, add_12_11,
   add_12_12, add_12_13, add_12_14, add_12_15,

   add_14_00, add_14_01, add_14_02, add_14_03,
   add_14_04, add_14_05, add_14_06, add_14_07,
   add_14_08, add_14_09, add_14_10, add_14_11,
   add_14_12, add_14_13, add_14_14, add_14_15,

   add_16_00, add_16_01, add_16_02, add_16_03,
   add_16_04, add_16_05, add_16_06, add_16_07,
   add_16_08, add_16_09, add_16_10, add_16_11,
   add_16_12, add_16_13, add_16_14, add_16_15,

   add_18_00, add_18_01, add_18_02, add_18_03,
   add_18_04, add_18_05, add_18_06, add_18_07,
   add_18_08, add_18_09, add_18_10, add_18_11,
   add_18_12, add_18_13, add_18_14, add_18_15,

   add_20_00, add_20_01, add_20_02, add_20_03,
   add_20_04, add_20_05, add_20_06, add_20_07,
   add_20_08, add_20_09, add_20_10, add_20_11,
   add_20_12, add_20_13, add_20_14, add_20_15,

   add_22_00, add_22_01, add_22_02, add_22_03,
   add_22_04, add_22_05, add_22_06, add_22_07,
   add_22_08, add_22_09, add_22_10, add_22_11,
   add_22_12, add_22_13, add_22_14, add_22_15,

   add_24_00, add_24_01, add_24_02, add_24_03,
   add_24_04, add_24_05, add_24_06, add_24_07,
   add_24_08, add_24_09, add_24_10, add_24_11,
   add_24_12, add_24_13, add_24_14, add_24_15,

   add_26_00, add_26_01, add_26_02, add_26_03,
   add_26_04, add_26_05, add_26_06, add_26_07,
   add_26_08, add_26_09, add_26_10, add_26_11,
   add_26_12, add_26_13, add_26_14, add_26_15,

   add_28_00, add_28_01, add_28_02, add_28_03,
   add_28_04, add_28_05, add_28_06, add_28_07,
   add_28_08, add_28_09, add_28_10, add_28_11,
   add_28_12, add_28_13, add_28_14, add_28_15,

   add_30_00, add_30_01, add_30_02, add_30_03,
   add_30_04, add_30_05, add_30_06, add_30_07,
   add_30_08, add_30_09, add_30_10, add_30_11,
   add_30_12, add_30_13, add_30_14, add_30_15,

   add_32_00, add_32_01, add_32_02, add_32_03,
   add_32_04, add_32_05, add_32_06, add_32_07,
   add_32_08, add_32_09, add_32_10, add_32_11,
   add_32_12, add_32_13, add_32_14, add_32_15,
  };



#if 0
/* No need currently for a subtraction routine */

static void
sub_00 (UINT4 *in, UINT4 decrement) {
  abort();
  return;
}



static void
sub_02_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_02_01 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 2);
  return;
}

static void
sub_02_02 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 4);
  return;
}

static void
sub_02_03 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 6);
  return;
}

static void
sub_02_04 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 8);
  return;
}

static void
sub_02_05 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 10);
  return;
}

static void
sub_02_06 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 12);
  return;
}

static void
sub_02_07 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 14);
  return;
}

static void
sub_02_08 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 16);
  return;
}

static void
sub_02_09 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 18);
  return;
}

static void
sub_02_10 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 20);
  return;
}

static void
sub_02_11 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 22);
  return;
}

static void
sub_02_12 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 24);
  return;
}

static void
sub_02_13 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 26);
  return;
}

static void
sub_02_14 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 28);
  return;
}

static void
sub_02_15 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 30);
  return;
}


static void
sub_04_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_04_01 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 4);
  return;
}

static void
sub_04_02 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 8);
  return;
}

static void
sub_04_03 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 12);
  return;
}

static void
sub_04_04 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 16);
  return;
}

static void
sub_04_05 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 20);
  return;
}

static void
sub_04_06 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 24);
  return;
}

static void
sub_04_07 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 28);
  return;
}

static void
sub_04_08 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_04_09 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_04_10 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_04_11 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 12);
  return;
}

static void
sub_04_12 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_04_13 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 20);
  return;
}

static void
sub_04_14 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 24);
  return;
}

static void
sub_04_15 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 28);
  return;
}


static void
sub_06_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_06_01 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 6);
  return;
}

static void
sub_06_02 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 12);
  return;
}

static void
sub_06_03 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 18);
  return;
}

static void
sub_06_04 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 24);
  return;
}

static void
sub_06_05 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 6 ) ;
  carry = ( (prev - decrement) >> (6 - 4) ) & 0x1;
  *in -= (decrement << 30);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 6 - 4 )) + carry;

  return;
}

static void
sub_06_06 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_06_07 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 10);
  return;
}

static void
sub_06_08 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_06_09 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 22);
  return;
}

static void
sub_06_10 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 6 ) ;
  carry = ( (prev - decrement) >> (6 - 2) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 6 - 2 )) + carry;

  return;
}

static void
sub_06_11 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 2);
  return;
}

static void
sub_06_12 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_06_13 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 14);
  return;
}

static void
sub_06_14 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 20);
  return;
}

static void
sub_06_15 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 26);
  return;
}


static void
sub_08_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_08_01 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 8);
  return;
}

static void
sub_08_02 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 16);
  return;
}

static void
sub_08_03 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 24);
  return;
}

static void
sub_08_04 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_08_05 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_08_06 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_08_07 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 24);
  return;
}

static void
sub_08_08 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_08_09 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_08_10 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_08_11 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 24);
  return;
}

static void
sub_08_12 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_08_13 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_08_14 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_08_15 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 24);
  return;
}



static void
sub_10_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_10_01 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 10);
  return;
}

static void
sub_10_02 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 20);
  return;
}

static void
sub_10_03 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 10 ) ;
  carry = ( (prev - decrement) >> (10 - 8) ) & 0x1;
  *in -= (decrement << 30);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 10 - 8 )) + carry;

  return;
}

static void
sub_10_04 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_10_05 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 18);
  return;
}

static void
sub_10_06 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;
  
  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 10 ) ;
  carry = ( (prev - decrement) >> (10 - 6) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 10 - 6 )) + carry;

  return;
}

static void
sub_10_07 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 6);
  return;
}

static void
sub_10_08 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_10_09 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 10 ) ;
  carry = ( (prev - decrement) >> (10 - 4) ) & 0x1;
  *in -= (decrement << 26);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 10 - 4 )) + carry;

  return;
}

static void
sub_10_10 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_10_11 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 14);
  return;
}

static void
sub_10_12 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 10 ) ;
  carry = ( (prev - decrement) >> (10 - 2) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 10 - 2 )) + carry;

  return;
}

static void
sub_10_13 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 2);
  return;
}

static void
sub_10_14 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 12);
  return;
}

static void
sub_10_15 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 22);
  return;
}



static void
sub_12_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_12_01 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 12);
  return;
}

static void
sub_12_02 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  carry = ( (prev - decrement) >> (12 - 4) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 12 - 4 )) + carry;

  return;
}

static void
sub_12_03 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_12_04 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_12_05 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  carry = ( (prev - decrement) >> (12 - 8)) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 12 - 8 )) + carry;

  return;
}

static void
sub_12_06 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_12_07 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 20);
  return;
}

static void
sub_12_08 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_12_09 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 12);
  return;
}

static void
sub_12_10 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  carry = ( (prev - decrement) >> (12 - 4) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 12 - 4 )) + carry;

  return;
}

static void
sub_12_11 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_12_12 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_12_13 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  carry = ( (prev - decrement) >> (12 - 8)) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 12 - 8 )) + carry;

  return;
}

static void
sub_12_14 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_12_15 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= (decrement << 20);
  return;
}


static void
sub_14_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_14_01 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 14);
  return;
}

static void
sub_14_02 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  28  )   % (1U << 14 ) ;
  carry = ( (prev - decrement) >> (14 - 10) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 14 - 10 )) + carry;

  return;
}

static void
sub_14_03 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 10);
  return;
}

static void
sub_14_04 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 14 ) ;
  carry = ( (prev - decrement) >> (14 - 6) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 14 - 6 )) + carry;

  return;
}

static void
sub_14_05 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 6);
  return;
}

static void
sub_14_06 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 14 ) ;
  carry = ( (prev - decrement) >> (14 - 2) ) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 14 - 2 )) + carry;

  return;
}

static void
sub_14_07 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 2);
  return;
}

static void
sub_14_08 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_14_09 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 14 ) ;
  carry = ( (prev - decrement) >> (14 - 12) ) & 0x1;
  *in -= (decrement << 30);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 14 - 12 )) + carry;

  return;
}

static void
sub_14_10 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 12);
  return;
}

static void
sub_14_11 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 14 ) ;
  carry = ( (prev - decrement) >> (14 - 8) ) & 0x1;
  *in -= (decrement << 26);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 14 - 8 )) + carry;

  return;
}

static void
sub_14_12 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_14_13 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 14 ) ;
  carry = ( (prev - decrement) >> (14 - 4) ) & 0x1;
  *in -= (decrement << 22);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 14 - 4 )) + carry;

  return;
}

static void
sub_14_14 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_14_15 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 18);
  return;
}


static void
sub_16_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_16_01 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 16);
  return;
}

static void
sub_16_02 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_16_03 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_16_04 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_16_05 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_16_06 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_16_07 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_16_08 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_16_09 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_16_10 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_16_11 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_16_12 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_16_13 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}

static void
sub_16_14 (UINT4 *in, UINT4 decrement) {
  in += 7 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_16_15 (UINT4 *in, UINT4 decrement) {
  in += 7 * WORD_INCR;
  *in -= (decrement << 16);
  return;
}


static void
sub_18_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_18_01 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  18  )   % (1U << 18 ) ;
  carry = ( (prev - decrement) >> (18 - 4) ) & 0x1;
  *in -= (decrement << 18);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 18 - 4 )) + carry;

  return;
}

static void
sub_18_02 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_18_03 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 18 ) ;
  carry = ( (prev - decrement) >> (18 - 8) ) & 0x1;
  *in -= (decrement << 22);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 18 - 8 )) + carry;

  return;
}

static void
sub_18_04 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_18_05 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 18 ) ;
  carry = ( (prev - decrement) >> (18 - 12) ) & 0x1;
  *in -= (decrement << 26);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 18 - 12 )) + carry;

  return;
}

static void
sub_18_06 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 12);
  return;
}

static void
sub_18_07 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 18 ) ;
  carry = ( (prev - decrement) >> (18 - 16) ) & 0x1;
  *in -= (decrement << 30);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 18 - 16 )) + carry;

  return;
}

static void
sub_18_08 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 18 ) ;
  carry = ( (prev - decrement) >> (18 - 2) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 18 - 2 )) + carry;

  return;
}

static void
sub_18_09 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= (decrement << 2);
  return;
}

static void
sub_18_10 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 18 ) ;
  carry = ( (prev - decrement) >> (18 - 6) ) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 18 - 6 )) + carry;

  return;
}

static void
sub_18_11 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 6);
  return;
}

static void
sub_18_12 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 18 ) ;
  carry = ( (prev - decrement) >> (18 - 10) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 18 - 10 )) + carry;

  return;
}

static void
sub_18_13 (UINT4 *in, UINT4 decrement) {
  in += 7 * WORD_INCR;
  *in -= (decrement << 10);
  return;
}

static void
sub_18_14 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 18 ) ;
  carry = ( (prev - decrement) >> (18 - 14) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 18 - 14 )) + carry;

  return;
}

static void
sub_18_15 (UINT4 *in, UINT4 decrement) {
  in += 8 * WORD_INCR;
  *in -= (decrement << 14);
  return;
}


static void
sub_20_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_20_01 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  carry = ( (prev - decrement) >> (20 - 8) ) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 20 - 8 )) + carry;

  return;
}

static void
sub_20_02 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_20_03 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  carry = ( (prev - decrement) >> (20 - 16) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 20 - 16 )) + carry;

  return;
}

static void
sub_20_04 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  carry = ( (prev - decrement) >> (20 - 4) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 20 - 4 )) + carry;

  return;
}

static void
sub_20_05 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_20_06 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  carry = ( (prev - decrement) >> (20 - 12) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 20 - 12 )) + carry;

  return;
}

static void
sub_20_07 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 12);
  return;
}

static void
sub_20_08 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_20_09 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  carry = ( (prev - decrement) >> (20 - 8) ) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 20 - 8 )) + carry;

  return;
}

static void
sub_20_10 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_20_11 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  carry = ( (prev - decrement) >> (20 - 16) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 20 - 16 )) + carry;

  return;
}

static void
sub_20_12 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  carry = ( (prev - decrement) >> (20 - 4) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 20 - 4 )) + carry;

  return;
}

static void
sub_20_13 (UINT4 *in, UINT4 decrement) {
  in += 8 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_20_14 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  carry = ( (prev - decrement) >> (20 - 12) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 20 - 12 )) + carry;

  return;
}

static void
sub_20_15 (UINT4 *in, UINT4 decrement) {
  in += 9 * WORD_INCR;
  *in -= (decrement << 12);
  return;
}


static void
sub_22_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_22_01 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  22  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 12) ) & 0x1;
  *in -= (decrement << 22);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 12 )) + carry;

  return;
}

static void
sub_22_02 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 2) ) & 0x1;
  *in -= (decrement << 12);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 2 )) + carry;

  return;
}

static void
sub_22_03 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 2);
  return;
}

static void
sub_22_04 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 14) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 14 )) + carry;

  return;
}

static void
sub_22_05 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 4) ) & 0x1;
  *in -= (decrement << 14);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 4 )) + carry;

  return;
}

static void
sub_22_06 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_22_07 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 16) ) & 0x1;
  *in -= (decrement << 26);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 16 )) + carry;

  return;
}

static void
sub_22_08 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 6) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 6 )) + carry;

  return;
}

static void
sub_22_09 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 6);
  return;
}

static void
sub_22_10 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 18) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 18 )) + carry;

  return;
}

static void
sub_22_11 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 8) ) & 0x1;
  *in -= (decrement << 18);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 8 )) + carry;

  return;
}

static void
sub_22_12 (UINT4 *in, UINT4 decrement) {
  in += 8 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_22_13 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 20) ) & 0x1;
  *in -= (decrement << 30);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 20 )) + carry;

  return;
}

static void
sub_22_14 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 22 ) ;
  carry = ( (prev - decrement) >> (22 - 10) ) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 22 - 10 )) + carry;

  return;
}

static void
sub_22_15 (UINT4 *in, UINT4 decrement) {
  in += 10 * WORD_INCR;
  *in -= (decrement << 10);
  return;
}



static void
sub_24_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_24_01 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev - decrement) >> (24 - 16) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 24 - 16 )) + carry;

  return;
}

static void
sub_24_02 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev - decrement) >> (24 - 8) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 24 - 8 )) + carry;

  return;
}

static void
sub_24_03 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_24_04 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_24_05 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev - decrement) >> (24 - 16) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 24 - 16 )) + carry;

  return;
}

static void
sub_24_06 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev - decrement) >> (24 - 8) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 24 - 8 )) + carry;

  return;
}

static void
sub_24_07 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_24_08 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_24_09 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev - decrement) >> (24 - 16) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 24 - 16 )) + carry;

  return;
}

static void
sub_24_10 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev - decrement) >> (24 - 8) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 24 - 8 )) + carry;

  return;
}

static void
sub_24_11 (UINT4 *in, UINT4 decrement) {
  in += 8 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}

static void
sub_24_12 (UINT4 *in, UINT4 decrement) {
  in += 9 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_24_13 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  carry = ( (prev - decrement) >> (24 - 16) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 24 - 16 )) + carry;

  return;
}

static void
sub_24_14 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  carry = ( (prev - decrement) >> (24 - 8) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 24 - 8 )) + carry;

  return;
}

static void
sub_24_15 (UINT4 *in, UINT4 decrement) {
  in += 11 * WORD_INCR;
  *in -= (decrement << 8);
  return;
}



static void
sub_26_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_26_01 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  26  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 20) ) & 0x1;
  *in -= (decrement << 26);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 20 )) + carry;

  return;
}

static void
sub_26_02 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 14) ) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 14 )) + carry;

  return;
}

static void
sub_26_03 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 8) ) & 0x1;
  *in -= (decrement << 14);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 8 )) + carry;

  return;
}

static void
sub_26_04 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 2) ) & 0x1;
  *in -= (decrement << 8);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 2 )) + carry;

  return;
}

static void
sub_26_05 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= (decrement << 2);
  return;
}

static void
sub_26_06 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 22) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 22 )) + carry;

  return;
}

static void
sub_26_07 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 16) ) & 0x1;
  *in -= (decrement << 22);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 16 )) + carry;

  return;
}

static void
sub_26_08 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 10) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 10 )) + carry;

  return;
}

static void
sub_26_09 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 4) ) & 0x1;
  *in -= (decrement << 10);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 4 )) + carry;

  return;
}

static void
sub_26_10 (UINT4 *in, UINT4 decrement) {
  in += 8 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_26_11 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  30  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 24) ) & 0x1;
  *in -= (decrement << 30);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 24 )) + carry;

  return;
}

static void
sub_26_12 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 18) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 18 )) + carry;

  return;
}

static void
sub_26_13 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 12) ) & 0x1;
  *in -= (decrement << 18);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 12 )) + carry;

  return;
}

static void
sub_26_14 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 26 ) ;
  carry = ( (prev - decrement) >> (26 - 6) ) & 0x1;
  *in -= (decrement << 12);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 26 - 6 )) + carry;

  return;
}

static void
sub_26_15 (UINT4 *in, UINT4 decrement) {
  in += 12 * WORD_INCR;
  *in -= (decrement << 6);
  return;
}


static void
sub_28_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_28_01 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 24) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 24 )) + carry;

  return;
}

static void
sub_28_02 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 20) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 20 )) + carry;

  return;
}

static void
sub_28_03 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 16) ) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 16 )) + carry;

  return;
}

static void
sub_28_04 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 12) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 12 )) + carry;

  return;
}

static void
sub_28_05 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 8) ) & 0x1;
  *in -= (decrement << 12);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 8 )) + carry;

  return;
}

static void
sub_28_06 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 4) ) & 0x1;
  *in -= (decrement << 8);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 4 )) + carry;

  return;
}

static void
sub_28_07 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}

static void
sub_28_08 (UINT4 *in, UINT4 decrement) {
  in += 7 * WORD_INCR;
  *in -= (decrement << 0);
  return;
}

static void
sub_28_09 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 24) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 24 )) + carry;

  return;
}

static void
sub_28_10 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 20) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 20 )) + carry;

  return;
}

static void
sub_28_11 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 16) ) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 16 )) + carry;

  return;
}

static void
sub_28_12 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 12) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 12 )) + carry;

  return;
}

static void
sub_28_13 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 8) ) & 0x1;
  *in -= (decrement << 12);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 8 )) + carry;

  return;
}

static void
sub_28_14 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 12 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  carry = ( (prev - decrement) >> (28 - 4) ) & 0x1;
  *in -= (decrement << 8);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 28 - 4 )) + carry;

  return;
}

static void
sub_28_15 (UINT4 *in, UINT4 decrement) {
  in += 13 * WORD_INCR;
  *in -= (decrement << 4);
  return;
}


static void
sub_30_00 (UINT4 *in, UINT4 decrement) {
  *in -= (decrement << 0);
  return;
}

static void
sub_30_01 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  prev = ( CONVERT(*in) >>  30  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 28) ) & 0x1;
  *in -= (decrement << 30);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 28 )) + carry;

  return;
}

static void
sub_30_02 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  28  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 26) ) & 0x1;
  *in -= (decrement << 28);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 26 )) + carry;

  return;
}

static void
sub_30_03 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  26  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 24) ) & 0x1;
  *in -= (decrement << 26);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 24 )) + carry;

  return;
}

static void
sub_30_04 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  24  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 22) ) & 0x1;
  *in -= (decrement << 24);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 22 )) + carry;

  return;
}

static void
sub_30_05 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  22  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 20) ) & 0x1;
  *in -= (decrement << 22);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 20 )) + carry;

  return;
}

static void
sub_30_06 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  20  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 18)) & 0x1;
  *in -= (decrement << 20);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 18 )) + carry;

  return;
}

static void
sub_30_07 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  18  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 16) ) & 0x1;
  *in -= (decrement << 18);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 16 )) + carry;

  return;
}

static void
sub_30_08 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 7 * WORD_INCR;
  prev = ( CONVERT(*in) >>  16  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 14) ) & 0x1;
  *in -= (decrement << 16);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 14 )) + carry;

  return;
}

static void
sub_30_09 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 8 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 12) ) & 0x1;
  *in -= (decrement << 14);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 12 )) + carry;

  return;
}

static void
sub_30_10 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 9 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 10) ) & 0x1;
  *in -= (decrement << 12);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 10 )) + carry;

  return;
}

static void
sub_30_11 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 10 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 8) ) & 0x1;
  *in -= (decrement << 10);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 8 )) + carry;

  return;
}

static void
sub_30_12 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 11 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 6) ) & 0x1;
  *in -= (decrement << 8);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 6 )) + carry;

  return;
}

static void
sub_30_13 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 12 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 4) ) & 0x1;
  *in -= (decrement << 6);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 4 )) + carry;

  return;
}

static void
sub_30_14 (UINT4 *in, UINT4 decrement) {
  UINT4 prev, carry;

  in += 13 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 30 ) ;
  carry = ( (prev - decrement) >> (30 - 2) ) & 0x1;
  *in -= (decrement << 4);

  in += 1 * WORD_INCR;
  *in -= (decrement >> ( 30 - 2 )) + carry;

  return;
}

static void
sub_30_15 (UINT4 *in, UINT4 decrement) {
  in += 14 * WORD_INCR;
  *in -= (decrement << 2);
  return;
}


static void
sub_32_00 (UINT4 *in, UINT4 decrement) {
  *in -= decrement;
  return;
}

static void
sub_32_01 (UINT4 *in, UINT4 decrement) {
  in += 1 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_02 (UINT4 *in, UINT4 decrement) {
  in += 2 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_03 (UINT4 *in, UINT4 decrement) {
  in += 3 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_04 (UINT4 *in, UINT4 decrement) {
  in += 4 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_05 (UINT4 *in, UINT4 decrement) {
  in += 5 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_06 (UINT4 *in, UINT4 decrement) {
  in += 6 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_07 (UINT4 *in, UINT4 decrement) {
  in += 7 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_08 (UINT4 *in, UINT4 decrement) {
  in += 8 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_09 (UINT4 *in, UINT4 decrement) {
  in += 9 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_10 (UINT4 *in, UINT4 decrement) {
  in += 10 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_11 (UINT4 *in, UINT4 decrement) {
  in += 11 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_12 (UINT4 *in, UINT4 decrement) {
  in += 12 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_13 (UINT4 *in, UINT4 decrement) {
  in += 13 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_14 (UINT4 *in, UINT4 decrement) {
  in += 14 * WORD_INCR;
  *in -= decrement;
  return;
}

static void
sub_32_15 (UINT4 *in, UINT4 decrement) {
  in += 15 * WORD_INCR;
  *in -= decrement;
  return;
}


typedef void (*Subtractor_T) (UINT4 *, UINT4);

static Subtractor_T subtractor_table[272] =
  {sub_00, sub_00, sub_00, sub_00,
   sub_00, sub_00, sub_00, sub_00,
   sub_00, sub_00, sub_00, sub_00,
   sub_00, sub_00, sub_00, sub_00,

   sub_02_00, sub_02_01, sub_02_02, sub_02_03,
   sub_02_04, sub_02_05, sub_02_06, sub_02_07,
   sub_02_08, sub_02_09, sub_02_10, sub_02_11,
   sub_02_12, sub_02_13, sub_02_14, sub_02_15,

   sub_04_00, sub_04_01, sub_04_02, sub_04_03,
   sub_04_04, sub_04_05, sub_04_06, sub_04_07,
   sub_04_08, sub_04_09, sub_04_10, sub_04_11,
   sub_04_12, sub_04_13, sub_04_14, sub_04_15,

   sub_06_00, sub_06_01, sub_06_02, sub_06_03,
   sub_06_04, sub_06_05, sub_06_06, sub_06_07,
   sub_06_08, sub_06_09, sub_06_10, sub_06_11,
   sub_06_12, sub_06_13, sub_06_14, sub_06_15,

   sub_08_00, sub_08_01, sub_08_02, sub_08_03,
   sub_08_04, sub_08_05, sub_08_06, sub_08_07,
   sub_08_08, sub_08_09, sub_08_10, sub_08_11,
   sub_08_12, sub_08_13, sub_08_14, sub_08_15,

   sub_10_00, sub_10_01, sub_10_02, sub_10_03,
   sub_10_04, sub_10_05, sub_10_06, sub_10_07,
   sub_10_08, sub_10_09, sub_10_10, sub_10_11,
   sub_10_12, sub_10_13, sub_10_14, sub_10_15,

   sub_12_00, sub_12_01, sub_12_02, sub_12_03,
   sub_12_04, sub_12_05, sub_12_06, sub_12_07,
   sub_12_08, sub_12_09, sub_12_10, sub_12_11,
   sub_12_12, sub_12_13, sub_12_14, sub_12_15,

   sub_14_00, sub_14_01, sub_14_02, sub_14_03,
   sub_14_04, sub_14_05, sub_14_06, sub_14_07,
   sub_14_08, sub_14_09, sub_14_10, sub_14_11,
   sub_14_12, sub_14_13, sub_14_14, sub_14_15,

   sub_16_00, sub_16_01, sub_16_02, sub_16_03,
   sub_16_04, sub_16_05, sub_16_06, sub_16_07,
   sub_16_08, sub_16_09, sub_16_10, sub_16_11,
   sub_16_12, sub_16_13, sub_16_14, sub_16_15,

   sub_18_00, sub_18_01, sub_18_02, sub_18_03,
   sub_18_04, sub_18_05, sub_18_06, sub_18_07,
   sub_18_08, sub_18_09, sub_18_10, sub_18_11,
   sub_18_12, sub_18_13, sub_18_14, sub_18_15,

   sub_20_00, sub_20_01, sub_20_02, sub_20_03,
   sub_20_04, sub_20_05, sub_20_06, sub_20_07,
   sub_20_08, sub_20_09, sub_20_10, sub_20_11,
   sub_20_12, sub_20_13, sub_20_14, sub_20_15,

   sub_22_00, sub_22_01, sub_22_02, sub_22_03,
   sub_22_04, sub_22_05, sub_22_06, sub_22_07,
   sub_22_08, sub_22_09, sub_22_10, sub_22_11,
   sub_22_12, sub_22_13, sub_22_14, sub_22_15,

   sub_24_00, sub_24_01, sub_24_02, sub_24_03,
   sub_24_04, sub_24_05, sub_24_06, sub_24_07,
   sub_24_08, sub_24_09, sub_24_10, sub_24_11,
   sub_24_12, sub_24_13, sub_24_14, sub_24_15,

   sub_26_00, sub_26_01, sub_26_02, sub_26_03,
   sub_26_04, sub_26_05, sub_26_06, sub_26_07,
   sub_26_08, sub_26_09, sub_26_10, sub_26_11,
   sub_26_12, sub_26_13, sub_26_14, sub_26_15,

   sub_28_00, sub_28_01, sub_28_02, sub_28_03,
   sub_28_04, sub_28_05, sub_28_06, sub_28_07,
   sub_28_08, sub_28_09, sub_28_10, sub_28_11,
   sub_28_12, sub_28_13, sub_28_14, sub_28_15,

   sub_30_00, sub_30_01, sub_30_02, sub_30_03,
   sub_30_04, sub_30_05, sub_30_06, sub_30_07,
   sub_30_08, sub_30_09, sub_30_10, sub_30_11,
   sub_30_12, sub_30_13, sub_30_14, sub_30_15,

   sub_32_00, sub_32_01, sub_32_02, sub_32_03,
   sub_32_04, sub_32_05, sub_32_06, sub_32_07,
   sub_32_08, sub_32_09, sub_32_10, sub_32_11,
   sub_32_12, sub_32_13, sub_32_14, sub_32_15,
  };
#endif





#define DIRECT_METAINFO_SIZE 1

/* Vertical */
UINT4
Bitpack64_incr (Oligospace_T oligo, UINT4 *ptrs, UINT4 *comp) {
  UINT4 *info, start;
  int nwritten, remainder;
  UINT4 *bitpack;
  int index, column;
#ifdef CHECK
  UINT4 prev, value;
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

#ifdef CHECK
  prev = Bitpack64_access(oligo,ptrs,comp);
  value = (incrementor_table[index])(&(bitpack[column]));
  if (value != prev || Bitpack64_access(oligo,ptrs,comp) != prev + 1) {
    packsize = nwritten*2;
    printf("Error at Bitpack64_incr with oligo %llu, packsize %d, remainder %d, column %d, index %d\n",
	   oligo,packsize,remainder,column,index);
    printf("bitpack:\n");
    for (i = 0; i < nwritten*4; i += 4) {
      printf("%08X %08X %08X %08X\n",bitpack[i],bitpack[i+1],bitpack[i+2],bitpack[i+3]);
    }
    printf("\n");
    abort();
  }
  return value;
#else
  return (incrementor_table[index])(&(bitpack[column]));
#endif

}


UINT4
Bitpack64_incr_bitpack (Oligospace_T oligo, int packsize, UINT4 *bitpack) {
  int nregisters, remainder;
  int index, column;
#ifdef CHECK
  UINT4 prev, value;
#endif

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*16 + remainder/4;
  column = remainder % 4;

#ifdef CHECK
  prev = Bitpack64_access_bitpack(oligo,packsize,bitpack);
  value = (incrementor_table[index])(&(bitpack[column]));
  if (value != prev || Bitpack64_access_bitpack(oligo,packsize,bitpack) != prev + 1) {
    printf("Error at Bitpack64_incr_bitpack with oligo %llu, packsize %d, remainder %d, column %d, index %d\n",
	   oligo,packsize,remainder,column,index);
    abort();
  }
  return value;
#else
  return (incrementor_table[index])(&(bitpack[column]));
#endif
}


void
Bitpack64_add (Oligospace_T oligo, UINT4 *ptrs, UINT4 *comp, UINT4 increment) {
  UINT4 *info, start;
  int nwritten, remainder;
  UINT4 *bitpack;
  int index, column;
#ifdef CHECK
  UINT4 prev;
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

#ifdef CHECK
  prev = Bitpack64_access(oligo,ptrs,comp);
  (adder_table[index])(&(bitpack[column]),increment);
  if (Bitpack64_access(oligo,ptrs,comp) != prev + increment) {
    packsize = nwritten*2;
    printf("Error at Bitpack64_add with oligo %llu, packsize %d, remainder %d, column %d, index %d\n",
	   oligo,packsize,remainder,column,index);
    printf("bitpack:\n");
    for (i = 0; i < nwritten*4; i += 4) {
      printf("%08X %08X %08X %08X\n",bitpack[i],bitpack[i+1],bitpack[i+2],bitpack[i+3]);
    }
    printf("\n");
    abort();
  }
#else
  (adder_table[index])(&(bitpack[column]),increment);
#endif

  return;
}



/* Assumes that increment > 0 */
void
Bitpack64_add_bitpack (Oligospace_T oligo, int packsize, UINT4 *bitpack, UINT4 increment) {
  int nregisters, remainder;
  int index, column;
#ifdef CHECK
  UINT4 prev;
#endif

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*16 + remainder/4;
  column = remainder % 4;

#ifdef CHECK
  prev = Bitpack64_access_bitpack(oligo,packsize,bitpack);
  (adder_table[index])(&(bitpack[column]),increment);
  if (Bitpack64_access_bitpack(oligo,packsize,bitpack) != prev + increment) {
    printf("Error at Bitpack64_add_bitpack with oligo %llu, packsize %d, remainder %d, column %d, index %d\n",
	   oligo,packsize,remainder,column,index);
    abort();
  }
#else
  (adder_table[index])(&(bitpack[column]),increment);
#endif

  return;
}




#if 0
/* No need currently for a subtraction routine */

void
Bitpack64_sub_bitpack (Oligospace_T oligo, int packsize, UINT4 *bitpack, UINT4 decrement) {
  int nregisters, remainder;
  int index, column;
#ifdef CHECK
  UINT4 prev;
#endif

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*16 + remainder/4;
  column = remainder % 4;

#ifdef CHECK
  prev = Bitpack64_access_bitpack(oligo,packsize,bitpack);
  (subtractor_table[index])(&(bitpack[column]),decrement);
  if (Bitpack64_access_bitpack(oligo,packsize,bitpack) != prev - decrement) {
    printf("Error at Bitpack64_sub_bitpack with oligo %llu, packsize %d, remainder %d, column %d, index %d\n",
	   oligo,packsize,remainder,column,index);
    abort();
  }
#else
  (subtractor_table[index])(&(bitpack[column]),decrement);
#endif

  return;
}
#endif


static void
transfer_00_02 (UINT4 *out, const UINT4 *in) {
  return;
}

static void
transfer_02_04 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 2 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 2 ) ;
  *out |= (value << 4);

  /* 02 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 2 ) ;
  *out |= (value << 8);

  /* 03 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 2 ) ;
  *out |= (value << 12);

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 2 ) ;
  *out |= (value << 16);

  /* 05 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 2 ) ;
  *out |= (value << 20);

  /* 06 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 2 ) ;
  *out |= (value << 24);

  /* 07 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 2 ) ;
  *out |= (value << 28);

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 2 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 09 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 2 ) ;
  *out |= (value << 4);

  /* 10 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 2 ) ;
  *out |= (value << 8);

  /* 11 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 2 ) ;
  *out |= (value << 12);

  /* 12 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 2 ) ;
  *out |= (value << 16);

  /* 13 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 2 ) ;
  *out |= (value << 20);

  /* 14 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 2 ) ;
  *out |= (value << 24);

  /* 15 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 2 ) ;
  *out |= (value << 28);

  return;
}

static void
transfer_04_06 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 4 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 4 ) ;
  *out |= (value << 6);

  /* 02 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 4 ) ;
  *out |= (value << 12);

  /* 03 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 4 ) ;
  *out |= (value << 18);

  /* 04 */
  value = ( CONVERT(*in) >>  16 )   % (1U << 4 ) ;
  *out |= (value << 24);

  /* 05 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 4 ) ;
  *out |= (value << 30);
  out += WORD_INCR;
  *out |= (value >> (6 - 4));

  /* 06 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 4 ) ;
  *out |= (value << 4);

  /* 07 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 4 ) ;
  *out |= (value << 10);

  /* 08 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0  )   % (1U << 4 ) ;
  *out |= (value << 16);

  /* 09 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 4 ) ;
  *out |= (value << 22);

  /* 10 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 4 ) ;
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (6 - 2));

  /* 11 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 4 ) ;
  *out |= (value << 2);

  /* 12 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 4 ) ;
  *out |= (value << 8);

  /* 13 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 4 ) ;
  *out |= (value << 14);

  /* 14 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 4 ) ;
  *out |= (value << 20);

  /* 15 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 4 ) ;
  *out |= (value << 26);

  return;
}


static void
transfer_06_08 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 6 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 6 ) ;
  *out |= (value << 8);

  /* 02 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;
  *out |= (value << 16);

  /* 03 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 6 ) ;
  *out |= (value << 24);

  /* 04 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 6 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 05 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 6 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 );
  *out |= (value << 8);

  /* 06 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 6 ) ;
  *out |= (value << 16);

  /* 07 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 6 ) ;
  *out |= (value << 24);

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 6 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 09 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 6 ) ;
  *out |= (value << 8);

  /* 10 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 6 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 );
  *out |= (value << 16);

  /* 11 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 6 ) ;
  *out |= (value << 24);

  /* 12 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 6 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 13 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;
  *out |= (value << 8);

  /* 14 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 6 ) ;
  *out |= (value << 16);

  /* 15 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 6 ) ;
  *out |= (value << 24);

  return;
}


static void
transfer_08_10 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *out |= (value << 10);

  /* 02 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
  *out |= (value << 20);

  /* 03 */
  value = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
  *out |= (value << 30);
  out += WORD_INCR;
  *out |= (value >> (10 - 8));

  /* 04 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *out |= (value << 8);

  /* 05 */
  value = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *out |= (value << 18);

  /* 06 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (10 - 6));

  /* 07 */
  value = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
  *out |= (value << 6);

  /* 08 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *out |= (value << 16);

  /* 09 */
  value = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *out |= (value << 26);
  out += WORD_INCR;
  *out |= (value >> (10 - 4));

  /* 10 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
  *out |= (value << 4);

  /* 11 */
  value = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
  *out |= (value << 14);

  /* 12 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (10 - 2));

  /* 13 */
  value = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *out |= (value << 2);

  /* 14 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 8 ) ;
  *out |= (value << 12);

  /* 15 */
  value = ( CONVERT(*in) >> 24 )   % (1U << 8 ) ;
  *out |= (value << 22);

  return;
}

static void
transfer_10_12 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 10 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;
  *out |= (value << 12);

  /* 02 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 10 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (12 - 4));

  /* 03 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 );
  *out |= (value << 4);

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;
  *out |= (value << 16);

  /* 05 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 10 ) ;
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (12 - 8));

  /* 06 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 );
  *out |= (value << 8);

  /* 07 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 10 ) ;
  *out |= (value << 20);

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 10 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 09 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 );
  *out |= (value << 12);

  /* 10 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 10 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (12 - 4));

  /* 11 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;
  *out |= (value << 4);

  /* 12 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 );
  *out |= (value << 16);

  /* 13 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 10 ) ;
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (12 - 8));

  /* 14 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;
  *out |= (value << 8);

  /* 15 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 10 ) ;
  *out |= (value << 20);

  return;
}

static void
transfer_12_14 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  *out |= (value << 14);

  /* 02 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (14 - 10));

  /* 03 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
  *out |= (value << 10);

  /* 04 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 12 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (14 - 6));

  /* 05 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  *out |= (value << 6);

  /* 06 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (14 - 2));

  /* 07 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 12 ) ;
  *out |= (value << 2);

  /* 08 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
  *out |= (value << 16);

  /* 09 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  *out |= (value << 30);
  out += WORD_INCR;
  *out |= (value >> (14 - 12));

  /* 10 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  *out |= (value << 12);

  /* 11 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
  *out |= (value << 26);
  out += WORD_INCR;
  *out |= (value >> (14 - 8));

  /* 12 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 12 ) ;
  *out |= (value << 8);

  /* 13 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  *out |= (value << 22);
  out += WORD_INCR;
  *out |= (value >> (14 - 4));

  /* 14 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  *out |= (value << 4);

  /* 15 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 12 ) ;
  *out |= (value << 18);

  return;
}


static void
transfer_14_16 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 14 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;
  *out |= (value << 16);

  /* 02 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 );
  out += WORD_INCR;
  *out |= (value << 0);

  /* 03 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;
  *out |= (value << 16);

  /* 04 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 );
  out += WORD_INCR;
  *out |= (value << 0);

  /* 05 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;
  *out |= (value << 16);

  /* 06 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 );
  out += WORD_INCR;
  *out |= (value << 0);

  /* 07 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 14 ) ;
  *out |= (value << 16);

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 14 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 09 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 );
  *out |= (value << 16);

  /* 10 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 11 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 );
  *out |= (value << 16);

  /* 12 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 13 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 );
  *out |= (value << 16);

  /* 14 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 15 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 14 ) ;
  *out |= (value << 16);

  return;
}

static void
transfer_16_18 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *out |= (value << 18);
  out += WORD_INCR;
  *out |= (value >> (18 - 4));

  /* 02 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *out |= (value << 4);

  /* 03 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *out |= (value << 22);
  out += WORD_INCR;
  *out |= (value >> (18 - 8));

  /* 04 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *out |= (value << 8);

  /* 05 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *out |= (value << 26);
  out += WORD_INCR;
  *out |= (value >> (18 - 12));

  /* 06 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *out |= (value << 12);

  /* 07 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *out |= (value << 30);
  out += WORD_INCR;
  *out |= (value >> (18 - 16));

  /* 08 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (18 - 2));

  /* 09 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *out |= (value << 2);

  /* 10 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (18 - 6));

  /* 11 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *out |= (value << 6);

  /* 12 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (18 - 10));

  /* 13 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *out |= (value << 10);

  /* 14 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 16 ) ;
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (18 - 14));

  /* 15 */
  value = ( CONVERT(*in) >> 16 )   % (1U << 16 ) ;
  *out |= (value << 14);

  return;
}


static void
transfer_18_20 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 18 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 18 - 4 );
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (20 - 8));

  /* 02 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 18 ) ;
  *out |= (value << 8);

  /* 03 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 18 - 8 );
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (20 - 16));

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 18 ) ;
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (20 - 4));

  /* 05 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 18 - 12 );
  *out |= (value << 4);

  /* 06 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 18 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (20 - 12));

  /* 07 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 18 - 16 );
  *out |= (value << 12);

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 18 - 2 );
  out += WORD_INCR;
  *out |= (value << 0);

  /* 09 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 18 ) ;
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (20 - 8));

  /* 10 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 18 - 6 );
  *out |= (value << 8);

  /* 11 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 18 ) ;
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (20 - 16));

  /* 12 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 18 - 10 );
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (20 - 4));

  /* 13 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 18 ) ;
  *out |= (value << 4);

  /* 14 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 18 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 14 ))<<( 18 - 14 );
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (20 - 12));

  /* 15 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 18 ) ;
  *out |= (value << 12);

  return;
}

static void
transfer_20_22 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 20 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 );
  *out |= (value << 22);
  out += WORD_INCR;
  *out |= (value >> (22 - 12));

  /* 02 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 20 ) ;
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (22 - 2));

  /* 03 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 );
  *out |= (value << 2);

  /* 04 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 );
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (22 - 14));

  /* 05 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 20 ) ;
  *out |= (value << 14);
  out += WORD_INCR;
  *out |= (value >> (22 - 4));

  /* 06 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 );
  *out |= (value << 4);

  /* 07 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 20 ) ;
  *out |= (value << 26);
  out += WORD_INCR;
  *out |= (value >> (22 - 16));

  /* 08 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0  )   % (1U << 20 ) ;
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (22 - 6));

  /* 09 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 20 - 8 );
  *out |= (value << 6);

  /* 10 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 20 ) ;
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (22 - 18));

  /* 11 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 20 - 16 );
  *out |= (value << 18);
  out += WORD_INCR;
  *out |= (value >> (22 - 8));

  /* 12 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 20 - 4 );
  *out |= (value << 8);

  /* 13 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 20 ) ;
  *out |= (value << 30);
  out += WORD_INCR;
  *out |= (value >> (22 - 20));

  /* 14 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 20 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 20 - 12 );
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (22 - 10));

  /* 15 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 20 ) ;
  *out |= (value << 10);

  return;
}

static void
transfer_22_24 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 22 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 22 - 12 );
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (24 - 16));

  /* 02 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 22 - 2 );
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (24 - 8));

  /* 03 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 22 ) ;
  *out |= (value << 8);

  /* 04 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 14 ))<<( 22 - 14 );
  out += WORD_INCR;
  *out |= (value << 0);

  /* 05 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 22 - 4 );
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (24 - 16));

  /* 06 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 22 ) ;
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (24 - 8));

  /* 07 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 22 - 16 );
  *out |= (value << 8);

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 22 - 6 );
  out += WORD_INCR;
  *out |= (value << 0);

  /* 09 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 22 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (24 - 16));

  /* 10 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 18 ))<<( 22 - 18 );
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (24 - 8));

  /* 11 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 22 - 8 );
  *out |= (value << 8);	/* was 16 */

  /* 12 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 22 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 13 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 22 - 20 );
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (24 - 16));

  /* 14 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 22 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 22 - 10 );
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (24 - 8));

  /* 15 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 22 ) ;
  *out |= (value << 8);

  return;
}

static void
transfer_24_26 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  *out |= (value << 26);
  out += WORD_INCR;
  *out |= (value >> (26 - 20));

  /* 02 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (26 - 14));

  /* 03 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
  *out |= (value << 14);
  out += WORD_INCR;
  *out |= (value >> (26 - 8));

  /* 04 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
  *out |= (value << 8);
  out += WORD_INCR;
  *out |= (value >> (26 - 2));

  /* 05 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  *out |= (value << 2);

  /* 06 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (26 - 22));

  /* 07 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
  *out |= (value << 22);
  out += WORD_INCR;
  *out |= (value >> (26 - 16));

  /* 08 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (26 - 10));

  /* 09 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  *out |= (value << 10);
  out += WORD_INCR;
  *out |= (value >> (26 - 4));

  /* 10 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  *out |= (value << 4);

  /* 11 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
  *out |= (value << 30);
  out += WORD_INCR;
  *out |= (value >> (26 - 24));

  /* 12 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0  )   % (1U << 24 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (26 - 18));

  /* 13 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 24 - 16 );
  *out |= (value << 18);
  out += WORD_INCR;
  *out |= (value >> (26 - 12));

  /* 14 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 24 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 24 - 8 );
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (26 - 6));

  /* 15 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 24 ) ;
  *out |= (value << 6);

  return;
}

static void
transfer_26_28 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 26 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 26 - 20 );
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (28 - 24));

  /* 02 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 14 ))<<( 26 - 14 );
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (28 - 20));

  /* 03 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 26 - 8 );
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (28 - 16));

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 26 - 2 );
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (28 - 12));	/* was (28 - 8) */

  /* 05 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 26 ) ;
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (28 - 8));

  /* 06 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 22 ))<<( 26 - 22 );
  *out |= (value << 8);
  out += WORD_INCR;
  *out |= (value >> (28 - 4));

  /* 07 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 26 - 16 );
  *out |= (value << 4);

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 26 - 10 );
  out += WORD_INCR;
  *out |= (value << 0);

  /* 09 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 26 - 4 );
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (28 - 24));

  /* 10 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 26 ) ;
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (28 - 20));

  /* 11 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 24 ))<<( 26 - 24 );
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (28 - 16));

  /* 12 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 18 ))<<( 26 - 18 );
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (28 - 12));

  /* 13 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 26 - 12 );
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (28 - 8));

  /* 14 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 26 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 26 - 6 );
  *out |= (value << 8);
  out += WORD_INCR;
  *out |= (value >> (28 - 4));

  /* 15 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 26 ) ;
  *out |= (value << 4);

  return;
}

static void
transfer_28_30 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 28 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 );
  *out |= (value << 30);
  out += WORD_INCR;
  *out |= (value >> (30 - 28));

  /* 02 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 );
  *out |= (value << 28);
  out += WORD_INCR;
  *out |= (value >> (30 - 26));

  /* 03 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 );
  *out |= (value << 26);
  out += WORD_INCR;
  *out |= (value >> (30 - 24));

  /* 04 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 );
  *out |= (value << 24);
  out += WORD_INCR;
  *out |= (value >> (30 - 22));

  /* 05 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 );
  *out |= (value << 22);
  out += WORD_INCR;
  *out |= (value >> (30 - 20));

  /* 06 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 );
  *out |= (value << 20);
  out += WORD_INCR;
  *out |= (value >> (30 - 18));

  /* 07 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 28 ) ;
  *out |= (value << 18);
  out += WORD_INCR;
  *out |= (value >> (30 - 16));

  /* 08 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0  )   % (1U << 28 ) ;
  *out |= (value << 16);
  out += WORD_INCR;
  *out |= (value >> (30 - 14));

  /* 09 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 24 ))<<( 28 - 24 );
  *out |= (value << 14);
  out += WORD_INCR;
  *out |= (value >> (30 - 12));

  /* 10 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 28 - 20 );
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (30 - 10));

  /* 11 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 28 - 16 );
  *out |= (value << 10);
  out += WORD_INCR;
  *out |= (value >> (30 - 8));

  /* 12 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 28 - 12 );
  *out |= (value << 8);
  out += WORD_INCR;
  *out |= (value >> (30 - 6));

  /* 13 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 28 - 8 );
  *out |= (value << 6);
  out += WORD_INCR;
  *out |= (value >> (30 - 4));

  /* 14 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 28 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 28 - 4 );
  *out |= (value << 4);
  out += WORD_INCR;
  *out |= (value >> (30 - 2));

  /* 15 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 28 ) ;
  *out |= (value << 2);

  return;
}


static void
transfer_30_32 (UINT4 *out, const UINT4 *in) {
  UINT4 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 30 ) ;
  *out |= value;

  /* 01 */
  value = ( CONVERT(*in) >>  30  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 28 ))<<( 30 - 28 );
  out += WORD_INCR;
  *out |= value;

  /* 02 */
  value = ( CONVERT(*in) >>  28  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 26 ))<<( 30 - 26 );
  out += WORD_INCR;
  *out |= value;

  /* 03 */
  value = ( CONVERT(*in) >>  26  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 24 ))<<( 30 - 24 );
  out += WORD_INCR;
  *out |= value;

  /* 04 */
  value = ( CONVERT(*in) >>  24  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 22 ))<<( 30 - 22 );
  out += WORD_INCR;
  *out |= value;

  /* 05 */
  value = ( CONVERT(*in) >>  22  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 20 ))<<( 30 - 20 );
  out += WORD_INCR;
  *out |= value;

  /* 06 */
  value = ( CONVERT(*in) >>  20  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 18 ))<<( 30 - 18 );
  out += WORD_INCR;
  *out |= value;

  /* 07 */
  value = ( CONVERT(*in) >>  18  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 16 ))<<( 30 - 16 );
  out += WORD_INCR;
  *out |= value;

  /* 08 */
  value = ( CONVERT(*in) >>  16  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 14 ))<<( 30 - 14 );
  out += WORD_INCR;
  *out |= value;

  /* 09 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 30 - 12 );
  out += WORD_INCR;
  *out |= value;

  /* 10 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 30 - 10 );
  out += WORD_INCR;
  *out |= value;

  /* 11 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 30 - 8 );
  out += WORD_INCR;
  *out |= value;

  /* 12 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 30 - 6 );
  out += WORD_INCR;
  *out |= value;

  /* 13 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 30 - 4 );
  out += WORD_INCR;
  *out |= value;

  /* 14 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 30 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 30 - 2 );
  out += WORD_INCR;
  *out |= value;

  /* 15 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 30 ) ;
  out += WORD_INCR;
  *out |= value;

  return;
}



typedef void (*Transfer_T) (UINT4 *, const UINT4 *);

static Transfer_T transfer_table[16] =
  {transfer_00_02, transfer_02_04, transfer_04_06, transfer_06_08,
   transfer_08_10, transfer_10_12, transfer_12_14, transfer_14_16,
   transfer_16_18, transfer_18_20, transfer_20_22, transfer_22_24,
   transfer_24_26, transfer_26_28, transfer_28_30, transfer_30_32};

UINT4 *
Bitpack64_realloc_one (int packsize, UINT4 *bitpack) {
  UINT4 *new;
  int nregisters;
#ifdef CHECK
  int i;
#endif

  new = (UINT4 *) CALLOC((packsize + 2) / 2 * 4,sizeof(UINT4));

  nregisters = packsize/2;
  (transfer_table[nregisters])(&(new[/*column*/0]),&(bitpack[/*column*/0]));
  (transfer_table[nregisters])(&(new[/*column*/1]),&(bitpack[/*column*/1]));
  (transfer_table[nregisters])(&(new[/*column*/2]),&(bitpack[/*column*/2]));
  (transfer_table[nregisters])(&(new[/*column*/3]),&(bitpack[/*column*/3]));

#ifdef CHECK
  for (i = 0; i < 64; i++) {
    if (Bitpack64_access_bitpack(i,packsize+2,new) != Bitpack64_access_bitpack(i,packsize,bitpack)) {
      fprintf(stderr,"Difference in packsize %d -> %d, index %d: %u != %u\n",
	      packsize,packsize+2,i,
	      Bitpack64_access_bitpack(i,packsize,bitpack),
	      Bitpack64_access_bitpack(i,packsize+2,new));
      abort();
    }
  }
#endif

  FREE(bitpack);
  return new;
}


UINT4 *
Bitpack64_realloc_multiple (int old_packsize, int new_packsize, UINT4 *bitpack) {
  UINT4 *new;
  int nregisters, packsize;
#ifdef CHECK
  int i;
#endif

  for (packsize = old_packsize; packsize < new_packsize; packsize += 2) {
    new = (UINT4 *) CALLOC((packsize + 2) / 2 * 4,sizeof(UINT4));

    nregisters = packsize/2;
    (transfer_table[nregisters])(&(new[/*column*/0]),&(bitpack[/*column*/0]));
    (transfer_table[nregisters])(&(new[/*column*/1]),&(bitpack[/*column*/1]));
    (transfer_table[nregisters])(&(new[/*column*/2]),&(bitpack[/*column*/2]));
    (transfer_table[nregisters])(&(new[/*column*/3]),&(bitpack[/*column*/3]));

#ifdef CHECK
    for (i = 0; i < 64; i++) {
      if (Bitpack64_access_bitpack(i,packsize+2,new) != Bitpack64_access_bitpack(i,packsize,bitpack)) {
	fprintf(stderr,"Difference in packsize %d -> %d, index %d: %u != %u\n",
		packsize,packsize+2,i,
		Bitpack64_access_bitpack(i,packsize,bitpack),
		Bitpack64_access_bitpack(i,packsize+2,new));
	abort();
      }
    }
#endif

    FREE(bitpack);
    bitpack = new;
  }

  return bitpack;
}

