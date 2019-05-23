static char rcsid[] = "$Id: epu16-bitpack64-incr.c 218163 2019-01-17 06:08:15Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "epu16-bitpack64-incr.h"
#include "epu16-bitpack64-access.h"

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#define CONVERT(x) Bigendian_convert_uint(x)
#else
#define CONVERT(x) x
#endif


/* #define CHECK 1 */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define BLOCKSIZE 64

/* Choice needs to match that in epu16-bitpack64-access */
#ifdef HORIZONTAL
#define WORD_INCR 1		/* 1 for horizontal; 8 for vertical */
#else
#define WORD_INCR 8
#endif

#define BIT0  0x1;
#define BIT2  0x4;
#define BIT4  0x10;
#define BIT6  0x40;
#define BIT8  0x100;
#define BIT10 0x400;
#define BIT12 0x1000;
#define BIT14 0x4000;

static UINT2
incr_00 (UINT2 *in) {
  abort();
  return 0U;
}

static UINT2
incr_02_00 (UINT2 *in) {
  UINT2 prev = ( CONVERT(*in) >> 0 )   % (1U << 2 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_02_01 (UINT2 *in) {
  UINT2 prev = ( CONVERT(*in) >> 2 )   % (1U << 2 ) ;
  *in += BIT2;
  return prev;
}

static UINT2
incr_02_02 (UINT2 *in) {
  UINT2 prev = ( CONVERT(*in) >> 4 )   % (1U << 2 ) ;
  *in += BIT4;
  return prev;
}

static UINT2
incr_02_03 (UINT2 *in) {
  UINT2 prev = ( CONVERT(*in) >> 6 )   % (1U << 2 ) ;
  *in += BIT6;
  return prev;
}

static UINT2
incr_02_04 (UINT2 *in) {
  UINT2 prev = ( CONVERT(*in) >> 8 )   % (1U << 2 ) ;
  *in += BIT8;
  return prev;
}

static UINT2
incr_02_05 (UINT2 *in) {
  UINT2 prev = ( CONVERT(*in) >> 10 )   % (1U << 2 ) ;
  *in += BIT10;
  return prev;
}

static UINT2
incr_02_06 (UINT2 *in) {
  UINT2 prev = ( CONVERT(*in) >> 12 )   % (1U << 2 ) ;
  *in += BIT12;
  return prev;
}

static UINT2
incr_02_07 (UINT2 *in) {
  UINT2 prev = ( CONVERT(*in) >> 14 )   % (1U << 2 ) ;
  *in += BIT14;
  return prev;
}



static UINT2
incr_04_00 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_04_01 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;
  *in += BIT4;
  return prev;
}

static UINT2
incr_04_02 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;
  *in += BIT8;
  return prev;
}

static UINT2
incr_04_03 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;
  *in += BIT12;
  return prev;
}

static UINT2
incr_04_04 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_04_05 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;
  *in += BIT4;
  return prev;
}

static UINT2
incr_04_06 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;
  *in += BIT8;
  return prev;
}

static UINT2
incr_04_07 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;
  *in += BIT12;
  return prev;
}


static UINT2
incr_06_00 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 6 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_06_01 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >>  6  )   % (1U << 6 ) ;
  *in += BIT6;
  return prev;
}

static UINT2
incr_06_02 (UINT2 *in) {
  UINT2 prev, carry;

  prev = ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;
  carry = ( (prev + 1) >> (6 - 2) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 );
  *in += carry;

  return prev;
}

static UINT2
incr_06_03 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 6 ) ;
  *in += BIT2;
  return prev;
}

static UINT2
incr_06_04 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 6 ) ;
  *in += BIT8;
  return prev;
}

static UINT2
incr_06_05 (UINT2 *in) {
  UINT2 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;
  carry = ( (prev + 1) >> (6 - 4));
  *in += BIT14;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 );
  *in += carry;

  return prev;
}

static UINT2
incr_06_06 (UINT2 *in) {
  UINT2 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 6 ) ;
  *in += BIT4;
  return prev;
}

static UINT2
incr_06_07 (UINT2 *in) {
  UINT2 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 6 ) ;
  *in += BIT10;
  return prev;
}


static UINT2
incr_08_00 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_08_01 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *in += BIT8;
  return prev;
}

static UINT2
incr_08_02 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_08_03 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *in += BIT8;
  return prev;
}

static UINT2
incr_08_04 (UINT2 *in) {
  UINT2 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_08_05 (UINT2 *in) {
  UINT2 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *in += BIT8;
  return prev;
}

static UINT2
incr_08_06 (UINT2 *in) {
  UINT2 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_08_07 (UINT2 *in) {
  UINT2 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *in += BIT8;
  return prev;
}



static UINT2
incr_10_00 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 10 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_10_01 (UINT2 *in) {
  UINT2 prev, carry;

  prev = ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;
  carry = ( (prev + 1) >> (10 - 4) );
  *in += BIT10;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 );
  *in += carry;

  return prev;
}

static UINT2
incr_10_02 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 10 ) ;
  *in += BIT4;
  return prev;
}

static UINT2
incr_10_03 (UINT2 *in) {
  UINT2 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;
  carry = ( (prev + 1) >> (10 - 8) );
  *in += BIT14;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 );
  *in += carry;

  return prev;
}

static UINT2
incr_10_04 (UINT2 *in) {
  UINT2 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;
  carry = ( (prev + 1) >> (10 - 2) );
  *in += BIT8;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 );
  *in += carry;

  return prev;
}

static UINT2
incr_10_05 (UINT2 *in) {
  UINT2 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 10 ) ;
  *in += BIT2;
  return prev;
}

static UINT2
incr_10_06 (UINT2 *in) {
  UINT2 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;
  carry = ( (prev + 1) >> (10 - 6) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 );
  *in += carry;

  return prev;
}

static UINT2
incr_10_07 (UINT2 *in) {
  UINT2 prev;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 10 ) ;
  *in += BIT6;
  return prev;
}


static UINT2
incr_12_00 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_12_01 (UINT2 *in) {
  UINT2 prev, carry;

  prev = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  carry = ( (prev + 1) >> (12 - 8) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  *in += carry;

  return prev;
}

static UINT2
incr_12_02 (UINT2 *in) {
  UINT2 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  carry = ( (prev + 1) >> (12 - 4) );
  *in += BIT8;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  *in += carry;

  return prev;
}

static UINT2
incr_12_03 (UINT2 *in) {
  UINT2 prev;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
  *in += BIT4;
  return prev;
}

static UINT2
incr_12_04 (UINT2 *in) {
  UINT2 prev;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_12_05 (UINT2 *in) {
  UINT2 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  carry = ( (prev + 1) >> (12 - 8) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  *in += carry;

  return prev;
}

static UINT2
incr_12_06 (UINT2 *in) {
  UINT2 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  carry = ( (prev + 1) >> (12 - 4) );
  *in += BIT8;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  *in += carry;

  return prev;
}

static UINT2
incr_12_07 (UINT2 *in) {
  UINT2 prev;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
  *in += BIT4;
  return prev;
}



static UINT2
incr_14_00 (UINT2 *in) {
  UINT2 prev;

  prev = ( CONVERT(*in) >>  0  )   % (1U << 14 ) ;
  *in += BIT0;
  return prev;
}

static UINT2
incr_14_01 (UINT2 *in) {
  UINT2 prev, carry;

  prev = ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 12) );
  *in += BIT14;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 );
  *in += carry;

  return prev;
}

static UINT2
incr_14_02 (UINT2 *in) {
  UINT2 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 10) );
  *in += BIT12;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 );
  *in += carry;

  return prev;
}

static UINT2
incr_14_03 (UINT2 *in) {
  UINT2 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 8) );
  *in += BIT10;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 );
  *in += carry;

  return prev;
}

static UINT2
incr_14_04 (UINT2 *in) {
  UINT2 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 6) );
  *in += BIT8;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 );
  *in += carry;

  return prev;
}

static UINT2
incr_14_05 (UINT2 *in) {
  UINT2 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 4) );
  *in += BIT6;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 );
  *in += carry;

  return prev;
}

static UINT2
incr_14_06 (UINT2 *in) {
  UINT2 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;
  carry = ( (prev + 1) >> (14 - 2) );
  *in += BIT4;

  in += 1 * WORD_INCR;
  prev |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 );
  *in += carry;

  return prev;
}

static UINT2
incr_14_07 (UINT2 *in) {
  UINT2 prev;

  in += 6 * WORD_INCR;
  prev = ( CONVERT(*in) >>  2  )   % (1U << 14 ) ;
  *in += BIT2;
  return prev;
}



static UINT2
incr_16_00 (UINT2 *in) {
  UINT2 prev;

  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT2
incr_16_01 (UINT2 *in) {
  UINT2 prev;

  in += 1 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT2
incr_16_02 (UINT2 *in) {
  UINT2 prev;

  in += 2 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT2
incr_16_03 (UINT2 *in) {
  UINT2 prev;

  in += 3 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT2
incr_16_04 (UINT2 *in) {
  UINT2 prev;

  in += 4 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT2
incr_16_05 (UINT2 *in) {
  UINT2 prev;

  in += 5 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}

static UINT2
incr_16_06 (UINT2 *in) {
  UINT2 prev;

  in += 6 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}


static UINT2
incr_16_07 (UINT2 *in) {
  UINT2 prev;

  in += 7 * WORD_INCR;
  prev = CONVERT(*in);
  *in += 1;
  return prev;
}


typedef UINT2 (*Incrementor_T) (UINT2 *);

static Incrementor_T incrementor_table[72] =
  {incr_00, incr_00, incr_00, incr_00,
   incr_00, incr_00, incr_00, incr_00,

   incr_02_00, incr_02_01, incr_02_02, incr_02_03,
   incr_02_04, incr_02_05, incr_02_06, incr_02_07,

   incr_04_00, incr_04_01, incr_04_02, incr_04_03,
   incr_04_04, incr_04_05, incr_04_06, incr_04_07,

   incr_06_00, incr_06_01, incr_06_02, incr_06_03,
   incr_06_04, incr_06_05, incr_06_06, incr_06_07,

   incr_08_00, incr_08_01, incr_08_02, incr_08_03,
   incr_08_04, incr_08_05, incr_08_06, incr_08_07,

   incr_10_00, incr_10_01, incr_10_02, incr_10_03,
   incr_10_04, incr_10_05, incr_10_06, incr_10_07,

   incr_12_00, incr_12_01, incr_12_02, incr_12_03,
   incr_12_04, incr_12_05, incr_12_06, incr_12_07,

   incr_14_00, incr_14_01, incr_14_02, incr_14_03,
   incr_14_04, incr_14_05, incr_14_06, incr_14_07,

   incr_16_00, incr_16_01, incr_16_02, incr_16_03,
   incr_16_04, incr_16_05, incr_16_06, incr_16_07,
  };


UINT2
Epu16_bitpack64_incr (Localspace_T oligo, int packsize, UINT2 *bitpack) {
  int nregisters, remainder;
  int index, column;
#ifdef CHECK
  UINT2 prev, value;
#endif

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*8 + remainder/8;
  column = remainder % 8;

#ifdef CHECK
  prev = Epu16_bitpack64_access(oligo,packsize,bitpack);
  value = (incrementor_table[index])(&(bitpack[column]));
  if (value != prev || Epu16_bitpack64_access(oligo,packsize,bitpack) != prev + 1) {
    printf("Error at Epu16_bitpack64_incr with oligo %llu, packsize %d, remainder %d, column %d, index %d\n",
	   oligo,packsize,remainder,column,index);
    abort();
  }
  return value;
#else
  return (incrementor_table[index])(&(bitpack[column]));
#endif
}



static void
add_00 (UINT2 *in, UINT2 increment) {
  abort();
  return;
}


static void
add_02_00 (UINT2 *in, UINT2 increment) {
  *in += (increment << 0);
  return;
}

static void
add_02_01 (UINT2 *in, UINT2 increment) {
  *in += (increment << 2);
  return;
}

static void
add_02_02 (UINT2 *in, UINT2 increment) {
  *in += (increment << 4);
  return;
}

static void
add_02_03 (UINT2 *in, UINT2 increment) {
  *in += (increment << 6);
  return;
}

static void
add_02_04 (UINT2 *in, UINT2 increment) {
  *in += (increment << 8);
  return;
}

static void
add_02_05 (UINT2 *in, UINT2 increment) {
  *in += (increment << 10);
  return;
}

static void
add_02_06 (UINT2 *in, UINT2 increment) {
  *in += (increment << 12);
  return;
}

static void
add_02_07 (UINT2 *in, UINT2 increment) {
  *in += (increment << 14);
  return;
}


static void
add_04_00 (UINT2 *in, UINT2 increment) {
  *in += (increment << 0);
  return;
}

static void
add_04_01 (UINT2 *in, UINT2 increment) {
  *in += (increment << 4);
  return;
}

static void
add_04_02 (UINT2 *in, UINT2 increment) {
  *in += (increment << 8);
  return;
}

static void
add_04_03 (UINT2 *in, UINT2 increment) {
  *in += (increment << 12);
  return;
}

static void
add_04_04 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_04_05 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_04_06 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_04_07 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 12);
  return;
}


static void
add_06_00 (UINT2 *in, UINT2 increment) {
  *in += (increment << 0);
  return;
}

static void
add_06_01 (UINT2 *in, UINT2 increment) {
  *in += (increment << 6);
  return;
}

static void
add_06_02 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  prev = ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;
  carry = ( (prev + increment) >> (6 - 2) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 ); */
  *in += carry;

  return;
}

static void
add_06_03 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 2);
  return;
}

static void
add_06_04 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_06_05 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;
  carry = ( (prev + increment) >> (6 - 4));
  *in += (increment << 14);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 ); */
  *in += carry;

  return;
}

static void
add_06_06 (UINT2 *in, UINT2 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_06_07 (UINT2 *in, UINT2 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 10);
  return;
}


static void
add_08_00 (UINT2 *in, UINT2 increment) {
  *in += (increment << 0);
  return;
}

static void
add_08_01 (UINT2 *in, UINT2 increment) {
  *in += (increment << 8);
  return;
}

static void
add_08_02 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_08_03 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_08_04 (UINT2 *in, UINT2 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_08_05 (UINT2 *in, UINT2 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 8);
  return;
}

static void
add_08_06 (UINT2 *in, UINT2 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_08_07 (UINT2 *in, UINT2 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 8);
  return;
}


static void
add_10_00 (UINT2 *in, UINT2 increment) {
  *in += (increment << 0);
  return;
}

static void
add_10_01 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  prev = ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;
  carry = ( (prev + increment) >> (10 - 4) );
  *in += (increment << 10);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 ); */
  *in += carry;

  return;
}

static void
add_10_02 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_10_03 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;
  carry = ( (prev + increment) >> (10 - 8) );
  *in += (increment << 14);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 ); */
  *in += carry;

  return;
}

static void
add_10_04 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;
  carry = ( (prev + increment) >> (10 - 2) );
  *in += (increment << 8);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 ); */
  *in += carry;

  return;
}

static void
add_10_05 (UINT2 *in, UINT2 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 2);
  return;
}

static void
add_10_06 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;
  carry = ( (prev + increment) >> (10 - 6) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 ); */
  *in += carry;

  return;
}

static void
add_10_07 (UINT2 *in, UINT2 increment) {
  in += 4 * WORD_INCR;
  *in += (increment << 6);
  return;
}


static void
add_12_00 (UINT2 *in, UINT2 increment) {
  *in += (increment << 0);
  return;
}

static void
add_12_01 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  prev = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  carry = ( (prev + increment) >> (12 - 8) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 ); */
  *in += carry;

  return;
}

static void
add_12_02 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  carry = ( (prev + increment) >> (12 - 4) );
  *in += (increment << 8);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 ); */
  *in += carry;

  return;
}

static void
add_12_03 (UINT2 *in, UINT2 increment) {
  in += 2 * WORD_INCR;
  *in += (increment << 4);
  return;
}

static void
add_12_04 (UINT2 *in, UINT2 increment) {
  in += 3 * WORD_INCR;
  *in += (increment << 0);
  return;
}

static void
add_12_05 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  carry = ( (prev + increment) >> (12 - 8) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 ); */
  *in += carry;

  return;
}

static void
add_12_06 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  carry = ( (prev + increment) >> (12 - 4) );
  *in += (increment << 8);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 ); */
  *in += carry;

  return;
}

static void
add_12_07 (UINT2 *in, UINT2 increment) {
  in += 5 * WORD_INCR;
  *in += (increment << 4);
  return;
}



static void
add_14_00 (UINT2 *in, UINT2 increment) {
  *in += (increment << 0);
  return;
}

static void
add_14_01 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  prev = ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 12) );
  *in += (increment << 14);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 ); */
  *in += carry;

  return;
}

static void
add_14_02 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 1 * WORD_INCR;
  prev = ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 10) );
  *in += (increment << 12);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 ); */
  *in += carry;

  return;
}

static void
add_14_03 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 2 * WORD_INCR;
  prev = ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 8) );
  *in += (increment << 10);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 ); */
  *in += carry;

  return;
}

static void
add_14_04 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 3 * WORD_INCR;
  prev = ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 6) );
  *in += (increment << 8);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 ); */
  *in += carry;

  return;
}

static void
add_14_05 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 4 * WORD_INCR;
  prev = ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 4) );
  *in += (increment << 6);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 ); */
  *in += carry;

  return;
}

static void
add_14_06 (UINT2 *in, UINT2 increment) {
  UINT2 prev, carry;

  in += 5 * WORD_INCR;
  prev = ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;
  carry = ( (prev + increment) >> (14 - 2) );
  *in += (increment << 4);

  in += 1 * WORD_INCR;
  /* prev |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 ); */
  *in += carry;

  return;
}

static void
add_14_07 (UINT2 *in, UINT2 increment) {
  in += 6 * WORD_INCR;
  *in += (increment << 2);
  return;
}



static void
add_16_00 (UINT2 *in, UINT2 increment) {
  *in += increment;
  return;
}

static void
add_16_01 (UINT2 *in, UINT2 increment) {
  in += 1 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_16_02 (UINT2 *in, UINT2 increment) {
  in += 2 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_16_03 (UINT2 *in, UINT2 increment) {
  in += 3 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_16_04 (UINT2 *in, UINT2 increment) {
  in += 4 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_16_05 (UINT2 *in, UINT2 increment) {
  in += 5 * WORD_INCR;
  *in += increment;
  return;
}

static void
add_16_06 (UINT2 *in, UINT2 increment) {
  in += 6 * WORD_INCR;
  *in += increment;
  return;
}


static void
add_16_07 (UINT2 *in, UINT2 increment) {
  in += 7 * WORD_INCR;
  *in += increment;
  return;
}


typedef void (*Adder_T) (UINT2 *, UINT2);

static Adder_T adder_table[72] =
  {add_00, add_00, add_00, add_00,
   add_00, add_00, add_00, add_00,

   add_02_00, add_02_01, add_02_02, add_02_03,
   add_02_04, add_02_05, add_02_06, add_02_07,

   add_04_00, add_04_01, add_04_02, add_04_03,
   add_04_04, add_04_05, add_04_06, add_04_07,

   add_06_00, add_06_01, add_06_02, add_06_03,
   add_06_04, add_06_05, add_06_06, add_06_07,

   add_08_00, add_08_01, add_08_02, add_08_03,
   add_08_04, add_08_05, add_08_06, add_08_07,

   add_10_00, add_10_01, add_10_02, add_10_03,
   add_10_04, add_10_05, add_10_06, add_10_07,

   add_12_00, add_12_01, add_12_02, add_12_03,
   add_12_04, add_12_05, add_12_06, add_12_07,

   add_14_00, add_14_01, add_14_02, add_14_03,
   add_14_04, add_14_05, add_14_06, add_14_07,

   add_16_00, add_16_01, add_16_02, add_16_03,
   add_16_04, add_16_05, add_16_06, add_16_07,
  };



/* Assumes that increment > 0 */
void
Epu16_bitpack64_add_bitpack (Localspace_T oligo, int packsize, UINT2 *bitpack, UINT2 increment) {
  int nregisters, remainder;
  int index, column;
#ifdef CHECK
  UINT2 prev;
#endif

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*8 + remainder/8;
  column = remainder % 8;

#ifdef CHECK
  prev = Epu16_bitpack64_access(oligo,packsize,bitpack);
  (adder_table[index])(&(bitpack[column]),increment);
  if (Epu16_bitpack64_access(oligo,packsize,bitpack) != prev + increment) {
    printf("Error at Epu16_bitpack64_add_bitpack with oligo %llu, packsize %d, remainder %d, column %d, index %d\n",
	   oligo,packsize,remainder,column,index);
    abort();
  }
#else
  (adder_table[index])(&(bitpack[column]),increment);
#endif

  return;
}



static void
transfer_00_02 (UINT2 *out, const UINT2 *in) {
  return;
}

static void
transfer_02_04 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

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
  out += WORD_INCR;
  *out |= (value << 0);

  /* 05 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 2 ) ;
  *out |= (value << 4);

  /* 06 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 2 ) ;
  *out |= (value << 8);

  /* 07 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 2 ) ;
  *out |= (value << 12);

  return;
}

static void
transfer_04_06 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 4 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 4 ) ;
  *out |= (value << 6);

  /* 02 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 4 ) ;
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (6 - 2));

  /* 03 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 4 ) ;
  *out |= (value << 2);

  /* 04 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0 )   % (1U << 4 ) ;
  *out |= (value << 8);

  /* 05 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 4 ) ;
  *out |= (value << 14);
  out += WORD_INCR;
  *out |= (value >> (6 - 4));

  /* 06 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 4 ) ;
  *out |= (value << 4);

  /* 07 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 4 ) ;
  *out |= (value << 10);

  return;
}


static void
transfer_06_08 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 6 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 6 ) ;
  *out |= (value << 8);


  /* 02 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 );
  out += WORD_INCR;
  *out |= (value << 0);


  /* 03 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 6 ) ;
  *out |= (value << 8);


  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 6 ) ;
  out += WORD_INCR;
  *out |= (value << 0);


  /* 05 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 );
  *out |= (value << 8);

  /* 06 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 6 ) ;
  out += WORD_INCR;
  *out |= (value << 0);

  /* 07 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 6 ) ;
  *out |= (value << 8);

  return;
}



static void
transfer_08_10 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  value = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *out |= (value << 10);
  out += WORD_INCR;
  *out |= (value >> (10 - 4));

  /* 02 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *out |= (value << 4);

  /* 03 */
  value = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *out |= (value << 14);
  out += WORD_INCR;
  *out |= (value >> (10 - 8));

  /* 04 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *out |= (value << 8);
  out += WORD_INCR;
  *out |= (value >> (10 - 2));

  /* 05 */
  value = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *out |= (value << 2);

  /* 06 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (10 - 6));

  /* 07 */
  value = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
  *out |= (value << 6);

  return;
}


static void
transfer_10_12 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 10 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 );
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (12 - 8));

  /* 02 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 10 ) ;
  *out |= (value << 8);
  out += WORD_INCR;
  *out |= (value >> (12 - 4));

  /* 03 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 );
  *out |= (value << 4);

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 );
  out += WORD_INCR;
  *out |= (value << 0);

  /* 05 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 10 ) ;
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (12 - 8));

  /* 06 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 );
  *out |= (value << 8);
  out += WORD_INCR;
  *out |= (value >> (12 - 4));

  /* 07 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 10 ) ;
  *out |= (value << 4);

  return;
}


static void
transfer_12_14 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
  *out |= (value << 0);

  /* 01 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  *out |= (value << 14);
  out += WORD_INCR;
  *out |= (value >> (14 - 12));

  /* 02 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  *out |= (value << 12);
  out += WORD_INCR;
  *out |= (value >> (14 - 10));

  /* 03 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
  *out |= (value << 10);
  out += WORD_INCR;
  *out |= (value >> (14 - 8));

  /* 04 */
  in += WORD_INCR;
  value = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
  *out |= (value << 8);
  out += WORD_INCR;
  *out |= (value >> (14 - 6));

  /* 05 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  *out |= (value << 6);
  out += WORD_INCR;
  *out |= (value >> (14 - 4));

  /* 06 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  *out |= (value << 4);
  out += WORD_INCR;
  *out |= (value >> (14 - 2));

  /* 07 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
  *out |= (value << 2);

  return;
}


static void
transfer_14_16 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  value = ( CONVERT(*in) >>  0  )   % (1U << 14 ) ;
  *out |= value;

  /* 01 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 );
  out += WORD_INCR;
  *out |= value;

  /* 02 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 );
  out += WORD_INCR;
  *out |= value;

  /* 03 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 );
  out += WORD_INCR;
  *out |= value;

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 8 );
  out += WORD_INCR;
  *out |= value;

  /* 05 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 );
  out += WORD_INCR;
  *out |= value;

  /* 06 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 );
  out += WORD_INCR;
  *out |= value;

  /* 07 */
  value = ( CONVERT(*in) >>  2  )   % (1U << 14 ) ;
  out += WORD_INCR;
  *out |= value;

  return;
}



typedef void (*Transfer_T) (UINT2 *, const UINT2 *);

static Transfer_T transfer_table[8] =
  {transfer_00_02, transfer_02_04, transfer_04_06, transfer_06_08,
   transfer_08_10, transfer_10_12, transfer_12_14, transfer_14_16};

UINT2 *
Epu16_bitpack64_realloc_one (int packsize, UINT2 *bitpack) {
  UINT2 *new;
  int nregisters;
#ifdef CHECK
  int i;
#endif

  new = (UINT2 *) CALLOC((packsize + 2) / 2 * 8,sizeof(UINT2));

  nregisters = packsize/2;
  (transfer_table[nregisters])(&(new[/*column*/0]),&(bitpack[/*column*/0]));
  (transfer_table[nregisters])(&(new[/*column*/1]),&(bitpack[/*column*/1]));
  (transfer_table[nregisters])(&(new[/*column*/2]),&(bitpack[/*column*/2]));
  (transfer_table[nregisters])(&(new[/*column*/3]),&(bitpack[/*column*/3]));
  (transfer_table[nregisters])(&(new[/*column*/4]),&(bitpack[/*column*/4]));
  (transfer_table[nregisters])(&(new[/*column*/5]),&(bitpack[/*column*/5]));
  (transfer_table[nregisters])(&(new[/*column*/6]),&(bitpack[/*column*/6]));
  (transfer_table[nregisters])(&(new[/*column*/7]),&(bitpack[/*column*/7]));

#ifdef CHECK
  for (i = 0; i < 64; i++) {
    if (Epu16_bitpack64_access(i,packsize+2,new) != Epu16_bitpack64_access(i,packsize,bitpack)) {
      fprintf(stderr,"Difference in packsize %d -> %d, index %d: %u != %u\n",
	      packsize,packsize+2,i,
	      Epu16_bitpack64_access(i,packsize,bitpack),
	      Epu16_bitpack64_access(i,packsize+2,new));
      abort();
    }
  }
#endif

  FREE(bitpack);
  return new;
}


UINT2 *
Epu16_bitpack64_realloc_multiple (int old_packsize, int new_packsize, UINT2 *bitpack) {
  UINT2 *new;
  int nregisters, packsize;
#ifdef CHECK
  int i;
#endif

  for (packsize = old_packsize; packsize < new_packsize; packsize += 2) {
    new = (UINT2 *) CALLOC((packsize + 2) / 2 * 8,sizeof(UINT2));

    nregisters = packsize/2;
    (transfer_table[nregisters])(&(new[/*column*/0]),&(bitpack[/*column*/0]));
    (transfer_table[nregisters])(&(new[/*column*/1]),&(bitpack[/*column*/1]));
    (transfer_table[nregisters])(&(new[/*column*/2]),&(bitpack[/*column*/2]));
    (transfer_table[nregisters])(&(new[/*column*/3]),&(bitpack[/*column*/3]));
    (transfer_table[nregisters])(&(new[/*column*/4]),&(bitpack[/*column*/4]));
    (transfer_table[nregisters])(&(new[/*column*/5]),&(bitpack[/*column*/5]));
    (transfer_table[nregisters])(&(new[/*column*/6]),&(bitpack[/*column*/6]));
    (transfer_table[nregisters])(&(new[/*column*/7]),&(bitpack[/*column*/7]));

#ifdef CHECK
    for (i = 0; i < 64; i++) {
      if (Epu16_bitpack64_access(i,packsize+2,new) != Epu16_bitpack64_access(i,packsize,bitpack)) {
	fprintf(stderr,"Difference in packsize %d -> %d, index %d: %u != %u\n",
		packsize,packsize+2,i,
		Epu16_bitpack64_access(i,packsize,bitpack),
		Epu16_bitpack64_access(i,packsize+2,new));
	abort();
      }
    }
#endif

    FREE(bitpack);
    bitpack = new;
  }

  return bitpack;
}

