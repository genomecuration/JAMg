static char rcsid[] = "$Id: epu16-bitpack64-access.c 214305 2018-03-19 23:40:43Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "epu16-bitpack64-access.h"

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

#ifdef HORIZONTAL
#define WORD_INCR 1		/* 1 for horizontal; 8 for vertical */
#else
#define WORD_INCR 8
#endif

static UINT2
access_00 (const UINT2 *in) {
  return 0U;
}

static UINT2
access_02_00 (const UINT2 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 2 ) ;
}

static UINT2
access_02_01 (const UINT2 *in) {
  return ( CONVERT(*in) >>  2  )   % (1U << 2 ) ;
}

static UINT2
access_02_02 (const UINT2 *in) {
  return ( CONVERT(*in) >>  4  )   % (1U << 2 ) ;
}

static UINT2
access_02_03 (const UINT2 *in) {
  return ( CONVERT(*in) >>  6  )   % (1U << 2 ) ;
}

static UINT2
access_02_04 (const UINT2 *in) {
  return ( CONVERT(*in) >>  8  )   % (1U << 2 ) ;
}

static UINT2
access_02_05 (const UINT2 *in) {
  return ( CONVERT(*in) >>  10  )   % (1U << 2 ) ;
}

static UINT2
access_02_06 (const UINT2 *in) {
  return ( CONVERT(*in) >>  12  )   % (1U << 2 ) ;
}

static UINT2
access_02_07 (const UINT2 *in) {
  return ( CONVERT(*in) >>  14  )   % (1U << 2 ) ;
}


static UINT2
access_04_00 (const UINT2 *in) {
  return ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;
}

static UINT2
access_04_01 (const UINT2 *in) {
  return ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;
}

static UINT2
access_04_02 (const UINT2 *in) {
  return ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;
}

static UINT2
access_04_03 (const UINT2 *in) {
  return ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;
}

static UINT2
access_04_04 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;
}

static UINT2
access_04_05 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;
}

static UINT2
access_04_06 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;
}

static UINT2
access_04_07 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;
}


static UINT2
access_06_00 (const UINT2 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 6 ) ;
}

static UINT2
access_06_01 (const UINT2 *in) {
  return ( CONVERT(*in) >>  6  )   % (1U << 6 ) ;
}

static UINT2
access_06_02 (const UINT2 *in) {
  UINT2 out;

  out = ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 );
  return out;
}

static UINT2
access_06_03 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 6 ) ;
}

static UINT2
access_06_04 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  8  )   % (1U << 6 ) ;
}

static UINT2
access_06_05 (const UINT2 *in) {
  UINT2 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 );
  return out;
}

static UINT2
access_06_06 (const UINT2 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 6 ) ;
}

static UINT2
access_06_07 (const UINT2 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  10  )   % (1U << 6 ) ;
}


static UINT2
access_08_00 (const UINT2 *in) {
  return ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
}

static UINT2
access_08_01 (const UINT2 *in) {
  return ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
}

static UINT2
access_08_02 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
}

static UINT2
access_08_03 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
}

static UINT2
access_08_04 (const UINT2 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
}

static UINT2
access_08_05 (const UINT2 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
}

static UINT2
access_08_06 (const UINT2 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;
}

static UINT2
access_08_07 (const UINT2 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;
}



static UINT2
access_10_00 (const UINT2 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 10 ) ;
}

static UINT2
access_10_01 (const UINT2 *in) {
  UINT2 out;

  out = ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 );
  return out;
}

static UINT2
access_10_02 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 10 ) ;
}

static UINT2
access_10_03 (const UINT2 *in) {
  UINT2 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 );
  return out;
}

static UINT2
access_10_04 (const UINT2 *in) {
  UINT2 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 );
  return out;
}

static UINT2
access_10_05 (const UINT2 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 10 ) ;
}

static UINT2
access_10_06 (const UINT2 *in) {
  UINT2 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 );
  return out;
}

static UINT2
access_10_07 (const UINT2 *in) {
  in += 4 * WORD_INCR;
  return ( CONVERT(*in) >>  6  )   % (1U << 10 ) ;
}



static UINT2
access_12_00 (const UINT2 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
}

static UINT2
access_12_01 (const UINT2 *in) {
  UINT2 out;

  out = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  return out;
}

static UINT2
access_12_02 (const UINT2 *in) {
  UINT2 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  return out;
}

static UINT2
access_12_03 (const UINT2 *in) {
  in += 2 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
}

static UINT2
access_12_04 (const UINT2 *in) {
  in += 3 * WORD_INCR;
  return ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;
}

static UINT2
access_12_05 (const UINT2 *in) {
  UINT2 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  return out;
}

static UINT2
access_12_06 (const UINT2 *in) {
  UINT2 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  return out;
}

static UINT2
access_12_07 (const UINT2 *in) {
  in += 5 * WORD_INCR;
  return ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;
}


static UINT2
access_14_00 (const UINT2 *in) {
  return ( CONVERT(*in) >>  0  )   % (1U << 14 ) ;
}

static UINT2
access_14_01 (const UINT2 *in) {
  UINT2 out;

  out = ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 );
  return out;
}

static UINT2
access_14_02 (const UINT2 *in) {
  UINT2 out;

  in += 1 * WORD_INCR;
  out = ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 );
  return out;
}

static UINT2
access_14_03 (const UINT2 *in) {
  UINT2 out;

  in += 2 * WORD_INCR;
  out = ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 );
  return out;
}

static UINT2
access_14_04 (const UINT2 *in) {
  UINT2 out;

  in += 3 * WORD_INCR;
  out = ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 );
  return out;
}

static UINT2
access_14_05 (const UINT2 *in) {
  UINT2 out;

  in += 4 * WORD_INCR;
  out = ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 );
  return out;
}

static UINT2
access_14_06 (const UINT2 *in) {
  UINT2 out;

  in += 5 * WORD_INCR;
  out = ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;
  in += 1 * WORD_INCR;
  out |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 );
  return out;
}

static UINT2
access_14_07 (const UINT2 *in) {
  in += 6 * WORD_INCR;
  return ( CONVERT(*in) >>  2  )   % (1U << 14 ) ;
}


static UINT2
access_16_00 (const UINT2 *in) {
  return CONVERT(*in);
}

static UINT2
access_16_01 (const UINT2 *in) {
  in += 1 * WORD_INCR;
  return CONVERT(*in);
}

static UINT2
access_16_02 (const UINT2 *in) {
  in += 2 * WORD_INCR;
  return CONVERT(*in);
}

static UINT2
access_16_03 (const UINT2 *in) {
  in += 3 * WORD_INCR;
  return CONVERT(*in);
}

static UINT2
access_16_04 (const UINT2 *in) {
  in += 4 * WORD_INCR;
  return CONVERT(*in);
}

static UINT2
access_16_05 (const UINT2 *in) {
  in += 5 * WORD_INCR;
  return CONVERT(*in);
}

static UINT2
access_16_06 (const UINT2 *in) {
  in += 6 * WORD_INCR;
  return CONVERT(*in);
}

static UINT2
access_16_07 (const UINT2 *in) {
  in += 7 * WORD_INCR;
  return CONVERT(*in);
}


typedef UINT2 (*Accessor_T) (const UINT2 *);

static Accessor_T accessor_table[72] =
  {access_00, access_00, access_00, access_00,
   access_00, access_00, access_00, access_00,

   access_02_00, access_02_01, access_02_02, access_02_03,
   access_02_04, access_02_05, access_02_06, access_02_07,

   access_04_00, access_04_01, access_04_02, access_04_03,
   access_04_04, access_04_05, access_04_06, access_04_07,

   access_06_00, access_06_01, access_06_02, access_06_03,
   access_06_04, access_06_05, access_06_06, access_06_07,

   access_08_00, access_08_01, access_08_02, access_08_03,
   access_08_04, access_08_05, access_08_06, access_08_07,

   access_10_00, access_10_01, access_10_02, access_10_03,
   access_10_04, access_10_05, access_10_06, access_10_07,

   access_12_00, access_12_01, access_12_02, access_12_03,
   access_12_04, access_12_05, access_12_06, access_12_07,

   access_14_00, access_14_01, access_14_02, access_14_03,
   access_14_04, access_14_05, access_14_06, access_14_07,

   access_16_00, access_16_01, access_16_02, access_16_03,
   access_16_04, access_16_05, access_16_06, access_16_07,

  };



UINT2
Epu16_bitpack64_access (Localspace_T oligo, int packsize, UINT2 *bitpack) {
  int nregisters, remainder;	/* nregisters is same as nwritten */
  int index, column;

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;

  index = nregisters*8 + remainder/8;
  column = remainder % 8;

  return (accessor_table[index])(&(bitpack[column]));
}

static UINT2 threshold[9] = 
  {/*0*/0U, /*2*/3U, /*4*/15U, /*6*/63U,
   /*8*/255U, /*10*/1023U, /*12*/4095U, /*14*/16383U,
   /*16*/65535U};

bool
Epu16_bitpack64_access_filledp (Localspace_T oligo, int packsize, UINT2 *bitpack) {
  int nregisters, remainder;	/* nregisters is same as nwritten */
  int index, column;
  UINT2 value;

  nregisters = packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*8 + remainder/8;
  column = remainder % 8;

  value = (accessor_table[index])(&(bitpack[column]));
  if (value == threshold[nregisters]) {
    return true;
  } else {
    return false;
  }
}

int
Epu16_bitpack64_access_new_packsize (Localspace_T oligo, int old_packsize, UINT2 *bitpack, int increment) {
  int new_packsize;
  int nregisters, remainder;	/* nregisters is same as nwritten */
  int index, column;
  UINT2 value;
  int firstbit;
#ifdef HAVE_BUILTIN_CLZ
#elif defined(HAVE_ASM_BSR)
  int msb;
#endif


  nregisters = old_packsize / 2;
  remainder = oligo % BLOCKSIZE;
  index = nregisters*8 + remainder/8;
  column = remainder % 8;

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



static void
extract_00 (UINT2 *out, const UINT2 *in) {
  /* 00 */
  out[0] = 0;

  /* 01 */
  out[8] = 0;

  /* 02 */
  out[16] = 0;

  /* 03 */
  out[24] = 0;

  /* 04 */
  out[32] = 0;

  /* 05 */
  out[40] = 0;

  /* 06 */
  out[48] = 0;

  /* 07 */
  out[56] = 0;

  return;
}


static void
extract_02 (UINT2 *out, const UINT2 *in) {

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 2 ) ;

  /* 01 */
  out[8] = ( CONVERT(*in) >>  2  )   % (1U << 2 ) ;

  /* 02 */
  out[16] = ( CONVERT(*in) >>  4  )   % (1U << 2 ) ;

  /* 03 */
  out[24] = ( CONVERT(*in) >>  6  )   % (1U << 2 ) ;

  /* 04 */
  out[32] = ( CONVERT(*in) >>  8 )   % (1U << 2 ) ;

  /* 05 */
  out[40] = ( CONVERT(*in) >>  10  )   % (1U << 2 ) ;

  /* 06 */
  out[48] = ( CONVERT(*in) >>  12  )   % (1U << 2 ) ;

  /* 07 */
  out[56] = ( CONVERT(*in) >>  14  )   % (1U << 2 ) ;

  return;
}


static void
extract_04 (UINT2 *out, const UINT2 *in) {

  /* 00 */
  out[0] = ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;

  /* 01 */
  out[8] = ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;

  /* 02 */
  out[16] = ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;

  /* 03 */
  out[24] = ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;

  /* 04 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >> 0 )   % (1U << 4 ) ;

  /* 05 */
  out[40] = ( CONVERT(*in) >> 4 )   % (1U << 4 ) ;

  /* 06 */
  out[48] = ( CONVERT(*in) >> 8 )   % (1U << 4 ) ;

  /* 07 */
  out[56] = ( CONVERT(*in) >> 12 )   % (1U << 4 ) ;

  return;
}


static void
extract_06 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 6 ) ;

  /* 01 */
  out[8] = ( CONVERT(*in) >>  6  )   % (1U << 6 ) ;

  /* 02 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 6 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 6 - 2 );
  out[16] = value;

  /* 03 */
  out[24] = ( CONVERT(*in) >>  2  )   % (1U << 6 ) ;

  /* 04 */
  out[32] = ( CONVERT(*in) >>  8  )   % (1U << 6 ) ;

  /* 05 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 6 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 6 - 4 );
  out[40] = value;

  /* 06 */
  out[48] = ( CONVERT(*in) >>  4  )   % (1U << 6 ) ;

  /* 07 */
  out[56] = ( CONVERT(*in) >>  10  )   % (1U << 6 ) ;

  return;
}


static void
extract_08 (UINT2 *out, const UINT2 *in) {

  /* 00 */
  out[0] = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;

  /* 01 */
  out[8] = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;

  /* 02 */
  in += WORD_INCR;
  out[16] = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;

  /* 03 */
  out[24] = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;

  /* 04 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;

  /* 05 */
  out[40] = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;

  /* 06 */
  in += WORD_INCR;
  out[48] = ( CONVERT(*in) >> 0 )   % (1U << 8 ) ;

  /* 07 */
  out[56] = ( CONVERT(*in) >> 8 )   % (1U << 8 ) ;

  return;
}


static void
extract_10 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 10 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 10 - 4 );
  out[8] = value;

  /* 02 */
  out[16] = ( CONVERT(*in) >>  4  )   % (1U << 10 ) ;

  /* 03 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 10 - 8 );
  out[24] = value;

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 10 - 2 );
  out[32] = value;

  /* 05 */
  out[40] = ( CONVERT(*in) >>  2  )   % (1U << 10 ) ;

  /* 06 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 10 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 10 - 6 );
  out[48] = value;

  /* 07 */
  out[56] = ( CONVERT(*in) >>  6  )   % (1U << 10 ) ;

  return;
}


static void
extract_12 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  out[8] = value;

  /* 02 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  out[16] = value;

  /* 03 */
  out[24] = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;

  /* 04 */
  in += WORD_INCR;
  out[32] = ( CONVERT(*in) >>  0  )   % (1U << 12 ) ;

  /* 05 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 12 - 8 );
  out[40] = value;

  /* 06 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 12 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 12 - 4 );
  out[48] = value;

  /* 07 */
  out[56] = ( CONVERT(*in) >>  4  )   % (1U << 12 ) ;

  return;
}


static void
extract_14 (UINT2 *out, const UINT2 *in) {
  UINT2 value;

  /* 00 */
  out[0] = ( CONVERT(*in) >>  0  )   % (1U << 14 ) ;

  /* 01 */
  value = ( CONVERT(*in) >>  14  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 12 ))<<( 14 - 12 );
  out[8] = value;

  /* 02 */
  value = ( CONVERT(*in) >>  12  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 10 ))<<( 14 - 10 );
  out[16] = value;

  /* 03 */
  value = ( CONVERT(*in) >>  10  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 8 ))<<( 14 - 8 );
  out[24] = value;

  /* 04 */
  value = ( CONVERT(*in) >>  8  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 6 ))<<( 14 - 6 );
  out[32] = value;

  /* 05 */
  value = ( CONVERT(*in) >>  6  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 4 ))<<( 14 - 4 );
  out[40] = value;

  /* 06 */
  value = ( CONVERT(*in) >>  4  )   % (1U << 14 ) ;
  in += WORD_INCR;
  value |= (CONVERT(*in) % (1U<< 2 ))<<( 14 - 2 );
  out[48] = value;

  /* 07 */
  out[56] = ( CONVERT(*in) >>  2  )   % (1U << 14 ) ;

  return;
}


static void
extract_16 (UINT2 *out, const UINT2 *in) {

  /* 00 */
  out[0] = CONVERT(*in);

  /* 01 */
  in += WORD_INCR;
  out[8] = CONVERT(*in);

  /* 02 */
  in += WORD_INCR;
  out[16] = CONVERT(*in);

  /* 03 */
  in += WORD_INCR;
  out[24] = CONVERT(*in);

  /* 04 */
  in += WORD_INCR;
  out[32] = CONVERT(*in);

  /* 05 */
  in += WORD_INCR;
  out[40] = CONVERT(*in);

  /* 06 */
  in += WORD_INCR;
  out[48] = CONVERT(*in);

  /* 07 */
  in += WORD_INCR;
  out[56] = CONVERT(*in);

  return;
}


typedef void (*Extractor_T) (UINT2 *, const UINT2 *);

static Extractor_T extractor_table[9] =
  {extract_00, extract_02, extract_04, extract_06,
   extract_08, extract_10, extract_12, extract_14,
   extract_16};

void
Epu16_bitpack64_extract (UINT2 *out, int packsize, UINT2 *bitpack) {
  int nregisters;

  nregisters = packsize / 2;
  (extractor_table[nregisters])(&(out[0]),&(bitpack[/*column*/0]));
  (extractor_table[nregisters])(&(out[1]),&(bitpack[/*column*/1]));
  (extractor_table[nregisters])(&(out[2]),&(bitpack[/*column*/2]));
  (extractor_table[nregisters])(&(out[3]),&(bitpack[/*column*/3]));
  (extractor_table[nregisters])(&(out[4]),&(bitpack[/*column*/4]));
  (extractor_table[nregisters])(&(out[5]),&(bitpack[/*column*/5]));
  (extractor_table[nregisters])(&(out[6]),&(bitpack[/*column*/6]));
  (extractor_table[nregisters])(&(out[7]),&(bitpack[/*column*/7]));

  return;
}
