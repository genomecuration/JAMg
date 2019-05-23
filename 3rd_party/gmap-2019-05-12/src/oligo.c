static char rcsid[] = "$Id: oligo.c 214305 2018-03-19 23:40:43Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "oligo.h" 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "assert.h"
#include "mem.h"
#include "types.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/************************************************************************
 *   Check
 ************************************************************************/

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

#if 0
/* For debugging purposes only.  Called by Kmer_search procedures */
char *
Oligo_one_nt (Oligospace_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Oligospace_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return nt;
}
#endif


static int oligosize;

/* Number of possible genestrands: 3.  Number of possible nucleotides: 4 */
static Oligospace_T forward_conv_5[3][4];
static Oligospace_T revcomp_conv_5[3][4];
static Oligospace_T forward_conv_3[3][4];
static Oligospace_T revcomp_conv_3[3][4];

static int char_to_int[128] =
  { -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 1..10 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 11..20 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 21..30 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 31..40 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 41..50 */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 51..60 */
    -1, -1, -1, -1,

  /* A,  B,  C,  D,  E,  F,  G,  H,  I,  J, */
     0, -1,  1, -1, -1, -1,  2, -1, -1, -1, /* 65..74 */

  /* K,  L,  M,  N,  O,  P,  Q,  R,  S,  T  */
    -1, -1, -1, -1, -1, -1, -1, -1, -1,  3, /* 75..84 */

  /* U,  V,  W,  X,  Y,  Z  */
    -1, -1, -1, -1, -1, -1,

    -1, -1, -1, -1, -1, -1,

  /* a,  b,  c,  d,  e,  f,  g,  h,  i,  j, */
     0, -1,  1, -1, -1, -1,  2, -1, -1, -1, /* 97..106 */

  /* k,  l,  m,  n,  o,  p,  q,  r,  s,  t  */
    -1, -1, -1, -1, -1, -1, -1, -1, -1,  3, /* 107..116 */

  /* u,  v,  w,  x,  y,  z  */
    -1, -1, -1, -1, -1, -1,

    -1, -1, -1, -1, -1};


static Oligostate_T
oligo_read_5 (int *querypos, Oligospace_T *forward, Oligospace_T *revcomp, 
	      Reader_T reader, int genestrand) {
  int count = 0;
  int c, d;


  *forward = *revcomp = 0U;
  while (count < oligosize && (c = Reader_getc_5(reader)) != '\0') {

#if 0
    switch (c) {
    case 'A': *forward <<= 2; *revcomp >>= 2; *revcomp |= LEFT_T; break;
    case 'C': *forward <<= 2; *forward |= RIGHT_C; 
      *revcomp >>= 2; *revcomp |= LEFT_G;  break;
    case 'G': *forward <<= 2; *forward |= RIGHT_G;
      *revcomp >>= 2; *revcomp |= LEFT_C; break;
    case 'T': *forward <<= 2; *forward |= RIGHT_T; *revcomp >>= 2; break;
    default: *forward = *revcomp = 0U; count = -1; 
      /* This counteracts count++ below */
    }

#else
    if ((d = char_to_int[c]) < 0) {
      *forward = *revcomp = 0U;
      count = -1;  /* This counteracts count++ below */
    } else {
      *forward <<= 2;
      *forward |= forward_conv_5[genestrand][d];

      *revcomp >>= 2;
      *revcomp |= revcomp_conv_5[genestrand][d];
    }
#endif

    count++;
    debug(printf("5' Read %c at %d, count = %d, oligo = %016lX, %016lX\n",
		 c,Reader_startpos(reader) - 1,count,*forward,*revcomp));
  }

  if (count < oligosize) {
    *forward = *revcomp = 0U;
    return DONE;
  } else {
    *querypos = Reader_startpos(reader) - oligosize;
    debug(printf("Read: Returning oligo %016lX for querypos %d, count = %d, and size = %d\n",
		 *forward,*querypos,count,oligosize));
    debug(printf("Setting querypos to be %u - %u = %u\n",
		 Reader_startpos(reader),oligosize,*querypos));
#if 0
    /* Holds only for STANDARD mode */
    assert((*forward & ~(~0U << 2*oligosize)) == Reader_check(reader,*querypos,oligosize));
#endif

    return VALID;
  }
}

static Oligostate_T
oligo_revise_5 (int *querypos, Oligospace_T *forward, Oligospace_T *revcomp, 
		Reader_T reader, int genestrand) {
  int c, d;

  if ((c = Reader_getc_5(reader)) == '\0') {
    *forward = *revcomp = 0U;
    debug(printf("5' Revision: read terminating char '0'\n"));
    return DONE;

  } else {
#if 0
    /* Old code */
    switch (c) {
    case 'A': *forward <<= 2; *revcomp >>= 2; *revcomp |= LEFT_T; break;
    case 'C': *forward <<= 2; *forward |= RIGHT_C; 
      *revcomp >>= 2; *revcomp |= LEFT_G; break;
    case 'G': *forward <<= 2; *forward |= RIGHT_G;
      *revcomp >>= 2; *revcomp |= LEFT_C; break;
    case 'T': *forward <<= 2; *forward |= RIGHT_T; *revcomp >>= 2; break;
    default: *forward = *revcomp = 0U; 
      *querypos = Reader_startpos(reader) - oligosize;
      return INVALID;
    }
#else
    if ((d = char_to_int[c]) < 0) {
      *forward = *revcomp = 0U;
      *querypos = Reader_startpos(reader) - oligosize;
      return INVALID;

    } else {
      *forward <<= 2;
      *forward |= forward_conv_5[genestrand][d];

      *revcomp >>= 2;
      *revcomp |= revcomp_conv_5[genestrand][d];
    }
#endif

    *querypos = Reader_startpos(reader) - oligosize;
    debug(printf("5' Revision: read char %c at %d, oligo = %016lX, %016lX at querypos %d\n",
		 c,Reader_startpos(reader) - 1,*forward,*revcomp,*querypos));
#if 0
    /* Holds only for STANDARD mode */
    assert((*forward & ~(~0U << 2*oligosize)) == Reader_check(reader,*querypos,oligosize));
#endif
    return VALID;
  }
}

Oligostate_T
Oligo_next_5 (Oligostate_T last_state, int *querypos, Oligospace_T *forward,
	      Oligospace_T *revcomp, Reader_T reader, int genestrand) {
  Oligostate_T state;

  if (last_state == DONE) {
    /* fprintf(stderr,"Called Oligo_next with last_state == DONE\n"); */
    return DONE;
  } else if (last_state != VALID) {
    /* INVALID or INIT */
    return oligo_read_5(&(*querypos),&(*forward),&(*revcomp),reader,genestrand);

  } else if ((state = oligo_revise_5(&(*querypos),&(*forward),&(*revcomp),reader,genestrand)) == DONE) {
    return DONE;
  } else if (state == VALID) {
    return VALID;
  } else {
    /* state == INVALID */
    return oligo_read_5(&(*querypos),&(*forward),&(*revcomp),reader,genestrand);
  }

  return state;
}


static Oligostate_T
oligo_read_3 (int *querypos, Oligospace_T *forward, Oligospace_T *revcomp, 
   	      Reader_T reader, int genestrand) {
  int count = 0;
  int c, d;

  *forward = *revcomp = 0U;
  while (count < oligosize && (c = Reader_getc_3(reader)) != '\0') {
#if 0
    /* Old code */
    switch (c) {
    case 'A': *forward >>= 2; *revcomp <<= 2; *revcomp |= RIGHT_T; break;
    case 'C': *forward >>= 2; *forward |= LEFT_C; 
      *revcomp <<= 2; *revcomp |= RIGHT_G;  break;
    case 'G': *forward >>= 2; *forward |= LEFT_G;
      *revcomp <<= 2; *revcomp |= RIGHT_C; break;
    case 'T': *forward >>= 2; *forward |= LEFT_T; *revcomp <<= 2; break;
    default: *forward = *revcomp = 0U; count = -1; 
      /* This counteracts count++ below */
    }
#else
    if ((d = char_to_int[c]) < 0) {
      *forward = *revcomp = 0U;
      count = -1;  /* This counteracts count++ below */
    } else {
      *forward >>= 2;
      *forward |= forward_conv_3[genestrand][d];

      *revcomp <<= 2;
      *revcomp |= revcomp_conv_3[genestrand][d];
    }
#endif

    count++;
    debug(printf("3' Read %c at %d, count = %d, oligo = %016lX, %016lX\n",
		 c,Reader_endpos(reader) + 1,count,*forward,*revcomp));
  }

  if (count < oligosize) {
    *forward = *revcomp = 0U;
    return DONE;
  } else {
    *querypos = Reader_endpos(reader) + 1;
    debug(printf("Read: Returning oligo %016lX for querypos %d, count = %d, and size = %d\n",
		 *forward,*querypos,count,oligosize));
    debug(printf("Setting querypos to be %u + 1 = %u\n",Reader_endpos(reader),*querypos));
#if 0
    /* Holds only for STANDARD mode */
    assert(((*forward >> (64 - 2*oligosize)) & ~(~0U << 2*oligosize)) == Reader_check(reader,*querypos,oligosize));
#endif
    return VALID;
  }
}


static Oligostate_T
oligo_revise_3 (int *querypos, Oligospace_T *forward, Oligospace_T *revcomp, 
		Reader_T reader, int genestrand) {
  int c, d;

  if ((c = Reader_getc_3(reader)) == '\0') {
    *forward = *revcomp = 0U;
    debug(printf("3' Revision: read terminating char '0'\n"));
    return DONE;

  } else {
#if 0
    /* Old code */
    switch (c) {
    case 'A': *forward >>= 2; *revcomp <<= 2; *revcomp |= RIGHT_T; break;
    case 'C': *forward >>= 2; *forward |= LEFT_C; 
      *revcomp <<= 2; *revcomp |= RIGHT_G; break;
    case 'G': *forward >>= 2; *forward |= LEFT_G;
      *revcomp <<= 2; *revcomp |= RIGHT_C; break;
    case 'T': *forward >>= 2; *forward |= LEFT_T; *revcomp <<= 2; break;
    default: *forward = *revcomp = 0U; 
      *querypos = Reader_endpos(reader) + 1;
      return INVALID;
    }
#else
    if ((d = char_to_int[c]) < 0) {
      *forward = *revcomp = 0U;
      *querypos = Reader_endpos(reader) + 1;
      return INVALID;
    } else {
      *forward >>= 2;
      *forward |= forward_conv_3[genestrand][d];

      *revcomp <<= 2;
      *revcomp |= revcomp_conv_3[genestrand][d];
    }
#endif

    *querypos = Reader_endpos(reader) + 1;
    debug(printf("3' Revision: read char %c at %d, oligo = %016lX, %016lX at querypos %d\n",
		 c,Reader_endpos(reader),*forward,*revcomp,*querypos));
#if 0
    /* Holds only for STANDARD mode */
    assert(((*forward >> (64 - 2*oligosize)) & ~(~0U << 2*oligosize)) == Reader_check(reader,*querypos,oligosize));
#endif
    return VALID;
  }
}


Oligostate_T
Oligo_next_3 (Oligostate_T last_state, int *querypos, Oligospace_T *forward,
	      Oligospace_T *revcomp, Reader_T reader, int genestrand) {
  Oligostate_T state;

  if (last_state == DONE) {
    /* fprintf(stderr,"Called Oligo_next with last_state == DONE\n"); */
    return DONE;
  } else if (last_state != VALID) {
    /* INVALID or INIT */
    return oligo_read_3(&(*querypos),&(*forward),&(*revcomp),reader,genestrand);

  } else if ((state = oligo_revise_3(&(*querypos),&(*forward),&(*revcomp),reader,genestrand)) == DONE) {
    return DONE;
  } else if (state == VALID) {
    return VALID;
  } else {
    /* state == INVALID */
    return oligo_read_3(&(*querypos),&(*forward),&(*revcomp),reader,genestrand);
  }
}


#ifndef GSNAP
/* For GMAP */
Oligostate_T
Oligo_skip_5 (Oligostate_T last_state, int *querypos, Oligospace_T *forward,
	      Oligospace_T *revcomp, Reader_T reader, int genestrand, int nskip) {
  int i = 0;

  while (i < nskip && last_state != DONE) {
    if (last_state == VALID) {
      last_state = oligo_revise_5(&(*querypos),&(*forward),&(*revcomp),reader,genestrand);
      i++;
    } else {			/* INVALID and INIT */
      last_state = oligo_read_5(&(*querypos),&(*forward),&(*revcomp),reader,genestrand);
      i += oligosize;
    }
  }
  return last_state;
}
#endif


#ifndef GSNAP
/* For GMAP */
Oligostate_T
Oligo_skip_3 (Oligostate_T last_state, int *querypos, Oligospace_T *forward,
	      Oligospace_T *revcomp, Reader_T reader, int genestrand, int nskip) {
  int i = 0;

  while (i < nskip && last_state != DONE) {
    if (last_state == VALID) {
      last_state = oligo_revise_3(&(*querypos),&(*forward),&(*revcomp),reader,genestrand);
      i++;
    } else {			/* INVALID and INIT */
      last_state = oligo_read_3(&(*querypos),&(*forward),&(*revcomp),reader,genestrand);
      i += oligosize;
    }
  }
  return last_state;
}
#endif


/* For mode STANDARD, CMET_STRANDED, ATOI_STRANDED, and TTOC_STRANDED, we set genestrand to be 0 */
/* For mode CMET_NONSTRANDED, ATOI_NONSTRANDED, and TTOC_NONSTRANDED, we set genestrand to be +1 for plus strand, and +2 for minus strand */
void
Oligo_setup (int index1part, int mode) {
  oligosize = index1part;

  if (mode == STANDARD) {
    memcpy(forward_conv_5[+0],((Oligospace_T[]){/*   */RIGHT_A, /*   */RIGHT_C, /*   */RIGHT_G, /*   */RIGHT_T}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_3[+0],((Oligospace_T[]){/*   */LEFT_A,  /*   */LEFT_C,  /*   */LEFT_G,  /*   */LEFT_T}),4*sizeof(Oligospace_T));

    memcpy(revcomp_conv_5[+0],((Oligospace_T[]){/*   */LEFT_T,  /*   */LEFT_G,  /*   */LEFT_C,  /*   */LEFT_A}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_3[+0],((Oligospace_T[]){/*   */RIGHT_T, /*   */RIGHT_G, /*   */RIGHT_C, /*   */RIGHT_A}),4*sizeof(Oligospace_T));

  } else if (mode == CMET_STRANDED) {
    memcpy(forward_conv_5[+0],((Oligospace_T[]){/*   */RIGHT_A, /*C>T*/RIGHT_T, /*   */RIGHT_G, /*   */RIGHT_T}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_3[+0],((Oligospace_T[]){/*   */LEFT_A,  /*C>T*/LEFT_T,  /*   */LEFT_G,  /*   */LEFT_T}),4*sizeof(Oligospace_T));

    memcpy(revcomp_conv_5[+0],((Oligospace_T[]){/*   */LEFT_T,  /*G>A*/LEFT_A,  /*   */LEFT_C,  /*   */LEFT_A}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_3[+0],((Oligospace_T[]){/*   */RIGHT_T, /*G>A*/RIGHT_A, /*   */RIGHT_C, /*   */RIGHT_A}),4*sizeof(Oligospace_T));


  } else if (mode == CMET_NONSTRANDED) {
    /* genestrand: +1 */
    /* Query from gplus strand: (C==T).  Matches genomicfwd kmer */
    /* Query from gminus strand: (C==T),  Then revcomp to match genomic genomicfwd kmer: (G==A) */

    /* genestrand: +2 (pcr step) */
    /* Query from gminus strand: (C==T). Then pcr: (G==A).  Matches genomicfwd kmer */
    /* Query from gplus strand: (C==T).  Then pcr: (G==A).  Then revcomp to match genomicfwd kmer: (C==T) */

    memcpy(forward_conv_5[+1],((Oligospace_T[]){/*   */RIGHT_A, /*C>T*/RIGHT_T, /*   */RIGHT_G, /*   */RIGHT_T}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_5[+2],((Oligospace_T[]){/*   */RIGHT_A, /*   */RIGHT_C, /*G>A*/RIGHT_A, /*   */RIGHT_T}),4*sizeof(Oligospace_T));

    memcpy(forward_conv_3[+1],((Oligospace_T[]){/*   */LEFT_A,  /*C>T*/LEFT_T,  /*   */LEFT_G,  /*   */LEFT_T}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_3[+2],((Oligospace_T[]){/*   */LEFT_A,  /*   */LEFT_C,  /*G>A*/LEFT_A,  /*   */LEFT_T}),4*sizeof(Oligospace_T));


    memcpy(revcomp_conv_5[+1],((Oligospace_T[]){/*   */LEFT_T,  /*G>A*/LEFT_A,  /*   */LEFT_C,  /*   */LEFT_A}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_5[+2],((Oligospace_T[]){/*   */LEFT_T,  /*   */LEFT_G,  /*C>T*/LEFT_T,  /*   */LEFT_A}),4*sizeof(Oligospace_T));

    memcpy(revcomp_conv_3[+1],((Oligospace_T[]){/*   */RIGHT_T, /*G>A*/RIGHT_A, /*   */RIGHT_C, /*   */RIGHT_A}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_3[+2],((Oligospace_T[]){/*   */RIGHT_T, /*   */RIGHT_G, /*C>T*/RIGHT_T, /*   */RIGHT_A}),4*sizeof(Oligospace_T));

  } else if (mode == ATOI_STRANDED) {
    memcpy(forward_conv_5[+0],((Oligospace_T[]){/*A>G*/RIGHT_G, /*   */RIGHT_C, /*   */RIGHT_G, /*   */RIGHT_T}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_3[+0],((Oligospace_T[]){/*A>G*/LEFT_G,  /*   */LEFT_C,  /*   */LEFT_G,  /*   */LEFT_T}),4*sizeof(Oligospace_T));

    memcpy(revcomp_conv_5[+0],((Oligospace_T[]){/*T>C*/LEFT_C,  /*   */LEFT_G,  /*   */LEFT_C,  /*   */LEFT_A}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_3[+0],((Oligospace_T[]){/*T>C*/RIGHT_C, /*   */RIGHT_G, /*   */RIGHT_C, /*   */RIGHT_A}),4*sizeof(Oligospace_T));

  } else if (mode == ATOI_NONSTRANDED) {
    memcpy(forward_conv_5[+1],((Oligospace_T[]){/*A>G*/RIGHT_G, /*   */RIGHT_C, /*   */RIGHT_G, /*   */RIGHT_T}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_5[+2],((Oligospace_T[]){/*   */RIGHT_A, /*   */RIGHT_C, /*   */RIGHT_G, /*T>C*/RIGHT_C}),4*sizeof(Oligospace_T));

    memcpy(forward_conv_3[+1],((Oligospace_T[]){/*A>G*/LEFT_G,  /*   */LEFT_C,  /*   */LEFT_G,  /*   */LEFT_T}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_3[+2],((Oligospace_T[]){/*   */LEFT_A,  /*   */LEFT_C,  /*   */LEFT_G,  /*T>C*/LEFT_C}),4*sizeof(Oligospace_T));


    memcpy(revcomp_conv_5[+1],((Oligospace_T[]){/*T>C*/LEFT_C,  /*   */LEFT_G,  /*   */LEFT_C,  /*   */LEFT_A}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_5[+2],((Oligospace_T[]){/*   */LEFT_T,  /*   */LEFT_G,  /*   */LEFT_C,  /*A>G*/LEFT_G}),4*sizeof(Oligospace_T));

    memcpy(revcomp_conv_3[+1],((Oligospace_T[]){/*T>C*/RIGHT_C, /*   */RIGHT_G, /*   */RIGHT_C, /*   */RIGHT_A}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_3[+2],((Oligospace_T[]){/*   */RIGHT_T, /*   */RIGHT_G, /*   */RIGHT_C, /*A>G*/RIGHT_G}),4*sizeof(Oligospace_T));

  } else if (mode == TTOC_STRANDED) {
    memcpy(forward_conv_5[+0],((Oligospace_T[]){/*   */RIGHT_A, /*   */RIGHT_C, /*   */RIGHT_G, /*T>C*/RIGHT_C}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_3[+0],((Oligospace_T[]){/*   */LEFT_A,  /*   */LEFT_C,  /*   */LEFT_G,  /*T>C*/LEFT_C}),4*sizeof(Oligospace_T));

    memcpy(revcomp_conv_5[+0],((Oligospace_T[]){/*   */LEFT_T,  /*   */LEFT_G,  /*   */LEFT_C,  /*A>G*/LEFT_G}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_3[+0],((Oligospace_T[]){/*   */RIGHT_T, /*   */RIGHT_G, /*   */RIGHT_C, /*A>G*/RIGHT_G}),4*sizeof(Oligospace_T));

  } else if (mode == TTOC_NONSTRANDED) {
    memcpy(forward_conv_5[+1],((Oligospace_T[]){/*   */RIGHT_A, /*   */RIGHT_C, /*   */RIGHT_G, /*T>C*/RIGHT_C}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_5[+2],((Oligospace_T[]){/*A>G*/RIGHT_G, /*   */RIGHT_C, /*   */RIGHT_G, /*   */RIGHT_T}),4*sizeof(Oligospace_T));

    memcpy(forward_conv_3[+1],((Oligospace_T[]){/*   */LEFT_A,  /*   */LEFT_C,  /*   */LEFT_G,  /*T>C*/LEFT_C}),4*sizeof(Oligospace_T));
    memcpy(forward_conv_3[+2],((Oligospace_T[]){/*A>G*/LEFT_G,  /*   */LEFT_C,  /*   */LEFT_G,  /*   */LEFT_T}),4*sizeof(Oligospace_T));


    memcpy(revcomp_conv_5[+1],((Oligospace_T[]){/*   */LEFT_T,  /*   */LEFT_G,  /*   */LEFT_C,  /*A>G*/LEFT_G}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_5[+2],((Oligospace_T[]){/*T>C*/LEFT_C,  /*   */LEFT_G,  /*   */LEFT_C,  /*   */LEFT_A}),4*sizeof(Oligospace_T));

    memcpy(revcomp_conv_3[+1],((Oligospace_T[]){/*   */RIGHT_T, /*   */RIGHT_G, /*   */RIGHT_C, /*A>G*/RIGHT_G}),4*sizeof(Oligospace_T));
    memcpy(revcomp_conv_3[+2],((Oligospace_T[]){/*T>C*/RIGHT_C, /*   */RIGHT_G, /*   */RIGHT_C, /*   */RIGHT_A}),4*sizeof(Oligospace_T));

  } else {
    fprintf(stderr,"Unexpected mode %d\n",mode);
    exit(9);
  }
}

