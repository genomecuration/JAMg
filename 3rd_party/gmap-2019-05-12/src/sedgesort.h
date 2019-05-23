/* $Id: sedgesort.h 216888 2018-10-07 16:37:57Z twu $ */
#ifndef SEDGESORT_INCLUDED
#define SEDGESORT_INCLUDED

#include "types.h"

extern void
Sedgesort_uint4 (register unsigned int array[], register int len);
#ifdef LARGE_GENOMES
extern void
Sedgesort_uint8 (register UINT8 array[], register int len);
#endif

extern int *
Sedgesort_order (register unsigned int array[], register int len);
extern int *
Sedgesort_order_int (register int array[], register int len);

#endif

