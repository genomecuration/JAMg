/* $Id: bitpack64-incr.h 212659 2018-01-20 00:58:14Z twu $ */
#ifndef BITPACK64_INCR_INCLUDED
#define BITPACK64_INCR_INCLUDED
#include "types.h"

/* For reading direct-coded bitstreams.  Horizontal format. */
extern UINT4
Bitpack64_incr (Oligospace_T oligo, UINT4 *ptrs, UINT4 *comp);

extern UINT4
Bitpack64_incr_bitpack (Oligospace_T oligo, int packsize, UINT4 *bitpack);

extern void
Bitpack64_add (Oligospace_T oligo, UINT4 *ptrs, UINT4 *comp, UINT4 increment);

extern void
Bitpack64_add_bitpack (Oligospace_T oligo, int packsize, UINT4 *bitpack, UINT4 increment);

#if 0
extern void
Bitpack64_sub_bitpack (Oligospace_T oligo, int packsize, UINT4 *bitpack, UINT4 decrement);
#endif

extern UINT4 *
Bitpack64_realloc_one (int packsize, UINT4 *bitpack);

extern UINT4 *
Bitpack64_realloc_multiple (int old_packsize, int new_packsize, UINT4 *bitpack);

#endif
