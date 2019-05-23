/* $Id: epu16-bitpack64-incr.h 214305 2018-03-19 23:40:43Z twu $ */
#ifndef EPU16_BITPACK64_INCR_INCLUDED
#define EPU16_BITPACK64_INCR_INCLUDED
#include "types.h"

extern UINT2
Epu16_bitpack64_incr (Localspace_T oligo, int packsize, UINT2 *bitpack);

extern UINT2 *
Epu16_bitpack64_realloc_one (int packsize, UINT2 *bitpack);

extern UINT2 *
Epu16_bitpack64_realloc_multiple (int old_packsize, int new_packsize, UINT2 *bitpack);

extern void
Epu16_bitpack64_add_bitpack (Localspace_T oligo, int packsize, UINT2 *bitpack, UINT2 increment);

#endif

