/* $Id: epu16-bitpack64-access.h 214305 2018-03-19 23:40:43Z twu $ */
#ifndef EPU16_BITPACK64_ACCESS_INCLUDED
#define EPU16_BITPACK64_ACCESS_INCLUDED
#include "types.h"
#include "bool.h"


extern UINT2
Epu16_bitpack64_access (Localspace_T oligo, int packsize, UINT2 *bitpack);

extern bool
Epu16_bitpack64_access_filledp (Localspace_T oligo, int packsize, UINT2 *bitpack);

extern int
Epu16_bitpack64_access_new_packsize (Localspace_T oligo, int old_packsize, UINT2 *bitpack, int increment);

extern void
Epu16_bitpack64_extract (UINT2 *out, int packsize, UINT2 *bitpack);

#endif

