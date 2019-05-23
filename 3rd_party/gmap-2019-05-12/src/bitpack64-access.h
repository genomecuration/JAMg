/* $Id: bitpack64-access.h 212659 2018-01-20 00:58:14Z twu $ */
#ifndef BITPACK64_ACCESS_INCLUDED
#define BITPACK64_ACCESS_INCLUDED
#include "types.h"
#include "bool.h"

/* For reading direct-coded bitstreams */
extern UINT4
Bitpack64_access (Oligospace_T oligo, UINT4 *ptrs, UINT4 *comp);

extern UINT4
Bitpack64_access_bitpack (Oligospace_T oligo, int packsize, UINT4 *bitpack);

extern bool
Bitpack64_access_filledp (Oligospace_T oligo, int packsize, UINT4 *bitpack);

extern int
Bitpack64_access_new_packsize (Oligospace_T oligo, int old_packsize, UINT4 *bitpack, int increment);

extern void
Bitpack64_extract_bitpack (UINT4 *out, int packsize, UINT4 *bitpack);

#endif
