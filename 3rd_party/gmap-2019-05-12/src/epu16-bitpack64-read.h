/* $Id: epu16-bitpack64-read.h 214305 2018-03-19 23:40:43Z twu $ */
#ifndef EPU16_BITPACK64_READ_INCLUDED
#define EPU16_BITPACK64_READ_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "types.h"


extern UINT4
Epu16_bitpack64_read_one (Localspace_T oligo, UINT2 *bitpackptrs, UINT2 *bitpackcomp);

extern void
Epu16_bitpack64_block_offsets (UINT4 *offsets, Localspace_T oligo,
			       UINT2 *bitpackptrs, UINT2 *bitpackcomp);

#endif

