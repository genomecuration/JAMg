/* $Id: epu16-bitpack64-write.h 214305 2018-03-19 23:40:43Z twu $ */
#ifndef EPU16_BITPACK64_WRITE_INCLUDED
#define EPU16_BITPACK64_WRITE_INCLUDED
#include <stdio.h>
#include "types.h"


extern UINT4
Epu16_bitpack64_append_differential (UINT4 *totalcount, FILE *ptrs_fp, FILE *comp_fp,
				     char *packsizes, UINT2 **bitpacks, Localspace_T n);


#endif
