/* $Id: bitpack64-write.h 180341 2015-12-07 18:29:40Z twu $ */
#ifndef BITPACK64_WRITE_INCLUDED
#define BITPACK64_WRITE_INCLUDED
#include <stdio.h>
#include "types.h"

extern int
Bitpack64_write_columnar (FILE *strm_fp, Positionsptr_T *strm_buffer, int strm_buffer_size, int strm_buffer_i,
			  const UINT4 *horizontal, int packsize);
extern int
Bitpack64_compute_q4_diffs_bidir (UINT4 *diffs, UINT4 *values);

/* Stores the $(n+1)$ values [0..n] */
extern void
Bitpack64_write_differential (char *ptrsfile, char *compfile, UINT4 *ascending, Oligospace_T n);
extern void
Bitpack64_write_differential_bitpacks (char *ptrsfile, char *compfile, char *packsizes, UINT4 **bitpacks,
				       Oligospace_T n);
extern void
Bitpack64_write_differential_huge (char *pagesfile, char *ptrsfile, char *compfile,
				   UINT8 *ascending, Oligospace_T n);
extern void
Bitpack64_write_differential_huge_bitpacks (char *pagesfile, char *ptrsfile, char *compfile,
					    char *packsizes, UINT4 **bitpacks, Oligospace_T n);

/* Stores the $n$ values [0..(n-1)] */
extern void
Bitpack64_write_direct (char *ptrsfile, char *compfile, UINT4 *direct, Oligospace_T n);

#endif
