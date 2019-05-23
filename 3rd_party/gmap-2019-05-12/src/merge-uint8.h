/* $Id: merge-uint8.h 216888 2018-10-07 16:37:57Z twu $ */
#ifndef MERGE_UINT8_INCLUDED
#define MERGE_UINT8_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "types.h"

#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
extern UINT8 *
Merge_uint8 (UINT8 *__restrict__ dest, UINT8 *__restrict__ A,
	     UINT8 *__restrict__ B, int nA, int nB);
#endif

#endif
