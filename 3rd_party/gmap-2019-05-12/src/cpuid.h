/* $Id: cpuid.h 219226 2019-05-12 22:32:05Z twu $ */
#ifndef CPUID_INCLUDED
#define CPUID_INCLUDED
#include "bool.h"

extern void
CPUID_support (bool *sse2_support_p, bool *ssse3_support_p, bool *sse4_1_support_p, bool *sse4_2_support_p,
	       bool *avx2_support_p, bool *avx512_support_p, bool *avx512bw_support_p);

#endif

