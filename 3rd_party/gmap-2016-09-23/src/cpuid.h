/* $Id: cpuid.h 171614 2015-08-10 23:27:29Z twu $ */
#ifndef CPUID_INCLUDED
#define CPUID_INCLUDED
#include "bool.h"

extern void
CPUID_support (bool *sse2_support_p, bool *ssse3_support_p, bool *sse4_1_support_p, bool *sse4_2_support_p, bool *avx2_support_p);

#endif

