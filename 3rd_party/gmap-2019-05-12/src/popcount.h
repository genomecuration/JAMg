/* $Id: popcount.h 207320 2017-06-14 19:37:19Z twu $ */
#ifndef POPCOUNT_INCLUDED
#define POPCOUNT_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_BUILTIN_CTZ, HAVE_BUILTIN_POPCOUNT, HAVE_BUILTIN_CLZ */
#endif

#if !defined(HAVE_SSE4_2) || (!defined(HAVE_TZCNT) && !defined(HAVE_BUILTIN_CTZ))
extern const int mod_37_bit_position[];
#endif

#if !defined(HAVE_SSE4_2) || (!defined(HAVE_POPCNT) && !defined(HAVE_MM_POPCNT) && !defined(HAVE_BUILTIN_POPCOUNT))
extern const int count_bits[];
#endif

#if !defined(HAVE_SSE4_2) || (!defined(HAVE_LZCNT) && !defined(HAVE_BUILTIN_CLZ))
extern const int clz_table[];
#endif

#endif

