/* $Id: genomicpos.h 214305 2018-03-19 23:40:43Z twu $ */
#ifndef GENOMICPOS_INCLUDED
#define GENOMICPOS_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include <stdlib.h>
#include "types.h"

/* A genomic position */
#ifdef LARGE_GENOMES
#include "uint8list.h"
typedef Uint8list_T Genomicposlist_T;
#else
#include "uintlist.h"
typedef Uintlist_T Genomicposlist_T;
#endif

/* A chromosomal position */
typedef UINT4 Chrpos_T;

extern char *
Genomicpos_commafmt (
#ifdef HAVE_64_BIT
		     UINT8 N
#else
		     UINT4 N
#endif
		     );
#ifdef MEMUSAGE
void
Genomicpos_commafmt_fill (char *string,
#ifdef HAVE_64_BIT
		     UINT8 N
#else
		     UINT4 N
#endif
		     );
#endif

extern int
UINT8_compare (const void *a, const void *b);
extern int
UINT4_compare (const void *a, const void *b);
extern int
UINT2_compare (const void *a, const void *b);
extern int
Univcoord_compare (const void *a, const void *b);
extern int
Chrpos_compare (const void *a, const void *b);

#endif
