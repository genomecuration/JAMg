#ifndef BITPACK64_READ_INCLUDED
#define BITPACK64_READ_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "types.h"

/* For reading differential-coded bitstreams */

extern UINT4
Bitpack64_read_one (Oligospace_T oligo,
		    UINT4 *bitpackptrs, UINT4 *bitpackcomp);

extern UINT8
Bitpack64_read_one_huge (Oligospace_T oligo, UINT4 *bitpackpages,
			 UINT4 *bitpackptrs, UINT4 *bitpackcomp);

extern void
Bitpack64_block_offsets (UINT4 *offsets, Oligospace_T oligo,
			 UINT4 *bitpackptrs, UINT4 *bitpackcomp);


#if defined(HAVE_64_BIT) && (defined(UTILITYP) || defined(LARGE_GENOMES))
extern void
Bitpack64_block_offsets_huge (UINT8 *offsets, Oligospace_T oligo,
			      UINT4 *bitpackpages, UINT4 *bitpackptrs, UINT4 *bitpackcomp);
#endif

#endif
