/* $Id: sarray-write.h 115133 2013-11-15 05:14:24Z twu $ */
#ifndef SARRAY_WRITE_INCLUDED
#define SARRAY_WRITE_INCLUDED
#include "types.h"
#include "genome.h"

extern void
Sarray_write_array (char *sarrayfile, Genome_T genomecomp, UINT4 genomelength);
extern void
Sarray_write_lcp (char *lcpptrsfile, char *lcpcompfile, char *sarrayfile,
		  Genome_T genomecomp, UINT4 genomelength);
extern void
Sarray_write_index (char *saindexfile, char *sarrayfile, char *lcpptrsfile, char *lcpcompfile,
		    Genome_T genomecomp, UINT4 genomelength);
extern void
Sarray_lcp_uncompress (char *lcpptrsfile, char *lcpcompfile, char *sarrayfile, UINT4 genomelength);

#endif

