/* $Id: sarray-write.h 132688 2014-04-08 17:34:25Z twu $ */
#ifndef SARRAY_WRITE_INCLUDED
#define SARRAY_WRITE_INCLUDED
#include "types.h"
#include "genome.h"

extern void
Sarray_write_array (char *sarrayfile, Genome_T genomecomp, UINT4 genomelength);

extern void
Sarray_write_index_separate (char *indexiptrsfile, char *indexicompfile, char *indexjptrsfile,char *indexjcompfile,
			     char *sarrayfile, Genome_T genomecomp, UINT4 genomelength, bool compressp);

extern void
Sarray_write_index_interleaved (char *indexptrsfile, char *indexcompfile,
				char *sarrayfile, Genome_T genomecomp, UINT4 genomelength, bool compressp);

extern UINT4 *
Sarray_compute_lcp (UINT4 *SA, UINT4 n);

extern unsigned char *
Sarray_discriminating_chars (UINT4 *nbytes, UINT4 *SA, Genome_T genome,
			     unsigned char *lcp_bytes, UINT4 *lcp_guide, UINT4 *lcp_exceptions, int guide_interval,
			     UINT4 n);

extern UINT4 *
Sarray_compute_child (unsigned char *lcp_bytes, UINT4 *lcp_guide, UINT4 *lcp_exceptions, UINT4 n);


extern void
Sarray_array_uncompress (Genome_T genomecomp, char *sarrayfile, char *plcpptrsfile, char *plcpcompfile,
			 UINT4 genomelength, UINT4 start, UINT4 end);

extern void
Sarray_child_uncompress (Genome_T genomecomp, unsigned char *lcpchilddc, UINT4 *lcp_guide, UINT4 *lcp_exceptions,
			 int n_lcp_exceptions, UINT4 *child_guide, UINT4 *child_exceptions, int n_child_exceptions,
			 UINT4 *SA, UINT4 genomelength, UINT4 start, UINT4 end);

extern void
Sarray_child_test (char *sarrayfile,
		   char *childbpfile, char *childs_pagesfile, char *childs_ptrsfile, char *childs_compfile,
		   char *childr_ptrsfile, char *childr_compfile, char *childx_ptrsfile, char *childx_compfile,
		   char *pioneerbpfile, char *pior_ptrsfile, char *pior_compfile,
		   char *piom_ptrsfile, char *piom_compfile, UINT4 genomelength, int sampling_interval);

#endif

