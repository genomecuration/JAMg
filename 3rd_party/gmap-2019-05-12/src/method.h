/* $Id: method.h 218675 2019-03-16 01:25:48Z twu $ */
#ifndef METHOD_INCLUDED
#define METHOD_INCLUDED

#include "filestring.h"


/* Keep in order of single_read */
typedef enum {TR, /* tr, from Kmer_search_transcriptome_{single|paired} */
	      KMER_EXACT, /* both ends, from Kmer_search_genome_ends_exact */
	      EXT, /* ext, from Extension_search */
	      KMER_APPROX, /* both ends, from Kmer_search_genome_ends_approx */
	      SEGMENT1, /* seg1, from Segment_search_genome */
	      DISTANT_RNA, /* distant, from Distant_rna_solve */
	      DISTANT_DNA, /* distant, from Distant_dna_solve */
	      TERMINAL, /* term, from Terminal_solve */
	      SEGMENT2, /* seg2, from Segment_search_genome */

	      KMER_ONE_END, /* one end, from Kmer_search_genome_one_end (3of4 algorithm) */
	      EXT_GMAP, /* ext-gmap, from Extension_search */

	      /* The rest are not currently used */
	      KMER_MIDDLE, /* middle, from Kmer_search_genome_one_kmer */
	      SEGMENT_GMAP, /* seg-gmap, from Path_solve_via_gmap */
	      KMER, /* kmer, from Kmer_search_genome_complete (not currently used) */
	      KMER_GMAP, /* kmer-gmap, from Path_solve_via_gmap (not currently used from KMER) */
	      REGION_GMAP /* not yet implemented */} Method_T;

extern char *
Method_string (Method_T method);
extern void
Method_samprint (Filestring_T fp, Method_T method);
extern void
Method_print (Filestring_T fp, Method_T method);

#endif

