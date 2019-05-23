/* $Id: segment-search.h 218473 2019-02-22 23:39:06Z twu $ */
#ifndef SEGMENT_SEARCH_INCLUDED
#define SEGMENT_SEARCH_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "bool.h"
#include "types.h"
#include "univcoord.h"
#include "mode.h"
#include "genomicpos.h"
#include "method.h"

#include "shortread.h"
#include "iit-read-univ.h"
#include "genome.h"

#include "stage1hr.h"
#include "resulthr.h"		/* For Pairtype_T */
#include "splice.h"

#ifdef PMAP
#include "oligoindex_pmap.h"
#else
/* #include "oligoindex_hr.h" */
#endif
#include "intlistpool.h"
#ifdef LARGE_GENOMES
#include "uint8listpool.h"
#else
#include "uintlistpool.h"
#endif
#include "listpool.h"
#include "univdiagpool.h"
#include "hitlistpool.h"
#include "record.h"


extern struct Record_T *
Segment_identify_lower (int *nrecords,
#ifdef LARGE_GENOMES
			unsigned char **positions_high, UINT4 **positions_low,
#else
			Univcoord_T **positions,
#endif		  
			int *npositions, bool *validp,
			
			Univcoord_T **stream_alloc, int *streamsize_alloc, int *diagterm_alloc,
			Chrpos_T max_pairlength, int querylength,
			Univcoord_T *ref_diagonals, int ref_ndiagonals);

extern struct Record_T *
Segment_identify_higher (int *nrecords,
#ifdef LARGE_GENOMES
			 unsigned char **positions_high, UINT4 **positions_low,
#else
			 Univcoord_T **positions,
#endif		  
			 int *npositions, bool *validp,

			 Univcoord_T **stream_alloc, int *streamsize_alloc, int *diagterm_alloc,
			 Chrpos_T max_pairlength, int querylength,
			 Univcoord_T *ref_diagonals, int ref_ndiagonals);

extern struct Record_T *
Segment_identify (int *nrecords,
#ifdef LARGE_GENOMES
		  unsigned char **positions_high,
#endif		  
		  UINT4 **positions, int *npositions, bool *validp,

#ifdef LARGE_GENOMES
		  unsigned char **stream_high_alloc, UINT4 **stream_low_alloc,
#else
		  Univcoord_T **stream_alloc,
#endif
		  int *streamsize_alloc, int *diagterm_alloc,
		  int querylength, int sizelimit);

extern void
Segment_search_all (int *found_score_overall, int *found_score_within_trims,
		    List_T *plus_hits, List_T *minus_hits,

		    struct Record_T *plus_records, int plus_nrecords, 
		    struct Record_T *minus_records, int minus_nrecords,

		    char *queryuc_ptr, char *queryrc, int querylength,
		    int *mismatch_positions_alloc, Spliceinfo_T spliceinfo,
		    Univcoord_T **stream_alloc, int *streamsize_alloc,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    int genestrand, bool paired_end_p, bool first_read_p,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
		    Method_T method, int level);

extern void
Segment_search_setup (int index1part_in, int index1interval_in,
		      int max_anchors_in, Univ_IIT_T chromosome_iit_in, int nchromosomes_in,
		      int circular_typeint_in, Mode_T mode_in,
		      Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		      Chrpos_T *splicedists_in, int nsplicesites_in,
		      int max_middle_deletions, Chrpos_T shortsplicedist_in);

extern void
Segment_search_cleanup ();

#endif

