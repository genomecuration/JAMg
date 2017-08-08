/* $Id: sarray-search.h 207324 2017-06-14 19:41:18Z twu $ */
#ifndef SARRAY_SEARCH_INCLUDED
#define SARRAY_SEARCH_INCLUDED
#include "access.h"
#include "bool.h"
#include "mode.h"
#include "genome.h"
#include "compress.h"
#include "genomicpos.h"
#include "splicetrie.h"
#include "iit-read-univ.h"
#include "sarray-read.h"


#define T Sarray_T

extern void
Sarray_search_setup (T sarray_fwd_in, T sarray_rev_in, Genome_T genome_in, Mode_T mode,
		     Univ_IIT_T chromosome_iit_in, int circular_typeint_in, bool *circularp_in,
		     Chrpos_T shortsplicedist_in, int splicing_penalty_in,
		     int min_intronlength_in, int max_deletionlength, int max_end_deletions,
		     int max_middle_insertions_in, int max_end_insertions,
		     Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		     Chrpos_T *splicedists_in, int nsplicesites_in);

extern List_T
Sarray_search_greedy (int *found_score, char *queryuc_ptr, char *queryrc, int querylength,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      int nmisses_allowed, int genestrand);

#undef T
#endif
