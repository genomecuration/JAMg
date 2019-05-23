/* $Id: distant-rna.h 218675 2019-03-16 01:25:48Z twu $ */
#ifndef DISTANT_RNA_INCLUDED
#define DISTANT_RNA_INCLUDED

#include "bool.h"
#include "list.h"
#include "univdiag.h"
#include "compress.h"
#include "genome.h"
#include "iit-read-univ.h"
#include "listpool.h"
#include "hitlistpool.h"


extern void
Distant_rna_solve (int *found_score_overall, int *found_score_within_trims,
		   List_T *hits_plus, List_T *hits_minus,

		   List_T *startfrags_plus, List_T *endfrags_plus,
		   List_T *startfrags_minus, List_T *endfrags_minus,

		   List_T queryfwd_plus_set, List_T queryfwd_minus_set,
		   List_T queryrev_plus_set, List_T queryrev_minus_set,

		   int *mismatch_positions_alloc, int *positions_alloc,
		   Compress_T query_compress_fwd, Compress_T query_compress_rev,
		   char *queryuc_ptr, char *queryrc, int querylength, int max_splice_mismatches,
		   int genestrand, bool first_read_p, Listpool_T listpool,
		   Hitlistpool_T hitlistpool, int level);

extern void
Distant_rna_setup (Univ_IIT_T chromosome_iit_in, int circular_typeint_in,
		   Genome_T genomebits_in, Genome_T genomebits_alt_in,
		   int index1part_in, int index1interval_in,
		   int subopt_levels_in, Chrpos_T shortsplicedist_in,
		   int min_distantsplicing_end_matches_in, int min_distantsplicing_identity_in,
		   int localsplicing_penalty_in, int distantsplicing_penalty_in);

#endif

