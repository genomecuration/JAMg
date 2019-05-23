/* $Id: concordance.h 218473 2019-02-22 23:39:06Z twu $ */
#ifndef CONCORDANCE_INCLUDED
#define CONCORDANCE_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "bool.h"
#include "list.h"
#include "types.h"
#include "genomicpos.h"
#include "ladder.h"
#include "stage3hr.h"
#include "record.h"

#include "listpool.h"
#include "hitlistpool.h"

#define T Stage3end_T


extern List_T
Concordance_pair_up_transcriptome (bool *abort_pairing_p, int *concordant_score_overall, List_T hitpairs,

				   Ladder_T ladder5_plus, Ladder_T ladder5_minus,
				   Ladder_T ladder3_plus, Ladder_T ladder3_minus,
#if 0
				   char *queryuc_ptr_5, char *queryuc_ptr_3,
				   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
				   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
				   Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				   Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
				   Listpool_T listpool, Hitlistpool_T hitlistpool,
				   int maxpairedpaths, int genestrand);

extern List_T
Concordance_pair_up_genome (bool *abort_pairing_p, int *adjacent_score, int *concordant_score_overall,
			    List_T *distant_hitpairs, List_T hitpairs, 

			    List_T hitlist5_gplus, List_T hitlist5_gminus,
			    List_T hitlist3_gplus, List_T hitlist3_gminus,
			    
			    Ladder_T ladder5_plus, Ladder_T ladder5_minus,
			    Ladder_T ladder3_plus, Ladder_T ladder3_minus,
			    int querylength5, int querylength3,

#if 0
			    char *queryuc_ptr_5, char *queryuc_ptr_3,
			    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			    Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
			    Listpool_T listpool, Hitlistpool_T hitlistpool,
			    int maxpairedpaths, int genestrand);

extern List_T
Concordance_pair_up_distant (bool *abort_pairing_p, int *concordant_score_overall,
			     List_T *samechr, List_T *conc_transloc, List_T hitpairs, 

			     List_T hitlist5_gplus, List_T hitlist5_gminus,
			     List_T hitlist3_gplus, List_T hitlist3_gminus,

			     Ladder_T ladder5_plus, Ladder_T ladder5_minus,
			     Ladder_T ladder3_plus, Ladder_T ladder3_minus,
			     int querylength5, int querylength3,
#if 0
			     char *queryuc_ptr_5, char *queryuc_ptr_3,
			     Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			     Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			     Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
			     Listpool_T listpool, Hitlistpool_T hitlistpool,
			     int maxpairedpaths, int genestrand);

extern void
Concordance_setup (int subopt_levels_in, bool novelsplicingp_in,
		   Chrpos_T pairmax_transcriptome_in,
		   Chrpos_T pairmax_linear_in, Chrpos_T pairmax_circular_in,
		   Chrpos_T expected_pairlength_in, Chrpos_T pairlength_deviation_in, 
		   bool *circularp_in, bool merge_samechr_p_in);

#undef T
#endif

