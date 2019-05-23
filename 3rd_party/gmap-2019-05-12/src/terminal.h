/* $Id: terminal.h 218473 2019-02-22 23:39:06Z twu $ */
#ifndef TERMINAL_INCLUDED
#define TERMINAL_INCLUDED

#include "bool.h"
#include "list.h"
#include "compress.h"
#include "genome.h"
#include "iit-read-univ.h"
#include "listpool.h"
#include "hitlistpool.h"

extern List_T
Terminal_solve_plus (int *found_score_overall, int *found_score_within_trims,
		     List_T queryfwd_plus_set, List_T queryrev_plus_set,

		     Compress_T query_compress_fwd, int querylength, 
		     int genestrand, Listpool_T listpool,
		     Hitlistpool_T hitlistpool, int level);

extern List_T
Terminal_solve_minus (int *found_score_overall, int *found_score_within_trims,
		      List_T queryfwd_minus_set, List_T queryrev_minus_set,

		      Compress_T query_compress_rev, int querylength, 
		      int genestrand, Listpool_T listpool,
		      Hitlistpool_T hitlistpool, int level);

extern void
Terminal_setup (Univ_IIT_T chromosome_iit_in, int circular_typeint_in,
		Genome_T genomebits_in, Genome_T genomebits_alt_in,
		int index1part_in, int index1interval_in, int subopt_levels_in);

#endif

