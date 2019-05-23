/* $Id: extension-search.h 218736 2019-03-25 16:53:44Z twu $ */
#ifndef EXTENSION_SEARCH_INCLUDED
#define EXTENSION_SEARCH_INCLUDED

#include "types.h"
#include "mode.h"
#include "list.h"

#include "iit-read-univ.h"
#include "genome.h"
#include "compress.h"

#include "stage1hr.h"
#include "indexdb.h"
#include "intlistpool.h"

#ifdef LARGE_GENOMES
#include "uint8listpool.h"
#else
#include "uintlistpool.h"
#endif

#include "listpool.h"
#include "univdiagpool.h"
#include "hitlistpool.h"


#define T Elt_T
typedef struct T *T;
struct T {
  int qstart;
  int qend;

  int nmatches;

  Univcoord_T *all_diagonals;
  int n_all_diagonals;

  Univcoord_T *diagonals;
  int ndiagonals;

  int lowi;
  int highi;
};


extern void
Extension_search_setup (Mode_T mode, Univ_IIT_T chromosome_iit_in, int circular_typeint_in, bool *circularp_in,
			Genome_T genomebits_in, Genome_T genomebits_alt_in, Indexdb_T indexdb_in, Indexdb_T indexdb2_in,
			Chrpos_T shortsplicedist_in, Chrpos_T shortsplicedist_novelend,
			int min_intronlength_in, int max_deletionlength, int max_end_deletions_in,
			int max_middle_insertions_in, int max_end_insertions,
			int index1part_in, int local1part_in);

extern void
Elt_gc (List_T *set);

extern void
Elt_set_queryfwd_extend_left (List_T list,
			      Compress_T query_compress, bool plusp, int genestrand);
extern void
Elt_set_queryrev_extend_right (List_T list, int querylength,
			       Compress_T query_compress, bool plusp, int genestrand);

extern void
Extension_search (int *found_score_overall, int *found_score_within_trims,
		  List_T *hits_gplus, List_T *hits_gminus, Stage1_T stage1,
		  char *queryuc_ptr, char *queryrc, int querylength,
		  Compress_T query_compress_fwd, Compress_T query_compress_rev, 
#if 0
		  Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		  Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		  Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
		  int nmismatches_allowed, int genestrand, bool paired_end_p, bool first_read_p,
		  int level);

#undef T
#endif


