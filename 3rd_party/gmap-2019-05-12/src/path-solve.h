/* $Id: path-solve.h 218692 2019-03-19 17:38:52Z twu $ */
#ifndef PATH_INCLUDED
#define PATH_INCLUDED

#include "types.h"
#include "univcoord.h"
#include "bool.h"
#include "list.h"
#include "univdiag.h"
#include "chrnum.h"
#include "compress.h"
#include "genome.h"
#include "localdb.h"
#include "stage3hr.h"
#include "splice.h"
#include "intlistpool.h"

#ifdef LARGE_GENOMES
#include "uint8listpool.h"
#else
#include "uintlistpool.h"
#endif

#include "listpool.h"
#include "hitlistpool.h"
#include "univdiagpool.h"


/* Used by extension-search and segment-search */
extern List_T
Path_solve_from_diagonals (bool *foundp, int *found_score_overall, int *found_score_within_trims, List_T hits,
			   Univcoord_T middle_diagonal_univdiagonal, int middle_diagonal_qstart, int middle_diagonal_qend,
			   List_T qend_diagonals, List_T qstart_diagonals, char *queryptr, int querylength,
			   int *mismatch_positions_alloc, Spliceinfo_T spliceinfo,
			   Univcoord_T **stream_alloc, int *streamsize_alloc, Compress_T query_compress, 
			   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			   bool plusp, int genestrand, int nmismatches_allowed, bool paired_end_p, bool first_read_p,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
			   Method_T method, int level);


#if 0
/* Used previously by segment-search.c */
List_T
Path_solve_from_path (bool *foundp, int *found_score, List_T hits, List_T complete_path,
		      Compress_T query_compress, int querylength, Genome_T genomebits,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      bool plusp, int genestrand, int max_mismatches_allowed, Method_T method, int level);
#endif


extern void
Path_solve_setup (Genome_T genomebits_in, Genome_T genomebits_alt_in, bool *circularp_in,
		  Localdb_T localdb_in, Localdb_T localdb2_in,
		  Chrpos_T shortsplicedist_novelend, int min_intronlength_in, int max_deletionlength,
		  int max_insertionlen_in, 
		  bool novelsplicingp_in, Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		  Chrpos_T *splicedists_in, int nsplicesites_in, 
		  int index1part_in, int index1interval_in, int local1part_in);

#endif

