/* $Id: stage1hr.h 218320 2019-01-31 19:56:35Z twu $ */
#ifndef STAGE1HR_INCLUDED
#define STAGE1HR_INCLUDED

typedef struct Stage1_T *Stage1_T;

#include "bool.h"
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "localdb.h"
#include "shortread.h"
#include "iit-read-univ.h"
#include "genome.h"
#include "resulthr.h"		/* For Pairtype_T */
#include "stage3hr.h"

#ifdef LARGE_GENOMES
#include "uint8listpool.h"
#else
#include "uintlistpool.h"
#endif

#include "intlistpool.h"
#include "listpool.h"
#include "univdiagpool.h"
#include "hitlistpool.h"
#include "splice.h"


#define T Stage1_T
struct T {
  /* Initialized by Stage1_fill_all_positions */
  bool *plus_validp;
  bool *minus_validp;

  Oligospace_T *plus_oligos;
  Oligospace_T *minus_oligos;

  bool *retrievedp_allocated;
  bool *plus_retrievedp;	/* points to above[index1interval-1] */
  bool *minus_retrievedp;	/* points to above[index1interval-1] */

#ifdef LARGE_GENOMES
  unsigned char **positions_high_allocated;
  unsigned char **plus_positions_high; /* points to above[index1interval-1] */
  unsigned char **minus_positions_high; /* points to above[index1interval-1] */
#endif

  UINT4 **positions_allocated;
  UINT4 **plus_positions; /* points to above[index1interval-1] */
  UINT4 **minus_positions; /* points to above[index1interval-1] */

  int *npositions_allocated;
  int *plus_npositions;		/* points to above[index1interval-1] */
  int *minus_npositions;	/* points to above[index1interval-1] */


  /* Memory allocated for mismatch_positions in Substring_new (Genome_mismatches_{left,right}),
     and for Distant_rna_solve */
  int *mismatch_positions_alloc;
  int *positions_alloc;

  /* Memory allocated for spliceinfo in kmer-search.c and
     path-solve.c, used by Splice_resolve_sense and
     Splice_resolve_antisense */
  Spliceinfo_T spliceinfo;

  /* Memory allocated for Segment_identify in segment-search.c, and
     Merge_diagonals in kmer-search.c (which needs four sets of
     arrays) */
#ifdef LARGE_GENOMES
  unsigned char **stream_high_alloc, **gplus_stream_high_array_5, **gminus_stream_high_array_5, **gplus_stream_high_array_3, **gminus_stream_high_array_3;
  UINT4 **stream_low_alloc, **gplus_stream_low_array_5, **gminus_stream_low_array_5, **gplus_stream_low_array_3, **gminus_stream_low_array_3;
#endif
  Univcoord_T **stream_alloc, **gplus_stream_array_5, **gminus_stream_array_5, **gplus_stream_array_3, **gminus_stream_array_3;
  Trcoord_T **tplus_stream_array, **tminus_stream_array;

  int *streamsize_alloc, *tplus_streamsize_array, *tminus_streamsize_array,
    *gplus_streamsize_array_5, *gminus_streamsize_array_5, *gplus_streamsize_array_3, *gminus_streamsize_array_3;
  int *querypos_diagterm_alloc, *tplus_diagterm_array, *tminus_diagterm_array,
    *gplus_diagterm_array_5, *gminus_diagterm_array_5, *gplus_diagterm_array_3, *gminus_diagterm_array_3;


  /* Intermediate calculations.  Need to save and resume from single_read to paired_read */
  /* Returned from Kmer_search_genome_exact and used to create stage1 for Kmer_search_genome_complete */
#ifdef LARGE_GENOMES
  unsigned char *plus_rawpositions_high_5[6];
  unsigned char *minus_rawpositions_high_5[6];
  unsigned char *plus_rawpositions_high_3[6];
  unsigned char *minus_rawpositions_high_3[6];
#endif
  UINT4 *plus_rawpositions_5[6];
  UINT4 *minus_rawpositions_5[6];
  UINT4 *plus_rawpositions_3[6];
  UINT4 *minus_rawpositions_3[6];

  int plus_nrawpositions_5[6];
  int minus_nrawpositions_5[6];
  int plus_nrawpositions_3[6];
  int minus_nrawpositions_3[6];

  int plus_diagterms_5[6];
  int minus_diagterms_5[6];
  int plus_diagterms_3[6];
  int minus_diagterms_3[6];


  /* Returned from Extension_search and used for distant splicing */
  List_T queryfwd_plus_set;
  List_T queryfwd_minus_set;
  List_T queryrev_plus_set;
  List_T queryrev_minus_set;
};


extern Stage3end_T *
Stage1_single_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    Shortread_T queryseq,
#if 0
		    Oligoindex_array_T oligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
#endif
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool);

extern Stage3pair_T *
Stage1_paired_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
		    Stage3end_T **stage3array5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
		    Stage3end_T **stage3array3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
		    Shortread_T queryseq5, Shortread_T queryseq3, Chrpos_T pairmax_linear,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool);

extern void
Stage1hr_setup (Univ_IIT_T transcript_iit_in, Transcriptome_T transcriptome_in, Genome_T transcriptomebits_in,
		bool use_only_transcriptome_p_in,

		Indexdb_T indexdb_fwd_in, Indexdb_T indexdb_rev_in, Indexdb_T indexdb_tr_in, Localdb_T localdb_in,

		int index1part_tr_in, int index1part_in, int index1interval_in, 
		double user_maxlevel_float_in, double user_mincoverage_float_in,

		Univ_IIT_T chromosome_iit_in, int nchromosomes_in,
		Genome_T genomecomp_in, Genome_T genomebits_in, Genome_T genomebits_alt_in,
		Mode_T mode_in, int maxpaths_search_in, int maxpaths_report_in,

		bool find_dna_chimeras_p_in, bool distances_observed_p_in, int subopt_levels_in,
		int max_middle_insertions, int max_middle_deletions,

		bool novelsplicingp_in, Chrpos_T shortsplicedist_in,
		Chrpos_T min_intronlength_in, Chrpos_T expected_pairlength_in, Chrpos_T pairlength_deviation_in,

		int distantsplicing_penalty_in,
		int min_distantsplicing_end_matches_in, int min_distantsplicing_identity_in);

#undef T
#endif

