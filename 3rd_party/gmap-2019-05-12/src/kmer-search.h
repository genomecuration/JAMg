/* $Id: kmer-search.h 218473 2019-02-22 23:39:06Z twu $ */
#ifndef KMER_SEARCH_INCLUDED
#define KMER_SEARCH_INCLUDED

#include "mode.h"
#include "list.h"
#include "indexdb.h"
#include "localdb.h"
#include "iit-read-univ.h"
#include "transcriptome.h"
#include "compress.h"
#include "genome.h"

#include "univdiag.h"
#include "stage1hr.h"
#include "stage3hr.h"
#include "listpool.h"
#include "hitlistpool.h"

typedef struct Path_T *Path_T;


extern List_T
Kmer_remap_transcriptome (char *remap_sequence, int remap_seqlength,
			  Chrnum_T chrnum, Chrpos_T lowbound, Chrpos_T highbound,
			  Univ_IIT_T transcript_iit, 
			  Genome_T transcriptomebits, Transcriptome_T transcriptome);

extern void
Kmer_search_transcriptome_single (int *found_score_overall, int *found_score_within_trims,
				  List_T *hits_gplus, List_T *hits_gminus,

				  char *queryuc_ptr, int querylength,
				  Trcoord_T **tplus_stream_array, int *tplus_streamsize_array, int *tplus_diagterm_array,
				  Trcoord_T **tminus_stream_array, int *tminus_streamsize_array, int *tminus_diagterm_array,
				  Compress_T query_compress_fwd, Compress_T query_compress_rev,
				  Univ_IIT_T transcript_iit, Transcriptome_T transcriptome, Genome_T transcriptomebits, 
				  int nmismatches_allowed,
				  Listpool_T listpool, Hitlistpool_T hitlistpool, int level);

#if 0
extern void
Kmer_search_transcriptome_paired (int *found_score_5, int *found_score_3,
				  List_T *hits5_gplus, List_T *hits5_gminus,
				  List_T *hits3_gplus, List_T *hits3_gminus,
				  
				  char *queryuc_ptr_5, int querylength5,
				  char *queryuc_ptr_3, int querylength3,
				  int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,

				  Trcoord_T **tplus_stream_array_5, int *tplus_streamsize_array_5, int *tplus_diagterm_array_5,
				  Trcoord_T **tminus_stream_array_5, int *tminus_streamsize_array_5, int *tminus_diagterm_array_5,
				  Trcoord_T **tplus_stream_array_3, int *tplus_streamsize_array_3, int *tplus_diagterm_array_3,
				  Trcoord_T **tminus_stream_array_3, int *tminus_streamsize_array_3, int *tminus_diagterm_array_3,

				  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
				  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

				  Univ_IIT_T transcript_iit,
				  Transcriptome_T transcriptome, Genome_T transcriptomebits, 
				  int nmismatches_allowed_5, int nmismatches_allowed_3,
				  int kmer_search_sizelimit, Listpool_T listpool, Hitlistpool_T hitlistpool, int level);
#endif

/* Does not take paired_end_p as a parameter.  ? Generates both sense and antisense */
extern void
Kmer_search_genome_ends_exact (bool *abort_exact_p, int *found_score_overall, int *found_score_within_trims,
			       List_T *hits_gplus, List_T *hits_gminus,
			       Stage1_T stage1, int querylength, int *mismatch_positions_alloc,
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       int genestrand, int nmismatches_allowed,
			       Listpool_T listpool, Hitlistpool_T hitlistpool, int level);

/* Takes sizelimit because merging can be expensive */
extern void
Kmer_search_genome_ends_approx (int *found_score_overall, int *found_score_within_trims,
				List_T *hits_gplus, List_T *hits_gminus, Stage1_T stage1,
				Compress_T query_compress_fwd, Compress_T query_compress_rev,
				int querylength, int genestrand, int nmismatches_allowed, int sizelimit,
				Listpool_T listpool, Hitlistpool_T hitlistpool, int level);

#if 0
extern void
Kmer_search_genome_one_end (int *found_score, List_T *hits_gplus, List_T *hits_gminus,

#ifdef LARGE_GENOMES
			    unsigned char *plus_rawpositions_high_5[], unsigned char *minus_rawpositions_high_5[],
			    unsigned char *plus_rawpositions_high_3[], unsigned char *minus_rawpositions_high_3[],
#endif
			    UINT4 *plus_rawpositions_5[], int plus_nrawpositions_5[], int plus_diagterms_5[],
			    UINT4 *minus_rawpositions_5[], int minus_nrawpositions_5[], int minus_diagterms_5[],
			    UINT4 *plus_rawpositions_3[], int plus_nrawpositions_3[], int plus_diagterms_3[],
			    UINT4 *minus_rawpositions_3[], int minus_nrawpositions_3[], int minus_diagterms_3[],

			    Compress_T query_compress_fwd, Compress_T query_compress_rev,
			    char *queryuc_ptr, char *queryrc, int querylength,
			    int nmismatches_allowed, bool paired_end_p, int sizelimit, Listpool_T listpool,
			    int level);
#endif

#if 0
extern void
Kmer_search_genome_one_kmer (int *found_score, List_T *hits_gplus, List_T *hits_gminus,
			     Stage1_T stage1, int *plus_specifici, int *minus_specifici,
			     Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev,

			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     char *queryuc_ptr, char *queryrc, int querylength, int genestrand,
			     int nmismatches_allowed, bool paired_end_p, int level);
#endif

extern void
Kmer_divide_remapped_hits (List_T *transcriptome_hits, List_T *genome_hits, List_T hits,
			   Genome_T genomecomp, Univ_IIT_T transcript_iit, Genome_T transcriptomebits, 
			   Transcriptome_T transcriptome, int genestrand);

extern void
Kmer_search_setup (Mode_T mode_in,
		   int index1part_tr_in, int index1part_in, int index1interval_in, int local1part_in,
		   int min_intronlength_in, int max_deletionlength,
		   Univ_IIT_T chromosome_iit_in, int circular_typeint_in, bool *circularp_in, 
		   Genome_T genomebits_in, Genome_T genomebits_alt_in,
		   Indexdb_T indexdb_in, Indexdb_T indexdb2_in, Indexdb_T indexdb_tr_in, Chrpos_T shortsplicedist_in,
		   bool novelsplicingp_in, Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		   Chrpos_T *splicedists_in, int nsplicesites_in);
#endif
