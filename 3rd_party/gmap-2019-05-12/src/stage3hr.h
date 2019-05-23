/* $Id: stage3hr.h 218675 2019-03-16 01:25:48Z twu $ */
#ifndef STAGE3HR_INCLUDED
#define STAGE3HR_INCLUDED

typedef struct Stage3end_T *Stage3end_T;
typedef struct Stage3pair_T *Stage3pair_T;

/* Should arrange in order of goodness, best to worst */
typedef enum {EXACT, SUB, SUBSTRINGS,
	      HALFSPLICE_DONOR, HALFSPLICE_ACCEPTOR, SPLICE, SAMECHR_SPLICE, TRANSLOC_SPLICE} Hittype_T;

#include <stdio.h>
#include "bool.h"
#include "sense.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "list.h"
#include "intlist.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "shortread.h"
#include "transcriptome.h"
#include "genome.h"
#include "compress.h"
#include "resulthr.h"
#include "substring.h"
#include "junction.h"
#include "method.h"
#include "simplepair.h"
#include "filestring.h"
#include "output.h"

#include "listpool.h"
#include "hitlistpool.h"


#define T Stage3end_T


extern void
Stage3hr_setup (bool transcriptomep_in, bool invert_first_p_in, bool invert_second_p_in,
		Genome_T genomecomp_in, Genome_T genomebits_in, Genome_T genomebits_alt_in,
		Univ_IIT_T chromosome_iit_in, int nchromosomes_in, int circular_typeint_in,

		IIT_T genes_iit_in, int *genes_divint_crosstable_in, int *genes_chrnum_crosstable_in, 
		Genome_T transcriptomebits_in, Transcriptome_T transcriptome_in, Univ_IIT_T transcript_iit_in,

		IIT_T tally_iit_in, int *tally_divint_crosstable_in,
		IIT_T runlength_iit_in, int *runlength_divint_crosstable_in,
		bool distances_observed_p,
		Chrpos_T pairmax_linear_in, Chrpos_T pairmax_circular_in,
		Chrpos_T expected_pairlength_in, Chrpos_T pairlength_deviation_in,
		int localsplicing_penalty_in, int indel_penalty_middle_in,
		int antistranded_penalty_in, bool favor_multiexon_p_in, int subopt_levels_in,
		int max_middle_insertions_in, int max_middle_deletions_in,
		Chrpos_T shortsplicedist_in, bool *circularp_in, bool *altlocp_in,
		Univcoord_T *alias_starts_in, Univcoord_T *alias_ends_in,
		bool ignore_trim_p, bool omit_concordant_uniq_p_in, bool omit_concordant_mult_p_in,
		char *failedinput_root_in, Outputtype_T output_type, bool merge_samechr_p_in,
		bool method_print_p_in, bool want_random_p_in);

extern char *
Stage3end_deletion_string (T this);

extern Hittype_T
Stage3end_hittype (T this);
extern char *
Stage3end_hittype_string (T this);
extern Method_T
Stage3end_method (T this);
extern int
Stage3end_genestrand (T this);
extern bool
Stage3end_transcriptomep (T this);
extern List_T
Stage3end_transcripts (T this);
extern void
Stage3end_set_transcripts (T this, List_T transcripts);
extern List_T
Stage3end_transcripts_other (T this);

extern void
Stage3end_transfer_transcripts (T dest, List_T sources);
extern void
Stage3end_transfer_transcripts_transloc (T dest, List_T sources);
extern bool
Stage3end_improved_by_gmap_p (T this);
extern void
Stage3end_set_improved_by_gmap (T this);
extern bool
Stage3end_distant_splice_p (T this);
extern Chrnum_T
Stage3end_chrnum (T this);
extern Chrnum_T
Stage3end_effective_chrnum (T this);
extern Chrnum_T
Stage3end_other_chrnum (T this);
extern Univcoord_T
Stage3end_chroffset (T this);
extern Univcoord_T
Stage3end_chrhigh (T this);
extern Chrpos_T
Stage3end_chrlength (T this);
extern Chrpos_T
Stage3end_chrpos_low (T this);
extern Chrpos_T
Stage3end_chrpos_high (T this);

extern Univcoord_T
Stage3end_genomicstart (T this);
extern Univcoord_T
Stage3end_genomicend (T this);
extern int
Stage3end_query_alignment_length (T this);
extern Chrpos_T
Stage3end_genomic_alignment_length (T this);
extern int
Stage3end_mapq_score (T this);
extern int
Stage3end_absmq_score (T this);
extern int
Stage3end_gmap_goodness (T this);
extern int
Stage3end_gmap_max_match_length (T this);
extern int
Stage3end_nmismatches_whole (T this);
extern int
Stage3end_nmismatches_bothdiff (T this);
extern int
Stage3end_nmismatches_refdiff (T this);
extern Endtype_T
Stage3end_start_endtype (T this);
extern Endtype_T
Stage3end_end_endtype (T this);
extern Endtype_T
Stage3end_gmap_start_endtype (T this);
extern Endtype_T
Stage3end_gmap_end_endtype (T this);
extern int
Stage3end_nindels (T this);
extern int
Stage3end_querylength (T this);
extern bool
Stage3end_plusp (T this);
extern bool
Stage3end_paired_usedp (T this);
extern int
Stage3end_max_trim (T this);
extern int
Stage3end_circularpos (T this);

extern Junction_T
Stage3end_junctionD (T this);
extern Junction_T
Stage3end_junctionA (T this);

extern List_T
Stage3end_substrings_LtoH (T this);
extern List_T
Stage3end_junctions_LtoH (T this);

extern Substring_T
Stage3end_substring1 (T this);
extern Substring_T
Stage3end_substringN (T this);
extern Substring_T
Stage3end_substring_for_concordance (T this, bool first_read_p);
extern Substring_T
Stage3end_substring_other (T this, bool first_read_p);
extern bool
Stage3end_donor_concordant_p (T this, bool first_read_p);

extern bool
Stage3end_donor_concordant_p (T this, bool first_read_p);
extern Substring_T
Stage3end_substring_donor (T this);
extern Substring_T
Stage3end_substring_acceptor (T this);
extern Substring_T
Stage3end_substringD (T this);
extern Substring_T
Stage3end_substringA (T this);
extern Substring_T
Stage3end_substringS (T this);

extern Substring_T
Stage3end_substring_low (T this, int hardclip_low);
extern Substring_T
Stage3end_substring_containing (T this, int querypos);
extern double
Stage3end_min_evalue (T this);
extern struct Simplepair_T *
Stage3end_pairarray (T this);
extern int
Stage3end_npairs (T this);
extern List_T
Stage3end_cigar_tokens (T this);
extern bool
Stage3end_gmap_intronp (T this);

extern Chrpos_T
Stage3end_distance (T this);
extern Chrpos_T
Stage3end_shortexonA_distance (T this);
extern Chrpos_T
Stage3end_shortexonD_distance (T this);
extern double
Stage3end_chimera_prob (T this);
extern double
Stage3end_shortexon_prob (T this);
extern Univcoord_T
Stage3end_chimera_segmenti_left (T this);
extern Univcoord_T
Stage3end_chimera_segmentj_left (T this);
extern int
Stage3end_chimera_segmenti_cmp (const void *a, const void *b);
extern int
Stage3end_chimera_segmentj_cmp (const void *a, const void *b);
extern int
Stage3end_shortexon_substringD_cmp (const void *a, const void *b);
extern int
Stage3end_shortexon_substringA_cmp (const void *a, const void *b);

extern int
Stage3end_sensedir (T this);
extern int
Stage3end_sensedir_distant_guess (T this);
extern int
Stage3end_cdna_direction (T this);

extern bool
Stage3end_start_has_alts_p (T this);
extern bool
Stage3end_end_has_alts_p (T this);
extern Univcoord_T *
Stage3end_start_alts_coords (T this);
extern Univcoord_T *
Stage3end_end_alts_coords (T this);
extern int
Stage3end_start_alts_ncoords (T this);
extern int
Stage3end_end_alts_ncoords (T this);


extern int
Stage3end_substrings_querystart (T this);
extern int
Stage3end_substrings_queryend (T this);
extern int
Stage3end_gmap_querystart (T this);
extern int
Stage3end_gmap_queryend (T this);
extern void
Stage3end_count_hits (int *npaths_primary, int *npaths_altloc, List_T hits);

extern void
Stage3end_free (T *old);
extern void
Stage3end_gc (List_T values);


extern bool
Stage3pair_distant_splice_p (Stage3pair_T this);
extern int
Stage3pair_genestrand (Stage3pair_T this);
extern Stage3end_T
Stage3pair_hit5 (Stage3pair_T this);
extern Stage3end_T
Stage3pair_hit3 (Stage3pair_T this);
extern int
Stage3pair_mapq_score (Stage3pair_T this);
extern int
Stage3pair_absmq_score (Stage3pair_T this);

extern List_T
Stage3pair_transcripts5 (Stage3pair_T this);
extern List_T
Stage3pair_transcripts3 (Stage3pair_T this);


extern Chrpos_T
Stage3pair_pairlength (Stage3pair_T this);
extern int
Stage3pair_relationship (Stage3pair_T this);
extern int
Stage3pair_total_trim (Stage3pair_T this);
extern int
Stage3pair_nmatches_to_trims (int *nmatches5, int *nmatches3, Stage3pair_T this);
extern bool
Stage3pair_concordantp (List_T hitpairs);
extern void
Stage3pair_count_hits (int *npaths_primary, int *npaths_altloc, List_T hitpairs);

extern List_T
Stage3pair_filter_nonconcordant (List_T hitpairs, Hitlistpool_T hitlistpool);
extern int
Stage3pair_overlap (int *hardclip5_low, int *hardclip5_high, int *hardclip3_low, int *hardclip3_high, Stage3pair_T this);

extern void
Stage3pair_free (Stage3pair_T *old);

extern char *
Stage3end_substrings_genomic_sequence (int *seqlength, T this, Genome_T genome);


extern T
Stage3end_new_terminal (int *found_score_overall, int *found_score_within_trims,
			Substring_T substring_in, int querylength,
			bool gplusp, int genestrand, int sensedir, Listpool_T listpool,
			Method_T method, int level);

extern T
Stage3end_new_precomputed (int *found_score_overall, int *found_score_within_trims, int nmismatches_bothdiff,
			   List_T substrings, List_T junctions, List_T transcripts, List_T transcripts_other,
			   int querylength, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			   bool gplusp, int genestrand, int sensedir, Listpool_T listpool, Method_T method, int level);

extern int
Stage3end_nmatches_substrings (Intlist_T endpoints, Univcoordlist_T lefts,
			       Intlist_T nmismatches_list, List_T junctions,
			       int querylength, Compress_T query_compress,
			       Substring_T qend_alts, Substring_T qstart_alts,
			       bool plusp, int genestrand,
			       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			       bool splice5p_in, bool splice3p_in, Listpool_T listpool);
extern T
Stage3end_new_substrings (int *found_score_overall, int *found_score_within_trims,
			  Intlist_T endpoints, Univcoordlist_T lefts,
			  Intlist_T nmismatches_list, List_T junctions,
			  int querylength, Compress_T query_compress,
			  Substring_T qend_alts, Substring_T qstart_alts,
			  bool plusp, int genestrand, int sensedir,
			  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			  bool splice5p_in, Splicetype_T splicetype5, double ambig_prob_5,
			  bool splice3p_in, Splicetype_T splicetype3, double ambig_prob_3,
			  Listpool_T listpool, Method_T method, int level);

extern T
Stage3end_new_substitution (int *found_score_overall, int *found_score_within_trims,
			    Univcoord_T left, int genomiclength,
			    int querylength, int *mismatch_positions_alloc, Compress_T query_compress,
			    bool plusp, int genestrand, int sensedir, int nmismatches_allowed,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			    Listpool_T listpool, Method_T method, int level);
extern T
Stage3end_new_splice (int *found_score_overall, int *found_score_within_trims,
		      Substring_T donor, Substring_T acceptor, 
		      Chrpos_T distance, bool shortdistancep, int querylength,
		      bool copy_donor_p, bool copy_acceptor_p,
		      bool first_read_p, int sensedir, Listpool_T listpool, Method_T method, int level);

extern T
Stage3end_new_distant (int *found_score_overall, int *found_score_within_trims,
		       Substring_T startfrag, Substring_T endfrag, int splice_pos,
		       int nmismatches1, int nmismatches2,
		       double prob1, double prob2, int sensedir_distant_guess,
		       Chrpos_T distance, bool shortdistancep, int querylength,
		       bool first_read_p, Listpool_T listpool, int level);

extern List_T
Stage3end_sort_bymatches (List_T hits);
extern List_T
Stage3end_sort_by_paired_seenp (List_T hits);

extern List_T
Stage3end_filter (List_T hits, Hitlistpool_T hitlistpool, int max_mismatches, int min_coverage);

extern Stage3end_T *
Stage3end_eval_and_sort (int npaths, int *first_absmq, int *second_absmq,
			 Stage3end_T *stage3array, char *queryuc_ptr,
			 char *quality_string, bool displayp);
extern Stage3end_T *
Stage3end_eval_and_sort_guided (int npaths, int *first_absmq, int *second_absmq, Stage3end_T guide,
				Stage3end_T *stage3array, char *queryuc_ptr,
				char *quality_string, bool displayp);
extern List_T
Stage3end_optimal_score (List_T hitlist, Hitlistpool_T hitlistpool, int querylength, bool finalp);
extern bool
Stage3pair_sense_consistent_p (List_T hitpairlist);
extern List_T
Stage3end_linearize_5 (List_T hitlist);
extern List_T
Stage3end_linearize_3 (List_T hitlist);
extern List_T
Stage3end_remove_circular_alias (List_T hitlist, Hitlistpool_T hitlistpool);
extern List_T
Stage3end_remove_duplicates (List_T hitlist, Hitlistpool_T hitlistpool);
extern T *
Stage3end_remove_duplicates_array (int *nunique, List_T *duplicates, T *hits, int n,
				   Hitlistpool_T hitlistpool);
extern List_T
Stage3end_merge_ambiguous (List_T hitlist);
extern List_T
Stage3end_merge_trim5 (List_T hitlist, Compress_T query_compress_fwd, Compress_T query_compress_rev);
extern List_T
Stage3end_merge_trim3 (List_T hitlist, Compress_T query_compress_fwd, Compress_T query_compress_rev);
extern int
Stage3end_hit_goodness_cmp (bool *equalp, Stage3end_T hit,
			    Stage3end_T best_hit, bool finalp);


extern List_T
Stage3end_remove_overlaps (List_T hitlist, Hitlistpool_T hitlistpool, bool finalp);
extern List_T
Stage3end_resolve_multimapping (List_T hitlist, Hitlistpool_T hitlistpool);
extern void
Stage3end_method_samprint (Filestring_T fp, T this);
extern void
Stage3end_method_print (Filestring_T fp, Method_T method);
extern Pairtype_T
Stage3_determine_pairtype (T hit5, T hit3, Stage3pair_T stage3pair);



/* If hit5 and hit3 are not NULL, then we know this is part of a pair */
extern void
Stage3end_print (Filestring_T fp, Stage3pair_T stage3pair, T this,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq,
		 Shortread_T headerseq, char *acc_suffix, bool invertp,
		 T hit5, T hit3, int pairedlength, int pairscore,
		 Pairtype_T pairtype, int mapq_score, bool first_read_p);

extern Pairtype_T
Stage3pair_determine_pairtype (Stage3pair_T this);
extern bool
Stage3pair_circularp (Stage3pair_T this);
extern bool
Stage3pair_altlocp (Stage3pair_T this);
extern void
Stage3end_print_transcript_info (Filestring_T fp, Intlist_T trnums, Intlist_T trstarts, Intlist_T trends,
				 bool invertp);
extern void
Stage3end_print_transcript_diff (Filestring_T fp, T hit, Intlist_T common_trnums, bool invertp);

extern void
Stage3pair_print_end (Filestring_T fp, Filestring_T fp_failedinput,
		      Result_T result, Resulttype_T resulttype,
		      char initchar, bool firstp, Univ_IIT_T chromosome_iit,
		      Shortread_T queryseq, Shortread_T headerseq1, Shortread_T headerseq2,
		      int maxpaths, bool quiet_if_excessive_p,
		      bool invertp, int quality_shift);

#ifdef RESOLVE_INSIDE_GENERAL
extern List_T
Stage3pair_resolve_insides (List_T hitpairlist, char *queryseq_ptr_5, char *queryseq_ptr_3,
			    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			    Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
			    Listpool_T listpool, Hitlistpool_T hitlistpool);
#endif

extern Stage3pair_T
Stage3pair_new (T hit5, T hit3, int genestrand, Pairtype_T pairtype,
#ifdef RESOLVE_INSIDE_GENERAL
		char *queryuc_ptr_5, char *queryuc_ptr_3,
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		Listpool_T listpool, bool expect_concordant_p, bool transcriptome_guided_p);

struct Simplepair_T *
Stage3pair_merge (int *npairs, int *querylength_merged, char **queryseq_merged, char **quality_merged,
		  Stage3pair_T this, Shortread_T queryseq5, Shortread_T queryseq3,
		  int querylength_5, int querylength_3, int clipdir,
		  int hardclip5_low, int hardclip5_high, int hardclip3_low, int hardclip3_high);

extern void
Stage3pair_privatize (Stage3pair_T *array, int npairs);

extern List_T
Stage3pair_sort_bymatches (List_T hits, Hitlistpool_T hitlistpool);

extern List_T
Stage3pair_remove_overlaps (List_T hitpairlist, Hitlistpool_T hitlistpool,
			    bool translocp, bool finalp);

extern List_T
Stage3pair_resolve_multimapping (List_T hitpairs, Hitlistpool_T hitlistpool);

extern List_T
Stage3pair_filter (List_T hits, Hitlistpool_T hitlistpool,
		   int max_mismatches_5, int max_mismatches_3,
		   int min_coverage_5, int min_coverage_3);

extern Stage3pair_T *
Stage3pair_eval_and_sort (int npaths, int *first_absmq, int *second_absmq,
			  Stage3pair_T *stage3pairarray,
			  char *queryuc_ptr_5, char *queryuc_ptr_3,
			  char *quality_string_5, char *quality_string_3);

extern List_T
Stage3pair_optimal_score (List_T hitpairlist, Hitlistpool_T hitlistpool,
			  int querylength5, int querylength3, bool finalp);

#if 0
extern List_T
Stage3end_unalias_circular (List_T hitlist);
#endif

extern List_T
Stage3pair_remove_circular_alias (List_T hitpairlist, Hitlistpool_T hitlistpool);

#undef T
#endif

