/* $Id: stage3.h 218191 2019-01-17 13:35:57Z twu $ */
#ifndef STAGE3_INCLUDED
#define STAGE3_INCLUDED

typedef struct Stage3middle_T *Stage3middle_T;
typedef struct Stage3_T *Stage3_T;

#include "bool.h"
#include "sense.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "types.h"
#include "list.h"
#include "sequence.h"
#include "genome.h"
#include "stage2.h"
#include "pairpool.h"
#include "diagpool.h"
#include "cellpool.h"
#include "splicetrie.h"
#include "splicetrie_build.h"	/* For Splicetype_T */
#include "dynprog.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "reader.h"		/* For cDNAEnd_T */
#include "chimera.h"
#include "stopwatch.h"
#ifdef PMAP
#include "oligoindex_pmap.h"
#else
#include "oligoindex_hr.h"
#endif
#include "filestring.h"
#include "output.h"		/* For Printtype_T */

#ifndef GSNAP
#include "gregion.h"
#endif

#define EXTRAQUERYGAP 20

/* POST_CANONICAL is the path_compute_final() step */
/* POST_TRIM is the path_trim() step */
typedef enum {NO_STAGE3DEBUG, POST_STAGE2, POST_SINGLES, POST_INTRONS,
	      POST_HMM, POST_SMOOTHING, POST_DUAL_INTRONS, POST_CYCLES, POST_DUAL_BREAKS,
	      POST_MIDDLE, POST_ENDS, POST_CANONICAL, POST_TRIM, POST_CHANGEPOINT, POST_DISTAL_MEDIAL} Stage3debug_T;

#define T Stage3_T

extern void
Stage3_setup (bool splicingp_in, bool novelsplicingp_in, bool require_splicedir_p_in,
	      Chrpos_T shortsplicedist_novelend, 
	      IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
	      int donor_typeint_in, int acceptor_typeint_in,
	      Univcoord_T *splicesites_in, bool *circularp_in, bool *altlocp_in,
	      Univcoord_T *alias_starts_in, Univcoord_T *alias_ends_in,
	      int min_intronlength_in, int max_deletionlength_in, int min_indel_end_matches_in,
	      int maxpeelback_distalmedial_in, int nullgap_in,
	      int extramaterial_end_in, int extramaterial_paired_in,
	      int extraband_single_in, int extraband_end_in, int extraband_paired_in,
	      int ngap_in, int maxintronlen_in, int maxintronlen_ends_in, int minendexon_in,
	      bool homopolymerp_in, Stage3debug_T stage3debug_in, Univcoord_T genome_totallength_in);

extern Chrnum_T
Stage3middle_chrnum (Stage3middle_T this);
extern Univcoord_T
Stage3middle_chroffset (Stage3middle_T this);
extern Univcoord_T
Stage3middle_chrhigh (Stage3middle_T this);
extern Chrpos_T
Stage3middle_chrlength (Stage3middle_T this);

extern bool
Stage3middle_watsonp (Stage3middle_T this);
extern int
Stage3middle_genestrand (Stage3middle_T this);

extern int
Stage3middle_goodness (Stage3middle_T this);
extern void
Stage3middle_free (Stage3middle_T *old);
extern int
Stage3middle_cmp (const void *a, const void *b);

extern bool
Stage3_chimera_left_p (T this);
extern bool
Stage3_chimera_right_p (T this);
extern bool
Stage3_watsonp (T this);
extern int
Stage3_genestrand (T this);
extern int
Stage3_cdna_direction (T this);
extern int
Stage3_sensedir (T this);
extern int
Stage3_straintype (T this);
extern int
Stage3_goodness (T this);
extern int
Stage3_absmq_score (T this);
extern int
Stage3_mapq_score (T this);
extern List_T
Stage3_pairs (T this);
extern struct Pair_T *
Stage3_pairarray (T this);
extern int
Stage3_npairs (T this);
extern int
Stage3_matches (T this);
extern int
Stage3_mismatches (T this);
extern int
Stage3_indels (T this);

extern int
Stage3_querystart (T this);
extern int
Stage3_queryend (T this);

extern bool
Stage3_joinable_left_p (T this);
extern bool
Stage3_joinable_right_p (T this);
extern void
Stage3_clear_joinable (T this);
extern void
Stage3_set_joinable_left (T this);
extern void
Stage3_set_joinable_right (T this);

extern void
Stage3_print_ends (T this);
extern Chrnum_T
Stage3_chrnum (T this);
extern Univcoord_T
Stage3_chroffset (T this);
extern Univcoord_T
Stage3_chrhigh (T this);
extern Chrpos_T
Stage3_chrlength (T this);
extern bool
Stage3_altloc_chr (Univcoord_T *alias_start, Univcoord_T *alias_end, T this);
extern Chrpos_T
Stage3_chrstart (T this);
extern Chrpos_T
Stage3_chrend (T this);
extern Univcoord_T
Stage3_genomicstart (T this);
extern Univcoord_T
Stage3_genomicend (T this);
extern void
Stage3_set_genomicend (T this, Univcoord_T genomicend);
extern int
Stage3_circularpos (T this);

extern int
Stage3_translation_start (T this);
extern int
Stage3_translation_end (T this);
extern int
Stage3_domain (T this);
extern int
Stage3_largemargin (int *newstart, int *newend, T this, int queryntlength);

extern double
Stage3_fracidentity (T this);
extern Univcoord_T
Stage3_genomicpos (T this, int querypos, bool headp);
extern int
Stage3_chimeric_goodness (int *matches1, int *matches2, T part1, T part2, int breakpoint);

extern bool
Stage3_passes_filter (T this, double min_trimmed_coverage, double min_identity);
extern bool
Stage3_passes_filter_chimera (Chimera_T chimera, double min_trimmed_coverage, double min_identity);
extern int
Stage3_cmp (const void *a, const void *b);
extern Chrpos_T
Stage3_genomiclength (T this);
extern int
Stage3_position_cmp (const void *a, const void *b);
extern int
Stage3_querystart_cmp (const void *a, const void *b);
extern int
Stage3_queryend_cmp (const void *a, const void *b);
extern int
Stage3_chrnum_cmp (const void *a, const void *b);
extern int
Stage3_chrnum_querystart_cmp (const void *a, const void *b);
extern int
Stage3_chrnum_queryend_cmp (const void *a, const void *b);
extern int
Stage3_identity_cmp (const void *a, const void *b);
extern bool
Stage3_overlap (T x, T y);

extern void
Stage3_compute_mapq (List_T stage3list);
extern void
Stage3_count_paths (int *npaths_primary, int *npaths_altloc, List_T stage3list);
extern void
Stage3_recompute_goodness (List_T stage3list);
extern void
Stage3_recompute_coverage (List_T stage3list, Sequence_T queryseq);
extern void
Stage3_free (T *old);

extern bool
Stage3_test_bounds (T this, int minpos, int maxpos);

#ifdef PMAP
extern void
Stage3_translate_cdna (T this, Sequence_T queryaaseq, bool strictp);
extern void
Stage3_backtranslate_cdna (T this);
#else
extern void
Stage3_translate_genomic (T this, int npairs, bool fulllengthp, int cds_startpos, int querylength,
			  bool truncatep, bool strictp);
#endif
extern void
Stage3_translate (T this,
#ifdef PMAP
		  Sequence_T queryseq,
#endif
		  int querylength, bool fulllengthp,
		  int cds_startpos, bool truncatep, bool strictp);
extern void
Stage3_translate_chimera (T this, T mate,
#ifdef PMAP
			  Sequence_T queryseq,
#endif
			  int querylength, bool fulllengthp,
			  int cds_startpos, bool truncatep, bool strictp);
extern void
Stage3_print_pathsummary (Filestring_T fp, T this, int pathnum, Univ_IIT_T chromosome_iit, Univ_IIT_T contig_iit, 
			  IIT_T altstrain_iit, Sequence_T queryseq, char *dbversion);
extern void
Stage3_print_pslformat_nt (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryseq);
#ifdef PMAP
extern void
Stage3_print_pslformat_pro (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryseq, bool strictp);
#endif
extern void
Stage3_print_gff3 (Filestring_T fp, T this, int pathnum, Univ_IIT_T chromosome_iit, Sequence_T usersegment,
		   Sequence_T queryseq, int querylength, Printtype_T printtype, char *sourcename);
#ifndef PMAP
extern void
Stage3_print_bedpe (Filestring_T fp, T this, Univ_IIT_T chromosome_iit);
extern void
Stage3_print_sam (Filestring_T fp, char *abbrev, T this, int pathnum, int npaths_primary, int npaths_altloc,
		  int absmq_score, int second_absmq, int mapq_score,
		  Univ_IIT_T chromosome_iit, Sequence_T usersegment,
		  Sequence_T queryseq, int chimera_part, Chimera_T chimera,
		  int quality_shift, bool sam_paired_p, char *sam_read_group_id);
#endif
extern void
Stage3_print_iit_map (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_iit_exon_map (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_splicesites (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_introns (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T queryseq);

extern void
Stage3_print_map (Filestring_T fp, T this, IIT_T map_iit, int *map_divint_crosstable, Univ_IIT_T chromosome_iit,
		  int pathnum, bool map_exons_p, bool map_bothstrands_p, int nflanking, bool print_comment_p);
extern void
Stage3_print_alignment (Filestring_T fp, T this, Genome_T genome,
			Univ_IIT_T chromosome_iit, Printtype_T printtype,
			bool continuousp, bool continuous_by_exon_p, bool genomefirstp,
			int invertmode, bool nointronlenp, int wraplength);

extern void
Stage3_print_coordinates (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, int invertmode);
extern void
Stage3_print_cdna (Filestring_T fp, T this, int wraplength);

extern void
Stage3_print_protein_genomic (Filestring_T fp, T this, int wraplength);

extern void
Stage3_print_compressed (Filestring_T fp, T this, Sequence_T queryseq, Univ_IIT_T chromosome_iit,
			 char *dbversion, Sequence_T usersegment, int pathnum, int npaths,
			 bool checksump, int chimerapos, int chimeraequivpos,
			 double donor_prob, double acceptor_prob, int chimera_cdna_direction);


extern T
Stage3_new (struct Pair_T *pairarray, List_T pairs, int npairs, int goodness,
	    int cdna_direction, int sensedir,
	    int matches, int unknowns, int mismatches, int qopens, int qindels,
	    int topens, int tindels, int ncanonical, int nsemicanonical, int nnoncanonical, 
	    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
	    bool watsonp, int genestrand, int querylength, int skiplength, int trimlength,
	    int straintype, char *strain, IIT_T altstrain_iit);

extern T
Stage3_new_from_pairs (List_T pairs, int cdna_direction, bool watsonp, int genestrand, int sensedir,
		       Pairpool_T pairpool, Sequence_T queryseq, int query_subseq_offset,
		       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength);

extern bool
Stage3_short_alignment_p (struct Pair_T *pairarray, int npairs, int querylength);

extern bool
Stage3_bad_stretch_p (struct Pair_T *pairarray, int npairs, int pos5, int pos3);

extern int
Stage3_good_part (struct Pair_T *pairarray, int npairs, int pos5, int pos3);

extern Stage3middle_T
Stage3_compute_middle (List_T stage2pairs, List_T all_stage2_starts, List_T all_stage2_ends,
#ifdef PMAP
		       char *queryaaseq_ptr,
#endif
		       char *queryseq_ptr, char *queryuc_ptr, int querylength,
		       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		       bool watsonp, int genestrand, bool jump_late_p, int maxpeelback,
#ifndef GSNAP
		       Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool, 
#endif
		       Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		       int sense_try);

extern struct Pair_T *
Stage3_compute_ends (int *cdna_direction, int *sensedir, List_T *finalpairs1, int *npairs1, int *goodness1,
		     int *matches1, int *nmatches_posttrim_1, int *max_match_length_1,
		     int *ambig_end_length_5_1, int *ambig_end_length_3_1,
		     Splicetype_T *ambig_splicetype_5_1, Splicetype_T *ambig_splicetype_3_1,
		     double *ambig_prob_5_1, double *ambig_prob_3_1,
		     int *unknowns1, int *mismatches1, int *qopens1, int *qindels1, int *topens1, int *tindels1,
		     int *ncanonical1, int *nsemicanonical1, int *nnoncanonical1, double *avg_splice_score_1,
#ifdef GSNAP
		     struct Pair_T **pairarray2, List_T *finalpairs2, int *npairs2, int *goodness2,
		     int *matches2, int *nmatches_posttrim_2, int *max_match_length_2,
		     int *ambig_end_length_5_2, int *ambig_end_length_3_2,
		     Splicetype_T *ambig_splicetype_5_2, Splicetype_T *ambig_splicetype_3_2,
		     double *ambig_prob_5_2, double *ambig_prob_3_2,
		     int *unknowns2, int *mismatches2, int *qopens2, int *qindels2, int *topens2, int *tindels2,
		     int *ncanonical2, int *nsemicanonical2, int *nnoncanonical2, double *avg_splice_score_2,
#endif

		     Stage3middle_T stage3middle,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, int querylength,
		     int skiplength, int query_subseq_offset,
		     Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		     int maxpeelback, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		     int sense_filter, Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool);

extern struct Pair_T *
Stage3_compute_one (int *cdna_direction, int *sensedir,
		    List_T *pairs1, int *npairs1, int *goodness1,
		    int *matches1, int *nmatches_posttrim_1, int *max_match_length_1,
		    int *ambig_end_length_5_1, int *ambig_end_length_3_1,
		    Splicetype_T *ambig_splicetype_5_1, Splicetype_T *ambig_splicetype_3_1,
		    double *ambig_prob_5_1, double *ambig_prob_3_1,
		    int *unknowns1, int *mismatches1, int *qopens1, int *qindels1, int *topens1, int *tindels1,
		    int *ncanonical1, int *nsemicanonical1, int *nnoncanonical1, double *avg_splice_score_1,
#ifdef GSNAP
		    struct Pair_T **pairarray2, List_T *pairs2, int *npairs2, int *goodness2,
		    int *matches2, int *nmatches_posttrim_2, int *max_match_length_2,
		    int *ambig_end_length_5_2, int *ambig_end_length_3_2,
		    Splicetype_T *ambig_splicetype_5_2, Splicetype_T *ambig_splicetype_3_2,
		    double *ambig_prob_5_2, double *ambig_prob_3_2,
		    int *unknowns2, int *mismatches2, int *qopens2, int *qindels2, int *topens2, int *tindels2,
		    int *ncanonical2, int *nsemicanonical2, int *nnoncanonical2, double *avg_splice_score_2,
#endif

		    List_T stage2pairs, List_T all_stage2_starts, List_T all_stage2_ends,
#ifdef PMAP
		    char *queryaaseq_ptr,
#endif
		    char *queryseq_ptr, char *queryuc_ptr, int querylength,
		    int skiplength, int query_subseq_offset,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		    bool watsonp, int genestrand, bool jump_late_p,
		    int maxpeelback,
#ifndef GSNAP
		    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool, 
#endif
		    Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    int sense_try, int sense_filter);

#ifndef GSNAP
extern T
Stage3_direct (Gregion_T gregion,
#ifdef PMAP
	       Sequence_T queryaaseq,
#endif
	       Sequence_T queryseq, Sequence_T queryuc, Pairpool_T pairpool, Genome_T genome,
	       Chrnum_T chrnum,  Univcoord_T chroffset, Chrpos_T chrpos, bool watsonp,
	       int ngap, Dynprog_T dynprogL, Dynprog_T dynprogR,
	       int extramaterial_end, int extraband_end);
#endif

extern bool
Stage3_mergeable (Stage3_T firstpart, Stage3_T secondpart, int exonexonpos, int queryntlength);

extern bool
Stage3_merge_chimera (T *new_left, T *new_right, T old_left, T old_right,
		      int minpos1, int maxpos1, int minpos2, int maxpos2,
		      Sequence_T queryseq,
		      char *queryseq_ptr, char *queryuc_ptr, Pairpool_T pairpool, 
		      Dynprog_T dynprogL, Dynprog_T dynprogR, int maxpeelback);
extern void
Stage3_extend_right (T this, int goal, int querylength,
		     char *queryseq_ptr, char *queryuc_ptr,
		     bool max_extend_p, Pairpool_T pairpool,
		     int genestrand, int maxpeelback);
extern void
Stage3_extend_left (T this, int goal,
		    char *queryseq_ptr, char *queryuc_ptr,
		    bool max_extend_p, Pairpool_T pairpool,
		    int genestrand, int maxpeelback);

extern void
Stage3_trim_right (T this, int goal, char *queryseq_ptr, Pairpool_T pairpool);
extern void
Stage3_trim_left (T this, int goal, char *queryseq_ptr, Pairpool_T pairpool);

extern T
Stage3_merge_local (T old_left, T old_right,
		    int minpos1, int maxpos1, int minpos2, int maxpos2,
		    Sequence_T queryseq,
#ifdef PMAP
		    char *queryaaseq_ptr,
#endif
		    char *queryseq_ptr, char *queryuc_ptr,
#ifndef GSNAP
		    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool, 
#endif
		    Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    int maxpeelback);

extern List_T
Stage3_split (T this, Sequence_T queryseq, Pairpool_T pairpool);

#ifndef PMAP
extern void
Stage3_guess_cdna_direction (T this);
#endif

#undef T
#endif
