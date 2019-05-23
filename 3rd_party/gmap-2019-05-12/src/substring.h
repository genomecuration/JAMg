/* $Id: substring.h 218675 2019-03-16 01:25:48Z twu $ */
#ifndef SUBSTRING_INCLUDED
#define SUBSTRING_INCLUDED

typedef enum {GMAP_NOT_APPLICABLE, GMAP_VIA_SUBSTRINGS, GMAP_VIA_SEGMENTS, GMAP_VIA_REGION} GMAP_source_T;
typedef enum {END, INS, DEL, FRAG, DON, ACC, AMB_DON, AMB_ACC, TERM} Endtype_T;

typedef struct Substring_T *Substring_T;

#include <stdio.h>
#include "mode.h"
#include "genomicpos.h"
#include "types.h"
#include "chrnum.h"
#include "shortread.h"
#include "genome.h"
#include "compress.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "bool.h"
#include "filestring.h"
#include "junction.h"
#include "intlist.h"
#include "list.h"
#ifdef LARGE_GENOMES
#include "uint8list.h"
#else
#include "uintlist.h"
#endif
#include "output.h"


#define T Substring_T

extern void
Substring_alias_circular (T this);
extern void
Substring_unalias_circular (T this);

extern void
Substring_free (T *old);
extern void
Substring_list_gc (List_T *old);

extern bool
Substring_contains_p (T this, int querypos);
extern int
Substring_compare (T substring1, T substring2, int alias1, int alias2, Chrpos_T chrlength1, Chrpos_T chrlength2);
extern bool
Substring_equal_p (T substring1, T substring2);

extern bool
Substring_overlap_p (T substring1, T substring2);
extern Chrpos_T
Substring_insert_length (int *pair_relationship, T substring5, T substring3);
extern bool
Substring_overlap_point_trimmed_p (T substring, Univcoord_T endpos);
extern Univcoord_T
Substring_overlap_segment_trimmed (T substring1, T substring2);

extern int
Substring_trim_qstart_nosplice (int *nmismatches, int *mismatch_positions_alloc,
				Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				bool plusp,int genestrand);
extern int
Substring_trim_qend_nosplice (int *nmismatches, int *mismatch_positions_alloc,
			      Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
			      bool plusp, int genestrand);

extern bool
Substring_trimmed_qstarts (int *result, Splicetype_T *splicetype, int **ambig_qstarts, double **ambig_probs_5,
			   Univcoord_T left, int qend,
			   bool plusp, int genestrand, int *mismatch_positions_alloc, Compress_T query_compress,
			   Univcoord_T chroffset, int sensedir);
extern bool
Substring_qstart_trim (int *trimpos, Splicetype_T *splicetype, double *ambig_prob_qstart,
		       Univcoord_T left, int qend,
		       bool plusp, int genestrand, int *mismatch_positions_alloc, Compress_T query_compress,
		       Univcoord_T chroffset, int sensedir);

extern bool
Substring_trimmed_qends (int *result, Splicetype_T *splicetype, int **ambig_qends, double **ambig_probs_3,
			 Univcoord_T left, int qstart, int querylength,
			 bool plusp, int genestrand, int *mismatch_positions_alloc, Compress_T query_compress,
			 Univcoord_T chroffset, int sensedir);

extern bool
Substring_qend_trim (int *trimpos, Splicetype_T *splicetype, double *ambig_prob_qend,
		     Univcoord_T left, int qstart, int querylength,
		     bool plusp, int genestrand, int *mismatch_positions_alloc, Compress_T query_compress,
		     Univcoord_T chroffset, int sensedir);

extern int
Substring_compute_nmatches (Univcoord_T left, int querystart, int queryend, int querylength,
			    bool plusp, int genestrand, Compress_T query_compress,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			    bool splice_querystart_p, bool splice_queryend_p, bool chrnum_fixed_p);

extern T
Substring_new (int nmismatches, Univcoord_T left, int querystart, int queryend, int querylength,
	       bool plusp, int genestrand, Compress_T query_compress,
	       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
	       bool splice_querystart_p, Splicetype_T splicetype_querystart, double ambig_prob_querystart, 
	       bool splice_queryend_p, Splicetype_T splicetype_queryend, double ambig_prob_queryend,
	       int sensedir);

extern void
Substring_set_querystart_pretrim (T this, int querystart_pretrim);
extern void
Substring_set_queryend_pretrim (T this, int queryend_pretrim);

extern T
Substring_new_simple (int nmismatches, Univcoord_T left, int querystart, int queryend, int querylength,
		      bool plusp, int genestrand, Compress_T query_compress,
		      Endtype_T start_endtype, Endtype_T end_endtype,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      int sensedir);

extern T
Substring_new_alts_D (int querystart, int queryend, int splice_pos, int querylength,
		      bool plusp, int genestrand, Compress_T query_compress,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      Univcoord_T *alts_coords, int *alts_knowni, int *alts_nmismatches, double *alts_probs,
		      int alts_ncoords, double alts_common_prob, bool substring1p);

extern T
Substring_new_alts_A (int querystart, int queryend, int splice_pos, int querylength,
		      bool plusp, int genestrand, Compress_T query_compress,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      Univcoord_T *alts_coords, int *alts_knowni, int *alts_nmismatches, double *alts_probs,
		      int alts_ncoords, double alts_common_prob, bool substring1p);

extern Univcoord_T
Substring_set_alt (double *donor_prob, double *acceptor_prob, Univcoord_T *genomicstart, Univcoord_T *genomicend,
		   T this, int bingoi);

extern float
Substring_compute_mapq (T this, char *quality_string);

extern int
Substring_display_prep (T this, char *queryuc_ptr, int querylength,
			int extraleft, int extraright, Genome_T genome);
extern char *
Substring_genomic_sequence (int *seqlength, T this, Genome_T genome);

extern Univcoord_T
Substring_left (T this);
extern Univcoord_T
Substring_splicecoord_D (T this);
extern Univcoord_T
Substring_splicecoord_A (T this);
extern char
Substring_chimera_strand (T this);
extern Chrpos_T
Substring_chr_splicecoord_D (T this, char donor_strand);
extern Chrpos_T
Substring_chr_splicecoord_A (T this, char acceptor_strand);
extern int
Substring_splicesitesD_knowni (T this);
extern int
Substring_splicesitesA_knowni (T this);

extern bool
Substring_plusp (T this);
extern int
Substring_sensedir (T this);
extern int
Substring_genestrand (T this);
extern char *
Substring_genomic_bothdiff (T this);
extern char *
Substring_genomic_refdiff (T this);
extern int
Substring_nmismatches_bothdiff (T this);
extern int
Substring_nmismatches_refdiff (T this);
extern int
Substring_nmatches_to_trims (T this);

/* Returns nmatches_plus_spliced_trims */
extern int
Substring_nmatches (T this);

extern void
Substring_set_nmismatches_terminal (T this, int nmismatches_whole, int nmismatches_bothdiff);
extern Endtype_T
Substring_start_endtype (T this);
extern Endtype_T
Substring_end_endtype (T this);
extern float
Substring_mapq_loglik (T this);
extern int
Substring_trim_querystart (T this);
extern int
Substring_trim_queryend (T this);
extern bool
Substring_trim_querystart_splicep (T this);
extern bool
Substring_trim_queryend_splicep (T this);
extern int
Substring_mandatory_trim_querystart (T this);
extern int
Substring_mandatory_trim_queryend (T this);

extern int
Substring_querystart (T this);
extern int
Substring_querystart_pretrim (T this);
extern int
Substring_queryend (T this);
extern int
Substring_querylength (T this);
extern int
Substring_match_length (T this);
extern int
Substring_match_length_pretrim (T this);
extern int
Substring_amb_length (T this);
extern int
Substring_start_amb_length (T this);
extern int
Substring_end_amb_length (T this);

extern Chrpos_T
Substring_genomic_alignment_length (T this);

extern Chrnum_T
Substring_chrnum (T this);
extern Univcoord_T
Substring_chroffset (T this);
extern Univcoord_T
Substring_chrhigh (T this);
extern Chrpos_T
Substring_chrlength (T this);
extern Chrpos_T
Substring_chrpos_low (T this);
extern Chrpos_T
Substring_chrpos_high (T this);

extern Chrpos_T
Substring_alignstart_trim_chr (T this);
extern Chrpos_T
Substring_alignend_trim_chr (T this);
extern Univcoord_T
Substring_alignstart_trim (T this);
extern Univcoord_T
Substring_alignend_trim (T this);
extern Univcoord_T
Substring_left_genomicseg (T this);
extern Univcoord_T
Substring_genomicstart (T this);
extern Univcoord_T
Substring_genomicend (T this);

extern double
Substring_amb_prob (T this);
extern double
Substring_amb_donor_prob (T this);
extern double
Substring_amb_acceptor_prob (T this);

extern double
Substring_siteD_prob (T this);
extern double
Substring_siteA_prob (T this);

extern int
Substring_siteD_pos (T this);
extern int
Substring_siteA_pos (T this);
extern int
Substring_siteN_pos (T this);

extern bool
Substring_ambiguous_p (T this);
extern bool
Substring_has_alts_p (T this);
extern int
Substring_alts_ncoords (T this);
extern Univcoord_T *
Substring_alts_coords (T this);
extern double
Substring_alts_common_prob (T this);
extern void
Substring_print_alts_coords (T this);
extern int *
Substring_alts_nmismatches (T this);

extern int
Substring_circularpos (T this);


extern T
Substring_copy (T old);

extern T
Substring_new_donor (int nmismatches, Univcoord_T donor_coord, int donor_knowni,
		     int querystart, int queryend, int sitepos, double donor_prob,
		     bool splice_querystart_p, Splicetype_T splicetype_querystart, double ambig_prob_querystart,
		     bool splice_queryend_p, Splicetype_T splicetype_queryend, double ambig_prob_queryend,

		     Univcoord_T left, Compress_T query_compress, int querylength, bool plusp, int genestrand, int sensedir,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength);
extern T
Substring_new_acceptor (int nmismatches, Univcoord_T acceptor_coord, int acceptor_knowni,
			int querystart, int queryend, int sitepos, double acceptor_prob,
			bool splice_querystart_p, Splicetype_T splicetype_querystart, double ambig_prob_querystart,
			bool splice_queryend_p, Splicetype_T splicetype_queryend, double ambig_prob_queryend,

			Univcoord_T left, Compress_T query_compress, int querylength, bool plusp, int genestrand, int sensedir,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength);

extern void
Substring_label_donor (T this, int splice_pos, double donor_prob, int sensedir_distant_guess);
extern void
Substring_label_acceptor (T this, int splice_pos, double acceptor_prob, int sensedir_distant_guess);

extern T
Substring_new_startfrag (int nmismatches, int querystart, int queryend,
			 Univcoord_T left, Compress_T query_compress, int querylength, bool plusp, int genestrand,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength);
extern T
Substring_new_endfrag (int nmismatches, int querystart, int queryend,
		       Univcoord_T left, Compress_T query_compress, int querylength, bool plusp, int genestrand,
		       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength);

extern T
Substring_trim_startfrag (int nmismatches, T old, int new_queryend);
extern T
Substring_trim_endfrag (int nmismatches, T old, int new_querystart);

extern List_T
Substring_sort_siteD_halves (List_T hitlist, Listpool_T listpool, bool ascendingp);
extern List_T
Substring_sort_siteA_halves (List_T hitlist, Listpool_T listpool, bool ascendingp);
extern List_T
Substring_sort_siteN_halves (List_T hitlist, Listpool_T listpool, bool ascendingp);


extern Chrpos_T
Substring_compute_chrpos (T this, int hardclip_low, bool hide_soft_clips_p);

extern double
Substring_evalue (T substring);

extern void
Substring_print_m8 (Filestring_T fp, T substring, Shortread_T headerseq, char *acc_suffix,
		    char *chr, bool invertp);
extern void
Substring_print_alignment (Filestring_T fp, Junction_T pre_junction, T substring, Junction_T post_junction,
			   Shortread_T queryseq, Genome_T genome, char *chr, bool invertp);

extern long int
Substring_tally (T this, IIT_T tally_iit, int *tally_divint_crosstable);

extern bool
Substring_runlength_p (T this, IIT_T runlength_iit, int *runlength_divint_crosstable);


extern int
Substring_count_mismatches_region (T this, int trim_querystart, int trim_queryend);
extern List_T
Substring_convert_to_pairs_out (List_T pairs, T substring, int querylength, Shortread_T queryseq,
				int hardclip_low, int hardclip_high, int queryseq_offset);

extern List_T
Substring_add_insertion_out (List_T pairs, T substringA, T substringB, int querylength,
			     int insertionlength, Shortread_T queryseq,
			     int hardclip_low, int hardclip_high, int queryseq_offset);
extern List_T
Substring_add_deletion_out (List_T pairs, T substringA, T substringB, int querylength,
			    char *deletion, int deletionlength,
			    int hardclip_low, int hardclip_high, int queryseq_offset);
extern List_T
Substring_add_intron_out (List_T pairs, T substringA, T substringB, int querylength,
			  int hardclip_low, int hardclip_high, int queryseq_offset);

extern void
Substring_setup (bool print_nsnpdiffs_p_in, bool print_snplabels_p_in,
		 bool show_refdiff_p_in, IIT_T snps_iit_in, int *snps_divint_crosstable_in,
		 Genome_T genomebits_in, Genome_T genomebits_alt_in,
		 Univ_IIT_T chromosome_iit_in, int nchromosomes_in, int circular_typeint_in,
		 IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		 IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
		 int donor_typeint_in, int acceptor_typeint_in,
		 bool novelsplicingp_in, bool knownsplicingp_in,
		 Outputtype_T output_type_in, Mode_T mode_in, Univcoord_T genomelength_in);


#undef T
#endif


