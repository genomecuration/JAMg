/* $Id: stage3hrdef.h 218689 2019-03-19 17:21:12Z twu $ */
#ifndef STAGE3HRDEF_INCLUDED
#define STAGE3HRDEF_INCLUDED

#include "bool.h"
#include "types.h"
#include "method.h"
#include "chrnum.h"
#include "iit-read.h"		/* For Overlap_T */
#include "stage3hr.h"		/* For Hittype_T */
#include "list.h"
#include "substring.h"
#include "resulthr.h"		/* For Pairtype_T */


/* Note: Substring_T has genomiclength, but not Stage3end_T */

/* TODO: Allow a Stage3end_T object to hold solutions for both
   sensedirs.  Then the pairing operations can select the sensedirs
   that are best */
#define T Stage3end_T
struct T {
  Hittype_T hittype;
  Method_T method;

  int level;

  int genestrand;
  bool improved_by_gmap_p;	/* true if GMAP alignment based on this hit is better */

  bool distant_splice_p;	/* Indicates a distant splice (same
				   chromosome or different ones).
				   Used for filtering (e.g., by
				   DISTANT_SPLICE_SPECIAL) */

  Chrnum_T chrnum; /* Needed for printing paired-end results.  A
		      chrnum of 0 indicates a translocation (two
		      different chromosomes), or an alignment that
		      needs to be printed as a translocation
		      (samechr_splice unless --merge_samechr is
		      selected).  Used for printing. */

  Chrnum_T effective_chrnum;	/* For determining concordance */
  Chrnum_T other_chrnum;	/* 0 for non-translocations, and other chrnum besides effective_chrnum for translocations */
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;

  int querylength;		/* Needed for overlap and pairlength calculations */
  int querylength_adj;		/* Adjusted for insertions */

  Univcoord_T genomicstart;
  Univcoord_T genomicend;
  bool plusp;

  Univcoord_T low;
  Univcoord_T high;

  Chrpos_T genomiclength;
  Chrpos_T guided_insertlength; /* Used only by Stage3end_eval_and_sort_guided */

  float mapq_loglik;
  int mapq_score;
  int absmq_score;		/* Absolute MAPQ, for XQ and X2 flags */

  int nsegments;

  /* score = querylength - nmatches - penalties */
  /* score_posttrim = querylength - nmatches_posttrim - penalties */
  /* score (best case) <= score_posttrim (worst case) */
  /* In computing concordance, can rank by score to get ambiguous
     alignments, but keep going until we reach score_posttrim to give
     complete alignments a chance */

  int score_overall;	      /* Over entire query, so trimming raises the score */
  int nmatches_to_trims;      /* Ignores any parts that are trimmed */

  int score_within_trims;	/* From one trim to the other */
  int nmatches_plus_spliced_trims;  /* Includes alts and ambiguous parts after good splice ends */

  double splice_score;		/* Used by various SPLICE types */

  /* if trim_querystart_splicep or trim_queryend_splicep is true, then trim is of type "unknown amb" */
  /* if trim_querystart_splicep or trim_queryend_splicep is false, then trim is of type "unknown" */

  int trim_querystart; /* Used by Stage3end_optimal_score for comparing terminals and non-terminals */
  int trim_queryend;
  bool trim_querystart_splicep;
  bool trim_queryend_splicep;

#if 0
  int penalties;		/* Indel penalties */
#endif
  int score_eventrim;		/* Temporary storage used by Stage3end_optimal_score */

  Overlap_T gene_overlap;
  long int tally;

  int nmismatches_bothdiff;
  int nmismatches_refdiff;	/* Set only for display */

  int nindels;			/* for indels */

  Chrpos_T distance;	/* for splicing or shortexon (sum of two distances) */

  /* For spliced alignments */
  int sensedir_for_concordance;
  int sensedir;			/* a private value */
  /* Possibilities:
     not spliced: sensedir_for_concordance is NULL.  sensedir_private is NULL.
     regular or transloc splice with certain sense: sensedir_for_concordance is {FORWARD,ANTI}.  sensedir(private) is the same {FORWARD,ANTI}.
     distant splice with uncertain sense: sensedir_for_concordance is NULL.  sensedir(private) is {FORWARD,ANTI}.
  */


  int nsplices;

#if 0
  bool start_ambiguous_p;
  bool end_ambiguous_p;
  int amb_length_donor;	/* For shortexon only */
  int amb_length_acceptor;	/* For shortexon only */
  double amb_prob_donor;	/* For shortexon */
  double amb_prob_acceptor;	/* For shortexon */
#endif

#if 0
  Univcoord_T *start_alts_coords;	/* Pointer to either alts_coords_donor or alts_coords_acceptor */
  Univcoord_T *end_alts_coords;	/* Pointer to either alts_coords_donor or alts_coords_acceptor */
  int start_alts_ncoords;		/* Equal to either nalts_coords_donor or nalts_coords_acceptor */
  int end_nalts_coords;           /* Equal to either nalts_coords_donor or nalts_coords_acceptor */
  Univcoord_T *alts_coords_donor;
  Univcoord_T *alts_coords_acceptor;
  int nalts_coords_donor;
  int nalts_coords_acceptor;

  int *start_alts_knowni;        /* Pointer to either alts_knowni_donor or alts_knowni_acceptor */
  int *end_alts_knowni;          /* Pointer to either alts_knowni_donor or alts_knowni_acceptor */
  int *alts_knowni_donor;
  int *alts_knowni_acceptor;

  int *start_alts_nmismatches;	/* Pointer to either alts_nmismatches_donor or alts_nmismatches_acceptor */
  int *end_alts_nmismatches;     /* Pointer to either alts_nmismatches_donor or alts_nmismatches_acceptor */
  double *start_alts_probs;	/* Pointer to either alts_probs_donor or alts_probs_acceptor */
  double *end_alts_probs;	/* Pointer to either alts_probs_donor or alts_probs_acceptor */

  int *alts_nmismatches_donor;
  int *alts_nmismatches_acceptor;
  double *alts_probs_donor;
  double *alts_probs_acceptor;
#endif


  /* For transcriptome alignment.  For fusions, transcripts corresponding to substring_for_concordance */
  List_T transcripts;

  /*  For fusions, transcripts corresponding to substring_other */
  List_T transcripts_other;

  List_T substrings_1toN;	/* query position 1 to N */
  List_T substrings_Nto1;	/* query position N to 1.  Keeps only pointers to the substrings. */
  List_T substrings_LtoH;	/* Chromosomal low-to-high.  Keeps only pointers to the substrings. */
  List_T substrings_HtoL;	/* Chromosomal high-to-low.  Keeps only pointers to the substrings. */

  Substring_T substring_for_concordance;
  Substring_T substring_other;

  List_T junctions_LtoH;
  /* List_T junctions_HtoL; */
  List_T junctions_1toN;
  List_T junctions_Nto1;

  bool paired_usedp;
  bool concordantp;    /* for paired-end.  set to true by Stage3_pair_up(). */

  int query_splicepos;		/* For splices.  Relative to querystart, so different from circularpos */

  int circularalias;			/* -1 if all below chrlength, 0 if straddles or NA (e.g., transloc), and +1 if above */
                                /* -2 if extends below beginning of circular chromosome, +2 if extends beyond end of second copy */
  int circularpos;		/* if circularalias == 0, then amount of queryseq below chrlength.  Defined relative to low */

  bool altlocp;
};


struct Stage3pair_T {
  Pairtype_T pairtype;
  int genestrand;

  T hit5;			/* Always a copy from the original */
  T hit3;			/* Always a copy from the original */

  Univcoord_T low;
  Univcoord_T high;
  Chrpos_T insertlength;
  int pair_relationship;
  int insertlength_expected_sign;	/* 1 if in (expected_pairlength_low, expected_pairlength_high),
					   0 if in (expected_pairlength_low, expected_pairlength_very_high), and
					   -1 if < expected_pairlength_low or > expected_pairlength_very_high */

  Chrpos_T outerlength;

  float mapq_loglik;
  int mapq_score;
  int absmq_score;

  int score_overall;
  int nmatches_to_trims;

  int score_within_trims;
  int nmatches_plus_spliced_trims;

  int nmismatches;		/* querylength - sum of nmatches */
  int score_eventrim;

  Overlap_T gene_overlap;
  long int tally;

#ifdef USE_ABSDIFFLENGTH
  Chrpos_T absdifflength;
#endif
#ifdef USE_BINGO
  bool absdifflength_bingo_p;
#endif
  int dir;			/* -1, 0, or +1 */
  bool sense_consistent_p;

  int nsplices;

  bool circularp;		/* If either hit5 or hit3 are circular */
  int alts_resolve_5; /* Resolution of ambiguous end for this particular pair */
  int alts_resolve_3; /* Resolution of ambiguous end for this particular pair */
  int alts_status_inside;

  /* After performing Transcript_concordance */
  /* Need to have separate intersection results, in case an end pairs with a different end. */
  List_T transcripts5;
  List_T transcripts3;
};

#undef T
#endif

