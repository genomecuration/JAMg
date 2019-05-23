static char rcsid[] = "$Id: stage3.c 218286 2019-01-23 16:46:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage3.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include <math.h>		/* For pow() */

#include "assert.h"
#include "mem.h"
#include "univcoord.h"

#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "comp.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "smooth.h"
#include "scores.h"
#include "intron.h"
#include "pbinom.h"
#include "changepoint.h"
#ifndef GSNAP
#include "translation.h"
#endif
#ifdef PMAP
#include "backtranslation.h"
#endif
#include "complement.h"
#include "iit-read.h"
#include "stage2.h"
#include "dynprog.h"		/* For Dynprog_consistent_p */
#include "dynprog_single.h"
#include "dynprog_genome.h"
#include "dynprog_cdna.h"
#include "dynprog_end.h"
#include "boyer-moore.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "fastlog.h"
#include "doublelist.h"


/* The following are the same as dividing by 2048 and 1024 */
#define goodness_intronlen(x) (x >> 11)
#define goodness_nonintronlen(x) (x >> 10)


#if 0
/* These iterations are very costly to compute */
#define MAXITER_CYCLES 5
#define MAXITER_SMOOTH_BY_SIZE 2
#define MAXITER_INTRONS 2
#define MAXITER_KNOWNSPLICE 1	/* Effectively turns off iteration */

#else
#define MAXITER_CYCLES 1
#define MAXITER_SMOOTH_BY_SIZE 1
#define MAXITER_INTRONS 1
#define MAXITER_KNOWNSPLICE 1	/* Effectively turns off iteration */
#endif


#define INTRON_PENALTY_INCONSISTENT 16
#define NONCANONICAL_PENALTY 12

#define SINGLESLEN 9	/* Should be same as MININTRONLEN */
#define MININTRONLEN 9		/* Determines when Dynprog_genome_gap gets called vs Dynprog_single_gap */
#define MININTRONLEN_FINAL 50	/* Determines when to perform final
				   pass to find canonical introns */
#define SUFF_MATCHES_KEEP 300

/* Old parameters for extrapeel */
/* #define SUFFCONSECUTIVE 5 */
/* #define MAXINCURSION 5 */


#define MINCOVERAGE 0.10  /* Not used anymore */

#define DYNPROGINDEX_MAJOR -1
#define DYNPROGINDEX_MINOR +1

#define DUAL_BREAK_PROB_THRESHOLD 0.90

/* If too small, e.g., 3, misses introns with a nearby mismatch.  If too large, e.g., 24, misses small exons */
#define MIN_STAGE2_FOR_DUALBREAK 6

#define MIN_MICROEXON_LENGTH 3
#define MAX_MICROEXON_LOOPS 10

#define THETA_SLACK 0.10
#define TRIM_END_PVALUE 1e-4

#define NEARBY_INDEL 6
#define INDEL_SPLICE_ENDLENGTH 12
#define NONCANONICAL_ACCEPT 15
#define NONCANONICAL_PERFECT_MATCHES 12

#define MAXPEELBACK_SCORE 5	/* For determining goodness of intron */
#define MAXPEELBACK_END 1000

#define DUALBREAK_QUERYJUMP_FACTOR 10

#define SCORE_SIGDIFF 5
#define PROB_SIGDIFF 0.5

#define END_SPLICESITE_SEARCH 10
#define END_SPLICESITE_PROB_MATCH 0.90
#define END_SPLICESITE_PROB_MISMATCH 0.95

#define MICROEXON_PROB_MATCH 0.50
#define MICROEXON_PROB_MISMATCH 0.80

#define END_MIN_EXONLENGTH 12
#define END_SUFFICIENT_EXONLENGTH_PCT 0.10 /* Defines length (in terms of percentage of query beyond which we can ignore maxintronlen_ends */

#if 0
/* No longer used.  Not sure why it was used before */
#define END_SPLICESITE_EXON_LENGTH 100  /* If shorter than this, then don't look for end splice site */
#endif


/* For Stage3_bad_stretch_p */
#define LOG_99 -0.01005033585
#define LOG_01 -4.605170186

#define LOG_9999 -0.000100005
#define LOG_90 -0.1053605
#define LOG_75 -0.2876821
#define LOG_25 -1.386294
#define LOG_10 -2.302585
#define LOG_0001 -9.21034

#if 0
/* Switches on 5 consecutive mismatches */
#define LOG_99_9999 -0.01015034085
#define LOG_99_0001 -9.220390708
#define LOG_25_0001 -10.59663473
#define LOG_25_9999 -1.386394366
#define LOG_01_9999 -4.605270191
#define LOG_01_0001 -13.81551056
#define LOG_75_0001 -9.498022444
#define LOG_75_9999 -0.2877820775
#endif

#if 1
#define LOG_99_999 -0.01105083619
#define LOG_99_001 -6.917805615
#define LOG_25_001 -8.29404964
#define LOG_25_999 -1.387294861
#define LOG_01_999 -4.606170686
#define LOG_01_001 -11.51292546
#define LOG_75_001 -7.195437351
#define LOG_75_999 -0.2886825728
#endif

#if 0
/* Switches on 4 consecutive mismatches */
#define LOG_99_99 -0.02010067171
#define LOG_99_01 -4.615220522
#define LOG_25_01 -5.991464547
#define LOG_25_99 -1.396344697
#define LOG_01_99 -4.615220522
#define LOG_01_01 -9.210340372
#define LOG_75_01 -4.892852258
#define LOG_75_99 -0.2977324083
#endif


/* #define EXTRACT_GENOMICSEG 1 */


static const Except_T gapcheck_error = {"Gap check failed"};
static const Except_T coordinate_error = {"Coordinate error"};

/* #define SHORTCUT 1 */		/* Skips re-solving introns if already canonical */
#define EXCESS_GAPHOLDERS 1
#define MAXITER 100		/* For peelback */

/* In debug mode, probably want to activate debug in pairpool.c and
   dynprog.c also */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#ifdef DEBUG0
#define debug0(x) x
#else 
#define debug0(x)
#endif

/* Pair dump */
#ifdef DEBUG1
#define debug1(x) x
#else 
#define debug1(x)
#endif

/* trim ends */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Fix adjacent indels */
#ifdef DEBUG4
#define debug4(x) x
#else 
#define debug4(x)
#endif

/* HMM */
#ifdef DEBUG5
#define debug5(x) x
#else 
#define debug5(x)
#endif

/* assign_gap_types and fill_in_gaps */
#ifdef DEBUG7
#define debug7(x) x
#else 
#define debug7(x)
#endif

/* stage3debug */
#ifdef DEBUG8
#define debug8(x) x
#else 
#define debug8(x)
#endif

/* bad_stretch_p */
#ifdef DEBUG9
#define debug9(x) x
#else 
#define debug9(x)
#endif

/* chimera */
#ifdef DEBUG10
#define debug10(x) x
#else 
#define debug10(x)
#endif

/* pick_cdna_direction */
#ifdef DEBUG11
#define debug11(x) x
#else 
#define debug11(x)
#endif

/* splicesitepos */
#ifdef DEBUG12
#define debug12(x) x
#else 
#define debug12(x)
#endif

/* trimming at novel splice sites at ends */
#ifdef DEBUG13
#define debug13(x) x
#else 
#define debug13(x)
#endif

/* build_dual_breaks */
#ifdef DEBUG14
#define debug14(x) x
#else 
#define debug14(x)
#endif

/* changepoint */
#ifdef DEBUG18
#define debug18(x) x
#else 
#define debug18(x)
#endif

/* mergeable */
#ifdef DEBUG20
#define debug20(x) x
#else 
#define debug20(x)
#endif

/* end_compare */
#ifdef DEBUG21
#define debug21(x) x
#else 
#define debug21(x)
#endif

#ifdef DEBUG99
#define debug99(x) x
#else 
#define debug99(x)
#endif




static bool splicingp;
static bool novelsplicingp;
static bool require_splicedir_p;
static Chrpos_T overall_end_distance_linear;
static Chrpos_T overall_end_distance_circular;

static IIT_T splicesites_iit;
static int *splicesites_divint_crosstable;

static int donor_typeint;
static int acceptor_typeint;

static Univcoord_T *splicesites;

static bool *circularp;
static bool *altlocp;
static Univcoord_T *alias_starts;
static Univcoord_T *alias_ends;

static int min_intronlength;
static int max_deletionlength;
static int min_indel_end_matches;

static int maxpeelback_distalmedial;
static int nullgap;
static int extramaterial_end;
static int extramaterial_paired;
static int extraband_single;
static int extraband_end;
static int extraband_paired;
static int ngap;
static int maxintronlen;	/* for middle */
static int maxintronlen_ends;
static int minendexon;

static bool maximize_coverage_p = false;
static Stage3debug_T stage3debug;
static Univcoord_T genome_totallength;

static bool homopolymerp;

void

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
	      bool homopolymerp_in, Stage3debug_T stage3debug_in, Univcoord_T genome_totallength_in) {
  splicingp = splicingp_in;
  novelsplicingp = novelsplicingp_in;
  require_splicedir_p = require_splicedir_p_in;
  if (shortsplicedist_novelend > (Chrpos_T) max_deletionlength_in) {
    overall_end_distance_linear = shortsplicedist_novelend;
  } else {
    overall_end_distance_linear = max_deletionlength_in;
  }
  overall_end_distance_circular = max_deletionlength_in;

  splicesites_iit = splicesites_iit_in;
  splicesites_divint_crosstable = splicesites_divint_crosstable_in;
  donor_typeint = donor_typeint_in;
  acceptor_typeint = acceptor_typeint_in;

  splicesites = splicesites_in;

  circularp = circularp_in;
  altlocp = altlocp_in;
  alias_starts = alias_starts_in;
  alias_ends = alias_ends_in;

  min_intronlength = min_intronlength_in;
  max_deletionlength = max_deletionlength_in;
  min_indel_end_matches = min_indel_end_matches_in;

  maxpeelback_distalmedial = maxpeelback_distalmedial_in;
  nullgap = nullgap_in;
  extramaterial_end = extramaterial_end_in;
  extramaterial_paired = extramaterial_paired_in;
  extraband_single = extraband_single_in;
  extraband_end = extraband_end_in;
  extraband_paired = extraband_paired_in;
  ngap = ngap_in;
  maxintronlen = maxintronlen_in;
  maxintronlen_ends = maxintronlen_ends_in;
  minendexon = minendexon_in;

  homopolymerp = homopolymerp_in;
  stage3debug = stage3debug_in;
  genome_totallength = genome_totallength_in;

  return;
}


/************************************************************************
 *   Stage 3 merges cDNA-genomic pairs from stage 2 (called "path")
 *   and from dynamic programming into a single list.  In this
 *   process, stage 3 may also have to pop a pair off the path, insert
 *   the dynamic programming results (called "gappairs") and then push
 *   the stored pair onto the list.  The relevant pointers are
 *   leftquerypos and leftgenomepos, which refer to the left (stored) pair;
 *   rightquerypos and rightgenomepos, which refer to the right (top) pair on
 *   the list (which may represent what the list should have for the
 *   purposes of dynamic programming); and querydp5, genomedp5,
 *   querydp3, and genomedp3, which refer to the dynamic programming
 *   indices, inclusive.
 *
 *   path has the end of the query sequence as its car.
 *   pairs has the beginning of the query sequence as its car.
 *
 *   Most procedures take the top of path and put it onto pairs:
 *
 *	 <- <- path  =====>  pairs -> ->
 *	   leftpair	     rightpair
 *
 *   For stage 1 and stage 2, the poly-A/T tails, if any, was stripped
 *   off, but for stage 3, we try to extend these tails if possible.
 *   Therefore, we use the full length and full sequence here, and add
 *   the offset to the path from stage 2.
 ************************************************************************/


struct Stage3middle_T {
  int goodness;
  double defect_rate_fwd;
  double defect_rate_rev;
  List_T pairs_fwd;
  List_T pairs_rev;

  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;
  bool watsonp;
  int genestrand;

  List_T all_stage2_starts;
  List_T all_stage2_ends;
};

Chrnum_T
Stage3middle_chrnum (Stage3middle_T this) {
  return this->chrnum;
}

Univcoord_T
Stage3middle_chroffset (Stage3middle_T this) {
  return this->chroffset;
}

Univcoord_T
Stage3middle_chrhigh (Stage3middle_T this) {
  return this->chrhigh;
}

Chrpos_T
Stage3middle_chrlength (Stage3middle_T this) {
  return this->chrlength;
}


bool
Stage3middle_watsonp (Stage3middle_T this) {
  return this->watsonp;
}

int
Stage3middle_genestrand (Stage3middle_T this) {
  return this->genestrand;
}

int
Stage3middle_goodness (Stage3middle_T this) {
  return this->goodness;
}

static Stage3middle_T
Stage3middle_new (int goodness, double defect_rate_fwd, double defect_rate_rev,
		  List_T pairs_fwd, List_T pairs_rev,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		  Chrpos_T chrlength, bool watsonp, int genestrand,
		  List_T all_stage2_starts, List_T all_stage2_ends) {
  Stage3middle_T new = (Stage3middle_T) MALLOC(sizeof(*new));


  new->goodness = goodness;
  new->defect_rate_fwd = defect_rate_fwd;
  new->defect_rate_rev = defect_rate_rev;
  new->pairs_fwd = pairs_fwd;
  new->pairs_rev = pairs_rev;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->watsonp = watsonp;
  new->genestrand = genestrand;

  new->all_stage2_starts = all_stage2_starts;
  new->all_stage2_ends = all_stage2_ends;

  return new;
}

void
Stage3middle_free (Stage3middle_T *old) {
  List_free(&(*old)->all_stage2_starts);
  List_free(&(*old)->all_stage2_ends);
  FREE(*old);
  return;
}

int
Stage3middle_cmp (const void *a, const void *b) {
  Stage3middle_T x = * (Stage3middle_T *) a;
  Stage3middle_T y = * (Stage3middle_T *) b;

  if (x->goodness > y->goodness) {
    return -1;
  } else if (y->goodness > x->goodness) {
    return +1;
  } else {
    return 0;
  }
}


#define T Stage3_T
struct T {
  struct Pair_T *pairarray;	/* The array version of pairs_fwd or pairs_rev, with the gaps substituted */
  bool pairarray_freeable_p;
  bool chimera_left_p;		/* Part of a chimera on its querystart end */
  bool chimera_right_p;		/* Part of a chimera on its queryend end */
  int npairs;

  List_T pairs;			/* Winning set of pairs */

  int straintype;
  char *strain;

  Chrnum_T chrnum;
  Univcoord_T chroffset;	/* Start of chromosome chrnum on genome */
  Univcoord_T chrhigh;	        /* End of chromosome chrnum on genome */
  Chrpos_T chrlength;
  int circularpos;

  Chrpos_T chrstart; /* Position on chromosome of start of genomicseg */
  Univcoord_T genomicstart;	/* Start of alignment */
  Univcoord_T genomicend;	/* End of alignment */
  int cdna_direction;
  int sensedir;
  bool watsonp;
  int genestrand;

  double trimmed_coverage;
  int matches;
  int unknowns;
  int mismatches;
  int qopens;
  int qindels;
  int topens;
  int tindels;
  int noncanonical;
  int goodness;
  int absmq_score;
  int mapq_score;

  int translation_start;
  int translation_end;
  int translation_length;

  int relaastart;
  int relaaend;

#if 0
  int stage2_source;
  int stage2_indexsize;
  double stage2_diag_runtime;
  double stage2_align_runtime;
  double stage2_mapfraction;
  int stage2_maxconsecutive;
  double stage3_runtime;
#endif

  bool joinable_left_p;
  bool joinable_right_p;
};


bool
Stage3_chimera_left_p (T this) {
  return this->chimera_left_p;
}

bool
Stage3_chimera_right_p (T this) {
  return this->chimera_right_p;
}

bool
Stage3_watsonp (T this) {
  return this->watsonp;
}

int
Stage3_genestrand (T this) {
  return this->genestrand;
}

int
Stage3_cdna_direction (T this) {
  return this->cdna_direction;
}

int
Stage3_sensedir (T this) {
  return this->sensedir;
}


int
Stage3_straintype (T this) {
  return this->straintype;
}

int
Stage3_goodness (T this) {
  debug(printf("Overall goodness:\n"));
  debug(printf("  %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels => goodness %d\n",
	       this->matches,this->mismatches,this->qopens,this->qindels,this->topens,this->tindels,this->goodness));

  return this->goodness;
}

int
Stage3_absmq_score (T this) {
  return this->absmq_score;
}

int
Stage3_mapq_score (T this) {
  return this->mapq_score;
}

List_T
Stage3_pairs (T this) {
  return this->pairs;
}

struct Pair_T *
Stage3_pairarray (T this) {
  return this->pairarray;
}

int
Stage3_npairs (T this) {
  return this->npairs;
}

int
Stage3_matches (T this) {
  return this->matches;
}

int
Stage3_mismatches (T this) {
  return this->mismatches;
}

int
Stage3_indels (T this) {
  /* This should be consistent with the output from Pair_print_pathsummary */
  return this->qindels + this->tindels;
}


int
Stage3_querystart (T this) {
  return Pair_querypos(&(this->pairarray[0]));
}

int
Stage3_queryend (T this) {
  return Pair_querypos(&(this->pairarray[this->npairs-1]));
}

bool
Stage3_joinable_left_p (T this) {
  return this->joinable_left_p;
}

bool
Stage3_joinable_right_p (T this) {
  return this->joinable_right_p;
}

void
Stage3_clear_joinable (T this) {
  this->joinable_left_p = false;
  this->joinable_right_p = false;
  return;
}

void
Stage3_set_joinable_left (T this) {
  this->joinable_left_p = true;
  return;
}

void
Stage3_set_joinable_right (T this) {
  this->joinable_right_p = true;
  return;
}


void
Stage3_print_ends (T this) {
  Pair_print_ends(this->pairs);
  printf(" chimera_left: %d, chimera_right: %d",this->chimera_left_p,this->chimera_right_p);
  printf(" goodness: %d",this->goodness);
  printf("\n");
  return;
}


Chrnum_T
Stage3_chrnum (T this) {
  return this->chrnum;
}

Univcoord_T
Stage3_chroffset (T this) {
  return this->chroffset;
}

Univcoord_T
Stage3_chrhigh (T this) {
  return this->chrhigh;
}

Chrpos_T
Stage3_chrlength (T this) {
  return this->chrlength;
}

bool
Stage3_altloc_chr (Univcoord_T *alias_start, Univcoord_T *alias_end, T this) {
#if 0
  *alias_start = alias_starts[this->chrnum];
  *alias_end = alias_ends[this->chrnum];
  return altlocp[this->chrnum];
#else
  return false;
#endif
}


Chrpos_T
Stage3_chrstart (T this) {
  return Pair_genomepos(&(this->pairarray[0]));
}

Chrpos_T
Stage3_chrend (T this) {
  return Pair_genomepos(&(this->pairarray[this->npairs-1]));
}

Univcoord_T
Stage3_genomicstart (T this) {
  /* Should be chroffset + Pair_genomepos(start) */
  return this->genomicstart;
}

Univcoord_T
Stage3_genomicend (T this) {
  /* Should be chroffset + Pair_genomepos(end) */
  return this->genomicend;
}

void
Stage3_set_genomicend (T this, Univcoord_T genomicend) {
  this->genomicend = genomicend;
  return;
}

int
Stage3_circularpos (T this) {
  return this->circularpos;
}



int
Stage3_translation_start (T this) {
  return this->translation_start;
}

int
Stage3_translation_end (T this) {
  return this->translation_end;
}


int
Stage3_domain (T this) {
  int querystart, queryend;

  querystart = Pair_querypos(&(this->pairarray[0]));
  queryend = Pair_querypos(&(this->pairarray[this->npairs-1]));

  return queryend - querystart + 1;
}


int
Stage3_largemargin (int *newstart, int *newend, T this, int queryntlength) {
  int leftmargin, rightmargin;
  int querystart, queryend;

  querystart = Pair_querypos(&(this->pairarray[0]));
  queryend = Pair_querypos(&(this->pairarray[this->npairs-1]));

  if ((leftmargin = querystart) < 0) {
    leftmargin = 0;
  }
  if ((rightmargin = queryntlength - queryend) < 0) {
    rightmargin = 0;
  }

  /* Return larger margin */
  *newstart = querystart;
  *newend = queryend + 1;
  if (leftmargin > rightmargin) {
    /* Trim left */
    return leftmargin;
  } else {
    return rightmargin;
  }
}


double
Stage3_fracidentity (T this) {
  int den;

  if ((den = this->matches + this->mismatches + this->qindels + this->tindels) == 0) {
    return 1.0;
  } else {
    return (double) this->matches/(double) den;
  }
}

Univcoord_T
Stage3_genomicpos (T this, int querypos, bool headp) {
  return this->chroffset + Pair_genomicpos(this->pairarray,this->npairs,querypos,headp);
}


int
Stage3_chimeric_goodness (int *matches1, int *matches2, T part1, T part2, int breakpoint) {
  int goodness1, goodness2, querystart, queryend;
  int unknowns1, mismatches1, qopens1, qindels1, topens1, tindels1, 
    ncanonical1, nsemicanonical1, nnoncanonical1;
  int unknowns2, mismatches2, qopens2, qindels2, topens2, tindels2, 
    ncanonical2, nsemicanonical2, nnoncanonical2;

  querystart = Pair_querypos(&(part1->pairarray[0]));
  debug10(printf("Chimeric goodness requested for part %d..%d\n",querystart+1,breakpoint));
  Pair_fracidentity_bounded(&(*matches1),&unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
			    &ncanonical1,&nsemicanonical1,&nnoncanonical1,
			    part1->pairarray,part1->npairs,part1->cdna_direction,
			    querystart,breakpoint);
  goodness1 = (*matches1) + MISMATCH*mismatches1 + QOPEN*qopens1 + QINDEL*qindels1 + TOPEN*topens1 + TINDEL*tindels1;
  debug10(printf("  %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels => %d\n",
		*matches1,mismatches1,qopens1,qindels1,topens1,tindels1,goodness1));

  queryend = Pair_querypos(&(part2->pairarray[part2->npairs-1]));
  debug10(printf("Chimeric goodness requested for part %d..%d\n",breakpoint+1,queryend+1));
  Pair_fracidentity_bounded(&(*matches2),&unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
			    &ncanonical2,&nsemicanonical2,&nnoncanonical2,
			    part2->pairarray,part2->npairs,part2->cdna_direction,
			    breakpoint,queryend);
  goodness2 = (*matches2) + MISMATCH*mismatches2 + QOPEN*qopens2 + QINDEL*qindels2 + TOPEN*topens2 + TINDEL*tindels2;
  debug10(printf("  %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels => %d\n",
		*matches2,mismatches2,qopens2,qindels2,topens2,tindels2,goodness2));

  return goodness1 + goodness2;
}


Chrpos_T
Stage3_genomiclength (T this) {
  if (this->genomicstart < this->genomicend) {
    return this->genomicend - this->genomicstart + 1U;
  } else {
    return this->genomicstart - this->genomicend + 1U;
  }
}


bool
Stage3_passes_filter (T this, double min_trimmed_coverage, double min_identity) {
  int den;

  if (this->trimmed_coverage < min_trimmed_coverage) {
    return false;
  } else if ((den = this->matches + this->mismatches + this->qindels + this->tindels) == 0) {
    /* fracidentity is 1.0 */
    return true;
  } else if ((double) this->matches/(double) den < min_identity) {
    return false;
  } else {
    return true;
  }
}

bool
Stage3_passes_filter_chimera (Chimera_T chimera, double min_trimmed_coverage, double min_identity) {
  int den;
  Stage3_T from, to;

  from = Chimera_left_part(chimera);
  to = Chimera_right_part(chimera);

  if (from->trimmed_coverage + to->trimmed_coverage < min_trimmed_coverage) {
    return false;
  } else if ((den = from->matches + from->mismatches + from->qindels + from->tindels +
	      to->matches + to->mismatches + to->qindels + to->tindels) == 0) {
    /* fracidentity is 1.0 */
    return true;
  } else if ((double) (from->matches + to->matches)/(double) den < min_identity) {
    return false;
  } else {
    return true;
  }
}



int
Stage3_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  Chrpos_T x_genomiclength, y_genomiclength;

  if (x->chimera_right_p == true && y->chimera_right_p == false) {
    return -1;
  } else if (y->chimera_right_p == true && x->chimera_right_p == false) {
    return +1;
  } else if (x->chimera_left_p == true && y->chimera_left_p == false) {
    return -1;
  } else if (y->chimera_left_p == true && x->chimera_left_p == false) {
    return +1;
  } else if (x->goodness > y->goodness) {
    return -1;
  } else if (y->goodness > x->goodness) {
    return +1;

    /* If we can achieve same goodness with fewer pairs, then it is a better alignment */
  } else if (x->npairs < y->npairs) {
    return -1;
  } else if (y->npairs < x->npairs) {
    return +1;

    /* If we can achieve same goodness with more matches, then it is a better alignment */
  } else if (x->matches > y->matches) {
    return -1;
  } else if (y->matches > x->matches) {
    return +1;

  } else if (x->straintype < y->straintype) {
    return -1;
  } else if (y->straintype < x->straintype) {
    return +1;
  } else {
    x_genomiclength = Stage3_genomiclength(x);
    y_genomiclength = Stage3_genomiclength(y);
    if (x_genomiclength < y_genomiclength) {
      return -1;
    } else if (y_genomiclength < x_genomiclength) {
      return +1;
    } else if (x->chrnum < y->chrnum) {
      return -1;
    } else if (y->chrnum < x->chrnum) {
      return +1;
    } else if (x->genomicstart < y->genomicstart) {
      return -1;
    } else if (y->genomicstart < x->genomicstart) {
      return +1;
    } else {
      return 0;
    }
  }
}


int
Stage3_position_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  int querypos1, querypos2;

  if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (y->genomicstart < x->genomicstart) {
    return +1;
  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (y->genomicend < x->genomicend) {
    return +1;
  } else {
    querypos1 = Pair_querypos(&(x->pairarray[0]));
    querypos2 = Pair_querypos(&(y->pairarray[0]));
    if (querypos1 < querypos2) {
      return -1;
    } else if (querypos2 < querypos1) {
      return +1;
    } else {
      querypos1 = Pair_querypos(&(x->pairarray[x->npairs-1]));
      querypos2 = Pair_querypos(&(y->pairarray[y->npairs-1]));
      if (querypos1 < querypos2) {
	return -1;
      } else if (querypos2 < querypos1) {
	return +1;
      } else {
	return 0;
      }
    }
  }
}
					  
int
Stage3_querystart_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  int x_querystart, y_querystart;

  x_querystart = Pair_querypos(&(x->pairarray[0]));
  y_querystart = Pair_querypos(&(y->pairarray[0]));

  if (x_querystart < y_querystart) {
    return -1;
  } else if (y_querystart < x_querystart) {
    return +1;
  } else {
    return 0;
  }
}

int
Stage3_queryend_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  int x_queryend, y_queryend;

  x_queryend = Pair_querypos(&(x->pairarray[x->npairs-1]));
  y_queryend = Pair_querypos(&(y->pairarray[y->npairs-1]));

  if (x_queryend < y_queryend) {
    return -1;
  } else if (y_queryend < x_queryend) {
    return +1;
  } else {
    return 0;
  }
}

int
Stage3_chrnum_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->chrnum < y->chrnum) {
    return -1;
  } else if (y->chrnum < x->chrnum) {
    return +1;
  } else {
    return 0;
  }
}


int
Stage3_chrnum_querystart_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  int x_querystart, y_querystart, x_length, y_length;

  if (x->chrnum < y->chrnum) {
    return -1;
  } else if (y->chrnum < x->chrnum) {
    return +1;
  } else {
    x_querystart = Pair_querypos(&(x->pairarray[0]));
    y_querystart = Pair_querypos(&(y->pairarray[0]));

    if (x_querystart < y_querystart) {
      return -1;
    } else if (y_querystart < x_querystart) {
      return +1;
    } else {
      /* Put longer segments at the end so they supersede earlier chimeric matches */
      x_length = Pair_querypos(&(x->pairarray[x->npairs-1])) - x_querystart;
      y_length = Pair_querypos(&(y->pairarray[y->npairs-1])) - y_querystart;

      if (x_length < y_length) {
	return -1;
      } else if (y_length < x_length) {
	return +1;
      } else {
	return 0;
      }
    }
  }
}

int
Stage3_chrnum_queryend_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  int x_queryend, y_queryend, x_length, y_length;

  if (x->chrnum < y->chrnum) {
    return -1;
  } else if (y->chrnum < x->chrnum) {
    return +1;
  } else {
    x_queryend = Pair_querypos(&(x->pairarray[x->npairs-1]));
    y_queryend = Pair_querypos(&(y->pairarray[y->npairs-1]));

    if (x_queryend < y_queryend) {
      return -1;
    } else if (y_queryend < x_queryend) {
      return +1;
    } else {
      /* Put longer segments at the end so they supersede earlier chimeric matches */
      x_length = x_queryend - Pair_querypos(&(x->pairarray[0]));
      y_length = y_queryend - Pair_querypos(&(y->pairarray[0]));

      if (x_length < y_length) {
	return -1;
      } else if (y_length < x_length) {
	return +1;
      } else {
	return 0;
      }
    }
  }
}


int
Stage3_identity_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x < y) {
    return -1;
  } else if (x > y) {
    return +1;
  } else {
    return 0;
  }
}


bool
Stage3_overlap (T x, T y) {

  if (x->straintype != y->straintype) {
    return false;
  } else if (x->watsonp != y->watsonp) {
    return false;
  } else if (x->watsonp) {
    if (x->genomicstart >= y->genomicstart && x->genomicstart <= y->genomicend) {
      return true;
    } else if (y->genomicstart >= x->genomicstart && y->genomicstart <= x->genomicend) {
      return true;
    } else {
      return false;
    }
  } else {
    if (x->genomicstart >= y->genomicend && x->genomicstart <= y->genomicstart) {
      return true;
    } else if (y->genomicstart >= x->genomicend && y->genomicstart <= x->genomicstart) {
      return true;
    } else {
      return false;
    }
  }
}

/************************************************************************
 *   Gaps
 ************************************************************************/

/* Note: In going through pairs and path, we have two methods:

   1.  pairptr = path; (to save the pointer)
       path = Pairpool_pop(path,&pair);

       pairs = List_push_existing(pairs,pairptr);

   2.  pair = (Pair_T) path->first;

       (refer to path->rest, since we haven't popped path yet)
       pairs = List_transfer_one(pairs,&path);  (combines a push and pop)

    In the code below, we sometimes mix these, using method 2 for speed,
    and method 1 for clarity.
*/


#if 0
static List_T
check_gaps (List_T pairs, Pairpool_T pairpool) {
  List_T path = NULL, pairptr;
  Pair_T pair, leftpair, rightpair;
  int queryjump, genomejump;

  debug(printf("\nBeginning check of gaps\n"));
  debug(printf("length = %d\n",List_length(pairs)));
  debug(Pair_dump_list(pairs,true));

  if (pairs == NULL) {
    return (List_T) NULL;
  }

  pairptr = pairs;
  pairs = Pairpool_pop(pairs,&pair);
  if (pair->gapp == true) {
    fprintf(stderr,"Gap check error: Unexpected gap at start of pairs\n");
    debug(printf("Gap check error: Unexpected gap at start of pairs\n"));
#ifndef DEBUG
    Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
  } else {
#ifdef WASTE
    path = Pairpool_push_existing(NULL,pairpool,pair);
#else
    path = List_push_existing(NULL,pairptr);
#endif
  }

  while (pairs != NULL) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true) {
      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Observed a gap at %d..%d with queryjump = %d, genomejump = %d\n",
		   leftpair->querypos,rightpair->querypos,pair->queryjump,pair->genomejump));

      queryjump = rightpair->querypos - leftpair->querypos - 1;
      genomejump = rightpair->genomepos - leftpair->genomepos - 1;
      /* if (leftpair->cdna == ' ') queryjump++; -- For old dynamic programming */
      /* if (leftpair->genome == ' ') genomejump++; -- For old dynamic programming */

      if (pair->queryjump != queryjump) {
	if (rightpair->querypos >= HALFLEN && leftpair->querypos < HALFLEN) {
	  debug(printf("Accept queryjump for gap at %d..%d as probable skiplength.  It's %d, should be %d\n",
		       leftpair->querypos,rightpair->querypos,pair->queryjump,queryjump));
	} else {
	  debug(printf("Gap check error: Wrong queryjump for gap at %d..%d.  It's %d, should be %d\n",
		       leftpair->querypos,rightpair->querypos,pair->queryjump,queryjump));
#ifndef DEBUG
	  Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
	}
      }
      if (pair->genomejump != genomejump) {
	debug(printf("Gap check error: Wrong genomejump for gap at %d..%d.  It's %d, should be %d\n",
		     leftpair->querypos,rightpair->querypos,pair->genomejump,genomejump));
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
      }
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_push_existing(path,pairptr);
#endif

      /* Process another pair after gap */
      if (pairs == NULL) {
	fprintf(stderr,"Gap check error: Unexpected gap at end of pairs\n");
	debug(printf("Gap check error: Unexpected gap at end of pairs\n"));
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
      }
      pairptr = pairs;
      pairs = Pairpool_pop(pairs,&pair);
      if (pair->gapp == true) {
	fprintf(stderr,"Gap check error: Unexpected gap after gap\n");
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
      }
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_push_existing(path,pairptr);
#endif

    } else {
      /* Not a gap */
      leftpair = path->first;
      queryjump = pair->querypos - leftpair->querypos - 1;
      genomejump = pair->genomepos - leftpair->genomepos - 1;
      /* if (leftpair->cdna == ' ') queryjump++; -- For old dynamic programming */
      /* if (leftpair->genome == ' ') genomejump++; -- For old dynamic programming */

      if (queryjump <= 0 && genomejump <= 0) {
#ifdef WASTE
	path = Pairpool_push_existing(path,pairpool,pair);
#else
	path = List_push_existing(path,pairptr);
#endif
      } else if (queryjump == 0 && genomejump == 0) {
#ifdef WASTE
	path = Pairpool_push_existing(path,pairpool,pair);
#else
	path = List_push_existing(path,pairptr);
#endif
      } else {
	fprintf(stderr,"Gap check error: Unexpected missing gap at %d..%d\n",leftpair->querypos,pair->querypos);
	debug(printf("Gap check error: Unexpected missing gap at %d..%d\n",leftpair->querypos,pair->querypos));
	debug(printf("Gap check error: Pushing a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		     leftpair->querypos,pair->querypos,queryjump,genomejump));
	/* One place we need accurate queryjump and genomejump */
	path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump,
				       /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#ifdef WASTE
	path = Pairpool_push_existing(path,pairpool,pair);
#else
	path = List_push_existing(path,pairptr);
#endif
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
      }
    }
  }	

  debug(printf("Done with check of gaps\n\n"));

  return path;
}
#endif


static char complCode[128] = COMPLEMENT_LC;

static char
get_genomic_nt (char *g_alt, int genomicpos, Univcoord_T chroffset, Univcoord_T chrhigh,
		bool watsonp) {
  char c2, c2_alt;
  Univcoord_T pos;

  if (watsonp) {
    if ((pos = chroffset + genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';
      
    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';
      
    } else {
      debug7(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char_blocks(&(*g_alt),pos)));
      return Genome_get_char_blocks(&(*g_alt),pos);
    }

  } else {
    if ((pos = chrhigh - genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';
      
    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';
      
    } else {
      c2 = Genome_get_char_blocks(&c2_alt,pos);
    }
    debug7(printf("At %u, genomicnt is %c\n",
		  genomicpos,complCode[(int) c2]));
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}

#if 0
static char
get_genomic_nt_genomicseg (char *g_alt, Chrpos_T genomicpos,
			   Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			   char *genomicseg_ptr, Genome_T genome,
			   bool use_genomicseg_p) {
  char g, c2, c2_alt;
  Univcoord_T pos;

  assert(use_genomicseg_p == false);
  /* Need to allow genomicpos < 0 and genomicpos >= genomiclength for iteration on finding chimeras */
  if (use_genomicseg_p == true && genomicpos < 0) {
    *g_alt = '*';
    return '*';

  } else if (use_genomicseg_p == true && genomicpos >= genomiclength) {
    *g_alt = '*';
    return '*';

  } else if (use_genomicseg_p) {
    debug7(printf("At %u, genomicnt is %c\n",genomicpos,genomicseg_ptr[genomicpos]));
    g = *g_alt = genomicseg_ptr[genomicpos];
    return g;

  } else if (watsonp) {
    if ((pos = chroffset + genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';
      
    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';
      
#if 0
    } else if (genome) {
      debug7(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char(genome,pos)));
      return Genome_get_char(genome,pos);
#endif

    } else {
      debug7(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char_blocks(pos)));
      return Genome_get_char_blocks(&(*g_alt),pos);
    }

  } else {
    if ((pos = chrhigh - genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      return '*';
      
    } else if (pos >= chrhigh) {
      return '*';
      
#if 0
    } else if (genome) {
      c2 = Genome_get_char(genome,pos);
#endif

    } else {
      c2 = Genome_get_char_blocks(&c2_alt,pos);
    }
    debug7(printf("At %u, genomicnt is %c\n",
		  genomicpos,complCode[(int) c2]));
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}
#endif


#if 0
static char *
get_genomic_seg (Chrpos_T genomicpos, Univcoord_T chroffset, Univcoord_T chrhigh,
		 int length, bool watsonp,
		 char *genomicseg_ptr, Genome_T genome,
		 bool use_genomicseg_p) {
  char *segment;

  if (use_genomicseg_p) {
    debug7(printf("At %u, genomicseg is %.*s\n",
		  genomicpos,length,genomicseg_ptr[genomicpos]));
    return &(genomicseg_ptr[genomicpos]);

  } else if (watsonp) {
    segment = (char *) CALLOC(length+1,sizeof(char));
    if (genome) {
      Genome_fill_buffer_simple(genome,chroffset + genomicpos,length,segment);
      debug7(printf("At %u, genomicseg is %s\n",genomicpos,segment));
      return segment;
    } else {
      Genome_fill_buffer_blocks(chroffset + genomicpos,length,segment);
      debug7(printf("At %u, genomicseg is %s\n",genomicpos,segment));
      return segment;
    }

  } else {
    if (genome) {
      Genome_fill_buffer_simple(genome,chrhigh - genomicpos,length,segment);
    } else {
      Genome_fill_buffer_blocks(chrhigh - genomicpos,length,segment);
    }
    make_complement_inplace(segment,length);
    debug7(printf("At %u, genomicnt is %s\n",genomicpos,segment));
    return segment;
  }
}
#endif


/* For use by stage3.c procedures */
static List_T
insert_gapholders (List_T pairs, char *queryseq_ptr, char *queryuc_ptr,
		   Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
		   int genestrand, Pairpool_T pairpool, bool finalp) {
  List_T path = NULL;
  Pair_T pair, leftpair, gappair = NULL;
  int queryjump, genomejump, i;
  bool firstp = true;
  char comp, c, g, g_alt;

  /* Remove all existing gaps */
  debug(printf("Beginning deletion/insertion of gaps\n"));

  /* Discard old gap(s) and indels */
  while (pairs != NULL) {
    /* pairptr = pairs; */
    /* pairs = Pairpool_pop(pairs,&pair); */
    pair = (Pair_T) pairs->first;
    if (pair->knowngapp == true) {
      /* Keep known introns */
      debug(printf("Keeping a known intron gap with queryjump = %d, genomejump = %d\n",
		   pair->queryjump,pair->genomejump));
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif

    } else if (pair->gapp == true) {
      debug(printf("Removing a gap with queryjump = %d, genomejump = %d\n",
		   pair->queryjump,pair->genomejump));
      pairs = Pairpool_pop(pairs,&pair);

    } else if (pair->cdna == ' ' || pair->genome == ' ') {
      /* New */
      debug(printf("Removing an indel\n"));
      pairs = Pairpool_pop(pairs,&pair);

    } else {
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif
    }
  }

  pairs = List_reverse(path);
  path = (List_T) NULL;

  if (pairs != NULL) {
    /* pairptr = pairs; */
    /* pairs = Pairpool_pop(pairs,&pair); */
    pair = (Pair_T) pairs->first;
#ifdef WASTE
    path = Pairpool_push_existing(path,pairpool,pair);
#else
    path = List_transfer_one(path,&pairs);
#endif
    leftpair = pair;
  }

  while (pairs != NULL) {
    /* pairptr = pairs; */
    /* pairs = Pairpool_pop(pairs,&pair); */

    /* queryjump = pair->querypos - leftpair->querypos - 1; */
    /* genomejump = pair->genomepos - leftpair->genomepos - 1; */
    /* if (leftpair->cdna == ' ') queryjump++; -- For old dynamic programming */
    /* if (leftpair->genome == ' ') genomejump++; -- For old dynamic programming */

#if 1
    pair = (Pair_T) pairs->first;
    queryjump = pair->querypos - leftpair->querypos - 1;
    genomejump = pair->genomepos - leftpair->genomepos - 1;
#else
    /* Needed when we did not remove indels */
    p = pairs;
    while (p != NULL && ((Pair_T) p->first)->cdna == ' ') {
      p = p->rest;
    }
    if (p == NULL) {
      queryjump = 0;
    } else {
      queryjump = ((Pair_T) p->first)->querypos - leftpair->querypos - 1;
    }

    p = pairs;
    while (p != NULL && ((Pair_T) p->first)->genome == ' ') {
      p = p->rest;
    }
    if (p == NULL) {
      genomejump = 0;
    } else {
      genomejump = ((Pair_T) p->first)->genomepos - leftpair->genomepos - 1;
    }
#endif

    pair = (Pair_T) pairs->first;
    if (pair->knowngapp == true) {
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif

    } else if (leftpair->knowngapp == true) {
      /* Ignore queryjump and genomejump information of gap pair */
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif
    
    } else if (queryjump <= 0 && genomejump <= 0) {
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif

    } else if (finalp == true && queryjump == genomejump) {
      /* Fill gap with nucleotides */
      debug(printf("Filling a gap with nucleotides at %d..%d because of queryjump %d == genomejump %d\n",
		   leftpair->querypos,pair->querypos,queryjump,genomejump));
      for (i = 1; i <= queryjump; i++) {
	c = queryuc_ptr[leftpair->querypos+i];
	g = get_genomic_nt(&g_alt,leftpair->genomepos+i,chroffset,chrhigh,watsonp);
	/* It is possible for a gap with c == g to occur in the middle of a repetitive oligo, such as poly-A */
	if (Dynprog_consistent_p(c,g,g_alt,genestrand) == true) {
	  comp = MATCH_COMP;
#ifdef PMAP
	} else if (Dynprog_consistent_p(c,g,g_alt) == true) {
	  comp = AMBIGUOUS_COMP;
#endif
	} else {
	  comp = MISMATCH_COMP;
	}
	debug(printf(" => query %c, genomic %c\n",queryseq_ptr[leftpair->querypos+i],g));

	path = Pairpool_push(path,pairpool,leftpair->querypos+i,leftpair->genomepos+i,queryseq_ptr[leftpair->querypos+i],
			     comp,g,g_alt,/*dynprogindex*/0);
      }
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif

    } else if (queryjump == 1 && genomejump == 1) {
      /* Handle a single mismatch by a simple fill */
      c = queryuc_ptr[leftpair->querypos+1];
      g = get_genomic_nt(&g_alt,leftpair->genomepos+1,chroffset,chrhigh,watsonp);
      /* It is possible for a gap with c == g to occur in the middle of a repetitive oligo, such as poly-A */
      if (Dynprog_consistent_p(c,g,g_alt,genestrand) == true) {
	comp = MATCH_COMP;
#ifdef PMAP
      } else if (Dynprog_consistent_p(c,g,g_alt) == true) {
	comp = AMBIGUOUS_COMP;
#endif
      } else {
	comp = MISMATCH_COMP;
      }
      debug(printf("Filling a gap at %d..%d because of queryjump = %d, genomejump = %d => query %c, genomic %c\n",
		   leftpair->querypos,pair->querypos,queryjump,genomejump,queryseq_ptr[leftpair->querypos+1],g));
      path = Pairpool_push(path,pairpool,leftpair->querypos+1,leftpair->genomepos+1,queryseq_ptr[leftpair->querypos+1],
			   comp,g,g_alt,/*dynprogindex*/0);
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif

    } else {
      /* Insert new gap.  Need accurate queryjump and genomejump */
      debug(printf("Inserting a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		   leftpair->querypos,pair->querypos,queryjump,genomejump));
      debug(printf("queryjump %d = pair->querypos %d - leftpair->querypos %d - 1\n",queryjump,pair->querypos,leftpair->querypos));
      debug(printf("genomejump %d = pair->genomepos %u - leftpair->genomepos %u - 1\n",genomejump,pair->genomepos,leftpair->genomepos));

      path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump,
				     /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      gappair = (Pair_T) path->first;
      if (firstp == true) {
	gappair->end_intron_p = true;
	firstp = false;
      }
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_transfer_one(path,&pairs);
#endif
    }

    leftpair = pair;
  }

  if (gappair != NULL) {
    gappair->end_intron_p = true;
  }
  debug(printf("Ending deletion/insertion of gaps\n"));

  /* debug(Pair_dump_list(path,true)); */
  return path;
}



/* Should call before peel_rightward and peel_leftward, so we don't
   run into gaps that are really indels */
static List_T
assign_gap_types (List_T path, int cdna_direction, bool watsonp, char *queryseq_ptr,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		  Pairpool_T pairpool) {
  List_T pairs = NULL, pairptr;
  Pair_T pair, leftpair, rightpair;
  Univcoord_T splicesitepos;
  int queryjump, genomejump, leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, curquerypos,
    introntype, intronlength, genomicpos;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt, c2, c2_alt;

  debug(printf("\n** Starting assign_gap_types\n"));
  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;
    if (pair->gapp == false) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pairs == NULL) {
      /* Discard initial gap */
      debug7(printf("Discard initial gap\n"));
      path = Pairpool_pop(path,&pair);

    } else if (path->rest == NULL) {
      /* Discard terminal gap */
      debug7(printf("Discard terminal gap\n"));
      path = Pairpool_pop(path,&pair);

    } else {
      queryjump = pair->queryjump;
      genomejump = pair->genomejump;
      debug7(printf("  Gap has queryjump %d, genomejump %d\n",queryjump,genomejump));

      if (queryjump == 0 && genomejump == 0) {
	debug7(printf("  Gap is a non-gap\n"));
	/* Discard the gap pair */
	path = Pairpool_pop(path,&pair);

      } else if (genomejump == 0) {
	debug7(printf("  Gap is a cDNA insertion, so replacing it with indels\n"));
	/* pair->comp = INDEL_COMP; */

	/* Discard the gap pair */
	path = Pairpool_pop(path,&pair);

	leftpair = path->first;
	rightpair = pairs->first;
	leftquerypos = leftpair->querypos;
	/* if (leftpair->cdna == ' ') leftquerypos--; -- For old dynamic programming */
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;

	debug7(printf("leftquerypos = %d, rightquerypos = %d\n",leftquerypos,rightquerypos));
	for (curquerypos = rightquerypos - 1; curquerypos > leftquerypos; --curquerypos) {
	  debug7(printf("  pushing indel at %d\n",curquerypos));
	  pairs = Pairpool_push(pairs,pairpool,curquerypos,rightgenomepos,
				queryseq_ptr[curquerypos],INDEL_COMP,/*genome*/' ',/*genomealt*/' ',
				/*dynprogindex*/0);
	}

      } else if (queryjump > 0 /* || stage3debug > NO_STAGE3DEBUG */) {
	debug7(printf("  Gap is a dual break\n"));
	pair->comp = DUALBREAK_COMP;
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_transfer_one(pairs,&path);
#endif

      } else {
	debug7(printf("Gap is an intron\n"));

	pairptr = path;		/* save */
	path = Pairpool_pop(path,&pair);

	leftpair = path->first;
	rightpair = pairs->first;

	leftquerypos = leftpair->querypos;
	leftgenomepos = leftpair->genomepos;
	/* if (leftpair->cdna == ' ') leftquerypos--; -- For old dynamic programming */
	/* if (leftpair->genome == ' ') leftgenomepos--; -- For old dynamic programming */
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;

	pair->queryjump = rightquerypos - leftquerypos - 1;
	pair->genomejump = rightgenomepos - leftgenomepos - 1;

	left1 = get_genomic_nt(&left1_alt,leftgenomepos+1,chroffset,chrhigh,watsonp);
	left2 = get_genomic_nt(&left2_alt,leftgenomepos+2,chroffset,chrhigh,watsonp);
	right2 = get_genomic_nt(&right2_alt,rightgenomepos-2,chroffset,chrhigh,watsonp);
	right1 = get_genomic_nt(&right1_alt,rightgenomepos-1,chroffset,chrhigh,watsonp);
	debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
	introntype = Intron_type(left1,left2,right2,right1,
				 left1_alt,left2_alt,right2_alt,right1_alt,
				 cdna_direction);
	debug7(printf("  Introntype at %u..%u is %s (cdna_direction %d)\n",
		      leftgenomepos,rightgenomepos,Intron_type_string(introntype),cdna_direction));

	intronlength = rightgenomepos - leftgenomepos - 1;
	if (intronlength < min_intronlength) {
	  debug7(printf("  Gap is too short to be an intron (intronlength %d).  Replacing with pairs from %d downto %d\n",
			intronlength,rightgenomepos-1,leftgenomepos+1));
	  for (genomicpos = rightgenomepos - 1; genomicpos > leftgenomepos; --genomicpos) {
	    c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
	    pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',/*comp*/SHORTGAP_COMP,c2,c2_alt,
				  /*dynprogindex*/0);
	  }
	  debug7(printf("  Gap is a short gap with queryjump %d, genomejump %d, so discarding the gap pair\n",queryjump,genomejump));
	  /* Discard the gap */

	} else if (cdna_direction > 0) {
	  pair->introntype = introntype;
	  switch (introntype) {
	  case GTAG_FWD: pair->comp = FWD_CANONICAL_INTRON_COMP; break;
	  case GCAG_FWD: pair->comp = FWD_GCAG_INTRON_COMP; break;
	  case ATAC_FWD: pair->comp = FWD_ATAC_INTRON_COMP; break;
	  case NONINTRON: pair->comp = NONINTRON_COMP; break;
	  default: 
	    printf("Unexpected intron type %d\n",introntype);
	    fprintf(stderr,"Unexpected intron type %d\n",introntype);
	    abort();
	  }
	  debug7(printf("  Gap is a fwd intron (intronlength %d), now of type %c\n",intronlength,pair->comp));

	  if (watsonp == true) {
	    splicesitepos = leftgenomepos + 1;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/+1)) {
	      debug12(printf("1. donor at splicesitepos %u is known\n",splicesitepos));
	      pair->donor_prob = 1.0;
	    } else {
	      pair->donor_prob = Maxent_hr_donor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("1. donor at splicesitepos %u has prob %f\n",splicesitepos,pair->donor_prob));
	    }

	    splicesitepos = rightgenomepos;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/+1)) {
	      debug12(printf("2. acceptor at splicesitepos %u is known\n",splicesitepos));
	      pair->acceptor_prob = 1.0;
	    } else {
	      pair->acceptor_prob = Maxent_hr_acceptor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("2. acceptor at splicesitepos %u has prob %f\n",splicesitepos,pair->acceptor_prob));
	    }

	  } else {
	    splicesitepos = (chrhigh - chroffset) - leftgenomepos;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/-1)) {
	      debug12(printf("3. antidonor at splicesitepos %u is known\n",splicesitepos));
	      pair->donor_prob = 1.0;
	    } else {
	      pair->donor_prob = Maxent_hr_antidonor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("3. antidonor at splicesitepos %u has prob %f\n",splicesitepos,pair->donor_prob));
	    }

	    splicesitepos = (chrhigh - chroffset) - rightgenomepos + 1;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/-1)) {
	      debug12(printf("4. antiacceptor at splicesitepos %u is known\n",splicesitepos));
	      pair->acceptor_prob = 1.0;
	    } else {
	      pair->acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("4. antiacceptor at splicesitepos %u has prob %f\n",splicesitepos,pair->acceptor_prob));
	    }
	  }

	  /* Push the gap back on */
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif
	  
#ifndef PMAP
	} else if (cdna_direction < 0) {
	  pair->introntype = introntype;
	  switch (introntype) {
	  case ATAC_REV: pair->comp = REV_ATAC_INTRON_COMP; break;
	  case GCAG_REV: pair->comp = REV_GCAG_INTRON_COMP; break;
	  case GTAG_REV: pair->comp = REV_CANONICAL_INTRON_COMP; break;
	  case NONINTRON: pair->comp = NONINTRON_COMP; break;
	  default: 
	    printf("Unexpected intron type %d\n",introntype);
	    fprintf(stderr,"Unexpected intron type %d\n",introntype);
	    abort();
	  }
	  debug7(printf("  Gap is a rev intron (intronlength %d), now of type %c\n",intronlength,pair->comp));

	  if (watsonp == true) {
	    splicesitepos = leftgenomepos + 1;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/-1)) {
	      debug12(printf("5. antiacceptor at splicesitepos %u is known\n",splicesitepos));
	      pair->acceptor_prob = 1.0;
	    } else {
	      pair->acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("5. antiacceptor at splicesitepos %u has prob %f\n",splicesitepos,pair->acceptor_prob));
	    }

	    splicesitepos = rightgenomepos;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/-1)) {
	      debug12(printf("6. antidonor at splicesitepos %u is known\n",splicesitepos));
	      pair->donor_prob = 1.0;
	    } else {
	      pair->donor_prob = Maxent_hr_antidonor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("6. antidonor at splicesitepos %u has prob %f\n",splicesitepos,pair->donor_prob));
	    }

	  } else {
	    splicesitepos = (chrhigh - chroffset) - leftgenomepos;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/+1)) {
	      debug12(printf("7. acceptor at splicesitepos %u is known\n",splicesitepos));
	      pair->acceptor_prob = 1.0;
	    } else {
	      pair->acceptor_prob = Maxent_hr_acceptor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("7. acceptor at splicesitepos %u has prob %f\n",splicesitepos,pair->acceptor_prob));
	    }

	    splicesitepos = (chrhigh - chroffset) - rightgenomepos + 1;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/+1)) {
	      debug12(printf("8. donor at splicesitepos %u is known\n",splicesitepos));
	      pair->donor_prob = 1.0;
	    } else {
	      pair->donor_prob = Maxent_hr_donor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("8. donor at splicesitepos %u has prob %f\n",splicesitepos,pair->donor_prob));
	    }
	  }

	  /* Push the gap back on */
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif

#endif	/* ifndef PMAP */

	} else {
	  /* cdna_direction == 0 */
	  pair->introntype = introntype;
	  switch (introntype) {
	  case GTAG_FWD: pair->comp = FWD_CANONICAL_INTRON_COMP; break;
	  case GCAG_FWD: pair->comp = FWD_GCAG_INTRON_COMP; break;
	  case ATAC_FWD: pair->comp = FWD_ATAC_INTRON_COMP; break;
#ifndef PMAP
	  case ATAC_REV: pair->comp = REV_ATAC_INTRON_COMP; break;
	  case GCAG_REV: pair->comp = REV_GCAG_INTRON_COMP; break;
	  case GTAG_REV: pair->comp = REV_CANONICAL_INTRON_COMP; break;
#endif
	  case NONINTRON: pair->comp = NONINTRON_COMP; break;
	  default: 
	    printf("Unexpected intron type %d\n",introntype);
	    fprintf(stderr,"Unexpected intron type %d\n",introntype);
	    abort();
	  }
	  pair->donor_prob = 0.0;
	  pair->acceptor_prob = 0.0;

	  /* Push the gap back on */
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif
	}
      }
    }
  }

  return pairs;
}



static List_T
assign_intron_probs (List_T path, int cdna_direction, bool watsonp,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     Pairpool_T pairpool) {
  List_T pairs = NULL, pairptr;
  Pair_T pair, leftpair, rightpair;
  Univcoord_T splicesitepos;
  int queryjump, genomejump, leftquerypos, leftgenomepos, rightquerypos, rightgenomepos,
    introntype, intronlength, genomicpos;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt, c2, c2_alt;

  debug(printf("\n** Starting assign_intron_probs with watsonp %d and cdna_direction %d\n",watsonp,cdna_direction));
  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;
    if (pair->gapp == false) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pairs == NULL) {
      /* Discard initial gap */
      path = Pairpool_pop(path,&pair);

    } else if (path->rest == NULL) {
      /* Discard terminal gap */
      path = Pairpool_pop(path,&pair);

    } else {
      queryjump = pair->queryjump;
      genomejump = pair->genomejump;

      if (queryjump == 0 && genomejump == 0) {
	debug7(printf("  Gap is a non-gap\n"));
	/* Discard the gap pair */
	path = Pairpool_pop(path,&pair);

      } else if (genomejump == 0) {
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_transfer_one(pairs,&path);
#endif

      } else if (queryjump > 0) {
	debug7(printf("  Gap is a dual break\n"));
	pair->comp = DUALBREAK_COMP;
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_transfer_one(pairs,&path);
#endif

      } else {
	debug7(printf("Gap is an intron\n"));

	pairptr = path;		/* save */
	path = Pairpool_pop(path,&pair);

	leftpair = path->first;
	rightpair = pairs->first;

	leftquerypos = leftpair->querypos;
	leftgenomepos = leftpair->genomepos;
	/* if (leftpair->cdna == ' ') leftquerypos--; -- For old dynamic programming */
	/* if (leftpair->genome == ' ') leftgenomepos--; -- For old dynamic programming */
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;

	pair->queryjump = rightquerypos - leftquerypos - 1;
	pair->genomejump = rightgenomepos - leftgenomepos - 1;

	left1 = get_genomic_nt(&left1_alt,leftgenomepos+1,chroffset,chrhigh,watsonp);
	left2 = get_genomic_nt(&left2_alt,leftgenomepos+2,chroffset,chrhigh,watsonp);
	right2 = get_genomic_nt(&right2_alt,rightgenomepos-2,chroffset,chrhigh,watsonp);
	right1 = get_genomic_nt(&right1_alt,rightgenomepos-1,chroffset,chrhigh,watsonp);
	debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
	introntype = Intron_type(left1,left2,right2,right1,
				 left1_alt,left2_alt,right2_alt,right1_alt,
				 cdna_direction);
	debug7(printf("  Introntype at %u..%u is %s (cdna_direction %d)\n",
		      leftgenomepos,rightgenomepos,Intron_type_string(introntype),cdna_direction));

	intronlength = rightgenomepos - leftgenomepos - 1;
	if (intronlength < min_intronlength) {
	  debug7(printf("  Gap is too short to be an intron (intronlength %d).  Replacing with pairs from %d downto %d\n",
			intronlength,rightgenomepos-1,leftgenomepos+1));
	  for (genomicpos = rightgenomepos - 1; genomicpos > leftgenomepos; --genomicpos) {
	    c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
	    pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',/*comp*/SHORTGAP_COMP,c2,c2_alt,
				  /*dynprogindex*/0);
	  }
	  debug7(printf("  Gap is a short gap with queryjump %d, genomejump %d, so discarding the gap pair\n",queryjump,genomejump));
	  /* Discard the gap */

	} else if (cdna_direction > 0) {
	  pair->introntype = introntype;
	  switch (introntype) {
	  case GTAG_FWD: pair->comp = FWD_CANONICAL_INTRON_COMP; break;
	  case GCAG_FWD: pair->comp = FWD_GCAG_INTRON_COMP; break;
	  case ATAC_FWD: pair->comp = FWD_ATAC_INTRON_COMP; break;
	  case NONINTRON: pair->comp = NONINTRON_COMP; break;
	  default: 
	    printf("Unexpected intron type %d\n",introntype);
	    fprintf(stderr,"Unexpected intron type %d\n",introntype);
	    abort();
	  }
	  debug7(printf("  Gap is a fwd intron (intronlength %d), now of type %c\n",intronlength,pair->comp));

	  if (watsonp == true) {
	    splicesitepos = leftgenomepos + 1;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/+1)) {
	      debug12(printf("1. donor at splicesitepos %u is known\n",splicesitepos));
	      pair->donor_prob = 1.0;
	    } else {
	      pair->donor_prob = Maxent_hr_donor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("1. donor at splicesitepos %u has prob %f\n",splicesitepos,pair->donor_prob));
	    }

	    splicesitepos = rightgenomepos;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/+1)) {
	      debug12(printf("2. acceptor at splicesitepos %u is known\n",splicesitepos));
	      pair->acceptor_prob = 1.0;
	    } else {
	      pair->acceptor_prob = Maxent_hr_acceptor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("2. acceptor at splicesitepos %u has prob %f\n",splicesitepos,pair->acceptor_prob));
	    }

	  } else {
	    splicesitepos = (chrhigh - chroffset) - leftgenomepos;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/-1)) {
	      debug12(printf("3. antidonor at splicesitepos %u is known\n",splicesitepos));
	      pair->donor_prob = 1.0;
	    } else {
	      pair->donor_prob = Maxent_hr_antidonor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("3. antidonor at splicesitepos %u has prob %f\n",splicesitepos,pair->donor_prob));
	    }

	    splicesitepos = (chrhigh - chroffset) - rightgenomepos + 1;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/-1)) {
	      debug12(printf("4. antiacceptor at splicesitepos %u is known\n",splicesitepos));
	      pair->acceptor_prob = 1.0;
	    } else {
	      pair->acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("4. antiacceptor at splicesitepos %u has prob %f\n",splicesitepos,pair->acceptor_prob));
	    }
	  }

	  /* Push the gap back on */
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif
	  
#ifndef PMAP
	} else if (cdna_direction < 0) {
	  pair->introntype = introntype;
	  switch (introntype) {
	  case ATAC_REV: pair->comp = REV_ATAC_INTRON_COMP; break;
	  case GCAG_REV: pair->comp = REV_GCAG_INTRON_COMP; break;
	  case GTAG_REV: pair->comp = REV_CANONICAL_INTRON_COMP; break;
	  case NONINTRON: pair->comp = NONINTRON_COMP; break;
	  default: 
	    printf("Unexpected intron type %d\n",introntype);
	    fprintf(stderr,"Unexpected intron type %d\n",introntype);
	    abort();
	  }
	  debug7(printf("  Gap is a rev intron (intronlength %d), now of type %c\n",intronlength,pair->comp));

	  if (watsonp == true) {
	    splicesitepos = leftgenomepos + 1;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/-1)) {
	      debug12(printf("5. antiacceptor at splicesitepos %u is known\n",splicesitepos));
	      pair->acceptor_prob = 1.0;
	    } else {
	      pair->acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("5. antiacceptor at splicesitepos %u has prob %f\n",splicesitepos,pair->acceptor_prob));
	    }

	    splicesitepos = rightgenomepos;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/-1)) {
	      debug12(printf("6. antidonor at splicesitepos %u is known\n",splicesitepos));
	      pair->donor_prob = 1.0;
	    } else {
	      pair->donor_prob = Maxent_hr_antidonor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("6. antidonor at splicesitepos %u has prob %f\n",splicesitepos,pair->donor_prob));
	    }

	  } else {
	    splicesitepos = (chrhigh - chroffset) - leftgenomepos;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/+1)) {
	      debug12(printf("7. acceptor at splicesitepos %u is known\n",splicesitepos));
	      pair->acceptor_prob = 1.0;
	    } else {
	      pair->acceptor_prob = Maxent_hr_acceptor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("7. acceptor at splicesitepos %u has prob %f\n",splicesitepos,pair->acceptor_prob));
	    }

	    splicesitepos = (chrhigh - chroffset) - rightgenomepos + 1;
	    if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
								      splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/+1)) {
	      debug12(printf("8. donor at splicesitepos %u is known\n",splicesitepos));
	      pair->donor_prob = 1.0;
	    } else {
	      pair->donor_prob = Maxent_hr_donor_prob(chroffset + splicesitepos,chroffset);
	      debug12(printf("8. donor at splicesitepos %u has prob %f\n",splicesitepos,pair->donor_prob));
	    }
	  }

	  /* Push the gap back on */
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif

#endif	/* ifndef PMAP */

	} else {
	  /* cdna_direction == 0 */
	  pair->introntype = introntype;
	  switch (introntype) {
	  case GTAG_FWD: pair->comp = FWD_CANONICAL_INTRON_COMP; break;
	  case GCAG_FWD: pair->comp = FWD_GCAG_INTRON_COMP; break;
	  case ATAC_FWD: pair->comp = FWD_ATAC_INTRON_COMP; break;
#ifndef PMAP
	  case ATAC_REV: pair->comp = REV_ATAC_INTRON_COMP; break;
	  case GCAG_REV: pair->comp = REV_GCAG_INTRON_COMP; break;
	  case GTAG_REV: pair->comp = REV_CANONICAL_INTRON_COMP; break;
#endif
	  case NONINTRON: pair->comp = NONINTRON_COMP; break;
	  default: 
	    printf("Unexpected intron type %d\n",introntype);
	    fprintf(stderr,"Unexpected intron type %d\n",introntype);
	    abort();
	  }
	  pair->donor_prob = 0.0;
	  pair->acceptor_prob = 0.0;

	  /* Push the gap back on */
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif
	}
      }
    }
  }

  return pairs;
}


/* Modeled after assign_gap_types */
static List_T
remove_indel_gaps (List_T path
#ifdef WASTE
		   , Pairpool_T pairpool
#endif
		   ) {
  List_T pairs = NULL, pairptr;
  Pair_T pair, leftpair, rightpair;
  int queryjump, genomejump, leftgenomepos, rightgenomepos, intronlength;

  debug(printf("\n** Starting remove_indel_gaps\n"));
  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;

    if (pair->gapp == false) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pairs == NULL) {
      /* Discard initial gap */
      path = Pairpool_pop(path,&pair);

    } else if (path->rest == NULL) {
      /* Discard terminal gap */
      path = Pairpool_pop(path,&pair);

    } else {
      queryjump = pair->queryjump;
      genomejump = pair->genomejump;

      if (queryjump == 0 && genomejump == 0) {
	debug7(printf("  Gap is a non-gap\n"));
	/* Discard the gap pair */
	path = Pairpool_pop(path,&pair);

      } else if (genomejump == 0) {
	debug7(printf("  Gap is a cDNA insertion\n"));
	/* pair->comp = INDEL_COMP; */
	/* Discard the gap pair */
	path = Pairpool_pop(path,&pair);

      } else if (queryjump > 0) {
	debug7(printf("  Gap is a dual break\n"));
	pair->comp = DUALBREAK_COMP;
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_transfer_one(pairs,&path);
#endif

      } else {
	debug7(printf("  Gap is an intron of type %c\n",pair->comp));

	pairptr = path;		/* save */
	path = Pairpool_pop(path,&pair);

	leftpair = path->first;
	rightpair = pairs->first;

	leftgenomepos = leftpair->genomepos;
	/* if (leftpair->genome == ' ') leftgenomepos--; -- For old dynamic programming */
	rightgenomepos = rightpair->genomepos;

	intronlength = rightgenomepos - leftgenomepos - 1;
	if (intronlength < min_intronlength) {
	  debug7(printf("  Gap is short (intronlength %d).  Adding pairs from %d downto %d\n",
			intronlength,rightgenomepos-1,leftgenomepos+1));
	  debug7(printf("  Gap is a short gap, so discarding the gap pair\n"));
	  /* Discard the gap pair */

	} else {
	  debug7(printf("  Gap is not short (intronlength %d)\n",intronlength));
	  /* Push the gap back on */
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif
	}
      }
    }
  }

  return pairs;
}



#ifdef PMAP
static List_T
undefine_nucleotides (char *queryseq_ptr, int querylength, List_T path, Pairpool_T pairpool, int width) {
  List_T pairs = NULL, pairptr;
  Pair_T pair, leftpair, rightpair;
  int leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, pos;

  debug(printf("\n** Starting undefine_nucleotides\n"));

  if (path != NULL) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
#ifdef WASTE
    pairs = Pairpool_push_existing(NULL,pairpool,pair);
#else
    pairs = List_push_existing(NULL,pairptr);
#endif
    rightquerypos = pair->querypos;
    rightgenomepos = pair->genomepos;
  }

  while (path != NULL) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == true) {
      leftpair = path->first;
      rightpair = pairs->first;

      leftquerypos = leftpair->querypos;
      leftgenomepos = leftpair->genomepos;
      /* if (leftpair->cdna == ' ') leftquerypos--; -- For old dynamic programming */
      /* if (leftpair->genome == ' ') leftgenomepos--; -- For old dynamic programming */

      rightquerypos = rightpair->querypos;
      rightgenomepos = rightpair->genomepos;

      debug(printf("Undefining around rightquerypos = %d and leftquerypos = %d\n",rightquerypos,leftquerypos));
      for (pos = rightquerypos; pos < rightquerypos + width && pos < querylength; pos++) {
	queryseq_ptr[pos] = BACKTRANSLATE_CHAR;
      }
      for (pos = leftquerypos; pos > leftquerypos - width && pos >= 0; --pos) {
	queryseq_ptr[pos] = BACKTRANSLATE_CHAR;
      }
    }
#ifdef WASTE
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
    pairs = List_push_existing(pairs,pairptr);
#endif
  }

  return pairs;
}  
#endif


static List_T
add_dualbreak (List_T pairs, char *queryseq_ptr, 
	       Univcoord_T chroffset, Univcoord_T chrhigh, int cdna_direction,
	       bool watsonp, Pair_T leftpair, Pair_T rightpair, Pairpool_T pairpool, int ngap) {
  int genomicpos, k;
  int leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, gapgenomepos, midpoint;
  int introntype;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  char c1, c2, c2_alt, comp;

  leftquerypos = leftpair->querypos;
  leftgenomepos = leftpair->genomepos;
  /* if (leftpair->cdna == ' ') leftquerypos--; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') leftgenomepos--; -- For old dynamic programming */
  rightquerypos = rightpair->querypos;
  rightgenomepos = rightpair->genomepos;

  /* Previously checked for genomicuc_ptr != NULL, but this does not
     work with second round of prepare_for_printing */
  left1 = get_genomic_nt(&left1_alt,leftgenomepos+1,chroffset,chrhigh,watsonp);
  left2 = get_genomic_nt(&left2_alt,leftgenomepos+2,chroffset,chrhigh,watsonp);
  right2 = get_genomic_nt(&right2_alt,rightgenomepos-2,chroffset,chrhigh,watsonp);
  right1 = get_genomic_nt(&right1_alt,rightgenomepos-1,chroffset,chrhigh,watsonp);
  
  debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
  introntype = Intron_type(left1,left2,right2,right1,
			   left1_alt,left2_alt,right2_alt,right1_alt,
			   cdna_direction);
  debug7(printf("  Introntype at %u..%u is %s (cdna_direction %d)\n",
		leftgenomepos,rightgenomepos,Intron_type_string(introntype),cdna_direction));
  switch (introntype) {
  case GTAG_FWD: comp = FWD_CANONICAL_INTRON_COMP; break;
  case GCAG_FWD: comp = FWD_GCAG_INTRON_COMP; break;
  case ATAC_FWD: comp = FWD_ATAC_INTRON_COMP; break;
#ifndef PMAP
  case ATAC_REV: comp = REV_ATAC_INTRON_COMP; break;
  case GCAG_REV: comp = REV_GCAG_INTRON_COMP; break;
  case GTAG_REV: comp = REV_CANONICAL_INTRON_COMP; break;
#endif
  case NONINTRON: comp = NONINTRON_COMP; break;
  default: 
    printf("Unexpected intron type %d\n",introntype);
    fprintf(stderr,"Unexpected intron type %d\n",introntype);
    abort();
  }
  /* End of check */


  /* queryjump = rightquerypos - leftquerypos - 1; */
  /* genomejump = rightgenomepos - leftgenomepos - 1; */

  if (rightgenomepos - leftgenomepos - 1 < ngap + ngap) {
    midpoint = (rightgenomepos + leftgenomepos) / 2;

    /* First insertion */
    for (genomicpos = rightgenomepos - 1; genomicpos >= midpoint; --genomicpos) {
      c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,
				     /*cdna*/' ',comp,introntype,c2,c2_alt,/*extraexonp*/true);
    }

    /* cDNA sequence */
    gapgenomepos = genomicpos + 1;
    for (k = rightquerypos - 1; k > leftquerypos; --k) {
#if 0				/* PMAP */
      c1 = Sequence_codon_char(queryaaseq_ptr[k/3],k%3);
#else
      c1 = queryseq_ptr[k];
#endif
      pairs = Pairpool_push_gapalign(pairs,pairpool,k,gapgenomepos,
				     c1,EXTRAEXON_COMP,/*introntype*/NONINTRON,c1,c1,
				     /*extraexonp*/true); /* Transfer cDNA char to genome */
    }

    /* Second insertion */
    for (genomicpos = midpoint - 1; genomicpos > leftgenomepos; --genomicpos) {
      c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push_gapalign(pairs,pairpool,leftquerypos,genomicpos,
				     /*cdna*/' ',comp,introntype,c2,c2_alt,/*extraexonp*/true);
    }

  } else {

    /* First insertion */
    for (k = 0, genomicpos = rightgenomepos - 1; k < ngap; k++, --genomicpos) {
      c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,
				     /*cdna*/' ',comp,introntype,c2,c2_alt,/*extraexonp*/true);
    }

    /* cDNA sequence */
    gapgenomepos = genomicpos + 1;
    for (k = rightquerypos - 1; k > leftquerypos; --k) {
#if 0				/* PMAP */
      c1 = Sequence_codon_char(queryaaseq_ptr[k/3],k%3);
#else
      c1 = queryseq_ptr[k];
#endif
      pairs = Pairpool_push_gapalign(pairs,pairpool,k,gapgenomepos,
				     c1,EXTRAEXON_COMP,/*introntype*/NONINTRON,c1,c1,
				     /*extraexonp*/true); /* Transfer cDNA char to genome */
    }

    /* Second insertion */
    genomicpos = leftgenomepos + ngap;
    for (k = 0; k < ngap; k++, --genomicpos) {
      c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push_gapalign(pairs,pairpool,leftquerypos,genomicpos,
				     /*cdna*/' ',comp,introntype,c2,c2_alt,/*extraexonp*/true);
    }
  }

  return pairs;
}


static List_T
add_intron (List_T pairs, Univcoord_T chroffset, Univcoord_T chrhigh,
	    Pair_T leftpair, Pair_T rightpair, char comp, int introntype, int ngap,
	    bool watsonp, Pairpool_T pairpool) {
  char c2, c2_alt;
  int intronlength, genomicpos;
  int leftgenomepos, rightquerypos, rightgenomepos, gapgenomepos;
  int i;

  leftgenomepos = leftpair->genomepos;
  /* if (leftpair->genome == ' ') leftgenomepos--; -- For old dynamic programming */
  rightquerypos = rightpair->querypos;
  rightgenomepos = rightpair->genomepos;

  intronlength = rightgenomepos - leftgenomepos - 1;

  debug7(printf("Adding gap of type %c of length %d\n",comp,intronlength));

#if 0
  /* Should not be necessary to fix introns at this point */
  if (cdna_direction >= 0) {
    switch (*comp) {
    case FWD_CANONICAL_INTRON_COMP: case FWD_GCAG_INTRON_COMP: case FWD_ATAC_INTRON_COMP: case NONINTRON: break;
    default: 
      debug7(printf("Unexpected intron comp %c.  Need to fix.\n",*comp));

      left1 = genomicuc_ptr[leftgenomepos+1];
      left2 = genomicuc_ptr[leftgenomepos+2];
      right2 = genomicuc_ptr[rightgenomepos-2];
      right1 = genomicuc_ptr[rightgenomepos-1];

      debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
      introntype = Intron_type(left1,left2,right2,right1,
			       left1_alt,left2_alt,right2_alt,right1_alt,
			       cdna_direction);
      debug7(printf("  Introntype at %u..%u is %s (cdna_direction %d)\n",
		    leftgenomepos,rightgenomepos,Intron_type_string(introntype),cdna_direction));
      switch (introntype) {
      case GTAG_FWD: *comp = FWD_CANONICAL_INTRON_COMP; break;
      case GCAG_FWD: *comp = FWD_GCAG_INTRON_COMP; break;
      case ATAC_FWD: *comp = FWD_ATAC_INTRON_COMP; break;
      case NONINTRON:
	intronlength = rightgenomepos - leftgenomepos - 1;
	if (intronlength < min_intronlength) {
	  *comp = SHORTGAP_COMP;	/* Will be printed as INDEL_COMP, but need to score as NONINTRON_COMP */
	} else {
	  *comp = NONINTRON_COMP;
	}
      }
    }
  } else {
    switch (*comp) {
    case REV_CANONICAL_INTRON_COMP: case REV_GCAG_INTRON_COMP: case REV_ATAC_INTRON_COMP: case NONINTRON: break;
    default: 
      debug7(printf("Unexpected intron comp %c.  Need to fix.\n",*comp));

      left1 = genomicuc_ptr[leftgenomepos+1];
      left2 = genomicuc_ptr[leftgenomepos+2];
      right2 = genomicuc_ptr[rightgenomepos-2];
      right1 = genomicuc_ptr[rightgenomepos-1];

      debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
      introntype = Intron_type(left1,left2,right2,right1,
			       left1_alt,left2_alt,right2_alt,right1_alt,
			       cdna_direction);
      debug7(printf("  Introntype at %u..%u is %s (cdna_direction %d)\n",
		    leftgenomepos,rightgenomepos,Intron_type_string(introntype),cdna_direction));
      switch (introntype) {
      case ATAC_REV: *comp = REV_ATAC_INTRON_COMP; break;
      case GCAG_REV: *comp = REV_GCAG_INTRON_COMP; break;
      case GTAG_REV: *comp = REV_CANONICAL_INTRON_COMP; break;
      case NONINTRON:
	intronlength = rightgenomepos - leftgenomepos - 1;
	if (intronlength < min_intronlength) {
	  *comp = SHORTGAP_COMP;	/* Will be printed as INDEL_COMP, but need to score as NONINTRON_COMP */
	} else {
	  *comp = NONINTRON_COMP;
	}
      }
    }
  }
#endif

  if (intronlength < ngap + ngap + 3) {
    for (i = 0, genomicpos = rightgenomepos - 1; i < intronlength; i++, --genomicpos) {
      c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,
				     /*cdna*/' ',comp,introntype,c2,c2_alt,/*extraexonp*/false);
    }
  } else {
    for (i = 0, genomicpos = rightgenomepos - 1; i < ngap; i++, --genomicpos) {
      c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,
				     /*cdna*/' ',comp,introntype,c2,c2_alt,/*extraexonp*/false);
      debug7(printf("Pushing %c at genomicpos %d\n",c2,genomicpos));
    }
    
    gapgenomepos = genomicpos + 1;
    pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,gapgenomepos,' ',INTRONGAP_COMP,/*introntype*/NONINTRON,
				   /*genome*/INTRONGAP_CHAR,/*genomealt*/INTRONGAP_CHAR,/*extraexonp*/false);
    pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,gapgenomepos,' ',INTRONGAP_COMP,/*introntype*/NONINTRON,
				   /*genome*/INTRONGAP_CHAR,/*genomealt*/INTRONGAP_CHAR,/*extraexonp*/false);
    pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,gapgenomepos,' ',INTRONGAP_COMP,/*introntype*/NONINTRON,
				   /*genome*/INTRONGAP_CHAR,/*genomealt*/INTRONGAP_CHAR,/*extraexonp*/false);

    genomicpos = leftgenomepos + ngap;
    for (i = ngap-1; i >= 0; --i, --genomicpos) {
      c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,
				     /*cdna*/' ',comp,introntype,c2,c2_alt,/*extraexonp*/false);
      debug7(printf("Pushing %c at genomicpos %d\n",c2,genomicpos));
    }
  }

  return pairs;
}



/************************************************************************
 *   Fix adjacent indels
 ************************************************************************/

/* Modeled after print_sam_forward in pair.c */
/* Handles indels next to gaps */
static List_T
fix_adjacent_indels (List_T pairs) {
  List_T path = NULL, pairptr;
  Pair_T this = NULL, pair;
  bool in_exon = false;
  int Mlength = 0, Ilength = 0, Dlength = 0;
  char last_token_type = ' ';
  int last_token_length = 0, i;

  debug4(printf("Starting fix_adjacent_indels: "));

  while (pairs != NULL) {
    this = (Pair_T) List_head(pairs);

    if (this->gapp) {
      if (in_exon == true) {

	if (Mlength > 0) {
	  last_token_type = 'M';
	  last_token_length = Mlength;
	  debug4(printf("%dM",Mlength));
	} else if (Ilength > 0) {
	  debug4(printf("%dI",Ilength));
	  if (last_token_type == 'I' || last_token_type == 'D') {
	    debug4(printf("fix_adjacent_indels found %d%c to %d%c\n",last_token_length,last_token_type,Ilength,'I'));
	    for (i = 0; i < last_token_length + Ilength; i++) {
	      path = Pairpool_pop(path,&pair);
	    }
	    last_token_type = 'I';
	    last_token_length = 0; /* Since we have already taken care of this */
	  } else {
	    last_token_type = 'I';
	    last_token_length = Ilength;
	  }
	} else if (Dlength > 0) {
	  debug4(printf("%dD",Dlength));
	  if (last_token_type == 'I' || last_token_type == 'D') {
	    debug4(printf("fix_adjacent_indels found %d%c to %d%c\n",last_token_length,last_token_type,Dlength,'D'));
	    for (i = 0; i < last_token_length + Dlength; i++) {
	      path = Pairpool_pop(path,&pair);
	    }
	    last_token_type = 'D';
	    last_token_length = 0; /* Since we have already taken care of this */
	  } else {
	    last_token_type = 'D';
	    last_token_length = Dlength;
	  }
	}

	Mlength = Ilength = Dlength = 0;
	in_exon = false;
      }

    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {

	if (last_token_type != ' ') {
	  /* Gap */
	  debug4(printf("?N"));
	  last_token_type = 'N';  /* Could potentially also be considered 'D' */
	  last_token_length = 0;

#if 0
	  query_gap = this->querypos - exon_queryend;
	  if (query_gap > 0) {
	    /* Dual gap.  Don't try to piece together.  */
	    debug4(printf("%dI",query_gap));
	    last_token_type = 'I';
	    last_token_length = query_gap;
	  }
#endif
	}

	in_exon = true;
      }

      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (this->genome == ' ') {
	  if (Mlength > 0) {
	    debug4(printf("%dM",Mlength));
	    last_token_type = 'M';
	    last_token_length = Mlength;
	    Mlength = 0;

	  } else if (Dlength > 0) {
	    /* unlikely */
	    debug4(printf("%dD",Dlength));
	    if (last_token_type == 'I' || last_token_type == 'D') {
	      debug4(printf("fix_adjacent_indels found %d%c to %d%c\n",last_token_length,last_token_type,Dlength,'D'));
	      for (i = 0; i < last_token_length + Dlength; i++) {
		path = Pairpool_pop(path,&pair);
	      }
	      last_token_type = 'D';
	      last_token_length = 0; /* Since we have already taken care of this */
	    } else {
	      last_token_type = 'D';
	      last_token_length = Dlength;
	      Dlength = 0;
	    }
	  }
	  Ilength++;

	} else if (this->cdna == ' ') {
	  if (Mlength > 0) {
	    debug4(printf("%dM",Mlength));
	    last_token_type = 'M';
	    last_token_length = Mlength;
	    Mlength = 0;

	  } else if (Ilength > 0) {
	    debug4(printf("%dI",Ilength));
	    if (last_token_type == 'I' || last_token_type == 'D') {
	      debug4(printf("fix_adjacent_indels found %d%c to %d%c\n",last_token_length,last_token_type,Ilength,'I'));
	      for (i = 0; i < last_token_length + Ilength; i++) {
		path = Pairpool_pop(path,&pair);
	      }
	      last_token_type = 'I';
	      last_token_length = 0; /* Since we have already taken care of this */
	    } else {
	      last_token_type = 'I';
	      last_token_length = Ilength;
	    }
	    Ilength = 0;
	  }
	  Dlength++;

	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	/* Count even if unknown base */

	if (Ilength > 0) {
	  debug4(printf("%dI",Ilength));
	  if (last_token_type == 'I' || last_token_type == 'D') {
	    debug4(printf("fix_adjacent_indels found %d%c to %d%c\n",last_token_length,last_token_type,Ilength,'I'));
	    for (i = 0; i < last_token_length + Ilength; i++) {
	      path = Pairpool_pop(path,&pair);
	    }
	    last_token_type = 'I';
	    last_token_length = 0; /* Since we have already taken care of this */
	  } else {
	    last_token_type = 'I';
	    last_token_length = Ilength;
	  }
	  Ilength = 0;

	} else if (Dlength > 0) {
	  debug4(printf("%dD",Dlength));
	  if (last_token_type == 'I' || last_token_type == 'D') {
	    debug4(printf("fix_adjacent_indels found %d%c to %d%c\n",last_token_length,last_token_type,Dlength,'D'));
	    for (i = 0; i < last_token_length + Dlength; i++) {
	      path = Pairpool_pop(path,&pair);
	    }
	    last_token_type = 'D';
	    last_token_length = 0; /* Since we have already taken care of this */
	  } else {
	    last_token_type = 'D';
	    last_token_length = Dlength;
	  }
	  Dlength = 0;
	}

	Mlength++;
      }
    }

    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
    path = Pairpool_push_existing(path,pairpool,pair);
#else
    path = List_push_existing(path,pairptr);
#endif
  }

  if (Mlength > 0) {
    debug4(printf("%dM",Mlength));
    /* last_token_type = 'M'; */
    /* last_token_length = Mlength; */
  } else if (Ilength > 0) {
    debug4(printf("%dI",Ilength));
    if (last_token_type == 'I' || last_token_type == 'D') {
      debug4(printf("fix_adjacent_indels found %d%c to %d%c\n",last_token_length,last_token_type,Ilength,'I'));
      for (i = 0; i < last_token_length + Ilength; i++) {
	path = Pairpool_pop(path,&pair);
      }
      /* last_token_type = 'I'; */
      /* last_token_length = 0; */ /* Since we have already taken care of this */
    } else {
      /* last_token_type = 'I'; */
      /* last_token_length = Ilength; */
    }
  } else if (Dlength > 0) {
    debug4(printf("%dD",Dlength));
    if (last_token_type == 'I' || last_token_type == 'D') {
      debug4(printf("fix_adjacent_indels found %d%c to %d%c\n",last_token_length,last_token_type,Dlength,'D'));
      for (i = 0; i < last_token_length + Dlength; i++) {
	path = Pairpool_pop(path,&pair);
      }
      /* last_token_type = 'D'; */
      /* last_token_length = 0; */ /* Since we have already taken care of this */
    } else {
      /* last_token_type = 'D'; */
      /* last_token_length = Dlength; */
    }
  }

  debug4(printf("\n"));

  return path;
}



#define NORMAL_STATE 0
#define INSERTION_STATE +1
#define DELETION_STATE -1


#if 0
static List_T
remove_adjacent_ins_del (bool *foundp, List_T pairs) {
  List_T path = NULL, pairptr;
  Pair_T this, pair;
  int state = NORMAL_STATE;

  *foundp = false;
  while (pairs != NULL) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&this);

    if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
      if (this->genome == ' ') {
	if (state == NORMAL_STATE) {
	  path = List_push_existing(path,pairptr);
	  state = INSERTION_STATE;

	} else if (state == INSERTION_STATE) {
	  /* Do nothing */
	  path = List_push_existing(path,pairptr);

	} else if (state == DELETION_STATE) {
	  /* Switch from insertion to deletion */
	  /* Remove past insertion */
	  while (path != NULL &&
		 (((Pair_T) path->first)->comp == INDEL_COMP || ((Pair_T) path->first)->comp == SHORTGAP_COMP) &&
		 ((Pair_T) path->first)->cdna == ' ') {
	    path = Pairpool_pop(path,&pair);
	  }
	  /* Remove future deletion */
	  while (pairs != NULL &&
		 (((Pair_T) pairs->first)->comp == INDEL_COMP || ((Pair_T) pairs->first)->comp == SHORTGAP_COMP) &&
		 ((Pair_T) pairs->first)->genome == ' ') {
	    pairs = Pairpool_pop(pairs,&pair);
	  }
	  *foundp = true;
	  state = NORMAL_STATE;

	} else {
	  abort();
	}

      } else if (this->cdna == ' ') {
	if (state == NORMAL_STATE) {
	  path = List_push_existing(path,pairptr);
	  state = DELETION_STATE;

	} else if (state == DELETION_STATE) {
	  /* Do nothing */
	  path = List_push_existing(path,pairptr);

	} else if (state == INSERTION_STATE) {
	  /* Switch from deletion to insertion */
	  /* Remove past deletion */
	  while (path != NULL &&
		 (((Pair_T) path->first)->comp == INDEL_COMP || ((Pair_T) path->first)->comp == SHORTGAP_COMP) &&
		 ((Pair_T) path->first)->genome == ' ') {
	    path = Pairpool_pop(path,&pair);
	  }
	  /* Remove future insertion */
	  while (pairs != NULL &&
		 (((Pair_T) pairs->first)->comp == INDEL_COMP || ((Pair_T) pairs->first)->comp == SHORTGAP_COMP) &&
		 ((Pair_T) pairs->first)->cdna == ' ') {
	    pairs = Pairpool_pop(pairs,&pair);
	  }
	  *foundp = true;
	  state = NORMAL_STATE;

	} else {
	  abort();
	}

      } else {
	abort();
      }

    } else {
      path = List_push_existing(path,pairptr);
      state = NORMAL_STATE;
    }
  }

  return path;
}
#endif



/************************************************************************
 *   Chop (trimming within end exons)
 ************************************************************************/

/* Called only by GMAP, because nucleotide matches in PMAP have several ambiguous matches. */
static List_T
clean_path_end3 (List_T path, int ambig_end_length_3) {
  Pair_T lastpair;

  debug(printf("Starting clean_path_end3\n"));
  debug(Pair_dump_list(path,true));
  if (ambig_end_length_3 == 0) {
    /* Remove any remaining nonmatches, gaps, or indels at 3' end */
    if (path != NULL) {
      lastpair = path->first;
      while (lastpair->gapp || (lastpair->comp != MATCH_COMP && lastpair->comp != DYNPROG_MATCH_COMP /*&& lastpair->comp != AMBIGUOUS_COMP*/)) {
	debug(printf("Removing nonmatch at 3' end: "));
	debug(Pair_dump_one(lastpair,/*zerobasedp*/true));
	debug(printf("\n"));
	path = Pairpool_pop(path,&lastpair);
	if (path == NULL) {
	  return NULL;
	} else {
	  lastpair = path->first;
	}
      }
    }

#ifdef PMAP
    while (path != NULL) {
      lastpair = path->first;
      if (lastpair->querypos % 3 == 2) {
	debug(printf("Ending clean_path_end3\n"));
	debug(Pair_dump_list(path,true));
	return path;
      } else {
	debug(printf("PMAP popping querypos %d to get to codon boundary\n",lastpair->querypos));
	path = Pairpool_pop(path,&lastpair);
      }
    }
#endif
  }

  debug(printf("Ending clean_path_end3\n"));
  debug(Pair_dump_list(path,true));
  return path;
}


static List_T
clean_pairs_end5 (List_T pairs, int ambig_end_length_5) {
  Pair_T firstpair;

  debug(printf("Starting clean_pairs_end5\n"));
  if (ambig_end_length_5 == 0) {
    /* Remove any remaining nonmatches, gaps, or indels at 5' end */
    if (pairs != NULL) {
      firstpair = pairs->first;
      while (firstpair->gapp || (firstpair->comp != MATCH_COMP && firstpair->comp != DYNPROG_MATCH_COMP /*&& firstpair->comp != AMBIGUOUS_COMP*/)) {
	debug(printf("Removing nonmatch at 5' end: "));
	debug(Pair_dump_one(firstpair,/*zerobasedp*/true));
	debug(printf("\n"));
	pairs = Pairpool_pop(pairs,&firstpair);
	if (pairs == NULL) {
	  return NULL;
	} else {
	  firstpair = pairs->first;
	}
      }
    }

#ifdef PMAP
    while (pairs != NULL) {
      firstpair = pairs->first;
      if (firstpair->querypos % 3 == 0) {
	return pairs;
      } else {
	debug(printf("PMAP popping querypos %d to get to codon boundary\n",firstpair->querypos));
	pairs = Pairpool_pop(pairs,&firstpair);
      }
    }
#endif
  }

  debug(printf("Ending clean_pairs_end5\n"));
  return pairs;
}


static List_T
clip_path_end3_chromosomal_bounds (List_T path, Univcoord_T chroffset, Univcoord_T chrhigh) {
  Pair_T pair;

  while (path != NULL && ((Pair_T) path->first)->genomepos >= chrhigh - chroffset) {
    /* Will clean gaps also, but that is okay */
    path = Pairpool_pop(path,&pair);
    debug(printf("Clipping a pair with %u\n",pair->genomepos));
  }

  return path;
}

static List_T
clip_pairs_end5_chromosomal_bounds (List_T pairs) {
  Pair_T pair;

  /* genomepos of 0 corresponds to chrhigh (probably should change this to be chrhigh - 1 in the future) */
  while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < 1) {
    /* Will clean gaps also, but that is okay */
    pairs = Pairpool_pop(pairs,&pair);
    debug(printf("Clipping a pair with %u\n",pair->genomepos));
  }

  return pairs;
}


/* Called only by GMAP, because nucleotide matches in PMAP have several ambiguous matches. */
static List_T
clean_path_end3_gap_indels (List_T path) {
  Pair_T lastpair;

  debug(printf("Starting clean_path_end3_gap_indels\n"));
  debug(Pair_dump_list(path,true));

  /* Remove any remaining gap/indels at 3' end, which can happen rarely */
  if (path != NULL) {
    lastpair = path->first;
    while (lastpair->gapp == true || lastpair->comp == INDEL_COMP || lastpair->comp == SHORTGAP_COMP) {
      debug(printf("Removing gap/indel at 3' end: "));
      debug(Pair_dump_one(lastpair,/*zerobasedp*/true));
      debug(printf("\n"));
      path = Pairpool_pop(path,&lastpair);
      if (path == NULL) {
	return NULL;
      } else {
	lastpair = path->first;
      }
    }
  }

#ifdef PMAP
  while (path != NULL) {
    lastpair = path->first;
    if (lastpair->querypos % 3 == 2) {
      debug(printf("Ending clean_path_end3_gap_indels\n"));
      debug(Pair_dump_list(path,true));
      return path;
    } else {
      debug(printf("PMAP popping querypos %d to get to codon boundary\n",lastpair->querypos));
      path = Pairpool_pop(path,&lastpair);
    }
  }
#endif

  debug(printf("Ending clean_path_end3_gap_indels\n"));
  debug(Pair_dump_list(path,true));
  return path;
}


static List_T
clean_pairs_end5_gap_indels (List_T pairs) {
  Pair_T firstpair;

  debug(printf("Starting clean_pairs_end5_gap_indels\n"));
  /* Remove any remaining gap/indels at 5' end, which can happen rarely */
  if (pairs != NULL) {
    firstpair = pairs->first;
    while (firstpair->gapp == true || firstpair->comp == INDEL_COMP || firstpair->comp == SHORTGAP_COMP) {
      debug(printf("Removing gap/indel at 5' end: "));
      debug(Pair_dump_one(firstpair,/*zerobasedp*/true));
      debug(printf("\n"));
      pairs = Pairpool_pop(pairs,&firstpair);
      if (pairs == NULL) {
	return NULL;
      } else {
	firstpair = pairs->first;
      }
    }
  }

#ifdef PMAP
  while (pairs != NULL) {
    firstpair = pairs->first;
    if (firstpair->querypos % 3 == 0) {
      return pairs;
    } else {
      debug(printf("PMAP popping querypos %d to get to codon boundary\n",firstpair->querypos));
      pairs = Pairpool_pop(pairs,&firstpair);
    }
  }
#endif

  debug(printf("Ending clean_pairs_end5_gap_indels\n"));
  return pairs;
}


/* Cleans to any gapp found within 20 bp of the breakpoint */
static List_T
clean_end_chimera (List_T end) {
  Pair_T lastpair;
  List_T peeled = NULL;
  int n = 0;

  debug10(printf("Starting clean_path_end_chimera\n"));
  while (end != NULL && n < 20) {
    lastpair = end->first;
    peeled = List_transfer_one(peeled,&end);
    if (lastpair->gapp == true) {
      debug10(printf("Cleaning end at a gapp\n"));
      peeled = (List_T) NULL;
    } else if (lastpair->comp == INDEL_COMP) {
      debug10(printf("Cleaning end at an indel\n"));
      peeled = (List_T) NULL;
    } else if (lastpair->comp == SHORTGAP_COMP) {
      debug10(printf("Cleaning end at an indel\n"));
      peeled = (List_T) NULL;
    } else {
      n++;
    }
  }

  end = Pairpool_transfer(end,peeled);
  debug10(printf("Ending clean_path_end_chimera\n"));
  return end;
}


#if 0
static List_T
chop_ends_by_changepoint (List_T pairs
#ifdef WASTE
			  , Pairpool_T pairpool
#endif
			  ) {
  List_T path;
  Pair_T pair;
  int *matchscores;
  int nmatches, ntotal, nmatches_left, ntotal_left, nmatches_right, ntotal_right;
  int left_edge, right_edge, length, i;
  int side;
  double theta;
  bool chop_left_p = false, chop_right_p = false;

  if (pairs == NULL) {
    return (List_T) NULL;
  } else {
    matchscores = Pair_matchscores_list(&nmatches,&ntotal,&length,pairs);
    debug18(printf("Overall, %d matches/%d total\n",nmatches,ntotal));
  }

  left_edge = Changepoint_left(&nmatches_left,&ntotal_left,matchscores,length);
  right_edge = Changepoint_right(&nmatches_right,&ntotal_right,matchscores,length);

  debug18(printf("At left edge %d (in 0..%d), %d matches/%d total\n",left_edge,List_length(pairs),nmatches_left,ntotal_left));
  debug18(printf("At right edge %d (in 0..%d), %d matches/%d total\n",right_edge,List_length(pairs),nmatches_right,ntotal_right));

  if (right_edge <= left_edge) {
    debug18(printf("Edges cross.  Need to select one.\n"));
    /* Need to select one side to chop. */
    if (ntotal_left == 0 || ntotal - ntotal_left <= 0) {
      side = +1;		/* chop right side */
    } else if (ntotal_right == 0 || ntotal - ntotal_right <= 0) {
      side = -1;		/* chop left side */
    } else {

#if 0
      theta = (double) (nmatches - nmatches_left)/(double) (ntotal - ntotal_left);
      /* Don't have artificially high expectations for theta, e.g., 1.00 */
      theta = theta - THETA_SLACK;
      if (theta < 0.10) {
	/* Protect against negative values */
	theta = 0.10;
      }
      debug18(printf("Testing on left: Pbinom(%d,%d,%f) = %g\n",
		    nmatches_left,ntotal_left,theta,Pbinom(nmatches_left,ntotal_left,theta)));
      pbinom_left = Pbinom(nmatches_left,ntotal_left,theta);

      theta = (double) (nmatches - nmatches_right)/(double) (ntotal - ntotal_right);
      /* Don't have artificially high expectations for theta, e.g., 1.00 */
      theta = theta - THETA_SLACK;
      if (theta < 0.10) {
	/* Protect against negative values */
	theta = 0.10;
      }
      
      debug18(printf("Testing on right: Pbinom(%d,%d,%f) = %g\n",
		    nmatches_right,ntotal_right,theta,Pbinom(nmatches_right,ntotal_right,theta)));
      pbinom_right = Pbinom(nmatches_right,ntotal_right,theta);
      
      if (pbinom_left < pbinom_right) {
	if (pbinom_left < TRIM_END_PVALUE) {
	  side = -1;		/* chop left side */
	} else {
	  side = 0;
	}
      } else if (pbinom_right < pbinom_left) {
	if (pbinom_right < TRIM_END_PVALUE) {
	  side = +1;		/* chop right side */
	} else {
	  side = 0;
	}
      } else {
	side = 0;
      }
#else
      /* Pick shortest side */
      if (ntotal_left < ntotal_right) {
	debug18(printf("left side is shorter\n"));
	side = -1;		/* chop left side */
      } else {
	debug18(printf("right side is shorter\n"));
	side = +1;		/* chop right side */
      }
#endif

    }

    if (side == -1) {
      debug18(printf("Chopping %d on left.\n",left_edge));
      for (i = 0; i < left_edge; i++) {
	pairs = Pairpool_pop(pairs,&pair);
      }
      chop_left_p = true;
    } else if (side == +1) {
      debug18(printf("Chopping %d - %d on right.\n",length,right_edge));
      path = List_reverse(pairs);
      for (i = 0; i < length - right_edge; i++) {
	path = Pairpool_pop(path,&pair);
      }

      pairs = (List_T) NULL;
#ifdef WASTE
      while (path) {
	path = Pairpool_pop(path,&pair);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }
#else
      pairs = Pairpool_transfer(pairs,path);
#endif
      chop_right_p = true;
    }

  } else {
    if (ntotal_left == 0) {
      path = List_reverse(pairs);
    } else if (ntotal - ntotal_left <= 0) {
      path = List_reverse(pairs);
    } else {
      theta = (double) (nmatches - nmatches_left)/(double) (ntotal - ntotal_left);
      /* Don't have artificially high expectations for theta, e.g., 1.00 */
      theta = theta - THETA_SLACK;
      if (theta < 0.10) {
	/* Protect against negative values */
	theta = 0.10;
      }

      debug18(printf("Testing on left: Pbinom(%d,%d,%f) = %g\n",
		    nmatches_left,ntotal_left,theta,Pbinom(nmatches_left,ntotal_left,theta)));
      if (Pbinom(nmatches_left,ntotal_left,theta) > TRIM_END_PVALUE) {
	path = List_reverse(pairs);
      } else {
	debug18(printf("Chopping %d on left.\n",left_edge));
	for (i = 0; i < left_edge; i++) {
	  pairs = Pairpool_pop(pairs,&pair);
	}
	path = (List_T) NULL;
#ifdef WASTE
	while (pairs) {
	  pairs = Pairpool_pop(pairs,&pair);
	  path = Pairpool_push_existing(path,pairpool,pair);
	}
#else
	path = Pairpool_transfer(path,pairs);
#endif
	debug18(printf("path is now length %d\n",List_length(path)));
	chop_left_p = true;
      }
    }

    if (ntotal_right == 0) {
      pairs = List_reverse(path);
    } else if (ntotal - ntotal_right <= 0) {
      pairs = List_reverse(path);
    } else {
      theta = (double) (nmatches - nmatches_right)/(double) (ntotal - ntotal_right);
      /* Don't have artificially high expectations for theta, e.g., 1.00 */
      theta = theta - THETA_SLACK;
      if (theta < 0.10) {
	/* Protect against negative values */
	theta = 0.10;
      }

      debug18(printf("Testing on right: Pbinom(%d,%d,%f) = %g\n",
		    nmatches_right,ntotal_right,theta,Pbinom(nmatches_right,ntotal_right,theta)));
      if (Pbinom(nmatches_right,ntotal_right,theta) > TRIM_END_PVALUE) {
	pairs = List_reverse(path);
      } else {
	debug18(printf("Chopping %d - %d on right.\n",length,right_edge));
	for (i = 0; i < length - right_edge; i++) {
	  path = Pairpool_pop(path,&pair);
	}
	pairs = (List_T) NULL;
#ifdef WASTE
	while (path) {
	  path = Pairpool_pop(path,&pair);
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
	}
#else
	pairs = Pairpool_transfer(pairs,path);
#endif
	chop_right_p = true;
      }
    }
  }
  
  FREE(matchscores);

  debug18(printf("Returning alignment of length %d\n",List_length(pairs)));

  return pairs;
}
#endif


#if 0
/* pairs -> pairs */
static List_T
trim_short_end_exons (bool *trim5p, bool *trim3p, List_T pairs, Pairpool_T pairpool, int minendexon) {
  List_T path, exon, pairptr;
  Pair_T pair;
  int exon_nmatches;

  debug18(printf("Starting trim_short_end5_exons\n"));
  debug18(Pair_dump_list(pairs,true));

  /* Handle first exon */
  if (pairs == NULL) {
    *trim5p = *trim3p = false;
    return (List_T) NULL;
  } else {
    pair = pairs->first;
  }

  exon = (List_T) NULL;
  exon_nmatches = 0;
  while (pairs != NULL && !pair->gapp) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
    if (pair->gapp == false && pair->comp != MISMATCH_COMP && pair->comp != INDEL_COMP && pair->comp != SHORTGAP_COMP) {
      exon_nmatches++;
    }
  }

  if (exon_nmatches >= minendexon) {
    debug18(printf("Keeping first exon of length %d\n",exon_nmatches));
    path = exon;		/* exon already has the gap */
    *trim5p = false;
  } else if (exon_nmatches == 0) {
    debug18(printf("Trimming first exon of length %d.  firstpair must be a gap.\n",exon_nmatches));
    pairs = Pairpool_pop(pairs,&pair); /* discard gap */
    assert(pair->gapp == true);
    path = (List_T) NULL;
    *trim5p = false;
  } else {
    debug18(printf("Trimming first exon of length %d\n",exon_nmatches));
    path = (List_T) NULL;
    *trim5p = true;
  }

#ifdef WASTE
  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    path = Pairpool_push_existing(path,pairpool,pair);

  }
#else
  path = Pairpool_transfer(path,pairs);
#endif

  /* Handle last exon */
  if (path == NULL) {
    *trim5p = *trim3p = false;
    return (List_T) NULL;
  } else {
    pair = path->first;
  }

  exon = (List_T) NULL;
  exon_nmatches = 0;
  while (path != NULL && !pair->gapp) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
    if (pair->gapp == false && pair->comp != MISMATCH_COMP && pair->comp != INDEL_COMP && pair->comp != SHORTGAP_COMP) {
      exon_nmatches++;
    }
  }

  if (exon_nmatches >= minendexon) {
    debug18(printf("Keeping last exon of length %d\n",exon_nmatches));
    pairs = exon;		/* exon already has the gap */
    *trim3p = false;
  } else if (exon_nmatches == 0) {
    debug18(printf("Trimming last exon of length %d.  firstpair must be a gap.\n",exon_nmatches));
    path = Pairpool_pop(path,&pair); /* discard gap */
    assert(pair->gapp == true);
    pairs = (List_T) NULL;
    *trim3p = false;
  } else {
    debug18(printf("Trimming last exon of length %d\n",exon_nmatches));
    pairs = (List_T) NULL;
    *trim3p = true;
  }

#ifdef WASTE
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
  }
#else
  pairs = Pairpool_transfer(pairs,path);
#endif

  debug18(printf("End of trim_short_end_exons: length = %d\n",List_length(pairs)));
  debug18(Pair_dump_list(pairs,true));
  return pairs;
}
#endif


#if 0
/* pairs -> path */
static List_T
trim_short_end5_exons (bool *trim5p, List_T pairs,
#ifdef WASTE
		       Pairpool_T pairpool,
#endif
		       int minendexon) {
  List_T path, exon, pairptr;
  Pair_T pair;
  int exon_nmatches, exon_nmismatches;

  debug18(printf("Starting trim_short_end5_exons\n"));
  debug18(Pair_dump_list(pairs,true));

  /* Handle first exon */
  if (pairs == NULL) {
    *trim5p = false;
    return (List_T) NULL;
  } else {
    pair = pairs->first;
  }

  exon = (List_T) NULL;
  exon_nmatches = exon_nmismatches = 0;
  while (pairs != NULL && !pair->gapp) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
    if (pair->gapp == true) {
      /* Skip */
    } else if (pair->comp == MISMATCH_COMP || pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
      exon_nmismatches++;
    } else {
      exon_nmatches++;
    }
  }

  if (exon_nmatches - exon_nmismatches >= minendexon) {
    debug18(printf("Keeping first exon of length %d\n",exon_nmatches));
    path = exon;		/* exon already has the gap */
    *trim5p = false;
  } else if (exon_nmatches == 0) {
    debug18(printf("Trimming first exon of length %d.  firstpair must be a gap.\n",exon_nmatches));
    pairs = Pairpool_pop(pairs,&pair); /* discard gap */
    path = (List_T) NULL;
    *trim5p = false;
  } else {
    debug18(printf("Trimming first exon of length %d\n",exon_nmatches));
    path = (List_T) NULL;
    *trim5p = true;
  }

#ifdef WASTE
  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    path = Pairpool_push_existing(path,pairpool,pair);

  }
#else
  path = Pairpool_transfer(path,pairs);
#endif

  debug18(printf("End of trim_short_end_exons: length = %d\n",List_length(pairs)));
  debug18(Pair_dump_list(pairs,true));
  return path;
}
#endif


#if 0
/* path -> pairs */
static List_T
trim_short_end3_exons (bool *trim3p, List_T path,
#ifdef WASTE
		       Pairpool_T pairpool,
#endif
		       int minendexon) {
  List_T pairs, exon, pairptr;
  Pair_T pair;
  int exon_nmatches, exon_nmismatches;

  debug18(printf("Starting trim_short_end3_exons\n"));
  debug18(Pair_dump_list(path,true));

  /* Handle last exon */
  if (path == NULL) {
    *trim3p = false;
    return (List_T) NULL;
  } else {
    pair = path->first;
  }

  exon = (List_T) NULL;
  exon_nmatches = exon_nmismatches = 0;
  while (path != NULL && !pair->gapp) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
    if (pair->gapp == true) {
      /* Skip */
    } else if (pair->comp == MISMATCH_COMP || pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
      exon_nmismatches++;
    } else {
      exon_nmatches++;
    }
  }

  if (exon_nmatches - exon_nmismatches >= minendexon) {
    debug18(printf("Keeping last exon of length %d\n",exon_nmatches));
    pairs = exon;		/* exon already has the gap */
    *trim3p = false;
  } else if (exon_nmatches == 0) {
    debug18(printf("Trimming last exon of length %d.  firstpair must be a gap.\n",exon_nmatches));
    path = Pairpool_pop(path,&pair); /* discard gap */
    pairs = (List_T) NULL;
    *trim3p = false;
  } else {
    debug18(printf("Trimming last exon of length %d\n",exon_nmatches));
    pairs = (List_T) NULL;
    *trim3p = true;
  }

#ifdef WASTE
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
  }
#else
  pairs = Pairpool_transfer(pairs,path);
#endif

  debug18(printf("End of trim_short_end_exons: length = %d\n",List_length(pairs)));
  debug18(Pair_dump_list(pairs,true));
  return pairs;
}
#endif



#ifdef STRICT
/* We want pair->queryjump == 0 or pair->genomejump == 0 */
static bool
dual_break_p (List_T pairs) {
  Pair_T pair;

  while (pairs != NULL) {
    pair = (Pair_T) pairs->first;
    /* This used to fail when we used UNKNOWNJUMP for gaps */
    if (pair->gapp == true && (pair->queryjump != 0 && pair->genomejump != 0)) {
      return true;
    }
    pairs = pairs->rest;
  }

  return false;
}
#endif


/* We want pair->queryjump >= 0 and pair->genomejump >= 0 */
static bool
negative_break_p (List_T pairs) {
  Pair_T pair;

  while (pairs != NULL) {
    pair = (Pair_T) pairs->first;
    /* This used to fail when we used UNKNOWNJUMP for gaps */
    if (pair->gapp == true && (pair->queryjump < 0 || pair->genomejump < 0)) {
      return true;
    }
    pairs = pairs->rest;
  }

  return false;
}



#if 0
/* We want pair->queryjump == 0 or pair->genomejump == 0 */
static int
dual_break_distance_from_end (int *npairs, int *totaljump, List_T pairs) {
  Pair_T pair;
  int nmatches, nmismatches;

  /* Handle first pair */
  if (pairs == NULL) {
    *totaljump = 0;
    return 0;
  } else {
    pair = pairs->first;
  }

  nmatches = nmismatches = 0;
  *npairs = 0;
  while (pairs != NULL && (pair->gapp == false || (pair->queryjump == 0 || pair->genomejump == 0))) {
    if (pair->gapp == true) {
      /* Skip */
    } else if (pair->comp == MISMATCH_COMP || pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
      nmismatches++;
    } else {
      nmatches++;
    }

    pairs = pairs->rest;
    if (pairs != NULL) {
      pair = (Pair_T) pairs->first;
    }
    *npairs += 1;
  }

  /* This used to fail when we used UNKNOWNJUMP for gaps */
  if (pair->gapp == true && (pair->queryjump != 0 || pair->genomejump != 0)) {
    *npairs += 1;		/* trim gap */
    *totaljump = DUALBREAK_QUERYJUMP_FACTOR * pair->queryjump;
  } else {
    *totaljump = 0;
  }

  return nmatches - nmismatches;
}
#endif


#ifdef STRICT
/* We want pair->queryjump == 0 or pair->genomejump == 0 */
static void
dual_break_count_npairs (int *npairs5, int *npairs3, List_T pairs) {
  Pair_T pair;
  int npairs = 0;

  assert(pairs != NULL);	/* We previously checked for this */
  pair = pairs->first;

  npairs = 0;
  while (pairs != NULL && (pair->gapp == false || (pair->queryjump == 0 || pair->genomejump == 0))) {
    pairs = pairs->rest;
    if (pairs != NULL) {
      pair = (Pair_T) pairs->first;
    }
    npairs += 1;
  }
  *npairs5 = npairs;

  /* Skip the gap */
  pairs = pairs->rest;
  npairs += 1;

  /* Count the rest */
  while (pairs != NULL) {
    pairs = pairs->rest;
    npairs += 1;
  }
  *npairs3 = npairs - (*npairs5) - 1 /* for the gap */;

  return;
}
#endif


/* We want pair->queryjump >= 0 and pair->genomejump >= 0 */
static void
negative_break_count_npairs (int *npairs5, int *npairs3, List_T pairs) {
  Pair_T pair;
  int npairs = 0;

  assert(pairs != NULL);	/* We previously checked for this */
  pair = pairs->first;

  npairs = 0;
  while (pairs != NULL && (pair->gapp == false || (pair->queryjump >= 0 && pair->genomejump >= 0))) {
    pairs = pairs->rest;
    if (pairs != NULL) {
      pair = (Pair_T) pairs->first;
    }
    npairs += 1;
  }
  *npairs5 = npairs;

  /* Skip the gap */
  pairs = pairs->rest;
  npairs += 1;

  /* Count the rest */
  while (pairs != NULL) {
    pairs = pairs->rest;
    npairs += 1;
  }
  *npairs3 = npairs - (*npairs5) - 1 /* for the gap */;

  return;
}


static List_T
trim_npairs (List_T pairs, int npairs) {
  int i;
  Pair_T pair;

  for (i = 0; i < npairs; i++) {
    pairs = Pairpool_pop(pairs,&pair);
  }
  return pairs;
}



#if 0
static bool
enough_matches (int matches, int genomejump
#if 0
		, double donor_prob, double acceptor_prob
#endif
		) {
#if 1
  if (genomejump > 100000) {
    return (matches >= 10) ? true : false;
  } else if (genomejump > 32000) {
    return (matches >= 9) ? true : false;
  } else if (genomejump > 8000) {
    return (matches >= 8) ? true : false;
  } else if (genomejump > 2000) {
    return (matches >= 7) ? true : false;
  } else {
    return (matches >= 6) ? true : false;
  }
#else
  double prob, prob_threshold;

  prob = 1 - pow(1.0-pow(4.0,(double) -matches),(double) genomejump);
  debug3(printf("Probability of exon of length %d with intron of length %d is %g\n",
		 matches,genomejump,prob));

#if 0
  prob_threshold = 1.0 - (1.0 - donor_prob)*(1.0 - acceptor_prob);
  debug3(printf("  Comparing with probability of splice %f and %f => %f\n",
		 donor_prob,acceptor_prob,prob_threshold));
#endif

  if (prob < 0.10) {
    return true;
  } else {
    return false;
  }
#endif
}
#endif


static bool
canonicalp (bool knowngapp, char comp, double donor_prob, double acceptor_prob, int cdna_direction) {

  if (knowngapp) {
    return true;

  } else if (donor_prob < 0.9 || acceptor_prob < 0.9) {
    return false;
    
  } else if (cdna_direction > 0) {
    if (comp == FWD_CANONICAL_INTRON_COMP || comp == FWD_GCAG_INTRON_COMP || comp == FWD_ATAC_INTRON_COMP) {
      return true;
    } else {
      return false;
    }
  } else if (cdna_direction < 0) {
    if (comp == REV_CANONICAL_INTRON_COMP || comp == REV_GCAG_INTRON_COMP || comp == REV_ATAC_INTRON_COMP) {
      return true;
    } else {
      return false;
    }
  } else {
#if 0
    /* Too much freedom.  Also, ambig_splicetypes depend on cdna_direction to be known. */
    if (comp == FWD_CANONICAL_INTRON_COMP || comp == FWD_GCAG_INTRON_COMP || comp == FWD_ATAC_INTRON_COMP ||
	comp == REV_CANONICAL_INTRON_COMP || comp == REV_GCAG_INTRON_COMP || comp == REV_ATAC_INTRON_COMP) {
      return true;
    } else {
      return false;
    }
#else
    return false;
#endif
  }

}


/* Copied from stage1hr.c.  A slightly different version is used in splice.c */
static int
sufficient_splice_prob_local (int nmatches, int nmismatches, double distal_spliceprob,
			      double medial_spliceprob) {
  debug3(printf("Checking for sufficient splice prob, based on %d matches, %d mismatches, support %d, distal spliceprob %f, and medial spliceprob %f\n",
		nmatches,nmismatches,support,distal_spliceprob,medial_spliceprob));
  nmatches -= 2*nmismatches;
  if (nmatches < 0) {
    return (int) false;
  } else if (nmatches < 7) {
    return (distal_spliceprob > 0.95 && medial_spliceprob > 0.90);
  } else if (nmatches < 11) {
    return (distal_spliceprob > 0.90 && medial_spliceprob > 0.85);
  } else if (nmatches < 15) {
    return (distal_spliceprob > 0.85 && medial_spliceprob > 0.80);
  } else if (nmatches < 19) {
    return (distal_spliceprob > 0.50 /*&& medial_spliceprob > 0.50*/);
  } else {
    return (int) true;
  }
}

static bool
sufficient_splice_prob_strict (double distal_spliceprob, double medial_spliceprob) {
  debug3(printf("Checking for sufficient splice prob, based on spliceprob %f, and medial spliceprob %f\n",
		distal_spliceprob,medial_spliceprob));
  if (distal_spliceprob > 0.95 && medial_spliceprob > 0.90) {
    return true;
  } else {
    return false;
  }
}



#if 0
static int
exon_length_5 (List_T pairs) {
  int exon_length = 0;
  List_T p;

  p = pairs;
  while (p != NULL && ((Pair_T) p->first)->gapp == false) {
    exon_length++;
    p = p->rest;
  }

  if (p == NULL) {
    debug13(printf("no intron found, so exon_length_5 = %d\n",END_SPLICESITE_EXON_LENGTH));
    return END_SPLICESITE_EXON_LENGTH;
  } else {
    debug13(printf("intron found, with exon_length_5 = %d\n",exon_length));
    return exon_length;
  }
}
#endif


#if 0
static int
exon_length_3 (List_T path) {
  int exon_length = 0;
  List_T p;

  p = path;
  while (p != NULL && ((Pair_T) p->first)->gapp == false) {
    exon_length++;
    p = p->rest;
  }

  if (p == NULL) {
    debug13(printf("no intron found, so exon_length_3 = %d\n",END_SPLICESITE_EXON_LENGTH));
    return END_SPLICESITE_EXON_LENGTH;
  } else {
    debug13(printf("intron found, with exon_length_3 = %d\n",exon_length));
    return exon_length;
  }
}
#endif


/* Also handles case where novelsplicingp == false */
/* pairs -> pairs */
static List_T
trim_end5_indels (List_T pairs, int ambig_end_length,
		  Dynprog_T dynprog, Univcoord_T chroffset, Univcoord_T chrhigh,
		  char *queryseq_ptr, char *queryuc_ptr,
		  bool watsonp, int genestrand, bool jump_late_p,
		  Pairpool_T pairpool, double defect_rate) {
  List_T path, exon, pairptr, p;
  Pair_T pair;
  int indel_score;
  int nindels;

  int finalscore, continuous_nmatches, continuous_nmismatches, continuous_nopens, continuous_nindels;
  int querydp3_medialgap, genomedp3_medialgap, queryjump, genomejump;
  List_T continuous_gappairs_medialgap;
  int dynprogindex_minor = 0;

  debug3(printf("Starting trim_end5_indels\n"));

  /* Handle first exon */
  if (pairs == NULL) {
    /* *trim5p = false; */
    return (List_T) NULL;
  } else if (ambig_end_length > 0) {
    /* Don't mess with ambiguous end */
    /* *trim5p = false; */
    return pairs;
  } else {
    pair = pairs->first;
    debug3(printf("querystart %d\n",pair->querypos));
  }

  exon = (List_T) NULL;
  while (pairs != NULL && pair->gapp == false && pair->comp != INDEL_COMP && pair->comp != SHORTGAP_COMP) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
  }

  while (pairs != NULL && pair->gapp == false && (((Pair_T) pairs->first)->comp == INDEL_COMP || ((Pair_T) pairs->first)->comp == SHORTGAP_COMP)) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
  }
  debug3(printf("5' End exon:\n"));
  debug3(Pair_dump_list(exon,true));


  if (exon == NULL) {
    /* *trim5p = false; */
    return pairs;

  } else {
    p = exon;
    nindels = 1;
    while (p != NULL && (((Pair_T) p->first)->comp == INDEL_COMP || ((Pair_T) p->first)->comp == SHORTGAP_COMP)) {
      p = List_next(p);
      nindels++;
    }

    indel_score = 0;
    /* Evaluate region distal to indel */
    while (p != NULL) {
      pair = (Pair_T) List_head(p);
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	indel_score += 1;
      } else {
	indel_score -= 3;
      }
      debug3(printf("5' querypos %d => indel_score %d\n",pair->querypos,indel_score));
      p = List_next(p);
    }

#if 0
    for ( i = 0; p != NULL && i < NEARBY_INDEL; p = List_next(p), i++) {
      medial = (Pair_T) p->first;
      if (medial->gapp) {
	debug3(printf("Saw splice medial to 5' end indel\n"));
	nearindelp = true;
      } else if (medial->comp == MATCH_COMP || medial->comp == DYNPROG_MATCH_COMP || medial->comp == AMBIGUOUS_COMP) {
	/* Skip */
      } else {
	debug3(printf("Saw mismatch %c medial to 5' end indel\n",medial->comp));
      }
    }
#endif

    if (pairs == NULL) {
      debug3(printf("No indel/gap\n"));
      path = exon;
      /* *trim5p = false; */

    } else if (exon == NULL) {
      debug3(printf("No 5' exon\n"));
      path = exon;
      /* *trim5p = false; */

#if 0
    } else if (nearindelp == true && max_nmatches < INDEL_SPLICE_ENDLENGTH) {
      debug3(printf("near indel with nmatches %d too low, so trimming it\n",max_nmatches));
      path = (List_T) NULL;
      /* *trim5p = true; */
#endif

    } else if (((Pair_T) pairs->first)->gapp == true) {
      debug3(printf("Peeled all the way to a gap, so not handling with this procedure\n"));
      path = exon;

    } else {
      querydp3_medialgap = ((Pair_T) pairs->first)->querypos - 1;
      genomedp3_medialgap = ((Pair_T) pairs->first)->genomepos - 1;
      queryjump = querydp3_medialgap + 1;
      genomejump = queryjump /*+ extramaterial_end*/;

      continuous_gappairs_medialgap = Dynprog_end5_gap(&dynprogindex_minor,&finalscore,
						       &continuous_nmatches,&continuous_nmismatches,&continuous_nopens,&continuous_nindels,
						       dynprog,&(queryseq_ptr[querydp3_medialgap]),&(queryuc_ptr[querydp3_medialgap]),
						       queryjump,genomejump,querydp3_medialgap,genomedp3_medialgap,
						       chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
						       extraband_end,defect_rate,/*endalign*/QUERYEND_NOGAPS,/*require_pos_score_p*/true);
      debug(printf("CONTINUOUS AT 5 (trim_end5_indels)?\n"));
      debug(printf("CONTINUOUS_GAPPAIRS_MEDIALGAP:\n"));
      debug(Pair_dump_list(continuous_gappairs_medialgap,true));
      debug3(printf("continuous finalscore %d\n",finalscore));

      if (finalscore > 0) {
	debug3(printf("Using continuous\n"));
        path = continuous_gappairs_medialgap;
	/* *trim5p = false; */
	
      } else if (indel_score < 0) {
	debug3(printf("Not enough matches, so trimming it\n"));
	path = (List_T) NULL;
	/* *trim5p = true; */

      } else {
	debug3(printf("Using indel, because indel_score %d > 0\n",indel_score));
	path = exon;		/* exon already has the indel */
	/* *trim5p = false; */
      }
    }

    path = Pairpool_transfer(path,pairs);

    pairs = List_reverse(path);
    pairs = clean_pairs_end5(pairs,ambig_end_length);

    debug3(printf("End of trim_end5_indels: length = %d\n",List_length(pairs)));
    debug3(Pair_dump_list(pairs,true));
    return pairs;
  }
}


/* Also handles case where novelsplicingp == false */
/* pairs -> pairs */
static List_T
trim_end5_exons (bool *indelp, bool *trim5p, int ambig_end_length, List_T pairs,
		 Dynprog_T dynprog, Univcoord_T chroffset, Univcoord_T chrhigh,
		 char *queryseq_ptr, char *queryuc_ptr, int querylength,
		 int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		 Pairpool_T pairpool, double defect_rate) {
  List_T path, exon, pairptr, p;
  Pair_T pair, splice = NULL, gappair;
  int max_nmatches, max_nmismatches, nmatches, nmismatches;
  int max_score, score;
  /* bool nearindelp = false; */
  double medial_prob;
  int peelback, last_indel;

  int finalscore, continuous_nmatches, continuous_nmismatches, continuous_nopens, continuous_nindels;
  int querydp3_medialgap, genomedp3_medialgap, queryjump, genomejump;
  List_T continuous_gappairs_medialgap;
  int dynprogindex_minor = 0;


  debug3(printf("Starting trim_end5_exons with ambig_end_length %d\n",ambig_end_length));
  debug3(Pair_dump_list(pairs,true));

  *indelp = false;

  /* Handle first exon */
  if (pairs == NULL) {
    *trim5p = false;
    debug3(printf("Ending trim_end5_exons because pairs is NULL\n"));
    return (List_T) NULL;

  } else if (0 && ambig_end_length > 0) {
    /* Previously didn't mess with ambiguous end, but this can lead to bad ends */
    *trim5p = false;
    debug3(printf("Ending trim_end5_exons because of ambiguous end\n"));
    return pairs;

  } else {
    pair = pairs->first;
    debug3(printf("querystart %d\n",pair->querypos));
    /* Normally expect pair->querypos to be 0, and want to start with -1 because of the gap */
#if 0
    if (pair->querypos <= ambig_end_length) {
      nmismatches = -1;
    } else {
      nmismatches = (pair->querypos - ambig_end_length) - 1;
    }
#endif
  }

  exon = (List_T) NULL;
  while (pairs != NULL && !pair->gapp /*&& pair->comp != INDEL_COMP && pair->comp != SHORTGAP_COMP */) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
  }


  if (exon == NULL) {
    *trim5p = false;
    return pairs;

  } else if (pair->gapp == false) {
    /* No intron */
    gappair = (Pair_T) NULL;

  } else {
    /* Look for an indel on the medial side */
    gappair = (Pair_T) List_head(exon);

    p = pairs;
    peelback = 0;
    last_indel = 0;
    while (p != NULL && peelback++ < 6) {
      if (((Pair_T) p->first)->comp == INDEL_COMP || ((Pair_T) p->first)->comp == SHORTGAP_COMP || ((Pair_T) p->first)->gapp) {
	last_indel = peelback;
      }
      p = p->rest;
    }
    if (last_indel == 6) {
      while (p != NULL && (((Pair_T) p->first)->comp == INDEL_COMP || ((Pair_T) p->first)->comp == SHORTGAP_COMP || ((Pair_T) p->first)->gapp)) {
	last_indel = peelback++;
	p = p->rest;
      }
    }
    
    if (last_indel > 0) {
      peelback = 0;
      while (pairs != NULL && peelback++ < last_indel) {
	debug3(printf("  Trimming indel on medial side\n"));
	pairptr = pairs;
	pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
	exon = Pairpool_push_existing(exon,pairpool,pair);
#else
	exon = List_push_existing(exon,pairptr);
#endif
      }
    }
  }

  debug3(printf("End exon:\n"));
  debug3(Pair_dump_list(exon,true));


  max_nmatches = max_nmismatches = 0;
  nmatches = nmismatches = 0;
  max_score = score = 0;
  for (p = exon; p != NULL; p = List_next(p)) {
    pair = (Pair_T) List_head(p);
    if (pair->gapp == true) {
      /* Skip the intron gap */
    } else if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
      score += 1;
      nmatches += 1;
    } else {
      score -= 3;
      nmismatches += 1;
    }
    if (score > max_score) {
      max_score = score;
      max_nmatches = nmatches;
      max_nmismatches = nmismatches;
    }
    debug3(printf("5' querypos %d => score %d, max_nmatches %d, max_nmismatches %d\n",
		  pair->querypos,score,max_nmatches,max_nmismatches));
  }


#if 0
  for (p = pairs, i = 0; p != NULL && i < NEARBY_INDEL; p = List_next(p), i++) {
    medial = (Pair_T) p->first;
    if (medial->comp == MATCH_COMP || medial->comp == DYNPROG_MATCH_COMP || medial->comp == AMBIGUOUS_COMP) {
      /* Skip */
    } else if (medial->comp == INDEL_COMP || medial->comp == SHORTGAP_COMP) {
      debug3(printf("Saw indel medial to 5' end intron\n"));
      nearindelp = true;
    } else {
      debug3(printf("Saw mismatch %c medial to 5' end intron\n",medial->comp));
    }
  }
#endif

  debug3(printf("Before end intron, nmatches %d, nmismatches %d\n",max_nmatches,max_nmismatches));
  if (pairs == NULL) {
    debug3(printf("No indel/gap\n"));
    path = exon;
    *trim5p = false;

  } else if (exon == NULL) {
    debug3(printf("No 5' exon\n"));
    path = exon;
    *trim5p = false;

  } else if (gappair == NULL) {
    debug3(printf("No gappair\n"));
    path = exon;
    *trim5p = false;

#if 0
  } else if (exon->rest != NULL && ((Pair_T) exon->rest->first)->disallowedp == true) {
      debug3(printf("Intron is disallowed, so trimming it\n"));
      path = (List_T) NULL;
      *trim5p = true;
#endif

#if 0
  } else if (List_length(exon) - 1 > List_length(pairs)) {
    /* Subtract 1 because gap is included in exon */
    debug3(printf("Exon is more than halfway across %d - 1 > %d, so keeping it\n",List_length(exon),List_length(pairs)));
    path = exon;		/* exon already has the gap */
    *trim5p = false;
#endif

#if 0
  } else if (nearindelp == true && max_nmatches < INDEL_SPLICE_ENDLENGTH) {
    debug3(printf("near indel with nmatches %d too low, so trimming it\n",max_nmatches));
    path = (List_T) NULL;
    *trim5p = true;
#endif

  } else {
    splice = gappair;
    if (nmatches < (int) (END_SUFFICIENT_EXONLENGTH_PCT * (double) querylength) &&
	splice->genomejump > maxintronlen_ends) {
      debug3(printf("End intron is too long, so trimming it\n"));
      path = (List_T) NULL;
      *trim5p = true;

    } else if (nmatches < minendexon) {
      debug3(printf("End exon is too short, so trimming it\n"));
      path = (List_T) NULL;
      *trim5p = true;

    } else if (splice->knowngapp == true && max_nmismatches == 0) {
      debug3(printf("Intron is known and no mismatches, so keeping it\n"));
      path = exon;		/* exon already has the gap */
      *trim5p = false;

#if 0
    } else if (splice->genomejump > maxintronlen) {
      debug3(printf("Intron length %d is too long, so trimming it\n",splice->genomejump));
      path = (List_T) NULL;
      *trim5p = true;
#endif

#if 0
    } else if (enough_matches(nmatches-nmismatches,splice->genomejump/*,splice->donor_prob,splice->acceptor_prob*/) == false) {
      debug3(printf("nmatches %d - nmismatches %d not enough for genomejump %d, so trimming it\n",
		    nmatches,nmismatches,splice->genomejump));
      path = (List_T) NULL;
      *trim5p = true;
#endif

#if 0
    } else if (max_score < 12) {
      /* This eliminates ambig end information */
      debug3(printf("max_score %d < 12, so trimming it\n",max_score));
      path = (List_T) NULL;
      *trim5p = true;
#endif

    } else if (sufficient_splice_prob_strict(/*distal_spliceprob*/cdna_direction >= 0 ? splice->donor_prob : splice->acceptor_prob,
					    /*medial_spliceprob*/cdna_direction >= 0 ? splice->acceptor_prob : splice->donor_prob)) {
      debug3(printf("Keeping first 5' exon with %d matches and %d mismatches\n",max_nmatches,max_nmismatches));
      path = exon;		/* exon already has the gap */
      *trim5p = false;

    } else {
      querydp3_medialgap = ((Pair_T) pairs->first)->querypos - 1;
      genomedp3_medialgap = ((Pair_T) pairs->first)->genomepos - 1;
      queryjump = querydp3_medialgap + 1;
      genomejump = queryjump + extramaterial_end;
      debug3(printf("head of pairs is "));
      debug3(Pair_dump_one((Pair_T) pairs->first,true));
      debug3(printf("\nmedial end starts at query %d, genome %d\n",querydp3_medialgap,genomedp3_medialgap));


      /* Set require_pos_score_p to be false for debugging purposes, but we don't use the result anyway when finalscore is <= 0 */
      continuous_gappairs_medialgap = Dynprog_end5_gap(&dynprogindex_minor,&finalscore,
						       &continuous_nmatches,&continuous_nmismatches,&continuous_nopens,&continuous_nindels,
						       dynprog,&(queryseq_ptr[querydp3_medialgap]),&(queryuc_ptr[querydp3_medialgap]),
						       queryjump,genomejump,querydp3_medialgap,genomedp3_medialgap,
						       chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
						       extraband_end,defect_rate,/*endalign*/QUERYEND_INDELS,/*require_pos_score_p*/false);
      debug(printf("CONTINUOUS AT 5 (trim_end5_exons)?\n"));
      debug(printf("CONTINUOUS_GAPPAIRS_MEDIALGAP:\n"));
      debug(Pair_dump_list(continuous_gappairs_medialgap,true));
      debug3(printf("continuous finalscore %d\n",finalscore));

      if (finalscore > 0) {
	debug3(printf("Using continuous\n"));
        path = continuous_gappairs_medialgap;

	if (continuous_nindels > 0) {
	  *trim5p = true;	/* So calling procedure iterates */
	  *indelp = true;	/* So calling procedure will call trim_end5_indels */
	} else {
	  *trim5p = false;
	}
	
      } else if (sufficient_splice_prob_local(max_nmatches,max_nmismatches,
	                                      /*distal_spliceprob*/cdna_direction >= 0 ? splice->donor_prob : splice->acceptor_prob,
					      /*medial_spliceprob*/cdna_direction >= 0 ? splice->acceptor_prob : splice->donor_prob)) {
        /* Want to keep for comparison of fwd and rev, even if probabilities are poor */
        debug3(printf("Keeping first 5' exon with %d matches and %d mismatches\n",max_nmatches,max_nmismatches));
	path = exon;		/* exon already has the gap */
	*trim5p = false;

      } else {
	/* TODO: Set ambig_end_length_5 here, so default output shows a donor or acceptor end type */
	debug3(printf("Fall through (bad probabilities %f and %f): trimming noncanonical 5' exon\n",splice->donor_prob,splice->acceptor_prob));
	
	medial_prob = (cdna_direction >= 0) ? splice->acceptor_prob : splice->donor_prob;
	if (canonicalp(splice->knowngapp,splice->comp,splice->donor_prob,splice->acceptor_prob,cdna_direction) == true &&
	    medial_prob > 0.95) {
	  *trim5p = false;		/* Not really, since we are trimming, but this stops further work */
	} else {
	  *trim5p = true;
	}
	path = (List_T) NULL;
      }
    }
  }

#ifdef WASTE
  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    path = Pairpool_push_existing(path,pairpool,pair);

  }
#else
  path = Pairpool_transfer(path,pairs);
#endif

  pairs = List_reverse(path);
  pairs = clean_pairs_end5(pairs,ambig_end_length);

  debug3(printf("End of trim_end5_exons: length = %d\n",List_length(pairs)));
  debug3(Pair_dump_list(pairs,true));

#ifdef CHECK_ASSERTIONS
  Pair_check_list_pairs(pairs);
#endif

  return pairs;
}


/* Also handles case where novelsplicingp == false */
/* path -> path */
static List_T
trim_end3_indels (List_T path, int ambig_end_length,
		  Dynprog_T dynprog, Univcoord_T chroffset, Univcoord_T chrhigh,
		  char *queryseq_ptr, char *queryuc_ptr, int querylength,
		  bool watsonp, int genestrand,bool jump_late_p,
		  Pairpool_T pairpool, double defect_rate) {
  List_T pairs, exon, pairptr, p;
  Pair_T pair;
  int indel_score;
  int nindels;

  int finalscore, continuous_nmatches, continuous_nmismatches, continuous_nopens, continuous_nindels;
  int querydp5_medialgap, genomedp5_medialgap, queryjump, genomejump;
  List_T continuous_gappairs_medialgap;
  int dynprogindex_minor = 0;

  debug3(printf("Starting trim_end3_indels\n"));

  /* Handle last exon */
  if (path == NULL) {
    /* *trim3p = false; */
    return (List_T) NULL;
  } else if (ambig_end_length > 0) {
    /* Don't mess with ambiguous end */
    /* *trim3p = false; */
    return path;
  } else {
    pair = path->first;
    debug3(printf("queryend %d\n",pair->querypos));
  }

  exon = (List_T) NULL;
  while (path != NULL && pair->gapp == false && pair->comp != INDEL_COMP && pair->comp != SHORTGAP_COMP) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
  }

  while (path != NULL && pair->gapp == false && (((Pair_T) path->first)->comp == INDEL_COMP || ((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
  }
  debug3(printf("3' End exon:\n"));
  debug3(Pair_dump_list(exon,true));


  if (exon == NULL) {
    /* *trim3p = false; */
    return path;

  } else {
    p = exon;
    nindels = 1;
    while (p != NULL && (((Pair_T) p->first)->comp == INDEL_COMP || ((Pair_T) p->first)->comp == SHORTGAP_COMP)) {
      p = List_next(p);
      nindels++;
    }

    indel_score = 0;
    /* Evaluate region distal to indel */
    while (p != NULL) {
      pair = (Pair_T) List_head(p);
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	indel_score += 1;
      } else {
	indel_score -= 3;
      }
      debug3(printf("3' querypos %d => indel_score %d\n",pair->querypos,indel_score));
      p = List_next(p);
    }

#if 0
    for ( i = 0; p != NULL && i < NEARBY_INDEL; p = List_next(p), i++) {
      medial = (Pair_T) p->first;
      if (medial->gapp) {
	debug3(printf("Saw splice medial to 3' end indeln"));
	nearindelp = true;
      } else if (medial->comp == MATCH_COMP || medial->comp == DYNPROG_MATCH_COMP || medial->comp == AMBIGUOUS_COMP) {
	/* Skip */
      } else {
	debug3(printf("Saw mismatch medial %c to 3' end indel\n",medial->comp));
      }
    }
#endif

    if (path == NULL) {
      debug3(printf("No indel/gap\n"));
      pairs = exon;
      /* *trim3p = false; */

    } else if (exon == NULL) {
      debug3(printf("No 3' exon\n"));
      pairs = exon;
      /* *trim3p = false; */

#if 0
    } else if (nearindelp == true && max_nmatches < INDEL_SPLICE_ENDLENGTH) {
      debug3(printf("near indel with nmatches %d too low, so trimming it\n",max_nmatches));
      pairs = (List_T) NULL;
      /* *trim3p = true; */
#endif
    
    } else if (((Pair_T) path->first)->gapp == true) {
      debug3(printf("Peeled all the way to a gap, so not handling with this procedure\n"));
      pairs = exon;

    } else {
      querydp5_medialgap = ((Pair_T) path->first)->querypos + 1;
      genomedp5_medialgap = ((Pair_T) path->first)->genomepos + 1;
      queryjump = querylength - querydp5_medialgap;
      genomejump = queryjump /*+ extramaterial_end*/;

      continuous_gappairs_medialgap = Dynprog_end3_gap(&dynprogindex_minor,&finalscore,
						       &continuous_nmatches,&continuous_nmismatches,&continuous_nopens,&continuous_nindels,
						       dynprog,&(queryseq_ptr[querydp5_medialgap]),&(queryuc_ptr[querydp5_medialgap]),
						       /*rlength*/queryjump,/*glength*/genomejump,
						       /*roffset*/querydp5_medialgap,/*goffset*/genomedp5_medialgap,
						       chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
						       extraband_end,defect_rate,/*endalign*/QUERYEND_NOGAPS,/*require_pos_score_p*/true);
      debug(printf("CONTINUOUS AT 3 (trim_end3_indels)?\n"));
      debug(printf("CONTINUOUS_GAPPAIRS_MEDIALGAP:\n"));
      debug(Pair_dump_list(continuous_gappairs_medialgap,true));
      debug3(printf("continuous finalscore %d\n",finalscore));

      if (finalscore > 0) {
	debug3(printf("Using continuous\n"));
	pairs = List_reverse(continuous_gappairs_medialgap);
	/* *trim3p = false; */
	
      } else if (indel_score < 0) {
	debug3(printf("Not enough matches, so trimming it\n"));
	pairs = (List_T) NULL;
	/* *trim3p = true; */

      } else {
	debug3(printf("Using indel, because indel_score %d > 0\n",indel_score));
	pairs = exon;
	/* *trim3p = false; */
      }
    }

    pairs = Pairpool_transfer(pairs,path);

    path = List_reverse(pairs);
    path = clean_path_end3(path,ambig_end_length);

    debug3(printf("End of trim_end3_indels: length = %d\n",List_length(path)));
    debug3(Pair_dump_list(path,true));
    return path;
  }
}


/* Also handles case where novelsplicingp == false */
/* path -> path */
static List_T
trim_end3_exons (bool *indelp, bool *trim3p, int ambig_end_length, List_T path,
		 Dynprog_T dynprog, Univcoord_T chroffset, Univcoord_T chrhigh,
		 char *queryseq_ptr, char *queryuc_ptr, int querylength,
		 int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		 Pairpool_T pairpool, double defect_rate) {
  List_T pairs, exon, pairptr, p;
  Pair_T pair, splice = NULL, gappair;
  int max_nmatches, max_nmismatches, nmatches, nmismatches;
  int max_score, score;
  /* bool nearindelp = false; */
  double medial_prob;
  int peelback, last_indel;

  int finalscore, continuous_nmatches, continuous_nmismatches, continuous_nopens, continuous_nindels;
  int querydp5_medialgap, genomedp5_medialgap, queryjump, genomejump;
  List_T continuous_gappairs_medialgap;
  int dynprogindex_minor = 0;

  debug3(printf("Starting trim_end3_exons with ambig_end_length %d\n",ambig_end_length));

  *indelp = false;

  /* Handle last exon */
  if (path == NULL) {
    *trim3p = false;
    debug3(printf("Ending trim_end3_exons because path is NULL\n"));
    return (List_T) NULL;

  } else if (0 && ambig_end_length > 0) {
    /* Previously didn't mess with ambiguous end, but this can lead to bad ends */
    *trim3p = false;
    debug3(printf("Ending trim_end3_exons because of ambiguous end\n"));
    return path;

  } else {
    pair = path->first;
    debug3(printf("queryend %d\n",pair->querypos));
#if 0
    /* Normally expect pair->querypos to be 0, and want to start with -1 because of the gap */
    if (pair->querypos >= (querylength - 1) - ambig_end_length) {
      nmismatches = -1;
    } else {
      nmismatches = (querylength - 1) - ambig_end_length - pair->querypos - 1;
    }
#endif
  }

  exon = (List_T) NULL;
  while (path != NULL && !pair->gapp) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
#ifdef WASTE
    exon = Pairpool_push_existing(exon,pairpool,pair);
#else
    exon = List_push_existing(exon,pairptr);
#endif
  }


  if (exon == NULL) {
    *trim3p = false;
    return path;

  } else if (pair->gapp == false) {
    /* No intron */
    gappair = (Pair_T) NULL;

  } else {
    /* Look for an indel on the medial side */
    gappair = (Pair_T) List_head(exon);

    p = path;
    peelback = 0;
    last_indel = 0;
    while (p != NULL && peelback++ < 6) {
      if (((Pair_T) p->first)->comp == INDEL_COMP || ((Pair_T) p->first)->comp == SHORTGAP_COMP || ((Pair_T) p->first)->gapp) {
	last_indel = peelback;
      }
      p = p->rest;
    }
    if (last_indel == 6) {
      while (p != NULL && (((Pair_T) p->first)->comp == INDEL_COMP || ((Pair_T) p->first)->comp == SHORTGAP_COMP || ((Pair_T) p->first)->gapp)) {
	last_indel = peelback++;
	p = p->rest;
      }
    }

    if (last_indel > 0) {
      peelback = 0;
      while (path != NULL && peelback++ < last_indel) {
	debug3(printf("  Trimming indel on medial side\n"));
	pairptr = path;
	path = Pairpool_pop(path,&pair);
#ifdef WASTE
	exon = Pairpool_push_existing(exon,pairpool,pair);
#else
	exon = List_push_existing(exon,pairptr);
#endif
      }
    }
  }

  debug3(printf("End exon:\n"));
  debug3(Pair_dump_list(exon,true));


  max_nmatches = max_nmismatches = 0;
  nmatches = nmismatches = 0;
  max_score = score = 0;
  for (p = exon; p != NULL; p = List_next(p)) {
    pair = (Pair_T) List_head(p);
    if (pair->gapp == true) {
      /* Skip the intron gap */
    } else if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
      score += 1;
      nmatches += 1;
    } else {
      score -= 3;
      nmismatches += 1;
    }
    if (score > max_score) {
      max_score = score;
      max_nmatches = nmatches;
      max_nmismatches = nmismatches;
    }
    debug3(printf("3' querypos %d => score %d, max_nmatches %d, max_nmismatches %d\n",
	          pair->querypos,score,max_nmatches,max_nmismatches));
  }


#if 0
  for (p = path, i = 0; p != NULL && i < NEARBY_INDEL; p = List_next(p), i++) {
    medial = (Pair_T) p->first;
    if (medial->comp == MATCH_COMP || medial->comp == DYNPROG_MATCH_COMP || medial->comp == AMBIGUOUS_COMP) {
      /* Skip */
    } else if (medial->comp == INDEL_COMP || medial->comp == SHORTGAP_COMP) {
      debug3(printf("Saw indel medial to 3' end intron\n"));
      nearindelp = true;
    } else {
      debug3(printf("Saw mismatch medial %c to 3' end intron\n",medial->comp));
    }
  }
#endif

  debug3(printf("Before end intron, nmatches %d, nmismatches %d\n",max_nmatches,max_nmismatches));
  if (path == NULL) {
    debug3(printf("No gap\n"));
    pairs = exon;
    *trim3p = false;

  } else if (exon == NULL) {
    debug3(printf("No 3' exon\n"));
    pairs = exon;
    *trim3p = false;

  } else if (gappair == NULL) {
    debug3(printf("No gappair\n"));
    pairs = exon;
    *trim3p = false;

#if 0
  } else if (exon->rest != NULL && ((Pair_T) exon->rest->first)->disallowedp == true) {
    debug3(printf("Intron is disallowed, so trimming it\n"));
    pairs = (List_T) NULL;
    *trim3p = true;
#endif

#if 0
  } else if (List_length(exon) - 1 > List_length(path)) {
    /* Subtract 1 because gap is included in exon */
    debug3(printf("Exon is more than halfway across %d - 1 > %d, so keeping it\n",List_length(exon),List_length(path)));
    pairs = exon;		/* exon already has the gap */
    *trim3p = false;
#endif

#if 0
  } else if (nearindelp == true && max_nmatches < INDEL_SPLICE_ENDLENGTH) {
    debug3(printf("near indel with nmatches %d too low, so trimming it\n",max_nmatches));
    pairs = (List_T) NULL;
    *trim3p = true;
#endif
    
  } else {
    splice = gappair;
    if (nmatches < (int) (END_SUFFICIENT_EXONLENGTH_PCT * (double) querylength) &&
	splice->genomejump > maxintronlen_ends) {
      debug3(printf("End intron is too long, so trimming it\n"));
      pairs = (List_T) NULL;
      *trim3p = true;

    } else if (nmatches < minendexon) {
      debug3(printf("End exon is too short, so trimming it\n"));
      pairs = (List_T) NULL;
      *trim3p = true;

    } else if (splice->knowngapp == true && max_nmismatches == 0) {
      debug3(printf("Intron is known and no mismatches, so keeping it\n"));
      pairs = exon;		/* exon already has the gap */
      *trim3p = false;

#if 0
    } else if (splice->genomejump > maxintronlen) {
      debug3(printf("Intron length %d is too long, so trimming it\n",pair->genomejump));
      pairs = (List_T) NULL;
      *trim3p = true;
#endif

#if 0
    } else if (enough_matches(nmatches-nmismatches,splice->genomejump/*,splice->donor_prob,splice->acceptor_prob*/) == false) {
      debug3(printf("nmatches %d - nmismatches %d not enough for genomejump %d, so trimming it\n",
		    nmatches,nmismatches,splice->genomejump));
      pairs = (List_T) NULL;
      *trim3p = true;
#endif

#if 0
    } else if (max_score < 12) {
      /* This eliminates ambig end information */
      debug3(printf("max_score %d < 12, so trimming it\n",max_score));
      pairs = (List_T) NULL;
      *trim3p = true;
#endif

    } else if (sufficient_splice_prob_strict(/*distal_spliceprob*/cdna_direction >= 0 ? splice->acceptor_prob : splice->donor_prob,
					     /*medial_spliceprob*/cdna_direction >= 0 ? splice->donor_prob : splice->acceptor_prob)) {
      debug3(printf("Keeping last 3' exon with %d matches and %d mismatches\n",max_nmatches,max_nmismatches));
      pairs = exon;		/* exon already has the gap */
      *trim3p = false;
	
    } else {
      querydp5_medialgap = ((Pair_T) path->first)->querypos + 1;
      genomedp5_medialgap = ((Pair_T) path->first)->genomepos + 1;
      queryjump = querylength - querydp5_medialgap;
      genomejump = queryjump + extramaterial_end;
      debug3(printf("head of path is "));
      debug3(Pair_dump_one((Pair_T) path->first,true));
      debug3(printf("\nmedial end starts at query %d, genome %d\n",querydp5_medialgap,genomedp5_medialgap));

      /* Set require_pos_score_p to be false for debugging purposes, but we don't use the result anyway when finalscore is <= 0 */
      continuous_gappairs_medialgap = Dynprog_end3_gap(&dynprogindex_minor,&finalscore,
						       &continuous_nmatches,&continuous_nmismatches,&continuous_nopens,&continuous_nindels,
						       dynprog,&(queryseq_ptr[querydp5_medialgap]),&(queryuc_ptr[querydp5_medialgap]),
						       queryjump,genomejump,querydp5_medialgap,genomedp5_medialgap,
						       chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
						       extraband_end,defect_rate,/*endalign*/QUERYEND_INDELS,/*require_pos_score_p*/false);
      debug(printf("CONTINUOUS AT 3 (trim_end3_exons)?\n"));
      debug(printf("CONTINUOUS_GAPPAIRS_MEDIALGAP:\n"));
      debug(Pair_dump_list(continuous_gappairs_medialgap,true));
      debug3(printf("continuous finalscore %d\n",finalscore));

      if (finalscore > 0) {
	debug3(printf("Using continuous\n"));
	pairs = List_reverse(continuous_gappairs_medialgap);
	if (continuous_nindels > 0) {
	  *trim3p = true;	/* So calling procedure iterates */
	  *indelp = true; /* So calling procedure will call trim_end3_indels */
	} else {
	  *trim3p = false;
	}
	
      } else if (sufficient_splice_prob_local(max_nmatches,max_nmismatches,
					      /*distal_spliceprob*/cdna_direction >= 0 ? splice->acceptor_prob : splice->donor_prob,
					      /*medial_spliceprob*/cdna_direction >= 0 ? splice->donor_prob : splice->acceptor_prob)) {
	/* Want to keep for comparison of fwd and rev, even if probabilities are poor */
	debug3(printf("Keeping last 3' exon with %d matches and %d mismatches\n",max_nmatches,max_nmismatches));
	pairs = exon;		/* exon already has the gap */
	*trim3p = false;

      } else {
	/* TODO: Set ambig_end_length_3 here, so default output shows a donor or acceptor end type */
	debug3(printf("Fall through (bad probabilities %f and %f): trimming noncanonical 3' exon\n",splice->donor_prob,splice->acceptor_prob));
	
	medial_prob = (cdna_direction >= 0) ? splice->donor_prob : splice->acceptor_prob;
	if (canonicalp(splice->knowngapp,splice->comp,splice->donor_prob,splice->acceptor_prob,cdna_direction) == true &&
	    medial_prob > 0.95) {
	  *trim3p = false;		/* Not really, since we are trimming, but this stops further work */
	} else {
	  *trim3p = true;
	}
	pairs = (List_T) NULL;
      }
    }
  }

#ifdef WASTE
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
  }
#else
  pairs = Pairpool_transfer(pairs,path);
#endif

  path = List_reverse(pairs);
  path = clean_path_end3(path,ambig_end_length);

  debug3(printf("End of trim_end3_exons: length = %d\n",List_length(path)));
  debug3(Pair_dump_list(path,true));

#ifdef CHECK_ASSERTIONS
  Pair_check_list_path(path);
#endif

  return path;
}



/* This procedure fills in introns and replaces non-canonical introns
   with deletions, so it should be called after all dynamic
   programming procedures */
static List_T 
fill_in_gaps (List_T path, Pairpool_T pairpool, char *queryseq_ptr, 
	      Univcoord_T chroffset, Univcoord_T chrhigh,
	      int cdna_direction, bool watsonp, int ngap) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;

  int leftquerypos, leftgenomepos, rightquerypos, rightgenomepos,
    introntype, intronlength, genomicpos, k;
  char comp;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt, c1, c2, c2_alt;
  bool intronp, introntype_found_p;


  if (path == NULL) {
    return (List_T) NULL;
  } else {
    pair = path->first;
  }

  while (path != NULL && ((Pair_T) path->first)->gapp == true) {
    /* Gap at beginning of alignment.  Can occur after smoothing. */
    debug7(printf("Gap %p at beginning of alignment\n",pair));
    path = Pairpool_pop(path,&pair);
  }

  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;

#ifdef PMAP
    if (pair->cdna == BACKTRANSLATE_CHAR) {
      pair->cdna = 'N';
    }
#endif
    if (pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->gapp == false) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (path->rest == NULL) {
      /* Gap at end of alignment.  Can occur after smoothing. */
      debug7(printf("Gap at end of alignment\n"));
      path = Pairpool_pop(path,&pair);

    } else if (pairs == NULL) {
      /* Gap at beginning of alignment.  Skip. */
      debug7(printf("Gap at beginning of alignment\n"));
      path = Pairpool_pop(path,&pair);

    } else {
      /* pairptr = path; -- save */
      path = Pairpool_pop(path,&pair);

      /* Discard gap; do not push */
      leftpair = path->first;
      rightpair = pairs->first;

      if (pair->comp == DUALBREAK_COMP) {
	pairs = add_dualbreak(pairs,queryseq_ptr,
			      chroffset,chrhigh,cdna_direction,watsonp,
			      leftpair,rightpair,pairpool,ngap);
      } else {

	leftquerypos = leftpair->querypos;
	leftgenomepos = leftpair->genomepos;
	/* if (leftpair->cdna == ' ') leftquerypos--; -- For old dynamic programming */
	/* if (leftpair->genome == ' ') leftgenomepos--; -- For old dynamic programming */
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;
	intronlength = rightgenomepos - leftgenomepos - 1;

	introntype_found_p = false;
	if (pair->knowngapp == true) {
	  debug7(printf("known intron is true, so an intron\n"));
	  intronp = true;
	} else if (splicingp == false) {
	  debug7(printf("splicingp is false, so not an intron, but an indel\n"));
	  intronp = false;
#if 0
	} else if (sensedir == SENSE_NULL) {
	  /* Can lead to very large deletions */
	  debug7(printf("sensedir == SENSE_NULL, so not an intron, but an indel\n"));
	  intronp = false;
#endif
	} else if (intronlength < min_intronlength) {
	  debug7(printf("intronlength %d < min_intronlength %d, so not an intron, but an indel\n",
			intronlength,min_intronlength));
	  intronp = false;
	} else if (intronlength >= max_deletionlength) {
	  debug7(printf("intronlength %d >= max_deletionlength %d, so an intron, not an indel\n",
			intronlength,max_deletionlength));
	  intronp = true;
	} else {
	  /* Possible short intron */
	  left1 = get_genomic_nt(&left1_alt,leftgenomepos+1,chroffset,chrhigh,watsonp);
	  left2 = get_genomic_nt(&left2_alt,leftgenomepos+2,chroffset,chrhigh,watsonp);
	  right2 = get_genomic_nt(&right2_alt,rightgenomepos-2,chroffset,chrhigh,watsonp);
	  right1 = get_genomic_nt(&right1_alt,rightgenomepos-1,chroffset,chrhigh,watsonp);
	  debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
	  introntype = Intron_type(left1,left2,right2,right1,
				   left1_alt,left2_alt,right2_alt,right1_alt,
				   cdna_direction);
	  introntype_found_p = true;
	  debug7(printf("  Introntype at %u..%u is %s (cdna_direction %d)\n",
			leftgenomepos,rightgenomepos,Intron_type_string(introntype),cdna_direction));

	  if ((cdna_direction >= 0 && introntype == GTAG_FWD)
#ifndef PMAP
	      || (cdna_direction <= 0 && introntype == GTAG_REV)
#endif
	      ) {
	    intronp = true;
	  } else {
	    intronp = false;
	  }
	}

	if (intronp == false) {
	  debug7(printf("  Gap is not an intron (intronlength %d).  Adding pairs from %d downto %d\n",
			intronlength,rightgenomepos-1,leftgenomepos+1));
	  for (k = rightquerypos - 1; k > leftquerypos; --k) {
#if 0				/* PMAP */
	    c1 = Sequence_codon_char(queryaaseq_ptr[k/3],k%3);
#else
	    c1 = queryseq_ptr[k];
#endif
	    pairs = Pairpool_push(pairs,pairpool,k,rightgenomepos,c1,INDEL_COMP,' ',' ',/*dynprogindex*/0);
	  }

	  for (genomicpos = rightgenomepos - 1; genomicpos > leftgenomepos; --genomicpos) {
	    c2 = get_genomic_nt(&c2_alt,genomicpos,chroffset,chrhigh,watsonp);
	    pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',/*comp*/SHORTGAP_COMP,c2,c2_alt,
				  /*dynprogindex*/0);
	  }
	} else {
	  if (introntype_found_p == false) {
	    left1 = get_genomic_nt(&left1_alt,leftgenomepos+1,chroffset,chrhigh,watsonp);
	    left2 = get_genomic_nt(&left2_alt,leftgenomepos+2,chroffset,chrhigh,watsonp);
	    right2 = get_genomic_nt(&right2_alt,rightgenomepos-2,chroffset,chrhigh,watsonp);
	    right1 = get_genomic_nt(&right1_alt,rightgenomepos-1,chroffset,chrhigh,watsonp);
	    debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
	    introntype = Intron_type(left1,left2,right2,right1,
				     left1_alt,left2_alt,right2_alt,right1_alt,
				     cdna_direction);
	    debug7(printf("  Introntype at %u..%u is %s (cdna_direction %d)\n",
			  leftgenomepos,rightgenomepos,Intron_type_string(introntype),cdna_direction));
	  }

	  if (cdna_direction == 0) {
	    /* cdna_direction of 0 should happen only from Stage3_merge_local_splice */
	    comp = NONINTRON_COMP;
	  } else {
	    comp = pair->comp;
	  }

	  debug7(printf("Adding an intron at %d..%d, currently of type %c, with introntype %d\n",
			leftpair->querypos,rightpair->querypos,comp,introntype));
	  pairs = add_intron(pairs,chroffset,chrhigh,leftpair,rightpair,comp,introntype,ngap,
			     watsonp,pairpool);
	}
      }
    }
  }

  debug7(printf("Final length: %d\n",List_length(pairs)));

  return pairs;
}

static List_T 
add_queryseq_offset (List_T path, int queryseq_offset
#ifdef WASTE
		     , Pairpool_T pairpool
#endif
		     ) {
  List_T pairs = NULL;
  Pair_T pair;

  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */

    /* Previously excluded cases where pair->gapp was true, but this failed on chimeric paths */
    /* Now excluding again, because we are running this before chimera detection and merging */
    pair = (Pair_T) path->first;
    if (pair->gapp == false) {
      pair->querypos += queryseq_offset;
    }
#ifdef WASTE
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
    pairs = List_transfer_one(pairs,&path);
#endif
  }

  return pairs;
}


static void
add_skiplength (List_T pairs, int skiplength) {
  List_T p;
  Pair_T pair;

  for (p = pairs; p != NULL; p = p->rest) {
    pair = (Pair_T) p->first;
    if (pair->gapp == true) {
      /* Skip */
    } else if (pair->querypos >= HALFLEN) {
      pair->querypos += skiplength;
    }
  }
  return;
}


static void
Stage3_free_pairarray (T *old) {
#if 0
  if ((*old)->pairarray_freeable_p == true) {
    FREE_OUT((*old)->pairarray);
    (*old)->pairarray_freeable_p = false;
  }
#else
  FREE_OUT((*old)->pairarray);
#endif
  return;
}


/* Does not alter pairs, except for adding subseq_offset to querypos,
   in case we need to re-compute alignment for chimera */
static struct Pair_T *
make_pairarray (int *npairs, List_T *pairs, int cdna_direction, bool watsonp,
		Pairpool_T pairpool, char *queryseq_ptr,
		Univcoord_T chroffset, Univcoord_T chrhigh,
		int ngap, int subseq_offset, int skiplength) {
  struct Pair_T *pairarray;
  List_T printpairs, printpath, path, p;
  Pair_T oldpair, newpair;

  
#if 0
  /* Better to use clip_pairs_end5_chromosomal_bound and
     clip_path_end3_chromosomal_bound instead */
  printpairs = Pairpool_copy_bounded(*pairs,pairpool,/*chrlength*/chrhigh - chroffset);
#else
  printpairs = Pairpool_copy(*pairs,pairpool);
#endif

  printpath = List_reverse(printpairs);
  printpairs = fill_in_gaps(printpath,pairpool,queryseq_ptr,
			    chroffset,chrhigh,cdna_direction,watsonp,ngap);

  if (subseq_offset != 0) {
    path = List_reverse(*pairs);
#ifdef WASTE
    *pairs = add_queryseq_offset(path,subseq_offset,pairpool);
#else
    *pairs = add_queryseq_offset(path,subseq_offset);
#endif
    
    printpath = List_reverse(printpairs);
#ifdef WASTE
    printpairs = add_queryseq_offset(printpath,subseq_offset,pairpool);
#else
    printpairs = add_queryseq_offset(printpath,subseq_offset);
#endif
  }

  if (skiplength != 0) {
    add_skiplength(*pairs,skiplength);
    add_skiplength(printpairs,skiplength);
  }

  if ((*npairs = List_length(printpairs)) == 0) {
    return (struct Pair_T *) NULL;
  } else {
    /* Used to be Pair_block_copy */
    newpair = pairarray = (struct Pair_T *) MALLOC_OUT(*npairs*sizeof(struct Pair_T));
    for (p = printpairs; p != NULL; p = p->rest) {
      oldpair = (Pair_T) p->first;
      memcpy(newpair++,oldpair,sizeof(struct Pair_T));
    }

    /* No need to free newpairs, since they belong to pairpool */
    Pair_set_genomepos(pairarray,*npairs,chroffset,chrhigh,watsonp);

    return pairarray;
  }
}


/* Does not alter pairs, except for adding subseq_offset to querypos,
   in case we need to re-compute alignment for chimera */
static bool
make_pairarray_merge (T this_left, int cdna_direction, bool watsonp,
		      Pairpool_T pairpool, char *queryseq_ptr,
		      Univcoord_T chroffset, Univcoord_T chrhigh,
		      int ngap, int subseq_offset, int skiplength, bool new_gap_p) {
  struct Pair_T *pairarray, *pairarray_save;
  List_T printpairs, printpath, path, p;
  Pair_T oldpair, newpair;
  int ncanonical, nsemicanonical;
  double min_splice_prob;

  pairarray_save = this_left->pairarray;

  if (new_gap_p == true) {
    path = List_reverse(this_left->pairs);
    this_left->pairs = assign_gap_types(path,cdna_direction,this_left->watsonp,queryseq_ptr,
					this_left->chrnum,this_left->chroffset,this_left->chrhigh,
					pairpool);
  }

  debug10(Pair_dump_list(this_left->pairs,true));
      
  this_left->cdna_direction = cdna_direction;
  
  printpairs = Pairpool_copy(this_left->pairs,pairpool);

  printpath = List_reverse(printpairs);
  printpairs = fill_in_gaps(printpath,pairpool,queryseq_ptr,
			    chroffset,chrhigh,cdna_direction,watsonp,ngap);

  if (List_length(printpairs) == 0) {
    this_left->pairarray = pairarray_save;
    /* this_left->pairarray_freeable_p = false; */
    return false;
    
  } else {
    if (subseq_offset != 0) {
      path = List_reverse(this_left->pairs);
#ifdef WASTE
      this_left->pairs = add_queryseq_offset(path,subseq_offset,pairpool);
#else
      this_left->pairs = add_queryseq_offset(path,subseq_offset);
#endif
    
      printpath = List_reverse(printpairs);
#ifdef WASTE
      printpairs = add_queryseq_offset(printpath,subseq_offset,pairpool);
#else
      printpairs = add_queryseq_offset(printpath,subseq_offset);
#endif
    }

    if (skiplength != 0) {
      add_skiplength(this_left->pairs,skiplength);
      add_skiplength(printpairs,skiplength);
    }

    Stage3_free_pairarray(&this_left);
    this_left->npairs = List_length(printpairs);

    /* Used to be Pair_block_copy */
    newpair = pairarray = (struct Pair_T *) MALLOC_OUT(this_left->npairs*sizeof(struct Pair_T));
    for (p = printpairs; p != NULL; p = p->rest) {
      oldpair = (Pair_T) p->first;
      memcpy(newpair++,oldpair,sizeof(struct Pair_T));
    }

    /* No need to free newpairs, since they belong to pairpool */
    Pair_set_genomepos(pairarray,this_left->npairs,chroffset,chrhigh,watsonp);
    this_left->pairarray = pairarray;
    this_left->pairarray_freeable_p = true;

    this_left->goodness =
      Pair_fracidentity_array(&this_left->matches,&this_left->unknowns,&this_left->mismatches,
			      &this_left->qopens,&this_left->qindels,&this_left->topens,&this_left->tindels,
			      &ncanonical,&nsemicanonical,&this_left->noncanonical,
			      &min_splice_prob,this_left->pairarray,this_left->npairs,this_left->cdna_direction);

    return true;
  }

}


static void
make_pairarrays_chimera (T this_left, T this_right,
			 char *queryseq_ptr, Pairpool_T pairpool, int gaplength, int ngap) {
  List_T printpairs_left, printpath_left, printpairs_right, printpath_right, p;
  Pair_T oldpair, newpair;
  int newnpairs;
  int ncanonical, nsemicanonical;
  double min_splice_prob;

  
  /* Revise statistics */
  Pair_fracidentity(&this_left->matches,&this_left->unknowns,&this_left->mismatches,
		    &this_left->qopens,&this_left->qindels,&this_left->topens,&this_left->tindels,
		    &ncanonical,&nsemicanonical,&this_left->noncanonical,
		    &min_splice_prob,this_left->pairs,this_left->cdna_direction);

  Pair_fracidentity(&this_right->matches,&this_right->unknowns,&this_right->mismatches,
		    &this_right->qopens,&this_right->qindels,&this_right->topens,&this_right->tindels,
		    &ncanonical,&nsemicanonical,&this_right->noncanonical,
		    &min_splice_prob,this_right->pairs,this_right->cdna_direction);


  printpairs_left = Pairpool_copy(this_left->pairs,pairpool);
  printpath_left = List_reverse(printpairs_left);
  printpairs_left = fill_in_gaps(printpath_left,pairpool,queryseq_ptr,
				 this_left->chroffset,this_left->chrhigh,
				 this_left->cdna_direction,this_left->watsonp,ngap);

  printpairs_right = Pairpool_copy(this_right->pairs,pairpool);
  printpath_right = List_reverse(printpairs_right);
  printpairs_right = fill_in_gaps(printpath_right,pairpool,queryseq_ptr,
				  this_right->chroffset,this_right->chrhigh,
				  this_right->cdna_direction,this_right->watsonp,ngap);


  /* Do not use subseq_offset or skiplength for chimeras, since we are
     working on the original queryseq, not querysubseq */

  this_left->npairs = List_length(printpairs_left);
  this_right->npairs = List_length(printpairs_right);

  Stage3_free_pairarray(&this_left);
  Stage3_free_pairarray(&this_right);
  if ((newnpairs = this_left->npairs + this_right->npairs) == 0) {
    this_left->pairarray = (struct Pair_T *) NULL;
    this_right->pairarray = (struct Pair_T *) NULL;
    this_left->pairarray_freeable_p = false;
    this_right->pairarray_freeable_p = false;

  } else {
    /* Need to have a single pairarray for this_left, so we can translate protein correctly */
    newpair = this_left->pairarray = (struct Pair_T *) MALLOC_OUT((newnpairs + gaplength)*sizeof(struct Pair_T));
    this_left->pairarray_freeable_p = true;

    for (p = printpairs_left; p != NULL; p = p->rest) {
      oldpair = (Pair_T) p->first;
      memcpy(newpair++,oldpair,sizeof(struct Pair_T));
    }
    Pair_set_genomepos(this_left->pairarray,this_left->npairs,this_left->chroffset,this_left->chrhigh,
		       this_left->watsonp);
    this_left->goodness =
      Pair_fracidentity_array(&this_left->matches,&this_left->unknowns,&this_left->mismatches,
			      &this_left->qopens,&this_left->qindels,&this_left->topens,&this_left->tindels,
			      &ncanonical,&nsemicanonical,&this_left->noncanonical,
			      &min_splice_prob,this_left->pairarray,this_left->npairs,this_left->cdna_direction);


    newpair = this_right->pairarray = &(this_left->pairarray[this_left->npairs + gaplength]);
    this_right->pairarray_freeable_p = false; /* This is the only case of setting pairarray_freeable_p to be false */

    for (p = printpairs_right; p != NULL; p = p->rest) {
      oldpair = (Pair_T) p->first;
      memcpy(newpair++,oldpair,sizeof(struct Pair_T));
    }
    Pair_set_genomepos(this_right->pairarray,this_right->npairs,this_right->chroffset,this_right->chrhigh,
		       this_right->watsonp);
    this_right->goodness =
      Pair_fracidentity_array(&this_right->matches,&this_right->unknowns,&this_right->mismatches,
			      &this_right->qopens,&this_right->qindels,&this_right->topens,&this_right->tindels,
			      &ncanonical,&nsemicanonical,&this_right->noncanonical,
			      &min_splice_prob,this_right->pairarray,this_right->npairs,this_right->cdna_direction);
  }

  return;
}


void
Stage3_count_paths (int *npaths_primary, int *npaths_altloc, List_T stage3list) {
  T this;

  *npaths_primary = *npaths_altloc = 0;

  while (stage3list != NULL) {
    this = (T) List_head(stage3list);
    if (altlocp[this->chrnum] == true) {
      (*npaths_altloc) += 1;
    } else {
      (*npaths_primary) += 1;
    }
    stage3list = List_next(stage3list);
  }

  return;
}


#define MAPQ_MAXIMUM_SCORE 40

void
Stage3_compute_mapq (List_T stage3list) {
  T this;
  List_T p;
  int best_absmq_score;
  float total = 0.0, q;

  if (stage3list != NULL) {
    /* Use the first entry to initialize best_absmq_score */
    p = stage3list;
    this = (T) List_head(p);
    best_absmq_score = this->absmq_score = this->matches - 10*this->mismatches;
    p = List_next(p);

    while (p != NULL) {
      this = (T) List_head(p);
      if ((this->absmq_score = this->matches - 10*this->mismatches) > best_absmq_score) {
	best_absmq_score = this->absmq_score;
      }
      p = List_next(p);
    }
  }

  for (p = stage3list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    this->absmq_score -= best_absmq_score;
    total += fasterexp(this->absmq_score);
  }

  for (p = stage3list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);

    if ((q = 1.0 - fasterexp(this->absmq_score) / total) < 1.0e-4 /* 10^-4.0 */) {
      this->mapq_score = 40;
    } else {
      this->mapq_score = rint(-10.0 * log10(q));
    }

    this->absmq_score += MAPQ_MAXIMUM_SCORE;
    if (this->absmq_score < 0) {
      this->absmq_score = 0;
    }

  }
}


void
Stage3_recompute_coverage (List_T stage3list, Sequence_T queryseq) {
  List_T p;
  T stage3;
  Pair_T start, end;
  int querypos1, querypos2;
  int trim_start, trim_end, skiplength;

  trim_start = Sequence_trim_start(queryseq);
  trim_end = Sequence_trim_end(queryseq);
  skiplength = Sequence_skiplength(queryseq);

  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (T) List_head(p);
    if (stage3->npairs == 0) {
      stage3->trimmed_coverage = 0.0;

    } else {
      start = &(stage3->pairarray[0]);
      end = &(stage3->pairarray[stage3->npairs - 1]);

      querypos1 = start->querypos;
      querypos2 = end->querypos;

#if 0
      if (querypos2 + 1 > trim_end) {
	effective_trim_end = querypos2 + 1;
      } else {
	effective_trim_end = trim_end;
      }
      if (querypos1 < trim_start) {
	effective_trim_start = querypos1;
      } else {
	effective_trim_start = trim_start;
      }
#endif

      stage3->trimmed_coverage = (double) (querypos2 - querypos1 + 1)/(double) (trim_end - trim_start + skiplength);
    }
  }

  return;
}


static List_T
pick_cdna_direction (int *winning_cdna_direction, int *sensedir,
		     List_T pairs_fwd, List_T pairs_rev, double defect_rate_fwd, double defect_rate_rev,
		     int nknown_fwd, int ncanonical_fwd, int nsemicanonical_fwd,
		     int nnoncanonical_fwd, int nbadintrons_fwd,
		     int nknown_rev, int ncanonical_rev, int nsemicanonical_rev,
		     int nnoncanonical_rev, int nbadintrons_rev,
		     double max_intron_score_fwd, double avg_donor_score_fwd, double avg_acceptor_score_fwd,
		     double max_intron_score_rev, double avg_donor_score_rev, double avg_acceptor_score_rev,
#ifdef COMPLEX_DIRECTION
		     int nmatches_fwd, int nmismatches_fwd, int nmatches_rev, int nmismatches_rev, int nindels_fwd, int nindels_rev,
		     int indel_alignment_score_fwd, int indel_alignment_score_rev,
#endif
		     int sense_filter) {
#if 0
  int canonical_score_fwd, canonical_score_rev;
#endif

  if (pairs_fwd) {
    /* canonical_score_fwd = ncanonical_fwd - nbadintrons_fwd + nsemicanonical_fwd - nnoncanonical_fwd; */
    debug11(printf("nknown_fwd %d, ncanonical_fwd %d, nbadintrons_fwd %d, nsemicanonical_fwd %d, nnoncanonical_fwd %d\n",
		   nknown_fwd,ncanonical_fwd,nbadintrons_fwd,nsemicanonical_fwd,nnoncanonical_fwd));
  }
  if (pairs_rev) {
    /* canonical_score_rev = ncanonical_rev - nbadintrons_rev + nsemicanonical_rev - nnoncanonical_rev; */
    debug11(printf("nknown_rev %d, ncanonical_rev %d, nbadintrons_rev %d, nsemicanonical_rev %d, nnoncanonical_rev %d\n",
		   nknown_rev,ncanonical_rev,nbadintrons_rev,nsemicanonical_rev,nnoncanonical_rev));
  }

  if (pairs_fwd == NULL && pairs_rev == NULL) {
    debug11(printf("pairs_fwd is NULL and pairs_rev is NULL\n"));
    *winning_cdna_direction = 0;
    *sensedir = SENSE_NULL;
    return (List_T) NULL;

  } else if (pairs_rev == NULL) {
    debug11(printf("pairs_rev is NULL, so fwd wins\n"));
    *winning_cdna_direction = +1;
    *sensedir = SENSE_FORWARD;
    return pairs_fwd;

  } else if (pairs_fwd == NULL) {
    debug11(printf("pairs_fwd is NULL, so rev wins\n"));
    *winning_cdna_direction = -1;
    *sensedir = SENSE_ANTI;
    return pairs_rev;

#if 0
  } else if (indel_alignment_score_fwd >= 0 && indel_alignment_score_rev < 0 && nbadintrons_rev > 0) {
    debug11(printf("indel_alignment_score_fwd %d positive and indel_alignment_score_rev %d negative for a bad intron, so fwd wins\n",
		   indel_alignment_score_fwd,indel_alignment_score_rev));
    *winning_cdna_direction = +1;

  } else if (indel_alignment_score_rev >= 0 && indel_alignment_score_fwd < 0 && nbadintrons_fwd > 0) {
    debug11(printf("indel_alignment_score_fwd %d negative for a bad intron and indel_alignment_score_rev %d positive, so rev wins\n",
		   indel_alignment_score_fwd,indel_alignment_score_rev));
    *winning_cdna_direction = -1;
#endif

#if 0
    /* Cannot use, because favors a terminal over a splice */
  } else if (nmismatches_fwd < nmismatches_rev) {
    debug11(printf("nmismatches fwd %d < nmismatches rev %d, so fwd wins\n",
		   nmismatches_fwd,nmismatches_rev));
    *winning_cdna_direction = +1;

  } else if (nmismatches_fwd > nmismatches_rev) {
    debug11(printf("nmismatches fwd %d > nmismatches rev %d, so rev wins\n",
		   nmismatches_fwd,nmismatches_rev));
    *winning_cdna_direction = -1;
#endif

  } else if (defect_rate_fwd > DEFECT_MEDQ && defect_rate_rev > DEFECT_MEDQ &&
	     avg_donor_score_fwd > 0.9 && avg_donor_score_rev < 0.5 &&
	     avg_acceptor_score_fwd > 0.9 && avg_acceptor_score_rev < 0.5) {
    debug11(printf("defect_rate %f, %f and intronscores fwd %f,%f > intronscores rev %f,%f, so fwd wins\n",
		   defect_rate_fwd,defect_rate_rev,avg_donor_score_fwd,avg_acceptor_score_fwd,
		   avg_donor_score_rev,avg_acceptor_score_rev));
    /* intronscores reveal a clear sensedir */
    *winning_cdna_direction = +1;

  } else if (defect_rate_fwd > DEFECT_MEDQ && defect_rate_rev > DEFECT_MEDQ &&
	     avg_donor_score_rev > 0.9 && avg_donor_score_fwd < 0.5 &&
	     avg_acceptor_score_rev > 0.9 && avg_acceptor_score_fwd < 0.5) {
    debug11(printf("defect_rate %f, %f and intronscores rev %f,%f > intronscores fwd %f,%f, so fwd wins\n",
		   defect_rate_fwd,defect_rate_rev,avg_donor_score_rev,avg_acceptor_score_rev,
		   avg_donor_score_fwd,avg_acceptor_score_fwd));
    /* intronscores reveal a clear sensedir */
    *winning_cdna_direction = -1;

  } else if (nknown_fwd > 0 && nknown_rev == 0) {
    debug11(printf("nknown_fwd %d && nknown_rev %d, so fwd wins\n",
		   nknown_fwd,nknown_rev));
    *winning_cdna_direction = +1;

  } else if (nknown_rev > 0 && nknown_fwd == 0) {
    debug11(printf("nknown_fwd %d && nknown_rev %d, so rev wins\n",
		   nknown_fwd,nknown_rev));
    *winning_cdna_direction = -1;

  } else if (ncanonical_fwd > 0 && ncanonical_rev == 0) {
    debug11(printf("ncanonical_fwd %d && ncanonical_rev %d, so fwd wins\n",
		   ncanonical_fwd,ncanonical_rev));
    *winning_cdna_direction = +1;

  } else if (ncanonical_rev > 0 && ncanonical_fwd == 0) {
    debug11(printf("ncanonical_fwd %d && ncanonical_rev %d, so rev wins\n",
		   ncanonical_fwd,ncanonical_rev));
    *winning_cdna_direction = -1;

#if 0
  } else if (canonical_score_fwd > canonical_score_rev + 1) {
    debug11(printf("canonical_score_fwd %d > canonical_score_rev %d + 1, so fwd wins\n",
		   canonical_score_fwd,canonical_score_rev));
    *winning_cdna_direction = +1;

  } else if (canonical_score_rev > canonical_score_fwd + 1) {
    debug11(printf("canonical_score_rev %d > canonical_score_fwd %d + 1, so rev wins\n",
		   canonical_score_rev,canonical_score_fwd));
    *winning_cdna_direction = -1;
#endif

#if 0
  } else if (alignment_score_fwd > alignment_score_rev + SCORE_SIGDIFF) {
    debug11(printf("alignment_score_fwd %d >> alignment_score_rev %d, so fwd wins\n",
		   alignment_score_fwd,alignment_score_rev));
    *winning_cdna_direction = +1;

  } else if (alignment_score_rev > alignment_score_fwd + SCORE_SIGDIFF) {
    debug11(printf("alignment_score_rev %d << alignment_score_fwd %d, so rev wins\n",
		   alignment_score_rev,alignment_score_fwd));
    *winning_cdna_direction = -1;
#endif

  } else if (nnoncanonical_fwd == 0 && nnoncanonical_rev > 0) {
    debug11(printf("nnoncanonical_fwd %d < nnoncanonical_rev %d, so fwd wins\n",
		   nnoncanonical_fwd,nnoncanonical_rev));
    *winning_cdna_direction = +1;

  } else if (nnoncanonical_rev == 0 && nnoncanonical_fwd > 0) {
    debug11(printf("nnoncanonical_rev %d < nnoncanonical_fwd %d, so rev wins\n",
		   nnoncanonical_rev,nnoncanonical_fwd));
    *winning_cdna_direction = -1;

  } else if (nbadintrons_fwd == 0 && nbadintrons_rev > 0) {
    debug11(printf("nbadintrons_fwd %d < nbadintrons_rev %d, so fwd wins\n",
		   nbadintrons_fwd,nbadintrons_rev));
    *winning_cdna_direction = +1;

  } else if (nbadintrons_rev == 0 && nbadintrons_fwd > 0) {
    debug11(printf("nbadintrons_rev %d < nbadintrons_fwd %d, so rev wins\n",
		   nbadintrons_rev,nbadintrons_fwd));
    *winning_cdna_direction = -1;

  } else if (avg_donor_score_fwd > avg_donor_score_rev + PROB_SIGDIFF && 
	     avg_acceptor_score_fwd > avg_acceptor_score_rev + PROB_SIGDIFF) {
    debug11(printf("intronscores fwd %f+%f > intronscores rev %f+%f, so fwd wins\n",
		   avg_donor_score_fwd,avg_acceptor_score_fwd,avg_donor_score_rev,avg_acceptor_score_rev));
    /* intronscores reveal a preferred sensedir */
    *winning_cdna_direction = +1;

  } else if (avg_donor_score_rev > avg_donor_score_fwd + PROB_SIGDIFF &&
	     avg_acceptor_score_rev > avg_acceptor_score_fwd + PROB_SIGDIFF) {
    debug11(printf("intronscores rev %f+%f > intronscores fwd %f+%f, so fwd wins\n",
		   avg_donor_score_rev,avg_acceptor_score_rev,avg_donor_score_fwd,avg_acceptor_score_fwd));
    /* intronscores reveal a preferred sensedir */
    *winning_cdna_direction = -1;

#if 0
  } else if (alignment_score_fwd > alignment_score_rev && alignment_score_fwd > 0) {
    debug11(printf("alignment_score_fwd %d > alignment_score_rev %d, so fwd wins\n",
		   alignment_score_fwd,alignment_score_rev));
    *winning_cdna_direction = +1;

  } else if (alignment_score_rev > alignment_score_fwd && alignment_score_rev > 0) {
    debug11(printf("alignment_score_rev %d < alignment_score_fwd %d, so rev wins\n",
		   alignment_score_rev,alignment_score_fwd));
    *winning_cdna_direction = -1;
#endif

  } else {
    debug11(printf("scores all equal, so fwd wins, but setting cdna_direction to be 0\n"));
    /* No clear intron direction, so allow under all sense_filters */
    *winning_cdna_direction = 0; 
    *sensedir = SENSE_NULL;
    return pairs_fwd;
  }

  debug11(printf("max_intron_score_fwd = %f, max_intron_score_rev = %f\n",max_intron_score_fwd,max_intron_score_rev));

  if (*winning_cdna_direction == +1) {
    if (ncanonical_fwd == 0 && nsemicanonical_fwd == 0 && nnoncanonical_fwd == 0) {
      *sensedir = SENSE_NULL;
    } else if (max_intron_score_fwd < 1.8) {
      *sensedir = SENSE_NULL;
    } else {
      *sensedir = SENSE_FORWARD;
    }
#ifndef PMAP
    if (sense_filter < 0) {
      return (List_T) NULL;
    }
#endif
    debug11(printf("winning_cdna_direction = %d, sensedir = %d\n",
		   *winning_cdna_direction,*sensedir));
    return pairs_fwd;

  } else if (*winning_cdna_direction == -1) {
    if (ncanonical_rev == 0 && nsemicanonical_rev == 0 && nnoncanonical_rev == 0) {
      *sensedir = SENSE_NULL;
    } else if (max_intron_score_rev < 1.8) {
      *sensedir = SENSE_NULL;
    } else {
      *sensedir = SENSE_ANTI;
    }
#ifndef PMAP
    if (sense_filter > 0) {
      return (List_T) NULL;
    }
#endif
    debug11(printf("winning_cdna_direction = %d, sensedir = %d\n",
		   *winning_cdna_direction,*sensedir));
    return pairs_rev;

  } else {
    fprintf(stderr,"Unexpected value %d for winning_cdna_direction\n",*winning_cdna_direction);
    abort();
  }
}


static int
initial_cdna_direction (List_T pairs_fwd, List_T pairs_rev,
			double avg_donor_score_fwd, double avg_acceptor_score_fwd,
			double avg_donor_score_rev, double avg_acceptor_score_rev
#ifdef COMPLEX_DIRECTION
			, int nmatches_fwd, int nmismatches_fwd, int nmatches_rev, int nmismatches_rev, int nindels_fwd, int nindels_rev,
			int indel_alignment_score_fwd, int indel_alignment_score_rev
#endif
			) {

  if (pairs_fwd == NULL && pairs_rev == NULL) {
    debug11(printf("pairs_fwd is NULL and pairs_rev is NULL\n"));
    return 0;

  } else if (pairs_rev == NULL) {
    debug11(printf("pairs_rev is NULL, so fwd wins\n"));
    return +1;

  } else if (pairs_fwd == NULL) {
    debug11(printf("pairs_fwd is NULL, so rev wins\n"));
    return -1;

  } else if (avg_donor_score_fwd > 0.9 && avg_acceptor_score_fwd > 0.9 &&
	     (avg_donor_score_rev < 0.5 || avg_acceptor_score_rev < 0.5)) {
    debug11(printf("intronscores fwd %f,%f > intronscores rev %f,%f, so fwd wins\n",
		   avg_donor_score_fwd,avg_acceptor_score_fwd,avg_donor_score_rev,avg_acceptor_score_rev));
    /* intronscores reveal a clear sensedir */
    return +1;

  } else if (avg_donor_score_rev > 0.9 && avg_acceptor_score_rev > 0.9 &&
	     (avg_donor_score_fwd < 0.5 || avg_acceptor_score_fwd < 0.5)) {
    debug11(printf("intronscores rev %f,%f > intronscores fwd %f,%f, so fwd wins\n",
		   avg_donor_score_rev,avg_acceptor_score_rev,avg_donor_score_fwd,avg_acceptor_score_fwd));
    /* intronscores reveal a clear sensedir */
    return -1;

#if 0
  } else if (alignment_score_fwd > alignment_score_rev + SCORE_SIGDIFF) {
    /* Cannot use alignment score until we do a full alignment */
    debug11(printf("alignment_score_fwd %d >> alignment_score_rev %d, so fwd wins\n",
		   alignment_score_fwd,alignment_score_rev));
    return +1;

  } else if (alignment_score_rev > alignment_score_fwd + SCORE_SIGDIFF) {
    debug11(printf("alignment_score_rev %d << alignment_score_fwd %d, so rev wins\n",
		   alignment_score_rev,alignment_score_fwd));
    return -1;
#endif

#if 0
  } else if (nnoncanonical_fwd < nnoncanonical_rev) {
    /* Not a good test until we do full dynamic programming */
    debug11(printf("nnoncanonical_fwd %d < nnoncanonical_rev %d, so fwd wins\n",
		   nnoncanonical_fwd,nnoncanonical_rev));
    return +1;

  } else if (nnoncanonical_rev < nnoncanonical_fwd) {
    debug11(printf("nnoncanonical_rev %d < nnoncanonical_fwd %d, so rev wins\n",
		   nnoncanonical_rev,nnoncanonical_fwd));
    return -1;
#endif

#if 0
  } else if (avg_donor_score_fwd + avg_acceptor_score_fwd > avg_donor_score_rev + avg_acceptor_score_rev + PROB_SIGDIFF) {
    /* Not a good test until we do full dynamic programming */
    debug11(printf("intronscores fwd %f+%f > intronscores rev %f+%f, so fwd wins\n",
		   avg_donor_score_fwd,avg_acceptor_score_fwd,avg_donor_score_rev,avg_acceptor_score_rev));
    /* intronscores reveal a preferred sensedir */
    return +1;

  } else if (avg_donor_score_rev + avg_acceptor_score_rev > avg_donor_score_fwd + avg_acceptor_score_fwd + PROB_SIGDIFF) {
    debug11(printf("intronscores rev %f+%f > intronscores fwd %f+%f, so fwd wins\n",
		   avg_donor_score_rev,avg_acceptor_score_rev,avg_donor_score_fwd,avg_acceptor_score_fwd));
    /* intronscores reveal a preferred sensedir */
    return -1;
#endif

  } else {
    return 0;
  }
}



T
Stage3_new (struct Pair_T *pairarray, List_T pairs, int npairs, int goodness, int cdna_direction, int sensedir,
	    int matches, int unknowns, int mismatches, int qopens, int qindels,
	    int topens, int tindels, int ncanonical, int nsemicanonical, int nnoncanonical,
	    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
	    bool watsonp, int genestrand, int querylength, int skiplength, int trimlength,
	    int straintype, char *strain, IIT_T altstrain_iit) {
  T new;
  Pair_T start, end;
  int *typematches, nmatches;
  int alias;

  List_T cigar_tokens;
  bool intronp;
  int hardclip_start, hardclip_end;


  /* pairs can be NULL from Stage3_split */
  if (pairs == NULL || npairs == 0) {
    return (T) NULL;
  }

  start = &(pairarray[0]);
  end = &(pairarray[npairs-1]);
  assert(start->gapp == false);
  assert(end->gapp == false);

  hardclip_start = start->querypos;
  hardclip_end = (querylength - 1) - end->querypos;

  cigar_tokens = Pair_compute_cigar(&intronp,&hardclip_start,&hardclip_end,pairarray,npairs,querylength,
				    watsonp,/*chimera_part*/0);
  if (Pair_cigar_length(cigar_tokens) + hardclip_start + hardclip_end != querylength) {
    fprintf(stderr,"Could not compute a valid cigar from the following alignment: %d + %d + %d != %d\n",
	    Pair_cigar_length(cigar_tokens),hardclip_start,hardclip_end,querylength);
    Pair_dump_array_stderr(pairarray,npairs,/*zerobasedp*/true);
#ifdef CHECK_ASSERTIONS
    abort();
#endif
    Pair_tokens_free(&cigar_tokens);
    return (T) NULL;

  } else {
    Pair_tokens_free(&cigar_tokens);

    new = (T) MALLOC_OUT(sizeof(*new)); /* Matches FREE_OUT in Stage3_free */
    debug99(printf("Creating %p\n",new));
  }

  new->pairarray = pairarray;
  new->pairarray_freeable_p = true;
  new->chimera_left_p = false;
  new->chimera_right_p = false;

  new->pairs = pairs;
  new->npairs = npairs;

  new->matches = matches;
  new->unknowns = unknowns;
  new->mismatches = mismatches;
  new->qopens = qopens;
  new->qindels = qindels;
  new->topens = topens;
  new->tindels = tindels;

  new->noncanonical = nsemicanonical + nnoncanonical;
  new->goodness = goodness;

#ifdef PMAP
  /* Should be +1 */
  new->cdna_direction = cdna_direction;
  new->sensedir = sensedir;
#else
  if (cdna_direction == 0 && require_splicedir_p == true) {
    new->cdna_direction = Pair_guess_cdna_direction_array(&new->sensedir,pairarray,npairs,/*invertedp*/false,
							  chroffset,watsonp);
  } else if (ncanonical == 0 && nsemicanonical == 0 /*&& nnoncanonical == 0*/) {
    new->cdna_direction = 0;
    new->sensedir = SENSE_NULL;	/* was sensedir, but this gives bad XS output */
  } else {
    new->cdna_direction = cdna_direction;
    new->sensedir = sensedir;
  }
#endif

#if 0
  nexons = Pair_nexons_approx(pairs);
  if (nexons > 2) {
    /* Favor spliced transcripts, but only if we're sure they're
       spliced (i.e., 3 or more exons).  A random intron shouldn't
       get credit. */
    new->goodness += nexons;
  }
#endif
    
  new->translation_start = 0;
  new->translation_end = 0;
  new->translation_length = 0;

  /* new->stage2_source = stage2_source; */
  /* new->stage2_indexsize = stage2_indexsize; */

  new->straintype = straintype;
  new->strain = strain;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->circularpos = Pair_circularpos(&alias,pairarray,npairs,chrlength,watsonp,querylength);

  new->watsonp = watsonp;
  new->genestrand = genestrand;

  new->genomicstart = chroffset + Pair_genomepos(start);
  new->genomicend = chroffset + Pair_genomepos(end);

  /* new->stage3_runtime = stage3_runtime; */

  new->trimmed_coverage = (double) (end->querypos - start->querypos + 1)/(double) (trimlength + skiplength);

  debug0(printf("Creating stage3 at chr %d:%u..%u, goodness %d, matches %d, npairs %d\n",
		chrnum,Stage3_chrstart(new),Stage3_chrend(new),new->goodness,new->matches,new->npairs));

  if (straintype == 0) {
    return new;
  } else {
    if (watsonp) {
      typematches = IIT_get_typed(&nmatches,altstrain_iit,/*divstring*/NULL,
				  new->genomicstart,new->genomicend,straintype,/*sortp*/false);
    } else {
      typematches = IIT_get_typed(&nmatches,altstrain_iit,/*divstring*/NULL,
				  new->genomicend,new->genomicstart,straintype,/*sortp*/false);
    }
    if (typematches == NULL) {
      Stage3_free(&new);
      return NULL;
    } else {
      FREE(typematches);
      return new;
    }
  }

  return new;
}


T
Stage3_new_from_pairs (List_T pairs, int cdna_direction, bool watsonp, int genestrand, int sensedir,
		       Pairpool_T pairpool, Sequence_T queryseq, int query_subseq_offset,
		       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength) {
  struct Pair_T *pairarray;
  int npairs;
  int goodness, matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double min_splice_prob;

  if (pairs == NULL) {
    return (T) NULL;
  } else {
    /* Use chroffset+chrlength instead of chrhigh, because we are calling this procedure after Pair_split_circular */
    pairarray = make_pairarray(&npairs,&pairs,cdna_direction,watsonp,
			       pairpool,/*queryseq_ptr*/Sequence_fullpointer(queryseq),
			       chroffset,/*chrhigh*/chroffset+chrlength,
			       ngap,query_subseq_offset,/*skiplength*/Sequence_skiplength(queryseq));

    goodness = Pair_fracidentity_array(&matches,&unknowns,&mismatches,
				       &qopens,&qindels,&topens,&tindels,
				       &ncanonical,&nsemicanonical,&nnoncanonical,
				       &min_splice_prob,pairarray,npairs,cdna_direction);
    
    return Stage3_new(pairarray,pairs,npairs,goodness,cdna_direction,sensedir,
		      matches,unknowns,mismatches,
		      qopens,qindels,topens,tindels,ncanonical,nsemicanonical,nnoncanonical,
		      chrnum,chroffset,chrhigh,chrlength,watsonp,genestrand,
		      /*querylength*/Sequence_fulllength(queryseq),
		      /*skiplength*/Sequence_skiplength(queryseq),
		      /*trimlength*/Sequence_trimlength(queryseq),
		      /*straintype*/0,/*strain*/NULL,/*altstrain_iit*/NULL);
  }
}


#if 0
void
Stage3_check_for_shared_pairs (T this, T that) {
  List_T p, q;
  Pair_T pair1, pair2;

  if (this == that) {
    return;
  } else {
    for (p = this->pairs; p != NULL; p = List_next(p)) {
      for (q = that->pairs; q != NULL; q = List_next(q)) {
	if (List_head(p) == List_head(q)) {
	  pair1 = (Pair_T) List_head(p);
	  pair2 = (Pair_T) List_head(q);
	  Pair_dump_one(pair1,true);
	  printf("\n");
	  Pair_dump_one(pair2,true);
	  printf("\n");
	  abort();
	}
      }
    }
  }

  return;
}


void
Stage3_list_check_for_shared_pairs (List_T stage3list) {
  List_T p, q, r, s;
  Pair_T pair1, pair2;

  for (r = stage3list; r != NULL; r = List_next(r)) {
    for (s = stage3list; s != NULL; s = List_next(s)) {
      Stage3_check_for_shared_pairs(List_head(r),List_head(s));
    }
  }

  return;
}
#endif


#if 0
static T
Stage3_copy (T old, Pairpool_T pairpool) {
  T new = (T) MALLOC_OUT(sizeof(*new));

  new->pairarray = (struct Pair_T *) MALLOC_OUT(old->npairs*sizeof(struct Pair_T));
  memcpy(new->pairarray,old->pairarray,old->npairs*sizeof(struct Pair_T));
  new->npairs = old->npairs;
  new->pairarray_freeable_p = true;

  new->chimera_left_p = false;
  new->chimera_right_p = false;

  new->pairs = Pairpool_copy(old->pairs,pairpool);
  new->matches = old->matches;
  new->unknowns = old->unknowns;
  new->mismatches = old->mismatches;
  new->qopens = old->qopens;
  new->qindels = old->qindels;
  new->topens = old->topens;
  new->tindels = old->tindels;

  new->noncanonical = old->noncanonical;
  new->goodness = old->goodness;

  new->cdna_direction = old->cdna_direction;
  new->sensedir = old->sensedir;

  new->translation_start = old->translation_start;
  new->translation_end = old->translation_end;
  new->translation_length = old->translation_length;

  new->straintype = old->straintype;
  new->strain = old->strain;

  new->chrnum = old->chrnum;
  new->chroffset = old->chroffset;
  new->chrhigh = old->chrhigh;
  new->chrlength = old->chrlength;
  new->circularpos = old->circularpos;

  new->watsonp = old->watsonp;
  new->genestrand = old->genestrand;

  new->genomicstart = old->genomicstart;
  new->genomicend = old->genomicend;

  /* new->stage3_runtime = stage3_runtime; */

  new->trimmed_coverage = old->trimmed_coverage;

  return new;
}
#endif


void
Stage3_free (T *old) {

  debug99(printf("Freeing %p\n",*old));
  if (*old) {
    /* Don't free strain.  Belongs to altstrain_iit. */
    if ((*old)->pairarray_freeable_p == true) {
      FREE_OUT((*old)->pairarray);
    }
    FREE_OUT(*old);
  }
  return;
}

#if 0
/* Needed for mutation analysis in align_relative */
void
Stage3_genomicbounds (Univcoord_T *genomicstart, Univcoord_T *genomiclength, T this) {
  *genomicstart = this->chroffset;
  *genomiclength = this->genomiclength;
  return;
}
#endif


bool
Stage3_test_bounds (T this, int minpos, int maxpos) {
  int nstart;

  if (Pairpool_count_bounded(&nstart,this->pairs,minpos,maxpos) >= 25) {
    return true;
  } else {
    return false;
  }
}


#ifdef PMAP
void
Stage3_translate_cdna (T this, Sequence_T queryaaseq, bool strictp) {
  Translation_via_cdna(&this->translation_start,&this->translation_end,&this->translation_length,
		       &this->relaastart,&this->relaaend,
		       this->pairarray,this->npairs,Sequence_fullpointer(queryaaseq),strictp);
  return;
}

void
Stage3_backtranslate_cdna (T this) {
  Backtranslation_cdna(this->pairarray,this->npairs,this->translation_start,this->translation_end);
  return;
}

#elif defined(GSNAP)

/* No need to perform translation */

#else

static void
truncate_fulllength (Stage3_T this, bool translatep, int cds_startpos, int querylength, bool strictp) {

  if (translatep == true) {
    if (this->cdna_direction < 0) {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,this->npairs,/*backwardsp*/true,/*revcompp*/true,/*fulllengthp*/true,
			      cds_startpos,querylength,strictp);
    } else {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,this->npairs,/*backwardsp*/false,/*revcompp*/false,/*fulllengthp*/true,
			      cds_startpos,querylength,strictp);
    }
  }

  this->npairs = Pair_clip_bounded_array(this->pairarray,this->npairs,
					 /*minpos*/Stage3_translation_start(this),
					 /*maxpos*/Stage3_translation_end(this));
  return;
}


void
Stage3_translate_genomic (T this, int npairs, bool fulllengthp, int cds_startpos, int querylength, bool truncatep, bool strictp) {

  if (this->cdna_direction < 0) {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairarray,npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp,
			    cds_startpos,querylength,strictp);
  } else {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairarray,npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp,
			    cds_startpos,querylength,strictp);
  }
  if (truncatep == true) {
    truncate_fulllength(this,/*translatep*/false,cds_startpos,querylength,strictp);
    if (this->cdna_direction < 0) {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp,
			      cds_startpos,querylength,strictp);
    } else {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp,
			      cds_startpos,querylength,strictp);
    }
  }

  return;
}
#endif


#if 0
/* Called previously by align_relative in gmap.c */
#ifndef GSNAP
void
Stage3_translate_cdna_via_reference (T this, T reference) {
  /* bool fixshiftp = !literalrefp; */

  if (this->watsonp == reference->watsonp) {
    if (reference->cdna_direction < 0) {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairarray,this->npairs,this->watsonp,/*backwardsp*/true,/*revcompp*/true,
				reference->pairarray,reference->npairs,reference->watsonp);
    } else {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairarray,this->npairs,this->watsonp,/*backwardsp*/false,/*revcompp*/false,
				reference->pairarray,reference->npairs,reference->watsonp);
    }
  } else {
    if (reference->cdna_direction < 0) {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairarray,this->npairs,this->watsonp,/*backwardsp*/false,/*revcompp*/false,
				reference->pairarray,reference->npairs,reference->watsonp);
    } else {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairarray,this->npairs,this->watsonp,/*backwardsp*/true,/*revcompp*/true,
				reference->pairarray,reference->npairs,reference->watsonp);
    }
  }

  return;
}
#endif
#endif


#if 0
/* Called previously by align_relative in gmap.c */
void
Stage3_fix_cdna_direction (T this, T reference) {
  if (this->cdna_direction == 0) {
    if (reference->cdna_direction > 0) {
      if (this->watsonp == reference->watsonp) {
	this->cdna_direction = +1;
      } else {
	this->cdna_direction = -1;
      }
    } else if (reference->cdna_direction < 0) {
      if (this->watsonp == reference->watsonp) {
	this->cdna_direction = -1;
      } else {
	this->cdna_direction = +1;
      }
    }
  }
  return;
}
#endif



#ifndef GSNAP
void
Stage3_translate (T this,
#ifdef PMAP
		  Sequence_T queryseq,
#endif
		  int querylength, bool fulllengthp,
		  int cds_startpos, bool truncatep, bool strictp) {

#ifdef PMAP
  Translation_via_cdna(&this->translation_start,&this->translation_end,&this->translation_length,
		       &this->relaastart,&this->relaaend,
		       this->pairarray,this->npairs,Sequence_fullpointer(queryseq),strictp);
  Backtranslation_cdna(this->pairarray,this->npairs,this->translation_start,this->translation_end);
#else
  if (this->cdna_direction < 0) {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairarray,this->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp,
			    cds_startpos,querylength,strictp);
  } else {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairarray,this->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp,
			    cds_startpos,querylength,strictp);
  }

  if (truncatep == true) {
    truncate_fulllength(this,/*translatep*/false,cds_startpos,querylength,strictp);
    if (this->cdna_direction < 0) {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,this->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp,
			      cds_startpos,querylength,strictp);
    } else {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,this->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp,
			      cds_startpos,querylength,strictp);
    }
  }
#endif

  return;
}
#endif


#ifndef GSNAP
void
Stage3_translate_chimera (T this, T mate,
#ifdef PMAP
			  Sequence_T queryseq,
#endif
			  int querylength, bool fulllengthp,
			  int cds_startpos, bool truncatep, bool strictp) {
  int npairs1, npairs2;
  int translation_start, translation_end, translation_length, relaastart, relaaend;

  npairs1 = this->npairs;
  npairs2 = mate->npairs;

#ifdef PMAP
  Translation_via_cdna(&translation_start,&translation_end,&translation_length,
		       &relaastart,&relaaend,
		       this->pairarray,npairs1 + npairs2,Sequence_fullpointer(queryseq),strictp);
  Backtranslation_cdna(this->pairarray,npairs1 + npairs2,translation_start,translation_end);
#else
  if (this->cdna_direction < 0) {
    Translation_via_genomic(&translation_start,&translation_end,&translation_length,
			    &relaastart,&relaaend,
			    this->pairarray,npairs1 + npairs2,/*backwardsp*/true,/*revcompp*/true,fulllengthp,
			    cds_startpos,querylength,strictp);
  } else {
    Translation_via_genomic(&translation_start,&translation_end,&translation_length,
			    &relaastart,&relaaend,
			    this->pairarray,npairs1 + npairs2,/*backwardsp*/false,/*revcompp*/false,fulllengthp,
			    cds_startpos,querylength,strictp);
  }

  if (truncatep == true) {
    truncate_fulllength(this,/*translatep*/false,cds_startpos,querylength,strictp);
    if (this->cdna_direction < 0) {
      Translation_via_genomic(&translation_start,&translation_end,&translation_length,
			      &relaastart,&relaaend,
			      this->pairarray,npairs1 + npairs2,/*backwardsp*/true,/*revcompp*/true,fulllengthp,
			      cds_startpos,querylength,strictp);
    } else {
      Translation_via_genomic(&translation_start,&translation_end,&translation_length,
			      &relaastart,&relaaend,
			      this->pairarray,npairs1 + npairs2,/*backwardsp*/false,/*revcompp*/false,fulllengthp,
			      cds_startpos,querylength,strictp);
    }
  }

#endif

  if (translation_start < npairs1) {
    this->translation_start = translation_start;
    mate->translation_start = 0;
  } else {
    this->translation_start = npairs1 - 1;
    mate->translation_start = translation_start - npairs1;
  }
  if (translation_end < npairs1) {
    this->translation_end = translation_end;
    mate->translation_end = 0;
  } else {
    this->translation_end = npairs1 - 1;
    mate->translation_end = translation_end - npairs1;
  }

  /* Additional checks to stay within array bounds */
  if (this->translation_end >= this->npairs) {
    this->translation_end = this->npairs - 1;
  }
  if (this->translation_start > this->translation_end) {
    this->translation_start = this->translation_end;
  }

  if (mate->translation_end >= mate->npairs) {
    mate->translation_end = mate->npairs - 1;
  }
  if (mate->translation_start > mate->translation_end) {
    mate->translation_start = mate->translation_end;
  }

  debug(printf("Converted translation %d..%d in %d+%d pairs to %d..%d and %d..%d\n",
	       translation_start,translation_end,this->npairs,mate->npairs,
	       this->translation_start,this->translation_end,mate->translation_start,mate->translation_end));

  this->translation_length = Pair_translation_length(this->pairarray,this->npairs);
  mate->translation_length = Pair_translation_length(mate->pairarray,mate->npairs);
  debug(printf("Original translation length %d => %d plus %d\n",
	       translation_length,this->translation_length,mate->translation_length));

  this->relaastart = this->pairarray[this->translation_start].aapos;
  this->relaaend = this->pairarray[this->translation_end].aapos;

  mate->relaastart = mate->pairarray[mate->translation_start].aapos;
  mate->relaaend = mate->pairarray[mate->translation_end].aapos;

  return;
}
#endif


#ifndef GSNAP
void
Stage3_print_pathsummary (Filestring_T fp, T this, int pathnum, Univ_IIT_T chromosome_iit, Univ_IIT_T contig_iit, 
			  IIT_T altstrain_iit, Sequence_T queryseq, char *dbversion) {
  Pair_T start, end;
  bool referencealignp;

  debug99(printf("Printing %p\n",this));
  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);
  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(fp,pathnum,start,end,this->chrnum,this->chroffset,
			 chromosome_iit,referencealignp,altstrain_iit,this->strain,contig_iit,
			 dbversion,Sequence_fulllength_given(queryseq),Sequence_skiplength(queryseq),
			 Sequence_trim_start(queryseq),Sequence_trim_end(queryseq),
			 Pair_nexons(this->pairarray,this->npairs),this->matches,this->unknowns,this->mismatches,
			 this->qopens,this->qindels,this->topens,this->tindels,
			 this->watsonp,this->cdna_direction,
			 this->translation_start,this->translation_end,this->translation_length,
			 /*relaastart*/0,/*relaaend*/0);
  Translation_print_comparison(fp,this->pairarray,this->npairs,this->relaastart,this->relaaend);
  FPRINTF(fp,"\n");

  return;
}
#endif


void
Stage3_print_pslformat_nt (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryaaseq) {
  Pair_T start, end;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_pslformat_nt(fp,this->pairarray,this->npairs,start,end,queryaaseq,this->chrnum,
			  chromosome_iit,usersegment,
			  /* Pair_nexons(this->pairarray,this->npairs), */
			  this->matches,this->unknowns,this->mismatches, 
			  this->watsonp);
  return;
}

#ifdef PMAP
void
Stage3_print_pslformat_pro (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryaaseq, bool strictp) {
  Pair_T start, end;

#if 0
  Stage3_translate_cdna(this,queryaaseq,strictp);
  Stage3_backtranslate_cdna(this);
#endif

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_pslformat_pro(fp,this->pairarray,this->npairs,start,end,queryaaseq,this->chrnum,
			   chromosome_iit,usersegment,
			   /* Pair_nexons(this->pairarray,this->npairs), */
			   this->watsonp,this->cdna_direction);
  return;
}
#endif


void
Stage3_print_gff3 (Filestring_T fp, T this, int pathnum, Univ_IIT_T chromosome_iit, Sequence_T usersegment,
		   Sequence_T queryseq, int querylength, Printtype_T printtype, char *sourcename) {
  Pair_T start, end;
  bool gff_gene_format_p, gff_estmatch_format_p;

  if (printtype == GFF3_GENE) {
    gff_gene_format_p = true;
    gff_estmatch_format_p = false;
  } else if (printtype == GFF3_MATCH_CDNA) {
    gff_gene_format_p = false;
    gff_estmatch_format_p = false;
  } else if (printtype == GFF3_MATCH_EST) {
    gff_gene_format_p = false;
    gff_estmatch_format_p = true;
  } else {
    fprintf(stderr,"Unexpected printtype %d\n",printtype);
    abort();
  }

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_gff3(fp,this->pairarray,this->npairs,pathnum,
		  Sequence_accession(queryseq),start,end,
		  this->chrnum,chromosome_iit,usersegment,
		  this->translation_end,querylength,Sequence_skiplength(queryseq),
		  this->matches,this->mismatches,this->qindels,this->tindels,this->unknowns,
		  this->watsonp,this->cdna_direction,gff_gene_format_p,gff_estmatch_format_p,
		  sourcename);
  return;
}


#ifndef PMAP
void
Stage3_print_bedpe (Filestring_T fp, T this, Univ_IIT_T chromosome_iit) {
  Pair_print_bedpe(fp,this->pairarray,this->npairs,this->chrnum,this->watsonp,chromosome_iit);
  return;
}
#endif



#ifndef GSNAP
#ifndef PMAP
/* Only for GMAP program */
void
Stage3_print_sam (Filestring_T fp, char *abbrev, T this, int pathnum, int npaths_primary, int npaths_altloc,
		  int absmq_score, int second_absmq, int mapq_score,
		  Univ_IIT_T chromosome_iit, Sequence_T usersegment,
		  Sequence_T queryseq, int chimera_part, Chimera_T chimera,
		  int quality_shift, bool sam_paired_p, char *sam_read_group_id) {
  int querylength;
  Chrpos_T chrpos;
  Pair_T pair;

  querylength = Sequence_fulllength_given(queryseq);
  if (this->watsonp == true) {
    pair = &(this->pairarray[0]);
    chrpos = pair->genomepos + 1U;
  } else {
    pair = &(this->pairarray[this->npairs-1]);
    chrpos = pair->genomepos + 1U;
  }

  if (this->circularpos > 0) {
    Pair_print_sam(fp,abbrev,this->pairarray,this->npairs,
		   Sequence_accession(queryseq),/*acc2*/NULL,this->chrnum,chromosome_iit,usersegment,
		   Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
		   /*hardclip5*/0,/*hardclip3*/querylength-this->circularpos,
		   querylength,this->watsonp,this->sensedir,chimera_part,chimera,
		   quality_shift,Sequence_firstp(queryseq),
		   pathnum,npaths_primary,npaths_altloc,absmq_score,second_absmq,chrpos,this->chrlength,
		   mapq_score,sam_paired_p,sam_read_group_id);
    Pair_print_sam(fp,abbrev,this->pairarray,this->npairs,
		   Sequence_accession(queryseq),/*acc2*/NULL,this->chrnum,chromosome_iit,usersegment,
		   Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
		   /*hardclip5*/this->circularpos,/*hardclip3*/0,
		   querylength,this->watsonp,this->sensedir,chimera_part,chimera,
		   quality_shift,Sequence_firstp(queryseq),
		   pathnum,npaths_primary,npaths_altloc,absmq_score,second_absmq,/*chrpos*/1,this->chrlength,
		   mapq_score,sam_paired_p,sam_read_group_id);
  } else {
    Pair_print_sam(fp,abbrev,this->pairarray,this->npairs,
		   Sequence_accession(queryseq),/*acc2*/NULL,this->chrnum,chromosome_iit,usersegment,
		   Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
		   /*hardclip5*/0,/*hardclip3*/0,querylength,
		   this->watsonp,this->sensedir,chimera_part,chimera,
		   quality_shift,Sequence_firstp(queryseq),
		   pathnum,npaths_primary,npaths_altloc,absmq_score,second_absmq,chrpos,this->chrlength,
		   mapq_score,sam_paired_p,sam_read_group_id);
  }

  return;
}
#endif
#endif


void
Stage3_print_iit_map (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_T start, end;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_iit_map(fp,queryseq,Sequence_accession(queryseq),start,end,
		     this->chrnum,chromosome_iit);
  return;
}

void
Stage3_print_iit_exon_map (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_T start, end;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_iit_exon_map(fp,this->pairarray,this->npairs,queryseq,Sequence_accession(queryseq),
			  start,end,this->chrnum,chromosome_iit);
  return;
}

void
Stage3_print_splicesites (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_print_splicesites(fp,this->pairarray,this->npairs,Sequence_accession(queryseq),
			 Pair_nexons(this->pairarray,this->npairs),this->chrnum,
			 chromosome_iit,this->watsonp);
  return;
}

void
Stage3_print_introns (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_print_introns(fp,this->pairarray,this->npairs,Sequence_accession(queryseq),
		     Pair_nexons(this->pairarray,this->npairs),this->chrnum,
		     chromosome_iit);
  return;
}


#if 0
/* Called previously by align_relative in gmap.c */
#ifndef GSNAP
void
Stage3_print_mutations (Filestring_T fp, T this, T reference, Univ_IIT_T chromosome_iit, Sequence_T queryseq,
			char *dbversion, bool showalignp, int invertmode, bool nointronlenp, int wraplength) {
  Pair_T start, end;
  bool referencealignp;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  /*  Pair_dump_array(this->pairarray,this->npairs,false); */

  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(fp,/*pathnum*/1,start,end,reference->chrnum,reference->chroffset,
			 chromosome_iit,referencealignp,/*altstrain_iit*/NULL,this->strain,/*contig_iit*/NULL,
			 dbversion,Sequence_fulllength_given(queryseq),Sequence_skiplength(queryseq),
			 Sequence_trim_start(queryseq),Sequence_trim_end(queryseq),
			 Pair_nexons(this->pairarray,this->npairs),this->matches,this->unknowns,this->mismatches,
			 this->qopens,this->qindels,this->topens,this->tindels,
			 this->watsonp,this->cdna_direction,0,0,0,this->relaastart,this->relaaend);
  Translation_print_comparison(fp,this->pairarray,this->npairs,this->relaastart,this->relaaend);
  FPRINTF(fp,"\n");

  if (showalignp == true) {
    Pair_print_alignment(fp,this->pairarray,this->npairs,reference->chrnum,reference->chroffset,
			 chromosome_iit,this->watsonp,invertmode,nointronlenp,wraplength);
  }
  debug1(Pair_dump_array(this->pairarray,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairarray,this->npairs));

  return;
}
#endif
#endif



static void
print_map (Filestring_T fp, T this, IIT_T map_iit, int *map_divint_crosstable,
	   Univ_IIT_T chromosome_iit, int pathnum, bool map_bothstrands_p,
	   int nflanking, bool print_comment_p) {
  int chrlow, chrhigh;
  Pair_T start, end;
  int chrpos1, chrpos2;
  int *iit_matches = NULL, nmatches, *leftflanks, nleftflanks, *rightflanks, nrightflanks;
  int divno, sign;
  char *chr;

  if ((divno = map_divint_crosstable[this->chrnum]) <= 0) {
    FPRINTF(fp,"  *Map hits for path %d (0):\n\n",pathnum);
    return;
  } else {
    chr = Chrnum_to_string(this->chrnum,chromosome_iit);
  }

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  if (this->watsonp) {
    chrlow = chrpos1 = Pair_genomepos(start);
    chrhigh = chrpos2 = Pair_genomepos(end);
    sign = +1;

  } else {
    chrhigh = chrpos1 = Pair_genomepos(start);
    chrlow = chrpos2 = Pair_genomepos(end);
    sign = -1;
  }

  if (map_bothstrands_p == true) {
    iit_matches = IIT_get_with_divno(&nmatches,map_iit,divno,chrlow,chrhigh,/*sortp*/false);
    if (nflanking > 0) {
      IIT_get_flanking_with_divno(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,map_iit,
				  divno,chrlow,chrhigh,nflanking,/*sign*/0);
    }
    if (nflanking > 0) {
      FPRINTF(fp,"  Map hits for path %d (%d|%d|%d):\n",pathnum,nleftflanks,nmatches,nrightflanks);
    } else {
      FPRINTF(fp,"  Map hits for path %d (%d):\n",pathnum,nmatches);
    }
    if (nflanking > 0) {
      IIT_print_header(fp,map_iit,leftflanks,nleftflanks,chr,
		       /*reversep*/true,/*relativep*/false,/*left*/0U,print_comment_p);
      FPRINTF(fp,"    ====================\n");
    }
    IIT_print_header(fp,map_iit,iit_matches,nmatches,chr,
		     /*reversep*/false,/*relativep*/false,/*left*/0U,print_comment_p);
    if (nflanking > 0) {
      FPRINTF(fp,"    ====================\n");
      IIT_print_header(fp,map_iit,rightflanks,nrightflanks,chr,
		       /*reversep*/false,/*relativep*/false,/*left*/0U,print_comment_p);
    }

  } else {
    iit_matches = IIT_get_signed_with_divno(&nmatches,map_iit,divno,chrlow,chrhigh,/*sortp*/true,sign);
    if (nflanking > 0) {
      IIT_get_flanking_with_divno(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,map_iit,
				  divno,chrlow,chrhigh,nflanking,sign);
    }
    if (nflanking > 0) {
      FPRINTF(fp,"  Map hits for path %d (%d|%d|%d):\n",pathnum,nleftflanks,nmatches,nrightflanks);
    } else {
      FPRINTF(fp,"  Map hits for path %d (%d):\n",pathnum,nmatches);
    }
    if (nflanking > 0) {
      IIT_print_header(fp,map_iit,leftflanks,nleftflanks,chr,
		       /*reversep*/true,/*relativep*/false,/*left*/0U,print_comment_p);
      FPRINTF(fp,"    ====================\n");
    }
    IIT_print_header(fp,map_iit,iit_matches,nmatches,chr,
		     /*reversep*/false,/*relativep*/false,/*left*/0U,print_comment_p);
    if (nflanking > 0) {
      FPRINTF(fp,"    ====================\n");
      IIT_print_header(fp,map_iit,rightflanks,nrightflanks,chr,
		       /*reversep*/false,/*relativep*/false,/*left*/0U,print_comment_p);
    }
  }
  FPRINTF(fp,"\n");

  if (nflanking > 0) {
    FREE(rightflanks);
    FREE(leftflanks);
  }

  FREE(iit_matches);
  FREE(chr);
  return;
}


/* Doesn't handle nflanking */
static void
print_exon_map (Filestring_T fp, T this, IIT_T map_iit, int *map_divint_crosstable, 
		Univ_IIT_T chromosome_iit, int pathnum, bool map_bothstrands_p, bool print_comment_p) {
  Uintlist_T exonbounds;
  Chrpos_T position1, position2;
  int *iit_matches = NULL, nmatches;
  int divno, exonno = 0;
  char *chr;

  if ((divno = map_divint_crosstable[this->chrnum]) <= 0) {
    FPRINTF(fp,"  *Map hits for path %d (0):\n\n",pathnum);
    return;
  } else {
    chr = Chrnum_to_string(this->chrnum,chromosome_iit);
  }

  exonbounds = Pair_exonbounds(this->pairarray,this->npairs);

  while (exonbounds != NULL) {
    exonbounds = Uintlist_pop(exonbounds,&position1);
    exonbounds = Uintlist_pop(exonbounds,&position2);

    if (map_bothstrands_p == true) {
      if (position1 < position2) {
	iit_matches = IIT_get(&nmatches,map_iit,chr,position1,position2,/*sortp*/true);
      } else {
	iit_matches = IIT_get(&nmatches,map_iit,chr,position2,position1,/*sortp*/true);
      }
      FPRINTF(fp,"  Map hits for path %d, exon %d (%d):\n",pathnum,++exonno,nmatches);
      IIT_print_header(fp,map_iit,iit_matches,nmatches,chr,
		       /*reversep*/false,/*relativep*/false,/*left*/0U,print_comment_p);

    } else {
      if (position1 < position2) {
	iit_matches = IIT_get_signed_with_divno(&nmatches,map_iit,divno,position1,position2,
						/*sortp*/true,/*sign*/+1);
      } else {
	iit_matches = IIT_get_signed_with_divno(&nmatches,map_iit,divno,position2,position1,
						/*sortp*/true,/*sign*/-1);
      }
      FPRINTF(fp,"  Map hits for path %d, exon %d (%d):\n",pathnum,++exonno,nmatches);
      IIT_print_header(fp,map_iit,iit_matches,nmatches,chr,
		       /*reversep*/false,/*relativep*/false,/*left*/0U,print_comment_p);
    }
    FPRINTF(fp,"\n");
    FREE(iit_matches);
  }

  return;
}

void
Stage3_print_map (Filestring_T fp, T this, IIT_T map_iit, int *map_divint_crosstable, Univ_IIT_T chromosome_iit,
		  int pathnum, bool map_exons_p, bool map_bothstrands_p, int nflanking,
		  bool print_comment_p) {
  if (map_exons_p == true) {
    print_exon_map(fp,this,map_iit,map_divint_crosstable,
		   chromosome_iit,pathnum,map_bothstrands_p,print_comment_p);
  } else {
    print_map(fp,this,map_iit,map_divint_crosstable,
	      chromosome_iit,pathnum,map_bothstrands_p,nflanking,print_comment_p);
  }
  return;
}



/* queryaaseq is used only by PMAP */
void
Stage3_print_alignment (Filestring_T fp, T this, Genome_T genome,
			Univ_IIT_T chromosome_iit, Printtype_T printtype,
			bool continuousp, bool continuous_by_exon_p, bool genomefirstp,
			int invertmode, bool nointronlenp, int wraplength) {
  if (continuous_by_exon_p == true) {
    Pair_print_exonsummary(fp,this->pairarray,this->npairs,this->chrnum,this->chroffset,
			   genome,chromosome_iit,this->watsonp,this->cdna_direction,genomefirstp,invertmode);
    Pair_print_continuous_byexon(fp,this->pairarray,this->npairs,this->watsonp,invertmode);

  } else if (continuousp == true) {
    Pair_print_continuous(fp,this->pairarray,this->npairs,this->watsonp,
			  genomefirstp,invertmode,nointronlenp);
  } else {
    /* Assumes Stage3_print_pathsummary already called */
    Pair_print_exonsummary(fp,this->pairarray,this->npairs,this->chrnum,this->chroffset,
			   genome,chromosome_iit,this->watsonp,this->cdna_direction,genomefirstp,invertmode);
    if (printtype == ALIGNMENT) {
      Pair_print_alignment(fp,this->pairarray,this->npairs,this->chrnum,this->chroffset,
			   chromosome_iit,this->watsonp,invertmode,nointronlenp,wraplength);
    }
  }
  debug1(Pair_dump_array(this->pairarray,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairarray,this->npairs));
  return;
}


void
Stage3_print_coordinates (Filestring_T fp, T this, Univ_IIT_T chromosome_iit, int invertmode) {
  Pair_print_coordinates(fp,this->pairarray,this->npairs,this->chrnum,this->chroffset,
			 chromosome_iit,this->watsonp,invertmode);
  return;
}


void
Stage3_print_cdna (Filestring_T fp, T this, int wraplength) {
#ifdef PMAP
  Pair_print_nucleotide_cdna(fp,this->pairarray,this->npairs,wraplength);
#else
  if (this->cdna_direction >= 0) {
    Pair_print_protein_cdna(fp,this->pairarray,this->npairs,wraplength,/*forwardp*/true);
  } else {
    Pair_print_protein_cdna(fp,&(this->pairarray[this->npairs-1]),this->npairs,wraplength,/*forwardp*/false);
  }
#endif
  return;
}

void
Stage3_print_protein_genomic (Filestring_T fp, T this, int wraplength) {
  if (this->cdna_direction >= 0) {
    Pair_print_protein_genomic(fp,this->pairarray,this->npairs,wraplength,/*forwardp*/true);
  } else {
    Pair_print_protein_genomic(fp,&(this->pairarray[this->npairs-1]),this->npairs,wraplength,/*forwardp*/false);
  }
  return;
}


void
Stage3_print_compressed (Filestring_T fp, T this, Sequence_T queryseq, Univ_IIT_T chromosome_iit,
			 char *dbversion, Sequence_T usersegment, int pathnum, int npaths,
			 bool checksump, int chimerapos, int chimeraequivpos,
			 double donor_prob, double acceptor_prob, int chimera_cdna_direction) {
  Pair_T start, end;

#if 0
#ifdef PMAP
  Stage3_translate_cdna(this,queryseq,strictp);
  Stage3_backtranslate_cdna(this);
#else
  if (truncatep == true) {
    truncate_fulllength(this,/*translatep*/true,/*cds_startpos*/-1,
			Sequence_fulllength_given(queryseq),strictp);
  }
#endif
#endif

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);
  Pair_print_compressed(fp,pathnum,npaths,start,end,queryseq,dbversion,usersegment,
			Pair_nexons(this->pairarray,this->npairs),
			Stage3_fracidentity(this),this->pairarray,this->npairs,
			this->chrnum,this->chroffset,chromosome_iit,
			Sequence_fulllength_given(queryseq),Sequence_skiplength(queryseq),
			Sequence_trim_start(queryseq),Sequence_trim_end(queryseq),
			checksump,chimerapos,chimeraequivpos,donor_prob,acceptor_prob,
			chimera_cdna_direction,this->strain,this->cdna_direction);
  return;
}



#if 0
static int
compute_introntype (char left1, char left2, char right2, char right1) {
  int leftdi, rightdi;

  if (left1 == 'G' && left2 == 'T') {
    leftdi = LEFT_GT;
  } else if (left1 == 'G' && left2 == 'C') {
    leftdi = LEFT_GC;
  } else if (left1 == 'A' && left2 == 'T') {
    leftdi = LEFT_AT;
#ifndef PMAP
  } else if (left1 == 'C' && left2 == 'T') {
    leftdi = LEFT_CT;
#endif
  } else {
    leftdi = 0x00;
  }

  if (right2 == 'A' && right1 == 'G') {
    rightdi = RIGHT_AG;
  } else if (right2 == 'A' && right1 == 'C') {
    rightdi = RIGHT_AC;
#ifndef PMAP
  } else if (right2 == 'G' && right1 == 'C') {
    rightdi = RIGHT_GC;
  } else if (right2 == 'A' && right1 == 'T') {
    rightdi = RIGHT_AT;
#endif
  } else {
    rightdi = 0x00;
  }

  return leftdi & rightdi;
}
#endif

#if 0
static char uppercaseCode[128] = UPPERCASE_U2T;
#endif

#if 0
static List_T
peel_leftward_old (bool *mismatchp, List_T *peeled_path, List_T path, int *querydp5, int *genomedp5, 
#ifdef WASTE
		   Pairpool_T pairpool,
#endif
		   int maxpeelback, bool throughmismatchp, bool quit_on_gap_p,
		   List_T *endgappairs, Pair_T *gappair, int *querydp5_medialgap, int *genomedp5_medialgap) {
  List_T peeled = NULL, rest = NULL, pairptr;
  Pair_T pair, nextpair, rightpair;
  int npeelback = 0, nconsecutive = 0, init_dynprogindex = DYNPROGINDEX_MINOR;
  bool stopp;
  int nmatches;
#if 0
  int nincursion = 0;
#endif

  *mismatchp = false;
  debug(printf("Peeling leftward:"));
  if (path == NULL) {
    debug(printf(" path is empty\n"));
  } else {
    pair = path->first;
    if (pair->gapp == true) {
      /* Throw away known gap */
      debug(printf(" Known_gap"));
      pairptr = path;
      path = Pairpool_pop(path,&pair);
#ifdef WASTE
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
#else
      peeled = List_push_existing(peeled,pairptr);
#endif
    }
    rest = path->rest;
    
    stopp = false;
    while (rest != NULL && stopp == false) {
      nextpair = rest->first;
      if (nextpair->gapp == true || nextpair->cdna == ' ' || nextpair->genome == ' ' || nextpair->protectedp == true) {
	stopp = true;
      } else {
	pairptr = path;
	path = Pairpool_pop(path,&pair);
#ifdef WASTE
	peeled = Pairpool_push_existing(peeled,pairpool,pair);
#else
	peeled = List_push_existing(peeled,pairptr);
#endif
	debug(printf(" Peel [");
	      Pair_dump_one(pair,/*zerobasedp*/true);
	      printf("]"));

	if (uppercaseCode[(int) pair->cdna] != uppercaseCode[(int) pair->genome]) {
	  *mismatchp = true;
	}

	if (++npeelback >= maxpeelback) {
	  stopp = true;
	}

	if (init_dynprogindex > 0 && pair->dynprogindex <= 0) {
	  init_dynprogindex = pair->dynprogindex;
	}

	rest = path->rest;
      }
    }

    /* Continue to peelback through little skips and mismatches */
    debug(printf("\n||"));

    stopp = false;
    while (rest != NULL && stopp == false) {
      nextpair = rest->first;
      if (nextpair->gapp == true) {
	/* Peel this one, but then stop at end of loop */
      } else if (nextpair->protectedp == true) {
	/* Stop because it's protected */
	stopp = true;
      } else if (nextpair->cdna != ' ' && nextpair->genome != ' ') {
	/* Stop because it looks okay */
	stopp = true;
      }

      pairptr = path;
      path = Pairpool_pop(path,&pair);
#ifdef WASTE
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
#else
      peeled = List_push_existing(peeled,pairptr);
#endif
      debug(printf(" Extrapeel [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
	
      if (uppercaseCode[(int) pair->cdna] != uppercaseCode[(int) pair->genome]) {
	*mismatchp = true;
      }

#if 0
      if (pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP || pair->comp == MISMATCH_COMP) {
	nconsecutive = 0;
      } else if (++nconsecutive >= SUFFCONSECUTIVE) {
	stopp = true;
      }
#endif
      
#if 0
      if (pair->dynprogindex != init_dynprogindex) {
	if (++nincursion >= MAXINCURSION) {
	  stopp = true;
	}
      }
#endif

      if (pair->gapp == true) {
	stopp = true;
      }

      rest = path->rest;
    }
  }

#ifdef PMAP
  /* Reverse process to codon boundary.  Cases:

     X X X | X -
     0 1 2   3 4

     X X X | - X
     0 1 2   3 3

     X X - X | X
     0 1 2 2   3

     Rule: nextpair->querypos % 3 == 0 */

  debug(printf("\n<<"));
  if (peeled != NULL) {
    rest = peeled->rest;
    stopp = false;
    while (rest != NULL && stopp == false) {
      pairptr = peeled;
      peeled = Pairpool_pop(peeled,&pair);
#ifdef WASTE
      path = Pairpool_push_existing(path,pairpool,pair);
#else
      path = List_push_existing(path,pairptr);
#endif
      debug(printf(" Mod3putback [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
      nextpair = rest->first;
      if (nextpair->querypos % 3 == 0) {
	stopp = true;
      }
      rest = peeled->rest;
    }
  }
#endif

  if (peeled == NULL) {
    /* Do not alter querydp5 or genomedp5 */
  } else {
    rightpair = peeled->first;
    while (peeled != NULL && (rightpair->gapp == true || rightpair->comp == INDEL_COMP || rightpair->comp == SHORTGAP_COMP)) {
      debug(printf("Ran into gap; undoing peel, case 1, rightpair gapp %d, comp %c\n",
		   rightpair->gapp,rightpair->comp));
      if (endgappairs != NULL) {
	path = Pairpool_transfer(path,*endgappairs);
	*endgappairs = (List_T) NULL;
      }

      if (quit_on_gap_p == true) {
	path = Pairpool_transfer(path,peeled);
	*peeled_path = (List_T) NULL;
	return path;

      } else {
	/* Put back 1 */
	/* if ((pairptr = peeled) != NULL) { */
	pairptr = peeled;
	peeled = Pairpool_pop(peeled,&pair);
	path = List_push_existing(path,pairptr);
	debug(printf(" Putback [");
	      Pair_dump_one(pair,/*zerobasedp*/true);
	      printf("]"));
	  /* } */

#if 0
	/* Put back 2 */
	if ((pairptr = peeled) != NULL) {
	  peeled = Pairpool_pop(peeled,&pair);
	  path = List_push_existing(path,pairptr);
	  debug(printf(" Putback [");
		Pair_dump_one(pair,/*zerobasedp*/true);
		printf("]"));
	}
#endif
      }

      rightpair = path->first;
    }

    if (path != NULL) {
      rightpair = path->first;
      *querydp5 = rightpair->querypos + 1;
      *genomedp5 = rightpair->genomepos + 1;
    } else if (peeled != NULL) {
      rightpair = peeled->first;
      *querydp5 = rightpair->querypos;
      *genomedp5 = rightpair->genomepos;
    } else {
      fprintf(stderr,"In peel_rightward, path and peeled are both NULL\n");
      abort();
    }
  }

  if (endgappairs != NULL) {
    if (path == NULL || (pair = path->first) == NULL || (pair->gapp == false && pair->comp != INDEL_COMP && pair->comp != SHORTGAP_COMP)) {
      *endgappairs = NULL;
      *querydp5_medialgap = *querydp5;
      *genomedp5_medialgap = *genomedp5;
    } else {
      pairptr = path;
      path = Pairpool_pop(path,&pair);
#ifdef WASTE
      *endgappairs = Pairpool_push_existing(NULL,pairpool,pair);
#else
      *endgappairs = List_push_existing(NULL,pairptr);
#endif
      debug(printf(" Peeling gap [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
      *gappair = pair;
      debug(printf(" gapcomp: '%c'",pair->comp));

      nmatches = 0;
      while (path != NULL && nmatches < 3) {
	pairptr = path;
	path = Pairpool_pop(path,&pair);
#ifdef WASTE
	*endgappairs = Pairpool_push_existing(*endgappairs,pairpool,pair);
#else
	*endgappairs = List_push_existing(*endgappairs,pairptr);
#endif
	debug(printf(" Peeling after gap [");
	      Pair_dump_one(pair,/*zerobasedp*/true);
	      printf("]"));
	if (uppercaseCode[(int) pair->cdna] == uppercaseCode[(int) pair->genome]) {
	  nmatches++;
	}
      }

      rightpair = (*endgappairs)->first;
      if (rightpair->gapp == true || rightpair->comp == INDEL_COMP || rightpair->comp == SHORTGAP_COMP) {
	debug(printf("Ran into gap; undoing peel, case 2\n"));
	path = Pairpool_transfer(path,*endgappairs);
	*endgappairs = (List_T) NULL;

	if (quit_on_gap_p == true) {
	  path = Pairpool_transfer(path,peeled);
	  *peeled_path = (List_T) NULL;
	  return path;

	} else {
	  /* Put back 1 */
	  if ((pairptr = peeled) != NULL) {
	    peeled = Pairpool_pop(peeled,&pair);
	    path = List_push_existing(path,pairptr);
	    debug(printf(" Putback [");
		  Pair_dump_one(pair,/*zerobasedp*/true);
		  printf("]"));
	  }

	  /* Put back 2 */
	  if ((pairptr = peeled) != NULL) {
	    peeled = Pairpool_pop(peeled,&pair);
	    path = List_push_existing(path,pairptr);
	    debug(printf(" Putback [");
		  Pair_dump_one(pair,/*zerobasedp*/true);
		  printf("]"));
	  }
	}
      }

      if (path != NULL) {
	rightpair = path->first;
	*querydp5_medialgap = rightpair->querypos + 1;
	*genomedp5_medialgap = rightpair->genomepos + 1;
      } else if (peeled != NULL) {
	rightpair = peeled->first;
	*querydp5_medialgap = rightpair->querypos;
	*genomedp5_medialgap = rightpair->genomepos;
      } else {
	fprintf(stderr,"In peel_rightward for medialgap, path and peeled are both NULL\n");
	abort();
      }
    }
  }

  /* assert(peeled == NULL || path == NULL || ((Pair_T) path->first)->comp != INDEL_COMP); */
  debug(
	if (path == NULL) {
	  printf(" => Top of path is NULL.");
	} else {
	  pair = path->first;
	  printf(" => Top of path is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp5 = %d, genomedp5 = %d\n",*querydp5,*genomedp5);
	);

  *peeled_path = peeled;
  return path;
}
#endif


static List_T
peel_leftward (int *n_peeled_indels, bool *protectedp, List_T *peeled_path, List_T path, int *querydp5, Chrpos_T *genomedp5, 
	       int maxpeelback, bool stop_at_indels_p) {
  List_T peeled = NULL;
  Pair_T pair, rightpair;
  int npeelback = 0, niter;
#if 0
  int nincursion = 0;
#endif

  *n_peeled_indels = 0;
  /* *protectedp = false; -- set by calling procedure */

  debug(printf("Peeling leftward with maxpeelback %d and stop_at_indels_p %d:",maxpeelback,stop_at_indels_p));

  /* Remove initial gaps */
  while (path != NULL &&
	 ( ((Pair_T) path->first)->gapp == true || 
	   ((Pair_T) path->first)->comp == INDEL_COMP ||
	   ((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    path = Pairpool_pop(path,&pair);
  }

  if (path == NULL) {
    debug(printf(" path is empty\n"));

  } else if (stop_at_indels_p == true) {
    pair = path->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&path);
    }

    /* Peel initial indels anyway */
    while (path != NULL && ( ((Pair_T) path->first)->comp == INDEL_COMP || ((Pair_T) path->first)->comp == SHORTGAP_COMP )) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
      peeled = List_transfer_one(peeled,&path);
    }

    while (npeelback < maxpeelback && path != NULL &&
	   ((Pair_T) path->first)->gapp == false &&
	   ((Pair_T) path->first)->comp != INDEL_COMP &&
	   ((Pair_T) path->first)->comp != SHORTGAP_COMP) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
      if (((Pair_T) path->first)->protectedp == true) {
	*protectedp = true;
      }
      peeled = List_transfer_one(peeled,&path);
      npeelback++;
    }
	   
  } else {
    /* Don't stop at indels, but do stop at gaps */
    pair = path->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&path);
    }

    niter = 0;
    while (npeelback < maxpeelback && niter < MAXITER && path != NULL &&
	   ((Pair_T) path->first)->gapp == false) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
      if (((Pair_T) path->first)->comp == MATCH_COMP || ((Pair_T) path->first)->comp == DYNPROG_MATCH_COMP || ((Pair_T) path->first)->comp == AMBIGUOUS_COMP) {
	npeelback++;
      } else if (((Pair_T) path->first)->comp == INDEL_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else if (((Pair_T) path->first)->comp == SHORTGAP_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else {
	npeelback--;
      }
      if (((Pair_T) path->first)->protectedp == true) {
	*protectedp = true;
      }
      niter++;
      peeled = List_transfer_one(peeled,&path);
    }
	   
    if (path != NULL && ((Pair_T) path->first)->gapp == true) {
      debug(printf(" Hit gap [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
    }
  }

  if (path != NULL &&
      ( ((Pair_T) path->first)->gapp == true || 
	((Pair_T) path->first)->comp == INDEL_COMP ||
	((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    /* Don't leave a gap or indel on the top of the path */
    while (peeled != NULL &&
	   ( ((Pair_T) peeled->first)->gapp == true ||
	     ((Pair_T) peeled->first)->comp == INDEL_COMP ||
	     ((Pair_T) peeled->first)->comp == SHORTGAP_COMP)) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      path = List_transfer_one(path,&peeled);
    }
    if (peeled != NULL) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      path = List_transfer_one(path,&peeled); /* This should be match or mismatch */
    }
  }

  if (path != NULL) {
    rightpair = path->first;
    *querydp5 = rightpair->querypos + 1;
    *genomedp5 = rightpair->genomepos + 1;
  } else if (peeled != NULL) {
    rightpair = peeled->first;
    *querydp5 = rightpair->querypos;
    *genomedp5 = rightpair->genomepos;
  } else {
    /* fprintf(stderr,"In peel_leftward, path and peeled are both NULL\n"); */
    /* abort(); */
  }

  debug(
	if (path == NULL) {
	  printf(" => Top of path is NULL.");
	} else {
	  pair = path->first;
	  printf(" => Top of path is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp5 = %d, genomedp5 = %d\n",*querydp5,*genomedp5);
	);

  *peeled_path = peeled;
  return path;
}
    

static List_T
peel_leftward_intron (int *n_peeled_indels, bool *protectedp, List_T *peeled_path, List_T path, int *querydp5, Chrpos_T *genomedp5, 
		      Chrpos_T genomedp3, bool stop_at_indels_p, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
		      int minpeelback, int min_mismatches) {
  List_T peeled = NULL;
  Pair_T pair, rightpair;
  int npeelback = 0, nmismatches = 0, niter;
  char cdna, intron_nt, intron_nt_alt;
#if 0
  int nincursion = 0;
#endif
  /* int maxpeelback = 12; */

  *n_peeled_indels = 0;
  /* *protectedp = false; -- set by calling procedure */

  debug(printf("Peeling leftward with genomedp3 %d and stop_at_indels_p %d:",genomedp3,stop_at_indels_p));

  /* Remove initial gaps */
  while (path != NULL &&
	 ( ((Pair_T) path->first)->gapp == true || 
	   ((Pair_T) path->first)->comp == INDEL_COMP ||
	   ((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    path = Pairpool_pop(path,&pair);
  }

  if (path == NULL) {
    debug(printf(" path is empty\n"));

  } else if (stop_at_indels_p == true) {
    pair = path->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&path);
    }

    /* Peel initial indels anyway */
    while (path != NULL && ( ((Pair_T) path->first)->comp == INDEL_COMP || ((Pair_T) path->first)->comp == SHORTGAP_COMP )) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
      peeled = List_transfer_one(peeled,&path);
    }

    while (/*npeelback < maxpeelback &&*/
	   (npeelback < minpeelback || nmismatches < min_mismatches) && path != NULL &&
	   ((Pair_T) path->first)->gapp == false &&
	   ((Pair_T) path->first)->comp != INDEL_COMP &&
	   ((Pair_T) path->first)->comp != SHORTGAP_COMP) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));

      intron_nt = get_genomic_nt(&intron_nt_alt,genomedp3--,chroffset,chrhigh,watsonp);
      if ((cdna = ((Pair_T) path->first)->cdna) != intron_nt && cdna != intron_nt_alt) {
	nmismatches++;
	debug(printf(" (1) Intron mismatch #%d: %c != %c or %c at %u\n",nmismatches,cdna,intron_nt,intron_nt_alt,genomedp3+1));
      }

      if (((Pair_T) path->first)->protectedp == true) {
	*protectedp = true;
      }
      peeled = List_transfer_one(peeled,&path);
      npeelback++;
    }
	   
  } else {
    /* Don't stop at indels, but do stop at gaps */
    pair = path->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&path);
    }

    niter = 0;
    while (/*npeelback < maxpeelback &&*/
	   (npeelback < minpeelback || nmismatches < min_mismatches) && niter < MAXITER && path != NULL &&
	   ((Pair_T) path->first)->gapp == false) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));

      intron_nt = get_genomic_nt(&intron_nt_alt,genomedp3--,chroffset,chrhigh,watsonp);
      if ((cdna = ((Pair_T) path->first)->cdna) != intron_nt && cdna != intron_nt_alt) {
	nmismatches++;
	debug(printf(" (2) Intron mismatch #%d: %c != %c or %c at %u\n",nmismatches,cdna,intron_nt,intron_nt_alt,genomedp3+1));
      }

      if (((Pair_T) path->first)->comp == MATCH_COMP || ((Pair_T) path->first)->comp == DYNPROG_MATCH_COMP || ((Pair_T) path->first)->comp == AMBIGUOUS_COMP) {
	npeelback++;
      } else if (((Pair_T) path->first)->comp == INDEL_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else if (((Pair_T) path->first)->comp == SHORTGAP_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else {
	npeelback--;
      }
      if (((Pair_T) path->first)->protectedp == true) {
	*protectedp = true;
      }
      niter++;
      peeled = List_transfer_one(peeled,&path);
    }
	   
    if (path != NULL && ((Pair_T) path->first)->gapp == true) {
      debug(printf(" Hit gap [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
    }
  }

  if (path != NULL &&
      ( ((Pair_T) path->first)->gapp == true || 
	((Pair_T) path->first)->comp == INDEL_COMP ||
	((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    /* Don't leave a gap or indel on the top of the path */
    while (peeled != NULL &&
	   ( ((Pair_T) peeled->first)->gapp == true ||
	     ((Pair_T) peeled->first)->comp == INDEL_COMP ||
	     ((Pair_T) peeled->first)->comp == SHORTGAP_COMP)) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      path = List_transfer_one(path,&peeled);
    }
    if (peeled != NULL) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      path = List_transfer_one(path,&peeled); /* This should be match or mismatch */
    }
  }

  if (path != NULL) {
    rightpair = path->first;
    *querydp5 = rightpair->querypos + 1;
    *genomedp5 = rightpair->genomepos + 1;
  } else if (peeled != NULL) {
    rightpair = peeled->first;
    *querydp5 = rightpair->querypos;
    *genomedp5 = rightpair->genomepos;
  } else {
    /* fprintf(stderr,"In peel_leftward, path and peeled are both NULL\n"); */
    /* abort(); */
  }

  debug(
	if (path == NULL) {
	  printf(" => Top of path is NULL.");
	} else {
	  pair = path->first;
	  printf(" => Top of path is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp5 = %d, genomedp5 = %d\n",*querydp5,*genomedp5);
	);

  *peeled_path = peeled;
  return path;
}


#if 0
/* Not sure if we need this, or if it causes GMAP to fail on some alignments */
static List_T
peel_leftward_contiguous (int *n_peeled_indels, bool *protectedp, List_T *peeled_path, List_T path, int *querydp5, Chrpos_T *genomedp5, 
			  int maxpeelback, bool stop_at_indels_p) {
  List_T peeled = NULL;
  Pair_T pair, rightpair;
  int npeelback = 0, niter;
#if 0
  int nincursion = 0;
#endif
  int last_querypos;
  Chrpos_T last_genomepos;


  *n_peeled_indels = 0;
  /* *protectedp = false; -- set by calling procedure */

  debug(printf("Peeling leftward with maxpeelback %d and stop_at_indels_p %d:",maxpeelback,stop_at_indels_p));

  /* Remove initial gaps */
  while (path != NULL &&
	 ( ((Pair_T) path->first)->gapp == true || 
	   ((Pair_T) path->first)->comp == INDEL_COMP ||
	   ((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    path = Pairpool_pop(path,&pair);
  }

  if (path == NULL) {
    debug(printf(" path is empty\n"));

  } else if (stop_at_indels_p == true) {
    pair = path->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&path);
    }

    /* Peel initial indels anyway */
    while (path != NULL && ( ((Pair_T) path->first)->comp == INDEL_COMP || ((Pair_T) path->first)->comp == SHORTGAP_COMP )) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
      peeled = List_transfer_one(peeled,&path);
    }

    if (path != NULL) {
      last_querypos = ((Pair_T) path->first)->querypos;
      last_genomepos = ((Pair_T) path->first)->genomepos;
    }
    while (npeelback < maxpeelback && path != NULL &&
	   ((Pair_T) path->first)->gapp == false &&
	   ((Pair_T) path->first)->comp != INDEL_COMP &&
	   ((Pair_T) path->first)->comp != SHORTGAP_COMP &&
	   ((Pair_T) path->first)->querypos + 1 >= last_querypos &&
	   ((Pair_T) path->first)->genomepos + 1 >= last_genomepos) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
      if (((Pair_T) path->first)->protectedp == true) {
	*protectedp = true;
      }
      last_querypos = ((Pair_T) path->first)->querypos;
      last_genomepos = ((Pair_T) path->first)->genomepos;
      peeled = List_transfer_one(peeled,&path);
      npeelback++;
    }
	   
  } else {
    /* Don't stop at indels, but do stop at gaps */
    pair = path->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&path);
    }

    niter = 0;
    if (path != NULL) {
      last_querypos = ((Pair_T) path->first)->querypos;
      last_genomepos = ((Pair_T) path->first)->genomepos;
    }
    while (npeelback < maxpeelback && niter < MAXITER && path != NULL &&
	   ((Pair_T) path->first)->gapp == false &&
	   ((Pair_T) path->first)->querypos + 1 >= last_querypos &&
	   ((Pair_T) path->first)->genomepos + 1 >= last_genomepos) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
      if (((Pair_T) path->first)->comp == MATCH_COMP || ((Pair_T) path->first)->comp == DYNPROG_MATCH_COMP || ((Pair_T) path->first)->comp == AMBIGUOUS_COMP) {
	npeelback++;
      } else if (((Pair_T) path->first)->comp == INDEL_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else if (((Pair_T) path->first)->comp == SHORTGAP_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else {
	npeelback--;
      }
      if (((Pair_T) path->first)->protectedp == true) {
	*protectedp = true;
      }
      niter++;
      last_querypos = ((Pair_T) path->first)->querypos;
      last_genomepos = ((Pair_T) path->first)->genomepos;
      peeled = List_transfer_one(peeled,&path);
    }
	   
    if (path != NULL && ((Pair_T) path->first)->gapp == true) {
      debug(printf(" Hit gap [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
    }
  }

  if (path != NULL &&
      ( ((Pair_T) path->first)->gapp == true || 
	((Pair_T) path->first)->comp == INDEL_COMP ||
	((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    /* Don't leave a gap or indel on the top of the path */
    while (peeled != NULL &&
	   ( ((Pair_T) peeled->first)->gapp == true ||
	     ((Pair_T) peeled->first)->comp == INDEL_COMP ||
	     ((Pair_T) peeled->first)->comp == SHORTGAP_COMP)) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      path = List_transfer_one(path,&peeled);
    }
    if (peeled != NULL) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      path = List_transfer_one(path,&peeled); /* This should be match or mismatch */
    }
  }

  if (path != NULL) {
    rightpair = path->first;
    *querydp5 = rightpair->querypos + 1;
    *genomedp5 = rightpair->genomepos + 1;
  } else if (peeled != NULL) {
    rightpair = peeled->first;
    *querydp5 = rightpair->querypos;
    *genomedp5 = rightpair->genomepos;
  } else {
    /* fprintf(stderr,"In peel_leftward, path and peeled are both NULL\n"); */
    /* abort(); */
  }

  debug(
	if (path == NULL) {
	  printf(" => Top of path is NULL.");
	} else {
	  pair = path->first;
	  printf(" => Top of path is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp5 = %d, genomedp5 = %d\n",*querydp5,*genomedp5);
	);

  *peeled_path = peeled;
  return path;
}
#endif
    

#if 0
/* Not sure if we need this, or if it causes GMAP to fail on some alignments */
static List_T
peel_leftward_intron_contiguous (int *n_peeled_indels, bool *protectedp, List_T *peeled_path, List_T path, int *querydp5, Chrpos_T *genomedp5, 
				 Chrpos_T genomedp3, bool stop_at_indels_p, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
				 int minpeelback, int min_mismatches) {
  List_T peeled = NULL;
  Pair_T pair, rightpair;
  int npeelback = 0, nmismatches = 0, niter;
  char cdna, intron_nt, intron_nt_alt;
#if 0
  int nincursion = 0;
#endif
  int maxpeelback = 12;
  int last_querypos;
  Chrpos_T last_genomepos;


  *n_peeled_indels = 0;
  /* *protectedp = false; -- set by calling procedure */

  debug(printf("Peeling leftward with genomedp3 %d and stop_at_indels_p %d:",genomedp3,stop_at_indels_p));

  /* Remove initial gaps */
  while (path != NULL &&
	 ( ((Pair_T) path->first)->gapp == true || 
	   ((Pair_T) path->first)->comp == INDEL_COMP ||
	   ((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    path = Pairpool_pop(path,&pair);
  }

  if (path == NULL) {
    debug(printf(" path is empty\n"));

  } else if (stop_at_indels_p == true) {
    pair = path->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&path);
    }

    /* Peel initial indels anyway */
    while (path != NULL && ( ((Pair_T) path->first)->comp == INDEL_COMP || ((Pair_T) path->first)->comp == SHORTGAP_COMP )) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
      peeled = List_transfer_one(peeled,&path);
    }

    if (path != NULL) {
      last_querypos = ((Pair_T) path->first)->querypos;
      last_genomepos = ((Pair_T) path->first)->genomepos;
    }
    while (/*npeelback < maxpeelback &&*/
	   (npeelback < minpeelback || nmismatches < min_mismatches) && path != NULL &&
	   ((Pair_T) path->first)->gapp == false &&
	   ((Pair_T) path->first)->comp != INDEL_COMP &&
	   ((Pair_T) path->first)->comp != SHORTGAP_COMP &&
	   ((Pair_T) path->first)->querypos + 1 >= last_querypos &&
	   ((Pair_T) path->first)->genomepos + 1 >= last_genomepos) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));

      intron_nt = get_genomic_nt(&intron_nt_alt,genomedp3--,chroffset,chrhigh,watsonp);
      if ((cdna = ((Pair_T) path->first)->cdna) != intron_nt && cdna != intron_nt_alt) {
	nmismatches++;
	debug(printf(" (1) Intron mismatch #%d: %c != %c or %c at %u\n",nmismatches,cdna,intron_nt,intron_nt_alt,genomedp3+1));
      }

      if (((Pair_T) path->first)->protectedp == true) {
	*protectedp = true;
      }
      last_querypos = ((Pair_T) path->first)->querypos;
      last_genomepos = ((Pair_T) path->first)->genomepos;
      peeled = List_transfer_one(peeled,&path);
      npeelback++;
    }
	   
  } else {
    /* Don't stop at indels, but do stop at gaps */
    pair = path->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&path);
    }

    niter = 0;
    if (path != NULL) {
      last_querypos = ((Pair_T) path->first)->querypos;
      last_genomepos = ((Pair_T) path->first)->genomepos;
    }
    while (/*npeelback < maxpeelback &&*/
	   (npeelback < minpeelback || nmismatches < min_mismatches) && niter < MAXITER && path != NULL &&
	   ((Pair_T) path->first)->gapp == false &&
	   ((Pair_T) path->first)->querypos + 1 >= last_querypos &&
	   ((Pair_T) path->first)->genomepos + 1 >= last_genomepos) {
      debug(printf(" Peel [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));

      intron_nt = get_genomic_nt(&intron_nt_alt,genomedp3--,chroffset,chrhigh,watsonp);
      if ((cdna = ((Pair_T) path->first)->cdna) != intron_nt && cdna != intron_nt_alt) {
	nmismatches++;
	debug(printf(" (2) Intron mismatch #%d: %c != %c or %c at %u\n",nmismatches,cdna,intron_nt,intron_nt_alt,genomedp3+1));
      }

      if (((Pair_T) path->first)->comp == MATCH_COMP || ((Pair_T) path->first)->comp == DYNPROG_MATCH_COMP || ((Pair_T) path->first)->comp == AMBIGUOUS_COMP) {
	npeelback++;
      } else if (((Pair_T) path->first)->comp == INDEL_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else if (((Pair_T) path->first)->comp == SHORTGAP_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else {
	npeelback--;
      }
      if (((Pair_T) path->first)->protectedp == true) {
	*protectedp = true;
      }
      niter++;
      last_querypos = ((Pair_T) path->first)->querypos;
      last_genomepos = ((Pair_T) path->first)->genomepos;
      peeled = List_transfer_one(peeled,&path);
    }
	   
    if (path != NULL && ((Pair_T) path->first)->gapp == true) {
      debug(printf(" Hit gap [");
	    Pair_dump_one(path->first,/*zerobasedp*/true);
	    printf("]"));
    }
  }

  if (path != NULL &&
      ( ((Pair_T) path->first)->gapp == true || 
	((Pair_T) path->first)->comp == INDEL_COMP ||
	((Pair_T) path->first)->comp == SHORTGAP_COMP)) {
    /* Don't leave a gap or indel on the top of the path */
    while (peeled != NULL &&
	   ( ((Pair_T) peeled->first)->gapp == true ||
	     ((Pair_T) peeled->first)->comp == INDEL_COMP ||
	     ((Pair_T) peeled->first)->comp == SHORTGAP_COMP)) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      path = List_transfer_one(path,&peeled);
    }
    if (peeled != NULL) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      path = List_transfer_one(path,&peeled); /* This should be match or mismatch */
    }
  }

  if (path != NULL) {
    rightpair = path->first;
    *querydp5 = rightpair->querypos + 1;
    *genomedp5 = rightpair->genomepos + 1;
  } else if (peeled != NULL) {
    rightpair = peeled->first;
    *querydp5 = rightpair->querypos;
    *genomedp5 = rightpair->genomepos;
  } else {
    /* fprintf(stderr,"In peel_leftward, path and peeled are both NULL\n"); */
    /* abort(); */
  }

  debug(
	if (path == NULL) {
	  printf(" => Top of path is NULL.");
	} else {
	  pair = path->first;
	  printf(" => Top of path is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp5 = %d, genomedp5 = %d\n",*querydp5,*genomedp5);
	);

  *peeled_path = peeled;
  return path;
}
#endif


#if 0
static List_T
peel_rightward_old (bool *mismatchp, List_T *peeled_pairs, List_T pairs, int *querydp3, int *genomedp3, 
#ifdef WASTE
		    Pairpool_T pairpool,
#endif
		    int maxpeelback, bool throughmismatchp, bool quit_on_gap_p,
		    List_T *endgappairs, Pair_T *gappair, int *querydp3_medialgap, int *genomedp3_medialgap) {
  List_T peeled = NULL, rest = NULL, pairptr;
  Pair_T pair, nextpair, leftpair;
  int npeelback = 0, nconsecutive = 0, init_dynprogindex = DYNPROGINDEX_MINOR;
  bool stopp;
  int nmatches;
#if 0
  int incursion = 0;
#endif

  *mismatchp = false;
  debug(printf("Peeling rightward:"));
  if (pairs == NULL) {
    debug(printf(" pairs is empty\n"));
  } else {
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Throw away known gap */
      debug(printf(" Known_gap"));
      pairptr = pairs;
      pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
#else
      peeled = List_push_existing(peeled,pairptr);
#endif
    }
    rest = pairs->rest;

    stopp = false;
    while (rest != NULL && stopp == false) {
      nextpair = rest->first;
      if (nextpair->gapp == true || nextpair->cdna == ' ' || nextpair->genome == ' ' || nextpair->protectedp == true) {
	stopp = true;
      } else {
	pairptr = pairs;
	pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
	peeled = Pairpool_push_existing(peeled,pairpool,pair);
#else
	peeled = List_push_existing(peeled,pairptr);
#endif
	debug(printf(" Peel [");
	      Pair_dump_one(pair,/*zerobasedp*/true);
	      printf("]"));
      
	if (uppercaseCode[(int) pair->cdna] != uppercaseCode[(int) pair->genome]) {
	  *mismatchp = true;
	}

	if (++npeelback >= maxpeelback) {
	  stopp = true;
	}

	if (init_dynprogindex > 0 && pair->dynprogindex <= 0) {
	  init_dynprogindex = pair->dynprogindex;
	}

	rest = pairs->rest;
      }
    }

    /* Continue to peelback through little skips and mismatches */
    debug(printf("\n||"));

    stopp = false;
    while (rest != NULL && stopp == false) {
      nextpair = rest->first;
      if (nextpair->gapp == true) {
	/* Peel this one, but then stop at end of loop */
      } else if (nextpair->protectedp == true) {
	/* Stop because it's protected */
	stopp = true;
      } else if (nextpair->cdna != ' ' && nextpair->genome != ' ') {
	/* Stop because it looks okay */
	stopp = true;
      }

      pairptr = pairs;
      pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
#else
      peeled = List_push_existing(peeled,pairptr);
#endif
      debug(printf(" Extrapeel [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
	
      if (uppercaseCode[(int) pair->cdna] != uppercaseCode[(int) pair->genome]) {
	*mismatchp = true;
      }

#if 0
      if (pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP || pair->comp == MISMATCH_COMP) {
	nconsecutive = 0;
      } else if (++nconsecutive >= SUFFCONSECUTIVE) {
	stopp = true;
      }
#endif
	
#if 0
      if (pair->dynprogindex != init_dynprogindex) {
	if (++nincursion >= MAXINCURSION) {
	  stopp = true;
	}
      }
#endif

      if (pair->gapp == true) {
	stopp = true;
      }

      rest = pairs->rest;
    }
  }

#ifdef PMAP
  /* Reverse process to codon boundary.  Cases:

     - X | X X X
     5 5   6 7 8

     X - | X X X
     5 6   6 7 8

     X | X - X X
     5   6 7 7 8

     Rule: pair->querypos % 3 == 0 */

  debug(printf("\n<<"));
  stopp = false;
  while (peeled != NULL && stopp == false) {
    pairptr = peeled;
    peeled = Pairpool_pop(peeled,&pair);
#ifdef WASTE
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
    pairs = List_push_existing(pairs,pairptr);
#endif
    debug(printf(" Mod3putback [");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	  printf("]"));
    if (pair->querypos % 3 == 0) {
      stopp = true;
    }
  }
#endif

  if (peeled == NULL) {
    /* Do not alter querydp3 or genomedp3 */
  } else {
    leftpair = peeled->first;
    while (peeled != NULL && (leftpair->gapp == true || leftpair->comp == INDEL_COMP || leftpair->comp == SHORTGAP_COMP)) {
      debug(printf("Ran into gap; undoing peel, case 3, leftpair gapp %d, comp %c\n",
		   leftpair->gapp,leftpair->comp));
      if (endgappairs != NULL) {
	pairs = Pairpool_transfer(pairs,*endgappairs);
	*endgappairs = (List_T) NULL;
      }

      if (quit_on_gap_p == true) {
	pairs = Pairpool_transfer(pairs,peeled);
	*peeled_pairs = (List_T) NULL;
	return pairs;

      } else {
	/* Put back 1 */
	/* if ((pairptr = peeled) != NULL) { */
	pairptr = peeled;
	peeled = Pairpool_pop(peeled,&pair);
	pairs = List_push_existing(pairs,pairptr);
	debug(printf(" Putback [");
	      Pair_dump_one(pair,/*zerobasedp*/true);
	      printf("]"));
	  /* } */

#if 0
	/* Put back 2 */
	if ((pairptr = peeled) != NULL) {
	  peeled = Pairpool_pop(peeled,&pair);
	  pairs = List_push_existing(pairs,pairptr);
	  debug(printf(" Putback [");
		Pair_dump_one(pair,/*zerobasedp*/true);
		printf("]"));
	}
#endif
      }

      leftpair = pairs->first;
    }

    if (pairs != NULL) {
      leftpair = pairs->first;
      *querydp3 = leftpair->querypos - 1;
      *genomedp3 = leftpair->genomepos - 1;
    } else if (peeled != NULL) {
      leftpair = peeled->first;
      *querydp3 = leftpair->querypos;
      *genomedp3 = leftpair->genomepos;
    } else {
      fprintf(stderr,"In peel_leftward, pairs and peeled are both NULL\n");
      abort();
    }
  }

  if (endgappairs != NULL) {
    if (pairs == NULL || (pair = pairs->first) == NULL || (pair->gapp == false && pair->comp != INDEL_COMP && pair->comp != SHORTGAP_COMP)) {
      *endgappairs = NULL;
      *querydp3_medialgap = *querydp3;
      *genomedp3_medialgap = *genomedp3;
    } else {
      pairptr = pairs;
      pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
      *endgappairs = Pairpool_push_existing(NULL,pairpool,pair);
#else
      *endgappairs = List_push_existing(NULL,pairptr);
#endif
      debug(printf(" Peeling gap [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
      *gappair = pair;
      debug(printf(" gapcomp: '%c'",pair->comp));

      nmatches = 0;
      while (pairs != NULL && nmatches < 3) {
	pairptr = pairs;
	pairs = Pairpool_pop(pairs,&pair);
#ifdef WASTE
	*endgappairs = Pairpool_push_existing(*endgappairs,pairpool,pair);
#else
	*endgappairs = List_push_existing(*endgappairs,pairptr);
#endif
	debug(printf(" Peeling after gap [");
	      Pair_dump_one(pair,/*zerobasedp*/true);
	      printf("]"));
	if (uppercaseCode[(int) pair->cdna] == uppercaseCode[(int) pair->genome]) {
	  nmatches++;
	}
      }

      leftpair = (*endgappairs)->first;
      if (leftpair->gapp == true || leftpair->comp == INDEL_COMP || leftpair->comp == SHORTGAP_COMP) {
	debug(printf("Ran into gap; undoing peel, case 4\n"));
	pairs = Pairpool_transfer(pairs,*endgappairs);
	*endgappairs = (List_T) NULL;

	if (quit_on_gap_p == true) {
	  pairs = Pairpool_transfer(pairs,peeled);
	  *peeled_pairs = (List_T) NULL;
	  return pairs;

	} else {
	  /* Put back 1 */
	  if ((pairptr = peeled) != NULL) {
	    peeled = Pairpool_pop(peeled,&pair);
	    pairs = List_push_existing(pairs,pairptr);
	    debug(printf(" Putback [");
		  Pair_dump_one(pair,/*zerobasedp*/true);
		  printf("]"));
	  }

	  /* Put back 2 */
	  if ((pairptr = peeled) != NULL) {
	    peeled = Pairpool_pop(peeled,&pair);
	    pairs = List_push_existing(pairs,pairptr);
	    debug(printf(" Putback [");
		  Pair_dump_one(pair,/*zerobasedp*/true);
		  printf("]"));
	  }
	}
      }

      if (pairs != NULL) {
	leftpair = pairs->first;
	*querydp3_medialgap = leftpair->querypos - 1;
	*genomedp3_medialgap = leftpair->genomepos - 1;
      } else if (peeled != NULL) {
	leftpair = peeled->first;
	*querydp3_medialgap = leftpair->querypos;
	*genomedp3_medialgap = leftpair->genomepos;
      } else {
	fprintf(stderr,"In peel_leftward for medialgap, pairs and peeled are both NULL\n");
	abort();
      }
    }
  }

  /* assert(peeled == NULL || pairs == NULL || ((Pair_T) pairs->first)->comp != INDEL_COMP); */
  debug(
	if (pairs == NULL) {
	  printf(" => Top of pairs is NULL.");
	} else {
	  pair = pairs->first;
	  printf(" => Top of pairs is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp3 = %d, genomedp3 = %d\n",*querydp3,*genomedp3);
	);

  *peeled_pairs = peeled;
  return pairs;
}
#endif


static List_T
peel_rightward (int *n_peeled_indels, bool *protectedp, List_T *peeled_pairs, List_T pairs, int *querydp3, Chrpos_T *genomedp3, 
		int maxpeelback, bool stop_at_indels_p) {
  List_T peeled = NULL;
  Pair_T pair, leftpair;
  int npeelback = 0, niter;
#if 0
  int incursion = 0;
#endif

  *n_peeled_indels = 0;
  /* *protectedp = false; -- set by calling procedure */

  debug(printf("Peeling rightward with maxpeelback %d and stop_at_indels_p %d:",maxpeelback,stop_at_indels_p));

  /* Remove initial gaps */
  while (pairs != NULL && 
	 ( ((Pair_T) pairs->first)->gapp == true ||
	   ((Pair_T) pairs->first)->comp == INDEL_COMP ||
	   ((Pair_T) pairs->first)->comp == SHORTGAP_COMP )) {
    pairs = Pairpool_pop(pairs,&pair);
  }

  if (pairs == NULL) {
    debug(printf(" pairs is empty\n"));

  } else if (stop_at_indels_p == true) {
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    /* Peel initial indels anyway */
    while (pairs != NULL && ( ((Pair_T) pairs->first)->comp == INDEL_COMP || ((Pair_T) pairs->first)->comp == INDEL_COMP )) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    while (npeelback < maxpeelback && pairs != NULL &&
	   ((Pair_T) pairs->first)->gapp == false &&
	   ((Pair_T) pairs->first)->comp != INDEL_COMP &&
	   ((Pair_T) pairs->first)->comp != SHORTGAP_COMP) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
      if (((Pair_T) pairs->first)->protectedp == true) {
	*protectedp = true;
      }
      peeled = List_transfer_one(peeled,&pairs);
      npeelback++;
    }

  } else {
    /* Don't stop at indels, but do stop at gaps */
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    niter = 0;
    while (npeelback < maxpeelback && niter < MAXITER && pairs != NULL &&
	   ((Pair_T) pairs->first)->gapp == false) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
      if (((Pair_T) pairs->first)->comp == MATCH_COMP || ((Pair_T) pairs->first)->comp == DYNPROG_MATCH_COMP || ((Pair_T) pairs->first)->comp == AMBIGUOUS_COMP) {
	npeelback++;
      } else if (((Pair_T) pairs->first)->comp == INDEL_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else if (((Pair_T) pairs->first)->comp == SHORTGAP_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else {
	npeelback--;
      }
      if (((Pair_T) pairs->first)->protectedp == true) {
	*protectedp = true;
      }
      niter++;
      peeled = List_transfer_one(peeled,&pairs);
    }

    if (pairs != NULL && ((Pair_T) pairs->first)->gapp == true) {
      debug(printf(" Hit gap [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
    }
  }
	   
  if (pairs != NULL &&
      ( ((Pair_T) pairs->first)->gapp == true || 
	((Pair_T) pairs->first)->comp == INDEL_COMP ||
	((Pair_T) pairs->first)->comp == SHORTGAP_COMP )) {
    /* Don't leave a gap or indel on the top of the pairs */
    while (peeled != NULL &&
	   ( ((Pair_T) peeled->first)->gapp == true ||
	     ((Pair_T) peeled->first)->comp == INDEL_COMP ||
	     ((Pair_T) peeled->first)->comp == SHORTGAP_COMP)) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      pairs = List_transfer_one(pairs,&peeled);
    }
    if (peeled != NULL) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      pairs = List_transfer_one(pairs,&peeled); /* This should be match or mismatch */
    }
  }

  if (pairs != NULL) {
    leftpair = pairs->first;
    *querydp3 = leftpair->querypos - 1;
    *genomedp3 = leftpair->genomepos - 1;
  } else if (peeled != NULL) {
    leftpair = peeled->first;
    *querydp3 = leftpair->querypos;
    *genomedp3 = leftpair->genomepos;
  } else {
    /* fprintf(stderr,"In peel_rightward, pairs and peeled are both NULL\n"); */
    /* abort(); */
  }

  debug(
	if (pairs == NULL) {
	  printf(" => Top of pairs is NULL.");
	} else {
	  pair = pairs->first;
	  printf(" => Top of pairs is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp3 = %d, genomedp3 = %d\n",*querydp3,*genomedp3);
	);

  *peeled_pairs = peeled;
  return pairs;
}


/* Instead of maxpeelback, follow the 5' intron until we get enough mismatches */
static List_T
peel_rightward_intron (int *n_peeled_indels, bool *protectedp, List_T *peeled_pairs, List_T pairs, int *querydp3, Chrpos_T *genomedp3, 
		       Chrpos_T genomedp5, bool stop_at_indels_p, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
		       int minpeelback, int min_mismatches) {
  List_T peeled = NULL;
  Pair_T pair, leftpair;
  int npeelback = 0, nmismatches = 0, niter;
  char cdna, intron_nt, intron_nt_alt;
#if 0
  int incursion = 0;
#endif
  /* int maxpeelback = 12; */

  *n_peeled_indels = 0;
  /* *protectedp = false; -- set by calling procedure */

  debug(printf("Peeling rightward with genomedp5 %d and stop_at_indels_p %d:",genomedp5,stop_at_indels_p));

  /* Remove initial gaps */
  while (pairs != NULL && 
	 ( ((Pair_T) pairs->first)->gapp == true ||
	   ((Pair_T) pairs->first)->comp == INDEL_COMP ||
	   ((Pair_T) pairs->first)->comp == SHORTGAP_COMP )) {
    pairs = Pairpool_pop(pairs,&pair);
  }

  if (pairs == NULL) {
    debug(printf(" pairs is empty\n"));

  } else if (stop_at_indels_p == true) {
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    /* Peel initial indels anyway */
    while (pairs != NULL && ( ((Pair_T) pairs->first)->comp == INDEL_COMP || ((Pair_T) pairs->first)->comp == INDEL_COMP )) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    while (/*npeelback < maxpeelback &&*/
	   (npeelback < minpeelback || nmismatches < min_mismatches) && pairs != NULL &&
	   ((Pair_T) pairs->first)->gapp == false &&
	   ((Pair_T) pairs->first)->comp != INDEL_COMP &&
	   ((Pair_T) pairs->first)->comp != SHORTGAP_COMP) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));

      intron_nt = get_genomic_nt(&intron_nt_alt,genomedp5++,chroffset,chrhigh,watsonp);
      if ((cdna = ((Pair_T) pairs->first)->cdna) != intron_nt && cdna != intron_nt_alt) {
	nmismatches++;
	debug(printf(" (3) Intron mismatch #%d: %c != %c or %c at %u\n",nmismatches,cdna,intron_nt,intron_nt_alt,genomedp5-1));
      }

      if (((Pair_T) pairs->first)->protectedp == true) {
	*protectedp = true;
      }
      peeled = List_transfer_one(peeled,&pairs);
      npeelback++;
    }

  } else {
    /* Don't stop at indels, but do stop at gaps */
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    niter = 0;
    while (/*npeelback < maxpeelback &&*/
	   (npeelback < minpeelback || nmismatches < min_mismatches) && niter < MAXITER && pairs != NULL &&
	   ((Pair_T) pairs->first)->gapp == false) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));

      intron_nt = get_genomic_nt(&intron_nt_alt,genomedp5++,chroffset,chrhigh,watsonp);
      if ((cdna = ((Pair_T) pairs->first)->cdna) != intron_nt && cdna != intron_nt_alt) {
	nmismatches++;
	debug(printf(" (4) Intron mismatch #%d: %c != %c or %c at %u\n",nmismatches,cdna,intron_nt,intron_nt_alt,genomedp5-1));
      }

      if (((Pair_T) pairs->first)->comp == MATCH_COMP || ((Pair_T) pairs->first)->comp == DYNPROG_MATCH_COMP || ((Pair_T) pairs->first)->comp == AMBIGUOUS_COMP) {
	npeelback++;
      } else if (((Pair_T) pairs->first)->comp == INDEL_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else if (((Pair_T) pairs->first)->comp == SHORTGAP_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else {
	npeelback--;
      }
      if (((Pair_T) pairs->first)->protectedp == true) {
	*protectedp = true;
      }
      niter++;
      peeled = List_transfer_one(peeled,&pairs);
    }

    if (pairs != NULL && ((Pair_T) pairs->first)->gapp == true) {
      debug(printf(" Hit gap [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
    }
  }
	   
  if (pairs != NULL &&
      ( ((Pair_T) pairs->first)->gapp == true || 
	((Pair_T) pairs->first)->comp == INDEL_COMP ||
	((Pair_T) pairs->first)->comp == SHORTGAP_COMP )) {
    /* Don't leave a gap or indel on the top of the pairs */
    while (peeled != NULL &&
	   ( ((Pair_T) peeled->first)->gapp == true ||
	     ((Pair_T) peeled->first)->comp == INDEL_COMP ||
	     ((Pair_T) peeled->first)->comp == SHORTGAP_COMP)) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      pairs = List_transfer_one(pairs,&peeled);
    }
    if (peeled != NULL) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      pairs = List_transfer_one(pairs,&peeled); /* This should be match or mismatch */
    }
  }

  if (pairs != NULL) {
    leftpair = pairs->first;
    *querydp3 = leftpair->querypos - 1;
    *genomedp3 = leftpair->genomepos - 1;
  } else if (peeled != NULL) {
    leftpair = peeled->first;
    *querydp3 = leftpair->querypos;
    *genomedp3 = leftpair->genomepos;
  } else {
    /* fprintf(stderr,"In peel_rightward, pairs and peeled are both NULL\n"); */
    /* abort(); */
  }

  debug(
	if (pairs == NULL) {
	  printf(" => Top of pairs is NULL.");
	} else {
	  pair = pairs->first;
	  printf(" => Top of pairs is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp3 = %d, genomedp3 = %d\n",*querydp3,*genomedp3);
	);

  *peeled_pairs = peeled;
  return pairs;
}


#if 0
/* Not sure if we need this, or if it causes GMAP to fail on some alignments */
static List_T
peel_rightward_contiguous (int *n_peeled_indels, bool *protectedp, List_T *peeled_pairs, List_T pairs, int *querydp3, Chrpos_T *genomedp3, 
			   int maxpeelback, bool stop_at_indels_p) {
  List_T peeled = NULL;
  Pair_T pair, leftpair;
  int npeelback = 0, niter;
#if 0
  int incursion = 0;
#endif
  int last_querypos;
  Chrpos_T last_genomepos;


  *n_peeled_indels = 0;
  /* *protectedp = false; -- set by calling procedure */

  debug(printf("Peeling rightward with maxpeelback %d and stop_at_indels_p %d:",maxpeelback,stop_at_indels_p));

  /* Remove initial gaps */
  while (pairs != NULL && 
	 ( ((Pair_T) pairs->first)->gapp == true ||
	   ((Pair_T) pairs->first)->comp == INDEL_COMP ||
	   ((Pair_T) pairs->first)->comp == SHORTGAP_COMP )) {
    pairs = Pairpool_pop(pairs,&pair);
  }

  if (pairs == NULL) {
    debug(printf(" pairs is empty\n"));

  } else if (stop_at_indels_p == true) {
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    /* Peel initial indels anyway */
    while (pairs != NULL && ( ((Pair_T) pairs->first)->comp == INDEL_COMP || ((Pair_T) pairs->first)->comp == INDEL_COMP )) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    if (pairs != NULL) {
      last_querypos = ((Pair_T) pairs->first)->querypos;
      last_genomepos = ((Pair_T) pairs->first)->genomepos;
    }
    while (npeelback < maxpeelback && pairs != NULL &&
	   ((Pair_T) pairs->first)->gapp == false &&
	   ((Pair_T) pairs->first)->comp != INDEL_COMP &&
	   ((Pair_T) pairs->first)->comp != SHORTGAP_COMP &&
	   ((Pair_T) pairs->first)->querypos <= last_querypos + 1 &&
	   ((Pair_T) pairs->first)->genomepos <= last_genomepos + 1) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
      if (((Pair_T) pairs->first)->protectedp == true) {
	*protectedp = true;
      }
      last_querypos = ((Pair_T) pairs->first)->querypos;
      last_genomepos = ((Pair_T) pairs->first)->genomepos;
      peeled = List_transfer_one(peeled,&pairs);
      npeelback++;
    }

  } else {
    /* Don't stop at indels, but do stop at gaps */
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    niter = 0;
    if (pairs != NULL) {
      last_querypos = ((Pair_T) pairs->first)->querypos;
      last_genomepos = ((Pair_T) pairs->first)->genomepos;
    }
    while (npeelback < maxpeelback && niter < MAXITER && pairs != NULL &&
	   ((Pair_T) pairs->first)->gapp == false &&
	   ((Pair_T) pairs->first)->querypos <= last_querypos + 1 &&
	   ((Pair_T) pairs->first)->genomepos <= last_genomepos + 1) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
      if (((Pair_T) pairs->first)->comp == MATCH_COMP || ((Pair_T) pairs->first)->comp == DYNPROG_MATCH_COMP || ((Pair_T) pairs->first)->comp == AMBIGUOUS_COMP) {
	npeelback++;
      } else if (((Pair_T) pairs->first)->comp == INDEL_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else if (((Pair_T) pairs->first)->comp == SHORTGAP_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else {
	npeelback--;
      }
      if (((Pair_T) pairs->first)->protectedp == true) {
	*protectedp = true;
      }
      niter++;
      last_querypos = ((Pair_T) pairs->first)->querypos;
      last_genomepos = ((Pair_T) pairs->first)->genomepos;
      peeled = List_transfer_one(peeled,&pairs);
    }

    if (pairs != NULL && ((Pair_T) pairs->first)->gapp == true) {
      debug(printf(" Hit gap [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
    }
  }
	   
  if (pairs != NULL &&
      ( ((Pair_T) pairs->first)->gapp == true || 
	((Pair_T) pairs->first)->comp == INDEL_COMP ||
	((Pair_T) pairs->first)->comp == SHORTGAP_COMP )) {
    /* Don't leave a gap or indel on the top of the pairs */
    while (peeled != NULL &&
	   ( ((Pair_T) peeled->first)->gapp == true ||
	     ((Pair_T) peeled->first)->comp == INDEL_COMP ||
	     ((Pair_T) peeled->first)->comp == SHORTGAP_COMP)) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      pairs = List_transfer_one(pairs,&peeled);
    }
    if (peeled != NULL) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      pairs = List_transfer_one(pairs,&peeled); /* This should be match or mismatch */
    }
  }

  if (pairs != NULL) {
    leftpair = pairs->first;
    *querydp3 = leftpair->querypos - 1;
    *genomedp3 = leftpair->genomepos - 1;
  } else if (peeled != NULL) {
    leftpair = peeled->first;
    *querydp3 = leftpair->querypos;
    *genomedp3 = leftpair->genomepos;
  } else {
    /* fprintf(stderr,"In peel_rightward, pairs and peeled are both NULL\n"); */
    /* abort(); */
  }

  debug(
	if (pairs == NULL) {
	  printf(" => Top of pairs is NULL.");
	} else {
	  pair = pairs->first;
	  printf(" => Top of pairs is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp3 = %d, genomedp3 = %d\n",*querydp3,*genomedp3);
	);

  *peeled_pairs = peeled;
  return pairs;
}
#endif


#if 0
/* Not sure if we need this, or if it causes GMAP to fail on some alignments */
/* Instead of maxpeelback, follow the 5' intron until we get enough mismatches */
static List_T
peel_rightward_intron_contiguous (int *n_peeled_indels, bool *protectedp, List_T *peeled_pairs, List_T pairs, int *querydp3, Chrpos_T *genomedp3, 
				  Chrpos_T genomedp5, bool stop_at_indels_p, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
				  int minpeelback, int min_mismatches) {
  List_T peeled = NULL;
  Pair_T pair, leftpair;
  int npeelback = 0, nmismatches = 0, niter;
  char cdna, intron_nt, intron_nt_alt;
#if 0
  int incursion = 0;
#endif
  int maxpeelback = 12;
  int last_querypos;
  Chrpos_T last_genomepos;


  *n_peeled_indels = 0;
  /* *protectedp = false; -- set by calling procedure */

  debug(printf("Peeling rightward with genomedp5 %d and stop_at_indels_p %d:",genomedp5,stop_at_indels_p));

  /* Remove initial gaps */
  while (pairs != NULL && 
	 ( ((Pair_T) pairs->first)->gapp == true ||
	   ((Pair_T) pairs->first)->comp == INDEL_COMP ||
	   ((Pair_T) pairs->first)->comp == SHORTGAP_COMP )) {
    pairs = Pairpool_pop(pairs,&pair);
  }

  if (pairs == NULL) {
    debug(printf(" pairs is empty\n"));

  } else if (stop_at_indels_p == true) {
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    /* Peel initial indels anyway */
    while (pairs != NULL && ( ((Pair_T) pairs->first)->comp == INDEL_COMP || ((Pair_T) pairs->first)->comp == INDEL_COMP )) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    if (pairs != NULL) {
      last_querypos = ((Pair_T) pairs->first)->querypos;
      last_genomepos = ((Pair_T) pairs->first)->genomepos;
    }
    while (/*npeelback < maxpeelback &&*/
	   (npeelback < minpeelback || nmismatches < min_mismatches) && pairs != NULL &&
	   ((Pair_T) pairs->first)->gapp == false &&
	   ((Pair_T) pairs->first)->comp != INDEL_COMP &&
	   ((Pair_T) pairs->first)->comp != SHORTGAP_COMP &&
	   ((Pair_T) pairs->first)->querypos <= last_querypos + 1 &&
	   ((Pair_T) pairs->first)->genomepos <= last_genomepos + 1) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));

      intron_nt = get_genomic_nt(&intron_nt_alt,genomedp5++,chroffset,chrhigh,watsonp);
      if ((cdna = ((Pair_T) pairs->first)->cdna) != intron_nt && cdna != intron_nt_alt) {
	nmismatches++;
	debug(printf(" (3) Intron mismatch #%d: %c != %c or %c at %u\n",nmismatches,cdna,intron_nt,intron_nt_alt,genomedp5-1));
      }

      if (((Pair_T) pairs->first)->protectedp == true) {
	*protectedp = true;
      }
      last_querypos = ((Pair_T) pairs->first)->querypos;
      last_genomepos = ((Pair_T) pairs->first)->genomepos;
      peeled = List_transfer_one(peeled,&pairs);
      npeelback++;
    }

  } else {
    /* Don't stop at indels, but do stop at gaps */
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Peel known gap */
      debug(printf(" Known_gap"));
      peeled = List_transfer_one(peeled,&pairs);
    }

    niter = 0;
    if (pairs != NULL) {
      last_querypos = ((Pair_T) pairs->first)->querypos;
      last_genomepos = ((Pair_T) pairs->first)->genomepos;
    }
    while (/*npeelback < maxpeelback &&*/
	   (npeelback < minpeelback || nmismatches < min_mismatches) && niter < MAXITER && pairs != NULL &&
	   ((Pair_T) pairs->first)->gapp == false &&
	   ((Pair_T) pairs->first)->querypos <= last_querypos + 1 &&
	   ((Pair_T) pairs->first)->genomepos <= last_genomepos + 1) {
      debug(printf(" Peel [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));

      intron_nt = get_genomic_nt(&intron_nt_alt,genomedp5++,chroffset,chrhigh,watsonp);
      if ((cdna = ((Pair_T) pairs->first)->cdna) != intron_nt && cdna != intron_nt_alt) {
	nmismatches++;
	debug(printf(" (4) Intron mismatch #%d: %c != %c or %c at %u\n",nmismatches,cdna,intron_nt,intron_nt_alt,genomedp5-1));
      }

      if (((Pair_T) pairs->first)->comp == MATCH_COMP || ((Pair_T) pairs->first)->comp == DYNPROG_MATCH_COMP || ((Pair_T) pairs->first)->comp == AMBIGUOUS_COMP) {
	npeelback++;
      } else if (((Pair_T) pairs->first)->comp == INDEL_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else if (((Pair_T) pairs->first)->comp == SHORTGAP_COMP) {
	*n_peeled_indels += 1;
	npeelback--;
      } else {
	npeelback--;
      }
      if (((Pair_T) pairs->first)->protectedp == true) {
	*protectedp = true;
      }
      niter++;
      last_querypos = ((Pair_T) pairs->first)->querypos;
      last_genomepos = ((Pair_T) pairs->first)->genomepos;
      peeled = List_transfer_one(peeled,&pairs);
    }

    if (pairs != NULL && ((Pair_T) pairs->first)->gapp == true) {
      debug(printf(" Hit gap [");
	    Pair_dump_one(pairs->first,/*zerobasedp*/true);
	    printf("]"));
    }
  }
	   
  if (pairs != NULL &&
      ( ((Pair_T) pairs->first)->gapp == true || 
	((Pair_T) pairs->first)->comp == INDEL_COMP ||
	((Pair_T) pairs->first)->comp == SHORTGAP_COMP )) {
    /* Don't leave a gap or indel on the top of the pairs */
    while (peeled != NULL &&
	   ( ((Pair_T) peeled->first)->gapp == true ||
	     ((Pair_T) peeled->first)->comp == INDEL_COMP ||
	     ((Pair_T) peeled->first)->comp == SHORTGAP_COMP)) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      pairs = List_transfer_one(pairs,&peeled);
    }
    if (peeled != NULL) {
      debug(printf(" Putback [");
	    Pair_dump_one(peeled->first,/*zerobasedp*/true);
	    printf("]"));
      pairs = List_transfer_one(pairs,&peeled); /* This should be match or mismatch */
    }
  }

  if (pairs != NULL) {
    leftpair = pairs->first;
    *querydp3 = leftpair->querypos - 1;
    *genomedp3 = leftpair->genomepos - 1;
  } else if (peeled != NULL) {
    leftpair = peeled->first;
    *querydp3 = leftpair->querypos;
    *genomedp3 = leftpair->genomepos;
  } else {
    /* fprintf(stderr,"In peel_rightward, pairs and peeled are both NULL\n"); */
    /* abort(); */
  }

  debug(
	if (pairs == NULL) {
	  printf(" => Top of pairs is NULL.");
	} else {
	  pair = pairs->first;
	  printf(" => Top of pairs is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp3 = %d, genomedp3 = %d\n",*querydp3,*genomedp3);
	);

  *peeled_pairs = peeled;
  return pairs;
}
#endif


/************************************************************************
 *  Traversal functions
 ************************************************************************/

/* For peel_rightward and peel_leftward, we set quit_on_gap_p = true,
   because we want to merge gaps in initial smoothing steps */

static List_T
traverse_single_gap (bool *filledp, int *dynprogindex, List_T pairs, List_T *path, 
		     Pair_T leftpair, Pair_T rightpair,
		     Univcoord_T chroffset, Univcoord_T chrhigh,
		     char *queryseq_ptr, char *queryuc_ptr, int querylength,
		     bool watsonp, int genestrand, bool jump_late_p,
		     Pairpool_T pairpool, Dynprog_T dynprog,
		     Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3,
		     int maxpeelback, double defect_rate, bool forcep, bool finalp) {
  List_T gappairs, peeled_pairs, peeled_path;
  int queryjump, genomejump;
  int querydp5, querydp3;
  Chrpos_T genomedp5, genomedp3;
  int nmatches, nmismatches, nopens, nindels;
  int unknowns, qopens, qindels, topens, tindels, ncanonical, nsemicanonical, nnoncanonical;
  int finalscore, origscore;
  bool protectedp;
  int n_peeled_indels;
  double min_splice_prob;
  /* int origqueryjump, origgenomejump; */

  debug(printf("\nTRAVERSE_SINGLE_GAP\n"));
  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  /* if (leftpair->cdna == ' ') querydp5--; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') genomedp5--; -- For old dynamic programming */
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

  /* origqueryjump = querydp3 - querydp5 + 1; */
  /* origgenomejump = genomedp3 - genomedp5 + 1; */

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
  protectedp = false;
  pairs = peel_rightward(&n_peeled_indels,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
			 maxpeelback,/*stop_at_indels_p*/true);
  *path = peel_leftward(&n_peeled_indels,&protectedp,&peeled_path,*path,&querydp5,&genomedp5,
			maxpeelback,/*stop_at_indels_p*/true);

  if (last_genomedp5 != NULL) {
    if (querydp5 < 0) {
      querydp5 = 0;
    }
    if (querydp3 >= querylength) {
      querydp3 = querylength - 1;
    }
    if (0 && finalp == false && genomedp5 == last_genomedp5[querydp5] && genomedp3 == last_genomedp3[querydp3]) {
      debug(printf("Already solved for %u..%u at %d..%d\n",genomedp5,genomedp3,querydp5,querydp3));

      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);

      *filledp = false;		/* This replaces the gap */
      return pairs;
    }
  }

  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;
  
  if (queryjump <= 0 || genomejump <= 0) {
    /* This prevents cases like queryjump 0, genomejump 1 from being solved */
    debug(printf("Unable to perform dynamic programming\n"));
    *filledp = false;

    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);

    return pairs;

  } else {
    gappairs = Dynprog_single_gap(&(*dynprogindex),&finalscore,
				  &nmatches,&nmismatches,&nopens,&nindels,dynprog,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  queryjump,genomejump,querydp5,genomedp5,
				  chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
				  extraband_single,defect_rate,/*widebandp*/true);
    if (protectedp == true) {
      debug(printf("Protecting gappairs\n"));
      Pair_protect_list(gappairs);
    }
    debug(Pair_dump_list(gappairs,true));
  }
  debug(printf("Gap pairs:\n"));
  debug(Pair_dump_list(gappairs,true));
  debug(printf("  Final score: %d\n",finalscore));

#if 0
  /* Old behavior: Depends on amount peeled */
  if (!forcep && nmismatches + nopens > nmatches) {
    /* Put back peeled pairs */
    debug(printf("Bad alignment, so undoing this solution\n"));
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    *filledp = false;
  } else {
    pairs = Pairpool_transfer(pairs,gappairs);
    *filledp = true;
  }
#else
  /* New behavior: Compares new score to orig score */
  if (forcep == true && gappairs != NULL) {
    /* Need to heed forcep, otherwise we get into bad problems with add_dual_break */
    /* But check that gappairs is not NULL, which would eliminate peeled pairs */
    /* Intended for build_dual_breaks */
    debug(printf("forcep is true, so transferring gappairs\n"));
    pairs = Pairpool_transfer(pairs,gappairs);
    *filledp = true;
  } else {
    Pair_fracidentity(&nmatches,&unknowns,&nmismatches,&qopens,&qindels, 
		      &topens,&tindels,&ncanonical,&nsemicanonical,&nnoncanonical,
		      &min_splice_prob,peeled_pairs,/*cdna_direction*/0);
    origscore = Dynprog_score(nmatches,nmismatches,qopens,qindels,topens,tindels,defect_rate);

    Pair_fracidentity(&nmatches,&unknowns,&nmismatches,&qopens,&qindels, 
		      &topens,&tindels,&ncanonical,&nsemicanonical,&nnoncanonical,
		      &min_splice_prob,peeled_path,/*cdna_direction*/0);
    origscore += Dynprog_score(nmatches,nmismatches,qopens,qindels,topens,tindels,defect_rate);
    debug(printf("  Orig score: %d, ",origscore));

    queryjump = (rightpair->querypos - leftpair->querypos - 1);
    if (queryjump > 0) {
      origscore += Dynprog_score(/*nmatches*/0,/*nmatches*/0,/*qopens*/1,/*qindels*/queryjump,
				 /*topens*/0,/*tindels*/0,defect_rate);
    }
    genomejump = (rightpair->genomepos - leftpair->genomepos - 1);
    if (genomejump > 0) {
      origscore += Dynprog_score(/*nmatches*/0,/*nmatches*/0,/*qopens*/0,/*qindels*/0,
				 /*topens*/1,/*tindels*/genomejump,defect_rate);
    }
    debug(printf("queryjump = %d, genomejump = %d, Orig score: %d\n",queryjump,genomejump,origscore));

    if (queryjump < 0 || genomejump < 0) {
      /* Cannot accept the previous queryjump or genomejump */
      debug(printf("Existing queryjump %d or genomejump %d is negative, so accepting this solution\n",
		   queryjump,genomejump));
      pairs = Pairpool_transfer(pairs,gappairs);
      *filledp = true;
#ifdef PMAP
    } else if (finalscore < 0 || finalscore < 2*origscore/3) { /* Allow for wobble */
#else
    } else if (finalscore < 0 || finalscore < origscore) {
#endif
      /* Put back peeled pairs */
      debug(printf("Forcep %d, and Bad alignment, so undoing this solution\n",forcep));
      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);
      *filledp = false;
    } else {
      debug(printf("Forcep %d, and Good alignment, so accepting this solution\n",forcep));
      pairs = Pairpool_transfer(pairs,gappairs);
      *filledp = true;
    }
  }

  if (last_genomedp5 != NULL && *filledp == true) {
    /* Save coordinates so we don't recompute this solution */
    last_genomedp5[querydp5] = genomedp5;
    last_genomedp3[querydp3] = genomedp3;
  }

#endif

  return pairs;
}

static List_T
traverse_cdna_gap (bool *filledp, bool *incompletep, int *dynprogindex_minor, int *dynprogindex_major,
		   List_T pairs, List_T *path, Pair_T leftpair, Pair_T rightpair,
		   Univcoord_T chroffset, Univcoord_T chrhigh,
		   char *queryseq_ptr, char *queryuc_ptr, int querylength,
		   bool watsonp, int genestrand, bool jump_late_p, Pairpool_T pairpool,
		   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, 
		   Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3,
		   int maxpeelback, double defect_rate, bool finalp) {
  List_T gappairs, peeled_pairs = NULL, peeled_path = NULL;
  int queryjump, genomejump;
  int querydp5, querydp3;
  Chrpos_T genomedp5, genomedp3;
  int finalscore;
  int nmatches, nmismatches, nopens, nindels;
  bool protectedp;
  int n_peeled_indels;

  debug(printf("\nTRAVERSE_CDNA_GAP\n"));
  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  /* if (leftpair->cdna == ' ') querydp5--; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') genomedp5--; -- For old dynamic programming */
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

#if 0
  if (leftpair->dynprogindex < 0 && leftpair->dynprogindex == rightpair->dynprogindex) {
    debug(printf("Re-peeling prior solution\n"));
    /* throughmismatchp = false; */
  } else {
    debug(printf("No prior solution\n"));
    /* throughmismatchp = true; */
  }
#endif

  protectedp = false;
  pairs = peel_rightward(&n_peeled_indels,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
			 maxpeelback,/*stop_at_indels_p*/true);
  *path = peel_leftward(&n_peeled_indels,&protectedp,&peeled_path,*path,&querydp5,&genomedp5,
			maxpeelback,/*stop_at_indels_p*/true);

  if (last_genomedp5 != NULL) {
    if (querydp5 < 0) {
      querydp5 = 0;
    }
    if (querydp3 >= querylength) {
      querydp3 = querylength - 1;
    }
    if (0 && finalp == false && genomedp5 == last_genomedp5[querydp5] && genomedp3 == last_genomedp3[querydp3]) {
      debug(printf("Already solved for %u..%u at %d..%d\n",genomedp5,genomedp3,querydp5,querydp3));

      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);

      *filledp = false;		/* This replaces the gap */
      return pairs;
    } else {
      last_genomedp5[querydp5] = genomedp5;
      last_genomedp3[querydp3] = genomedp3;
    }
  }


#if 0
  if (peeled_pairs == NULL || peeled_path == NULL) {
    debug(printf("Skipping this because unable to peel\n"));
    *filledp = false;
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    return pairs;
  }
#endif

  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;

  if (queryjump > 0 && queryjump <= genomejump + MININTRONLEN) {
    debug(printf("Really a single gap, not a cDNA gap, since queryjump %d <= genomejump %d + minintronlen %d\n",
		 queryjump,genomejump,MININTRONLEN));
    gappairs = Dynprog_single_gap(&(*dynprogindex_minor),&finalscore,
				  &nmatches,&nmismatches,&nopens,&nindels,dynprogM,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  queryjump,genomejump,querydp5,genomedp5,
				  chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
				  extraband_single,defect_rate,/*widebandp*/true);
    debug(Pair_dump_list(gappairs,true));
    debug(printf("  Score: %d\n",finalscore));
    pairs = Pairpool_transfer(pairs,gappairs);
    *filledp = true;

  } else {
    /* Set queryjump approximately equal to genomejump to have square
       dynamic programming matrices */
    queryjump = genomejump + extramaterial_paired;
    gappairs = Dynprog_cdna_gap(&(*dynprogindex_major),&finalscore,&(*incompletep),dynprogL,dynprogR,
				&(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				&(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
#if 0
				&(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
#endif
				/*length1L*/queryjump,/*length1R*/queryjump,/*length2*/genomejump,
				/*offset1L*/querydp5,/*revoffset1R*/querydp3,/*offset2*/genomedp5,
				chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
				extraband_paired,defect_rate);
    debug(Pair_dump_list(gappairs,true));
    *filledp = true;
    if (gappairs == NULL) {
      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,
				      /*leftpair*/(*path)->first,/*rightpair*/pairs->first,/*knownp*/false);
    } else {
      pairs = Pairpool_transfer(pairs,gappairs);
    }
  }

  return pairs;
}


/* genome_gap is usually an intron */
/* Do not set shiftp to false */
static List_T
traverse_genome_gap (bool *filledp, bool *shiftp, int *dynprogindex_minor, int *dynprogindex_major,
		     List_T pairs, List_T *path, Pair_T leftpair, Pair_T rightpair,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     char *queryseq_ptr, char *queryuc_ptr, int querylength,
		     int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		     Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3,
		     int maxpeelback, double defect_rate, bool finalp, bool simplep) {
  List_T gappairs, peeled_pairs = NULL, peeled_path = NULL, p;
  Pair_T pair;
  int queryjump, genomejump;
  int querydp5, querydp3;
  Chrpos_T genomedp5, genomedp3, orig_genomedp5, orig_genomedp3;
  int minpeelback, min_mismatches;
  int new_leftgenomepos, new_rightgenomepos;
  double left_prob, right_prob;
  int finalscore, nmatches, nmismatches, nopens, nindels, exonhead, introntype;
  int acceptable_nmismatches;
  bool stop_at_indels_p, protectedp;
  int n_peeled_indels_rightward, n_peeled_indels_leftward;
#ifndef GSNAP
  double prob2, prob3;
#endif
#ifndef PMAP
  List_T micropairs;
  int microintrontype;
#endif

#ifdef SHORTCUT
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
#endif

  debug(printf("\nTRAVERSE_GENOME_GAP\n"));

  stop_at_indels_p = false;	/* ? true when finalp == true */

  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  /* if (leftpair->cdna == ' ') querydp5--; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') genomedp5--; -- For old dynamic programming */
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

#if 0
  if (leftpair->dynprogindex < 0 && leftpair->dynprogindex == rightpair->dynprogindex) {
    debug(printf("Re-peeling prior solution\n"));
    /* throughmismatchp = false; */
  } else {
    debug(printf("No prior solution\n"));
    /* throughmismatchp = true; */
  }
#endif

#ifdef SHORTCUT
  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;

  if (querydp5 != querydp3 + 1) {
    protectedp = false;
    pairs = peel_rightward(&n_peeled_indels_rightward,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
			   maxpeelback,stop_at_indels_p);
    *path = peel_leftward(&n_peeled_indels_leftward,&protectedp,&peeled_path,*path,&querydp5,&genomedp5,
			  maxpeelback,stop_at_indels_p);

  } else {

#if 0
    left1 = genomicuc_ptr[genomedp5];
    left2 = genomicuc_ptr[genomedp5+1];
    right2 = genomicuc_ptr[genomedp3-1];
    right1 = genomicuc_ptr[genomedp3];
#else
    left1 = get_genomic_nt(&left1_alt,genomedp5,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,genomedp5+1,chroffset,chrhigh,watsonp);
    right2 = get_genomic_nt(&right2_alt,genomedp3-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,genomedp3,chroffset,chrhigh,watsonp);
#endif

    introntype = Intron_type(left1,left2,right2,right1,
			     left1_alt,left2_alt,right2_alt,right1_alt,
			     cdna_direction);
    debug(printf("Introntype at %u..%u is %s\n",genomedp5-1,genomedp3+1,Intron_type_string(introntype)));

    protectedp = false;
    pairs = peel_rightward(&n_peeled_indels_rightward,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
			   maxpeelback,stop_at_indels_p);
    *path = peel_leftward(&n_peeled_indels_leftward,&protectedp,&peeled_path,*path,&querydp5,&genomedp5,
			  maxpeelback,stop_at_indels_p);

    if (finalp == false && novelsplicingp == true /* && mismatch_rightward_p == false && mismatch_leftward_p == false */) {
      debug(printf("No mismatches seen\n"));
      if ((cdna_direction > 0 && introntype == GTAG_FWD)
#ifndef PMAP
	  || (cdna_direction < 0 && introntype == GTAG_REV)
#endif
	  ) {
	debug(printf("Skipping because intron is already canonical\n"));
	*filledp = false;		/* Calling procedure will replace the gap */
	pairs = Pairpool_transfer(pairs,peeled_pairs);
	*path = Pairpool_transfer(*path,peeled_path);
	return pairs;
      }
    }
  }

#else  /* SHORTCUT */

  orig_genomedp5 = genomedp5;
  orig_genomedp3 = genomedp3;

  if (defect_rate < DEFECT_HIGHQ) {
    minpeelback = 6;
    min_mismatches = 2;
  } else if (defect_rate < DEFECT_MEDQ) {
    minpeelback = 8;
    min_mismatches = 3;
  } else {
    minpeelback = 10;
    min_mismatches = 4;
  }

  protectedp = false;
  pairs = peel_rightward_intron(&n_peeled_indels_rightward,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
				orig_genomedp5,stop_at_indels_p,chroffset,chrhigh,watsonp,minpeelback,min_mismatches);
  *path = peel_leftward_intron(&n_peeled_indels_leftward,&protectedp,&peeled_path,*path,&querydp5,&genomedp5,
			       orig_genomedp3,stop_at_indels_p,chroffset,chrhigh,watsonp,minpeelback,min_mismatches);

  if (last_genomedp5 != NULL) {
    if (querydp5 < 0) {
      querydp5 = 0;
    }
    if (querydp3 >= querylength) {
      querydp3 = querylength - 1;
    }
    if (0 && finalp == false && genomedp5 == last_genomedp5[querydp5] && genomedp3 == last_genomedp3[querydp3]) {
      debug(printf("Already solved for %u..%u at %d..%d\n",genomedp5,genomedp3,querydp5,querydp3));

      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);

      *filledp = false;		/* This replaces the gap */
      return pairs;
    } else {
      last_genomedp5[querydp5] = genomedp5;
      last_genomedp3[querydp3] = genomedp3;
    }
  }

#endif


  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;

  /* genomedp5 + genomejump - 1 >= genomedp3 - genomejump + 1) ?  but doesn't work on AA669154, chr1*/
  if (queryjump > 0 && genomejump <= queryjump + MININTRONLEN) {
    debug(printf("Really a single gap, not an intron\n"));
    gappairs = Dynprog_single_gap(&(*dynprogindex_minor),&finalscore,
				  &nmatches,&nmismatches,&nopens,&nindels,dynprogM,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  queryjump,genomejump,querydp5,genomedp5,
				  chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
				  extraband_single,defect_rate,/*widebandp*/true);
    if (protectedp == true) {
      debug(printf("Protecting gappairs\n"));
      Pair_protect_list(gappairs);
    }
    debug(Pair_dump_list(gappairs,true));
    debug(printf("  Score: %d\n",finalscore));

    pairs = Pairpool_transfer(pairs,gappairs);
    *filledp = true;

  } else {
    /* Set genomejump approximately equal to queryjump to have square
       dynamic programming matrices */
    /* The canonical reward for finalp == true is too high */
    genomejump = queryjump + extramaterial_paired;
    gappairs = Dynprog_genome_gap(&(*dynprogindex_major),&finalscore,&new_leftgenomepos,&new_rightgenomepos,
				  &left_prob,&right_prob,&nmatches,&nmismatches,&nopens,&nindels,
				  &exonhead,&introntype,dynprogL,dynprogR,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  queryjump,genomejump,genomejump,querydp5,/*goffsetL*/genomedp5,/*rev_goffsetR*/genomedp3,
				  chrnum,chroffset,chrhigh,
				  cdna_direction,watsonp,genestrand,jump_late_p,pairpool,
				  extraband_paired + n_peeled_indels_leftward + n_peeled_indels_rightward,
				  defect_rate,maxpeelback,/*halfp*/false,/*finalp*/false);
    if (protectedp == true) {
      debug(printf("Protecting gappairs\n"));
      Pair_protect_list(gappairs);
    }

    if (gappairs == NULL) {
      if (simplep == true) {
	*shiftp = true;
	debug(printf("Shift, since gappairs is NULL and simplep is true\n"));
      } else {
	/* *shiftp = false; */
	debug(printf("No shift, since gappairs is NULL: intron is disallowed?\n"));
      }

    } else if (new_leftgenomepos != (int) leftpair->genomepos || new_rightgenomepos != (int) rightpair->genomepos) {
      *shiftp = true;
      debug(printf("Shift in intron location from %d..%d to %d..%d\n",
		   leftpair->genomepos,rightpair->genomepos,new_leftgenomepos,new_rightgenomepos));
    } else {
      /* *shiftp = false; */
      debug(printf("No shift in intron location\n"));
    }
    debug(Pair_dump_list(gappairs,true));
    debug(printf("  gappairs score (%d..%d, %u..%u, dir %d): %d\n",
		 querydp5,querydp3,genomedp5,genomedp3,cdna_direction,finalscore));
    debug(fprintf(stderr,"  gappairs score (%d..%d, %u..%u, dir %d): %d\n",
		  querydp5,querydp3,genomedp5,genomedp3,cdna_direction,finalscore));

#if 0
    /* prob = 1.0 - (1.0 - left_prob)*(1.0 - right_prob); */
    if (finalp == true && novelsplicingp == true && (left_prob < 0.90 || right_prob < 0.90)) {
      /* Bad intron.  See if alternative with indel is better.  Check
	 only on finalp, because earlier steps may need to iterate. */
      debug(printf("Checking alternative because found a bad intron with probs %f and %f\n",
		   left_prob,right_prob));
      gappairs_alt = Dynprog_genome_gap(&(*dynprogindex_major),&finalscore_alt,&new_leftgenomepos,&new_rightgenomepos,
					&left_prob_alt,&right_prob_alt,&nmatches,&nmismatches_alt,&nopens,&nindels,
					&exonhead,&introntype,dynprogL,dynprogR,
					&(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
					queryjump,genomejump,genomejump,querydp5,genomedp5,genomedp3,
					chrnum,chroffset,chrhigh,cdna_direction,watsonp,genestrand,jump_late_p,pairpool,
					extraband_paired + n_peeled_indels_leftward + n_peeled_indels_rightward,
					defect_rate,maxpeelback,/*halfp*/false,/*finalp*/true);
      if (protectedp == true) {
	debug(printf("Protecting gappairs_alt\n"));
	Pair_protect_list(gappairs_alt);
      }

      debug(Pair_dump_list(gappairs_alt,true));
      debug(printf("  gappairs_alt score: %d, left prob %f, right prob %f\n",finalscore_alt,left_prob_alt,right_prob_alt));
      debug(fprintf(stderr,"  gappairs_alt score: %d\n",finalscore_alt));
      if (gappairs_alt != NULL && left_prob_alt > left_prob && right_prob_alt > right_prob) {
	debug(printf(" switching to alt\n"));
	gappairs = gappairs_alt;
      } else {
	debug(printf(" keeping original\n"));
      }
    }
#endif

    if (defect_rate < DEFECT_HIGHQ) {
      acceptable_nmismatches = 2;
    } else if (defect_rate < DEFECT_MEDQ) {
      acceptable_nmismatches = 2;
    } else {
      acceptable_nmismatches = 3;
    }

    debug(printf("nmismatches = %d, nopens = %d, nindels = %d.  acceptable nmismatches = %d\n",
		 nmismatches,nopens,nindels,acceptable_nmismatches));
    
    if (gappairs == NULL) {
      *filledp = false;
      if (simplep == true) {
	/* Put back peeled pairs */
	debug(printf("gappairs is false, but simple, so allowed\n"));
	for (p = peeled_pairs; p != NULL; p = p->rest) {
	  pair = (Pair_T) p->first;
	}
	for (p = peeled_path; p != NULL; p = p->rest) {
	  pair = (Pair_T) p->first;
	}
      } else {
	/* Put back peeled pairs, but mark pairs as disallowed */
	debug(printf("gappairs is false, so disallowed\n"));
	for (p = peeled_pairs; p != NULL; p = p->rest) {
	  pair = (Pair_T) p->first;
	  pair->disallowedp = true;
	}
	for (p = peeled_path; p != NULL; p = p->rest) {
	  pair = (Pair_T) p->first;
	  pair->disallowedp = true;
	}
      }

      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);
      introntype = NONINTRON;

#if 0
    } else if (!finalp && finalscore < 0) {
      *filledp = false;
      /* Put back peeled pairs */
      debug(printf("Not forced and finalscore is negative\n"));
      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);
      introntype = NONINTRON;
#endif

#if 0
    } else if (defect_rate > DEFECT_MEDQ) {
      /* Should look for them, especially for GSNAP short reads */
      /* Don't look for microexons in low-quality sequences */
      debug(printf("Don't look for microexon in low-quality sequence\n"));
      *filledp = true;
      pairs = Pairpool_transfer(pairs,gappairs);
#endif

    } else if (introntype != NONINTRON && nmismatches <= acceptable_nmismatches && nopens <= 1 && nindels <= 3) {
      debug(printf("introntype != NONINTRON and nmismatches, nopens, nindels low\n"));
      *filledp = true;
      pairs = Pairpool_transfer(pairs,gappairs);

    } else {
#ifdef PMAP
      *filledp = true;
      pairs = Pairpool_transfer(pairs,gappairs);
#else      
      *filledp = true;

#ifdef GSNAP
      micropairs = (List_T) NULL;
      microintrontype = NONINTRON;
#else
      /* Expensive call, because it requires BoyerMoore */
      debug(printf("Calling microexon because introntype == %d or nmismatches %d > acceptable %d or nopens %d > 1 or nindels %d > 3\n",
		   introntype,nmismatches,acceptable_nmismatches,nopens,nindels));
      micropairs = Dynprog_microexon_int(&prob2,&prob3,&(*dynprogindex_major),&microintrontype,
					 /*sequence1*/&(queryseq_ptr[querydp5]),
					 /*sequenceuc1*/&(queryuc_ptr[querydp5]),
					 /*length1*/queryjump,/*length2L:genomejump,*//*length2R:genomejump,*/
					 /*offset1*/querydp5,/*offset2L*/genomedp5,/*revoffset2R*/genomedp3,
					 cdna_direction,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,
					 watsonp,genestrand,pairpool);
#endif
      
      if (micropairs == NULL) {
	debug(printf("No microexon found\n"));
	pairs = Pairpool_transfer(pairs,gappairs);
	/* *shiftp = false; */
      } else {
	debug(Pair_dump_list(micropairs,/*zerobasedp*/true));
	debug(printf("\n"));

#if 0
	if (1 || (nindels == 0 && nmismatches < 4)) {
	  /* Have a higher standard */
	  if (prob2 >= 0.95 && prob3 >= 0.95) {
	    debug(printf("Transferring microexon pairs\n"));
	    pairs = Pairpool_transfer(pairs,micropairs);
	    introntype = microintrontype;
	    *shiftp = true;
	  } else {
	    pairs = Pairpool_transfer(pairs,gappairs);
	  }
	} else {
	  /* Have a lower standard */
	  if (prob2 >= 0.90 || prob3 >= 0.90) {
	    debug(printf("Transferring microexon pairs\n"));
	    pairs = Pairpool_transfer(pairs,micropairs);
	    introntype = microintrontype;
	    *shiftp = true;
	  } else {
	    pairs = Pairpool_transfer(pairs,gappairs);
	  }
	}
#else
	/* Just transfer */
	debug(printf("Transferring microexon pairs\n"));
	pairs = Pairpool_transfer(pairs,micropairs);
	introntype = microintrontype;
	*shiftp = true;
#endif

      }
#endif
    }
  }

  return pairs;
}


static List_T
traverse_dual_genome_gap (int *dynprogindex, List_T pairs, List_T *path, 
			  Pair_T leftpair, Pair_T rightpair, bool left_end_intron_p, bool right_end_intron_p,
			  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			  int midquerypos, int midgenomepos, 
			  char *queryseq_ptr, char *queryuc_ptr, int querylength,
			  int cdna_direction, bool watsonp, int genestrand,
			  bool jump_late_p, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR, 
			  Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3,
			  int maxpeelback, double defect_rate, bool finalp) {
  List_T single_gappairs, dual_gappairs_1 = NULL, dual_gappairs_2 = NULL,
    right_gappairs = NULL, left_gappairs = NULL, peeled_pairs, peeled_path;
  int queryjump, genomejump;
  int querydp5, querydp3;
  Chrpos_T genomedp5, genomedp3, orig_genomedp5, orig_genomedp3;
  int minpeelback, min_mismatches;
  int new_leftgenomepos, new_rightgenomepos;
  double single_left_prob, single_right_prob, dual_left_prob_1, dual_right_prob_1, dual_left_prob_2, dual_right_prob_2;
  int querydp5_dual, querydp3_dual, genomedp5_dual, genomedp3_dual;
  int querydp5_left, querydp3_left, genomedp5_left, genomedp3_left;
  int querydp5_right, querydp3_right, genomedp5_right, genomedp3_right;
  int single_nmatches = 0, dual_nmatches_1 = 0, dual_nmatches_2 = 0, left_nmatches = 0, right_nmatches = 0;
  int single_score, dual_score_1, dual_score_2, single_goodness, dual_goodness, 
    nmismatches, nopens, nindels, exonhead, right_exonhead, left_exonhead;
  int left_score, right_score, left_goodness = 0, right_goodness = 0;
  int middle_exonlength;
  int single_introntype, dual_introntype_1, dual_introntype_2, left_introntype, right_introntype;
#ifdef DEBUG
  int interexon_region;
  double middle_exonprob;
#endif
  bool singlep = false, single_canonical_p, dual_canonical_p, protectedp;
  int n_peeled_indels;

  debug(printf("\nTRAVERSE_DUAL_GENOME_GAP: left_end_intron_p %d, right_end_intron_p %d\n",
	       left_end_intron_p,right_end_intron_p));
#if 0
  if (cdna_direction > 0) {
    canonical_introntype = GTAG_FWD;
    semicanonical_introntype_1 = ATAC_FWD;
    semicanonical_introntype_2 = GCAG_FWD;
#ifndef PMAP
  } else {
    canonical_introntype = GTAG_REV;
    semicanonical_introntype_1 = ATAC_REV;
    semicanonical_introntype_2 = GCAG_REV;
#endif
  }
#endif

  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  /* if (leftpair->cdna == ' ') querydp5--; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') genomedp5--; -- For old dynamic programming */
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;


  orig_genomedp5 = genomedp5;
  orig_genomedp3 = genomedp3;

  if (defect_rate < DEFECT_HIGHQ) {
    minpeelback = 6;
    min_mismatches = 2;
  } else if (defect_rate < DEFECT_MEDQ) {
    minpeelback = 8;
    min_mismatches = 3;
  } else {
    minpeelback = 10;
    min_mismatches = 4;
  }

  protectedp = false;
  pairs = peel_rightward_intron(&n_peeled_indels,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
				orig_genomedp5,/*stop_at_indels_p*/false,chroffset,chrhigh,watsonp,minpeelback,min_mismatches);
  *path = peel_leftward_intron(&n_peeled_indels,&protectedp,&peeled_path,*path,&querydp5,&genomedp5,
			       orig_genomedp3,/*stop_at_indels_p*/false,chroffset,chrhigh,watsonp,minpeelback,min_mismatches);

  if (last_genomedp5 != NULL) {
    if (querydp5 < 0) {
      querydp5 = 0;
    }
    if (querydp3 >= querylength) {
      querydp3 = querylength - 1;
    }
    if (0 && finalp == false && genomedp5 == last_genomedp5[querydp5] && genomedp3 == last_genomedp3[querydp3]) {
      /* Don't want to abort this procedure early */
      debug(printf("Already solved for %u..%u at %d..%d\n",genomedp5,genomedp3,querydp5,querydp3));

      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);

      return pairs;
    } else {
      last_genomedp5[querydp5] = genomedp5;
      last_genomedp3[querydp3] = genomedp3;
    }
  }

  queryjump = querydp3 - querydp5 + 1;
  genomejump = queryjump + extramaterial_paired;

  if (queryjump > nullgap) {
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,
				    /*leftpair*/(*path)->first,/*rightpair*/pairs->first,/*knownp*/false);

    return pairs;
  }

  if (genomedp5 + genomejump - 1 >= genomedp3 - genomejump + 1) {
    debug(printf("Bounds don't make sense for dual intron gap: %d + %d - 1 >= %d - %d + 1\n\n",
		 genomedp5,genomejump,genomedp3,genomejump));
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,
				    /*leftpair*/(*path)->first,/*rightpair*/pairs->first,/*knownp*/false);

    return pairs;
  }
  
  /* Want finalp == true to get best chance for canonical splice site */
  single_gappairs = Dynprog_genome_gap(&(*dynprogindex),&single_score,&new_leftgenomepos,&new_rightgenomepos,
				       &single_left_prob,&single_right_prob,&single_nmatches,&nmismatches,&nopens,&nindels,
				       &exonhead,&single_introntype,dynprogL,dynprogR,
				       &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				       queryjump,genomejump,genomejump,querydp5,genomedp5,genomedp3,
				       chrnum,chroffset,chrhigh,
				       cdna_direction,watsonp,genestrand,jump_late_p,pairpool,extraband_paired,
				       defect_rate,maxpeelback,/*halfp*/false,/*finalp*/true);

  debug(Pair_check_list_pairs(single_gappairs));

  /* Okay to have one indel, because may need to shift an island */
  if (nopens <= 1) {
    single_goodness = (single_nmatches + nindels) + MISMATCH*nmismatches;
  } else {
    single_goodness = single_nmatches + MISMATCH*nmismatches + QOPEN*(nopens-1) + QINDEL*nindels;
  }

#if 0
  if (single_gappairs == NULL) {
    single_canonical_p = false;
  } else if (single_introntype == canonical_introntype) {
    single_canonical_p = true;
  } else if (single_introntype == semicanonical_introntype_1 ||
	     single_introntype == semicanonical_introntype_2) {
    single_canonical_p = false;
  } else {
    single_canonical_p = false;
  }
#else
  if (single_left_prob > 0.9 && single_right_prob > 0.9) {
    single_canonical_p = true;
  } else {
    single_canonical_p = false;
  }
#endif


  /* Right of short exon */
  querydp5_dual = midquerypos;
  genomedp5_dual = midgenomepos;
  querydp3_dual = querydp3;	/* From peel_rightward_intron */
  genomedp3_dual = genomedp3;	/* From peel_rightward_intron */

  queryjump = querydp3_dual - querydp5_dual + 1;
  genomejump = queryjump + extramaterial_paired;

  if (genomedp5_dual + genomejump - 1 >= genomedp3_dual) {
    /* Bounds don't make sense */
    debug(printf("Bounds don't make sense on right of dual intron gap: %d + %d - 1 >= %d\n\n",
		 genomedp5_dual,genomejump,genomedp3_dual));
    dual_gappairs_2 = NULL;

  } else {
    /* Want finalp == true to get best chance for canonical splice site */
    dual_gappairs_2 = Dynprog_genome_gap(&(*dynprogindex),&dual_score_2,&new_leftgenomepos,&new_rightgenomepos,
					 &dual_left_prob_2,&dual_right_prob_2,&dual_nmatches_2,&nmismatches,&nopens,&nindels,
					 &right_exonhead,&dual_introntype_2,dynprogL,dynprogR,
					 &(queryseq_ptr[querydp5_dual]),&(queryuc_ptr[querydp5_dual]),
					 queryjump,genomejump,genomejump,
					 querydp5_dual,genomedp5_dual,genomedp3_dual,
					 chrnum,chroffset,chrhigh,
					 cdna_direction,watsonp,genestrand,jump_late_p,pairpool,extraband_paired,
					 defect_rate,maxpeelback,/*halfp*/true,/*finalp*/true);

    dual_goodness = dual_nmatches_2 + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

    /* Left of short exon */
    querydp5_dual = querydp5;	/* From peel_leftward_intron */
    genomedp5_dual = genomedp5;	/* From peel_leftward_intron */
    querydp3_dual = midquerypos-1;
    genomedp3_dual = midgenomepos-1;

    queryjump = querydp3_dual - querydp5_dual + 1;
    genomejump = queryjump + extramaterial_paired;

    if (genomedp5_dual + genomejump - 1 >= genomedp3_dual) {
      /* Bounds don't make sense */
      debug(printf("Bounds don't make sense on left of dual intron gap: %d + %d - 1 >= %d\n\n",
		   genomedp5_dual,genomejump,genomedp3_dual));
      dual_gappairs_1 = NULL;

    } else {
      /* Want finalp == true to get best chance for canonical splice site */
      dual_gappairs_1 = Dynprog_genome_gap(&(*dynprogindex),&dual_score_1,&new_leftgenomepos,&new_rightgenomepos,
					   &dual_left_prob_1,&dual_right_prob_1,&dual_nmatches_1,&nmismatches,&nopens,&nindels,
					   &left_exonhead,&dual_introntype_1,dynprogL,dynprogR,
					   &(queryseq_ptr[querydp5_dual]),&(queryuc_ptr[querydp5_dual]),
					   queryjump,genomejump,genomejump,
					   querydp5_dual,genomedp5_dual,genomedp3_dual,
					   chrnum,chroffset,chrhigh,
					   cdna_direction,watsonp,genestrand,jump_late_p,pairpool,extraband_paired,
					   defect_rate,maxpeelback,/*halfp*/true,/*finalp*/true);

      dual_goodness += dual_nmatches_1 + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

      if (dual_gappairs_1 == NULL || dual_gappairs_2 == NULL) {
	dual_canonical_p = false;
#if 0
      } else if (dual_introntype_1 == canonical_introntype && dual_introntype_2 == canonical_introntype) {
	dual_canonical_p = true;
#endif
#if 0
      } else if (dual_left_prob_1 > 0.9 && dual_right_prob_1 > 0.9 &&
		 dual_left_prob_2 > 0.9 && dual_right_prob_2 > 0.9) {
	dual_canonical_p = true;
#endif
      } else if (dual_left_prob_1 > 0.9 || dual_right_prob_1 > 0.9 ||
		 dual_left_prob_2 > 0.9 || dual_right_prob_2 > 0.9) {
	dual_canonical_p = true;
      } else {
	dual_canonical_p = false;
      }
    }
  }

  if (dual_gappairs_2 == NULL || dual_gappairs_1 == NULL) {
    debug(printf("Single score wins because dual_guappairs_2 is NULL or dual_gappairs_1 is NULL\n"));
    debug(printf("Loser: dual_gappairs\n"));
    debug(Pair_dump_list(dual_gappairs_2,true));
    debug(Pair_dump_list(dual_gappairs_1,true));
    debug(printf("Winner: single gap pairs\n"));
    debug(Pair_dump_list(single_gappairs,true));
    /* pairs = Pairpool_transfer(pairs,single_gappairs); -- Wait until we check for left_goodness and right_goodness */
    singlep = true;

  } else {
    middle_exonlength = right_exonhead-left_exonhead;
    debug(printf("Middle exon is %d - %d = %d bp in interexon region of %d bp\n",
		 right_exonhead,left_exonhead,right_exonhead-left_exonhead,new_rightgenomepos-new_leftgenomepos));
    if (middle_exonlength <= 0) {
      debug(middle_exonprob = 0.0);
    } else {
      debug(interexon_region = new_rightgenomepos - new_leftgenomepos);
      
#if 0
      if (dual_introntype_2 == canonical_introntype) {
	middle_exonlength += DUAL_HALFCANONICAL_POINTS;
	debug(printf("Add canonical credit of %d for right intron\n",DUAL_HALFCANONICAL_POINTS));
      }
      if (dual_introntype_1 == canonical_introntype) {
	middle_exonlength += DUAL_HALFCANONICAL_POINTS;
	debug(printf("Add canonical credit of %d for left intron\n",DUAL_HALFCANONICAL_POINTS));
      }
#else
      if (dual_left_prob_2 > 0.9 && dual_right_prob_2 > 0.9) {
	middle_exonlength += DUAL_HALFCANONICAL_POINTS;
	debug(printf("Add canonical credit of %d for right intron\n",DUAL_HALFCANONICAL_POINTS));
      }
      if (dual_left_prob_1 > 0.9 && dual_right_prob_1 > 0.9) {
	middle_exonlength += DUAL_HALFCANONICAL_POINTS;
	debug(printf("Add canonical credit of %d for left intron\n",DUAL_HALFCANONICAL_POINTS));
      }
#endif

#ifdef DEBUG      
      middle_exonprob = 1.0-pow(1.0-pow(4.0,-(double) middle_exonlength),(double) interexon_region);
#endif
      
      debug(printf("Single score = %d (%d matches).  Single canonical: %d.  Dual score = %d & %d (%d & %d matches).  Dual canonical: %d.  ",
		   single_score,single_nmatches,single_canonical_p,
		   dual_score_1,dual_score_2,dual_nmatches_1,dual_nmatches_2,
		   dual_canonical_p));
      debug(printf("Single goodness = %d.  Dual goodness = %d.  ",
		   single_goodness,dual_goodness));
      debug(printf("Probability is %g (ignoring).  ",middle_exonprob));
    }

    /* Want high threshold for accepting dual intron */
    if (dual_canonical_p == true && /*middle_exonprob < 0.001 &&*/
	single_canonical_p == false && single_goodness <= dual_goodness) {
      debug(printf("Dual scores win\n"));
      debug(printf("Loser: single_gappairs\n"));
      debug(Pair_dump_list(single_gappairs,true));
      debug(printf("Winner: Transferring dual_gappairs_2 onto pairs\n"));
      Pair_protect_list(dual_gappairs_2);
      debug(Pair_dump_list(dual_gappairs_2,true));
      pairs = Pairpool_transfer(pairs,dual_gappairs_2);
      debug(printf("Winner: Transferring dual_gappairs_1 onto pairs\n"));
      Pair_protect_list(dual_gappairs_1);
      debug(Pair_dump_list(dual_gappairs_1,true));
      pairs = Pairpool_transfer(pairs,dual_gappairs_1);
    } else {
      debug(printf("Single score wins\n"));
      debug(printf("Loser: dual_gappairs\n"));
      debug(Pair_dump_list(dual_gappairs_2,true));
      debug(Pair_dump_list(dual_gappairs_1,true));
      debug(printf("Winner: single gappairs\n"));
      debug(Pair_dump_list(single_gappairs,true));
      /* pairs = Pairpool_transfer(pairs,single_gappairs); -- Wait until we check for left_goodness and right_goodness */
      singlep = true;
    }
  }

  if (singlep == true) {
    if (single_gappairs == NULL) {
      /* Need to handle the possibility that Dynprog_genome_gap artificially returns NULL for single_gappairs */
      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);

    } else {
      if (right_end_intron_p == true) {
	/* Keep left intron only and extend right from short exon */
	querydp5_right = querydp5;
	genomedp5_right = genomedp5;
	querydp3_right = midquerypos;
	genomedp3_right = midgenomepos;

	queryjump = querydp3_right - querydp5_right + 1;
	genomejump = queryjump + extramaterial_paired;

	if (genomedp5_right + genomejump - 1 >= genomedp3_right) {
	  /* Bounds don't make sense */
	  debug(printf("Bounds don't make sense if we omit right intron: %d + %d - 1 >= %d\n\n",
		       genomedp5_right,genomejump,genomedp3_right));
	  right_gappairs = NULL;

	} else {
	  right_gappairs = Dynprog_genome_gap(&(*dynprogindex),&right_score,&new_leftgenomepos,&new_rightgenomepos,
					      &single_left_prob,&single_right_prob,&right_nmatches,&nmismatches,&nopens,&nindels,
					      &right_exonhead,&right_introntype,dynprogL,dynprogR,
					      &(queryseq_ptr[querydp5_right]),&(queryuc_ptr[querydp5_right]),
					      queryjump,genomejump,genomejump,
					      querydp5_right,genomedp5_right,genomedp3_right,
					      chrnum,chroffset,chrhigh,
					      cdna_direction,watsonp,genestrand,jump_late_p,pairpool,extraband_paired,
					      defect_rate,maxpeelback,/*halfp*/false,/*finalp*/false);

	  right_goodness = right_nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;
	  debug(printf("Right goodness (keeping left intron only) = %d\n",right_goodness));

	  if (right_goodness > single_goodness) {
	    debug(printf("New winner: right gappairs\n"));
	    debug(Pair_dump_list(right_gappairs,true));
	    single_gappairs = right_gappairs;
	    single_goodness = right_goodness;
	  }
	}
      }

      if (left_end_intron_p == true) {
	/* Keep right intron only and extend left from short exon */
	querydp5_left = midquerypos;
	genomedp5_left = midgenomepos;
	querydp3_left = querydp3;
	genomedp3_left = genomedp3;

	queryjump = querydp3_left - querydp5_left + 1;
	genomejump = queryjump + extramaterial_paired;

	if (genomedp5_left + genomejump - 1 >= genomedp3_left) {
	  /* Bounds don't make sense */
	  debug(printf("Bounds don't make sense if we omit left intron: %d + %d - 1 >= %d\n\n",
		       genomedp5_left,genomejump,genomedp3_left));
	  left_gappairs = NULL;

	} else {
	  left_gappairs = Dynprog_genome_gap(&(*dynprogindex),&left_score,&new_leftgenomepos,&new_rightgenomepos,
					     &single_left_prob,&single_right_prob,&left_nmatches,&nmismatches,&nopens,&nindels,
					     &left_exonhead,&left_introntype,dynprogL,dynprogR,
					     &(queryseq_ptr[querydp5_left]),&(queryuc_ptr[querydp5_left]),
					     queryjump,genomejump,genomejump,
					     querydp5_left,genomedp5_left,genomedp3_left,
					     chrnum,chroffset,chrhigh,
					     cdna_direction,watsonp,genestrand,jump_late_p,pairpool,extraband_paired,
					     defect_rate,maxpeelback,/*halfp*/false,/*finalp*/false);
      
	  left_goodness = left_nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;
	  debug(printf("Left goodness (keeping right intron only) = %d\n",left_goodness));

	  if (left_goodness > single_goodness) {
	    debug(printf("New winner: left gappairs\n"));
	    debug(Pair_dump_list(left_gappairs,true));
	    single_gappairs = left_gappairs;
	    single_goodness = left_goodness;
	  }
	}
      }

      /* Finally transfer best single result */
      if (single_gappairs == right_gappairs) {
	pairs = Pairpool_transfer(pairs,peeled_pairs);
      }
      pairs = Pairpool_transfer(pairs,single_gappairs);
      if (single_gappairs == left_gappairs) {
	*path = Pairpool_transfer(*path,peeled_path);
      }
    }

  }

  return pairs;
}


#if 0
static bool
good_end_intron_p (Pair_T gappair, int cdna_direction) {
  if (gappair->knowngapp == true) {
    return true;

  } else if (cdna_direction > 0) {
    if (gappair->comp == FWD_CANONICAL_INTRON_COMP || (gappair->donor_prob >= 0.90 && gappair->acceptor_prob >= 0.90)) {
      return true;
    } else {
      return false;
    }
  } else if (cdna_direction < 0) {
    if (gappair->comp == REV_CANONICAL_INTRON_COMP || (gappair->donor_prob >= 0.90 && gappair->acceptor_prob >= 0.90)) {
      return true;
    } else {
      return false;
    }
  } else {
    if (gappair->comp == FWD_CANONICAL_INTRON_COMP || gappair->comp == REV_CANONICAL_INTRON_COMP || 
	(gappair->donor_prob >= 0.90 && gappair->acceptor_prob >= 0.90)) {
      return true;
    } else {
      return false;
    }
  }
}
#endif


/* Note on QUERYEND_INDELS.  Profiling shows that using
   QUERYEND_INDELS caused compute_scores_lookup to be called too
   often, slowing program down. */


/* to_queryend_p must be true for distalmedial_ending, since we are
   comparing alternatives.  But it can be false for extend_ending,
   which just tries to improve the ends. */

static List_T
distalmedial_ending5 (bool *knownsplicep, bool *chop_exon_p, int *dynprogindex_minor,
		      int *finalscore, int *ambig_end_length, double *ambig_prob,
		      List_T *pairs, int leftquerypos, Pair_T rightpair,
		      Univcoord_T chroffset, Univcoord_T chrhigh,
		      char *queryseq_ptr, char *queryuc_ptr,
		      bool watsonp, int genestrand, bool jump_late_p, Pairpool_T pairpool,
		      Dynprog_T dynprog, int maxpeelback, double defect_rate) {
  List_T peeled_pairs, continuous_gappairs_medialgap = NULL;
  int queryjump, genomejump;
  int querydp5, querydp3_distalgap, querydp3_medialgap;
  Chrpos_T genomedp3_distalgap, genomedp3_medialgap;
  int continuous_goodness_distalgap = 0, continuous_goodness_medialgap = 0,
    nmatches, nmismatches, nopens, nindels;
  bool protectedp;
  int n_peeled_indels;
  bool knownsplice_medial_p = false;

  debug(printf("\nDISTALMEDIAL_ENDING5\n"));

  querydp5 = leftquerypos + 1;
#if 0
  genomedp5 = leftgenomepos + 1; /* 0 */
#endif
  querydp3_distalgap = querydp3_medialgap = rightpair->querypos - 1;
  genomedp3_distalgap = genomedp3_medialgap = rightpair->genomepos - 1;

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
  protectedp = false;
  *pairs = peel_rightward(&n_peeled_indels,&protectedp,&peeled_pairs,*pairs,&querydp3_distalgap,&genomedp3_distalgap,
			  maxpeelback,/*stop_at_indels_p*/true);

  continuous_goodness_distalgap = Pair_fracidentity_score(peeled_pairs);
  /* continuous_goodness_distalgap += Pair_fracidentity_score(endgappairs); */
  debug(printf("continuous_goodness_distalgap (%d pairs) is %d\n",
	       List_length(peeled_pairs),continuous_goodness_distalgap));

#if 0
  /* gappair wasn't initialized */
  if (good_end_intron_p(gappair,cdna_direction) == false) {
    debug(printf("Subtracting points from continuous distal because noncanonical\n"));
    continuous_goodness_distalgap -= CANONICAL_POINTS;
  } else if (gappair->comp == DUALBREAK_COMP) {
    debug(printf("Subtracting points from continuous distal because of dual break\n"));
    continuous_goodness_distalgap -= (CANONICAL_POINTS + CANONICAL_POINTS);
  }
#endif

  /* Solve if gap were not present */
  queryjump = querydp3_medialgap - querydp5 + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

#ifdef EXTRACT_GENOMICSEG
  genomedp5 = genomedp3_medialgap - genomejump + 1;
  /* Make sure we don't go past the beginning */
  if (genomedp5 < 0) {
    genomedp5 = 0;
    genomejump = genomedp3_medialgap - genomedp5 + 1;
  }
#endif

  debug(printf("Stage 3: traverse_ending5: Dynamic programming at 5' end (medial to gap): querydp5 = %d, querydp3 = %d, genomedp3 = %d\n",
	       querydp5,querydp3_medialgap,genomedp3_medialgap));

  continuous_gappairs_medialgap = Dynprog_end5_gap(&(*dynprogindex_minor),&(*finalscore),
						   &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						   &(queryseq_ptr[querydp3_medialgap]),&(queryuc_ptr[querydp3_medialgap]),
						   queryjump,genomejump,querydp3_medialgap,genomedp3_medialgap,
						   chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
						   extraband_end,defect_rate,/*endalign*/QUERYEND_INDELS,/*require_pos_score_p*/false);
  *ambig_end_length = 0;
  *ambig_prob = 0.0;

  continuous_goodness_medialgap = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;
  debug(printf("Continuous_goodness_medialgap %d = %d + %d*%d + %d*%d + %d*%d\n",
	       continuous_goodness_medialgap,nmatches,MISMATCH,nmismatches,QOPEN,nopens,QINDEL,nindels));
  
  if (continuous_goodness_distalgap > continuous_goodness_medialgap) {
    debug(printf("Continuous distal wins: %d > %d\n",continuous_goodness_distalgap,continuous_goodness_medialgap));
    *ambig_end_length = 0;
    *ambig_prob = 0.0;

    /* *pairs = Pairpool_transfer(*pairs,endgappairs); */
    *chop_exon_p = false;
    /* Let previous value of knownsplicep stand */
    debug(printf("Returning peeled pairs:\n"));
    debug(Pair_dump_list(peeled_pairs,true));
    debug(printf("\n"));
    return peeled_pairs;

  } else {
    debug(printf("Continuous medial wins: %d > %d\n",
		 continuous_goodness_medialgap,continuous_goodness_distalgap));
    *chop_exon_p = true;
    *knownsplicep = knownsplice_medial_p;
    return continuous_gappairs_medialgap;
  }
}


static List_T
extend_ending5 (bool *knownsplicep, int *dynprogindex_minor,
		int *finalscore, int *ambig_end_length, Splicetype_T *ambig_splicetype, double *ambig_prob,
		List_T *pairs, int leftquerypos, Pair_T rightpair,
		Univcoord_T chroffset, Univcoord_T chrhigh,
		Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		char *queryseq_ptr, char *queryuc_ptr,
		int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		Pairpool_T pairpool, Dynprog_T dynprog, int maxpeelback,
		double defect_rate, Endalign_T endalign, bool forcep) {
  List_T continuous_gappairs_distalgap = NULL, peeled_pairs = NULL;
  int queryjump, genomejump;
  int querydp5, querydp3_distalgap;
  Chrpos_T genomedp3_distalgap;
  int nmatches, nmismatches, nopens, nindels;
  bool protectedp = false;
  int n_peeled_indels = 0;
  Pair_T firstpair;


  debug(printf("\nEXTEND_ENDING5 with endalign %s and maxpeelback %d\n",
	       Dynprog_endalign_string(endalign),maxpeelback));

  querydp5 = leftquerypos + 1;
#if 0
  genomedp5 = leftgenomepos + 1; /* 0 */
#endif
  querydp3_distalgap = rightpair->querypos - 1;
  genomedp3_distalgap = rightpair->genomepos - 1;

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */

  if (endalign == QUERYEND_NOGAPS) {
    /* Don't peelback on extension */
  } else if (maxpeelback == 0) {
    /* Actually, we should peelback after trim_ends, because indel placement could be wrong */
    /* Don't peelback on BEST_LOCAL after trim_ends */
  } else {
    protectedp = false;
    *pairs = peel_rightward(&n_peeled_indels,&protectedp,&peeled_pairs,*pairs,&querydp3_distalgap,&genomedp3_distalgap,
			    maxpeelback,/*stop_at_indels_p*/true);
  }
  
  queryjump = querydp3_distalgap - querydp5 + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

#if 0
  genomedp5 = genomedp3_distalgap - genomejump + 1;
#endif
#ifdef EXTRACT_GENOMICSEG
  /* Make sure we don't go past the beginning */
  if (genomedp5 < 0) {
    genomedp5 = 0;
    genomejump = genomedp3_distalgap - genomedp5 + 1;
  }
#endif

  debug(printf("Stage 3 (dir %d), extend_ending5: Dynamic programming at 5' end (distal to gap): querydp5 = %d, querydp3 = %d, genomedp3 = %d\n",
	       cdna_direction,querydp5,querydp3_distalgap,genomedp3_distalgap));


  if (endalign == QUERYEND_GAP && splicesites != NULL) {
    continuous_gappairs_distalgap = Dynprog_end5_known(&(*knownsplicep),&(*dynprogindex_minor),&(*finalscore),
						       &(*ambig_end_length),&(*ambig_splicetype),
						       &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						       &(queryseq_ptr[querydp3_distalgap]),&(queryuc_ptr[querydp3_distalgap]),
						       queryjump,genomejump,querydp3_distalgap,genomedp3_distalgap,
						       chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
						       cdna_direction,watsonp,genestrand,jump_late_p,pairpool,
						       extraband_end,defect_rate);
    if (*ambig_end_length > 0) {
      *ambig_prob = 2.0;
    }
  } else {
    continuous_gappairs_distalgap = Dynprog_end5_gap(&(*dynprogindex_minor),&(*finalscore),
						     &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						     &(queryseq_ptr[querydp3_distalgap]),&(queryuc_ptr[querydp3_distalgap]),
						     /*rlength*/queryjump,/*glength*/genomejump,
						     /*rev_roffset*/querydp3_distalgap,/*rev_goffset*/genomedp3_distalgap,
						     chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
						     extraband_end,defect_rate,endalign,/*require_pos_score_p*/false);
    *ambig_end_length = 0;
    *ambig_prob = 0.0;
    *knownsplicep = false;
  }

  debug(printf("  finalscore: %d\n",*finalscore));
  if (continuous_gappairs_distalgap == NULL) {
    return peeled_pairs;
  } else {
    firstpair = List_head(continuous_gappairs_distalgap);
    if (0 && firstpair->querypos != querydp3_distalgap) {
      /* Not a good test anymore, since we are halting peelbacks at gaps */
      /* Must have an indel between the gappairs and the rest of the read */
      debug(printf("Detected indel between gappairs %d and the rest of the read %d\n",
		   firstpair->querypos,querydp3_distalgap));
      return peeled_pairs;

    } else if (forcep == true) {
      /* For example, needed for extending middles of chimeras */
      return continuous_gappairs_distalgap;
    } else if (*finalscore <= 0) {
      /* Need this branch to eliminate bad exons at ends */
      *knownsplicep = false;
      return peeled_pairs;
    } else {
      return continuous_gappairs_distalgap;
    }
  }
}


static List_T
distalmedial_ending3 (bool *knownsplicep, bool *chop_exon_p, int *dynprogindex_minor,
		      int *finalscore, int *ambig_end_length, double *ambig_prob,
		      List_T *path, Pair_T leftpair, int rightquerypos,
		      Univcoord_T chroffset, Univcoord_T chrhigh,
		      char *queryseq_ptr, char *queryuc_ptr,
		      bool watsonp, int genestrand, bool jump_late_p,
		      Pairpool_T pairpool, Dynprog_T dynprog, int maxpeelback, double defect_rate) {
  List_T peeled_path, continuous_gappairs_medialgap = NULL;
  int queryjump, genomejump;
  int querydp5_distalgap, querydp3, querydp5_medialgap;
  Chrpos_T genomedp5_distalgap, genomedp5_medialgap;
  int continuous_goodness_distalgap = 0, continuous_goodness_medialgap = 0,
    nmatches, nmismatches, nopens, nindels;
  bool protectedp;
  int n_peeled_indels;
  bool knownsplice_medial_p = false;


  debug(printf("\nDISTALMEDIAL_ENDING3\n"));
  
  querydp5_distalgap = leftpair->querypos + 1;
  genomedp5_distalgap = leftpair->genomepos + 1;
  /* if (leftpair->cdna == ' ') querydp5_distalgap--; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') genomedp5_distalgap--; -- For old dynamic programming */

  querydp5_medialgap = querydp5_distalgap;
  genomedp5_medialgap = genomedp5_distalgap;
  querydp3 = rightquerypos - 1;
  /* genomedp3 = rightgenomepos - 1; */

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
  protectedp = false;
  *path = peel_leftward(&n_peeled_indels,&protectedp,&peeled_path,*path,&querydp5_distalgap,&genomedp5_distalgap,
			maxpeelback,/*stop_at_indels_p*/true);
  
  continuous_goodness_distalgap = Pair_fracidentity_score(peeled_path);
  /* continuous_goodness_distalgap += Pair_fracidentity_score(endgappairs); */
  debug(printf("continuous_goodness_distalgap (%d pairs) is %d\n",
	       List_length(peeled_path),continuous_goodness_distalgap));

#if 0
  /* gappair wasn't initialized */
  if (good_end_intron_p(gappair,cdna_direction) == false) {
    debug(printf("Subtracting points from continuous distal because noncanonical\n"));
    continuous_goodness_distalgap -= CANONICAL_POINTS;
  } else if (gappair->comp == DUALBREAK_COMP) {
    debug(printf("Subtracting points from continuous distal because of dual break\n"));
    continuous_goodness_distalgap -= (CANONICAL_POINTS + CANONICAL_POINTS);
  }
#endif

  /* Solve if gap were not present */
  queryjump = querydp3 - querydp5_medialgap + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

#ifdef EXTRACT_GENOMICSEG
  genomedp3 = genomedp5_medialgap + genomejump - 1;
  /* Make sure we don't go past the end */
  if (genomedp3 > genomiclength - 1) {
    genomedp3 = genomiclength - 1;
    genomejump = genomedp3 - genomedp5_medialgap + 1;
  }
#endif
    
  debug(printf("Stage 3: distalmedial_ending3: Dynamic programming at 3' end (medial to gap): querydp5 = %d, querydp3 = %d, genomedp5 = %u\n",
	       querydp5_medialgap,querydp3,genomedp5_medialgap));

  debug(printf("Before solving the 3' end, here is the path:\n"));
  debug(Pair_dump_list(*path,true));
  debug(printf("\n"));

  continuous_gappairs_medialgap = Dynprog_end3_gap(&(*dynprogindex_minor),&(*finalscore),
						   &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						   &(queryseq_ptr[querydp5_medialgap]),&(queryuc_ptr[querydp5_medialgap]),
						   queryjump,genomejump,querydp5_medialgap,genomedp5_medialgap,
						   chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
						   extraband_end,defect_rate,/*endalign*/QUERYEND_INDELS,
						   /*require_pos_score_p*/false);
  *ambig_end_length = 0;
  *ambig_prob = 0.0;

  debug(printf("Medial gap\n"));
  debug(Pair_dump_list(continuous_gappairs_medialgap,true));

  continuous_goodness_medialgap = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;
  debug(printf("Continuous_goodness_medialgap %d = %d + %d*%d + %d*%d + %d*%d\n",
	       continuous_goodness_medialgap,nmatches,MISMATCH,nmismatches,QOPEN,nopens,QINDEL,nindels));

  if (continuous_goodness_distalgap > continuous_goodness_medialgap) {
    debug(printf("Continuous distal wins: %d > %d\n",continuous_goodness_distalgap,continuous_goodness_medialgap));
    *ambig_end_length = 0;
    *ambig_prob = 0.0;

    /* *path = Pairpool_transfer(*path,endgappairs); */
    *chop_exon_p = false;
    /* Let previous value of knownsplicep stand */
    debug(printf("Returning peeled path:\n"));
    debug(Pair_dump_list(peeled_path,true));
    debug(printf("\n"));
    return peeled_path;

  } else {
    debug(printf("Continuous medial wins: %d > %d\n",continuous_goodness_medialgap,continuous_goodness_distalgap));
    *chop_exon_p = true;
    *knownsplicep = knownsplice_medial_p;
    return List_reverse(continuous_gappairs_medialgap);
  }
}


static List_T
extend_ending3 (bool *knownsplicep, int *dynprogindex_minor, int *finalscore,
		int *ambig_end_length, Splicetype_T *ambig_splicetype, double *ambig_prob,
		List_T *path, Pair_T leftpair, int rightquerypos,
		int querylength, Univcoord_T chroffset, Univcoord_T chrhigh,
		Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		char *queryseq_ptr, char *queryuc_ptr,
		int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		Pairpool_T pairpool, Dynprog_T dynprog, int maxpeelback,
		double defect_rate, Endalign_T endalign, bool forcep) {
  List_T continuous_gappairs_distalgap = NULL, peeled_path = NULL;
  int queryjump, genomejump;
  int querydp5_distalgap, querydp3;
  Chrpos_T genomedp5_distalgap;
  int nmatches, nmismatches, nopens, nindels;
  bool protectedp = false;
  int n_peeled_indels = 0;
  Pair_T firstpair;

  debug(printf("\nEXTEND_ENDING3 with endalign %s and maxpeelback %d\n",
	       Dynprog_endalign_string(endalign),maxpeelback));
  
  querydp5_distalgap = leftpair->querypos + 1;
  genomedp5_distalgap = leftpair->genomepos + 1;
  /* if (leftpair->cdna == ' ') querydp5_distalgap--; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') genomedp5_distalgap--; -- For old dynamic programming */
  querydp3 = rightquerypos - 1;
  /* genomedp3 = rightgenomepos - 1; */
  debug(printf("Set dynprog 3 end to be querydp3 = %d\n",querydp3));

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
  if (endalign == QUERYEND_NOGAPS) {
    /* Don't peelback on extension */
  } else if (maxpeelback == 0) {
    /* Actually, we should peelback after trim_ends, because indel placement could be wrong */
    /* Don't peelback on BEST_LOCAL after trim_ends */
  } else {
    protectedp = false;
    *path = peel_leftward(&n_peeled_indels,&protectedp,&peeled_path,*path,&querydp5_distalgap,&genomedp5_distalgap,
			  maxpeelback,/*stop_at_indels_p*/true);
  }

  queryjump = querydp3 - querydp5_distalgap + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

  /* genomedp3 = genomedp5_distalgap + genomejump - 1; */
#ifdef EXTRACT_GENOMICSEG
  /* Make sure we don't go past the end */
  if (genomedp3 > genomiclength - 1) {
    genomedp3 = genomiclength - 1;
    genomejump = genomedp3 - genomedp5_distalgap + 1;
  }
#endif

  debug(printf("Stage 3 (dir %d), extend_ending3: Dynamic programming at 3' end (distal to gap): querydp5 = %d, querydp3 = %d, genomedp5 = %d\n",
	       cdna_direction,querydp5_distalgap,querydp3,genomedp5_distalgap));
  
  if (endalign == QUERYEND_GAP && splicesites != NULL) {
    continuous_gappairs_distalgap = Dynprog_end3_known(&(*knownsplicep),&(*dynprogindex_minor),&(*finalscore),
						       &(*ambig_end_length),&(*ambig_splicetype),
						       &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						       &(queryseq_ptr[querydp5_distalgap]),&(queryuc_ptr[querydp5_distalgap]),
						       queryjump,genomejump,querydp5_distalgap,genomedp5_distalgap,
						       querylength,chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
						       cdna_direction,watsonp,genestrand,jump_late_p,pairpool,
						       extraband_end,defect_rate);
    if (*ambig_end_length > 0) {
      *ambig_prob = 2.0;
    }
  } else {
    continuous_gappairs_distalgap = Dynprog_end3_gap(&(*dynprogindex_minor),&(*finalscore),
						     &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						     &(queryseq_ptr[querydp5_distalgap]),&(queryuc_ptr[querydp5_distalgap]),
						     queryjump,genomejump,querydp5_distalgap,genomedp5_distalgap,
						     chroffset,chrhigh,watsonp,genestrand,jump_late_p,pairpool,
						     extraband_end,defect_rate,endalign,/*require_pos_score_p*/false);
    *ambig_end_length = 0;
    *ambig_prob = 0.0;
    *knownsplicep = false;
  }

  debug(printf("  finalscore: %d\n",*finalscore));
  if (continuous_gappairs_distalgap == NULL) {
    return peeled_path;
  } else {
    continuous_gappairs_distalgap = List_reverse(continuous_gappairs_distalgap);
    firstpair = List_head(continuous_gappairs_distalgap);
    if (0 && firstpair->querypos != querydp5_distalgap) {
      /* Not a good test anymore, since we are halting peelbacks at gaps */
      /* Must have an indel between the gappairs and the rest of the read */
      debug(printf("Detected indel between gappairs %d and the rest of the read %d\n",
		   firstpair->querypos,querydp5_distalgap));
      return peeled_path;
      
    } else if (forcep == true) {
      /* For example, needed for extending middle of chimeras */
      return continuous_gappairs_distalgap;
    } else if (*finalscore <= 0) {
      /* Need this branch to eliminate bad exons at ends */
      *knownsplicep = false;
      return peeled_path;
    } else {
      return continuous_gappairs_distalgap;
    }
  }
}


#ifndef GSNAP
/* Modified from trim_novel_spliceends.  Note that code for 5' end
   here is taken from the 3' end of trim_novel_spliceends, and vice
   versa */
static void
find_dual_break_spliceends (List_T path, List_T pairs, 
			    Doublelist_T *spliceprobs5, Doublelist_T *spliceprobs3,
			    Univcoordlist_T *splice_positions_5, Univcoordlist_T *splice_positions_3,
			    int cdna_direction, bool watsonp,
			    Univcoord_T chroffset, Univcoord_T chrhigh) {
  List_T p;
  int exondist, i;
  int querypos;

  Pair_T pair;
  Univcoord_T genomicpos, start_genomicpos, middle_genomicpos, end_genomicpos;
  Univcoord_T start, middle, end; /* start to middle has mismatches, while middle to end has none */
  double donor_prob, acceptor_prob;


  debug13(printf("\nEntered find_dual_break_spliceends with cdna_direction %d\n",cdna_direction));
  *splice_positions_5 = *splice_positions_3 = (Univcoordlist_T) NULL;
  *spliceprobs5 = *spliceprobs3 = (Doublelist_T) NULL;


  /* 5' intron end */
  if (path != NULL) {
    p = path;
    while (p != NULL && ((Pair_T) p->first)->gapp == true) {
      p = List_next(p);
    }

    if (p != NULL) {
      pair = (Pair_T) List_head(p);
      querypos = pair->querypos + 1; /* Because start_genomicpos = start + 1 */
      start = middle = end = pair->genomepos;
      debug13(printf("Initializing start and end to be %u\n",start));

      i = 0;
      while (i < END_SPLICESITE_SEARCH) {
	if ((p = List_next(p)) == NULL) {
	  break;
	} else if (pair->gapp == true) {
	  break;
	} else {
	  end = pair->genomepos;
	  debug13(printf("Resetting end to be %u\n",end));
	}
	pair = (Pair_T) List_head(p);
	i++;
      }

      start = middle + 5;
      querypos += 5;

      /* Find distance from end to intron, if any */
      exondist = 0;
      while (p != NULL && ((Pair_T) List_head(p))->gapp == false &&
	     exondist < END_MIN_EXONLENGTH) {
	p = List_next(p);
	exondist++;
      }
      debug13(printf("exondist is %d\n",exondist));

      if (cdna_direction > 0) {
	if (watsonp) {
	  /* splicetype5 = splicetype5_mm = DONOR; */

	  start_genomicpos = start + 1;
	  middle_genomicpos = middle + 1;
	  end_genomicpos = end + 1;

	  /* assert(start_genomicpos >= end_genomicpos); */
	  genomicpos = start_genomicpos;
	  while (genomicpos >= middle_genomicpos &&
		 genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	    donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 1 */
	    debug13(printf("5', watson, sense anti %d %u %u %f mm",querypos,chroffset+genomicpos,genomicpos,donor_prob));
	    if (donor_prob > 0.5) {
	      *splice_positions_5 = Univcoordlist_push(*splice_positions_5,genomicpos - 1);
	      *spliceprobs5 = Doublelist_push(*spliceprobs5,donor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos--;
	    querypos--;
	  }
	  while (genomicpos >= end_genomicpos &&
		 genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	    donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 1 */
	    debug13(printf("5', watson, sense anti %d %u %u %f",querypos,chroffset+genomicpos,genomicpos,donor_prob));
	    if (donor_prob > 0.5) {
	      *splice_positions_5 = Univcoordlist_push(*splice_positions_5,genomicpos - 1);
	      *spliceprobs5 = Doublelist_push(*spliceprobs5,donor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos--;
	    querypos--;
	  }
	  debug13(printf("\n"));

	} else {
	  /* splicetype5 = splicetype5_mm = ANTIDONOR; */

	  start_genomicpos = (start > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - start;
	  middle_genomicpos = (middle > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - middle;
	  end_genomicpos = (end > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - end;

	  /* assert(start_genomicpos <= end_genomicpos); */
	  genomicpos = start_genomicpos;
	  while (genomicpos <= middle_genomicpos &&
		 genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	    donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 3 */
	    debug13(printf("5', crick, sense forward %d %u %u %f mm",querypos,chroffset+genomicpos,genomicpos,donor_prob));
	    if (donor_prob > 0.5) {
	      *splice_positions_5 = Univcoordlist_push(*splice_positions_5,(chrhigh - chroffset) - genomicpos);
	      *spliceprobs5 = Doublelist_push(*spliceprobs5,donor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos++;
	    querypos--;
	  }
	  while (genomicpos <= end_genomicpos &&
		 genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	    donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 3 */
	    debug13(printf("5', crick, sense forward %d %u %u %f",querypos,chroffset+genomicpos,genomicpos,donor_prob));
	    if (donor_prob > 0.5) {
	      *splice_positions_5 = Univcoordlist_push(*splice_positions_5,(chrhigh - chroffset) - genomicpos);
	      *spliceprobs5 = Doublelist_push(*spliceprobs5,donor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos++;
	    querypos--;
	  }
	  debug13(printf("\n"));
	}

      } else if (cdna_direction < 0) {
	if (watsonp) {
	  /* splicetype5 = splicetype5_mm = ANTIACCEPTOR; */

	  start_genomicpos = start + 1;
	  middle_genomicpos = middle + 1;
	  end_genomicpos = end + 1;

	  /* assert(start_genomicpos >= end_genomicpos); */
	  genomicpos = start_genomicpos;
	  while (genomicpos >= middle_genomicpos &&
		 genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	    acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 5 */
	    debug13(printf("5', watson, sense forward %d %u %u %f mm",querypos,chroffset+genomicpos,genomicpos,acceptor_prob));
	    if (acceptor_prob > 0.5) {
	      *splice_positions_5 = Univcoordlist_push(*splice_positions_5,genomicpos - 1);
	      *spliceprobs5 = Doublelist_push(*spliceprobs5,acceptor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos--;
	    querypos--;
	  }
	  while (genomicpos >= end_genomicpos &&
		 genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	    acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 5 */
	    debug13(printf("5', watson, sense forward %d %u %u %f",querypos,chroffset+genomicpos,genomicpos,acceptor_prob));
	    if (acceptor_prob > 0.5) {
	      *splice_positions_5 = Univcoordlist_push(*splice_positions_5,genomicpos - 1);
	      *spliceprobs5 = Doublelist_push(*spliceprobs5,acceptor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos--;
	    querypos--;
	  }
	  debug13(printf("\n"));

	} else {
	  /* splicetype5 = splicetype5_mm = ACCEPTOR; */

	  start_genomicpos = (start > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - start;
	  middle_genomicpos = (middle > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - middle;
	  end_genomicpos = (end > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - end;

	  /* assert(start_genomicpos <= end_genomicpos); */
	  genomicpos = start_genomicpos;
	  while (genomicpos <= middle_genomicpos &&
		 genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	    acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 7 */
	    debug13(printf("5', crick, sense anti %d %u %u %f mm",querypos,chroffset+genomicpos,genomicpos,acceptor_prob));
	    if (acceptor_prob > 0.5) {
	      *splice_positions_5 = Univcoordlist_push(*splice_positions_5,(chrhigh - chroffset) - genomicpos);
	      *spliceprobs5 = Doublelist_push(*spliceprobs5,acceptor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos++;
	    querypos--;
	  }
	  while (genomicpos <= end_genomicpos &&
		 genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	    acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 7 */
	    debug13(printf("5', crick, sense anti %d %u %u %f",querypos,chroffset+genomicpos,genomicpos,acceptor_prob));
	    if (acceptor_prob > 0.5) {
	      *splice_positions_5 = Univcoordlist_push(*splice_positions_5,(chrhigh - chroffset) - genomicpos);
	      *spliceprobs5 = Doublelist_push(*spliceprobs5,acceptor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos++;
	    querypos--;
	  }
	  debug13(printf("\n"));
	}
      
      } else {
	fprintf(stderr,"Not expecting cdna_direction to be 0\n");
	abort();
      }
    }
  }

  /* 3' intron end */
  if (pairs != NULL) {
    p = pairs;
    while (p != NULL && ((Pair_T) p->first)->gapp == true) {
      p = List_next(p);
    }

    if (p != NULL) {
      pair = (Pair_T) List_head(p);
      querypos = pair->querypos;
      start = middle = end = pair->genomepos;
      debug13(printf("Initializing start and end to be %u\n",start));

      i = 0;
      while (i < END_SPLICESITE_SEARCH) {
	if ((p = List_next(p)) == NULL) {
	  break;
	} else if (pair->gapp == true) {
	  break;
	} else {
	  end = pair->genomepos;
	  debug13(printf("Resetting end to be %u\n",end));
	}
	pair = (Pair_T) List_head(p);
	i++;
      }

      start = middle - 5;
      querypos -= 5;

      /* Find distance from end to intron, if any */
      exondist = 0;
      while (p != NULL && ((Pair_T) List_head(p))->gapp == false &&
	     exondist < END_MIN_EXONLENGTH) {
	p = List_next(p);
	exondist++;
      }
      debug13(printf("exondist is %d\n",exondist));

      if (cdna_direction > 0) {
	if (watsonp) {
	  /* splicetype3 = splicetype3_mm = ACCEPTOR; */

	  start_genomicpos = start;
	  middle_genomicpos = middle;
	  end_genomicpos = end;

	  /* assert(start_genomicpos <= end_genomicpos); */
	  genomicpos = start_genomicpos;
	  while (genomicpos <= middle_genomicpos &&
		 genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	    acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 2 */
	    debug13(printf("3', watson, sense forward %d %u %u %f mm",querypos,chroffset+genomicpos,genomicpos,acceptor_prob));
	    if (acceptor_prob > 0.5) {
	      *splice_positions_3 = Univcoordlist_push(*splice_positions_3,genomicpos);
	      *spliceprobs3 = Doublelist_push(*spliceprobs3,acceptor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos++;
	    querypos++;
	  }
	  while (genomicpos <= end_genomicpos &&
		 genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	    acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 2 */
	    debug13(printf("3', watson, sense forward %d %u %u %f",querypos,chroffset+genomicpos,genomicpos,acceptor_prob));
	    if (acceptor_prob > 0.5) {
	      *splice_positions_3 = Univcoordlist_push(*splice_positions_3,genomicpos);
	      *spliceprobs3 = Doublelist_push(*spliceprobs3,acceptor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos++;
	    querypos++;
	  }
	  debug13(printf("\n"));

	} else {
	  /* splicetype3 = splicetype3_mm = ANTIACCEPTOR; */

	  start_genomicpos = (start > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - start + 1;
	  middle_genomicpos = (middle > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - middle + 1;
	  end_genomicpos = (end > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - end + 1;

	  /* assert(start_genomicpos >= end_genomicpos); */
	  genomicpos = start_genomicpos;
	  while (genomicpos >= middle_genomicpos &&
		 genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	    acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 4 */
	    debug13(printf("3', crick, sense anti %d %u %u %f mm",querypos,chroffset+genomicpos,genomicpos,acceptor_prob));
	    if (acceptor_prob > 0.5) {
	      *splice_positions_3 = Univcoordlist_push(*splice_positions_3,(chrhigh - chroffset) - genomicpos + 1);
	      *spliceprobs3 = Doublelist_push(*spliceprobs3,acceptor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos--;
	    querypos++;
	  }
	  while (genomicpos >= end_genomicpos &&
		 genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	    acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 4 */
	    debug13(printf("3', crick, sense anti %d %u %u %f",querypos,chroffset+genomicpos,genomicpos,acceptor_prob));
	    if (acceptor_prob > 0.5) {
	      *splice_positions_3 = Univcoordlist_push(*splice_positions_3,(chrhigh - chroffset) - genomicpos + 1);
	      *spliceprobs3 = Doublelist_push(*spliceprobs3,acceptor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos--;
	    querypos++;
	  }
	  debug13(printf("\n"));
	}

      } else if (cdna_direction < 0) {
	if (watsonp) {
	  /* splicetype3 = splicetype3_mm = ANTIDONOR; */

	  start_genomicpos = start;
	  middle_genomicpos = middle;
	  end_genomicpos = end;
	
	  /* assert(start_genomicpos <= end_genomicpos); */
	  genomicpos = start_genomicpos;
	  while (genomicpos <= middle_genomicpos &&
		 genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	    donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 6 */
	    debug13(printf("3', watson, sense anti %d %u %u %f mm",querypos,chroffset+genomicpos,genomicpos,donor_prob));
	    if (donor_prob > 0.5) {
	      *splice_positions_3 = Univcoordlist_push(*splice_positions_3,genomicpos);
	      *spliceprobs3 = Doublelist_push(*spliceprobs3,donor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos++;
	    querypos++;
	  }
	  while (genomicpos <= end_genomicpos &&
		 genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	    donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 6 */
	    debug13(printf("3', watson, sense anti %d %u %u %f",querypos,chroffset+genomicpos,genomicpos,donor_prob));
	    if (donor_prob > 0.5) {
	      *splice_positions_3 = Univcoordlist_push(*splice_positions_3,genomicpos);
	      *spliceprobs3 = Doublelist_push(*spliceprobs3,donor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos++;
	    querypos++;
	  }
	  debug13(printf("\n"));

	} else {
	  /* splicetype3 = splicetype3_mm = DONOR; */

	  start_genomicpos = (start > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - start + 1;
	  middle_genomicpos = (middle > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - middle + 1;
	  end_genomicpos = (end > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - end + 1;

	  /* assert(start_genomicpos >= end_genomicpos); */
	  genomicpos = start_genomicpos;
	  while (genomicpos >= middle_genomicpos &&
		 genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	    donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 8 */
	    debug13(printf("3', crick, sense forward %d %u %u %f mm",querypos,chroffset+genomicpos,genomicpos,donor_prob));
	    if (donor_prob > 0.5) {
	      *splice_positions_3 = Univcoordlist_push(*splice_positions_3,(chrhigh - chroffset) - genomicpos + 1);
	      *spliceprobs3 = Doublelist_push(*spliceprobs3,donor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos--;
	    querypos++;
	  }
	  while (genomicpos >= end_genomicpos &&
		 genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	    donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 8 */
	    debug13(printf("3', crick, sense forward %d %u %u %f",querypos,chroffset+genomicpos,genomicpos,donor_prob));
	    if (donor_prob > 0.5) {
	      *splice_positions_3 = Univcoordlist_push(*splice_positions_3,(chrhigh - chroffset) - genomicpos + 1);
	      *spliceprobs3 = Doublelist_push(*spliceprobs3,donor_prob);
	      debug13(printf(" **"));
	    }
	    debug13(printf("\n"));
	    genomicpos--;
	    querypos++;
	  }
	  debug13(printf("\n"));
	}
      
      } else {
	fprintf(stderr,"Not expecting cdna_direction to be 0\n");
	abort();
      }
    }
  }

  return;
}
#endif


#ifndef GSNAP
/* Modified from make_microexon_pairs_double in dynprog_single.c */
static List_T
add_microexon_pairs (List_T pairs, int querydpM, Chrpos_T genomedpM, int lengthM,
		     char *queryseq_ptr, char *queryuc_ptr,
		     Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
		     int genestrand, Pairpool_T pairpool, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int i;

  /* Microexon */
  for (i = 0; i < lengthM; i++) {
    c1 = queryseq_ptr[querydpM+i];
    c1_uc = queryuc_ptr[querydpM+i];

    c2 = get_genomic_nt(&c2_alt,genomedpM+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[genomedpM+i]);
#endif

    if (c1_uc == c2 || c1_uc == c2_alt) {
      pairs = Pairpool_push(pairs,pairpool,querydpM+i,genomedpM+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (Dynprog_consistent_p(c1_uc,/*g*/c2,/*g_alt*/c2_alt,genestrand) == true) {
      pairs = Pairpool_push(pairs,pairpool,querydpM+i,genomedpM+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,querydpM+i,genomedpM+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  return pairs;
}
#endif


#ifndef GSNAP
static List_T
traverse_dual_break (List_T pairs, List_T *path, Pair_T leftpair, Pair_T rightpair,
		     Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, int querylength,
		     int cdna_direction, bool watsonp, int genestrand, Pairpool_T pairpool, int maxpeelback,
		     Oligoindex_array_T oligoindices_minor,
		     Diagpool_T diagpool, Cellpool_T cellpool, int *dynprogindex) {
  List_T gappairs, peeled_pairs = NULL, peeled_path = NULL, q, r;
  int querydp5, querydp3;
  Chrpos_T genomedp5, genomedp3;
  Univcoord_T best_splicepos5, best_splicepos3, best_splicepos5_with_canonical, best_splicepos3_with_canonical,
    splicepos5, splicepos3;
  bool protectedp;
  int n_peeled_indels;
  Pair_T firstpair, lastpair;
  Chrpos_T chrstart, chrend;

  Univcoordlist_T splice_positions_5, splice_positions_3, a, b;
  Doublelist_T spliceprobs5, spliceprobs3, x, y;
  Intlist_T hits = NULL, p;
  int best_middlelength, best_middlelength_with_canonical, middlelength;
  int best_candidate, best_candidate_with_canonical, candidate;
  Chrpos_T splicesitepos;
  double bestprob = 0.0, bestprob_with_canonical = 0.0, prob1, prob2, prob3, prob4;
  char intron1, intron2, intron3, intron4;
  char c1_alt, c2_alt, c3_alt, c4_alt;
  int nloops;

  debug14(printf("\nTRAVERSE_DUAL_BREAK\n"));

  /* First, try to find a microexon */
  find_dual_break_spliceends(*path,pairs,&spliceprobs5,&spliceprobs3,&splice_positions_5,&splice_positions_3,
			     cdna_direction,watsonp,chroffset,chrhigh);

  if (cdna_direction > 0) {
    intron1 = 'G'; intron2 = 'T'; intron3 = 'A'; intron4 = 'G';
  } else if (cdna_direction < 0) {
    intron1 = 'C'; intron2 = 'T'; intron3 = 'A'; intron4 = 'C';
  }


  debug(printf("Number of 5' splice positions: %d\n",Univcoordlist_length(splice_positions_5)));
  debug(printf("Number of 3' splice positions: %d\n",Univcoordlist_length(splice_positions_3)));
  nloops = Univcoordlist_length(splice_positions_5) * Univcoordlist_length(splice_positions_3);

  if (nloops < MAX_MICROEXON_LOOPS) {
    for (a = splice_positions_5, x = spliceprobs5; a != NULL;
	 a = Univcoordlist_next(a), x = Doublelist_next(x)) {
      prob1 = Doublelist_head(x);
      splicepos5 = Univcoordlist_head(a);

      q = *path;
      while (q != NULL && ((Pair_T) q->first)->genomepos > splicepos5) {
	q = List_next(q);
      }
      if (q == NULL) {
	leftpair = (Pair_T) NULL;
      } else {
	leftpair = (Pair_T) q->first;
	querydp5 = leftpair->querypos + 1;
	genomedp5 = leftpair->genomepos + 1;
      }
      
      for (b = splice_positions_3, y = spliceprobs3; b != NULL;
	   b = Univcoordlist_next(b), y = Doublelist_next(y)) {
	prob4 = Doublelist_head(y);
	splicepos3 = Univcoordlist_head(b);
	
	r = pairs;
	while (r != NULL && ((Pair_T) r->first)->genomepos < splicepos3) {
	  r = List_next(r);
	}
	if (r == NULL) {
	  rightpair = (Pair_T) NULL;
	} else {
	  rightpair = (Pair_T) r->first;
	  querydp3 = rightpair->querypos - 1;
	  genomedp3 = rightpair->genomepos - 1;
	}
	
	if (leftpair != NULL && rightpair != NULL &&
	    (middlelength = querydp3 - querydp5 + 1) > MIN_MICROEXON_LENGTH) {
	  debug(printf("middlelength %d = %d - %d + 1\n",middlelength,querydp3,querydp5));
	  hits = BoyerMoore_nt(&(queryuc_ptr[querydp5]),/*querylen*/middlelength,
			       /*textleft*/genomedp5,/*textlen*/genomedp3 - genomedp5 + 1,
			       chroffset,chrhigh,watsonp);
	  for (p = hits; p != NULL; p = Intlist_next(p)) {
	    candidate = genomedp5 + Intlist_head(p);
	    
	    /* Not handling known splice sites yet */
	    if (watsonp == true) {
	      if (cdna_direction > 0) {
		splicesitepos = chroffset + (candidate-1) + 1;
		prob2 = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
		splicesitepos = chroffset + candidate+middlelength;
		prob3 = Maxent_hr_donor_prob(splicesitepos,chroffset);
	      } else {
		splicesitepos = chroffset + (candidate-1) + 1;
		prob2 = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
		splicesitepos = chroffset + candidate+middlelength;
		prob3 = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	      }
	    } else {
	      if (cdna_direction > 0) {
		splicesitepos = chrhigh - (candidate-1);
		prob2 = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
		splicesitepos = chrhigh - (candidate+middlelength) + 1;
		prob3 = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	      } else {
		splicesitepos = chrhigh - (candidate-1);
		prob2 = Maxent_hr_donor_prob(splicesitepos,chroffset);
		splicesitepos = chrhigh - (candidate+middlelength) + 1;
		prob3 = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	      }
	    }
	    
	    debug13(printf("end probabilities: prob1 = %f, prob4 = %f, microexon probabilities: prob2 = %f, prob3 = %f\n",prob1,prob4,prob2,prob3));
	    if (prob1 + prob2 + prob3 + prob4 > bestprob) {
	      best_splicepos5 = splicepos5;
	      best_candidate = candidate;
	      best_middlelength = middlelength;
	      best_splicepos3 = splicepos3;
	      bestprob = prob1 + prob2 + prob3 + prob4;
	    }
	    
	    debug(printf("candidate: at %u\n",candidate));
	    debug(printf("intron3 %c\n",get_genomic_nt(&c3_alt,candidate-2,chroffset,chrhigh,watsonp)));
	    debug(printf("intron4 %c\n",get_genomic_nt(&c4_alt,candidate-1,chroffset,chrhigh,watsonp)));
	    debug(printf("intron1 %c\n",get_genomic_nt(&c1_alt,candidate+middlelength,chroffset,chrhigh,watsonp)));
	    debug(printf("intron2 %c\n",get_genomic_nt(&c2_alt,candidate+middlelength+1,chroffset,chrhigh,watsonp)));
	    
	    if (/*genomicuc[candidate - 2]*/ get_genomic_nt(&c3_alt,candidate-2,chroffset,chrhigh,watsonp) == intron3 &&
		/*genomicuc[candidate - 1]*/ get_genomic_nt(&c4_alt,candidate-1,chroffset,chrhigh,watsonp)  == intron4 &&
		/*genomicuc[candidate + middlelength]*/ get_genomic_nt(&c1_alt,candidate+middlelength,chroffset,chrhigh,watsonp) == intron1 &&
		/*genomicuc[candidate + middlelength + 1]*/ get_genomic_nt(&c2_alt,candidate+middlelength+1,chroffset,chrhigh,watsonp) == intron2) {
	      debug(printf("  Canonical microexon at %d >>> %d..%d >>> %d\n",genomedp5,candidate,candidate+middlelength,genomedp3));
	      if (prob1 + prob2 + prob3 + prob4 > bestprob_with_canonical) {
		best_splicepos5_with_canonical = splicepos5;
		best_candidate_with_canonical = candidate;
		best_middlelength_with_canonical = middlelength;
		best_splicepos3_with_canonical = splicepos3;
		bestprob_with_canonical = prob1 + prob2 + prob3 + prob4;
	      }
	    }
	  }
	  
	  Intlist_free(&hits);
	}
      }
    }
  }

  debug13(printf("best prob is %f\n",bestprob));
  if (bestprob > 3.0) {
    while ((*path) != NULL && ((Pair_T) (*path)->first)->genomepos > best_splicepos5) {
      *path = List_next(*path);
    }
    while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < best_splicepos3) {
      pairs = List_next(pairs);
    }
    leftpair = (Pair_T) (*path)->first;
    querydp5 = leftpair->querypos + 1;

    debug13(printf("Making microexon pairs with splicepos5 %u, candidate %u, middlelength %d, splicepos3 %u\n",
		   best_splicepos5,best_candidate,best_middlelength,best_splicepos3));
    *path = add_microexon_pairs(*path,/*querydpM*/querydp5,/*genomedpM*/best_candidate,
				/*lengthM*/best_middlelength,queryseq_ptr,queryuc_ptr,
				chroffset,chrhigh,watsonp,genestrand,pairpool,*dynprogindex);
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);

    Univcoordlist_free(&splice_positions_5);
    Univcoordlist_free(&splice_positions_3);
    Doublelist_free(&spliceprobs3);
    Doublelist_free(&spliceprobs5);

    return pairs;

  } else if (bestprob_with_canonical > 0.0) {
    while ((*path) != NULL && ((Pair_T) (*path)->first)->genomepos > best_splicepos5_with_canonical) {
      *path = List_next(*path);
    }
    while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < best_splicepos3_with_canonical) {
      pairs = List_next(pairs);
    }
    leftpair = (Pair_T) (*path)->first;
    querydp5 = leftpair->querypos + 1;

    debug13(printf("Making microexon pairs with splicepos5 %u, candidate %u, middlelength %d, splicepos3 %u\n",
		   best_splicepos5_with_canonical,best_candidate_with_canonical,best_middlelength_with_canonical,best_splicepos3_with_canonical));
    *path = add_microexon_pairs(*path,/*querydpM*/querydp5,/*genomedpM*/best_candidate_with_canonical,
				/*lengthM*/best_middlelength_with_canonical,queryseq_ptr,queryuc_ptr,
				chroffset,chrhigh,watsonp,genestrand,pairpool,*dynprogindex);
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);

    Univcoordlist_free(&splice_positions_5);
    Univcoordlist_free(&splice_positions_3);
    Doublelist_free(&spliceprobs3);
    Doublelist_free(&spliceprobs5);

    return pairs;

  } else {
    Univcoordlist_free(&splice_positions_5);
    Univcoordlist_free(&splice_positions_3);
    Doublelist_free(&spliceprobs3);
    Doublelist_free(&spliceprobs5);
  }



  /* Try to solve without a microexon */
  if (*path == NULL) {
    leftpair = (Pair_T) NULL;
  } else {
    leftpair = (Pair_T) (*path)->first;
  }

  if (pairs == NULL) {
    rightpair = (Pair_T) NULL;
  } else {
    rightpair = (Pair_T) pairs->first;
  }

  if (leftpair != NULL && rightpair != NULL) {
    querydp5 = leftpair->querypos + 1;
    genomedp5 = leftpair->genomepos + 1;
    /* if (leftpair->cdna == ' ') querydp5--; -- For old dynamic programming */
    /* if (leftpair->genome == ' ') genomedp5--; -- For old dynamic programming */
    querydp3 = rightpair->querypos - 1;
    genomedp3 = rightpair->genomepos - 1;
  } else if (leftpair == NULL) {
    querydp5 = 0;
    genomedp5 = rightpair->genomepos - rightpair->querypos - 100;
    querydp3 = rightpair->querypos - 1;
    genomedp3 = rightpair->genomepos - 1;
  } else if (rightpair == NULL) {
    querydp5 = leftpair->querypos + 1;
    genomedp5 = leftpair->genomepos + 1;
    /* if (leftpair->cdna == ' ') querydp5--; -- For old dynamic programming */
    /* if (leftpair->genome == ' ') genomedp5--; -- For old dynamic programming */
    querydp3 = querylength - 1;
    genomedp3 = leftpair->genomepos + (querylength - leftpair->querypos) + 100;
  }

  /* Previously skipped this, but need to do at least a little
     peelback to avoid gaps at either end */
  protectedp = false;
  pairs = peel_rightward(&n_peeled_indels,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
			 maxpeelback,/*stop_at_indels_p*/true);
  *path = peel_leftward(&n_peeled_indels,&protectedp,&peeled_path,*path,&querydp5,&genomedp5,
			maxpeelback,/*stop_at_indels_p*/true);
#ifdef PMAP
  querydp3 /= 3;
  querydp5 /= 3;
#endif

  debug14(printf("genome %d..%d, query %d..%d\n",genomedp5,genomedp3,querydp5,querydp3));
#ifdef EXTRACT_GENOMICSEG
  debug14(printf("genome %.*s\n",genomedp3-genomedp5+1,&(genomicseg_ptr[genomedp5])));
#endif
#ifdef PMAP
  debug14(printf("query %.*s\n",querydp3-querydp5+1,&(queryaaseq_ptr[querydp5])));
#else
  debug14(printf("query %.*s\n",querydp3-querydp5+1,&(queryseq_ptr[querydp5])));
#endif

  if (watsonp) {
    chrstart = genomedp5;
    chrend = genomedp3;
  } else {
    chrstart = (chrhigh - chroffset) - genomedp3;
    chrend = (chrhigh - chroffset) - genomedp5;
  }

  /* chrend can be equal to chrstart - 1 if peelback fails on both ends */
  assert(chrend + 1 >= chrstart);
  debug14(printf("Starting stage2 with chrstart %u, chrend %u, watsonp %d\n",
		 chrstart,chrend,watsonp));
  gappairs = Stage2_compute_one(
#ifdef PMAP
				  &(queryaaseq_ptr[querydp5]),&(queryaaseq_ptr[querydp5]),
				  /*querylength*/querydp3-querydp5+1,/*query_offset*/querydp5*3,
#else
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  /*querylength*/querydp3-querydp5+1,/*query_offset*/querydp5,
#endif
				  chrstart,chrend,chroffset,chrhigh,/*plusp*/watsonp,genestrand,

				  oligoindices_minor,pairpool,diagpool,cellpool,
				  /*localp should be false*/true,/*skip_repetitive_p*/false,
				  /*favor_right_p*/false,/*debug_graphic_p*/false);
    
  debug14(printf("Internal stage2 result:\n"));
  debug14(Pair_dump_list(gappairs,true));

  if (gappairs == NULL) {
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    if (*path != NULL && pairs != NULL) {
      /* Do not put a gap at the end of the alignment */
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,
				      /*leftpair*/(*path)->first,/*rightpair*/pairs->first,/*knownp*/false);
    }

  } else {
    lastpair = (Pair_T) gappairs->first;
    firstpair = (Pair_T) List_last_value(gappairs);
    debug14(printf("gappairs goes from %d to %d\n",firstpair->querypos,lastpair->querypos));
    if (1 || (firstpair->querypos == querydp5 && lastpair->querypos == querydp3)) {
      /* Note: We have to take this branch, otherwise we get unexpected comp errors */
      /* fprintf(stderr,"%d..%d .. %d..%d\n",querydp5,firstpair->querypos,lastpair->querypos,querydp3); */
      debug14(printf("  => entire query sequence bridged or not, but taking it regardless\n"));
      pairs = Pairpool_transfer(pairs,gappairs);

    } else {
      debug14(printf("  => entire query sequence not bridged, so abort\n"));
      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);
      if (*path != NULL && pairs != NULL) {
	/* Do not put a gap at the end of the alignment */
	pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,
					/*leftpair*/(*path)->first,/*rightpair*/pairs->first,/*knownp*/false);
      }
    }
  }

  return pairs;
}
#endif


/************************************************************************
 *   End of traversal functions
 ************************************************************************/

static List_T
build_dual_breaks (bool *dual_break_p, int *dynprogindex_minor, int *dynprogindex_major, List_T path,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		   char *queryaaseq_ptr,
#endif
		   char *queryseq_ptr, char *queryuc_ptr, int querylength,
		   int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		   Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3, int maxpeelback,
#ifndef GSNAP
		   Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		   double defect_rate, bool finalp, bool simplep) {

  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;
  bool filledp, shiftp;

  debug(printf("Entered build_dual_breaks\n"));
  *dual_break_p = false;
  debug(Pair_dump_list(path,true));

  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;

    /* Cannot rely on previous procedures to assign pair->comp value */
    if (pair->gapp == false /*|| pair->comp != DUALBREAK_COMP*/) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif
    } else if (path->rest == NULL || pairs == NULL) {
      debug(printf("Observed a gap at the end of the alignment, case 1\n"));
      path = Pairpool_pop(path,&pair);

    } else {
      /* pairptr = path; -- save */
      path = Pairpool_pop(path,&pair);

      leftpair = path->first;
      rightpair = pairs->first;
      if (leftpair->querypos < 0 || rightpair->querypos < 0) {
	debug(printf("Observed a gap at the end of the alignment, case 2\n"));
      } else {
	debug(printf("Observed a gap at %d..%d with queryjump = %d, genomejump = %d\n",
		     leftpair->querypos,rightpair->querypos,pair->queryjump,pair->genomejump));
	
	if (0 && finalp == true) {
	  /* If genomejump is too large, this causes problem with allocation in Dynprog_T objects */
	  debug(printf("  Final: solve as a single gap\n"));
	  pairs = traverse_single_gap(&filledp,&(*dynprogindex_minor),pairs,&path,leftpair,rightpair,
				      chroffset,chrhigh,
				      queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
				      jump_late_p,pairpool,dynprogM,last_genomedp5,last_genomedp3,
				      maxpeelback,defect_rate,/*forcep*/true,/*finalp*/false);

	} else if (pair->genomejump - pair->queryjump < SINGLESLEN &&
		   pair->queryjump - pair->genomejump < SINGLESLEN) {
	  debug(printf("  Can be solved as a single gap\n"));
	  pairs = traverse_single_gap(&filledp,&(*dynprogindex_minor),pairs,&path,leftpair,rightpair,
				      chroffset,chrhigh,
				      queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
				      jump_late_p,pairpool,dynprogM,last_genomedp5,last_genomedp3,
				      maxpeelback,defect_rate,/*forcep*/true,/*finalp*/false);

#ifdef GSNAP
	} else {
#else
	} else if (pair->queryjump < MIN_STAGE2_FOR_DUALBREAK) {
#endif
	  debug(printf("  Too small for a dual break\n"));
	  pairs = traverse_genome_gap(&filledp,&shiftp,&(*dynprogindex_minor),&(*dynprogindex_major),
				      pairs,&path,leftpair,rightpair,chrnum,chroffset,chrhigh,
				      queryseq_ptr,queryuc_ptr,querylength,
				      cdna_direction,watsonp,genestrand,jump_late_p,
				      pairpool,dynprogL,dynprogM,dynprogR,last_genomedp5,last_genomedp3,
				      maxpeelback,defect_rate,/*finalp*/false,simplep);

#ifndef GSNAP
	} else {
	  debug(printf("  Solving as a dual break\n"));
	  *dual_break_p = true;
	  pairs = traverse_dual_break(pairs,&path,leftpair,rightpair,chroffset,chrhigh,
#ifdef PMAP
				      queryaaseq_ptr,
#endif
				      queryseq_ptr,queryuc_ptr,querylength,cdna_direction,watsonp,genestrand,
				      pairpool,maxpeelback,oligoindices_minor,
				      diagpool,cellpool,&(*dynprogindex_major));
#endif
	}
      }
    }
  }

  debug(printf("After build_dual_breaks:\n"));
  debug(Pair_dump_list(pairs,true));

  return pairs;
}


  
/* Note: querypos is actually indexsize nt to the left of the last nt match.
   
	||||||||********   X  X	 XX X	X
	       ^	 <- queryjump->	 ^
	    leftquerypos		 rightquerypos

		<-     querydpspan     ->
*/
static List_T
build_path_end3 (bool *knownsplicep, int *ambig_end_length_3, Splicetype_T *ambig_splicetype_3, double *ambig_prob_3,
		 bool *chop_exon_p, int *dynprogindex_minor,
		 List_T path, Univcoord_T chroffset, Univcoord_T chrhigh, int querylength,
		 Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		 char *queryseq_ptr, char *queryuc_ptr,
		 int cdna_direction, bool watsonp, int genestrand, bool jump_late_p, int maxpeelback,
		 double defect_rate, Pairpool_T pairpool, Dynprog_T dynprogL,
		 bool extendp, Endalign_T endalign, bool forcep) {
  List_T gappairs;
  Pair_T leftpair;
  /* int genomejump */
  int queryjump, rightquerypos;
  int finalscore;

#if 0
  if (*ambig_end_length_3 > 0) {
    debug(printf("ambig_end_length_3 is %d, so returning path\n",*ambig_end_length_3));
    return path;
  } else {
    path = clean_path_end3_gap_indels(path);
  }
#else
  /* Always want to clean indels at end */
  path = clean_path_end3_gap_indels(path);
  if (*ambig_end_length_3 > 0) {
    debug(printf("ambig_end_length_3 is %d, so returning path\n",*ambig_end_length_3));
    return path;
  }
#endif


  *knownsplicep = false;
  if (path == NULL) {
    *ambig_end_length_3 = 0;
    *ambig_prob_3 = 0.0;
    return (List_T) NULL;
  } else {
    leftpair = path->first;
  }
  debug(printf("Stage 3 (dir %d): 3' end: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d\n",
	       cdna_direction,leftpair->querypos,querylength,leftpair->genomepos));
  if (leftpair->querypos < 0) {
    *ambig_end_length_3 = 0;
    *ambig_prob_3 = 0.0;
    return (List_T) NULL;
    /* abort(); */
  }

  queryjump = querylength - leftpair->querypos - 1;
  /* genomejump = genomiclength - leftpair->genomepos - 1; */
  /* if (leftpair->cdna == ' ') queryjump++; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') genomejump++; */

  /* Note difference with 5' case.  We use queryjump+1 here instead of queryjump and genomejump */
  /* Do use nullgap here.  Truncating back to entire exon can slow algorithm down significantly. */
  if (/* 0 && */ queryjump+1 > nullgap) {
    rightquerypos = leftpair->querypos + nullgap + 1;
    debug(printf("Since queryjump+1 %d > nullgap %d, setting rightquerypos %d = %d + %d + 1\n",
		 queryjump+1,nullgap,rightquerypos,leftpair->querypos,nullgap));
  } else {
    rightquerypos = querylength;
  }

  if (extendp == true) {
    debug(printf("Running extend_ending3\n"));
    *chop_exon_p = false;
    gappairs = extend_ending3(&(*knownsplicep),&(*dynprogindex_minor),&finalscore,
			      &(*ambig_end_length_3),&(*ambig_splicetype_3),&(*ambig_prob_3),
			      &path,leftpair,rightquerypos,querylength,
			      chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
			      queryseq_ptr,queryuc_ptr,
			      cdna_direction,watsonp,genestrand,jump_late_p,pairpool,dynprogL,maxpeelback,
			      defect_rate,endalign,forcep);
  } else {
    /* Looks like we aren't calling this anymore */
    abort();
    debug(printf("Running distalmedial_ending3\n"));
    gappairs = distalmedial_ending3(&(*knownsplicep),&(*chop_exon_p),&(*dynprogindex_minor),
				    &finalscore,&(*ambig_end_length_3),&(*ambig_prob_3),
				    &path,leftpair,rightquerypos,chroffset,chrhigh,
				    queryseq_ptr,queryuc_ptr,
				    watsonp,genestrand,jump_late_p,pairpool,dynprogL,
				    maxpeelback,defect_rate);
  }

  debug(printf("Gappairs from build_path_end3:\n"));
  debug(Pair_dump_list(gappairs,true));

  path = Pairpool_transfer(path,gappairs);

  debug(printf("Final result of build_path_end3:\n"));
  debug(Pair_dump_list(path,true));
  debug(printf("\n"));

  debug(printf("ambig_end_length_3 is %d\n",*ambig_end_length_3));

  return path;
}


/* Schematic:
   
       <- queryjump ->
       X  X  XX X   X ********||||||||
      ^		      ^
   leftquerypos	      rightquerypos

       <-    querydpspan    ->

*/
static List_T
build_pairs_end5 (bool *knownsplicep, int *ambig_end_length_5, Splicetype_T *ambig_splicetype_5, double *ambig_prob_5,
		  bool *chop_exon_p, int *dynprogindex_minor, List_T pairs,
		  Univcoord_T chroffset, Univcoord_T chrhigh,
		  Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		  char *queryseq_ptr, char *queryuc_ptr,
		  int cdna_direction, bool watsonp, int genestrand, bool jump_late_p, int maxpeelback,
		  double defect_rate, Pairpool_T pairpool, Dynprog_T dynprogR,
		  bool extendp, Endalign_T endalign, bool forcep) {
  List_T gappairs;
  Pair_T rightpair;
  int queryjump, leftquerypos;
  int finalscore;
  /* int genomejump */

#if 0
  if (*ambig_end_length_5 > 0) {
    debug(printf("ambig_end_length_5 is %d, so returning pairs\n",*ambig_end_length_5));
    return pairs;
  } else {
    pairs = clean_pairs_end5_gap_indels(pairs);
  }
#else
  /* Always want to clean indels at end */
  pairs = clean_pairs_end5_gap_indels(pairs);
  if (*ambig_end_length_5 > 0) {
    debug(printf("ambig_end_length_5 is %d, so returning pairs\n",*ambig_end_length_5));
    return pairs;
  }
#endif


  *knownsplicep = false;
  if (pairs == NULL) {
    *ambig_end_length_5 = 0;
    *ambig_prob_5 = 0.0;
    return (List_T) NULL;
  } else {
    rightpair = pairs->first;
  }
  debug(printf("Stage 3 (dir %d): 5' end: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d\n",
	       cdna_direction,-1,rightpair->querypos,-1));
  if (rightpair->querypos < 0) {
    *ambig_end_length_5 = 0;
    *ambig_prob_5 = 0.0;
    return (List_T) NULL;
    /* abort(); */
  }

  queryjump = rightpair->querypos; /* - leftquerypos (-1) - 1 */
  /* genomejump = rightpair->genomepos; */  /* - leftgenomepos (-1) - 1 */

  /* Note difference with 3' case.  We use queryjump here instead of queryjump+1 */
  /* Do use nullgap here.  Truncating back to entire exon can slow algorithm significantly. */
  if (/*0 && */ queryjump > nullgap) {
    leftquerypos = rightpair->querypos - nullgap - 1;
  } else {
    leftquerypos = -1;
  }

  if (extendp == true) {
    debug(printf("Running extend_ending5\n"));
    *chop_exon_p = false;
    gappairs = extend_ending5(&(*knownsplicep),&(*dynprogindex_minor),
			      &finalscore,&(*ambig_end_length_5),&(*ambig_splicetype_5),&(*ambig_prob_5),
			      &pairs,leftquerypos,rightpair,
			      chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
			      queryseq_ptr,queryuc_ptr,
			      cdna_direction,watsonp,genestrand,jump_late_p,pairpool,dynprogR,
			      maxpeelback,defect_rate,endalign,forcep);
  } else {
    /* Looks like we aren't calling this anymore */
    abort();
    debug(printf("Running distalmedial_ending5\n"));
    gappairs = distalmedial_ending5(&(*knownsplicep),&(*chop_exon_p),&(*dynprogindex_minor),
				    &finalscore,&(*ambig_end_length_5),&(*ambig_prob_5),
				    &pairs,leftquerypos,rightpair,chroffset,chrhigh,
				    queryseq_ptr,queryuc_ptr,
				    watsonp,genestrand,jump_late_p,pairpool,dynprogR,
				    maxpeelback,defect_rate);
  }

  debug(printf("Gappairs from build_pairs_end5:\n"));
  debug(Pair_dump_list(gappairs,true));

  pairs = Pairpool_transfer(pairs,gappairs);

  debug(printf("Final result of build_pairs_end5:\n"));
  debug(Pair_dump_list(pairs,true));
  debug(printf("\n"));

  debug(printf("ambig_end_length_5 is %d\n",*ambig_end_length_5));

  return pairs;
}


/* maxsize can be either 3 or nullgap */
static List_T
build_pairs_singles (int *dynprogindex, List_T path, int maxsize,
		     Univcoord_T chroffset, Univcoord_T chrhigh,
		     char *queryseq_ptr, char *queryuc_ptr, int querylength,
		     bool watsonp, int genestrand, bool jump_late_p, int maxpeelback, double defect_rate,
		     Pairpool_T pairpool, Dynprog_T dynprogM, 
		     Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3, bool forcep, bool finalp) {
  List_T pairs = NULL, pairptr;
  Pair_T pair, leftpair, rightpair;
  bool filledp;

  debug(printf("\n** Starting build_pairs_singles with maxsize %d\n",maxsize));

  /* Remove gaps at beginning */
  while (path != NULL && ((Pair_T) path->first)->gapp == true) {
    path = Pairpool_pop(path,&pair);
  }

  while (path != NULL && path->rest != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;
    if (pair->gapp == false) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->queryjump > maxsize) {
      /* Large gap.  Do nothing */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      /* cDNA insertion.  Do nothing */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->genomejump > pair->queryjump + SINGLESLEN) {
      /* Intron.  Do nothing */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (path->rest == NULL || pairs == NULL) {
      fprintf(stderr,"Single gap at beginning or end of alignment\n");
      abort();

    } else {
      /* Guarantees: queryjump <= nullgap && genomejump < queryjump - EXTRAQUERYGAP &&
	 genomejump <= queryjump + MININTRONLEN, meaning that score matrix is nearly square */
      pairptr = path;		/* save */
      path = Pairpool_pop(path,&pair);

      leftpair = path->first;
      rightpair = pairs->first;
	
      debug(printf("Stage 3: Traversing single gap small: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d, queryjump %d, genomejump %d\n",
		   leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos,pair->queryjump,pair->genomejump));
      pairs = traverse_single_gap(&filledp,&(*dynprogindex),pairs,&path,leftpair,rightpair,
				  chroffset,chrhigh,queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
				  jump_late_p,pairpool,dynprogM,last_genomedp5,last_genomedp3,
				  maxpeelback,defect_rate,forcep,finalp);
      /* (old comment:) forcep needs to be true here to avoid subsequent anomalies in building dualintrons, e.g., XM_376610.2_mRNA on 7:127885572..127888991 */
      if (filledp == true) {
	/* Discard the gap */
	debug(printf("Discarding gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
      } else {
	/* Replace the gap */
	debug(printf("Replacing gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_push_existing(pairs,pairptr);
#endif

      }
    }
  }

  /* Handle last entry if not a gap */
  if (path != NULL && ((Pair_T) path->first)->gapp == false) {
    pair = (Pair_T) path->first;
    pairs = List_transfer_one(pairs,&path);
  }

  return pairs;
}


static List_T
build_pairs_dualintrons (int *dynprogindex, List_T path,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			 char *queryseq_ptr, char *queryuc_ptr, int querylength,
			 int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
			 int maxpeelback, double defect_rate,
			 Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR,
			 Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3) {
  List_T pairs = NULL, midexon_pairs = NULL, pairptr;
  Pair_T pair, leftpair, midleftpair, midpair, midrightpair, rightpair;
  int midquerypos, midgenomepos;
  bool left_end_intron_p = false, right_end_intron_p, exonp;

  debug(printf("\n** Starting build_pairs_dualintrons\n"));
  debug(Pair_dump_list(path,true));

  /* Remove gaps at beginning */
  while (path != NULL && ((Pair_T) path->first)->gapp == true) {
    path = Pairpool_pop(path,&pair);
  }

  while (path != NULL && path->rest != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;
    
    if (pair->gapp == false) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->queryjump > nullgap) {
      /* Large gap.  Do nothing */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      /* cDNA insertion.  Do nothing */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->genomejump <= pair->queryjump + MININTRONLEN) {
      /* Single gap.  Do nothing */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else {
      pairptr = path;		/* save */
      path = Pairpool_pop(path,&pair);

      midpair = path->first;
      if (midpair->shortexonp == false) {
	/* Long exon; do nothing */
	debug(printf("I see a long exon at %d...do nothing\n",midpair->querypos));
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_push_existing(pairs,pairptr);
#endif

      } else {
	/* Short exon */
	debug(printf("I see a short exon at %d...crossing\n",midpair->querypos));
	right_end_intron_p = pair->end_intron_p;
#ifdef WASTE
	path = Pairpool_pop(path,&midpair);
	midexon_pairs = Pairpool_push_existing(NULL,pairpool,midpair);
#else
	midpair = path->first;
	midexon_pairs = List_transfer_one(NULL,&path);
#endif
	midrightpair = midpair;

	exonp = true;
	while (path != NULL && exonp) {
#ifdef WASTE
	  path = Pairpool_pop(path,&midpair);
#else
	  midpair = path->first;
#endif
	  if (midpair->gapp == true) {
	    left_end_intron_p = midpair->end_intron_p;
	    exonp = false;
#ifdef WASTE
#else
	    path = path->rest;
#endif
	  } else {
#ifdef WASTE
	    midexon_pairs = Pairpool_push_existing(midexon_pairs,pairpool,midpair);
#else
	    midexon_pairs = List_transfer_one(midexon_pairs,&path);
#endif
	  }
	}
	debug(printf("Finished crossing a short exon\n"));

	if (path == NULL) {
	  /* Short exon is the first one.  Process it. */
	  debug(printf("path is NULL so Pairpool_push_existing\n"));
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair); /* initial gap */
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif
	  pairs = Pairpool_transfer(pairs,List_reverse(midexon_pairs));

	} else {
	  /* Perform dual intron gap */
	  debug(printf("path is not NULL, so doing dual intron gap\n"));
	  midleftpair = midexon_pairs->first;
	  midgenomepos = (midleftpair->genomepos + midrightpair->genomepos)/2;
	  midquerypos = midrightpair->querypos - (midrightpair->genomepos - midgenomepos);
	  leftpair = path->first;
	  rightpair = pairs->first;
	  if (midquerypos <= leftpair->querypos || midquerypos >= rightpair->querypos) {
	    /* Skip */

	  } else {	  
	    debug(printf("Stage 3 (dir %d): Traversing dual intron gap: leftquerypos = %d, midquerypos = %d, rightquerypos = %d, leftgenomepos = %d, midgenomepos = %d, rightgenomepos = %d\n",
			 cdna_direction,leftpair->querypos,midquerypos,rightpair->querypos,
			 leftpair->genomepos,midgenomepos,rightpair->genomepos));
	  
	    pairs = traverse_dual_genome_gap(&(*dynprogindex),pairs,&path,leftpair,rightpair,
					     left_end_intron_p,right_end_intron_p,
					     chrnum,chroffset,chrhigh,midquerypos,midgenomepos,
					     queryseq_ptr,queryuc_ptr,querylength,cdna_direction,watsonp,genestrand,
					     jump_late_p,pairpool,dynprogL,dynprogR,last_genomedp5,last_genomedp3,
					     maxpeelback,defect_rate,/*finalp*/false);
	  }
	}
      }
    }
  }

  /* Handle last entry if not a gap */
  if (path != NULL && ((Pair_T) path->first)->gapp == false) {
    pair = (Pair_T) path->first;
    pairs = List_transfer_one(pairs,&path);
  }

  return pairs;
}  


static List_T
build_pairs_introns (bool *shiftp, bool *incompletep, 
		     int *dynprogindex_minor, int *dynprogindex_major, List_T path,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, int querylength,
		     int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		     int maxpeelback, double defect_rate,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
#ifndef GSNAP
		     Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		     Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3, bool finalp, bool simplep) {
  List_T pairs = NULL, pairptr;
  Pair_T pair, leftpair, rightpair;
  bool filledp;
  int minintronlen;

  debug(printf("\n** Starting build_pairs_introns\n"));
  debug(Pair_dump_list(path,true));

  if (finalp == true) {
    minintronlen = MININTRONLEN_FINAL;
  } else {
    minintronlen = MININTRONLEN;
  }

  /* Remove gaps at beginning */
  while (path != NULL && ((Pair_T) path->first)->gapp == true) {
    path = Pairpool_pop(path,&pair);
  }

  *shiftp = *incompletep = false;
  while (path != NULL && path->rest != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;
    if (pair->gapp == false) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->queryjump > nullgap) {
      if (pair->genomejump < 16) {
	/* Not enough genome material to run stage 2 */
	pairptr = path;		/* save */
	path = Pairpool_pop(path,&pair);

	leftpair = path->first;
	rightpair = pairs->first;
	debug(printf("Stage 3 (dir %d): Traversing cDNA gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d, queryjump %d, genomejump %d\n",
		     cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos,pair->queryjump,pair->genomejump));
	pairs = traverse_cdna_gap(&filledp,&(*incompletep),&(*dynprogindex_minor),&(*dynprogindex_major),
				  pairs,&path,leftpair,rightpair,
				  chroffset,chrhigh,queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
				  jump_late_p,pairpool,dynprogL,dynprogM,dynprogR,
				  last_genomedp5,last_genomedp3,maxpeelback,defect_rate,/*finalp*/true);

	if (filledp == true) {
	  /* Discard gap */
	  debug(printf("Discarding gap ");
		Pair_dump_one(pair,true);
		printf("\n"));
	} else {
	  /* Replace the gap */
	  debug(printf("Replacing gap ");
		Pair_dump_one(pair,true);
		printf("\n"));
#ifdef WASTE
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	  pairs = List_push_existing(pairs,pairptr);
#endif
	}

      } else {
#ifndef GSNAP
	/* Solve as dual break.  Should have already been done by build_dual_breaks */
	/* pairptr = path; */		/* save */
	/* path = Pairpool_pop(path,&pair); */
	leftpair = path->first;
	rightpair = pairs->first;
	pairs = traverse_dual_break(pairs,&path,leftpair,rightpair,chroffset,chrhigh,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,cdna_direction,watsonp,genestrand,
				    pairpool,maxpeelback,oligoindices_minor,
				    diagpool,cellpool,&(*dynprogindex_major));
#endif
      }

    } else if (finalp == false && pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      /* If finalp is true, then will need to solve as singles */
      pairptr = path;		/* save */
      path = Pairpool_pop(path,&pair);

      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3 (dir %d): Traversing cDNA gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d, queryjump %d, genomejump %d\n",
		   cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos,pair->queryjump,pair->genomejump));
      pairs = traverse_cdna_gap(&filledp,&(*incompletep),&(*dynprogindex_minor),&(*dynprogindex_major),
				pairs,&path,leftpair,rightpair,
				chroffset,chrhigh,queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
				jump_late_p,pairpool,dynprogL,dynprogM,dynprogR,
				last_genomedp5,last_genomedp3,maxpeelback,defect_rate,/*finalp*/true);

      if (filledp == true) {
	/* Discard gap */
	debug(printf("Discarding gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
      } else {
	/* Replace the gap */
	debug(printf("Replacing gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_push_existing(pairs,pairptr);
#endif
      }

    } else if (pair->genomejump > pair->queryjump + minintronlen) {
      /* Previously was 2*MININTRONLEN, and comment said needed space for two introns */
      /* We will make the score matrices nearly square */
      pairptr = path;		/* save */
      path = Pairpool_pop(path,&pair);

      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3 (dir %d): Traversing paired gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d, queryjump %d, genomejump %d\n",
		   cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos,pair->queryjump,pair->genomejump));
      /* fprintf(stderr,"donor prob %f, acceptor prob %f\n",pair->donor_prob,pair->acceptor_prob); */
      pairs = traverse_genome_gap(&filledp,&(*shiftp),&(*dynprogindex_minor),&(*dynprogindex_major),
				  pairs,&path,leftpair,rightpair,chrnum,chroffset,chrhigh,
				  queryseq_ptr,queryuc_ptr,querylength,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  pairpool,dynprogL,dynprogM,dynprogR,last_genomedp5,last_genomedp3,
				  maxpeelback,defect_rate,finalp,simplep);
      /* Previously had forcep == true, because previously thought that adding large gap is not a good solution */

      if (filledp == true) {
	/* Discard the gap */
	debug(printf("Discarding gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
      } else {
	/* Replace the gap */
	debug(printf("Replacing gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_push_existing(pairs,pairptr);
#endif
      }

    } else if (pair->genomejump > pair->queryjump + SINGLESLEN) {
      /* Intron length shorter than MININTRONLEN_FINAL.  Just replace the gap */
      debug(printf("Short intron; not candidate for final calculation.  Replacing gap ");
	    Pair_dump_one(pair,true);
	    printf("\n"));
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else {
      /* Single gap; force fill */
      pairptr = path;		/* save */
      path = Pairpool_pop(path,&pair);

      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3 (dir %d): Traversing single gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d.  queryjump = %d, genomejump = %d\n",
		   cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos,
		   pair->queryjump,pair->genomejump));
      pairs = traverse_single_gap(&filledp,&(*dynprogindex_minor),pairs,&path,leftpair,rightpair,
				  chroffset,chrhigh,
				  queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
				  jump_late_p,pairpool,dynprogM,last_genomedp5,last_genomedp3,
				  maxpeelback,defect_rate,/*forcep*/false,finalp);

      if (filledp == true) {
	/* Discard the gap */
	debug(printf("Discarding gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
      } else {
	/* Replace the gap */
	debug(printf("Replacing gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_push_existing(pairs,pairptr);
#endif
      }
    }
  }

  /* Handle last entry if not a gap */
  if (path != NULL && ((Pair_T) path->first)->gapp == false) {
    pair = (Pair_T) path->first;
    pairs = List_transfer_one(pairs,&path);
  }

  debug(printf("\n** Finishing build_pairs_introns\n"));
  return pairs;
}  


static int
score_alignment (int *nmatches, int *nmismatches, int *nindels,
#ifdef COMPLEX_DIRECTION
		 int *indel_alignment_score,
#endif
		 int *nsemicanonical, int *nnoncanonical, List_T pairs, int cdna_direction) {
  int ncanonical;		/* Do not return this; use score_introns instead */
  int nunknowns, qopens, qindels, topens, tindels;
  double min_splice_prob;

  Pair_fracidentity(&(*nmatches),&nunknowns,&(*nmismatches),&qopens,&qindels,&topens,&tindels,
		    &ncanonical,&(*nsemicanonical),&(*nnoncanonical),&min_splice_prob,
		    pairs,cdna_direction);
  debug11(printf("%d matches, %d nmismatches, %d+%d qgaps, %d+%d tgaps => alignment_score is %d\n",
		 *nmatches,*nmismatches,qopens,qindels,topens,tindels,
		 MATCH*(*nmatches) + MISMATCH*(*nmismatches) + QOPEN*(qopens + qindels) + TOPEN*(topens + tindels)));

  debug(printf("%d matches, %d nmismatches, %d+%d qgaps, %d+%d tgaps => alignment_score is %d\n",
	       *nmatches,*nmismatches,qopens,qindels,topens,tindels,
	       MATCH*(*nmatches) + MISMATCH*(*nmismatches) + QOPEN*(qopens + qindels) + TOPEN*(topens + tindels)));

#ifdef COMPLEX_DIRECTION
  *indel_alignment_score = QOPEN*(qopens + qindels) + TOPEN*(topens + tindels);
#endif

  *nindels = qindels + tindels;
  return MATCH*(*nmatches) + MISMATCH*(*nmismatches) + QOPEN*(qopens + qindels) + TOPEN*(topens + tindels);
}  


static List_T
score_introns (double *max_intron_score, double *avg_donor_score, double *avg_acceptor_score,
	       int *nknown, int *ncanonical, int *nbadintrons, List_T path, int cdna_direction, bool watsonp,
	       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh
#ifdef WASTE	       
	       , Pairpool_T pairpool
#endif
	       ) {
  List_T pairs = NULL, pairptr, p;
  Pair_T pair, leftpair, rightpair;
  Univcoord_T splicesitepos;
  int minintronlen;
  double donor_score, acceptor_score;
  int nintrons = 0;
  int total_matches, total_denominator;
  int max_neighborhood_score, neighborhood_score, neighborhood_length;
#if 0
  char gbuffer1[MAXENT_MAXLENGTH];
#endif

  debug11(printf("\n** Starting score_introns with cdna_direction %d\n",cdna_direction));
  debug11(Pair_dump_list(path,true));

  minintronlen = MININTRONLEN_FINAL;

  *max_intron_score = *avg_donor_score = *avg_acceptor_score = 0.0;
  *nknown = *ncanonical = *nbadintrons = 0;

  total_matches = total_denominator = 0;
  for (p = path; p != NULL; p = p->rest) {
    pair = (Pair_T) p->first;
    if (pair->gapp == true) {
      /* Skip */
    } else {
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	total_matches++;
      }
      total_denominator++;
    }
  }


  while (path != NULL) {
    /* pairptr = path; */
    /* path = Pairpool_pop(path,&pair); */
    pair = (Pair_T) path->first;

    if (pair->gapp == false) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->queryjump > nullgap) {
      debug11(printf("pair->queryjump %d > nullgap %d\n",pair->queryjump,nullgap));
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      debug11(printf("pair->queryjump %d > pair->genomejump %d + EXTRAQUERYGAP %d\n",
		     pair->queryjump,pair->genomejump,EXTRAQUERYGAP));
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else if (pair->genomejump > pair->queryjump + minintronlen) {
      debug11(printf("pair->genomejump %d > pair->queryjump %d + minintronlen %d\n",
		     pair->genomejump,pair->queryjump,minintronlen));
      pairptr = path;	/* save */
      path = Pairpool_pop(path,&pair);

      /* Look at right neighborhood */
      max_neighborhood_score = neighborhood_score = 0;
      neighborhood_length = 0;
      for (p = pairs; p != NULL && neighborhood_length < 25 && ((Pair_T) (p->first))->gapp == false; p = p->rest) {
	rightpair = p->first;
	if (rightpair->comp == MATCH_COMP || rightpair->comp == DYNPROG_MATCH_COMP || rightpair->comp == AMBIGUOUS_COMP) {
	  neighborhood_score += 1;
	} else {
	  neighborhood_score -= 3;
	}
	if (neighborhood_score > max_neighborhood_score) {
	  max_neighborhood_score = neighborhood_score;
	}
	neighborhood_length += 1;
      }

      debug11(printf("right neighborhood: max_neighborhood_score %d, neighborhood_length %d\n",
		     max_neighborhood_score,neighborhood_length));
      if (max_neighborhood_score >= 6 ||
	  (neighborhood_length < 10 && max_neighborhood_score > neighborhood_length - 1)) {
	/* Alignment in right neighborhood okay.  Look at left neighborhood */
	max_neighborhood_score = neighborhood_score = 0;
	neighborhood_length = 0;
	for (p = path; p != NULL && neighborhood_length < 25 && ((Pair_T) (p->first))->gapp == false; p = p->rest) {
	  leftpair = p->first;
	  if (leftpair->comp == MATCH_COMP || leftpair->comp == DYNPROG_MATCH_COMP || leftpair->comp == AMBIGUOUS_COMP) {
	    neighborhood_score += 1;
	  } else {
	    neighborhood_score -= 3;
	  }
	  if (neighborhood_score > max_neighborhood_score) {
	    max_neighborhood_score = neighborhood_score;
	  }
	  neighborhood_length += 1;
	}

	debug11(printf("left neighborhood: max_neighborhood_score %d, neighborhood_length %d\n",
		       max_neighborhood_score,neighborhood_length));
	if (max_neighborhood_score >= 6 ||
	    (neighborhood_length < 10 && max_neighborhood_score > neighborhood_length - 1)) {
	  /* Alignment in left neighborhood okay */
	  leftpair = path->first;
	  rightpair = pairs->first;

	  debug11(printf("pair->comp = %c\n",pair->comp));

	  if (cdna_direction == +1) {
	    if (watsonp) {
	      splicesitepos = leftpair->genomepos + 1;
	      debug11(printf("1. looking up splicesites_iit for donor at #%d:%u..%u\n",chrnum,splicesitepos,splicesitepos+1));
	      if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
									splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/+1)) {
		debug11(printf(" => known\n"));
		pair->knowngapp = true;
		donor_score = 1.0;
	      } else {
		donor_score = Maxent_hr_donor_prob(chroffset + splicesitepos,chroffset);
	      }

	      splicesitepos = rightpair->genomepos;
	      debug11(printf("2. looking up splicesites_iit for acceptor at #%d:%u..%u\n",chrnum,splicesitepos,splicesitepos+1));
	      if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
									splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/+1)) {
		debug11(printf(" => known\n"));
		pair->knowngapp = true;
		acceptor_score = 1.0;
	      } else {
		acceptor_score = Maxent_hr_acceptor_prob(chroffset + splicesitepos,chroffset);
	      }

	    } else {
	      splicesitepos = (chrhigh - chroffset) - leftpair->genomepos;
	      debug11(printf("3. looking up splicesites_iit for donor at #%d:%u..%u\n",chrnum,splicesitepos,splicesitepos+1));
	      if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
									splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/-1)) {
		debug11(printf(" => known\n"));
		pair->knowngapp = true;
		donor_score = 1.0;
	      } else {
		donor_score = Maxent_hr_antidonor_prob(chroffset + splicesitepos,chroffset);
	      }

	      splicesitepos = (chrhigh - chroffset) - rightpair->genomepos + 1;
	      debug11(printf("4. looking up splicesites_iit for acceptor at #%d:%u..%u\n",chrnum,splicesitepos,splicesitepos+1));
	      if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
									splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/-1)) {
		debug11(printf(" => known\n"));
		pair->knowngapp = true;
		acceptor_score = 1.0;
	      } else {
		acceptor_score = Maxent_hr_antiacceptor_prob(chroffset + splicesitepos,chroffset);
	      }
	    }
	    debug11(printf("donor score at %u is %f, watson %d, cdna_direction %d, comp %c\n",
			   leftpair->genomepos,donor_score,watsonp,cdna_direction,pair->comp));
	    debug11(printf("acceptor score at %u is %f, watson %d, cdna_direction %d, comp %c\n",
			   rightpair->genomepos,acceptor_score,watsonp,cdna_direction,pair->comp));
	    nintrons += 1;
	    if (pair->knowngapp == true) {
	      *nknown += 1;
	      *ncanonical += 1;
	    } else if (pair->comp == FWD_CANONICAL_INTRON_COMP) {
	      *ncanonical += 1;
	    } else if (donor_score < 0.9 && acceptor_score < 0.9) {
	      *nbadintrons = 1;
	    }
	    *avg_donor_score += donor_score;
	    *avg_acceptor_score += acceptor_score;
	    if (donor_score + acceptor_score > *max_intron_score) {
	      *max_intron_score = donor_score + acceptor_score;
	    }

	  } else if (cdna_direction == -1) {

#if 0
	    make_complement_buffered(gbuffer1,&(genomicuc_ptr[leftpair->genomepos - ACCEPTOR_MODEL_RIGHT_MARGIN]),
				     ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1);
	    acceptor_score = Maxent_acceptor_prob(gbuffer1);
	    make_complement_buffered(gbuffer1,&(genomicuc_ptr[rightpair->genomepos - DONOR_MODEL_RIGHT_MARGIN - 1]),
				     DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1);
	    donor_score = Maxent_donor_prob(gbuffer1);
#endif

	    if (watsonp) {
	      splicesitepos = leftpair->genomepos + 1;
	      debug11(printf("5. looking up splicesites_iit for acceptor at #%d:%u..%u\n",chrnum,splicesitepos,splicesitepos+1));
	      if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
									splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/-1)) {
		debug11(printf(" => known\n"));
		pair->knowngapp = true;
		acceptor_score = 1.0;
	      } else {
		acceptor_score = Maxent_hr_antiacceptor_prob(chroffset + splicesitepos,chroffset);
	      }


	      splicesitepos = rightpair->genomepos;
	      debug11(printf("6. looking up splicesites_iit for donor at #%d:%u..%u\n",chrnum,splicesitepos,splicesitepos+1));
	      if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
									splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/-1)) {
		debug11(printf(" => known\n"));
		pair->knowngapp = true;
		donor_score = 1.0;
	      } else {
		donor_score = Maxent_hr_antidonor_prob(chroffset + splicesitepos,chroffset);
	      }

	    } else {
	      splicesitepos = (chrhigh - chroffset) - leftpair->genomepos;
	      debug11(printf("7. looking up splicesites_iit for acceptor at #%d:%u..%u\n",chrnum,splicesitepos,splicesitepos+1));
	      if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
									splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/+1)) {
		debug11(printf(" => known\n"));
		pair->knowngapp = true;
		acceptor_score = 1.0;
	      } else {
		acceptor_score = Maxent_hr_acceptor_prob(chroffset + splicesitepos,chroffset);
	      }

	      splicesitepos = (chrhigh - chroffset) - rightpair->genomepos + 1;
	      debug11(printf("8. looking up splicesites_iit for donor at #%d:%u..%u\n",chrnum,splicesitepos,splicesitepos+1));
	      if (splicesites_iit && IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
									splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/+1)) {
		debug11(printf(" => known\n"));
		pair->knowngapp = true;
		donor_score = 1.0;
	      } else {
		donor_score = Maxent_hr_donor_prob(chroffset + splicesitepos,chroffset);
	      }
	    }
	    debug11(printf("donor score at %u is %f, watson %d, cdna_direction %d, comp %c\n",
			   leftpair->genomepos,donor_score,watsonp,cdna_direction,pair->comp));
	    debug11(printf("acceptor score at %u is %f, watson %d, cdna_direction %d, comp %c\n",
			   rightpair->genomepos,acceptor_score,watsonp,cdna_direction,pair->comp));
	    nintrons += 1;
	    if (pair->knowngapp == true) {
	      *nknown += 1;
	      *ncanonical += 1;
	    } else if (pair->comp == REV_CANONICAL_INTRON_COMP) {
	      *ncanonical += 1;
	    } else if (donor_score < 0.9 && acceptor_score < 0.9) {
	      *nbadintrons += 1;
	    }
	    *avg_donor_score += donor_score;
	    *avg_acceptor_score += acceptor_score;
	    if (donor_score + acceptor_score > *max_intron_score) {
	      *max_intron_score = donor_score + acceptor_score;
	    }

	  }
	  debug11(printf("\n"));
	}
      }

#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_push_existing(pairs,pairptr);
#endif

    } else if (pair->genomejump > pair->queryjump + SINGLESLEN) {
      /* Intron length shorter than MININTRONLEN_FINAL.  Just replace the gap */
      debug11(printf("pair->genomejump %d > pair->queryjump %d + SINGLESLEN %d\n",
		     pair->genomejump,pair->queryjump,SINGLESLEN));

#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif

    } else {
      /* Single gap; force fill */
      debug11(printf("pair->queryjump %d, pair->genomejump %d => single gap\n",
		     pair->queryjump,pair->genomejump));
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_transfer_one(pairs,&path);
#endif
    }
  }

  /* Want average scores */
  if (nintrons > 0) {
    *avg_donor_score /= (double) nintrons;
    *avg_acceptor_score /= (double) nintrons;
  }

  debug11(printf("max_intron_score = %f, avg_donor_score = %f, avg_acceptor_score = %f\n",
		 *max_intron_score,*avg_donor_score,*avg_acceptor_score));
  return pairs;
}  


static int
end_compare (List_T x, List_T y, int cdna_direction, bool watsonp,
	     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
	     bool pairsp) {
  List_T pairs1, pairs2, path1, path2;
  double max_intron_score;

  int nknown_1, ncanonical_1, nsemicanonical_1, nnoncanonical_1, nbadintrons_1;
  int nknown_2, ncanonical_2, nsemicanonical_2, nnoncanonical_2, nbadintrons_2;
  double avg_donor_score_1, avg_acceptor_score_1;
  double avg_donor_score_2, avg_acceptor_score_2;
  int alignment_score_1, alignment_score_2;
  int nmatches_1, nmismatches_1, nmatches_2, nmismatches_2, nindels_1, nindels_2;
#ifdef COMPLEX_DIRECTION
  int indel_alignment_score_1, indel_alignment_score_2;
#endif


  if (pairsp == true) {
    pairs1 = x;
    pairs2 = y;

    path1 = List_reverse(pairs1);
    debug11(printf("Calling score_introns for end_compare on path1\n"));
    pairs1 = score_introns(&max_intron_score,&avg_donor_score_1,&avg_acceptor_score_1,&nknown_1,&ncanonical_1,
			   &nbadintrons_1,path1,cdna_direction,watsonp,chrnum,chroffset,chrhigh);
    alignment_score_1 = score_alignment(&nmatches_1,&nmismatches_1,&nindels_1,
#ifdef COMPLEX_DIRECTION
					&indel_alignment_score_1,
#endif
					&nsemicanonical_1,&nnoncanonical_1,pairs1,cdna_direction);

    path2 = List_reverse(pairs2);
    debug11(printf("Calling score_introns for end_compare on path2\n"));
    pairs2 = score_introns(&max_intron_score,&avg_donor_score_2,&avg_acceptor_score_2,&nknown_2,&ncanonical_2,
			   &nbadintrons_2,path2,cdna_direction,watsonp,chrnum,chroffset,chrhigh);
    alignment_score_2 = score_alignment(&nmatches_2,&nmismatches_2,&nindels_2,
#ifdef COMPLEX_DIRECTION
					&indel_alignment_score_2,
#endif
					&nsemicanonical_2,&nnoncanonical_2,pairs2,cdna_direction);

  } else {
    path1 = x;
    path2 = y;

    debug11(printf("Calling score_introns for end_compare on path1\n"));
    pairs1 = score_introns(&max_intron_score,&avg_donor_score_1,&avg_acceptor_score_1,&nknown_1,&ncanonical_1,
			   &nbadintrons_1,path1,cdna_direction,watsonp,chrnum,chroffset,chrhigh);
    alignment_score_1 = score_alignment(&nmatches_1,&nmismatches_1,&nindels_1,
#ifdef COMPLEX_DIRECTION
					&indel_alignment_score_1,
#endif
					&nsemicanonical_1,&nnoncanonical_1,pairs1,cdna_direction);

    path1 = List_reverse(pairs1);
    debug11(printf("Calling score_introns for end_compare on path2\n"));
    pairs2 = score_introns(&max_intron_score,&avg_donor_score_2,&avg_acceptor_score_2,&nknown_2,&ncanonical_2,
			   &nbadintrons_2,path2,cdna_direction,watsonp,chrnum,chroffset,chrhigh);
    alignment_score_2 = score_alignment(&nmatches_2,&nmismatches_2,&nindels_2,
#ifdef COMPLEX_DIRECTION
					&indel_alignment_score_2,
#endif
					&nsemicanonical_2,&nnoncanonical_2,pairs2,cdna_direction);
    path2 = List_reverse(pairs2);
  }

  if (avg_donor_score_1 > 0.9 && avg_acceptor_score_1 > 0.9 &&
      (avg_donor_score_2 < 0.5 || avg_acceptor_score_2 < 0.5)) {
    debug21(printf("intronscores orig %f,%f > intronscores new %f,%f, so original wins\n",
		   avg_donor_score_1,avg_acceptor_score_1,avg_donor_score_2,avg_acceptor_score_2));
    /* intronscores reveal a clear sensedir */
    return -1;

  } else if (avg_donor_score_2 > 0.9 && avg_acceptor_score_2 > 0.9 &&
	     (avg_donor_score_1 < 0.5 || avg_acceptor_score_1 < 0.5)) {
    debug21(printf("intronscores new %f,%f > intronscores orig %f,%f, so new one wins\n",
		   avg_donor_score_2,avg_acceptor_score_2,avg_donor_score_1,avg_acceptor_score_1));
    /* intronscores reveal a clear sensedir */
    return +1;

  } else if (alignment_score_1 > alignment_score_2 + SCORE_SIGDIFF) {
    debug21(printf("alignment_score_1 %d >> alignment_score_2 %d, so original wins\n",
		   alignment_score_1,alignment_score_2));
    return -1;

  } else if (alignment_score_2 > alignment_score_1 + SCORE_SIGDIFF) {
    debug21(printf("alignment_score_2 %d << alignment_score_1 %d, so new one wins\n",
		   alignment_score_2,alignment_score_1));
    return +1;

  } else if (nnoncanonical_1 < nnoncanonical_2) {
    debug21(printf("nnoncanonical_1 %d < nnoncanonical_2 %d, so original wins\n",
		   nnoncanonical_1,nnoncanonical_2));
    return -1;

  } else if (nnoncanonical_2 < nnoncanonical_1) {
    debug21(printf("nnoncanonical_2 %d < nnoncanonical_1 %d, so new one wins\n",
		   nnoncanonical_2,nnoncanonical_1));
    return +1;

  } else if (avg_donor_score_1 + avg_acceptor_score_1 > avg_donor_score_2 + avg_acceptor_score_2 + PROB_SIGDIFF) {
    debug21(printf("intronscores orig %f+%f > intronscores new %f+%f, so original wins\n",
		   avg_donor_score_1,avg_acceptor_score_1,avg_donor_score_2,avg_acceptor_score_2));
    /* intronscores reveal a preferred sensedir */
    return -1;

  } else if (avg_donor_score_2 + avg_acceptor_score_2 > avg_donor_score_1 + avg_acceptor_score_1 + PROB_SIGDIFF) {
    debug21(printf("intronscores new %f+%f > intronscores orig %f+%f, so new one wins\n",
		   avg_donor_score_2,avg_acceptor_score_2,avg_donor_score_1,avg_acceptor_score_1));
    /* intronscores reveal a preferred sensedir */
    return +1;

  } else if (alignment_score_1 > alignment_score_2) {
    debug21(printf("alignment_score_1 %d > alignment_score_2 %d, so original wins\n",
		   alignment_score_1,alignment_score_2));
    return -1;

  } else if (alignment_score_2 > alignment_score_1) {
    debug21(printf("alignment_score_2 %d < alignment_score_1 %d, so new one wins\n",
		   alignment_score_2,alignment_score_1));
    return +1;

  } else {
    debug21(printf("scores all equal\n"));
    return 0;
  }
}


#if 0
static List_T
filter_goodness_hmm (bool *filterp, List_T pairs, double defect_rate) {
  Pair_T pair;
  List_T path, p;
  float prev_vprob_good = 0.0, prev_vprob_bad = 0.0, vprob_good, vprob_bad;
  float good_incr_prob, bad_incr_prob;
  float emission_prob;
  State_T state;

  if (defect_rate == 0.0) {
    defect_rate = 0.001;
  }

  debug5(printf("Beginning filter_goodness_hmm with defect rate %f\n",defect_rate));

  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = (Pair_T) List_head(p);
    debug5(printf("hmm querypos %d (%c %c %c): ",pair->querypos,pair->genome,pair->comp,pair->cdna));

    /* state: GOOD */
    if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
      emission_prob = 1.0 - defect_rate;
    } else {
      emission_prob = defect_rate;
    }

#ifdef COMPUTE_LOG
    good_incr_prob = log(emission_prob) + log(/*transition_prob*/0.99); /* Prob(prev good state -> good state) */
    bad_incr_prob = log(emission_prob) + log(/*transition_prob*/0.10);  /* Prob(prev bad state -> good state) */
#else
    good_incr_prob = fasterlog(emission_prob) + LOG_99; /* Prob(prev good state -> good state) */
    bad_incr_prob = fasterlog(emission_prob) + LOG_10;  /* Prob(prev bad state -> good state) */
#endif

    debug5(printf("state GOOD: %f+%f %f+%f ",prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob));
    if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
      vprob_good = prev_vprob_good + good_incr_prob;
      pair->vstate_good = GOOD;
      debug5(printf(" =>GOOD.  "));
    } else {
      vprob_good = prev_vprob_bad + bad_incr_prob;
      pair->vstate_good = BAD;
      debug5(printf(" =>BAD.   "));
    }

    /* state: BAD */
    if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
      emission_prob = 0.25;
    } else {
      emission_prob = 0.75;
    }

#ifdef COMPUTE_LOG
    good_incr_prob = log(emission_prob) + log(/*transition_prob*/0.01);  /* Prob(prev good state -> bad state) */
    bad_incr_prob = log(emission_prob) + log(/*transition_prob*/0.90);  /* Prob(prev bad state -> bad state) */
#else
    good_incr_prob = fasterlog(emission_prob) + LOG_01; /* Prob(prev good state -> bad state) */
    bad_incr_prob = fasterlog(emission_prob) + LOG_90; /* Prob(prev bad state -> bad state) */
#endif

    debug5(printf("state BAD: %f+%f %f+%f ",prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob));
    if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
      vprob_bad = prev_vprob_good + good_incr_prob;
      pair->vstate_bad = GOOD;
      debug5(printf(" =>GOOD.\n"));
    } else {
      vprob_bad = prev_vprob_bad + bad_incr_prob;
      pair->vstate_bad = BAD;
      debug5(printf(" =>BAD.\n"));
    }

    prev_vprob_good = vprob_good;
    prev_vprob_bad = vprob_bad;
  }

  if (prev_vprob_good > prev_vprob_bad) {
    state = GOOD;
  } else {
    state = BAD;
  }

  path = List_reverse(pairs);
  pairs = (List_T) NULL;

  *filterp = false;
  while (path != NULL) {
    pair = path->first;
    pair->state = state;

#ifdef DEBUG5
    Pair_dump_one(pair,/*zerobasedp*/false);
    printf("\n");
#endif

    if (state == GOOD) {
      pairs = List_transfer_one(pairs,&path);
      state = pair->vstate_good;
    } else {
      *filterp = true;
      path = path->rest;
      state = pair->vstate_bad;
    }
  }
  
  return pairs;
}
#endif


#if 0
static List_T
filter_indels_hmm (bool *filterp, List_T pairs) {
  Pair_T pair;
  List_T path, p;
  float prev_vprob_good = 0.0, prev_vprob_bad = 0.0, vprob_good, vprob_bad;
  float good_incr_prob, bad_incr_prob;
  float emission_prob;
  State_T state;

  debug5(printf("Beginning filter_indels_hmm\n"));

  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = (Pair_T) List_head(p);
    debug5(printf("indels querypos %d (%c %c %c): ",pair->querypos,pair->genome,pair->comp,pair->cdna));

    /* state: GOOD */
    /* These emission probs should add to 1.0 */
    if (pair->comp != INDEL_COMP) {
      emission_prob = 0.9999;	/* Prob(good state -> match/mismatch) */
    } else if (pair->comp != SHORTGAP_COMP) {
      emission_prob = 0.9999;	/* Prob(good state -> match/mismatch) */
    } else {
      emission_prob = 0.0001;	/* Prob(good state -> indel) */
    }

    /* These transition probs should complement those for state BAD */
#ifdef COMPUTE_LOG
    good_incr_prob = log(emission_prob) + log(/*transition_prob*/0.9999);  /* Prob(prev good state -> good state) */
    bad_incr_prob = log(emission_prob) + log(/*transition_prob*/0.25);   /* Prob(prev bad state -> good state) */
#else
    good_incr_prob = fasterlog(emission_prob) + LOG_9999;  /* Prob(prev good state -> good state) */
    bad_incr_prob = fasterlog(emission_prob) + LOG_25;   /* Prob(prev bad state -> good state) */
#endif

    debug5(printf("state GOOD: %f+%f %f+%f ",prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob));
    if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
      vprob_good = prev_vprob_good + good_incr_prob;
      pair->vstate_good = GOOD;
      debug5(printf(" =>GOOD.  "));
    } else {
      vprob_good = prev_vprob_bad + bad_incr_prob;
      pair->vstate_good = BAD;
      debug5(printf(" =>BAD.   "));
    }

    /* state: BAD */
    /* These emission probs should add to 1.0 */
    if (pair->comp != INDEL_COMP && pair->comp != SHORTGAP_COMP) {
      emission_prob = 0.5; 	/* Prob(bad state -> match/mismatch) */
    } else {
      emission_prob = 0.5;	/* Prob(bad state -> indel) */
    }

#ifdef COMPUTE_LOG
    good_incr_prob = log(emission_prob) + log(/*transition_prob*/0.0001);   /* Prob(prev good state -> bad state) */
    bad_incr_prob = log(emission_prob) + log(/*transition_prob*/0.75);  /* Prob(prev bad state -> bad state) */
#else
    good_incr_prob = fasterlog(emission_prob) + LOG_0001;   /* Prob(prev good state -> bad state) */
    bad_incr_prob = fasterlog(emission_prob) + LOG_75;  /* Prob(prev bad state -> bad state) */
#endif

    debug5(printf("state BAD: %f+%f %f+%f ",prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob));
    if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
      vprob_bad = prev_vprob_good + good_incr_prob;
      pair->vstate_bad = GOOD;
      debug5(printf(" =>GOOD.\n"));
    } else {
      vprob_bad = prev_vprob_bad + bad_incr_prob;
      pair->vstate_bad = BAD;
      debug5(printf(" =>BAD.\n"));
    }

    prev_vprob_good = vprob_good;
    prev_vprob_bad = vprob_bad;
  }

  if (prev_vprob_good > prev_vprob_bad) {
    state = GOOD;
  } else {
    state = BAD;
  }

  path = List_reverse(pairs);
  pairs = (List_T) NULL;

  *filterp = false;
  while (path != NULL) {
    pair = path->first;
    pair->state = state;

#ifdef DEBUG5
    Pair_dump_one(pair,/*zerobasedp*/false);
    printf("\n");
#endif

    if (state == GOOD) {
      pairs = List_transfer_one(pairs,&path);
      state = pair->vstate_good;
    } else {
      *filterp = true;
      path = path->rest;
      state = pair->vstate_bad;
    }
  }

  return pairs;
}
#endif



bool
Stage3_short_alignment_p (struct Pair_T *pairarray, int npairs, int querylength) {
  int querystart, queryend;

  if (npairs == 0) {
    return true;
  } else {
    querystart = pairarray[0].querypos;
    queryend = pairarray[npairs-1].querypos;
    if (queryend - querystart + 1 < querylength/3) {
      return true;
    } else {
      return false;
    }
  }
}


#if 0
/* Uses hmm */
/* Modified from Substring_bad_stretch_p */
bool
Stage3_bad_stretch_p (struct Pair_T *pairarray, int npairs, int pos5, int pos3) {
  Pair_T pair;
  int i;
  double vprob_good, vprob_bad, prev_vprob_good, prev_vprob_bad, good_incr_prob, bad_incr_prob;
  bool indelp = false, mismatchp, matchp;

  /* Initialize priors */
#ifdef COMPUTE_LOG
  prev_vprob_good = log(0.99);
  prev_vprob_bad = log(0.01);
#else
  prev_vprob_good = LOG_99;
  prev_vprob_bad = LOG_01;
#endif

  for (i = 0; i < npairs; i++) {
    pair = &(pairarray[i]);
    if (pair->querypos < pos5) {
      /* Skip */
      matchp = mismatchp = false;
    } else if (pair->querypos >= pos3) {
      /* Skip */
      matchp = mismatchp = false;
    } else if (pair->gapp == true) {
      /* Skip */
      matchp = mismatchp = false;
      indelp = false;
    } else if (pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
      if (indelp == true) {
	/* Skip, because we count each indel just once  */
	matchp = mismatchp = false;
      } else {
	/* Count each gap as a mismatch */
	mismatchp = true;
	matchp = false;
	indelp = true;
      }
    } else if (pair->comp == MISMATCH_COMP) {
      mismatchp = true;
      matchp = false;
      indelp = false;
    } else {
      mismatchp = false;
      matchp = true;
      indelp = false;
    }
      
    if (mismatchp == true) {
      debug9(printf("querypos %d (mismatch): ",pair->querypos));

      /* state: GOOD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.001) + log(/*transition_prob*/0.999);
      bad_incr_prob = log(/*emission_prob*/0.001) + log(/*transition_prob*/0.001);
#else
      good_incr_prob = LOG_01_999;
      bad_incr_prob = LOG_01_001;
#endif

      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_good = prev_vprob_good + good_incr_prob;
      } else {
	/* vprob_good = prev_vprob_bad + bad_incr_prob; */
	return true;
      }

      /* state: BAD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.75) + log(/*transition_prob*/0.001);
      bad_incr_prob = log(/*emission_prob*/0.75) + log(/*transition_prob*/0.999);
#else
      good_incr_prob = LOG_75_001;
      bad_incr_prob = LOG_75_999;
#endif
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_bad = prev_vprob_good + good_incr_prob;
      } else {
	vprob_bad = prev_vprob_bad + bad_incr_prob;
      }

      debug9(printf("vprob_good %f, vprob_bad %f",vprob_good,vprob_bad));
      prev_vprob_good = vprob_good;
      prev_vprob_bad = vprob_bad;
      debug9(printf("\n"));

    } else if (matchp == true) {
      debug9(printf("querypos %d (match): ",pair->querypos));

      /* state: GOOD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.999) + log(/*transition_prob*/0.999);
      bad_incr_prob = log(/*emission_prob*/0.999) + log(/*transition_prob*/0.001);
#else
      good_incr_prob = LOG_99_999;
      bad_incr_prob = LOG_99_001;
#endif

      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_good = prev_vprob_good + good_incr_prob;
      } else {
	/* vprob_good = prev_vprob_bad + bad_incr_prob; */
	return true;
      }

      /* state: BAD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.25) + log(/*transition_prob*/0.001);
      bad_incr_prob = log(/*emission_prob*/0.25) + log(/*transition_prob*/0.999);
#else
      good_incr_prob = LOG_25_001;
      bad_incr_prob = LOG_25_999;
#endif
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_bad = prev_vprob_good + good_incr_prob;
      } else {
	vprob_bad = prev_vprob_bad + bad_incr_prob;
      }

      debug9(printf("vprob_good %f, vprob_bad %f",vprob_good,vprob_bad));
      prev_vprob_good = vprob_good;
      prev_vprob_bad = vprob_bad;
      debug9(printf("\n"));
    }
    
  }

  return false;
}
#else

/* Uses a window */
bool
Stage3_bad_stretch_p (struct Pair_T *pairarray, int npairs, int pos5, int pos3) {
  int *nindels;
  struct Pair_T *ptr;
  Pair_T pair;
  int i;

  nindels = (int *) CALLOCA(npairs,sizeof(int));

  i = 0;
  ptr = pairarray;
  while (i < npairs) {
    pair = ptr++;
    i++;

    if (pair->querypos < pos5) {
      /* Skip */
    } else if (pair->querypos >= pos3) {
      /* Skip */
    } else if (pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
      if (pair->genome == ' ') {
	nindels[i] = 1;
	while (i < npairs && pair->genome == ' ') {
	  pair = ptr++;
	  i++;
	}
	ptr--;
	i--;

      } else {
	nindels[i] = 1;
	while (i < npairs && pair->cdna == ' ') {
	  pair = ptr++;
	  i++;
	}
	ptr--;
	i--;

      }
    }
  }

  /* Compute cumulative count of indel openings */
  for (i = 1; i < npairs; i++) {
    nindels[i] += nindels[i-1];
  }
   
  /* Look for more than 3 indel openings in a span of 25 pairs */
  for (i = 0; i < npairs - 25; i++) {
    if (nindels[i+25] - nindels[i] > 3) {
      FREEA(nindels);
      return true;
    }
  }

  FREEA(nindels);

  return false;
}

#endif



int
Stage3_good_part (struct Pair_T *pairarray, int npairs, int pos5, int pos3) {
  Pair_T pair;
  int i;
  double vprob_good, vprob_bad, prev_vprob_good, prev_vprob_bad, good_incr_prob, bad_incr_prob;
  bool indelp = false, mismatchp, matchp, stopp;
  int ngoodleft, ngoodright;

  /* Initialize priors */
#ifdef COMPUTE_LOG
  prev_vprob_good = log(0.99);
  prev_vprob_bad = log(0.01);
#else
  prev_vprob_good = LOG_99;
  prev_vprob_bad = LOG_01;
#endif

  ngoodleft = 0;
  stopp = false;
  for (i = 0; i < npairs && stopp == false; i++) {
    pair = &(pairarray[i]);
    if (pair->querypos < pos5) {
      /* Skip */
      matchp = mismatchp = false;
    } else if (pair->querypos >= pos3) {
      /* Skip */
      matchp = mismatchp = false;
    } else if (pair->gapp == true) {
      /* Skip */
      matchp = mismatchp = false;
      indelp = false;
    } else if (pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
      if (indelp == true) {
	/* Skip, because we count each indel just once  */
	matchp = mismatchp = false;
      } else {
	/* Count each gap as a mismatch */
	mismatchp = true;
	matchp = false;
	indelp = true;
      }
    } else if (pair->comp == MISMATCH_COMP) {
      mismatchp = true;
      matchp = false;
      indelp = false;
    } else {
      mismatchp = false;
      matchp = true;
      indelp = false;
    }
      
    if (mismatchp == true) {
      debug9(printf("querypos %d (mismatch): ",pair->querypos));

      /* state: GOOD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.001) + log(/*transition_prob*/0.999);
      bad_incr_prob = log(/*emission_prob*/0.001) + log(/*transition_prob*/0.001);
#else
      good_incr_prob = LOG_01_999;
      bad_incr_prob = LOG_01_001;
#endif

      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_good = prev_vprob_good + good_incr_prob;
	ngoodleft++;
      } else {
	/* vprob_good = prev_vprob_bad + bad_incr_prob; */
	stopp = true;
      }

      /* state: BAD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.75) + log(/*transition_prob*/0.001);
      bad_incr_prob = log(/*emission_prob*/0.75) + log(/*transition_prob*/0.999);
#else
      good_incr_prob = LOG_75_001;
      bad_incr_prob = LOG_75_999;
#endif
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_bad = prev_vprob_good + good_incr_prob;
      } else {
	vprob_bad = prev_vprob_bad + bad_incr_prob;
      }

      debug9(printf("vprob_good %f, vprob_bad %f",vprob_good,vprob_bad));
      prev_vprob_good = vprob_good;
      prev_vprob_bad = vprob_bad;
      debug9(printf("\n"));

    } else if (matchp == true) {
      debug9(printf("querypos %d (match): ",pair->querypos));

      /* state: GOOD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.999) + log(/*transition_prob*/0.999);
      bad_incr_prob = log(/*emission_prob*/0.999) + log(/*transition_prob*/0.001);
#else
      good_incr_prob = LOG_99_999;
      bad_incr_prob = LOG_99_001;
#endif

      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_good = prev_vprob_good + good_incr_prob;
	ngoodleft++;
      } else {
	/* vprob_good = prev_vprob_bad + bad_incr_prob; */
	stopp = true;
      }

      /* state: BAD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.25) + log(/*transition_prob*/0.001);
      bad_incr_prob = log(/*emission_prob*/0.25) + log(/*transition_prob*/0.999);
#else
      good_incr_prob = LOG_25_001;
      bad_incr_prob = LOG_25_999;
#endif
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_bad = prev_vprob_good + good_incr_prob;
      } else {
	vprob_bad = prev_vprob_bad + bad_incr_prob;
      }

      debug9(printf("vprob_good %f, vprob_bad %f",vprob_good,vprob_bad));
      prev_vprob_good = vprob_good;
      prev_vprob_bad = vprob_bad;
      debug9(printf("\n"));
    }
    
  }


  /* Initialize priors */
#ifdef COMPUTE_LOG
  prev_vprob_good = log(0.99);
  prev_vprob_bad = log(0.01);
#else
  prev_vprob_good = LOG_99;
  prev_vprob_bad = LOG_01;
#endif

  ngoodright = 0;
  stopp = false;
  for (i = npairs - 1; i >= 0 && stopp == false; i--) {
    pair = &(pairarray[i]);
    if (pair->querypos < pos5) {
      /* Skip */
      matchp = mismatchp = false;
    } else if (pair->querypos >= pos3) {
      /* Skip */
      matchp = mismatchp = false;
    } else if (pair->gapp == true) {
      /* Skip */
      matchp = mismatchp = false;
      indelp = false;
    } else if (pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
      if (indelp == true) {
	/* Skip, because we count each indel just once  */
	matchp = mismatchp = false;
      } else {
	/* Count each gap as a mismatch */
	mismatchp = true;
	matchp = false;
	indelp = true;
      }
    } else if (pair->comp == MISMATCH_COMP) {
      mismatchp = true;
      matchp = false;
      indelp = false;
    } else {
      mismatchp = false;
      matchp = true;
      indelp = false;
    }
      
    if (mismatchp == true) {
      debug9(printf("querypos %d (mismatch): ",pair->querypos));

      /* state: GOOD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.001) + log(/*transition_prob*/0.999);
      bad_incr_prob = log(/*emission_prob*/0.001) + log(/*transition_prob*/0.001);
#else
      good_incr_prob = LOG_01_999;
      bad_incr_prob = LOG_01_001;
#endif

      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_good = prev_vprob_good + good_incr_prob;
	ngoodright++;
      } else {
	/* vprob_good = prev_vprob_bad + bad_incr_prob; */
	stopp = true;
      }

      /* state: BAD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.75) + log(/*transition_prob*/0.001);
      bad_incr_prob = log(/*emission_prob*/0.75) + log(/*transition_prob*/0.999);
#else
      good_incr_prob = LOG_75_001;
      bad_incr_prob = LOG_75_999;
#endif
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_bad = prev_vprob_good + good_incr_prob;
      } else {
	vprob_bad = prev_vprob_bad + bad_incr_prob;
      }

      debug9(printf("vprob_good %f, vprob_bad %f",vprob_good,vprob_bad));
      prev_vprob_good = vprob_good;
      prev_vprob_bad = vprob_bad;
      debug9(printf("\n"));

    } else if (matchp == true) {
      debug9(printf("querypos %d (match): ",pair->querypos));

      /* state: GOOD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.999) + log(/*transition_prob*/0.999);
      bad_incr_prob = log(/*emission_prob*/0.999) + log(/*transition_prob*/0.001);
#else
      good_incr_prob = LOG_99_999;
      bad_incr_prob = LOG_99_001;
#endif

      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_good = prev_vprob_good + good_incr_prob;
	ngoodright++;
      } else {
	/* vprob_good = prev_vprob_bad + bad_incr_prob; */
	stopp = true;
      }

      /* state: BAD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.25) + log(/*transition_prob*/0.001);
      bad_incr_prob = log(/*emission_prob*/0.25) + log(/*transition_prob*/0.999);
#else
      good_incr_prob = LOG_25_001;
      bad_incr_prob = LOG_25_999;
#endif
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_bad = prev_vprob_good + good_incr_prob;
      } else {
	vprob_bad = prev_vprob_bad + bad_incr_prob;
      }

      debug9(printf("vprob_good %f, vprob_bad %f",vprob_good,vprob_bad));
      prev_vprob_good = vprob_good;
      prev_vprob_bad = vprob_bad;
      debug9(printf("\n"));
    }
    
  }

  debug9(printf("ngoodleft = %d, ngoodright = %d\n",ngoodleft,ngoodright));

  if (ngoodleft > ngoodright) {
    return ngoodleft;
  } else {
    return ngoodright;
  }
}




#if 0
static int
score_nconsecutive (List_T pairs) {
  int score = 0, nconsecutive;
  Pair_T pair;
  List_T path, p;
  float prev_vprob_good = 0.0, prev_vprob_bad = 0.0, vprob_good, vprob_bad;
  float defect_rate = 0.001, good_incr_prob, bad_incr_prob;
  float emission_prob;
  State_T state;


  debug(printf("Beginning score_nconsecutive\n"));

  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = (Pair_T) List_head(p);
    debug5(printf("hmm querypos %d (%c %c %c): ",pair->querypos,pair->genome,pair->comp,pair->cdna));

    /* state: GOOD */
    if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
      emission_prob = 1.0 - defect_rate;
    } else {
      emission_prob = defect_rate;
    }

#ifdef COMPUTE_LOG
    good_incr_prob = log(emission_prob) + log(/*transition_prob*/0.99); /* Prob(prev good state -> good state) */
    bad_incr_prob = log(emission_prob) + log(/*transition_prob*/0.10);  /* Prob(prev bad state -> good state) */
#else
    good_incr_prob = fasterlog(emission_prob) + LOG_99; /* Prob(prev good state -> good state) */
    bad_incr_prob = fasterlog(emission_prob) + LOG_10;  /* Prob(prev bad state -> good state) */
#endif

    debug5(printf("state GOOD: %f+%f %f+%f ",prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob));
    if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
      vprob_good = prev_vprob_good + good_incr_prob;
      pair->vstate_good = GOOD;
      debug5(printf(" =>GOOD.  "));
    } else {
      vprob_good = prev_vprob_bad + bad_incr_prob;
      pair->vstate_good = BAD;
      debug5(printf(" =>BAD.   "));
    }

    /* state: BAD */
    if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
      emission_prob = 0.25;
    } else {
      emission_prob = 0.75;
    }

#ifdef COMPUTE_LOG
    good_incr_prob = log(emission_prob) + log(/*transition_prob*/0.01);  /* Prob(prev good state -> bad state) */
    bad_incr_prob = log(emission_prob) + log(/*transition_prob*/0.90);  /* Prob(prev bad state -> bad state) */
#else
    good_incr_prob = fasterlog(emission_prob) + LOG_01;  /* Prob(prev good state -> bad state) */
    bad_incr_prob = fasterlog(emission_prob) + LOG_90;  /* Prob(prev bad state -> bad state) */
#endif

    debug5(printf("state BAD: %f+%f %f+%f ",prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob));
    if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
      vprob_bad = prev_vprob_good + good_incr_prob;
      pair->vstate_bad = GOOD;
      debug5(printf(" =>GOOD.\n"));
    } else {
      vprob_bad = prev_vprob_bad + bad_incr_prob;
      pair->vstate_bad = BAD;
      debug5(printf(" =>BAD.\n"));
    }

    prev_vprob_good = vprob_good;
    prev_vprob_bad = vprob_bad;
  }

  if (prev_vprob_good > prev_vprob_bad) {
    state = GOOD;
  } else {
    state = BAD;
  }

  nconsecutive = 0;
  path = List_reverse(pairs);

  for (p = path; p != NULL; p = p->rest) {
    pair = p->first;
    debug5(printf("hmm querypos %d (%c %c %c): ",pair->querypos,pair->genome,pair->comp,pair->cdna));
    if (state == GOOD) {
      nconsecutive++;
      debug5(printf("state is GOOD, nconsecutive %d\n",nconsecutive));
      state = pair->vstate_good;
    } else {
      if (nconsecutive > score) {
	score = nconsecutive;
      }
      nconsecutive = 0;
      debug5(printf("state is BAD, score %d\n",score));
      state = pair->vstate_bad;
    }
  }
  
  if (nconsecutive > score) {
    score = nconsecutive;
  }
  debug5(printf("At end, score %d\n",score));

  pairs = List_reverse(path);

  return score;
}
#endif



static List_T
path_compute_dir (double *defect_rate, List_T pairs,
		  int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
#ifdef PMAP
		  char *queryaaseq_ptr,
#endif
		  char *queryseq_ptr, char *queryuc_ptr, int querylength,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		  int maxpeelback,
#ifndef GSNAP
		  Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool, 
#endif
		  Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		  Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3) {
  List_T path = NULL;
  int dynprogindex_minor = DYNPROGINDEX_MINOR, dynprogindex_major = DYNPROGINDEX_MAJOR;
  int iter0, iter1, iter2;
  bool shiftp, incompletep;
  bool shortp, badp, deletep, dual_break_p;
  int matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double min_splice_prob;

  bool indelp, trim5p, trim3p;

  
  iter0 = 0;
  dual_break_p = true;

  while (dual_break_p == true && iter0 < MAXITER_CYCLES) {
    /* path = List_reverse(pairs); */
    /* Need to insert gapholders after Pairpool_join_end5 and Pairpool_join_end3 */
    path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			     pairpool,/*finalp*/false);

#ifdef PMAP
#if 0
    /* Pass 1b: undefine nucleotides around gaps.  path --> path */
    pairs = undefine_nucleotides(queryseq_ptr,querylength,path,pairpool,/*width*/6);
    path = List_reverse(pairs);
#endif
#endif

    /* Pass 2A: solve straight gaps (small).  path --> pairs (for defect rate) */
    debug(printf("\n*** Pass 2A (dir %d): Solve straight gaps (small).  Iteration0 %d\n",
		 cdna_direction,iter0));
    debug(Pair_dump_list(path,true));
    pairs = build_pairs_singles(&dynprogindex_minor,path,/*maxsize*/3,
				chroffset,chrhigh,queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
				jump_late_p,maxpeelback,/*defect_rate*/0.0,pairpool,dynprogM,
				last_genomedp5,last_genomedp3,/*forcep*/false,/*finalp*/false);
#ifdef DEBUG8
    if (stage3debug == POST_SINGLES) {
      path = List_reverse(pairs);
      return path;
    }
#endif

    if (homopolymerp == true) {
      /* Pass 2B: fix adjacent indels */
      /* >>pairs */
#if 0
      /* gapholders shouldn't be necessary before fix_adjacent_indels,
	 but is necessary afterward for build_pairs_singles */
      path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			       pairpool,/*finalp*/false);
      pairs = List_reverse(path);
#endif

      debug(printf("\n*** Pass 2B (dir %d): Fix adjacent indels.  Iteration0 %d\n",
		   cdna_direction,iter0));
      path = fix_adjacent_indels(pairs);
      pairs = List_reverse(path);
      path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			       pairpool,/*finalp*/false);


      /* Pass 2C: solve straight gaps again.  path --> pairs (for defect rate) */
      debug(printf("\n*** Pass 2C (dir %d): Solve straight gaps again (large).  Iteration0 %d\n",
		   cdna_direction,iter0));
      pairs = build_pairs_singles(&dynprogindex_minor,path,/*maxsize*/nullgap,
				  chroffset,chrhigh,queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
				  jump_late_p,maxpeelback,/*defect_rate*/0.0,pairpool,dynprogM,
				  last_genomedp5,last_genomedp3,/*forcep*/false,/*finalp*/false);
      /* <<pairs */
    }

    /* Compute defect rate here */
    Pair_fracidentity(&matches,&unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
		      &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
		      pairs,/*cdna_direction*/0);
    if (matches + mismatches == 0) {
      *defect_rate = 0.0;
    } else {
      *defect_rate = (double) mismatches/(double) (matches + mismatches);
    }
    debug(printf("defect_rate = %f (%d matches, %d mismatches)\n",*defect_rate,matches,mismatches));
    

    /* Pass 3: Smoothing */
    debug(printf("*** Pass 3 (dir %d): Smooth\n",cdna_direction));

    /* Smoothing by probability */
    path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			     pairpool,/*finalp*/false);
    pairs = assign_intron_probs(path,cdna_direction,watsonp,chrnum,chroffset,chrhigh,pairpool);
    Smooth_reset(pairs);
    pairs = Smooth_pairs_by_intronprobs(&badp,pairs,pairpool);

#if 0
    /* Smoothing by netgap.  Can crash or stall, and generally doesn't do anything except for very low-identity alignments. */
    debug(printf("\n*** Pass 1 (dir %d): Initial smoothing by net gap.  Iteration1 %d\n",
		 cdna_direction,iter1));
    pairs = Smooth_pairs_by_netgap(&smoothp,pairs,pairpool);
#endif
      
    /* Smoothing by size: This can undo the short exons found by traverse_dual_genome, so we use protectedp in traverse_dual_genome  */
    debug(printf("*** Pass 3a (dir %d): Smoothing by size.  Iteration0 %d\n",
		 cdna_direction,iter0));
    path = List_reverse(pairs);
    pairs = remove_indel_gaps(path);
    pairs = Smooth_pairs_by_size(&shortp,&deletep,pairs,pairpool,/*stage2_indexsize*/6);
    debug(printf("  => Result of Pass 3a (smoothing): shortp is %d, deletep is %d\n",shortp,deletep));
    debug(Pair_dump_list(pairs,/*zerobasedp*/true));
      
#ifdef DEBUG8
    if (stage3debug == POST_SMOOTHING) {
      path = List_reverse(pairs);
      return path;
    }
#endif


    /* Pass 5: introns */
    /* >>pairs */
    debug(printf("\n*** Pass 5 (dir %d): Smooth and solve dual introns iteratively.  Iteration0 %d\n",
		 cdna_direction,iter0));
    iter1 = 0;
    badp = true;
    while (badp == true && iter1 < MAXITER_SMOOTH_BY_SIZE) {
      /* Pass 5a: single introns */
      debug(printf("*** Pass 5a (dir %d): Solve introns.  Iteration0 %d, iteration1 %d\n",
		   cdna_direction,iter0,iter1));

      iter2 = 0;
      shiftp = true;
      while ((shiftp == true || incompletep == true) && iter2++ < MAXITER_INTRONS) {
	path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
				 pairpool,/*finalp*/false);
	pairs = build_pairs_introns(&shiftp,&incompletep,
				    &dynprogindex_minor,&dynprogindex_major,path,
				    chrnum,chroffset,chrhigh,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,
				    cdna_direction,watsonp,genestrand,jump_late_p,
				    maxpeelback,*defect_rate,pairpool,dynprogL,dynprogM,dynprogR,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool,
#endif
				    last_genomedp5,last_genomedp3,/*finalp*/false,/*simplep*/true);
	debug(printf("  => Result of Pass 5 (introns):\n"));
	debug(Pair_dump_list(pairs,/*zerobasedp*/true));
      }

#ifdef DEBUG8
      if (stage3debug == POST_INTRONS) {
	path = List_reverse(pairs);
	return path;
      }
#endif

      /* Re-evaluate any small exons inserted by build_dual_breaks */
      path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			       pairpool,/*finalp*/false);
      pairs = assign_intron_probs(path,cdna_direction,watsonp,chrnum,chroffset,chrhigh,pairpool);
      Smooth_reset(pairs);
      pairs = Smooth_pairs_by_intronprobs(&badp,pairs,pairpool);

      debug(printf("*** Pass 5b (dir %d): Solve dual introns.  Iteration0 %d, Iteration1 %d\n",
		   cdna_direction,iter0,iter1));
      if (badp == false) {
	debug(printf("  no badp, so do nothing\n"));
      } else {
	debug(printf("  badp is true, so running build_pairs_dualintrons\n"));
	path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
				 pairpool,/*finalp*/false);
	/* XX */
	/* pairs = assign_gap_types(path,cdna_direction,watsonp,queryseq_ptr,
	   chrnum,chroffset,chrhigh,pairpool); */
	
	/* Pass 3b: dual introns.  pairs --> pairs */
	pairs = build_pairs_dualintrons(&dynprogindex_major,path,chrnum,chroffset,chrhigh,
					queryseq_ptr,queryuc_ptr,querylength,
					cdna_direction,watsonp,genestrand,jump_late_p,
					maxpeelback,*defect_rate,pairpool,dynprogL,dynprogR,
					last_genomedp5,last_genomedp3);
	debug(printf("  => Result of Pass 5b (dual introns):\n"));
	debug(Pair_dump_list(pairs,/*zerobasedp*/true));

	/* Trimming needed for path_compute_end5 and path_compute_end3 to add new start and end exons */
	pairs = trim_end5_exons(&indelp,&trim5p,/*ambig_end_length*/0,pairs,dynprogR,chroffset,chrhigh,
				queryseq_ptr,queryuc_ptr,querylength,
				cdna_direction,watsonp,genestrand,jump_late_p,pairpool,*defect_rate);
	path = List_reverse(pairs);
	path = trim_end3_exons(&indelp,&trim3p,/*ambig_end_length_3*/0,path,dynprogL,chroffset,chrhigh,
			       queryseq_ptr,queryuc_ptr,querylength,
			       cdna_direction,watsonp,genestrand,jump_late_p,pairpool,*defect_rate);
	pairs = List_reverse(path);
      }

#ifdef DEBUG8
      if (stage3debug == POST_DUAL_INTRONS) {
	path = List_reverse(pairs);
	return path;
      }
#endif

      iter1++;
      debug(printf("At end of inner loop: iter1 %d, badp %d\n",iter1,badp));
    }

#ifdef DEBUG8
    if (stage3debug == POST_CYCLES) {
      path = List_reverse(pairs);
      return path;
    }
#endif


    /* Pass 99: Fix dual breaks */
    /* >>pairs */
    debug(printf("\n*** Pass 99 (dir %d): Fix dual breaks.  Iteration0 %d\n",cdna_direction,iter0));
    /* pairs = remove_indel_gaps(path); */
    path = List_reverse(pairs);

    pairs = build_dual_breaks(&dual_break_p,&dynprogindex_minor,&dynprogindex_major,path,
			      chrnum,chroffset,chrhigh,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      queryseq_ptr,queryuc_ptr,querylength,
			      cdna_direction,watsonp,genestrand,jump_late_p,pairpool,
			      dynprogL,dynprogM,dynprogR,last_genomedp5,last_genomedp3,maxpeelback,
#ifndef GSNAP
			      oligoindices_minor,diagpool,cellpool,
#endif
			      *defect_rate,/*finalp*/false,/*simplep*/true);


#if 0
    path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			     pairpool,/*finalp*/false);
    debug(Pair_dump_list(path,/*zerobasedp*/true));
    pairs = List_reverse(path);
    debug14(printf("Result of build_dual_breaks\n"));
    debug14(Pair_dump_list(pairs,true));
    debug(printf("Result of build_dual_breaks\n"));
    debug(Pair_dump_list(pairs,true));
#endif
      
#ifdef GSNAP
    /* Too expensive to loop */
    dual_break_p = false;
    /* filterp = false; */
#endif
    iter0++;
    debug(printf("At end of outer loop: dual_break_p %d\n",dual_break_p));
  }

  path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			   pairpool,/*finalp*/false);

  return path;
}


#ifdef STRICT
/* Written for pairs, but works for path also */
static List_T
remove_dual_breaks (List_T pairs) {
  int npairs5, npairs3;
  List_T path;
  bool donep;


  donep = false;
  while (donep == false) {
    if (pairs == NULL) {
      debug(printf("pairs is NULL\n"));
      donep = true;
    } else if (dual_break_p(pairs) == false) {
      debug(printf("dual break not found\n"));
      donep = true;
    } else {
      dual_break_count_npairs(&npairs5,&npairs3,pairs);
      debug(printf("npairs5 %d, npairs3 %d, total npairs %d\n",
		   npairs5,npairs3,List_length(pairs)));
      if (npairs5 < npairs3) {
	/* Remove from 5' end */
	debug(printf("trimming dual break from 5' end, %d pairs\n",npairs5));
	pairs = trim_npairs(pairs,npairs5 + 1 /*for the gap*/);
      } else {
	/* Remove from 3' end */
	path = List_reverse(pairs);
	path = trim_npairs(path,npairs3 + 1 /*for the gap*/);
	pairs = List_reverse(path);
      }
      debug(printf("Now alignment length is %d\n",List_length(pairs)));
    }
  }

  return pairs;
}
#endif


/* Written for pairs, but works for path also */
static List_T
remove_negative_breaks (List_T pairs) {
  int npairs5, npairs3;
  List_T path;
  bool donep;


  donep = false;
  while (donep == false) {
    if (pairs == NULL) {
      debug(printf("pairs is NULL\n"));
      donep = true;
    } else if (negative_break_p(pairs) == false) {
      debug(printf("negative break not found\n"));
      donep = true;
    } else {
      negative_break_count_npairs(&npairs5,&npairs3,pairs);
      debug(printf("npairs5 %d, npairs3 %d, total npairs %d\n",
		   npairs5,npairs3,List_length(pairs)));
      if (npairs5 < npairs3) {
	/* Remove from 5' end */
	debug(printf("trimming negative break from 5' end, %d pairs\n",npairs5));
	pairs = trim_npairs(pairs,npairs5 + 1 /*for the gap*/);
      } else {
	/* Remove from 3' end */
	path = List_reverse(pairs);
	path = trim_npairs(path,npairs3 + 1 /*for the gap*/);
	pairs = List_reverse(path);
      }
      debug(printf("Now alignment length is %d\n",List_length(pairs)));
    }
  }

  return pairs;
}


static List_T
path_compute_end5 (int *ambig_end_length_5, Splicetype_T *ambig_splicetype_5, double *ambig_prob_5,
		   double defect_rate, List_T pairs, int cdna_direction,
		   bool watsonp, int genestrand, bool jump_late_p,
		   char *queryseq_ptr, char *queryuc_ptr, int querylength,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		   Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		   int maxpeelback, Pairpool_T pairpool, Dynprog_T dynprogR) {
  List_T path = NULL;
  int iter1;
  int dynprogindex_minor = DYNPROGINDEX_MINOR;
  int nmatches, nunknowns, nmismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double min_splice_prob;
  bool knownsplice5p, chop_exon_p;
  bool trim5p, indelp;

#if 0
  int distance5, totaljump5;
  int npairs5;

  Pair_T start;
  int querypos;
  Chrpos_T chrstart, chrend;
  List_T all_stage2_starts;
#endif


  *ambig_end_length_5 = 0;
  *ambig_prob_5 = 0.0;

  /* Pass 7: Remove dual or negative breaks */
  /* >>pairs */
  debug(printf("\n*** Pass 7a (dir %d): Remove dual breaks at ends.  Initially alignment length is %d\n",
	       cdna_direction,List_length(path)));
#ifdef STRICT
  pairs = remove_dual_breaks(pairs);
#else
  pairs = remove_negative_breaks(pairs);
#endif
  debug(printf("Final alignment length is %d\n",List_length(pairs)));
  /* <<pairs */


  /* Extend to query end, so we get an accurate count of matches and mismatches */
  /* This is the first extension */
  /* >>pairs */
  debug(printf("\n*** Pass 8 (dir %d): Extend to 5' end and determine distalmedial\n",cdna_direction));
  pairs = build_pairs_end5(&knownsplice5p,&(*ambig_end_length_5),&(*ambig_splicetype_5),&(*ambig_prob_5),
			   &chop_exon_p,&dynprogindex_minor,pairs,
			   chroffset,chrhigh,
			   knownsplice_limit_low,knownsplice_limit_high,
			   queryseq_ptr,queryuc_ptr,
			   cdna_direction,watsonp,genestrand,jump_late_p,
			   maxpeelback,defect_rate,pairpool,dynprogR,
			   /*extendp*/true,/*endalign*/QUERYEND_GAP,/*forcep*/false);

  /* Necessary to insert gaps and assign gap types (fills in cDNA
     insertions, so they don't get trimmed), in case an insertion was
     introduced at ends */
  path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			   pairpool,/*finalp*/false);
  pairs = assign_gap_types(path,cdna_direction,watsonp,queryseq_ptr,
			   chrnum,chroffset,chrhigh,pairpool);

  Pair_fracidentity(&nmatches,&nunknowns,&nmismatches,&qopens,&qindels,&topens,&tindels,
		    &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
		    pairs,cdna_direction);
  if (ncanonical > 0 && nnoncanonical > 0) {
    /* Pass 10: Remove noncanonical end exons: pairs -> pairs */
    debug(printf("\n*** Pass 10 (dir %d): Remove noncanonical end exons\n",cdna_direction));

    if (maximize_coverage_p == true) {
      trim5p = false;
    } else if (knownsplice5p == true) {
      /* Don't trim at known splice sites */
      trim5p = false;
    } else {
      trim5p = true;
    }
    debug(printf("trim5p = %d\n",trim5p));

    /* Using iter1 to avoid the possibility of an infinite loop */
    iter1 = 0;
    while (iter1 < 5 && trim5p == true) {
      pairs = trim_end5_exons(&indelp,&trim5p,*ambig_end_length_5,pairs,dynprogR,chroffset,chrhigh,
			      queryseq_ptr,queryuc_ptr,querylength,
			      cdna_direction,watsonp,genestrand,jump_late_p,pairpool,defect_rate);
      if (indelp == true) {
	pairs = trim_end5_indels(pairs,*ambig_end_length_5,dynprogR,chroffset,chrhigh,
				 queryseq_ptr,queryuc_ptr,
				 watsonp,genestrand,jump_late_p,pairpool,defect_rate);
      }
      if (trim5p == true) {
	pairs = build_pairs_end5(&knownsplice5p,&(*ambig_end_length_5),&(*ambig_splicetype_5),&(*ambig_prob_5),
				 &chop_exon_p,&dynprogindex_minor,pairs,
				 chroffset,chrhigh,
				 knownsplice_limit_low,knownsplice_limit_high,
				 queryseq_ptr,queryuc_ptr,
				 cdna_direction,watsonp,genestrand,jump_late_p,
				 maxpeelback,defect_rate,pairpool,dynprogR,/*extendp*/true,
				 /*endalign*/BEST_LOCAL,/*forcep*/false);
	debug3(printf("AFTER 5' REBUILD\n"));
	debug3(Pair_dump_list(pairs,true));
      }

      /* Stop trimming at known splice sites */
      if (knownsplice5p == true) {
	trim5p = false;
      }
      iter1++;
    }

#if 0
    pairs = build_pairs_end5(&knownsplice5p,&(*ambig_end_length_5),&(*ambig_splicetype_5),&(*ambig_prob_5),
			     &chop_exon_p,&dynprogindex_minor,pairs,
			     chroffset,chrhigh,
			     knownsplice_limit_low,knownsplice_limit_high,
			     queryseq_ptr,queryuc_ptr,
			     cdna_direction,watsonp,genestrand,jump_late_p,
			     maxpeelback,defect_rate,pairpool,dynprogR,/*extendp*/true,
			     /*endalign*/QUERYEND_INDELS,/*forcep*/false);
#endif
  }


#if 0
  debug(printf("Pass 11 (dir %d): Final extension, end5\n",cdna_direction));
  /* This is the second extension */
  /* Perform final extension without gaps so we can compare fwd and rev properly */
  pairs = build_pairs_end5(&knownsplice5p,&(*ambig_end_length_5),&(*ambig_splicetype_5),&(*ambig_prob_5),
			   &chop_exon_p,&dynprogindex_minor,pairs,
			   chroffset,chrhigh,
			   knownsplice_limit_low,knownsplice_limit_high,
			   queryseq_ptr,queryuc_ptr,
			   cdna_direction,watsonp,genestrand,jump_late_p,
			   maxpeelback,defect_rate,pairpool,dynprogR,
			   /*extendp*/true,/*endalign*/QUERYEND_NOGAPS,/*forcep*/false);
#endif

  debug(Pair_dump_list(pairs,true));
  debug(printf("End of path_compute_end5\n"));

  return pairs;
}


static List_T
path_compute_end3 (int *ambig_end_length_3, Splicetype_T *ambig_splicetype_3, double *ambig_prob_3,
		   double defect_rate, List_T path, int cdna_direction,
		   bool watsonp, int genestrand, bool jump_late_p, int querylength,
		   char *queryseq_ptr, char *queryuc_ptr,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		   Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		   int maxpeelback, Pairpool_T pairpool, Dynprog_T dynprogL) {
  List_T pairs = NULL;
  int iter1;
  int dynprogindex_minor = DYNPROGINDEX_MINOR;
  int nmatches, nunknowns, nmismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double min_splice_prob;
  bool knownsplice3p, chop_exon_p;
  bool trim3p, indelp;

#if 0
  int distance3, totaljump3;
  int npairs3;

  Pair_T end;
  int querypos;
  Chrpos_T chrstart, chrend;
  List_T all_stage2_ends;
#endif


  *ambig_end_length_3 = 0;
  *ambig_prob_3 = 0.0;

  /* Pass 7: Remove dual or negative breaks */
  /* >>path */
  debug(printf("\n*** Pass 7b (dir %d): Remove dual breaks at ends.  Initially alignment length is %d\n",
	       cdna_direction,List_length(path)));
#ifdef STRICT
  path = remove_dual_breaks(path);
#else
  path = remove_negative_breaks(path);
#endif
  debug(printf("Final alignment length is %d\n",List_length(path)));
  /* <<path */


  /* Extend to ends */
  /* This is the first extension */
  /* >>path */
  debug(printf("\n*** Pass 8 (dir %d): Extend to 3' end and determine distalmedial\n",cdna_direction));
  path = build_path_end3(&knownsplice3p,&(*ambig_end_length_3),&(*ambig_splicetype_3),&(*ambig_prob_3),
			 &chop_exon_p,&dynprogindex_minor,path,
			 chroffset,chrhigh,querylength,
			 knownsplice_limit_low,knownsplice_limit_high,
			 queryseq_ptr,queryuc_ptr,
			 cdna_direction,watsonp,genestrand,jump_late_p,
			 maxpeelback,defect_rate,pairpool,dynprogL,
			 /*extendp*/true,/*endalign*/QUERYEND_GAP,/*forcep*/false);

  /* Necessary to insert gaps and assign gap types (fills in cDNA
     insertions, so they don't get trimmed), in case an insertion was
     introduced at ends */
  pairs = List_reverse(path);
  path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			   pairpool,/*finalp*/false);
  pairs = assign_gap_types(path,cdna_direction,watsonp,queryseq_ptr,
			   chrnum,chroffset,chrhigh,pairpool);
  path = List_reverse(pairs);

  Pair_fracidentity(&nmatches,&nunknowns,&nmismatches,&qopens,&qindels,&topens,&tindels,
		    &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
		    pairs,cdna_direction);
  if (ncanonical > 0 && nnoncanonical > 0) {
    /* Pass 10: Remove noncanonical end exons: pairs -> pairs */
    debug(printf("\n*** Pass 10 (dir %d): Remove noncanonical end exons\n",cdna_direction));

    if (maximize_coverage_p == true) {
      trim3p = false;
    } else if (knownsplice3p == true) {
      trim3p = false;
    } else {
      trim3p = true;
    }
    debug(printf("trim3p = %d\n",trim3p));

    /* Using iter1 to avoid the possibility of an infinite loop */
    iter1 = 0;
    while (iter1 < 5 && trim3p == true) {
      path = trim_end3_exons(&indelp,&trim3p,*ambig_end_length_3,path,dynprogL,chroffset,chrhigh,
			     queryseq_ptr,queryuc_ptr,querylength,
			     cdna_direction,watsonp,genestrand,jump_late_p,pairpool,defect_rate);
      if (indelp == true) {
	path = trim_end3_indels(path,*ambig_end_length_3,dynprogL,chroffset,chrhigh,
				queryseq_ptr,queryuc_ptr,querylength,
				watsonp,genestrand,jump_late_p,pairpool,defect_rate);
      }
      if (trim3p == true) {
	path = build_path_end3(&knownsplice3p,&(*ambig_end_length_3),&(*ambig_splicetype_3),&(*ambig_prob_3),
			       &chop_exon_p,&dynprogindex_minor,path,
			       chroffset,chrhigh,querylength,
			       knownsplice_limit_low,knownsplice_limit_high,
			       queryseq_ptr,queryuc_ptr,
			       cdna_direction,watsonp,genestrand,jump_late_p,
			       maxpeelback,defect_rate,pairpool,dynprogL,/*extendp*/true,
			       /*endalign*/BEST_LOCAL,/*forcep*/false);
	debug3(printf("AFTER 3' REBUILD\n"));
	debug3(Pair_dump_list(path,true));
      }

      if (knownsplice3p == true) {
	trim3p = false;
      }
      iter1++;
    }

#if 0
    path = build_path_end3(&knownsplice3p,&(*ambig_end_length_3),&(*ambig_splicetype_3),&(*ambig_prob_3),
			   &chop_exon_p,&dynprogindex_minor,path,
			   chroffset,chrhigh,querylength,
			   knownsplice_limit_low,knownsplice_limit_high,
			   queryseq_ptr,queryuc_ptr,
			   cdna_direction,watsonp,genestrand,jump_late_p,
			   maxpeelback,defect_rate,pairpool,dynprogL,/*extendp*/true,
			   /*endalign*/QUERYEND_NOGAPS,/*forcep*/false);
#endif
  }


#if 0
  debug(printf("Pass 11 (dir %d): Final extension, end3\n",cdna_direction));
  /* This is the second extension */
  /* Perform final extension without gaps so we can compare fwd and rev properly */
  path = build_path_end3(&knownsplice3p,&(*ambig_end_length_3),&(*ambig_splicetype_3),&(*ambig_prob_3),
			 &chop_exon_p,&dynprogindex_minor,path,
			 chroffset,chrhigh,querylength,
			 knownsplice_limit_low,knownsplice_limit_high,
			 queryseq_ptr,queryuc_ptr,
			 cdna_direction,watsonp,genestrand,jump_late_p,
			 maxpeelback,defect_rate,pairpool,dynprogL,
			 /*extendp*/true,/*endalign*/QUERYEND_NOGAPS,/*forcep*/false);
#endif

  debug(Pair_dump_list(path,true));
  debug(printf("End of path_compute_end3\n"));
  return path;
}


static List_T
path_compute_final (double defect_rate, List_T pairs, int cdna_direction, bool watsonp, int genestrand,
		    bool jump_late_p, int querylength,
#ifdef PMAP
		    char *queryaaseq_ptr,
#endif
		    char *queryseq_ptr, char *queryuc_ptr,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    int maxpeelback,
#ifndef GSNAP
		    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool, 
#endif
		    Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    Chrpos_T *last_genomedp5, Chrpos_T *last_genomedp3) {
  List_T path = NULL;
  int dynprogindex_minor = DYNPROGINDEX_MINOR, dynprogindex_major = DYNPROGINDEX_MAJOR;
  bool dual_break_p = true;
  bool shiftp, incompletep;

  /* List_T all_stage2_starts, all_stage2_ends; */
  /* Univcoord_T position; */


  debug(printf("Entering path_compute_final\n"));

#if 0
  if (circularp[chrnum] == true) {
    overall_end_distance = overall_end_distance_circular;
  } else {
    overall_end_distance = overall_end_distance_linear;
  }
#endif

  path = List_reverse(pairs);
  pairs = build_pairs_singles(&dynprogindex_minor,path,/*maxsize*/nullgap,
			      chroffset,chrhigh,queryseq_ptr,queryuc_ptr,querylength,watsonp,genestrand,
			      jump_late_p,maxpeelback,defect_rate,pairpool,dynprogM,
			      last_genomedp5,last_genomedp3,/*forcep*/true,/*finalp*/true);

#if 1
  /* Okay to use finalp == true, as long as Dynprog_genome_gap is called with finalp == false */
  debug(printf("\n*** Pass 999 (dir %d): Final pass to find canonical introns\n",cdna_direction));
  debug(Pair_dump_list(pairs,true));
  path = List_reverse(pairs);	/* ? insert_gapholders() */
  pairs = build_pairs_introns(&shiftp,&incompletep,
			      &dynprogindex_minor,&dynprogindex_major,path,
			      chrnum,chroffset,chrhigh,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      queryseq_ptr,queryuc_ptr,querylength,
			      cdna_direction,watsonp,genestrand,jump_late_p,
			      maxpeelback,defect_rate,pairpool,dynprogL,dynprogM,dynprogR,
#ifndef GSNAP
			      oligoindices_minor,diagpool,cellpool,
#endif
			      last_genomedp5,last_genomedp3,/*finalp*/true,/*simplep*/true);
#endif

  path = List_reverse(pairs);
  pairs = build_dual_breaks(&dual_break_p,&dynprogindex_minor,&dynprogindex_major,path,
			    chrnum,chroffset,chrhigh,
#ifdef PMAP
			    queryaaseq_ptr,
#endif
			    queryseq_ptr,queryuc_ptr,querylength,
			    cdna_direction,watsonp,genestrand,jump_late_p,pairpool,
			    dynprogL,dynprogM,dynprogR,last_genomedp5,last_genomedp3,maxpeelback,
#ifndef GSNAP
			    oligoindices_minor,diagpool,cellpool,
#endif
			    defect_rate,/*finalp*/true,/*simplep*/true);

  path = insert_gapholders(pairs,queryseq_ptr,queryuc_ptr,chroffset,chrhigh,watsonp,genestrand,
			   pairpool,/*finalp*/true);
  pairs = assign_gap_types(path,cdna_direction,watsonp,queryseq_ptr,
			   chrnum,chroffset,chrhigh,pairpool);

  return pairs;
}



#ifdef GSNAP
/* I believe this function never gets called with SENSE_NULL */
static List_T
trim_novel_spliceends (int *new_sensedir, List_T pairs, 
		       int *ambig_end_length_5, int *ambig_end_length_3,
		       Splicetype_T *ambig_splicetype_5, Splicetype_T *ambig_splicetype_3,
		       double *ambig_prob_5, double *ambig_prob_3,
		       int orig_sensedir, bool watsonp, int querylength,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       bool knownsplice5p, bool knownsplice3p) {
  List_T path, p;
  int i;

  Pair_T pair;
  int trim5, trim3, exondist;
  Univcoord_T genomicpos, start_genomicpos, middle_genomicpos, end_genomicpos;
  Univcoord_T splice_genomepos_5, splice_genomepos_3, splice_genomepos_5_mm, splice_genomepos_3_mm;
  Univcoord_T start, middle, end; /* start to middle has mismatches, while middle to end has none */
  double donor_prob, acceptor_prob;
  double max_prob_5 = 0.0, max_prob_3 = 0.0,
    max_prob_sense_forward_5 = 0.0, max_prob_sense_anti_5 = 0.0,
    max_prob_sense_forward_3 = 0.0, max_prob_sense_anti_3 = 0.0;
  double max_prob_5_mm = 0.0, max_prob_3_mm = 0.0,
    max_prob_sense_forward_5_mm = 0.0, max_prob_sense_anti_5_mm = 0.0,
    max_prob_sense_forward_3_mm = 0.0, max_prob_sense_anti_3_mm = 0.0;
  Splicetype_T splicetype5, splicetype3, splicetype5_mm, splicetype3_mm;
  /* int splice_cdna_direction_5, splice_cdna_direction_3; */
  int splice_sensedir_5, splice_sensedir_3;
  /* int splice_cdna_direction_5_mm, splice_cdna_direction_3_mm; */
  int splice_sensedir_5_mm, splice_sensedir_3_mm;
  int nmismatches, *scorei;
  bool mismatchp;

  debug13(printf("\nEntered gmap_trim_novel_spliceends with orig_sensedir %d, ambig_end_lengths %d and %d\n",
		 orig_sensedir,*ambig_end_length_5,*ambig_end_length_3));
  *new_sensedir = SENSE_NULL;

  Pair_trim_distances(&trim5,&trim3,pairs);
  debug13(printf("Trim distances are %d and %d\n",trim5,trim3));
  if (trim5 > trim3) {
    scorei = (int *) MALLOC((trim5 + 1) * sizeof(int));
  } else {
    scorei = (int *) MALLOC((trim3 + 1) * sizeof(int));
  }


  path = List_reverse(pairs);
  if (path != NULL && knownsplice3p == false && *ambig_end_length_3 == 0
      /* && exon_length_3(path) >= END_SPLICESITE_EXON_LENGTH */) {
    /* See if there is a good splice site at the 3' end */
    /* debug13(Pair_dump_list(path,true)); */

    pair = (Pair_T) List_head(p = path);
    start = middle = end = pair->genomepos;
    debug13(printf("Initializing start and end to be %u\n",start));

    if (pair->querypos != querylength - 1) {
      mismatchp = true;
    } else {
      mismatchp = false;
    }

    i = 0;
    nmismatches = 0;
    while (i < trim3) {
      if ((p = List_next(p)) == NULL) {
	break;
      } else if (pair->gapp == true) {
	break;
      } else if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	middle = pair->genomepos;
	scorei[i] = nmismatches;
	debug13(printf("Resetting middle to be %u\n",middle));
      } else {
	middle = pair->genomepos;
	scorei[i] = ++nmismatches;
	mismatchp = true;
	debug13(printf("Resetting middle to be %u\n",middle));
      }
      pair = (Pair_T) List_head(p);
      i++;
    }
    scorei[i] = ++nmismatches;

    while (i < trim3 + END_SPLICESITE_SEARCH) {
      if ((p = List_next(p)) == NULL) {
	break;
      } else if (pair->gapp == true) {
	break;
      } else {
	end = pair->genomepos;
	debug13(printf("Resetting end to be %u\n",end));
      }
      pair = (Pair_T) List_head(p);
      i++;
    }

    /* Find distance from end to intron, if any */
    exondist = 0;
    while (p != NULL && ((Pair_T) List_head(p))->gapp == false &&
	   exondist < END_MIN_EXONLENGTH) {
      p = List_next(p);
      exondist++;
    }
    debug13(printf("exondist is %d\n",exondist));

    if (mismatchp == false) {
      /* Allow perfect overhangs into intron */
      /* Note: pairs may not extend all the way to the end, which is why we look at end pair to initialize mismatchp */
      debug13(printf("Allowing perfect overhang into potential intron\n"));

    } else if (orig_sensedir == SENSE_FORWARD) {
      if (watsonp) {
	splicetype3 = splicetype3_mm = DONOR;

	start_genomicpos = start + 1;
	middle_genomicpos = middle + 1;
	end_genomicpos = end + 1;

	/* assert(start_genomicpos >= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos >= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("3', watson, sense anti %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos--;
	  i++;
	}
	while (genomicpos >= middle_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 1 */
	  debug13(printf("3', watson, sense anti %u %u %f mm %d\n",chroffset+genomicpos,genomicpos,donor_prob,nmismatches - scorei[i]));
	  if (donor_prob > max_prob_3_mm) {
	    max_prob_3_mm = donor_prob;
	    splice_genomepos_3_mm = genomicpos - 1;
	  }
	  genomicpos--;
	  debug13(i++);
	}
	while (genomicpos >= end_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 1 */
	  debug13(printf("3', watson, sense anti %u %u %f\n",chroffset+genomicpos,genomicpos,donor_prob));
	  if (donor_prob > max_prob_3) {
	    max_prob_3 = donor_prob;
	    splice_genomepos_3 = genomicpos - 1;
	  }
	  genomicpos--;
	}
	debug13(printf("\n"));

      } else {
	splicetype3 = splicetype3_mm = ANTIDONOR;

	start_genomicpos = (start > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - start;
	middle_genomicpos = (middle > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - middle;
	end_genomicpos = (end > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - end;

	/* assert(start_genomicpos <= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos <= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("3', crick, sense forward %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos++;
	  i++;
	}
	while (genomicpos <= middle_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 3 */
	  debug13(printf("3', crick, sense forward %u %u %f mm %d\n",chroffset+genomicpos,genomicpos,donor_prob,nmismatches - scorei[i]));
	  if (donor_prob > max_prob_3_mm) {
	    max_prob_3_mm = donor_prob;
	    splice_genomepos_3_mm = (chrhigh - chroffset) - genomicpos;
	  }
	  genomicpos++;
	  debug13(i++);
	}
	while (genomicpos <= end_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 3 */
	  debug13(printf("3', crick, sense forward %u %u %f\n",chroffset+genomicpos,genomicpos,donor_prob));
	  if (donor_prob > max_prob_3) {
	    max_prob_3 = donor_prob;
	    splice_genomepos_3 = (chrhigh - chroffset) - genomicpos;
	  }
	  genomicpos++;
	}
	debug13(printf("\n"));
      }

    } else if (orig_sensedir == SENSE_ANTI) {
      if (watsonp) {
	splicetype3 = splicetype3_mm = ANTIACCEPTOR;

	start_genomicpos = start + 1;
	middle_genomicpos = middle + 1;
	end_genomicpos = end + 1;

	/* assert(start_genomicpos >= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos >= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("3', watson, sense forward %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos--;
	  i++;
	}
	while (genomicpos >= middle_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 5 */
	  debug13(printf("3', watson, sense forward %u %u %f mm %d\n",chroffset+genomicpos,genomicpos,acceptor_prob,nmismatches - scorei[i]));
	  if (acceptor_prob > max_prob_3_mm) {
	    max_prob_3_mm = acceptor_prob;
	    splice_genomepos_3_mm = genomicpos - 1;
	  }
	  genomicpos--;
	  debug13(i++);
	}
	while (genomicpos >= end_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 5 */
	  debug13(printf("3', watson, sense forward %u %u %f\n",chroffset+genomicpos,genomicpos,acceptor_prob));
	  if (acceptor_prob > max_prob_3) {
	    max_prob_3 = acceptor_prob;
	    splice_genomepos_3 = genomicpos - 1;
	  }
	  genomicpos--;
	}
	debug13(printf("\n"));

      } else {
	splicetype3 = splicetype3_mm = ACCEPTOR;

	start_genomicpos = (start > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - start;
	middle_genomicpos = (middle > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - middle;
	end_genomicpos = (end > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - end;

	/* assert(start_genomicpos <= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos <= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("3', crick, sense anti %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos++;
	  i++;
	}
	while (genomicpos <= middle_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 7 */
	  debug13(printf("3', crick, sense anti %u %u %f mm %d\n",chroffset+genomicpos,genomicpos,acceptor_prob,nmismatches - scorei[i]));
	  if (acceptor_prob > max_prob_3_mm) {
	    max_prob_3_mm = acceptor_prob;
	    splice_genomepos_3_mm = (chrhigh - chroffset) - genomicpos;
	  }
	  genomicpos++;
	  debug13(i++);
	}
	while (genomicpos <= end_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 7 */
	  debug13(printf("3', crick, sense anti %u %u %f\n",chroffset+genomicpos,genomicpos,acceptor_prob));
	  if (acceptor_prob > max_prob_3) {
	    max_prob_3 = acceptor_prob;
	    splice_genomepos_3 = (chrhigh - chroffset) - genomicpos;
	  }
	  genomicpos++;
	}
	debug13(printf("\n"));
      }
      
    } else {
      if (watsonp) {
	start_genomicpos = start + 1;
	middle_genomicpos = middle + 1;
	end_genomicpos = end + 1;

	/* assert(start_genomicpos >= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos >= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("3', watson, sense null %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos--;
	  i++;
	}
	while (genomicpos >= middle_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 1 */
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 5 */
	  debug13(printf("3', watson, sense null %u %u %f %f mm %d\n",chroffset+genomicpos,genomicpos,donor_prob,acceptor_prob,nmismatches - scorei[i]));
	  if (donor_prob > max_prob_sense_forward_3_mm) {
	    max_prob_sense_forward_3_mm = donor_prob;
	    if (donor_prob > max_prob_3_mm) {
	      max_prob_3_mm = donor_prob;
	      splice_genomepos_3_mm = genomicpos - 1;
	      /* splice_cdna_direction_3_mm = +1; */
	      splice_sensedir_3_mm = SENSE_FORWARD;
	      splicetype3_mm = DONOR;
	    }
	  }
	  if (acceptor_prob > max_prob_sense_anti_3_mm) {
	    max_prob_sense_anti_3_mm = acceptor_prob;
	    if (acceptor_prob > max_prob_3_mm) {
	      max_prob_3_mm = acceptor_prob;
	      splice_genomepos_3_mm = genomicpos - 1;
	      /* splice_cdna_direction_3_mm = -1; */
	      splice_sensedir_3_mm = SENSE_ANTI;
	      splicetype3_mm = ANTIACCEPTOR;
	    }
	  }
	  genomicpos--;
	  debug13(i++);
	}
	while (genomicpos >= end_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 1 */
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 5 */
	  debug13(printf("3', watson, sense null %u %u %f %f\n",chroffset+genomicpos,genomicpos,donor_prob,acceptor_prob));
	  if (donor_prob > max_prob_sense_forward_3) {
	    max_prob_sense_forward_3 = donor_prob;
	    if (donor_prob > max_prob_3) {
	      max_prob_3 = donor_prob;
	      splice_genomepos_3 = genomicpos - 1;
	      /* splice_cdna_direction_3 = +1; */
	      splice_sensedir_3 = SENSE_FORWARD;
	      splicetype3 = DONOR;
	    }
	  }
	  if (acceptor_prob > max_prob_sense_anti_3) {
	    max_prob_sense_anti_3 = acceptor_prob;
	    if (acceptor_prob > max_prob_3) {
	      max_prob_3 = acceptor_prob;
	      splice_genomepos_3 = genomicpos - 1;
	      /* splice_cdna_direction_3 = -1; */
	      splice_sensedir_3 = SENSE_ANTI;
	      splicetype3 = ANTIACCEPTOR;
	    }
	  }
	  genomicpos--;
	}
	debug13(printf("\n"));

      } else {
	start_genomicpos = (start > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - start;
	middle_genomicpos = (middle > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - middle;
	end_genomicpos = (end > chrhigh - chroffset) ? 0 : (chrhigh - chroffset) - end;

	/* assert(start_genomicpos <= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos <= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("3', crick, sense null %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos++;
	  i++;
	}
	while (genomicpos <= middle_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 3 */
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 7 */
	  debug13(printf("3', crick, sense null %u %u %f %f mm %d\n",chroffset+genomicpos,genomicpos,donor_prob,acceptor_prob,nmismatches - scorei[i]));
	  if (donor_prob > max_prob_sense_forward_3_mm) {
	    max_prob_sense_forward_3_mm = donor_prob;
	    if (donor_prob > max_prob_3_mm) {
	      max_prob_3_mm = donor_prob;
	      splice_genomepos_3_mm = (chrhigh - chroffset) - genomicpos;
	      /* splice_cdna_direction_3_mm = +1; */
	      splice_sensedir_3_mm = SENSE_FORWARD;
	      splicetype3_mm = ANTIDONOR;
	    }
	  }
	  if (acceptor_prob > max_prob_sense_anti_3_mm) {
	    max_prob_sense_anti_3_mm = acceptor_prob;
	    if (acceptor_prob > max_prob_3_mm) {
	      max_prob_3_mm = acceptor_prob;
	      splice_genomepos_3_mm = (chrhigh - chroffset) - genomicpos;
	      /* splice_cdna_direction_3_mm = -1; */
	      splice_sensedir_3_mm = SENSE_ANTI;
	      splicetype3_mm = ACCEPTOR;
	    }
	  }
	  genomicpos++;
	  debug13(i++);
	}
	while (genomicpos <= end_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 3 */
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 7 */
	  debug13(printf("3', crick, sense null %u %u %f %f\n",chroffset+genomicpos,genomicpos,donor_prob,acceptor_prob));
	  if (donor_prob > max_prob_sense_forward_3) {
	    max_prob_sense_forward_3 = donor_prob;
	    if (donor_prob > max_prob_3) {
	      max_prob_3 = donor_prob;
	      splice_genomepos_3 = (chrhigh - chroffset) - genomicpos;
	      /* splice_cdna_direction_3 = +1; */
	      splice_sensedir_3 = SENSE_FORWARD;
	      splicetype3 = ANTIDONOR;
	    }
	  }
	  if (acceptor_prob > max_prob_sense_anti_3) {
	    max_prob_sense_anti_3 = acceptor_prob;
	    if (acceptor_prob > max_prob_3) {
	      max_prob_3 = acceptor_prob;
	      splice_genomepos_3 = (chrhigh - chroffset) - genomicpos;
	      /* splice_cdna_direction_3 = -1; */
	      splice_sensedir_3 = SENSE_ANTI;
	      splicetype3 = ACCEPTOR;
	    }
	  }
	  genomicpos++;
	}
	debug13(printf("\n"));
      }
    }

    if (orig_sensedir != SENSE_NULL) {
      if (max_prob_3 > END_SPLICESITE_PROB_MATCH) {
	debug13(printf("Found good splice %s on 3' end at %u with probability %f\n",
		       Splicetype_string(splicetype3),splice_genomepos_3,max_prob_3));
	while (path != NULL && ((Pair_T) path->first)->genomepos > splice_genomepos_3) {
	  path = Pairpool_pop(path,&pair);
	}
	/* path = clean_path_end3(path); -- Gives wrong end */
	if (path != NULL) {
	  *ambig_end_length_3 = (querylength - 1) - ((Pair_T) path->first)->querypos;
	  *ambig_splicetype_3 = splicetype3;
	  *ambig_prob_3 = max_prob_3;
	  debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	}

      } else if (max_prob_3_mm > END_SPLICESITE_PROB_MISMATCH) {
	debug13(printf("Found good mismatch splice %s on 3' end at %u with probability %f\n",
		       Splicetype_string(splicetype3_mm),splice_genomepos_3_mm,max_prob_3_mm));
	while (path != NULL && ((Pair_T) path->first)->genomepos > splice_genomepos_3_mm) {
	  path = Pairpool_pop(path,&pair);
	}
	/* path = clean_path_end3(path); -- Gives wrong end */
	if (path != NULL) {
	  *ambig_end_length_3 = (querylength - 1) - ((Pair_T) path->first)->querypos;
	  *ambig_splicetype_3 = splicetype3_mm;
	  *ambig_prob_3 = max_prob_3_mm;
	  debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	}
      }
    }
  }

  /* 5' end */

  pairs = List_reverse(path);
  if (pairs != NULL && knownsplice5p == false && *ambig_end_length_5 == 0
      /* && exon_length_5(pairs) >= END_SPLICESITE_EXON_LENGTH */) {
    /* See if there is a good splice site at the 5' end */
    /* debug13(Pair_dump_list(pairs,true)); */

    pair = (Pair_T) List_head(p = pairs);
    start = middle = end = pair->genomepos;
    debug13(printf("Initializing start and end to be %u\n",start));

    if (pair->querypos != 0) {
      mismatchp = true;
    } else {
      mismatchp = false;
    }

    i = 0;
    nmismatches = 0;
    while (i < trim5) {
      if ((p = List_next(p)) == NULL) {
	break;
      } else if (pair->gapp == true) {
	break;
      } else if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	middle = pair->genomepos;
	scorei[i] = nmismatches;
	debug13(printf("Resetting middle to be %u\n",middle));
      } else {
	middle = pair->genomepos;
	scorei[i] = ++nmismatches;
	mismatchp = true;
	debug13(printf("Resetting middle to be %u\n",middle));
      }
      pair = (Pair_T) List_head(p);
      i++;
    }
    scorei[i] = nmismatches;

    while (i < trim5 + END_SPLICESITE_SEARCH) {
      if ((p = List_next(p)) == NULL) {
	break;
      } else if (pair->gapp == true) {
	break;
      } else {
	end = pair->genomepos;
	debug13(printf("Resetting end to be %u\n",end));
      }
      pair = (Pair_T) List_head(p);
      i++;
    }

    /* Find distance from end to intron, if any */
    exondist = 0;
    while (p != NULL && ((Pair_T) List_head(p))->gapp == false &&
	   exondist < END_MIN_EXONLENGTH) {
      p = List_next(p);
      exondist++;
    }
    debug13(printf("exondist is %d\n",exondist));

    if (mismatchp == false) {
      /* Allow perfect overhangs into intron */
      /* Note: pairs may not extend all the way to the end, which is why we look at end pair to initialize mismatchp */
      debug13(printf("Allowing perfect overhang into potential intron\n"));

    } else if (orig_sensedir == SENSE_FORWARD) {
      if (watsonp) {
	splicetype5 = splicetype5_mm = ACCEPTOR;

	start_genomicpos = start;
	middle_genomicpos = middle;
	end_genomicpos = end;

	/* assert(start_genomicpos <= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos <= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("5', watson, sense forward %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos++;
	  i++;
	}
	while (genomicpos <= middle_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 2 */
	  debug13(printf("5', watson, sense forward %u %u %f mm %d\n",chroffset+genomicpos,genomicpos,acceptor_prob,nmismatches - scorei[i]));
	  if (acceptor_prob > max_prob_5_mm) {
	    max_prob_5_mm = acceptor_prob;
	    splice_genomepos_5_mm = genomicpos;
	  }
	  genomicpos++;
	  debug13(i++);
	}
	while (genomicpos <= end_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 2 */
	  debug13(printf("5', watson, sense forward %u %u %f\n",chroffset+genomicpos,genomicpos,acceptor_prob));
	  if (acceptor_prob > max_prob_5) {
	    max_prob_5 = acceptor_prob;
	    splice_genomepos_5 = genomicpos;
	  }
	  genomicpos++;
	}
	debug13(printf("\n"));

      } else {
	splicetype5 = splicetype5_mm = ANTIACCEPTOR;

	start_genomicpos = (start > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - start + 1;
	middle_genomicpos = (middle > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - middle + 1;
	end_genomicpos = (end > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - end + 1;

	/* assert(start_genomicpos >= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos >= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("5', crick, sense anti %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos--;
	  i++;
	}
	while (genomicpos >= middle_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 4 */
	  debug13(printf("5', crick, sense anti %u %u %f mm %d\n",chroffset+genomicpos,genomicpos,acceptor_prob,nmismatches - scorei[i]));
	  if (acceptor_prob > max_prob_5_mm) {
	    max_prob_5_mm = acceptor_prob;
	    splice_genomepos_5_mm = (chrhigh - chroffset) - genomicpos + 1;
	  }
	  genomicpos--;
	  debug13(i++);
	}
	while (genomicpos >= end_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 4 */
	  debug13(printf("5', crick, sense anti %u %u %f\n",chroffset+genomicpos,genomicpos,acceptor_prob));
	  if (acceptor_prob > max_prob_5) {
	    max_prob_5 = acceptor_prob;
	    splice_genomepos_5 = (chrhigh - chroffset) - genomicpos + 1;
	  }
	  genomicpos--;
	}
	debug13(printf("\n"));
      }

    } else if (orig_sensedir == SENSE_ANTI) {
      if (watsonp) {
	splicetype5 = splicetype5_mm = ANTIDONOR;

	start_genomicpos = start;
	middle_genomicpos = middle;
	end_genomicpos = end;
	
	/* assert(start_genomicpos <= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos <= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("5', watson, sense anti %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos++;
	  i++;
	}
	while (genomicpos <= middle_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 6 */
	  debug13(printf("5', watson, sense anti %u %u %f mm %d\n",chroffset+genomicpos,genomicpos,donor_prob,nmismatches - scorei[i]));
	  if (donor_prob > max_prob_5_mm) {
	    max_prob_5_mm = donor_prob;
	    splice_genomepos_5_mm = genomicpos;
	  }
	  genomicpos++;
	  debug13(i++);
	}
	while (genomicpos <= end_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 6 */
	  debug13(printf("5', watson, sense anti %u %u %f\n",chroffset+genomicpos,genomicpos,donor_prob));
	  if (donor_prob > max_prob_5) {
	    max_prob_5 = donor_prob;
	    splice_genomepos_5 = genomicpos;
	  }
	  genomicpos++;
	}
	debug13(printf("\n"));

      } else {
	splicetype5 = splicetype5_mm = DONOR;

	start_genomicpos = (start > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - start + 1;
	middle_genomicpos = (middle > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - middle + 1;
	end_genomicpos = (end > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - end + 1;

	/* assert(start_genomicpos >= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos >= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("5', crick, sense forward %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos--;
	  i++;
	}
	while (genomicpos >= middle_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 8 */
	  debug13(printf("5', crick, sense forward %u %u %f mm %d\n",chroffset+genomicpos,genomicpos,donor_prob,nmismatches - scorei[i]));
	  if (donor_prob > max_prob_5_mm) {
	    max_prob_5_mm = donor_prob;
	    splice_genomepos_5_mm = (chrhigh - chroffset) - genomicpos + 1;
	  }
	  genomicpos--;
	  debug13(i++);
	}
	while (genomicpos >= end_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 8 */
	  debug13(printf("5', crick, sense forward %u %u %f\n",chroffset+genomicpos,genomicpos,donor_prob));
	  if (donor_prob > max_prob_5) {
	    max_prob_5 = donor_prob;
	    splice_genomepos_5 = (chrhigh - chroffset) - genomicpos + 1;
	  }
	  genomicpos--;
	}
	debug13(printf("\n"));
      }
      
    } else {
      if (watsonp) {
	start_genomicpos = start;
	middle_genomicpos = middle;
	end_genomicpos = end;

	/* assert(start_genomicpos <= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos <= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("5', watson, sense null %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos++;
	  debug13(i++);
	}
	while (genomicpos <= middle_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 2 */
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 6 */
	  debug13(printf("5', watson, sense null %u %u %f %f mm %d\n",chroffset+genomicpos,genomicpos,donor_prob,acceptor_prob,nmismatches - scorei[i]));
	  if (acceptor_prob > max_prob_sense_forward_5_mm) {
	    max_prob_sense_forward_5_mm = acceptor_prob;
	    if (acceptor_prob > max_prob_5_mm) {
	      max_prob_5_mm = acceptor_prob;
	      splice_genomepos_5_mm = genomicpos;
	      /* splice_cdna_direction_5_mm = +1; */
	      splice_sensedir_5_mm = SENSE_FORWARD;
	      splicetype5_mm = ACCEPTOR;
	    }
	  }
	  if (donor_prob > max_prob_sense_anti_5_mm) {
	    max_prob_sense_anti_5_mm = donor_prob;
	    if (donor_prob > max_prob_5_mm) {
	      max_prob_5_mm = donor_prob;
	      splice_genomepos_5_mm = genomicpos;
	      /* splice_cdna_direction_5_mm = -1; */
	      splice_sensedir_5_mm = SENSE_ANTI;
	      splicetype5_mm = ANTIDONOR;
	    }
	  }
	  genomicpos++;
	  debug13(i++);
	}
	while (genomicpos <= end_genomicpos &&
	       genomicpos <= end_genomicpos + exondist - END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_acceptor_prob(chroffset+genomicpos,chroffset); /* Case 2 */
	  donor_prob = Maxent_hr_antidonor_prob(chroffset+genomicpos,chroffset); /* Case 6 */
	  debug13(printf("5', watson, sense null %u %u %f %f\n",chroffset+genomicpos,genomicpos,donor_prob,acceptor_prob));
	  if (acceptor_prob > max_prob_sense_forward_5) {
	    max_prob_sense_forward_5 = acceptor_prob;
	    if (acceptor_prob > max_prob_5) {
	      max_prob_5 = acceptor_prob;
	      splice_genomepos_5 = genomicpos;
	      /* splice_cdna_direction_5 = +1; */
	      splice_sensedir_5 = SENSE_FORWARD;
	      splicetype5 = ACCEPTOR;
	    }
	  }
	  if (donor_prob > max_prob_sense_anti_5) {
	    max_prob_sense_anti_5 = donor_prob;
	    if (donor_prob > max_prob_5) {
	      max_prob_5 = donor_prob;
	      splice_genomepos_5 = genomicpos;
	      /* splice_cdna_direction_5 = -1; */
	      splice_sensedir_5 = SENSE_ANTI;
	      splicetype5 = ANTIDONOR;
	    }
	  }
	  genomicpos++;
	}
	debug13(printf("\n"));

      } else {
	start_genomicpos = (start > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - start + 1;
	middle_genomicpos = (middle > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - middle + 1;
	end_genomicpos = (end > chrhigh - chroffset) ? 1 : (chrhigh - chroffset) - end + 1;

	/* assert(start_genomicpos >= end_genomicpos); */
	genomicpos = start_genomicpos;
	i = 0;
	while (genomicpos >= middle_genomicpos && scorei[i] < nmismatches - 2) {
	  debug13(printf("5', crick, sense null %u %u score %d\n",chroffset+genomicpos,genomicpos,nmismatches - scorei[i]));
	  genomicpos--;
	  i++;
	}
	while (genomicpos >= middle_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 4 */
	  donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 8 */
	  debug13(printf("5', crick, sense null %u %u %f %f mm %d\n",chroffset+genomicpos,genomicpos,donor_prob,acceptor_prob,nmismatches - scorei[i]));
	  if (acceptor_prob > max_prob_sense_forward_5_mm) {
	    max_prob_sense_forward_5_mm = acceptor_prob;
	    if (acceptor_prob > max_prob_5_mm) {
	      max_prob_5_mm = acceptor_prob;
	      splice_genomepos_5_mm = (chrhigh - chroffset) - genomicpos + 1;
	      /* splice_cdna_direction_5_mm = +1; */
	      splice_sensedir_5_mm = SENSE_FORWARD;
	      splicetype5_mm = ANTIACCEPTOR;
	    }
	  }
	  if (donor_prob > max_prob_sense_anti_5_mm) {
	    max_prob_sense_anti_5_mm = donor_prob;
	    if (donor_prob > max_prob_5_mm) {
	      max_prob_5_mm = donor_prob;
	      splice_genomepos_5_mm = (chrhigh - chroffset) - genomicpos + 1;
	      /* splice_cdna_direction_5_mm = -1; */
	      splice_sensedir_5_mm = SENSE_ANTI;
	      splicetype5_mm = DONOR;
	    }
	  }
	  genomicpos--;
	  debug13(i++);
	}
	while (genomicpos >= end_genomicpos &&
	       genomicpos >= end_genomicpos - exondist + END_MIN_EXONLENGTH) {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(chroffset+genomicpos,chroffset); /* Case 4 */
	  donor_prob = Maxent_hr_donor_prob(chroffset+genomicpos,chroffset); /* Case 8 */
	  debug13(printf("5', crick, sense null %u %u %f %f\n",chroffset+genomicpos,genomicpos,donor_prob,acceptor_prob));
	  if (acceptor_prob > max_prob_sense_forward_5) {
	    max_prob_sense_forward_5 = acceptor_prob;
	    if (acceptor_prob > max_prob_5) {
	      max_prob_5 = acceptor_prob;
	      splice_genomepos_5 = (chrhigh - chroffset) - genomicpos + 1;
	      /* splice_cdna_direction_5 = +1; */
	      splice_sensedir_5 = SENSE_FORWARD;
	      splicetype5 = ANTIACCEPTOR;
	    }
	  }
	  if (donor_prob > max_prob_sense_anti_5) {
	    max_prob_sense_anti_5 = donor_prob;
	    if (donor_prob > max_prob_5) {
	      max_prob_5 = donor_prob;
	      splice_genomepos_5 = (chrhigh - chroffset) - genomicpos + 1;
	      /* splice_cdna_direction_5 = -1; */
	      splice_sensedir_5 = SENSE_ANTI;
	      splicetype5 = DONOR;
	    }
	  }
	  genomicpos--;
	}
	debug13(printf("\n"));
      }
    }

    if (orig_sensedir != SENSE_NULL) {
      if (max_prob_5 > END_SPLICESITE_PROB_MATCH) {
	debug13(printf("Found good splice %s on 5' end at %u with probability %f\n",
		       Splicetype_string(splicetype5),splice_genomepos_5,max_prob_5));
	while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < splice_genomepos_5) {
	  pairs = Pairpool_pop(pairs,&pair);
	}
	/* pairs = clean_pairs_end5(pairs); -- gives wrong end */
	if (pairs != NULL) {
	  *ambig_end_length_5 = ((Pair_T) pairs->first)->querypos;
	  *ambig_splicetype_5 = splicetype5;
	  *ambig_prob_5 = max_prob_5;
	  debug13(printf("Set ambig_end_length_5 to be %d\n",*ambig_end_length_5));
	}
      } else if (max_prob_5_mm > END_SPLICESITE_PROB_MISMATCH) {
	debug13(printf("Found good mismatch splice %s on 5' end at %u with probability %f\n",
		       Splicetype_string(splicetype5_mm),splice_genomepos_5_mm,max_prob_5_mm));
	while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < splice_genomepos_5_mm) {
	  pairs = Pairpool_pop(pairs,&pair);
	}
	/* pairs = clean_pairs_end5(pairs); -- gives wrong end */
	if (pairs != NULL) {
	  *ambig_end_length_5 = ((Pair_T) pairs->first)->querypos;
	  *ambig_splicetype_5 = splicetype5_mm;
	  *ambig_prob_5 = max_prob_5_mm;
	  debug13(printf("Set ambig_end_length_5 to be %d\n",*ambig_end_length_5));
	}
      }
    }
  }

  if (orig_sensedir == SENSE_NULL) {
    if (max_prob_3 >= END_SPLICESITE_PROB_MATCH || max_prob_5 >= END_SPLICESITE_PROB_MATCH) {
      if (max_prob_3 >= END_SPLICESITE_PROB_MATCH && max_prob_5 >= END_SPLICESITE_PROB_MATCH
	  && max_prob_sense_forward_3 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_anti_3 < END_SPLICESITE_PROB_MATCH
	  && max_prob_sense_forward_5 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_anti_5 < END_SPLICESITE_PROB_MATCH) {
	/* Forward sense wins on both sides */

	debug13(printf("Found good splice %s on 5' end at %u with probability %f\n",
		       Splicetype_string(splicetype5),splice_genomepos_5,max_prob_5));
	while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < splice_genomepos_5) {
	  pairs = Pairpool_pop(pairs,&pair);
	}
	/* pairs = clean_pairs_end5(pairs); -- gives wrong end */
	if (pairs != NULL) {
	  *ambig_end_length_5 = ((Pair_T) pairs->first)->querypos;
	  *ambig_splicetype_5 = splicetype5;
	  *ambig_prob_5 = max_prob_5;
	  /* *cdna_direction = splice_cdna_direction_5; */
	  debug13(printf("Set ambig_end_length_5 to be %d\n",*ambig_end_length_5));
	}

	debug13(printf("Found good splice %s on 3' end at %u with probability %f\n",
		       Splicetype_string(splicetype3),splice_genomepos_3,max_prob_3));
	path = List_reverse(pairs);
	while (path != NULL && ((Pair_T) path->first)->genomepos > splice_genomepos_3) {
	  path = Pairpool_pop(path,&pair);
	}
	/* path = clean_path_end3(path); -- gives wrong end */
	if (path != NULL) {
	  *ambig_end_length_3 = (querylength - 1) - ((Pair_T) path->first)->querypos;
	  *ambig_splicetype_3 = splicetype3;
	  *ambig_prob_3 = max_prob_3;
	  /* *cdna_direction = splice_cdna_direction_3; */
	  debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	}
	*new_sensedir = SENSE_FORWARD;
	pairs = List_reverse(path);

      } else if (max_prob_3 >= END_SPLICESITE_PROB_MATCH && max_prob_5 >= END_SPLICESITE_PROB_MATCH
		 && max_prob_sense_anti_3 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_forward_3 < END_SPLICESITE_PROB_MATCH
		 && max_prob_sense_anti_5 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_forward_5 < END_SPLICESITE_PROB_MATCH) {

	/* Anti sense wins on both sides */
	debug13(printf("Found good splice %s on 5' end at %u with probability %f\n",
		       Splicetype_string(splicetype5),splice_genomepos_5,max_prob_5));
	while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < splice_genomepos_5) {
	  pairs = Pairpool_pop(pairs,&pair);
	}
	/* pairs = clean_pairs_end5(pairs); -- gives wrong end */
	if (pairs != NULL) {
	  *ambig_end_length_5 = ((Pair_T) pairs->first)->querypos;
	  *ambig_splicetype_5 = splicetype5;
	  *ambig_prob_5 = max_prob_5;
	  /* *cdna_direction = splice_cdna_direction_5; */
	  debug13(printf("Set ambig_end_length_5 to be %d\n",*ambig_end_length_5));
	}

	debug13(printf("Found good splice %s on 3' end at %u with probability %f\n",
		       Splicetype_string(splicetype3),splice_genomepos_3,max_prob_3));
	path = List_reverse(pairs);
	while (path != NULL && ((Pair_T) path->first)->genomepos > splice_genomepos_3) {
	  path = Pairpool_pop(path,&pair);
	}
	/* path = clean_path_end3(path); -- gives wrong end */
	if (path != NULL) {
	  *ambig_end_length_3 = (querylength - 1) - ((Pair_T) path->first)->querypos;
	  *ambig_splicetype_3 = splicetype3;
	  *ambig_prob_3 = max_prob_3;
	  /* *cdna_direction = splice_cdna_direction_3; */
	  debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	}
	*new_sensedir = SENSE_ANTI;
	pairs = List_reverse(path);

      } else if (max_prob_3 > max_prob_5) {
	/* Consider just 3' end */
	debug13(printf("Found good splice %s on 3' end at %u with probability %f\n",
		       Splicetype_string(splicetype3),splice_genomepos_3,max_prob_3));
	path = List_reverse(pairs);
	while (path != NULL && ((Pair_T) path->first)->genomepos > splice_genomepos_3) {
	  path = Pairpool_pop(path,&pair);
	}
	/* path = clean_path_end3(path); -- gives wrong end */
	if (path != NULL) {
	  *ambig_end_length_3 = (querylength - 1) - ((Pair_T) path->first)->querypos;
	  *ambig_splicetype_3 = splicetype3;
	  *ambig_prob_3 = max_prob_3;
	  /* *cdna_direction = splice_cdna_direction_3; */
	  debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	  if (max_prob_sense_forward_3 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_anti_3 < END_SPLICESITE_PROB_MATCH
	      && max_prob_sense_anti_5 < END_SPLICESITE_PROB_MATCH) {
	    *new_sensedir = splice_sensedir_3;
	  } else if (max_prob_sense_anti_3 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_forward_3 < END_SPLICESITE_PROB_MATCH
		     && max_prob_sense_forward_5 < END_SPLICESITE_PROB_MATCH) {
	    *new_sensedir = splice_sensedir_3;
	  } else {
	    /* Not enough evidence to set sensedir */
	  }
	}
	pairs = List_reverse(path);

      } else {
	/* Consider just 5' end */
	debug13(printf("Found good splice %s on 5' end at %u with probability %f\n",
		       Splicetype_string(splicetype5),splice_genomepos_5,max_prob_5));
	while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < splice_genomepos_5) {
	  pairs = Pairpool_pop(pairs,&pair);
	}
	/* pairs = clean_pairs_end5(pairs); -- gives wrong end */
	if (pairs != NULL) {
	  *ambig_end_length_5 = ((Pair_T) pairs->first)->querypos;
	  *ambig_splicetype_5 = splicetype5;
	  *ambig_prob_5 = max_prob_5;
	  /* *cdna_direction = splice_cdna_direction_5; */
	  debug13(printf("Set ambig_end_length_5 to be %d\n",*ambig_end_length_5));
	  if (max_prob_sense_forward_5 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_anti_5 < END_SPLICESITE_PROB_MATCH
	      && max_prob_sense_anti_3 < END_SPLICESITE_PROB_MATCH) {
	    *new_sensedir = splice_sensedir_5;
	  } else if (max_prob_sense_anti_5 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_forward_5 < END_SPLICESITE_PROB_MATCH
		     && max_prob_sense_forward_3 < END_SPLICESITE_PROB_MATCH) {
	    *new_sensedir = splice_sensedir_5;
	  } else {
	    /* Not enough evidence to set sensedir */
	  }
	}
      }

    } else if (max_prob_3_mm >= END_SPLICESITE_PROB_MISMATCH || max_prob_5_mm >= END_SPLICESITE_PROB_MISMATCH) {
      if (max_prob_3_mm > max_prob_5_mm) {
	debug13(printf("Found good mismatch splice %s on 3' end at %u with probability %f\n",
		       Splicetype_string(splicetype3_mm),splice_genomepos_3_mm,max_prob_3_mm));
	path = List_reverse(pairs);
	while (path != NULL && ((Pair_T) path->first)->genomepos > splice_genomepos_3_mm) {
	  path = Pairpool_pop(path,&pair);
	}
	/* path = clean_path_end3(path); -- gives wrong end */
	if (path != NULL) {
	  *ambig_end_length_3 = (querylength - 1) - ((Pair_T) path->first)->querypos;
	  *ambig_splicetype_3 = splicetype3_mm;
	  *ambig_prob_3 = max_prob_3_mm;
	  /* *cdna_direction = splice_cdna_direction_3_mm; */
	  debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	  if (max_prob_sense_forward_3_mm >= END_SPLICESITE_PROB_MISMATCH && max_prob_sense_anti_3_mm < END_SPLICESITE_PROB_MISMATCH
	      && max_prob_sense_anti_5_mm < END_SPLICESITE_PROB_MISMATCH) {
	    *new_sensedir = splice_sensedir_3_mm;
	  } else if (max_prob_sense_anti_3_mm >= END_SPLICESITE_PROB_MISMATCH && max_prob_sense_forward_3_mm < END_SPLICESITE_PROB_MISMATCH
		     && max_prob_sense_forward_5_mm < END_SPLICESITE_PROB_MISMATCH) {
	    *new_sensedir = splice_sensedir_3_mm;
	  } else {
	    /* Not enough evidence to set sensedir */
	  }
	}
	pairs = List_reverse(path);
      } else {
	debug13(printf("Found good mismatch splice %s on 5' end at %u with probability %f\n",
		       Splicetype_string(splicetype5_mm),splice_genomepos_5_mm,max_prob_5_mm));
	while (pairs != NULL && ((Pair_T) pairs->first)->genomepos < splice_genomepos_5_mm) {
	  pairs = Pairpool_pop(pairs,&pair);
	}
	/* pairs = clean_pairs_end5(pairs); -- gives wrong end */
	if (pairs != NULL) {
	  *ambig_end_length_5 = ((Pair_T) pairs->first)->querypos;
	  *ambig_splicetype_5 = splicetype5_mm;
	  *ambig_prob_5 = max_prob_5_mm;
	  /* *cdna_direction = splice_cdna_direction_5_mm; */
	  debug13(printf("Set ambig_end_length_5 to be %d\n",*ambig_end_length_5));
	  if (max_prob_sense_forward_5_mm >= END_SPLICESITE_PROB_MISMATCH && max_prob_sense_anti_5_mm < END_SPLICESITE_PROB_MISMATCH
	      && max_prob_sense_anti_3_mm < END_SPLICESITE_PROB_MISMATCH) {
	    *new_sensedir = splice_sensedir_5_mm;
	  } else if (max_prob_sense_anti_5_mm >= END_SPLICESITE_PROB_MISMATCH && max_prob_sense_forward_5_mm < END_SPLICESITE_PROB_MISMATCH
		     && max_prob_sense_forward_3_mm < END_SPLICESITE_PROB_MISMATCH) {
	    *new_sensedir = splice_sensedir_5_mm;
	  } else {
	    /* Not enough evidence to set sensedir */
	  }
	}
      }
    }
  }

  FREE(scorei);
  return pairs;
}
#endif



#if 0
/* Still somewhat buggy */
static List_T
trim_novel_spliceends_new (List_T pairs, 
			   int *ambig_end_length_5, int *ambig_end_length_3,
			   Splicetype_T *ambig_splicetype_5, Splicetype_T *ambig_splicetype_3,
			   double *ambig_prob_5, double *ambig_prob_3,
			   int *sensedir, bool watsonp, int querylength,
			   Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool knownsplice5p, bool knownsplice3p) {
  List_T path, p;
  int i;

  Pair_T pair, prev;
  int trim5, trim3, exondist5, exondist3;

 /* start to middle has mismatches, while middle to end has none */
  Univcoord_T start5, middle5, end5, start3, middle3, end3;
  Univcoord_T genomicstart5, genomicend3;
  bool solve5p, solve3p, mismatchp;


  debug13(printf("\nEntered gmap_trim_novel_spliceends with sensedir %d, ambig_end_lengths %d and %d\n",
		 *sensedir,*ambig_end_length_5,*ambig_end_length_3));

  debug13(Pair_dump_list(pairs,true));

  Pair_trim_distances(&trim5,&trim3,pairs);
  debug13(printf("Trim distances (where we would trim theoretically) are %d and %d\n",trim5,trim3));

  if (pairs != NULL && knownsplice5p == false && *ambig_end_length_5 == 0
      /* && exon_length_5(pairs) >= END_SPLICESITE_EXON_LENGTH */) {
    /* See if there is a good splice site at the 5' end */

    pair = (Pair_T) List_head(p = pairs);
    middle5 = (Univcoord_T) pair->genomepos - 1;

    if (pair->querypos != 0) {
      mismatchp = true;
    } else {
      mismatchp = false;
    }

    i = 0;
    while (i < trim5) {
      if ((p = List_next(p)) == NULL) {
	break;
      } else if (pair->gapp == true) {
	break;
      } else if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	middle5 = pair->genomepos - 1;
	debug13(printf("Resetting middle to be %u\n",middle5));
      } else {
	middle5 = pair->genomepos - 1;
	mismatchp = true;
	debug13(printf("Resetting middle to be %u\n",middle5));
      }
      prev = pair;
      pair = (Pair_T) List_head(p);
      i++;
    }

    end5 = middle5;
    while (i < trim5 + END_SPLICESITE_SEARCH) {
      if ((p = List_next(p)) == NULL) {
	break;
      } else if (pair->gapp == true) {
	break;
      } else {
	end5 = pair->genomepos - 1;
	debug13(printf("Resetting end to be %u\n",end5));
      }
      prev = pair;
      pair = (Pair_T) List_head(p);
      i++;
    }

    /* Determine genomicstart5 after indels hve been skipped */
    if (pair->gapp == true) {
      pair = prev;
    }
    if (watsonp) {
      debug13(printf("Plus: Setting genomicstart5 to be genomepos %u - querypos %d\n",
		     pair->genomepos,pair->querypos));
      genomicstart5 = pair->genomepos - pair->querypos;
    } else {
      debug13(printf("Minus: Setting genomicstart5 to be genomepos %u - querypos %d - 1\n",
		     pair->genomepos,pair->querypos));
      genomicstart5 = pair->genomepos - pair->querypos - 1;
    }
    if ((start5 = middle5 - END_SPLICESITE_SEARCH) < genomicstart5) {
      start5 = genomicstart5;
    }


    /* Find distance from end to intron, if any */
    exondist5 = 0;
    while (p != NULL && ((Pair_T) List_head(p))->gapp == false &&
	   exondist5 < END_MIN_EXONLENGTH) {
      p = List_next(p);
      exondist5++;
    }
    debug13(printf("exondist5 is %d\n",exondist5));
  }

  if (mismatchp == false) {
    solve5p = false;
    genomicstart5 = 0;

  } else if (watsonp) {
    solve5p = true;
    genomicstart5 = chroffset + genomicstart5;

    start5 = chroffset + start5;
    middle5 = chroffset + middle5;
    end5 = chroffset + end5;
    debug13(printf("\n2 Set end points for 5' trim to be %u..%u..%u, plusp %d\n",
		   start5 - chroffset,middle5 - chroffset,end5 - chroffset,watsonp));

  } else {
    solve5p = true;
    genomicstart5 = chrhigh - genomicstart5;

    start5 = chrhigh - start5;
    middle5 = chrhigh - middle5;
    end5 = chrhigh - end5;
    debug13(printf("\n2 Set end points for 5' trim to be %u..%u..%u, plusp %d\n",
		   start5 - chroffset,middle5 - chroffset,end5 - chroffset,watsonp));
  }


  path = List_reverse(pairs);
  if (path != NULL && knownsplice3p == false && *ambig_end_length_3 == 0
      /* && exon_length_3(path) >= END_SPLICESITE_EXON_LENGTH */) {
    /* See if there is a good splice site at the 3' end */
    /* debug13(Pair_dump_list(path,true)); */

    pair = (Pair_T) List_head(p = path);
    middle3 = (Univcoord_T) pair->genomepos + 1;

    if (pair->querypos != querylength - 1) {
      mismatchp = true;
    } else {
      mismatchp = false;
    }

    i = 0;
    while (i < trim3) {
      if ((p = List_next(p)) == NULL) {
	break;
      } else if (pair->gapp == true) {
	break;
      } else if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	middle3 = pair->genomepos + 1;
	debug13(printf("Resetting middle to be %u\n",middle3));
      } else {
	middle3 = pair->genomepos + 1;
	mismatchp = true;
	debug13(printf("Resetting middle to be %u\n",middle3));
      }
      prev = pair;
      pair = (Pair_T) List_head(p);
      i++;
    }

    end3 = middle3;
    while (i < trim3 + END_SPLICESITE_SEARCH) {
      if ((p = List_next(p)) == NULL) {
	break;
      } else if (pair->gapp == true) {
	break;
      } else {
	end3 = pair->genomepos + 1;
	debug13(printf("Resetting end to be %u\n",end3));
      }
      prev = pair;
      pair = (Pair_T) List_head(p);
      i++;
    }

    /* Determine genomicend3 after indels hve been skipped */
    if (pair->gapp == true) {
      pair = prev;
    }
    if (watsonp) {
      debug13(printf("Plus: Setting genomicend3 to be genomepos %u + (querylength %d - querypos %d)\n",
		     pair->genomepos,querylength,pair->querypos));
      genomicend3 = pair->genomepos + (querylength - pair->querypos);
    } else {
      debug13(printf("Minus: Setting genomicend3 to be genomepos %u + (querylength %d - 1 - querypos %d)\n",
		     pair->genomepos,querylength,pair->querypos));
      genomicend3 = pair->genomepos + (querylength - 1 - pair->querypos);
    }
    if ((start3 = middle3 + END_SPLICESITE_SEARCH) > genomicend3) {
      start3 = genomicend3;
    }


    /* Find distance from end to intron, if any */
    exondist3 = 0;
    while (p != NULL && ((Pair_T) List_head(p))->gapp == false &&
	   exondist3 < END_MIN_EXONLENGTH) {
      p = List_next(p);
      exondist3++;
    }
    debug13(printf("exondist3 is %d\n",exondist3));
  }

  if (mismatchp == false) {
    solve3p = false;
    genomicend3 = 0;

  } else if (watsonp) {
    solve3p = true;
    genomicend3 = chroffset + genomicend3;

    start3 = chroffset + start3;
    middle3 = chroffset + middle3;
    end3 = chroffset + end3;

  } else {
    solve3p = true;
    genomicend3 = chrhigh - genomicend3;
    start3 = chrhigh - start3;
    middle3 = chrhigh - middle3;
    end3 = chrhigh - end3;
  }
  debug13(printf("\n2 Set end points for 3' trim to be %u..%u..%u, plusp %d\n",
		 start3 - chroffset,middle3 - chroffset,end3 - chroffset,watsonp));

  Splice_trim_novel_spliceends(&(*ambig_end_length_5),&(*ambig_end_length_3),
			       &(*ambig_splicetype_5),&(*ambig_splicetype_3),
			       &(*ambig_prob_5),&(*ambig_prob_3),/*orig_sensedir*/*sensedir,
			       start5,middle5,end5,solve5p,start3,middle3,end3,solve3p,
			       genomicstart5,genomicend3,chroffset,/*plusp*/watsonp);

  while (path != NULL && ((Pair_T) path->first)->querypos > (querylength - 1) - *ambig_end_length_3) {
    path = Pairpool_pop(path,&pair);
  }

  pairs = List_reverse(path);
  while (pairs != NULL && ((Pair_T) pairs->first)->querypos < *ambig_end_length_5) {
    pairs = Pairpool_pop(pairs,&pair);
  }

  debug13(printf("Returning ambig_end_length_5 %d and ambig_end_length_3 %d, probs %f and %f\n",
		 *ambig_end_length_5,*ambig_end_length_3,*ambig_prob_5,*ambig_prob_3));
  debug13(Pair_dump_list(pairs,true));

  return pairs;
}
#endif



static List_T
path_trim (double defect_rate, int *ambig_end_length_5, int *ambig_end_length_3,
	   Splicetype_T *ambig_splicetype_5, Splicetype_T *ambig_splicetype_3,
	   double *ambig_prob_5, double *ambig_prob_3,
	   List_T pairs, int *cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
	   int querylength,
#ifdef GSNAP
	   int orig_sensedir,
#endif
	   char *queryseq_ptr, char *queryuc_ptr,
	   Univcoord_T chroffset, Univcoord_T chrhigh,
	   Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
	   int maxpeelback, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR) {
  List_T path = NULL;
  int dynprogindex_minor = DYNPROGINDEX_MINOR;
  bool chop_exon_p;
  bool knownsplice5p = false, knownsplice3p = false;
  bool trimp, trim5p, trim3p, trim5p_ignore, trim3p_ignore, indelp;
  int iter = 0;
#ifdef GSNAP
  int new_sensedir;
#endif

#ifdef GSNAP
  debug(printf("Entering path_trim with cdna_direction %d and sensedir %d\n",*cdna_direction,orig_sensedir));
  debug3(printf("Entering path_trim with cdna_direction %d and sensedir %d\n",*cdna_direction,orig_sensedir));
#else
  debug(printf("Entering path_trim with cdna_direction %d\n",*cdna_direction));
  debug3(printf("Entering path_trim with cdna_direction %d\n",*cdna_direction));
#endif

  debug3(Pair_dump_list(pairs,true));

#ifdef GSNAP
  if (novelsplicingp == true) {
    pairs = trim_novel_spliceends(&new_sensedir,pairs,&(*ambig_end_length_5),&(*ambig_end_length_3),
				  &(*ambig_splicetype_5),&(*ambig_splicetype_3),
				  &(*ambig_prob_5),&(*ambig_prob_3),
				  orig_sensedir,watsonp,querylength,
				  chroffset,chrhigh,knownsplice5p,knownsplice3p);
  }
#endif
  debug13(printf("After trim_novel_spliceends\n"));
  debug13(Pair_dump_list(pairs,true));


  if (pairs == NULL) {
    return (List_T) NULL;
  } else if (maximize_coverage_p == true) {
    /* Don't trim ends */
  } else {
    knownsplice5p = knownsplice3p = false;
    debug3(printf("Before Pair_trim_ends\n"));
    debug3(Pair_dump_list(pairs,true));
    debug3(printf("\n"));

    /* Done anyway within loop below */
    /* pairs = Pair_trim_ends(&trim5p,&trim3p,pairs); */
    trimp = trim5p = trim3p = true;
    
    debug3(printf("After Pair_trim_ends: trim5p %d, trim3p %d\n",trim5p,trim3p));
    debug3(Pair_dump_list(pairs,true));
    debug3(printf("\n"));

    while (iter++ < 3 && trimp == true) {
      trimp = false;
      /* Revised: Using QUERYEND_NOGAPS combined with Pair_trim_ends */
      /* Old: Extend with BEST_LOCAL to get right local (not global) answer,
	 and with maxpeelback == 0 to ensure we perform no peelback */
      if (trim5p == true) {
	debug3(printf("Extending at 5'\n"));
	/* This is the third and final extension */
	pairs = build_pairs_end5(&knownsplice5p,&(*ambig_end_length_5),&(*ambig_splicetype_5),&(*ambig_prob_5),
				 &chop_exon_p,&dynprogindex_minor,pairs,
				 chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
				 queryseq_ptr,queryuc_ptr,
				 *cdna_direction,watsonp,genestrand,jump_late_p,
				 maxpeelback,defect_rate,pairpool,dynprogR,
				 /*extendp*/true,/*endalign*/QUERYEND_NOGAPS,/*forcep*/false);
	pairs = trim_end5_exons(&indelp,&trim5p,*ambig_end_length_5,pairs,dynprogR,chroffset,chrhigh,
				queryseq_ptr,queryuc_ptr,querylength,*cdna_direction,watsonp,genestrand,
				jump_late_p,pairpool,defect_rate);
	if (indelp == true) {
	  pairs = trim_end5_indels(pairs,*ambig_end_length_5,dynprogR,chroffset,chrhigh,
				   queryseq_ptr,queryuc_ptr,
				   watsonp,genestrand,jump_late_p,pairpool,defect_rate);
	}
	if (trim5p == true) {
	  trimp = true;
	}
      }

      if (trim3p == true) {
	debug3(printf("Extending at 3'\n"));
	/* This is the third and final extension */
	path = List_reverse(pairs);
	path = build_path_end3(&knownsplice3p,&(*ambig_end_length_3),&(*ambig_splicetype_3),&(*ambig_prob_3),
			       &chop_exon_p,&dynprogindex_minor,path,
			       chroffset,chrhigh,querylength,
			       knownsplice_limit_low,knownsplice_limit_high,
			       queryseq_ptr,queryuc_ptr,
			       *cdna_direction,watsonp,genestrand,jump_late_p,
			       maxpeelback,defect_rate,pairpool,dynprogL,
			       /*extendp*/true,/*endalign*/QUERYEND_NOGAPS,/*forcep*/false);
	path = trim_end3_exons(&indelp,&trim3p,*ambig_end_length_3,path,dynprogL,chroffset,chrhigh,
			       queryseq_ptr,queryuc_ptr,querylength,
			       *cdna_direction,watsonp,genestrand,jump_late_p,pairpool,defect_rate);
	if (indelp == true) {
	  path = trim_end3_indels(path,*ambig_end_length_3,dynprogL,chroffset,chrhigh,
				  queryseq_ptr,queryuc_ptr,querylength,
				  watsonp,genestrand,jump_late_p,pairpool,defect_rate);
	}
	if (trim3p == true) {
	  trimp = true;
	}
	pairs = List_reverse(path);
      }

      /* Important to end the alignment with Pair_trim_ends, or else trimming will be faulty */
      /* Also, doing trimming within each loop yields better results in a small number of cases */
      pairs = Pair_trim_ends(&trim5p_ignore,&trim3p_ignore,pairs,*ambig_end_length_5,*ambig_end_length_3);
    }

    debug3(printf("After trim ends:\n"));
    debug3(Pair_dump_list(pairs,true));
  }

  if (watsonp) {
    path = List_reverse(pairs);
    path = clip_path_end3_chromosomal_bounds(path,chroffset,chrhigh);
    path = clean_path_end3_gap_indels(path);
    pairs = List_reverse(path);
  } else {
    pairs = clip_pairs_end5_chromosomal_bounds(pairs);
    pairs = clean_pairs_end5_gap_indels(pairs);
  }

  /* Cannot put trim_novel_spliceends here, which can generate an infinite loop in calling procedures */

  debug3(printf("Final result of path_trim: chroffset = %u, cdna_direction %d, new_sensedir %d\n",
		chroffset,*cdna_direction,new_sensedir));
  debug3(Pair_dump_list(pairs,true));
  debug3(printf("\n"));

  return pairs;
}


/* Using alloca for last_genomedp5 and last_genomedp3 can cause stack overflow */
Stage3middle_T
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
		       int sense_try) {
  int goodness_fwd, goodness_rev;
  Chrpos_T *last_genomedp5_fwd = NULL, *last_genomedp3_fwd = NULL, *last_genomedp5_rev = NULL, *last_genomedp3_rev = NULL;
  List_T pairs_fwd, pairs_rev, path_fwd, path_rev;
  double defect_rate_fwd, defect_rate_rev;
  

  /* stage2pairs = Stage2_middle(stage2); */
#if defined(DEBUG0) || defined(DEBUG11)
  if (watsonp == true) {
    printf("Stage 3: *** Starting stage 3 middle at chrnum #%d, chrstart %u, chrend %u, sense try %d)\n",
	   chrnum,((Pair_T) stage2pairs->first)->genomepos,
	   ((Pair_T) List_last_value(stage2pairs))->genomepos,sense_try);
  } else {
    printf("Stage 3: *** Starting stage 3 middle at chrnum #%d, chrstart %u, chrend %u, sense try %d)\n",
	   chrnum,(chrhigh - chroffset) - ((Pair_T) stage2pairs->first)->genomepos,
	   (chrhigh - chroffset) - ((Pair_T) List_last_value(stage2pairs))->genomepos,sense_try);
  }
#endif

#ifdef DEBUG
  if (watsonp == true) {
    printf("Stage 3: *** Stage3_compute_middle: Starting stage 3 at chrnum #%d, chrstart %u, chrend %u, sense try %d)\n",
	   chrnum,((Pair_T) stage2pairs->first)->genomepos,
	   ((Pair_T) List_last_value(stage2pairs))->genomepos,sense_try);
  } else {
    printf("Stage 3: *** Stage3_compute_middle: Starting stage 3 at chrnum #%d, chrstart %u, chrend %u, sense try %d)\n",
	   chrnum,(chrhigh - chroffset) - ((Pair_T) stage2pairs->first)->genomepos,
	   (chrhigh - chroffset) - ((Pair_T) List_last_value(stage2pairs))->genomepos,sense_try);
  }

  printf("stage2pairs\n");
  Pair_dump_list(stage2pairs,true);
#endif


#ifdef PMAP
  pairs_fwd = stage2pairs;
  pairs_rev = (List_T) NULL;
  /* do_final_p = true; */
#else
  if (splicingp == false) {
    pairs_fwd = stage2pairs;
    pairs_rev = (List_T) NULL;
  } else if (sense_try > 0) {
    pairs_fwd = stage2pairs;
    pairs_rev = (List_T) NULL;
  } else if (sense_try < 0) {
    pairs_fwd = (List_T) NULL;
    pairs_rev = stage2pairs;
  } else {
    /* sense_try == 0: Should try both even if no introns (cf, AA011563) */
    pairs_fwd = stage2pairs;
    pairs_rev = Pairpool_copy(stage2pairs,pairpool);
  }
#endif


  /* 1.  Middle */
  if (pairs_fwd == NULL) {
    path_fwd = (List_T) NULL;
  } else {
    last_genomedp5_fwd = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
    last_genomedp3_fwd = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
    debug(printf("*** Solve path_fwd\n"));
    path_fwd = path_compute_dir(&defect_rate_fwd,pairs_fwd,/*cdna_direction*/+1,
				watsonp,genestrand,jump_late_p,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				maxpeelback,
#ifndef GSNAP
				oligoindices_minor,diagpool,cellpool, 
#endif
				pairpool,dynprogL,dynprogM,dynprogR,
				last_genomedp5_fwd,last_genomedp3_fwd/*,clean:ends_p:true*/);
    FREE(last_genomedp3_fwd);
    FREE(last_genomedp5_fwd);
  }

  if (pairs_rev == NULL) {
    path_rev = (List_T) NULL;
  } else {
    last_genomedp5_rev = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
    last_genomedp3_rev = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
    debug(printf("*** Solve path_rev\n"));
    path_rev = path_compute_dir(&defect_rate_rev,pairs_rev,/*cdna_direction*/-1,
				watsonp,genestrand,jump_late_p,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				maxpeelback,
#ifndef GSNAP
				oligoindices_minor,diagpool,cellpool, 
#endif
				pairpool,dynprogL,dynprogM,dynprogR,
				last_genomedp5_rev,last_genomedp3_rev/*,clean_ends_p:true*/);
    FREE(last_genomedp5_rev);
    FREE(last_genomedp3_rev);
  }

  pairs_fwd = assign_gap_types(path_fwd,/*cdna_direction*/+1,watsonp,queryseq_ptr,
			       chrnum,chroffset,chrhigh,pairpool);
  pairs_rev = assign_gap_types(path_rev,/*cdna_direction*/-1,watsonp,queryseq_ptr,
			       chrnum,chroffset,chrhigh,pairpool);

  goodness_fwd = Pair_goodness_simple(pairs_fwd);
  goodness_rev = Pair_goodness_simple(pairs_rev);
  debug0(printf("goodness %d fwd, %d rev\n",goodness_fwd,goodness_rev));
  if (goodness_fwd >= goodness_rev) {
    return Stage3middle_new(goodness_fwd,defect_rate_fwd,defect_rate_rev,pairs_fwd,pairs_rev,
			    chrnum,chroffset,chrhigh,chrlength,watsonp,genestrand,
			    all_stage2_starts,all_stage2_ends);
  } else {
    return Stage3middle_new(goodness_rev,defect_rate_fwd,defect_rate_rev,pairs_fwd,pairs_rev,
			    chrnum,chroffset,chrhigh,chrlength,watsonp,genestrand,
			    all_stage2_starts,all_stage2_ends);
  }
}


/* Using alloca for last_genomedp5 and last_genomedp3 can cause stack overflow */
struct Pair_T *
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
		     int sense_filter, Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  struct Pair_T *pairarray1;
  List_T all_stage2_starts, all_stage2_ends, p;
  double defect_rate_fwd, defect_rate_rev, defect_rate_temp, defect_rate;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  /* Chrpos_T chrlength; */

  bool watsonp, jump_late_p;
  int genestrand;

  Univcoord_T position;
  Chrpos_T chrstart, chrend, overall_end_distance;
  Pair_T start, end;
  int querypos;

  Chrpos_T *last_genomedp5_fwd = NULL, *last_genomedp3_fwd = NULL, *last_genomedp5_rev = NULL, *last_genomedp3_rev = NULL;
  List_T pairs_pretrim, pairs_fwd, pairs_rev, best_pairs, temp_pairs, path_fwd, path_rev, best_path, temp_path;
  List_T copy;
  List_T joined_ends, joined_starts;
  int nknown_fwd, ncanonical_fwd, nsemicanonical_fwd, nnoncanonical_fwd,
    nknown_rev, ncanonical_rev, nsemicanonical_rev, nnoncanonical_rev;
  int nbadintrons_fwd, nbadintrons_rev;
  double min_splice_prob_1;
  double max_intron_score_fwd = 0.0, max_intron_score_rev = 0.0,
    avg_donor_score_fwd = 0.0, avg_acceptor_score_fwd = 0.0,
    avg_donor_score_rev = 0.0, avg_acceptor_score_rev = 0.0;
  int nmatches_fwd, nmismatches_fwd, nmatches_rev, nmismatches_rev, nindels_fwd, nindels_rev;
  int fwd_ambig_end_length_5 = 0, fwd_ambig_end_length_3 = 0, rev_ambig_end_length_5 = 0, rev_ambig_end_length_3 = 0, temp_ambig_end_length;
  Splicetype_T fwd_ambig_splicetype_5, fwd_ambig_splicetype_3, rev_ambig_splicetype_5, rev_ambig_splicetype_3, temp_ambig_splicetype;
  double fwd_ambig_prob_5, fwd_ambig_prob_3, rev_ambig_prob_5, rev_ambig_prob_3, temp_ambig_prob;
#ifdef GSNAP
  List_T pairs_fwd_copy, pairs_rev_copy;
  double min_splice_prob_2;
#endif


  last_genomedp5_fwd = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
  last_genomedp3_fwd = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
  last_genomedp5_rev = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
  last_genomedp3_rev = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));

  chrnum = stage3middle->chrnum;
  chroffset = stage3middle->chroffset;
  chrhigh = stage3middle->chrhigh;
  /* chrlength = stage3middle->chrlength; */

  watsonp = stage3middle->watsonp;
  genestrand = stage3middle->genestrand;
  jump_late_p = watsonp ? false : true;

  defect_rate_fwd = stage3middle->defect_rate_fwd;
  defect_rate_rev = stage3middle->defect_rate_rev;

  pairs_fwd = stage3middle->pairs_fwd;
  pairs_rev = stage3middle->pairs_rev;
  path_fwd = List_reverse(pairs_fwd);
  path_rev = List_reverse(pairs_rev);

  if (path_fwd != NULL && path_rev != NULL) {
    debug11(printf("Calling score_introns for path_fwd\n"));
    pairs_fwd = score_introns(&max_intron_score_fwd,&avg_donor_score_fwd,&avg_acceptor_score_fwd,
			      &nknown_fwd,&ncanonical_fwd,&nbadintrons_fwd,path_fwd,/*cdna_direction*/+1,watsonp,
			      chrnum,chroffset,chrhigh
#ifdef WASTE
			      ,pairpool
#endif
			      );

    debug11(printf("Calling score_introns for path_rev\n"));
    pairs_rev = score_introns(&max_intron_score_rev,&avg_donor_score_rev,&avg_acceptor_score_rev,
			      &nknown_rev,&ncanonical_rev,&nbadintrons_rev,path_rev,/*cdna_direction*/-1,watsonp,
			      chrnum,chroffset,chrhigh
#ifdef WASTE
			      ,pairpool
#endif			      
			      );
    
    if ((*cdna_direction = initial_cdna_direction(pairs_fwd,pairs_rev,
						  avg_donor_score_fwd,avg_acceptor_score_fwd,
						  avg_donor_score_rev,avg_acceptor_score_rev)) > 0) {
      debug11(printf("Initial cdna direction is %d\n",*cdna_direction));
      path_fwd = List_reverse(pairs_fwd);
      path_rev = (List_T) NULL;

    } else if (*cdna_direction < 0) {
      debug11(printf("Initial cdna direction is %d\n",*cdna_direction));
      path_fwd = (List_T) NULL;
      path_rev = List_reverse(pairs_rev);

    } else {
      debug11(printf("Initial cdna direction is %d\n",*cdna_direction));
      path_fwd = List_reverse(pairs_fwd);
      path_rev = List_reverse(pairs_rev);
    }
  }



  if (circularp[chrnum] == true) {
    overall_end_distance = overall_end_distance_circular;
  } else {
    overall_end_distance = overall_end_distance_linear;
  }

  /* 2.  3' and 5' ends (possibly multiple) */
  if (path_fwd == NULL) {
    pairs_fwd = (List_T) NULL;
  } else {
    /* 3' end */
    end = (Pair_T) List_head(path_fwd);
    if ((querypos = end->querypos) >= querylength - 10) {
      all_stage2_ends = (List_T) NULL;

    } else {
      if (watsonp == true) {
	chrstart = end->genomepos;
	
	/* We actually don't care if alignment goes past chrhigh,
	   since that can be fixed later by trimming.  But we do care
	   about going past the genome bounds */
	if ((position = chroffset + chrstart + overall_end_distance) > genome_totallength) {
	  chrend = (Chrpos_T) (genome_totallength - chroffset);
	} else {
	  chrend = (Chrpos_T) (position - chroffset);
	}
	
      } else {
	chrend = (chrhigh - chroffset) - end->genomepos;
      
	/* This check is necessary, because we cannot have negative values for chrstart */
	if (chrend < overall_end_distance) {
	  chrstart = 0;
	} else {
	  chrstart = chrend - overall_end_distance;
	}
      }
      
      debug(printf("QUERYPOS %d, CHRSTART %u, CHREND %u\n",querypos,chrstart,chrend));
      all_stage2_ends = Stage2_compute_ends(&(queryseq_ptr[querypos]),&(queryuc_ptr[querypos]),
					    /*querylength*/querylength - querypos,/*query_offset*/querypos,
					    chrstart,chrend,
					    chroffset,chrhigh,/*plusp*/watsonp,genestrand,
					    oligoindices_minor,pairpool,diagpool,cellpool,/*localp*/false,
					    /*skip_repetitive_p*/false,/*favor_right_p*/false,/*max_nalignments*/1,
					    /*debug_graphic_p*/false);
#ifdef DEBUG
      printf("fwd stage2ends\n");
      for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
	Pair_dump_list((List_T) List_head(p),true);
      }
#endif
    }

    if (all_stage2_ends == NULL) {
      best_path = path_compute_end3(&fwd_ambig_end_length_3,&fwd_ambig_splicetype_3,&fwd_ambig_prob_3,
				    defect_rate_fwd,path_fwd,/*cdna_direction*/+1,watsonp,genestrand,
				    jump_late_p,querylength,
				    queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				    knownsplice_limit_low,knownsplice_limit_high,
				    maxpeelback,pairpool,dynprogL);
    } else {
      best_path = Pairpool_remove_gapholders(path_fwd); /* Pairpool_join cannot handle gapholders */
      joined_ends = (List_T) NULL;
      for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
#ifdef PMAP      
        copy = Pairpool_join_end3(/*path*/path_fwd,/*end3_pairs*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
#else
	if (path_rev == NULL) {
	  /* Won't need ends anymore */
	  copy = Pairpool_join_end3(/*path*/path_fwd,/*end3_pairs*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
	} else {
	  copy = Pairpool_join_end3(/*path*/path_fwd,/*end3_pairs*/(List_T) List_head(p),pairpool,/*copy_end_p*/true);
	}
#endif
	joined_ends = List_push(joined_ends,(void *) copy);
      }

      for (p = joined_ends; p != NULL; p = List_next(p)) {
        copy = (List_T) List_head(p);
	debug(printf("*** Solve path_fwd joined end\n"));
	path_fwd = path_compute_dir(&defect_rate_temp,/*pairs*/List_reverse(copy),/*cdna_direction*/+1,
				    watsonp,genestrand,jump_late_p,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				    maxpeelback,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool, 
#endif
				    pairpool,dynprogL,dynprogM,dynprogR,
				    last_genomedp5_fwd,last_genomedp3_fwd/*,clean_ends_p:false*/);

	temp_path = path_compute_end3(&temp_ambig_end_length,&temp_ambig_splicetype,&temp_ambig_prob,
				      defect_rate_temp,path_fwd,/*cdna_direction*/+1,watsonp,genestrand,
				      jump_late_p,querylength,
				      queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				      knownsplice_limit_low,knownsplice_limit_high,
				      maxpeelback,pairpool,dynprogL);

	if (temp_path != NULL && end_compare(best_path,temp_path,/*cdna_direction*/+1,watsonp,
					     chrnum,chroffset,chrhigh,/*pairsp*/false) > 0) {
	  best_path = temp_path;
	  fwd_ambig_end_length_3 = temp_ambig_end_length;
	  fwd_ambig_splicetype_3 = temp_ambig_splicetype;
	  fwd_ambig_prob_3 = temp_ambig_prob;
	  defect_rate_fwd = defect_rate_temp;
	  debug21(printf("New best path:\n"));
	  debug21(Pair_dump_list(best_path,true));
	}
      }

      List_free(&joined_ends);
      List_free(&all_stage2_ends);
    }

    /* 5' end */
    pairs_fwd = List_reverse(best_path);

    start = (Pair_T) List_head(pairs_fwd);
    if ((querypos = start->querypos) <= 10) {
      all_stage2_starts = (List_T) NULL;

    } else {
      if (watsonp == true) {
	chrend = start->genomepos;
	/* This check is necessary, because we cannot have negative values for chrstart */
	if (chrend < overall_end_distance) {
	  chrstart = 0;
	} else {
	  chrstart = chrend - overall_end_distance;
	}

      } else {
	chrstart = (chrhigh - chroffset) - start->genomepos;

	/* We actually don't care if alignment goes past chrhigh,
	   since that can be fixed later by trimming.  But we do care
	   about going past the genome bounds */
	if ((position = chroffset + chrstart + overall_end_distance) > genome_totallength) {
	  chrend = (Chrpos_T) (genome_totallength - chroffset);
	} else {
	  chrend = (Chrpos_T) (position - chroffset);
	}
      }
      debug(printf("QUERYPOS %d, CHRSTART %u, CHREND %u\n",querypos,chrstart,chrend));
      all_stage2_starts = Stage2_compute_starts(&(queryseq_ptr[0]),&(queryuc_ptr[0]),querypos,/*query_offset*/0,
						chrstart,chrend,
						chroffset,chrhigh,/*plusp*/watsonp,genestrand,
						oligoindices_minor,pairpool,diagpool,cellpool,/*localp*/false,
						/*skip_repetitive_p*/false,/*favor_right_p*/false,/*max_nalignments*/1,
						/*debug_graphic_p*/false);
#ifdef DEBUG
      printf("fwd stage2starts\n");
      for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
	Pair_dump_list((List_T) List_head(p),true);
      }
#endif
    }

    if (all_stage2_starts == NULL) {
      best_pairs = path_compute_end5(&fwd_ambig_end_length_5,&fwd_ambig_splicetype_5,&fwd_ambig_prob_5,
				     defect_rate_fwd,pairs_fwd,/*cdna_direction*/+1,
				     watsonp,genestrand,jump_late_p,
				     queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				     knownsplice_limit_low,knownsplice_limit_high,
				     maxpeelback,pairpool,dynprogR);
    } else {
      best_pairs = Pairpool_remove_gapholders(pairs_fwd); /* Pairpool_join cannot handle gapholders */
      joined_starts = (List_T) NULL;
      for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
#ifdef PMAP
        copy = Pairpool_join_end5(/*pairs*/pairs_fwd,/*end5_path*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
#else
	if (path_rev == NULL) {
	  /* Won't need ends anymore */
	  copy = Pairpool_join_end5(/*pairs*/pairs_fwd,/*end5_path*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
	} else {
	  copy = Pairpool_join_end5(/*pairs*/pairs_fwd,/*end5_path*/(List_T) List_head(p),pairpool,/*copy_end_p*/true);
	}
#endif
	joined_starts = List_push(joined_starts,(void *) copy);
      }

      for (p = joined_starts; p != NULL; p = List_next(p)) {
	copy = (List_T) List_head(p);
	debug(printf("*** Solve path_fwd joined start\n"));
	path_fwd = path_compute_dir(&defect_rate_temp,/*pairs*/copy,/*cdna_direction*/+1,
				    watsonp,genestrand,jump_late_p,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				    maxpeelback,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool, 
#endif
				    pairpool,dynprogL,dynprogM,dynprogR,
				    last_genomedp5_fwd,last_genomedp3_fwd/*,clean_ends_p:false*/);

	temp_pairs = path_compute_end5(&temp_ambig_end_length,&temp_ambig_splicetype,&temp_ambig_prob,
				       defect_rate_temp,/*pairs*/List_reverse(path_fwd),
				       /*cdna_direction*/+1,watsonp,genestrand,jump_late_p,
				       queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				       knownsplice_limit_low,knownsplice_limit_high,
				       maxpeelback,pairpool,dynprogR);

	if (temp_pairs != NULL && end_compare(best_pairs,temp_pairs,/*cdna_direction*/+1,watsonp,
					      chrnum,chroffset,chrhigh,/*pairsp*/true) > 0) {
	  best_pairs = temp_pairs;
	  fwd_ambig_end_length_5 = temp_ambig_end_length;
	  fwd_ambig_splicetype_5 = temp_ambig_splicetype;
	  fwd_ambig_prob_5 = temp_ambig_prob;
	  defect_rate_fwd = defect_rate_temp;
	  debug21(printf("New best pairs:\n"));
	  debug21(Pair_dump_list(best_pairs,true));
	}
      }

      List_free(&joined_starts);
      List_free(&all_stage2_starts);
    }

    pairs_fwd = best_pairs;
  }


#ifndef PMAP
  if (path_rev == NULL) {
    pairs_rev = (List_T) NULL;
  } else {
    /* 3' end */
    end = (Pair_T) List_head(path_rev);
    if ((querypos = end->querypos) >= querylength - 10) {
      all_stage2_ends = (List_T) NULL;

    } else {
      if (watsonp == true) {
	chrstart = end->genomepos;
	
	/* We actually don't care if alignment goes past chrhigh,
	   since that can be fixed later by trimming.  But we do care
	   about going past the genome bounds */
	if ((position = chroffset + chrstart + overall_end_distance) > genome_totallength) {
	  chrend = (Chrpos_T) (genome_totallength - chroffset);
	} else {
	  chrend = (Chrpos_T) (position - chroffset);
	}
	
      } else {
	chrend = (chrhigh - chroffset) - end->genomepos;
      
	/* This check is necessary, because we cannot have negative values for chrstart */
	if (chrend < overall_end_distance) {
	  chrstart = 0;
	} else {
	  chrstart = chrend - overall_end_distance;
	}
      }
      
      debug(printf("QUERYPOS %d, CHRSTART %u, CHREND %u\n",querypos,chrstart,chrend));
      all_stage2_ends = Stage2_compute_ends(&(queryseq_ptr[querypos]),&(queryuc_ptr[querypos]),
					    /*querylength*/querylength - querypos,/*query_offset*/querypos,
					    chrstart,chrend,
					    chroffset,chrhigh,/*plusp*/watsonp,genestrand,
					    oligoindices_minor,pairpool,diagpool,cellpool,/*localp*/false,
					    /*skip_repetitive_p*/false,/*favor_right_p*/false,/*max_nalignments*/1,
					    /*debug_graphic_p*/false);
#ifdef DEBUG
      printf("rev stage2ends\n");
      for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
	Pair_dump_list((List_T) List_head(p),true);
      }
#endif
    }

    if (all_stage2_ends == NULL) {
      best_path = path_compute_end3(&rev_ambig_end_length_3,&rev_ambig_splicetype_3,&rev_ambig_prob_3,
				    defect_rate_rev,path_rev,/*cdna_direction*/-1,watsonp,genestrand,
				    jump_late_p,querylength,
				    queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				    knownsplice_limit_low,knownsplice_limit_high,
				    maxpeelback,pairpool,dynprogL);

    } else {
      best_path = Pairpool_remove_gapholders(path_rev); /* Pairpool_join cannot handle gapholders */
      joined_ends = (List_T) NULL;
      for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
	copy = Pairpool_join_end3(/*path*/path_rev,/*end3_pairs*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
	joined_ends = List_push(joined_ends,(void *) copy);
      }

      for (p = joined_ends; p != NULL; p = List_next(p)) {
        copy = (List_T) List_head(p);
	debug(printf("*** Solve path_rev joined end\n"));
	path_rev = path_compute_dir(&defect_rate_temp,/*pairs*/List_reverse(copy),/*cdna_direction*/-1,
				    watsonp,genestrand,jump_late_p,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				    maxpeelback,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool, 
#endif
				    pairpool,dynprogL,dynprogM,dynprogR,
				    last_genomedp5_rev,last_genomedp3_rev/*,clean_ends_p:false*/);
	
	temp_path = path_compute_end3(&temp_ambig_end_length,&temp_ambig_splicetype,&temp_ambig_prob,
				      defect_rate_temp,path_rev,/*cdna_direction*/-1,watsonp,genestrand,
				      jump_late_p,querylength,
				      queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				      knownsplice_limit_low,knownsplice_limit_high,
				      maxpeelback,pairpool,dynprogL);
	
	if (temp_path != NULL && end_compare(best_path,temp_path,/*cdna_direction*/-1,watsonp,
					     chrnum,chroffset,chrhigh,/*pairsp*/false) > 0) {
	  best_path = temp_path;
	  rev_ambig_end_length_3 = temp_ambig_end_length;
	  rev_ambig_splicetype_3 = temp_ambig_splicetype;
	  rev_ambig_prob_3 = temp_ambig_prob;
	  defect_rate_rev = defect_rate_temp;
	  debug21(printf("New best path:\n"));
	  debug21(Pair_dump_list(best_path,true));
	}
      }

      List_free(&joined_ends);
      List_free(&all_stage2_ends);
    }

    /* 5' end */
    pairs_rev = List_reverse(best_path);

    start = (Pair_T) List_head(pairs_rev);
    if ((querypos = start->querypos) <= 10) {
      all_stage2_starts = (List_T) NULL;

    } else {
      if (watsonp == true) {
	chrend = start->genomepos;
	/* This check is necessary, because we cannot have negative values for chrstart */
	if (chrend < overall_end_distance) {
	  chrstart = 0;
	} else {
	  chrstart = chrend - overall_end_distance;
	}

      } else {
	chrstart = (chrhigh - chroffset) - start->genomepos;

	/* We actually don't care if alignment goes past chrhigh,
	   since that can be fixed later by trimming.  But we do care
	   about going past the genome bounds */
	if ((position = chroffset + chrstart + overall_end_distance) > genome_totallength) {
	  chrend = (Chrpos_T) (genome_totallength - chroffset);
	} else {
	  chrend = (Chrpos_T) (position - chroffset);
	}
      }
      debug(printf("QUERYPOS %d, CHRSTART %u, CHREND %u\n",querypos,chrstart,chrend));
      all_stage2_starts = Stage2_compute_starts(&(queryseq_ptr[0]),&(queryuc_ptr[0]),querypos,/*query_offset*/0,
						chrstart,chrend,
						chroffset,chrhigh,/*plusp*/watsonp,genestrand,
						oligoindices_minor,pairpool,diagpool,cellpool,/*localp*/false,
						/*skip_repetitive_p*/false,/*favor_right_p*/false,/*max_nalignments*/1,
						/*debug_graphic_p*/false);
#ifdef DEBUG
      printf("rev stage2starts\n");
      for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
	Pair_dump_list((List_T) List_head(p),true);
      }
#endif
    }

    if (all_stage2_starts == NULL) {
      best_pairs = path_compute_end5(&rev_ambig_end_length_5,&rev_ambig_splicetype_5,&rev_ambig_prob_5,
				     defect_rate_rev,pairs_rev,/*cdna_direction*/-1,
				     watsonp,genestrand,jump_late_p,
				     queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				     knownsplice_limit_low,knownsplice_limit_high,
				     maxpeelback,pairpool,dynprogR);

    } else {
      best_pairs = Pairpool_remove_gapholders(pairs_rev); /* Pairpool_join cannot handle gapholders */
      joined_starts = (List_T) NULL;
      for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
	copy = Pairpool_join_end5(/*pairs*/pairs_rev,/*end5_path*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
	joined_starts = List_push(joined_starts,(void *) copy);
      }

      for (p = joined_starts; p != NULL; p = List_next(p)) {
	copy = (List_T) List_head(p);
	debug(printf("*** Solve path_rev joined start\n"));
	path_rev = path_compute_dir(&defect_rate_temp,/*pairs*/copy,/*cdna_direction*/-1,
				    watsonp,genestrand,jump_late_p,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				    maxpeelback,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool, 
#endif
				    pairpool,dynprogL,dynprogM,dynprogR,
				    last_genomedp5_rev,last_genomedp3_rev/*,clean_ends_p:false*/);
      
	temp_pairs = path_compute_end5(&temp_ambig_end_length,&temp_ambig_splicetype,&temp_ambig_prob,
				       defect_rate_temp,/*pairs*/List_reverse(path_rev),
				       /*cdna_direction*/-1,watsonp,genestrand,jump_late_p,
				       queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				       knownsplice_limit_low,knownsplice_limit_high,
				       maxpeelback,pairpool,dynprogR);
	if (temp_pairs != NULL && end_compare(best_pairs,temp_pairs,/*cdna_direction*/-1,watsonp,
					      chrnum,chroffset,chrhigh,/*pairsp*/true) > 0) {
	  best_pairs = temp_pairs;
	  rev_ambig_end_length_5 = temp_ambig_end_length;
	  rev_ambig_splicetype_5 = temp_ambig_splicetype;
	  rev_ambig_prob_5 = temp_ambig_prob;
	  defect_rate_rev = defect_rate_temp;
	  debug21(printf("New best pairs:\n"));
	  debug21(Pair_dump_list(best_pairs,true));
	}
      }

      List_free(&joined_starts);
      List_free(&all_stage2_starts);
    }

    pairs_rev = best_pairs;
  }
#endif


  pairs_fwd = path_compute_final(defect_rate_fwd,pairs_fwd,/*cdna_direction*/+1,
				 watsonp,genestrand,jump_late_p,querylength,
#ifdef PMAP
				 queryaaseq_ptr,
#endif
				 queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				 maxpeelback,
#ifndef GSNAP
				 oligoindices_minor,diagpool,cellpool, 
#endif
				 pairpool,dynprogL,dynprogM,dynprogR,
				 last_genomedp5_fwd,last_genomedp3_fwd);

  pairs_rev = path_compute_final(defect_rate_rev,pairs_rev,/*cdna_direction*/-1,
				 watsonp,genestrand,jump_late_p,querylength,
#ifdef PMAP
				 queryaaseq_ptr,
#endif
				 queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				 maxpeelback,
#ifndef GSNAP
				 oligoindices_minor,diagpool,cellpool, 
#endif
				 pairpool,dynprogL,dynprogM,dynprogR,
				 last_genomedp5_rev,last_genomedp3_rev);

  FREE(last_genomedp3_rev);
  FREE(last_genomedp5_rev);
  FREE(last_genomedp3_fwd);
  FREE(last_genomedp5_fwd);


  debug(printf("Forward:\n"));
  debug(Pair_dump_list(pairs_fwd,true));
  debug(printf("\n"));

  debug(printf("Reverse:\n"));
  debug(Pair_dump_list(pairs_rev,true));
  debug(printf("\n"));

  debug11(printf("Forward:\n"));
  debug11(Pair_dump_list(pairs_fwd,true));
  debug11(printf("\n"));

  debug11(printf("Reverse:\n"));
  debug11(Pair_dump_list(pairs_rev,true));
  debug11(printf("\n"));

  debug(printf("Intronscores: %f,%f fwd, %f,%f rev\n",
	       avg_donor_score_fwd,avg_acceptor_score_fwd,avg_donor_score_rev,avg_acceptor_score_rev));
  if (pairs_rev == NULL) {
    pairs_pretrim = pairs_fwd;
    *cdna_direction = +1;
    *sensedir = SENSE_FORWARD;

  } else if (pairs_fwd == NULL) {
    pairs_pretrim = pairs_rev;
    *cdna_direction = -1;
    *sensedir = SENSE_ANTI;

  } else {
    path_fwd = List_reverse(pairs_fwd);
    debug11(printf("Calling score_introns for path_fwd before path_trim\n"));
    pairs_fwd = score_introns(&max_intron_score_fwd,&avg_donor_score_fwd,&avg_acceptor_score_fwd,
			      &nknown_fwd,&ncanonical_fwd,&nbadintrons_fwd,path_fwd,/*cdna_direction*/+1,watsonp,
			      chrnum,chroffset,chrhigh
#ifdef WASTE
			      ,pairpool
#endif
			      );
    /* alignment_score_fwd = */ score_alignment(&nmatches_fwd,&nmismatches_fwd,&nindels_fwd,
#ifdef COMPLEX_DIRECTION
						&indel_alignment_score_fwd,
#endif
						&nsemicanonical_fwd,&nnoncanonical_fwd,
						pairs_fwd,/*cdna_direction*/+1);
    
    path_rev = List_reverse(pairs_rev);
    debug11(printf("Calling score_introns for path_rev before path_trim\n"));
    pairs_rev = score_introns(&max_intron_score_rev,&avg_donor_score_rev,&avg_acceptor_score_rev,
			      &nknown_rev,&ncanonical_rev,&nbadintrons_rev,path_rev,/*cdna_direction*/-1,watsonp,
			      chrnum,chroffset,chrhigh
#ifdef WASTE
			      ,pairpool
#endif			      
			      );
    /* alignment_score_rev = */ score_alignment(&nmatches_rev,&nmismatches_rev,&nindels_rev,
#ifdef COMPLEX_DIRECTION
						&indel_alignment_score_rev,
#endif
						&nsemicanonical_rev,&nnoncanonical_rev,
						pairs_rev,/*cdna_direction*/-1);

    pairs_pretrim = pick_cdna_direction(&(*cdna_direction),&(*sensedir),pairs_fwd,pairs_rev,
					defect_rate_fwd,defect_rate_rev,
					nknown_fwd,ncanonical_fwd,nsemicanonical_fwd,nnoncanonical_fwd,nbadintrons_fwd,
					nknown_rev,ncanonical_rev,nsemicanonical_rev,nnoncanonical_rev,nbadintrons_rev,
					max_intron_score_fwd,avg_donor_score_fwd,avg_acceptor_score_fwd,
					max_intron_score_rev,avg_donor_score_rev,avg_acceptor_score_rev,
#ifdef COMPLEX_DIRECTION
					nmatches_fwd,nmismatches_fwd,nmatches_rev,nmismatches_rev,nindels_fwd,nindels_rev,
					indel_alignment_score_fwd,indel_alignment_score_rev,
#endif
					sense_filter);
  }


  if (pairs_pretrim == NULL) {
#if 0
    *npairs1 = 0;
    *goodness1 = 0;
    *nmatches_posttrim_1 = 0;
    *ambig_end_length_5_1 = *ambig_end_length_3_1 = 0;
    *ambig_prob_5_1 = *ambig_prob_3_1 = 0.0;
#endif
    return (struct Pair_T *) NULL;
  }

  if (splicingp == false) {
    *sensedir = SENSE_NULL;
  }

#ifdef GSNAP
  if (*cdna_direction == 0) {
    /* If both pairarrays are returned, then first one is fwd and second one is rev */
    debug11(printf("Initial cdna_direction is 0\n"));
    *ambig_end_length_5_1 = fwd_ambig_end_length_5;
    *ambig_end_length_3_1 = fwd_ambig_end_length_3;
    *ambig_splicetype_5_1 = fwd_ambig_splicetype_5;
    *ambig_splicetype_3_1 = fwd_ambig_splicetype_3;
    *ambig_prob_5_1 = fwd_ambig_prob_5;
    *ambig_prob_3_1 = fwd_ambig_prob_3;

    *cdna_direction = +1;

    /* path_trim alters pairs_fwd, so make a copy in case we use it for pairs_pretrim */
    pairs_fwd_copy = Pairpool_copy(pairs_fwd,pairpool);
    *finalpairs1 = path_trim(defect_rate_fwd,&(*ambig_end_length_5_1),&(*ambig_end_length_3_1),
			     &(*ambig_splicetype_5_1),&(*ambig_splicetype_3_1),
			     &(*ambig_prob_5_1),&(*ambig_prob_3_1),
			     pairs_fwd_copy,&(*cdna_direction),watsonp,genestrand,
			     jump_late_p,querylength,
#ifdef GSNAP
			     /*orig_sensedir*/SENSE_FORWARD,
#endif
			     queryseq_ptr,queryuc_ptr,
			     chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
			     maxpeelback,pairpool,dynprogL,dynprogR);

    *ambig_end_length_5_2 = rev_ambig_end_length_5;
    *ambig_end_length_3_2 = rev_ambig_end_length_3;
    *ambig_splicetype_5_2 = rev_ambig_splicetype_5;
    *ambig_splicetype_3_2 = rev_ambig_splicetype_3;
    *ambig_prob_5_2 = rev_ambig_prob_5;
    *ambig_prob_3_2 = rev_ambig_prob_3;

    *cdna_direction = -1;

    pairs_rev_copy = Pairpool_copy(pairs_rev,pairpool);
    *finalpairs2 = path_trim(defect_rate_rev,&(*ambig_end_length_5_2),&(*ambig_end_length_3_2),
			     &(*ambig_splicetype_5_2),&(*ambig_splicetype_3_2),
			     &(*ambig_prob_5_2),&(*ambig_prob_3_2),
			     pairs_rev_copy,&(*cdna_direction),watsonp,genestrand,
			     jump_late_p,querylength,
#ifdef GSNAP
			     /*orig_sensedir*/SENSE_ANTI,
#endif
			     queryseq_ptr,queryuc_ptr,
			     chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
			     maxpeelback,pairpool,dynprogL,dynprogR);

    if (*finalpairs1 != NULL && *finalpairs2 != NULL) {
      debug11(printf("Both directions are non-null, so returning both\n"));
      /* Pairarray 1 (cdna_direction +1): */
      *nmatches_posttrim_1 = Pair_nmatches_posttrim(&(*max_match_length_1),*finalpairs1,/*pos5*/*ambig_end_length_5_1,
						    /*pos3*/querylength - (*ambig_end_length_3_1));
      pairarray1 = make_pairarray(&(*npairs1),&(*finalpairs1),/*cdna_direction*/+1,watsonp,
				  pairpool,queryseq_ptr,chroffset,chrhigh,
				  ngap,query_subseq_offset,skiplength);
      *goodness1 = Pair_fracidentity_array(&(*matches1),&(*unknowns1),&(*mismatches1),
					   &(*qopens1),&(*qindels1),&(*topens1),&(*tindels1),
					   &(*ncanonical1),&(*nsemicanonical1),&(*nnoncanonical1),
					   &min_splice_prob_1,pairarray1,*npairs1,/*cdna_direction*/+1);
      *avg_splice_score_1 = avg_donor_score_fwd + avg_acceptor_score_fwd;


      debug0(printf("Result 1 (%d pairs): %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels, splice score %f\n",
		    *npairs1,*matches1,*mismatches1,*qopens1,*qindels1,*topens1,*tindels1,*avg_splice_score_1));
      debug0(Pair_dump_array(pairarray1,*npairs1,/*zerobasedp*/true));

      /* Note avg_donor_score_fwd and so on do not include evaluations
	 of the end splice junctions.  So if cdna_direction == 0,
	 callers should assume that the sensedir is not known */

      if (0 /*&& Pair_identical_p(*finalpairs1,*finalpairs2) == true*/) {
	/* This causes misses in resolve-inside procedures */
	debug0(printf("Result 2 is identical to Result 1, so not returning it\n"));
	*pairarray2 = (struct Pair_T *) NULL;

      } else {
	/* Pairarray 2 (cdna_direction -1): */
	*nmatches_posttrim_2 = Pair_nmatches_posttrim(&(*max_match_length_2),*finalpairs2,/*pos5*/*ambig_end_length_5_2,
						      /*pos3*/querylength - (*ambig_end_length_3_2));
	*pairarray2 = make_pairarray(&(*npairs2),&(*finalpairs2),/*cdna_direction*/-1,watsonp,
				     pairpool,queryseq_ptr,chroffset,chrhigh,
				     ngap,query_subseq_offset,skiplength);
	*goodness2 = Pair_fracidentity_array(&(*matches2),&(*unknowns2),&(*mismatches2),
					     &(*qopens2),&(*qindels2),&(*topens2),&(*tindels2),
					     &(*ncanonical2),&(*nsemicanonical2),&(*nnoncanonical2),
					     &min_splice_prob_2,*pairarray2,*npairs2,/*cdna_direction*/-1);
	*avg_splice_score_2 = avg_donor_score_rev + avg_acceptor_score_rev;
	
	debug0(printf("Result 2 (%d pairs): %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels, splice score %f\n",
		      *npairs2,*matches2,*mismatches2,*qopens2,*qindels2,*topens2,*tindels2,*avg_splice_score_2));
	debug0(Pair_dump_array(*pairarray2,*npairs2,/*zerobasedp*/true));
      }

      *cdna_direction = 0;
      *sensedir = SENSE_NULL;
      return pairarray1;

    } else if (*finalpairs1 != NULL) {
      debug11(printf("Only forward direction is non-null, so retrying...\n"));
      pairs_pretrim = pairs_fwd;
      *cdna_direction = +1;
      /* Continue below */

    } else if (*finalpairs2 != NULL) {
      debug11(printf("Only reverse direction is non-null, so retrying...\n"));
      pairs_pretrim = pairs_rev;
      *cdna_direction = -1;
      /* Continue below */

    } else {
      return (struct Pair_T *) NULL;
    }
  }
#endif

  if (*cdna_direction > 0) {
    debug11(printf("Solving for forward direction\n"));
    *ambig_end_length_5_1 = fwd_ambig_end_length_5;
    *ambig_end_length_3_1 = fwd_ambig_end_length_3;
    *ambig_splicetype_5_1 = fwd_ambig_splicetype_5;
    *ambig_splicetype_3_1 = fwd_ambig_splicetype_3;
    *ambig_prob_5_1 = fwd_ambig_prob_5;
    *ambig_prob_3_1 = fwd_ambig_prob_3;
    *sensedir = SENSE_FORWARD;
    *avg_splice_score_1 = avg_donor_score_fwd + avg_acceptor_score_fwd;
    defect_rate = defect_rate_fwd;

  } else if (*cdna_direction < 0) {
    debug11(printf("Solving for reverse direction\n"));
    *ambig_end_length_5_1 = rev_ambig_end_length_5;
    *ambig_end_length_3_1 = rev_ambig_end_length_3;
    *ambig_splicetype_5_1 = rev_ambig_splicetype_5;
    *ambig_splicetype_3_1 = rev_ambig_splicetype_3;
    *ambig_prob_5_1 = rev_ambig_prob_5;
    *ambig_prob_3_1 = rev_ambig_prob_3;
    *sensedir = SENSE_ANTI;
    *avg_splice_score_1 = avg_donor_score_rev + avg_acceptor_score_rev;
    defect_rate = defect_rate_rev;

  } else {
#ifdef GSNAP
    abort();
#else
    debug11(printf("Solving for unknown (forward) direction\n"));
    *ambig_end_length_5_1 = fwd_ambig_end_length_5;
    *ambig_end_length_3_1 = fwd_ambig_end_length_3;
    *ambig_splicetype_5_1 = fwd_ambig_splicetype_5;
    *ambig_splicetype_3_1 = fwd_ambig_splicetype_3;
    *ambig_prob_5_1 = fwd_ambig_prob_5;
    *ambig_prob_3_1 = fwd_ambig_prob_3;
    *sensedir = SENSE_FORWARD;
    *avg_splice_score_1 = 0.0;
    defect_rate = defect_rate_fwd;
#endif
  }

  /* Okay for path_trim to alter pairs_pretrim */
  *finalpairs1 = path_trim(defect_rate,&(*ambig_end_length_5_1),&(*ambig_end_length_3_1),
			   &(*ambig_splicetype_5_1),&(*ambig_splicetype_3_1),
			   &(*ambig_prob_5_1),&(*ambig_prob_3_1),
			   pairs_pretrim,&(*cdna_direction),watsonp,genestrand,
			   jump_late_p,querylength,
#ifdef GSNAP
			   /*orig_sensedir*/*sensedir,
#endif
			   queryseq_ptr,queryuc_ptr,
			   chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
			   maxpeelback,pairpool,dynprogL,dynprogR);

  *nmatches_posttrim_1 = Pair_nmatches_posttrim(&(*max_match_length_1),*finalpairs1,/*pos5*/*ambig_end_length_5_1,
						/*pos3*/querylength - (*ambig_end_length_3_1));
  pairarray1 = make_pairarray(&(*npairs1),&(*finalpairs1),*cdna_direction,watsonp,
			      pairpool,queryseq_ptr,chroffset,chrhigh,
			      ngap,query_subseq_offset,skiplength);
  *goodness1 = Pair_fracidentity_array(&(*matches1),&(*unknowns1),&(*mismatches1),
				       &(*qopens1),&(*qindels1),&(*topens1),&(*tindels1),
				       &(*ncanonical1),&(*nsemicanonical1),&(*nnoncanonical1),
				       &min_splice_prob_1,pairarray1,*npairs1,*cdna_direction);
  /* *avg_splice_score_1 assigned above */

  
  debug0(printf("Result (%d pairs): %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels, splice score %f\n",
		*npairs1,*matches1,*mismatches1,*qopens1,*qindels1,*topens1,*tindels1,*avg_splice_score_1));
  debug0(Pair_dump_array(pairarray1,*npairs1,/*zerobasedp*/true));
  
#ifdef GSNAP
  *pairarray2 = (struct Pair_T *) NULL;
  *npairs2 = 0;
#endif

  debug11(printf("Final cdna direction is %d\n",*cdna_direction));
  debug11(printf("Final sensedir is %d\n",*sensedir));

  return pairarray1;
}



/* Using alloca for last_genomedp5 and last_genomedp3 can cause stack overflow */
/* Combines Stage3_compute_middle and Stage3_compute_ends */
struct Pair_T *
Stage3_compute_one (int *cdna_direction, int *sensedir, List_T *finalpairs1, int *npairs1, int *goodness1,
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
		    int sense_try, int sense_filter) {
  struct Pair_T *pairarray1;
  List_T p;
  Chrpos_T *last_genomedp5_fwd = NULL, *last_genomedp3_fwd = NULL, *last_genomedp5_rev = NULL, *last_genomedp3_rev = NULL;
  List_T pairs_pretrim, pairs_fwd, pairs_rev, best_pairs, temp_pairs, path_fwd, path_rev, best_path, temp_path;
  List_T copy;
  List_T joined_ends, joined_starts;
  int nknown_fwd, ncanonical_fwd, nsemicanonical_fwd, nnoncanonical_fwd,
    nknown_rev, ncanonical_rev, nsemicanonical_rev, nnoncanonical_rev;
  int nbadintrons_fwd, nbadintrons_rev;
  double min_splice_prob_1;
  double max_intron_score_fwd = 0.0, max_intron_score_rev = 0.0,
    avg_donor_score_fwd = 0.0, avg_acceptor_score_fwd = 0.0,
    avg_donor_score_rev = 0.0, avg_acceptor_score_rev = 0.0;
  double defect_rate_fwd, defect_rate_rev, defect_rate_temp, defect_rate;
  int nmatches_fwd, nmismatches_fwd, nmatches_rev, nmismatches_rev, nindels_fwd, nindels_rev;
  int fwd_ambig_end_length_5 = 0, fwd_ambig_end_length_3 = 0, rev_ambig_end_length_5 = 0, rev_ambig_end_length_3 = 0, temp_ambig_end_length;
  Splicetype_T fwd_ambig_splicetype_5, fwd_ambig_splicetype_3, rev_ambig_splicetype_5, rev_ambig_splicetype_3, temp_ambig_splicetype;
  double fwd_ambig_prob_5, fwd_ambig_prob_3, rev_ambig_prob_5, rev_ambig_prob_3, temp_ambig_prob;
#ifdef GSNAP
  List_T pairs_fwd_copy, pairs_rev_copy;
  double min_splice_prob_2;
#endif
  

#ifdef COMPLEX_DIRECTION
  int indel_alignment_score_fwd, indel_alignment_score_rev;
#endif

  /* stage2pairs = Stage2_middle(stage2); */
#if defined(DEBUG0) || defined(DEBUG11)
  if (watsonp == true) {
    printf("Stage 3: *** Starting stage 3 at chrnum #%d, chrstart %u, chrend %u, query_subseq_offset %d, sense try %d)\n",
	   chrnum,((Pair_T) stage2pairs->first)->genomepos,
	   ((Pair_T) List_last_value(stage2pairs))->genomepos,query_subseq_offset,sense_try);
  } else {
    printf("Stage 3: *** Starting stage 3 at chrnum #%d, chrstart %u, chrend %u, query_subseq_offset %d, sense try %d)\n",
	   chrnum,(chrhigh - chroffset) - ((Pair_T) stage2pairs->first)->genomepos,
	   (chrhigh - chroffset) - ((Pair_T) List_last_value(stage2pairs))->genomepos,query_subseq_offset,sense_try);
  }
#endif

#ifdef DEBUG
  if (watsonp == true) {
    printf("Stage 3: *** Stage3_compute_one: Starting stage 3 at chrnum #%d, chrstart %u, chrend %u, query_subseq_offset %d, sense try %d)\n",
	   chrnum,((Pair_T) stage2pairs->first)->genomepos,
	   ((Pair_T) List_last_value(stage2pairs))->genomepos,query_subseq_offset,sense_try);
  } else {
    printf("Stage 3: *** Stage3_compute_one: Starting stage 3 at chrnum #%d, chrstart %u, chrend %u, query_subseq_offset %d, sense try %d)\n",
	   chrnum,(chrhigh - chroffset) - ((Pair_T) stage2pairs->first)->genomepos,
	   (chrhigh - chroffset) - ((Pair_T) List_last_value(stage2pairs))->genomepos,query_subseq_offset,sense_try);
  }

  printf("stage2pairs\n");
  Pair_dump_list(stage2pairs,true);

  printf("stage2starts\n");
  for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
    Pair_dump_list((List_T) List_head(p),true);
  }

  printf("stage2ends\n");
  for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
    Pair_dump_list((List_T) List_head(p),true);
  }
#endif


#ifdef PMAP
  pairs_fwd = stage2pairs;
  pairs_rev = (List_T) NULL;
  /* do_final_p = true; */
#else
  if (splicingp == false) {
    pairs_fwd = stage2pairs;
    pairs_rev = (List_T) NULL;
  } else if (sense_try > 0) {
    pairs_fwd = stage2pairs;
    pairs_rev = (List_T) NULL;
  } else if (sense_try < 0) {
    pairs_fwd = (List_T) NULL;
    pairs_rev = stage2pairs;
  } else {
    /* sense_try == 0: Should try both even if no introns (cf, AA011563) */
    pairs_fwd = stage2pairs;
    pairs_rev = Pairpool_copy(stage2pairs,pairpool);
  }
#endif


  /* 1.  Middle */
  if (pairs_fwd == NULL) {
    path_fwd = (List_T) NULL;
  } else {
    last_genomedp5_fwd = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
    last_genomedp3_fwd = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
    debug(printf("*** Solve path_fwd\n"));
    path_fwd = path_compute_dir(&defect_rate_fwd,pairs_fwd,/*cdna_direction*/+1,
				watsonp,genestrand,jump_late_p,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				maxpeelback,
#ifndef GSNAP
				oligoindices_minor,diagpool,cellpool, 
#endif
				pairpool,dynprogL,dynprogM,dynprogR,
				last_genomedp5_fwd,last_genomedp3_fwd/*,clean_ends_p:true*/);
    /* FREE(last_genomedp3_fwd); -- Do not free here, but at end */
    /* FREE(last_genomedp5_fwd); -- Do not free here, but at end */
  }

  if (pairs_rev == NULL) {
    path_rev = (List_T) NULL;
  } else {
    last_genomedp5_rev = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
    last_genomedp3_rev = (Chrpos_T *) CALLOC(querylength,sizeof(Chrpos_T));
    debug(printf("*** Solve path_rev\n"));
    path_rev = path_compute_dir(&defect_rate_rev,pairs_rev,/*cdna_direction*/-1,
				watsonp,genestrand,jump_late_p,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				maxpeelback,
#ifndef GSNAP
				oligoindices_minor,diagpool,cellpool, 
#endif
				pairpool,dynprogL,dynprogM,dynprogR,
				last_genomedp5_rev,last_genomedp3_rev/*,clean_ends_p:true*/);
    /* FREE(last_genomedp5_rev); -- Do not free here, but at end */
    /* FREE(last_genomedp3_rev); -- Do not free here, but at end */
  }


#ifdef GSNAP
  if (path_fwd != NULL && path_rev != NULL) {
    /* Pick cdna_direction based on initial alignment to avoid unnecessary computation */
    pairs_fwd = assign_gap_types(path_fwd,/*cdna_direction*/+1,watsonp,queryseq_ptr,
				 chrnum,chroffset,chrhigh,pairpool);
    path_fwd = List_reverse(pairs_fwd);
    debug11(printf("Calling score_introns for path_fwd after path_compute_dir\n"));
    pairs_fwd = score_introns(&max_intron_score_fwd,&avg_donor_score_fwd,&avg_acceptor_score_fwd,
			      &nknown_fwd,&ncanonical_fwd,&nbadintrons_fwd,path_fwd,/*cdna_direction*/+1,watsonp,
			      chrnum,chroffset,chrhigh
#ifdef WASTE
			      ,pairpool
#endif
			      );

    pairs_rev = assign_gap_types(path_rev,/*cdna_direction*/-1,watsonp,queryseq_ptr,
				 chrnum,chroffset,chrhigh,pairpool);
    path_rev = List_reverse(pairs_rev);
    debug11(printf("Calling score_introns for path_rev after path_compute_dir\n"));
    pairs_rev = score_introns(&max_intron_score_rev,&avg_donor_score_rev,&avg_acceptor_score_rev,
			      &nknown_rev,&ncanonical_rev,&nbadintrons_rev,path_rev,/*cdna_direction*/-1,watsonp,
			      chrnum,chroffset,chrhigh
#ifdef WASTE
			      ,pairpool
#endif			      
			      );
    
    if ((*cdna_direction = initial_cdna_direction(pairs_fwd,pairs_rev,
						  avg_donor_score_fwd,avg_acceptor_score_fwd,
						  avg_donor_score_rev,avg_acceptor_score_rev)) > 0) {
      debug(printf("Initial cdna direction is %d\n",*cdna_direction));
      path_fwd = List_reverse(pairs_fwd);
      path_rev = (List_T) NULL;

    } else if (*cdna_direction < 0) {
      debug(printf("Initial cdna direction is %d\n",*cdna_direction));
      path_fwd = (List_T) NULL;
      path_rev = List_reverse(pairs_rev);

    } else {
      debug(printf("Initial cdna direction is %d\n",*cdna_direction));
      path_fwd = List_reverse(pairs_fwd);
      path_rev = List_reverse(pairs_rev);
    }
  }
#endif


  /* 2.  3' and 5' ends (possibly multiple) */
  debug(printf("Stage2 has %d starts and %d ends\n",List_length(all_stage2_starts),List_length(all_stage2_ends)));
  if (path_fwd == NULL) {
    pairs_fwd = (List_T) NULL;
  } else {
    /* 3' end */
    if (all_stage2_ends == NULL) {
      best_path = path_compute_end3(&fwd_ambig_end_length_3,&fwd_ambig_splicetype_3,&fwd_ambig_prob_3,
				    defect_rate_fwd,path_fwd,/*cdna_direction*/+1,watsonp,genestrand,
				    jump_late_p,querylength,
				    queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				    knownsplice_limit_low,knownsplice_limit_high,
				    maxpeelback,pairpool,dynprogL);
    } else {
      best_path = Pairpool_remove_gapholders(path_fwd); /* Pairpool_join cannot handle gapholders */
      joined_ends = (List_T) NULL;
      for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
#ifdef PMAP      
        copy = Pairpool_join_end3(/*path*/path_fwd,/*end3_pairs*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
#else
	if (path_rev == NULL) {
	  /* Won't need ends anymore */
	  copy = Pairpool_join_end3(/*path*/path_fwd,/*end3_pairs*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
	} else {
	  copy = Pairpool_join_end3(/*path*/path_fwd,/*end3_pairs*/(List_T) List_head(p),pairpool,/*copy_end_p*/true);
	}
#endif
	joined_ends = List_push(joined_ends,(void *) copy);
      }

      for (p = joined_ends; p != NULL; p = List_next(p)) {
        copy = (List_T) List_head(p);
	debug(printf("*** Solve path_fwd joined end\n"));
	path_fwd = path_compute_dir(&defect_rate_temp,/*pairs*/List_reverse(copy),/*cdna_direction*/+1,
				    watsonp,genestrand,jump_late_p,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				    maxpeelback,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool, 
#endif
				    pairpool,dynprogL,dynprogM,dynprogR,
				    last_genomedp5_fwd,last_genomedp3_fwd/*,clean_ends_p:false*/);

	temp_path = path_compute_end3(&temp_ambig_end_length,&temp_ambig_splicetype,&temp_ambig_prob,
				      defect_rate_temp,path_fwd,/*cdna_direction*/+1,watsonp,genestrand,
				      jump_late_p,querylength,
				      queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				      knownsplice_limit_low,knownsplice_limit_high,
				      maxpeelback,pairpool,dynprogL);

	if (temp_path != NULL && end_compare(best_path,temp_path,/*cdna_direction*/+1,watsonp,
					     chrnum,chroffset,chrhigh,/*pairsp*/false) > 0) {
	  best_path = temp_path;
	  fwd_ambig_end_length_3 = temp_ambig_end_length;
	  fwd_ambig_splicetype_3 = temp_ambig_splicetype;
	  fwd_ambig_prob_3 = temp_ambig_prob;
	  defect_rate_fwd = defect_rate_temp;
	  debug21(printf("New best path:\n"));
	  debug21(Pair_dump_list(best_path,true));
	}
      }

      List_free(&joined_ends);
    }

    /* 5' end */
    pairs_fwd = List_reverse(best_path);
    if (all_stage2_starts == NULL) {
      best_pairs = path_compute_end5(&fwd_ambig_end_length_5,&fwd_ambig_splicetype_5,&fwd_ambig_prob_5,
				     defect_rate_fwd,pairs_fwd,/*cdna_direction*/+1,
				     watsonp,genestrand,jump_late_p,
				     queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				     knownsplice_limit_low,knownsplice_limit_high,
				     maxpeelback,pairpool,dynprogR);
    } else {
      best_pairs = Pairpool_remove_gapholders(pairs_fwd); /* Pairpool_join cannot handle gapholders */
      joined_starts = (List_T) NULL;
      for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
#ifdef PMAP
        copy = Pairpool_join_end5(/*pairs*/pairs_fwd,/*end5_path*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
#else
	if (path_rev == NULL) {
	  /* Won't need ends anymore */
	  copy = Pairpool_join_end5(/*pairs*/pairs_fwd,/*end5_path*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
	} else {
	  copy = Pairpool_join_end5(/*pairs*/pairs_fwd,/*end5_path*/(List_T) List_head(p),pairpool,/*copy_end_p*/true);
	}
#endif
	joined_starts = List_push(joined_starts,(void *) copy);
      }

      for (p = joined_starts; p != NULL; p = List_next(p)) {
	copy = (List_T) List_head(p);
	debug(printf("*** Solve path_fwd joined start\n"));
	path_fwd = path_compute_dir(&defect_rate_temp,/*pairs*/copy,/*cdna_direction*/+1,
				    watsonp,genestrand,jump_late_p,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				    maxpeelback,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool, 
#endif
				    pairpool,dynprogL,dynprogM,dynprogR,
				    last_genomedp5_fwd,last_genomedp3_fwd/*,clean_ends_p:false*/);

	temp_pairs = path_compute_end5(&temp_ambig_end_length,&temp_ambig_splicetype,&temp_ambig_prob,
				       defect_rate_temp,/*pairs*/List_reverse(path_fwd),
				       /*cdna_direction*/+1,watsonp,genestrand,jump_late_p,
				       queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				       knownsplice_limit_low,knownsplice_limit_high,
				       maxpeelback,pairpool,dynprogR);

	if (temp_pairs != NULL && end_compare(best_pairs,temp_pairs,/*cdna_direction*/+1,watsonp,
					      chrnum,chroffset,chrhigh,/*pairsp*/true) > 0) {
	  best_pairs = temp_pairs;
	  fwd_ambig_end_length_5 = temp_ambig_end_length;
	  fwd_ambig_splicetype_5 = temp_ambig_splicetype;
	  fwd_ambig_prob_5 = temp_ambig_prob;
	  defect_rate_fwd = defect_rate_temp;
	  debug21(printf("New best pairs:\n"));
	  debug21(Pair_dump_list(best_pairs,true));
	}
      }

      List_free(&joined_starts);
    }

    pairs_fwd = best_pairs;
  }


#ifndef PMAP
  if (path_rev == NULL) {
    pairs_rev = (List_T) NULL;
  } else {
    /* 3' end */
    if (all_stage2_ends == NULL) {
      best_path = path_compute_end3(&rev_ambig_end_length_3,&rev_ambig_splicetype_3,&rev_ambig_prob_3,
				    defect_rate_rev,path_rev,/*cdna_direction*/-1,watsonp,genestrand,
				    jump_late_p,querylength,
				    queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				    knownsplice_limit_low,knownsplice_limit_high,
				    maxpeelback,pairpool,dynprogL);

    } else {
      best_path = Pairpool_remove_gapholders(path_rev); /* Pairpool_join cannot handle gapholders */
      joined_ends = (List_T) NULL;
      for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
	copy = Pairpool_join_end3(/*path*/path_rev,/*end3_pairs*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
	joined_ends = List_push(joined_ends,(void *) copy);
      }

      for (p = joined_ends; p != NULL; p = List_next(p)) {
        copy = (List_T) List_head(p);
	debug(printf("*** Solve path_rev joined end\n"));
	path_rev = path_compute_dir(&defect_rate_temp,/*pairs*/List_reverse(copy),/*cdna_direction*/-1,
				    watsonp,genestrand,jump_late_p,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				    maxpeelback,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool, 
#endif
				    pairpool,dynprogL,dynprogM,dynprogR,
				    last_genomedp5_rev,last_genomedp3_rev/*,clean_ends_p:false*/);
	
	temp_path = path_compute_end3(&temp_ambig_end_length,&temp_ambig_splicetype,&temp_ambig_prob,
				      defect_rate_temp,path_rev,/*cdna_direction*/-1,watsonp,genestrand,
				      jump_late_p,querylength,
				      queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				      knownsplice_limit_low,knownsplice_limit_high,
				      maxpeelback,pairpool,dynprogL);
	
	if (temp_path != NULL && end_compare(best_path,temp_path,/*cdna_direction*/-1,watsonp,
					     chrnum,chroffset,chrhigh,/*pairsp*/false) > 0) {
	  best_path = temp_path;
	  rev_ambig_end_length_3 = temp_ambig_end_length;
	  rev_ambig_splicetype_3 = temp_ambig_splicetype;
	  rev_ambig_prob_3 = temp_ambig_prob;
	  defect_rate_rev = defect_rate_temp;
	  debug21(printf("New best path:\n"));
	  debug21(Pair_dump_list(best_path,true));
	}
      }

      List_free(&joined_ends);
    }

    /* 5' end */
    pairs_rev = List_reverse(best_path);
    if (all_stage2_starts == NULL) {
      best_pairs = path_compute_end5(&rev_ambig_end_length_5,&rev_ambig_splicetype_5,&rev_ambig_prob_5,
				     defect_rate_rev,pairs_rev,/*cdna_direction*/-1,
				     watsonp,genestrand,jump_late_p,
				     queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				     knownsplice_limit_low,knownsplice_limit_high,
				     maxpeelback,pairpool,dynprogR);

    } else {
      best_pairs = Pairpool_remove_gapholders(pairs_rev); /* Pairpool_join cannot handle gapholders */
      joined_starts = (List_T) NULL;
      for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
	copy = Pairpool_join_end5(/*pairs*/pairs_rev,/*end5_path*/(List_T) List_head(p),pairpool,/*copy_end_p*/false);
	joined_starts = List_push(joined_starts,(void *) copy);
      }

      for (p = joined_starts; p != NULL; p = List_next(p)) {
	copy = (List_T) List_head(p);
	debug(printf("*** Solve path_rev joined start\n"));
	path_rev = path_compute_dir(&defect_rate_temp,/*pairs*/copy,/*cdna_direction*/-1,
				    watsonp,genestrand,jump_late_p,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				    maxpeelback,
#ifndef GSNAP
				    oligoindices_minor,diagpool,cellpool, 
#endif
				    pairpool,dynprogL,dynprogM,dynprogR,
				    last_genomedp5_rev,last_genomedp3_rev/*,clean_ends_p:false*/);
      
	temp_pairs = path_compute_end5(&temp_ambig_end_length,&temp_ambig_splicetype,&temp_ambig_prob,
				       defect_rate_temp,/*pairs*/List_reverse(path_rev),
				       /*cdna_direction*/-1,watsonp,genestrand,jump_late_p,
				       queryseq_ptr,queryuc_ptr,querylength,chrnum,chroffset,chrhigh,
				       knownsplice_limit_low,knownsplice_limit_high,
				       maxpeelback,pairpool,dynprogR);
	if (temp_pairs != NULL && end_compare(best_pairs,temp_pairs,/*cdna_direction*/-1,watsonp,
					      chrnum,chroffset,chrhigh,/*pairsp*/true) > 0) {
	  best_pairs = temp_pairs;
	  rev_ambig_end_length_5 = temp_ambig_end_length;
	  rev_ambig_splicetype_5 = temp_ambig_splicetype;
	  rev_ambig_prob_5 = temp_ambig_prob;
	  defect_rate_rev = defect_rate_temp;
	  debug21(printf("New best pairs:\n"));
	  debug21(Pair_dump_list(best_pairs,true));
	}
      }

      List_free(&joined_starts);
    }

    pairs_rev = best_pairs;
  }
#endif

  pairs_fwd = path_compute_final(defect_rate_fwd,pairs_fwd,/*cdna_direction*/+1,
				 watsonp,genestrand,jump_late_p,querylength,
#ifdef PMAP
				 queryaaseq_ptr,
#endif
				 queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				 maxpeelback,
#ifndef GSNAP
				 oligoindices_minor,diagpool,cellpool, 
#endif
				 pairpool,dynprogL,dynprogM,dynprogR,
				 last_genomedp5_fwd,last_genomedp3_fwd);

  pairs_rev = path_compute_final(defect_rate_rev,pairs_rev,/*cdna_direction*/-1,
				 watsonp,genestrand,jump_late_p,querylength,
#ifdef PMAP
				 queryaaseq_ptr,
#endif
				 queryseq_ptr,queryuc_ptr,chrnum,chroffset,chrhigh,
				 maxpeelback,
#ifndef GSNAP
				 oligoindices_minor,diagpool,cellpool, 
#endif
				 pairpool,dynprogL,dynprogM,dynprogR,
				 last_genomedp5_rev,last_genomedp3_rev);

  FREE(last_genomedp3_rev);
  FREE(last_genomedp5_rev);
  FREE(last_genomedp3_fwd);
  FREE(last_genomedp5_fwd);


  debug(printf("Forward:\n"));
  debug(Pair_dump_list(pairs_fwd,true));
  debug(printf("\n"));

  debug(printf("Reverse:\n"));
  debug(Pair_dump_list(pairs_rev,true));
  debug(printf("\n"));

  debug11(printf("Forward:\n"));
  debug11(Pair_dump_list(pairs_fwd,true));
  debug11(printf("\n"));

  debug11(printf("Reverse:\n"));
  debug11(Pair_dump_list(pairs_rev,true));
  debug11(printf("\n"));

  debug(printf("Intronscores: %f,%f fwd, %f,%f rev\n",
	       avg_donor_score_fwd,avg_acceptor_score_fwd,avg_donor_score_rev,avg_acceptor_score_rev));
  if (pairs_rev == NULL) {
    pairs_pretrim = pairs_fwd;
    *cdna_direction = +1;
    *sensedir = SENSE_FORWARD;

  } else if (pairs_fwd == NULL) {
    pairs_pretrim = pairs_rev;
    *cdna_direction = -1;
    *sensedir = SENSE_ANTI;

  } else {
    path_fwd = List_reverse(pairs_fwd);
    debug11(printf("Calling score_introns for path_fwd before path_trim\n"));
    pairs_fwd = score_introns(&max_intron_score_fwd,&avg_donor_score_fwd,&avg_acceptor_score_fwd,
			      &nknown_fwd,&ncanonical_fwd,&nbadintrons_fwd,path_fwd,/*cdna_direction*/+1,watsonp,
			      chrnum,chroffset,chrhigh
#ifdef WASTE
			      ,pairpool
#endif
			      );
    /* alignment_score_fwd = */ score_alignment(&nmatches_fwd,&nmismatches_fwd,&nindels_fwd,
#ifdef COMPLEX_DIRECTION
						&indel_alignment_score_fwd,
#endif
						&nsemicanonical_fwd,&nnoncanonical_fwd,
						pairs_fwd,/*cdna_direction*/+1);
    
    path_rev = List_reverse(pairs_rev);
    debug11(printf("Calling score_introns for path_rev before path_trim\n"));
    pairs_rev = score_introns(&max_intron_score_rev,&avg_donor_score_rev,&avg_acceptor_score_rev,
			      &nknown_rev,&ncanonical_rev,&nbadintrons_rev,path_rev,/*cdna_direction*/-1,watsonp,
			      chrnum,chroffset,chrhigh
#ifdef WASTE
			      ,pairpool
#endif			      
			      );
    /* alignment_score_rev = */ score_alignment(&nmatches_rev,&nmismatches_rev,&nindels_rev,
#ifdef COMPLEX_DIRECTION
						&indel_alignment_score_rev,
#endif
						&nsemicanonical_rev,&nnoncanonical_rev,
						pairs_rev,/*cdna_direction*/-1);

    pairs_pretrim = pick_cdna_direction(&(*cdna_direction),&(*sensedir),pairs_fwd,pairs_rev,
					defect_rate_fwd,defect_rate_rev,
					nknown_fwd,ncanonical_fwd,nsemicanonical_fwd,nnoncanonical_fwd,nbadintrons_fwd,
					nknown_rev,ncanonical_rev,nsemicanonical_rev,nnoncanonical_rev,nbadintrons_rev,
					max_intron_score_fwd,avg_donor_score_fwd,avg_acceptor_score_fwd,
					max_intron_score_rev,avg_donor_score_rev,avg_acceptor_score_rev,
#ifdef COMPLEX_DIRECTION
					nmatches_fwd,nmismatches_fwd,nmatches_rev,nmismatches_rev,nindels_fwd,nindels_rev,
					indel_alignment_score_fwd,indel_alignment_score_rev,
#endif
					sense_filter);
  }


  if (pairs_pretrim == NULL) {
#if 0
    *npairs1 = 0;
    *goodness1 = 0;
    *nmatches_posttrim_1 = 0;
    *ambig_end_length_5_1 = *ambig_end_length_3_1 = 0;
    *ambig_prob_5_1 = *ambig_prob_3_1 = 0.0;
#endif
    return (struct Pair_T *) NULL;
  }

  if (splicingp == false) {
    *sensedir = SENSE_NULL;
  }

#ifdef GSNAP
  if (*cdna_direction == 0) {
    /* If both pairarrays are returned, then first one is fwd and second one is rev */
    debug11(printf("Initial cdna_direction is 0\n"));
    *ambig_end_length_5_1 = fwd_ambig_end_length_5;
    *ambig_end_length_3_1 = fwd_ambig_end_length_3;
    *ambig_splicetype_5_1 = fwd_ambig_splicetype_5;
    *ambig_splicetype_3_1 = fwd_ambig_splicetype_3;
    *ambig_prob_5_1 = fwd_ambig_prob_5;
    *ambig_prob_3_1 = fwd_ambig_prob_3;

    *cdna_direction = +1;

    /* path_trim alters pairs_fwd, so make a copy in case we use it for pairs_pretrim */
    pairs_fwd_copy = Pairpool_copy(pairs_fwd,pairpool);
    *finalpairs1 = path_trim(defect_rate_fwd,&(*ambig_end_length_5_1),&(*ambig_end_length_3_1),
			     &(*ambig_splicetype_5_1),&(*ambig_splicetype_3_1),
			     &(*ambig_prob_5_1),&(*ambig_prob_3_1),
			     pairs_fwd_copy,&(*cdna_direction),watsonp,genestrand,
			     jump_late_p,querylength,
#ifdef GSNAP
			     /*orig_sensedir*/SENSE_FORWARD,
#endif
			     queryseq_ptr,queryuc_ptr,
			     chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
			     maxpeelback,pairpool,dynprogL,dynprogR);

    *ambig_end_length_5_2 = rev_ambig_end_length_5;
    *ambig_end_length_3_2 = rev_ambig_end_length_3;
    *ambig_splicetype_5_2 = rev_ambig_splicetype_5;
    *ambig_splicetype_3_2 = rev_ambig_splicetype_3;
    *ambig_prob_5_2 = rev_ambig_prob_5;
    *ambig_prob_3_2 = rev_ambig_prob_3;

    *cdna_direction = -1;

    pairs_rev_copy = Pairpool_copy(pairs_rev,pairpool);
    *finalpairs2 = path_trim(defect_rate_rev,&(*ambig_end_length_5_2),&(*ambig_end_length_3_2),
			     &(*ambig_splicetype_5_2),&(*ambig_splicetype_3_2),
			     &(*ambig_prob_5_2),&(*ambig_prob_3_2),
			     pairs_rev_copy,&(*cdna_direction),watsonp,genestrand,
			     jump_late_p,querylength,
#ifdef GSNAP
			     /*orig_sensedir*/SENSE_ANTI,
#endif
			     queryseq_ptr,queryuc_ptr,
			     chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
			     maxpeelback,pairpool,dynprogL,dynprogR);

    if (*finalpairs1 != NULL && *finalpairs2 != NULL) {
      debug11(printf("Both directions are non-null, so returning both\n"));
      /* Pairarray 1 (cdna_direction +1): */
      *nmatches_posttrim_1 = Pair_nmatches_posttrim(&(*max_match_length_1),*finalpairs1,/*pos5*/*ambig_end_length_5_1,
						    /*pos3*/querylength - (*ambig_end_length_3_1));
      pairarray1 = make_pairarray(&(*npairs1),&(*finalpairs1),/*cdna_direction*/+1,watsonp,
				  pairpool,queryseq_ptr,chroffset,chrhigh,
				  ngap,query_subseq_offset,skiplength);
      *goodness1 = Pair_fracidentity_array(&(*matches1),&(*unknowns1),&(*mismatches1),
					   &(*qopens1),&(*qindels1),&(*topens1),&(*tindels1),
					   &(*ncanonical1),&(*nsemicanonical1),&(*nnoncanonical1),
					   &min_splice_prob_1,pairarray1,*npairs1,/*cdna_direction*/+1);
      *avg_splice_score_1 = avg_donor_score_fwd + avg_acceptor_score_fwd;


      debug0(printf("Result 1 (%d pairs): %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels, splice score %f\n",
		    *npairs1,*matches1,*mismatches1,*qopens1,*qindels1,*topens1,*tindels1,*avg_splice_score_1));
      debug0(Pair_dump_array(pairarray1,*npairs1,/*zerobasedp*/true));

      /* Note avg_donor_score_fwd and so on do not include evaluations
	 of the end splice junctions.  So if cdna_direction == 0,
	 callers should assume that the sensedir is not known */

      if (0 /*&& Pair_identical_p(*finalpairs1,*finalpairs2) == true*/) {
	/* This causes misses in resolve-inside procedures */
	debug0(printf("Result 2 is identical to Result 1, so not returning it\n"));
	*pairarray2 = (struct Pair_T *) NULL;

      } else {
	/* Pairarray 2 (cdna_direction -1): */
	*nmatches_posttrim_2 = Pair_nmatches_posttrim(&(*max_match_length_2),*finalpairs2,/*pos5*/*ambig_end_length_5_2,
						      /*pos3*/querylength - (*ambig_end_length_3_2));
	*pairarray2 = make_pairarray(&(*npairs2),&(*finalpairs2),/*cdna_direction*/-1,watsonp,
				     pairpool,queryseq_ptr,chroffset,chrhigh,
				     ngap,query_subseq_offset,skiplength);
	*goodness2 = Pair_fracidentity_array(&(*matches2),&(*unknowns2),&(*mismatches2),
					     &(*qopens2),&(*qindels2),&(*topens2),&(*tindels2),
					     &(*ncanonical2),&(*nsemicanonical2),&(*nnoncanonical2),
					     &min_splice_prob_2,*pairarray2,*npairs2,/*cdna_direction*/-1);
	*avg_splice_score_2 = avg_donor_score_rev + avg_acceptor_score_rev;
	
	debug0(printf("Result 2 (%d pairs): %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels, splice score %f\n",
		      *npairs2,*matches2,*mismatches2,*qopens2,*qindels2,*topens2,*tindels2,*avg_splice_score_2));
	debug0(Pair_dump_array(*pairarray2,*npairs2,/*zerobasedp*/true));
      }

      *cdna_direction = 0;
      *sensedir = SENSE_NULL;
      return pairarray1;

    } else if (*finalpairs1 != NULL) {
      debug11(printf("Only forward direction is non-null, so retrying...\n"));
      pairs_pretrim = pairs_fwd;
      *cdna_direction = +1;
      /* Continue below */

    } else if (*finalpairs2 != NULL) {
      debug11(printf("Only reverse direction is non-null, so retrying...\n"));
      pairs_pretrim = pairs_rev;
      *cdna_direction = -1;
      /* Continue below */

    } else {
      return (struct Pair_T *) NULL;
    }
  }
#endif

  if (*cdna_direction > 0) {
    debug11(printf("Solving for forward direction\n"));
    *ambig_end_length_5_1 = fwd_ambig_end_length_5;
    *ambig_end_length_3_1 = fwd_ambig_end_length_3;
    *ambig_splicetype_5_1 = fwd_ambig_splicetype_5;
    *ambig_splicetype_3_1 = fwd_ambig_splicetype_3;
    *ambig_prob_5_1 = fwd_ambig_prob_5;
    *ambig_prob_3_1 = fwd_ambig_prob_3;
    *sensedir = SENSE_FORWARD;
    *avg_splice_score_1 = avg_donor_score_fwd + avg_acceptor_score_fwd;
    defect_rate = defect_rate_fwd;

  } else if (*cdna_direction < 0) {
    debug11(printf("Solving for reverse direction\n"));
    *ambig_end_length_5_1 = rev_ambig_end_length_5;
    *ambig_end_length_3_1 = rev_ambig_end_length_3;
    *ambig_splicetype_5_1 = rev_ambig_splicetype_5;
    *ambig_splicetype_3_1 = rev_ambig_splicetype_3;
    *ambig_prob_5_1 = rev_ambig_prob_5;
    *ambig_prob_3_1 = rev_ambig_prob_3;
    *sensedir = SENSE_ANTI;
    *avg_splice_score_1 = avg_donor_score_rev + avg_acceptor_score_rev;
    defect_rate = defect_rate_rev;

  } else {
#ifdef GSNAP
    abort();
#else
    debug11(printf("Solving for unknown (forward) direction\n"));
    *ambig_end_length_5_1 = fwd_ambig_end_length_5;
    *ambig_end_length_3_1 = fwd_ambig_end_length_3;
    *ambig_splicetype_5_1 = fwd_ambig_splicetype_5;
    *ambig_splicetype_3_1 = fwd_ambig_splicetype_3;
    *ambig_prob_5_1 = fwd_ambig_prob_5;
    *ambig_prob_3_1 = fwd_ambig_prob_3;
    *sensedir = SENSE_FORWARD;
    *avg_splice_score_1 = 0.0;
    defect_rate = defect_rate_fwd;
#endif
  }

  /* Okay for path_trim to alter pairs_pretrim */
  *finalpairs1 = path_trim(defect_rate,&(*ambig_end_length_5_1),&(*ambig_end_length_3_1),
			   &(*ambig_splicetype_5_1),&(*ambig_splicetype_3_1),
			   &(*ambig_prob_5_1),&(*ambig_prob_3_1),
			   pairs_pretrim,&(*cdna_direction),watsonp,genestrand,
			   jump_late_p,querylength,
#ifdef GSNAP
			   /*orig_sensedir*/*sensedir,
#endif
			   queryseq_ptr,queryuc_ptr,
			   chroffset,chrhigh,knownsplice_limit_low,knownsplice_limit_high,
			   maxpeelback,pairpool,dynprogL,dynprogR);

  *nmatches_posttrim_1 = Pair_nmatches_posttrim(&(*max_match_length_1),*finalpairs1,/*pos5*/*ambig_end_length_5_1,
						/*pos3*/querylength - (*ambig_end_length_3_1));
  pairarray1 = make_pairarray(&(*npairs1),&(*finalpairs1),*cdna_direction,watsonp,
			      pairpool,queryseq_ptr,chroffset,chrhigh,
			      ngap,query_subseq_offset,skiplength);
  *goodness1 = Pair_fracidentity_array(&(*matches1),&(*unknowns1),&(*mismatches1),
				       &(*qopens1),&(*qindels1),&(*topens1),&(*tindels1),
				       &(*ncanonical1),&(*nsemicanonical1),&(*nnoncanonical1),
				       &min_splice_prob_1,pairarray1,*npairs1,*cdna_direction);
  /* *avg_splice_score_1 assigned above */

  
  debug0(printf("Result (%d pairs): %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels, splice score %f\n",
		*npairs1,*matches1,*mismatches1,*qopens1,*qindels1,*topens1,*tindels1,*avg_splice_score_1));
  debug0(Pair_dump_array(pairarray1,*npairs1,/*zerobasedp*/true));
  
#ifdef GSNAP
  *pairarray2 = (struct Pair_T *) NULL;
  *npairs2 = 0;
#endif

  debug11(printf("Final cdna direction is %d\n",*cdna_direction));
  debug11(printf("Final sensedir is %d\n",*sensedir));

  return pairarray1;
}


/************************************************************************
 *  Merging
 ************************************************************************/

bool
Stage3_mergeable (Stage3_T firstpart, Stage3_T secondpart,
		  int breakpoint, int queryntlength) {
  Pair_T end1, start2;
  bool watsonp, connectablep = false;
  Chrpos_T endchrpos1, startchrpos2;
  int npairs_left, npairs_right, nstart;
  int cdna_direction_1, cdna_direction_2;


  assert(breakpoint > 0);	/* Caller should check for this */
  assert(firstpart->pairs != NULL);
  assert(secondpart->pairs != NULL);

  debug20(printf("Stage3_mergeable called with breakpoint %d, watsonp %d and %d, and cdna_directions %d and %d\n",
		 breakpoint,firstpart->watsonp,secondpart->watsonp,firstpart->cdna_direction,secondpart->cdna_direction));
  debug10(Stage3_print_ends(firstpart));
  debug10(Stage3_print_ends(secondpart));

  if (firstpart->chrnum != secondpart->chrnum) {
    debug20(printf("not mergeable: chrnum %d != chrnum %d\n",firstpart->chrnum,secondpart->chrnum));
    return false;

  } else if (firstpart->watsonp != secondpart->watsonp) {
    debug20(printf("not mergeable: watsonp %d != watsonp %d\n",firstpart->watsonp,secondpart->watsonp));
    return false;

#if 0
  } else if (firstpart->sensedir != secondpart->sensedir &&
	     firstpart->sensedir != SENSE_NULL && secondpart->sensedir != SENSE_NULL) {
    /* Could be mergeable if an intron is trimmed during the merge */
    debug20(printf("not mergeable: sensedir %d != sensedir %d\n",
		   firstpart->sensedir,secondpart->sensedir));
    return false;
#endif

  } else {
    /* Find end pairs. Ignore cdna directions for now. */
    end1 = Pair_end_bound(&cdna_direction_1,firstpart->pairs,breakpoint);
    start2 = Pair_start_bound(&cdna_direction_2,secondpart->pairs,breakpoint+1);

    if ((watsonp = firstpart->watsonp) == true) {
      endchrpos1 = end1->genomepos;
      startchrpos2 = start2->genomepos;

      debug20(printf("? connectable, watson: endchrpos1 %d at querypos %d versus startchrpos2 %d at querypos %d\n",
		     endchrpos1,end1->querypos,startchrpos2,start2->querypos));
      if (endchrpos1 < startchrpos2) {
	/* Deletion */
	/* *genomejump = startchrpos2 - endchrpos1 - 1; */
	debug20(printf("endchrpos1 < startchrpos2, so deletion of length %u\n",startchrpos2 - endchrpos1 - 1));
	if (startchrpos2 < endchrpos1 + maxintronlen) {
	  connectablep = true;
	}

      } else if (startchrpos2 + (end1->querypos - start2->querypos) + 100 >= endchrpos1) {
	/* Insertion */
	debug20(printf("startchrpos2 + (%d - %d) + %d >= endchrpos2, so insertion\n",end1->querypos,start2->querypos,20));
	/* *genomejump = 0; */
	connectablep = true;
      }


    } else {
      /* These are genomicpos, not really chrpos.  If go to chrpos, need to rewrite logic. */
      endchrpos1 = firstpart->chrhigh - end1->genomepos;
      startchrpos2 = secondpart->chrhigh - start2->genomepos;

      debug20(printf("? connectable, crick: startchrpos2 %u at querypos %d versus endchrpos1 %u at querypos %d\n",
		     startchrpos2,start2->querypos,endchrpos1,end1->querypos));
      if (startchrpos2 < endchrpos1) {
	/* Deletion */
	/* *genomejump = endchrpos1 - startchrpos2 - 1; */
	debug20(printf("startchrpos2 < endchrpos1, so deletion of length %u\n",endchrpos1 - startchrpos2 - 1));
	if (endchrpos1 < startchrpos2 + maxintronlen) {
	  connectablep = true;
	}
      } else if (endchrpos1 + (end1->querypos - start2->querypos) + 100 >= startchrpos2) {
	/* Insertion */
	debug20(printf("endchrpos1 + (%d - %d) + %d >= endchrpos1, so insertion\n",end1->querypos,start2->querypos,20));
	/* *genomejump = 0; */
	connectablep = true;
      }
    }

    if (connectablep == false) {
      debug20(printf("result: not mergeable\n\n"));
      return false;
    } else {
      npairs_left = Pairpool_count_bounded(&nstart,firstpart->pairs,0,breakpoint);
      npairs_right = Pairpool_count_bounded(&nstart,secondpart->pairs,breakpoint,queryntlength);
      debug20(printf("Predicted after splicing: npairs_left %d, npairs_right %d\n",npairs_left,npairs_right));
      if (npairs_left < 25 || npairs_right < 25) {
	return false;
      } else {
	/* *queryjump = start2->querypos - end1->querypos - 1; */
	debug20(printf("result: mergeable: queryjump = %d - %d - 1 = %d\n\n",
		       start2->querypos,end1->querypos,start2->querypos - end1->querypos - 1));
		       
	return true;
      }
    }
  }
}


bool
Stage3_merge_chimera (T *new_left, T *new_right, T old_left, T old_right,
		      int minpos1, int maxpos1, int minpos2, int maxpos2,
		      Sequence_T queryseq, char *queryseq_ptr, char *queryuc_ptr, Pairpool_T pairpool, 
		      Dynprog_T dynprogL, Dynprog_T dynprogR, int maxpeelback) {
  List_T left_pairs, right_pairs;
  List_T path, peeled_path, pairs, peeled_pairs;
  bool knownsplicep, chop_exon_p;
  int ambig_end_length_5 = 0, ambig_end_length_3 = 0;	/* Need to be set for build_pairs_end5 and build_path_end3 */
  double ambig_prob_5, ambig_prob_3;
  int dynprogindex_minor = 0;
  Splicetype_T ambig_splicetype;

  Pair_T endpair;
  int querydp5, querydp3, n_peeled_indels;
  Chrpos_T genomedp5, genomedp3;
  bool protectedp;

  debug10(printf("Entered stage3_merge_chimera with minpos1 %d, maxpos1 %d, minpos2 %d, maxpos2 %d\n",
		minpos1,maxpos1,minpos2,maxpos2));

#ifdef DEBUG10
  printf("ORIG LEFT PAIRS:\n");
  Pair_dump_list(old_left->pairs,true);
  printf("ORIG RIGHT PAIRS:\n");
  Pair_dump_list(old_right->pairs,true);
#endif

  left_pairs = Pairpool_copy(old_left->pairs,pairpool);
  left_pairs = Pair_clip_bounded_list_5(left_pairs,minpos1,maxpos1);
  right_pairs = Pairpool_copy(old_right->pairs,pairpool);
  right_pairs = Pair_clip_bounded_list_3(right_pairs,minpos2,maxpos2);

#ifdef DEBUG10
  printf("CLIPPED LEFT PAIRS:\n");
  Pair_dump_list(left_pairs,true);
  printf("CLIPPED RIGHT PAIRS:\n");
  Pair_dump_list(right_pairs,true);
#endif

  if (left_pairs == NULL || right_pairs == NULL) {
    *new_left = *new_right = (T) NULL;
    return false;

  } else {
    *new_left = Stage3_new_from_pairs(left_pairs,old_left->cdna_direction,old_left->watsonp,
				      old_left->genestrand,old_left->sensedir,pairpool,queryseq,
				      /*query_subseq_offset*/Sequence_subseq_offset(queryseq),
				      old_left->chrnum,old_left->chroffset,old_left->chrhigh,old_left->chrlength);
    *new_right = Stage3_new_from_pairs(right_pairs,old_right->cdna_direction,old_right->watsonp,
				       old_right->genestrand,old_right->sensedir,pairpool,queryseq,
				       /*query_subseq_offset*/Sequence_subseq_offset(queryseq),
				       old_right->chrnum,old_right->chroffset,old_right->chrhigh,old_right->chrlength);

    path = List_reverse((*new_left)->pairs);

    /* To avoid indels at chimeric join, need to peelback, clean ends, extend with nogaps, and then re-clip*/
    endpair = (Pair_T) path->first;
    querydp5 = endpair->querypos + 1;
    genomedp5 = endpair->genomepos + 1;
    protectedp = false;
    path = peel_leftward(&n_peeled_indels,&protectedp,&peeled_path,path,&querydp5,&genomedp5,
			 maxpeelback,/*stop_at_indels_p*/true);
    path = clean_path_end3_gap_indels(path);

    /* Have to use forcep == true, because we cannot put back the
       pairs removed from peel_leftward and clean_path_end3_gap_indels */
    path = build_path_end3(&knownsplicep,&ambig_end_length_3,&ambig_splicetype,&ambig_prob_3,
			   &chop_exon_p,&dynprogindex_minor,path,
			   (*new_left)->chroffset,(*new_left)->chrhigh,/*querylength*/maxpos1+1,
			   /*knownsplice_limit_low*/-1,/*knownsplice_limit_high*/0,
			   queryseq_ptr,queryuc_ptr,
			   (*new_left)->cdna_direction,(*new_left)->watsonp,(*new_left)->genestrand,
			   /*jump_late_p*/(*new_left)->watsonp ? false : true,
			   maxpeelback,/*defect_rate*/0.0,pairpool,dynprogL,
			   /*extendp*/true,/*endalign*/QUERYEND_NOGAPS,/*forcep*/true);

    (*new_left)->pairs = List_reverse(path);
    (*new_left)->pairs = Pair_clip_bounded_list_5((*new_left)->pairs,minpos1,maxpos1);
    debug10(printf("CLEANED LEFT PAIRS:\n"));
    debug10(Pair_dump_list((*new_left)->pairs,true));

    /* To avoid indels at chimeric join, need to peelback, clean ends, extend with nogaps, and then re-clip*/
    pairs = (*new_right)->pairs;

    endpair = (Pair_T) pairs->first;
    querydp3 = endpair->querypos - 1;
    genomedp3 = endpair->genomepos - 1;
    protectedp = false;
    pairs = peel_rightward(&n_peeled_indels,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
			   maxpeelback,/*stop_at_indels_p*/true);
    pairs = clean_pairs_end5_gap_indels(pairs);

    /* Have to use forcep == true, because we cannot put back the
       pairs removed from peel_rightward and clean_pairs_end5_gap_indels */
    pairs = build_pairs_end5(&knownsplicep,&ambig_end_length_5,&ambig_splicetype,&ambig_prob_5,
			     &chop_exon_p,&dynprogindex_minor,pairs,
			     (*new_right)->chroffset,(*new_right)->chrhigh,
			     /*knownsplice_limit_low*/-1,/*knownsplice_limit_high*/0,
			     queryseq_ptr,queryuc_ptr,
			     (*new_right)->cdna_direction,(*new_right)->watsonp,(*new_right)->genestrand,
			     /*jump_late_p*/(*new_right)->watsonp ? false : true,
			     maxpeelback,/*defect_rate*/0.0,pairpool,dynprogR,
			     /*extendp*/true,/*endalign*/QUERYEND_NOGAPS,/*forcep*/true);
    (*new_right)->pairs = Pair_clip_bounded_list_3(pairs,minpos2,maxpos2);
    debug10(printf("CLEANED RIGHT PAIRS:\n"));
    debug10(Pair_dump_list((*new_right)->pairs,true));

    if ((*new_left)->pairs == NULL || (*new_right)->pairs == NULL) {
      Stage3_free(&(*new_left));
      Stage3_free(&(*new_right));
      *new_left = *new_right = (T) NULL;
      return false;
    } else {
      make_pairarrays_chimera(*new_left,*new_right,queryseq_ptr,pairpool,/*gaplength*/0,ngap);

      (*new_left)->chimera_right_p = true;
      (*new_right)->chimera_left_p = true;
      return true;
    }
  }
}



void
Stage3_extend_right (T this, int goal, int querylength,
		     char *queryseq_ptr, char *queryuc_ptr,
		     bool max_extend_p, Pairpool_T pairpool,
		     int genestrand, int maxpeelback) {
  List_T path, peeled_path;
  Pair_T leftpair;

  int nconsecutive_mismatches;
  int querypos, querydp5;
  Chrpos_T genomedp5;
  int genomepos;
  char c, c_upper, g, g_alt, comp;
  bool protectedp;
  int n_peeled_indels;

  int ncanonical, nsemicanonical;
  double min_splice_prob;


  debug10(printf("Entered Stage3_extend_right with goal %d\n",goal));
  debug10(printf("LEFT BEFORE FILL\n"));
  debug10(Pair_dump_list(this->pairs,true));
  debug10(printf("END_LEFT BEFORE FILL\n"));


  path = List_reverse(this->pairs);
  leftpair = (Pair_T) path->first;

  debug(printf("\nEXTEND_RIGHT\n"));
  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  /* if (leftpair->cdna == ' ') querydp5--; -- For old dynamic programming */
  /* if (leftpair->genome == ' ') genomedp5--; -- For old dynamic programming */

  protectedp = false;
  path = peel_leftward(&n_peeled_indels,&protectedp,&peeled_path,path,&querydp5,&genomedp5,
		       maxpeelback,/*stop_at_indels_p*/true);
  if (path == NULL) {
    querypos = querydp5 - 1;
    genomepos = genomedp5 - 1;
  } else {
    path = clean_path_end3_gap_indels(path);
    leftpair = (Pair_T) path->first;
    querypos = leftpair->querypos;
    genomepos = leftpair->genomepos;
  }

  if (this->watsonp == true) {
    /* pos = this->chroffset + genomepos; */
    debug10(printf("watsonp on left is true.  pos is %u.  goal querypos is %d\n",genomepos,goal));

    querypos++;
    genomepos++;
    /* pos++; */
    while (querypos < goal /* && pos <= this->chrhigh */) {
      c = queryseq_ptr[querypos];
      c_upper = queryuc_ptr[querypos];
      /* g = Genome_get_char(genome,pos); */
      g = get_genomic_nt(&g_alt,genomepos,this->chroffset,this->chrhigh,/*watsonp*/true);
      if (g != '*') {
	if (Dynprog_consistent_p(c_upper,g,g_alt,genestrand) == true) {
	  comp = MATCH_COMP;
#ifdef PMAP
	} else if (Dynprog_consistent_p(c_upper,g,g_alt) == true) {
	  comp = AMBIGUOUS_COMP;
#endif
	} else {
	  comp = MISMATCH_COMP;
	}
	debug10(printf("At querypos %d and pos %u on left, have %c and %c\n",
		       querypos,genomepos,c,g));
	path = Pairpool_push(path,pairpool,querypos,genomepos,c,comp,g,g_alt,/*dynprogindex*/0);
      }
      querypos++;
      genomepos++;
      /* pos++; */
    }

    if (max_extend_p == true) {
      debug10(printf("\nGoal achieved.  Now looking for consecutive mismatches\n"));
      nconsecutive_mismatches = 0;
      while (querypos < querylength /* && pos <= this->chrhigh */ && nconsecutive_mismatches < 3) {
	c = queryseq_ptr[querypos];
	c_upper = queryuc_ptr[querypos];
	/* g = Genome_get_char(genome,pos); */
	g = get_genomic_nt(&g_alt,genomepos,this->chroffset,this->chrhigh,/*watsonp*/true);
	if (g != '*') {
	  if (Dynprog_consistent_p(c_upper,g,g_alt,genestrand) == true) {
	    comp = MATCH_COMP;
	    nconsecutive_mismatches = 0;
#ifdef PMAP
	  } else if (Dynprog_consistent_p(c_upper,g,g_alt) == true) {
	    comp = AMBIGUOUS_COMP;
	    /* Don't consider as match or mismatch */
#endif
	  } else {
	    comp = MISMATCH_COMP;
	    nconsecutive_mismatches += 1;
	  }
	  debug10(printf("At querypos %d and pos %u on left, have %c and %c\n",
			 querypos,genomepos,c,g));
	  path = Pairpool_push(path,pairpool,querypos,genomepos,c,comp,g,g_alt,/*dynprogindex*/0);
	}
	querypos++;
	genomepos++;
	/* pos++; */
      }
    }

  } else {
    /* pos = this->chrhigh - genomepos; */
    debug10(printf("watsonp on left is false.  pos is %u.  querypos is %d.  want to go up to goal querypos %d\n",
		   genomepos,querypos,goal));

    querypos++;
    genomepos++;
    /* pos--; */
    while (querypos < goal /* && pos != this->chroffset - 1U */) {
      c = queryseq_ptr[querypos];
      c_upper = queryuc_ptr[querypos];
      /* g = complCode[(int) Genome_get_char(genome,pos)]; */
      g = get_genomic_nt(&g_alt,genomepos,this->chroffset,this->chrhigh,/*watsonp*/false);
      if (g != '*') {
	if (Dynprog_consistent_p(c_upper,g,g_alt,genestrand) == true) {
	  comp = MATCH_COMP;
#ifdef PMAP
	} else if (Dynprog_consistent_p(c_upper,g,g_alt) == true) {
	  comp = AMBIGUOUS_COMP;
#endif
	} else {
	  comp = MISMATCH_COMP;
	}
	debug10(printf("At querypos %d and pos %u on left, have %c and %c\n",
		       querypos,genomepos,c,g));
	path = Pairpool_push(path,pairpool,querypos,genomepos,c,comp,g,g_alt,/*dynprogindex*/0);
      }
      querypos++;
      genomepos++;
      /* pos--; */
    }

    if (max_extend_p == true) {
      debug10(printf("\nGoal achieved.  Now looking for consecutive mismatches\n"));
      nconsecutive_mismatches = 0;
      while (querypos < querylength /* && pos != this->chroffset - 1U */ && nconsecutive_mismatches < 3) {
	c = queryseq_ptr[querypos];
	c_upper = queryuc_ptr[querypos];
	/* g = complCode[(int) Genome_get_char(genome,pos)]; */
	g = get_genomic_nt(&g_alt,genomepos,this->chroffset,this->chrhigh,/*watsonp*/false);
	if (g != '*') {
	  if (Dynprog_consistent_p(c_upper,g,g_alt,genestrand) == true) {
	    comp = MATCH_COMP;
	    nconsecutive_mismatches = 0;
#ifdef PMAP
	  } else if (Dynprog_consistent_p(c_upper,g,g_alt) == true) {
	    comp = AMBIGUOUS_COMP;
	    /* Don't count as either match or mismatch */
#endif
	  } else {
	    comp = MISMATCH_COMP;
	    nconsecutive_mismatches += 1;
	  }
	  debug10(printf("At querypos %d and pos %u on left, have %c and %c\n",
			 querypos,genomepos,c,g));
	  path = Pairpool_push(path,pairpool,querypos,genomepos,c,comp,g,g_alt,/*dynprogindex*/0);
	}
	querypos++;
	genomepos++;
	/* pos--; */
      }
    }
  }

  this->pairs = List_reverse(path);

  debug10(printf("LEFT AFTER FILL\n"));
  debug10(Pair_dump_list(this->pairs,true));
  debug10(printf("END_LEFT AFTER FILL\n"));

  Stage3_free_pairarray(&this);
  this->pairarray = make_pairarray(&this->npairs,&this->pairs,this->cdna_direction,
				   this->watsonp,pairpool,queryseq_ptr,
				   this->chroffset,this->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0);
  this->goodness = Pair_fracidentity_array(&this->matches,&this->unknowns,&this->mismatches,
					   &this->qopens,&this->qindels,&this->topens,&this->tindels,
					   &ncanonical,&nsemicanonical,&this->noncanonical,
					   &min_splice_prob,this->pairarray,this->npairs,this->cdna_direction);

  if (this->pairarray == NULL) {
    this->pairarray_freeable_p = false;
  } else {
    this->pairarray_freeable_p = true;
  }

  return;
}


void
Stage3_extend_left (T this, int goal,
		    char *queryseq_ptr, char *queryuc_ptr,
		    bool max_extend_p, Pairpool_T pairpool,
		    int genestrand, int maxpeelback) {
  List_T pairs, peeled_pairs;
  Pair_T rightpair;

  int nconsecutive_mismatches;
  int querypos, querydp3;
  Chrpos_T genomedp3;
  int genomepos;
  char c, c_upper, g, g_alt, comp;
  bool protectedp;
  int n_peeled_indels;

  int ncanonical, nsemicanonical;
  double min_splice_prob;


  debug10(printf("Entered Stage3_extend_left with goal %d\n",goal));
  debug10(printf("RIGHT BEFORE FILL\n"));
  debug10(Pair_dump_list(this->pairs,true));
  debug10(printf("END_RIGHT BEFORE FILL\n"));


  /* Do not call insert_gapholders */
  pairs = this->pairs;
  rightpair = (Pair_T) pairs->first;

  debug(printf("\nEXTEND_LEFT\n"));
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

  protectedp = false;
  pairs = peel_rightward(&n_peeled_indels,&protectedp,&peeled_pairs,pairs,&querydp3,&genomedp3,
			 maxpeelback,/*stop_at_indels_p*/true);
  if (pairs == NULL) {
    querypos = querydp3 + 1;
    genomepos = genomedp3 + 1;
  } else {
    pairs = clean_pairs_end5_gap_indels(pairs);
    rightpair = (Pair_T) pairs->first;
    querypos = rightpair->querypos;
    genomepos = rightpair->genomepos;
  }
  
  if (this->watsonp == true) {
    /* pos = this->chroffset + genomepos; */
    debug10(printf("watsonp on right is true.  pos is %u.  goal querypos is %d\n",genomepos,goal));

    querypos--;
    genomepos--;
    /* pos--; */
    while (querypos >= goal /* && pos != this->chroffset - 1U */) {
      c = queryseq_ptr[querypos];
      c_upper = queryuc_ptr[querypos];
      /* g = Genome_get_char(genome,pos); */
      g = get_genomic_nt(&g_alt,genomepos,this->chroffset,this->chrhigh,/*watsonp*/true);
      if (g != '*') {
	if (Dynprog_consistent_p(c_upper,g,g_alt,genestrand) == true) {
	  comp = MATCH_COMP;
#ifdef PMAP
	} else if (Dynprog_consistent_p(c_upper,g,g_alt) == true) {
	  comp = AMBIGUOUS_COMP;
#endif
	} else {
	  comp = MISMATCH_COMP;
	}
	debug10(printf("At querypos %d and pos %u on right, have %c and %c\n",
		       querypos,genomepos,c,g));
	pairs = Pairpool_push(pairs,pairpool,querypos,genomepos,c,comp,g,g_alt,/*dynprogindex*/0);
      }
      querypos--;
      genomepos--;
      /* pos--; */
    }

    if (max_extend_p == true) {
      debug10(printf("\nGoal achieved.  Now looking for consecutive mismatches\n"));
      nconsecutive_mismatches = 0;
      while (querypos >= 0 /* && pos != this->chroffset - 1U */ && nconsecutive_mismatches < 3) {
	c = queryseq_ptr[querypos];
	c_upper = queryuc_ptr[querypos];
	/* g = Genome_get_char(genome,pos); */
	g = get_genomic_nt(&g_alt,genomepos,this->chroffset,this->chrhigh,/*watsonp*/true);
	if (g != '*') {
	  if (Dynprog_consistent_p(c_upper,g,g_alt,genestrand) == true) {
	    comp = MATCH_COMP;
	    nconsecutive_mismatches = 0;
#ifdef PMAP
	  } else if (Dynprog_consistent_p(c_upper,g,g_alt) == true) {
	    comp = AMBIGUOUS_COMP;
	    /* Don't count as either match or mismatch */
#endif
	  } else {
	    comp = MISMATCH_COMP;
	    nconsecutive_mismatches += 1;
	  }
	  debug10(printf("At querypos %d and pos %u on right, have %c and %c\n",
			 querypos,genomepos,c,g));
	  pairs = Pairpool_push(pairs,pairpool,querypos,genomepos,c,comp,g,g_alt,/*dynprogindex*/0);
	}
	querypos--;
	genomepos--;
	/* pos--; */
      }
    }

  } else {
    /* pos = this->chrhigh - genomepos; */
    debug10(printf("watsonp on right is false.  pos is %u.  goal querypos is %d\n",genomepos,goal));

    querypos--;
    genomepos--;
    /* pos++; */
    while (querypos >= goal /* && pos <= this->chrhigh */) {
      c = queryseq_ptr[querypos];
      c_upper = queryuc_ptr[querypos];
      /* g = complCode[(int) Genome_get_char(genome,pos)]; */
      g = get_genomic_nt(&g_alt,genomepos,this->chroffset,this->chrhigh,/*watsonp*/false);
      if (g != '*') {
	if (Dynprog_consistent_p(c_upper,g,g_alt,genestrand) == true) {
	  comp = MATCH_COMP;
#ifdef PMAP
	} else if (Dynprog_consistent_p(c_upper,g,g_alt) == true) {
	  comp = AMBIGUOUS_COMP;
#endif
	} else {
	  comp = MISMATCH_COMP;
	}
	debug10(printf("At querypos %d and pos %u on right, have %c and %c\n",
		       querypos,genomepos,c,g));
	pairs = Pairpool_push(pairs,pairpool,querypos,genomepos,c,comp,g,g_alt,/*dynprogindex*/0);
      }
      querypos--;
      genomepos--;
      /* pos++; */
    }

    if (max_extend_p == true) {
      debug10(printf("\nGoal achieved.  Now looking for consecutive mismatches\n"));
      nconsecutive_mismatches = 0;
      while (querypos >= 0 /* && pos <= this->chrhigh */ && nconsecutive_mismatches > 3) {
	c = queryseq_ptr[querypos];
	c_upper = queryuc_ptr[querypos];
	/* g = complCode[(int) Genome_get_char(genome,pos)]; */
	g = get_genomic_nt(&g_alt,genomepos,this->chroffset,this->chrhigh,/*watsonp*/false);
	if (g != '*') {
	  if (Dynprog_consistent_p(c_upper,g,g_alt,genestrand) == true) {
	    comp = MATCH_COMP;
	    nconsecutive_mismatches = 0;
#ifdef PMAP
	  } else if (Dynprog_consistent_p(c_upper,g,g_alt) == true) {
	    comp = AMBIGUOUS_COMP;
	    /* Don't count as either match or mismatch */
#endif
	  } else {
	    comp = MISMATCH_COMP;
	    nconsecutive_mismatches += 1;
	  }
	  debug10(printf("At querypos %d and pos %u on right, have %c and %c\n",
			 querypos,genomepos,c,g));
	  pairs = Pairpool_push(pairs,pairpool,querypos,genomepos,c,comp,g,g_alt,/*dynprogindex*/0);
	}
	querypos--;
	genomepos--;
	/* pos++; */
      }
    }
  }

  this->pairs = pairs;

  debug10(printf("RIGHT AFTER FILL\n"));
  debug10(Pair_dump_list(this->pairs,true));
  debug10(printf("END_RIGHT AFTER FILL\n"));

  Stage3_free_pairarray(&this);
  this->pairarray = make_pairarray(&this->npairs,&this->pairs,this->cdna_direction,
				   this->watsonp,pairpool,queryseq_ptr,
				   this->chroffset,this->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0);
  this->goodness = Pair_fracidentity_array(&this->matches,&this->unknowns,&this->mismatches,
					   &this->qopens,&this->qindels,&this->topens,&this->tindels,
					   &ncanonical,&nsemicanonical,&this->noncanonical,
					   &min_splice_prob,this->pairarray,this->npairs,this->cdna_direction);

  if (this->pairarray == NULL) {
    this->pairarray_freeable_p = false;
  } else {
    this->pairarray_freeable_p = true;
  }

  return;
}


void
Stage3_trim_right (T this, int goal, char *queryseq_ptr, Pairpool_T pairpool) {
  List_T path;
  Pair_T pair;

  int ncanonical, nsemicanonical;
  double min_splice_prob;


  debug10(printf("Entered Stage3_trim_right with goal %d\n",goal));
  debug10(printf("LEFT BEFORE TRIM\n"));
  debug10(Pair_dump_list(this->pairs,true));
  debug10(printf("END_LEFT BEFORE TRIM\n"));

  path = List_reverse(this->pairs);

  while (((Pair_T) path->first)->querypos > goal /* && pos <= this->chrhigh */) {
    path = Pairpool_pop(path,&pair);
  }

  this->pairs = List_reverse(path);

  debug10(printf("LEFT AFTER TRIM\n"));
  debug10(Pair_dump_list(this->pairs,true));
  debug10(printf("END_LEFT AFTER TRIM\n"));

  Stage3_free_pairarray(&this);
  this->pairarray = make_pairarray(&this->npairs,&this->pairs,this->cdna_direction,
				   this->watsonp,pairpool,queryseq_ptr,
				   this->chroffset,this->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0);
  this->goodness = Pair_fracidentity_array(&this->matches,&this->unknowns,&this->mismatches,
					   &this->qopens,&this->qindels,&this->topens,&this->tindels,
					   &ncanonical,&nsemicanonical,&this->noncanonical,
					   &min_splice_prob,this->pairarray,this->npairs,this->cdna_direction);

  if (this->pairarray == NULL) {
    this->pairarray_freeable_p = false;
  } else {
    this->pairarray_freeable_p = true;
  }

  return;
}


void
Stage3_trim_left (T this, int goal, char *queryseq_ptr, Pairpool_T pairpool) {

  List_T pairs;
  Pair_T pair;

  int ncanonical, nsemicanonical;
  double min_splice_prob;


  debug10(printf("Entered Stage3_trim_left with goal %d\n",goal));
  debug10(printf("RIGHT BEFORE TRIM\n"));
  debug10(Pair_dump_list(this->pairs,true));
  debug10(printf("END_RIGHT BEFORE TRIM\n"));


  /* Do not call insert_gapholders */
  pairs = this->pairs;
  
  while (((Pair_T) pairs->first)->querypos < goal) {
    pairs = Pairpool_pop(pairs,&pair);
  }

  this->pairs = pairs;

  debug10(printf("RIGHT AFTER TRIM\n"));
  debug10(Pair_dump_list(this->pairs,true));
  debug10(printf("END_RIGHT AFTER TRIM\n"));

  Stage3_free_pairarray(&this);
  this->pairarray = make_pairarray(&this->npairs,&this->pairs,this->cdna_direction,
				   this->watsonp,pairpool,queryseq_ptr,
				   this->chroffset,this->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0);
  this->goodness = Pair_fracidentity_array(&this->matches,&this->unknowns,&this->mismatches,
					   &this->qopens,&this->qindels,&this->topens,&this->tindels,
					   &ncanonical,&nsemicanonical,&this->noncanonical,
					   &min_splice_prob,this->pairarray,this->npairs,this->cdna_direction);

  if (this->pairarray == NULL) {
    this->pairarray_freeable_p = false;
  } else {
    this->pairarray_freeable_p = true;
  }

  return;
}



#if 0
static void
adjust_genomepos (T this, int delta) {
  Pair_T pair;
  List_T p;

  for (p = this->pairs; p != NULL; p = List_next(p)) {
    pair = (Pair_T) List_head(p);
    pair->genomepos += delta;
  }

  return;
}
#endif


static bool
merge_local_single (T this_left, T this_right,
		    int minpos1, int maxpos1, int minpos2, int maxpos2,
		    char *queryseq_ptr, char *queryuc_ptr,
		    Pairpool_T pairpool, Dynprog_T dynprogM,
		    int maxpeelback) {
  bool successp;
  Pair_T leftpair, rightpair;
  List_T path, orig_left_pairs, orig_right_pairs;
  bool watsonp, filledp;

  int ncanonical, nsemicanonical;
  double min_splice_prob;


#ifdef EXTRACT_GENOMICSEG
  char *genomicseg_ptr = NULL;
#endif
  int dynprogindex_minor = 0;


  orig_left_pairs = Pairpool_copy(this_left->pairs,pairpool);
  orig_right_pairs = Pairpool_copy(this_right->pairs,pairpool);

  this_left->pairs = Pair_clip_bounded_list_5(this_left->pairs,minpos1,maxpos1);
  this_right->pairs = Pair_clip_bounded_list_3(this_right->pairs,minpos2,maxpos2);

  /* Stage3_free_pairarray(&this_left); */
  /* Stage3_free_pairarray(&this_right); */

  if (this_left->pairs == NULL && this_right->pairs == NULL) {
    this_left->pairs = orig_left_pairs;
    this_right->pairs = orig_right_pairs;

#if 0
    this_left->pairarray = (struct Pair_T *) NULL;
    this_right->pairarray = (struct Pair_T *) NULL;
    this_left->pairarray_freeable_p = false;
    this_right->pairarray_freeable_p = false;
#endif
    return false;

  } else if ((watsonp = this_left->watsonp) == true) {
    debug10(printf("watsonp %d\n",watsonp));

#if 0
    /* Has no effect on plus strand */
    Pair_set_genomepos_list(this_left->pairs,chroffset,chrhigh,/*watsonp*/true);
    Pair_set_genomepos_list(this_right->pairs,chroffset,chrhigh,/*watsonp*/true);
#endif

    debug10(printf("LEFT\n"));
    debug10(Pair_dump_list(this_left->pairs,true));
    debug10(printf("END LEFT\n"));

    debug10(printf("RIGHT\n"));
    debug10(Pair_dump_list(this_right->pairs,true));
    debug10(printf("END RIGHT\n"));

#ifdef EXTRACT_GENOMICSEG
    firstpair = (Pair_T) List_head(this_left->pairs);
    lastpair = (Pair_T) List_last_value(this_right->pairs);
    firstpos = firstpair->genomepos;
    lastpos = lastpair->genomepos;
    left = this_left->chroffset + firstpos;
    genomicseg = Genome_get_segment(genome,left,genomiclength,/*chromosome_iit*/NULL,/*revcomp*/false);
    genomicseg_ptr = genomicuc_ptr = Sequence_fullpointer(genomicseg);
#endif

#if 0
    /* Has no effect on plus strand */
    Pair_set_genomepos_list(this_left->pairs,chroffset,chrhigh,/*watsonp*/true);
    Pair_set_genomepos_list(this_right->pairs,chroffset,chrhigh,/*watsonp*/true);
#endif

    debug10(printf("LEFT\n"));
    debug10(Pair_dump_list(this_left->pairs,true));
    debug10(printf("END LEFT\n"));

    debug10(printf("RIGHT\n"));
    debug10(Pair_dump_list(this_right->pairs,true));
    debug10(printf("END RIGHT\n"));

    path = List_reverse(this_left->pairs);

    leftpair = (Pair_T) path->first;
    rightpair = (Pair_T) this_right->pairs->first;

    debug10(printf("Running traverse_single_gap\n"));
    if ((this_right->pairs = traverse_single_gap(&filledp,&dynprogindex_minor,this_right->pairs,&path,leftpair,rightpair,
						 this_right->chroffset,this_right->chrhigh,
						 queryseq_ptr,queryuc_ptr,/*querylength*/0,watsonp,this_right->genestrand,
						 /*jump_late_p*/watsonp ? false : true,pairpool,dynprogM,
						 /*last_genomedp5*/NULL,/*last_genomedp3*/NULL,
						 maxpeelback,/*defect_rate*/0,/*forcep*/false,/*finalp*/true)) == NULL) {
      debug10(printf(" => failed\n"));
      successp = false;
    } else if (filledp == false) {
      debug10(printf(" => failed\n"));
      successp = false;
    } else {
      debug10(printf(" => succeeded\n"));
      successp = true;
    }
    this_left->pairs = List_reverse(path);

  } else {
    debug10(printf("watsonp %d\n",watsonp));

#if 0
    /* Do not change list, just pairarray */
    Pair_set_genomepos_list(this_left->pairs,chroffset,chrhigh,/*watsonp*/false);
    Pair_set_genomepos_list(this_right->pairs,chroffset,chrhigh,/*watsonp*/false);
#endif

    debug10(printf("LEFT\n"));
    debug10(Pair_dump_list(this_left->pairs,true));
    debug10(printf("END LEFT\n"));

    debug10(printf("RIGHT\n"));
    debug10(Pair_dump_list(this_right->pairs,true));
    debug10(printf("END RIGHT\n"));

#ifdef EXTRACT_GENOMICSEG
    firstpair = (Pair_T) List_head(this_left->pairs);
    lastpair = (Pair_T) List_last_value(this_right->pairs);
    firstpos = firstpair->genomepos;
    lastpos = lastpair->genomepos;
    left = this_right->chroffset + lastpos;
    genomicseg = Genome_get_segment(genome,left,genomiclength,/*chromosome_iit*/NULL,/*revcomp*/true);
    genomicseg_ptr = genomicuc_ptr = Sequence_fullpointer(genomicseg);
#endif

#if 0
    /* Do not change list, just pairarray */
    Pair_set_genomepos_list(this_left->pairs,chroffset,chrhigh,/*watsonp*/false);
    Pair_set_genomepos_list(this_right->pairs,chroffset,chrhigh,/*watsonp*/false);
#endif

    debug10(printf("LEFT\n"));
    debug10(Pair_dump_list(this_left->pairs,true));
    debug10(printf("END LEFT\n"));

    debug10(printf("RIGHT\n"));
    debug10(Pair_dump_list(this_right->pairs,true));
    debug10(printf("END RIGHT\n"));

    path = List_reverse(this_left->pairs);

    leftpair = (Pair_T) path->first;
    rightpair = (Pair_T) this_right->pairs->first;

    if ((this_right->pairs = traverse_single_gap(&filledp,&dynprogindex_minor,this_right->pairs,&path,leftpair,rightpair,
						 this_right->chroffset,this_right->chrhigh,
						 queryseq_ptr,queryuc_ptr,/*querylength*/0,watsonp,this_right->genestrand,
						 /*jump_late_p*/watsonp ? false : true,pairpool,dynprogM,
						 /*last_genomedp5*/NULL,/*last_genomedp3*/NULL,
						 maxpeelback,/*defect_rate*/0,/*forcep*/false,/*finalp*/true)) == NULL) {
      debug10(printf(" => failed\n"));
      successp = false;
    } else if (filledp == false) {
      debug10(printf(" => failed\n"));
      successp = false;
    } else {
      debug10(printf(" => succeeded\n"));
      successp = true;
    }
    this_left->pairs = List_reverse(path);
  }

  if (successp == false) {
    this_left->pairs = orig_left_pairs;
    this_left->pairarray = make_pairarray(&this_left->npairs,&this_left->pairs,this_left->cdna_direction,
					  this_left->watsonp,pairpool,queryseq_ptr,
					  this_left->chroffset,this_left->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0);
    this_left->goodness = Pair_fracidentity_array(&this_left->matches,&this_left->unknowns,&this_left->mismatches,
						  &this_left->qopens,&this_left->qindels,&this_left->topens,&this_left->tindels,
						  &ncanonical,&nsemicanonical,&this_left->noncanonical,
						  &min_splice_prob,this_left->pairarray,this_left->npairs,this_left->cdna_direction);

    this_right->pairs = orig_right_pairs;
    this_right->pairarray = make_pairarray(&this_right->npairs,&this_right->pairs,this_right->cdna_direction,
					   this_right->watsonp,pairpool,queryseq_ptr,
					   this_right->chroffset,this_right->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0);
    this_right->goodness = Pair_fracidentity_array(&this_right->matches,&this_right->unknowns,&this_right->mismatches,
						   &this_right->qopens,&this_right->qindels,&this_right->topens,&this_right->tindels,
						   &ncanonical,&nsemicanonical,&this_right->noncanonical,
						   &min_splice_prob,this_right->pairarray,this_right->npairs,this_right->cdna_direction);

  } else {
    this_left->pairs = List_append(this_left->pairs,this_right->pairs);
    this_right->pairs = (List_T) NULL;
  }


  debug10(printf(" => returning successp %d\n",successp));
  return successp;
}


static List_T
recompute_for_cdna_direction (int *cdna_direction, List_T pairs, int genestrand, bool watsonp,
#ifdef PMAP
			      char *queryaaseq_ptr,
#endif
			      char *queryseq_ptr, char *queryuc_ptr,
			      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
#ifndef GSNAP
			      Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool, 
#endif
			      Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			      int maxpeelback) {
  List_T pairs_fwd, path_fwd, pairs_rev, path_rev, copy;
  double max_intron_score_fwd = 0.0, max_intron_score_rev = 0.0,
    avg_donor_score_fwd = 0.0, avg_acceptor_score_fwd = 0.0,
    avg_donor_score_rev = 0.0, avg_acceptor_score_rev = 0.0;
  int nmatches_fwd, nmismatches_fwd, nindels_fwd,
    nknown_fwd, ncanonical_fwd, nsemicanonical_fwd, nnoncanonical_fwd, nbadintrons_fwd,
    nmatches_rev, nmismatches_rev, nindels_rev,
    nknown_rev, ncanonical_rev, nsemicanonical_rev, nnoncanonical_rev, nbadintrons_rev;
  int sensedir;

  double defect_rate_fwd, defect_rate_rev;

  copy = Pairpool_copy(pairs,pairpool);

  /* Compute fwd */
  path_fwd = path_compute_dir(&defect_rate_fwd,/*pairs*/copy,/*cdna_direction*/+1,watsonp,
			      genestrand,/*jump_late_p*/watsonp ? false : true,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      queryseq_ptr,queryuc_ptr,/*querylength*/0,chrnum,chroffset,chrhigh,
			      maxpeelback,
#ifndef GSNAP
			      oligoindices_minor,diagpool,cellpool, 
#endif
			      pairpool,dynprogL,dynprogM,dynprogR,
			      /*last_genomedp5*/NULL,/*last_genomedp3*/NULL/*,clean_ends_p:true*/);
  pairs_fwd = score_introns(&max_intron_score_fwd,&avg_donor_score_fwd,&avg_acceptor_score_fwd,
			    &nknown_fwd,&ncanonical_fwd,&nbadintrons_fwd,path_fwd,/*cdna_direction*/+1,watsonp,
			    chrnum,chroffset,chrhigh
#ifdef WASTE
			    ,pairpool
#endif
			    );
  /* alignment_score_fwd = */ score_alignment(&nmatches_fwd,&nmismatches_fwd,&nindels_fwd,
#ifdef COMPLEX_DIRECTION
					      &indel_alignment_score_fwd,
#endif
					      &nsemicanonical_fwd,&nnoncanonical_fwd,
					      pairs_fwd,/*cdna_direction*/+1);
    

  /* Compute rev */
  path_rev = path_compute_dir(&defect_rate_rev,/*pairs*/pairs,/*cdna_direction*/-1,watsonp,
			      genestrand,/*jump_late_p*/watsonp ? false : true,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      queryseq_ptr,queryuc_ptr,/*querylength*/0,chrnum,chroffset,chrhigh,
			      maxpeelback,
#ifndef GSNAP
			      oligoindices_minor,diagpool,cellpool, 
#endif
			      pairpool,dynprogL,dynprogM,dynprogR,
			      /*last_genomedp5*/NULL,/*last_genomedp3*/NULL/*,clean_ends_p:true*/);
  pairs_rev = score_introns(&max_intron_score_rev,&avg_donor_score_rev,&avg_acceptor_score_rev,
			    &nknown_rev,&ncanonical_rev,&nbadintrons_rev,path_rev,/*cdna_direction*/-1,watsonp,
			    chrnum,chroffset,chrhigh
#ifdef WASTE
			    ,pairpool
#endif			      
			    );
  /* alignment_score_rev = */ score_alignment(&nmatches_rev,&nmismatches_rev,&nindels_rev,
#ifdef COMPLEX_DIRECTION
					      &indel_alignment_score_rev,
#endif
					      &nsemicanonical_rev,&nnoncanonical_rev,
					      pairs_rev,/*cdna_direction*/-1);

  pairs = pick_cdna_direction(&(*cdna_direction),&sensedir,pairs_fwd,pairs_rev,
			      defect_rate_fwd,defect_rate_rev,
			      nknown_fwd,ncanonical_fwd,nsemicanonical_fwd,nnoncanonical_fwd,nbadintrons_fwd,
			      nknown_rev,ncanonical_rev,nsemicanonical_rev,nnoncanonical_rev,nbadintrons_rev,
			      max_intron_score_fwd,avg_donor_score_fwd,avg_acceptor_score_fwd,
			      max_intron_score_rev,avg_donor_score_rev,avg_acceptor_score_rev,
#ifdef COMPLEX_DIRECTION
			      nmatches_fwd,nmismatches_fwd,nmatches_rev,nmismatches_rev,nindels_fwd,nindels_rev,
			      indel_alignment_score_fwd,indel_alignment_score_rev,
#endif
			      /*sense_filter*/0);

  /* Don't know if we need to call path_compute_final */

  return pairs;
}


T
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
		    int maxpeelback) {
  T new_left, new_right;
  Pair_T end1, start2, leftpair, rightpair;
  List_T left_pairs, right_pairs, path;
  bool watsonp, filledp, shiftp, incompletep;
  int cdna_direction, cdna_direction_1, cdna_direction_2;
  bool make_dir_consistent_p;

  int intronlength, queryjump, genomejump;

  int dynprogindex_minor = 0, dynprogindex_major = 0;


  debug10(printf("Entered Stage3_merge_local with bounds1 %d..%d and bounds2 %d..%d\n",
		 minpos1,maxpos1,minpos2,maxpos2));

#ifdef DEBUG10
  printf("ORIG LEFT PAIRS:\n");
  Pair_dump_list(old_left->pairs,true);
  printf("ORIG RIGHT PAIRS:\n");
  Pair_dump_list(old_right->pairs,true);
#endif

  left_pairs = Pairpool_copy(old_left->pairs,pairpool);
  left_pairs = Pair_clip_bounded_list_5(left_pairs,minpos1,maxpos1);
  right_pairs = Pairpool_copy(old_right->pairs,pairpool);
  right_pairs = Pair_clip_bounded_list_3(right_pairs,minpos2,maxpos2);

#ifdef DEBUG10
  printf("CLIPPED LEFT PAIRS:\n");
  Pair_dump_list(left_pairs,true);
  printf("CLIPPED RIGHT PAIRS:\n");
  Pair_dump_list(right_pairs,true);
#endif


  path = clean_end_chimera(List_reverse(left_pairs));
  right_pairs = clean_end_chimera(right_pairs);

  Pairpool_clean_join(&path,&right_pairs);

  if (path == NULL || right_pairs == NULL) {
    /* Do not attach copies of pairs to this_left or this_right */
    return (T) NULL;

  } else {
    new_left = Stage3_new_from_pairs(left_pairs,old_left->cdna_direction,old_left->watsonp,
				     old_left->genestrand,old_left->sensedir,pairpool,queryseq,
				     /*query_subseq_offset*/Sequence_subseq_offset(queryseq),
				     old_left->chrnum,old_left->chroffset,old_left->chrhigh,old_left->chrlength);
    new_right = Stage3_new_from_pairs(right_pairs,old_right->cdna_direction,old_right->watsonp,
				      old_right->genestrand,old_right->sensedir,pairpool,queryseq,
				     /*query_subseq_offset*/Sequence_subseq_offset(queryseq),
				      old_right->chrnum,old_right->chroffset,old_right->chrhigh,old_right->chrlength);

    new_left->pairs = List_reverse(path);
    new_right->pairs = right_pairs;
  }

#ifdef DEBUG10
  printf("JOINED LEFT PAIRS:\n");
  Pair_dump_list(new_left->pairs,true);
  printf("JOINED RIGHT PAIRS:\n");
  Pair_dump_list(new_right->pairs,true);
#endif

  Stage3_free_pairarray(&new_left);
  Stage3_free_pairarray(&new_right);


  watsonp = new_left->watsonp;
#if 0
  if (watsonp == true) {
    debug10(printf("watsonp true\n"));
    
    firstpair = (Pair_T) List_head(new_left->pairs);
    lastpair = (Pair_T) List_last_value(new_right->pairs);
    firstpos = firstpair->genomepos;
    lastpos = lastpair->genomepos;
    left = new_left->chroffset + firstpos;

  } else {
    debug10(printf("watsonp false\n"));
    
    firstpair = (Pair_T) List_head(new_left->pairs);
    lastpair = (Pair_T) List_last_value(new_right->pairs);
    firstpos = firstpair->genomepos;
    lastpos = lastpair->genomepos;
    left = new_right->chroffset + lastpos;
  }
#endif

  /* Determine if need to make cdna_direction consistent */
  end1 = Pair_end_bound(&cdna_direction_1,new_left->pairs,/*breakpoint*/maxpos1);
  start2 = Pair_start_bound(&cdna_direction_2,new_right->pairs,/*breakpoint+1*/minpos2);
  debug10(printf("cdna_directions up to breakpoint are %d and %d\n",cdna_direction_1,cdna_direction_2));
  assert(end1 != NULL);
  assert(start2 != NULL);

  if (cdna_direction_1 > 0 && cdna_direction_2 < 0) {
    make_dir_consistent_p = true;

  } else if (cdna_direction_1 < 0 && cdna_direction_2 > 0) {
    make_dir_consistent_p = true;

  } else {
    make_dir_consistent_p = false;
    if (cdna_direction_1 == 0) {
      cdna_direction = cdna_direction_2;
    } else if (cdna_direction_2 == 0) {
      cdna_direction = cdna_direction_1;
    } else {
      cdna_direction = cdna_direction_1;
    }
    debug10(printf("cdna_direction is %d\n",cdna_direction));
  }


  /* Determine if the gap is an intron or not */
  end1 = (Pair_T) List_last_value(new_left->pairs);
  start2 = (Pair_T) List_head(new_right->pairs);
  queryjump = start2->querypos - end1->querypos - 1;
  genomejump = start2->genomepos - end1->genomepos - 1;
  intronlength = genomejump - queryjump;

  debug10(printf("intronlength %d = (start2->genomepos %d - end1->genomepos %d - 1) - (start2->querypos %d - end1->querypos %d - 1)\n",
		 intronlength,start2->genomepos,end1->genomepos,start2->querypos,end1->querypos));

  if (intronlength >= min_intronlength && splicingp == true) {
    debug10(printf("intronlength %d >= min_intronlength %d, so an intron\n",
		   intronlength,min_intronlength));
    /* Intron */
    path = List_reverse(new_left->pairs);
    leftpair = (Pair_T) path->first;
    rightpair = (Pair_T) new_right->pairs->first;

    if (make_dir_consistent_p == true) {
      /* Solve intron when re-computing for cdna_direction */
      debug10(printf("intron, but make dir consistent\n"));

      new_right->pairs = Pairpool_push_gapholder(new_right->pairs,pairpool,queryjump,genomejump,
						  /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      new_left->pairs = List_reverse(path);
      new_left->pairs = List_append(new_left->pairs,new_right->pairs);
      new_right->pairs = (List_T) NULL;

      new_left->pairs =
	recompute_for_cdna_direction(&cdna_direction,new_left->pairs,new_left->genestrand,watsonp,
#ifdef PMAP
				     queryaaseq_ptr,
#endif
				     queryseq_ptr,queryuc_ptr,
				     new_left->chrnum,new_left->chroffset,new_left->chrhigh,
#ifndef GSNAP
				     oligoindices_minor,diagpool,cellpool,
#endif
				     pairpool,dynprogL,dynprogM,dynprogR,maxpeelback);

    } else {
      debug10(printf("traverse_genome_gap with cdna_direction %d...",cdna_direction));
      new_right->pairs = traverse_genome_gap(&filledp,&shiftp,&dynprogindex_minor,&dynprogindex_major,
					      new_right->pairs,&path,leftpair,rightpair,
					      new_left->chrnum,new_left->chroffset,new_left->chrhigh,
					      queryseq_ptr,queryuc_ptr,/*querylength*/0,cdna_direction,
					      watsonp,new_right->genestrand,
					      /*jump_late_p*/watsonp ? false : true,pairpool,
					      dynprogL,dynprogM,dynprogR,/*last_genomedp5*/NULL,/*last_genomedp3*/NULL,
					      maxpeelback,/*defect_rate*/0,/*finalp*/true,/*simplep*/false);
      debug10(printf("done"));
      
      if (filledp == false) {
	new_right->pairs = Pairpool_push_gapholder(new_right->pairs,pairpool,queryjump,genomejump,
						    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      }
    
      new_left->pairs = List_reverse(path);
      new_left->pairs = List_append(new_left->pairs,new_right->pairs);
      new_right->pairs = (List_T) NULL;
    }

    if (make_pairarray_merge(new_left,cdna_direction,new_left->watsonp,pairpool,queryseq_ptr,
			     new_left->chroffset,new_left->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0,
			     /*new_gap_p*/true) == false) {
      Stage3_free(&new_left);
      Stage3_free(&new_right);
      return (T) NULL;
    }

  } else if (intronlength < 0) { /* Was intronlength < -EXTRAQUERYGAP, but this missed some short insertions */
    /* If traverse_cdna_gap fails, causes seg faults later on */
    /* cDNA gap */
    debug10(printf("cDNA gap, but make dir consistent\n"));

    path = List_reverse(new_left->pairs);
    leftpair = (Pair_T) path->first;
    rightpair = (Pair_T) new_right->pairs->first;

    if (make_dir_consistent_p == true) {
      /* Solve cDNA gap when re-computing for cdna_direction */
      new_right->pairs = Pairpool_push_gapholder(new_right->pairs,pairpool,queryjump,genomejump,
						  /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      new_left->pairs = List_reverse(path);
      new_left->pairs = List_append(new_left->pairs,new_right->pairs);
      new_right->pairs = (List_T) NULL;

      new_left->pairs =
	recompute_for_cdna_direction(&cdna_direction,new_left->pairs,new_left->genestrand,watsonp,
#ifdef PMAP
				     queryaaseq_ptr,
#endif
				     queryseq_ptr,queryuc_ptr,
				     new_left->chrnum,new_left->chroffset,new_left->chrhigh,
#ifndef GSNAP
				     oligoindices_minor,diagpool,cellpool,
#endif
				     pairpool,dynprogL,dynprogM,dynprogR,maxpeelback);

    } else {
      debug10(printf("traverse_cdna_gap..."));
      new_right->pairs = traverse_cdna_gap(&filledp,&incompletep,&dynprogindex_minor,&dynprogindex_major,
					    new_right->pairs,&path,leftpair,rightpair,
					    new_left->chroffset,new_left->chrhigh,
					    queryseq_ptr,queryuc_ptr,/*querylength*/0,watsonp,new_right->genestrand,
					    /*jump_late_p*/watsonp ? false : true,pairpool,
					    dynprogL,dynprogM,dynprogR,/*last_genomedp5*/NULL,/*last_genomedp3*/NULL,
					    maxpeelback,/*defect_rate*/0,/*finalp*/true);
      debug10(printf("done"));

      if (filledp == false) {
	new_right->pairs = Pairpool_push_gapholder(new_right->pairs,pairpool,queryjump,genomejump,
						    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      }

      new_left->pairs = List_reverse(path);
      new_left->pairs = List_append(new_left->pairs,new_right->pairs);
      new_right->pairs = (List_T) NULL;
    }

    if (make_pairarray_merge(new_left,cdna_direction,new_left->watsonp,pairpool,queryseq_ptr,
			     new_left->chroffset,new_left->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0,
			     /*new_gap_p*/true) == false) {
      Stage3_free(&new_left);
      Stage3_free(&new_right);
      return (T) NULL;
    }

  } else {
    /* Single gap */
    debug10(printf("intronlength %d, so a single gap\n",intronlength));
    debug10(printf("Before\n"));
    debug10(Pair_dump_list(new_left->pairs,true));
    debug10(Pair_dump_list(new_right->pairs,true));

    end1 = (Pair_T) List_last_value(new_left->pairs);
    start2 = (Pair_T) List_head(new_right->pairs);

    debug10(printf("Running merge_local_single\n"));
    if (merge_local_single(new_left,new_right,
			   minpos1,/*maxpos1*/end1->querypos,
			   /*minpos2*/start2->querypos,maxpos2,
			   queryseq_ptr,queryuc_ptr,
			   pairpool,dynprogM,maxpeelback) == false) {
      Stage3_free(&new_left);
      Stage3_free(&new_right);
      return (T) NULL;

    } else if (make_dir_consistent_p == true) {
      debug10(printf("Need to make dir consistent\n"));
      new_left->pairs =
	recompute_for_cdna_direction(&cdna_direction,new_left->pairs,new_left->genestrand,watsonp,
#ifdef PMAP
				     queryaaseq_ptr,
#endif
				     queryseq_ptr,queryuc_ptr,
				     new_left->chrnum,new_left->chroffset,new_left->chrhigh,
#ifndef GSNAP
				     oligoindices_minor,diagpool,cellpool,
#endif
				     pairpool,dynprogL,dynprogM,dynprogR,maxpeelback);
    }

    if (make_pairarray_merge(new_left,cdna_direction,new_left->watsonp,pairpool,queryseq_ptr,
			     new_left->chroffset,new_left->chrhigh,ngap,/*subseq_offset*/0,/*skiplength*/0,
			     /*new_gap_p*/false) == false) {
      Stage3_free(&new_left);
      Stage3_free(&new_right);
      return (T) NULL;
    }
    
    debug10(printf("After\n"));
    debug10(Pair_dump_list(new_left->pairs,true));
  }
      
  new_left->cdna_direction = cdna_direction;
  
  Stage3_free(&new_right);
  return new_left;
}


List_T
Stage3_split (T this, Sequence_T queryseq, Pairpool_T pairpool) {
  List_T split_objects = NULL;
  Stage3_T stage3;
  struct Pair_T *pairarray;
  int npairs;
  int goodness, matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  double min_splice_prob;

  List_T path, pairs, p;
  Pair_T pair;
  bool large_intron_p;

  debug(printf("\n** Starting Stage3_split with watsonp %d and cdna_direction %d\n",this->watsonp,this->cdna_direction));

  large_intron_p = false;
  for (p = this->pairs; p != NULL; p = p->rest) {
    pair = (Pair_T) p->first;
    if (pair->gapp == true && pair->genomejump > maxintronlen) {
      debug(printf("Found large intron of size %d\n",pair->genomejump));
      large_intron_p = true;
    }
  }

  if (large_intron_p == false) {
    return (List_T) NULL;

  } else {
    pairs = (List_T) NULL;
    path = List_reverse(this->pairs);
    while (path != NULL) {
      /* pairptr = path; */
      /* path = Pairpool_pop(path,&pair); */
      pair = (Pair_T) path->first;
      if (pair->gapp == false) {
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_transfer_one(pairs,&path);
#endif

      } else if (pair->genomejump <= maxintronlen) {
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_transfer_one(pairs,&path);
#endif

      } else {
	/* Start a new path */
	/* Pair_dump_list(pairs,true); */

	pairarray = make_pairarray(&npairs,&pairs,this->cdna_direction,this->watsonp,
				   pairpool,/*queryseq_ptr*/Sequence_fullpointer(queryseq),
				   this->chroffset,this->chrhigh,
				   ngap,/*query_subseq_offset*/Sequence_subseq_offset(queryseq),
				   /*skiplength*/Sequence_skiplength(queryseq));

	goodness = Pair_fracidentity_array(&matches,&unknowns,&mismatches,
					   &qopens,&qindels,&topens,&tindels,
					   &ncanonical,&nsemicanonical,&nnoncanonical,
					   &min_splice_prob,pairarray,npairs,this->cdna_direction);

	if ((stage3 = Stage3_new(pairarray,pairs,npairs,goodness,this->cdna_direction,this->sensedir,
				 matches,unknowns,mismatches,
				 qopens,qindels,topens,tindels,ncanonical,nsemicanonical,nnoncanonical,
				 this->chrnum,this->chroffset,this->chrhigh,this->chrlength,
				 this->watsonp,this->genestrand,
				 /*querylength*/Sequence_fulllength(queryseq),
				 /*skiplength*/Sequence_skiplength(queryseq),
				 /*trimlength*/Sequence_trimlength(queryseq),
				 this->straintype,this->strain,/*altstrain_iit*/NULL)) != NULL) {
	  split_objects = List_push(split_objects,(void *) stage3);
	}

	pairs = (List_T) NULL;
	path = path->rest;	/* Discard gap */
      }
    }

    /* Handle final path */
    /* Pair_dump_list(pairs,true); */

    pairarray = make_pairarray(&npairs,&pairs,this->cdna_direction,this->watsonp,
			       pairpool,/*queryseq_ptr*/Sequence_fullpointer(queryseq),
			       this->chroffset,this->chrhigh,
			       ngap,/*query_subseq_offset*/Sequence_subseq_offset(queryseq),
			       /*skiplength*/Sequence_skiplength(queryseq));

    goodness = Pair_fracidentity_array(&matches,&unknowns,&mismatches,
				       &qopens,&qindels,&topens,&tindels,
				       &ncanonical,&nsemicanonical,&nnoncanonical,
				       &min_splice_prob,pairarray,npairs,this->cdna_direction);

    if ((stage3 = Stage3_new(pairarray,pairs,npairs,goodness,this->cdna_direction,this->sensedir,
			     matches,unknowns,mismatches,
			     qopens,qindels,topens,tindels,ncanonical,nsemicanonical,nnoncanonical,
			     this->chrnum,this->chroffset,this->chrhigh,this->chrlength,
			     this->watsonp,this->genestrand,
			     /*querylength*/Sequence_fulllength(queryseq),
			     /*skiplength*/Sequence_skiplength(queryseq),
			     /*trimlength*/Sequence_trimlength(queryseq),
			     this->straintype,this->strain,/*altstrain_iit*/NULL)) != NULL) {
      split_objects = List_push(split_objects,(void *) stage3);
    }
  }

  return split_objects;
}


#ifndef PMAP
void
Stage3_guess_cdna_direction (T this) {
  this->cdna_direction = Pair_guess_cdna_direction_array(&this->sensedir,this->pairarray,this->npairs,
							 /*invertedp*/false,this->chroffset,this->watsonp);
  Pair_fix_cdna_direction_array(this->pairarray,this->npairs,this->cdna_direction);
  return;
}
#endif
