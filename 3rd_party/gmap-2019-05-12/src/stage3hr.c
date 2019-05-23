static char rcsid[] = "$Id: stage3hr.c 218689 2019-03-19 17:21:12Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage3hr.h"
#include "stage3hrdef.h"

#include <stdlib.h>		/* For qsort */
#include <string.h>
#include <strings.h>
#include <ctype.h>		/* For islower */
#include <math.h>		/* For exp() and log10() */
#include "assert.h"
#include "mem.h"
#include "univcoord.h"

#include "chrnum.h"
#include "complement.h"
#include "interval.h"
#include "univdiag.h"
#include "univdiagdef.h"
#include "substring.h"
#include "junction.h"
#include "genome128_hr.h"
#include "mapq.h"
#include "cigar.h"
#include "comp.h"		/* For Stage3end_run_gmap */
#include "maxent_hr.h"
#include "fastlog.h"
#include "transcript.h"
#include "kmer-search.h"



/* Scores for alts_status_inside */
#define ALTS_RESOLVED_BYLENGTH 0
#define ALTS_NOT_AMBIGUOUS 1


/* Eliminates distant splices if short splices are found */
/* #define DISTANT_SPLICE_SPECIAL 1 */

#define CONCORDANT_TEXT "concordant"
#define PAIRED_TEXT "paired"
#define UNPAIRED_TEXT "unpaired"

#ifdef USE_TALLY_RATIO
#define TALLY_RATIO 2.0
#endif

#define SUBSUMPTION_SLOP 10	/* Should allow for short insert lengths */
#define NMATCHES_SLOP 6
#define NMATCHES_TO_TRIMS_SLOP 9 /* Looser to allow for different splice options */
#define OUTERLENGTH_SLOP 1000
#define AMB_PENALTY 8		/* For long ambiguous ends */

#define SCORE_EVENTRIM_SLOP 2
#define SCORE_INDELS_EVENTRIM 1 /* Needed to compare genomic positions with and without indels */
#define EVENTRIM_BADINTRON_PENALTY 2
#define DO_FINAL 1


#ifdef CHECK_ASSERTIONS
#define CHECK_NMISMATCHES 1
#endif


#if 0
/* This is a bad idea.  Better to use nconcordant as a guide to stopping */
#define MAX_HITS 100		/* For evaluating concordance */
#endif

/* #define USE_ALLOCA_FOR_HITS 1 -- can lead to stack overflow */


/* Stage3end_new */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* transcript-guided alignment */
/* May want to turn on debug2 in transcript.c */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Stage3end_T comparisons.  Need to modify calls from path-solve.c */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Stage3pair_T comparisons */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif


/* Stage3end_nmatches_substrings */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif


/* Resolving insides */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* insert length calculation */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Stage3end_new_gmap check for chromosomal bounds */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif

/* circular chromosomes */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* substring_gmap */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif

/* Stage3_determine_pairtype */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

/* Stage3pair_overlap */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif



#define MAPQ_MAXIMUM_SCORE 40

static bool omit_concordant_uniq_p = false;
static bool omit_concordant_mult_p = false;
static bool ignore_trim_p = false;

/* static int kmer_search_sizelimit = 100; */

static int subopt_levels;

static bool want_random_p;
static bool transcriptomep;
static bool invert_first_p;
static bool invert_second_p;
static Genome_T genomecomp;
static Genome_T genomebits;
static Genome_T genomebits_alt;

static Univ_IIT_T chromosome_iit;
static int nchromosomes;
static int circular_typeint;

static Genome_T transcriptomebits;
static Transcriptome_T transcriptome;
static Univ_IIT_T transcript_iit;
static bool remap_transcriptome_p = false;

static IIT_T genes_iit;
static int *genes_divint_crosstable;
static int *genes_chrnum_crosstable;
static IIT_T tally_iit;
static int *tally_divint_crosstable;
static IIT_T runlength_iit;
static int *runlength_divint_crosstable;

static Chrpos_T pairmax_linear;
static Chrpos_T pairmax_circular;

static Chrpos_T expected_pairlength;
static Chrpos_T pairlength_deviation;

static Chrpos_T expected_pairlength_low;
static Chrpos_T expected_pairlength_high;
static Chrpos_T expected_pairlength_very_high;

static int localsplicing_penalty;
static int indel_penalty_middle;
static int antistranded_penalty;
static bool favor_multiexon_p;

static int ambig_end_interval;	/* For penalizing large ambiguous ends
				   in GMAP alignments, since such ends
				   should have been found */

static Chrpos_T overall_max_distance;
static int max_middle_insertions_default; /* If negative, then compute querylength - 2*min_indel_end_matches */
static int max_middle_deletions;

static bool *circularp;
static bool *altlocp;
static Univcoord_T *alias_starts;
static Univcoord_T *alias_ends;

static char *failedinput_root;
static Outputtype_T output_type;
static bool merge_samechr_p;
static bool method_print_p = false;


/* Probably not good to use in certain genomic regions, unless we also
   use known splicesites with distance information. */
/* But sometimes need to use to get correct mapping */
static bool favor_ambiguous_p;


void
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
		bool ignore_trim_p_in, bool omit_concordant_uniq_p_in, bool omit_concordant_mult_p_in,
		char *failedinput_root_in, Outputtype_T output_type_in, bool merge_samechr_p_in,
		bool method_print_p_in, bool want_random_p_in) {

  transcriptomep = transcriptomep_in;
  invert_first_p = invert_first_p_in;
  invert_second_p = invert_second_p_in;
  genomecomp = genomecomp_in;
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

  chromosome_iit = chromosome_iit_in;
  nchromosomes = nchromosomes_in;
  circular_typeint = circular_typeint_in;
  genes_iit = genes_iit_in;
  genes_divint_crosstable = genes_divint_crosstable_in;
  genes_chrnum_crosstable = genes_chrnum_crosstable_in;

  transcriptomebits = transcriptomebits_in;
  transcriptome = transcriptome_in;
  transcript_iit = transcript_iit_in;

  tally_iit = tally_iit_in;
  tally_divint_crosstable = tally_divint_crosstable_in;
  runlength_iit = runlength_iit_in;
  runlength_divint_crosstable = runlength_divint_crosstable_in;
  localsplicing_penalty = localsplicing_penalty_in;
  indel_penalty_middle = indel_penalty_middle_in;
  antistranded_penalty = antistranded_penalty_in;
  favor_multiexon_p = favor_multiexon_p_in;

  pairmax_linear = pairmax_linear_in;
  pairmax_circular = pairmax_circular_in;
  expected_pairlength = expected_pairlength_in;
  pairlength_deviation = pairlength_deviation_in;

  if (pairlength_deviation > expected_pairlength) {
    expected_pairlength_low = 0;
  } else {
    expected_pairlength_low = expected_pairlength - pairlength_deviation;
  }
  expected_pairlength_high = expected_pairlength + pairlength_deviation;
  expected_pairlength_very_high = expected_pairlength + 10*pairlength_deviation;

  if (distances_observed_p == true) {
    favor_ambiguous_p = false;
  } else {
    favor_ambiguous_p = true;
  }

#if 0
  ambig_end_interval = index1part + (index1interval - 1);
#else
  ambig_end_interval = 8;	/* Since GMAP uses 8-mers */
#endif

  subopt_levels = subopt_levels_in;

  max_middle_insertions_default = max_middle_insertions_in;
  max_middle_deletions = max_middle_deletions_in;

  overall_max_distance = shortsplicedist_in;
  if (max_middle_deletions > (int) overall_max_distance) {
    overall_max_distance = max_middle_deletions;
  }
  if (max_middle_insertions_default > (int) overall_max_distance) {
    overall_max_distance = max_middle_insertions_default;
  }

  circularp = circularp_in;
  altlocp = altlocp_in;
  alias_starts = alias_starts_in;
  alias_ends = alias_ends_in;

  failedinput_root = failedinput_root_in;

  ignore_trim_p = ignore_trim_p_in;
  omit_concordant_uniq_p = omit_concordant_uniq_p_in;
  omit_concordant_mult_p = omit_concordant_mult_p_in;

  output_type = output_type_in;
  merge_samechr_p = merge_samechr_p_in;
  method_print_p = method_print_p_in;
  want_random_p = want_random_p_in;

  return;
}



#define T Stage3end_T

Hittype_T
Stage3end_hittype (T this) {
  return this->hittype;
}

static char *
hittype_string (Hittype_T hittype) {
  switch (hittype) {
  case EXACT: return "exact";
  case SUB: return "sub";
  case HALFSPLICE_DONOR: return "donor";
  case HALFSPLICE_ACCEPTOR: return "acceptor";
  case SPLICE: return "splice";
  case SAMECHR_SPLICE: return "samechr_splice";
  case TRANSLOC_SPLICE: return "transloc_splice";
  case SUBSTRINGS: return "substrings";
  default: abort();
  }
}

char *
Stage3end_hittype_string (T this) {
  return hittype_string(this->hittype);
}

Method_T
Stage3end_method (T this) {
  return this->method;
}


int
Stage3end_genestrand (T this) {
  return this->genestrand;
}

bool
Stage3end_transcriptomep (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return false;
  } else if (this->method == TR) {
    return true;
  } else {
    return false;
  }
}

List_T
Stage3end_transcripts (T this) {
  return this->transcripts;
}

void
Stage3end_set_transcripts (T this, List_T transcripts) {
  List_free(&this->transcripts);
  this->transcripts = transcripts;
  return;
}

List_T
Stage3end_transcripts_other (T this) {
  return this->transcripts_other;
}


#if 0
void
Stage3end_transfer_transcripts (T dest, List_T sources) {
  List_T p, q;
  T source;
  Transcript_T transcript;

  for (p = sources; p != NULL; p = List_next(p)) {
    source = (T) List_head(p);
    debug2(printf("Transferring %d transcripts from %s to %s\n",
		  List_length(source->transcripts),hittype_string(source->hittype),hittype_string(dest->hittype)));
    for (q = source->transcripts; q != NULL; q = List_next(q)) {
      transcript = (Transcript_T) List_head(q);
      if (Transcript_in_list_p(transcript,dest->transcripts) == true) {
	Transcript_free(&transcript);
      } else {
	printf("Pushing onto transcripts %p,",dest->transcripts);
	dest->transcripts = List_push(dest->transcripts,(void *) transcript);
	printf(" now %p\n",dest->transcripts);
      }
    }
    List_free(&source->transcripts);
    debug2(Transcript_print_nums(dest->transcripts));
    debug2(printf("\n"));

    Stage3end_free(&source);
  }

  return;
}
#endif

#if 0
static void
Stage3end_transfer_transcripts_other (T dest, List_T sources) {
  List_T p, q;
  T source;
  Transcript_T transcript;

  for (p = sources; p != NULL; p = List_next(p)) {
    source = (T) List_head(p);
    for (q = source->transcripts; q != NULL; q = List_next(q)) {
      transcript = (Transcript_T) List_head(q);
      if (Transcript_in_list_p(transcript,dest->transcripts_other) == true) {
	Transcript_free(&transcript);
      } else {
	printf("Pushing onto transcripts %p,",dest->transcripts);
	dest->transcripts_other = List_push(dest->transcripts_other,(void *) transcript);
	printf(" now %p\n",dest->transcripts);
      }
    }
    List_free(&source->transcripts);
    Stage3end_free(&source);
  }

  return;
}
#endif


static void
Stage3end_transfer_transcripts_one (T dest, T source) {
  List_T q;
  Transcript_T transcript;

  debug2(printf("Transferring %d transcripts from %s to %s\n",
		List_length(source->transcripts),hittype_string(source->hittype),hittype_string(dest->hittype)));

  for (q = source->transcripts; q != NULL; q = List_next(q)) {
    transcript = (Transcript_T) List_head(q);
    if (Transcript_in_list_p(transcript,dest->transcripts) == true) {
      Transcript_free(&transcript);
    } else {
      dest->transcripts = List_push(dest->transcripts,(void *) transcript);
    }
  }
  List_free(&source->transcripts);

  for (q = source->transcripts_other; q != NULL; q = List_next(q)) {
    transcript = (Transcript_T) List_head(q);
    if (Transcript_in_list_p(transcript,dest->transcripts_other) == true) {
      Transcript_free(&transcript);
    } else {
      dest->transcripts_other = List_push(dest->transcripts_other,(void *) transcript);
    }
  }
  List_free(&source->transcripts_other);

  debug2(Transcript_print_nums(dest->transcripts));
  debug2(printf("\n"));

  return;
}

static void
Stage3pair_transfer_transcripts_one (Stage3pair_T dest, Stage3pair_T source) {
  List_T q;
  Transcript_T transcript;

  for (q = source->transcripts5; q != NULL; q = List_next(q)) {
    transcript = (Transcript_T) List_head(q);
    if (Transcript_in_list_p(transcript,dest->transcripts5) == true) {
      Transcript_free(&transcript);
    } else {
      dest->transcripts5 = List_push(dest->transcripts5,(void *) transcript);
    }
  }
  List_free(&source->transcripts5);

  for (q = source->transcripts3; q != NULL; q = List_next(q)) {
    transcript = (Transcript_T) List_head(q);
    if (Transcript_in_list_p(transcript,dest->transcripts3) == true) {
      Transcript_free(&transcript);
    } else {
      dest->transcripts3 = List_push(dest->transcripts3,(void *) transcript);
    }
  }
  List_free(&source->transcripts3);


  Stage3end_transfer_transcripts_one(dest->hit5,source->hit5);
  Stage3end_transfer_transcripts_one(dest->hit3,source->hit3);

  return;
}



bool
Stage3end_distant_splice_p (T this) {
  if (this->distant_splice_p == true) {
    return true;
  } else {
    return false;
  }
}


Chrnum_T
Stage3end_chrnum (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return 0;
  } else {
    return this->chrnum;
  }
}

Chrnum_T
Stage3end_effective_chrnum (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return 0;
  } else {
    return this->effective_chrnum;
  }
}

Chrnum_T
Stage3end_other_chrnum (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return 0;
  } else {
    return this->other_chrnum;
  }
}

Univcoord_T
Stage3end_chroffset (T this) {
  return this->chroffset;
}

Univcoord_T
Stage3end_chrhigh (T this) {
  return this->chrhigh;
}

Chrpos_T
Stage3end_chrlength (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return 0;
  } else {
    return this->chrlength;
  }
}

Chrpos_T
Stage3end_chrpos_low (T this) {
  return this->low - this->chroffset;
}

Chrpos_T
Stage3end_chrpos_high (T this) {
  return this->high - this->chroffset;
}


Univcoord_T
Stage3end_genomicstart (T this) {
  return this->genomicstart;
}

Univcoord_T
Stage3end_genomicend (T this) {
  return this->genomicend;
}

/* For Goby */
int
Stage3end_query_alignment_length (T this) {
  int length = 0;
  List_T p;
  Substring_T substring;
  Junction_T junction;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    length += Substring_match_length(substring);
  }
  for (p = this->junctions_LtoH; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == INS_JUNCTION) {
      length += Junction_nindels(junction);
    }
  }

  return length;
}

Chrpos_T
Stage3end_genomic_alignment_length (T this) {
  Chrpos_T length = 0;
  List_T p;
  Substring_T substring;
  Junction_T junction;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    length += Substring_genomic_alignment_length(substring);
  }
  for (p = this->junctions_LtoH; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == DEL_JUNCTION) {
      length += (Chrpos_T) Junction_nindels(junction);
    }
  }

  return length;
}


static Chrpos_T
Stage3end_chrpos_low_trim (T this) {
  Substring_T substring_low;
  List_T p;

  p = this->substrings_LtoH;
  substring_low = (Substring_T) List_head(p);
  if (Substring_has_alts_p(substring_low) == true) {
    p = List_next(p);
    substring_low = (Substring_T) List_head(p);
  }

  return Substring_chrpos_low(substring_low);
}

static Chrpos_T
Stage3end_chrpos_high_trim (T this) {
  Substring_T substring_high;
  List_T p;

  p = this->substrings_HtoL;
  substring_high = (Substring_T) List_head(p);
  if (Substring_has_alts_p(substring_high) == true) {
    p = List_next(p);
    substring_high = (Substring_T) List_head(p);
  }

  return Substring_chrpos_high(substring_high);
}

int
Stage3end_mapq_score (T this) {
  return this->mapq_score;
}

int
Stage3end_absmq_score (T this) {
  return this->absmq_score;
}

int
Stage3end_nmismatches_bothdiff (T this) {
  return this->nmismatches_bothdiff;
}

int
Stage3end_nmismatches_refdiff (T this) {
  return this->nmismatches_refdiff;
}

/* Called only for terminals */
Endtype_T
Stage3end_start_endtype (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_LtoH);
  return Substring_start_endtype(substring);
}

/* Called only for terminals */
Endtype_T
Stage3end_end_endtype (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_LtoH);
  return Substring_end_endtype(substring);
}

int
Stage3end_nindels (T this) {
  return this->nindels;
}

int
Stage3end_querylength (T this) {
  return this->querylength;
}

bool
Stage3end_plusp (T this) {
  return this->plusp;
}

bool
Stage3end_paired_usedp (T this) {
  return this->paired_usedp;
}

int
Stage3end_max_trim (T this) {
  if (this->trim_querystart > this->trim_queryend) {
    return this->trim_querystart;
  } else {
    return this->trim_queryend;
  }
}


static int
start_amb_length (T this) {
  return Substring_start_amb_length((Substring_T) List_head(this->substrings_1toN));
}

static int
end_amb_length (T this) {
  return Substring_end_amb_length((Substring_T) List_head(this->substrings_Nto1));
}

#if 0
static int
n_amb_ends (T this) {
  int n = 0;

  if (start_amb_length(this) > 0) {
    n++;
  }
  if (end_amb_length(this) > 0) {
    n++;
  }

  return n;
}
#endif


#ifdef DEBUG8
static int
amb_length (T this) {
  return Substring_start_amb_length((Substring_T) List_head(this->substrings_1toN)) +
    Substring_end_amb_length((Substring_T) List_head(this->substrings_Nto1));
}
#endif


/* Two types of ambiguity: known amb (mapped to >1 genomic place) and unknown amb (splice site seen) */
static bool
known_ambiguous_p (T this) {
  if (Substring_ambiguous_p((Substring_T) List_head(this->substrings_1toN))) {
    return true;
  } else if (Substring_ambiguous_p((Substring_T) List_head(this->substrings_Nto1))) {
    return true;
  } else {
    return false;
  }
}


/* Includes amb and non-amb */
int
Stage3end_total_trim (T this) {
  return this->trim_querystart + this->trim_queryend;
}


int
Stage3end_circularpos (T this) {
  return this->circularpos;
}


Junction_T
Stage3end_junctionD (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Junction_T) List_head(this->junctions_Nto1);
  } else {
    return (Junction_T) List_head(this->junctions_1toN);
  }
}

Junction_T
Stage3end_junctionA (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Junction_T) List_head(this->junctions_1toN);
  } else {
    return (Junction_T) List_head(this->junctions_Nto1);
  }
}

List_T
Stage3end_substrings_LtoH (T this) {
  return this->substrings_LtoH;
}

List_T
Stage3end_junctions_LtoH (T this) {
  return this->junctions_LtoH;
}


/* Called only by samprint currently */
Substring_T
Stage3end_substring1 (T this) {
  return (Substring_T) List_head(this->substrings_1toN);
}

/* Called only by samprint currently */
Substring_T
Stage3end_substringN (T this) {
  return (Substring_T) List_head(this->substrings_Nto1);
}


Substring_T
Stage3end_substring_for_concordance (T this, bool first_read_p) {
  if (first_read_p == true) {
    return (Substring_T) List_head(this->substrings_Nto1);
  } else {
    return (Substring_T) List_head(this->substrings_1toN);
  }
}

Substring_T
Stage3end_substring_other (T this, bool first_read_p) {
  if (first_read_p == true) {
    return (Substring_T) List_head(this->substrings_1toN);
  } else {
    return (Substring_T) List_head(this->substrings_Nto1);
  }
}


bool
Stage3end_donor_concordant_p (T this, bool first_read_p) {
  if (this->sensedir != SENSE_ANTI) {
    if (first_read_p == true) {
      return false;
    } else {
      return true;
    }
  } else {
    if (first_read_p == true) {
      return true;
    } else {
      return false;
    }
  }
}


Substring_T
Stage3end_substring_donor (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Substring_T) List_head(this->substrings_Nto1);
  } else if (this->sensedir == SENSE_FORWARD) {
    return (Substring_T) List_head(this->substrings_1toN);
  } else {
    fprintf(stderr,"sensedir is SENSE_NULL in Stage3end_substring_donor\n");
    abort();
  }
}

Substring_T
Stage3end_substring_acceptor (T this) {
  if (this->sensedir == SENSE_ANTI) { 
    return (Substring_T) List_head(this->substrings_1toN);
  } else if (this->sensedir == SENSE_FORWARD) {
    return (Substring_T) List_head(this->substrings_Nto1);
  } else {
    fprintf(stderr,"sensedir is SENSE_NULL in Stage3end_substring_acceptor\n");
    abort();
  }
}

/* Now same as Stage3end_substring_donor */
Substring_T
Stage3end_substringD (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Substring_T) List_head(this->substrings_Nto1);
  } else {
    return (Substring_T) List_head(this->substrings_1toN);
  }
}

/* Now same as Stage3end_substring_acceptor */
Substring_T
Stage3end_substringA (T this) {
  if (this->sensedir == SENSE_ANTI) {
    return (Substring_T) List_head(this->substrings_1toN);
  } else {
    return (Substring_T) List_head(this->substrings_Nto1);
  }
}


Substring_T
Stage3end_substringS (T this) {
  return (Substring_T) List_head(List_next(this->substrings_1toN));
}



/* Same logic as in print_substrings in samprint.c to get the first substring for CIGAR or MD string */
Substring_T
Stage3end_substring_low (T this, int hardclip_low) {
  List_T p;

  debug15(printf("Entered Stage3end_substring_low\n"));

  if (this == NULL) {
    return (Substring_T) NULL;

  } else if (this->plusp == true) {
    p = this->substrings_LtoH;
    if (Substring_has_alts_p((Substring_T) List_head(p)) == true) {
      p = List_next(p);
    }
    while (p != NULL && Substring_queryend((Substring_T) List_head(p)) <= hardclip_low) {
      debug15(printf("Plus: Skipping substring %d..%d against hardclip_low %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     hardclip_low));
      p = List_next(p);
    }

    if (p == NULL) {
      return (Substring_T) NULL;
    } else {
      debug15(printf("Plus: Returning substring %d..%d against hardclip_low %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     hardclip_low));
      return (Substring_T) List_head(p);
    }

  } else {
#ifdef DEBUG15
    for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
      printf("LtoH: %d..%d\n",
	     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)));
    }
#endif

    p = this->substrings_LtoH;
    if (Substring_has_alts_p((Substring_T) List_head(p)) == true) {
      p = List_next(p);
    }

    while (p != NULL && Substring_querystart((Substring_T) List_head(p)) >= this->querylength - hardclip_low) {
      debug15(printf("Minus: Skipping substring %d..%d against %d = querylength %d - hardclip_low %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     this->querylength - hardclip_low,this->querylength,hardclip_low));
      p = List_next(p);
    }

    if (p == NULL) {
      return (Substring_T) NULL;
    } else {
      debug15(printf("Minus: Returning substring %d..%d against %d = querylength %d - hardclip_low %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     this->querylength - hardclip_low,this->querylength,hardclip_low));
      return (Substring_T) List_head(p);
    }
  }
}


#if 0
/* Modified from Stage3end_substring_low */
Substring_T
Stage3end_substring_high (T this, int hardclip_high) {
  List_T p;

  debug15(printf("Entered Stage3end_substring_high\n"));

  if (this == NULL) {
    return (Substring_T) NULL;

  } else if (this->plusp == true) {
    p = this->substrings_HtoL;
    if (Substring_has_alts_p((Substring_T) List_head(p)) == true) {
      p = List_next(p);
    }

    while (p != NULL && Substring_querystart((Substring_T) List_head(p)) >= this->querylength - hardclip_high) {
      debug15(printf("Plus: Skipping substring %d..%d against %d = querylength %d - hardclip_high %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     this->querylength - hardclip_high,this->querylength,hardclip_high));
      p = List_next(p);
    }

    if (p == NULL) {
      return (Substring_T) NULL;
    } else {
      debug15(printf("Plus: Returning substring %d..%d against %d = querylength %d - hardclip_high %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     this->querylength - hardclip_high,this->querylength,hardclip_high));
      return (Substring_T) List_head(p);
    }

  } else {
#ifdef DEBUG15
    for (p = this->substrings_HtoL; p != NULL; p = List_next(p)) {
      printf("HtoL: %d..%d\n",
	     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)));
    }
#endif

    p = this->substrings_HtoL;
    if (Substring_has_alts_p((Substring_T) List_head(p)) == true) {
      p = List_next(p);
    }

    while (p != NULL && Substring_queryend((Substring_T) List_head(p)) <= hardclip_high) {
      debug15(printf("Minus: Skipping substring %d..%d against hardclip_high %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     hardclip_high));
      p = List_next(p);
    }

    if (p == NULL) {
      return (Substring_T) NULL;
    } else {
      debug15(printf("Minus: Returning substring %d..%d against hardclip_high %d\n",
		     Substring_querystart((Substring_T) List_head(p)),Substring_queryend((Substring_T) List_head(p)),
		     hardclip_high));
      return (Substring_T) List_head(p);
    }
  }
}
#endif



Substring_T
Stage3end_substring_containing (T this, int querypos) {
  Substring_T substring;
  List_T p;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_contains_p(substring,querypos) == true) {
      return substring;
    }
  }
  return (Substring_T) NULL;
}


double
Stage3end_min_evalue (T this) {
  double min_evalue = 1000.0, evalue;
  Substring_T substring;
  List_T p;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if ((evalue = Substring_evalue(substring)) < min_evalue) {
      min_evalue = evalue;
    }
  }

  return min_evalue;
}


Chrpos_T
Stage3end_distance (T this) {
  return this->distance;
}

double
Stage3end_chimera_prob (T this) {
  List_T p;
  Junction_T junction;

  for (p = this->junctions_1toN; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == CHIMERA_JUNCTION) {
      return Junction_prob(junction);
    }
  }

  return 0.0;
}

static double
Stage3end_prob (T this) {
  double prob = 0.0;
  List_T p;
  Junction_T junction;

  for (p = this->junctions_LtoH; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    prob += Junction_prob(junction);
  }

  return prob;
}


/* Should eventually look for substrings adjacent to the chimeric junction */
Univcoord_T
Stage3end_chimera_segmenti_left (T this) {
  Univcoord_T x_segmenti, x_segmentj;
  Substring_T substring_donor, substring_acceptor;

  if (this->sensedir == SENSE_ANTI) {
    substring_donor = (Substring_T) List_head(this->substrings_Nto1);
    substring_acceptor = (Substring_T) List_head(this->substrings_1toN);
  } else {
    substring_donor = (Substring_T) List_head(this->substrings_1toN);
    substring_acceptor = (Substring_T) List_head(this->substrings_Nto1);
  }

  x_segmenti = Substring_left_genomicseg(substring_donor);
  x_segmentj = Substring_left_genomicseg(substring_acceptor);
  if (x_segmenti < x_segmentj) {
    return x_segmenti;
  } else {
    return x_segmentj;
  }
}  

/* Should eventually look for substrings adjacent to the chimeric junction */
Univcoord_T
Stage3end_chimera_segmentj_left (T this) {
  Univcoord_T x_segmenti, x_segmentj;
  Substring_T substring_donor, substring_acceptor;

  if (this->sensedir == SENSE_ANTI) {
    substring_donor = (Substring_T) List_head(this->substrings_Nto1);
    substring_acceptor = (Substring_T) List_head(this->substrings_1toN);
  } else {
    substring_donor = (Substring_T) List_head(this->substrings_1toN);
    substring_acceptor = (Substring_T) List_head(this->substrings_Nto1);
  }

  x_segmenti = Substring_left_genomicseg(substring_donor);
  x_segmentj = Substring_left_genomicseg(substring_acceptor);
  if (x_segmenti > x_segmentj) {
    return x_segmenti;
  } else {
    return x_segmentj;
  }
}  


int
Stage3end_chimera_segmenti_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoord_T x_segmenti, x_segmentj, y_segmenti, y_segmentj, temp;
  Substring_T x_substring_donor, x_substring_acceptor,
    y_substring_donor, y_substring_acceptor;

  if (x->sensedir == SENSE_ANTI) {
    x_substring_donor = (Substring_T) List_head(x->substrings_Nto1);
    x_substring_acceptor = (Substring_T) List_head(x->substrings_1toN);
  } else {
    x_substring_donor = (Substring_T) List_head(x->substrings_1toN);
    x_substring_acceptor = (Substring_T) List_head(x->substrings_Nto1);
  }

  if (y->sensedir == SENSE_ANTI) {
    y_substring_donor = (Substring_T) List_head(y->substrings_Nto1);
    y_substring_acceptor = (Substring_T) List_head(y->substrings_1toN);
  } else {
    y_substring_donor = (Substring_T) List_head(y->substrings_1toN);
    y_substring_acceptor = (Substring_T) List_head(y->substrings_Nto1);
  }

  x_segmenti = Substring_left_genomicseg(x_substring_donor);
  x_segmentj = Substring_left_genomicseg(x_substring_acceptor);
  if (x_segmentj < x_segmenti) {
    temp = x_segmentj;
    x_segmentj = x_segmenti;
    x_segmenti = temp;
  }

  y_segmenti = Substring_left_genomicseg(y_substring_donor);
  y_segmentj = Substring_left_genomicseg(y_substring_acceptor);
  if (y_segmentj < y_segmenti) {
    temp = y_segmentj;
    y_segmentj = y_segmenti;
    y_segmenti = temp;
  }

  if (x_segmenti < y_segmenti) {
    return -1;
  } else if (y_segmenti < x_segmenti) {
    return +1;
  } else if (x_segmentj > y_segmentj) {
    return -1;
  } else if (y_segmentj > x_segmentj) {
    return +1;
  } else {
    return 0;
  }
}



int
Stage3end_chimera_segmentj_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoord_T x_segmenti, x_segmentj, y_segmenti, y_segmentj, temp;
  Substring_T x_substring_donor, x_substring_acceptor,
    y_substring_donor, y_substring_acceptor;

  if (x->sensedir == SENSE_ANTI) {
    x_substring_donor = (Substring_T) List_head(x->substrings_Nto1);
    x_substring_acceptor = (Substring_T) List_head(x->substrings_1toN);
  } else {
    x_substring_donor = (Substring_T) List_head(x->substrings_1toN);
    x_substring_acceptor = (Substring_T) List_head(x->substrings_Nto1);
  }

  if (y->sensedir == SENSE_ANTI) {
    y_substring_donor = (Substring_T) List_head(y->substrings_Nto1);
    y_substring_acceptor = (Substring_T) List_head(y->substrings_1toN);
  } else {
    y_substring_donor = (Substring_T) List_head(y->substrings_1toN);
    y_substring_acceptor = (Substring_T) List_head(y->substrings_Nto1);
  }


  x_segmenti = Substring_left_genomicseg(x_substring_donor);
  x_segmentj = Substring_left_genomicseg(x_substring_acceptor);
  if (x_segmentj < x_segmenti) {
    temp = x_segmentj;
    x_segmentj = x_segmenti;
    x_segmenti = temp;
  }

  y_segmenti = Substring_left_genomicseg(y_substring_donor);
  y_segmentj = Substring_left_genomicseg(y_substring_acceptor);
  if (y_segmentj < y_segmenti) {
    temp = y_segmentj;
    y_segmentj = y_segmenti;
    y_segmenti = temp;
  }

  if (x_segmentj < y_segmentj) {
    return -1;
  } else if (y_segmentj < x_segmentj) {
    return +1;
  } else if (x_segmenti > y_segmenti) {
    return -1;
  } else if (y_segmenti > x_segmenti) {
    return +1;
  } else {
    return 0;
  }
}


int
Stage3end_sensedir (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return SENSE_NULL;
  } else {
    return this->sensedir;
  }
}

#if 0
int
Stage3end_cdna_direction (T this) {
  if (this == NULL) {
    return SENSE_NULL;
  } else if (this->sensedir == SENSE_FORWARD) {
    return +1;
  } else if (this->sensedir == SENSE_ANTI) {
    return -1;
  } else {
    return SENSE_NULL;
  }
}
#endif

#if 0
bool
Stage3end_start_ambiguous_p (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  return Substring_ambiguous_p(substring);
}
#endif

#if 0
bool
Stage3end_end_ambiguous_p (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  return Substring_ambiguous_p(substring);
}
#endif

bool
Stage3end_start_has_alts_p (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  return Substring_has_alts_p(substring);
}

bool
Stage3end_end_has_alts_p (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  return Substring_has_alts_p(substring);
}


Univcoord_T *
Stage3end_start_alts_coords (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  if (Substring_has_alts_p(substring) == false) {
    return (Univcoord_T *) NULL;
  } else {
    return Substring_alts_coords(substring);
  }
}

Univcoord_T *
Stage3end_end_alts_coords (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  if (Substring_has_alts_p(substring) == false) {
    return (Univcoord_T *) NULL;
  } else {
    return Substring_alts_coords(substring);
  }
}

int
Stage3end_start_alts_ncoords (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  if (Substring_has_alts_p(substring) == false) {
    return 0;
  } else {
    return Substring_alts_ncoords(substring);
  }
}

int
Stage3end_end_alts_ncoords (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  if (Substring_has_alts_p(substring) == false) {
    return 0;
  } else {
    return Substring_alts_ncoords(substring);
  }
}


int
Stage3end_substrings_querystart (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  return Substring_querystart(substring);
}

int
Stage3end_substrings_queryend (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  return Substring_queryend(substring);
}

#ifdef RESOLVE_INSIDE_GENERAL
static int
Stage3end_querystart (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_1toN);
  return Substring_querystart(substring);
}
#endif

#ifdef RESOLVE_INSIDE_GENERAL
static int
Stage3end_queryend (T this) {
  Substring_T substring;

  substring = (Substring_T) List_head(this->substrings_Nto1);
  return Substring_queryend(substring);
}
#endif

int
Stage3end_trimlength (T this) {
  return this->trim_querystart + this->trim_queryend;
}


void
Stage3end_count_hits (int *npaths_primary, int *npaths_altloc, List_T hits) {
  T hit;

  *npaths_primary = *npaths_altloc = 0;

  while (hits != NULL) {
    hit = (T) List_head(hits);
    if (altlocp[hit->chrnum] == true) {
      *npaths_altloc += 1;
    } else {
      *npaths_primary += 1;
    }
    hits = List_next(hits);
  }

  return;
}

static long int
Stage3end_compute_tally (T this) {
  long int tally = 0L;
  List_T p;
  Substring_T substring;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    tally += Substring_tally(substring,tally_iit,tally_divint_crosstable);
  }

  return tally;
}

static bool
Stage3end_runlength_p (T this) {
  List_T p;
  Substring_T substring;

  for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    if (Substring_runlength_p(substring,runlength_iit,runlength_divint_crosstable) == true) {
      return true;
    }
  }

  return false;
}


void
Stage3end_free (T *old) {
  List_T p;
  Substring_T substring;
  Junction_T junction;


  if (*old != NULL) {
    debug0(printf("Freeing Stage3end %p from method %s\n",*old,Method_string((*old)->method)));

    if ((*old)->transcripts_other != NULL) {
      Transcript_gc(&(*old)->transcripts_other);
    }
    if ((*old)->transcripts != NULL) {
      Transcript_gc(&(*old)->transcripts);
    }

    for (p = (*old)->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      Substring_free(&substring);
    }
    /* List_free(&(*old)->substrings_1toN); -- allocated by Listpool_push */
    /* List_free(&(*old)->substrings_Nto1); -- allocated by Listpool_push */
    /* List_free(&(*old)->substrings_LtoH); -- allocated by Listpool_push */
    /* List_free(&(*old)->substrings_HtoL); -- allocated by Listpool_push */

    for (p = (*old)->junctions_1toN; p != NULL; p = List_next(p)) {
      junction = (Junction_T) List_head(p);
      Junction_free(&junction);
    }
    /* List_free(&(*old)->junctions_1toN); -- allocated by Listpool_push */
    /* List_free(&(*old)->junctions_Nto1); -- allocated by Listpool_push */
    /* List_free(&(*old)->junctions_LtoH); -- allocated by Listpool_push */
    /* List_free(&(*old)->junctions_HtoL); */

    FREE_OUT(*old);
  }

  return;
}


/* Used for freeing list contents in Concordance_pair_up procedures */
/* Do not free the list itself, though, which was previously freed in
   stage1hr.c, and now allocated by Hitlistpool_T */
void
Stage3end_gc (List_T values) {
  List_T p;
  T hit;

  for (p = values; p != NULL; p = p->rest) {
    if ((hit = (T) p->first) != NULL) {
      Stage3end_free(&hit);
    }
  }
  Hitlist_free(&values);
  return;
}



bool
Stage3pair_distant_splice_p (Stage3pair_T this) {
  if (this->hit5 != NULL && this->hit5->distant_splice_p == true) {
    return true;
  } else if (this->hit3 != NULL && this->hit3->distant_splice_p == true) {
    return true;
  } else {
    return false;
  }
}


int
Stage3pair_genestrand (Stage3pair_T this) {
  return this->genestrand;
}

Stage3end_T
Stage3pair_hit5 (Stage3pair_T this) {
  return this->hit5;
}

Stage3end_T
Stage3pair_hit3 (Stage3pair_T this) {
  return this->hit3;
}

int
Stage3pair_mapq_score (Stage3pair_T this) {
  return this->mapq_score;
}

int
Stage3pair_absmq_score (Stage3pair_T this) {
  return this->absmq_score;
}

List_T
Stage3pair_transcripts5 (Stage3pair_T this) {
  return this->transcripts5;
}

List_T
Stage3pair_transcripts3 (Stage3pair_T this) {
  return this->transcripts3;
}

Chrpos_T
Stage3pair_pairlength (Stage3pair_T this) {
  return this->insertlength;
}

int
Stage3pair_relationship (Stage3pair_T this) {
  return this->pair_relationship;
}

int
Stage3pair_total_trim (Stage3pair_T this) {
  return Stage3end_total_trim(this->hit5) + Stage3end_total_trim(this->hit3);
}

int
Stage3pair_max_trim (Stage3pair_T this) {
  int trim5, trim3;
  T hit;

#if 0
  /* Don't want ambiguous ends for purpose of defining concordant terminals */
  trim5 = Stage3end_total_trim(this->hit5);
  trim3 = Stage3end_total_trim(this->hit3);
#else
  hit = this->hit5;
  trim5 = hit->trim_querystart + hit->trim_queryend;
  hit = this->hit3;
  trim3 = hit->trim_querystart + hit->trim_queryend;
#endif

  if (trim5 > trim3) {
    return trim5;
  } else {
    return trim3;
  }
}

int
Stage3pair_nmatches_to_trims (int *nmatches5, int *nmatches3, Stage3pair_T this) {
  *nmatches5 = this->hit5->nmatches_to_trims;
  *nmatches3 = this->hit3->nmatches_to_trims;
  return this->nmatches_to_trims;
}


bool
Stage3pair_concordantp (List_T hitpairs) {
  List_T p;
  Stage3pair_T hitpair;

  for (p = hitpairs; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
#if 0
    /* Not necessary, since we are getting the result after GMAP align pair */
    if (Stage3_determine_pairtype(hitpair->hit5,hitpair->hit3,hitpair) == CONCORDANT) {
      return true;
    }
#else
    if (hitpair->pairtype == CONCORDANT) {
      return true;
    }
#endif
  }
  return false;
}

void
Stage3pair_count_hits (int *npaths_primary, int *npaths_altloc, List_T hitpairs) {
  Stage3pair_T hitpair;

  *npaths_primary = *npaths_altloc = 0;

  while (hitpairs != NULL) {
    hitpair = (Stage3pair_T) List_head(hitpairs);
    if (altlocp[hitpair->hit5->chrnum] == true) {
      *npaths_altloc += 1;
    } else if (altlocp[hitpair->hit3->chrnum] == true) {
      *npaths_altloc += 1;
    } else {
      *npaths_primary += 1;
    }
    hitpairs = List_next(hitpairs);
  }

  return;
}

List_T
Stage3pair_filter_nonconcordant (List_T hitpairs, Hitlistpool_T hitlistpool) {
  List_T filtered = NULL, p;
  Stage3pair_T hitpair;

  for (p = hitpairs; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    if (hitpair->pairtype != CONCORDANT) {
      Stage3pair_free(&hitpair);
    } else {
      filtered = Hitlist_push(filtered,hitlistpool,(void *) hitpair);
    }
  }
  Hitlist_free(&hitpairs);
  return filtered;
}


/* Returns true if ilengths are valid */
static bool
find_ilengths (int *ilength_low, int *ilength_high, Stage3end_T hit, Univcoord_T common_genomicpos) {
  List_T p, q;
  Substring_T substring;
  Junction_T junction;


  debug15(printf("Finding ilengths for common_genomicpos %u\n",(Chrpos_T) (common_genomicpos - chroffset)));
  if (hit->plusp == true) {
#ifdef DEBUG15
    printf("plus.  Checking common genomicpos %llu against\n",common_genomicpos - hit->chroffset);
    for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      printf("substring %p: %u..%u, trim %d..%d\n",
	     substring,Substring_alignstart_trim(substring) - hit->chroffset,
	     Substring_alignend_trim(substring) - 1U - hit->chroffset,
	     Substring_trim_querystart(substring),Substring_trim_queryend(substring));
    }
    printf("\n");
#endif
    /* Plus: Subtract 1 from alignend */
    *ilength_low = 0;
    for (p = hit->substrings_1toN, q = hit->junctions_1toN; p != NULL; p = List_next(p), q = List_next(q)) {
      substring = (Substring_T) List_head(p);
      debug15(printf("substring %p: %u..%u, trim %d..%d\n",substring,
		     Substring_alignstart_trim(substring) - hit->chroffset,
		     Substring_alignend_trim(substring) - 1U - hit->chroffset,
		     Substring_trim_querystart(substring),Substring_trim_queryend(substring)));
      if (Substring_overlap_point_trimmed_p(substring,common_genomicpos) == false) {
	*ilength_low += Substring_genomic_alignment_length(substring);
	if (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_low += Junction_nindels(junction);
	  }
	}

      } else {
	*ilength_low += (common_genomicpos - Substring_alignstart_trim(substring) + 1);
	*ilength_high = ((Substring_alignend_trim(substring) - 1) - common_genomicpos + 1);
	p = List_next(p);
	while (p != NULL) {
	  substring = (Substring_T) List_head(p);
	  *ilength_high += Substring_genomic_alignment_length(substring);
	  p = List_next(p);
	}
	while (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_high += Junction_nindels(junction);
	  }
	  q = List_next(q);
	}
	debug15(printf("Plus: Have ilength_low %d and ilength_high %d\n",*ilength_low,*ilength_high));
	return true;
      }
    }
  } else {
#ifdef DEBUG15
    printf("minus.  Checking common genomicpos %llu against\n",common_genomicpos - hit->chroffset);
    for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      printf("substring %p: %u..%u, trim %d..%d\n",
	     substring,Substring_alignstart_trim(substring) - hit->chroffset,
	     Substring_alignend_trim(substring) - 1U - hit->chroffset,
	     Substring_trim_querystart(substring),Substring_trim_queryend(substring));
    }
    printf("\n");
#endif
    /* Minus: Subtract 1 from alignstart */
    *ilength_high = 0;
    for (p = hit->substrings_1toN, q = hit->junctions_1toN; p != NULL; p = List_next(p), q = List_next(q)) {
      substring = (Substring_T) List_head(p);
      debug15(printf("substring: %u..%u\n",
		     Substring_alignstart_trim(substring) - 1U - hit->chroffset,
		     Substring_alignend_trim(substring) - hit->chroffset));
      if (Substring_overlap_point_trimmed_p(substring,common_genomicpos) == false) {
	*ilength_high += Substring_genomic_alignment_length(substring);
	if (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_high += Junction_nindels(junction);
	  }
	}

      } else {
	*ilength_high += ((Substring_alignstart_trim(substring) - 1) - common_genomicpos + 1);
	*ilength_low = (common_genomicpos - (Substring_alignend_trim(substring) /*+ 1*/) + 1);
	p = List_next(p);
	while (p != NULL) {
	  substring = (Substring_T) List_head(p);
	  *ilength_low += Substring_genomic_alignment_length(substring);
	  p = List_next(p);
	}
	while (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_low += Junction_nindels(junction);
	  }
	  q = List_next(q);
	}
	debug15(printf("Minus: Have ilength_low %d and ilength_high %d\n",*ilength_low,*ilength_high));
	return true;
      }
    }
  }

  return false;
}



/* Needed to compute overlap properly.  Based on pair_insert_length below, plus code for handling GMAP. */
static Univcoord_T
pair_common_genomicpos (Stage3end_T hit5, Stage3end_T hit3) {
  Univcoord_T common_genomicpos;
  Univcoord_T start5, end5, start3, end3;
  List_T p, q;
  Substring_T substring, substring5, substring3;

  if (hit5->plusp == true && hit3->plusp == true) {
    /* plus/plus */
    debug15(printf("Computing overlap using substrings plus/plus\n"));

    start5 = hit5->genomicstart + hit5->trim_querystart + start_amb_length(hit5);
    end5 = (hit5->genomicend - 1) - hit5->trim_queryend - end_amb_length(hit5);
    start3 = hit3->genomicstart + hit3->trim_querystart + start_amb_length(hit3);
    end3 = (hit3->genomicend - 1) - hit3->trim_queryend - end_amb_length(hit3);
    debug15(printf("hit5 endpoints are %u..%u.  hit3 endpoints are %u..%u\n",
		   start5-hit5->chroffset,end5-hit5->chroffset,start3-hit3->chroffset,end3-hit3->chroffset));

    if (end3 < start5) {
      /* Case 1 */
      return false;
    } else if (end5 < start3) {
      /* Case 6 */
      return false;
    } else if (start3 < start5) {
      if (end3 < end5) {
	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus/plus case 2a: start5 %u\n",start5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start5)) {
	    return start5;
	  }
	}

	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus/plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	for (p = hit5->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end3)) {
	    return end3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus/plus case 3\n"));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (end3 < end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus/plus case 4\n"));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus/plus case 5a\n"));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus/plus case 5b\n"));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus/plus general\n"));
    for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
      substring3 = (Substring_T) List_head(p);
      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring5 = (Substring_T) List_head(q);
	if ((common_genomicpos = Substring_overlap_segment_trimmed(substring5,substring3)) != 0) {
	  return common_genomicpos;
	}
      }
    }

    return 0;

  } else if (hit5->plusp == true && hit3->plusp == false) {
    /* plus/minus */
    debug15(printf("Computing overlap using substrings plus/minus\n"));
    return 0;

#if 0
    start5 = hit5->genomicstart + hit5->trim_querystart + start_amb_length(hit5);
    end5 = hit5->genomicend - hit5->trim_queryend - end_amb_length(hit5);
    start3 = hit3->genomicstart - hit3->trim_querystart - start_amb_length(hit3);
    end3 = hit3->genomicend + hit3->trim_queryend + end_amb_length(hit3);

    if (start3 < start5) {
      /* Case 1 */
      return 0;
    } else if (end5 < end3) {
      /* Case 6 */
      return 0;
    } else if (end3 < start5) {
      if (start3 < end5) {
	/* Case 2: Tails overlap.  Go from start5 to start3 */
	debug15(printf("plus case 2a: start5 %u\n",start5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	}

	/* Case 2: Tails overlap.  Go from start5 to start3 */
	debug15(printf("plus case 2b: start3 %u\n",start3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus case 3\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (start3 < end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus case 4\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus case 5a\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus case 5b\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus general: hit3->substring1\n"));
    if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring2 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring0 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring1)) != 0) {
      return common_genomicpos;
    }

    if (hit3->substring2 != NULL) {
      debug15(printf("plus general: hit3->substring2\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring2)) != 0) {
	return common_genomicpos;
      }
    }

    if (hit3->substring0 != NULL) {
      debug15(printf("plus general: hit3->substring0\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring0)) != 0) {
	return common_genomicpos;
      }
    }

    return 0U;
#endif

  } else if (hit5->plusp == false && hit3->plusp == true) {
    /* minus/plus */
    debug15(printf("Computing overlap using substrings minus/plus\n"));
    return 0;

#if 0
    start5 = hit5->genomicstart - hit5->trim_querystart - start_amb_length(hit5);
    end5 = hit5->genomicend + hit5->trim_queryend + end_amb_length(hit5);
    start3 = hit3->genomicstart + hit3->trim_querystart + start_amb_length(hit3);
    end3 = hit3->genomicend - hit3->trim_queryend - end_amb_length(hit3);

    if (end3 < end5) {
      /* Case 1 */
      return 0;
    } else if (start5 < start3) {
      /* Case 6 */
      return 0;
    } else if (start3 < end5) {
      if (end3 < start5) {
	/* Case 2: Tails overlap.  Go from end5 to end3 */
	debug15(printf("plus case 2a: end5 %u\n",end5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	}

	/* Case 2: Tails overlap.  Go from end5 to end3 */
	debug15(printf("plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus case 3\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (end3 < start5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus case 4\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus case 5a\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus case 5b\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus general: hit3->substring1\n"));
    if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring2 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring0 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring1)) != 0) {
      return common_genomicpos;
    }

    if (hit3->substring2 != NULL) {
      debug15(printf("plus general: hit3->substring2\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring2)) != 0) {
	return common_genomicpos;
      }
    }

    if (hit3->substring0 != NULL) {
      debug15(printf("plus general: hit3->substring0\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring0)) != 0) {
	return common_genomicpos;
      }
    }

    return 0;
#endif

  } else if (hit5->plusp == false && hit3->plusp == false) {
    /* minus/minus */
    debug15(printf("Computing overlap using substrings minus/minus\n"));

    start5 = (hit5->genomicstart - 1) - hit5->trim_querystart /*- start_amb_length(hit5)*/;
    end5 = hit5->genomicend + hit5->trim_queryend /*+ end_amb_length(hit5)*/;
    start3 = (hit3->genomicstart - 1) - hit3->trim_querystart /*- start_amb_length(hit3)*/;
    end3 = hit3->genomicend + hit3->trim_queryend /*+ end_amb_length(hit3)*/;
    debug15(printf("hit5 endpoints are %u..%u.  hit3 endpoints are %u..%u\n",
		   start5-hit5->chroffset,end5-hit5->chroffset,start3-hit3->chroffset,end3-hit3->chroffset));

    if (end3 > start5) {
      /* Case 1 */
      return 0;
    } else if (end5 > start3) {
      /* Case 6 */
      return 0;
    } else if (start3 > start5) {
      if (end3 > end5) {
	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("minus/minus case 2a: start5 %llu (%u)\n",start5,start5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start5)) {
	    return start5;
	  }
	}

	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	for (p = hit5->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end3)) {
	    return end3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("minus/minus case 3: end5 %u\n",end5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}

	/* Fall through to general algorithm */
      }

    } else {
      if (end3 > end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("minus/minus case 4: start3 %u\n",(Chrpos_T) (start3 - hit3->chroffset)));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("minus case 5a: start3 %u\n",start3 - hit3->chroffset));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("minus case 5b: end5 %u\n",end5 - hit5->chroffset));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("minus/minus general\n"));
    for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
      substring3 = (Substring_T) List_head(p);
      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring5 = (Substring_T) List_head(q);
	if ((common_genomicpos = Substring_overlap_segment_trimmed(substring5,substring3)) != 0) {
	  return common_genomicpos;
	}
      }
    }

    return 0;

  } else {
    abort();
    return 0;
  }
}


static bool
test_hardclips (Univcoord_T *common_genomicpos, int hardclip_low, Stage3end_T hit_low,
		int hardclip_high, Stage3end_T hit_high, Univcoord_T chroffset) {
  Substring_T low_substring, high_substring;
  int low_querypos, high_querypos;
  int low_querylength, high_querylength;
  bool plusp;

  low_querylength = hit_low->querylength;
  high_querylength = hit_high->querylength;

  debug15(printf("Entering test_hardclips with hardclip_low %d, hardclip_high %d\n",
		 hardclip_low,hardclip_high));
  debug15(printf("querylength_low %d, querylength_high %d\n",low_querylength,high_querylength));

  plusp = Stage3end_plusp(hit_low);

  if (plusp == true) {
    low_querypos = hardclip_low;
    high_querypos = high_querylength /*- 1*/ - hardclip_high;
    debug15(printf("Both substrings, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
      debug15(printf("Fails because low_querypos %d gives a NULL substring\n",low_querypos));
      return false;
    } else if (Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring) {
      debug15(printf("Fails because low_querypos %d - 1 gives substring %p\n",
		     low_querypos,Stage3end_substring_containing(hit_low,low_querypos-1)));
      return false;
    } else if (Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring) {
      debug15(printf("Fails because low_querypos %d + 1 gives substring %p\n",
		     low_querypos,Stage3end_substring_containing(hit_low,low_querypos+1)));
      return false;
    } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
      debug15(printf("Fails because high_querypos %d gives a NULL substring\n",high_querypos));
      return false;
    } else if (Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring) {
      debug15(printf("Fails because high_querypos %d - 1 gives substring %p\n",
		     high_querypos,Stage3end_substring_containing(hit_high,high_querypos-1)));
      return false;
    } else if (Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring) {
      debug15(printf("Fails because high_querypos %d + 1 gives substring %p\n",
		     high_querypos,Stage3end_substring_containing(hit_high,high_querypos+1)));
      return false;
    } else if (Substring_genomicstart(low_substring) + low_querypos - chroffset != Substring_genomicstart(high_substring) + high_querypos - chroffset) {
      debug15(printf("Fails because low chrpos %u != high chrpos %u\n",
		     Substring_genomicstart(low_substring) + low_querypos - chroffset,
		     Substring_genomicstart(high_substring) + high_querypos - chroffset));
      return false;
    } else {
      *common_genomicpos = Substring_genomicstart(low_substring) + low_querypos; /* Want univcoord */
      debug15(printf("Succeeds with common point %u\n",*common_genomicpos - chroffset));
      return true;
    }

  } else {
    low_querypos = low_querylength /*- 1*/ - hardclip_low;
    high_querypos = hardclip_high;
    debug15(printf("Both substrings, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
      debug15(printf("Fails because low_querypos %d gives a NULL substring\n",low_querypos));
      return false;
    } else if (Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring) {
      debug15(printf("Fails because low_querypos %d - 1 gives substring %p\n",
		     low_querypos,Stage3end_substring_containing(hit_low,low_querypos-1)));
      return false;
    } else if (Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring) {
      debug15(printf("Fails because low_querypos %d + 1 gives substring %p\n",
		     low_querypos,Stage3end_substring_containing(hit_low,low_querypos+1)));
      return false;
    } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
      debug15(printf("Fails because high_querypos %d gives a NULL substring\n",high_querypos));
      return false;
    } else if (Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring) {
      debug15(printf("Fails because high_querypos %d - 1 gives substring %p\n",
		     high_querypos,Stage3end_substring_containing(hit_high,high_querypos-1)));
      return false;
    } else if (Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring) {
      debug15(printf("Fails because high_querypos %d + 1 gives substring %p\n",
		     high_querypos,Stage3end_substring_containing(hit_high,high_querypos+1)));
      return false;
    }  else if ((Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset) {
      debug15(printf("Fails because low chrpos %u != high chrpos %u\n",
		     (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset,
		     (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset));
      return false;
    } else {
      *common_genomicpos = (Substring_genomicstart(low_substring) - 1) - low_querypos; /* Want univcoord */
      debug15(printf("Succeeds with common point %u\n",*common_genomicpos - chroffset));
      return true;
    }
  }
}



/* Replaces adjust_hardclips in samprint.c */
static Univcoord_T
adjust_hardclips_right (int *shift, int hardclip_low, Stage3end_T hit_low,
			int hardclip_high, Stage3end_T hit_high, Univcoord_T chroffset) {
  Substring_T low_substring, high_substring;
  int low_querypos, high_querypos;
  int low_querylength, high_querylength;
  Chrpos_T low_chrpos, high_chrpos;
  bool plusp;


  low_querylength = hit_low->querylength;
  high_querylength = hit_high->querylength;

  debug15(printf("Entering adjust_hardclips_right with hardclip_low %d, hardclip_high %d\n",
		 hardclip_low,hardclip_high));
  *shift = 1;			/* Making an initial move before each while loop */
  plusp = Stage3end_plusp(hit_low);

  if (plusp == true) {
    low_querypos = hardclip_low;
    high_querypos = high_querylength /*- 1*/ - hardclip_high;
    debug15(printf("Both substrings, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    low_querypos++;
    high_querypos++;
    debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
    while ((low_querypos + 1) < low_querylength && (high_querypos + 1) < high_querylength &&
	   ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	    Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	    Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	    (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	    Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	    Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	    Substring_genomicstart(low_substring) + low_querypos - chroffset != Substring_genomicstart(high_substring) + high_querypos - chroffset)) {
      (*shift) += 1;
      if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	low_querypos++;
      } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	high_querypos++;
      } else {
	low_chrpos = Substring_genomicstart(low_substring) + low_querypos - chroffset;
	high_chrpos = Substring_genomicstart(high_substring) + high_querypos - chroffset;
	if (low_chrpos < high_chrpos) {
	  debug15(printf("low_chrpos %u < high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	  low_querypos++;
	} else if (high_chrpos < low_chrpos) {
	  debug15(printf("high_chrpos %u < low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	  high_querypos++;
	} else {
	  low_querypos++;
	  high_querypos++;
	}
      }
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
    }

    if ((low_querypos + 1) >= low_querylength ||
	(high_querypos + 1) >= high_querylength ||
	(low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
      *shift = 0;
      return 0;
    } else {
      debug15(printf("Returning %u + %d\n",Substring_genomicstart(low_substring) - chroffset,
		     low_querypos));
      assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
      assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
      assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
      assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
      assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
      assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
      return Substring_genomicstart(low_substring) + low_querypos; /* Want univcoord */
    }

  } else {
    low_querypos = low_querylength /*- 1*/ - hardclip_low;
    high_querypos = hardclip_high;
    debug15(printf("Both substrings, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    low_querypos--;
    high_querypos--;
    debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
    while ((low_querypos - 1) >= 0 && (high_querypos - 1) >= 0 &&
	   ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	    Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	    Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	    (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	    Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	    Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	    (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset)) {
      (*shift) += 1;
      if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	low_querypos--;
      } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	high_querypos--;
      } else {
	low_chrpos = (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset;
	high_chrpos = (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset;
	if (low_chrpos < high_chrpos) {
	  debug15(printf("low_chrpos %u < high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	  low_querypos--;
	} else if (high_chrpos < low_chrpos) {
	  debug15(printf("high_chrpos %u < low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	  high_querypos--;
	} else {
	  low_querypos--;
	  high_querypos--;
	}
      }
      debug15(printf("right shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
    }

    if ((low_querypos - 1) < 0 ||
	(high_querypos - 1) < 0 ||
	(low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
      *shift = 0;
      return 0;
    } else {
      debug15(printf("Returning %u - %d\n",Substring_genomicstart(low_substring) - chroffset,
		     low_querypos));
      assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
      assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
      assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
      assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
      assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
      assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
      return (Substring_genomicstart(low_substring) - 1) - low_querypos; /* Want univcoord */
    }
  }
}


/* Replaces adjust_hardclips in samprint.c */
static Univcoord_T
adjust_hardclips_left (int *shift, int hardclip_low, Stage3end_T hit_low,
		       int hardclip_high, Stage3end_T hit_high, Univcoord_T chroffset) {
  Substring_T low_substring, high_substring;
  int low_querypos, high_querypos;
  int low_querylength, high_querylength;
  Chrpos_T low_chrpos, high_chrpos;
  bool plusp;


  low_querylength = hit_low->querylength;
  high_querylength = hit_high->querylength;

  debug15(printf("Entering adjust_hardclips_left with hardclip_low %d, hardclip_high %d\n",
		 hardclip_low,hardclip_high));
  *shift = 1;			/* Making an initial move before each while loop */
  plusp = Stage3end_plusp(hit_low);

  if (plusp == true) {
    low_querypos = hardclip_low;
    high_querypos = high_querylength /*- 1*/ - hardclip_high;
    debug15(printf("Both substrings, plus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    low_querypos--;
    high_querypos--;
    debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
    while ((low_querypos - 1) >= 0 && (high_querypos - 1) >= 0 &&
	   ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	    Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	    Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	    (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	    Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	    Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	    Substring_genomicstart(low_substring) + low_querypos - chroffset != Substring_genomicstart(high_substring) + high_querypos - chroffset)) {
      (*shift) += 1;
      if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	low_querypos--;
      } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	high_querypos--;
      } else {
	low_chrpos = Substring_genomicstart(low_substring) + low_querypos - chroffset;
	high_chrpos = Substring_genomicstart(high_substring) + high_querypos - chroffset;
	if (low_chrpos > high_chrpos) {
	  debug15(printf("low_chrpos %u > high_chrpos %u, so decreasing low_querypos\n",low_chrpos,high_chrpos));
	  low_querypos--;
	} else if (high_chrpos > low_chrpos) {
	  debug15(printf("high_chrpos %u > low_chrpos %u, so decreasing high_querypos\n",high_chrpos,low_chrpos));
	  high_querypos--;
	} else {
	  low_querypos--;
	  high_querypos--;
	}
      }
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
    }

    if ((low_querypos - 1) < 0 || (high_querypos - 1) < 0 ||
	(low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
      *shift = 0;
      return 0;
    } else {
      debug15(printf("Returning %u + %d\n",Substring_genomicstart(low_substring) - chroffset,
		     low_querypos));
      assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
      assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
      assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
      assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
      assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
      assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
      return Substring_genomicstart(low_substring) + low_querypos; /* Want univcoord */
    }

  } else {
    low_querypos = low_querylength /*- 1*/ - hardclip_low;
    high_querypos = hardclip_high;
    debug15(printf("Both substrings, minus.  low_querypos %d, high_querypos %d\n",low_querypos,high_querypos));

    low_querypos++;
    high_querypos++;
    debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
    while ((low_querypos + 1) < low_querylength && (high_querypos + 1) < high_querylength &&
	   ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	    Stage3end_substring_containing(hit_low,low_querypos-1) != low_substring ||
	    Stage3end_substring_containing(hit_low,low_querypos+1) != low_substring ||
	    (high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL ||
	    Stage3end_substring_containing(hit_high,high_querypos-1) != high_substring ||
	    Stage3end_substring_containing(hit_high,high_querypos+1) != high_substring ||
	    (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset != (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset)) {
      (*shift) += 1;
      if ((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL) {
	low_querypos++;
      } else if ((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) == NULL) {
	high_querypos++;
      } else {
	low_chrpos = (Substring_genomicstart(low_substring) - 1) - low_querypos - chroffset;
	high_chrpos = (Substring_genomicstart(high_substring) - 1) - high_querypos - chroffset;
	if (low_chrpos > high_chrpos) {
	  debug15(printf("low_chrpos %u > high_chrpos %u, so advancing low_querypos\n",low_chrpos,high_chrpos));
	  low_querypos++;
	} else if (high_chrpos > low_chrpos) {
	  debug15(printf("high_chrpos %u > low_chrpos %u, so advancing high_querypos\n",high_chrpos,low_chrpos));
	  high_querypos++;
	} else {
	  low_querypos++;
	  high_querypos++;
	}
      }
      debug15(printf("left shift %d: Advancing to low_querypos %d and high_querypos %d\n",*shift,low_querypos,high_querypos));
    }

    if ((low_querypos + 1) >= low_querylength || (high_querypos + 1) >= high_querylength ||
	(low_substring = Stage3end_substring_containing(hit_low,low_querypos)) == NULL ||
	Stage3end_substring_containing(hit_high,high_querypos) == NULL) {
      *shift = 0;
      return 0;
    } else {
      debug15(printf("Returning %u - %d\n",Substring_genomicstart(low_substring) - chroffset,
		     low_querypos));
      assert((low_substring = Stage3end_substring_containing(hit_low,low_querypos)) != NULL);
      assert((high_substring = Stage3end_substring_containing(hit_high,high_querypos)) != NULL);
      assert(Stage3end_substring_containing(hit_low,low_querypos-1) == low_substring);
      assert(Stage3end_substring_containing(hit_low,low_querypos+1) == low_substring);
      assert(Stage3end_substring_containing(hit_high,high_querypos-1) == high_substring);
      assert(Stage3end_substring_containing(hit_high,high_querypos+1) == high_substring);
      return (Substring_genomicstart(low_substring) - 1) - low_querypos; /* Want univcoord */
    }
  }
}



/* Note: Do not alter this->insertlength, which is used for SAM
   output.  The insertlength computed here is used only for performing
   --clip-overlap or --merge-overlap */
int
Stage3pair_overlap (int *hardclip5_low, int *hardclip5_high, int *hardclip3_low, int *hardclip3_high, Stage3pair_T this) {
  Stage3end_T hit5, hit3;
  int clipdir;
  int ilength53, ilength35, ilength5_low, ilength5_high, ilength3_low, ilength3_high;
  int common_shift, common_left, common_right;
  Univcoord_T common_genomicpos, common_genomicpos_right, common_genomicpos_left;
  int shift_right, shift_left;
#ifdef DEBUG15
  int overlap;
#endif


  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;

  hit5 = this->hit5;
  hit3 = this->hit3;

  debug15(printf("Entered Stage3pair_overlap with hittype %s and %s\n",
		 hittype_string(hit5->hittype),hittype_string(hit3->hittype)));
  if (hit5->hittype == SAMECHR_SPLICE || hit5->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (hit3->hittype == SAMECHR_SPLICE || hit3->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (hit5->plusp != hit3->plusp) {
    debug15(printf("The two ends are not on the same strand, so returning 0\n"));
    return 0;
  } else {
    debug15(printf("hit5 trim_querystart %d + amb_start %d, trim_queryend %d + amb_end %d, hit3 trim_querystart %d + amb_start %d, trim_queryend %d + amb_end %d\n",
		   hit5->trim_querystart,start_amb_length(hit5),hit5->trim_queryend,end_amb_length(hit5),
		   hit3->trim_querystart,start_amb_length(hit3),hit3->trim_queryend,end_amb_length(hit3)));
    if (hit5->plusp == true) {
      /* plus */
#if 0
      hit5_trimmed_length = hit5->querylength - hit5->trim_querystart - hit5->trim_queryend - start_amb_length(hit5) - end_amb_length(hit5);
      hit3_trimmed_length = hit3->querylength - hit3->trim_querystart - hit3->trim_queryend - start_amb_length(hit3) - end_amb_length(hit3);
      totallength = hit5_trimmed_length + hit3_trimmed_length;
      debug15(printf("totallength = %d, hit5 trimmed length = %d, hit3 trimmed length = %d\n",
		     totallength,hit5_trimmed_length,hit3_trimmed_length));
      debug15(printf("original insertlength: %d, trim+amb5: %d..%d, trim+amb3: %d..%d\n",
		     this->insertlength,hit5->trim_querystart + start_amb_length(hit5),
		     hit5->trim_queryend + end_amb_length(hit5),hit3->trim_querystart + start_amb_length(hit3),
		     hit3->trim_queryend + end_amb_length(hit3)));
#endif

      if ((common_genomicpos = pair_common_genomicpos(hit5,hit3)) == 0) {
	debug15(printf("Cannot determine a common point, so returning 0\n"));
	return 0;

      } else if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos) == false ||
		 find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos) == false) {
	debug15(printf("Cannot determine ilengths, so returning 0\n"));
	return 0;

      } else {
	debug15(printf("Inclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	debug15(printf("ilength53 is %d, ilength 35 is %d\n",ilength5_low + ilength3_high - 1,ilength3_low + ilength5_high - 1));

	common_left = (ilength5_low < ilength3_low) ? ilength5_low : ilength3_low;
	common_right = (ilength5_high < ilength3_high) ? ilength5_high : ilength3_high;
	if (common_right > common_left) {
	  common_shift = common_right/2 - (common_left - 1)/2;
	  debug15(printf("Common shift is %d = common_right %d/2 - (common_left %d - 1)/2\n",
			 common_shift,common_right,common_left));
	  assert(ilength5_low > 0);
	  assert(ilength3_low > 0);
	  ilength5_low -= 1;
	  ilength3_low -= 1;
	} else {
	  common_shift = (common_right - 1)/2 - common_left/2;
	  debug15(printf("Common shift is %d = (common_right %d - 1)/2 - common_left %d/2\n",
			 common_shift,common_right,common_left));
	  assert(ilength5_high > 0);
	  assert(ilength3_high > 0);
	  ilength5_high -= 1;
	  ilength3_high -= 1;
	}
	debug15(printf("Exclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));


	if ((ilength53 = ilength5_low + ilength3_high) >= (ilength35 = ilength3_low + ilength5_high)) {
	  /* Use >=, not >, so we favor clipping heads over clipping tails in case of a tie */
	  debug15(printf("plus, ilength53 is longer.  Clipping heads.\n"));
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 common_left+common_right-1,common_left,common_right));
	  clipdir = +1;

	  /* Want to clip 5 high and 3 low */
	  *hardclip5_high = ilength5_high - common_shift;
	  *hardclip3_low = ilength3_low + common_shift;
	  debug15(printf("Overlap clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_high += hit5->trim_queryend /*+ end_amb_length(hit5)*/;
	  *hardclip3_low += hit3->trim_querystart /*+ start_amb_length(hit3)*/;
	  debug15(printf("Ambig clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_high = ilength5_high /*- common_shift*/;
		*hardclip3_low = ilength3_low /*+ common_shift*/;
		*hardclip5_high += hit5->trim_queryend /*+ end_amb_length(hit5)*/;
		*hardclip3_low += hit3->trim_querystart /*+ start_amb_length(hit3)*/;
		if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - hit3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength3_low > ilength5_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength3_low > 0);
	      ilength3_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength5_high > 0);
	      ilength5_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_high = ilength5_high /*- common_shift*/;
	    *hardclip3_low = ilength3_low /*+ common_shift*/;
	    debug15(printf("Initial computation of clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	    *hardclip5_high += hit5->trim_queryend /*+ end_amb_length(hit5)*/;
	    *hardclip3_low += hit3->trim_querystart /*+ start_amb_length(hit3)*/;
	    debug15(printf("Recomputed clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

#if 0
	  if (*hardclip5_high < 0) {
	    *hardclip5_high = 0;
	  }
	  if (*hardclip3_low < 0) {
	    *hardclip3_low = 0;
	  }
	  debug15(printf("Positive clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif

	} else {
	  debug15(printf("plus, ilength35 is longer.  Clipping tails.\n"));
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 common_left+common_right-1,common_left,common_right));
	  clipdir = -1;

	  /* Want to clip 5 low and 3 high */
	  *hardclip5_low = ilength5_low + common_shift;
	  *hardclip3_high = ilength3_high - common_shift;
	  debug15(printf("Overlap clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_low += hit5->trim_querystart /*+ start_amb_length(hit5)*/;
	  *hardclip3_high += hit3->trim_queryend /*+ end_amb_length(hit3)*/;
	  debug15(printf("Ambig clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_low = ilength5_low /*+ common_shift*/;
		*hardclip3_high = ilength3_high /*- common_shift*/;
		*hardclip5_low += hit5->trim_querystart /*+ start_amb_length(hit5)*/;
		*hardclip3_high += hit3->trim_queryend /*+ end_amb_length(hit3)*/;
		if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - hit3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength5_low > ilength3_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength5_low > 0);
	      ilength5_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength3_high > 0);
	      ilength3_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_low = ilength5_low /*+ common_shift*/;
	    *hardclip3_high = ilength3_high /*- common_shift*/;
	    debug15(printf("Initial computation of clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	    *hardclip5_low += hit5->trim_querystart /*+ start_amb_length(hit5)*/;
	    *hardclip3_high += hit3->trim_queryend /*+ end_amb_length(hit3)*/;
	    debug15(printf("Recomputed clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

#if 0
	  if (*hardclip5_low < 0) {
	    *hardclip5_low = 0;
	  }
	  if (*hardclip3_high < 0) {
	    *hardclip3_high = 0;
	  }
	  debug15(printf("Positive clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif
	}

	debug15(printf("returning clipdir %d\n",clipdir));
	return clipdir;
      }

    } else {
      /* minus */
#if 0
      hit5_trimmed_length = hit5->querylength - hit5->trim_querystart - hit5->trim_queryend - start_amb_length(hit5) - end_amb_length(hit5);
      hit3_trimmed_length = hit3->querylength - hit3->trim_querystart - hit3->trim_queryend - start_amb_length(hit3) - end_amb_length(hit3);
      totallength = hit5_trimmed_length + hit3_trimmed_length;
      debug15(printf("totallength = %d, hit5 trimmed length = %d, hit3 trimmed length = %d\n",
		     totallength,hit5_trimmed_length,hit3_trimmed_length));
      debug15(printf("original insertlength: %d, trim+amb5: %d..%d, trim+amb3: %d..%d\n",
		     this->insertlength,hit5->trim_querystart + start_amb_length(hit5),
		     hit5->trim_queryend + hit5->end_amb_length,hit3->trim_querystart + start_amb_length(hit3),
		     hit3->trim_queryend + hit3->end_amb_length));
#endif

      if ((common_genomicpos = pair_common_genomicpos(hit5,hit3)) == 0) {
	debug15(printf("Cannot determine a common point, so returning 0\n"));
	return 0;

      } else if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos) == false ||
		 find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos) == false) {
	debug15(printf("Cannot determine ilengths, so returning 0\n"));
	return 0;

      } else {
	debug15(printf("Inclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	debug15(printf("ilength53lh is %d, ilength35lh is %d\n",ilength5_low + ilength3_high - 1,ilength3_low + ilength5_high - 1));

	common_left = (ilength5_low < ilength3_low) ? ilength5_low : ilength3_low;
	common_right = (ilength5_high < ilength3_high) ? ilength5_high : ilength3_high;
	if (common_right > common_left) {
	  common_shift = common_right/2 - (common_left - 1)/2;
	  debug15(printf("Common shift is %d = common_right %d/2 - (common_left %d - 1)/2\n",
			 common_shift,common_right,common_left));
	  assert(ilength5_low > 0);
	  assert(ilength3_low > 0);
	  ilength5_low -= 1;
	  ilength3_low -= 1;
	} else {
	  common_shift = (common_right - 1)/2 - common_left/2;
	  debug15(printf("Common shift is %d = (common_right %d - 1)/2 - common_left %d/2\n",
			 common_shift,common_right,common_left));
	  assert(ilength5_high > 0);
	  assert(ilength3_high > 0);
	  ilength5_high -= 1;
	  ilength3_high -= 1;
	}
	debug15(printf("Exclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	if ((ilength53 = ilength5_low + ilength3_high) > (ilength35 = ilength3_low + ilength5_high)) {
	  /* Use >, not >=, so we favor clipping heads over clipping tails in case of a tie */
	  debug15(printf("minus, ilength53 is longer.  Clipping tails.\n"));
	  debug15(overlap = common_left + common_right - 1);
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 overlap,common_left,common_right));
	  clipdir = +1;


	  /* Want to clip 5 high and 3 low */
	  *hardclip5_high = ilength5_high - common_shift;
	  *hardclip3_low = ilength3_low + common_shift;
	  debug15(printf("Overlap clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_high += hit5->trim_querystart /*+ start_amb_length(hit5)*/;
	  *hardclip3_low += hit3->trim_queryend /*+ end_amb_length(hit3)*/;
	  debug15(printf("Ambig clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_high = ilength5_high /*- common_shift*/;
		*hardclip3_low = ilength3_low /*+ common_shift*/;
		*hardclip5_high += hit5->trim_querystart /*+ start_amb_length(hit5)*/;
		*hardclip3_low += hit3->trim_queryend /*+ end_amb_length(hit3)*/;
		if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - hit3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength3_low > ilength5_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength3_low > 0);
	      ilength3_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength5_high > 0);
	      ilength5_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_high = ilength5_high /*- common_shift*/;
	    *hardclip3_low = ilength3_low /*+ common_shift*/;
	    debug15(printf("Initial computation of clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	    *hardclip5_high += hit5->trim_querystart /*+ start_amb_length(hit5)*/;
	    *hardclip3_low += hit3->trim_queryend /*+ end_amb_length(hit3)*/;
	    debug15(printf("Recomputed clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

#if 0
	  if (*hardclip5_high < 0) {
	    *hardclip5_high = 0;
	  }
	  if (*hardclip3_low < 0) {
	    *hardclip3_low = 0;
	  }
	  debug15(printf("Positive clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif

	} else {
	  debug15(printf("minus, ilength35 is longer.  Clipping heads.\n"));
	  debug15(overlap = common_left + common_right - 1);
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 overlap,common_left,common_right));
	  clipdir = -1;

	  /* Want to clip 5 low and 3 high */
	  *hardclip5_low = ilength5_low + common_shift;
	  *hardclip3_high = ilength3_high - common_shift;
	  debug15(printf("Overlap clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_low += hit5->trim_queryend /*+ end_amb_length(hit5)*/;
	  *hardclip3_high += hit3->trim_querystart /*+ start_amb_length(hit3)*/;
	  debug15(printf("Ambig clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip5_low,hit5,*hardclip3_high,hit3,hit3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_low = ilength5_low /*+ common_shift*/;
		*hardclip3_high = ilength3_high /*- common_shift*/;
		*hardclip5_low += hit5->trim_queryend /*+ end_amb_length(hit5)*/;
		*hardclip3_high += hit3->trim_querystart /*+ start_amb_length(hit3)*/;
		if (test_hardclips(&common_genomicpos,*hardclip3_low,hit3,*hardclip5_high,hit5,hit3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - hit3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,hit5,common_genomicpos) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,hit3,common_genomicpos) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength5_low > ilength3_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength5_low > 0);
	      ilength5_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength3_high > 0);
	      ilength3_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_low = ilength5_low /*+ common_shift*/;
	    *hardclip3_high = ilength3_high /*- common_shift*/;
	    debug15(printf("Initial computation of clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	    *hardclip5_low += hit5->trim_queryend /*+ end_amb_length(hit5)*/;
	    *hardclip3_high += hit3->trim_querystart /*+ start_amb_length(hit3)*/;
	    debug15(printf("Recomputed clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

#if 0
	  if (*hardclip5_low < 0) {
	    *hardclip5_low = 0;
	  }
	  if (*hardclip3_high < 0) {
	    *hardclip3_high = 0;
	  }
	  debug15(printf("Positive clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif
	}
      }

      debug15(printf("returning clipdir %d\n",clipdir));
      return clipdir;
    }
  }
}


void
Stage3pair_free (Stage3pair_T *old) {
  debug0(printf("Freeing pair %p with hits %p and %p\n",*old,(*old)->hit5,(*old)->hit3));
  assert((*old)->hit3 != NULL);
  debug0(printf("Freeing end3 at %p\n",(*old)->hit3));
  Stage3end_free(&(*old)->hit3);

  assert((*old)->hit5 != NULL);
  debug0(printf("Freeing end5 at %p\n",(*old)->hit5));
  Stage3end_free(&(*old)->hit5);

  if ((*old)->transcripts5 != NULL) {
    Transcript_gc(&(*old)->transcripts5);
  }
  if ((*old)->transcripts3 != NULL) {
    Transcript_gc(&(*old)->transcripts3);
  }

  FREE_OUT(*old);
  return;
}



#if 0
static long int
Stage3pair_tally (Stage3pair_T this) {

  if (tally_iit == NULL) {
    return 0L;
  } else if (this->tally >= 0) {
    return this->tally;
  } else {
    this->tally = Stage3end_compute_tally(this->hit5) + Stage3end_compute_tally(this->hit3);
    return this->tally;
  }
}
#endif


static char complCode[128] = COMPLEMENT_LC;

#if 0
static char *
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC_OUT(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}
#endif

static char *
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return sequence;
}

char *
Stage3end_substrings_genomic_sequence (int *seqlength, T this, Genome_T genome) {
  char *gbuffer;
  List_T p, q;
  Substring_T substring;
  Junction_T junction;
  int querypos, querystart, queryend, querylength, substring_length;

  *seqlength = 0;
  for (p = this->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
#ifdef NO_SOFT_CLIPS
    querystart = Substring_querystart_orig(substring);
    queryend = Substring_queryend_orig(substring);
#else
    querystart = Substring_querystart(substring);
    queryend = Substring_queryend(substring);
#endif
    *seqlength += queryend - querystart;
  }
  for (p = this->junctions_1toN; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == DEL_JUNCTION) {
      *seqlength += Junction_nindels(junction);
    }
  }

  gbuffer = (char *) MALLOC((*seqlength+1) * sizeof(char));
  if (this->plusp == true) {
    /* Build from querystart to queryend, so we don't wipe out sequence with terminating \0 character */
    querypos = 0;
    for (p = this->substrings_1toN, q = this->junctions_1toN; p != NULL; p = List_next(p), q = List_next(q)) {
      substring = (Substring_T) List_head(p);
#ifdef NO_SOFT_CLIPS
      querystart = Substring_querystart_orig(substring);
      queryend = Substring_queryend_orig(substring);
#else
      querystart = Substring_querystart(substring);
      queryend = Substring_queryend(substring);
#endif
      substring_length = queryend - querystart;
      Genome_fill_buffer_simple(genome,Substring_left(substring) + querystart,
				substring_length,&(gbuffer[querypos]));
      querypos += substring_length;

      if (q != NULL) {
	junction = (Junction_T) List_head(q);
	if (Junction_type(junction) == DEL_JUNCTION) {
	  substring_length = Junction_nindels(junction);
	  Genome_fill_buffer_simple(genome,Junction_deletionpos(junction),
				    substring_length,&(gbuffer[querypos]));
	  querypos += substring_length;
	}
      }
    }

    return gbuffer;

  } else {
    /* Build from queryend to querystart, so we don't wipe out sequence with terminating \0 character */
    querypos = 0;
    querylength = this->querylength;
    for (p = this->substrings_Nto1, q = this->junctions_Nto1; p != NULL; p = List_next(p), q = List_next(q)) {
      substring = (Substring_T) List_head(p);
#ifdef NO_SOFT_CLIPS
      querystart = Substring_querystart_orig(substring);
      queryend = Substring_queryend_orig(substring);
#else
      querystart = Substring_querystart(substring);
      queryend = Substring_queryend(substring);
#endif
      substring_length = queryend - querystart;
      Genome_fill_buffer_simple(genome,Substring_left(substring) + (querylength - queryend),
				substring_length,&(gbuffer[querypos]));
      querypos += substring_length;

      if (q != NULL) {
	junction = (Junction_T) List_head(q);
	if (Junction_type(junction) == DEL_JUNCTION) {
	  substring_length = Junction_nindels(junction);
	  Genome_fill_buffer_simple(genome,Junction_deletionpos(junction),
				    substring_length,&(gbuffer[querypos]));
	  querypos += substring_length;
	}
      }
    }

    return make_complement_inplace(gbuffer,*seqlength);
  }
}


const Except_T Copy_Substring = { "Substring invalid during copy" };

static T
Stage3end_copy (T old, Listpool_T listpool) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  List_T p;
  Substring_T old_substring, new_substring;
  Junction_T old_junction, new_junction;

  debug0(printf("Copying Stage3end %p -> %p of type %s\n",
		old,new,hittype_string(old->hittype)));

  new->hittype = old->hittype;
  new->method = old->method;
  new->level = old->level;
  new->genestrand = old->genestrand;

  new->distant_splice_p = old->distant_splice_p;
  new->chrnum = old->chrnum;
  new->effective_chrnum = old->effective_chrnum;
  new->other_chrnum = old->other_chrnum;
  new->chroffset = old->chroffset;
  new->chrhigh = old->chrhigh;
  new->chrlength = old->chrlength;

  new->querylength = old->querylength;
  new->querylength_adj = old->querylength_adj;

  new->genomicstart = old->genomicstart;
  new->genomicend = old->genomicend;
  new->plusp = old->plusp;

  new->low = old->low;
  new->high = old->high;
  new->genomiclength = old->genomiclength;
  new->guided_insertlength = old->guided_insertlength;

  new->mapq_loglik = old->mapq_loglik;
  new->mapq_score = old->mapq_score;
  new->absmq_score = old->absmq_score;

  new->nsegments = old->nsegments;

  new->score_overall = old->score_overall;
  new->nmatches_to_trims = old->nmatches_to_trims;

  new->score_within_trims = old->score_within_trims;
  new->nmatches_plus_spliced_trims = old->nmatches_plus_spliced_trims;
  /* new->nmatches_amb = old->nmatches_amb; */
  new->splice_score = old->splice_score;

  new->trim_querystart = old->trim_querystart;
  new->trim_queryend = old->trim_queryend;
  new->trim_querystart_splicep = old->trim_querystart_splicep;
  new->trim_queryend_splicep = old->trim_queryend_splicep;

  /* new->penalties = old->penalties; */
  new->score_eventrim = old->score_eventrim;

  new->gene_overlap = old->gene_overlap;
  new->tally = old->tally;

  new->nmismatches_bothdiff = old->nmismatches_bothdiff;
  new->nmismatches_refdiff = old->nmismatches_refdiff;

  new->nindels = old->nindels;

  new->distance = old->distance;

  new->sensedir_for_concordance = old->sensedir_for_concordance;
  new->sensedir = old->sensedir;

  new->nsplices = old->nsplices;

  new->substrings_1toN = (List_T) NULL;
  new->substrings_Nto1 = (List_T) NULL;
  new->substrings_LtoH = (List_T) NULL;
  new->substrings_HtoL = (List_T) NULL;

  new->junctions_1toN = (List_T) NULL;
  new->junctions_Nto1 = (List_T) NULL;
  new->junctions_LtoH = (List_T) NULL;
  /* new->junctions_HtoL = (List_T) NULL; */

  new->transcripts = Transcript_copy_list(old->transcripts);
  new->transcripts_other = Transcript_copy_list(old->transcripts_other);

  for (p = old->substrings_1toN; p != NULL; p = List_next(p)) {
    old_substring = (Substring_T) List_head(p);
    new_substring = Substring_copy(old_substring);
    new->substrings_1toN = Listpool_push(new->substrings_1toN,listpool,(void *) new_substring);
  }

  for (p = old->junctions_1toN; p != NULL; p = List_next(p)) {
    old_junction = (Junction_T) List_head(p);
    new_junction = Junction_copy(old_junction);
    new->junctions_1toN = Listpool_push(new->junctions_1toN,listpool,(void *) new_junction);
  }

  new->substrings_Nto1 = Listpool_copy(new->substrings_1toN,listpool); /* Before reversal of 1toN */
  new->junctions_Nto1 = Listpool_copy(new->junctions_1toN,listpool);   /* Before reversal of 1toN */

  /* Reversals to handle builds of 1toN */
  new->substrings_1toN = List_reverse(new->substrings_1toN);
  new->junctions_1toN = List_reverse(new->junctions_1toN);

  if (old->chrnum == 0) {
    /* Translocation */
    if (old->sensedir == SENSE_FORWARD) {
      new->substrings_LtoH = Listpool_copy(new->substrings_1toN,listpool);
      new->junctions_LtoH = Listpool_copy(new->junctions_1toN,listpool);
      new->substrings_HtoL = Listpool_copy(new->substrings_Nto1,listpool);
      /* new->junctions_HtoL = Listpool_copy(new->junctions_Nto1,listpool); */
    } else if (old->sensedir == SENSE_ANTI) {
      new->substrings_LtoH = Listpool_copy(new->substrings_Nto1,listpool);
      new->junctions_LtoH = Listpool_copy(new->junctions_Nto1,listpool);
      new->substrings_HtoL = Listpool_copy(new->substrings_1toN,listpool);
      /* new->junctions_HtoL = Listpool_copy(new->junctions_1toN,listpool); */
    } else {
      /* SENSE_NULL: Treat as SENSE_FORWARD */
      new->substrings_LtoH = Listpool_copy(new->substrings_1toN,listpool);
      new->junctions_LtoH = Listpool_copy(new->junctions_1toN,listpool);
      new->substrings_HtoL = Listpool_copy(new->substrings_Nto1,listpool);
      /* new->junctions_HtoL = Listpool_copy(new->junctions_Nto1,listpool); */
    }

  } else {
    if (old->plusp == true) {
      new->substrings_LtoH = Listpool_copy(new->substrings_1toN,listpool);
      new->junctions_LtoH = Listpool_copy(new->junctions_1toN,listpool);
      new->substrings_HtoL = Listpool_copy(new->substrings_Nto1,listpool);
      /* new->junctions_HtoL = Listpool_copy(new->junctions_Nto1,listpool); */
    } else {
      new->substrings_LtoH = Listpool_copy(new->substrings_Nto1,listpool);
      new->junctions_LtoH = Listpool_copy(new->junctions_Nto1,listpool);
      new->substrings_HtoL = Listpool_copy(new->substrings_1toN,listpool);
      /* new->junctions_HtoL = Listpool_copy(new->junctions_1toN,listpool); */
    }
  }
  /* Actually, the assertion is excluded only for the JOIN hittype */
  assert(new->hittype == SPLICE || Substring_querystart(List_head(new->substrings_1toN)) <= Substring_querystart(List_head(new->substrings_Nto1)));

  new->paired_usedp = old->paired_usedp;
  new->concordantp = old->concordantp;

  new->query_splicepos = old->query_splicepos;
  new->circularalias = old->circularalias;
  new->circularpos = old->circularpos;
  debug12(printf("Copying circularpos of %d from hit %p to hit %p\n",new->circularpos,old,new));

  new->altlocp = old->altlocp;

  return new;
}


static int
compute_circularpos (int *circularalias, T hit) {
  int circularpos;
  List_T p;
  Substring_T substring;


  debug12(printf("Computing circularpos on hit at %u..%u, plusp %d, with trim left %d and trim right %d\n",
		 hit->genomicstart - hit->chroffset,hit->genomicend - hit->chroffset,
		 hit->plusp,hit->trim_querystart,hit->trim_queryend));
  if (circularp[hit->chrnum] == false) {
    debug12(printf("Chromosome #%d is not circular\n",hit->chrnum));
    /* This also handles hit->chrnum == 0, where translocation cannot be circular */
    *circularalias = 0;
    return -1;

  } else if (Stage3end_chrpos_low_trim(hit) >= hit->chrlength) {
    /* All of read after trimming is in high part.  Previously
       checked hit->high against hit->chrhigh, for circularalias of
       +2, but that should be fixed now */

    debug12(printf("All of read after trimming %u..%u is in high part\n",
		   Stage3end_chrpos_low_trim(hit),Stage3end_chrpos_high_trim(hit)));
    *circularalias = +1;		/* All of read is in second copy */
    debug12(printf("For hit %p, pair circularpos is -1, circularalias is %d\n",hit,*circularalias));
    return -1;

  } else if (Stage3end_chrpos_high_trim(hit) < hit->chrlength) {
    /* All of read after trimming is in low part.  Previously
       checked hit->low against hit->chroffset for circularalias of
       -2, but that should be fixed now */

    debug12(printf("All of read after trimming %u..%u is in low part\n",
		   Stage3end_chrpos_low_trim(hit),Stage3end_chrpos_high_trim(hit)));
    *circularalias = -1;		/* All of read is in first copy */
    debug12(printf("For hit %p, pair circularpos is -1, circularalias is %d\n",hit,*circularalias));
    return -1;

  } else {
    *circularalias = 0;	/* Straddling middle */
    for (p = hit->substrings_LtoH; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      if ((circularpos = Substring_circularpos(substring)) > 0) {
	debug12(printf("For hit %p, returning circularpos %d from substring (plus)\n",hit,circularpos));
	return circularpos;
      }
    }      
    debug12(printf("For hit %p, pair circularpos is -1, circularalias is %d\n",hit,*circularalias));
    return -1;
  }
}


/* Modified from Stage3end_new_precomputed for a single substring */
T
Stage3end_new_terminal (int *found_score_overall, int *found_score_within_trims,
			Substring_T substring_in, int querylength,
			bool gplusp, int genestrand, int sensedir, Listpool_T listpool,
			Method_T method, int level) {
  T new;

  Substring_T substring;
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;

  Univcoord_T genomicstart, genomicend;
  List_T substrings;
  List_T p;
  int adj = 0;


  substring = Substring_copy(substring_in);
  chrnum = Substring_chrnum(substring);
  chroffset = Substring_chroffset(substring);
  chrhigh = Substring_chrhigh(substring);
  chrlength = Substring_chrlength(substring);

  debug0(printf("Entered Stage3end_new_terminal, method %s, with chrnum %d, query %d..%d\n",
		Method_string(method),chrnum,Substring_querystart(substring),Substring_queryend(substring)));

  new = (T) MALLOC_OUT(sizeof(*new));
  new->hittype = SUBSTRINGS;
  new->method = method;
  new->level = level;

  new->querylength = querylength;
  new->querylength_adj = querylength + adj;

  /* Caller must not free these lists */
  new->transcripts = (List_T) NULL;
  new->transcripts_other = (List_T) NULL;

  /* Unlike Stage3end_new_substrings, where substrings and junctions
     are in opposite orders, substrings and junctions here are in the
     same order. */

  substrings = Listpool_push(NULL,listpool,(void *) substring);
  new->substrings_1toN = substrings;
  new->substrings_Nto1 = Listpool_copy(substrings,listpool);
  new->substrings_LtoH = Listpool_copy(substrings,listpool);
  new->substrings_HtoL = Listpool_copy(substrings,listpool);
  /* Do not use substrings after this */

  new->junctions_1toN = new->junctions_Nto1 = new->junctions_LtoH = /* new->junctions_HtoL = */ (List_T) NULL;
  /* There is no junctions_HtoL field */
  /* Do not use junctions after this */

#if 0
  /* No need to reverse for a single substring */
  if (gplusp == true) {
    /* Substrings, head to tail, are query low to high and genome low to high */
    new->substrings_HtoL = List_reverse(new->substrings_HtoL);
  } else {
    /* Substrings, head to tail, are query low to high and genome high to low */
    new->substrings_LtoH = List_reverse(new->substrings_LtoH);
    new->junctions_LtoH = List_reverse(new->junctions_LtoH);
  }
#endif

#ifdef DEBUG0
  printf("NEW SUBSTRING\n");
  printf("%d..%d\t%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\n",Substring_querystart(substring),Substring_queryend(substring),
    Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
    Substring_nmatches_to_trims(substring),Substring_amb_length(substring));
  printf("\n");
#endif
  

  genomicstart = Substring_genomicstart(substring);
  genomicend = Substring_genomicend(substring);
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;

  new->genestrand = genestrand;
  new->splice_score = 0.0;

  new->distant_splice_p = false;
  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->plusp = gplusp;
  new->sensedir_for_concordance = new->sensedir = sensedir;

  new->nindels = 0;
  new->nmismatches_refdiff = 0;	/* Set later */
  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(substring); /* Trimmed */
  /* new->nmismatches_refdiff = 0; */
  new->nsegments = List_length(new->substrings_1toN);

  new->trim_querystart = Substring_trim_querystart(substring);
  new->trim_queryend = Substring_trim_queryend(substring);
  new->trim_querystart_splicep = Substring_trim_querystart_splicep(substring);
  new->trim_queryend_splicep = Substring_trim_queryend_splicep(substring);
  debug0(printf("substrings trim_querystart %d, trim_queryend %d\n",new->trim_querystart,new->trim_queryend));


  new->nmatches_to_trims = 0;
  /* Note: Cannot use substrings variable here.  Need to use new->substrings_1toN */
  for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    new->nmatches_to_trims += Substring_nmatches_to_trims(substring);
  }
  debug0(printf("**Setting nmatches_to_trims to be %d\n",new->nmatches_to_trims));

  new->nmatches_plus_spliced_trims = new->nmatches_to_trims + Substring_amb_length(substring);
  assert(new->nmatches_plus_spliced_trims <= querylength);

  /* Used for global comparisons */
  new->score_within_trims = querylength - new->nmatches_plus_spliced_trims + Substring_amb_length(substring)/AMB_PENALTY;
  new->score_overall = querylength - new->nmatches_to_trims;

  /* found_score_overall does not compensate for spliced ends, so gives motivation to find distant splicing */
  if (new->score_overall < *found_score_overall) {
    *found_score_overall = new->score_overall;
  }
  /* found_score_within_trims does compensate for spliced trims, and guides how much further alignment is necessary */
  if (new->score_within_trims < *found_score_within_trims) {
    *found_score_within_trims = new->score_within_trims;
  }


  /* new->penalties = 0; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  new->nsplices = 0;
  new->distance = 0U;

  new->paired_usedp = false;
  new->concordantp = false;

  new->query_splicepos = -1;
  new->circularpos = compute_circularpos(&new->circularalias,new);

  if ((new->altlocp = altlocp[chrnum]) == false) {
    debug0(printf("Returning primary %p from Stage3end_new_terminal with found_score %d\n",new,*found_score_within_trims));
    debug0(printf("Method %s: Stage3end_new_terminal returning %p at %u..%u\n\n",
		  Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset));
    return new;

  } else {
    debug0(printf("Returning altloc %p from Stage3end_new_terminal with found_score %d\n",new,*found_score_within_trims));
    debug0(printf("Method %s: Stage3end_new_terminal returning %p at %u..%u\n\n",
		  Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset));
    return new;
  }
}



T
Stage3end_new_precomputed (int *found_score_overall, int *found_score_within_trims, int nmismatches_bothdiff,
			   List_T substrings, List_T junctions, List_T transcripts, List_T transcripts_other,
			   int querylength, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			   bool gplusp, int genestrand, int sensedir, Listpool_T listpool, Method_T method, int level) {
  T new;

  Univcoord_T genomicstart, genomicend;
  Substring_T substring, substring1, substringN;
  Junction_T junction;
  List_T p;
  int adj = 0;
  int nsites;
  double prob_total;


#ifdef DEBUG0
  printf("Entered Stage3end_new_precomputed, method %s, with gplusp %d\n",Method_string(method),gplusp);
  printf("%d substrings\n",List_length(substrings));
  printf("%d junctions\n",List_length(junctions));
#endif
  assert(List_length(substrings) == List_length(junctions) + 1);

  new = (T) MALLOC_OUT(sizeof(*new));
  new->hittype = SUBSTRINGS;
  new->method = method;
  new->level = level;

  new->querylength = querylength;
  new->querylength_adj = querylength + adj;

  /* Caller must not free these lists */
  new->transcripts = transcripts;
  new->transcripts_other = transcripts_other;

  /* Unlike Stage3end_new_substrings, where substrings and junctions
     are in opposite orders, substrings and junctions here are in the
     same order. */

  new->substrings_1toN = substrings;
  new->substrings_Nto1 = List_reverse(Listpool_copy(substrings,listpool));
  new->junctions_1toN = junctions;
  new->junctions_Nto1 = List_reverse(Listpool_copy(junctions,listpool));


  /* There is no junctions_HtoL field */

  if (gplusp == true) {
    /* Substrings, head to tail, are query low to high and genome low to high */
    new->substrings_LtoH = Listpool_copy(new->substrings_1toN,listpool);
    new->substrings_HtoL = Listpool_copy(new->substrings_Nto1,listpool);
    new->junctions_LtoH = Listpool_copy(new->junctions_1toN,listpool);
    /* new->junctions_HtoL = Listpool_copy(new->junctions_Nto1,listpool); */
  } else {
    /* Substrings, head to tail, are query low to high and genome high to low */
    new->substrings_LtoH = Listpool_copy(new->substrings_Nto1,listpool);
    new->substrings_HtoL = Listpool_copy(new->substrings_1toN,listpool);
    new->junctions_LtoH = Listpool_copy(new->junctions_Nto1,listpool);
    /* new->junctions_HtoL = Listpool_copy(new->junctions_1toN,listpool); */
  }
  /* Do not use substrings after this */
  /* Do not use junctions after this */



#ifdef DEBUG0
  printf("NEW SUBSTRINGS (query order)\n");
  for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = List_head(p);
    if (Substring_ambiguous_p(substring) == true) {
      printf("%d..%d\t%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\tprobs:%f and %f\n",Substring_querystart(substring),Substring_queryend(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring),Substring_amb_donor_prob(substring),Substring_amb_acceptor_prob(substring)
);
    } else {
      printf("%d..%d\t%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\n",Substring_querystart(substring),Substring_queryend(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring));
    }
  }
  printf("\n");

  printf("NEW SUBSTRINGS (genome order)\n");
  for (p = new->substrings_LtoH; p != NULL; p = List_next(p)) {
    substring = List_head(p);
    if (Substring_ambiguous_p(substring) == true) {
      printf("%d..%d\t%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\tprobs:%f and %f\n",Substring_querystart(substring),Substring_queryend(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring),Substring_amb_donor_prob(substring),Substring_amb_acceptor_prob(substring));
    } else {
      printf("%d..%d\t%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\n",Substring_querystart(substring),Substring_queryend(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring));
    }
  }
  printf("\n");

  printf("NEW JUNCTIONS (query order)\n");
  for (p = new->junctions_1toN; p != NULL; p = List_next(p)) {
    junction = List_head(p);
    printf("splice distance %u, nindels %d\n",Junction_splice_distance(junction),Junction_nindels(junction));
  }
  printf("\n");

  printf("NEW JUNCTIONS (genome order)\n");
  for (p = new->junctions_LtoH; p != NULL; p = List_next(p)) {
    junction = List_head(p);
    printf("splice distance %u, nindels %d\n",Junction_splice_distance(junction),Junction_nindels(junction));
  }
  printf("\n");
#endif
  

  substring1 = (Substring_T) List_head(new->substrings_1toN);
  substringN = (Substring_T) List_head(new->substrings_Nto1);

  genomicstart = Substring_genomicstart(substring1);
  genomicend = Substring_genomicend(substringN); /* DOESN'T WORK FOR AMBIGUOUS */
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;

  new->genestrand = genestrand;

  new->distant_splice_p = false;
  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->plusp = gplusp;

  new->nindels = 0;
  new->nmismatches_bothdiff = nmismatches_bothdiff; /* Trimmed */
  /* new->nmismatches_refdiff = 0; */
  new->nsegments = List_length(new->substrings_1toN);

  new->trim_querystart = Substring_trim_querystart(substring1);
  new->trim_queryend = Substring_trim_queryend(substringN);
  new->trim_querystart_splicep = Substring_trim_querystart_splicep(substring1);
  new->trim_queryend_splicep = Substring_trim_queryend_splicep(substringN);
  debug0(printf("substrings trim_querystart %d, trim_queryend %d\n",new->trim_querystart,new->trim_queryend));

  /* new->nmatches_to_trims = querylength_trimmed - nmismatches_whole; */
  new->nmatches_to_trims = 0;
  /* Note: Cannot use substrings variable here.  Need to use new->substrings_1toN */
  for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    new->nmatches_to_trims += Substring_nmatches_to_trims(substring);
  }
  debug0(printf("**Setting nmatches_to_trims to be %d\n",new->nmatches_to_trims));

  new->nmatches_plus_spliced_trims = new->nmatches_to_trims + Substring_start_amb_length(substring) + Substring_end_amb_length(substring);
  for (p = new->junctions_1toN; p != NULL; p = List_next(p)) {
    junction = List_head(p);
    new->nmatches_plus_spliced_trims += Junction_ninserts(junction);
  }
  assert(new->nmatches_plus_spliced_trims <= querylength);

  new->score_within_trims = querylength - new->nmatches_plus_spliced_trims +
    Substring_start_amb_length(substring1)/AMB_PENALTY + Substring_end_amb_length(substringN)/AMB_PENALTY;

  new->score_overall = querylength - new->nmatches_to_trims;

  /* found_score_overall does not compensate for spliced ends, so gives motivation to find distant splicing */
  if (new->score_overall < *found_score_overall) {
    *found_score_overall = new->score_overall;
  }
  /* found_score_within_trims does compensate for spliced trims, and guides how much further alignment is necessary */
  if (new->score_within_trims < *found_score_within_trims) {
    *found_score_within_trims = new->score_within_trims;
  }


  /* new->penalties = 0; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  prob_total = 0.0;
  nsites = 0;
  new->nsplices = 0;
  for (p = junctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == SPLICE_JUNCTION) {
      prob_total += Junction_splice_score(junction);
      nsites += 2;
      new->nsplices += 1;
    }
  }
  if (nsites == 0) {
    new->splice_score = 0.0;
  } else {
    new->splice_score = prob_total / (double) nsites;
  }
  debug0(printf("SPLICE SCORE: %f\n",new->splice_score));

  new->sensedir_for_concordance = new->sensedir = sensedir;

  new->distance = 0U;

  new->paired_usedp = false;
  new->concordantp = false;

  new->query_splicepos = -1;
  new->circularpos = compute_circularpos(&new->circularalias,new);

  if ((new->altlocp = altlocp[chrnum]) == false) {
    debug0(printf("Returning primary %p from Stage3end_new_precomputed with found_score %d\n",new,*found_score_within_trims));
    debug0(printf("Method %s: Stage3end_new_precomputed returning %p at %u..%u with splice_score %f\n\n",
		  Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset,new->splice_score));
    return new;

  } else {
    debug0(printf("Returning altloc %p from Stage3end_new_precomputed with found_score %d\n",new,*found_score_within_trims));
    debug0(printf("Method %s: Stage3end_new_precomputed returning %p at %u..%u with splice_score %f\n\n",
		  Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset,new->splice_score));
    return new;
  }
}


int
Stage3end_nmatches_substrings (Intlist_T endpoints, Univcoordlist_T lefts,
			       Intlist_T nmismatches_list, List_T junctions,
			       int querylength, Compress_T query_compress,
			       Substring_T qend_alts, Substring_T qstart_alts,
			       bool plusp, int genestrand,
			       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			       bool splice5p_in, bool splice3p_in, Listpool_T listpool) {
  int nmatches = 0, substring_nmatches;
  int qstart, qend;
  Univcoord_T left;
  Intlist_T r, x;
  Univcoordlist_T q;
  Junction_T junction;
#ifdef MAKE_JUNCTION
  Junction_T qstart_junction = NULL, qend_junction = NULL;
  double donor_prob, acceptor_prob;
#endif
  List_T newjunctions, p, j;
  bool splice5p, splice3p;
  int adj0;			/* deletions - insertions */
  int nmismatches, indel_score = 0, nindels = 0;
  int nindelbreaks, n_large_indels;
  /* double donor_prob, acceptor_prob; */


  debug7(printf("Entered Stage3end_nmatches_substrings with %s, plusp %d, splice5p %d, splice3p %d\n",
		Intlist_to_string(endpoints),plusp,splice5p_in,splice3p_in));

#ifdef DEBUG0
  printf("Entered Stage3end_nmatches_substrings, at left %u [%u], with chrnum #%d, plusp %d, and endpoints %s\n",
	 Univcoordlist_head(lefts),Univcoordlist_head(lefts) - chroffset,chrnum,plusp,Intlist_to_string(endpoints));
  printf("There are %d endpoints, %d lefts, %d nmismatches, and %d junctions\n",
	 Intlist_length(endpoints),Univcoordlist_length(lefts),Intlist_length(nmismatches_list),List_length(junctions));
  if (qstart_alts != NULL) {
    printf("qstart_alts at %d..%d\n",Substring_querystart(qstart_alts),Substring_queryend(qstart_alts));
  }
  if (qend_alts != NULL) {
    printf("qend_alts at %d..%d\n",Substring_querystart(qend_alts),Substring_queryend(qend_alts));
  }
  printf("Endpoints: %s\n",Intlist_to_string(endpoints));
  printf("Lefts: %s\n",Univcoordlist_to_string_offset(lefts,chroffset));
  printf("Mismatches: %s\n",Intlist_to_string(nmismatches_list));
#endif

  assert(Univcoordlist_length(lefts) == Intlist_length(endpoints) - 1);
  assert(Intlist_length(nmismatches_list) == Intlist_length(endpoints) - 1);
  assert(List_length(junctions) == Intlist_length(endpoints) - 2);

  
  newjunctions = Listpool_copy(junctions,listpool);


#ifdef DEBUG0
  for (p = junctions; p != NULL; p = List_next(p)) {
    Junction_print((Junction_T) List_head(p));
  }
  printf("\n");
#endif


  qstart = Intlist_head(endpoints);
  nmismatches = Intlist_head(nmismatches_list);

  if (plusp == true) {
    j = newjunctions;		/* Put here before we handle querystart_alts */
    if (qstart_alts != NULL) {
      debug7(printf("Adding %d matches for qstart_alts\n",Substring_nmatches(qstart_alts)));
      nmatches += Substring_nmatches(qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
#ifdef MAKE_JUNCTION
      donor_prob = Substring_amb_donor_prob(qstart_alts);
      acceptor_prob = Substring_amb_acceptor_prob(qstart_alts);
      qstart_junction = Junction_new_splice(/*distance*/0,orig_sensedir,donor_prob,acceptor_prob);
      newjunctions = Listpool_push(newjunctions,listpool,(void *) qstart_junction);
#else
      newjunctions = Listpool_push(newjunctions,listpool,(void *) NULL);
#endif
      splice5p = false;
    } else {
      splice5p = splice5p_in;
    }

    /* Add qpos to get alignstart/alignend */
    for (q = lefts, x = nmismatches_list, r = Intlist_next(endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
      qend = Intlist_head(r);
      nmismatches = Intlist_head(x);
      left = Univcoordlist_head(q);
      debug0(printf("Working on qstart %d..qend %d at left %u\n",qstart,qend,left));

      /* genomicstart = left; */
      /* genomicend = left + querylength; */
      /* alignstart = genomicstart + qstart; */
      /* alignend = genomicstart + qend; */

      if (nmismatches >= 0) {
	debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",left - chroffset,qstart,qend));
	debug7(printf("%d vs %d\n",nmismatches,
		      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							/*pos5*/qstart,/*pos3*/qend,/*plusp*/true,genestrand)));
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
								/*pos5*/qstart,/*pos3*/qend,/*plusp*/true,genestrand));
#endif
      } else {
	nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							/*pos5*/qstart,/*pos3*/qend,/*plusp*/true,genestrand);
	Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	debug7(printf("nmismatches %d from genome over querypos %d..%d\n",nmismatches,qstart,qend));
      }
      if (Univcoordlist_next(q) != NULL || qend_alts != NULL) {
	splice3p = false;
      } else {
	splice3p = splice3p_in;
      }

      if (splice5p == false && splice3p == false) {
	/* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
	debug7(printf("Shortcut computes matches of %d = (%d - %d) - nmismatches %d\n",
		      (qend-qstart)-nmismatches,qend,qstart,nmismatches));
	nmatches += (qend - qstart) - nmismatches;
      } else if ((substring_nmatches = Substring_compute_nmatches(left,/*querystart*/qstart,/*queryend*/qend,querylength,
								  /*plusp*/true,genestrand,query_compress,
								  chrnum,chroffset,chrhigh,chrlength,
								  /*splice_querystart_p*/splice5p,
								  /*splice_queryend_p*/splice3p,/*chrnum_fixed_p*/true)) < 0) {
	/* Don't know how to fix junctions */
	debug0(printf("Poor substring (plus) for %d..%d, so returning -1 from Stage3end_nmatches_substrings\n",
		      qstart,qend));
	return -1;
      } else {
	debug7(printf("Substring_compute_nmatches returns nmatches %d over querypos %d..%d\n",
		      substring_nmatches,qstart,qend));
	nmatches += substring_nmatches;
      }

      /* Prepare for next iteration */
      qstart = qend;
      if (j != NULL) {
	if ((junction = (Junction_T) List_head(j)) == NULL) {
	  /* qstart_junction */
	} else if ((adj0 = Junction_adj(junction)) != 0) {
	  /* adj += adj0; */
	  indel_score += indel_penalty_middle;
	  nindels += Junction_nindels(junction);
	  if (adj0 < 0) {
	    debug7(printf("Adjusting qstart %d up by %d\n",qstart,-adj0));
	    qstart -= adj0;	/* Insertion */
	  }
	}
      }
      splice5p = false;
    }

  } else {
    j = newjunctions;		/* Put here before we handle querystart_alts */
    if (qstart_alts != NULL) {
      debug7(printf("Adding %d matches for qstart_alts\n",Substring_nmatches(qstart_alts)));
      nmatches += Substring_nmatches(qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
#ifdef MAKE_JUNCTION
      donor_prob = Substring_amb_donor_prob(qstart_alts);
      acceptor_prob = Substring_amb_acceptor_prob(qstart_alts);
      qstart_junction = Junction_new_splice(/*distance*/0,orig_sensedir,donor_prob,acceptor_prob);
      /* printf("Creating junction with donor_prob %f and acceptor_prob %f\n",donor_prob,acceptor_prob); */
      newjunctions = Listpool_push(newjunctions,listpool,(void *) qstart_junction);
#else
      newjunctions = Listpool_push(newjunctions,listpool,(void *) NULL);
#endif
      splice5p = false;
    } else {
      splice5p = splice5p_in;
    }

    /* Subtract querypos to get alignstart/alignend */
    for (q = lefts, x = nmismatches_list, r = Intlist_next(endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
      qend = Intlist_head(r);
      nmismatches = Intlist_head(x);
      left = Univcoordlist_head(q);
      debug0(printf("Working on qstart %d..qend %d at left %u\n",qstart,qend,left));

      /* genomicend = left; */
      /* genomicstart = left + querylength; */
      /* genomicend_adj = genomicend - adj; */
      /* genomicstart_adj = genomicend - adj; */
      /* alignstart = genomicstart - (querylength - qend); */
      /* alignend = genomicstart - (querylength - qstart); */

      if (nmismatches >= 0) {
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
								/*pos5*/qstart,/*pos3*/qend,/*plusp*/false,genestrand));
#endif
      } else {
	nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							/*pos5*/qstart,/*pos3*/qend,/*plusp*/false,genestrand);
	Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	debug7(printf("nmismatches %d from genome over querypos %d..%d\n",nmismatches,querylength - qend,querylength - qstart));
      }
      if (Univcoordlist_next(q) != NULL || qend_alts != NULL) {
	splice3p = false;
      } else {
	splice3p = splice3p_in;
      }

      if (splice5p == false && splice3p == false) {
	/* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
	debug7(printf("Shortcut computes matches of %d = (%d - %d) - nmismatches %d\n",
		      (qend-qstart)-nmismatches,querylength - qstart,querylength - qend,nmismatches));
	nmatches += (qend - qstart) - nmismatches;
      } else if ((substring_nmatches = Substring_compute_nmatches(left,/*querystart*/querylength - qend,
								  /*queryend*/querylength - qstart,querylength,
								  /*plusp*/false,genestrand,query_compress,
								  chrnum,chroffset,chrhigh,chrlength,
								  /*splice_querystart_p*/splice3p,
								  /*splice_queryend_p*/splice5p,/*chrnum_fixed_p*/true)) < 0) {
	/* Don't know how to fix junctions */
	debug0(printf("Poor substring (minus) for querypos %d..%d, so returning -1 from Stage3end_new_substrings\n",
		      querylength - qend,querylength - qstart));
	return -1;
      } else {
	debug7(printf("Substring_compute_nmatches returns nmatches %d over querypos %d..%d\n",
		      substring_nmatches,querylength - qend,querylength - qstart));
	nmatches += substring_nmatches;
      }

      /* Prepare for next iteration */
      qstart = qend;
      if (j != NULL) {
	if ((junction = (Junction_T) List_head(j)) == NULL) {
	  /* qstart_junction */
	} else if ((adj0 = Junction_adj(junction)) != 0) {
	  /* adj += adj0; */
	  indel_score += indel_penalty_middle;
	  nindels += Junction_nindels(junction);
	  if (adj0 < 0) {
	    debug7(printf("Adjusting qstart %d up by %d\n",qstart,-adj0));
	    qstart -= adj0;	/* Insertion */
	  }
	}
      }
      splice5p = false;
    }
  }

  if (qend_alts != NULL) {
    debug7(printf("Adding %d matches for qend_alts\n",Substring_nmatches(qend_alts)));
    nmatches += Substring_nmatches(qend_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
#ifdef MAKE_JUNCTION
    newjunctions = List_reverse(newjunctions);
    donor_prob = Substring_amb_donor_prob(qend_alts);
    acceptor_prob = Substring_amb_acceptor_prob(qend_alts);
    qend_junction = Junction_new_splice(/*distance*/0,orig_sensedir,donor_prob,acceptor_prob);
    /* printf("Creating junction with donor_prob %f and acceptor_prob %f\n",donor_prob,acceptor_prob); */
    newjunctions = Listpool_push(newjunctions,listpool,(void *) qend_junction);
    newjunctions = List_reverse(newjunctions);
#endif
  }

    
  nindelbreaks = 0;
  n_large_indels = 0;

  for (p = newjunctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    /* CHIMERA_JUNCTION not possible */
    if (junction == NULL) {
      /* qstart_junction */
    } else if (Junction_type(junction) == SPLICE_JUNCTION) {
      /* No indel breaks.  ? Add penalty for bad splice probs */

    } else if (Junction_type(junction) == INS_JUNCTION) {
      nindelbreaks++;
      if (Junction_nindels(junction) > 6) {
	n_large_indels++;
      }
    } else if (Junction_type(junction) == DEL_JUNCTION) {
      nindelbreaks++;
      if (Junction_nindels(junction) > 6) {
	n_large_indels++;
      }
    }
  }

#if 0
  nmatches = nmatches - nindelbreaks*indel_penalty_middle - n_large_indels*3;
  for (p = newjunctions; p != NULL; p = List_next(p)) {
    if ((junction = List_head(p)) != NULL) {
      nmatches += Junction_ninserts(junction);
    }
  }
#endif


#ifdef MAKE_JUNCTION
  Junction_free(&qstart_junction);
  Junction_free(&qend_junction);
#endif
  
  debug7(printf("Stage3end_nmatches_substrings returning %d matches\n",nmatches));
  /* List_free(&newjunctions); -- allocated by Listpool_push */

  assert(nmatches <= querylength);
  return nmatches;
}



/* endpoints are all in qstart/qend convention.  Need to convert to
   querystart and queryend when creating Substring_T objects */
/* Three actions at each end: extend, chop, or compute_trim */
T
Stage3end_new_substrings (int *found_score_overall, int *found_score_within_trims,
			  Intlist_T endpoints, Univcoordlist_T lefts,
			  Intlist_T nmismatches_list, List_T junctions,
			  int querylength, Compress_T query_compress,
			  Substring_T qend_alts, Substring_T qstart_alts,
			  bool plusp, int genestrand, int sensedir,
			  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			  bool splice5p_in, Splicetype_T splicetype5, double ambig_prob_5,
			  bool splice3p_in, Splicetype_T splicetype3, double ambig_prob_3,
			  Listpool_T listpool, Method_T method, int level) {
  T new;

  Univcoord_T genomicstart, genomicend;
  int querylength_trimmed = 0;
  int qstart, qend;
  Univcoord_T left;
  Intlist_T r, x;
  Univcoordlist_T q;
  Substring_T substring, substring1, substringN;
  Junction_T junction;
  List_T substrings = NULL, p, j;
  List_T newjunctions;
  bool splice5p, splice3p;
  int adj = 0, adj0;			/* deletions - insertions */
  int nmismatches, indel_score = 0, nindels = 0;
  int nmismatches_bothdiff = 0;
  int new_sensedir;
  bool contradictionp;
  int nsites, nindelbreaks, n_large_indels;
  double prob_total, donor_prob, acceptor_prob;


  debug7(printf("Entered Stage3end_new_substrings with %s, plusp %d, splice5p %d, splice3p %d\n",
		Intlist_to_string(endpoints),plusp,splice5p_in,splice3p_in));

#ifdef DEBUG0
  printf("Entered Stage3end_new_substrings, at left %u [%u], with chrnum #%d, plusp %d, sensedir %d, and endpoints %s\n",
	 Univcoordlist_head(lefts),Univcoordlist_head(lefts) - chroffset,chrnum,plusp,sensedir,Intlist_to_string(endpoints));
  printf("There are %d endpoints, %d lefts, %d nmismatches, and %d junctions\n",
	 Intlist_length(endpoints),Univcoordlist_length(lefts),Intlist_length(nmismatches_list),List_length(junctions));
  if (qstart_alts != NULL) {
    printf("qstart_alts at %d..%d.  ",Substring_querystart(qstart_alts),Substring_queryend(qstart_alts));
    Substring_print_alts_coords(qstart_alts);
    printf("\n");
  }
  if (qend_alts != NULL) {
    printf("qend_alts at %d..%d.  ",Substring_querystart(qend_alts),Substring_queryend(qend_alts));
    Substring_print_alts_coords(qend_alts);
    printf("\n");
  }
  printf("Endpoints: %s\n",Intlist_to_string(endpoints));
  printf("Lefts: %s\n",Univcoordlist_to_string_offset(lefts,chroffset));
  printf("Mismatches: %s\n",Intlist_to_string(nmismatches_list));
#endif

  assert(Univcoordlist_length(lefts) == Intlist_length(endpoints) - 1);
  assert(Intlist_length(nmismatches_list) == Intlist_length(endpoints) - 1);
  assert(List_length(junctions) == Intlist_length(endpoints) - 2);

  
  newjunctions = Junction_copy_list(junctions,listpool);

#ifdef DEBUG0
  for (p = newjunctions; p != NULL; p = List_next(p)) {
    Junction_print((Junction_T) List_head(p));
  }
  printf("\n");
#endif

  qstart = Intlist_head(endpoints);
  nmismatches = Intlist_head(nmismatches_list);

  if (plusp == true) {
    j = newjunctions;		/* Put here before we handle qstart_alts */
    if (qstart_alts != NULL) {
      substrings = Listpool_push(substrings,listpool,(void *) Substring_copy(qstart_alts));
      donor_prob = Substring_amb_donor_prob(qstart_alts);
      acceptor_prob = Substring_amb_acceptor_prob(qstart_alts);
      junction = Junction_new_splice(/*distance*/0,sensedir,donor_prob,acceptor_prob);
      newjunctions = Listpool_push(newjunctions,listpool,(void *) junction);
      splice5p = false;
    } else {
      splice5p = splice5p_in;
    }

    /* Add qpos to get alignstart/alignend */
    for (q = lefts, x = nmismatches_list, r = Intlist_next(endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
      qend = Intlist_head(r);
      nmismatches = Intlist_head(x);
      left = Univcoordlist_head(q);
      debug0(printf("Working on qstart %d..qend %d at left %u\n",qstart,qend,left));

      /* genomicstart = left; */
      /* genomicend = left + querylength; */
      /* alignstart = genomicstart + qstart; */
      /* alignend = genomicstart + queryend; */

      if (nmismatches >= 0) {
	debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",left - chroffset,qstart,qend));
	debug7(printf("%d vs %d\n",nmismatches,
		      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							/*pos5*/qstart,/*pos3*/qend,/*plusp*/true,genestrand)));
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
								/*pos5*/qstart,/*pos3*/qend,/*plusp*/true,genestrand));
#endif
      } else {
	nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							/*pos5*/qstart,/*pos3*/qend,/*plusp*/true,genestrand);
      }
      if (Univcoordlist_next(q) != NULL || qend_alts != NULL) {
	splice3p = false;
      } else {
	splice3p = splice3p_in;
      }

      if ((substring = Substring_new(nmismatches,left,/*querystart*/qstart,/*queryend*/qend,querylength,
				     /*plusp*/true,genestrand,query_compress,
				     chrnum,chroffset,chrhigh,chrlength,
				     /*splice_querystart_p*/splice5p,/*splicetype_querystart*/splicetype5,
				     /*ambig_prob_querystart*/ambig_prob_5,
				     /*splice_queryend_p*/splice3p,/*splicetype_queryend*/splicetype3,
				     /*ambig_prob_queryend*/ambig_prob_3,sensedir)) == NULL) {
	/* Don't know how to fix junctions */
	debug0(printf("Poor substring (plus) for %d..%d, so returning NULL from Stage3end_new_substrings\n",
		      qstart,qend));
	for (p = substrings; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (substring == qstart_alts) {
	    /* qstart_alts freed by calling procedure.  Need to free junction created for querystart_alts. */
	    /* junctions = List_pop(junctions,(void **) &junction); */
	    /* Junction_free(&junction); */
	  } else {
	    Substring_free(&substring);
	  }
	}
	/* List_free(&substrings); -- allocated by Listpool_push */
	debug0(printf("Stage3end_new_substrings returning NULL\n"));
	Junction_gc(&newjunctions);
	return (T) NULL;

      } else {
	debug7(printf("Substring_new returns nmismatches %d, nmatches %d, ambp %d, amb %d over querypos %d..%d\n",
		      Substring_nmismatches_bothdiff(substring),
		      Substring_nmatches(substring),Substring_ambiguous_p(substring),
		      Substring_amb_length(substring),Substring_querystart(substring),Substring_queryend(substring)));

	debug0(printf("Substring_new returns nmismatches %d, nmatches %d, ambp %d, amb %d over querypos %d..%d\n",
		      Substring_nmismatches_bothdiff(substring),
		      Substring_nmatches(substring),Substring_ambiguous_p(substring),
		      Substring_amb_length(substring),Substring_querystart(substring),Substring_queryend(substring)));
	substrings = Listpool_push(substrings,listpool,(void *) substring);
	nmismatches_bothdiff += Substring_nmismatches_bothdiff(substring);
	querylength_trimmed += Substring_match_length(substring);
      }

      /* Prepare for next iteration */
      qstart = qend;
      if (j != NULL) {
	junction = (Junction_T) List_head(j);
	if ((adj0 = Junction_adj(junction)) != 0) {
	  adj += adj0;
	  indel_score += indel_penalty_middle;
	  nindels += Junction_nindels(junction);
	  if (adj0 < 0) {
	    qstart -= adj0;	/* Insertion */
	  }
	}
      }
      splice5p = false;
    }

  } else {
    j = newjunctions;		/* Put here before we handle querystart_alts */
    if (qstart_alts != NULL) {
      substrings = Listpool_push(substrings,listpool,(void *) Substring_copy(qstart_alts));
      donor_prob = Substring_amb_donor_prob(qstart_alts);
      acceptor_prob = Substring_amb_acceptor_prob(qstart_alts);
      junction = Junction_new_splice(/*distance*/0,sensedir,donor_prob,acceptor_prob);
      /* printf("Creating junction with donor_prob %f and acceptor_prob %f\n",donor_prob,acceptor_prob); */
      newjunctions = Listpool_push(newjunctions,listpool,(void *) junction);
      splice5p = false;
    } else {
      splice5p = splice5p_in;
    }

    /* Subtract qpos to get alignstart/alignend */
    for (q = lefts, x = nmismatches_list, r = Intlist_next(endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
      qend = Intlist_head(r);
      nmismatches = Intlist_head(x);
      left = Univcoordlist_head(q);
      debug0(printf("Working on qstart %d..qend %d at left %u\n",qstart,qend,left));

      /* genomicend = left; */
      /* genomicstart = left + querylength; */
      /* genomicend_adj = genomicend - adj; */
      /* genomicstart_adj = genomicend - adj; */
      /* alignstart = genomicstart - (querylength - qend); */
      /* alignend = genomicstart - (querylength - qstart); */

      if (nmismatches >= 0) {
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
								/*pos5*/qstart,/*pos3*/qend,/*plusp*/false,genestrand));
#endif
      } else {
	nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							/*pos5*/qstart,/*pos3*/qend,/*plusp*/false,genestrand);
      }
      if (Univcoordlist_next(q) != NULL || qend_alts != NULL) {
	splice3p = false;
      } else {
	splice3p = splice3p_in;
      }

      if ((substring = Substring_new(nmismatches,left,/*querystart*/querylength - qend,
				     /*queryend*/querylength - qstart,querylength,
				     /*plusp*/false,genestrand,query_compress,
				     chrnum,chroffset,chrhigh,chrlength,
				     /*splice_querystart_p*/splice3p,/*splicetype_querystart*/splicetype3,
				     /*ambig_prob_querystart*/ambig_prob_3,
				     /*splice_queryend_p*/splice5p,/*splicetype_queryend*/splicetype5,
				     /*ambig_prob_queryend*/ambig_prob_5,sensedir)) == NULL) {
	/* Don't know how to fix junctions */
	debug0(printf("Poor substring (minus) for %d..%d, so returning NULL from Stage3end_new_substrings\n",
		      querylength - qend,querylength - qstart));
	for (p = substrings; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (substring == qstart_alts) {
	    /* querystart_alts freed by calling procedure.  Need to free junction created for querystart_alts. */
	    /* junctions = List_pop(junctions,(void **) &junction); */
	    /* Junction_free(&junction); */
	  } else {
	    Substring_free(&substring);
	  }
	}
	/* List_free(&substrings); -- allocated by Listpool_push */

	debug0(printf("Stage3end_new_substrings returning NULL\n"));
	Junction_gc(&newjunctions);
	return (T) NULL;

      } else {
	debug7(printf("Substring_new returns nmismatches %d, nmatches %d, ambp %d, amb %d over querypos %d..%d\n",
		      Substring_nmismatches_bothdiff(substring),
		      Substring_nmatches(substring),Substring_ambiguous_p(substring),
		      Substring_amb_length(substring),Substring_querystart(substring),Substring_queryend(substring)));

	debug0(printf("Substring_new returns nmismatches %d, nmatches %d, ambp %d, amb %d over querypos %d..%d\n",
		      Substring_nmismatches_bothdiff(substring),
		      Substring_nmatches(substring),Substring_ambiguous_p(substring),
		      Substring_amb_length(substring),Substring_querystart(substring),Substring_queryend(substring)));
	substrings = Listpool_push(substrings,listpool,(void *) substring);
	nmismatches_bothdiff += Substring_nmismatches_bothdiff(substring);
	querylength_trimmed += Substring_match_length(substring);
      }

      /* Prepare for next iteration */
      qstart = qend;
      if (j != NULL) {
	junction = (Junction_T) List_head(j);
	if ((adj0 = Junction_adj(junction)) != 0) {
	  adj += adj0;
	  indel_score += indel_penalty_middle;
	  nindels += Junction_nindels(junction);
	  if (adj0 < 0) {
	    qstart -= adj0;	/* Insertion */
	  }
	}
      }
      splice5p = false;
    }
  }

  if (qend_alts != NULL) {
    substrings = Listpool_push(substrings,listpool,(void *) Substring_copy(qend_alts));
    newjunctions = List_reverse(newjunctions);
    donor_prob = Substring_amb_donor_prob(qend_alts);
    acceptor_prob = Substring_amb_acceptor_prob(qend_alts);
    junction = Junction_new_splice(/*distance*/0,sensedir,donor_prob,acceptor_prob);
    /* printf("Creating junction with donor_prob %f and acceptor_prob %f\n",donor_prob,acceptor_prob); */
    newjunctions = Listpool_push(newjunctions,listpool,(void *) junction);
    newjunctions = List_reverse(newjunctions);
  }

#ifdef DEBUG0
  printf("NEW JUNCTIONS\n");
  for (p = newjunctions; p != NULL; p = List_next(p)) {
    Junction_print(List_head(p));
  }
  printf("\n");
#endif


  if (plusp == true) {
    substring1 = List_last_value(substrings);
    substringN = List_head(substrings);
  } else {
    substring1 = List_head(substrings);
    substringN = List_last_value(substrings);
  }

  debug0(printf("Trim left: %d.  Trim right: %d\n",
		Substring_trim_querystart(substring1),Substring_trim_queryend(substringN)));
  if (Substring_chrnum(substring1) != Substring_chrnum(substringN)) {
    debug0(printf("ABORTING BECAUSE SUBSTRINGS HAVE DIFFERENT CHRNUMS: %d AND %d\n",
		  Substring_chrnum(substring1),Substring_chrnum(substringN)));
    for (p = substrings; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      if (substring == qstart_alts || substring == qend_alts) {
	/* qstart_alts and qend_alts freed by calling procedure */
      } else {
	Substring_free(&substring);
      }
    }
    /* List_free(&substrings); -- allocated by Listpool_push */

    debug0(printf("Stage3end_new_substrings returning NULL\n"));
    Junction_gc(&newjunctions);
    return (T) NULL;
  }

  new = (T) MALLOC_OUT(sizeof(*new));
  new->hittype = SUBSTRINGS;
  new->method = method;
  new->level = level;

  new->transcripts = (List_T) NULL;
  new->transcripts_other = (List_T) NULL;

  new->querylength = querylength;
  new->querylength_adj = querylength + adj;

  /* Note differences between substrings and junctions.  Substrings
     were pushed onto lists above, and junctions were created by the
     caller, so they are originally in opposite orders */
  new->substrings_HtoL = substrings;
  new->substrings_LtoH = List_reverse(Listpool_copy(substrings,listpool));
  new->junctions_LtoH = newjunctions;
  /* new->junctions_HtoL = List_reverse(Listpool_copy(newjunctions,listpool)); */

  if (plusp == true) {
    new->substrings_1toN = Listpool_copy(new->substrings_LtoH,listpool);
    new->substrings_Nto1 = Listpool_copy(new->substrings_HtoL,listpool);

    new->junctions_1toN = Listpool_copy(new->junctions_LtoH,listpool);
    new->junctions_Nto1 = List_reverse(Listpool_copy(new->junctions_LtoH,listpool));

  } else {
    new->substrings_1toN = Listpool_copy(new->substrings_HtoL,listpool);
    new->substrings_Nto1 = Listpool_copy(new->substrings_LtoH,listpool);

    new->junctions_1toN = List_reverse(Listpool_copy(new->junctions_LtoH,listpool));
    new->junctions_Nto1 = Listpool_copy(new->junctions_LtoH,listpool);
  }


#ifdef DEBUG0
  printf("NEW SUBSTRINGS\n");
  for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = List_head(p);
    if (Substring_has_alts_p(substring) == true) {
      printf("%d..%d\t#%d\talts\tmatches_to_trims: %d\tamb:%d\t%d common_prob:%f alts:",
	     Substring_querystart(substring),Substring_queryend(substring),Substring_chrnum(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring),
	     Substring_alts_ncoords(substring),Substring_alts_common_prob(substring));
      Substring_print_alts_coords(substring);
      printf("\n");

    } else if (Substring_ambiguous_p(substring) == true) {
      printf("%d..%d\t#%d\t%u..%u\tambig\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\tprobs:%f and %f\n",
	     Substring_querystart(substring),Substring_queryend(substring),Substring_chrnum(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),
	     Substring_nmismatches_bothdiff(substring),Substring_nmatches_to_trims(substring),Substring_amb_length(substring),
	     Substring_amb_donor_prob(substring),Substring_amb_acceptor_prob(substring));
    } else {
      printf("%d..%d\t#%d\t%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\n",
	     Substring_querystart(substring),Substring_queryend(substring),Substring_chrnum(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),
	     Substring_nmismatches_bothdiff(substring),Substring_nmatches_to_trims(substring),
	     Substring_amb_length(substring));
    }
  }
  printf("\n");
#endif
  

  substring1 = (Substring_T) List_head(new->substrings_1toN);
  substringN = (Substring_T) List_head(new->substrings_Nto1);

  genomicstart = Substring_genomicstart(substring1);
  genomicend = Substring_genomicend(substringN); /* DOESN'T WORK FOR AMBIGUOUS */
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;

  new->genestrand = genestrand;

  new->distant_splice_p = false;
  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->plusp = plusp;

  if (sensedir != SENSE_NULL) {
    debug0(printf("sensedir is %d (original)\n",sensedir));
    new->sensedir = sensedir;
  } else {
    new->sensedir = SENSE_NULL;
    contradictionp = false;
    for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      debug0(printf("substring has sensedir %d\n",Substring_sensedir(substring)));
      if (Substring_sensedir(substring) == SENSE_NULL) {
	/* Ignore */
      } else if (new_sensedir == SENSE_NULL) {
	new_sensedir = Substring_sensedir(substring);
      } else if (Substring_sensedir(substring) != new_sensedir) {
	contradictionp = true;
      }
    }
  
    for (p = new->junctions_1toN; p != NULL; p = List_next(p)) {
      junction = (Junction_T) List_head(p);
      debug0(printf("junction has sensedir %d\n",Junction_sensedir(junction)));
      if (Junction_sensedir(junction) == SENSE_NULL) {
	/* Ignore.  Probably an indel. */
      } else if (new_sensedir == SENSE_NULL) {
	new_sensedir = Junction_sensedir(junction);
      } else if (Junction_sensedir(junction) != new_sensedir) {
	contradictionp = true;
      }
    }

    if (contradictionp == true) {
      debug0(printf("CONTRADICTION IN SENSEDIR\n"));
      new->sensedir = SENSE_NULL;
    } else {
      debug0(printf("sensedir is %d\n",new_sensedir));
      new->sensedir = new_sensedir;
    }
  }
  new->sensedir_for_concordance = new->sensedir;

  prob_total = 0.0;
  nsites = 0;
  if (splice5p_in == true) {
    prob_total += ambig_prob_5;
    nsites++;
  }
  if (splice3p_in == true) {
    prob_total += ambig_prob_3;
    nsites++;
  }

  new->nsplices = 0;
  for (p = newjunctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == SPLICE_JUNCTION) {
      prob_total += Junction_splice_score(junction);
      nsites += 2;
      new->nsplices += 1;
    }
  }
  if (nsites == 0) {
    new->splice_score = 0.0;
  } else {
    new->splice_score = prob_total / (double) nsites;
  }
  debug0(printf("SPLICE SCORE: %f\n",new->splice_score));


  nindelbreaks = 0;
  n_large_indels = 0;
  for (p = newjunctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    /* CHIMERA_JUNCTION not possible */
    if (Junction_type(junction) == INS_JUNCTION) {
      nindelbreaks++;
      if (Junction_nindels(junction) > 6) {
	n_large_indels++;
      }
    } else if (Junction_type(junction) == DEL_JUNCTION) {
      nindelbreaks++;
      if (Junction_nindels(junction) > 6) {
	n_large_indels++;
      }
    }
  }


  /* nmismatches_bothdiff is computed after trimming */
  new->nindels = nindels;
  new->nmismatches_bothdiff = nmismatches_bothdiff; /* Trimmed */
  /* new->nmismatches_refdiff = 0; */
  new->nsegments = List_length(new->substrings_1toN);

  new->trim_querystart = Substring_trim_querystart(substring1);
  new->trim_queryend = Substring_trim_queryend(substringN);
  new->trim_querystart_splicep = Substring_trim_querystart_splicep(substring1);
  new->trim_queryend_splicep = Substring_trim_queryend_splicep(substringN);
  debug0(printf("substrings trim_querystart %d, trim_queryend %d\n",new->trim_querystart,new->trim_queryend));


  new->nmatches_to_trims = 0;
  /* Note: Cannot use substrings variable here.  Need to use new->substrings_1toN */
  for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    new->nmatches_to_trims += Substring_nmatches_to_trims(substring);
  }
  debug0(printf("Setting nmatches_to_trims to be %d\n",new->nmatches_to_trims));

  new->nmatches_plus_spliced_trims = new->nmatches_to_trims + Substring_start_amb_length(substring1) + Substring_end_amb_length(substringN);
  debug0(printf("Setting nmatches_plus_spliced_trims to be %d = %d + %d + %d\n",
		new->nmatches_plus_spliced_trims,new->nmatches_to_trims,
		Substring_start_amb_length(substring1),Substring_end_amb_length(substringN)));

  for (p = new->junctions_1toN; p != NULL; p = List_next(p)) {
    junction = List_head(p);
    new->nmatches_plus_spliced_trims += Junction_ninserts(junction);
  }
  assert(new->nmatches_plus_spliced_trims >= 0);
  assert(new->nmatches_plus_spliced_trims <= querylength);

  new->score_overall = querylength - new->nmatches_to_trims;
  new->score_within_trims = querylength - new->nmatches_plus_spliced_trims + 
    Substring_start_amb_length(substring1)/AMB_PENALTY + Substring_end_amb_length(substringN)/AMB_PENALTY;

  /* found_score_overall does not compensate for spliced ends, so gives motivation to find distant splicing */
  if (new->score_overall < *found_score_overall) {
    *found_score_overall = new->score_overall;
  }
  /* found_score_within_trims does compensate for spliced trims, and guides how much further alignment is necessary */
  if (new->score_within_trims < *found_score_within_trims) {
    *found_score_within_trims = new->score_within_trims;
  }


  /* new->penalties = 0; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  new->distance = 0U;

  new->paired_usedp = false;
  new->concordantp = false;

  new->query_splicepos = -1;
  new->circularpos = compute_circularpos(&new->circularalias,new);

  debug0(printf("%d substrings\n",List_length(new->substrings_1toN)));
  debug0(printf("%d junctions\n",List_length(new->junctions_1toN)));
  assert(List_length(new->substrings_1toN) == List_length(new->junctions_1toN) + 1);


  /* Previously checked for (new->circularalias == +2 || new->circularalias == -2) */

  debug7(printf("Stage3end_new_substrings returning %d matches_plus_spliced_trims\n",new->nmatches_plus_spliced_trims));

  if (new->circularpos >= 0) {
    new->altlocp = false;
    debug0(printf("*****Returning circular %p from Stage3end_new_substrings with score %d (found_score %d), nmatches %d\n\n",
		  new,new->score_within_trims,*found_score_within_trims,new->nmatches_plus_spliced_trims));
    debug0(printf("Method %s: Stage3end_new_substrings returning %p at %u..%u with %d nmismatches_bothdiff\n\n",
		  Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset,new->nmismatches_bothdiff));
    return new;
    
  } else if ((new->altlocp = altlocp[chrnum]) == false) {
    debug0(printf("*****Returning primary %p from Stage3end_new_substrings with score %d within trims, %d overall (found_score %d), nmatches %d\n\n",
		  new,new->score_within_trims,new->score_overall,*found_score_within_trims,new->nmatches_plus_spliced_trims));
    debug0(printf("Method %s: Stage3end_new_substrings returning %p at %u..%u with splice score %f, nmatches %d\n\n",
		  Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset,new->splice_score,
		  new->nmatches_plus_spliced_trims));
    return new;

  } else {
    debug0(printf("*****Returning altloc %p from Stage3end_new_substrings with score %d within trims, %d overall (found_score %d), nmatches %d\n\n",
		  new,new->score_within_trims,new->score_overall,*found_score_within_trims,new->nmatches_plus_spliced_trims));
    debug0(printf("Method %s: Stage3end_new_substrings returning %p at %u..%u with splice score %f, nmatches %d\n\n",
		  Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset,new->splice_score,
		  new->nmatches_plus_spliced_trims));
    return new;
  }
}


#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


T
Stage3end_new_substitution (int *found_score_overall, int *found_score_within_trims,
			    Univcoord_T left, int genomiclength, int querylength,
			    int *mismatch_positions_alloc, Compress_T query_compress,
			    bool plusp, int genestrand, int sensedir, int nmismatches_allowed,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			    Chrpos_T chrlength, Listpool_T listpool, Method_T method, int level) {
  T new;
  Substring_T substring;
  int qstart, qend, nmismatches;
  bool splice_querystart_p, splice_queryend_p;
  Splicetype_T splicetype_querystart, splicetype_queryend;
  double ambig_prob_querystart, ambig_prob_queryend;


  debug0(printf("Entered Stage3end_new_substitution, sensedir %d, method %s, at left %u [%u] and chrhigh %u\n",
		sensedir,Method_string(method),left,left - chroffset,chrhigh));

  if (plusp == true) {
    splice_querystart_p = Substring_qstart_trim(&qstart,&splicetype_querystart,&ambig_prob_querystart,
						left,/*qstart:0,*//*qend*/genomiclength,/*querylength:genomiclength,*/
						plusp,genestrand,mismatch_positions_alloc,query_compress,chroffset,sensedir);
    splice_queryend_p = Substring_qend_trim(&qend,&splicetype_queryend,&ambig_prob_queryend,
					    left,/*qstart*/0,/*qend:genomiclength,*//*querylength*/genomiclength,
					    plusp,genestrand,mismatch_positions_alloc,query_compress,chroffset,sensedir);

    debug0(printf("Trimming querystart yields splicep %d, qstart %d, prob %f\n",splice_querystart_p,qstart,ambig_prob_querystart));
    debug0(printf("Trimming queryend yields splicep %d, qend %d, prob %f\n",splice_queryend_p,qend,ambig_prob_queryend));

    if (qstart < 0 || qend < 0) {
      debug0(printf("Returning NULL\n"));
      return (T) NULL;

    } else if (qend <= qstart) {
      /* Otherwise, calling Genome_count_mismatches_substring will not be defined */
      debug0(printf("Returning NULL\n"));
      return (T) NULL;

    } else if ((nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
								/*pos5*/qstart,/*pos3*/qend,
								/*plusp*/true,genestrand)) > nmismatches_allowed) {
      debug0(printf("Returning NULL\n"));
      return (T) NULL;

    } else if ((substring = Substring_new(nmismatches,left,/*querystart*/qstart,/*queryend*/qend,
					  /*querylength*/genomiclength,/*plusp*/true,genestrand,query_compress,
					  chrnum,chroffset,chrhigh,chrlength,
					  splice_querystart_p,splicetype_querystart,ambig_prob_querystart,
					  splice_queryend_p,splicetype_queryend,ambig_prob_queryend,
					  sensedir)) == NULL) {
      debug0(printf("Returning NULL\n"));
      return (T) NULL;
    }

  } else {
    /* trim_querystart and trim_queryend Genome_count_mismatches_substring are flipped, but not for Substring_new */
    splice_querystart_p = Substring_qend_trim(&qend,&splicetype_querystart,&ambig_prob_querystart,
					      left,/*querystart*/0,/*qend:genomiclength,*//*querylength*/genomiclength,
					      plusp,genestrand,mismatch_positions_alloc,query_compress,chroffset,sensedir);
    splice_queryend_p = Substring_qstart_trim(&qstart,&splicetype_queryend,&ambig_prob_queryend,
					      left,/*qstart:0,*//*qend*/genomiclength,/*querylength:genomiclength,*/
					      plusp,genestrand,mismatch_positions_alloc,query_compress,chroffset,sensedir);

    debug0(printf("Trimming querystart yields splicep %d, qstart %d, prob %f\n",splice_querystart_p,qstart,ambig_prob_querystart));
    debug0(printf("Trimming queryend yields splicep %d, qend %d, prob %f\n",splice_queryend_p,qend,ambig_prob_queryend));

    if (qstart < 0 || qend < 0) {
      debug0(printf("Returning NULL\n"));
      return (T) NULL;

    } else if (qend <= qstart) {
      /* Otherwise, calling Genome_count_mismatches_substring will not be defined */
      debug0(printf("Returning NULL\n"));
      return (T) NULL;

    } else if ((nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
								/*pos5*/qstart,/*pos3*/qend,
								/*plusp*/false,genestrand)) > nmismatches_allowed) {
      debug0(printf("Returning NULL\n"));
      return (T) NULL;

    } else if ((substring = Substring_new(nmismatches,left,/*querystart*/querylength - qend,/*queryend*/querylength - qstart,
					  /*querylength*/genomiclength,/*plusp*/false,genestrand,query_compress,
					  chrnum,chroffset,chrhigh,chrlength,
					  splice_querystart_p,splicetype_querystart,ambig_prob_querystart,
					  splice_queryend_p,splicetype_queryend,ambig_prob_queryend,
					  sensedir)) == NULL) {
      debug0(printf("Returning NULL\n"));
      return (T) NULL;
    }
  }

  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_substitution %p: left %llu, chrnum %d, nmismatches %d\n",
		new,(unsigned long long) left,Substring_chrnum(substring),nmismatches));
  
  new->substrings_LtoH = Listpool_push(NULL,listpool,(void *) substring);
  new->substrings_HtoL = Listpool_push(NULL,listpool,(void *) substring);
  new->substrings_1toN = Listpool_push(NULL,listpool,(void *) substring);
  new->substrings_Nto1 = Listpool_push(NULL,listpool,(void *) substring);
  
  new->junctions_LtoH = (List_T) NULL;
  /* new->junctions_HtoL = (List_T) NULL; */
  new->junctions_1toN = (List_T) NULL;
  new->junctions_Nto1 = (List_T) NULL;
  
  new->transcripts = (List_T) NULL;
  new->transcripts_other = (List_T) NULL;
  
  new->querylength_adj = new->querylength = genomiclength;
  if (plusp == true) {
    new->low = new->genomicstart = left;
    new->high = new->genomicend = left + genomiclength;
    
  } else {
    new->low = new->genomicend = left;
    new->high = new->genomicstart = left + genomiclength;
  }
  
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;
  
#if 0
  if (nmismatches == 0) {
    /* Proper hittype needed so we can eliminate identical hits */
    new->hittype = EXACT;
  } else {
    new->hittype = SUB;
  }
#else
  new->hittype = SUB;
#endif
  new->method = method;
  new->level = level;
  
  new->genestrand = genestrand;
  
  /* Note: It is possible that Substring_new has assigned a new chrnum, different from the one given */
  new->distant_splice_p = false;
  new->chrnum = new->effective_chrnum = Substring_chrnum(substring);
  new->other_chrnum = 0;
  new->chroffset = Substring_chroffset(substring);
  new->chrhigh = Substring_chrhigh(substring);
  new->chrlength = Substring_chrlength(substring);
  new->plusp = plusp;
  new->sensedir_for_concordance = new->sensedir = sensedir;
  
#if 0
  new->mapq_loglik = Substring_mapq_loglik(substring);
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif
  
  new->nindels = 0;
  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(substring);
  new->nsegments = 1;
  
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1); */
  
  if (plusp) {
    new->trim_querystart = qstart;
    new->trim_queryend = querylength - qend;
  } else {
    new->trim_querystart = querylength - qend;
    new->trim_queryend = qstart;
  }
  new->trim_querystart_splicep = splice_querystart_p;
  new->trim_queryend_splicep = splice_queryend_p;
  debug0(printf("  trim on left: %d (splicep %d)\n",new->trim_querystart,new->trim_querystart_splicep));
  debug0(printf("  trim on right: %d (splicep %d)\n",new->trim_queryend,new->trim_queryend_splicep));

  
  new->nmatches_to_trims = Substring_nmatches_to_trims(substring);
  new->nmatches_plus_spliced_trims = new->nmatches_to_trims + Substring_amb_length(substring);
  assert(new->nmatches_plus_spliced_trims <= querylength);

  new->score_within_trims = querylength - new->nmatches_plus_spliced_trims + Substring_amb_length(substring)/AMB_PENALTY;
  new->score_overall = querylength - new->nmatches_to_trims;
  
  /* found_score_overall does not compensate for spliced ends, so gives motivation to find distant splicing */
  if (new->score_overall < *found_score_overall) {
    *found_score_overall = new->score_overall;
  }
  /* found_score_within_trims does compensate for spliced trims, and guides how much further alignment is necessary */
  if (new->score_within_trims < *found_score_within_trims) {
    *found_score_within_trims = new->score_within_trims;
  }
  
  
  /* new->penalties = 0; */
  
  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;
  
  new->nsplices = 0;
  if (splice_querystart_p == true && splice_queryend_p == true) {
    new->splice_score = (ambig_prob_querystart + ambig_prob_queryend)/2.0;
  } else if (splice_querystart_p == true) {
    new->splice_score = ambig_prob_querystart;
  } else if (splice_queryend_p == true) {
    new->splice_score = ambig_prob_queryend;
  } else {
    new->splice_score = 0.0;
  }
  
  new->distance = 0U;
  
  new->paired_usedp = false;
  new->concordantp = false;
  
  new->query_splicepos = -1;
  new->circularpos = compute_circularpos(&new->circularalias,new);
  
  debug0(printf("Method %s: Stage3end_new_substitution returning %p at %u..%u with nmatches_to_trims %d and amb length %d+%d\n",
		Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset,new->nmatches_to_trims,
		start_amb_length(new),end_amb_length(new)));
  
  /* Previously checked for (new->circularalias == +2 || new->circularalias == -2) */
  
  if (new->circularpos >= 0) {
    new->altlocp = false;
    return new;
    
  } else if ((new->altlocp = altlocp[chrnum]) == false) {
    return new;
    
  } else {
    return new;
  }
}



/* Previously allowed donor or acceptor to be NULL, when we performed Splice_group_by_segment */
/* Previously new->substring1 was donor and new->substring2 was acceptor */
/* TODO: Modify a Stage3end_new_splice to take two Stage3end_T parts, somewhat like a Stage3pair_T */
T
Stage3end_new_splice (int *found_score_overall, int *found_score_within_trims,
		      Substring_T donor, Substring_T acceptor,
		      Chrpos_T distance, bool shortdistancep, int querylength,
		      bool copy_donor_p, bool copy_acceptor_p, bool first_read_p, int orig_sensedir,
		      Listpool_T listpool, Method_T method, int level) {
  T new;
  Substring_T substring_for_concordance; /* always the inner substring */
  Substring_T substring_other;		 /* the outer substring */
  Substring_T substring, substring1, substringN;
  Junction_T junction;

  List_T transcripts;
  char *remap_sequence;
  int remap_seqlength;
  double donor_prob, acceptor_prob;
#ifdef DEBUG0
  List_T p;
#endif


  if (Substring_nmatches_to_trims(donor) < 15 || 
      Substring_nmatches_to_trims(acceptor) < 15) {
    /* Not enough evidence to find each end of the translocation */
    return (T) NULL;
  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
  }

  donor_prob = Substring_siteD_prob(donor);
  acceptor_prob = Substring_siteA_prob(acceptor);

  debug0(printf("Stage3end_new_splice, method %s: %p with first_read_p %d, sensedir %d, donor substring %p and acceptor substring %p, donor_prob %f and acceptor_prob %f\n",
		Method_string(method),new,first_read_p,orig_sensedir,donor,acceptor,donor_prob,acceptor_prob));

#if 0
  assert(Substring_match_length_orig(donor) + Substring_match_length_orig(acceptor) + amb_length == querylength);
#endif

  new->querylength_adj = new->querylength = querylength;

  new->nindels = 0;

  new->transcripts = (List_T) NULL;
  new->transcripts_other = (List_T) NULL;

  new->splice_score = donor_prob + acceptor_prob;

  new->method = method;
  new->level = level;

  if (shortdistancep == true) {
    new->distant_splice_p = false;

    new->hittype = SPLICE;
    new->genestrand = Substring_genestrand(donor);
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->chrhigh = Substring_chrhigh(donor);
    new->chrlength = Substring_chrlength(donor);

    assert(Substring_plusp(donor) == Substring_plusp(acceptor));
    assert(SENSE_CONSISTENT_P(Substring_sensedir(donor),Substring_sensedir(acceptor)));

  } else {
    new->distant_splice_p = true;

    if (Substring_chrnum(donor) == Substring_chrnum(acceptor) &&
	Substring_plusp(donor) == Substring_plusp(acceptor) &&
	SENSE_CONSISTENT_P(Substring_sensedir(donor),Substring_sensedir(acceptor))) {
      new->genestrand = Substring_genestrand(donor);
      new->hittype = SAMECHR_SPLICE;
      new->chrnum = Substring_chrnum(donor);
      new->chroffset = Substring_chroffset(donor);
      new->chrhigh = Substring_chrhigh(donor);
      new->chrlength = Substring_chrlength(donor);
    } else {
      new->hittype = TRANSLOC_SPLICE;
      new->genestrand = 0;
      new->chrnum = 0;
      new->chroffset = 0;
      new->chrhigh = 0;
      new->chrlength = 0;
    }
  }

  /* printf("Making splice with shortdistancep = %d, donor chrnum %d, and acceptor chrnum %d => chrnum %d\n",
     shortdistancep,Substring_chrnum(donor),Substring_chrnum(acceptor),new->chrnum); */

  new->distance = distance;
  new->guided_insertlength = 0U;
  new->nsegments = 2;
  new->nsplices = 1;

  /* Define substrings and junctions */
  if (new->chrnum != 0) {
    new->sensedir = orig_sensedir;
    junction = Junction_new_splice(distance,orig_sensedir,donor_prob,acceptor_prob);

  } else if (Substring_querystart(donor) < Substring_querystart(acceptor)) {
    /* Translocation, sense */
    new->sensedir = SENSE_FORWARD;
    junction = Junction_new_chimera(/*sensedir*/SENSE_FORWARD,donor_prob,acceptor_prob);

  } else {
    /* Translocation, antisense */
    new->sensedir = SENSE_ANTI;
    junction = Junction_new_chimera(/*sensedir*/SENSE_ANTI,donor_prob,acceptor_prob);
  }
  new->sensedir_for_concordance = new->sensedir;

  debug0(printf("donor querypos %d..%d\n",Substring_querystart(donor),Substring_queryend(donor)));
  debug0(printf("acceptor querypos %d..%d\n",Substring_querystart(acceptor),Substring_queryend(acceptor)));
  debug0(printf("sensedir %d\n",new->sensedir));


  new->junctions_LtoH = Listpool_push(NULL,listpool,(void *) junction);
  /* new->junctions_HtoL = Listpool_push(NULL,listpool,(void *) junction); */
  new->junctions_1toN = Listpool_push(NULL,listpool,(void *) junction);
  new->junctions_Nto1 = Listpool_push(NULL,listpool,(void *) junction);

  donor = copy_donor_p ? Substring_copy(donor) : donor;
  acceptor = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
  if (new->sensedir != SENSE_ANTI) {
    /* SENSE_FORWARD or SENSE_NULL */
    /* Order is donor (substring1), acceptor (substring2) */
    new->substrings_1toN = Listpool_push(NULL,listpool,(void *) acceptor);
    new->substrings_1toN = Listpool_push(new->substrings_1toN,listpool,(void *) donor);
  } else {
    /* SENSE_ANTI */
    /* Order is acceptor (substring1), donor (substring2) */
    new->substrings_1toN = Listpool_push(NULL,listpool,(void *) donor);
    new->substrings_1toN = Listpool_push(new->substrings_1toN,listpool,(void *) acceptor);
  }
  new->substrings_Nto1 = List_reverse(Listpool_copy(new->substrings_1toN,listpool));
  assert(Substring_querystart(List_head(new->substrings_1toN)) < Substring_querystart(List_head(new->substrings_Nto1)));
  /* Done assigning substrings */


  if (new->chrnum != 0) {
    /* Ordinary splice.  No need to distinguish effective_chrnum and other_chrnum */
    substring_for_concordance = substring_other = (Substring_T) NULL;
    new->effective_chrnum = new->chrnum;
    new->other_chrnum = 0;

    /* Define coordinates as usual */
    substring1 = (Substring_T) List_head(new->substrings_1toN);
    substringN = (Substring_T) List_head(new->substrings_Nto1);
    new->genomicstart = Substring_genomicstart(substring1);
    new->genomicend = Substring_genomicend(substringN);
    new->plusp = Substring_plusp(substring1);

  } else {
    /* Translocation.  Concordant substring is the inner one */
    if (first_read_p == true) {
      substring_for_concordance = (Substring_T) List_head(new->substrings_Nto1);
      substring_other = (Substring_T) List_head(new->substrings_1toN);
      debug0(printf("Since first read, substring for concordance is at chr %d\n",Substring_chrnum(substring_for_concordance)));
    } else {
      substring_for_concordance = (Substring_T) List_head(new->substrings_1toN);
      substring_other = (Substring_T) List_head(new->substrings_Nto1);
      debug0(printf("Since second read, substring for concordance is at chr %d\n",Substring_chrnum(substring_for_concordance)));
    }
      
    new->effective_chrnum = Substring_chrnum(substring_for_concordance);
    new->other_chrnum = Substring_chrnum(substring_other);

    /* Define coordinates based on substring for concordance */
    new->genomicstart = Substring_genomicstart(substring_for_concordance);
    new->genomicend = Substring_genomicend(substring_for_concordance);
    new->plusp = Substring_plusp(substring_for_concordance);
      

  }

#ifdef DEBUG0
  printf("NEW SUBSTRINGS (query order)\n");
  for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = List_head(p);
    if (Substring_ambiguous_p(substring) == true) {
      printf("%d..%d\t%d:%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\tprobs:%f and %f\n",
	     Substring_querystart(substring),Substring_queryend(substring),Substring_chrnum(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring),
	     Substring_amb_donor_prob(substring),Substring_amb_acceptor_prob(substring));
    } else {
      printf("%d..%d\t%d:%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\n",
	     Substring_querystart(substring),Substring_queryend(substring),Substring_chrnum(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring));
    }
  }
  printf("\n");
#endif


  /* This plusp is somewhat artificial, based on substring_for_concordance,
     but it defines order of substrings_LtoH */
  if (new->plusp == true) {
    new->substrings_LtoH = Listpool_copy(new->substrings_1toN,listpool);
    new->substrings_HtoL = Listpool_copy(new->substrings_Nto1,listpool);

    new->low = new->genomicstart;
    new->high = new->genomicend;

  } else {
    new->substrings_LtoH = Listpool_copy(new->substrings_Nto1,listpool);
    new->substrings_HtoL = Listpool_copy(new->substrings_1toN,listpool);

    new->low = new->genomicend;
    new->high = new->genomicstart;
  }
  new->genomiclength = new->high - new->low;

  debug0(printf("  hittype is %s, plusp %d, genomicpos %u..%u\n",
		hittype_string(new->hittype),new->plusp,new->genomicstart - new->chroffset,new->genomicend - new->chroffset));

  substring = (Substring_T) List_head(new->substrings_1toN);
  new->trim_querystart = Substring_trim_querystart(substring);
  new->trim_querystart_splicep = Substring_trim_querystart_splicep(substring);

  substring = (Substring_T) List_head(new->substrings_Nto1);
  new->trim_queryend = Substring_trim_queryend(substring);
  new->trim_queryend_splicep = Substring_trim_queryend_splicep(substring);
  

  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(donor) + Substring_nmismatches_bothdiff(acceptor);
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + Substring_nmismatches_refdiff(acceptor); */

  new->nmatches_to_trims = Substring_nmatches_to_trims(donor) + Substring_nmatches_to_trims(acceptor);
  new->nmatches_plus_spliced_trims = new->nmatches_to_trims;
  assert(new->nmatches_plus_spliced_trims <= querylength);

  new->score_within_trims = querylength - new->nmatches_plus_spliced_trims;
  new->score_overall = querylength - new->nmatches_to_trims;

  /* found_score_overall does not compensate for spliced ends, so gives motivation to find distant splicing */
  if (new->score_overall < *found_score_overall) {
    *found_score_overall = new->score_overall;
  }
  /* found_score_within_trims does compensate for spliced trims, and guides how much further alignment is necessary */
  if (new->score_within_trims < *found_score_within_trims) {
    *found_score_within_trims = new->score_within_trims;
  }

  debug0(printf("New splice has donor %d + acceptor %d matches, sensedir %d\n",
		Substring_nmatches(donor),Substring_nmatches(acceptor),new->sensedir));

  /* new->penalties = splicing_penalty; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

#if 0
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  new->paired_usedp = false;
  new->concordantp = false;

  if (new->sensedir != SENSE_ANTI) {
    assert(Substring_queryend(donor) == Substring_querystart(acceptor));
    new->query_splicepos = Substring_queryend(donor);
  } else {
    assert(Substring_queryend(acceptor) == Substring_querystart(donor));
    new->query_splicepos = Substring_queryend(acceptor);
  }
  assert(new->query_splicepos > 0 && new->query_splicepos < querylength - 1);


  new->circularpos = compute_circularpos(&new->circularalias,new);
  /* Previously checked for (new->circularalias == +2 || new->circularalias == -2) */

  if (new->circularpos >= 0) {
    new->altlocp = false;
  } else if ((new->altlocp = altlocp[new->chrnum]) == false) {
  } else {
  }
  
  if (transcriptomep == true && remap_transcriptome_p == true && substring_for_concordance != NULL) {
    /* Remap substring_for_concordance */
    remap_sequence = Substring_genomic_sequence(&remap_seqlength,substring_for_concordance,genomecomp);
    if ((transcripts = Kmer_remap_transcriptome(remap_sequence,remap_seqlength,new->effective_chrnum,
						Substring_chrpos_low(substring_for_concordance),
						Substring_chrpos_high(substring_for_concordance),
						transcript_iit,transcriptomebits,transcriptome)) != NULL) {
      new->transcripts = transcripts;
    }
    FREE(remap_sequence);
    
    /* Remap substring_other */
    remap_sequence = Substring_genomic_sequence(&remap_seqlength,substring_other,genomecomp);
    if ((transcripts = Kmer_remap_transcriptome(remap_sequence,remap_seqlength,new->other_chrnum,
						Substring_chrpos_low(substring_other),
						Substring_chrpos_high(substring_other),
						transcript_iit,transcriptomebits,transcriptome)) != NULL) {
      new->transcripts_other = transcripts;
    }
    FREE(remap_sequence);
  }
  
  debug0(printf("Returning new splice %p at genomic %u..%u, donor %p (%u => %u), acceptor %p (%u => %u), score %d\n",
		new,new->genomicstart - new->chroffset,new->genomicend - new->chroffset,donor,
		donor == NULL ? 0 : Substring_left_genomicseg(donor),
		donor == NULL ? 0 : Substring_splicecoord_D(donor),
		acceptor,acceptor == NULL ? 0 : Substring_left_genomicseg(acceptor),
		acceptor == NULL ? 0 : Substring_splicecoord_A(acceptor),new->score_within_trims));
  debug0(printf("sensedir %d\n",new->sensedir));
  return new;
}


T
Stage3end_new_distant (int *found_score_overall, int *found_score_within_trims,
		       Substring_T startfrag, Substring_T endfrag, int splice_pos,
		       int nmismatches1, int nmismatches2,
		       double prob1, double prob2, int sensedir_distant_guess,
		       Chrpos_T distance, bool shortdistancep, int querylength,
		       bool first_read_p, Listpool_T listpool, int level) {
  T new;
  Substring_T substring_for_concordance; /* always the inner substring */
  Substring_T substring_other;		 /* the outer substring */
  Substring_T substring;
  Substring_T donor, acceptor;
  Junction_T junction;

  List_T transcripts;
  char *remap_sequence;
  int remap_seqlength;
#ifdef DEBUG0
  List_T p;
#endif

  new = (T) MALLOC_OUT(sizeof(*new));

  debug0(printf("Stage3end_new_distant: %p with first_read_p %d\n",new,first_read_p));

  new->querylength_adj = new->querylength = querylength;

  new->nindels = 0;

  new->transcripts = (List_T) NULL;
  new->transcripts_other = (List_T) NULL;

  new->splice_score = 0.0;

  new->method = DISTANT_DNA;
  new->level = level;

  if (shortdistancep == true) {
    new->distant_splice_p = false;

    new->hittype = SPLICE;
    new->genestrand = Substring_genestrand(startfrag);
    new->chrnum = Substring_chrnum(startfrag);
    new->chroffset = Substring_chroffset(startfrag);
    new->chrhigh = Substring_chrhigh(startfrag);
    new->chrlength = Substring_chrlength(startfrag);

    assert(Substring_plusp(startfrag) == Substring_plusp(endfrag));
    assert(SENSE_CONSISTENT_P(Substring_sensedir(startfrag),Substring_sensedir(endfrag)));

  } else {
    new->distant_splice_p = true;

    new->hittype = TRANSLOC_SPLICE;
    new->genestrand = 0;
    new->chrnum = 0;
    new->chroffset = 0;
    new->chrhigh = 0;
    new->chrlength = 0;
  }

  /* printf("Making splice with shortdistancep = %d, startfrag chrnum %d, and endfrag chrnum %d => chrnum %d\n",
     shortdistancep,Substring_chrnum(startfrag),Substring_chrnum(endfrag),new->chrnum); */

  new->distance = distance;
  new->guided_insertlength = 0U;
  new->nsegments = 2;
  new->nsplices = 1;

  /* Trim startfrag and endfrag at splice_pos */
  startfrag = Substring_trim_startfrag(nmismatches1,/*old*/startfrag,/*new_queryend*/splice_pos);
  endfrag = Substring_trim_endfrag(nmismatches2,/*old*/endfrag,/*new_querystart*/splice_pos);

  /* Define substrings and junctions */
  new->sensedir_for_concordance = SENSE_NULL; /* Because we are uncertain about splice */
  new->sensedir = sensedir_distant_guess;
  if (sensedir_distant_guess != SENSE_ANTI) {
    /* Order is donor (substring1), acceptor (substring2) */
    donor = startfrag;
    Substring_label_donor(donor,splice_pos,prob1,sensedir_distant_guess);

    acceptor = endfrag;
    Substring_label_acceptor(acceptor,splice_pos,prob2,sensedir_distant_guess);

    new->substrings_1toN = Listpool_push(NULL,listpool,(void *) acceptor);
    new->substrings_1toN = Listpool_push(new->substrings_1toN,listpool,(void *) donor);

  } else {
    /* Order is acceptor (substring1), donor (substring2) */
    acceptor = startfrag;
    Substring_label_acceptor(acceptor,splice_pos,prob1,sensedir_distant_guess);

    donor = endfrag;
    Substring_label_donor(donor,splice_pos,prob2,sensedir_distant_guess);

    new->substrings_1toN = Listpool_push(NULL,listpool,(void *) donor);
    new->substrings_1toN = Listpool_push(new->substrings_1toN,listpool,(void *) acceptor);
  }

  if (shortdistancep == true) {
    junction = Junction_new_splice(distance,sensedir_distant_guess,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor));
  } else {
    junction = Junction_new_chimera(sensedir_distant_guess,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor));
  }

  new->junctions_LtoH = Listpool_push(NULL,listpool,(void *) junction);
  /* new->junctions_HtoL = Listpool_push(NULL,listpool,(void *) junction); */
  new->junctions_1toN = Listpool_push(NULL,listpool,(void *) junction);
  new->junctions_Nto1 = Listpool_push(NULL,listpool,(void *) junction);

  new->substrings_Nto1 = List_reverse(Listpool_copy(new->substrings_1toN,listpool));
  assert(Substring_querystart(List_head(new->substrings_1toN)) < Substring_querystart(List_head(new->substrings_Nto1)));
  /* Done assigning substrings */


  /* Translocation.  Concordant substring is the inner one */
  if (first_read_p == true) {
    substring_for_concordance = (Substring_T) List_head(new->substrings_Nto1);
    substring_other = (Substring_T) List_head(new->substrings_1toN);
    debug0(printf("Since first read, substring for concordance is at chr %d\n",Substring_chrnum(substring_for_concordance)));
  } else {
    substring_for_concordance = (Substring_T) List_head(new->substrings_1toN);
    substring_other = (Substring_T) List_head(new->substrings_Nto1);
    debug0(printf("Since second read, substring for concordance is at chr %d\n",Substring_chrnum(substring_for_concordance)));
  }
      
  new->effective_chrnum = Substring_chrnum(substring_for_concordance);
  new->other_chrnum = Substring_chrnum(substring_other);

  /* Define coordinates based on substring for concordance */
  new->genomicstart = Substring_genomicstart(substring_for_concordance);
  new->genomicend = Substring_genomicend(substring_for_concordance);
  new->plusp = Substring_plusp(substring_for_concordance);
      
#ifdef DEBUG0
  printf("NEW SUBSTRINGS (query order)\n");
  for (p = new->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = List_head(p);
    if (Substring_ambiguous_p(substring) == true) {
      printf("%d..%d\t%d:%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\tprobs:%f and %f\n",
	     Substring_querystart(substring),Substring_queryend(substring),Substring_chrnum(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring),
	     Substring_amb_donor_prob(substring),Substring_amb_acceptor_prob(substring));
    } else {
      printf("%d..%d\t%d:%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\n",
	     Substring_querystart(substring),Substring_queryend(substring),Substring_chrnum(substring),
	     Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	     Substring_nmatches_to_trims(substring),Substring_amb_length(substring));
    }
  }
  printf("\n");
#endif


  /* This plusp is somewhat artificial, based on substring_for_concordance,
     but it defines order of substrings_LtoH */
  if (new->plusp == true) {
    new->substrings_LtoH = Listpool_copy(new->substrings_1toN,listpool);
    new->substrings_HtoL = Listpool_copy(new->substrings_Nto1,listpool);

    new->low = new->genomicstart;
    new->high = new->genomicend;

  } else {
    new->substrings_LtoH = Listpool_copy(new->substrings_Nto1,listpool);
    new->substrings_HtoL = Listpool_copy(new->substrings_1toN,listpool);

    new->low = new->genomicend;
    new->high = new->genomicstart;
  }
  new->genomiclength = new->high - new->low;

  debug0(printf("  hittype is %s, plusp %d, genomicpos %u..%u\n",
		hittype_string(new->hittype),new->plusp,new->genomicstart - new->chroffset,new->genomicend - new->chroffset));

  substring = (Substring_T) List_head(new->substrings_1toN);
  new->trim_querystart = Substring_trim_querystart(substring);
  new->trim_querystart_splicep = Substring_trim_querystart_splicep(substring);

  substring = (Substring_T) List_head(new->substrings_Nto1);
  new->trim_queryend = Substring_trim_queryend(substring);
  new->trim_queryend_splicep = Substring_trim_queryend_splicep(substring);
  

  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(startfrag) + Substring_nmismatches_bothdiff(endfrag);
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(startfrag) + Substring_nmismatches_refdiff(endfrag); */

  new->nmatches_to_trims = Substring_nmatches_to_trims(startfrag) + Substring_nmatches_to_trims(endfrag);
  new->nmatches_plus_spliced_trims = new->nmatches_to_trims;
  assert(new->nmatches_plus_spliced_trims <= querylength);

  new->score_within_trims = querylength - new->nmatches_plus_spliced_trims;
  new->score_overall = querylength - new->nmatches_to_trims;

  /* found_score_overall does not compensate for spliced ends, so gives motivation to find distant splicing */
  if (new->score_overall < *found_score_overall) {
    *found_score_overall = new->score_overall;
  }
  /* found_score_within_trims does compensate for spliced trims, and guides how much further alignment is necessary */
  if (new->score_within_trims < *found_score_within_trims) {
    *found_score_within_trims = new->score_within_trims;
  }

  debug0(printf("New distant has startfrag %d + endfrag %d matches, sensedir %d\n",
		Substring_nmatches(startfrag),Substring_nmatches(endfrag),new->sensedir));

  /* new->penalties = splicing_penalty; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

#if 0
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  new->paired_usedp = false;
  new->concordantp = false;
  new->query_splicepos = splice_pos;

  new->circularpos = compute_circularpos(&new->circularalias,new);
  /* Previously checked for (new->circularalias == +2 || new->circularalias == -2) */

  if (new->circularpos >= 0) {
    new->altlocp = false;
  } else if ((new->altlocp = altlocp[new->chrnum]) == false) {
  } else {
  }
  
  if (transcriptomep == true && remap_transcriptome_p == true && substring_for_concordance != NULL) {
    /* Remap substring_for_concordance */
    remap_sequence = Substring_genomic_sequence(&remap_seqlength,substring_for_concordance,genomecomp);
    if ((transcripts = Kmer_remap_transcriptome(remap_sequence,remap_seqlength,new->effective_chrnum,
						Substring_chrpos_low(substring_for_concordance),
						Substring_chrpos_high(substring_for_concordance),
						transcript_iit,transcriptomebits,transcriptome)) != NULL) {
      new->transcripts = transcripts;
    }
    FREE(remap_sequence);
    
    /* Remap substring_other */
    remap_sequence = Substring_genomic_sequence(&remap_seqlength,substring_other,genomecomp);
    if ((transcripts = Kmer_remap_transcriptome(remap_sequence,remap_seqlength,new->other_chrnum,
						Substring_chrpos_low(substring_other),
						Substring_chrpos_high(substring_other),
						transcript_iit,transcriptomebits,transcriptome)) != NULL) {
      new->transcripts_other = transcripts;
    }
    FREE(remap_sequence);
  }
  
  debug0(printf("Returning new distant %p at genomic %u..%u, startfrag %p (%u => ), endfrag %p (%u => ), score %d\n",
		new,new->genomicstart - new->chroffset,new->genomicend - new->chroffset,
		startfrag,Substring_left_genomicseg(startfrag),endfrag,Substring_left_genomicseg(endfrag),
		new->score_within_trims));
  return new;
}


#if 0
T
Stage3end_new_gmap (int *found_score, int nmatches_to_trims, int max_match_length,
		    int ambig_end_length_5, int ambig_end_length_3,
		    double avg_splice_score, struct Simplepair_T *pairarray, int npairs,
		    Univcoord_T left, bool plusp, int genestrand,
		    int querylength, Compress_T query_compress,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		    int cdna_direction, int orig_sensedir, Listpool_T listpool, Method_T method, int level) {
  List_T substrings, junctions;
#ifdef DEBUG0
  Substring_T substring;
  Junction_T junction;
  List_T p;
#endif

  T new;
  /* double prob5_sense_forward, prob5_sense_anti, prob3_sense_forward, prob3_sense_anti; */
  int nmismatches;

  Univcoord_T new_chroffset;
  Chrpos_T adjustment;
  Univcoord_T alignstart, alignend;
  Chrpos_T chrpos_start, chrpos_end;
  int outofbounds_start, outofbounds_end;
  int trimmed_start, trimmed_length, n1, n2;


  debug0(printf("Entered Stage3end_new_gmap, method %s, with max_match_length %d, orig_sensedir %d, ambig_end_length_5 %d, ambig_end_length_3 %d, and method %s\n",
		Method_string(method),max_match_length,orig_sensedir,ambig_end_length_5,ambig_end_length_3,Method_string(method)));
  assert(orig_sensedir == SENSE_NULL || orig_sensedir == SENSE_ANTI || orig_sensedir == SENSE_FORWARD);

  if (max_match_length < gmap_min_nconsecutive) {
    debug0(printf("  Bad GMAP: max_match_length %d < %d, so returning NULL\n",max_match_length,gmap_min_nconsecutive));
    FREE_OUT(pairarray);
    return (T) NULL;
  }

  /* Compute chromosomal bounds and trimming */
  trimmed_start = 0;
  trimmed_length = npairs;

  /* Use original chroffset to compute alignstart and alignend */
  Pairarray_chrpos_bounds(&chrpos_start,&chrpos_end,pairarray,npairs);
  alignstart = chroffset + chrpos_start;
  alignend = chroffset + chrpos_end;
  debug11(printf("alignstart %u, alignend %u\n",alignstart,alignend));

  /* Logic modeled after Substring_new */
  new_chroffset = chroffset;
  if (plusp) {
    if (alignstart >= chrhigh) {
      /* Need to recompute chromosome bounds */
      debug11(printf("Plus: recomputing chromosomal bounds\n"));
      chrnum = Univ_IIT_get_one(chromosome_iit,alignstart,alignstart);
      Univ_IIT_interval_bounds(&new_chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
    }

    if (alignend <= chrhigh) {
      /* Alignment is within the chromosome */
      debug11(printf("alignend <= chrhigh, so alignment is within the chromosome\n"));
      outofbounds_start = outofbounds_end = 0;
    } else {
      /* Alignment straddles the high bound of the chromosome */
      outofbounds_start = chrhigh - alignstart;
      outofbounds_end = alignend - chrhigh;
      if (outofbounds_start >= outofbounds_end) {
	/* Keep in existing chromosome */
	outofbounds_start = 0;
      } else if (++chrnum > nchromosomes) {
	/* Went past the end of the genome */
	FREE_OUT(pairarray);
	return (T) NULL;
      } else {
	/* Move to next chromosome */
	debug11(printf("Moving to next chromosome\n"));
	Univ_IIT_interval_bounds(&new_chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	if (alignend <= chrhigh) {
	  /* Alignment is within the new chromosome */
	  debug11(printf("alignend <= new chrhigh %u, so alignment is within the new chromosome\n",chrhigh));
	  outofbounds_end = 0;
	} else {
	  outofbounds_end = alignend - chrhigh;
	}
      }
    }

#ifdef DEBUG11
    if (outofbounds_start != 0 || outofbounds_end != 0) {
      printf("outofbounds_start %d, outofbounds_end %d\n",outofbounds_start,outofbounds_end);
      printf("Original plus\n");
      Pair_dump_array(pairarray,npairs,true);
    }
#endif

    if (outofbounds_start == 0) {
      n1 = npairs;
    } else {
      debug11(printf("Counting ge from end for chrbound %u = %u - %u\n",
		     new_chroffset - chroffset,new_chroffset,chroffset));
      n1 = Pair_count_ge_fromend(pairarray,npairs,/*chrbound*/new_chroffset - chroffset);
    }
    if (outofbounds_end == 0) {
      n2 = npairs;
    } else {
      debug11(printf("Counting lt from start for chrbound %u = %u - %u\n",
		     chrhigh - chroffset,chrhigh,chroffset));
      n2 = Pair_count_lt_fromstart(pairarray,npairs,/*chrbound*/chrhigh - chroffset);
    }

    if ((trimmed_length = n1 + n2 - npairs) <= 0) {
      debug11(printf("Everything got trimmed! (unusual)\n"));
      FREE_OUT(pairarray);
      return (T) NULL;
    } else {
      trimmed_start = npairs - n1;
    }

  } else {
    if (alignend >= chrhigh) {
      /* Need to recompute chromosome bounds */
      debug11(printf("Minus: Recomputing chromosomal bounds\n"));
      chrnum = Univ_IIT_get_one(chromosome_iit,alignend,alignend);
      Univ_IIT_interval_bounds(&new_chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
    }

    if (alignstart <= chrhigh) {
      /* Alignment is within the chromosome */
      debug11(printf("alignstart <= chrhigh, so alignment is within the chromosome\n"));
      outofbounds_start = outofbounds_end = 0;
    } else {
      /* Alignment straddles the high bound of the chromosome */
      outofbounds_start = alignstart - chrhigh;
      outofbounds_end = chrhigh - alignend;
      if (outofbounds_end >= outofbounds_start) {
	/* Keep in existing chromosome */
	outofbounds_end = 0;
      } else if (++chrnum > nchromosomes) {
	/* Went past the end of the genome */
	FREE_OUT(pairarray);
	return (T) NULL;
      } else {
	/* Move to next chromosome */
	debug11(printf("Moving to next chromosome\n"));
	Univ_IIT_interval_bounds(&new_chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	if (alignstart <= chrhigh) {
	  debug11(printf("alignstart <= new chrhigh %u, so alignment is within the new chromosome\n",chrhigh));
	  outofbounds_start = 0;
	} else {
	  outofbounds_start = alignstart - chrhigh;
	}
      }
    }

#ifdef DEBUG11
    if (outofbounds_start != 0 || outofbounds_end != 0) {
      printf("outofbounds_start %d, outofbounds_end %d\n",outofbounds_start,outofbounds_end);
      printf("Original minus\n");
      Pair_dump_array(pairarray,npairs,true);
    }
#endif

    if (outofbounds_start == 0) {
      n1 = npairs;
    } else {
      debug11(printf("Counting lt from end for chrbound %u = %u - %u\n",
		     chrhigh - chroffset,chrhigh,chroffset));
      n1 = Pair_count_lt_fromend(pairarray,npairs,/*chrbound*/chrhigh - chroffset);
    }
    if (outofbounds_end == 0) {
      n2 = npairs;
    } else {
      debug11(printf("Counting ge from start for chrbound %u = %u - %u\n",
		     new_chroffset - chroffset,new_chroffset,chroffset));
      n2 = Pair_count_ge_fromstart(pairarray,npairs,/*chrbound*/new_chroffset - chroffset);
    }

    if ((trimmed_length = n1 + n2 - npairs) <= 0) {
      debug11(printf("Everything got trimmed! (unusual)\n"));
      FREE_OUT(pairarray);
      return (T) NULL;
    } else {
      trimmed_start = npairs - n1;
    }
  }

  debug11(printf("n1 %d, n2 %d, npairs %d\n",n1,n2,npairs));
  if ((adjustment = new_chroffset - chroffset) > 0) {
    debug11(printf("Adjustment %u\n",adjustment));
    Pair_subtract_genomepos(pairarray,npairs,adjustment);
  }

  debug0(printf("Stage3end_new_gmap: left %llu, chrhigh %llu, chrnum %d, nmatches_to_trims %d, cdna_direction %d, orig_sensedir %d, avg_splice_score %f, max_match_length %d\n",
		(unsigned long long) left,(unsigned long long) chrhigh,chrnum,nmatches_to_trims,cdna_direction,orig_sensedir,avg_splice_score,max_match_length));
  debug0(printf("  ambig_end_length_5 %d, ambig_end_length_3 %d\n",ambig_end_length_5,ambig_end_length_3));
  debug0(Pair_dump_genome_array(pairarray,npairs));
  debug0(Pair_dump_comp_array(pairarray,npairs));

  if ((substrings = Pairarray_convert_to_substrings(&junctions,&nmismatches,&(pairarray[trimmed_start]),trimmed_length,
						    querylength,/*watsonp*/plusp,genestrand,query_compress,orig_sensedir,
						    chrnum,chroffset,chrhigh,chrlength,listpool)) == NULL) {
    FREE_OUT(pairarray);
    return (T) NULL;

  } else {
#ifdef DEBUG0
    printf("nmismatches %d\n",nmismatches);
    printf("NEW SUBSTRINGS FROM GMAP (query order)\n");
    for (p = substrings; p != NULL; p = List_next(p)) {
      substring = List_head(p);
      if (Substring_ambiguous_p(substring) == true) {
	printf("%d..%d\t%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\tprobs:%f and %f\n",Substring_querystart(substring),Substring_queryend(substring),
	       Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	       Substring_nmatches_to_trims(substring),Substring_amb_length(substring),
	       Substring_amb_donor_prob(substring),Substring_amb_acceptor_prob(substring));
      } else {
	printf("%d..%d\t%u..%u\tmismatches:%d\tmatches_to_trims:%d\tamb:%d\n",Substring_querystart(substring),Substring_queryend(substring),
	       Substring_alignstart_trim_chr(substring),Substring_alignend_trim_chr(substring),Substring_nmismatches_bothdiff(substring),
	       Substring_nmatches_to_trims(substring),Substring_amb_length(substring));
      }
    }
    printf("\n");

    printf("NEW JUNCTIONS FROM GMAP (query order)\n");
    for (p = junctions; p != NULL; p = List_next(p)) {
      junction = List_head(p);
      printf("splice distance %u, nindels %d, probs %f and %f\n",
	     Junction_splice_distance(junction),Junction_nindels(junction),
	     Junction_donor_prob(junction),Junction_acceptor_prob(junction));
    }
    printf("\n");
#endif
  
    FREE_OUT(pairarray);
  }

  /* Modified from resolve_end_diagonals in path-solve.c */
#if 0
  /* From Dynprog_end procedures, need alts_coords, knowni, nmismatchesi, probi, and probj */
  if (ambig_end_length_5 > 0) {
    if (plusp == true) {
      if (orig_sensedir != SENSE_ANTI) {
      }
    }
  }
#endif


  new = Stage3end_new_precomputed(&(*found_score),/*nmismatches_bothdiff*/nmismatches,
				  substrings,junctions,/*transcripts*/NULL,/*transcripts_other*/NULL,
				  querylength,chrnum,chroffset,chrhigh,chrlength,
				  plusp,genestrand,orig_sensedir,listpool,method,level);

  if (new->nmatches_to_trims < querylength/2) {
    debug0(printf("  Bad GMAP: nmatches %d < querylength %d/2, so returning NULL\n",
		  new->nmatches_to_trims,querylength));
    Stage3end_free(&new);
    return (T) NULL;
  } else {
    debug0(printf("Method %s: Stage3end_new_gmap returning %p at %u..%u with splice score %f, nmatches %d\n",
		  Method_string(method),new,new->genomicstart - chroffset,new->genomicend - chroffset,new->splice_score,new->nmatches));
    return new;
  }
}
#endif



static int
Stage3end_output_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->nmatches_plus_spliced_trims > y->nmatches_plus_spliced_trims) {
    return -1;
  } else if (y->nmatches_plus_spliced_trims > x->nmatches_plus_spliced_trims) {
    return +1;
  } else if (x->mapq_loglik > y->mapq_loglik) {
    return -1;
  } else if (y->mapq_loglik > x->mapq_loglik) {
    return +1;
  } else if (x->distant_splice_p == false && y->distant_splice_p == true) {
    return -1;
  } else if (y->distant_splice_p == false && x->distant_splice_p == true) {
    return +1;
  } else if (x->guided_insertlength > 0 && y->guided_insertlength == 0) {
    return -1;
  } else if (y->guided_insertlength > 0 && x->guided_insertlength == 0) {
    return +1;
  } else if (x->guided_insertlength < y->guided_insertlength) {
    return -1;
  } else if (y->guided_insertlength < x->guided_insertlength) {
    return +1;
#if 0
  } else if (x->nmatches_to_trims > y->nmatches_to_trims) {
    return -1;
  } else if (y->nmatches_to_trims > x->nmatches_to_trims) {
    return +1;
#endif
  } else if (x->score_within_trims < y->score_within_trims) {
    return -1;
  } else if (y->score_within_trims < x->score_within_trims) {
    return +1;

    /* This genomic ordering will be undone if want_random_p is true */
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (y->genomicstart < x->genomicstart) {
    return +1;

  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (y->genomicend < x->genomicend) {
    return +1;

  } else if (x->plusp == true && y->plusp == false) {
    return -1;
  } else if (x->plusp == false && y->plusp == true) {
    return +1;

  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (y->hittype < x->hittype) {
    return +1;

  } else {
    return 0;
  }
}


static int
Stage3pair_output_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

#ifdef USE_BINGO
  if (x->absdifflength_bingo_p == true && y->absdifflength_bingo_p == false) {
    return -1;
  } else if (y->absdifflength_bingo_p == true && x->absdifflength_bingo_p == false) {
    return +1;
  }
#endif

  if (x->nmatches_plus_spliced_trims > y->nmatches_plus_spliced_trims) {
    return -1;
  } else if (y->nmatches_plus_spliced_trims > x->nmatches_plus_spliced_trims) {
    return +1;
  } else if (x->mapq_loglik > y->mapq_loglik) {
    return -1;
  } else if (y->mapq_loglik > x->mapq_loglik) {
    return +1;
  } else if (x->insertlength > 0 && y->insertlength == 0) {
    return -1;
  } else if (y->insertlength > 0 && x->insertlength == 0) {
    return +1;
  } else if (x->insertlength < y->insertlength) {
    return -1;
  } else if (y->insertlength < x->insertlength) {
    return +1;
#if 0
  } else if (x->nmatches_to_trims > y->nmatches_to_trims) {
    return -1;
  } else if (y->nmatches_to_trims > x->nmatches_to_trims) {
    return +1;
#endif
  } else if (x->score_within_trims < y->score_within_trims) {
    return -1;
  } else if (y->score_within_trims < x->score_within_trims) {
    return +1;

    /* This genomic ordering will be undone if want_random_p is true */
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;

  } else if (x->high < y->high) {
    return -1;
  } else if (y->high < x->high) {
    return +1;

  } else {
    return 0;
  }
}



static float
Stage3end_compute_mapq (Stage3end_T this, char *quality_string) {
  List_T p;
  Substring_T substring;

  if (this == NULL) {
    return 0.0;

  } else {
    this->mapq_loglik = 0.0;
    for (p = this->substrings_LtoH; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      this->mapq_loglik += Substring_compute_mapq(substring,quality_string);
    }
  }

  return this->mapq_loglik;
}



static void
Stage3end_display_prep (Stage3end_T this, char *queryuc_ptr, int alts_resolve, bool first_read_p) {
  List_T p, q;
  Substring_T substring, anchor;
  Junction_T pre_junction, post_junction, junction;
  Junctiontype_T type;
  int extraleft, extraright;
  Univcoord_T left, ignore;
  double donor_prob, acceptor_prob;
  bool sam_print_xt_p = false;
  /* int type; */
  /* int extralow, extrahigh; */

  if (this != NULL) {
    if (output_type == SAM_OUTPUT) {
      if (this->hittype == TRANSLOC_SPLICE ||
	  (this->hittype == SAMECHR_SPLICE && merge_samechr_p == false)) {
	  /* This is the condition in samprint to print the XT field, which needs the splice information */
	sam_print_xt_p = true;
      }
    }

    debug0(printf("Doing a display prep of end %p\n",this));
    /* Resolve ambiguous end */
    if (alts_resolve >= 0) {
      if (first_read_p == true) {
	substring = (Substring_T) List_head(this->substrings_Nto1);
	anchor = (Substring_T) List_head(List_next(this->substrings_Nto1));
	junction = (Junction_T) List_head(this->junctions_Nto1);
	left = Substring_set_alt(&donor_prob,&acceptor_prob,&ignore,&this->genomicend,substring,alts_resolve);
	if (this->plusp == true) {
	  Junction_set_unambiguous(junction,left - Substring_left(anchor),donor_prob,acceptor_prob);
	} else {
	  Junction_set_unambiguous(junction,Substring_left(anchor) - left,donor_prob,acceptor_prob);
	}

      } else {
	substring = (Substring_T) List_head(this->substrings_1toN);
	anchor = (Substring_T) List_head(List_next(this->substrings_1toN));
	junction = (Junction_T) List_head(this->junctions_1toN);
	left = Substring_set_alt(&donor_prob,&acceptor_prob,&this->genomicstart,&ignore,substring,alts_resolve);
	if (this->plusp == true) {
	  Junction_set_unambiguous(junction,Substring_left(anchor) - left,donor_prob,acceptor_prob);
	} else {
	  Junction_set_unambiguous(junction,left - Substring_left(anchor),donor_prob,acceptor_prob);
	}
      }
    }

    this->nmismatches_refdiff = 0;

    /* First segments */
    /* For operations on substrings, proceed in 1toN order, not LtoH order */
    substring = (Substring_T) List_head(this->substrings_1toN);
    if (output_type == STD_OUTPUT) {
      extraleft = Substring_querystart(substring); /* terminal start */
    } else {
      extraleft = 0;
    }

    if (List_length(this->substrings_1toN) == 1) {
      post_junction = (Junction_T) NULL;
      if (output_type == STD_OUTPUT) {
	extraright = this->querylength - Substring_queryend(substring); /* terminal end */
      } else {
	extraright = 0;
      }
    } else {
      post_junction = (Junction_T) List_head(this->junctions_1toN);
      /* Junction_print(post_junction); */

      if (output_type == M8_OUTPUT) {
	extraright = 0;
      } else if ((type = Junction_type(post_junction)) == CHIMERA_JUNCTION || sam_print_xt_p == true) {
	extraright = 2;
      } else if (output_type == SAM_OUTPUT) {
	extraright = 0;
      } else if (type == SPLICE_JUNCTION) {
	extraright = 2;
      } else if (first_read_p == true && type == DEL_JUNCTION) {
	extraright = Junction_nindels(post_junction);
      } else {
	extraright = 0;
      }
    }
      
    if (Substring_has_alts_p(substring) == true) {
      /* Skip */
    } else {
      this->nmismatches_refdiff += 
	Substring_display_prep(substring,queryuc_ptr,this->querylength,
			       extraleft,extraright,genomecomp);
    }

    assert(List_length(this->substrings_1toN) == List_length(this->junctions_1toN) + 1);
    if ((p = List_next(this->substrings_1toN)) == NULL) {
      /* No middle segments */
    } else {
      for (q = List_next(this->junctions_1toN); q != NULL; p = List_next(p), q = List_next(q)) {
	/* Middle segments */
	pre_junction = post_junction;
	post_junction = List_head(q);

	/* Junction_print(pre_junction); */
	/* Junction_print(post_junction); */

	if (output_type == M8_OUTPUT) {
	  extraleft = 0;
	} else if ((type = Junction_type(pre_junction)) == CHIMERA_JUNCTION || sam_print_xt_p == true) {
	  extraleft = 2;
	} else if (output_type == SAM_OUTPUT) {
	  extraleft = 0;
	} else if (type == SPLICE_JUNCTION) {
	  extraleft = 2;
	} else if (first_read_p == false && type == DEL_JUNCTION) {
	  extraleft = Junction_nindels(pre_junction);
	} else {
	  extraleft = 0;
	}

	if (output_type == M8_OUTPUT) {
	  extraright = 0;
	} else if ((type = Junction_type(post_junction)) == CHIMERA_JUNCTION || sam_print_xt_p == true) {
	  extraright = 2;
	} else if (output_type == SAM_OUTPUT) {
	  extraright = 0;
	} else if (type == SPLICE_JUNCTION) {
	  extraright = 2;
	} else if (first_read_p == true && type == DEL_JUNCTION) {
	  extraright = Junction_nindels(post_junction);
	} else {
	  extraright = 0;
	}

	substring = (Substring_T) List_head(p);
	if (Substring_has_alts_p(substring) == true) {
	  /* Skip */
	} else {
	  this->nmismatches_refdiff += 
	    Substring_display_prep(substring,queryuc_ptr,this->querylength,
				   extraleft,extraright,genomecomp);
	}
      }

      /* Last segment */
      pre_junction = post_junction;
      /* Junction_print(pre_junction); */

      if (output_type == M8_OUTPUT) {
	extraleft = 0;
      } else if ((type = Junction_type(pre_junction)) == CHIMERA_JUNCTION || sam_print_xt_p == true) {
	extraleft = 2;
      } else if (output_type == SAM_OUTPUT) {
	extraleft = 0;
      } else if (type == SPLICE_JUNCTION) {
	extraleft = 2;
      } else if (first_read_p == false && type == DEL_JUNCTION) {
	extraleft = Junction_nindels(pre_junction);
      } else {
	extraleft = 0;
      }

      substring = (Substring_T) List_head(p);
      if (output_type == STD_OUTPUT) {
	extraright = this->querylength - Substring_queryend(substring);
      } else {
	extraright = 0;
      }
	
      if (Substring_has_alts_p(substring) == true) {
	/* Skip */
      } else {
	this->nmismatches_refdiff += 
	  Substring_display_prep(substring,queryuc_ptr,this->querylength,
				 extraleft,extraright,genomecomp);
      }
    }
  }

  return;
}


List_T
Stage3end_filter (List_T hits, Hitlistpool_T hitlistpool, int max_mismatches, int min_coverage) {
  List_T newhits = NULL, p;
  Stage3end_T hit;

  if (ignore_trim_p == false) {
    for (p = hits; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
#if 0    
      printf("Comparing score %d against max_mismatches %d, and coverage %d against min_coverage %d\n",
	     hit->score,max_mismatches,hit->querylength - hit->trim_querystart - hit->trim_queryend,min_coverage);
#endif
      if (hit->score_overall > max_mismatches) {
	Stage3end_free(&hit);
      } else if (hit->querylength - hit->trim_querystart - hit->trim_queryend < min_coverage) {
	Stage3end_free(&hit);
      } else {
	newhits = Hitlist_push(newhits,hitlistpool,(void *) hit);
      }
    }

  } else {
    for (p = hits; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
#if 0    
      printf("Comparing score_ignore_trim %d against max_mismatches %d, and coverage %d against min_coverage %d\n",
	     hit->score_ignore_trim,max_mismatches,hit->querylength - hit->trim_querystart - hit->trim_queryend,min_coverage);
#endif
      if (hit->score_within_trims > max_mismatches) {
	Stage3end_free(&hit);
      } else if (hit->querylength - hit->trim_querystart - hit->trim_queryend < min_coverage) {
	Stage3end_free(&hit);
      } else {
	newhits = Hitlist_push(newhits,hitlistpool,(void *) hit);
      }
    }
  }

  Hitlist_free(&hits);
  return newhits;
}




Stage3end_T *
Stage3end_eval_and_sort (int npaths, int *first_absmq, int *second_absmq,
			 Stage3end_T *stage3array, char *queryuc_ptr,
			 char *quality_string, bool displayp) {
  float maxlik, loglik;
  float total, q;		/* For Bayesian mapq calculation */
  int compute_npaths;

  int randomi, i;
  Stage3end_T temp, hit;

  if (npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (npaths == 1) {
    hit = stage3array[0];
    hit->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    hit->mapq_score = MAPQ_max_quality_score(quality_string,hit->querylength);
    hit->absmq_score = MAPQ_MAXIMUM_SCORE;

    if (displayp == true) {
      Stage3end_display_prep(hit,queryuc_ptr,/*alts_resolve*/-1,/*first_read_p*/true);
    }
    *first_absmq = hit->absmq_score;
    *second_absmq = 0;

  } else {
    /* Compute mapq_loglik */
    for (i = 0; i < npaths; i++) {
      Stage3end_compute_mapq(stage3array[i],quality_string);
    }

    /* Sort by nmatches, then mapq */
    qsort(stage3array,npaths,sizeof(Stage3end_T),Stage3end_output_cmp);

    if (want_random_p) {
      /* Randomize among best alignments */
      i = 1;
      while (i < npaths && Stage3end_output_cmp(&(stage3array[i]),&(stage3array[0])) == 0) {
	i++;
      }
      if (i > 1) {		/* i is number of ties */
	/* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
	randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) i));
	/* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
	temp = stage3array[0];
	stage3array[0] = stage3array[randomi];
	stage3array[randomi] = temp;
      }
    }

    /* Enforce monotonicity */
    for (i = npaths - 1; i > 0; i--) {
      if (stage3array[i-1]->mapq_loglik < stage3array[i]->mapq_loglik) {
	stage3array[i-1]->mapq_loglik = stage3array[i]->mapq_loglik;
      }
    }
    maxlik = stage3array[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < npaths; i++) {
      stage3array[i]->mapq_loglik -= maxlik;
    }

#if 0
    /* Save on computation if possible */
    /* Not possible, since we are going to select randomly from among all npaths */
    if (npaths < maxpaths) {
      compute_npaths = npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }
#else
    compute_npaths = npaths;
#endif

    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3array[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3array[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = stage3array[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < npaths; i++) {
      total += (stage3array[i]->mapq_loglik = fasterexp(stage3array[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3array[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3array[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3array[i]->mapq_score = 96;
      } else {
	stage3array[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

    if (displayp == true) {
      /* Prepare for display */
      for (i = 0; i < compute_npaths; i++) {
	Stage3end_display_prep(stage3array[i],queryuc_ptr,/*alts_resolve*/-1,/*first_read_p*/true);
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3array[0]->mapq_score >= mapq_unique_score &&
	stage3array[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3end_free(&(stage3array[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3array;
}


static int
insertlength_expected (Chrpos_T insertlength) {
  if (insertlength < expected_pairlength_low) {
    return -1;
  } else if (insertlength > expected_pairlength_very_high) {
    return -1;
  } else if (insertlength > expected_pairlength_high) {
    return 0;
  } else {
    return +1;
  }
}


/* For concordant ends */
static Chrpos_T
pair_insert_length (int *pair_relationship, Stage3end_T hit5, Stage3end_T hit3) {
  List_T p, q;
  Substring_T substring5, substring3;

  if (hit5->plusp != hit3->plusp) {
    debug10(printf("pair_insert_length: hit5->plusp %d != hit3->plusp %d, so returning 0\n",
		   hit5->plusp,hit3->plusp));
    *pair_relationship = 0;
    return 0;
  }

  if (hit5->chrnum != 0 && hit3->chrnum != 0) {
    for (q = hit3->substrings_1toN; q != NULL; q = List_next(q)) {
      substring3 = (Substring_T) List_head(q);
      for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	substring5 = (Substring_T) List_head(p);
	if (Substring_overlap_p(substring5,substring3)) {
	  debug10(printf("Calling Substring_insert_length on %d..%d and %d..%d\n",
			 Substring_querystart(substring5),Substring_queryend(substring5),
			 Substring_querystart(substring3),Substring_queryend(substring3)));
	  return Substring_insert_length(&(*pair_relationship),substring5,substring3);
	}
      }
    }
  }

  /* No overlap found between any combination of substrings */
  if (hit5->plusp == true) {
    if (hit5->genomicend > hit3->genomicstart + hit5->querylength + hit3->querylength) {
      debug10(printf("pair_insert_length: no overlap found, and %u - %u + %d + %d < 0, so returning 0\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicend - hit5->chroffset,
		     hit5->querylength,hit3->querylength));
      *pair_relationship = 0;
      return 0;
    } else {
      debug10(printf("pair_insert_length: no overlap found, so returning %u - %u + %d + %d\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicend - hit5->chroffset,
		     hit5->querylength,hit3->querylength));
    }
    *pair_relationship = +1;
    return hit3->genomicstart - hit5->genomicend + hit5->querylength + hit3->querylength;

  } else {
    if (hit3->genomicstart > hit5->genomicend + hit5->querylength + hit3->querylength) {
      debug10(printf("pair_insert_length: no overlap found, and %u - %u + %d + %d < 0, so returning 0\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicstart - hit3->chroffset,
		     hit5->querylength,hit3->querylength));
      *pair_relationship = 0;
      return 0;
    } else {
      debug10(printf("pair_insert_length: no overlap found, so returning %u - %u + %d + %d\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicstart - hit3->chroffset,
		     hit5->querylength,hit3->querylength));
      *pair_relationship = -1;
      return hit5->genomicend - hit3->genomicstart + hit5->querylength + hit3->querylength;
    }
  }
}



/* For unpaired ends */
static Chrpos_T
pair_insert_length_unpaired (Stage3end_T hit5, Stage3end_T hit3) {

  if (hit5->effective_chrnum != hit3->effective_chrnum) {
    debug10(printf("pair_insert_length: hit5->plusp %d != hit3->plusp %d, so returning 0\n",
		   hit5->plusp,hit3->plusp));
    return 0;
  } else if (hit5->distant_splice_p == true) {
    return 0;
  } else if (hit3->distant_splice_p == true) {
    return 0;
  } else if (hit5->high < hit3->low) {
    return hit3->low - hit5->high + hit5->querylength + hit3->querylength;
  } else if (hit3->high < hit5->low) {
    return hit5->low - hit3->high + hit5->querylength + hit3->querylength;
  } else {
    return hit5->querylength + hit3->querylength;
  }
}


Stage3end_T *
Stage3end_eval_and_sort_guided (int npaths, int *first_absmq, int *second_absmq, Stage3end_T guide,
				Stage3end_T *stage3array, char *queryuc_ptr,
				char *quality_string, bool displayp) {
  float maxlik, loglik;
  float total, q;		/* For Bayesian mapq calculation */
  int compute_npaths;

  int randomi, i;
  Stage3end_T temp, hit;

  if (npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (npaths == 1) {
    hit = stage3array[0];
    hit->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    hit->mapq_score = MAPQ_max_quality_score(quality_string,hit->querylength);
    hit->absmq_score = MAPQ_MAXIMUM_SCORE;

    if (displayp == true) {
      Stage3end_display_prep(hit,queryuc_ptr,/*alts_resolve*/-1,/*first_read_p*/true);
    }
    *first_absmq = hit->absmq_score;
    *second_absmq = 0;

  } else {
    /* Compute mapq_loglik */
    for (i = 0; i < npaths; i++) {
      Stage3end_compute_mapq(stage3array[i],quality_string);
    }

    /* Compute insert_length relative to guide.  This is the only change from the unguided procedure. */
    for (i = 0; i < npaths; i++) {
      stage3array[i]->guided_insertlength = pair_insert_length_unpaired(stage3array[i],guide);
    }

    /* Sort by nmatches, then mapq */
    qsort(stage3array,npaths,sizeof(Stage3end_T),Stage3end_output_cmp);

    if (want_random_p) {
      /* Randomize among best alignments */
      i = 1;
      while (i < npaths && Stage3end_output_cmp(&(stage3array[i]),&(stage3array[0])) == 0) {
	i++;
      }
      if (i > 1) {		/* i is number of ties */
	/* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
	randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) i));
	/* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
	temp = stage3array[0];
	stage3array[0] = stage3array[randomi];
	stage3array[randomi] = temp;
      }
    }

    /* Enforce monotonicity */
    for (i = npaths - 1; i > 0; i--) {
      if (stage3array[i-1]->mapq_loglik < stage3array[i]->mapq_loglik) {
	stage3array[i-1]->mapq_loglik = stage3array[i]->mapq_loglik;
      }
    }
    maxlik = stage3array[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < npaths; i++) {
      stage3array[i]->mapq_loglik -= maxlik;
    }

#if 0
    /* Save on computation if possible */
    /* Not possible, since we are going to select randomly from among all paths */
    if (npaths < maxpaths) {
      compute_npaths = npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }
#else
    compute_npaths = npaths;
#endif

    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3array[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3array[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = stage3array[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < npaths; i++) {
      total += (stage3array[i]->mapq_loglik = fasterexp(stage3array[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3array[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3array[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3array[i]->mapq_score = 96;
      } else {
	stage3array[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

    if (displayp == true) {
      /* Prepare for display */
      for (i = 0; i < compute_npaths; i++) {
	Stage3end_display_prep(stage3array[i],queryuc_ptr,/*alts_resolve*/-1,/*first_read_p*/true);
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3array[0]->mapq_score >= mapq_unique_score &&
	stage3array[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3end_free(&(stage3array[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3array;
}


/* Note: single-end terminals can be present with non-terminals when
   paired-end reads are searched for concordance, which can accumulate
   terminal alignments */

/* Pre-final: max (max-terminal, min-other)
   Final: max (min-terminal, max-GMAP, min-other) */


static List_T
Stage3end_optimal_score_prefinal (bool *eliminatedp, List_T hitlist, 
				  Hitlistpool_T hitlistpool, int querylength) {
  List_T optimal = NULL, p, q;
  T hit;
  Substring_T substring;
  Junction_T junction;
  int n;
  int cutoff_level;
  int minscore = querylength;
  int trim_querystart = 0, trim_queryend = 0, trim_querystart_0, trim_queryend_0;


#ifdef DISTANT_SPLICE_SPECIAL
  bool shortdistance_p = false;
#endif


  *eliminatedp = false;
  n = List_length(hitlist);
  debug4(printf("\nEntered Stage3end_optimal_score with %d hits\n",n));

  if (n <= 1) {
    return hitlist;
  }

  /* Use eventrim for comparing alignments.  Previously picked
     smallest trims, but now picking largest ones */
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    debug4(printf("hit %u..%u method %s, nsegments %d, nindels %d, trim_querystart: %d%s, trim_queryend %d%s, start_ambig %d, end_ambig %d.  sensedir %d\n",
		  hit->genomicstart - hit->chroffset,hit->genomicend - hit->chroffset,Method_string(hit->method),
		  hit->nsegments,hit->nindels,hit->trim_querystart,hit->trim_querystart_splicep ? " (splice)" : "",
		  hit->trim_queryend,hit->trim_queryend_splicep ? " (splice)" : "",
		  start_amb_length(hit),end_amb_length(hit),hit->sensedir));

    if (hit->trim_querystart_splicep == true) {
      /* Skip */
    } else if (hit->trim_querystart > trim_querystart) {
      trim_querystart = hit->trim_querystart;
    }
    if (hit->trim_queryend_splicep == true) {
      /* Skip */
    } else if (hit->trim_queryend > trim_queryend) {
      trim_queryend = hit->trim_queryend;
    }
  }

  if (trim_querystart == querylength) {
    trim_querystart = 0;
  }
  if (trim_queryend == querylength) {
    trim_queryend = 0;
  }
  debug4(printf("trim_querystart: %d, trim_queryend %d\n",trim_querystart,trim_queryend));

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

#ifdef CONSIDER_ENDS_IN_EVAL
    hit->score_eventrim = hit->trim_querystart / 8 + hit->trim_queryend / 8;
#else
    hit->score_eventrim = 0;
#endif

    debug4(printf("score OTHER:"));

    if (trim_querystart + trim_queryend >= querylength) {
      for (q = hit->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	hit->score_eventrim += Substring_nmismatches_bothdiff(substring);
      }
	
    } else {
      for (q = hit->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	trim_querystart_0 = trim_querystart;
	trim_queryend_0 = trim_queryend;
	if (Substring_mandatory_trim_querystart(substring) > trim_querystart_0) {
	  trim_querystart_0 = Substring_mandatory_trim_querystart(substring);
	}
	if (Substring_mandatory_trim_queryend(substring) > trim_queryend_0) {
	  trim_queryend_0 = Substring_mandatory_trim_queryend(substring);
	}
	hit->score_eventrim += Substring_count_mismatches_region(substring,trim_querystart_0,trim_queryend_0);
	debug4(printf("  substring (%d..%d) %d.",trim_querystart,trim_queryend,
		      Substring_count_mismatches_region(substring,trim_querystart_0,trim_queryend_0)));
      }
    }

    for (q = hit->junctions_1toN; q != NULL; q = List_next(q)) {
      junction = (Junction_T) List_head(q);
      if (Junction_nindels(junction) > 0) {
	hit->score_eventrim += indel_penalty_middle;
	debug4(printf(" => add %d.",indel_penalty_middle));
      }
    }


#if 0
    /* Accept a single indel */
#ifdef SCORE_INDELS_EVENTRIM
    if (hit->hittype == INSERTION || hit->hittype == DELETION) {
      debugee(printf("  indel at %d",hit->indel_pos));
      if (hit->indel_pos > trim_querystart && hit->indel_pos < querylength - trim_queryend) {
	hit->score_eventrim += indel_penalty_middle;
	debug4(printf(" => add %d.",indel_penalty_middle));
      }
    }
#endif
#endif
    debug4(printf("  RESULT: %d\n",hit->score_eventrim));

    if (hit->score_eventrim < minscore) {
      minscore = hit->score_eventrim;
    }
  }
  debug4(printf("MINSCORE: %d\n",minscore));


  /* Prefinal: Use score_eventrim */
  debug4(printf("Stage3end_optimal_score over %d hits: minscore = %d + subopt:%d\n",
		n,minscore,subopt_levels));
  minscore += subopt_levels;
  cutoff_level = minscore;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    if (hit->score_eventrim > cutoff_level + SCORE_EVENTRIM_SLOP) {
      debug4(printf("Prefinal: Eliminating hit %p at %u..%u with score_eventrim %d > cutoff_level %d\n",
		    hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
		    hit->score_eventrim,cutoff_level));
      Stage3end_free(&hit);
      *eliminatedp = true;

    } else {
      debug4(printf("Prefinal: Keeping hit %p at %u..%u with score_eventrim %d <= cutoff_level %d\n",
		    hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
		    hit->score_eventrim,cutoff_level));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hit);
    }
  }
  Hitlist_free(&hitlist);


#if 0
  /* Filter on nsegments */
  if (finalp == true && optimal != NULL) {
    hitlist = optimal;
    optimal = (List_T) NULL;

    hit = (T) hitlist->first;
    best_nsegments = hit->nsegments;

    for (p = hitlist; p != NULL; p = p->rest) {
      hit = (T) p->first;
      if (hit->nsegments < best_nsegments) {
	best_nsegments = hit->nsegments;
      }
    }

    for (p = hitlist; p != NULL; p = p->rest) {
      hit = (T) p->first;
      if (hit->nsegments > best_nsegments + 2) {
	debug4(printf("Eliminating a hit with nsegments %d\n",hit->nsegments));
	Stage3end_free(&hit);
	*eliminatedp = true;
      } else {
	debug4(printf("Keeping a hit with nsegments %d, nindels %d\n",hit->nsegments,hit->nindels));
	optimal = Hitlist_push(optimal,hitlitpool,(void *) hit);
      }
    }

    Hitlist_free(&hitlist);
  }
#endif

  debug4(printf("hitlist now has %d entries\n",List_length(optimal)));
  return optimal;
}


static int
hit_position_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  
  if (x->plusp < y->plusp) {
    return -1;
  } else if (y->plusp < x->plusp) {
    return -1;
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return +1;
  } else if (y->high < x->high) {
    return -1;
  } else {
    return 0;
  }
}

static bool
hit_overlap_p (T x, T y) {
  if (x->plusp != y->plusp) {
    return false;		/* Different strands */
  } else if (x->high < y->low) {
    return false;
  } else if (x->low > y->high) {
    return false;
  } else {
    return true;
  }
}

static List_T
Stage3end_optimal_score_final (bool *eliminatedp, List_T hitlist, Hitlistpool_T hitlistpool,
			       int querylength) {
  List_T optimal = NULL, p;
  T *hits, hit;
  int n, i, j, k;
  int best_nsegments;
  int best_nmatches_to_trims;
  double max_splice_score;
  int max_nmatches = 0, cutoff_level;
  /* int trim_querystart, trim_queryend, min_trim; */
  bool *eliminate, keptp;

  /* Relies on Path_solve_from_diagonals to maximize nsegments at each locus */

  *eliminatedp = false;
  n = List_length(hitlist);
  debug4(printf("\nEntered Stage3end_optimal_score with %d hits\n",n));

  if (n <= 1) {
    return hitlist;
  }

#ifdef DEBUG4
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (Stage3end_T) p->first;
    printf("%p %u..%u method %s, score_eventrim %d, nmatches %d (%d to_trims)\n",
		  hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
	   Method_string(hit->method),hit->score_eventrim,hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims);
  }
#endif

  /* Prune based on nmatches (to get the splice ends) */
  max_nmatches = 0;
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (Stage3end_T) p->first;
    if (hit->nmatches_plus_spliced_trims > max_nmatches) {
      max_nmatches = hit->nmatches_plus_spliced_trims;
      assert(max_nmatches <= querylength);
    }
  }

  cutoff_level = max_nmatches - subopt_levels;
  debug4(printf("cutoff level %d = max_nmatches %d\n",cutoff_level,max_nmatches));

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) p->first;

    if (hit->nmatches_plus_spliced_trims < cutoff_level - NMATCHES_SLOP) {
      debug4(printf("Final (nmatches %d < %d): Eliminating hit %p at %u..%u with nmatches %d (%d to_trims) < cutoff_level %d\n",
		    hit->nmatches_plus_spliced_trims,cutoff_level - NMATCHES_SLOP,
		    hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
		    hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,cutoff_level));
      Stage3end_free(&hit);
      *eliminatedp = true;

    } else {
      debug4(printf("Final (nmatches %d >= %d): Keeping hit %p at %u..%u with nmatches %d (%d to_trims) >= cutoff_level %d\n",
		    hit->nmatches_plus_spliced_trims,cutoff_level - NMATCHES_SLOP,
		    hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
		    hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,cutoff_level));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hit);
    }
  }
  Hitlist_free(&hitlist);


  /* Prune based on nmatches_to_trims */
  hitlist = optimal;
  optimal = (List_T) NULL;

  best_nmatches_to_trims = 0;
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (Stage3end_T) p->first;
    if (hit->nmatches_to_trims > best_nmatches_to_trims) {
      best_nmatches_to_trims = hit->nmatches_to_trims;
      assert(best_nmatches_to_trims <= querylength);
    }
  }

  cutoff_level = best_nmatches_to_trims - subopt_levels;
  debug4(printf("cutoff level %d = best_nmatches_to_trims %d\n",cutoff_level,best_nmatches_to_trims));

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) p->first;

    /* Do not allow slop for final */
    if (hit->nmatches_to_trims < cutoff_level /*- NMATCHES_TO_TRIMS_SLOP*/) {
      debug4(printf("Final (nmatches_to_trims %d < %d): Eliminating hit %p at %u..%u with nmatches_to_trims %d (%d to_trims) < cutoff_level %d\n",
		    hit->nmatches_to_trims,cutoff_level - NMATCHES_TO_TRIMS_SLOP,
		    hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
		    hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,cutoff_level));
      Stage3end_free(&hit);
      *eliminatedp = true;

    } else {
      debug4(printf("Final (nmatches_to_trims %d >= %d): Keeping hit %p at %u..%u with nmatches_to_trims %d (%d to_trims) >= cutoff_level %d\n",
		    hit->nmatches_to_trims,cutoff_level /*- NMATCHES_TO_TRIMS_SLOP*/,
		    hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
		    hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,cutoff_level));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hit);
    }
  }
  Hitlist_free(&hitlist);



  /* Eliminate within loci (1) */
  hitlist = optimal;
  optimal = (List_T) NULL;

  keptp = false;
  hits = (T *) List_to_array_n(&n,hitlist);
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  qsort(hits,n,sizeof(T),hit_position_cmp);
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hit_overlap_p(hits[j],hits[i]) == true) {
      j++;
    }
    if (j - i > 1) {
      best_nmatches_to_trims = 0;
      for (k = i; k < j; k++) {
	hit = hits[k];
	if (hit->nmatches_to_trims > best_nmatches_to_trims) {
	  best_nmatches_to_trims = hit->nmatches_to_trims;
	}
      }
      debug4(printf("best_nmatches_to_trims %d\n",best_nmatches_to_trims));

      for (k = i; k < j; k++) {
	hit = hits[k];
	/* Do not allow slop for final */
	if (hit->nmatches_to_trims < best_nmatches_to_trims /*- NMATCHES_TO_TRIMS_SLOP*/) {
	  debug4(printf("Within loci (nmatches_to_trims): Marking hit %p for elimination at %u..%u with nsegments %d, nmatches %d (%d to_trims), splice_score %f\n",
			hit,hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nsegments,
			hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->splice_score));
	  eliminate[k] = true;
	} else {
	  keptp = true;
	}
      }
    }

    i = j;
  }
  
  if (keptp == false) {
    optimal = hitlist;
  } else {
    for (k = 0; k < n; k++) {
      hit = hits[k];
      if (eliminate[k] == true) {
	debug4(printf("Within loci: Eliminating hit %p at %u..%u with nsegments %d, nmatches %d (%d to_trims), splice_score %f\n",
		      hit,hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nsegments,
		      hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->splice_score));
	Stage3end_free(&hit);
	*eliminatedp = true;
      } else {
	optimal = Hitlist_push(optimal,hitlistpool,(void *) hit);
      }
    }
    Hitlist_free(&hitlist);
  }
  FREE(hits);
  FREE(eliminate);


  /* Eliminate within loci (2) */
  hitlist = optimal;
  optimal = (List_T) NULL;

  keptp = false;
  hits = (T *) List_to_array_n(&n,hitlist);
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  qsort(hits,n,sizeof(T),hit_position_cmp);
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hit_overlap_p(hits[j],hits[i]) == true) {
      j++;
    }
    if (j - i > 1) {
      best_nsegments = 0;
      max_splice_score = 0.0;
      for (k = i; k < j; k++) {
	hit = hits[k];
	if (hit->nsegments > best_nsegments) {
	  best_nsegments = hit->nsegments;
	  max_splice_score = hit->splice_score;

	} else if (hit->nsegments == best_nsegments) {
	  if (hit->splice_score > max_splice_score) {
	    max_splice_score = hit->splice_score;
	  }
	}
      }
      debug8(printf("best_nsegments %d, max_splice_score %f\n",
		    best_nsegments,max_splice_score));

      for (k = i; k < j; k++) {
	hit = hits[k];
	if (hit->nsegments < best_nsegments) {
	debug4(printf("Within loci (nsegments %d < %d): Marking hit %p for elimination at %u..%u with nsegments %d, nmatches %d (%d to_trims), splice_score %f\n",
		      hit->nsegments,best_nsegments,
		      hit,hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nsegments,
		      hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->splice_score));
	  eliminate[k] = true;

	} else if (hit->splice_score < max_splice_score) {
	debug4(printf("Within loci (splice score %f < %f): Marking hit %p for elimination at %u..%u with nsegments %d, nmatches %d (%d to_trims), splice_score %f\n",
		      hit->splice_score,max_splice_score,
		      hit,hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nsegments,
		      hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->splice_score));
	  eliminate[k] = true;

	} else {
	  keptp = true;
	}
      }
    }

    i = j;
  }
  
  if (keptp == false) {
    optimal = hitlist;
  } else {
    for (k = 0; k < n; k++) {
      hit = hits[k];
      if (eliminate[k] == true) {
	debug4(printf("Within loci: Eliminating hit %p at %u..%u with nsegments %d, nmatches %d (%d to_trims), splice_score %f\n",
		      hit,hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nsegments,
		      hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->splice_score));
	Stage3end_free(&hit);
	*eliminatedp = true;
      } else {
	optimal = Hitlist_push(optimal,hitlistpool,(void *) hit);
      }
    }
    Hitlist_free(&hitlist);
  }
  FREE(hits);
  FREE(eliminate);

#if 0
  /* Filter on trim amount */
  hitlist = optimal;
  optimal = (List_T) NULL;
  min_trim = querylength;
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->trim_querystart_splicep == true) {
      /* Skip */
      trim_querystart = 0;
    } else {
      trim_querystart = hit->trim_querystart;
    }
    if (hit->trim_queryend_splicep == true) {
      /* Skip */
      trim_queryend = 0;
    } else {
      trim_queryend = hit->trim_queryend;
    }
      
    if (trim_querystart + trim_queryend < min_trim) {
      min_trim = trim_querystart + trim_queryend;
    }
  }
    
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->trim_querystart_splicep == true) {
      /* Skip */
      trim_querystart = 0;
    } else {
      trim_querystart = hit->trim_querystart;
    }
    if (hit->trim_queryend_splicep == true) {
      /* Skip */
      trim_queryend = 0;
    } else {
      trim_queryend = hit->trim_queryend;
    }
      
    if (trim_querystart + trim_queryend > min_trim) {
      debug4(printf("Final: Eliminating hit %p at %u..%u with trim %d + %d > min_trim %d\n",
		    hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
		    trim_querystart,trim_queryend,min_trim));
      Stage3end_free(&hit);
      *eliminatedp = true;
	
    } else {
      debug4(printf("Final: Keeping hit %p at %u..%u with trim %d + %d == min_trim %d\n",
		    hit,hit->low - hit->chroffset,hit->high - hit->chroffset,
		    trim_querystart,trim_queryend,min_trim));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hit);
    }
  }
  Hitlist_free(&hitlist);
#endif


  debug4(printf("Exiting Stage3end_optimal_score_final with %d hits\n",List_length(optimal)));
  return optimal;
}



List_T
Stage3end_optimal_score (List_T hitlist, Hitlistpool_T hitlistpool, int querylength, bool finalp) {
  List_T optimal;
  bool eliminatedp;

  if (finalp == false) {
    optimal = Stage3end_optimal_score_prefinal(&eliminatedp,hitlist,hitlistpool,querylength);
    while (eliminatedp == true) {
      optimal = Stage3end_optimal_score_prefinal(&eliminatedp,optimal,hitlistpool,querylength);
    }

  } else {
    optimal = Stage3end_optimal_score_final(&eliminatedp,hitlist,hitlistpool,querylength);
    while (eliminatedp == true) {
      optimal = Stage3end_optimal_score_final(&eliminatedp,optimal,hitlistpool,querylength);
    }
  }

  return optimal;
}


static void
unalias_circular (T hit) {
  Chrpos_T chrlength = hit->chrlength;
  List_T p;
  Substring_T substring;

  assert(hit->circularalias == +1);
  debug12(printf("Calling unalias_circular on substrings\n"));
  for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    Substring_unalias_circular(substring);
  }

  /* Doesn't fix hitpair->low and hitpair->high */
  hit->genomicstart -= chrlength;
  hit->genomicend -= chrlength;
  hit->low -= chrlength;
  hit->high -= chrlength;

  hit->circularalias = -1;

  return;
}


#if 0
List_T
Stage3end_unalias_circular (List_T hitlist) {
  List_T p;
  T hit;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->circularalias == +1) {
      unalias_circular(hit);
    }
  }

  return hitlist;
}
#endif

List_T
Stage3end_remove_circular_alias (List_T hitlist, Hitlistpool_T hitlistpool) {
  List_T newlist = NULL, p;
  T hit;

  debug12(printf("Calling Stage3end_remove_circular_alias on %d hits\n",List_length(hitlist)));
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    if (hit->circularalias == +1) {
      /* First, try to salvage alias +1 */
      unalias_circular(hit);
    }

    if (hit->chrnum == 0) {
      /* Translocation */
      newlist = Hitlist_push(newlist,hitlistpool,(void *) hit);

    } else if (Stage3end_chrpos_low_trim(hit) >= hit->chrlength) {
      /* All in circular alias */
      debug12(printf("Freeing hit because all is in circular alias\n"));
      Stage3end_free(&hit);

    } else {
      newlist = Hitlist_push(newlist,hitlistpool,(void *) hit);
    }
  }

  Hitlist_free(&hitlist);
  return newlist;
}


#if 0
int
Stage3end_noptimal (List_T hitlist, int querylength) {
  int noptimal;
  List_T p;
  T hit;
  int minscore = querylength;

  noptimal = 0;
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->score < minscore) {
      minscore = hit->score;
      noptimal = 0;
    }
    if (hit->score == minscore) {
      noptimal++;
    }
  }

  return noptimal;
}
#endif


static Univcoord_T
normalize_coord (Univcoord_T orig, int circularalias, Chrpos_T chrlength) {
  if (circularalias == +1) {
    return orig - chrlength;
  } else {
    return orig;
  }
}



static int
duplicate_sort_cmp (const void *a, const void *b) {
  int cmp;
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoord_T x_genomicstart, y_genomicstart;
  Univcoord_T x_genomicend, y_genomicend;
  List_T p, q;
  Substring_T x_substring, y_substring;

  if (altlocp[x->chrnum] == true && altlocp[y->chrnum] == true) {
    if (alias_ends[y->chrnum] >= alias_starts[x->chrnum] &&
	alias_starts[y->chrnum] <= alias_ends[x->chrnum]) {
      /* The primary regions overlap */
      return 0;
    } else if (alias_starts[x->chrnum] < alias_starts[y->chrnum]) {
      return -1;
    } else if (alias_starts[y->chrnum] < alias_starts[x->chrnum]) {
      return +1;
    } else if (alias_ends[x->chrnum] < alias_ends[y->chrnum]) {
      return -1;
    } else if (alias_ends[y->chrnum] < alias_ends[x->chrnum]) {
      return +1;
    } else {
      return 0;
    }

  } else if (altlocp[x->chrnum] == true) {
    if (y->genomicend >= alias_starts[x->chrnum] &&
	y->genomicstart <= alias_ends[x->chrnum]) {
      /* y overlaps with the primary region for x */
      return +1;		/* Put primary region first */
    }
    /* Don't overlap, so fall through to rest of procedure */

  } else if (altlocp[y->chrnum] == true) {
    if (alias_ends[y->chrnum] >= x->genomicstart &&
	alias_starts[y->chrnum] <= x->genomicend) {
      /* x overlaps with the primary region for y */
      return -1;		/* Put primary region first */
    }
    /* Don't overlap, so fall through to rest of procedure */
  }


  x_genomicstart = normalize_coord(x->genomicstart,x->circularalias,x->chrlength);
  x_genomicend = normalize_coord(x->genomicend,x->circularalias,x->chrlength);

  y_genomicstart = normalize_coord(y->genomicstart,y->circularalias,y->chrlength);
  y_genomicend = normalize_coord(y->genomicend,y->circularalias,y->chrlength);

  
  if (x_genomicstart < y_genomicstart) {
    return -1;
  } else if (x_genomicstart > y_genomicstart) {
    return +1;
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (x->hittype > y->hittype) {
    return +1;
  } else if (x_genomicend < y_genomicend) {
    return -1;
  } else if (x_genomicend > y_genomicend) {
    return +1;

    /* sensedir is relevant for transcriptome-guided alignment, with overlapping genes */
  } else if (x->sensedir > y->sensedir) {
    return -1;
  } else if (y->sensedir > x->sensedir) {
    return +1;

  } else {
    for (p = x->substrings_1toN, q = y->substrings_1toN; p != NULL && q != NULL; p = List_next(p), q = List_next(q)) {
      x_substring = (Substring_T) List_head(p);
      y_substring = (Substring_T) List_head(q);
      if ((cmp = Substring_compare(x_substring,y_substring,x->circularalias,y->circularalias,x->chrlength,y->chrlength)) != 0) {
	return cmp;
      }
    }
    if (p == NULL && q != NULL) {
      return -1;
    } else if (p != NULL && q == NULL) {
      return +1;
    }

#if 0
    /* Need to change to search on junctions */
    if (x->indel_low < y->indel_low) {
      return -1;
    } else if (y->indel_low < x->indel_low) {
      return +1;
    }
#endif

    return 0;
  }
}

/* Same as duplicate_sort_cmp, except for indel_low */
static int
duplicate_equiv_cmp (const void *a, const void *b) {
  int cmp;
  T x = * (T *) a;
  T y = * (T *) b;
  List_T p, q;
  Substring_T x_substring, y_substring;

  Univcoord_T x_genomicstart, x_genomicend, y_genomicstart, y_genomicend;

  if (altlocp[x->chrnum] == true && altlocp[y->chrnum] == true) {
    if (alias_ends[y->chrnum] >= alias_starts[x->chrnum] &&
	alias_starts[y->chrnum] <= alias_ends[x->chrnum]) {
      /* The primary regions overlap */
      return 0;
    }

  } else if (altlocp[x->chrnum] == true) {
    if (y->genomicend >= alias_starts[x->chrnum] &&
	y->genomicstart <= alias_ends[x->chrnum]) {
      /* y overlaps with the primary region for x */
      return 0;
    }

  } else if (altlocp[y->chrnum] == true) {
    if (alias_ends[y->chrnum] >= x->genomicstart &&
	alias_starts[y->chrnum] <= x->genomicend) {
      /* x overlaps with the primary region for y */
      return 0;
    }
  }

  x_genomicstart = normalize_coord(x->genomicstart,x->circularalias,x->chrlength);
  x_genomicend = normalize_coord(x->genomicend,x->circularalias,x->chrlength);

  y_genomicstart = normalize_coord(y->genomicstart,y->circularalias,y->chrlength);
  y_genomicend = normalize_coord(y->genomicend,y->circularalias,y->chrlength);

  if (x_genomicstart < y_genomicstart) {
    return -1;
  } else if (x_genomicstart > y_genomicstart) {
    return +1;
#if 0
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (x->hittype > y->hittype) {
    return +1;
#endif
  } else if (x_genomicend < y_genomicend) {
    return -1;
  } else if (x_genomicend > y_genomicend) {
    return +1;

    /* sensedir is relevant for transcriptome-guided alignment, with overlapping genes */
  } else if (x->sensedir > y->sensedir) {
    return -1;
  } else if (y->sensedir > x->sensedir) {
    return +1;

  } else {
    for (p = x->substrings_1toN, q = y->substrings_1toN; p != NULL && q != NULL; p = List_next(p), q = List_next(q)) {
      x_substring = (Substring_T) List_head(p);
      y_substring = (Substring_T) List_head(q);
      if ((cmp = Substring_compare(x_substring,y_substring,x->circularalias,y->circularalias,x->chrlength,y->chrlength)) != 0) {
	return cmp;
      }
    }
    if (p == NULL && q != NULL) {
      return -1;
    } else if (p != NULL && q == NULL) {
      return +1;
    } else {
      return 0;
    }
  }
}


#if defined(DEBUG0) || defined(DEBUG4)
static void
Stage3end_print_substrings (Stage3end_T hit) {
  List_T p;
  Substring_T substring;

  for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
    if ((substring = (Substring_T) List_head(p)) == NULL) {
      printf("NA ");
    } else {
      printf("#%d:%llu..%llu ",
	     Substring_chrnum(substring),
	     (unsigned long long) Substring_alignstart_trim(substring),
	     (unsigned long long) Substring_alignend_trim(substring));
    }
  }
  return;
}
#endif


const Except_T Duplicate_Pairing = { "Duplicates both seen in pairing" };

List_T
Stage3end_remove_duplicates (List_T hitlist, Hitlistpool_T hitlistpool) {
#ifdef DEBUG4
  List_T p;
#endif
  T x, y, *hits;
  int n, usedi, i, j, k;
  bool *eliminate, eliminatep;

  debug4(printf("Entered Stage3end_remove_duplicates with %d hits\n",List_length(hitlist)));
  if ((n = List_length(hitlist)) == 0) {
    return (List_T) NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    hits = (T *) MALLOCA(n * sizeof(T));    
    List_fill_array((void **) hits,hitlist); /* hitlist is a return value */
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hits = (T *) List_to_array(hitlist,NULL);
#endif
  }


  /* By equivalence */
  debug4(printf("Stage3end_remove_duplicates: checking %d hits by equivalence class\n",n));
  qsort(hits,n,sizeof(T),duplicate_sort_cmp);

  debug4(
	 for (i = 0; i < n; i++) {
	   x = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, circularalias %d, nmatches %d (%d to_trims), score %d, sense %d ",
		  i,Method_string(x->method),x,x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
		  x->circularalias,x->nmatches_plus_spliced_trims,x->nmatches_to_trims,x->score_within_trims,x->sensedir);
	   Stage3end_print_substrings(x);
	   if (x->transcripts != NULL) {
	     Transcript_print_list(x->transcripts);
	   }
	   printf("\n");
	 }
	 );

  eliminatep = false;
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && duplicate_equiv_cmp(&(hits[j]),&(hits[i])) == 0) {
      j++;
    }

    if (j > i+1) {
      debug4(printf("Equivalence class #%d through #%d.  ",i,j-1));

      x = hits[i];
      if (x->paired_usedp == true) {
	usedi = i;
      } else {
	usedi = -1;
      }
      
      for (k = i+1; k < j; k++) {
	y = hits[k];
	if (y->paired_usedp == true) {
	  if (usedi >= 0) {
	    debug4(printf("  #%d equivalent to #%d and both used (%p and %p)\n",k,usedi,hits[k],hits[usedi]));
#if 0
	    /* This doesn't matter anymore.  Example from NM_001033853:
	       TTGCCCTTGGTCACCCCGATGACGTCGATCATCTCATCCTGCCCAAACACTTGGTTCACAGGTACCTGCTGCTCA
	       AGTGATGAATCCAAGAGGCGTTTCTATAAGAATTGGCATAAATCTAAGAAGAAGGCCCACCTGATGGAGATCCAG */
	    fprintf(stderr,"Duplicates of Stage3end_T both seen\n");
#if 0
	    /* No longer providing queryseq1 and queryseq2 */
	    Shortread_print_query_pairedend_fasta(stderr,queryseq1,queryseq2,
						  /*invert_first_p*/false,/*invert_second_p*/true);
#endif
	    Except_raise(&Duplicate_Pairing, __FILE__, __LINE__);
#endif
	  } else {
	    usedi = k;
	  }
	}
      }

      if (usedi < 0) {
	debug4(printf("None used yet so eliminating #%d through #%d\n",i+1,j-1));
	for (k = i+1; k < j; k++) {
	  y = hits[k];
	  if (y->transcripts != NULL) {
	    x->transcripts = List_append(y->transcripts,x->transcripts);
	    y->transcripts = (List_T) NULL;
	  }
	  eliminate[k] = true;
	  eliminatep = true;
	}
      } else {
	debug4(printf("One used already so eliminating all but #%d\n",usedi));
	for (k = i; k < j; k++) {
	  if (k != usedi) {
	    y = hits[k];
	    if (y->transcripts != NULL) {
	      x->transcripts = List_append(y->transcripts,x->transcripts);
	      y->transcripts = (List_T) NULL;
	    }
	    eliminate[k] = true;
	    eliminatep = true;
	  }
	}
      }
    }

    i = j;
  }
    

#if 0
  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
  }
#endif

  if (eliminatep == false) {
    debug4(printf("No eliminations, so hitlist is unchanged\n"));
  } else {
    Hitlist_free(&hitlist);
    for (i = n-1; i >= 0; i--) {
      x = hits[i];
      if (eliminate[i] == false) {
#ifdef DEBUG4
	printf("  Keeping #%d:%u..%u, score %d, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d, sensedir = %d) ",
	       x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
	       x->score_within_trims,x->nmatches_plus_spliced_trims,x->nindels,x->distance,x->chrnum,x->plusp,x->sensedir);
	Stage3end_print_substrings(x);
	if (x->transcripts != NULL) {
	  Transcript_print_nums(x->transcripts);
	}
	printf("\n");
#endif
	hitlist = Hitlist_push(hitlist,hitlistpool,(void *) x);

      } else {
#ifdef DEBUG4
	printf("  Eliminating #%d:%u..%u, score %d, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d, sensedir = %d) ",
	       x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
	       x->score_within_trims,x->nmatches_plus_spliced_trims,x->nindels,x->distance,x->chrnum,x->plusp,x->sensedir);
	Stage3end_print_substrings(x);
	if (x->transcripts != NULL) {
	  Transcript_print_nums(x->transcripts);
	}
	printf("\n");
#endif
	Stage3end_free(&x);
      }
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hits);
  FREEA(eliminate);
#else
  FREE(hits);
  FREE(eliminate);
#endif

#ifdef DEBUG4
  for (p = hitlist, i = 0; p != NULL; p = p->rest, i++) {
    x = (T) p->first;
    printf("  Final %d: #%d:%u..%u (plusp = %d, sensedir = %d) ",
	   i,x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,x->plusp,x->sensedir);
    Stage3end_print_substrings(x);
    if (x->transcripts != NULL) {
      Transcript_print_nums(x->transcripts);
    }
    printf("\n");
  }
#endif

  debug4(printf("Exited Stage3end_remove_duplicates with %d hits\n",List_length(hitlist)));
  return hitlist;
}



T *
Stage3end_remove_duplicates_array (int *nunique, List_T *duplicates, T *hits, int nhits,
				   Hitlistpool_T hitlistpool) {
  T *unique, *out, x, y;
  int usedi, i, j, k;
  bool *eliminate, eliminatep;

  debug4(printf("Entered Stage3end_remove_duplicates_array with %d hits\n",nhits));
  if (nhits == 0) {
    *nunique = 0;
    return (T *) NULL;

  } else {
    eliminate = (bool *) CALLOC(nhits,sizeof(bool));
  }


  /* By equivalence */
  debug4(printf("Stage3end_remove_duplicates_array: checking %d hits by equivalence class\n",nhits));
  qsort(hits,nhits,sizeof(T),duplicate_sort_cmp);

  debug4(
	 for (i = 0; i < nhits; i++) {
	   x = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, circularalias %d, nmatches %d (%d to_trims), score %d, sense %d ",
		  i,Method_string(x->method),x,x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
		  x->circularalias,x->nmatches_plus_spliced_trims,x->nmatches_to_trims,x->score_within_trims,x->sensedir);
	   Stage3end_print_substrings(x);
	   if (x->transcripts != NULL) {
	     Transcript_print_list(x->transcripts);
	   }
	   printf("\n");
	 }
	 );

  eliminatep = false;
  i = 0;
  while (i < nhits) {
    j = i+1;
    while (j < nhits && duplicate_equiv_cmp(&(hits[j]),&(hits[i])) == 0) {
      j++;
    }

    if (j > i+1) {
      debug4(printf("Equivalence class #%d through #%d.  ",i,j-1));

      x = hits[i];
      if (x->paired_usedp == true) {
	usedi = i;
      } else {
	usedi = -1;
      }
      
      for (k = i+1; k < j; k++) {
	y = hits[k];
	if (y->paired_usedp == true) {
	  if (usedi >= 0) {
	    debug4(printf("  #%d equivalent to #%d and both used (%p and %p)\n",k,usedi,hits[k],hits[usedi]));
#if 0
	    /* This doesn't matter anymore.  Example from NM_001033853:
	       TTGCCCTTGGTCACCCCGATGACGTCGATCATCTCATCCTGCCCAAACACTTGGTTCACAGGTACCTGCTGCTCA
	       AGTGATGAATCCAAGAGGCGTTTCTATAAGAATTGGCATAAATCTAAGAAGAAGGCCCACCTGATGGAGATCCAG */
	    fprintf(stderr,"Duplicates of Stage3end_T both seen\n");
#if 0
	    /* No longer providing queryseq1 and queryseq2 */
	    Shortread_print_query_pairedend_fasta(stderr,queryseq1,queryseq2,
						  /*invert_first_p*/false,/*invert_second_p*/true);
#endif
	    Except_raise(&Duplicate_Pairing, __FILE__, __LINE__);
#endif
	  } else {
	    usedi = k;
	  }
	}
      }

      if (usedi < 0) {
	debug4(printf("None used yet so eliminating #%d through #%d\n",i+1,j-1));
	for (k = i+1; k < j; k++) {
	  y = hits[k];
	  if (y->transcripts != NULL) {
	    x->transcripts = List_append(y->transcripts,x->transcripts);
	    y->transcripts = (List_T) NULL;
	  }
	  eliminate[k] = true;
	  eliminatep = true;
	}
      } else {
	debug4(printf("One used already so eliminating all but #%d\n",usedi));
	for (k = i; k < j; k++) {
	  if (k != usedi) {
	    y = hits[k];
	    if (y->transcripts != NULL) {
	      x->transcripts = List_append(y->transcripts,x->transcripts);
	      y->transcripts = (List_T) NULL;
	    }
	    eliminate[k] = true;
	    eliminatep = true;
	  }
	}
      }
    }

    i = j;
  }
    

#if 0
  nkept = 0;
  for (i = 0; i < nhits; i++) {
    if (eliminate[i] == false) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
  }
#endif

  if (eliminatep == false) {
    debug4(printf("No eliminations, so hits are unchanged\n"));
    unique = hits;
    *nunique = nhits;

  } else {
    /* Caller needs (*nunique)+1, but since we are guaranteed to have one elimination, nhits will suffice */
    out = unique = (T *) MALLOC(nhits*sizeof(T));

    for (i = nhits-1; i >= 0; i--) {
      x = hits[i];
      if (eliminate[i] == false) {
#ifdef DEBUG4
	printf("  Keeping #%d:%u..%u, score %d, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d, sensedir = %d) ",
	       x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
	       x->score_within_trims,x->nmatches_plus_spliced_trims,x->nindels,x->distance,x->chrnum,x->plusp,x->sensedir);
	Stage3end_print_substrings(x);
	if (x->transcripts != NULL) {
	  Transcript_print_nums(x->transcripts);
	}
	printf("\n");
#endif
	*out++ = x;

      } else {
#ifdef DEBUG4
	printf("  Eliminating #%d:%u..%u, score %d, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d, sensedir = %d) ",
	       x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,
	       x->score_within_trims,x->nmatches_plus_spliced_trims,x->nindels,x->distance,x->chrnum,x->plusp,x->sensedir);
	Stage3end_print_substrings(x);
	if (x->transcripts != NULL) {
	  Transcript_print_nums(x->transcripts);
	}
	printf("\n");
#endif
	/* Stage3end_free(&x); -- Cannot free, because newladder and ladder might share this hit */
	*duplicates = Hitlist_push(*duplicates,hitlistpool,(void *) x);
      }
    }

    *nunique = out - unique;
    FREE(hits);
  }

  FREE(eliminate);

#ifdef DEBUG4
  for (i = 0; i < *nunique; i++) {
    x = unique[i];
    printf("  Final %d: #%d:%u..%u (plusp = %d, sensedir = %d) ",
	   i,x->chrnum,x->genomicstart - x->chroffset,x->genomicend - x->chroffset,x->plusp,x->sensedir);
    Stage3end_print_substrings(x);
    if (x->transcripts != NULL) {
      Transcript_print_nums(x->transcripts);
    }
    printf("\n");
  }
#endif

  debug4(printf("Exited Stage3end_remove_duplicates_array with %d hits\n",*nunique));
  return unique;
}



#if 0
static bool
extra_ambiguous_ends_p (List_T substrings) {
  int nambiguous;
  List_T p;

  p = substrings;
  nambiguous = 0;
  while (Substring_ambiguous_p((Substring_T) List_head(p)) == true) {
    p = List_next(p);
    nambiguous += 1;
  }
  if (nambiguous > 1) {
    return true;
  }

  substrings = List_reverse(substrings);

  p = substrings;
  nambiguous = 0;
  while (Substring_ambiguous_p((Substring_T) List_head(p)) == true) {
    p = List_next(p);
    nambiguous += 1;
  }

  substrings = List_reverse(substrings);

  if (nambiguous > 1) {
    return true;
  } else {
    return false;
  }
}
#endif




#if 0
static int
recompute_nmismatches (T hit) {
  int nmismatches = 0;
  List_T p;
  Substring_T substring;

  assert(hit->method == TR);

  for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    nmismatches += Substring_nmismatches_bothdiff(substring);
  }

  return nmismatches;
}
#endif



#if 0
static void
zap_before (List_T *substrings, List_T *junctions, Substring_T substring) {
  void *ignore;

#if 0
  List_T p;
  printf("Zap before looking for %p in",substring);
  for (p = *substrings; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");
#endif
  
  while ((Substring_T) List_head(*substrings) != substring) {
    *substrings = Listpool_pop(*substrings,&ignore);
    *junctions = Listpool_pop(*junctions,&ignore);
  }
  return;
}

static void
zap_before_1 (List_T *substrings, Substring_T substring) {
  void *ignore;

#if 0
  List_T p;
  printf("Zap_before_1 looking for %p in",substring);
  for (p = *substrings; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");
#endif

  while ((Substring_T) List_head(*substrings) != substring) {
    *substrings = Listpool_pop(*substrings,&ignore);
  }
  return;
}
#endif


#if 0
static void
zap_before_free (List_T *substrings, List_T *junctions, Substring_T substring) {
  Substring_T unwanted1;
  Junction_T unwanted2;

#if 0
  List_T p;
  printf("Zap_before_free looking for %p in",substring);
  for (p = *substrings; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");
#endif

  while ((Substring_T) List_head(*substrings) != substring) {
    *substrings = Listpool_pop(*substrings,(void **) &unwanted1);
    Substring_free(&unwanted1);
    *junctions = Listpool_pop(*junctions,(void **) &unwanted2);
    Junction_free(&unwanted2);
  }
  return;
}

static void
zap_after_1 (List_T substrings, Substring_T substring) {
  List_T save;

#if 0
  List_T p;
  printf("Zap_after_1 looking for %p in",substring);
  for (p = substrings; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");
#endif

  while ((Substring_T) List_head(substrings) != substring) {
    substrings = List_next(substrings);
  }

  save = substrings;
  substrings = List_next(substrings);
  List_tail_set(save,(List_T) NULL);
  /* List_free(&substrings); -- allocated by Listpool_push */

  return;
}
#endif


/* For zap_after functions, junctions needs to lag behind substrings
   by one link, and could potentially be changed to NULL */
#if 0
static void
zap_after (List_T substrings, List_T *junctions, Substring_T substring) {
  List_T ptr, save;

#if 0
  List_T p;
  printf("Zap_after looking for %p in",substring);
  for (p = substrings; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");
#endif

  if ((Substring_T) List_head(substrings) == substring) {
    save = substrings;
    substrings = List_next(substrings);
    List_tail_set(save,(List_T) NULL);
    /* List_free(&substrings); -- allocated by Listpool_push */
    
    /* List_free(&(*junctions)); -- allocated by Listpool_push */
    *junctions = NULL;

  } else {
    substrings = List_next(substrings);
    ptr = *junctions;
    while ((Substring_T) List_head(substrings) != substring) {
      substrings = List_next(substrings);
      ptr = List_next(ptr);
    }

    save = substrings;
    substrings = List_next(substrings);
    List_tail_set(save,(List_T) NULL);
    /* List_free(&substrings); -- allocated by Listpool_push */

    save = ptr;
    ptr = List_next(ptr);
    List_tail_set(save,(List_T) NULL);
    /* List_free(&ptr); -- allocated by Listpool_push */
    /* *junctions is the same */
  }

  return;
}
#endif


#if 0
static void
zap_after_free (List_T substrings, List_T *junctions, Substring_T substring) {
  List_T ptr, save;
  Substring_T unwanted1;
  Junction_T unwanted2;

#if 0
  List_T p;
  printf("Zap_after_free looking for %p in",substring);
  for (p = substrings; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");
#endif

  if ((Substring_T) List_head(substrings) == substring) {
    save = substrings;
    substrings = List_next(substrings);
    List_tail_set(save,(List_T) NULL);
    while (substrings != NULL) {
      substrings = Listpool_pop(substrings,(void **) &unwanted1);
      Substring_free(&unwanted1);
    }
    
    while (*junctions != NULL) {
      *junctions = List_pop(*junctions,(void **) &unwanted2);
      Junction_free(&unwanted2);
    }
    /* *junctions == NULL */

  } else {
    substrings = List_next(substrings);
    ptr = *junctions;
    while ((Substring_T) List_head(substrings) != substring) {
      substrings = List_next(substrings);
      ptr = List_next(ptr);
    }

    save = substrings;
    substrings = List_next(substrings);
    List_tail_set(save,(List_T) NULL);
    while (substrings != NULL) {
      substrings = Listpool_pop(substrings,(void **) &unwanted1);
      Substring_free(&unwanted1);
    }

    save = ptr;
    ptr = List_next(ptr);
    List_tail_set(save,(List_T) NULL);
    while (ptr != NULL) {
      ptr = List_pop(ptr,(void **) &unwanted2);
      Junction_free(&unwanted2);
    }
    /* *junctions is the same */
  }

  return;
}
#endif




#if 0
List_T
Stage3end_reject_trimlengths (List_T hits, Hitlistpool_T hitlistpool) {
  List_T filtered = NULL, p;
  T hit;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->trim_querystart + hit->trim_queryend >= reject_trimlength) {
      Stage3end_free(&hit);
    } else {
      filtered = Hitlist_push(filtered,hitlistpool,(void *) hit);
    }
  }

  Hitlist_free(&hits);
  return filtered;
}
#endif


/* Used for eliminating exact duplicates.  Also sorts secondarily by hittype. */
static int
hit_sort_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  debug4(printf("Comparing %s: #%d:%u..%u, circularalias %d, nmatches %d (%d to_trims), score %d with %s: #%d:%u..%u, circularalias %d, nmatches %d (%d to_trims), score %d\n",
		Method_string(x->method),x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,
		x->circularalias,x->nmatches_plus_spliced_trims,x->nmatches_to_trims,x->score_within_trims,
		Method_string(y->method),y->chrnum,y->genomicstart-y->chroffset,y->genomicend-y->chroffset,
		y->circularalias,y->nmatches_plus_spliced_trims,x->nmatches_to_trims,y->score_within_trims));

  if (altlocp[x->chrnum] == true && altlocp[y->chrnum] == true) {
    if (alias_ends[y->chrnum] >= alias_starts[x->chrnum] &&
	alias_starts[y->chrnum] <= alias_ends[x->chrnum]) {
      /* The primary regions overlap */
      return 0;
    } else if (alias_starts[x->chrnum] < alias_starts[y->chrnum]) {
      return -1;
    } else if (alias_starts[y->chrnum] < alias_starts[x->chrnum]) {
      return +1;
    } else if (alias_ends[x->chrnum] < alias_ends[y->chrnum]) {
      return -1;
    } else if (alias_ends[y->chrnum] < alias_ends[x->chrnum]) {
      return +1;
    } else {
      return 0;
    }

  } else if (altlocp[x->chrnum] == true) {
    if (y->genomicend >= alias_starts[x->chrnum] &&
	y->genomicstart <= alias_ends[x->chrnum]) {
      /* y overlaps with the primary region for x */
      return +1;		/* Put primary region first */
    }
    /* Don't overlap, so fall through to rest of procedure */

  } else if (altlocp[y->chrnum] == true) {
    if (alias_ends[y->chrnum] >= x->genomicstart &&
	alias_starts[y->chrnum] <= x->genomicend) {
      /* x overlaps with the primary region for y */
      return -1;		/* Put primary region first */
    }
    /* Don't overlap, so fall through to rest of procedure */
  }


  if (x->plusp > y->plusp) {
    return -1;
  } else if (y->plusp > x->plusp) {
    return +1;

#if 0
  } else if (x->high < y->low) {
    return -1;
  } else if (y->high < x->low) {
    return +1;
#else
  } else if (x->low < y->low) {
    debug4(printf("Returning -1 for low\n"));
    return -1;
  } else if (y->low < x->low) {
    debug4(printf("Returning +1 for low\n"));
    return +1;

  } else if (x->high < y->high) {
    debug4(printf("Returning -1 for high\n"));
    return -1;
  } else if (y->high < x->high) {
    debug4(printf("Returning +1 for high\n"));
    return +1;
#endif

  } else if (x->score_within_trims < y->score_within_trims) {
    return -1;
  } else if (y->score_within_trims < x->score_within_trims) {
    return +1;
  } else if (x->nmatches_plus_spliced_trims > y->nmatches_plus_spliced_trims) {
    return -1;
  } else if (y->nmatches_plus_spliced_trims > x->nmatches_plus_spliced_trims) {
    return +1;
#if 0
  } else if (x->nmatches_to_trims > y->nmatches_plus_to_trims) {
    return -1;
  } else if (y->nmatches_to_trims > x->nmatches_plus_to_trims) {
    return +1;
#endif

#if 0
  } else if (x->nsplices < y->nsplices) {
    return -1;
  } else if (y->nsplices < x->nsplices) {
    return +1;
#endif

#if 0
    /* Hittype is irrelevant */
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (y->hittype < x->hittype) {
    return +1;
#endif

    /* Prioritize last method used */
  } else if (x->method > y->method) {
    return -1;
  } else if (y->method > x->method) {
    return +1;

#if 0
  } else if (y->start_amb_length + y->end_amb_length == 0 &&
	     x->start_amb_length + x->end_amb_length > 0) {
    return -1;
  } else if (x->start_amb_length + x->end_amb_length == 0 &&
	     y->start_amb_length + y->end_amb_length > 0) {
    return +1;
#endif

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

  } else if (x->altlocp < y->altlocp) {
    return -1;
  } else if (y->altlocp < x->altlocp) {
    return +1;

  } else if (x->sensedir != 0 && y->sensedir == 0) {
    return -1;
  } else if (y->sensedir != 0 && x->sensedir == 0) {
    return +1;

  } else {
    debug4(printf("Returning 0 for equivalent\n"));
    return 0;
  }
}

/* Same as hit_sort_cmp, except for hittype, nmatches_to_trims, and indel_low */
static int
hit_equiv_cmp (Stage3end_T x, Stage3end_T y) {

  if (altlocp[x->chrnum] == true && altlocp[y->chrnum] == true) {
    if (alias_ends[y->chrnum] >= alias_starts[x->chrnum] &&
	alias_starts[y->chrnum] <= alias_ends[x->chrnum]) {
      /* The primary regions overlap */
      return 0;
    }

  } else if (altlocp[x->chrnum] == true) {
    if (y->genomicend >= alias_starts[x->chrnum] &&
	y->genomicstart <= alias_ends[x->chrnum]) {
      /* y overlaps with the primary region for x */
      return 0;
    }

  } else if (altlocp[y->chrnum] == true) {
    if (alias_ends[y->chrnum] >= x->genomicstart &&
	alias_starts[y->chrnum] <= x->genomicend) {
      /* x overlaps with the primary region for y */
      return 0;
    }
  }

  if (x->plusp > y->plusp) {
    return -1;
  } else if (y->plusp > x->plusp) {
    return +1;
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return +1;
  } else if (y->high < x->high) {
    return -1;

  } else if (x->score_within_trims < y->score_within_trims) {
    return -1;
  } else if (y->score_within_trims < x->score_within_trims) {
    return +1;
  } else if (x->nmatches_plus_spliced_trims > y->nmatches_plus_spliced_trims) {
    return -1;
  } else if (y->nmatches_plus_spliced_trims > x->nmatches_plus_spliced_trims) {
    return +1;
#if 0
  } else if (x->nmatches_to_trims > y->nmatches_to_trims) {
    return -1;
  } else if (y->nmatches_to_trims > x->nmatches_to_trims) {
    return +1;
#endif

#if 0
    /* Causes hits to not be recognized as equivalent */
  } else if (x->nsplices < y->nsplices) {
    return -1;
  } else if (y->nsplices < x->nsplices) {
    return +1;
#endif

#if 0
  } else if (y->start_amb_length + y->end_amb_length == 0 &&
	     x->start_amb_length + x->end_amb_length > 0) {
    return -1;
  } else if (x->start_amb_length + x->end_amb_length == 0 &&
	     y->start_amb_length + y->end_amb_length > 0) {
    return +1;
#endif

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

#if 0
    /* Used for sorting but not equiv */
  } else if (x->sensedir != 0 && y->sensedir == 0) {
    return -1;
  } else if (y->sensedir != 0 && x->sensedir == 0) {
    return +1;
#endif

#if 0
  } else if (x->sensedir == y->sensedir) {
    return 0;
  } else if (x->sensedir > y->sensedir) {
    return +1;
  } else if (y->sensedir > x->sensedir) {
    return -1;
#endif

  } else if (x->splice_score > y->splice_score) {
    debug4(printf(" => loses by splice score\n"));
    return -1;

  } else if (y->splice_score > x->splice_score) {
    debug4(printf(" => wins by splice score\n"));
    return +1;

  } else {
    debug4(printf(" => identical for sorting purposes\n"));
    return 0;
  }
}


int
Stage3end_hit_goodness_cmp (bool *equalp, Stage3end_T hit,
			    Stage3end_T best_hit, bool finalp) {
  double prob1, prob2;

#ifdef PRE_RESOLVE_MULTIMAPPING
  if (Stage3end_tally(x) > TALLY_RATIO*Stage3end_tally(y)) {
    debug4(printf("  #%d overlaps #%d and tally %ld > %f*%ld, so marking %d for elimination\n",
		  i,j,x->tally,TALLY_RATIO,y->tally,j));
    eliminate[j] = true;
  } else if (Stage3end_tally(y) > TALLY_RATIO*Stage3end_tally(x)) {
    debug4(printf("  #%d overlaps #%d and tally %f*%ld < %ld, so marking %d for elimination\n",
		  i,j,TALLY_RATIO,x->tally,y->tally,i));
    eliminate[i] = true;
  }
#endif

  *equalp = false;

  /* Favors definitive splices over ambiguous ones (by using nmatches_to_trims) */
  if (known_ambiguous_p(hit) == true && known_ambiguous_p(best_hit) == false) {
    return -1;

  } else if (known_ambiguous_p(hit) == false && known_ambiguous_p(best_hit) == true) {
    return +1;

  } else if (hit->nsegments > best_hit->nsegments) {
    if (hit->nmatches_plus_spliced_trims > best_hit->nmatches_plus_spliced_trims) {
      /* More segments and strictly more matches */
      debug4(printf("More segments and strictly more matches (to_trims)\n"));
      return +1;
    } else {
      /* More segments, but don't add anything */
      debug4(printf("More segments but don't add anything\n"));
      return -1;
    }

  } else if (hit->nsegments < best_hit->nsegments) {
    if (hit->nmatches_plus_spliced_trims >= best_hit->nmatches_plus_spliced_trims) {
      /* Fewer segments, but same or more matches */
      debug4(printf("Fewer segments and same or more matches (to_trims)\n"));
      return +1;
    } else {
      debug4(printf("Fewer segments and don't add anything\n"));
      /* Fewer segments, and don't add anything */
      return -1;
    }

#if 0
  } else if (hit->nmatches_to_trims < best_hit->nmatches_to_trims) {
    /* Favors longer alignments to potentially wrong splice sites */
    debug4(printf("  => loses by nmatches_to_trims\n"));
    return -1;

  } else if (hit->nmatches_to_trims > best_hit->nmatches_to_trims) {
    debug4(printf("  => wins by nmatches_to_trims\n"));
    return +1;
#endif

#if 0
  } else if (hit->nsplices > best_hit->nsplices) {
    debug4(printf("  => loses by nsplices: %d > %d in best\n",hit->nsplices,best_hit->nsplices));
    return -1;
  } else if (hit->nsplices < best_hit->nsplices) {
    debug4(printf("  => wins by nsplices: %d < %d in best\n",hit->nsplices,best_hit->nsplices));
    return +1;
#endif

  } else if (hit->hittype > best_hit->hittype) {
    debug4(printf("  => loses by hittype\n"));
    return -1;
  } else if (hit->hittype < best_hit->hittype) {
    debug4(printf("  => wins by hittype\n"));
    return +1;

#if 0
  } else if (start_amb_length(hit) + end_amb_length(hit) > 0 &&
	     start_amb_length(best_hit) + end_amb_length(best_hit) == 0) {
    debug4(printf("  => loses by ambiguity\n"));
    return -1;
  } else if (start_amb_length(hit) + end_amb_length(hit) == 0 &&
	     start_amb_length(best_hit) + end_amb_length(best_hit) > 0) {
    debug4(printf("  => wins by ambiguity\n"));
    return +1;
#endif

  } else if (hit->nindels > best_hit->nindels) {
    debug4(printf("  => loses by nindels\n"));
    return -1;
  } else if (hit->nindels < best_hit->nindels) {
    debug4(printf("  => wins by nindels\n"));
    return +1;

  } else if (hit->distant_splice_p == true && best_hit->distant_splice_p == false) {
    debug4(printf("  => loses because distant splice\n"));
    return -1;
  } else if (hit->distant_splice_p == false && best_hit->distant_splice_p == true) {
    debug4(printf("  => wins because not distant splice\n"));
    return +1;

  } else if (finalp == false) {
    debug4(printf("  => indistinguishable\n"));
    return 0;

  } else if (hit->hittype == TRANSLOC_SPLICE && best_hit->hittype == TRANSLOC_SPLICE) {
    prob1 = hit->splice_score;
    prob2 = best_hit->splice_score;

    if (prob1 < prob2) {
      debug4(printf("  => loses by TRANSLOC_SPLICE splice prob %f vs %f\n",prob1,prob2));
      return -1;
    } else if (prob1 > prob2) {
      debug4(printf("  => wins by TRANSLOC_SPLICE splice prob %f vs %f\n",prob1,prob2));
      return +1;
    } else {
      debug4(printf("  => equal\n"));
      *equalp = true;
      return 0;
    }

  } else {
    prob1 = Stage3end_prob(hit);
    prob2 = Stage3end_prob(best_hit);
    if (prob1 < prob2) {
      debug4(printf("  => loses by splice prob %f vs %f\n",prob1,prob2));
      return -1;
    } else if (prob1 > prob2) {
      debug4(printf("  => wins by splice prob %f vs %f\n",prob1,prob2));
      return +1;
    }

    if (hit->genomiclength > best_hit->genomiclength) {
      debug4(printf("  => loses by genomiclength\n"));
      return -1;
    } else if (hit->genomiclength < best_hit->genomiclength) {
      debug4(printf("  => wins by genomiclength\n"));
      return +1;

    } else {
      debug4(printf("  => equal\n"));
      *equalp = true;
      return 0;
    }
  }
}


/* Not clear how to handle altloc */
static bool
hit_subsumption (Stage3end_T x, Stage3end_T y) {
#if 0
  /* Should not be necessary if we compute chrnum correctly */
  if (x->chrnum != y->chrnum) {
    return false;		/* Different chromosomes, possibly due to straddles */
  }
#endif
    
  if (x->plusp != y->plusp) {
    return false;		/* Different strands */
  } else if (x->low <= y->low && x->high >= y->high) {
    return true;
  } else if (y->low <= x->low && y->high >= x->high) {
    return true;
  } else {
    return false;
  }
}

/* Not clear how to handle altloc */
static bool
hit_endpoint_equivp (Stage3end_T x, Stage3end_T y) {
  if (x->plusp != y->plusp) {
    return false;		/* Different strands */
  } else if (x->low != y->low) {
    return false;
  } else if (x->high != y->high) {
    return false;
  } else {
    return true;
  }
}


static bool
hit_bad_superstretch_p (Stage3end_T hit_k, Stage3end_T *hits, int k, int j, bool finalp) {
  int a;
  bool equalp;

  for (a = k+1; a <= j; a++) {
    if (hit_subsumption(hit_k,hits[a]) == true) {
      debug4(printf("Testing %d because stretches over %d",k,a));
      if (Stage3end_hit_goodness_cmp(&equalp,hits[a],hit_k,finalp) > 0 || equalp == true) {
	debug4(printf(" => eliminating\n"));
	return true;
      }
      debug4(printf("\n"));
    }
  }
  return false;
}


static List_T
remove_overlaps_distant (List_T hitlist, Hitlistpool_T hitlistpool) {
  List_T unique = NULL;
  T best_hit, hit, *hits;
  int cmp;
  int n, i, j, k, besti;
  bool *eliminate, equalp;
#ifdef PRE_RESOLVE_MULTIMAPPING
  long int best_tally;
#endif

  if ((n = List_length(hitlist)) == 0) {
    return (List_T) NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    hits = (T *) MALLOCA(n * sizeof(T));
    List_fill_array((void **) hits,hitlist);
    Hitlist_free(&hitlist);
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hits = (T *) List_to_array(hitlist,NULL);
    Hitlist_free(&hitlist);
#endif
  }

  debug4(printf("Step 0.  Checking for duplicates among distant\n"));
  qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp);

  /* Find clusters from left */
  i = 0;
  while (i < n) {
    j = i;
    while (j+1 < n && hit_endpoint_equivp(hits[i],hits[j+1]) == true) {
      j = j+1;
    }

    if (j > i) {
      debug4(printf("Cluster from %d up through %d\n",i,j));

      best_hit = hits[i];
      besti = i;
      debug4(printf("Assume best is %d\n",besti));

      for (k = i+1; k <= j; k++) {
	cmp = Stage3end_hit_goodness_cmp(&equalp,hits[k],best_hit,/*finalp*/true);
	debug4(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	if (cmp > 0) {
	  best_hit = hits[k];
	  besti = k;
	  debug4(printf("Best is now %d\n",besti));
	}
      }

      for (k = i; k <= j; k++) {
	if (k == besti) {
	  /* Skip */
	} else if (Stage3end_hit_goodness_cmp(&equalp,hits[k],best_hit,/*finalp*/true) < 0 || equalp == true) {
	  debug4(printf("  Eliminating hit %d from left, because beaten by %d\n",k,besti));
	  eliminate[k] = true;
	}
      }
    }
      
    i = j+1;
  }

  for (i = n-1; i >= 0; i--) {
    hit = hits[i];
    if (eliminate[i] == false) {
      unique = Hitlist_push(unique,hitlistpool,(void *) hit);
    } else if (hit->paired_usedp == true) {
      unique = Hitlist_push(unique,hitlistpool,(void *) hit);
    } else {
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hits);
  FREEA(eliminate);
#else
  FREE(hits);
  FREE(eliminate);
#endif

  debug4(printf("Returning %d unique distant splices\n",List_length(unique)));
  return unique;
}




List_T
Stage3end_remove_overlaps (List_T hitlist, Hitlistpool_T hitlistpool, bool finalp) {
  List_T unique = NULL, distant = NULL, local = NULL, p;
  T best_hit, hit, parent, *hits, *prev;
  int cmp;
  int nkept, n, i, j, k, besti;
  bool *eliminate, equalp;
  int *parenti;
#ifdef PRE_RESOLVE_MULTIMAPPING
  long int best_tally;
#endif

  
  debug4(printf("Entered Stage3end_remove_overlaps with %d hits: %s\n",
		List_length(hitlist),finalp == true ? "FINAL" : "not final"));

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (T) List_head(p);
    if (hit->distant_splice_p == false) {
      local = Hitlist_push(local,hitlistpool,(void *) hit);
    } else {
      distant = Hitlist_push(distant,hitlistpool,(void *) hit);
    }
  }
  Hitlist_free(&hitlist);
  
  distant = remove_overlaps_distant(distant,hitlistpool);

  if ((n = List_length(local)) == 0) {
    return distant;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    hits = (T *) MALLOCA(n * sizeof(T));
    List_fill_array((void **) hits,local);
    Hitlist_free(&local);
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hits = (T *) List_to_array(local,NULL);
    Hitlist_free(&local);
#endif
  }


  /* Step 1.  Check for exact duplicates */
  /* Probably don't want to eliminate aliases at this point */
  debug4(printf("Step 1.  Checking for exact duplicates\n"));
  qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp);

  debug4(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, circularalias %d, nmatches %d (%d to_trims), score %d",
		  i,Method_string(hit->method),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->circularalias,hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->score_within_trims);
	   if (hit->transcripts != NULL) {
	     Transcript_print_list(hit->transcripts);
	   }
	   printf("\n");
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    debug4(printf(" %d,%d",i,j));
    while (j < n && hit_equiv_cmp(hits[j],hits[i]) == 0) {
      debug4(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
  }
  debug4(printf("\n"));


  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
#ifdef USE_ALLOCA_FOR_HITS
  hits = (Stage3end_T *) MALLOCA(nkept * sizeof(Stage3end_T));
#else
  hits = (Stage3end_T *) MALLOC(nkept * sizeof(Stage3end_T));
#endif

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug4(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      best_hit = hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug4(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      hits[j++] = hit;
    } else {
      debug4(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      Stage3end_transfer_transcripts_one(best_hit,hit);
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(prev);
#else
  FREE(prev);
#endif


  /* Step 2: Check for superstretches */
  n = nkept;
  debug4(printf("Step 2.  Checking for superstretches among %d hits within subsumption clusters\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }

  debug4(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches %d (%d to_trims), score %d",
		  i,Method_string(hit->method),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->score_within_trims);
	   if (hit->transcripts != NULL) {
	     Transcript_print_list(hit->transcripts);
	   }
	   printf("\n");
	 }
	 );

  /* Find clusters */
  i = 0;
  while (i < n) {
    j = i;
    /* Previously checked if (hits[i]->distant_splice_p == false) */
    while (j+1 < n && hit_subsumption(hits[i],hits[j+1]) == true) {
      j = j+1;
    }

    if (j > i) {
      debug4(printf("Cluster from %d up through %d\n",i,j));

      /* Find bad superstretches */
      for (k = i; k <= j; k++) {
	/* Previously checked if (hits[i]->distant_splice_p == false) */
	if (hit_bad_superstretch_p(hits[k],hits,k,j,finalp) == true) {
	  eliminate[k] = true;
	  /* parenti[k] = j; */
	}
      }
    }

    i = j+1;
  }

  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
#ifdef USE_ALLOCA_FOR_HITS
  hits = (Stage3end_T *) MALLOCA(nkept * sizeof(Stage3end_T));
#else
  hits = (Stage3end_T *) MALLOC(nkept * sizeof(Stage3end_T));
#endif

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug4(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug4(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      hits[j++] = hit;
    } else {
      debug4(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      /* parent = prev[parenti[i]]; */
      /* Stage3end_transfer_transcripts_one(parent,hit); */
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(prev);
#else
  FREE(prev);
#endif


  /* Step 3: Check for best within subsumption clusters */
  n = nkept;
  debug4(printf("Checking for best among %d hits within subsumption clusters\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp); -- No need since original order was kept */
  
  debug4(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches %d (%d to_trims), score %d",
		  i,Method_string(hit->method),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->score_within_trims);
	   if (hit->transcripts != NULL) {
	     Transcript_print_list(hit->transcripts);
	   }
	   printf("\n");
	 }
	 );

  /* Find clusters from left */
  i = 0;
  while (i < n) {
    j = i;
    /* Previously checked if (hits[i]->distant_splice_p == false) */
    while (j+1 < n && hit_subsumption(hits[i],hits[j+1]) == true) {
      j = j+1;
    }

    if (j > i) {
      debug4(printf("Cluster from %d up through %d\n",i,j));

      best_hit = hits[i];
      besti = i;
      debug4(printf("Assume best is %d\n",besti));

      for (k = i+1; k <= j; k++) {
	/* Previously checked if (hits[i]->distant_splice_p == false) */
	cmp = Stage3end_hit_goodness_cmp(&equalp,hits[k],best_hit,finalp);
	debug4(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	if (cmp > 0) {
	  best_hit = hits[k];
	  besti = k;
	  debug4(printf("Best is now %d\n",besti));
	}
      }

      for (k = i; k <= j; k++) {
	if (k == besti) {
	  /* Skip */
	  /* Previously checked if (hits[i]->distant_splice_p == false) */
	} else if (Stage3end_hit_goodness_cmp(&equalp,hits[k],best_hit,finalp) < 0 || equalp == true) {
	  debug4(printf("  Eliminating hit %d from left, because beaten by %d\n",k,besti));
	  eliminate[k] = true;
	  /* parenti[k] = i; */
	}
      }
    }
      
    i = j+1;
  }


  /* Find clusters starting from right */
  j = n - 1;
  while (j >= 0) {
    i = j;
    /* Previously checked if (hits[i]->distant_splice_p == false) */
    while (i-1 >= 0 && hit_subsumption(hits[j],hits[i-1]) == true) {
      i = i-1;
    }

    if (i < j) {
      debug4(printf("Cluster from %d down through %d\n",j,i));
      best_hit = hits[i];
      besti = i;
      debug4(printf("Assume best is %d\n",besti));

      for (k = i+1; k <= j; k++) {
	/* Previously checked if (hits[i]->distant_splice_p == false) */
	cmp = Stage3end_hit_goodness_cmp(&equalp,hits[k],best_hit,finalp);
	debug4(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	if (cmp > 0) {
	  best_hit = hits[k];
	  besti = k;
	  debug4(printf("Best is now %d\n",besti));
	}
      }

      for (k = i; k <= j; k++) {
	if (k == besti) {
	  /* Skip */
	  /* Previously checked if (hits[i]->distant_splice_p == false) */
	} else if (Stage3end_hit_goodness_cmp(&equalp,hits[k],best_hit,finalp) < 0 || equalp == true) {
	  debug4(printf("  Eliminating hit %d from right, because beaten by %d\n",k,besti));
	  eliminate[k] = true;
	  /* parenti[k] = i; */
	}
      }
    }
      
    j = i-1;
  }


  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
#ifdef USE_ALLOCA_FOR_HITS
  hits = (Stage3end_T *) MALLOCA(nkept * sizeof(Stage3end_T));
#else
  hits = (Stage3end_T *) MALLOC(nkept * sizeof(Stage3end_T));
#endif

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug4(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug4(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      hits[j++] = hit;
    } else {
      debug4(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d, sensedir = %d)\n",
		    hit->low - hit->chroffset,hit->high - hit->chroffset,hit->nmatches_plus_spliced_trims,
		    hit->plusp,hit->sensedir));
      /* parent = prev[parenti[i]]; */
      /* Stage3end_transfer_transcripts_one(parent,hit); */
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(prev);
  parenti = (int *) CALLOCA(nkept,sizeof(int));
#else
  FREE(prev);
  parenti = (int *) CALLOC(nkept,sizeof(int));
#endif


  /* Step 4: Check for identity */
  n = nkept;
  debug4(printf("Checking for duplicates among %d hits by identity\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp); -- No need since original order was kept */

  debug4(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches %d (%d to_trims), score %d",
		  i,Method_string(hit->method),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches_plus_spliced_trims,hit->nmatches_to_trims,hit->score_within_trims);
	   if (hit->transcripts != NULL) {
	     Transcript_print_list(hit->transcripts);
	   }
	   printf("\n");
	 }
	 );

  i = 0;
  while (i < n) {
    debug4(printf("Looking at %d with score %d\n",i,hits[i]->score_within_trims));
    j = i+1;
    while (j < n && hit_equiv_cmp(hits[j],hits[i]) == 0) {
      debug4(printf("  %d equal to %d\n",j,i));
      eliminate[j] = true;
      parenti[j] = i;
      j++;
    }

    i = j;
  }

  for (i = n-1; i >= 0; i--) {
    hit = hits[i];
    if (eliminate[i] == false) {
      unique = Hitlist_push(unique,hitlistpool,(void *) hit);
    } else if (hit->paired_usedp == true) {
      unique = Hitlist_push(unique,hitlistpool,(void *) hit);
    } else {
      parent = hits[parenti[i]]; /* Not prev, since we are using hits instead */
      Stage3end_transfer_transcripts_one(parent,hit);
      Stage3end_free(&hit);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hits);
  FREEA(eliminate);
  FREEA(parenti);
#else
  FREE(hits);
  FREE(eliminate);
  FREE(parenti);
#endif


#ifdef PRE_RESOLVE_MULTIMAPPING
  if (use_tally_p == true && tally_iit != NULL) {
    if ((n = List_length(unique)) > 1) {
#ifdef USE_ALLOCA_FOR_HITS
      hits = (T *) MALLOCA(n * sizeof(T));
      List_fill_array((void **) hits,unique);
      Hitlist_free(&unique);
#else
      hits = (T *) List_to_array(unique,NULL);
      Hitlist_free(&unique);
#endif

      best_tally = 0;
      for (i = 0; i < n; i++) {
	if (hits[i]->tally < 0) {
	  hits[i]->tally = Stage3end_compute_tally(hits[i]);
	}
	if (hits[i]->tally > best_tally) {
	  best_tally = hits[i]->tally;
	}
      }

      unique = (List_T) NULL;
      for (i = 0; i < n; i++) {
	if (hits[i]->tally < best_tally) {
	  /* Stage3end_free(&(hits[i])); */
	} else {
	  unique = Hitlist_push(unique,hitlistpool,(void *) hits[i]);
	}
      }

#ifdef USE_ALLOCA_FOR_HITS
      FREEA(hits);
#else
      FREE(hits);
#endif
    }
  }
#endif

  unique = List_append(unique,distant);
  debug4(printf("Exited Stage3end_remove_overlaps with %d hits\n",List_length(unique)));
  return unique;
}


List_T
Stage3end_resolve_multimapping (List_T hits, Hitlistpool_T hitlistpool) {
  List_T resolve1, resolve2, resolve3, p;
  Stage3end_T hit;

  /* Overlap_T best_overlap; */
  long int best_tally;
  double tally_threshold;
  bool runlengthp;

  if (List_length(hits) <= 1) {
    return hits;
  }

  resolve1 = hits;

  if (tally_iit == NULL) {
    resolve2 = resolve1;
  } else {
    best_tally = 0L;
    for (p = resolve1; p != NULL; p = p->rest) {
      hit = (Stage3end_T) p->first;
      if ((hit->tally = Stage3end_compute_tally(hit)) > best_tally) {
	best_tally = hit->tally;
      }
    }
    if (best_tally == 0L) {
      resolve2 = resolve1;
    } else {
      resolve2 = (List_T) NULL;
#ifdef USE_TALLY_RATIO
      tally_threshold = (double) best_tally / TALLY_RATIO;
#else
      tally_threshold = 1.0;
#endif
      for (p = resolve1; p != NULL; p = p->rest) {
	hit = (Stage3end_T) p->first;
	if ((double) hit->tally < tally_threshold) {
	  Stage3end_free(&hit);
	} else {
	  resolve2 = Hitlist_push(resolve2,hitlistpool,(void *) hit);
	}
      }
      Hitlist_free(&resolve1);
    }
  }


  if (List_length(resolve2) <= 1) {
    return resolve2;
  }

  if (runlength_iit == NULL) {
    resolve3 = resolve2;
  } else {
    runlengthp = false;
    for (p = resolve2; p != NULL; p = p->rest) {
      hit = (Stage3end_T) p->first;
      if (Stage3end_runlength_p(hit) == true) {
	runlengthp = true;
      }
    }
    if (runlengthp == false) {
      resolve3 = resolve2;
    } else {
      resolve3 = (List_T) NULL;
      for (p = resolve2; p != NULL; p = p->rest) {
	hit = (Stage3end_T) p->first;
	if (Stage3end_runlength_p(hit) == false) {
	  Stage3end_free(&hit);
	} else {
	  resolve3 = Hitlist_push(resolve3,hitlistpool,(void *) hit);
	}
      }
      Hitlist_free(&resolve2);
    }
  }


  return resolve3;
}


Pairtype_T
Stage3_determine_pairtype (T hit5, T hit3, Stage3pair_T stage3pair) {
  int pairmax;

  debug14(printf("Entered Stage3_determine_pairtype\n"));
  if (hit5->effective_chrnum != hit3->effective_chrnum) {
    debug14(printf("Returning unpaired\n"));
    return UNPAIRED;
  } else if (hit5->plusp != hit3->plusp) {
    debug14(printf("Returning paired_inversion\n"));
    return PAIRED_INVERSION;
  } else if (hit5->plusp == true) {
    if (hit3->genomicend < hit5->genomicstart) {
      debug14(printf("Returning paired_scramble\n"));
      return PAIRED_SCRAMBLE;
    } else {
      if (circularp[hit5->effective_chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }
      if (stage3pair != NULL && stage3pair->transcripts5 != NULL && stage3pair->transcripts3 != NULL) {
	debug14(printf("Returning concordant based on transcriptome\n"));
	return CONCORDANT;
      } else if (hit3->genomicstart > hit5->genomicend + pairmax) {
	debug14(printf("Returning paired_toolong\n"));
	return PAIRED_TOOLONG;
      } else {
	debug14(printf("Returning concordant\n"));
	return CONCORDANT;
      }
    }
  } else {
    if (hit3->genomicend > hit5->genomicstart) {
      debug14(printf("Returning paired_scramble\n"));
      return PAIRED_SCRAMBLE;
    } else {
      if (circularp[hit3->effective_chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }
      if (stage3pair != NULL && stage3pair->transcripts5 != NULL && stage3pair->transcripts3 != NULL) {
	debug14(printf("Returning concordant based on transcriptome\n"));
	return CONCORDANT;
      } else if (hit3->genomicstart + pairmax < hit5->genomicend) {
	debug14(printf("Returning paired_toolong\n"));
	return PAIRED_TOOLONG;
      } else {
	debug14(printf("Returning concordant\n"));
	return CONCORDANT;
      }
    }
  }
}


#if 0
/* Previously, samprint.c called this, but it can lead to incorrect answers when transcripts are added later */
Pairtype_T
Stage3pair_pairtype (Stage3pair_T this) {
  return this->pairtype;
}
#else
Pairtype_T
Stage3pair_determine_pairtype (Stage3pair_T this) {
  return Stage3_determine_pairtype(this->hit5,this->hit3,this);
}
#endif

bool
Stage3pair_circularp (Stage3pair_T this) {
  return this->circularp;
}

bool
Stage3pair_altlocp (Stage3pair_T this) {
  if (altlocp[this->hit5->chrnum] == true) {
    return true;
  } else if (altlocp[this->hit3->chrnum] == true) {
    return true;
  } else {
    return false;
  }
}


#if 0
static char *
unpaired_type_text (T hit5, T hit3) {
  if (hit5->chrnum != hit3->chrnum) {
    return UNPAIRED_INTERCHROM_TEXT;
  } else if (hit5->plusp != hit3->plusp) {
    return PAIRED_INVERSION_TEXT;
  } else if (hit5->plusp == true) {
    if (hit3->genomicstart < hit5->genomicstart) {
      return PAIRED_SCRAMBLE_TEXT;
    } else {
      return UNPAIRED_TOOLONG_TEXT;
    }
  } else {
    if (hit5->genomicstart < hit3->genomicstart) {
      return PAIRED_SCRAMBLE_TEXT;
    } else {
      return UNPAIRED_TOOLONG_TEXT;
    }
  }
}
#endif




/* Has a copy in pair.c */
static void
print_pair_info (Filestring_T fp, T hit5, T hit3, int insertlength, int pairscore,
		 Pairtype_T pairtype) {

  assert(hit5->effective_chrnum == hit3->effective_chrnum); /* Same chromosomes */

#if 0
  /* Doesn't hold for paired (inversion) */
  assert(hit5->plusp == hit3->plusp);	/* Same direction */
#endif

#ifndef NO_COMPARE
  FPRINTF(fp,"pair_score:%d",pairscore);
  FPRINTF(fp,",insert_length:%d",insertlength);
#endif

  switch (pairtype) {
  case CONCORDANT: break;
  case PAIRED_SCRAMBLE: FPRINTF(fp,",pairtype:scramble"); break;
  case PAIRED_INVERSION: FPRINTF(fp,",pairtype:inversion"); break;
  case PAIRED_TOOLONG: FPRINTF(fp,",pairtype:toolong"); break;
  case CONCORDANT_TRANSLOCATIONS: break;
  case PAIRED_UNSPECIFIED: abort();
  case UNPAIRED: abort();
  case UNSPECIFIED: abort();
  }

  return;
}




static void
print_substrings (Filestring_T fp, Stage3pair_T stage3pair, T this,
		  int score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
		  Shortread_T headerseq, char *acc_suffix, bool invertp, T hit5, T hit3, int insertlength,
		  int pairscore, Pairtype_T pairtype, int mapq_score, bool first_read_p) {
  char *single_chr, *chr;
  bool allocp, alloc1p, pairinfo_printed_p = false;
  List_T substrings, junctions, p, q;
  Substring_T substring;
  Junction_T pre_junction, post_junction;
  int nblocks;

  if (this->chrnum == 0) {
    single_chr = (char *) NULL;
    alloc1p = false;
  } else {
    single_chr = Univ_IIT_label(chromosome_iit,this->chrnum,&alloc1p);
  }
  if (invertp == true) {
    substrings = this->substrings_Nto1;
    junctions = this->junctions_Nto1;
  } else {
    substrings = this->substrings_1toN;
    junctions = this->junctions_1toN;
  }

  if (output_type == M8_OUTPUT) {
    for (p = substrings; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      if (Substring_has_alts_p(substring) == true) {
	/* Skip */
      } else {
	if ((chr = single_chr) == NULL) {
	  chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
	}
	Substring_print_m8(fp,substring,headerseq,acc_suffix,chr,invertp);
	if (single_chr == NULL && allocp == true) {
	  FREE(chr);
	}
      }
    }

  } else {
    if ((nblocks = List_length(substrings)) == 1) {
      post_junction = (Junction_T) NULL;
    } else {
      post_junction = (Junction_T) List_head(junctions);
    }
    substring = (Substring_T) List_head(substrings);
    if (Substring_has_alts_p(substring) == true) {
      nblocks -= 1;
    }
    substring = (Substring_T) List_last_value(substrings);
    if (Substring_has_alts_p(substring) == true) {
      nblocks -= 1;
    }


    /* First line */
    substring = (Substring_T) List_head(substrings);
    if (Substring_has_alts_p(substring) == true) {
      /* Skip */
    } else {
      if ((chr = single_chr) == NULL) {
	chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
      }
      FPRINTF(fp," ");
      Substring_print_alignment(fp,/*pre_junction*/NULL,substring,post_junction,queryseq,genomecomp,chr,invertp);
      if (single_chr == NULL && allocp == true) {
	FREE(chr);
      }

      /* Alignment info */
#ifndef NO_COMPARE      
      FPRINTF(fp,"\tsegs:%d,align_score:%d,mapq:%d",nblocks,score,mapq_score);
#endif
      if (method_print_p == true) {      
	Method_print(fp,this->method);
      }
    
      /* Transcriptome info */
      if (stage3pair != NULL) {
	if (first_read_p == true) {
	  if (stage3pair->transcripts5 != NULL) {
	    FPRINTF(fp,"\t");
	    FPRINTF(fp,"Transcripts:");
	    Transcript_print_info(fp,stage3pair->transcripts5,transcript_iit,invertp);
	    if (this->transcripts_other != NULL) {
	      FPRINTF(fp,"\tOther:");
	      Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
	    }
	  }
	} else {
	  if (stage3pair->transcripts3 != NULL) {
	    FPRINTF(fp,"\t");
	    FPRINTF(fp,"Transcripts:");
	    Transcript_print_info(fp,stage3pair->transcripts3,transcript_iit,invertp);
	    if (this->transcripts_other != NULL) {
	      FPRINTF(fp,"\tOther:");
	      Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
	    }
	  }
	}

      } else if (this->method == TR || this->transcripts != NULL) {
	FPRINTF(fp,"\t");
	FPRINTF(fp,"Transcripts:");
	Transcript_print_info(fp,this->transcripts,transcript_iit,invertp);
	if (this->transcripts_other != NULL) {
	  FPRINTF(fp,"\tOther:");
	  Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
	}
      }

      /* Pairing info */
      if (hit5 != NULL && hit3 != NULL) {
	FPRINTF(fp,"\t");
	print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
      pairinfo_printed_p = true;

      FPRINTF(fp,"\n");
    }

    if ((p = List_next(substrings)) == NULL) {
      /* Done */
    } else {
      /* Middle lines */
      for (q = List_next(junctions); q != NULL; p = List_next(p), q = List_next(q)) {
	pre_junction = post_junction;
	post_junction = List_head(q);

	substring = (Substring_T) List_head(p);
	if (Substring_has_alts_p(substring) == true) {
	  /* Skip */
	} else {
	  if (pairinfo_printed_p == true) {
	    FPRINTF(fp,",");
	  } else {
	    FPRINTF(fp," ");
	  }
	  if ((chr = single_chr) == NULL) {
	    chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
	  }
	  Substring_print_alignment(fp,pre_junction,substring,post_junction,queryseq,genomecomp,chr,invertp);
	  if (single_chr == NULL && allocp == true) {
	    FREE(chr);
	  }

	  if (pairinfo_printed_p == false) {
	    /* Alignment info */
#ifndef NO_COMPARE
	    FPRINTF(fp,"\tsegs:%d,align_score:%d,mapq:%d",nblocks,score,mapq_score);
#endif
	    if (method_print_p == true) {
	      Method_print(fp,this->method);
	    }
    
	    /* Transcriptome info */
	      if (stage3pair != NULL) {
		if (first_read_p == true) {
		  if (stage3pair->transcripts5 != NULL) {
		    FPRINTF(fp,"\t");
		    FPRINTF(fp,"Transcripts:");
		    Transcript_print_info(fp,stage3pair->transcripts5,transcript_iit,invertp);
		    if (this->transcripts_other != NULL) {
		      FPRINTF(fp,"\tOther:");
		      Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
		    }
		  }
		} else {
		  if (stage3pair->transcripts3 != NULL) {
		    FPRINTF(fp,"\t");
		    FPRINTF(fp,"Transcripts:");
		    Transcript_print_info(fp,stage3pair->transcripts3,transcript_iit,invertp);
		    if (this->transcripts_other != NULL) {
		      FPRINTF(fp,"\tOther:");
		      Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
		    }
		  }
		}

	      } else if (this->method == TR || this->transcripts != NULL) {
		FPRINTF(fp,"\t");
		FPRINTF(fp,"Transcripts:");
		Transcript_print_info(fp,this->transcripts,transcript_iit,invertp);
		if (this->transcripts_other != NULL) {
		  FPRINTF(fp,"\tOther:");
		  Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
		}
	      }

	    /* Pairing info */
	    if (hit5 != NULL && hit3 != NULL) {
	      FPRINTF(fp,"\t");
	      print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
	    }
	    pairinfo_printed_p = true;
	  }
	  FPRINTF(fp,"\n");
	}
      }

      /* Last line */
      pre_junction = post_junction;

      substring = (Substring_T) List_head(p);
      if (Substring_has_alts_p(substring) == true) {
	/* Skip */
      } else {
	if (pairinfo_printed_p == true) {
	  FPRINTF(fp,",");
	} else {
	  FPRINTF(fp," ");
	}
	if ((chr = single_chr) == NULL) {
	  chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
	}
	Substring_print_alignment(fp,pre_junction,substring,/*post_junction*/NULL,queryseq,genomecomp,chr,invertp);
	if (single_chr == NULL && allocp == true) {
	  FREE(chr);
	}

	if (pairinfo_printed_p == false) {
	  /* Alignment info */
#ifndef NO_COMPARE
	  FPRINTF(fp,"\tsegs:%d,align_score:%d,mapq:%d",nblocks,score,mapq_score);
#endif
	  if (method_print_p == true) {
	    Method_print(fp,this->method);
	  }
	  
	  /* Transcriptome info */
	  if (stage3pair != NULL) {
	    if (first_read_p == true) {
	      if (stage3pair->transcripts5 != NULL) {
		FPRINTF(fp,"\t");
		FPRINTF(fp,"Transcripts:");
		Transcript_print_info(fp,stage3pair->transcripts5,transcript_iit,invertp);
		if (this->transcripts_other != NULL) {
		  FPRINTF(fp,"\tOther:");
		  Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
		}
	      }
	    } else {
	      if (stage3pair->transcripts3 != NULL) {
		FPRINTF(fp,"\t");
		FPRINTF(fp,"Transcripts:");
		Transcript_print_info(fp,stage3pair->transcripts3,transcript_iit,invertp);
		if (this->transcripts_other != NULL) {
		  FPRINTF(fp,"\tOther:");
		  Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
		}
	      }
	    }
	    
	  } else if (this->method == TR || this->transcripts != NULL) {
	    FPRINTF(fp,"\t");
	    FPRINTF(fp,"Transcripts:");
	    Transcript_print_info(fp,this->transcripts,transcript_iit,invertp);
	    if (this->transcripts_other != NULL) {
	      FPRINTF(fp,"\tOther:");
	      Transcript_print_info(fp,this->transcripts_other,transcript_iit,invertp);
	    }
	  }
	  
	  /* Pairing info */
	  if (hit5 != NULL && hit3 != NULL) {
	    FPRINTF(fp,"\t");
	    print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
	  }
	  /* pairinfo_printed_p = true; */
	}
	FPRINTF(fp,"\n");
      }
    }
  }

  if (alloc1p == true) {
    FREE(single_chr);
  }
}



/* May substitute paired-end loglik for single-end loglik */
void
Stage3end_print (Filestring_T fp, Stage3pair_T stage3pair, T this,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, Shortread_T headerseq,
		 char *acc_suffix, bool invertp, T hit5, T hit3, int insertlength,
		 int pairscore, Pairtype_T pairtype, int mapq_score, bool first_read_p) {

  /* TODO: Instead of score_within_trims, which contains penalties for
     ambiguous lengths, use (querylength - this->nmatches_plus_spliced_trims) instead */
  print_substrings(fp,stage3pair,this,this->score_within_trims,
		   chromosome_iit,queryseq,headerseq,acc_suffix,invertp,
		   hit5,hit3,insertlength,pairscore,pairtype,mapq_score,first_read_p);

  return;
}


static void
print_query_header (Filestring_T fp, char initchar, Shortread_T queryseq, bool invertp) {
  FPRINTF(fp,"%c",initchar);
  if (invertp == false) {
    Shortread_print_oneline(fp,queryseq);
  } else {
    Shortread_print_oneline_revcomp(fp,queryseq);
  }

  return;
}



static void
print_barcode_and_quality (Filestring_T fp, Shortread_T queryseq, bool invertp, int quality_shift) {
  char *barcode;

  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    FPRINTF(fp,"\tbarcode:%s",barcode);
  }

  if (Shortread_quality_string(queryseq) != NULL) {
    FPRINTF(fp,"\t");
    if (invertp == false) {
      Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			      quality_shift,/*show_chopped_p*/true);
    } else {
      Shortread_print_quality_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				      quality_shift,/*show_chopped_p*/true);
    }
  }

  return;
}


void
Stage3pair_print_end (Filestring_T fp, Filestring_T fp_failedinput,
		      Result_T result, Resulttype_T resulttype,
		      char initchar, bool firstp, Univ_IIT_T chromosome_iit,
		      Shortread_T queryseq, Shortread_T headerseq1, Shortread_T headerseq2,
		      int maxpaths, bool quiet_if_excessive_p,
		      bool invertp, int quality_shift) {
  Stage3pair_T *stage3pairarray, stage3pair;
  T *stage3array, *stage3array_mate, this, hit5, hit3;
  int npaths_primary, npaths_altloc, npaths_mate_primary, npaths_mate_altloc, pathnum;
  int first_absmq, second_absmq;
  bool excessivep, translocationp;

  if (resulttype == PAIREDEND_NOMAPPING) {
    if (output_type != M8_OUTPUT) {
      Filestring_set_split_output(fp,OUTPUT_NM);
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t0 %s",UNPAIRED_TEXT);

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
    
      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
      FPRINTF(fp,"\n");
    }
    /* If failedinput_root != NULL, then this case is handled by calling procedure */

  } else if (resulttype == CONCORDANT_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    stage3pair = stage3pairarray[0];
    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;

    if (stage3pair->circularp == true) {
      Filestring_set_split_output(fp,OUTPUT_CC);
    } else {
      Filestring_set_split_output(fp,OUTPUT_CU);
    }

    if (omit_concordant_uniq_p == true && stage3pair->circularp == false) {
      /* Skip printing */
      Filestring_set_split_output(fp,OUTPUT_NONE);

    } else {
      if (output_type != M8_OUTPUT) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t1 %s",CONCORDANT_TEXT);
    
	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);
      }
    
      if (firstp == true) {
	Stage3end_print(fp,stage3pair,hit5,
			chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			invertp,hit5,hit3,stage3pair->insertlength,
			stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
			/*first_read_p*/true);
      } else {
	Stage3end_print(fp,stage3pair,hit3,
			chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			invertp,hit5,hit3,stage3pair->insertlength,
			stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
			/*first_read_p*/false);
      }

      if (output_type != M8_OUTPUT) {
	FPRINTF(fp,"\n");
      }
    }

  } else if (resulttype == CONCORDANT_TRANSLOC) {
    Filestring_set_split_output(fp,OUTPUT_CT);
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      if (output_type != M8_OUTPUT) {
	/* No xs category for transloc, so ignore quiet-if-excessive_p */
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,CONCORDANT_TEXT);
	FPRINTF(fp," (transloc)");
	
	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
      
	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);

	/* No further output */
	FPRINTF(fp,"\n");
      }

      if (failedinput_root != NULL) {
	Shortread_print_query_singleend(fp_failedinput,queryseq,headerseq1);
      }

    } else {
      if (output_type != M8_OUTPUT) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,CONCORDANT_TEXT);
	FPRINTF(fp," (transloc)");

	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
      
	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);
      }

      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (firstp == true) {
	  Stage3end_print(fp,stage3pair,hit5,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
			  /*first_read_p*/true);
	} else {
	  Stage3end_print(fp,stage3pair,hit3,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
			  /*first_read_p*/false);
	}
      }

      if (output_type != M8_OUTPUT) {
	FPRINTF(fp,"\n");
      }
    }


  } else if (resulttype == CONCORDANT_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (omit_concordant_mult_p == true) {
      /* Skip printing */
      Filestring_set_split_output(fp,OUTPUT_NONE);

    } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      Filestring_set_split_output(fp,OUTPUT_CX);
      if (output_type != M8_OUTPUT) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,CONCORDANT_TEXT);
	
	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);

	/* No further output */
	FPRINTF(fp,"\n");
	
	if (failedinput_root != NULL) {
	  Shortread_print_query_singleend(fp_failedinput,queryseq,headerseq1);
	}
      }

    } else {
      Filestring_set_split_output(fp,OUTPUT_CM);
      if (output_type != M8_OUTPUT) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,CONCORDANT_TEXT);
	
	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
	
	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);
      }

      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;
	
	if (firstp == true) {
	  Stage3end_print(fp,stage3pair,hit5,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
			  /*first_read_p*/true);
	} else {
	  Stage3end_print(fp,stage3pair,hit3,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
			  /*first_read_p*/false);
	}
      }
      
      if (output_type != M8_OUTPUT) {
	FPRINTF(fp,"\n");
      }
    }

  } else if (resulttype == PAIRED_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    stage3pair = stage3pairarray[0];

    if (stage3pair->circularp == true) {
      Filestring_set_split_output(fp,OUTPUT_PC);
    } else if (stage3pair->pairtype == PAIRED_INVERSION) {
      Filestring_set_split_output(fp,OUTPUT_PI);
    } else if (stage3pair->pairtype == PAIRED_SCRAMBLE) {
      Filestring_set_split_output(fp,OUTPUT_PS);
    } else if (stage3pair->pairtype == PAIRED_TOOLONG) {
      Filestring_set_split_output(fp,OUTPUT_PL);
    } else {
      fprintf(stderr,"Unexpected pairtype %d\n",stage3pair->pairtype);
      abort();
    }
    
    if (output_type != M8_OUTPUT) {
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t1 %s",PAIRED_TEXT);

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
    }

    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;
    
    if (firstp == true) {
      Stage3end_print(fp,stage3pair,hit5,
		      chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
		      invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
		      /*first_read_p*/true);
    } else {
      Stage3end_print(fp,stage3pair,hit3,
		      chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
		      invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
		      /*first_read_p*/false);
    }

    if (output_type != M8_OUTPUT) {
      FPRINTF(fp,"\n");
    }

  } else if (resulttype == PAIRED_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      Filestring_set_split_output(fp,OUTPUT_PX);
      if (output_type != M8_OUTPUT) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,PAIRED_TEXT);
	
	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);

	/* No further output */
	FPRINTF(fp,"\n");

	if (failedinput_root != NULL) {
	  Shortread_print_query_singleend(fp_failedinput,queryseq,headerseq1);
	}
      }

    } else {
      Filestring_set_split_output(fp,OUTPUT_PM);
      if (output_type != M8_OUTPUT) {
	print_query_header(fp,initchar,queryseq,invertp);
	FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,PAIRED_TEXT);

	print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

	FPRINTF(fp,"\t");
	Shortread_print_header(fp,headerseq1,headerseq2);
      }

      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (firstp == true) {
	  Stage3end_print(fp,stage3pair,hit5,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
			  /*first_read_p*/true);
	} else {
	  Stage3end_print(fp,stage3pair,hit3,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			  invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score_within_trims,stage3pair->pairtype,stage3pair->mapq_score,
			  /*first_read_p*/false);
	}
      }

      if (output_type != M8_OUTPUT) {
	FPRINTF(fp,"\n");
      }
    }


  } else {
    /* Print as singles */
    if (firstp == true) {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      stage3array_mate = (T *) Result_array2(&npaths_mate_primary,&npaths_mate_altloc,&first_absmq,&second_absmq,result);
      stage3array = (T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    } else {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      stage3array_mate = (T *) Result_array(&npaths_mate_primary,&npaths_mate_altloc,&first_absmq,&second_absmq,result);
      stage3array = (T *) Result_array2(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    }

    excessivep = false;
    translocationp = false;
    if (resulttype == HALFMAPPING_UNIQ) {
      if (npaths_primary + npaths_altloc > 0 && Stage3end_circularpos(stage3array[0]) > 0) {
	Filestring_set_split_output(fp,OUTPUT_HC);
      } else if (npaths_mate_primary + npaths_mate_altloc > 0 && Stage3end_circularpos(stage3array_mate[0]) > 0) {
	Filestring_set_split_output(fp,OUTPUT_HC);
      } else {
	Filestring_set_split_output(fp,OUTPUT_HU);
      }

    } else if (resulttype == HALFMAPPING_TRANSLOC) {
      Filestring_set_split_output(fp,OUTPUT_HT);
      translocationp = true;

    } else if (resulttype == HALFMAPPING_MULT) {
      if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
	Filestring_set_split_output(fp,OUTPUT_HX);
	excessivep = true;
      } else {
	Filestring_set_split_output(fp,OUTPUT_HM);
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      if (npaths_primary + npaths_altloc > 0 && Stage3end_circularpos(stage3array[0]) > 0) {
	Filestring_set_split_output(fp,OUTPUT_UC);
      } else if (npaths_mate_primary + npaths_mate_altloc > 0 && Stage3end_circularpos(stage3array_mate[0]) > 0) {
	Filestring_set_split_output(fp,OUTPUT_UC);
      } else {
	Filestring_set_split_output(fp,OUTPUT_UU);
      }

    } else if (resulttype == UNPAIRED_TRANSLOC) {
      Filestring_set_split_output(fp,OUTPUT_UT);
      translocationp = true;

    } else if (resulttype == UNPAIRED_MULT) {
      if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
	Filestring_set_split_output(fp,OUTPUT_UX);
	excessivep = true;
      } else {
	Filestring_set_split_output(fp,OUTPUT_UM);
      }

    } else {
      fprintf(stderr,"Resulttype is %s\n",Resulttype_string(resulttype));
      abort();
    }

    if (output_type != M8_OUTPUT) {
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,UNPAIRED_TEXT);
      if (translocationp == true) {
	FPRINTF(fp," (transloc)");
      }

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
    }

    if (excessivep == true) {
      /* No output */
      if (failedinput_root != NULL) {
	Shortread_print_query_singleend(fp_failedinput,queryseq,headerseq1);
      }
					      
    } else {
      if (firstp == true) {
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	  this = stage3array[pathnum-1];
	  Stage3end_print(fp,/*stage3pair*/NULL,this,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/1",
			  invertp,/*hit5*/(T) NULL,/*hit3*/(T) NULL,
			  /*insertlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,this->mapq_score,
			  /*first_read_p*/true);
	}
      } else {
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	  this = stage3array[pathnum-1];
	  Stage3end_print(fp,/*stage3pair*/NULL,this,
			  chromosome_iit,queryseq,headerseq1,/*acc_suffix*/"/2",
			  invertp,/*hit5*/(T) NULL,/*hit3*/(T) NULL,
			  /*insertlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,this->mapq_score,
			  /*first_read_p*/false);
	}
      }
    }

    if (output_type != M8_OUTPUT) {
      FPRINTF(fp,"\n");
    }
  }

  return;
}


#ifdef RESOLVE_INSIDE_GENERAL
/* Used for resolve_inside_general_splice_plus and resolve_inside_general_splice_minus */
static List_T
Stage3end_convert_to_pairs (List_T pairs, T hit, char *queryuc_ptr,
			    Chrpos_T chrlength, Pairpool_T pairpool) {
  List_T p, q;
  /* Chrpos_T genomicpos1, genomicpos2; */
  Substring_T substring, prev_substring;
  Junction_T junction;
  Junctiontype_T type;
  char *deletion_string;

  if (hit->hittype == TRANSLOC_SPLICE) {
    /* Cannot handle translocations within a single GMAP alignment */
    abort();
    return NULL;
    
  } else {
    p = hit->substrings_1toN;
    prev_substring = (Substring_T) List_head(p);
    debug9(printf("Converting substring\n"));
    /* Normally done during Stage3pair_eval_and_sort */
    Substring_display_prep(prev_substring,queryuc_ptr,hit->querylength,/*extraleft*/0,/*extraright*/0,
			   genomecomp);
    pairs = Substring_convert_to_pairs(pairs,prev_substring,queryuc_ptr,chrlength,pairpool);

    for (q = hit->junctions_1toN, p = List_next(p); p != NULL; q = List_next(q), p = List_next(p)) {
      junction = (Junction_T) List_head(q);
      substring = (Substring_T) List_head(p);
    
      if ((type = Junction_type(junction)) == INS_JUNCTION) {
	debug9(printf("Converting insertion\n"));
	pairs = Substring_add_insertion(pairs,prev_substring,substring,
					/*insertionlength*/Junction_nindels(junction),
					queryuc_ptr,pairpool);
      } else if (type == DEL_JUNCTION) {
	debug9(printf("Converting deletion\n"));
	deletion_string = Junction_deletion_string(junction,genomecomp,hit->plusp);
	pairs = Substring_add_deletion(pairs,prev_substring,substring,
				       deletion_string,/*deletionlength*/Junction_nindels(junction),
				       pairpool);
	FREE(deletion_string);
      } else if (type == SPLICE_JUNCTION) {
	/* Causes problems with bad comps.  Stage3_compute_one will insert gaps anyway */
	debug9(printf("(Not converting splice)\n"));
	/* pairs = Substring_add_intron(pairs,prev_substring,substring,pairpool); */
	
      } else {
	abort();
      }
    
      debug9(printf("Converting substring\n"));
      /* Normally done during Stage3pair_eval_and_sort */
      Substring_display_prep(substring,queryuc_ptr,hit->querylength,/*extraleft*/0,/*extraright*/0,
			     genomecomp);
      pairs = Substring_convert_to_pairs(pairs,substring,queryuc_ptr,chrlength,pairpool);
      prev_substring = substring;
    }

    debug9(Simplepair_dump_list(pairs,true));
    return pairs;
  }
}
#endif


/* Used only for --merge-overlap features, so obey hardclip and not querystart/queryend */
/* If use querylength_adj, ss.bug.4 fails.  If use querylength, ss.bug.3 fails */
static List_T
Stage3end_convert_to_pairs_out (List_T pairs, T hit, Shortread_T queryseq,
				int hardclip_low, int hardclip_high, int queryseq_offset) {
  List_T p, q;
  /* Chrpos_T genomicpos1, genomicpos2; */
  Substring_T substring, prev_substring;
  Junction_T junction;
  Junctiontype_T type;
  char *deletion_string;

  if (hit->hittype == TRANSLOC_SPLICE) {
    /* Cannot handle translocations within a single GMAP alignment */
    abort();
    return NULL;
    
  } else {
    p = hit->substrings_1toN;
    prev_substring = (Substring_T) List_head(p);
    pairs = Substring_convert_to_pairs_out(pairs,prev_substring,hit->querylength,
					   queryseq,hardclip_low,hardclip_high,queryseq_offset);

    for (q = hit->junctions_1toN, p = List_next(p); p != NULL; q = List_next(q), p = List_next(p)) {
      junction = (Junction_T) List_head(q);
      substring = (Substring_T) List_head(p);
    
      if ((type = Junction_type(junction)) == INS_JUNCTION) {
	pairs = Substring_add_insertion_out(pairs,prev_substring,substring,hit->querylength,
					    /*insertionlength*/Junction_nindels(junction),queryseq,
					    hardclip_low,hardclip_high,queryseq_offset);
      } else if (type == DEL_JUNCTION) {
	deletion_string = Junction_deletion_string(junction,genomecomp,hit->plusp);
	pairs = Substring_add_deletion_out(pairs,prev_substring,substring,hit->querylength,
					   deletion_string,/*deletionlength*/Junction_nindels(junction),
					   hardclip_low,hardclip_high,queryseq_offset);
      } else if (type == SPLICE_JUNCTION) {
	pairs = Substring_add_intron_out(pairs,prev_substring,substring,hit->querylength,
					 hardclip_low,hardclip_high,queryseq_offset);
	
      } else {
	abort();
      }
    
      pairs = Substring_convert_to_pairs_out(pairs,substring,hit->querylength,
					     queryseq,hardclip_low,hardclip_high,queryseq_offset);
      prev_substring = substring;
    }

    debug15(Simplepair_dump_list(pairs,true));
    return pairs;
  }
}


/* Don't want querylength_adj */
struct Simplepair_T *
Stage3pair_merge (int *npairs, int *querylength_merged, char **queryseq_merged, char **quality_merged,
		  Stage3pair_T this, Shortread_T queryseq5, Shortread_T queryseq3,
		  int querylength5, int querylength3, int clipdir,
		  int hardclip5_low, int hardclip5_high, int hardclip3_low, int hardclip3_high) {
  struct Simplepair_T *pairarray, *newpair;
  Simplepair_T oldpair;
  List_T pairs, pairs5, pairs3, p;
  T hit5, hit3;
  int querylengthA, querylengthB;
  char *queryseq_ptr_5, *queryseq_ptr_3, *quality_ptr_5, *quality_ptr_3;
#ifdef CHECK_ASSERTIONS
  Chrpos_T genomicpos1, genomicpos2;
#endif

  hit5 = this->hit5;
  hit3 = this->hit3;
  queryseq_ptr_5 = Shortread_fullpointer_uc(queryseq5);
  queryseq_ptr_3 = Shortread_fullpointer_uc(queryseq3);
  quality_ptr_5 = Shortread_quality_string(queryseq5);
  quality_ptr_3 = Shortread_quality_string(queryseq3);

  if (hit5->plusp == true) {
    if (clipdir > 0) {
      pairs5 = Stage3end_convert_to_pairs_out(NULL,hit5,queryseq5,hardclip5_low,hardclip5_high,/*queryseq_offset*/0);
      pairs5 = Simplepair_strip_gaps_at_head(pairs5);

      pairs3 = Stage3end_convert_to_pairs_out(NULL,hit3,queryseq3,hardclip3_low,hardclip3_high,
					      /*queryseq_offset*/querylength5-hardclip5_low-hardclip5_high-hardclip3_low-hardclip3_high);
      pairs3 = Simplepair_strip_gaps_at_tail(pairs3);

#ifdef CHECK_ASSERTIONS
      genomicpos1 = Simplepair_head_genomepos(pairs5);
      genomicpos2 = Simplepair_last_genomepos(pairs3);
      if (genomicpos2 != genomicpos1 + 1U) {
	printf("Accession %s, plus\n",Shortread_accession(queryseq5));
	printf("Expected genomicpos2 %u == genomicpos1 %u + 1\n",genomicpos2,genomicpos1);
	Simplepair_dump_list(pairs5,true);
	Simplepair_dump_list(pairs3,true);
	abort();
      }
#endif      

      pairs = List_append(pairs3,pairs5);

      querylengthA = querylength5 - hardclip5_low - hardclip5_high;
      querylengthB = querylength3 - hardclip3_low - hardclip3_high;
      *querylength_merged = querylengthA + querylengthB;

      *queryseq_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
      strncpy(*queryseq_merged,queryseq_ptr_5,querylengthA);
      strncpy(&((*queryseq_merged)[querylengthA]),&(queryseq_ptr_3[querylength3 - querylengthB]),querylengthB);
      (*queryseq_merged)[querylengthA+querylengthB] = '\0';

      if (quality_ptr_5 == NULL || quality_ptr_3 == NULL) {
	*quality_merged = (char *) NULL;
      } else {
	*quality_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
	strncpy(*quality_merged,quality_ptr_5,querylengthA);
	strncpy(&((*quality_merged)[querylengthA]),&(quality_ptr_3[querylength3 - querylengthB]),querylengthB);
	(*quality_merged)[querylengthA+querylengthB] = '\0';
      }

    } else if (clipdir < 0) {
      pairs3 = Stage3end_convert_to_pairs_out(NULL,hit3,queryseq3,hardclip3_low,hardclip3_high,/*queryseq_offset*/0);
      pairs3 = Simplepair_strip_gaps_at_head(pairs3);

      pairs5 = Stage3end_convert_to_pairs_out(NULL,hit5,queryseq5,hardclip5_low,hardclip5_high,
					      /*queryseq_offset*/querylength3-hardclip3_low-hardclip3_high-hardclip5_low-hardclip5_high);
      pairs5 = Simplepair_strip_gaps_at_tail(pairs5);

#ifdef CHECK_ASSERTIONS
      genomicpos1 = Simplepair_head_genomepos(pairs3);
      genomicpos2 = Simplepair_last_genomepos(pairs5);
      if (genomicpos2 != genomicpos1 + 1U) {
	printf("Accession %s, plus, clipdir %d\n",Shortread_accession(queryseq5),clipdir);
	printf("Expected genomicpos2 %u == genomicpos1 %u + 1\n",genomicpos2,genomicpos1);
	printf("Begin of pairs3\n");
	Simplepair_dump_list(pairs3,true);
	printf("Begin of pairs5\n");
	Simplepair_dump_list(pairs5,true);
	abort();
      }
#endif      

      pairs = List_append(pairs5,pairs3);

      querylengthA = querylength3 - hardclip3_low - hardclip3_high;
      querylengthB = querylength5 - hardclip5_low - hardclip5_high;
      *querylength_merged = querylengthA + querylengthB;

      *queryseq_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
      strncpy(*queryseq_merged,queryseq_ptr_3,querylengthA);
      strncpy(&((*queryseq_merged)[querylengthA]),&(queryseq_ptr_5[querylength5 - querylengthB]),querylengthB);
      (*queryseq_merged)[querylengthA+querylengthB] = '\0';

      if (quality_ptr_5 == NULL || quality_ptr_3 == NULL) {
	*quality_merged = (char *) NULL;
      } else {
	*quality_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
	strncpy(*quality_merged,quality_ptr_3,querylengthA);
	strncpy(&((*quality_merged)[querylengthA]),&(quality_ptr_5[querylength5 - querylengthB]),querylengthB);
	(*quality_merged)[querylengthA+querylengthB] = '\0';
      }

    } else {
      abort();
    }

  } else {
    if (clipdir > 0) {
      pairs3 = Stage3end_convert_to_pairs_out(NULL,hit3,queryseq3,hardclip3_low,hardclip3_high,/*queryseq_offset*/0);
      pairs3 = Simplepair_strip_gaps_at_head(pairs3);

      pairs5 = Stage3end_convert_to_pairs_out(NULL,hit5,queryseq5,hardclip5_low,hardclip5_high,
					      /*queryseq_offset*/querylength3-hardclip3_low-hardclip3_high-hardclip5_low-hardclip5_high);
      pairs5 = Simplepair_strip_gaps_at_tail(pairs5);

#ifdef CHECK_ASSERTIONS
      genomicpos1 = Simplepair_head_genomepos(pairs3);
      genomicpos2 = Simplepair_last_genomepos(pairs5);
      if (genomicpos2 != genomicpos1 - 1U) {
	printf("Accession %s, minus\n",Shortread_accession(queryseq5));
	printf("Expected genomicpos2 %u == genomicpos1 %u - 1\n",genomicpos2,genomicpos1);
	Simplepair_dump_list(pairs3,true);
	Simplepair_dump_list(pairs5,true);
	abort();
      }
#endif      

      pairs = List_append(pairs5,pairs3);

      querylengthA = querylength3 - hardclip3_low - hardclip3_high;
      querylengthB = querylength5 - hardclip5_low - hardclip5_high;
      *querylength_merged = querylengthA + querylengthB;

      *queryseq_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
      strncpy(*queryseq_merged,queryseq_ptr_3,querylengthA);
      strncpy(&((*queryseq_merged)[querylengthA]),&(queryseq_ptr_5[querylength5 - querylengthB]),querylengthB);
      (*queryseq_merged)[querylengthA+querylengthB] = '\0';

      if (quality_ptr_5 == NULL || quality_ptr_3 == NULL) {
	*quality_merged = (char *) NULL;
      } else {
	*quality_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
	strncpy(*quality_merged,quality_ptr_3,querylengthA);
	strncpy(&((*quality_merged)[querylengthA]),&(quality_ptr_5[querylength5 - querylengthB]),querylengthB);
	(*quality_merged)[querylengthA+querylengthB] = '\0';
      }

    } else if (clipdir < 0) {
      pairs5 = Stage3end_convert_to_pairs_out(NULL,hit5,queryseq5,hardclip5_low,hardclip5_high,/*queryseq_offset*/0);
      pairs5 = Simplepair_strip_gaps_at_head(pairs5);

      pairs3 = Stage3end_convert_to_pairs_out(NULL,hit3,queryseq3,hardclip3_low,hardclip3_high,
					      /*queryseq_offset*/querylength5-hardclip5_low-hardclip5_high-hardclip3_low-hardclip3_high);
      pairs3 = Simplepair_strip_gaps_at_tail(pairs3);

#ifdef CHECK_ASSERTIONS
      genomicpos1 = Simplepair_head_genomepos(pairs5);
      genomicpos2 = Simplepair_last_genomepos(pairs3);
      if (genomicpos2 != genomicpos1 - 1U) {
	printf("Accession %s, minus\n",Shortread_accession(queryseq5));
	printf("Expected genomicpos2 %u == genomicpos1 %u - 1\n",genomicpos2,genomicpos1);
	Simplepair_dump_list(pairs5,true);
	Simplepair_dump_list(pairs3,true);
	abort();
      }
#endif      

      pairs = List_append(pairs3,pairs5);

      querylengthA = querylength5 - hardclip5_low - hardclip5_high;
      querylengthB = querylength3 - hardclip3_low - hardclip3_high;
      *querylength_merged = querylengthA + querylengthB;

      *queryseq_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
      strncpy(*queryseq_merged,queryseq_ptr_5,querylengthA);
      strncpy(&((*queryseq_merged)[querylengthA]),&(queryseq_ptr_3[querylength3 - querylengthB]),querylengthB);
      (*queryseq_merged)[querylengthA+querylengthB] = '\0';

      if (quality_ptr_5 == NULL || quality_ptr_3 == NULL) {
	*quality_merged = (char *) NULL;
      } else {
	*quality_merged = (char *) MALLOC_OUT((querylengthA+querylengthB+1) * sizeof(char));
	strncpy(*quality_merged,quality_ptr_5,querylengthA);
	strncpy(&((*quality_merged)[querylengthA]),&(quality_ptr_3[querylength3 - querylengthB]),querylengthB);
	(*quality_merged)[querylengthA+querylengthB] = '\0';
      }

    } else {
      abort();
    }
  }

  pairs = List_reverse(pairs);
  /* Simplepair_dump_list(pairs,true); */

  *npairs = List_length(pairs);
  newpair = pairarray = (struct Simplepair_T *) MALLOC_OUT((*npairs)*sizeof(struct Simplepair_T));
  for (p = pairs; p != NULL; p = p->rest) {
    oldpair = (Simplepair_T) p->first;
    memcpy(newpair++,oldpair,sizeof(struct Simplepair_T));
    Simplepair_free_out(&oldpair);
  }
  List_free_out(&pairs);

  return pairarray;
}



#ifdef RESOLVE_INSIDE_GENERAL
/* Previously had private5p and private3p as parameters, but now that
   we are always copying, can assume they are private */
static bool
resolve_inside_general_splice_plus (T *oldhit5, T *oldhit3, int *alts_resolve_5, int *alts_resolve_3,
				    Compress_T query5_compress_fwd, Compress_T query3_compress_fwd, 
				    char *queryuc_ptr_5, char *queryuc_ptr_3, int querylength5, int querylength3,
				    int genestrand, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
				    Listpool_T listpool) {
  bool changep = false;
  T newhit5 = NULL, newhit3 = NULL, hit5 = *oldhit5, hit3 = *oldhit3;

#ifdef DEBUG9
  List_T p;
#endif
  List_T stage2pairs, all_stage2_starts, all_stage2_ends;
  int queryend, endlength;
  Chrpos_T chrstart, chrend;
  struct Simplepair_T *pairarray1, *pairarray2;
  List_T pairs1, pairs2;
  int found_score_5 = querylength5, found_score_3 = querylength3;

  int cdna_direction, sensedir, sense_try;
  int npairs1, goodness1, matches1, nmatches_to_trims_1,
    max_match_length_1, ambig_end_length_5_1, ambig_end_length_3_1,
    unknowns1, mismatches1, qopens1, qindels1, topens1, tindels1,
    ncanonical1, nsemicanonical1, nnoncanonical1;
  int npairs2, goodness2, matches2, nmatches_to_trims_2,
    max_match_length_2, ambig_end_length_5_2, ambig_end_length_3_2,
    unknowns2, mismatches2, qopens2, qindels2, topens2, tindels2,
    ncanonical2, nsemicanonical2, nnoncanonical2;
  double ambig_prob_5_1, ambig_prob_3_1, avg_splice_score_1;
  double ambig_prob_5_2, ambig_prob_3_2, avg_splice_score_2;
  Splicetype_T ambig_splicetype_5_1, ambig_splicetype_3_1;
  Splicetype_T ambig_splicetype_5_2, ambig_splicetype_3_2;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;

  Univcoord_T start;


  debug9(printf("Entered resolve_inside_general_splice_plus with hittypes %s (%d..%d) and %s (%d..%d)\n",
		hittype_string(hit5->hittype),Stage3end_querystart(hit5),Stage3end_queryend(hit5),
		hittype_string(hit3->hittype),Stage3end_querystart(hit3),Stage3end_queryend(hit3)));

  if (hit5->genomicstart > hit3->genomicstart || hit5->genomicend > hit3->genomicend) {
    /* Scramble, which could occur with circular chromosomes */
    debug9(printf("Scramble, possibly from a circular chromosome.  Not solving at this time\n"));

#if 0
  } else if (hit5->querylength - 1 - Stage3end_queryend(hit5) > 10 && Stage3end_querystart(hit3) > 10) {
    /* Now trying to solve both insides at this time.  If inconsistent, then skip. */
    /* Both insides need to be resolved.  Not solving at this time */
    debug9(printf("Dual to be resolved on inside.  Not solving at this time\n"));
#endif

  } else {
    if (hit5->hittype == SPLICE) {
      /* Cannot convert a SPLICE to pairs, and a SPLICE is currently
	 defined to have only two substrings */
      debug9(printf("5' end has hittype of SPLICE.  Not solving at this time\n"));
      
    } else if (hit5->chrnum != 0 && (endlength = hit5->querylength - 1 - Stage3end_queryend(hit5)) > 10 &&
	       hit3->genomicstart > hit5->chroffset) {
      chrend = hit3->genomicstart - hit5->chroffset; /* Use hit5->chroffset in case hit3 is a transloc */
      chrstart = subtract_bounded(chrend,(Chrpos_T) expected_pairlength + pairlength_deviation + endlength,0);
      if (chrstart < hit5->genomicend - hit5->chroffset) {
	debug9(printf("Revising chrstart\n"));
	chrstart = hit5->genomicend - hit5->chroffset;
      }
      queryend = Stage3end_queryend(hit5) + 1;
      debug9(printf("Resolve plus 5': For ends, chrstart %u, chrend %u\n",chrstart,chrend));
      if (chrstart < chrend &&
	  (all_stage2_ends = Stage2_compute_ends(
#ifdef PMAP
						 &(queryaaseq_ptr[queryend]),&(queryaaseq_ptr[queryend]),
						 /*querylength*/endlength,/*query_offset*/0*3,
#else
						 &(queryuc_ptr_5[queryend]),&(queryuc_ptr_5[queryend]),
						 /*querylength*/endlength,/*query_offset*/queryend,
#endif
						 chrstart,chrend,hit5->chroffset,hit5->chrhigh,/*plusp*/true,genestrand,
						 
						 oligoindices_minor,pairpool,diagpool,cellpool,
						 /*localp should be false*/true,/*skip_repetitive_p*/false,
						 /*favor_right_p*/false,/*max_nalignments*/2,/*debug_graphic_p*/false)) != NULL) {
	
	debug9(printf("Got %d ends\n",List_length(all_stage2_ends)));
	debug9(printf("5' end to be resolved on inside\n"));
#ifdef DEBUG9
	for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
	  Simplepair_dump_list(List_head(p),true);
	}
#endif
	stage2pairs = Stage3end_convert_to_pairs(/*pairs*/NULL,hit5,queryuc_ptr_5,
						 /*chrlength*/hit5->chrhigh - hit5->chroffset,pairpool);
	debug9(Simplepair_dump_list(stage2pairs,true));
	
	knownsplice_limit_high = ((Simplepair_T) stage2pairs->first)->genomepos + hit5->chroffset;
	stage2pairs = List_reverse(stage2pairs);
	knownsplice_limit_low = ((Simplepair_T) stage2pairs->first)->genomepos + hit5->chroffset;
	
	if ((sensedir = Stage3end_sensedir(hit3)) == SENSE_FORWARD) {
	  sense_try = +1;
	} else if (sensedir == SENSE_ANTI) {
	  sense_try = -1;
	} else {
	  sense_try = 0;
	}
	
	if ((pairarray1 = Stage3_compute_one(&cdna_direction,&sensedir,&pairs1,&npairs1,&goodness1,
					     &matches1,&nmatches_to_trims_1,&max_match_length_1,
					     &ambig_end_length_5_1,&ambig_end_length_3_1,
					     &ambig_splicetype_5_1,&ambig_splicetype_3_1,
					     &ambig_prob_5_1,&ambig_prob_3_1,
					     &unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
					     &ncanonical1,&nsemicanonical1,&nnoncanonical1,&avg_splice_score_1,
					     
					     &pairarray2,&pairs2,&npairs2,&goodness2,
					     &matches2,&nmatches_to_trims_2,&max_match_length_2,
					     &ambig_end_length_5_2,&ambig_end_length_3_2,
					     &ambig_splicetype_5_2,&ambig_splicetype_3_2,
					     &ambig_prob_5_2,&ambig_prob_3_2,
					     &unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
					     &ncanonical2,&nsemicanonical2,&nnoncanonical2,&avg_splice_score_2,
					     
					     stage2pairs,/*all_stage2_starts*/NULL,all_stage2_ends,
#ifdef END_KNOWNSPLICING_SHORTCUT
					     cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
					     watsonp ? query_compress_fwd : query_compress_rev,
#endif
					     /*queryseq_ptr*/queryuc_ptr_5,/*queryuc_ptr*/queryuc_ptr_5,
					     querylength5,/*skiplength*/0,/*query_subseq_offset*/0,
					     hit5->chrnum,hit5->chroffset,hit5->chrhigh,
					     knownsplice_limit_low,knownsplice_limit_high,/*plusp*/true,genestrand,
					     /*jump_late_p*/false,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
					     sense_try,/*sense_filter*/0)) == NULL) {
	  
	} else if (pairarray2 != NULL) {
	  if (avg_splice_score_1 > avg_splice_score_2) {
	    FREE_OUT(pairarray2);
	    start = subtract_bounded(hit5->chroffset + Pair_genomepos(&(pairarray1[0])),
				     /*minusterm*/Pair_querypos(&(pairarray1[0])),hit5->chroffset);
#if 0
	    end = add_bounded(hit5->chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			      /*plusterm*/querylength5 - 1 - Pair_querypos(&(pairarray1[npairs1-1])),hit5->chrhigh);
#endif
	    newhit5 = Stage3end_new_gmap(&found_score_5,nmatches_to_trims_1,max_match_length_1,
					 ambig_end_length_5_1,ambig_end_length_3_1,
					 avg_splice_score_1,pairarray1,npairs1,
					 /*left*/start,/*plusp*/true,genestrand,
					 querylength5,query5_compress_fwd,
					 hit5->chrnum,hit5->chroffset,hit5->chrhigh,hit5->chrlength,
					 /*cdna_direction*/+1,/*sensedir*/SENSE_FORWARD,listpool,
					 hit5->method,hit5->level);
	    
	  } else {
	    FREE_OUT(pairarray1);
	    start = subtract_bounded(hit5->chroffset + Pair_genomepos(&(pairarray2[0])),
				     /*minusterm*/Pair_querypos(&(pairarray2[0])),hit5->chroffset);
#if 0
	    end = add_bounded(hit5->chroffset + Pair_genomepos(&(pairarray2[npairs2-1])),
			      /*plusterm*/querylength5 - 1 - Pair_querypos(&(pairarray2[npairs2-1])),hit5->chrhigh);
#endif
	    newhit5 = Stage3end_new_gmap(&found_score_5,nmatches_to_trims_2,max_match_length_2,
					 ambig_end_length_5_2,ambig_end_length_3_2,
					 avg_splice_score_2,pairarray2,npairs2,
					 /*left*/start,/*plusp*/true,genestrand,
					 querylength5,query5_compress_fwd,
					 hit5->chrnum,hit5->chroffset,hit5->chrhigh,hit5->chrlength,
					 /*cdna_direction*/-1,/*sensedir*/SENSE_ANTI,listpool,
					 hit5->method,hit5->level);
	  }

	} else {
	  start = subtract_bounded(hit5->chroffset + Pair_genomepos(&(pairarray1[0])),
				   /*minusterm*/Pair_querypos(&(pairarray1[0])),hit5->chroffset);
#if 0
	  end = add_bounded(hit5->chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			    /*plusterm*/querylength5 - 1 - Pair_querypos(&(pairarray1[npairs1-1])),hit5->chrhigh);
#endif	  
	  newhit5 = Stage3end_new_gmap(&found_score_5,nmatches_to_trims_1,max_match_length_1,
				       ambig_end_length_5_1,ambig_end_length_3_1,
				       avg_splice_score_1,pairarray1,npairs1,
				       /*left*/start,/*plusp*/true,genestrand,
				       querylength5,query5_compress_fwd,
				       hit5->chrnum,hit5->chroffset,hit5->chrhigh,hit5->chrlength,
				       cdna_direction,sensedir,listpool,hit5->method,hit5->level);
	}
	
	List_free(&all_stage2_ends);
      }
    }
      

    if (hit3->hittype == SPLICE) {
      /* Cannot convert a SPLICE to pairs, and a SPLICE is currently
	 defined to have only two substrings */
      debug9(printf("3' end has hittype of SPLICE.  Not solving at this time\n"));
      
    } else if (hit3->chrnum != 0 && (endlength = Stage3end_querystart(hit3)) > 10 &&
	       hit5->genomicend > hit3->chroffset) {
      chrstart = hit5->genomicend - hit3->chroffset; /* Use hit3->chroffset in case hit5 is a transloc */
      chrend = add_bounded(chrstart,(Chrpos_T) expected_pairlength + pairlength_deviation + endlength,hit5->chrhigh);
      if (chrend > hit3->genomicstart - hit3->chroffset) {
	debug9(printf("Revising chrend\n"));
	chrend = hit3->genomicstart - hit3->chroffset;
      }
      debug9(printf("Resolve plus 3': For starts, chrstart %u, chrend %u\n",chrstart,chrend));
      if (chrstart < chrend && 
	  (all_stage2_starts = Stage2_compute_starts(
#ifdef PMAP
						     &(queryaaseq_ptr[0]),&(queryaaseq_ptr[0]),
						     /*querylength*/endlength,/*query_offset*/0*3,
#else
						     /*queryseq_ptr*/&(queryuc_ptr_3[0]),
						     /*queryuc_ptr*/&(queryuc_ptr_3[0]),
						     /*querylength*/endlength,/*query_offset*/0,
#endif
						     chrstart,chrend,hit3->chroffset,hit3->chrhigh,/*plusp*/true,genestrand,
						     
						     oligoindices_minor,pairpool,diagpool,cellpool,
						     /*localp should be false*/true,/*skip_repetitive_p*/false,
						     /*favor_right_p*/true,/*max_nalignments*/2,/*debug_graphic_p*/false)) != NULL) {
	
	debug9(printf("Got %d starts\n",List_length(all_stage2_starts)));
	debug9(printf("3' start to be resolved on inside\n"));
#ifdef DEBUG9
	for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
	  Simplepair_dump_list(List_head(p),true);
	}
#endif
	stage2pairs = Stage3end_convert_to_pairs(/*pairs*/NULL,hit3,queryuc_ptr_3,
						 /*chrlength*/hit3->chrhigh - hit3->chroffset,pairpool);
	debug9(Simplepair_dump_list(stage2pairs,true));
	
	knownsplice_limit_high = ((Simplepair_T) stage2pairs->first)->genomepos + hit3->chroffset;
	stage2pairs = List_reverse(stage2pairs);
	knownsplice_limit_low = ((Simplepair_T) stage2pairs->first)->genomepos + hit3->chroffset;
	
	if ((sensedir = Stage3end_sensedir(hit5)) == SENSE_FORWARD) {
	  sense_try = +1;
	} else if (sensedir == SENSE_ANTI) {
	  sense_try = -1;
	} else {
	  sense_try = 0;
	}
	
	if ((pairarray1 = Stage3_compute_one(&cdna_direction,&sensedir,&pairs1,&npairs1,&goodness1,
					     &matches1,&nmatches_to_trims_1,&max_match_length_1,
					     &ambig_end_length_5_1,&ambig_end_length_3_1,
					     &ambig_splicetype_5_1,&ambig_splicetype_3_1,
					     &ambig_prob_5_1,&ambig_prob_3_1,
					     &unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
					     &ncanonical1,&nsemicanonical1,&nnoncanonical1,&avg_splice_score_1,
					     
					     &pairarray2,&pairs2,&npairs2,&goodness2,
					     &matches2,&nmatches_to_trims_2,&max_match_length_2,
					     &ambig_end_length_5_2,&ambig_end_length_3_2,
					     &ambig_splicetype_5_2,&ambig_splicetype_3_2,
					     &ambig_prob_5_2,&ambig_prob_3_2,
					     &unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
					     &ncanonical2,&nsemicanonical2,&nnoncanonical2,&avg_splice_score_2,
					     
					     stage2pairs,all_stage2_starts,/*all_stage2_ends*/NULL,
#ifdef END_KNOWNSPLICING_SHORTCUT
					     cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
					     watsonp ? query_compress_fwd : query_compress_rev,
#endif
					     /*queryseq_ptr*/queryuc_ptr_3,/*queryuc_ptr*/queryuc_ptr_3,
					     querylength3,/*skiplength*/0,/*query_subseq_offset*/0,
					     hit3->chrnum,hit3->chroffset,hit3->chrhigh,
					     knownsplice_limit_low,knownsplice_limit_high,/*plusp*/true,genestrand,
					     /*jump_late_p*/false,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
					     sense_try,/*sense_filter*/0)) == NULL) {

	} else if (pairarray2 != NULL) {
	  if (avg_splice_score_1 > avg_splice_score_2) {
	    FREE_OUT(pairarray2);
	    start = subtract_bounded(hit3->chroffset + Pair_genomepos(&(pairarray1[0])),
				     /*minusterm*/Pair_querypos(&(pairarray1[0])),hit3->chroffset);
#if 0
	    end = add_bounded(hit3->chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			      /*plusterm*/querylength3 - 1 - Pair_querypos(&(pairarray1[npairs1-1])),hit3->chrhigh);
#endif
	    newhit3 = Stage3end_new_gmap(&found_score_3,nmatches_to_trims_1,max_match_length_1,
					 ambig_end_length_5_1,ambig_end_length_3_1,
					 avg_splice_score_1,pairarray1,npairs1,
					 /*left*/start,/*plusp*/true,genestrand,
					 querylength3,query3_compress_fwd,
					 hit3->chrnum,hit3->chroffset,hit3->chrhigh,hit3->chrlength,
					 /*cdna_direction*/+1,/*sensedir*/SENSE_FORWARD,listpool,
					 hit3->method,hit3->level);

	  } else {
	    FREE_OUT(pairarray1);
	    start = subtract_bounded(hit3->chroffset + Pair_genomepos(&(pairarray2[0])),
				     /*minusterm*/Pair_querypos(&(pairarray2[0])),hit3->chroffset);
#if 0
	    end = add_bounded(hit3->chroffset + Pair_genomepos(&(pairarray2[npairs2-1])),
			      /*plusterm*/querylength3 - 1 - Pair_querypos(&(pairarray2[npairs2-1])),hit3->chrhigh);
#endif
	    newhit3 = Stage3end_new_gmap(&found_score_3,nmatches_to_trims_2,max_match_length_2,
					 ambig_end_length_5_2,ambig_end_length_3_2,
					 avg_splice_score_2,pairarray2,npairs2,
					 /*left*/start,/*plusp*/true,genestrand,
					 querylength3,query3_compress_fwd,
					 hit3->chrnum,hit3->chroffset,hit3->chrhigh,hit3->chrlength,
					 /*cdna_direction*/-1,/*sensedir*/SENSE_ANTI,listpool,
					 hit3->method,hit3->level);
	  }

	} else {
	  start = subtract_bounded(hit3->chroffset + Pair_genomepos(&(pairarray1[0])),
				   /*minusterm*/Pair_querypos(&(pairarray1[0])),hit3->chroffset);
#if 0
	  end = add_bounded(hit3->chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			    /*plusterm*/querylength3 - 1 - Pair_querypos(&(pairarray1[npairs1-1])),hit3->chrhigh);
#endif
	  newhit3 = Stage3end_new_gmap(&found_score_3,nmatches_to_trims_1,max_match_length_1,
				       ambig_end_length_5_1,ambig_end_length_3_1,
				       avg_splice_score_1,pairarray1,npairs1,
				       /*left*/start,/*plusp*/true,genestrand,
				       querylength3,query3_compress_fwd,
				       hit3->chrnum,hit3->chroffset,hit3->chrhigh,hit3->chrlength,
				       cdna_direction,sensedir,listpool,hit3->method,hit3->level);
	}

	List_free(&all_stage2_starts);
      }
    }


    if (newhit5 != NULL && newhit3 != NULL) {
      /* Check for consistency */
      if (hit5->genomicend >= hit3->genomicstart) {
	debug9(printf("INCONSISTENT\n"));
	Stage3end_free(&newhit5);
	Stage3end_free(&newhit3);

      } else {
	Stage3end_free(&(*oldhit5));
	Stage3end_free(&(*oldhit3));
	debug9(printf("Both 5' and 3' resolved on inside\n"));
	*oldhit5 = newhit5;
	*oldhit3 = newhit3;
	*alts_resolve_5 = *alts_resolve_3 = -1;
	changep = true;
      }

    } else if (newhit5 != NULL) {
      Stage3end_free(&(*oldhit5));
      debug9(printf("5' resolved on inside\n"));
      *oldhit5 = newhit5;
      *alts_resolve_5 = -1;
      changep = true;

    } else if (newhit3 != NULL) {
      Stage3end_free(&(*oldhit3));
      debug9(printf("3' resolved on inside\n"));
      *oldhit3 = newhit3;
      *alts_resolve_3 = -1;
      changep = true;
    }
  }

  return changep;
}
#endif


#ifdef RESOLVE_INSIDE_GENERAL
/* Previously had private5p and private3p as parameters, but now that
   we are always copying, can assume they are private */
static bool
resolve_inside_general_splice_minus (T *oldhit5, T *oldhit3, int *alts_resolve_5, int *alts_resolve_3,
				     Compress_T query5_compress_rev, Compress_T query3_compress_rev,
				     char *queryuc_ptr_5, char *queryuc_ptr_3, int querylength5, int querylength3,
				     int genestrand, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				     Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
				     Listpool_T listpool) {
  bool changep = false;
  T newhit5 = NULL, newhit3 = NULL, hit5 = *oldhit5, hit3 = *oldhit3;

#ifdef DEBUG9
  List_T p;
#endif
  List_T stage2pairs, all_stage2_starts, all_stage2_ends;
  int queryend, endlength;
  Chrpos_T chrstart, chrend;
  struct Simplepair_T *pairarray1, *pairarray2;
  List_T pairs1, pairs2;
  int found_score_5 = querylength5, found_score_3 = querylength3;

  int cdna_direction, sensedir, sense_try;
  int npairs1, goodness1, matches1, nmatches_to_trims_1,
    max_match_length_1, ambig_end_length_5_1, ambig_end_length_3_1,
    unknowns1, mismatches1, qopens1, qindels1, topens1, tindels1,
    ncanonical1, nsemicanonical1, nnoncanonical1;
  int npairs2, goodness2, matches2, nmatches_to_trims_2,
    max_match_length_2, ambig_end_length_5_2, ambig_end_length_3_2,
    unknowns2, mismatches2, qopens2, qindels2, topens2, tindels2,
    ncanonical2, nsemicanonical2, nnoncanonical2;
  double ambig_prob_5_1, ambig_prob_3_1, avg_splice_score_1;
  double ambig_prob_5_2, ambig_prob_3_2, avg_splice_score_2;
  Splicetype_T ambig_splicetype_5_1, ambig_splicetype_3_1;
  Splicetype_T ambig_splicetype_5_2, ambig_splicetype_3_2;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;

  Univcoord_T end;


  debug9(printf("Entered resolve_inside_general_splice_minus with hittypes %s (%d..%d) and %s (%d..%d)\n",
		hittype_string(hit5->hittype),Stage3end_querystart(hit5),Stage3end_queryend(hit5),
		hittype_string(hit3->hittype),Stage3end_querystart(hit3),Stage3end_queryend(hit3)));

  if (hit5->genomicstart < hit3->genomicstart || hit5->genomicend < hit3->genomicend) {
    /* Scramble, which could occur with circular chromosomes */
    debug9(printf("Scramble, possibly from a circular chromosome.  Not solving at this time\n"));

#if 0
  } else if (hit5->querylength - 1 - Stage3end_queryend(hit5) > 10 && Stage3end_querystart(hit3) > 10) {
    /* Now trying to solve both insides at this time.  If inconsistent, then skip. */
    /* Both insides need to be resolved.  Not solving at this time */
    debug9(printf("Dual to be resolved on inside.  Not solving at this time\n"));
#endif

  } else {
    if (hit5->hittype == SPLICE) {
      /* Cannot convert a SPLICE to pairs, and a SPLICE is currently
	 defined to have only two substrings */
      debug9(printf("5' end has hittype of SPLICE.  Not solving at this time\n"));

    } else if (hit5->chrnum != 0 && (endlength = hit5->querylength - 1 - Stage3end_queryend(hit5)) > 10 &&
	       hit3->genomicstart > hit5->chroffset) {
      chrstart = hit3->genomicstart - hit5->chroffset; /* Use hit5->chroffset in case hit3 is a transloc */
      chrend = add_bounded(chrstart,(Chrpos_T) expected_pairlength + pairlength_deviation + endlength,hit3->chrhigh);
      if (chrend > hit5->genomicend - hit5->chroffset) {
	debug9(printf("Revising chrend\n"));
	chrend = hit5->genomicend - hit5->chroffset;
      }
      queryend = Stage3end_queryend(hit5) + 1;
      debug9(printf("Resolve minus 5': For ends, chrstart %u, chrend %u\n",chrstart,chrend));
      if (chrstart < chrend && 
	  (all_stage2_ends = Stage2_compute_ends(
#ifdef PMAP
						 &(queryaaseq_ptr[queryend]),&(queryaaseq_ptr[queryend]),
						 /*querylength*/endlength,/*query_offset*/0*3,
#else
						 /*queryseq_ptr*/&(queryuc_ptr_5[queryend]),
						 /*queryuc_ptr*/&(queryuc_ptr_5[queryend]),
						 /*querylength*/endlength,/*query_offset*/queryend,
#endif
						 chrstart,chrend,hit5->chroffset,hit5->chrhigh,/*plusp*/false,genestrand,
					    
						 oligoindices_minor,pairpool,diagpool,cellpool,
						 /*localp should be false*/true,/*skip_repetitive_p*/false,
						 /*favor_right_p*/false,/*max_nalignments*/2,/*debug_graphic_p*/false)) != NULL) {

	debug9(printf("Got %d ends\n",List_length(all_stage2_ends)));
	debug9(printf("5' end to be resolved on inside\n"));
#ifdef DEBUG9
	for (p = all_stage2_ends; p != NULL; p = List_next(p)) {
	  Simplepair_dump_list(List_head(p),true);
	}
#endif
	stage2pairs = Stage3end_convert_to_pairs(/*pairs*/NULL,hit5,queryuc_ptr_5,
						 /*chrlength*/hit5->chrhigh - hit5->chroffset,pairpool);
	debug9(Simplepair_dump_list(stage2pairs,true));

	knownsplice_limit_low = ((Simplepair_T) stage2pairs->first)->genomepos + hit5->chroffset;
	stage2pairs = List_reverse(stage2pairs);
	knownsplice_limit_high = ((Simplepair_T) stage2pairs->first)->genomepos + hit5->chroffset;

	if ((sensedir = Stage3end_sensedir(hit3)) == SENSE_FORWARD) {
	  sense_try = +1;
	} else if (sensedir == SENSE_ANTI) {
	  sense_try = -1;
	} else {
	  sense_try = 0;
	}

	if ((pairarray1 = Stage3_compute_one(&cdna_direction,&sensedir,&pairs1,&npairs1,&goodness1,
					     &matches1,&nmatches_to_trims_1,&max_match_length_1,
					     &ambig_end_length_5_1,&ambig_end_length_3_1,
					     &ambig_splicetype_5_1,&ambig_splicetype_3_1,
					     &ambig_prob_5_1,&ambig_prob_3_1,
					     &unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
					     &ncanonical1,&nsemicanonical1,&nnoncanonical1,&avg_splice_score_1,
				       
					     &pairarray2,&pairs2,&npairs2,&goodness2,
					     &matches2,&nmatches_to_trims_2,&max_match_length_2,
					     &ambig_end_length_5_2,&ambig_end_length_3_2,
					     &ambig_splicetype_5_2,&ambig_splicetype_3_2,
					     &ambig_prob_5_2,&ambig_prob_3_2,
					     &unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
					     &ncanonical2,&nsemicanonical2,&nnoncanonical2,&avg_splice_score_2,

					     stage2pairs,/*all_stage2_starts*/NULL,all_stage2_ends,
#ifdef END_KNOWNSPLICING_SHORTCUT
					     cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
					     watsonp ? query_compress_fwd : query_compress_rev,
#endif
					     /*queryseq_ptr*/queryuc_ptr_5,/*queryuc_ptr*/queryuc_ptr_5,
					     querylength5,/*skiplength*/0,/*query_subseq_offset*/0,
					     hit5->chrnum,hit5->chroffset,hit5->chrhigh,
					     knownsplice_limit_low,knownsplice_limit_high,/*plusp*/false,genestrand,
					     /*jump_late_p*/true,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
					     sense_try,/*sense_filter*/0)) == NULL) {

	} else if (pairarray2 != NULL) {
	  if (avg_splice_score_1 > avg_splice_score_2) {
	    FREE_OUT(pairarray2);
#if 0
	    start = add_bounded(hit5->chroffset + Pair_genomepos(&(pairarray1[0])),
				/*plusterm*/Pair_querypos(&(pairarray1[0])),hit5->chrhigh);
#endif
	    end = subtract_bounded(hit5->chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
				   /*minusterm*/querylength5 - 1 - Pair_querypos(&(pairarray1[npairs1-1])),hit5->chroffset);
	    newhit5 = Stage3end_new_gmap(&found_score_5,nmatches_to_trims_1,max_match_length_1,
					 ambig_end_length_5_1,ambig_end_length_3_1,
					 avg_splice_score_1,pairarray1,npairs1,
					 /*left*/end,/*plusp*/false,genestrand,
					 querylength5,query5_compress_rev,
					 hit5->chrnum,hit5->chroffset,hit5->chrhigh,hit5->chrlength,
					 /*cdna_direction*/+1,/*sensedir*/SENSE_FORWARD,listpool,
					 hit5->method,hit5->level);

	  } else {
	    FREE_OUT(pairarray1);
#if 0
	    start = add_bounded(hit5->chroffset + Pair_genomepos(&(pairarray2[0])),
				/*plusterm*/Pair_querypos(&(pairarray2[0])),hit5->chrhigh);
#endif
	    end = subtract_bounded(hit5->chroffset + Pair_genomepos(&(pairarray2[npairs2-1])),
				   /*minusterm*/querylength5 - 1 - Pair_querypos(&(pairarray2[npairs2-1])),hit5->chroffset);
	    newhit5 = Stage3end_new_gmap(&found_score_5,nmatches_to_trims_2,max_match_length_2,
					 ambig_end_length_5_2,ambig_end_length_3_2,
					 avg_splice_score_2,pairarray2,npairs2,
					 /*left*/end,/*plusp*/false,genestrand,
					 querylength5,query5_compress_rev,
					 hit5->chrnum,hit5->chroffset,hit5->chrhigh,hit5->chrlength,
					 /*cdna_direction*/-1,/*sensedir*/SENSE_ANTI,listpool,
					 hit5->method,hit5->level);
	  }

	} else {
#if 0
	  start = add_bounded(hit5->chroffset + Pair_genomepos(&(pairarray1[0])),
			      /*plusterm*/Pair_querypos(&(pairarray1[0])),hit5->chrhigh);
#endif
	  end = subtract_bounded(hit5->chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
				 /*minusterm*/querylength5 - 1 - Pair_querypos(&(pairarray1[npairs1-1])),hit5->chroffset);
	  newhit5 = Stage3end_new_gmap(&found_score_5,nmatches_to_trims_1,max_match_length_1,
				       ambig_end_length_5_1,ambig_end_length_3_1,
				       avg_splice_score_1,pairarray1,npairs1,
				       /*left*/end,/*plusp*/false,genestrand,
				       querylength5,query5_compress_rev,
				       hit5->chrnum,hit5->chroffset,hit5->chrhigh,hit5->chrlength,
				       cdna_direction,sensedir,listpool,hit5->method,hit5->level);
	}

	List_free(&all_stage2_ends);
      }
    }

    if (hit3->hittype == SPLICE) {
      /* Cannot convert a SPLICE to pairs, and a SPLICE is currently
	 defined to have only two substrings */
      debug9(printf("3' end has hittype of SPLICE.  Not solving at this time\n"));

    } else if (hit3->chrnum != 0 && (endlength = Stage3end_querystart(hit3)) > 10 &&
	       hit5->genomicend > hit3->chroffset) {
      chrend = hit5->genomicend - hit3->chroffset; /* Use hit3->chroffset in case hit5 is a transloc */
      chrstart = subtract_bounded(chrend,(Chrpos_T) expected_pairlength + pairlength_deviation + endlength,0);
      if (chrstart < hit3->genomicstart - hit3->chroffset) {
	debug9(printf("Revising chrstart\n"));
	chrstart = hit3->genomicstart - hit3->chroffset;
      }
      debug9(printf("Resolve minus 3': For starts, chrstart %u, chrend %u\n",chrstart,chrend));
      if (chrstart < chrend && 
	  (all_stage2_starts = Stage2_compute_starts(
#ifdef PMAP
						     &(queryaaseq_ptr[0]),&(queryaaseq_ptr[0]),
						     /*querylength*/endlength,/*query_offset*/0*3,
#else
						     /*queryseq_ptr*/&(queryuc_ptr_3[0]),
						     /*queryuc_ptr*/&(queryuc_ptr_3[0]),
						     /*querylength*/endlength,/*query_offset*/0,
#endif
						     chrstart,chrend,hit3->chroffset,hit3->chrhigh,/*plusp*/false,genestrand,
					    
						     oligoindices_minor,pairpool,diagpool,cellpool,
						     /*localp should be false*/true,/*skip_repetitive_p*/false,
						     /*favor_right_p*/true,/*max_nalignments*/2,/*debug_graphic_p*/false)) != NULL) {

	debug9(printf("Got %d starts\n",List_length(all_stage2_starts)));
	debug9(printf("3' start to be resolved on inside\n"));
#ifdef DEBUG9
	for (p = all_stage2_starts; p != NULL; p = List_next(p)) {
	  Simplepair_dump_list(List_head(p),true);
	}
#endif
	stage2pairs = Stage3end_convert_to_pairs(/*pairs*/NULL,hit3,queryuc_ptr_3,
						 /*chrlength*/hit3->chrhigh - hit3->chroffset,pairpool);
	debug9(Simplepair_dump_list(stage2pairs,true));

	knownsplice_limit_low = ((Simplepair_T) stage2pairs->first)->genomepos + hit3->chroffset;
	stage2pairs = List_reverse(stage2pairs);
	knownsplice_limit_high = ((Simplepair_T) stage2pairs->first)->genomepos + hit3->chroffset;

	if ((sensedir = Stage3end_sensedir(hit5)) == SENSE_FORWARD) {
	  sense_try = +1;
	} else if (sensedir == SENSE_ANTI) {
	  sense_try = -1;
	} else {
	  sense_try = 0;
	}

	if ((pairarray1 = Stage3_compute_one(&cdna_direction,&sensedir,&pairs1,&npairs1,&goodness1,
					     &matches1,&nmatches_to_trims_1,&max_match_length_1,
					     &ambig_end_length_5_1,&ambig_end_length_3_1,
					     &ambig_splicetype_5_1,&ambig_splicetype_3_1,
					     &ambig_prob_5_1,&ambig_prob_3_1,
					     &unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
					     &ncanonical1,&nsemicanonical1,&nnoncanonical1,&avg_splice_score_1,

					     &pairarray2,&pairs2,&npairs2,&goodness2,
					     &matches2,&nmatches_to_trims_2,&max_match_length_2,
					     &ambig_end_length_5_2,&ambig_end_length_3_2,
					     &ambig_splicetype_5_2,&ambig_splicetype_3_2,
					     &ambig_prob_5_2,&ambig_prob_3_2,
					     &unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
					     &ncanonical2,&nsemicanonical2,&nnoncanonical2,&avg_splice_score_2,

					     stage2pairs,all_stage2_starts,/*all_stage2_ends*/NULL,
#ifdef END_KNOWNSPLICING_SHORTCUT
					     cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
					     watsonp ? query_compress_fwd : query_compress_rev,
#endif
					     /*queryseq_ptr*/queryuc_ptr_3,/*queryuc_ptr*/queryuc_ptr_3,
					     querylength3,/*skiplength*/0,/*query_subseq_offset*/0,
					     hit3->chrnum,hit3->chroffset,hit3->chrhigh,
					     knownsplice_limit_low,knownsplice_limit_high,/*plusp*/false,genestrand,
					     /*jump_late_p*/true,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
					     sense_try,/*sense_filter*/0)) == NULL) {

	} else if (pairarray2 != NULL) {
	  if (avg_splice_score_1 > avg_splice_score_2) {
	    FREE_OUT(pairarray2);
#if 0
	    start = add_bounded(hit3->chroffset + Pair_genomepos(&(pairarray1[0])),
				/*plusterm*/Pair_querypos(&(pairarray1[0])),hit3->chrhigh);
#endif
	    end = subtract_bounded(hit3->chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
				   /*minusterm*/querylength3 - 1 - Pair_querypos(&(pairarray1[npairs1-1])),hit3->chroffset);
	    newhit3 = Stage3end_new_gmap(&found_score_3,nmatches_to_trims_1,max_match_length_1,
					 ambig_end_length_5_1,ambig_end_length_3_1,
					 avg_splice_score_1,pairarray1,npairs1,
					 /*left*/end,/*plusp*/false,genestrand,
					 querylength3,query3_compress_rev,
					 hit3->chrnum,hit3->chroffset,hit3->chrhigh,hit3->chrlength,
					 /*cdna_direction*/+1,/*sensedir*/SENSE_FORWARD,listpool,
					 hit3->method,hit3->level);

	  } else {
	    FREE_OUT(pairarray1);
#if 0
	    start = add_bounded(hit3->chroffset + Pair_genomepos(&(pairarray2[0])),
				/*plusterm*/Pair_querypos(&(pairarray2[0])),hit3->chrhigh);
#endif
	    end = subtract_bounded(hit3->chroffset + Pair_genomepos(&(pairarray2[npairs2-1])),
				   /*minusterm*/querylength3 - 1 - Pair_querypos(&(pairarray2[npairs2-1])),hit3->chroffset);
	    newhit3 = Stage3end_new_gmap(&found_score_3,nmatches_to_trims_2,max_match_length_2,
					 ambig_end_length_5_2,ambig_end_length_3_2,
					 avg_splice_score_2,pairarray2,npairs2,
					 /*left*/end,/*plusp*/false,genestrand,
					 querylength3,query3_compress_rev,
					 hit3->chrnum,hit3->chroffset,hit3->chrhigh,hit3->chrlength,
					 /*cdna_direction*/-1,/*sensedir*/SENSE_ANTI,listpool,
					 hit3->method,hit3->level);
	  }

	} else {
#if 0
	  start = add_bounded(hit3->chroffset + Pair_genomepos(&(pairarray1[0])),
			      /*plusterm*/Pair_querypos(&(pairarray1[0])),hit3->chrhigh);
#endif
	  end = subtract_bounded(hit3->chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
				 /*minusterm*/querylength3 - 1 - Pair_querypos(&(pairarray1[npairs1-1])),hit3->chroffset);
	  newhit3 = Stage3end_new_gmap(&found_score_3,nmatches_to_trims_1,max_match_length_1,
				       ambig_end_length_5_1,ambig_end_length_3_1,
				       avg_splice_score_1,pairarray1,npairs1,
				       /*left*/end,/*plusp*/false,genestrand,
				       querylength3,query3_compress_rev,
				       hit3->chrnum,hit3->chroffset,hit3->chrhigh,hit3->chrlength,
				       cdna_direction,sensedir,listpool,hit3->method,hit3->level);
	}

	List_free(&all_stage2_starts);
      }
    }


    if (newhit5 != NULL && newhit3 != NULL) {
      /* Check for consistency */
      if (hit5->genomicend <= hit3->genomicstart) {
	debug9(printf("INCONSISTENT\n"));
	Stage3end_free(&newhit5);
	Stage3end_free(&newhit3);

      } else {
	Stage3end_free(&(*oldhit5));
	Stage3end_free(&(*oldhit3));
	debug9(printf("Both 5' and 3' resolved on inside\n"));
	*oldhit5 = newhit5;
	*oldhit3 = newhit3;
	*alts_resolve_5 = *alts_resolve_3 = -1;
	changep = true;
      }

    } else if (newhit5 != NULL) {
      Stage3end_free(&(*oldhit5));
      debug9(printf("5' resolved on inside\n"));
      *oldhit5 = newhit5;
      *alts_resolve_5 = -1;
      changep = true;

    } else if (newhit3 != NULL) {
      Stage3end_free(&(*oldhit3));
      debug9(printf("3' resolved on inside\n"));
      *oldhit3 = newhit3;
      *alts_resolve_3 = -1;
      changep = true;
    }
  }

  return changep;
}
#endif


static int
compute_insertlength (int *pair_relationship, Stage3pair_T this) {
  T hit5, hit3;
  int querylength5, querylength3;

  hit5 = this->hit5;
  hit3 = this->hit3;
  querylength5 = hit5->querylength;
  querylength3 = hit3->querylength;

  debug10(printf("Computing insertlength on %u..%u to %u..%u\n",
		 hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
		 hit3->genomicend - hit3->chroffset,hit3->genomicstart - hit3->chroffset));

  if (hit5->plusp == true && hit3->plusp == false) {
    /* Have 5-start..end and 3-end..start */
    /*   or 3-end..start and 5-start..end */

    *pair_relationship = 0;
    if (hit5->genomicend < hit3->genomicend) {
      return (hit3->genomicend - hit5->genomicend) + querylength5 + querylength3;
    } else if (hit3->genomicstart < hit5->genomicstart) {
      return (hit5->genomicstart - hit3->genomicstart) + querylength5 + querylength3;
    } else {
      return pair_insert_length_unpaired(hit5,hit3);
    }

  } else if (hit5->plusp == false && hit3->plusp == true) {
    /* Have 5-end..start and 3-start..end */
    /*   or 3-start..end and 5-end..start */

    *pair_relationship = 0;
    if (hit5->genomicstart < hit3->genomicstart) {
      return (hit3->genomicstart - hit5->genomicstart) + querylength5 + querylength3;
    } else if (hit3->genomicend < hit5->genomicend) {
      return (hit5->genomicend - hit3->genomicend) + querylength5 + querylength3;
    } else {
      return pair_insert_length_unpaired(hit5,hit3);
    }

  } else if (hit5->plusp == true) {
    /* Concordant directions on same chromosome (plus) */
    debug10(printf("Concordant on plus strand\n"));
    /* Have 5-start..end and 3-start..end */
    if (hit5->genomicend < hit3->genomicstart) {
      /* No overlap */
      *pair_relationship = +1;
      return (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
    } else {
      return pair_insert_length(&(*pair_relationship),hit5,hit3);
    }


  } else {
    /* Concordant directions on same chromosome (minus) */
    debug10(printf("Concordant on minus strand\n"));
    /* Have 3-end..start and 5-end..start */
    if (hit3->genomicstart < hit5->genomicend) {
      /* No overlap */
      *pair_relationship = -1;
      return (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
    } else {
      return pair_insert_length(&(*pair_relationship),hit5,hit3);
    }
  }
}


#ifdef RESOLVE_INSIDE_GENERAL
static void
resolve_insides (Stage3pair_T stage3pair, char *queryuc_ptr_5, char *queryuc_ptr_3,
		 Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		 Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		 Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		 Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
		 Listpool_T listpool) {
  T hit5, hit3;
  int querylength5, querylength3;
  bool changep;

  return;

  hit5 = stage3pair->hit5;
  hit3 = stage3pair->hit3;
  querylength5 = hit5->querylength;
  querylength3 = hit3->querylength;
    
  debug9(printf("Before resolve, hits are %p and %p\n",hit5,hit3));
  if (hit5->plusp == true && hit3->plusp == true) {
    changep = resolve_inside_general_splice_plus(&hit5,&hit3,&stage3pair->alts_resolve_5,&stage3pair->alts_resolve_3,
						 query5_compress_fwd,query3_compress_fwd,
						 queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
						 stage3pair->genestrand,pairpool,dynprogL,dynprogM,dynprogR,
						 oligoindices_minor,diagpool,cellpool,listpool);
  } else if (hit5->plusp == false && hit3->plusp == false) {
    changep = resolve_inside_general_splice_minus(&hit5,&hit3,&stage3pair->alts_resolve_5,&stage3pair->alts_resolve_3,
						  query5_compress_rev,query3_compress_rev,
						  queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
						  stage3pair->genestrand,pairpool,dynprogL,dynprogM,dynprogR,
						  oligoindices_minor,diagpool,cellpool,listpool);
  } else {
    changep = false;
  }
  
  if (changep == true) {
    debug9(printf("After resolve, hits are %p and %p\n",hit5,hit3));
    stage3pair->hit5 = hit5;
    stage3pair->hit3 = hit3;
    stage3pair->insertlength = compute_insertlength(&stage3pair->pair_relationship,stage3pair);
    
    /* Rest of this code is taken from the bottom of Stage3pair_new */
    
    stage3pair->score_overall = hit5->score_overall + hit3->score_overall;
    stage3pair->score_within_trim = hit5->score_within_trims + hit3->score_within_trims;
    
    stage3pair->nmatches_to_trims = hit5->nmatches_to_trims + hit3->nmatches_to_trims;
    stage3pair->nmatches_plus_spliced_trims = hit5->nmatches_plus_spliced_trims + hit3->nmatches_plus_spliced_trims;
    /* stage3pair->overlap_known_gene_p = false; -- initialized later when resolving multimappers */
    stage3pair->tally = -1L;
    
    stage3pair->low = (hit5->low < hit3->low) ? hit5->low : hit3->low;
    stage3pair->high = (hit5->high > hit3->high) ? hit5->high : hit3->high;
    
#if 0
    if (stage3pair->low > stage3pair->high) {
      fprintf(stderr,"stage3pair->low %u > stage3pair->high %u, hit5->chrnum %d\n",
	      stage3pair->low - stage3pair->chroffset,stage3pair->high - stage3pair->chroffset,hit5->chrnum);
      abort();
    }
#endif
    
    if (hit5->chrnum == 0 || hit3->chrnum == 0) {
      stage3pair->outerlength = querylength5 + querylength3;
    } else {
      stage3pair->outerlength = stage3pair->high - stage3pair->low;
    }
    
    stage3pair->nsplices = hit5->nsplices + hit3->nsplices;
    
    debug0(printf("Revised new pair %p from %p and %p\n",stage3pair,hit5,hit3));
    debug0(printf("  methods %s and %s\n",Method_string(hit5->method),Method_string(hit3->method)));
    debug0(printf("  sensedirs %d and %d\n",hit5->sensedir,hit3->sensedir));
    debug0(printf("  chrpos %u..%u and %u..%u\n",
		  hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
		  hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset));
    
    if (hit5->circularpos < 0 && hit3->circularpos < 0) {
      stage3pair->circularp = false;
    } else {
      stage3pair->circularp = true;
    }
    
    /* Fixing insertlength for circular pairs */
    if (stage3pair->insertlength > hit5->chrlength) {
      stage3pair->insertlength -= hit5->chrlength;
    }
    
    /* Note: the new hit5 or hit3 is guaranteed to have private5p or private3p set to true, respectively */
    if (hit5->circularalias == +1) {
      debug0(printf("Unaliasing 5' end\n"));
      unalias_circular(stage3pair->hit5);
    }
    
    if (hit3->circularalias == +1) {
      debug0(printf("Unaliasing 3' end\n"));
      unalias_circular(stage3pair->hit3);
    }
  }

  return;
}
#endif


#ifdef RESOLVE_INSIDE_GENERAL
List_T
Stage3pair_resolve_insides (List_T hitpairlist, char *queryuc_ptr_5, char *queryuc_ptr_3,
			    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			    Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
			    Listpool_T listpool, Hitlistpool_T hitlistpool) {
  List_T result = NULL, p;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  int querylength5, querylength3;
  bool changep;

  return hitpairlist;

  for (p = hitpairlist; p != NULL; p = List_next(p)) {
    stage3pair = (Stage3pair_T) List_head(p);
    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;
    querylength5 = hit5->querylength;
    querylength3 = hit3->querylength;
    
    debug9(printf("Before resolve, hits are %p and %p\n",hit5,hit3));
    if (hit5->plusp == true && hit3->plusp == true) {
      changep = resolve_inside_general_splice_plus(&hit5,&hit3,&stage3pair->alts_resolve_5,&stage3pair->alts_resolve_3,
						   query5_compress_fwd,query3_compress_fwd,
						   queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
						   stage3pair->genestrand,pairpool,dynprogL,dynprogM,dynprogR,
						   oligoindices_minor,diagpool,cellpool,listpool);
    } else if (hit5->plusp == false && hit3->plusp == false) {
      changep = resolve_inside_general_splice_minus(&hit5,&hit3,&stage3pair->alts_resolve_5,&stage3pair->alts_resolve_3,
						    query5_compress_rev,query3_compress_rev,
						    queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
						    stage3pair->genestrand,pairpool,dynprogL,dynprogM,dynprogR,
						    oligoindices_minor,diagpool,cellpool,listpool);
    } else {
      changep = false;
    }

    if (changep == true) {
      debug9(printf("After resolve, hits are %p and %p\n",hit5,hit3));
      stage3pair->hit5 = hit5;
      stage3pair->hit3 = hit3;
      stage3pair->insertlength = compute_insertlength(&stage3pair->pair_relationship,stage3pair);

      /* Rest of this code is taken from the bottom of Stage3pair_new */

      stage3pair->score_overall = hit5->score_overall + hit3->score_overall;
      stage3pair->score_within_trims = hit5->score_within_trims + hit3->score_within_trims;

      stage3pair->nmatches_to_trims = hit5->nmatches_to_trims + hit3->nmatches_to_trims;
      stage3pair->nmatches_plus_spliced_trims = hit5->nmatches_plus_spliced_trims + hit3->nmatches_plus_spliced_trims;

      /* stage3pair->overlap_known_gene_p = false; -- initialized later when resolving multimappers */
      stage3pair->tally = -1L;

      stage3pair->low = (hit5->low < hit3->low) ? hit5->low : hit3->low;
      stage3pair->high = (hit5->high > hit3->high) ? hit5->high : hit3->high;

#if 0
      if (stage3pair->low > stage3pair->high) {
	fprintf(stderr,"stage3pair->low %u > stage3pair->high %u, hit5->chrnum %d\n",
		stage3pair->low - stage3pair->chroffset,stage3pair->high - stage3pair->chroffset,hit5->chrnum);
	abort();
      }
#endif

      if (hit5->chrnum == 0 || hit3->chrnum == 0) {
	stage3pair->outerlength = querylength5 + querylength3;
      } else {
	stage3pair->outerlength = stage3pair->high - stage3pair->low;
      }

      stage3pair->nsplices = hit5->nsplices + hit3->nsplices;

      debug0(printf("Revised new pair %p from %p and %p\n",stage3pair,hit5,hit3));
      debug0(printf("  methods %s and %s\n",Method_string(hit5->method),Method_string(hit3->method)));
      debug0(printf("  sensedirs %d and %d\n",hit5->sensedir,hit3->sensedir));
      debug0(printf("  chrpos %u..%u and %u..%u\n",
		    hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
		    hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset));

      if (hit5->circularpos < 0 && hit3->circularpos < 0) {
	stage3pair->circularp = false;
      } else {
	stage3pair->circularp = true;
      }

      /* Fixing insertlength for circular pairs */
      if (stage3pair->insertlength > hit5->chrlength) {
	stage3pair->insertlength -= hit5->chrlength;
      }
      
      /* Note: the new hit5 or hit3 is guaranteed to have private5p or private3p set to true, respectively */
      if (hit5->circularalias == +1) {
	debug0(printf("Unaliasing 5' end\n"));
	unalias_circular(stage3pair->hit5);
      }

      if (hit3->circularalias == +1) {
	debug0(printf("Unaliasing 3' end\n"));
	unalias_circular(stage3pair->hit3);
      }
    }

    result = Hitlist_push(result,hitlistpool,(void *) stage3pair);
  }
    
  Hitlist_free(&hitpairlist);
  return result;
}
#endif



/* Should not set ambiguous flag in substrings, because resolution of
   an ambiguity depends on a particular pair of ends */

static void
resolve_inside_alts_splice_plus (int *alts_resolve_5, int *alts_resolve_3,
				 int *alts_status_inside, T hit5, T hit3, int querylength5, int querylength3) {
  Chrpos_T best_insertlength, insertlength;
  Univcoord_T genomicstart, genomicend;
  int besti5 = -1, besti3 = -1, i, j;
  int best_nmismatches, nmismatches;

  Substring_T substring5, substring3;
  Univcoord_T *end_alts_coords, *start_alts_coords;
  int *end_alts_nmismatches, *start_alts_nmismatches;
  int end_amb_length_5, start_amb_length_3;


  debug9(printf("resolve plus: hit5 %p (%s) and hit3 %p (%s)\n",
		hit5,Method_string(hit5->method),hit3,Method_string(hit3->method)));

  substring5 = (Substring_T) List_head(hit5->substrings_Nto1); /* the substring for concordance */
  debug9(printf("Testing substring5 %p %d..%d alts_p %d\n",
		substring5,Stage3end_substrings_querystart(hit5),Stage3end_substrings_queryend(hit5),
		Substring_has_alts_p(substring5)));

  substring3 = (Substring_T) List_head(hit3->substrings_1toN); /* the substring for concordance (was Nto1) */
  debug9(printf("Testing substring3 %p %d..%d alts_p %d\n",
		substring3,Stage3end_substrings_querystart(hit3),Stage3end_substrings_queryend(hit3),
		Substring_has_alts_p(substring3)));

  if (substring5 != NULL && Substring_has_alts_p(substring5) == true && 
      substring3 != NULL && Substring_has_alts_p(substring3) == true) {
    debug9(printf("Resolve plus case 1: Got alts at 5' and alts at 3':"));
    end_alts_coords = Substring_alts_coords(substring5);
    end_alts_nmismatches = Substring_alts_nmismatches(substring5);
    start_alts_coords = Substring_alts_coords(substring3);
    start_alts_nmismatches = Substring_alts_nmismatches(substring3);
    end_amb_length_5 = end_amb_length(hit5);
    start_amb_length_3 = start_amb_length(hit3);
    
    best_insertlength = (Chrpos_T) -1;
    best_nmismatches = querylength5 + querylength3;
    for (i = 0; i < Substring_alts_ncoords(substring5); i++) {
      genomicend = end_alts_coords[i] + end_amb_length_5;
      for (j = 0; j < Substring_alts_ncoords(substring3); j++) {
	genomicstart = start_alts_coords[j] - start_amb_length_3;
	debug9(printf(" %u,%u",(Chrpos_T) (genomicend - hit5->chroffset),(Chrpos_T) (genomicstart - hit3->chroffset)));
	if (genomicend < genomicstart) {
	  /* Look for valid insertlength */
	  insertlength = genomicstart - genomicend + querylength5 + querylength3;
	  debug9(printf(" (insertlength %u)",insertlength));

	  if (insertlength < best_insertlength) {
	    besti5 = i;
	    besti3 = j;
	    best_insertlength = insertlength;
	    best_nmismatches = end_alts_nmismatches[i] + start_alts_nmismatches[j];
	    debug9(printf("*"));
	  } else if (insertlength == best_insertlength &&
		     (nmismatches = end_alts_nmismatches[i] + start_alts_nmismatches[j]) < best_nmismatches) {
	    besti5 = i;
	    besti3 = j;
	    best_nmismatches = nmismatches;
	    debug9(printf("*"));
	  } else if (nmismatches == best_nmismatches) {
	    debug9(printf("tie"));
	  }
	}
      }
    }

    if (besti5 >= 0 && besti3 >= 0) {
      debug9(printf("\nBEST HAS INSERTLENGTH %u AND NMISMATCHES %d\n",best_insertlength,best_nmismatches));
      *alts_resolve_5 = besti5;
      *alts_resolve_3 = besti3;
      *alts_status_inside = ALTS_RESOLVED_BYLENGTH;
      hit5->genomicend = end_alts_coords[besti5] + end_amb_length_5;
      hit3->genomicstart = start_alts_coords[besti3] - start_amb_length_3;
    }
    debug9(printf("\n"));

  } else if (substring5 != NULL && Substring_has_alts_p(substring5) == true) {
    debug9(printf("Resolve plus case 2: Got alts at 5':"));
    end_alts_coords = Substring_alts_coords(substring5);
    end_alts_nmismatches = Substring_alts_nmismatches(substring5);
    end_amb_length_5 = end_amb_length(hit5);

    best_insertlength = (Chrpos_T) -1;
    best_nmismatches = querylength5;
    for (i = 0; i < Substring_alts_ncoords(substring5); i++) {
      genomicend = end_alts_coords[i] + end_amb_length_5;
      debug9(printf(" %u",(Chrpos_T) (genomicend - hit5->chroffset)));
      if (genomicend < hit3->genomicstart /*allow overlap*/+ querylength3) {
	/* Look for valid insertlength */
	insertlength = hit3->genomicstart - genomicend + querylength5 + querylength3;
	debug9(printf(" (insertlength %u)",insertlength));

	if (insertlength < best_insertlength) {
	  besti5 = i;
	  best_insertlength = insertlength;
	  best_nmismatches = end_alts_nmismatches[i];
	  debug9(printf("*"));
	} else if (insertlength == best_insertlength &&
		   (nmismatches = end_alts_nmismatches[i]) < best_nmismatches) {
	  besti5 = i;
	  best_nmismatches = nmismatches;
	  debug9(printf("*"));
	} else if (nmismatches == best_nmismatches) {
	  debug9(printf("tie"));
	}
      }
    }

    if (besti5 >= 0) {
      debug9(printf("\nBEST HAS INSERTLENGTH %u WITH NMISMATCHES %d\n",best_insertlength,best_nmismatches));
      *alts_resolve_5 = besti5;
      *alts_status_inside = ALTS_RESOLVED_BYLENGTH;
      hit5->genomicend = end_alts_coords[besti5] + end_amb_length_5;
    }
    debug9(printf("\n"));

  } else if (substring3 != NULL && Substring_has_alts_p(substring3) == true) {
    debug9(printf("Resolve plus case 3: Got alts at 3':"));
    start_alts_coords = Substring_alts_coords(substring3);
    start_alts_nmismatches = Substring_alts_nmismatches(substring3);
    start_amb_length_3 = start_amb_length(hit3);
    
    best_insertlength = (Chrpos_T) -1;
    best_nmismatches = querylength3;
    for (j = 0; j < Substring_alts_ncoords(substring3); j++) {
      genomicstart = start_alts_coords[j] - start_amb_length_3;
      debug9(printf(" %u",(Chrpos_T) (genomicstart - hit3->chroffset)));
      if (hit5->genomicend < genomicstart /*allow overlap*/+ querylength5) {
	/* Look for valid insertlength */
	insertlength = genomicstart - hit5->genomicend + querylength5 + querylength3;
	debug9(printf(" (insertlength %u)",insertlength));

	if (insertlength < best_insertlength) {
	  besti3 = j;
	  best_insertlength = insertlength;
	  best_nmismatches = start_alts_nmismatches[j];
	  debug9(printf("*"));
	} else if (insertlength == best_insertlength &&
		   (nmismatches = start_alts_nmismatches[j]) < best_nmismatches) {
	  besti3 = j;
	  best_nmismatches = nmismatches;
	  debug9(printf("*"));
	} else if (nmismatches == best_nmismatches) {
	  debug9(printf("tie"));
	}
      }
    }

    if (besti3 >= 0) {
      debug9(printf("\nBEST HAS INSERTLENGTH %u WITH NMISMATCHES %d\n",best_insertlength,best_nmismatches));
      *alts_resolve_3 = besti3;
      *alts_status_inside = ALTS_RESOLVED_BYLENGTH;
      hit3->genomicstart = start_alts_coords[besti3] - start_amb_length_3;
    }
    debug9(printf("\n"));
  }

  return;
}


static void
resolve_inside_alts_splice_minus (int *alts_resolve_5, int *alts_resolve_3,
				  int *alts_status_inside, T hit5, T hit3, int querylength5, int querylength3) {
  Chrpos_T best_insertlength, insertlength;
  Univcoord_T genomicstart, genomicend;
  int besti5 = -1, besti3 = -1, i, j;
  int best_nmismatches, nmismatches;

  Substring_T substring5, substring3;
  Univcoord_T *end_alts_coords, *start_alts_coords;
  int *end_alts_nmismatches, *start_alts_nmismatches;
  int end_amb_length_5, start_amb_length_3;


  debug9(printf("resolve minus: hit5 %p (%s) and hit3 %p (%s)\n",
		hit5,Method_string(hit5->method),hit3,Method_string(hit3->method)));

  substring5 = (Substring_T) List_head(hit5->substrings_Nto1); /* the substring for concordance */
  debug9(printf("Testing substring5 %p %d..%d alts_p %d\n",
		substring5,Stage3end_substrings_querystart(hit5),Stage3end_substrings_queryend(hit5),
		Substring_has_alts_p(substring5)));

  substring3 = (Substring_T) List_head(hit3->substrings_1toN); /* the substring for concordance */
  debug9(printf("Testing substring3 %p %d..%d alts_p %d\n",
		substring3,Stage3end_substrings_querystart(hit3),Stage3end_substrings_queryend(hit3),
		Substring_has_alts_p(substring3)));

  if (substring5 != NULL && Substring_has_alts_p(substring5) == true &&
      substring3 != NULL && Substring_has_alts_p(substring3) == true) {
    debug9(printf("Resolve minus case 1: Got alts at 5' and alts at 3':"));
    end_alts_coords = Substring_alts_coords(substring5);
    end_alts_nmismatches = Substring_alts_nmismatches(substring5);
    start_alts_coords = Substring_alts_coords(substring3);
    start_alts_nmismatches = Substring_alts_nmismatches(substring3);
    end_amb_length_5 = end_amb_length(hit5);
    start_amb_length_3 = start_amb_length(hit3);

    best_insertlength = (Chrpos_T) -1;
    best_nmismatches = querylength5 + querylength3;
    for (i = 0; i < Substring_alts_ncoords(substring5); i++) {
      genomicend = end_alts_coords[i] - end_amb_length_5;
      for (j = 0; j < Substring_alts_ncoords(substring3); j++) {
	genomicstart = start_alts_coords[j] + start_amb_length_3;
	debug9(printf(" %u,%u",(Chrpos_T) (genomicend - hit5->chroffset),(Chrpos_T) (genomicstart - hit3->chroffset)));
	if (genomicstart < genomicend) {
	  /* Look for valid insertlength */
	  insertlength = genomicend - genomicstart + querylength5 + querylength3;
	  debug9(printf(" (insertlength %u)",insertlength));

	  if (insertlength < best_insertlength) {
	    besti5 = i;
	    besti3 = j;
	    best_insertlength = insertlength;
	    best_nmismatches = end_alts_nmismatches[i] + start_alts_nmismatches[j];
	    debug9(printf("*"));
	  } else if (insertlength == best_insertlength &&
		     (nmismatches = end_alts_nmismatches[i] + start_alts_nmismatches[j]) < best_nmismatches) {
	    besti5 = i;
	    besti3 = j;
	    best_nmismatches = nmismatches;
	    debug9(printf("*"));
	  } else if (nmismatches == best_nmismatches) {
	    debug9(printf("tie"));
	  }
	}
      }
    }

    if (besti5 >= 0 && besti3 >= 0) {
      debug9(printf("\nBEST HAS INSERTLENGTH %u AND NMISMATCHES %d\n",best_insertlength,best_nmismatches));
      *alts_resolve_5 = besti5;
      *alts_resolve_3 = besti3;
      *alts_status_inside = ALTS_RESOLVED_BYLENGTH;
      hit5->genomicend = end_alts_coords[besti5] - end_amb_length_5;
      hit3->genomicstart = start_alts_coords[besti3] + start_amb_length_3;
    }
    debug9(printf("\n"));

  } else if (substring5 != NULL && Substring_has_alts_p(substring5) == true) {
    debug9(printf("Resolve minus case 2: Got alts at 5':"));
    end_alts_coords = Substring_alts_coords(substring5);
    end_alts_nmismatches = Substring_alts_nmismatches(substring5);
    end_amb_length_5 = end_amb_length(hit5);

    best_insertlength = (Chrpos_T) -1;
    best_nmismatches = querylength5;
    for (i = 0; i < Substring_alts_ncoords(substring5); i++) {
      genomicend = end_alts_coords[i] - end_amb_length_5;
      debug9(printf(" %u",(Chrpos_T) (genomicend - hit5->chroffset)));
      debug9(printf(" (%u <? %u + %d)",hit3->genomicstart,genomicend,querylength3));
      if (hit3->genomicstart < genomicend /*allow overlap*/+ querylength3) {
	/* Look for valid insertlength */
	insertlength = genomicend - hit3->genomicstart + querylength5 + querylength3;
	debug9(printf(" (insertlength %u)",insertlength));

	if (insertlength < best_insertlength) {
	  besti5 = i;
	  best_insertlength = insertlength;
	  best_nmismatches = end_alts_nmismatches[i];
	  debug9(printf("*"));
	} else if (insertlength == best_insertlength &&
		   (nmismatches = end_alts_nmismatches[i]) < best_nmismatches) {
	  besti5 = i;
	  best_nmismatches = nmismatches;
	  debug9(printf("*"));
	} else if (nmismatches == best_nmismatches) {
	  debug9(printf("tie"));
	}
      }
    }

    if (besti5 >= 0) {
      debug9(printf("\nBEST HAS INSERTLENGTH %u WITH NMISMATCHES %d\n",best_insertlength,best_nmismatches));
      *alts_resolve_5 = besti5;
      *alts_status_inside = ALTS_RESOLVED_BYLENGTH;
      hit5->genomicend = end_alts_coords[besti5] - end_amb_length_5;
    }
    debug9(printf("\n"));

  } else if (substring3 != NULL && Substring_has_alts_p(substring3) == true) {
    debug9(printf("Resolve minus case 3: Got alts at 3':"));
    start_alts_coords = Substring_alts_coords(substring3);
    start_alts_nmismatches = Substring_alts_nmismatches(substring3);
    start_amb_length_3 = start_amb_length(hit3);

    best_insertlength = (Chrpos_T) -1;
    best_nmismatches = querylength3;
    for (j = 0; j < Substring_alts_ncoords(substring3); j++) {
      genomicstart = start_alts_coords[j] + start_amb_length_3;
      debug9(printf(" %u",(Chrpos_T) (genomicstart - hit3->chroffset)));
      if (genomicstart < hit5->genomicend /*allow overlap*/+ querylength5) {
	/* Look for valid insertlength */
	insertlength = hit5->genomicend - genomicstart + querylength5 + querylength3;
	debug9(printf(" (insertlength %u)",insertlength));

	if (insertlength < best_insertlength) {
	  besti3 = j;
	  best_insertlength = insertlength;
	  best_nmismatches = start_alts_nmismatches[j];
	  debug9(printf("*"));
	} else if (insertlength == best_insertlength &&
		   (nmismatches = start_alts_nmismatches[j]) < best_nmismatches) {
	  besti3 = j;
	  best_nmismatches = nmismatches;
	  debug9(printf("*"));
	} else if (nmismatches == best_nmismatches) {
	  debug9(printf("tie"));
	}
      }
    }

    if (besti3 >= 0) {
      debug9(printf("\nBEST HAS INSERTLENGTH %u WITH NMISMATCHES %d\n",best_insertlength,best_nmismatches));
      *alts_resolve_3 = besti3;
      *alts_status_inside = ALTS_RESOLVED_BYLENGTH;
      hit3->genomicstart = start_alts_coords[besti3] + start_amb_length_3;
    }
    debug9(printf("\n"));
  }

  return;
}



static void
alias_circular (T hit) {
  Chrpos_T chrlength = hit->chrlength;
  List_T p;
  Substring_T substring;

  assert(hit->circularalias == -1);
  for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
    substring = (Substring_T) List_head(p);
    Substring_alias_circular(substring);
  }

  /* Doesn't fix hitpair->low and hitpair->high */
  hit->genomicstart += chrlength;
  hit->genomicend += chrlength;
  hit->low += chrlength;
  hit->high += chrlength;

  hit->circularalias = +1;

  return;
}


/* Previously allowed for private5p or private3p to be true.  But now
   always copying (because concordance procedure can delete hits), and
   so private5p and private3p are essentially true. */
Stage3pair_T
Stage3pair_new (T hit5, T hit3, int genestrand, Pairtype_T pairtype,
#ifdef RESOLVE_INSIDE_GENERAL
		char *queryuc_ptr_5, char *queryuc_ptr_3,
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		Listpool_T listpool, bool expect_concordant_p, bool transcriptome_guided_p) {
  Stage3pair_T new;
  /* Stage3end_T copy; */
  Substring_T substring1, substringN;

  /* int found_score = 0; */
  bool overreach5p, overreach3p;
  Chrpos_T pairmax;

  int querylength5 = hit5->querylength;
  int querylength3 = hit3->querylength;

  char *remap_sequence;
  int remap_seqlength;
  List_T transcripts;


  debug0(printf("\nStage3pair_new called with pairtype %s and chrnum %d, %d (effective %d, %d), expect_concordant_p %d\n",
		Pairtype_string(pairtype),hit5->chrnum,hit3->chrnum,hit5->effective_chrnum,hit3->effective_chrnum,expect_concordant_p));

  /* Always make a copy, because concordance procedure might delete the hit */
  hit5 = Stage3end_copy(hit5,listpool);
  hit3 = Stage3end_copy(hit3,listpool);

  new = (Stage3pair_T) MALLOC_OUT(sizeof(*new));

  if (pairtype == PAIRED_UNSPECIFIED || pairtype == UNSPECIFIED) {
    /* Can get here from running GMAP improvement on a paired result */
    pairtype = Stage3_determine_pairtype(hit5,hit3,/*stage3pair*/NULL);
    debug10(printf("  Changing pairtype to %s\n",Pairtype_string(pairtype)));
    if (pairtype == CONCORDANT) {
      expect_concordant_p = true;
    }
  }
  new->pairtype = pairtype;
  new->genestrand = genestrand;

  new->alts_resolve_5 = -1;
  new->alts_resolve_3 = -1;
  new->alts_status_inside = ALTS_NOT_AMBIGUOUS;


#if 0
  new->mapq_loglik = hit5->mapq_loglik + hit3->mapq_loglik;
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  if (hit5->plusp == true && hit3->plusp == false) {
    debug10(printf("plus/minus\n"));
    new->dir = 0;
    
    /* Have 5-start..end and 3-end..start */
    /*   or 3-end..start and 5-start..end */

    new->pair_relationship = 0;
    if (hit5->genomicend < hit3->genomicend) {
      new->insertlength = (hit3->genomicend - hit5->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else if (hit3->genomicstart < hit5->genomicstart) {
      new->insertlength = (hit5->genomicstart - hit3->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else {
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->plusp == false && hit3->plusp == true) {
    debug10(printf("minus/plus\n"));
    new->dir = 0;
    
    /* Have 5-end..start and 3-start..end */
    /*   or 3-start..end and 5-end..start */

    new->pair_relationship = 0;
    if (hit5->genomicstart < hit3->genomicstart) {
      new->insertlength = (hit3->genomicstart - hit5->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else if (hit3->genomicend < hit5->genomicend) {
      new->insertlength = (hit5->genomicend - hit3->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else {
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->plusp == true) {
    /* Concordant directions on same chromosome (plus) */
    debug10(printf("*Concordant on plus strand\n"));
    new->dir = +1;

    if (expect_concordant_p == true) {
      overreach5p = overreach3p = false;
      if (hit5->hittype == SPLICE) {

	substringN = (Substring_T) List_head(hit5->substrings_Nto1);
	if (Substring_alignstart_trim(substringN) > hit3->genomicend) {
	  substring1 = (Substring_T) List_head(hit5->substrings_1toN);
	  if (Substring_alignend_trim(substring1) < hit3->genomicstart) {
	    overreach5p = true;
	  }
	}
      }
      if (hit3->hittype == SPLICE) {
	substring1 = (Substring_T) List_head(hit3->substrings_1toN);
	if (Substring_alignend_trim(substring1) < hit5->genomicstart) {
	  substringN = (Substring_T) List_head(hit3->substrings_Nto1);
	  if (Substring_alignstart_trim(substringN) > hit5->genomicend) {
	    overreach3p = true;
	  }
	}
      }

      if (overreach5p == true || overreach3p == true) {
	/* Either overreach */
	debug0(printf("  Returning NULL because of dual overreach\n"));
	Stage3end_free(&hit5);	/* This was the copy */
	Stage3end_free(&hit3);	/* This was the copy */
	FREE_OUT(new);
	return (Stage3pair_T) NULL;

#if 0
      } else if (overreach5p == true) {
	/* Overreach of hit5 */
	debug9(printf("Overreach of hit5 of type SPLICE.  Removing substring2\n"));
	if (hit5->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_bothdiff(hit5->substring1),
				      /*nmismatches_acceptor*/0,/*donor*/hit5->substring1,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*alts_coords_donor*/NULL,/*alts_coords_acceptor*/NULL,
				      /*alts_nmismatches_donor*/NULL,/*alts_nmismatches_acceptor*/NULL,
				      /*alts_probs_donor*/NULL,/*alts_probs_acceptor*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				      /*sensedir*/hit5->sensedir,listpool,hit5->method,hit5->level);
	} else if (hit5->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_bothdiff(hit5->substring1),/*donor*/NULL,
				      /*acceptor*/hit5->substring1,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*alts_coords_donor*/NULL,/*alts_coords_acceptor*/NULL,
				      /*alts_nmismatches_donor*/NULL,/*alts_nmismatches_acceptor*/NULL,
				      /*alts_probs_donor*/NULL,/*alts_probs_acceptor*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				      /*sensedir*/hit5->sensedir,listpool,hit5->method,hit5->level);
	} else {
	  abort();
	}
	Stage3end_free(&hit5);	/* This was the copy */
	hit5 = copy;

      } else if (overreach3p == true) {
	/* Overreach of hit3 */
	debug9(printf("Overreach of hit3 of type SPLICE.  Removing substring1\n"));
	if (hit3->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_bothdiff(hit3->substring2),/*donor*/NULL,
				      /*acceptor*/hit3->substring2,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*alts_coords_donor*/NULL,/*alts_coords_acceptor*/NULL,
				      /*alts_nmismatches_donor*/NULL,/*alts_nmismatches_acceptor*/NULL,
				      /*alts_probs_donor*/NULL,/*alts_probs_acceptor*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				      /*sensedir*/hit3->sensedir,listpool,hit3->method,hit3->level);
	} else if (hit3->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_bothdiff(hit3->substring2),
				      /*nmismatches_acceptor*/0,/*donor*/hit3->substring2,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*alts_coords_donor*/NULL,/*alts_coords_acceptor*/NULL,
				      /*alts_nmismatches_donor*/NULL,/*alts_nmismatches_acceptor*/NULL,
				      /*alts_probs_donor*/NULL,/*alts_probs_acceptor*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				      /*sensedir*/hit3->sensedir,listpool,hit3->method,hit3->level);
	} else {
	  abort();
	}
	Stage3end_free(&hit3);	/* This was the copy */
	hit3 = copy;
#endif
      }

      /* Try to resolve ambiguity on inside of concordant ends */
      debug9(printf("Calling resolve_inside_alts_splice_plus\n"));
      resolve_inside_alts_splice_plus(&new->alts_resolve_5,&new->alts_resolve_3,
				      &new->alts_status_inside,hit5,hit3,querylength5,querylength3);
      debug9(printf("For pair %p (%p and %p), set alts_resolve_5 to be %d and alts_resolve_3 to be %d\n",
		    new,hit5,hit3,new->alts_resolve_5,new->alts_resolve_3));
    }

    /* Have 5-start..end and 3-start..end */
    if (hit5->genomicend < hit3->genomicstart) {
      /* No overlap */
      new->pair_relationship = +1;
      new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		     new->insertlength,hit3->genomicstart - hit3->chroffset,
		     hit5->genomicend - hit5->chroffset,querylength5,querylength3));
#if 0
    } else if (hit5->genomicend > hit3->genomicend + SUBSUMPTION_SLOP) {
      /* hit5 subsumes hit3 */
      debug10(printf("plus, subsumption %u > %u\n",
		     hit5->genomicend - hit5->chroffset,hit3->genomicend - hit3->chroffset));
      new->pair_relationship = 0;
      new->insertlength = 0;
      new->insertlength_expected_sign = false;
#endif
    } else {
      new->insertlength = pair_insert_length(&new->pair_relationship,hit5,hit3);
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    }


  } else {
    /* Concordant directions on same chromosome (minus) */
    debug10(printf("*Concordant on minus strand\n"));
    new->dir = -1;

    if (expect_concordant_p == true) {
      overreach5p = overreach3p = false;
      if (hit5->hittype == SPLICE) {
	debug10(printf("Have splice on 5' end\n"));
	substringN = (Substring_T) List_head(hit5->substrings_Nto1);
	if (Substring_alignstart_trim(substringN) < hit3->genomicend) {
	  substring1 = (Substring_T) List_head(hit5->substrings_1toN);
	  if (Substring_alignend_trim(substring1) > hit3->genomicstart) {
	    overreach5p = true;
	  }
	}
      }
      if (hit3->hittype == SPLICE) {
	debug10(printf("Have splice on 3' end\n"));
	substring1 = (Substring_T) List_head(hit3->substrings_1toN);
	if (Substring_alignend_trim(substring1) > hit5->genomicstart) {
	  substringN = (Substring_T) List_head(hit3->substrings_Nto1);
	  if (Substring_alignstart_trim(substringN) < hit5->genomicend) {
	    overreach3p = true;
	  }
	}
      }

      if (overreach5p == true || overreach3p == true) {
	/* Either overreach */
	debug0(printf("  Returning NULL because of dual overreach\n"));
	Stage3end_free(&hit5); /* This was the copy */
	Stage3end_free(&hit3); /* This was the copy */
	FREE_OUT(new);
	return (Stage3pair_T) NULL;

#if 0
      } else if (overreach5p == true) {
	/* Overreach of hit5 */
	debug9(printf("Overreach of hit5 of type SPLICE.  Removing substring2\n"));
	if (hit5->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_bothdiff(hit5->substring1),
				      /*nmismatches_acceptor*/0,/*donor*/hit5->substring1,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*alts_coords_donor*/NULL,/*alts_coords_acceptor*/NULL,
				      /*alts_nmismatches_donor*/NULL,/*alts_nmismatches_acceptor*/NULL,
				      /*alts_probs_donor*/NULL,/*alts_probs_acceptor*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				      /*sensedir*/hit5->sensedir,listpool,hit5->method,hit5->level);
	} else if (hit5->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_bothdiff(hit5->substring1),/*donor*/NULL,
				      /*acceptor*/hit5->substring1,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*alts_coords_donor*/NULL,/*alts_coords_acceptor*/NULL,
				      /*alts_nmismatches_donor*/NULL,/*alts_nmismatches_acceptor*/NULL,
				      /*alts_probs_donor*/NULL,/*alts_probs_acceptor*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				      /*sensedir*/hit5->sensedir,listpool,hit5->method,hit5->level);
	} else {
	  abort();
	}
	Stage3end_free(&hit5);	/* This was the copy */
	hit5 = copy;

      } else if (overreach3p == true) {
	/* Overreach of hit3 */
	debug9(printf("Overreach of hit3 of type SPLICE.  Removing substring1\n"));
	if (hit3->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_bothdiff(hit3->substring2),/*donor*/NULL,
				      /*acceptor*/hit3->substring2,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*alts_coords_donor*/NULL,/*alts_coords_acceptor*/NULL,
				      /*alts_nmismatches_donor*/NULL,/*alts_nmismatches_acceptor*/NULL,
				      /*alts_probs_donor*/NULL,/*alts_probs_acceptor*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				      /*sensedir*/hit3->sensedir,listpool,hit3->method,hit3->level);
	} else if (hit3->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_bothdiff(hit3->substring2),
				      /*nmismatches_acceptor*/0,/*donor*/hit3->substring2,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,/*amb_length*/0,/*amb_prob*/0.0,
				      /*alts_coords_donor*/NULL,/*alts_coords_acceptor*/NULL,
				      /*alts_nmismatches_donor*/NULL,/*alts_nmismatches_acceptor*/NULL,
				      /*alts_probs_donor*/NULL,/*alts_probs_acceptor*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				      /*sensedir*/hit3->sensedir,listpool,hit3->method,hit3->level);
	} else {
	  abort();
	}
	Stage3end_free(&hit3);	/* This was the copy */
	hit3 = copy;
#endif
      }

      /* Try to resolve ambiguity on inside of concordant ends */
      debug9(printf("Calling resolve_inside_alts_splice_minus\n"));
      resolve_inside_alts_splice_minus(&new->alts_resolve_5,&new->alts_resolve_3,
				       &new->alts_status_inside,hit5,hit3,querylength5,querylength3);
      debug9(printf("For pair %p (%p and %p), set alts_resolve_5 to be %d and alts_resolve_3 to be %d\n",
		    new,hit5,hit3,new->alts_resolve_5,new->alts_resolve_3));
    }

    /* Have 3-end..start and 5-end..start */
    if (hit3->genomicstart < hit5->genomicend) {
      /* No overlap */
      new->pair_relationship = -1;
      new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		     new->insertlength,hit5->genomicend - hit5->chroffset,
		     hit3->genomicstart - hit3->chroffset,querylength5,querylength3));
#if 0
    } else if (hit3->genomicstart > hit5->genomicstart + SUBSUMPTION_SLOP) {
      /* hit3 subsumes hit5 */
      debug10(printf("minus, subsumption %u > %u\n",
		     hit3->genomicstart - hit3->chroffset,hit5->genomicstart - hit5->chroffset));
      new->pair_relationship = 0;
      new->insertlength = 0;
      new->insertlength_expected_sign = false;
#endif
    } else {
      new->insertlength = pair_insert_length(&new->pair_relationship,hit5,hit3);
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    }
  }

  debug10(printf("\nGot initial insertlength of %d\n",new->insertlength));

  new->hit5 = hit5;
  new->hit3 = hit3;

#ifdef RESOLVE_INSIDE_GENERAL
  if (expect_concordant_p == true && transcriptome_guided_p == false) {
    /* Previously performed resolve_insides to final hitlist, but
       results are better if it is done here */
    resolve_insides(new,queryuc_ptr_5,queryuc_ptr_3,
		    query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
		    pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,listpool);
    hit5 = new->hit5;		/* May have changed in resolve_insides */
    hit3 = new->hit3;		/* May have changed in resolve_insides */
    debug10(printf("\nAfter resolve_insides, got insertlength of %d\n",new->insertlength));
  }
#endif

  /* Was new->insertlength <= 0, but this eliminates legitimate overlaps */
  /* Was new->insertlength < -pairmax, but this allows overreach */
  if (new->insertlength <= 0) {	/* Not possible, since insertlength is unsigned */
    /* Not concordant */
#ifdef USE_BINGO
    new->absdifflength_bingo_p = false;
#endif
#ifdef USE_ABSDIFFLENGTH
    new->absdifflength = (Chrpos_T) -1;
#endif

    if (expect_concordant_p == true) {
      debug0(printf("  Returning NULL, because insertlength %u, so not concordant\n",new->insertlength));
      Stage3end_free(&hit5);	/* This was the copy */
      Stage3end_free(&hit3);	/* This was the copy */
      FREE_OUT(new);
      return (Stage3pair_T) NULL;
    }

  } else {
    if (transcriptome_guided_p == true) {
      pairmax = (Chrpos_T) -1;
    } else if (circularp[hit5->effective_chrnum] == true) {
      pairmax = pairmax_circular;
    } else {
      pairmax = pairmax_linear;
    }
    if (new->insertlength > pairmax && expect_concordant_p == true) {
      debug0(printf("  Returning NULL because insertlength %u > pairmax %d\n",new->insertlength,pairmax));
      Stage3end_free(&hit5);	/* This was the copy */
      Stage3end_free(&hit3);	/* This was the copy */
      FREE_OUT(new);
      return (Stage3pair_T) NULL;
      
    } else {
#ifdef USE_ABSDIFFLENGTH
      if (new->insertlength < expected_pairlength) {
	new->absdifflength = expected_pairlength - new->insertlength;
      } else {
	new->absdifflength = new->insertlength - expected_pairlength;
      }
#endif
#ifdef USE_BINGO
      if (new->absdifflength <= pairlength_deviation) {
	new->absdifflength_bingo_p = true;
      } else {
	new->absdifflength_bingo_p = false;
      }
#endif
    }
  }

  if (SENSE_CONSISTENT_P(hit5->sensedir_for_concordance,hit3->sensedir_for_concordance)) {
    debug0(printf("senses %d and %d are consistent\n",hit5->sensedir_for_concordance,hit3->sensedir_for_concordance));
    new->sense_consistent_p = true;

  } else if (expect_concordant_p == true) {
    debug0(printf("  Returning NULL, because senses are not consistent\n"));
    Stage3end_free(&hit5); 	/* This was the copy */
    Stage3end_free(&hit3);	/* This was the copy */
    FREE_OUT(new);
    return (Stage3pair_T) NULL;

  } else {
    debug0(printf("senses are inconsistent, but allowable\n"));
    new->sense_consistent_p = false;
  }

  /* Do not alter score, so the alignmnent terminates at the known splice site  */
  new->score_overall = hit5->score_overall + hit3->score_overall;
  new->score_within_trims = hit5->score_within_trims + hit3->score_within_trims;

  new->nmatches_to_trims = hit5->nmatches_to_trims + hit3->nmatches_to_trims;
  new->nmatches_plus_spliced_trims = hit5->nmatches_plus_spliced_trims + hit3->nmatches_plus_spliced_trims;

  /* new->overlap_known_gene_p = false; -- initialized later when resolving multimappers */
  new->tally = -1L;

  new->low = (hit5->low < hit3->low) ? hit5->low : hit3->low;
  new->high = (hit5->high > hit3->high) ? hit5->high : hit3->high;

#if 0
  if (new->low > new->high) {
    fprintf(stderr,"new->low %u > new->high %u, hit5->chrnum %d\n",
	    new->low - new->chroffset,new->high - new->chroffset,hit5->chrnum);
    abort();
  }
#endif

  if (hit5->chrnum == 0 || hit3->chrnum == 0) {
    new->outerlength = querylength5 + querylength3;
  } else {
    new->outerlength = new->high - new->low;
  }

  if (expect_concordant_p == true) {
    hit5->paired_usedp = true;
    hit3->paired_usedp = true;
  }

  new->nsplices = hit5->nsplices + hit3->nsplices;

  debug0(printf("Created new pair %p from %p and %p (nmatches_to_trims %d+%d)\n",
		new,hit5,hit3,hit5->nmatches_to_trims,hit3->nmatches_to_trims));
  debug0(printf("  methods %s and %s\n",Method_string(hit5->method),Method_string(hit3->method)));
  debug0(printf("  sensedirs %d and %d\n",hit5->sensedir,hit3->sensedir));
  debug0(printf("  chrpos %u..%u and %u..%u\n",
		hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
		hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset));

  if (hit5->circularpos < 0 && hit3->circularpos < 0) {
    new->circularp = false;
  } else {
    new->circularp = true;
  }

  /* Fixing insertlength for circular pairs */
  if (new->insertlength > hit5->chrlength) {
    new->insertlength -= hit5->chrlength;
  }

  if (hit5->circularalias == +1) {
    debug0(printf("Unaliasing 5' end\n"));
    unalias_circular(hit5);
  }

  if (hit3->circularalias == +1) {
    debug0(printf("Unaliasing 3' end\n"));
    unalias_circular(hit3);
  }

  if (remap_transcriptome_p == false) {
    /* Do not remap */

  } else if (hit5->transcripts != NULL && hit3->transcripts != NULL) {
    /* No need to remap */

  } else if (hit5->transcripts != NULL && hit3->transcripts == NULL) {
    debug0(printf("Remapping 3' end to transcriptome to match 5' end at %d:%u..%u\n",
		  hit5->chrnum,hit5->low - hit5->chroffset,hit5->high - hit5->chroffset));
    remap_sequence = Stage3end_substrings_genomic_sequence(&remap_seqlength,hit3,genomecomp);
    debug0(printf("%s\n",remap_sequence));

    if ((transcripts = Kmer_remap_transcriptome(remap_sequence,remap_seqlength,hit3->chrnum,
						/*lowbound*/hit3->low - hit3->chroffset,
						/*highbound*/hit3->high - hit3->chroffset,
						transcript_iit,transcriptomebits,transcriptome)) != NULL) {
      hit3->transcripts = transcripts;
    }
    FREE(remap_sequence);

  } else if (hit5->transcripts == NULL && hit3->transcripts != NULL) {
    debug0(printf("Remapping 5' end to transcriptome to match 3' end at %d:%u..%u\n",
		  hit3->chrnum,hit3->low - hit3->chroffset,hit3->high - hit3->chroffset));

    remap_sequence = Stage3end_substrings_genomic_sequence(&remap_seqlength,hit5,genomecomp);
    debug0(printf("%s\n",remap_sequence));
    if ((transcripts = Kmer_remap_transcriptome(remap_sequence,remap_seqlength,hit5->chrnum,
						/*lowbound*/hit5->low - hit5->chroffset,
						/*highbound*/hit5->high - hit5->chroffset,
						transcript_iit,transcriptomebits,transcriptome)) != NULL) {
      hit5->transcripts = transcripts;
    }
    FREE(remap_sequence);
  }

  /* Need this in addition to Stage3end_filter_concordant_tr, to
     eliminate any inconsistent transcripts */
  Transcript_concordance(&new->transcripts5,&new->transcripts3,hit5->transcripts,hit3->transcripts);
  pairtype = Stage3_determine_pairtype(hit5,hit3,/*stage3pair*/new);
  debug0(printf("%d transcripts5, %d transcripts3\n",List_length(new->transcripts5),List_length(new->transcripts3)));

  /* assert((int) new->insertlength >= 0); */
  return new;
}


/* Used for eliminating exact duplicates.  Also sorts secondarily by hittype. */
static int
hitpair_sort_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;
  
  Univcoord_T x_hit5_high, x_hit5_low, y_hit5_high, y_hit5_low;
  Univcoord_T x_hit3_high, x_hit3_low, y_hit3_high, y_hit3_low;
  Univcoord_T x_low, x_high, y_low, y_high;
  
  debug8(printf("  Comparing (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), circularalias %d|%d, nmatches: %d (%d to_trims), amb_lengths %d and %d, sensedirs %d-%d, score %f+%f\n",
		Pairtype_string(x->pairtype),Method_string(x->hit5->method),
		Method_string(x->hit3->method),x,
		x->hit5->low - x->hit5->chroffset,x->hit5->high - x->hit5->chroffset,
		x->hit3->low - x->hit3->chroffset,x->hit3->high - x->hit3->chroffset,
		x->dir,x->hit5->circularalias,x->hit3->circularalias,x->nmatches_plus_spliced_trims,x->nmatches_to_trims,
		amb_length(x->hit5),amb_length(x->hit3),x->hit5->sensedir,x->hit3->sensedir,
		x->hit5->splice_score,x->hit3->splice_score));

  debug8(printf("       with (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), circularalias %d|%d, nmatches: %d (%d to_trims), amb_lengths %d and %d, sensedirs %d-%d, score %f+%f\n",
		Pairtype_string(y->pairtype),Method_string(y->hit5->method),
		Method_string(y->hit3->method),y,
		y->hit5->low - y->hit5->chroffset,y->hit5->high - y->hit5->chroffset,
		y->hit3->low - y->hit3->chroffset,y->hit3->high - y->hit3->chroffset,
		y->dir,y->hit5->circularalias,y->hit3->circularalias,y->nmatches_plus_spliced_trims,y->nmatches_to_trims,
		amb_length(y->hit5),amb_length(y->hit3),y->hit5->sensedir,y->hit3->sensedir,
		y->hit5->splice_score,y->hit3->splice_score));

  x_hit5_low = normalize_coord(x->hit5->low,x->hit5->circularalias,x->hit5->chrlength);
  x_hit5_high = normalize_coord(x->hit5->high,x->hit5->circularalias,x->hit5->chrlength);

  x_hit3_low = normalize_coord(x->hit3->low,x->hit3->circularalias,x->hit3->chrlength);
  x_hit3_high = normalize_coord(x->hit3->high,x->hit3->circularalias,x->hit3->chrlength);

  x_low = (x_hit5_low < x_hit3_low) ? x_hit5_low : x_hit3_low;
  x_high = (x_hit5_high > x_hit3_high) ? x_hit5_high : x_hit3_high;


  y_hit5_low = normalize_coord(y->hit5->low,y->hit5->circularalias,y->hit5->chrlength);
  y_hit5_high = normalize_coord(y->hit5->high,y->hit5->circularalias,y->hit5->chrlength);

  y_hit3_low = normalize_coord(y->hit3->low,y->hit3->circularalias,y->hit3->chrlength);
  y_hit3_high = normalize_coord(y->hit3->high,y->hit3->circularalias,y->hit3->chrlength);

  y_low = (y_hit5_low < y_hit3_low) ? y_hit5_low : y_hit3_low;
  y_high = (y_hit5_high > y_hit3_high) ? y_hit5_high : y_hit3_high;


  if (x->dir != 0 && y->dir == 0) {
    return -1;
  } else if (x->dir == 0 && y->dir != 0) {
    return +1;
  } else if (x->dir > 0 && y->dir < 0) {
    return -1;
  } else if (x->dir < 0 && y->dir > 0) {
    return +1;

#if 0
  } else if (x->high < y->low) {
    return -1;
  } else if (y->high < x->low) {
    return +1;

  } else if (x->hit5->high < y->hit5->low) {
    return -1;
  } else if (y->hit5->high < x->hit5->low) {
    return +1;

  } else if (x->hit3->high < y->hit3->low) {
    return -1;
  } else if (y->hit3->high < x->hit3->low) {
    return +1;
#else
    /* low to high pattern needed for finding overlaps */
  } else if (x_low < y_low) {
    debug8(printf("Returning -1 for low\n"));
    return -1;
  } else if (y_low < x_low) {
    debug8(printf("Returning +1 for low\n"));
    return +1;

  } else if (x_high > y_high) {
    debug8(printf("Returning -1 for high\n"));
    return -1;
  } else if (y_high > x_high) {
    debug8(printf("Returning +1 for high\n"));
    return +1;

    /* Need to check inside ends to avoid declaring unequal hitpairs equal */
  } else if (x_hit5_low < y_hit5_low) {
    return -1;
  } else if (y_hit5_low < x_hit5_low) {
    return +1;

  } else if (x_hit5_high < y_hit5_high) {
    return -1;
  } else if (y_hit5_high < x_hit5_high) {
    return +1;

  } else if (x_hit3_low < y_hit3_low) {
    return -1;
  } else if (y_hit3_low < x_hit3_low) {
    return +1;

  } else if (x_hit3_high < y_hit3_high) {
    return -1;
  } else if (y_hit3_high < x_hit3_high) {
    return +1;
#endif


  } else if (x->score_within_trims < y->score_within_trims) {
    return -1;
  } else if (y->score_within_trims < x->score_within_trims) {
    return +1;
  } else if (x->nmatches_plus_spliced_trims > y->nmatches_plus_spliced_trims) {
    return -1;
  } else if (y->nmatches_plus_spliced_trims > x->nmatches_plus_spliced_trims) {
    return +1;
#if 0
  } else if (x->nmatches_to_trims > y->nmatches_to_trims) {
    return -1;
  } else if (y->nmatches_to_trims > x->nmatches_to_trims) {
    return +1;
#endif

#if 0
  } else if (x->nsplices < y->nsplices) {
    return -1;
  } else if (y->nsplices < x->nsplices) {
    return +1;
#endif

  } else if (x->alts_status_inside < y->alts_status_inside) {
    return -1;
  } else if (y->alts_status_inside < x->alts_status_inside) {
    return +1;

#if 0
    /* Hittype should be irrelevant */
  } else if (x->hit5->hittype < y->hit5->hittype) {
    return -1;
  } else if (y->hit5->hittype < x->hit5->hittype) {
    return +1;
  } else if (x->hit3->hittype < y->hit3->hittype) {
    return -1;
  } else if (y->hit3->hittype < x->hit3->hittype) {
    return +1;
#endif

#if 0
    /* Hittype should be irrelevant, except for transcriptome getting priority */
  } else if (x->hit5->hittype == TRANSCRIPTOME && x->hit3->hittype == TRANSCRIPTOME &&
	     (y->hit5->hittype != TRANSCRIPTOME || y->hit3->hittype != TRANSCRIPTOME)) {
    return -1;
  } else if (y->hit5->hittype == TRANSCRIPTOME && y->hit3->hittype == TRANSCRIPTOME &&
	     (x->hit5->hittype != TRANSCRIPTOME || x->hit3->hittype != TRANSCRIPTOME)) {
    return +1;
#endif

#if 0
  } else if ((x->alts_resolve_5 != -1 && x->alts_resolve_3 != -1) &&
	     (y->alts_resolve_5 == -1 || y->alts_resolve_3 == -1)) {
    /* x is resolved, y is ambiguous.  x wins */
    return -1;
  } else if ((y->alts_resolve_5 != -1 && y->alts_resolve_3 != -1) &&
	     (x->alts_resolve_5 == -1 || x->alts_resolve_3 == -1)) {
    /* y is resolved, x is ambiguous.  y wins */
    return +1;
#endif

#if 0
  } else if (x->hit5->start_amb_length + x->hit5->end_amb_length +
	     x->hit3->start_amb_length + x->hit3->end_amb_length == 0 &&
	     y->hit5->start_amb_length + y->hit5->end_amb_length +
	     y->hit3->start_amb_length + y->hit3->end_amb_length > 0) {
    /* x is resolved, y is ambiguous.  x wins */
    return -1;
  } else if (y->hit5->start_amb_length + y->hit5->end_amb_length +
	     y->hit3->start_amb_length + y->hit3->end_amb_length == 0 &&
	     x->hit5->start_amb_length + x->hit5->end_amb_length +
	     x->hit3->start_amb_length + x->hit3->end_amb_length > 0) {
    /* y is resolved, x is ambiguous.  y wins */
    return +1;
#endif

  } else if (x->sense_consistent_p == true && y->sense_consistent_p == false) {
    debug8(printf(" => loses by sense_consistent_p\n"));
    return -1;
  } else if (x->sense_consistent_p == false && y->sense_consistent_p == true) {
    debug8(printf(" => wins by sense_consistent_p\n"));
    return +1;

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

#if 0
  } else if (x->sense_consistent_p == true) {
    if ((x->hit5->sensedir_for_concordance != 0 || x->hit3->sensedir_for_concordance != 0) &&
	(y->hit5->sensedir_for_concordance == 0 && y->hit3->sensedir_for_concordance == 0)) {
      return -1;
    } else if ((y->hit5->sensedir_for_concordance != 0 || y->hit3->sensedir_for_concordance != 0) &&
	       (x->hit5->sensedir_for_concordance == 0 && x->hit3->sensedir_for_concordance == 0)) {
      return +1;
    } else {
      return 0;
    }
#endif

  } else if (x->hit5->splice_score + x->hit3->splice_score >
	     y->hit5->splice_score + y->hit3->splice_score) {
    debug8(printf(" => loses by splice score\n"));
    return -1;

  } else if (y->hit5->splice_score + y->hit3->splice_score >
	     x->hit5->splice_score + x->hit3->splice_score) {
    debug8(printf(" => wins by splice score\n"));
    return +1;

  } else {
    debug8(printf(" => identical for sorting purposes\n"));
    return 0;
  }
}


/* Same as hitpair_sort_cmp, except for hittype, nmatches_to_trims, and indel_low */
static int
hitpair_equiv_cmp (Stage3pair_T x, Stage3pair_T y) {
  Univcoord_T x_hit5_high, x_hit5_low, y_hit5_high, y_hit5_low;
  Univcoord_T x_hit3_high, x_hit3_low, y_hit3_high, y_hit3_low;
  Univcoord_T x_low, x_high, y_low, y_high;
  
  x_hit5_low = normalize_coord(x->hit5->low,x->hit5->circularalias,x->hit5->chrlength);
  x_hit5_high = normalize_coord(x->hit5->high,x->hit5->circularalias,x->hit5->chrlength);

  x_hit3_low = normalize_coord(x->hit3->low,x->hit3->circularalias,x->hit3->chrlength);
  x_hit3_high = normalize_coord(x->hit3->high,x->hit3->circularalias,x->hit3->chrlength);

  x_low = (x_hit5_low < x_hit3_low) ? x_hit5_low : x_hit3_low;
  x_high = (x_hit5_high > x_hit3_high) ? x_hit5_high : x_hit3_high;


  y_hit5_low = normalize_coord(y->hit5->low,y->hit5->circularalias,y->hit5->chrlength);
  y_hit5_high = normalize_coord(y->hit5->high,y->hit5->circularalias,y->hit5->chrlength);

  y_hit3_low = normalize_coord(y->hit3->low,y->hit3->circularalias,y->hit3->chrlength);
  y_hit3_high = normalize_coord(y->hit3->high,y->hit3->circularalias,y->hit3->chrlength);

  y_low = (y_hit5_low < y_hit3_low) ? y_hit5_low : y_hit3_low;
  y_high = (y_hit5_high > y_hit3_high) ? y_hit5_high : y_hit3_high;


  if (x->dir != 0 && y->dir == 0) {
    return -1;
  } else if (x->dir == 0 && y->dir != 0) {
    return +1;
  } else if (x->dir > 0 && y->dir < 0) {
    return -1;
  } else if (x->dir < 0 && y->dir > 0) {
    return +1;
  } else if (x_low < y_low) {
    return -1;
  } else if (y_low < x_low) {
    return +1;
  } else if (x_high < y_high) {
    return -1;
  } else if (y_high < x_high) {
    return +1;

  } else if (x_hit5_low < y_hit5_low) {
    return -1;
  } else if (y_hit5_low < x_hit5_low) {
    return +1;
  } else if (x_hit5_high < y_hit5_high) {
    return -1;
  } else if (y_hit5_high < x_hit5_high) {
    return +1;

  } else if (x_hit3_low < y_hit3_low) {
    return -1;
  } else if (y_hit3_low < x_hit3_low) {
    return +1;
  } else if (x_hit3_high < y_hit3_high) {
    return -1;
  } else if (y_hit3_high < x_hit3_high) {
    return +1;

#if 0
  } else if (x->score_within_trims < y->score_within_trims) {
    return -1;
  } else if (y->score_within_trims < x->score_within_trims) {
    return +1;
#endif

  } else if (x->nmatches_plus_spliced_trims > y->nmatches_plus_spliced_trims) {
    return -1;
  } else if (y->nmatches_plus_spliced_trims > x->nmatches_plus_spliced_trims) {
    return +1;

#if 0
  } else if (x->nmatches_to_trims > y->nmatches_to_trims) {
    return -1;
  } else if (y->nmatches_to_trims > x->nmatches_to_trims) {
    return +1;
#endif

#if 0
    /* Causes hits to not be recognized as equivalent */
  } else if (x->nsplices < y->nsplices) {
    return -1;
  } else if (y->nsplices < x->nsplices) {
    return +1;
#endif

  } else if (x->alts_status_inside < y->alts_status_inside) {
    return -1;
  } else if (y->alts_status_inside < x->alts_status_inside) {
    return +1;

#if 0
  } else if (x->hit5->start_amb_length + x->hit5->end_amb_length +
	     x->hit3->start_amb_length + x->hit3->end_amb_length > 0 &&
	     y->hit5->start_amb_length + y->hit5->end_amb_length +
	     y->hit3->start_amb_length + y->hit3->end_amb_length == 0) {
    return -1;
  } else if (y->hit5->start_amb_length + y->hit5->end_amb_length +
	     y->hit3->start_amb_length + y->hit3->end_amb_length > 0 &&
	     x->hit5->start_amb_length + x->hit5->end_amb_length +
	     x->hit3->start_amb_length + x->hit3->end_amb_length == 0) {
    return +1;
#endif

  } else if (x->sense_consistent_p == true && y->sense_consistent_p == false) {
    return -1;
  } else if (x->sense_consistent_p == false && y->sense_consistent_p == true) {
    return +1;

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

#if 0
  } else if (x->sense_consistent_p == true) {
    /* Used for sorting, but not equiv */
    if ((x->hit5->sensedir_for_concordance != 0 || x->hit3->sensedir_for_concordance != 0) &&
	(y->hit5->sensedir_for_concordance == 0 && y->hit3->sensedir_for_concordance == 0)) {
      return -1;
    } else if ((y->hit5->sensedir_for_concordance != 0 || y->hit3->sensedir_for_concordance != 0) &&
	       (x->hit5->sensedir_for_concordance == 0 && x->hit3->sensedir_for_concordance == 0)) {
      return +1;
    } else {
      return 0;
    }
#endif

#if 0
  } else if (x->hit5->sensedir_for_concordance == y->hit5->sensedir_for_concordance &&
	     x->hit3->sensedir_for_concordance == y->hit3->sensedir_for_concordance) {
    return 0;
  } else if (x->hit5->sensedir_for_concordance > y->hit5->sensedir_for_concordance) {
    return +1;
  } else if (y->hit5->sensedir_for_concordance > x->hit5->sensedir_for_concordance) {
    return -1;
  } else if (x->hit3->sensedir_for_concordance > y->hit3->sensedir_for_concordance) {
    return +1;
  } else if (y->hit3->sensedir_for_concordance > x->hit3->sensedir_for_concordance) {
    return -1;
#endif

  } else {
    return 0;
  }
}


static int
hitpair_position_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;
  
  if (x->dir < y->dir) {
    return -1;
  } else if (y->dir < x->dir) {
    return -1;
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return +1;
  } else if (y->high < x->high) {
    return -1;
  } else {
    return 0;
  }
}


#if 0
static bool
hitpair_equal (Stage3pair_T x, Stage3pair_T y) {
  if (x->dir != y->dir) {
    debug8(printf("=>F "));
    return false;		/* Different strands */
  } else if (x->hit5->low != y->hit5->low) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit5->high != y->hit5->high) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit3->low != y->hit3->low) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit3->high != y->hit3->high) {
    debug8(printf("=>F "));
    return false;
  } else {
    debug8(printf("=>T "));
    return true;
  }
}
#endif


static bool
hitpair_overlap_p (Stage3pair_T x, Stage3pair_T y) {
  /* printf("Checking for overlap of %u..%u and %u..%u ",x->low,x->high,y->low,y->high); */
  if (x->dir != y->dir) {
    /* printf("=> false\n"); */
    return false;		/* Different strands */
  } else if (x->high < y->low) {
    /* printf("=> false\n"); */
    return false;
  } else if (x->low > y->high) {
    /* printf("=> false\n"); */
    return false;
  } else {
    /* printf("=> true\n"); */
    return true;
  }
}


static bool
hitpair_subsumption (Stage3pair_T x, Stage3pair_T y) {
  if (x->dir != y->dir) {
    return false;		/* Different strands */

  } else if (x->low <= y->low && x->high >= y->high) {
    return true;
  } else if (y->low <= x->low && y->high >= x->high) {
    return true;
    
    /* Test each end of the pair.  Example: 1586..1512 and 1400..1468 should subsume 1586..1512 and 1564..1617 */
  } else if (x->hit5->low <= y->hit5->low && x->hit5->high >= y->hit5->high) {
    return true;
  } else if (y->hit5->low <= x->hit5->low && y->hit5->high >= x->hit5->high) {
    return true;

  } else if (x->hit3->low <= y->hit3->low && x->hit3->high >= y->hit3->high) {
    return true;
  } else if (y->hit3->low <= x->hit3->low && y->hit3->high >= x->hit3->high) {
    return true;

  } else {
    return false;
  }
}


static int
pair_matches_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

  if (x->nmatches_plus_spliced_trims > y->nmatches_plus_spliced_trims) {
    return -1;
  } else if (x->nmatches_plus_spliced_trims < y->nmatches_plus_spliced_trims) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Stage3pair_sort_bymatches (List_T hits, Hitlistpool_T hitlistpool) {
  List_T sorted = NULL;
  Stage3pair_T *array;
  int n, i;

  
  if ((n = List_length(hits)) == 0) {
    return (List_T) NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    array = (Stage3pair_T *) MALLOCA(n * sizeof(Stage3pair_T));
    List_fill_array((void **) array,hits);
    Hitlist_free(&hits);
#else
    array = (Stage3pair_T *) List_to_array(hits,NULL);
    Hitlist_free(&hits);
#endif

    qsort(array,n,sizeof(Stage3pair_T),pair_matches_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = Hitlist_push(sorted,hitlistpool,(void *) array[i]);
    }
#ifdef USE_ALLOCA_FOR_HITS
    FREEA(array);
#else
    FREE(array);
#endif

    return sorted;
  }
}



#if 0
List_T
Stage3pair_remove_duplicates_exact (List_T hitpairlist) {
  List_T unique = NULL;
  Stage3pair_T hitpair, *hitpairs;
  int n, i, j;
  bool *eliminate;

  debug8(printf("Entered Stage3pair_remove_duplicates_exact with %d pairs\n",n));
  if ((n = List_length(hitpairlist)) == 0) {
    return NULL;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) MALLOCA(n * sizeof(Stage3pair_T));
    List_fill_array((void **) hitpairs,hitpairlist);
    Hitlist_free(&hitpairlist);
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    Hitlist_free(&hitpairlist);
#endif
  }

  debug8(printf("Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), circularalias %d|%d, nmatches: %d (%d to_trims)\n",
		  i,Pairtype_string(hitpair->pairtype),Method_string(hitpair->hit5->method),
		  Method_string(hitpair->hit3->method),hitpair,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->hit5->circularalias,hitpair->hit3->circularalias,
		  hitpair->nmatches_plus_spliced_trims,hitpair->nmatches_to_trims);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hitpair_equiv_cmp(hitpairs[j],hitpairs[i]) == 0) {
      debug8(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
  }

  for (i = n-1; i >= 0; i--) {
    hitpair = hitpairs[i];
    if (eliminate[i] == false) {
      unique = Hitlist_push(unique,hitlistpool,(void *) hitpair);
    } else {
      Stage3pair_free(&hitpair);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hitpairs);
  FREEA(eliminate);
#else
  FREE(hitpairs);
  FREE(eliminate);
#endif

  debug8(printf("Exited Stage3pair_remove_duplicates_exact with %d pairs\n",List_length(unique)));
  return unique;
}
#endif


static int
hitpair_goodness_cmp (bool *equalp, Stage3pair_T hitpair,
		      Stage3pair_T best_hitpair, bool finalp) {
  double prob1, prob2;
  /* Chrpos_T total_querylength, best_total_querylength; */
  double zscore, best_zscore;

#if 0
  int hitpair_nmatches, best_hitpair_nmatches;
  int max_trim_querystart, max_trim_queryend;
  Stage3end_T hit5, besthit5, hit3, besthit3;

  if (hitpair->absdifflength_bingo_p < best_hitpair->absdifflength_bingo_p) {
    /* k is worse */
    debug8(printf(" => loses by absdifflength (bingo)\n"));
    return -1;
  } else if (hitpair->absdifflength_bingo_p > best_hitpair->absdifflength_bingo_p) {
    /* k is better */
    debug8(printf(" => wins by absdifflength (bingo)\n"));
    return +1;
  }
#endif

#ifdef PRE_RESOLVE_MULTIMAPPING
  if (TALLY_RATIO*Stage3pair_tally(hitpair) < Stage3pair_tally(best_hitpair)) {
    /* k is worse */
    debug8(printf(" => loses by tally\n"));
    return -1;
  } else if (Stage3pair_tally(hitpair) > TALLY_RATIO*Stage3pair_tally(best_hitpair)) {
    /* k is better */
    debug8(printf(" => wins by tally\n"));
    return +1;
  }
#endif

  *equalp = false;


  /* Previously, we favored ambiguous splices over definitive ones, but
     now that we are generating Stage3end_T objects with and without the
     end exons, we prefer definitive splices */
  if (known_ambiguous_p(hitpair->hit5) == true && known_ambiguous_p(best_hitpair->hit5) == false &&
      known_ambiguous_p(hitpair->hit3) == known_ambiguous_p(best_hitpair->hit3) &&
      hitpair->insertlength <= best_hitpair->insertlength) {
    debug8(printf("Case 1\n"));
    return -1;

  } else if (known_ambiguous_p(hitpair->hit5) == false && known_ambiguous_p(best_hitpair->hit5) == true &&
	     known_ambiguous_p(hitpair->hit3) == known_ambiguous_p(best_hitpair->hit3) &&
	     hitpair->insertlength >= best_hitpair->insertlength) {
    debug8(printf("Case 2\n"));
    return +1;

  } else if (known_ambiguous_p(hitpair->hit3) == true && known_ambiguous_p(best_hitpair->hit3) == false &&
	     known_ambiguous_p(hitpair->hit5) == known_ambiguous_p(best_hitpair->hit5) &&
	     hitpair->insertlength <= best_hitpair->insertlength) {
    debug8(printf("Case 3\n"));
    return -1;

  } else if (known_ambiguous_p(hitpair->hit3) == false && known_ambiguous_p(best_hitpair->hit3) == true &&
	     known_ambiguous_p(hitpair->hit5) == known_ambiguous_p(best_hitpair->hit5) &&
	     hitpair->insertlength > best_hitpair->insertlength) {
    debug8(printf("Case 4\n"));
    return +1;

  } else if (hitpair->hit5->nsegments + hitpair->hit3->nsegments > best_hitpair->hit5->nsegments + best_hitpair->hit3->nsegments) {
    if (hitpair->nmatches_plus_spliced_trims > best_hitpair->nmatches_plus_spliced_trims + NMATCHES_SLOP) {
      /* More segments and strictly more matches */
      debug8(printf("More segments and strictly more matches (to_trims)\n"));
      return +1;
    } else if (hitpair->nmatches_plus_spliced_trims < best_hitpair->nmatches_plus_spliced_trims - NMATCHES_SLOP) {
      /* More segments and don't add anything */
      debug8(printf("More segments but don't add anything\n"));
      return -1;
    } else {
      /* More segments, but same nmatches */
      debug8(printf("More segments but same nmatches\n"));
      return 0;
    }

  } else if (hitpair->hit5->nsegments + hitpair->hit3->nsegments < best_hitpair->hit5->nsegments + best_hitpair->hit3->nsegments) {
    if (hitpair->nmatches_plus_spliced_trims > best_hitpair->nmatches_plus_spliced_trims + NMATCHES_SLOP) {
      /* Fewer segments, but same or more matches */
      debug8(printf("Fewer segments and same or more matches (to_trims)\n"));
      return +1;
    } else if (hitpair->nmatches_plus_spliced_trims < best_hitpair->nmatches_plus_spliced_trims - NMATCHES_SLOP) {
      /* Fewer segments and don't add anything */
      debug8(printf("Fewer segments and don't add anything\n"));
      return -1;
    } else {
      /* Fewer segments but same nmatches */
      debug8(printf("Fewer segments but same nmatches\n"));
      return 0;
    }

#if 0
  } else if ((hitpair->hit5->hittype != TRANSCRIPTOME || hitpair->hit3->hittype != TRANSCRIPTOME) &&
	     (best_hitpair->hit5->hittype == TRANSCRIPTOME || best_hitpair->hit3->hittype == TRANSCRIPTOME)) {
    /* k is worse */
    debug8(printf(" => loses by transcriptome\n"));
    return -1;

  } else if ((hitpair->hit5->hittype == TRANSCRIPTOME || hitpair->hit3->hittype == TRANSCRIPTOME) &&
	     (best_hitpair->hit5->hittype != TRANSCRIPTOME || best_hitpair->hit3->hittype != TRANSCRIPTOME)) {
    /* k is better */
    debug8(printf(" => wins by transcriptome\n"));
    return +1;
#endif

  } else if (hitpair->nmatches_plus_spliced_trims < best_hitpair->nmatches_plus_spliced_trims - NMATCHES_SLOP) {
    /* k is worse */
    debug8(printf(" => loses by nmatches\n"));
    return -1;
  } else if (hitpair->nmatches_plus_spliced_trims > best_hitpair->nmatches_plus_spliced_trims + NMATCHES_SLOP) {
    /* k is better */
    debug8(printf(" => wins by nmatches\n"));
    return +1;

#if 0
  } else if (hitpair->nsplices > best_hitpair->nsplices) {
    /* k is worse */
    debug8(printf(" => loses by nsplices: %d > %d in best\n",hitpair->nsplices,best_hitpair->nsplices));
    return -1;
  } else if (hitpair->nsplices < best_hitpair->nsplices) {
    /* k is better */
    debug8(printf(" => wins by nsplices: %d < %d in best\n",hitpair->nsplices,best_hitpair->nsplices));
    return +1;
#endif

  } else if (hitpair->alts_status_inside > best_hitpair->alts_status_inside) {
    /* k is worse */
    debug8(printf(" => loses by alts_status_inside\n"));
    return -1;
  } else if (hitpair->alts_status_inside < best_hitpair->alts_status_inside) {
    /* k is better */
    debug8(printf(" => wins by alts_status_inside\n"));
    return +1;


  } else if (hitpair->hit5->hittype > best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype >= best_hitpair->hit3->hittype) {
    /* k is worse */
    debug8(printf(" => loses by hittype\n"));
    return -1;

  } else if (hitpair->hit5->hittype >= best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype > best_hitpair->hit3->hittype) {
    /* k is worse */
    debug8(printf(" => loses by hittype\n"));
    return -1;

  } else if (hitpair->hit5->hittype < best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype <= best_hitpair->hit3->hittype) {
    /* k is better */
    debug8(printf(" => wins by hittype\n"));
    return +1;

  } else if (hitpair->hit5->hittype <= best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype < best_hitpair->hit3->hittype) {
    /* k is better */
    debug8(printf(" => wins by hittype\n"));
    return +1;

#if 0
  } else if ((hitpair->alts_resolve_5 == -1 || hitpair->alts_resolve_3 == -1) &&
	     (best_hitpair->alts_resolve_5 != -1 && best_hitpair->alts_resolve_3 != -1)) {
    /* best_hitpair is resolved, hitpair is ambiguous.  best_hitpair wins */
    debug8(printf(" => loses by resolve_inside\n"));
    return -1;

  } else if ((hitpair->alts_resolve_5 != -1 && hitpair->alts_resolve_3 != -1) &&
	     (best_hitpair->alts_resolve_5 == -1 || best_hitpair->alts_resolve_3 == -1)) {
    /* hitpair is resolved, best_hitpair is ambiguous.  hitpair wins */
    debug8(printf(" => wins by resolve_inside: %d, %d, %d, %d\n",
		  hitpair->alts_resolve_5,hitpair->alts_resolve_3,
		  best_hitpair->alts_resolve_5,best_hitpair->alts_resolve_3));
    return +1;
#endif

#if 0
  } else if (n_amb_ends(hitpair->hit5) + n_amb_ends(hitpair->hit3) >
	     n_amb_ends(best_hitpair->hit5) + n_amb_ends(best_hitpair->hit3)) {
    /* k is worse */
    debug8(printf(" => loses by ambiguity\n"));
    return -1;

  } else if (n_amb_ends(hitpair->hit5) + n_amb_ends(hitpair->hit3) <
	     n_amb_ends(best_hitpair->hit5) + n_amb_ends(best_hitpair->hit3)) {
    /* k is better */
    debug8(printf(" => wins by ambiguity\n"));
    return +1;
#endif

  } else if (hitpair->hit5->splice_score + hitpair->hit3->splice_score >
	     best_hitpair->hit5->splice_score + best_hitpair->hit3->splice_score) {
    /* k is worse */
    debug8(printf(" => loses by splice score\n"));
    return -1;

  } else if (hitpair->hit5->splice_score + hitpair->hit3->splice_score >
	     best_hitpair->hit5->splice_score + best_hitpair->hit3->splice_score) {
    /* k is better */
    debug8(printf(" => wins by splice score\n"));
    return +1;

#if 0
  } else if (hitpair->absdifflength < best_hitpair->absdifflength) {
    /* k is worse */
    debug8(printf(" => loses by absdifflength\n"));
    return -1;
  } else if (hitpair->absdifflength > best_hitpair->absdifflength) {
    /* k is better */
    debug8(printf(" => wins by absdifflength\n"));
    return +1;
#endif

  } else if (finalp == false) {
    debug8(printf("  => indistinguishable\n"));
    return 0;

#ifdef USE_ABSDIFFLENGTH
    /* If insert length is within deviation of expected pairlength, favor it */
  } else if (best_hitpair->absdifflength <= (Chrpos_T) pairlength_deviation &&
	     hitpair->absdifflength > (Chrpos_T) pairlength_deviation) {
    /* k is worse */
    debug8(printf(" => loses by absdifflength within deviation %d\n",pairlength_deviation));
    return -1;
  } else if (hitpair->absdifflength <= (Chrpos_T) pairlength_deviation &&
	     best_hitpair->absdifflength > (Chrpos_T) pairlength_deviation) {
    /* k is better */
    debug8(printf(" => wins by absdifflength within deviation %d\n",pairlength_deviation));
    return +1;
#endif

#if 0
    /* Previously favored longer insert lengths to give more compact
       splices.  However, we now accept splices first that give
       expected pairlength */
  } else if (hitpair->insertlength_expected_sign == -1 && best_hitpair->insertlength_expected_sign == +1) {
    /* k is worse */
    debug8(printf(" => loses by insertlength_expected_sign\n"));
    return -1;
  } else if (hitpair->insertlength_expected_sign == +1 && best_hitpair->insertlength_expected_sign == -1) {
    /* k is better */
    debug8(printf(" => wins by insertlength_expected_sign\n"));
    return +1;
#endif

    /* Next we look at splice probability */
  } else {
    debug8(printf(" => prob"));
    prob1 = Stage3end_prob(hitpair->hit5) + Stage3end_prob(hitpair->hit3);
    prob2 = Stage3end_prob(best_hitpair->hit5) + Stage3end_prob(best_hitpair->hit3);
    if (prob1 + 0.3 < prob2) {
      /* k is worse */
      debug8(printf(" => loses by dual splice prob %f vs %f\n",prob1,prob2));
      return -1;
    } else if (prob1 > prob2 + 0.3) {
      /* k is better */
      debug8(printf(" => wins by dual splice prob %f vs %f\n",prob1,prob2));
      return +1;
    } else {
      debug8(printf(" => neither wins\n"));
    }


#if 0
    /* Overlapping ends worse than separate ends */
    total_querylength = (Chrpos_T) (hitpair->hit5->querylength + hitpair->hit3->querylength);
    best_total_querylength = (Chrpos_T) (best_hitpair->hit5->querylength + best_hitpair->hit3->querylength);

    if (hitpair->insertlength <= total_querylength && best_hitpair->insertlength > best_total_querylength) {
      debug8(printf(" => loses by being overlapping\n"));
      return -1;
    } else if (hitpair->insertlength > total_querylength && best_hitpair->insertlength <= best_total_querylength) {
      debug8(printf(" => wins by being separate\n"));
      return +1;

      /* Next, favor shorter outerlengths to give more compact splices or closer pairs */
    } else if (hitpair->outerlength > best_hitpair->outerlength + OUTERLENGTH_SLOP) {
      /* k is worse */
      debug8(printf(" => loses by outerlength\n"));
      return -1;
    } else if (hitpair->outerlength + OUTERLENGTH_SLOP < best_hitpair->outerlength) {
      /* k is better */
      debug8(printf(" => wins by outerlength\n"));
      return +1;
      
    } else {
#if 0
      if (hitpair->insertlength_expected_sign >= 0 && best_hitpair->insertlength_expected_sign >= 0) {
	/* Both insert lengths are short, so favor shorter insert length */
	debug8(printf(" => short insertlengths"));
	/* Favor shorter insert lengths */
	if (hitpair->insertlength > best_hitpair->insertlength) {
	  /* k is worse */
	  debug8(printf(" => loses by insertlength\n"));
	  return -1;
	} else if (hitpair->insertlength < best_hitpair->insertlength) {
	  /* k is better */
	  debug8(printf(" => wins by insertlength\n"));
	  return +1;
	}
      }
#endif

      /* Both insert lengths are long, so favor longer insert length to give more compact splices */
      debug8(printf(" => long insertlengths"));
      if (hitpair->insertlength < best_hitpair->insertlength) {
	/* k is worse */
	debug8(printf(" => loses by insertlength\n"));
	return -1;
      } else if (hitpair->insertlength > best_hitpair->insertlength) {
	/* k is better */
	debug8(printf(" => wins by insertlength\n"));
	return +1;
      }

      debug8(printf("  => equal\n"));
      *equalp = true;
      return 0;
    }
#endif

    /* Look at expected pairlength and pairlength deviation */
    if (hitpair->insertlength < expected_pairlength) {
      zscore = (double) (expected_pairlength - (Chrpos_T) hitpair->insertlength) / (double) pairlength_deviation;
    } else {
      zscore = (double) ((Chrpos_T) hitpair->insertlength - expected_pairlength) / (double) pairlength_deviation;
    }
    if (best_hitpair->insertlength < expected_pairlength) {
      best_zscore = (double) (expected_pairlength - (Chrpos_T) best_hitpair->insertlength) / (double) pairlength_deviation;
    } else {
      best_zscore = (double) ((Chrpos_T) best_hitpair->insertlength - expected_pairlength) / (double) pairlength_deviation;
    }
    debug8(printf("expected_pairlength %u, pairlength_deviation %u\n",expected_pairlength,pairlength_deviation));
    debug8(printf("Comparing insertlength %d (z score %f) with best_insertlength %d (zscore %f)\n",
		  hitpair->insertlength,zscore,best_hitpair->insertlength,best_zscore));

    if (zscore > best_zscore + 1.0) {
      /* k is worse */
      debug8(printf(" => loses by insertlength and zscore\n"));
      return -1;
    } else if (best_zscore > zscore + 1.0) {
      /* k is better */
      debug8(printf(" => wins by insertlength and zscore\n"));
      return +1;
    }
    
    debug8(printf("  => equal\n"));
    *equalp = true;
    return 0;
  }
}


#if 0
static bool
hitpair_bad_superstretch_p (Stage3pair_T hitpair_k, Stage3pair_T *hitpairs, int k, int j,
			    bool finalp) {
  int a;
  bool equalp;

  for (a = k+1; a <= j; a++) {
    if (hitpair_subsumption(hitpair_k,hitpairs[a]) == true) {
      debug8(printf("Testing %d because stretches over %d",k,a));
      if (hitpair_goodness_cmp(&equalp,hitpairs[a],
			       hitpair_k,finalp) > 0 || equalp == true) {
	debug8(printf(" => eliminating\n"));
	return true;
      }
      debug8(printf("\n"));
    }
  }
  return false;
}
#endif


/* Recursive, list-based approach */
static List_T
pair_remove_bad_superstretches (bool *keep_p, Stage3pair_T superstretch, List_T list,
				Hitlistpool_T hitlistpool, bool finalp) {
  List_T result = NULL, better, equal, p, q, r;
  Stage3pair_T stage3pair, hitpair;
  bool equalp, this_kept_p;
  int cmp;

  *keep_p = true;

  p = list;
  while (p != NULL) {
    stage3pair = (Stage3pair_T) List_head(p);

    q = List_next(p);
    while (q != NULL && hitpair_subsumption(stage3pair,(Stage3pair_T) List_head(q)) == true) {
#ifdef DEBUG8
      printf("  This (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), nmatches: %d (%d to_trims), insertlength %d, alts_status_inside %d, amb_lengths %d and %d\n",
	     Pairtype_string(stage3pair->pairtype),Method_string(stage3pair->hit5->method),
	     Method_string(stage3pair->hit3->method),stage3pair,
	     stage3pair->hit5->low - stage3pair->hit5->chroffset,stage3pair->hit5->high - stage3pair->hit5->chroffset,
	     stage3pair->hit3->low - stage3pair->hit3->chroffset,stage3pair->hit3->high - stage3pair->hit3->chroffset,
	     stage3pair->dir,stage3pair->nmatches_plus_spliced_trims,stage3pair->nmatches_to_trims,
	     stage3pair->insertlength,stage3pair->alts_status_inside,amb_length(stage3pair->hit5),amb_length(stage3pair->hit3));

      hitpair = (Stage3pair_T) List_head(q);
      printf("subsumes that (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), nmatches: %d (%d to_trims), insertlength %d, alts_status_inside %d, amb_lengths %d and %d\n",
	     Pairtype_string(hitpair->pairtype),Method_string(hitpair->hit5->method),
	     Method_string(hitpair->hit3->method),hitpair,
	     hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
	     hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
	     hitpair->dir,hitpair->nmatches_plus_spliced_trims,hitpair->nmatches_to_trims,
	     hitpair->insertlength,hitpair->alts_status_inside,amb_length(hitpair->hit5),amb_length(hitpair->hit3));
#endif
      q = List_next(q);
    }

    if (q == p) {
      result = Hitlist_push(result,hitlistpool,(void *) stage3pair);
      if (superstretch != NULL && 
	  (hitpair_goodness_cmp(&equalp,stage3pair,superstretch,finalp) > 0 || equalp == true)) {
	*keep_p = false;
      }
      p = List_next(q);

    } else {
      /* Cluster */
      debug8(printf("Processing cluster of size %d - %d\n",List_length(p),List_length(q)));
      better = equal = (List_T) NULL;
      for (r = List_next(p); r != q; r = List_next(r)) {
	debug8(printf("Calling hitpair_goodness_cmp\n"));
	cmp = hitpair_goodness_cmp(&equalp,(Stage3pair_T) List_head(r),stage3pair,finalp);
	debug8(printf("cmp = %d, equalp = %d\n",cmp,equalp));
	if (cmp > 0) {
	  better = Hitlist_push(better,hitlistpool,(void *) List_head(r));
	} else if (cmp < 0) {
	  hitpair = (Stage3pair_T) List_head(r);
	  Stage3pair_free(&hitpair);
	} else {
	  equal = Hitlist_push(equal,hitlistpool,(void *) List_head(r));
	}
      }

      debug8(printf("Found %d better, %d equal\n",List_length(better),List_length(equal)));

      if (better == NULL) {
	/* All children are equal to parent */
	debug8(printf("All children are equivalent, so keeping parent and all (equal) children\n"));
	result = Hitlist_push(result,hitlistpool,(void *) stage3pair);
	equal = List_reverse(equal); /* Keep original order */
	for (r = equal; r != NULL; r = List_next(r)) {
	  hitpair = (Stage3pair_T) List_head(r);
	  result = Hitlist_push(result,hitlistpool,(void *) hitpair);
	}
	Hitlist_free(&equal);

      } else {
	/* Exists a child better than parent */
	debug8(printf("Exists a child better than parent, so deleting parent and equal and calling recursively among all (better) children\n"));
	Stage3pair_free(&stage3pair);
	for (r = equal; r != NULL; r = List_next(r)) {
	  hitpair = (Stage3pair_T) List_head(r);
	  Stage3pair_free(&hitpair);
	}
	Hitlist_free(&equal);

	if (List_length(better) == 1) {
	  hitpair = (Stage3pair_T) List_head(better);
	  result = Hitlist_push(result,hitlistpool,(void *) hitpair);
	  Hitlist_free(&better);

	} else {
	  /* Don't call List_reverse(better) */
	  result = List_append(result,pair_remove_bad_superstretches(&this_kept_p,/*superstretch*/NULL,better,
								     hitlistpool,finalp));
#if 0
	  /* Already deleted parent */
	  if (this_kept_p == false) {
	    Stage3pair_free(&stage3pair);
	  } else {
	    /* Compare stage3pair against the current parent */
	    result = Hitlist_push(result,hitlistpool,(void *) stage3pair);
	    if (superstretch != NULL && 
		(hitpair_goodness_cmp(&equalp,stage3pair,superstretch,finalp) >= 0 || equalp == true)) {
	      *keep_p = false;
	    }
	  }
#endif
	}
      }

      p = q;
    }
  }

  Hitlist_free(&list);

  debug8(printf("Returning result of length %d\n",List_length(result)));
  return List_reverse(result);
}


static List_T
pair_remove_overlaps (List_T hitpairlist, Hitlistpool_T hitlistpool,
		      bool translocp, bool finalp) {
  List_T unique = NULL;
  Stage3pair_T hitpair, parent, *hitpairs;
  int nkept, n, i, j;
  bool *eliminate;
  int *parenti;
  bool keep_p;

  n = List_length(hitpairlist);
  debug8(printf("  Entering pair_remove_overlaps with %d pairs: %s\n",
		n,finalp == true ? "FINAL" : "not final"));

  if (n < 2) {
    debug8(printf("  Exiting pair_remove_overlaps with %d < 2 pairs\n",n));
    return hitpairlist;
  } else {
#ifdef USE_ALLOCA_FOR_HITS
    eliminate = (bool *) CALLOCA(n,sizeof(bool));
    parenti = (int *) CALLOCA(n,sizeof(int));
    hitpairs = (Stage3pair_T *) MALLOCA(n * sizeof(Stage3pair_T));
    List_fill_array((void **) hitpairs,hitpairlist);
    Hitlist_free(&hitpairlist);
#else
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    parenti = (int *) CALLOC(n,sizeof(int));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    Hitlist_free(&hitpairlist);
#endif
  }

  /* Step 1.  Check for exact duplicates */
  debug8(printf("  Step 1.  Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %p, %u..%u|%u..%u (dir = %d), circularalias %d|%d, nmatches: %d (%d to_trims), amb_lengths %d and %d, sensedirs %d and %d.",
		  i,Pairtype_string(hitpair->pairtype),Method_string(hitpair->hit5->method),
		  Method_string(hitpair->hit3->method),hitpair,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->hit5->circularalias,hitpair->hit3->circularalias,hitpair->nmatches_plus_spliced_trims,hitpair->nmatches_to_trims,
		  amb_length(hitpair->hit5),amb_length(hitpair->hit3),hitpair->hit5->sensedir,hitpair->hit3->sensedir);
	   if (hitpair->hit5->hittype == TRANSLOC_SPLICE) {
	     printf("  5' TRANSLOC splice probs %f",hitpair->hit5->splice_score);
	   }
	   if (hitpair->hit3->hittype == TRANSLOC_SPLICE) {
	     printf("  3' TRANSLOC splice probs %f",hitpair->hit3->splice_score);
	   }
	   printf("\n");
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    debug8(printf(" %d,%d",i,j));
    while (j < n && hitpair_equiv_cmp(hitpairs[j],hitpairs[i]) == 0) {
      debug8(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      parenti[j] = i;
      j++;
    }
    i = j;
  }
  debug8(printf("\n"));

  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    }
  }
  debug8(printf("nkept = %d\n",nkept));

  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    debug8(printf("All entries eliminate one another, so keep the first one\n"));
    eliminate[0] = false;
    nkept = 1;
  }

  for (i = n - 1; i >= 0; --i) {
    hitpair = hitpairs[i];
    if (eliminate[i] == false) {
      debug8(printf("  Keeping %s|%s %u..%u|%u..%u, nmatches (trimmed) %d, score %d, (dir = %d)\n",
		    Method_string(hitpair->hit5->method),Method_string(hitpair->hit3->method),
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches_plus_spliced_trims,hitpair->score_overall,hitpair->dir));
      unique = Hitlist_push(unique,hitlistpool,(void *) hitpair);

    } else {
      debug8(printf("  Eliminating %s|%s %u..%u|%u..%u, nmatches (trimmed) %d, score %d, (dir = %d)\n",
		    Method_string(hitpair->hit5->method),Method_string(hitpair->hit3->method),
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches_plus_spliced_trims,hitpair->score_overall,hitpair->dir));

      parent = hitpairs[parenti[i]];
      Stage3pair_transfer_transcripts_one(parent,hitpair);
      Stage3pair_free(&hitpair);
    }
  }

#ifdef USE_ALLOCA_FOR_HITS
  FREEA(hitpairs);
  FREEA(eliminate);
  FREEA(parenti);
#else
  FREE(hitpairs);
  FREE(eliminate);
  FREE(parenti);
#endif

  debug8(printf("  Step 2.  Checking for bad superstretches\n"));
  if (0 && translocp == true) {
    return unique;
  } else {
    return pair_remove_bad_superstretches(&keep_p,/*superstretch*/NULL,unique,
					  hitlistpool,finalp);
  }
}


List_T
Stage3pair_remove_overlaps (List_T hitpairlist, Hitlistpool_T hitlistpool,
			    bool translocp, bool finalp) {
  List_T unique_separate, unique_overlapping,
    separate = NULL, overlapping = NULL, p;
  Stage3pair_T hitpair_separate, hitpair_overlapping, hitpair;

  Stage3pair_T *array_separate, *array_overlapping;
  Univcoord_T low, high;
  bool subsumedp, equalp;
  int n_separate, n_overlapping, i, j;


  debug8(printf("Entered Stage3pair_remove_overlaps with %d hitpairs\n",List_length(hitpairlist)));
  for (p = hitpairlist; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    if (hitpair->insertlength <= (Chrpos_T) (hitpair->hit5->querylength + hitpair->hit3->querylength)) {
      overlapping = Hitlist_push(overlapping,hitlistpool,(void *) hitpair);
    } else {
      separate = Hitlist_push(separate,hitlistpool,(void *) hitpair);
    }
  }
  Hitlist_free(&hitpairlist);

  debug8(printf("Calling Stage3pair_remove_overlaps for separate pair ends\n"));
  unique_separate = pair_remove_overlaps(separate,hitlistpool,translocp,finalp);

  debug8(printf("Calling Stage3pair_remove_overlaps for overlapping pair ends\n"));
  unique_overlapping = pair_remove_overlaps(overlapping,hitlistpool,translocp,finalp);
  
  if (unique_overlapping == NULL) {
    debug8(printf("Unique overlapping is NULL\n"));
    return unique_separate;
  } else if (unique_separate == NULL) {
    debug8(printf("Unique separate is NULL\n"));
    return unique_overlapping;
  } else {
    debug8(printf("Have both overlapping and separate\n"));
    n_overlapping = List_length(unique_overlapping);
#ifdef USE_ALLOCA_FOR_HITS
    array_overlapping = (Stage3pair_T *) MALLOCA(n_overlapping * sizeof(Stage3pair_T));
    List_fill_array((void **) array_overlapping,unique_overlapping);
#else
    array_overlapping = (Stage3pair_T *) List_to_array(unique_overlapping,NULL);
#endif

    n_separate = List_length(unique_separate);
#ifdef USE_ALLOCA_FOR_HITS
    array_separate = (Stage3pair_T *) MALLOCA(n_separate * sizeof(Stage3pair_T));
    List_fill_array((void **) array_separate,unique_separate);
#else
    array_separate = (Stage3pair_T *) List_to_array(unique_separate,NULL);
#endif

    qsort(array_overlapping,n_overlapping,sizeof(Stage3pair_T),hitpair_position_cmp);
    qsort(array_separate,n_separate,sizeof(Stage3pair_T),hitpair_position_cmp);

    /* 1.  First, favor overlapping (with smaller insertlengths) */
    /* Keep unique_overlapping and filter unique_separate into indep_separate */
    Hitlist_free(&unique_separate);
    unique_separate = (List_T) NULL;

    i = j = 0;
    for (i = 0; i < n_separate; i++) {
      hitpair_separate = array_separate[i];
      low = hitpair_separate->low;
      high = hitpair_separate->high;
      while (j >= 0 && array_overlapping[j]->high >= low) {
	j--;
      }
      j += 1;

      subsumedp = false;
      while (j < n_overlapping && subsumedp == false && array_overlapping[j]->low <= high) {
	if (hitpair_goodness_cmp(&equalp,array_overlapping[j],
				 hitpair_separate,finalp) > 0) {
	  debug8(printf("overlapping pair %d better than separate pair %d\n",j,i));
	  subsumedp = hitpair_subsumption(hitpair_separate,array_overlapping[j]);
	  debug8(printf("  checking if separate pair %d subsumes overlapping pair %d => %d\n",
			i,j,subsumedp));
	}
	j++;
      }
      j -= 1;

      if (subsumedp == true) {
	Stage3pair_free(&hitpair_separate);
      } else {
        unique_separate = Hitlist_push(unique_separate,hitlistpool,(void *) hitpair_separate);
      }
    }

#ifdef USE_ALLOCA_FOR_HITS
    FREEA(array_separate);
#else
    FREE(array_separate);
#endif

    if ((n_separate = List_length(unique_separate)) == 0) {
#ifdef USE_ALLOCA_FOR_HITS
      FREEA(array_overlapping);
#else
      FREE(array_overlapping);
#endif
      return unique_overlapping;

    } else {
#ifdef USE_ALLOCA_FOR_HITS
      array_separate = (Stage3pair_T *) MALLOCA(n_separate * sizeof(Stage3pair_T));
      List_fill_array((void **) array_separate,unique_separate);
#else
      array_separate = (Stage3pair_T *) List_to_array(unique_separate,NULL);
#endif

      /* 2.  Second, favor separate (with larger insertlengths) */
      /* Keep indep_separate and filter unique_overlapping into indep_overlapping */
      Hitlist_free(&unique_overlapping);
      unique_overlapping = (List_T) NULL;

      i = j = 0;
      for (i = 0; i < n_overlapping; i++) {
	hitpair_overlapping = array_overlapping[i];
	low = hitpair_overlapping->low;
	high = hitpair_overlapping->high;
	while (j >= 0 && array_separate[j]->high >= low) {
	  j--;
	}
	j += 1;

	subsumedp = false;
	while (j < n_separate && subsumedp == false && array_separate[j]->low <= high) {
	  if (hitpair_goodness_cmp(&equalp,array_separate[j],
				   hitpair_overlapping,finalp) > 0) {
	    debug8(printf("separate pair %d better than overlapping pair %d\n",j,i));
	    subsumedp = hitpair_subsumption(array_separate[j],hitpair_overlapping);
	    debug8(printf("  checking if separate pair %d subsumes overlapping pair %d => %d\n",
			  j,i,subsumedp));
	  }
	  j++;
	}
	j -= 1;
	
	if (subsumedp == true) {
	  Stage3pair_free(&hitpair_overlapping);
	} else {
	  unique_overlapping = Hitlist_push(unique_overlapping,hitlistpool,(void *) hitpair_overlapping);
	}
      }
    }

#ifdef USE_ALLOCA_FOR_HITS
    FREEA(array_separate);
    FREEA(array_overlapping);
#else
    FREE(array_separate);
    FREE(array_overlapping);
#endif

    return List_append(unique_overlapping,unique_separate);
  }
}


List_T
Stage3pair_resolve_multimapping (List_T hitpairs, Hitlistpool_T hitlistpool) {
  List_T resolve1, resolve2, resolve3, p;
  Stage3pair_T hitpair;

  long int best_tally;
  double tally_threshold;
  bool runlengthp;


  if (List_length(hitpairs) <= 1) {
    return hitpairs;
  }

#if 0
  if (genes_iit == NULL) {
    resolve1 = hitpairs;
  } else {
    best_overlap = NO_KNOWN_GENE;
    for (p = hitpairs; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((hitpair->gene_overlap = Stage3pair_gene_overlap(hitpair)) > best_overlap) {
	best_overlap = hitpair->gene_overlap;
      }
    }
    if (best_overlap == NO_KNOWN_GENE) {
      resolve1 = hitpairs;
    } else {
      resolve1 = (List_T) NULL;
      for (p = hitpairs; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if (hitpair->gene_overlap < best_overlap) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve1 = Hitlist_push(resolve1,hitlistpool,(void *) hitpair);
	}
      }
      Hitlist_free(&hitpairs);
    }
  }
      
  if (List_length(resolve1) <= 1) {
    return resolve1;
  }
#else
  resolve1 = hitpairs;
#endif

  if (tally_iit == NULL) {
    resolve2 = resolve1;
  } else {
    best_tally = 0L;
    for (p = resolve1; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((hitpair->tally = Stage3end_compute_tally(hitpair->hit5) + Stage3end_compute_tally(hitpair->hit3)) > best_tally) {
	best_tally = hitpair->tally;
      }
    }
    if (best_tally == 0L) {
      resolve2 = resolve1;
    } else {
      resolve2 = (List_T) NULL;
#ifdef USE_TALLY_RATIO
      tally_threshold = (double) best_tally / TALLY_RATIO;
#else
      tally_threshold = 1.0;
#endif
      for (p = resolve1; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if ((double) hitpair->tally < tally_threshold) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve2 = Hitlist_push(resolve2,hitlistpool,(void *) hitpair);
	}
      }
      Hitlist_free(&resolve1);
    }
  }

  if (List_length(resolve2) <= 1) {
    return resolve2;
  }

  if (runlength_iit == NULL) {
    resolve3 = resolve2;
  } else {
    runlengthp = false;
    for (p = resolve2; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (Stage3end_runlength_p(hitpair->hit5) == true || Stage3end_runlength_p(hitpair->hit3) == true) {
	runlengthp = true;
      }
    }
    if (runlengthp == false) {
      resolve3 = resolve2;
    } else {
      resolve3 = (List_T) NULL;
      for (p = resolve2; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if (Stage3end_runlength_p(hitpair->hit5) == false && Stage3end_runlength_p(hitpair->hit3) == false) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve3 = Hitlist_push(resolve3,hitlistpool,(void *) hitpair);
	}
      }
      Hitlist_free(&resolve2);
    }
  }


  return resolve3;
}


List_T
Stage3pair_filter (List_T hits, Hitlistpool_T hitlistpool,
		   int max_mismatches_5, int max_mismatches_3,
		   int min_coverage_5, int min_coverage_3) {
  List_T newhits = NULL, p;
  Stage3end_T hit5, hit3;
  Stage3pair_T hitpair;

  if (ignore_trim_p == false) {
    for (p = hits; p != NULL; p = List_next(p)) {
      hitpair = (Stage3pair_T) List_head(p);
      hit5 = hitpair->hit5;
      hit3 = hitpair->hit3;
      if (hit5->score_overall > max_mismatches_5 || hit3->score_overall > max_mismatches_3) {
	Stage3pair_free(&hitpair);
      } else if (hit5->querylength - hit5->trim_querystart - hit5->trim_queryend < min_coverage_5 &&
		 hit3->querylength - hit3->trim_querystart - hit3->trim_queryend < min_coverage_3) {
	Stage3pair_free(&hitpair);
      } else {
	newhits = Hitlist_push(newhits,hitlistpool,(void *) hitpair);
      }
    }
  } else {
    for (p = hits; p != NULL; p = List_next(p)) {
      hitpair = (Stage3pair_T) List_head(p);
      hit5 = hitpair->hit5;
      hit3 = hitpair->hit3;
      if (hit5->score_within_trims > max_mismatches_5 || hit3->score_within_trims > max_mismatches_3) {
	Stage3pair_free(&hitpair);
      } else if (hit5->querylength - hit5->trim_querystart - hit5->trim_queryend < min_coverage_5 &&
		 hit3->querylength - hit3->trim_querystart - hit3->trim_queryend < min_coverage_3) {
	Stage3pair_free(&hitpair);
      } else {
	newhits = Hitlist_push(newhits,hitlistpool,(void *) hitpair);
      }
    }
  }

  Hitlist_free(&hits);
  return newhits;
}


Stage3pair_T *
Stage3pair_eval_and_sort (int npaths, int *first_absmq, int *second_absmq,
			  Stage3pair_T *stage3pairarray,
			  char *queryuc_ptr_5, char *queryuc_ptr_3,
			  char *quality_string_5, char *quality_string_3) {
  float maxlik, loglik;

  float total, q;
  int mapq_score;

  int compute_npaths;
  int randomi, i;
  Stage3pair_T temp, hitpair;

  if (npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (npaths == 1) {
    hitpair = stage3pairarray[0];
    hitpair->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    hitpair->mapq_score = MAPQ_max_quality_score(quality_string_5,hitpair->hit5->querylength);
    if ((mapq_score = MAPQ_max_quality_score(quality_string_3,hitpair->hit3->querylength)) > stage3pairarray[0]->mapq_score) {
      hitpair->mapq_score = mapq_score;
    }
    hitpair->absmq_score = MAPQ_MAXIMUM_SCORE;

    Stage3end_display_prep(hitpair->hit5,queryuc_ptr_5,hitpair->alts_resolve_5,/*first_read_p*/true);
    Stage3end_display_prep(hitpair->hit3,queryuc_ptr_3,hitpair->alts_resolve_3,/*first_read_p*/false);
    if (hitpair->alts_resolve_5 >= 0 || hitpair->alts_resolve_3 >= 0) {
      hitpair->insertlength = compute_insertlength(&hitpair->pair_relationship,hitpair);
      /* assert((int) hitpair->insertlength > 0);  -- can fail for bad resolutions */
    }

    *first_absmq = hitpair->absmq_score;
    *second_absmq = 0;

  } else {

    /* Resolve ambiguities, needed for computing mapq */
    for (i = 0; i < npaths; i++) {
      hitpair = stage3pairarray[i];
      Stage3end_display_prep(hitpair->hit5,queryuc_ptr_5,hitpair->alts_resolve_5,/*first_read_p*/true);
      Stage3end_display_prep(hitpair->hit3,queryuc_ptr_3,hitpair->alts_resolve_3,/*first_read_p*/false);
      if (hitpair->alts_resolve_5 >= 0 || hitpair->alts_resolve_3 >= 0) {
	hitpair->insertlength = compute_insertlength(&hitpair->pair_relationship,hitpair);
      }
    }


    /* Compute mapq_loglik */
    for (i = 0; i < npaths; i++) {
      hitpair = stage3pairarray[i];
      hitpair->mapq_loglik =
	Stage3end_compute_mapq(hitpair->hit5,quality_string_5);
      hitpair->mapq_loglik +=
	Stage3end_compute_mapq(hitpair->hit3,quality_string_3);
    }

    /* Sort by nmatches, then mapq, and then insert length */
    qsort(stage3pairarray,npaths,sizeof(Stage3pair_T),Stage3pair_output_cmp);

    if (want_random_p) {
      /* Randomize among best alignments */
      i = 1;
      while (i < npaths && Stage3pair_output_cmp(&(stage3pairarray[i]),&(stage3pairarray[0])) == 0) {
	i++;
      }
      if (i > 1) {		/* i is number of ties */
	/* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
	randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) i));
	/* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
	temp = stage3pairarray[0];
	stage3pairarray[0] = stage3pairarray[randomi];
	stage3pairarray[randomi] = temp;
      }
    }

    /* Enforce monotonicity */
    for (i = npaths - 1; i > 0; i--) {
      if (stage3pairarray[i-1]->mapq_loglik < stage3pairarray[i]->mapq_loglik) {
	stage3pairarray[i-1]->mapq_loglik = stage3pairarray[i]->mapq_loglik;
      }
    }
    maxlik = stage3pairarray[0]->mapq_loglik;

    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < npaths; i++) {
      stage3pairarray[i]->mapq_loglik -= maxlik;
    }

#if 0
    /* Save on computation if possible */
    /* Doesn't work */
    if (npaths < maxpaths) {
      compute_npaths = npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }
#else
    compute_npaths = npaths;
#endif


    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3pairarray[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3pairarray[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3pairarray[0]->absmq_score;
    *second_absmq = stage3pairarray[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < npaths; i++) {
      total += (stage3pairarray[i]->mapq_loglik = fasterexp(stage3pairarray[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3pairarray[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3pairarray[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3pairarray[i]->mapq_score = 96;
      } else {
	stage3pairarray[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3pairarray[0]->mapq_score >= mapq_unique_score &&
	stage3pairarray[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3pair_free(&(stage3pairarray[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3pairarray;
}


static List_T
Stage3pair_optimal_score_prefinal (bool *eliminatedp, List_T hitpairlist,
				   Hitlistpool_T hitlistpool, int querylength5, int querylength3) {
  List_T optimal = NULL, p, q;
  Stage3pair_T hitpair;
  T hit5, hit3;
  Substring_T substring;
  Junction_T junction;
  int cutoff_level_5, cutoff_level_3;
  int n;
  int minscore5 = querylength5, minscore3 = querylength3, minscore = querylength5 + querylength3;
#ifdef USE_OPTIMAL_SCORE_BINGO
  int minscore_bingo = querylength5 + querylength3;
#endif
  int trim_querystart_5 = 0, trim_queryend_5 = 0, trim_querystart_3 = 0, trim_queryend_3 = 0,
    trim_querystart_0, trim_queryend_0;


#if 0 /* DISTANT_SPLICE_SPECIAL */
  bool shortdistance_p = false;
#endif


  *eliminatedp = false;
  n = List_length(hitpairlist);
  debug8(printf("\nEntered Stage3pair_optimal_score_prefinal with %d hitpairs\n",n));
  
  if (n <= 1) {
    return hitpairlist;
  }


  /* Use eventrim for comparing alignments.  Previously picked
     smallest trims, but now picking largest ones */
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;

    debug8(printf("hit5 %u..%u method %s, nsegments %d, nindels %d, trim_querystart: %d%s, trim_queryend %d%s, start_ambig %d, end_ambig %d.  hit3 %u..%u method %s, nsegments %d, nindels %d, trim_querystart %d%s, trim_queryend %d%s, start_ambig %d, end_ambig %d, sensedirs %d and %d, splice scores %f and %f.\n",
		  hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,Method_string(hit5->method),
		  hit5->nsegments,hit5->nindels,hit5->trim_querystart,hit5->trim_querystart_splicep ? " (splice)" : "",
		  hit5->trim_queryend,hit5->trim_queryend_splicep ? " (splice)" : "",
		  start_amb_length(hit5),end_amb_length(hit5),
		  hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset,Method_string(hit3->method),
		  hit3->nsegments,hit3->nindels,hit3->trim_querystart,hit3->trim_querystart_splicep ? " (splice)" : "",
		  hit3->trim_queryend,hit3->trim_queryend_splicep ? " (splice)" : "",
		  start_amb_length(hit3),end_amb_length(hit3),hit5->sensedir,hit3->sensedir,hit5->splice_score,hit3->splice_score));

    if (hit5->trim_querystart_splicep == true) {
      /* Skip */
    } else if (hit5->trim_querystart > trim_querystart_5) {
      trim_querystart_5 = hit5->trim_querystart;
    }
    if (hit5->trim_queryend_splicep == true) {
      /* Skip */
    } else if (hit5->trim_queryend > trim_queryend_5) {
      trim_queryend_5 = hit5->trim_queryend;
    }

    if (hit3->trim_querystart_splicep == true) {
      /* Skip */
    } else if (hit3->trim_querystart > trim_querystart_3) {
      trim_querystart_3 = hit3->trim_querystart;
    }
    if (hit3->trim_queryend_splicep == true) {
      /* Skip */
    } else if (hit3->trim_queryend > trim_queryend_3) {
      trim_queryend_3 = hit3->trim_queryend;
    }
  }

  if (trim_querystart_5 == querylength5) {
    trim_querystart_5 = 0;
  }
  if (trim_queryend_5 == querylength5) {
    trim_queryend_5 = 0;
  }
  if (trim_querystart_3 == querylength3) {
    trim_querystart_3 = 0;
  }
  if (trim_queryend_3 == querylength3) {
    trim_queryend_3 = 0;
  }

  debug8(printf("overall 5': trim_querystart %d, trim_queryend %d\n",trim_querystart_5,trim_queryend_5));
  debug8(printf("overall 3': trim_querystart %d, trim_queryend %d\n",trim_querystart_3,trim_queryend_3));


  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;

#ifdef CONSIDER_ENDS_IN_EVAL
    hit5->score_eventrim = hit5->trim_querystart / 8 + hit5->trim_queryend / 8;
#else
    hit5->score_eventrim = 0;
#endif

    debug8(printf("score 5' OTHER:"));

    if (trim_querystart_5 + trim_queryend_5 >= querylength5) {
      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	hit5->score_eventrim += Substring_nmismatches_bothdiff(substring);
      }

    } else {
      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	trim_querystart_0 = trim_querystart_5;
	trim_queryend_0 = trim_queryend_5;
	if (Substring_mandatory_trim_querystart(substring) > trim_querystart_0) {
	  trim_querystart_0 = Substring_mandatory_trim_querystart(substring);
	}
	if (Substring_mandatory_trim_queryend(substring) > trim_queryend_0) {
	  trim_queryend_0 = Substring_mandatory_trim_queryend(substring);
	}
	hit5->score_eventrim += Substring_count_mismatches_region(substring,trim_querystart_0,trim_queryend_0);
	debug8(printf("  substring (%d..%d) %d.",trim_querystart_5,trim_queryend_5,
		      Substring_count_mismatches_region(substring,trim_querystart_0,trim_queryend_0)));
      }
    }

    for (q = hit5->junctions_1toN; q != NULL; q = List_next(q)) {
      junction = (Junction_T) List_head(q);
      if (Junction_nindels(junction) > 0) {
	hit5->score_eventrim += indel_penalty_middle;
	debug8(printf(" => add %d.",indel_penalty_middle));
      }
    }


#if 0
    /* Accept a single indel */
#ifdef SCORE_INDELS_EVENTRIM
    if (hit5->hittype == INSERTION || hit5->hittype == DELETION) {
      debug8(printf("  indel at %d",hit5->indel_pos));
      if (hit5->indel_pos > trim_querystart_5 && hit5->indel_pos < querylength5 - trim_queryend_5) {
	hit5->score_eventrim += indel_penalty_middle;
	debug8(printf(" => add %d.",indel_penalty_middle));
      }
    }
#endif
#endif
    debug8(printf("  RESULT: %d\n",hit5->score_eventrim));

    if (hitpair->hit5->score_eventrim < minscore5) {
      minscore5 = hitpair->hit5->score_eventrim;
    }


#ifdef CONSIDER_ENDS_IN_EVAL
    hit3->score_eventrim = hit3->trim_querystart / 8 + hit3->trim_queryend / 8;
#else
    hit3->score_eventrim = 0;
#endif

    debug8(printf("score 3' OTHER:"));

    if (trim_querystart_3 + trim_queryend_3 >= querylength3) {
      for (q = hit3->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	hit3->score_eventrim += Substring_nmismatches_bothdiff(substring);
      }

    } else {
      for (q = hit3->substrings_1toN; q != NULL; q = List_next(q)) {
	substring = (Substring_T) List_head(q);
	trim_querystart_0 = trim_querystart_3;
	trim_queryend_0 = trim_queryend_3;
	if (Substring_mandatory_trim_querystart(substring) > trim_querystart_0) {
	  trim_querystart_0 = Substring_mandatory_trim_querystart(substring);
	}
	if (Substring_mandatory_trim_queryend(substring) > trim_queryend_0) {
	  trim_queryend_0 = Substring_mandatory_trim_queryend(substring);
	}
	hit3->score_eventrim += Substring_count_mismatches_region(substring,trim_querystart_0,trim_queryend_0);
	debug8(printf("  substring (%d..%d) %d.",trim_querystart_3,trim_queryend_3,
		      Substring_count_mismatches_region(substring,trim_querystart_0,trim_queryend_0)));
      }
    }

    for (q = hit3->junctions_1toN; q != NULL; q = List_next(q)) {
      junction = (Junction_T) List_head(q);
      if (Junction_nindels(junction) > 0) {
	hit3->score_eventrim += indel_penalty_middle;
	debug8(printf(" => add %d.",indel_penalty_middle));
      }
    }

#if 0
    /* Accept a single indel */
#ifdef SCORE_INDELS_EVENTRIM
    if (hit3->hittype == INSERTION || hit3->hittype == DELETION) {
      debug8(printf("  indel at %d",hit3->indel_pos));
      if (hit3->indel_pos > trim_querystart_3 && hit3->indel_pos < querylength3 - trim_queryend_3) {
	hit3->score_eventrim += indel_penalty_middle;
	debug8(printf(" => add %d.",indel_penalty_middle));
      }
    }
#endif
#endif
    debug8(printf("  RESULT: %d\n",hit3->score_eventrim));

    if (hitpair->hit3->score_eventrim < minscore3) {
      minscore3 = hitpair->hit3->score_eventrim;
    }


    /* Compute for hitpair */
    debug8(printf("hitpair score_eventrim %d = %d + %d\n",
		  hit5->score_eventrim + hit3->score_eventrim,
		  hit5->score_eventrim,hit3->score_eventrim));
    hitpair->score_eventrim = hit5->score_eventrim + hit3->score_eventrim;
    if (hitpair->score_eventrim < minscore) {
      minscore = hitpair->score_eventrim;
    }

  }
  debug8(printf("MINSCORE: %d\n",minscore));


  /* Prefinal: Use score_eventrim */
  debug8(printf("Stage3pair_optimal_score_prefinal over %d pairs: minscore = %d and %d + subopt:%d\n",
		n,minscore5,minscore3,subopt_levels));

  /* finalp == false.  Add suboptimal_mismatches to each end. */
  minscore5 += subopt_levels;
  minscore3 += subopt_levels;
  cutoff_level_5 = minscore5;
  cutoff_level_3 = minscore3;

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if (hitpair->hit5->score_eventrim > cutoff_level_5 + SCORE_EVENTRIM_SLOP && hitpair->hit3->score_eventrim > cutoff_level_3 + SCORE_EVENTRIM_SLOP) {
      debug8(printf("Prefinal: Eliminating hit pair %p at %u..%u|%u..%u with score_eventrim_5 %d > cutoff_level_5 %d and score_eventrim_3 %d > cutoff_level_3 %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->hit5->score_eventrim,cutoff_level_5,hitpair->hit3->score_eventrim,cutoff_level_3,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      Stage3pair_free(&hitpair);
      *eliminatedp = true;

    } else {
      debug8(printf("Prefinal: Keeping hit pair %p at %u..%u|%u..%u with score_eventrim_5 %d > cutoff_level_5 %d and score_eventrim_3 %d <= cutoff_level_3 %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->hit5->score_eventrim,cutoff_level_5,hitpair->hit3->score_eventrim,cutoff_level_3,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
    }
  }
  Hitlist_free(&hitpairlist);


#if 0
  /* Filter on nsegments */
  if (finalp == true && optimal != NULL) {
    hitpairlist = optimal;
    optimal = (List_T) NULL;

    hitpair = (Stage3pair_T) hitpairlist->first;
    best_nsegments = hitpair->hit5->nsegments + hitpair->hit3->nsegments;
    best_nsegments_5 = hitpair->hit5->nsegments;
    best_nsegments_3 = hitpair->hit3->nsegments;

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->hit5->nsegments + hitpair->hit3->nsegments < best_nsegments) {
	best_nsegments = hitpair->hit5->nsegments + hitpair->hit3->nsegments;
      }
      if (hitpair->hit5->nsegments < best_nsegments_5) {
	best_nsegments_5 = hitpair->hit5->nsegments;
      }
      if (hitpair->hit3->nsegments < best_nsegments_3) {
	best_nsegments_3 = hitpair->hit3->nsegments;
      }
    }

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->hit5->nsegments + hitpair->hit3->nsegments > best_nsegments + 2) {
	debug8(printf("Eliminating hit pair %p with nsegments %d+%d, sensedirs %d and %d, splice scores %f and %f\n",
		      hitpair,hitpair->hit5->nsegments,hitpair->hit3->nsegments,
		      hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
	Stage3pair_free(&hitpair);
	*eliminatedp = true;
      } else {
	debug8(printf("Keeping hit pair %p with nsegments %d+%d, sensedirs %d and %d, splice scores %f and %f\n",
		      hitpair,hitpair->hit5->nsegments,hitpair->hit3->nsegments,
		      hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
	optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
      }
    }

    Hitlist_free(&hitpairlist);
  }
#endif


#if 0
  /* Filter on pairlength */
  if (optimal != NULL) {
    hitpairlist = optimal;
    optimal = (List_T) NULL;

    hitpair = (Stage3pair_T) hitpairlist->first;
    best_absdifflength = hitpair->absdifflength;
    best_outerlength = hitpair->outerlength;

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->absdifflength < best_absdifflength) {
	best_absdifflength = hitpair->absdifflength;
	best_outerlength = hitpair->outerlength;
      } else if (hitpair->absdifflength > best_absdifflength) {
	/* Skip */
      } else if (hitpair->outerlength < best_outerlength) {
	best_outerlength = hitpair->outerlength;
      }
    }

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->absdifflength > best_absdifflength) {
	debug8(printf("Eliminating hit pair %p with absdifflength %d\n",hitpair,hitpair->absdifflength));
	Stage3pair_free(&hitpair);
	*eliminatedp = true;
      } else if (hitpair->outerlength > best_outerlength + OUTERLENGTH_SLOP) {
	debug8(printf("Eliminating hit pair %p with outerlength %u\n",hitpair,hitpair->outerlength));
	Stage3pair_free(&hitpair);
	*eliminatedp = true;
      } else {
	debug8(printf("Keeping hit pair %p with absdifflength %d and outerlength %d\n",
		      hitpair,hitpair->absdifflength,hitpair->outerlength));
	optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
      }
    }

    Hitlist_free(&hitpairlist);
  }
#endif

  debug8(printf("Exiting Stage3pair_optimal_score_prefinal with %d hits\n",List_length(optimal)));
  return optimal;
}


static int
calc_insertlength_score (Chrpos_T insertlength) {
  if (insertlength > 80000) {
    return 2;
  } else if (insertlength > 1000) {
    return 1;
  } else {
    return 0;
  }
}


/* Desired criteria: (A) within locus: (A.1) nsegments within locus,
   to get most complete alignment; (A.2) insertlength; and (A.3)
   splice_score, to get the correct sensedir.  (B) between loci:
   nmatches (and not nmatches_to_trims), to end alignments at the
   splice site */

static List_T
Stage3pair_optimal_score_final (bool *eliminatedp, List_T hitpairlist,
				Hitlistpool_T hitlistpool, int querylength5, int querylength3) {
  List_T optimal = NULL, p;
  Stage3pair_T *hitpairs, hitpair;
  int n, i, j, k;
  int best_nsegments, nsegments;
  int best_insertlength_score, insertlength_score;
  int best_nmatches_to_trims, nmatches_to_trims;
  double max_splice_score, splice_score;
  int max_nmatches = 0, cutoff_level;
  Chrpos_T best_insertlength, best_outerlength;
  /* int trim5_left, trim5_right, trim3_left, trim3_right, min_trim; */
  bool *eliminate, keptp;

  /* Relies on Path_solve_from_diagonals to maximize the number of segments at each locus */

  *eliminatedp = false;
  n = List_length(hitpairlist);
  debug8(printf("\nEntered Stage3pair_optimal_score_final with %d hitpairs\n",n));
  
  if (n <= 1) {
    return hitpairlist;
  }

#ifdef DEBUG8
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    printf("%p %p %u..%u|%u..%u methods %s and %s, nsegments %d+%d, nmatches %d+%d (%d+%d to_trims), pairlength %u, outerlength %u, sensedirs %d and %d, splice scores %f and %f\n",
	   hitpair->hit5,hitpair->hit3,
	   hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
	   hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
	   Method_string(hitpair->hit5->method),Method_string(hitpair->hit3->method),
	   hitpair->hit5->nsegments,hitpair->hit3->nsegments,
	   hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,
	   hitpair->insertlength,hitpair->outerlength,
	   hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score);
  }
#endif


  /* Prune based on nmatches (to get the splice ends) */
  max_nmatches = 0;
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit5->nmatches_plus_spliced_trims + hitpair->hit3->nmatches_plus_spliced_trims > max_nmatches) {
      max_nmatches = hitpair->hit5->nmatches_plus_spliced_trims + hitpair->hit3->nmatches_plus_spliced_trims;
      assert(max_nmatches <= querylength5 + querylength3);
    }
  }

  /* May not want to be greedy on cutoff level here.  Might want to raise subopt_levels */
  cutoff_level = max_nmatches - subopt_levels;
  debug8(printf("cutoff level %d = max_nmatches %d\n",cutoff_level,max_nmatches));

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if (hitpair->hit5->nmatches_plus_spliced_trims + hitpair->hit3->nmatches_plus_spliced_trims < cutoff_level - NMATCHES_SLOP) {
      debug8(printf("Final (nmatches %d < %d): Eliminating hit pair %p at %u..%u|%u..%u with nmatches %d (%d+%d) (%d+%d to_trims) < cutoff_level %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair->hit5->nmatches_plus_spliced_trims + hitpair->hit3->nmatches_plus_spliced_trims,cutoff_level - NMATCHES_SLOP,
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,
		    hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,cutoff_level,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      Stage3pair_free(&hitpair);
      *eliminatedp = true;

    } else {
      debug8(printf("Final (nmatches %d >= %d): Keeping hit pair %p at %u..%u|%u..%u (%d+%d substrings) with nmatches %d (%d+%d) (%d+%d to_trims) >= cutoff_level %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair->hit5->nmatches_plus_spliced_trims + hitpair->hit3->nmatches_plus_spliced_trims,cutoff_level - NMATCHES_SLOP,
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    List_length(hitpair->hit5->substrings_1toN),List_length(hitpair->hit3->substrings_1toN),
		    hitpair->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,
		    hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,cutoff_level,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
    }
  }
  Hitlist_free(&hitpairlist);
  hitpairlist = optimal;


  /* Prune based on nmatches_to_trims */
  optimal = (List_T) NULL;

  best_nmatches_to_trims = 0;
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit5->nmatches_to_trims + hitpair->hit3->nmatches_to_trims > best_nmatches_to_trims) {
      best_nmatches_to_trims = hitpair->hit5->nmatches_to_trims + hitpair->hit3->nmatches_to_trims;
      assert(best_nmatches_to_trims <= querylength5 + querylength3);
    }
  }

  cutoff_level = best_nmatches_to_trims - subopt_levels;
  debug8(printf("cutoff level %d = best_nmatches_to_trims %d\n",cutoff_level,best_nmatches_to_trims));

  /* Do not allow slop for final */
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if (hitpair->hit5->nmatches_to_trims + hitpair->hit3->nmatches_to_trims < cutoff_level /*- NMATCHES_TO_TRIMS_SLOP*/) {
      debug8(printf("Final (nmatches_to_trims %d < %d): Eliminating hit pair %p at %u..%u|%u..%u with nmatches %d (%d+%d) (%d+%d to_trims) < cutoff_level %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair->hit5->nmatches_to_trims + hitpair->hit3->nmatches_to_trims,cutoff_level /*- NMATCHES_TO_TRIMS_SLOP*/,
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,
		    hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,cutoff_level,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      Stage3pair_free(&hitpair);
      *eliminatedp = true;

    } else {
      debug8(printf("Final (nmatches %d_to_trims >= %d): Keeping hit pair %p at %u..%u|%u..%u (%d+%d substrings) with nmatches %d (%d+%d) (%d+%d to_trims) >= cutoff_level %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair->hit5->nmatches_to_trims + hitpair->hit3->nmatches_to_trims,cutoff_level /*- NMATCHES_TO_TRIMS_SLOP*/,
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    List_length(hitpair->hit5->substrings_1toN),List_length(hitpair->hit3->substrings_1toN),
		    hitpair->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,
		    hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,cutoff_level,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
    }
  }
  Hitlist_free(&hitpairlist);



  /* Eliminate within loci (1) */
  hitpairlist = optimal;
  optimal = (List_T) NULL;

  keptp = false;
  hitpairs = (Stage3pair_T *) List_to_array_n(&n,hitpairlist);
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp);
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hitpair_overlap_p(hitpairs[j],hitpairs[i]) == true) {
      j++;
    }
    if (j - i > 1) {
      debug8(printf("Found a group from %d to %d\n",i,j));
      best_nmatches_to_trims = 0;
      for (k = i; k < j; k++) {
	hitpair = hitpairs[k];
	if ((nmatches_to_trims = hitpair->hit5->nmatches_to_trims + hitpair->hit3->nmatches_to_trims) > best_nmatches_to_trims) {
	  best_nmatches_to_trims = nmatches_to_trims;
	}
      }
      debug8(printf("best_nmatches_to_trims %d\n",best_nmatches_to_trims));
      
      for (k = i; k < j; k++) {
	hitpair = hitpairs[k];
	/* Do not allow slop for final */
	if ((nmatches_to_trims = hitpair->hit5->nmatches_to_trims + hitpair->hit3->nmatches_to_trims) < best_nmatches_to_trims /*- NMATCHES_TO_TRIMS_SLOP*/) {
	  debug8(printf("Within loci (nmatches_to_trims %d < %d): Marking hit pair %p for elimination at %u..%u|%u..%u with nsegments %d+%d, pairlength %u, nmatches %d+%d (%d+%d to_trims), sensedirs %d and %d, splice scores %f and %f\n",
			nmatches_to_trims,best_nmatches_to_trims - NMATCHES_TO_TRIMS_SLOP,
			hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
			hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
			hitpair->hit5->nsegments,hitpair->hit3->nsegments,hitpair->insertlength,
			hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,
			hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
	  eliminate[k] = true;
	} else {
	  keptp = true;
	}
      }
    }

    i = j;
  }
  
  if (keptp == false) {
    optimal = hitpairlist;
  } else {
    for (k = 0; k < n; k++) {
      hitpair = hitpairs[k];
      if (eliminate[k] == true) {
	debug8(printf("Within loci: Eliminating hit pair %p at %u..%u|%u..%u with nsegments %d+%d, pairlength %u, nmatches %d+%d (%d+%d to_trims), sensedirs %d and %d, splice scores %f and %f\n",
		      hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      hitpair->hit5->nsegments,hitpair->hit3->nsegments,hitpair->insertlength,
		      hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,
		      hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
	Stage3pair_free(&hitpair);
	*eliminatedp = true;
      } else {
	optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
      }
    }
    Hitlist_free(&hitpairlist);
  }
  FREE(hitpairs);
  FREE(eliminate);


  /* Eliminate within loci (2) */
  hitpairlist = optimal;
  optimal = (List_T) NULL;

  keptp = false;
  hitpairs = (Stage3pair_T *) List_to_array_n(&n,hitpairlist);
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp);
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hitpair_overlap_p(hitpairs[j],hitpairs[i]) == true) {
      j++;
    }
    if (j - i > 1) {
      debug8(printf("Found a group from %d to %d\n",i,j));
      best_nsegments = 0;
      best_insertlength_score = 99;
      max_splice_score = 0.0;
      for (k = i; k < j; k++) {
	hitpair = hitpairs[k];
	if ((nsegments = hitpair->hit5->nsegments + hitpair->hit3->nsegments) > best_nsegments) {
	  best_nsegments = nsegments;
	  best_insertlength_score = calc_insertlength_score(hitpair->insertlength);
	  max_splice_score = hitpair->hit5->splice_score + hitpair->hit3->splice_score;

	} else if (nsegments == best_nsegments) {
	  if ((insertlength_score = calc_insertlength_score(hitpair->insertlength)) < best_insertlength_score) {
	    best_insertlength_score = insertlength_score;
	    max_splice_score = hitpair->hit5->splice_score + hitpair->hit3->splice_score;

	  } else if (insertlength_score == best_insertlength_score) {
	    if ((splice_score = hitpair->hit5->splice_score + hitpair->hit3->splice_score) > max_splice_score) {
	      max_splice_score = splice_score;
	    }
	  }
	}
      }
      debug8(printf("best_nsegments %d, best_insertlength_score %d, max_splice_score %f\n",
		    best_nsegments,best_insertlength_score,max_splice_score));
      
      for (k = i; k < j; k++) {
	hitpair = hitpairs[k];
	if ((nsegments = hitpair->hit5->nsegments + hitpair->hit3->nsegments) < best_nsegments) {
	  debug8(printf("Within loci (nsegments %d < %d): Marking hit pair %p for elimination at %u..%u|%u..%u with nsegments %d+%d, pairlength %u, nmatches %d+%d (%d+%d to_trims), sensedirs %d and %d, splice scores %f and %f\n",
			nsegments,best_nsegments,
			hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
			hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
			hitpair->hit5->nsegments,hitpair->hit3->nsegments,hitpair->insertlength,
			hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,
			hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
	  eliminate[k] = true;

	} else if (calc_insertlength_score(hitpair->insertlength) > best_insertlength_score) {
	  debug8(printf("Within loci (insertlength score %d > %d): Marking hit pair %p for elimination at %u..%u|%u..%u with nsegments %d+%d, pairlength %u, nmatches %d+%d (%d+%d to_trims), sensedirs %d and %d, splice scores %f and %f\n",
			calc_insertlength_score(hitpair->insertlength),best_insertlength_score,
			hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
			hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
			hitpair->hit5->nsegments,hitpair->hit3->nsegments,hitpair->insertlength,
			hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,
			hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
	  eliminate[k] = true;

	} else if (hitpair->hit5->splice_score + hitpair->hit3->splice_score < max_splice_score) {
	  debug8(printf("Within loci (splice_score %f < %f): Marking hit pair %p for elimination at %u..%u|%u..%u with nsegments %d+%d, pairlength %u, nmatches %d+%d (%d+%d to_trims), sensedirs %d and %d, splice scores %f and %f\n",
			hitpair->hit5->splice_score + hitpair->hit3->splice_score,max_splice_score,
			hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
			hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
			hitpair->hit5->nsegments,hitpair->hit3->nsegments,hitpair->insertlength,
			hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,
			hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
	  eliminate[k] = true;

	} else {
	  keptp = true;
	}
      }
    }

    i = j;
  }
  
  if (keptp == false) {
    optimal = hitpairlist;
  } else {
    for (k = 0; k < n; k++) {
      hitpair = hitpairs[k];
      if (eliminate[k] == true) {
	debug8(printf("Within loci: Eliminating hit pair %p at %u..%u|%u..%u with nsegments %d+%d, pairlength %u, nmatches %d+%d (%d+%d to_trims), sensedirs %d and %d, splice scores %f and %f\n",
		      hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      hitpair->hit5->nsegments,hitpair->hit3->nsegments,hitpair->insertlength,
		      hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,
		      hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
	Stage3pair_free(&hitpair);
	*eliminatedp = true;
      } else {
	optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
      }
    }
    Hitlist_free(&hitpairlist);
  }
  FREE(hitpairs);
  FREE(eliminate);



#if 0
  /* Filter on trim amount */
  hitpairlist = optimal;
  optimal = (List_T) NULL;
  min_trim = querylength5 + querylength3;
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if (hitpair->hit5->trim_querystart_splicep == true) {
      /* Skip */
      trim5_left = 0;
    } else {
      trim5_left = hitpair->hit5->trim_querystart;
    }
    if (hitpair->hit5->trim_queryend_splicep == true) {
      /* Skip */
      trim5_right = 0;
    } else {
      trim5_right = hitpair->hit5->trim_queryend;
    }

    if (hitpair->hit3->trim_querystart_splicep == true) {
      /* Skip */
      trim3_left = 0;
    } else {
      trim3_left = hitpair->hit3->trim_querystart;
    }
    if (hitpair->hit3->trim_queryend_splicep == true) {
      /* Skip */
      trim3_right = 0;
    } else {
      trim3_right = hitpair->hit3->trim_queryend;
    }

    if (trim5_left + trim5_right + trim3_left + trim3_right < min_trim) {
      min_trim = trim5_left + trim5_right + trim3_left + trim3_right;
    }
  }

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if (hitpair->hit5->trim_querystart_splicep == true) {
      /* Skip */
      trim5_left = 0;
    } else {
      trim5_left = hitpair->hit5->trim_querystart;
    }
    if (hitpair->hit5->trim_queryend_splicep == true) {
      /* Skip */
      trim5_right = 0;
    } else {
      trim5_right = hitpair->hit5->trim_queryend;
    }

    if (hitpair->hit3->trim_querystart_splicep == true) {
      /* Skip */
      trim3_left = 0;
    } else {
      trim3_left = hitpair->hit3->trim_querystart;
    }
    if (hitpair->hit3->trim_queryend_splicep == true) {
      /* Skip */
      trim3_right = 0;
    } else {
      trim3_right = hitpair->hit3->trim_queryend;
    }

    if (trim5_left + trim5_right + trim3_left + trim3_right > min_trim) {
      debug8(printf("Final (trim): Eliminating hit pair %p at %u..%u|%u..%u for trim %d+%d+%d+%d\n",
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    trim5_left,trim5_right,trim3_left,trim3_right));
      Stage3pair_free(&hitpair);
      *eliminatedp = true;

    } else {
      debug8(printf("Final (trim): Keeping hit pair %p at %u..%u|%u..%u (%d+%d substrings) for trim %d+%d+%d+%d\n",
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    List_length(hitpair->hit5->substrings_1toN),List_length(hitpair->hit3->substrings_1toN),
		    trim5_left,trim5_right,trim3_left,trim3_right));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
    }
  }
  Hitlist_free(&hitpairlist);
#endif



  /* Then find smallest insert length and outerlength across loci */
  hitpairlist = optimal;
  optimal = (List_T) NULL;

  best_insertlength = (Chrpos_T) -1;
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->insertlength < best_insertlength) {
      best_insertlength = hitpair->insertlength;
    }
  }

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if (hitpair->insertlength > best_insertlength) {
      debug8(printf("Final (insertlength %u > %u): Eliminating hit pair %p at %u..%u|%u..%u with nmatches %d (%d+%d) (%d+%d to_trims) < cutoff_level %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair->insertlength,best_insertlength,
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,
		    hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,cutoff_level,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      Stage3pair_free(&hitpair);
      *eliminatedp = true;

    } else {
      debug8(printf("Final (insertlength %u): Keeping hit pair %p at %u..%u|%u..%u (%d+%d substrings) with nmatches %d (%d+%d) (%d+%d to_trims) >= cutoff_level %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair->insertlength,
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    List_length(hitpair->hit5->substrings_1toN),List_length(hitpair->hit3->substrings_1toN),
		    hitpair->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,
		    hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,cutoff_level,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
    }
  }
  Hitlist_free(&hitpairlist);


  /* Finally find smallest outerlength across loci */
  hitpairlist = optimal;
  optimal = (List_T) NULL;

  best_outerlength = (Chrpos_T) -1;
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->outerlength < best_outerlength) {
      best_outerlength = hitpair->outerlength;
    }
  }

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if (hitpair->outerlength > best_outerlength + OUTERLENGTH_SLOP) {
      debug8(printf("Final (outerlength %u > %u): Eliminating hit pair %p at %u..%u|%u..%u with nmatches %d (%d+%d) (%d+%d to_trims) < cutoff_level %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair->outerlength,best_outerlength + OUTERLENGTH_SLOP,
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,
		    hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,cutoff_level,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      Stage3pair_free(&hitpair);
      *eliminatedp = true;

    } else {
      debug8(printf("Final (outerlength %u): Keeping hit pair %p at %u..%u|%u..%u (%d+%d substrings) with nmatches %d (%d+%d) (%d+%d to_trims) >= cutoff_level %d, sensedirs %d and %d, splice scores %f and %f\n",
		    hitpair->outerlength,
		    hitpair,hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    List_length(hitpair->hit5->substrings_1toN),List_length(hitpair->hit3->substrings_1toN),
		    hitpair->nmatches_plus_spliced_trims,hitpair->hit5->nmatches_plus_spliced_trims,hitpair->hit3->nmatches_plus_spliced_trims,
		    hitpair->hit5->nmatches_to_trims,hitpair->hit3->nmatches_to_trims,cutoff_level,
		    hitpair->hit5->sensedir,hitpair->hit3->sensedir,hitpair->hit5->splice_score,hitpair->hit3->splice_score));
      optimal = Hitlist_push(optimal,hitlistpool,(void *) hitpair);
    }
  }
  Hitlist_free(&hitpairlist);



  debug8(printf("Exiting Stage3pair_optimal_score_final with %d hits\n",List_length(optimal)));
  return optimal;
}




List_T
Stage3pair_optimal_score (List_T hitpairlist, Hitlistpool_T hitlistpool,
			  int querylength5, int querylength3, bool finalp) {
  List_T optimal;
  bool eliminatedp;

  if (finalp == false) {
    optimal = Stage3pair_optimal_score_prefinal(&eliminatedp,hitpairlist,hitlistpool,
						querylength5,querylength3);
    while (eliminatedp == true) {
      optimal = Stage3pair_optimal_score_prefinal(&eliminatedp,optimal,hitlistpool,
						  querylength5,querylength3);
    }

  } else {
    optimal = Stage3pair_optimal_score_final(&eliminatedp,hitpairlist,hitlistpool,
					     querylength5,querylength3);
    while (eliminatedp == true) {
      optimal = Stage3pair_optimal_score_final(&eliminatedp,optimal,hitlistpool,
					       querylength5,querylength3);
    }
  }

  return optimal;
}


#if 0
/* Called when computing GMAP alignment in stage1hr.c */
bool
Stage3pair_sense_consistent_p (List_T hitpairlist) {
  Stage3pair_T hitpair;
  T hit5, hit3;
  List_T p;

  for (p = hitpairlist; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;
    if (hit5->sensedir_for_concordance == hit3->sensedir_for_concordance) {
      return true;
    }
  }
  return false;
}
#endif


/* Want to unalias plus and alias minus */
List_T
Stage3end_linearize_5 (List_T hitlist) {
  T hit;
  List_T p;
#ifdef DEBUG12
  Chrpos_T chrlength;
#endif

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (T) List_head(p);
    debug12(chrlength = hit->chrlength);
    debug12(printf("Looking at 5' end %u..%u against chrlength %u\n",
		   hit->genomicstart - hit->chroffset,hit->genomicend - hit->chroffset,chrlength));

    if (hit->circularalias == 0) {
      /* Skip */

    } else if (hit->circularalias == +1) {
      if (hit->plusp == true) {
	unalias_circular(hit);
      }

    } else if (hit->circularalias == -1) {
      if (hit->plusp == false) {
	alias_circular(hit);
      }
    }
  }

  return hitlist;
}


/* Want to alias plus and unalias minus */
List_T
Stage3end_linearize_3 (List_T hitlist) {
  T hit;
  List_T p;
#ifdef DEBUG12
  Chrpos_T chrlength;
#endif

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (T) List_head(p);
    debug12(chrlength = hit->chrlength);
    debug12(printf("Looking at 3' end %u..%u against chrlength %u\n",
		   hit->genomicstart - hit->chroffset,hit->genomicend - hit->chroffset,chrlength));

    if (hit->circularalias == 0) {
      /* Skip */

    } else if (hit->circularalias == -1) {
      if (hit->plusp == true) {
	alias_circular(hit);
      }

    } else if (hit->circularalias == +1) {
      if (hit->plusp == false) {
	unalias_circular(hit);
      }
    }
  }

  return hitlist;
}



List_T
Stage3pair_remove_circular_alias (List_T hitpairlist, Hitlistpool_T hitlistpool) {
  List_T newlist = NULL, p;
  Stage3pair_T hitpair;
  int trim;

  debug12(printf("Stage3pair_remove_circular_alias called with %d hitpairs\n",
		 List_length(hitpairlist)));
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

#if 0
    /* Not sure if this is necessary */
    if (hitpair->hit5->circularalias == +1 && hitpair->hit3->circularalias == +1) {
      /* First, try to salvage alias +1 */
      unalias_circular(hitpair->hit5);
      unalias_circular(hitpair->hit3);
    }
#endif

    if (hitpair->hit5->plusp == true) {
      trim = hitpair->hit5->trim_querystart;
    } else {
      trim = hitpair->hit3->trim_queryend;
    }

    if (hitpair->low + trim >= hitpair->hit5->chroffset + hitpair->hit5->chrlength) {
      /* Both ends in circular alias */
      debug12(printf("Both ends in circular alias\n"));
      Stage3pair_free(&hitpair);

    } else {
      newlist = Hitlist_push(newlist,hitlistpool,(void *) hitpair);
    }
  }

  Hitlist_free(&hitpairlist);
  return newlist;
}


