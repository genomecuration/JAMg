static char rcsid[] = "$Id: substring.c 219220 2019-05-12 22:28:07Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "substring.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>		/* For log and exp */

#include "assert.h"
#include "mem.h"
#include "maxent_hr.h"
#include "complement.h"
#include "genome128_hr.h"
#include "mapq.h"
#include "comp.h"
#include "splice.h"
#include "sense.h"
#include "simplepair.h"		/* For Simplepair_new_out */


/* Designed to allow 1 match to offset 1 mismatch.  To handle 2 matches vs 2 mismatches, penalize for multiple mismatches */
#define TRIM_MATCH_SCORE 1
#define TRIM_MISMATCH_SCORE_LAST -1 /* Requires 1 match to compensate */
#define TRIM_MISMATCH_SCORE_MULT -4 /* Requires 4 matches to compensate */

#define KNOWN_SPLICESITE_EDGE 1	/* Cannot be 0; otherwise will hit found splice sites */

#define CONSISTENT_TEXT "consistent"
#define TRANSLOCATION_TEXT "translocation"
#define INVERSION_TEXT "inversion"
#define SCRAMBLE_TEXT "scramble"

#define MIN_EXON_LENGTH 20
#define END_SPLICESITE_PROB_MATCH 0.90
#define END_SPLICESITE_PROB_MISMATCH 0.95


#ifdef CHECK_ASSERTIONS
#define CHECK_NMISMATCHES 1
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* mark_mismatches */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Substring_new */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Substring_overlap_p and Substring_insert_length */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


/* splice site probs */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


/* Substring_convert_to_pairs */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif


/* contains known splicesite */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif


/* trimming.  may also want to turn on DEBUG8 in pair.c */
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


/* Substring_set_alt */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Substring_circularpos */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif


#define LOG_99 -0.01005033585
#define LOG_01 -4.605170186

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

#if 0
#define LOG_99_999 -0.01105083619
#define LOG_99_001 -6.917805615
#define LOG_25_001 -8.29404964
#define LOG_25_999 -1.387294861
#define LOG_01_999 -4.606170686
#define LOG_01_001 -11.51292546
#define LOG_75_001 -7.195437351
#define LOG_75_999 -0.2886825728
#endif

/* Switches on 4 consecutive mismatches */
#define LOG_99_99 -0.02010067171
#define LOG_99_01 -4.615220522
#define LOG_25_01 -5.991464547
#define LOG_25_99 -1.396344697
#define LOG_01_99 -4.615220522
#define LOG_01_01 -9.210340372
#define LOG_75_01 -4.892852258
#define LOG_75_99 -0.2977324083


static bool print_nsnpdiffs_p;
static bool print_snplabels_p;
static bool show_refdiff_p;

static IIT_T snps_iit;
static int *snps_divint_crosstable;

static Genome_T genomebits;
static Genome_T genomebits_alt;

static Univ_IIT_T chromosome_iit;
static int nchromosomes;
static int circular_typeint;

static IIT_T genes_iit;
static int *genes_divint_crosstable;

static IIT_T splicesites_iit;
static int *splicesites_divint_crosstable;

static int donor_typeint;
static int acceptor_typeint;

static bool novelsplicingp;
static bool knownsplicingp;
static Outputtype_T output_type;

static Mode_T mode;

static double genomelength;	/* For BLAST E-value */


static char complCode[128] = COMPLEMENT_LC;

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



#define T Substring_T

struct T {
  /* Maintaining this pointer so we don't have to pass in the values
     for Substring_display_prep based on plusp.  Caller needs to keep
     these values until all computation is finished. */
  Compress_T query_compress; /* Pointer to the associated query_compress */
  

  int nmismatches_bothdiff;	/* Over region left after trimming */
  int nmismatches_refdiff;	/* Over region left after trimming */
  /* nsnpdiffs = nmismatches_bothdiff - nmismatches_refdiff */

  int nmatches_plus_spliced_trims;
  int nmatches_to_trims;

  /* int trim_querystart; */
  /* int trim_queryend; */
  int mandatory_trim_querystart;	/* Resulting from unalias, plus SOFT_CLIPS_AVOID_CIRCULARIZATION */
  int mandatory_trim_queryend;	/* Resulting from alias, plus SOFT_CLIPS_AVOID_CIRCULARIZATION */
  int start_amb_length;
  int end_amb_length;

  bool trim_querystart_splicep;
  bool trim_queryend_splicep;

  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;

  Univcoord_T left; /* for plus: alignstart - querystart(orig).  for
		       minus: alignend - (querylength -
		       queryend(orig)).  Set when substring is created
		       or made unambiguous, and remains constant */

  Univcoord_T genomicstart;	/* For region corresponding to entire querylength (if extrapolated) */
  Univcoord_T genomicend;

  Endtype_T start_endtype;
  Endtype_T end_endtype;

  int querystart_pretrim;	/* From original invocation, before trimming */
  int queryend_pretrim;		/* From original invocation, before trimming */
  int querystart;		/* For part that aligns to genome, post-trim */
  int queryend;
  int amb_splice_pos;		/* Used for ambiguous substrings.  Based on left, not querystart */
  int querylength;

  Univcoord_T alignstart_trim;	/* For part that aligns to genome, excluding part that is trimmed (post-trim) */
  Univcoord_T alignend_trim;

  bool plusp;
  int genestrand;

  char *genomic_bothdiff; /* In same direction as query.  NULL if same
  			     as query.  Has dashes outside of
  			     querystart..(queryend-1) for indels and
  			     splices.  Has lowercase to indicate
  			     trimmed regions, terminal regions,
  			     deletions, splice dinucleotides, and
  			     mismatches from ref and alt.  Use for
  			     MAPQ computations and scoring. */

  char *genomic_refdiff;  /* Same as above, but lowercase for
			     mismatches from ref only.  For
			     non-SNP-tolerant alignment, this is just
			     a pointer to genomic_bothdiff.  Use for
			     NM and MD computations.  */

  float mapq_loglik;

  int sensedir;

  Univcoord_T splicecoord_D;
  int splicesitesD_knowni;	/* Needed for intragenic_splice_p in stage1hr.c */

  int siteD_pos;
  double siteD_prob;

  /* for shortexon (always use *_1 for acceptor and *_2 for donor) */
  /* for donor/acceptor: the ambiguous position */
  Univcoord_T splicecoord_A;
  int splicesitesA_knowni;

  int siteA_pos;
  double siteA_prob;

  /* Note: For DNA fusions, use both splicecoord_D and splicecoord_A */
  int siteN_pos;


  int alts_ncoords;
  Univcoord_T *alts_coords;
  int *alts_knowni;
  int *alts_nmismatches;
  double *alts_probs;

  Endtype_T amb_type;		/* Ambiguous DONs or ACCs */
};


void
Substring_alias_circular (T this) {
  Chrpos_T chrlength;
  int i;

  if (this == NULL) {
    /* Skip */

  } else if (this->alts_ncoords > 0) {
    /* Cannot rely on this->left, which is 0.  Let mandatory_trim_querystart and mandatory_trim_queryend be 0 */
    chrlength = this->chrlength;
    for (i = 0; i < this->alts_ncoords; i++) {
      this->alts_coords[i] += chrlength;
    }

  } else {
    chrlength = this->chrlength;

    if (this->left + this->querylength > this->chroffset + chrlength) {
      debug2(printf("For alias, this->left %u + this->querylength %d > offset %u + chrlength %u\n",
		    this->left,this->querylength,this->chroffset,chrlength));
      if (this->plusp == true) {
	this->mandatory_trim_queryend = (this->left + this->querylength) - (this->chroffset + chrlength);
	debug2(printf("For alias, plusp true, setting mandatory_trim_queryend to be %d\n",this->mandatory_trim_queryend));
	assert(this->mandatory_trim_queryend >= 0);
      } else {
	this->mandatory_trim_querystart = (this->left + this->querylength) - (this->chroffset + chrlength);
	debug2(printf("For alias, plusp false, setting mandatory_trim_querystart to be %d\n",this->mandatory_trim_querystart));
	assert(this->mandatory_trim_querystart >= 0);
      }
    }

    this->left += chrlength;
    this->genomicstart += chrlength;
    this->genomicend += chrlength;
    this->alignstart_trim += chrlength;
    this->alignend_trim += chrlength;
    this->splicecoord_D += chrlength;
    this->splicecoord_A += chrlength;
  }

  return;
}


void
Substring_unalias_circular (T this) {
  Chrpos_T chrlength;
  int i;

  if (this == NULL) {
    /* Skip */

  } else if (this->alts_ncoords > 0) {
    /* Cannot rely on this->left, which is 0.  Let mandatory_trim_querystart and mandatory_trim_queryend be 0 */
    chrlength = this->chrlength;
    for (i = 0; i < this->alts_ncoords; i++) {
      this->alts_coords[i] -= chrlength;
    }

  } else {
    chrlength = this->chrlength;

    if (this->left < this->chroffset + chrlength) {
      debug2(printf("For unalias, this->left %u < chroffset %u + chrlength %d\n",
		    this->left,this->chroffset,chrlength));
      if (this->plusp == true) {
	this->mandatory_trim_querystart = (this->chroffset + chrlength) - this->left;
	debug2(printf("For unalias, plusp true, setting mandatory_trim_querystart to be %d\n",this->mandatory_trim_querystart));
	assert(this->mandatory_trim_querystart >= 0);
      } else {
	this->mandatory_trim_queryend = (this->chroffset + chrlength) - this->left;
	debug2(printf("For unalias, plusp false, setting mandatory_trim_queryend to be %d\n",this->mandatory_trim_queryend));
	assert(this->mandatory_trim_queryend >= 0);
      }
    }

    this->left -= chrlength;
    this->genomicstart -= chrlength;
    this->genomicend -= chrlength;
    this->alignstart_trim -= chrlength;
    this->alignend_trim -= chrlength;
    this->splicecoord_D -= chrlength;
    this->splicecoord_A -= chrlength;
  }

  return;
}



static void
fill_w_dashes (char *string, int start, int end) {
  int i;

  for (i = start; i < end; i++) {
    string[i] = '-';
  }
  return;
}

static void
fill_w_stars (char *string, int start, int end) {
  int i;

  for (i = start; i < end; i++) {
    string[i] = '*';
  }
  return;
}



void
Substring_free (T *old) {
  if (*old) {
    debug2(printf("Freeing substring %p\n",*old));
    if ((*old)->alts_ncoords > 0) {
      FREE_OUT((*old)->alts_coords);
      FREE_OUT((*old)->alts_knowni);
      FREE_OUT((*old)->alts_nmismatches);
      FREE_OUT((*old)->alts_probs);
    }
    if ((*old)->genomic_bothdiff != NULL) {
      if ((*old)->genomic_refdiff != (*old)->genomic_bothdiff) {
	FREE_OUT((*old)->genomic_refdiff);
      }
      FREE_OUT((*old)->genomic_bothdiff);
    }

    FREE_OUT(*old);
  }
  return;
}


void
Substring_list_gc (List_T *old) {
  List_T p;
  T substring;

  for (p = *old; p != NULL; p = p->rest) {
    substring = (Substring_T) p->first;
    Substring_free(&substring);
  }
  /* List_free(&(*old)); -- allocated by Listpool_push */
  return;
}



bool
Substring_contains_p (T this, int querypos) {
  if (querypos >= this->querystart && querypos < this->queryend) {
    return true;
  } else {
    return false;
  }
}


int
Substring_compare (T substring1, T substring2, int alias1, int alias2, Chrpos_T chrlength1, Chrpos_T chrlength2) {
  Univcoord_T alignstart1, alignend1, alignstart2, alignend2;

  if (substring1 == NULL && substring2 == NULL) {
    return 0;
  } else if (substring1 == NULL) {
    return -1;
  } else if (substring2 == NULL) {
    return +1;
  } else {
    alignstart1 = substring1->alignstart_trim;
    alignend1 = substring1->alignend_trim;
    if (alias1 < 0) {
      alignstart1 += chrlength1;
      alignend1 += chrlength1;
    }

    alignstart2 = substring2->alignstart_trim;
    alignend2 = substring2->alignend_trim;
    if (alias2 < 0) {
      alignstart2 += chrlength2;
      alignend2 += chrlength2;
    }

    if (alignstart1 < alignstart2) {
      return -1;
    } else if (alignstart1 > alignstart2) {
      return +1;
    } else if (alignend1 < alignend2) {
      return -1;
    } else if (alignend1 > alignend2) {
      return +1;
    } else {
      return 0;
    }
  }
}


bool
Substring_equal_p (T substring1, T substring2) {
  Univcoord_T low1, high1, low2, high2;

  if (substring1->plusp == true) {
    low1 = substring1->alignstart_trim;
    high1 = substring1->alignend_trim;
  } else {
    low1 = substring1->alignend_trim;
    high1 = substring1->alignstart_trim;
  }
  if (high1 > 0) {
    high1 -= 1;
  }

  if (substring2->plusp == true) {
    low2 = substring2->alignstart_trim;
    high2 = substring2->alignend_trim;
  } else {
    low2 = substring2->alignend_trim;
    high2 = substring2->alignstart_trim;
  }
  if (high2 > 0) {
    high2 -= 1;
  }

  debug3(printf("Checking equivalence between %u..%u and %u..%u",low1,high1,low2,high2));

  if (low1 == low2 && high1 == high2) {
    return true;
  } else {
    return false;
  }
}


bool
Substring_overlap_p (T substring1, T substring2) {
  Univcoord_T low1, high1, low2, high2;

  if (substring1->plusp == true) {
    low1 = substring1->alignstart_trim;
    high1 = substring1->alignend_trim;
  } else {
    low1 = substring1->alignend_trim;
    high1 = substring1->alignstart_trim;
  }
  if (high1 > 0) {
    high1 -= 1;
  }

  if (substring2->plusp == true) {
    low2 = substring2->alignstart_trim;
    high2 = substring2->alignend_trim;
  } else {
    low2 = substring2->alignend_trim;
    high2 = substring2->alignstart_trim;
  }
  if (high2 > 0) {
    high2 -= 1;
  }

  debug3(printf("Checking overlap between %u..%u and %u..%u",low1,high1,low2,high2));

  if (high2 < low1) {
    debug3(printf(" => no because %u < %u\n",high2,low1));
    return false;
  } else if (low2 > high1) {
    debug3(printf(" => no because %u > %u\n",low2,high1));
    return false;
  } else {
    debug3(printf(" => yes\n"));
    return true;
  }
}


Chrpos_T
Substring_insert_length (int *pair_relationship, T substring5, T substring3) {
  Univcoord_T low, high, low5, high5, low3, high3;

  if (substring5->plusp == true) {
    low5 = substring5->genomicstart;
    high5 = substring5->genomicend;
  } else {
    low5 = substring5->genomicend;
    high5 = substring5->genomicstart;
  }

  if (substring3->plusp == true) {
    low3 = substring3->genomicstart;
    high3 = substring3->genomicend;
  } else {
    low3 = substring3->genomicend;
    high3 = substring3->genomicstart;
  }

  if (low5 < low3) {
    low = low5;
  } else {
    low = low3;
  }

  if (high5 > high3) {
    high = high5;
  } else {
    high = high3;
  }

  if (low5 < low3 && high5 < high3) {
    *pair_relationship = +1;
  } else if (low3 < low5 && high3 < high5) {
    *pair_relationship = -1;
  } else if (substring5->plusp == true && substring3->plusp == true) {
    *pair_relationship = +1;
  } else if (substring5->plusp == false && substring3->plusp == false) {
    *pair_relationship = -1;
  } else {
    *pair_relationship = 0;
  }

  debug3(printf("Returning %u - %u.  pair_relationship %d\n",high,low,*pair_relationship));
  return high - low;
}


bool
Substring_overlap_point_trimmed_p (T substring, Univcoord_T endpos) {
  Univcoord_T low, high;

  if (substring == NULL) {
    return false;
  }

  if (substring->plusp == true) {
    low = substring->alignstart_trim;
    high = substring->alignend_trim;
    if (high > 0) {
      high -= 1;
    }
    debug3(printf("Checking overlap between plus %u..%u and %u",low,high,endpos));
  } else {
    low = substring->alignend_trim;
    high = substring->alignstart_trim;
    if (high > 0) {
      high -= 1;
    }
    debug3(printf("Checking overlap between minus %u..%u and %u",low,high,endpos));
  }


  if (endpos < low) {
    debug3(printf(" => no because %u < %u\n",endpos,low));
    return false;
  } else if (endpos > high) {
    debug3(printf(" => no because %u > %u\n",endpos,high));
    return false;
  } else {
    debug3(printf(" => yes\n"));
    return true;
  }
}


Univcoord_T
Substring_overlap_segment_trimmed (T substring1, T substring2) {
  Univcoord_T maxlow, minhigh;
  Univcoord_T low1, high1, low2, high2;

  if (substring1->plusp == true) {
    low1 = substring1->alignstart_trim;
    high1 = substring1->alignend_trim;
  } else {
    low1 = substring1->alignend_trim;
    high1 = substring1->alignstart_trim;
  }
  if (high1 > 0) {
    high1 -= 1;
  }

  if (substring2->plusp == true) {
    low2 = substring2->alignstart_trim;
    high2 = substring2->alignend_trim;
  } else {
    low2 = substring2->alignend_trim;
    high2 = substring2->alignstart_trim;
  }
  if (high2 > 0) {
    high2 -= 1;
  }

  debug3(printf("Checking overlap between %u..%u and %u..%u",low1,high1,low2,high2));

  if (high2 < low1) {
    debug3(printf(" => no because %u < %u\n",high2,low1));
    return 0;
  } else if (low2 > high1) {
    debug3(printf(" => no because %u > %u\n",low2,high1));
    return 0;
  } else {
    maxlow = (low1 > low2) ? low1 : low2;
    minhigh = (high1 < high2) ? high1 : high2;
    debug3(printf(" => yes.  maxlow %llu, minhigh %llu.  returning %llu\n",
		  maxlow,minhigh,maxlow + (minhigh - maxlow)/2));
    return maxlow + (minhigh - maxlow)/2;
  }
}


static void
mark_mismatches_cmet_gsnap (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("\n"));
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'G' && query[i] == 'A') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'C' && query[i] == 'T') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }

  return;
}


static void
mark_mismatches_cmet_sam (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'G' && query[i] == 'A') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'A';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'C' && query[i] == 'T') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'T';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }
  
  return;
}



static void
mark_mismatches_atoi_gsnap (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'T' && query[i] == 'C') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'A' && query[i] == 'G') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }
  
  return;
}



static void
mark_mismatches_atoi_sam (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'T' && query[i] == 'C') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'C';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'A' && query[i] == 'G') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'G';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }

  return;
}


static void
mark_mismatches_ttoc_gsnap (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'A' && query[i] == 'G') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'T' && query[i] == 'C') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }
  
  return;
}



static void
mark_mismatches_ttoc_sam (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'A' && query[i] == 'G') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'G';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'T' && query[i] == 'C') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'C';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }

  return;
}



/* The result variable should be an existing genomic_bothdiff or genomic_refdiff, and then reassigned to that */
static char *
embellish_genomic (char *result, char *genomic_diff, char *query, int querystart, int queryend, int querylength,
		   int mandatory_trim_querystart_in, int mandatory_trim_queryend_in, int extraleft, int extraright,
		   bool plusp, int genestrand) {
  int i, j, k;
  int mandatory_trim_querystart, mandatory_trim_queryend;

  debug1(printf("Entered embellish_genomic with querystart %d, queryend %d, querylength %d, mandatory_trim_querystart %d, mandatory_trim_queryend %d\n",
		querystart,queryend,querylength,mandatory_trim_querystart_in,mandatory_trim_queryend_in));
  debug1(printf("genomic_diff %s\n",genomic_diff));

  assert(mandatory_trim_querystart_in >= 0);
  assert(mandatory_trim_queryend_in >= 0);

  if (plusp == true) {
    mandatory_trim_querystart = mandatory_trim_querystart_in;
    mandatory_trim_queryend = mandatory_trim_queryend_in;
  } else {
    mandatory_trim_querystart = mandatory_trim_queryend_in;
    mandatory_trim_queryend = mandatory_trim_querystart_in;
  }

  if (result == NULL) {
    result = (char *) MALLOC_OUT((querylength+1) * sizeof(char));
  }
  result[querylength] = '\0';
#if 0
  for (i = 0; i < querylength; i++) {
    result[i] = '?';
  }
#endif

  /* Add aligned region with lower-case diffs, surrounded by dashes */
  fill_w_dashes(result,0,querystart);
  fill_w_stars(result,0,mandatory_trim_querystart);

  /* Don't need to know adj anymore, because each substring has its own left */
  debug1(printf("Copying from genomic_diff[%d] to result[%d] for a length of %d - %d\n",querystart,querystart,queryend,querystart));
  strncpy(&(result[querystart]),&(genomic_diff[querystart]),queryend-querystart);
  debug1(printf("A g1: %s (%d..%d) extraleft:%d extraright:%d\n",result,querystart,queryend,extraleft,extraright));

  if (mode == STANDARD) {
    /* Skip */
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    mark_mismatches_cmet_gsnap(result,query,querystart,queryend,genestrand);
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    mark_mismatches_atoi_gsnap(result,query,querystart,queryend,genestrand);
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    mark_mismatches_ttoc_gsnap(result,query,querystart,queryend,genestrand);
  } else {
    abort();
  }

  fill_w_dashes(result,queryend,querylength);
  fill_w_stars(result,querylength - mandatory_trim_queryend,querylength);
  debug1(printf("B g1: %s\n",result));

  /* Add terminal ends as lower-case.  Previously stopped at i >= 0 */
  for (k = 0, i = j = querystart - 1; k < extraleft && i >= mandatory_trim_querystart /*&& j >= 0*/; k++, i--, j--) {
    result[i] = (char) tolower(genomic_diff[j]);
    /* printf("k=%d i=%d result[i]=%c\n",k,i,result[i]); */
    assert(result[i] == 'a' || result[i] == 'c' || result[i] == 'g' || result[i] == 't' || result[i] == 'n' || result[i] == '*');
  }
  
  /* Add terminal ends as lower-case.  Previously stopped at i < querylength */
  for (k = 0, i = j = queryend; k < extraright && i < querylength - mandatory_trim_queryend /*&& j < genomiclength*/; k++, i++, j++) {
    result[i] = (char) tolower(genomic_diff[j]);
    /* printf("k=%d i=%d result[i]=%c\n",k,i,result[i]); */
  }
  debug1(printf("C g1: %s\n",result));

  return result;
}


static char *
embellish_genomic_sam (char *result, char *genomic_diff, char *query, int querystart, int queryend, int querylength,
		       int mandatory_trim_querystart_in, int mandatory_trim_queryend_in, int extraleft, int extraright,
		       bool plusp, int genestrand) {
  int i, j, k;
  int start, end;
  int mandatory_trim_querystart, mandatory_trim_queryend;

  assert(mandatory_trim_querystart_in >= 0);
  assert(mandatory_trim_queryend_in >= 0);

  if (plusp == true) {
    mandatory_trim_querystart = mandatory_trim_querystart_in;
    mandatory_trim_queryend = mandatory_trim_queryend_in;
  } else {
    mandatory_trim_querystart = mandatory_trim_queryend_in;
    mandatory_trim_queryend = mandatory_trim_querystart_in;
  }


  if (result == NULL) {
    result = (char *) MALLOC_OUT((querylength+1) * sizeof(char));
  }
  result[querylength] = '\0';

  if ((start = querystart - extraleft) < 0) {
    start = 0;
  }
  if ((end = queryend + extraright) > querylength) {
    end = querylength;
  }

  fill_w_dashes(result,0,start);
  fill_w_stars(result,0,mandatory_trim_querystart);
  strncpy(&(result[start]),&(genomic_diff[start]),end-start);
  fill_w_dashes(result,end,querylength);
  fill_w_stars(result,querylength - mandatory_trim_queryend,querylength);

  if (mode == STANDARD) {
    /* Skip */
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    mark_mismatches_cmet_sam(result,query,querystart,queryend,genestrand);
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    mark_mismatches_atoi_sam(result,query,querystart,queryend,genestrand);
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    mark_mismatches_ttoc_sam(result,query,querystart,queryend,genestrand);
  } else {
    abort();
  }

  /* Add terminal ends as lower-case.  Previously stopped at i >= 0 */
  for (k = 0, i = querystart-1, j = querystart-1; i >= mandatory_trim_querystart /*&& j >= 0*/; k++, i--, j--) {
    if (query[i] == genomic_diff[j]) {
      result[i] = genomic_diff[j];
    } else {
      result[i] = (char) tolower(genomic_diff[j]);
    }
    /* printf("k=%d i=%d j=%d result[i]=%c\n",k,i,j,result[i]); */
  }

  /* Add terminal ends as lower-case.  Previously stopped at i < querylength */
  for (k = 0, i = queryend, j = queryend; i < querylength - mandatory_trim_queryend /*&& j < genomiclength*/; k++, i++, j++) {
    if (query[i] == genomic_diff[j]) {
      result[i] = genomic_diff[j];
    } else {
      result[i] = (char) tolower(genomic_diff[j]);
    }
  }

  return result;
}


int
Substring_trim_qstart_nosplice (int *nmismatches, int *mismatch_positions_alloc,
				Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				bool plusp,int genestrand) {
  int score = 0;
  int trimpos = pos5, pos, prevpos, i;
  int *mismatch_positions = mismatch_positions_alloc;
  int alignlength = pos3 - pos5;


  debug8(printf("Entered Substring_trim_qstart_nosplice with pos5 %d, pos3 %d, mismatch_scores %d/%d, match_score %d\n",
		pos5,pos3,TRIM_MISMATCH_SCORE_LAST,TRIM_MISMATCH_SCORE_MULT,TRIM_MATCH_SCORE));
  debug8(printf("Calling Genome_mismatches_right_trim with left %u, pos5 %d, pos3 %d\n",left,pos5,pos3));
  *nmismatches = Genome_mismatches_right_trim(mismatch_positions,/*max_mismatches*/alignlength,
					     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
					     left,pos5,pos3,plusp,genestrand);
  debug8(printf("%d mismatches:",*nmismatches));
  debug8(
	 for (i = 0; i <= *nmismatches; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");
	 );

  if (*nmismatches == 0) {
    return pos5;
  }

  /* Advance until we get to a good region */
  pos = prevpos = pos3;
  i = 0;
  while (i < (*nmismatches) - 1 && score <= 0) {
    pos = mismatch_positions[i]; /* position at the mismatch */
    score += (prevpos - pos - 1)*TRIM_MATCH_SCORE + TRIM_MISMATCH_SCORE_MULT;
    if (score < 0) {
      /* Distal matches did not compensate for the mismatch, so advance to here */
      debug8(printf("Advance qstart pos %d, score %d => 0, pos %d\n",pos,score,pos));
      score = 0;
    } else {
      debug8(printf("Advance qstart pos %d, score %d, pos %d\n",pos,score,pos));
    }
    prevpos = pos;			   /* On the mismatch */
    i++;
  }

  /* Now trim */
  trimpos = prevpos = pos;
  score = 0;
  debug8(printf("Trim qend pos %d, score %d, trimpos %d\n",pos,score,trimpos));
  while (i < *nmismatches) {
    pos = mismatch_positions[i];
    score += (prevpos - pos - 1)*TRIM_MATCH_SCORE + TRIM_MISMATCH_SCORE_MULT;
    if (score >= 0) {
      trimpos = pos;
      score = 0;
    }
    debug8(printf("Trim qstart pos %d, score %d, trimpos %d\n",pos,score,trimpos));
    prevpos = pos;
    i++;
  }
  
  score += (prevpos - pos5 - 1)*TRIM_MATCH_SCORE + TRIM_MISMATCH_SCORE_LAST;
  if (score >= 0) {
    trimpos = pos5 - 1;		/* At the "mismatch" */
    /* score = 0; */
  }
  debug8(printf("Final qstart pos %d, score %d, trimpos %d\n",pos,score,trimpos));

  return trimpos + 1;		/* One position after the mismatch */
}


int
Substring_trim_qend_nosplice (int *nmismatches, int *mismatch_positions_alloc,
			      Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
			      bool plusp, int genestrand) {
  int score = 0;
  int trimpos = pos3, pos, prevpos, i;
  int *mismatch_positions = mismatch_positions_alloc;
  int alignlength = pos3 - pos5;


  debug8(printf("Entered Substring_trim_qend_nosplice with pos5 %d, pos3 %d, mismatch_scores %d/%d, match_score %d\n",
		pos5,pos3,TRIM_MISMATCH_SCORE_LAST,TRIM_MISMATCH_SCORE_MULT,TRIM_MATCH_SCORE));
  debug8(printf("Calling Genome_mismatches_left_trim with left %u, pos5 %d, pos3 %d\n",left,pos5,pos3));
  *nmismatches = Genome_mismatches_left_trim(mismatch_positions,/*max_mismatches*/alignlength,
					     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
					     left,pos5,pos3,plusp,genestrand);

  debug8(printf("%d mismatches:",*nmismatches));
  debug8(
	 for (i = 0; i <= *nmismatches; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");
	 );

  if (*nmismatches == 0) {
    return pos3;
  }

  /* Advance until we get to a good region */
  pos = prevpos = pos5;
  i = 0;
  while (i < (*nmismatches) - 1 && score <= 0) {
    pos = mismatch_positions[i]; /* position at the mismatch */
    score += (pos - prevpos - 1)*TRIM_MATCH_SCORE + TRIM_MISMATCH_SCORE_MULT;
    if (score < 0) {
      /* Distal matches did not compensate for the mismatch, so advance to here */
      debug8(printf("Advance qend pos %d, score %d => 0, pos %d\n",pos,score,pos));
      score = 0;
    } else {
      debug8(printf("Advance qend pos %d, score %d, pos %d\n",pos,score,pos));
    }
    prevpos = pos;			   /* On the mismatch */
    i++;
  }
  
  /* Now trim */
  trimpos = prevpos = pos;
  score = 0;
  while (i < *nmismatches) {
    pos = mismatch_positions[i];
    score += (pos - prevpos - 1)*TRIM_MATCH_SCORE + TRIM_MISMATCH_SCORE_MULT;
    if (score >= 0) {
      trimpos = pos;
      score = 0;
    }
    debug8(printf("Trim qend pos %d, score %d, trimpos %d\n",pos,score,trimpos));
    prevpos = pos;
    i++;
  }
  
  score += (pos3 - prevpos - 1)*TRIM_MATCH_SCORE + TRIM_MISMATCH_SCORE_LAST;
  if (score >= 0) {
    trimpos = pos3;
    /* score = 0; */
  }
  debug8(printf("Final qend pos %d, score %d, trimpos %d\n",pos,score,trimpos));
    
  return trimpos;		/* qend is outside the region */
}


/* Returns true if splice was found, false otherwise.  If true, result
   is the number of spliceends.  If false, result is the single trim
   position */
bool
Substring_trimmed_qstarts (int *result, Splicetype_T *splicetype, int **ambig_qstarts, double **ambig_probs_5,
			   Univcoord_T left, int qend,
			   bool plusp, int genestrand, int *mismatch_positions_alloc, Compress_T query_compress,
			   Univcoord_T chroffset, int sensedir) {
  int trimpos;
  int nspliceends;
  int nmismatches;

  debug8(printf("\n***Entered Substring_trimmed_qstarts with left %u, genome %u..%u, qend %d..%d, plusp %d, sensedir %d\n",
		left-chroffset,left+0-chroffset,left+qend-chroffset,0,qend,plusp,sensedir));

  trimpos = Substring_trim_qstart_nosplice(&nmismatches,mismatch_positions_alloc,query_compress,left,
					   /*pos5*/0,/*pos3*/qend,plusp,genestrand);
  debug8(printf("trimpos %d (relative to %d)\n",trimpos,/*qstart*/0));

  if (trimpos >= qend) {
    debug8(printf("trimpos %d >= qend %d, so returning -1\n",trimpos,qend));
    *result = -1;
    return false;

  } else if (novelsplicingp == false || trimpos == 0) {
    *result = trimpos;
    return false;

  } else if ((nspliceends = Splice_trim_novel_spliceends_5(&(*splicetype),&(*ambig_qstarts),&(*ambig_probs_5),
							   left,/*qstart*/trimpos,qend,
							   mismatch_positions_alloc,nmismatches,
							   chroffset,plusp,sensedir)) == 0) {
    *result = trimpos;
    return false;
      
  } else {
    *result = nspliceends;
    return true;
  }
}


bool
Substring_qstart_trim (int *trimpos, Splicetype_T *splicetype, double *ambig_prob_qstart,
		       Univcoord_T left, int qend,
		       bool plusp, int genestrand, int *mismatch_positions_alloc, Compress_T query_compress,
		       Univcoord_T chroffset, int sensedir) {
  int nspliceends;
  int nmismatches;
  int *ambig_qstarts;
  double *ambig_probs_5;

  debug8(printf("\n***Entered Substring_qstart_trim with left %u, genome %u..%u, qend %d..%d, plusp %d, sensedir %d\n",
		left-chroffset,left+0-chroffset,left+qend-chroffset,0,qend,plusp,sensedir));

  *trimpos = Substring_trim_qstart_nosplice(&nmismatches,mismatch_positions_alloc,query_compress,left,
					    /*pos5*/0,/*pos3*/qend,plusp,genestrand);
  debug8(printf("trimpos %d (relative to %d)\n",*trimpos,/*qstart*/0));

  if (*trimpos >= qend) {
    debug8(printf("trimpos %d >= qend %d, so returning -1\n",*trimpos,qend));
    *trimpos = -1;
    return false;

  } else if (novelsplicingp == false || *trimpos == 0) {
    /* Keep given trim */
    *ambig_prob_qstart = 0.0;
    return false;

  } else if ((nspliceends = Splice_trim_novel_spliceends_5(&(*splicetype),&ambig_qstarts,&ambig_probs_5,
							   left,/*qstart*/(*trimpos),qend,
							   mismatch_positions_alloc,nmismatches,
							   chroffset,plusp,sensedir)) == 0) {
    /* Keep given trim */
    *ambig_prob_qstart = 0.0;
    return false;

  } else {
    /* TODO: Make sure this is the farthest one */
    *trimpos = ambig_qstarts[0];
    *ambig_prob_qstart = ambig_probs_5[0];
    FREE(ambig_qstarts);
    FREE(ambig_probs_5);
    return true;
  }
}


/* Returns true if splice was found, false otherwise.  If true, result
   is the number of spliceends.  If false, result is the single trim
   position */
bool
Substring_trimmed_qends (int *result, Splicetype_T *splicetype, int **ambig_qends, double **ambig_probs_3,
			 Univcoord_T left, int qstart, int querylength,
			 bool plusp, int genestrand, int *mismatch_positions_alloc, Compress_T query_compress,
			 Univcoord_T chroffset, int sensedir) {
  int trimpos;
  int nspliceends;
  int nmismatches;

  debug8(printf("\n***Entered Substring_trimmed_qends with left %u, genome %u..%u, qstart %d..%d, plusp %d, sensedir %d\n",
		left-chroffset,left+qstart-chroffset,left+querylength-chroffset,qstart,querylength,plusp,sensedir));

  trimpos = Substring_trim_qend_nosplice(&nmismatches,mismatch_positions_alloc,query_compress,left,
					 /*pos5*/qstart,/*pos3*/querylength,plusp,genestrand);
  debug8(printf("trimpos %d (relative to %d)\n",trimpos,/*qend*/querylength));

  if (trimpos <= qstart) {
    debug8(printf("trimpos %d <= qstart %d, so returning -1\n",trimpos,qstart));
    *result = -1;
    return false;

  } else if (novelsplicingp == false || trimpos == querylength) {
    *result = trimpos;
    return false;

  } else if ((nspliceends = Splice_trim_novel_spliceends_3(&(*splicetype),&(*ambig_qends),&(*ambig_probs_3),
							   left,qstart,/*qend*/trimpos,querylength,
							   mismatch_positions_alloc,nmismatches,
							   chroffset,plusp,sensedir)) == 0) {
    *result = trimpos;
    return false;

  } else {
    *result = nspliceends;
    return true;
  }
}

bool
Substring_qend_trim (int *trimpos, Splicetype_T *splicetype, double *ambig_prob_qend,
		     Univcoord_T left, int qstart, int querylength,
		     bool plusp, int genestrand, int *mismatch_positions_alloc, Compress_T query_compress,
		     Univcoord_T chroffset, int sensedir) {
  int nspliceends;
  int nmismatches;
  int *ambig_qends;
  double *ambig_probs_3;

  debug8(printf("\n***Entered Substring_qend_trim with left %u, genome %u..%u, qstart %d..%d, plusp %d, sensedir %d\n",
		left-chroffset,left+qstart-chroffset,left+querylength-chroffset,qstart,querylength,plusp,sensedir));

  *trimpos = Substring_trim_qend_nosplice(&nmismatches,mismatch_positions_alloc,query_compress,left,
					  /*pos5*/qstart,/*pos3*/querylength,plusp,genestrand);
  debug8(printf("trimpos %d (relative to %d)\n",*trimpos,/*qend*/querylength));

  if (*trimpos <= qstart) {
    debug8(printf("trimpos %d <= qstart %d, so returning -1\n",*trimpos,qstart));
    *trimpos = -1;
    return false;

  } else if (novelsplicingp == false || *trimpos == querylength) {
    /* Keep given trim */
    *ambig_prob_qend = 0.0;
    return false;

  } else if ((nspliceends = Splice_trim_novel_spliceends_3(&(*splicetype),&ambig_qends,&ambig_probs_3,
							   left,qstart,/*qend*/*trimpos,querylength,
							   mismatch_positions_alloc,nmismatches,
							   chroffset,plusp,sensedir)) == 0) {
    /* Keep given trim */
    *ambig_prob_qend = 0.0;
    return false;

  } else {
    /* TODO: Make sure this is the farthest one */
    *trimpos = ambig_qends[0];
    *ambig_prob_qend = ambig_probs_3[0];
    FREE(ambig_qends);
    FREE(ambig_probs_3);
    return true;
  }
}



/* Modified from Substring_new to return nmatches_plus_spliced_trims only */
int
Substring_compute_nmatches (Univcoord_T left, int querystart, int queryend, int querylength,
			    bool plusp, int genestrand, Compress_T query_compress,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			    bool splice_querystart_p, bool splice_queryend_p, bool chrnum_fixed_p) {
  int nmatches_plus_spliced_trims, nmismatches_bothdiff, amb_length;
  Univcoord_T alignstart, alignend;
  int outofbounds_start, outofbounds_end;
  int trim_at_querystart, trim_at_queryend;

  
  debug2(printf("\n***Entered Substring_compute_nmatches with left %u, chr %d, genome %u..%u, query %d..%d, plusp %d, splice_querystart_p %d, splice_queryend_p %d\n",
		left,chrnum,left+querystart-chroffset,left+queryend-chroffset,querystart,queryend,plusp,splice_querystart_p,splice_queryend_p));

#if 0
  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
  } else {
    genomicstart = left + querylength;
    genomicend = left;
  }
#endif

  trim_at_querystart = querystart;
  trim_at_queryend = querylength - queryend;

  if (chrnum_fixed_p == true) {
    /* Chromosomal information determined by the middle diagonal */
    if (plusp) {
      alignstart = left + trim_at_querystart;
      alignend = left + (querylength - trim_at_queryend);
      
      if (alignstart >= chroffset) {
	outofbounds_start = 0;
      } else {
	outofbounds_start = chroffset - alignstart;
      }
      if (alignend <= chrhigh) {
	outofbounds_end = 0;
      } else {
	outofbounds_end = alignend - chrhigh;
      }

    } else {
      alignstart = left + (querylength - trim_at_querystart);
      alignend = left + trim_at_queryend;

      if (alignend >= chroffset) {
	outofbounds_end = 0;
      } else {
	outofbounds_end = chroffset - alignend;
      }
      if (alignstart <= chrhigh) {
	outofbounds_start = 0;
      } else {
	outofbounds_start = alignstart - chrhigh;
      }
    }

  } else {
    /* Determine chromosomal information */
    if (plusp) {
      alignstart = left + trim_at_querystart;
      alignend = left + (querylength - trim_at_queryend);

      if (alignend < chroffset) {
	/* Need to recompute chromosome bounds (unexpected since chroffset should be based on left) */
	debug2(printf("Plus: recomputing chromosomal bounds\n"));
	chrnum = Univ_IIT_get_one(chromosome_iit,alignstart,alignstart);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

      } else if (alignstart >= chrhigh) {
	/* Need to recompute chromosome bounds */
	debug2(printf("Plus: recomputing chromosomal bounds\n"));
	chrnum = Univ_IIT_get_one(chromosome_iit,alignstart,alignstart);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
      }

      if (alignend <= chrhigh) {
	/* Alignment is within the chromosome */
	debug2(printf("alignend <= chrhigh, so alignment is within the chromosome\n"));
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
	  return -1;
	} else {
	  /* Move to next chromosome */
	  debug2(printf("Moving to next chromosome\n"));
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  if (alignend <= chrhigh) {
	    /* Alignment is within the new chromosome */
	    debug2(printf("alignend <= new chrhigh %u, so alignment is within the new chromosome\n",chrhigh));
	    outofbounds_end = 0;
	  } else {
	    outofbounds_end = alignend - chrhigh;
	  }
	}
      }

    } else {
      alignstart = left + (querylength - trim_at_querystart);
      alignend = left + trim_at_queryend;

      if (alignstart < chroffset) {
	/* Need to recompute chromosome bounds (unexpected since chroffset should be based on left) */
	debug2(printf("Minus: Recomputing chromosomal bounds\n"));
	chrnum = Univ_IIT_get_one(chromosome_iit,alignend,alignend);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

      } else if (alignend >= chrhigh) {
	/* Need to recompute chromosome bounds */
	debug2(printf("Minus: Recomputing chromosomal bounds\n"));
	chrnum = Univ_IIT_get_one(chromosome_iit,alignend,alignend);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
      }

      if (alignstart <= chrhigh) {
	/* Alignment is within the chromosome */
	debug2(printf("alignstart <= chrhigh, so alignment is within the chromosome\n"));
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
	  return -1;
	} else {
	  /* Move to next chromosome */
	  debug2(printf("Moving to next chromosome\n"));
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  if (alignstart <= chrhigh) {
	    debug2(printf("alignstart <= new chrhigh %u, so alignment is within the new chromosome\n",chrhigh));
	    outofbounds_start = 0;
	  } else {
	    outofbounds_start = alignstart - chrhigh;
	  }
	}
      }
    }
  }

  /* outofbounds values are relative to preliminary alignstart and alignend */
  trim_at_querystart += outofbounds_start;
  trim_at_queryend += outofbounds_end;

  querystart = trim_at_querystart;
  queryend = querylength - trim_at_queryend;
  if (querystart >= queryend) {
    /* Can happen if multiple mismatches result in trimming down to 0 */
    return -1;

  } else if (plusp == true) {
    debug2(printf("Counting mismatches from querystart %d to queryend %d\n",querystart,queryend));
    nmismatches_bothdiff =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
					/*pos5*/querystart,/*pos3*/queryend,/*plusp*/true,genestrand);
  } else {
    debug2(printf("Counting mismatches from querystart %d to queryend %d\n",querystart,queryend));
    nmismatches_bothdiff =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
					/*pos5*/querylength - queryend,/*pos3*/querylength - querystart,
					/*plusp*/false,genestrand);
  }

  /* Needs to be identical to code in Substring_nmatches.  Don't need to handle alts_substring, though. */
  amb_length = 0;
  if (splice_querystart_p == true) {
    amb_length += trim_at_querystart;
  }
  if (splice_queryend_p == true) {
    amb_length += trim_at_queryend;
  }
  nmatches_plus_spliced_trims = amb_length + (queryend - querystart) - nmismatches_bothdiff;
  
  debug2(printf("Substring_compute_nmatches returning %d matches\n",nmatches_plus_spliced_trims));
  return nmatches_plus_spliced_trims;
}


/* Want querylength and not querylength_adj */
/* If orig_nmismatches < 0, then this procedure (re-)computes the value */

/* For trimming: querystart and queryend should be used as a basis for
   trim_novel_spliceends (which search to both sides of the point, but
   not for trim_querystart and trim_queryend, which search only
   inwards.  For those procedures, should use querystart of 0 and
   queryend of querylength, and the results are absolute, not relative
   to the given querystart and queryend */

/* Assumes that chrnum is determined correctly, even for straddles, from middle diagonal */
T
Substring_new (int nmismatches, Univcoord_T left, int querystart, int queryend, int querylength,
	       bool plusp, int genestrand, Compress_T query_compress,
	       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
	       bool splice_querystart_p, Splicetype_T splicetype_querystart, double ambig_prob_querystart, 
	       bool splice_queryend_p, Splicetype_T splicetype_queryend, double ambig_prob_queryend,
	       int sensedir) {
  T new;
  int trim_querystart, trim_queryend;
  Univcoord_T alignstart, alignend;
  int outofbounds_start, outofbounds_end;

  
  debug2(printf("\n***Entered Substring_new with left %u, chr %d, genome %d..%d, query %d..%d, plusp %d, sensedir %d, splice_querystart_p %d, splice_queryend_p %d\n",
		left,chrnum,(int) (left+querystart-chroffset),(int) (left+queryend-chroffset),querystart,queryend,plusp,sensedir,splice_querystart_p,splice_queryend_p));

  /* assert(queryend > querystart); -- Assertion does not hold.  Sometimes queryend == querystart */

  new = (T) MALLOC_OUT(sizeof(*new));
  debug2(printf("substring %p:\n",new));

  /* Required values if we abort and invoke Substring_free early */
  new->alts_ncoords = 0;
  new->genomic_bothdiff = (char *) NULL;
  new->genomic_refdiff = (char *) NULL;

  new->querylength = querylength;

  /* These values are required by substring_trim_novel_spliceends */
  new->left = left;
  new->plusp = plusp;
  new->genestrand = genestrand;
  new->query_compress = query_compress;

  if (plusp == true) {
    new->genomicstart = left;
    new->genomicend = left + querylength;
  } else {
    new->genomicstart = left + querylength;
    new->genomicend = left;
  }

  new->trim_querystart_splicep = splice_querystart_p;
  new->trim_queryend_splicep = splice_queryend_p;

  new->querystart_pretrim = querystart;
  new->queryend_pretrim = queryend;

  trim_querystart = querystart;
  trim_queryend = querylength - queryend;
      
  new->sensedir = sensedir;

  new->amb_splice_pos = 0;
  new->splicecoord_D = new->splicecoord_A = 0;
  new->siteD_pos = new->siteA_pos = new->siteN_pos = 0;
  new->siteD_prob = new->siteA_prob = 0.0;

  if (splice_querystart_p == false) {
    new->start_amb_length = 0;
    new->amb_type = END;
    new->start_endtype = END;

  } else if (splicetype_querystart == DONOR || splicetype_querystart == ANTIDONOR) {
    new->start_amb_length = querystart;
    new->amb_type = DON;
    new->start_endtype = DON;
    new->siteD_prob = ambig_prob_querystart;

  } else {
    new->start_amb_length = querystart;
    new->amb_type = ACC;
    new->start_endtype = ACC;
    new->siteA_prob = ambig_prob_querystart;
  }

  if (splice_queryend_p == false) {
    new->end_amb_length = 0;
    new->amb_type = END;
    new->end_endtype = END;

  } else if (splicetype_queryend == DONOR || splicetype_queryend == ANTIDONOR) {
    new->end_amb_length = querylength - queryend;
    new->amb_type = DON;
    new->end_endtype = DON;
    new->siteD_prob = ambig_prob_queryend;

  } else {
    new->end_amb_length = querylength - queryend;
    new->amb_type = ACC;
    new->end_endtype = ACC;
    new->siteA_prob = ambig_prob_queryend;
  }


  /* Chromosomal information determined by the middle diagonal */
  if (plusp) {
    alignstart = left + trim_querystart;
    alignend = left + (querylength - trim_queryend);
    
#if 1
    /* Needed to protect against substrings on the wrong chromosome */
    if (alignend < chroffset || alignstart >= chrhigh) {
      debug2(printf("Substring fails because alignend %u < chroffset %u or alignstart %u >= chrhigh %u\n",
		    alignend,chroffset,alignstart,chrhigh));
      Substring_free(&new);
      return (T) NULL;
    }
#endif

    if (alignstart >= chroffset) {
      outofbounds_start = 0;
    } else {
      outofbounds_start = chroffset - alignstart;
    }
    if (alignend <= chrhigh) {
      outofbounds_end = 0;
    } else {
      outofbounds_end = alignend - chrhigh;
    }
    
  } else {
    alignstart = left + (querylength - trim_querystart);
    alignend = left + trim_queryend;
    
#if 1
    /* Needed to protect against substrings on the wrong chromosome */
    if (alignstart < chroffset || alignend >= chrhigh) {
      debug2(printf("Substring fails because alignstart %u < chroffset %u or alignend %u >= chrhigh %u\n",
		    alignstart,chroffset,alignend,chrhigh));
      Substring_free(&new);
      return (T) NULL;
    }
#endif
    
    if (alignend >= chroffset) {
      outofbounds_end = 0;
    } else {
      outofbounds_end = chroffset - alignend;
    }
    if (alignstart <= chrhigh) {
      outofbounds_start = 0;
    } else {
      outofbounds_start = alignstart - chrhigh;
    }
  }

  /* outofbounds values are relative to preliminary alignstart and alignend */
  trim_querystart += outofbounds_start;
  trim_queryend += outofbounds_end;

#if 0
  if (trim_querystart + trim_queryend >= querylength) {
    debug2(printf("Substring fails because trim_querystart %d + trim_queryend %d >= querylength %d\n",
		  trim_querystart,trim_queryend,querylength));
    Substring_free(&new);
    return (T) NULL;
  }
#endif

#ifdef DEBUG2
  printf("Outofbounds start %d, outofbounds end %d\n",outofbounds_start,outofbounds_end);
  if (outofbounds_start > 0) {
    printf("Out of bounds: Revising trim_querystart to be %d\n",trim_querystart);
  }
  if (outofbounds_end > 0) {
    printf("Out of bounds: Revising trim_queryend to be %d\n",trim_queryend);
  }
  printf("Got trims of %d and %d\n",trim_querystart,trim_queryend);
#endif


  debug2(printf("\n***chrnum %d (chroffset %u, chrhigh %u), plusp %d\n",chrnum,chroffset,chrhigh,plusp));
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;

  /* mandatory_trim_querystart and mandatory_trim_queryend also set during aliasing for circular chromosomes */
  if (left < new->chroffset) {
    new->mandatory_trim_querystart = new->chroffset - left;
    assert(new->mandatory_trim_querystart >= 0);
  } else {
    new->mandatory_trim_querystart = 0;
  }
  if (left + querylength >= new->chrhigh) {
    new->mandatory_trim_queryend = (left + querylength) - new->chrhigh;
    assert(new->mandatory_trim_queryend >= 0);
  } else {
    new->mandatory_trim_queryend = 0;
  }
  debug2(printf("mandatory_trim_querystart %d, mandatory_trim_queryend %d\n",new->mandatory_trim_querystart,new->mandatory_trim_queryend));


  if (querystart != trim_querystart || queryend != querylength - trim_queryend) {
    nmismatches = -1;		/* Need to recalculate, because of change in querystart */
  }
  new->querystart = trim_querystart;
  new->queryend = querylength - trim_queryend;
  debug2(printf("querystart %d, queryend %d\n",new->querystart,new->queryend));
     
  assert(new->querystart >= 0);
  assert(new->querystart <= querylength);
  assert(new->queryend >= 0);
  assert(new->queryend <= querylength);

  /* Compute coordinates */
  if (new->queryend <= new->querystart) {
    /* Can happen if multiple mismatches result in trimming down to 0 */
    Substring_free(&new);
    return (T) NULL;

  } else if (plusp == true) {
    new->alignstart_trim = left + new->querystart;
    new->alignend_trim = left + new->queryend;
    
    debug2(printf("left is %u\n",new->left));
    debug2(printf("genomicstart is %u, genomicend is %u\n",new->genomicstart,new->genomicend));
    debug2(printf("querylength is %d, alignstart is %u, alignend is %u\n",querylength,new->alignstart_trim,new->alignend_trim));
    assert(new->alignstart_trim >= chroffset);
    assert(new->alignend_trim <= chrhigh);
    
    debug2(printf("Counting mismatches from querystart %d to queryend %d\n",new->querystart,new->queryend));
    if (nmismatches >= 0) {
      debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",left - chroffset,new->querystart,new->queryend));
      debug7(printf("%d vs %d\n",nmismatches,
		    Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
						      /*pos5*/new->querystart,/*pos3*/new->queryend,
						      /*plusp*/true,genestrand)));
#ifdef CHECK_NMISMATCHES
      assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							      /*pos5*/new->querystart,/*pos3*/new->queryend,
							      /*plusp*/true,genestrand));
#endif
    } else {
      nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
						      /*pos5*/new->querystart,/*pos3*/new->queryend,
						      /*plusp*/true,genestrand);
    }

  } else {
    new->alignstart_trim = left + (querylength - new->querystart);
    new->alignend_trim = left + (querylength - new->queryend);
    
    debug2(printf("left is %u\n",new->left));
    debug2(printf("genomicstart is %u, genomicend is %u\n",new->genomicstart,new->genomicend));
    debug2(printf("querylength is %d, alignstart is %u, alignend is %u\n",querylength,new->alignstart_trim,new->alignend_trim));
    assert(new->alignstart_trim <= chrhigh);
    assert(new->alignend_trim >= chroffset);

    if (nmismatches >= 0) {
      debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",left - chroffset,new->querystart,new->queryend));
      debug7(printf("%d vs %d\n",nmismatches,
		    Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
						      /*pos5*/querylength - new->queryend,/*pos3*/querylength - new->querystart,
						      /*plusp*/false,genestrand)));
#ifdef CHECK_NMISMATCHES
      assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							      /*pos5*/querylength - new->queryend,/*pos3*/querylength - new->querystart,
							      /*plusp*/false,genestrand));
#endif
    } else {
      nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
						      /*pos5*/querylength - new->queryend,/*pos3*/querylength - new->querystart,
						      /*plusp*/false,genestrand);
    }
  }
  
  new->nmismatches_bothdiff = nmismatches;
  new->nmismatches_refdiff = new->nmismatches_bothdiff; /* Will be recalculated later */

  new->nmatches_to_trims = (new->queryend - new->querystart) - new->nmismatches_bothdiff;
  new->nmatches_plus_spliced_trims = new->nmatches_to_trims + new->start_amb_length + new->end_amb_length;

  debug2(printf("nmatches_to_trims %d = queryend %d - querystart %d - nmismatches_bothdiff %d\n",
		new->nmatches_to_trims,new->queryend,new->querystart,new->nmismatches_bothdiff));
  assert(new->nmatches_to_trims >= 0);

  new->alts_ncoords = 0;
  new->alts_coords = (Univcoord_T *) NULL;
  new->alts_knowni = (int *) NULL;
  new->alts_nmismatches = (int *) NULL;
  new->alts_probs = (double *) NULL;

  debug2(printf("Substring_new returning %d matches\n",new->nmatches_plus_spliced_trims));
#if 0
  debug2(printf("** Returning substring %p, query %d..%d, trim %d..%d, nmatches_plus_spliced_trims %d, nmismatches_refdiff %d, nmismatches_bothdiff %d, amb_lengths %d and %d\n",
		new,new->querystart,new->queryend,trim_querystart,trim_queryend,nmatches_plus_spliced_trims,
		new->nmismatches_refdiff,new->nmismatches_bothdiff,Substring_start_amb_length(new),Substring_end_amb_length(new)));
#endif

  assert(new->nmismatches_bothdiff >= 0);
  return new;
}



#if 0
T
Substring_new_compute_chrnum (int nmismatches, Univcoord_T left, int querystart, int queryend, int querylength,
			      bool plusp, int genestrand, Compress_T query_compress,
			      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			      bool splice_querystart_p, Splicetype_T splicetype_querystart, double ambig_prob_querystart, 
			      bool splice_queryend_p, Splicetype_T splicetype_queryend, double ambig_prob_queryend,
			      int sensedir, bool chrnum_fixed_p) {
  T new;
  int trim_querystart, trim_queryend;

  Univcoord_T alignstart, alignend;
  int outofbounds_start, outofbounds_end;

  
  debug2(printf("\n***Entered Substring_new with left %u, chr %d, genome %d..%d, query %d..%d, plusp %d, sensedir %d, splice_querystart_p %d, splice_queryend_p %d\n",
		left,chrnum,(int) (left+querystart-chroffset),(int) (left+queryend-chroffset),querystart,queryend,plusp,sensedir,splice_querystart_p,splice_queryend_p));

  /* assert(queryend > querystart); -- Assertion does not hold.  Sometimes queryend == querystart */

  new = (T) MALLOC_OUT(sizeof(*new));
  debug2(printf("substring %p:\n",new));

  /* Required values if we abort and invoke Substring_free early */
  new->alts_ncoords = 0;
  new->genomic_bothdiff = (char *) NULL;
  new->genomic_refdiff = (char *) NULL;

  new->querylength = querylength;

  /* These values are required by substring_trim_novel_spliceends */
  new->left = left;
  new->plusp = plusp;
  new->genestrand = genestrand;
  new->query_compress = query_compress;

  if (plusp == true) {
    new->genomicstart = left;
    new->genomicend = left + querylength;
  } else {
    new->genomicstart = left + querylength;
    new->genomicend = left;
  }

  new->trim_querystart_splicep = splice_querystart_p;
  new->trim_queryend_splicep = splice_queryend_p;

  new->querystart_pretrim = querystart;
  new->queryend_pretrim = queryend;

  trim_querystart = querystart;
  trim_queryend = querylength - queryend;
      
  new->sensedir = sensedir;

  new->amb_splice_pos = 0;
  new->splicecoord_D = new->splicecoord_A = 0;
  new->siteD_pos = new->siteA_pos = new->siteN_pos = 0;
  new->siteD_prob = new->siteA_prob = 0.0;

  if (splice_querystart_p == false) {
    new->start_amb_length = 0;
    new->amb_type = END;
    new->start_endtype = END;

  } else if (splicetype_querystart == DONOR || splicetype_querystart == ANTIDONOR) {
    new->start_amb_length = querystart;
    new->amb_type = DON;
    new->start_endtype = DON;
    new->siteD_prob = ambig_prob_querystart;

  } else {
    new->start_amb_length = querystart;
    new->amb_type = ACC;
    new->start_endtype = ACC;
    new->siteA_prob = ambig_prob_querystart;
  }

  if (splice_queryend_p == false) {
    new->end_amb_length = 0;
    new->amb_type = END;
    new->end_endtype = END;

  } else if (splicetype_queryend == DONOR || splicetype_queryend == ANTIDONOR) {
    new->end_amb_length = querylength - queryend;
    new->amb_type = DON;
    new->end_endtype = DON;
    new->siteD_prob = ambig_prob_queryend;

  } else {
    new->end_amb_length = querylength - queryend;
    new->amb_type = ACC;
    new->end_endtype = ACC;
    new->siteA_prob = ambig_prob_queryend;
  }


  if (chrnum_fixed_p == true) {
    /* Chromosomal information determined by the middle diagonal */
    if (plusp) {
      alignstart = left + trim_querystart;
      alignend = left + (querylength - trim_queryend);
      
#if 0
      if (alignend < chroffset || alignstart >= chrhigh) {
	debug2(printf("Substring fails because alignend %u < chroffset %u or alignstart %u >= chrhigh %u\n",
		      alignend,chroffset,alignstart,chrhigh));
	Substring_free(&new);
	return (T) NULL;
      }
#endif

      if (alignstart >= chroffset) {
	outofbounds_start = 0;
      } else {
	outofbounds_start = chroffset - alignstart;
      }
      if (alignend <= chrhigh) {
	outofbounds_end = 0;
      } else {
	outofbounds_end = alignend - chrhigh;
      }

    } else {
      alignstart = left + (querylength - trim_querystart);
      alignend = left + trim_queryend;

#if 0
      if (alignstart < chroffset || alignend >= chrhigh) {
	debug2(printf("Substring fails because alignstart %u < chroffset %u or alignend %u >= chrhigh %u\n",
		      alignstart,chroffset,alignend,chrhigh));
	Substring_free(&new);
	return (T) NULL;
      }
#endif

      if (alignend >= chroffset) {
	outofbounds_end = 0;
      } else {
	outofbounds_end = chroffset - alignend;
      }
      if (alignstart <= chrhigh) {
	outofbounds_start = 0;
      } else {
	outofbounds_start = alignstart - chrhigh;
      }
    }

  } else {
    /* Determine chromosomal information */
    if (plusp) {
      alignstart = left + trim_querystart;
      alignend = left + (querylength - trim_queryend);

      if (alignend < chroffset) {
	/* Need to recompute chromosome bounds (unexpected since chroffset should be based on left) */
	debug2(printf("Plus: recomputing chromosomal bounds\n"));
	chrnum = Univ_IIT_get_one(chromosome_iit,alignstart,alignstart);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

      } else if (alignstart >= chrhigh) {
	/* Need to recompute chromosome bounds */
	debug2(printf("Plus: recomputing chromosomal bounds\n"));
	chrnum = Univ_IIT_get_one(chromosome_iit,alignstart,alignstart);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
      }

      if (alignend <= chrhigh) {
	/* Alignment is within the chromosome */
	debug2(printf("alignend <= chrhigh, so alignment is within the chromosome\n"));
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
	  Substring_free(&new);
	  return (T) NULL;
	} else {
	  /* Move to next chromosome */
	  debug2(printf("Moving to next chromosome\n"));
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  if (alignend <= chrhigh) {
	    /* Alignment is within the new chromosome */
	    debug2(printf("alignend <= new chrhigh %u, so alignment is within the new chromosome\n",chrhigh));
	    outofbounds_end = 0;
	  } else {
	    outofbounds_end = alignend - chrhigh;
	  }
	}
      }

    } else {
      alignstart = left + (querylength - trim_querystart);
      alignend = left + trim_queryend;

      if (alignstart < chroffset) {
	/* Need to recompute chromosome bounds (unexpected since chroffset should be based on left) */
	debug2(printf("Minus: Recomputing chromosomal bounds\n"));
	chrnum = Univ_IIT_get_one(chromosome_iit,alignend,alignend);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

      } else if (alignend >= chrhigh) {
	/* Need to recompute chromosome bounds */
	debug2(printf("Minus: Recomputing chromosomal bounds\n"));
	chrnum = Univ_IIT_get_one(chromosome_iit,alignend,alignend);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
      }

      if (alignstart <= chrhigh) {
	/* Alignment is within the chromosome */
	debug2(printf("alignstart <= chrhigh, so alignment is within the chromosome\n"));
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
	  Substring_free(&new);
	  return (T) NULL;
	} else {
	  /* Move to next chromosome */
	  debug2(printf("Moving to next chromosome\n"));
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  if (alignstart <= chrhigh) {
	    debug2(printf("alignstart <= new chrhigh %u, so alignment is within the new chromosome\n",chrhigh));
	    outofbounds_start = 0;
	  } else {
	    outofbounds_start = alignstart - chrhigh;
	  }
	}
      }
    }
  }

  /* outofbounds values are relative to preliminary alignstart and alignend */
  trim_querystart += outofbounds_start;
  trim_queryend += outofbounds_end;

#if 0
  if (trim_querystart + trim_queryend >= querylength) {
    debug2(printf("Substring fails because trim_querystart %d + trim_queryend %d >= querylength %d\n",
		  trim_querystart,trim_queryend,querylength));
    Substring_free(&new);
    return (T) NULL;
  }
#endif

#ifdef DEBUG2
  printf("Outofbounds start %d, outofbounds end %d\n",outofbounds_start,outofbounds_end);
  if (outofbounds_start > 0) {
    printf("Out of bounds: Revising trim_querystart to be %d\n",trim_querystart);
  }
  if (outofbounds_end > 0) {
    printf("Out of bounds: Revising trim_queryend to be %d\n",trim_queryend);
  }
  printf("Got trims of %d and %d\n",trim_querystart,trim_queryend);
#endif


  debug2(printf("\n***chrnum %d (chroffset %u, chrhigh %u), plusp %d\n",chrnum,chroffset,chrhigh,plusp));
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;

  /* mandatory_trim_querystart and mandatory_trim_queryend also set during aliasing for circular chromosomes */
  if (left < new->chroffset) {
    new->mandatory_trim_querystart = new->chroffset - left;
    assert(new->mandatory_trim_querystart >= 0);
  } else {
    new->mandatory_trim_querystart = 0;
  }
  if (left + querylength >= new->chrhigh) {
    new->mandatory_trim_queryend = (left + querylength) - new->chrhigh;
    assert(new->mandatory_trim_queryend >= 0);
  } else {
    new->mandatory_trim_queryend = 0;
  }
  debug2(printf("mandatory_trim_querystart %d, mandatory_trim_queryend %d\n",new->mandatory_trim_querystart,new->mandatory_trim_queryend));


  if (querystart != trim_querystart || queryend != querylength - trim_queryend) {
    nmismatches = -1;		/* Need to recalculate, because of change in querystart */
  }
  new->querystart = trim_querystart;
  new->queryend = querylength - trim_queryend;
  debug2(printf("querystart %d, queryend %d\n",new->querystart,new->queryend));
     
  assert(new->querystart >= 0);
  assert(new->querystart <= querylength);
  assert(new->queryend >= 0);
  assert(new->queryend <= querylength);

  /* Compute coordinates */
  if (new->queryend <= new->querystart) {
    /* Can happen if multiple mismatches result in trimming down to 0 */
    Substring_free(&new);
    return (T) NULL;

  } else if (plusp == true) {
    new->alignstart_trim = left + new->querystart;
    new->alignend_trim = left + new->queryend;
    
    debug2(printf("left is %u\n",new->left));
    debug2(printf("genomicstart is %u, genomicend is %u\n",new->genomicstart,new->genomicend));
    debug2(printf("querylength is %d, alignstart is %u, alignend is %u\n",querylength,new->alignstart_trim,new->alignend_trim));
    assert(new->alignstart_trim >= chroffset);
    assert(new->alignend_trim <= chrhigh);
    
    debug2(printf("Counting mismatches from querystart %d to queryend %d\n",new->querystart,new->queryend));
    if (nmismatches >= 0) {
      debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",left - chroffset,new->querystart,new->queryend));
      debug7(printf("%d vs %d\n",nmismatches,
		    Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
						      /*pos5*/new->querystart,/*pos3*/new->queryend,
						      /*plusp*/true,genestrand)));
#ifdef CHECK_NMISMATCHES
      assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							      /*pos5*/new->querystart,/*pos3*/new->queryend,
							      /*plusp*/true,genestrand));
#endif
    } else {
      nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
						      /*pos5*/new->querystart,/*pos3*/new->queryend,
						      /*plusp*/true,genestrand);
    }

  } else {
    new->alignstart_trim = left + (querylength - new->querystart);
    new->alignend_trim = left + (querylength - new->queryend);
    
    debug2(printf("left is %u\n",new->left));
    debug2(printf("genomicstart is %u, genomicend is %u\n",new->genomicstart,new->genomicend));
    debug2(printf("querylength is %d, alignstart is %u, alignend is %u\n",querylength,new->alignstart_trim,new->alignend_trim));
    assert(new->alignstart_trim <= chrhigh);
    assert(new->alignend_trim >= chroffset);

    if (nmismatches >= 0) {
      debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",left - chroffset,new->querystart,new->queryend));
      debug7(printf("%d vs %d\n",nmismatches,
		    Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
						      /*pos5*/querylength - new->queryend,
						      /*pos3*/querylength - new->querystart,
						      /*plusp*/false,genestrand)));
#ifdef CHECK_NMISMATCHES
      assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							      /*pos5*/querylength - new->queryend,
							      /*pos3*/querylength - new->querystart,
							      /*plusp*/false,genestrand));
#endif
    } else {
      nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
						      /*pos5*/querylength - new->queryend,
						      /*pos3*/querylength - new->querystart,
						      /*plusp*/false,genestrand);
    }
  }
  
  new->nmismatches_bothdiff = nmismatches;
  new->nmismatches_refdiff = new->nmismatches_bothdiff; /* Will be recalculated later */

  new->nmatches_to_trims = (new->queryend - new->querystart) - new->nmismatches_bothdiff;
  new->nmatches_plus_spliced_trims = new->nmatches_to_trims + new->start_amb_length + new->end_amb_length;

  debug2(printf("nmatches_plus_spliced_trims %d = queryend %d - querystart %d - nmismatches_bothdiff %d\n",
		new->nmatches_plus_spliced_trims,new->queryend,new->querystart,new->nmismatches_bothdiff));
  assert(new->nmatches_to_trims >= 0);

  new->alts_ncoords = 0;
  new->alts_coords = (Univcoord_T *) NULL;
  new->alts_knowni = (int *) NULL;
  new->alts_nmismatches = (int *) NULL;
  new->alts_probs = (double *) NULL;

  debug2(printf("Substring_new returning %d matches\n",new->nmatches_plus_spliced_trims));
#if 0
  debug2(printf("** Returning substring %p, query %d..%d, trim %d..%d, nmatches_plus_spliced_trims %d, nmismatches_refdiff %d, nmismatches_bothdiff %d, amb_lengths %d and %d\n",
		new,new->querystart,new->queryend,trim_querystart,trim_queryend,nmatches_plus_spliced_trims,
		new->nmismatches_refdiff,new->nmismatches_bothdiff,Substring_start_amb_length(new),Substring_end_amb_length(new)));
#endif

  assert(new->nmismatches_bothdiff >= 0);
  return new;
}
#endif


void
Substring_set_querystart_pretrim (T this, int querystart_pretrim) {
  this->querystart_pretrim = querystart_pretrim;
  return;
}

void
Substring_set_queryend_pretrim (T this, int queryend_pretrim) {
  this->queryend_pretrim = queryend_pretrim;
  return;
}


/* Used when converting pairs to substrings.  Assumes nmismatches and
   sensedir are correct, and no trimming is needed */
T
Substring_new_simple (int nmismatches, Univcoord_T left, int querystart, int queryend, int querylength,
		      bool plusp, int genestrand, Compress_T query_compress,
		      Endtype_T start_endtype, Endtype_T end_endtype,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      int sensedir) {
  T new;

  
  debug2(printf("\n***Entered Substring_new_simple with left %u, chr %d, genome %u..%u, query %d..%d, plusp %d\n",
		left,chrnum,left+querystart-chroffset,left+queryend-chroffset,querystart,queryend,plusp));

  /* assert(queryend > querystart); -- Assertion does not hold.  Sometimes queryend == querystart */
  new = (T) MALLOC_OUT(sizeof(*new));
  debug2(printf("substring %p:\n",new));


  /* Required values if we abort and invoke Substring_free early */
  new->alts_ncoords = 0;
  new->genomic_bothdiff = (char *) NULL;
  new->genomic_refdiff = (char *) NULL;

  new->start_endtype = start_endtype;
  new->end_endtype = end_endtype;

  /* Revised later for the first and last substrings */
  new->querystart_pretrim = querystart;
  new->queryend_pretrim = queryend;

  new->querylength = querylength;

  /* These values are required by substring_trim_novel_spliceends */
  new->left = left;
  new->plusp = plusp;
  new->genestrand = genestrand;
  new->query_compress = query_compress;

  if (plusp == true) {
    new->genomicstart = left;
    new->genomicend = left + querylength;
  } else {
    new->genomicstart = left + querylength;
    new->genomicend = left;
  }

  new->amb_splice_pos = 0;

  new->splicecoord_D = new->splicecoord_A = 0;
  new->siteD_pos = new->siteA_pos = new->siteN_pos = 0;

  new->siteD_prob = new->siteA_prob = 0.0;

  new->trim_querystart_splicep = new->trim_queryend_splicep = false;
  new->start_amb_length = new->end_amb_length = 0;
  new->sensedir = sensedir;

  /* Assume coordinates are all within the given chromosome */
  debug2(printf("\n***chrnum %d (chroffset %u, chrhigh %u), plusp %d\n",chrnum,chroffset,chrhigh,plusp));
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;

  /* mandatory_trim_querystart and mandatory_trim_queryend also set during aliasing for circular chromosomes */
  if (left < new->chroffset) {
    new->mandatory_trim_querystart = new->chroffset - left;
    assert(new->mandatory_trim_querystart >= 0);
  } else {
    new->mandatory_trim_querystart = 0;
  }
  if (left + querylength >= new->chrhigh) {
    new->mandatory_trim_queryend = (left + querylength) - new->chrhigh;
    assert(new->mandatory_trim_queryend >= 0);
  } else {
    new->mandatory_trim_queryend = 0;
  }
  debug2(printf("mandatory_trim_querystart %d, mandatory_trim_queryend %d\n",new->mandatory_trim_querystart,new->mandatory_trim_queryend));

  new->querystart = querystart;
  new->queryend = queryend;
  debug2(printf("querystart %d, queryend %d\n",new->querystart,new->queryend));
     
  assert(new->querystart >= 0);
  assert(new->querystart <= querylength);
  assert(new->queryend >= 0);
  assert(new->queryend <= querylength);

  /* Compute coordinates */
  if (plusp == true) {
    new->alignstart_trim = left + new->querystart;
    new->alignend_trim = left + new->queryend;

    debug2(printf("left is %u\n",new->left));
    debug2(printf("genomicstart is %u, genomicend is %u\n",new->genomicstart,new->genomicend));
    debug2(printf("querylength is %d, alignstart is %u, alignend is %u\n",querylength,new->alignstart_trim,new->alignend_trim));
    assert(new->alignstart_trim >= chroffset);
    assert(new->alignend_trim <= chrhigh);
    
  } else {
    new->alignstart_trim = left + (querylength - new->querystart);
    new->alignend_trim = left + (querylength - new->queryend);

    debug2(printf("left is %u\n",new->left));
    debug2(printf("genomicstart is %u, genomicend is %u\n",new->genomicstart,new->genomicend));
    debug2(printf("querylength is %d, alignstart is %u, alignend is %u\n",querylength,new->alignstart_trim,new->alignend_trim));
    assert(new->alignstart_trim <= chrhigh);
    assert(new->alignend_trim >= chroffset);
  }


  new->nmismatches_bothdiff = nmismatches;
  /* new->nmatches_plus_spliced_trims = (new->queryend - new->querystart) - new->nmismatches_bothdiff; */

  new->nmatches_plus_spliced_trims = new->nmatches_to_trims = (new->queryend - new->querystart) - new->nmismatches_bothdiff;

  new->alts_ncoords = 0;
  new->alts_coords = (Univcoord_T *) NULL;
  new->alts_knowni = (int *) NULL;
  new->alts_nmismatches = (int *) NULL;
  new->alts_probs = (double *) NULL;
  new->amb_type = END;

  debug2(printf("** Returning substring %p, query %d..%d, trim %d..%d, nmismatches %d\n",
		new,new->querystart,new->queryend,querystart,querylength - queryend,new->nmismatches_bothdiff));

  assert(new->nmismatches_bothdiff >= 0);

  if (new->querystart >= new->queryend) {
    /* Can happen if multiple mismatches result in trimming down to 0 */
    Substring_free(&new);
    return (T) NULL;
  } else {
    return new;
  }
}


static Univcoord_T
Univcoord_max (Univcoord_T *list, int n) {
  Univcoord_T m = list[0];
  int i;

  for (i = 1; i < n; i++) {
    if (list[i] > m) {
      m = list[i];
    }
  }

  return m;
}

static Univcoord_T
Univcoord_min (Univcoord_T *list, int n) {
  Univcoord_T m = list[0];
  int i;

  for (i = 1; i < n; i++) {
    if (list[i] < m) {
      m = list[i];
    }
  }

  return m;
}

static int
int_min (int *list, int n) {
  int m = list[0];
  int i;

  for (i = 1; i < n; i++) {
    if (list[i] < m) {
      m = list[i];
    }
  }

  return m;
}


T
Substring_new_alts_D (int querystart, int queryend, int splice_pos, int querylength,
		      bool plusp, int genestrand, Compress_T query_compress,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      Univcoord_T *alts_coords, int *alts_knowni, int *alts_nmismatches, double *alts_probs,
		      int alts_ncoords, double alts_common_prob, bool substring1p) {
  T new = (T) MALLOC_OUT(sizeof(*new));
#ifdef DEBUG2
  int i;
#endif

  debug2(printf("Entered Substring_new_alts_D with chrnum %d (chroffset %u, chrhigh %u), %d..%d, splice_pos %d, querylength %d, plusp %d\n",
		chrnum,chroffset,chrhigh,querystart,queryend,splice_pos,querylength,plusp));
  debug2(printf("alts_common_prob %f\n",alts_common_prob));

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;

  new->left = 0;
  if (plusp == true) {
    new->genomicstart = Univcoord_max(alts_coords,alts_ncoords);
    new->genomicend = Univcoord_min(alts_coords,alts_ncoords);
  } else {
    new->genomicstart = Univcoord_min(alts_coords,alts_ncoords);
    new->genomicend = Univcoord_max(alts_coords,alts_ncoords);
  }

  new->start_endtype = END;
  new->end_endtype = END;

  new->querystart = new->querystart_pretrim = querystart;
  new->queryend = new->queryend_pretrim = queryend;
  new->querylength = querylength;

  new->alignstart_trim = 0;
  new->alignend_trim = 0;

  new->plusp = plusp;
  new->genestrand = genestrand;
  new->query_compress = query_compress;

  new->siteD_prob = 0.0;
  new->siteA_prob = alts_common_prob;

  new->nmismatches_bothdiff = int_min(alts_nmismatches,alts_ncoords);
  debug2(printf("nmismatches_bothdiff due to alts_nmismatches is %d\n",new->nmismatches_bothdiff));

  /* Works because new->querystart == querystart and new->queryend == queryend */
  /* new->nmatches_plus_spliced_trims = (queryend - querystart) - new->nmismatches_bothdiff; */

  new->genomic_bothdiff = (char *) NULL;
  new->genomic_refdiff = (char *) NULL;

  if (substring1p == true) {
    debug2(printf("substring1p is true, so setting amb_lengths to be %d and %d\n",queryend,0));
    if (alts_common_prob >= 0.9) {
      /* Treat as a substring */
      new->nmatches_to_trims = new->nmatches_plus_spliced_trims = queryend - new->nmismatches_bothdiff;
      new->start_amb_length = new->end_amb_length = 0;
    } else {
      /* Treat as ambiguous splice site */
      new->nmatches_to_trims = new->nmatches_plus_spliced_trims = 0;
      new->start_amb_length = queryend - new->nmismatches_bothdiff; /* This is the amb length for the entire Stage3end_T object */
      new->end_amb_length = 0;
    }
  } else {
    debug2(printf("substring1p is false, so setting amb_lengths to be %d and %d\n",0,querylength - querystart));
    if (alts_common_prob >= 0.9) {
      /* Treat as a substring */
      new->nmatches_to_trims = new->nmatches_plus_spliced_trims = querylength - querystart - new->nmismatches_bothdiff;
      new->start_amb_length = new->end_amb_length = 0;
    } else {
      /* Treat as ambiguous splice site */
      new->nmatches_to_trims = new->nmatches_plus_spliced_trims = 0;
      new->start_amb_length = 0;
      new->end_amb_length = querylength - querystart - new->nmismatches_bothdiff; /* This is the amb length for the entire Stage3end_T object */
    }
  }
  new->mandatory_trim_querystart = 0;
  new->mandatory_trim_queryend = 0;
  new->trim_querystart_splicep = new->trim_queryend_splicep = false;
  /* new->start_amb_length = new->end_amb_length = 0; -- Set above */

#ifdef DEBUG2
  for (i = 0; i < alts_ncoords; i++) {
    printf("%u %f\n",alts_coords[i],alts_probs[i]);
  }
#endif
  new->alts_ncoords = alts_ncoords;
  new->alts_coords = alts_coords;
  new->alts_knowni = alts_knowni;
  new->alts_nmismatches = alts_nmismatches;
  new->alts_probs = alts_probs;
  new->amb_splice_pos = splice_pos;
  new->amb_type = DON;

  assert(new->nmismatches_bothdiff >= 0);

  if (new->querystart >= new->queryend) {
    /* Can happen if multiple mismatches result in trimming down to 0 */
    Substring_free(&new);
    return (T) NULL;
  } else {
    return new;
  }
}

T
Substring_new_alts_A (int querystart, int queryend, int splice_pos, int querylength,
		      bool plusp, int genestrand, Compress_T query_compress,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      Univcoord_T *alts_coords, int *alts_knowni, int *alts_nmismatches, double *alts_probs,
		      int alts_ncoords, double alts_common_prob, bool substring1p) {
  T new = (T) MALLOC_OUT(sizeof(*new));
#ifdef DEBUG2
  int i;
#endif
  
  debug2(printf("Entered Substring_new_alts_A with chrnum %d (chroffset %u, chrhigh %u), %d..%d, splice_pos %d, querylength %d, plusp %d\n",
		chrnum,chroffset,chrhigh,querystart,queryend,splice_pos,querylength,plusp));
  debug2(printf("alts_common_prob %f\n",alts_common_prob));

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;

  new->left = 0;
  if (plusp == true) {
    new->genomicstart = Univcoord_max(alts_coords,alts_ncoords);
    new->genomicend = Univcoord_min(alts_coords,alts_ncoords);
  } else {
    new->genomicstart = Univcoord_min(alts_coords,alts_ncoords);
    new->genomicend = Univcoord_max(alts_coords,alts_ncoords);
  }

  new->start_endtype = END;
  new->end_endtype = END;

  new->querystart = new->querystart_pretrim = querystart;
  new->queryend = new->queryend_pretrim = queryend;
  new->querylength = querylength;

  new->alignstart_trim = 0;
  new->alignend_trim = 0;

  new->plusp = plusp;
  new->genestrand = genestrand;
  new->query_compress = query_compress;

  new->siteA_prob = 0.0;
  new->siteD_prob = alts_common_prob;

  new->nmismatches_bothdiff = int_min(alts_nmismatches,alts_ncoords);
  debug2(printf("nmismatches_bothdiff due to alts_nmismatches is %d\n",new->nmismatches_bothdiff));

  /* Works because new->querystart == querystart and new->queryend == queryend */
  /* new->nmatches_plus_spliced_trims = (queryend - querystart) - new->nmismatches_bothdiff; */

  new->genomic_bothdiff = (char *) NULL;
  new->genomic_refdiff = (char *) NULL;

  if (substring1p == true) {
    debug2(printf("substring1p is true, so setting trims to be %d and %d\n",querystart,0));
    if (alts_common_prob >= 0.9) {
      /* Treat as a substring */
      new->nmatches_to_trims = new->nmatches_plus_spliced_trims = queryend - new->nmismatches_bothdiff;
      new->start_amb_length = new->end_amb_length = 0;
    } else {
      /* Treat as ambiguous splice site */
      new->nmatches_to_trims = new->nmatches_plus_spliced_trims = 0;
      new->start_amb_length = queryend - new->nmismatches_bothdiff; /* This is the amb length for the entire Stage3end_T object */
      new->end_amb_length = 0;
    }
  } else {
    debug2(printf("substring1p is false, so setting trims to be %d and %d\n",0,querylength - queryend));
    new->nmatches_to_trims = querylength - querystart - new->nmismatches_bothdiff;
    if (alts_common_prob >= 0.9) {
      /* Treat as a substring */
      new->nmatches_to_trims = new->nmatches_plus_spliced_trims = querylength - querystart - new->nmismatches_bothdiff;
      new->start_amb_length = new->end_amb_length = 0;
    } else {
      /* Treat as ambiguous splice site */
      new->nmatches_to_trims = new->nmatches_plus_spliced_trims = 0;
      new->start_amb_length = 0;
      new->end_amb_length = querylength - querystart - new->nmismatches_bothdiff; /* This is the amb length for the entire Stage3end_T object */
    }
  }
  new->mandatory_trim_querystart = 0;
  new->mandatory_trim_queryend = 0;
  new->trim_querystart_splicep = new->trim_queryend_splicep = false;

#ifdef DEBUG2
  for (i = 0; i < alts_ncoords; i++) {
    printf("%u %f\n",alts_coords[i],alts_probs[i]);
  }
#endif

  new->alts_ncoords = alts_ncoords;
  new->alts_coords = alts_coords;
  new->alts_knowni = alts_knowni;
  new->alts_nmismatches = alts_nmismatches;
  new->alts_probs = alts_probs;
  new->amb_splice_pos = splice_pos;
  new->amb_type = ACC;

  assert(new->nmismatches_bothdiff >= 0);

  if (new->querystart >= new->queryend) {
    /* Can happen if multiple mismatches result in trimming down to 0 */
    Substring_free(&new);
    return (T) NULL;
  } else {
    return new;
  }
}


Univcoord_T
Substring_set_alt (double *donor_prob, double *acceptor_prob, Univcoord_T *genomicstart, Univcoord_T *genomicend,
		   T this, int bingoi) {

#ifdef DEBUG10
  printf("Entered Substring_set_alt for %d..%d.  plusp %d, ",this->querystart,this->queryend,this->plusp);
  if (this->amb_type == DON) {
    printf("type DON\n");
  } else {
    printf("type ACC\n");
  }
#endif

  assert(this->alts_probs != NULL);
  
  this->nmismatches_refdiff = this->nmismatches_bothdiff = this->alts_nmismatches[bingoi];

  if (this->plusp == true) {
    if (this->amb_type == DON) {
      *acceptor_prob = this->siteA_prob;
      *donor_prob = this->siteD_prob = this->alts_probs[bingoi];
      this->splicecoord_D = this->alts_coords[bingoi];
      this->splicesitesD_knowni = this->alts_knowni[bingoi];
      this->left = this->splicecoord_D - this->amb_splice_pos;
    } else {
      *donor_prob = this->siteD_prob;
      *acceptor_prob = this->siteA_prob = this->alts_probs[bingoi];
      this->splicecoord_A = this->alts_coords[bingoi];
      this->splicesitesA_knowni = this->alts_knowni[bingoi];
      this->left = this->splicecoord_A - this->amb_splice_pos;
    }

    debug10(printf("left %u [%u]\n",this->left,this->left - this->chroffset));
    *genomicstart = this->genomicstart = this->left;
    *genomicend = this->genomicend = this->left + this->querylength;
    debug10(printf("genomicstart %u [%u]\n",*genomicstart,*genomicstart - this->chroffset));
    debug10(printf("genomicend %u [%u]\n",*genomicend,*genomicend - this->chroffset));

    this->alignstart_trim = this->genomicstart + this->querystart;
    this->alignend_trim =  this->genomicstart + this->queryend;
    /* this->nmatches_plus_spliced_trims = (this->alignend_trim - this->alignstart_trim) - this->nmismatches_bothdiff; */

    debug10(printf("querypos %d..%d, alignstart is %u (%u), alignend is %u (%u), genomicstart is %u, genomicend is %u\n",
		  this->querystart,this->queryend,this->alignstart_trim,this->alignstart_trim - this->chroffset,
		  this->alignend_trim,this->alignend_trim - this->chroffset,this->genomicstart,this->genomicend));

  } else {
    if (this->amb_type == DON) {
      *acceptor_prob = this->siteA_prob;
      *donor_prob = this->siteD_prob = this->alts_probs[bingoi];
      this->splicecoord_D = this->alts_coords[bingoi];
      this->splicesitesD_knowni = this->alts_knowni[bingoi];
      /* this->left = this->splicecoord_D - (this->querylength - this->amb_splice_pos); */
      this->left = this->splicecoord_D - this->amb_splice_pos; /* amb_splice_pos is based on left */
      debug10(printf("left %u [%u] = %u - (%d - %d)\n",
		     this->left,this->left - this->chroffset,this->splicecoord_D - this->chroffset,
		     this->querylength,this->amb_splice_pos));
    } else {
      *donor_prob = this->siteD_prob;
      *acceptor_prob = this->siteA_prob = this->alts_probs[bingoi];
      this->splicecoord_A = this->alts_coords[bingoi];
      this->splicesitesA_knowni = this->alts_knowni[bingoi];
      /* this->left = this->splicecoord_A - (this->querylength - this->amb_splice_pos); */
      this->left = this->splicecoord_A - this->amb_splice_pos; /* amb_splice_pos is based on left */
      debug10(printf("left %u [%u] = %u - (%d - %d)\n",
		     this->left,this->left - this->chroffset,this->splicecoord_A - this->chroffset,
		     this->querylength,this->amb_splice_pos));
    }

    *genomicend = this->genomicend = this->left;
    *genomicstart = this->genomicstart = this->left + this->querylength;
    debug10(printf("genomicstart %u [%u]\n",*genomicstart,*genomicstart - this->chroffset));
    debug10(printf("genomicend %u [%u]\n",*genomicend,*genomicend - this->chroffset));

    this->alignend_trim = this->genomicstart - this->queryend;
    this->alignstart_trim = this->genomicstart - this->querystart;
    /* this->nmatches_plus_spliced_trims = (this->alignstart_trim - this->alignend_trim) - this->nmismatches_bothdiff; */

    debug10(printf("querypos %d..%d, alignstart is %u (%u), alignend is %u (%u), genomicstart is %u, genomicend is %u\n",
		  this->querystart,this->queryend,this->alignstart_trim,this->alignstart_trim - this->chroffset,
		  this->alignend_trim,this->alignend_trim - this->chroffset,this->genomicstart,this->genomicend));
  }

  this->start_amb_length = this->end_amb_length = 0;
  this->amb_type = END;

  this->nmatches_plus_spliced_trims = this->nmatches_to_trims = this->queryend - this->querystart - this->nmismatches_bothdiff;

  this->alts_ncoords = 0;	/* To signify that alts have been resolved */
  FREE_OUT(this->alts_coords);
  FREE_OUT(this->alts_knowni);
  FREE_OUT(this->alts_nmismatches);
  FREE_OUT(this->alts_probs);

  return this->left;
}


float
Substring_compute_mapq (T this, char *quality_string) {
  int mapq_start, mapq_end;
  float best_loglik, loglik;
  Univcoord_T left, splicecoord;
  int i;

  /* mapq */
  mapq_start = this->querystart;
  mapq_end = this->queryend;

  /* It appears from simulated reads that it is better not to trim in
     computing MAPQ.  The correct mapping then tends to be selected
     with a higher MAPQ score. */
  /* But if all ends are terminals, then terminal parts should not be
     included in MAPQ scoring */

  if (this->nmismatches_bothdiff == 0) {  /* was this->exactp == true */
    /* this->mapq_loglik = MAPQ_loglik_exact(quality_string,0,querylength); */
    this->mapq_loglik = 0.0;

  } else if (this->alts_ncoords > 0) {
    if (this->plusp == true) {
      splicecoord = this->alts_coords[0];
      left = splicecoord - this->amb_splice_pos;
      best_loglik = MAPQ_loglik(/*ome*/genomebits,/*ome_alt*/genomebits_alt,
				this->query_compress,left,mapq_start,mapq_end,
				this->querylength,quality_string,/*plusp*/true,this->genestrand);
      for (i = 1; i < this->alts_ncoords; i++) {
	splicecoord = this->alts_coords[i];
	left = splicecoord - this->amb_splice_pos;
	if ((loglik = MAPQ_loglik(/*ome*/genomebits,/*ome_alt*/genomebits_alt,
				  this->query_compress,left,mapq_start,mapq_end,
				  this->querylength,quality_string,/*plusp*/true,this->genestrand)) > best_loglik) {
	  best_loglik = loglik;
	}
      }
    } else {
      splicecoord = this->alts_coords[0];
      left = splicecoord - (this->querylength - this->amb_splice_pos);
      best_loglik = MAPQ_loglik(/*ome*/genomebits,/*ome_alt*/genomebits_alt,
				this->query_compress,left,mapq_start,mapq_end,
				this->querylength,quality_string,/*plusp*/false,this->genestrand);
      for (i = 1; i < this->alts_ncoords; i++) {
	splicecoord = this->alts_coords[i];
	left = splicecoord - (this->querylength - this->amb_splice_pos);
	if ((loglik = MAPQ_loglik(/*ome*/genomebits,/*ome_alt*/genomebits_alt,
				  this->query_compress,left,mapq_start,mapq_end,
				  this->querylength,quality_string,/*plusp*/false,this->genestrand)) > best_loglik) {
	  best_loglik = loglik;
	}
      }
    }

    this->mapq_loglik = best_loglik;

  } else {
    debug2(printf("trim_querystart %d, trim_queryend %d, mapq_start = %d, mapq_end = %d\n",
		  this->querystart,this->querylength - this->queryend,mapq_start,mapq_end));
    this->mapq_loglik = MAPQ_loglik(/*ome*/genomebits,/*ome_alt*/genomebits_alt,
				    this->query_compress,this->left,mapq_start,mapq_end,
				    this->querylength,quality_string,this->plusp,this->genestrand);
    debug2(printf("Substring %u..%u gets loglik %f\n",this->genomicstart - this->chroffset,
		  this->genomicend - this->chroffset,this->mapq_loglik));
  }

  return this->mapq_loglik;
}


/* Sets genomic_bothdiff, genomic_refdiff, and nmismatches_refdiff */
int
Substring_display_prep (T this, char *queryuc_ptr, int querylength,
			int extraleft, int extraright, Genome_T genome) {

#if defined(HAVE_ALLOCA)
  char *genomic_diff, *gbuffer;
#else
  char *genomic_diff, *gbuffer;
#endif

  /* printf("left:%u nmismatches:%d extraleft:%d extraright:%d\n",
     this->left - this->chroffset,this->nmismatches_bothdiff,extraleft,extraright); */

  if (output_type == M8_OUTPUT && snps_iit == NULL) {
    /* Print procedures do not refer to genomic sequence */
    this->genomic_refdiff = this->genomic_bothdiff = (char *) NULL;
    this->nmismatches_refdiff = this->nmismatches_bothdiff;

  } else if (this->nmismatches_bothdiff == 0 && extraleft == 0 && extraright == 0 && snps_iit == NULL) {
    /* Print procedures will refer to query sequence */
    this->genomic_refdiff = this->genomic_bothdiff = (char *) NULL;
    this->nmismatches_refdiff = 0;

  } else if (this->plusp == true) {
    /* Used to be this->genomiclength, but doesn't work for large insertions */
#if defined(HAVE_ALLOCA)
    gbuffer = (char *) ALLOCA((querylength+1) * sizeof(char));
#else
    gbuffer = (char *) MALLOC((querylength+1) * sizeof(char));
#endif

    debug1(printf("Obtaining genomic_diff from left %u (%u) for querylength %d\n",
		  this->left,this->left - this->chroffset,querylength));
    Genome_fill_buffer_simple(genome,this->left,querylength,gbuffer);
    genomic_diff = gbuffer;

    Genome_mark_mismatches(genomic_diff,querylength,this->query_compress,
			   this->left,/*pos5*/this->querystart,/*pos3*/this->queryend,
			   /*plusp*/true,this->genestrand);

    /* Need to perform embellish to put dashes in */
    this->genomic_bothdiff = 
      embellish_genomic(this->genomic_bothdiff,genomic_diff,queryuc_ptr,this->querystart,this->queryend,
			querylength,this->mandatory_trim_querystart,this->mandatory_trim_queryend,
			extraleft,extraright,/*plusp*/true,this->genestrand);

    if (snps_iit == NULL) {
      this->genomic_refdiff = this->genomic_bothdiff;
      this->nmismatches_refdiff = this->nmismatches_bothdiff;

    } else {
      this->nmismatches_refdiff =
	Genome_count_mismatches_substring_ref(genomebits,this->query_compress,this->left,
					      /*pos5*/this->alignstart_trim - this->left,
					      /*pos3*/this->alignend_trim - this->left,
					      /*plusp*/true,this->genestrand);
	
      Genome_mark_mismatches_ref(genomic_diff,querylength,this->query_compress,this->left,
				 /*pos5*/this->querystart,/*pos3*/this->queryend,
				 /*plusp*/true,this->genestrand);
      if (output_type == STD_OUTPUT) {
	this->genomic_refdiff =
	  embellish_genomic(this->genomic_refdiff,genomic_diff,queryuc_ptr,this->querystart,this->queryend,
			    querylength,this->mandatory_trim_querystart,this->mandatory_trim_queryend,
			    extraleft,extraright,/*plusp*/true,this->genestrand);
      }
    }

    if (output_type == SAM_OUTPUT) {
      this->genomic_refdiff =
	embellish_genomic_sam(this->genomic_refdiff,genomic_diff,queryuc_ptr,this->querystart,this->queryend,
			      querylength,this->mandatory_trim_querystart,this->mandatory_trim_queryend,
			      extraleft,extraright,/*plusp*/true,this->genestrand);
    }

#if defined(HAVE_ALLOCA)
    FREEA(gbuffer);
#else
    FREE(gbuffer);
#endif

  } else {
    /* Used to be this->genomiclength, but doesn't work for large insertions */
#if defined(HAVE_ALLOCA)
    gbuffer = (char *) ALLOCA((querylength+1) * sizeof(char));
#else
    gbuffer = (char *) MALLOC((querylength+1) * sizeof(char));
#endif

    debug1(printf("Obtaining genomic_diff from left %u (%u) for querylength %d, and complemented\n",
		  this->left,this->left - this->chroffset,querylength));
    Genome_fill_buffer_simple(genome,this->left,querylength,gbuffer);
    genomic_diff = make_complement_inplace(gbuffer,querylength);

    Genome_mark_mismatches(genomic_diff,querylength,this->query_compress,
			   this->left,/*pos5*/querylength - this->queryend,
			   /*pos3*/querylength - this->querystart,
			   /*plusp*/false,this->genestrand);

    /* Need to perform embellish to put dashes in */
    this->genomic_bothdiff =
      embellish_genomic(this->genomic_bothdiff,genomic_diff,/*not queryrc*/queryuc_ptr,this->querystart,this->queryend,
			querylength,this->mandatory_trim_querystart,this->mandatory_trim_queryend,
			extraleft,extraright,/*plusp*/false,this->genestrand);

    if (snps_iit == NULL) {
      this->genomic_refdiff = this->genomic_bothdiff;
      this->nmismatches_refdiff = this->nmismatches_bothdiff;
      
    } else {
      this->nmismatches_refdiff = 
	Genome_count_mismatches_substring_ref(genomebits,this->query_compress,this->left,
					      /*pos5*/this->alignend_trim - this->left,
					      /*pos3*/this->alignstart_trim - this->left,/*plusp*/false,
					      this->genestrand);
      
      Genome_mark_mismatches_ref(genomic_diff,querylength,this->query_compress,this->left,
				 /*pos5*/querylength - this->queryend,
				 /*pos3*/querylength - this->querystart,
				 /*plusp*/false,this->genestrand);
      
      if (output_type == STD_OUTPUT) {
	this->genomic_refdiff =
	  embellish_genomic(this->genomic_refdiff,genomic_diff,/*not queryrc*/queryuc_ptr,this->querystart,this->queryend,
			    querylength,this->mandatory_trim_querystart,this->mandatory_trim_queryend,
			    extraleft,extraright,/*plusp*/false,this->genestrand);
      }
    }

    if (output_type == SAM_OUTPUT) {
      this->genomic_refdiff =
	embellish_genomic_sam(this->genomic_refdiff,genomic_diff,/*not queryrc*/queryuc_ptr,this->querystart,this->queryend,
			      querylength,this->mandatory_trim_querystart,this->mandatory_trim_queryend,
			      extraleft,extraright,/*plusp*/false,this->genestrand);
    }

#if defined(HAVE_ALLOCA)
    FREEA(gbuffer);
#else
    FREE(gbuffer);
#endif
  }

  return this->nmismatches_refdiff;
}


char *
Substring_genomic_sequence (int *seqlength, T this, Genome_T genome) {
  char *gbuffer;

  *seqlength = this->queryend - this->querystart;
  gbuffer = (char *) MALLOC((*seqlength+1) * sizeof(char));
  if (this->plusp == true) {
    Genome_fill_buffer_simple(genome,this->left + this->querystart,*seqlength,gbuffer);
    return gbuffer;
  } else {
    Genome_fill_buffer_simple(genome,this->left + (this->querylength - this->queryend),*seqlength,gbuffer);
    return make_complement_inplace(gbuffer,*seqlength);
  }
}


Univcoord_T
Substring_left (T this) {
  return this->left;
}


Univcoord_T
Substring_splicecoord_D (T this) {
  return this->splicecoord_D;
}

Univcoord_T
Substring_splicecoord_A (T this) {
  return this->splicecoord_A;
}


/* Called only by samprint, for donor or acceptor substrings */
char
Substring_chimera_strand (T this) {
  int sensedir;

  if ((sensedir = this->sensedir) != SENSE_ANTI) {
    if (this->plusp == true) {
      return '+';
    } else {
      return '-';
    }
  } else {
    if (this->plusp == true) {
      return '-';
    } else {
      return '+';
    }
  }
}


/* Called only by samprint */
Chrpos_T
Substring_chr_splicecoord_D (T this, char donor_strand) {
  if (donor_strand == '+') {
    return (Chrpos_T) (this->splicecoord_D - this->chroffset);
  } else if (donor_strand == '-') {
    return (Chrpos_T) (this->splicecoord_D - this->chroffset + 1);
  } else {
    abort();
  }
}

/* Called only by samprint */
Chrpos_T
Substring_chr_splicecoord_A (T this, char acceptor_strand) {
  if (acceptor_strand == '+') {
    return (Chrpos_T) (this->splicecoord_A - this->chroffset + 1);
  } else if (acceptor_strand == '-') {
    return (Chrpos_T) (this->splicecoord_A - this->chroffset);
  } else {
    abort();
  }
}

int
Substring_splicesitesD_knowni (T this) {
  return this->splicesitesD_knowni;
}

int
Substring_splicesitesA_knowni (T this) {
  return this->splicesitesA_knowni;
}

bool
Substring_plusp (T this) {
  return this->plusp;
}

int
Substring_sensedir (T this) {
  return this->sensedir;
}

int
Substring_genestrand (T this) {
  return this->genestrand;
}

char *
Substring_genomic_bothdiff (T this) {
  return this->genomic_bothdiff;
}

char *
Substring_genomic_refdiff (T this) {
  return this->genomic_refdiff;
}

int
Substring_nmismatches_bothdiff (T this) {
  return this->nmismatches_bothdiff;
}

int
Substring_nmismatches_refdiff (T this) {
  return this->nmismatches_refdiff;
}

int
Substring_nmatches_to_trims (T this) {
#if 0
  if (this->alts_ncoords > 0) {
    /* The entire substring is "trimmed" */
    return 0;
  } else {
    return this->queryend - this->querystart - this->nmismatches_bothdiff;
  }
#else
  return this->nmatches_to_trims;
#endif
}

/* nmatches_to_trims plus amb_length */
int
Substring_nmatches (T this) {
#if 0
  int amb_length;

  if (this->alts_ncoords > 0) {
    return this->start_amb_length + this->end_amb_length;
  } else {
    amb_length = 0;
    if (this->trim_querystart_splicep == true) {
      amb_length += this->querystart;
    }
    if (this->trim_queryend_splicep == true) {
      amb_length += this->querylength - this->queryend;
    }
    return amb_length + /*nmatches_to_trims*/ (this->queryend - this->querystart - this->nmismatches_bothdiff);
  }
#else
  return this->nmatches_plus_spliced_trims;
#endif
}


Endtype_T
Substring_start_endtype (T this) {
  return this->start_endtype;
}

Endtype_T
Substring_end_endtype (T this) {
  return this->end_endtype;
}

float
Substring_mapq_loglik (T this) {
  return this->mapq_loglik;
}

int
Substring_trim_querystart (T this) {
  return this->querystart;
}

int
Substring_trim_queryend (T this) {
  return this->querylength - this->queryend;
}

bool
Substring_trim_querystart_splicep (T this) {
  return this->trim_querystart_splicep;
}

bool
Substring_trim_queryend_splicep (T this) {
  return this->trim_queryend_splicep;
}

int
Substring_mandatory_trim_querystart (T this) {
  return this->mandatory_trim_querystart;
}

int
Substring_mandatory_trim_queryend (T this) {
  return this->mandatory_trim_queryend;
}


int
Substring_querystart (T this) {
  return this->querystart;
}

int
Substring_querystart_pretrim (T this) {
  return this->querystart_pretrim;
}

int
Substring_queryend (T this) {
  return this->queryend;
}

int
Substring_querylength (T this) {
  return this->querylength;
}

int
Substring_match_length (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return this->queryend - this->querystart;
  }
}

int
Substring_match_length_pretrim (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return this->queryend_pretrim - this->querystart_pretrim;
  }
}


/* Mapped and unmapped */
int
Substring_amb_length (T this) {
  return this->start_amb_length + this->end_amb_length;
}


/* Mapped and unmapped */
int
Substring_start_amb_length (T this) {
  return this->start_amb_length;
}


/* Mapped and unmapped */
int
Substring_end_amb_length (T this) {
  return this->end_amb_length;
}


Chrpos_T
Substring_genomic_alignment_length (T this) {
  /* Don't check for alts, as long as we return 0 */
  /* assert(this->alts_ncoords == 0); */

  if (this == NULL) {
    return (Chrpos_T) 0;
  } else if (this->alts_ncoords > 0) {
    return (Chrpos_T) 0;
  } else if (this->alignend_trim > this->alignstart_trim) {
    return this->alignend_trim - this->alignstart_trim;
  } else {
    return this->alignstart_trim - this->alignend_trim;
  }
}


Chrnum_T
Substring_chrnum (T this) {
  return this->chrnum;
}

Univcoord_T
Substring_chroffset (T this) {
  return this->chroffset;
}

Univcoord_T
Substring_chrhigh (T this) {
  return this->chrhigh;
}

Chrpos_T
Substring_chrlength (T this) {
  return this->chrlength;
}

Chrpos_T
Substring_chrpos_low (T this) {
  /* Need to check for alts, because then alignstart_trim and alignend_trim are set to 0 */
  assert(this->alts_ncoords == 0);

  if (this->plusp == true) {
    return this->alignstart_trim - this->chroffset;
  } else {
    return this->alignend_trim - this->chroffset;
  }
}

Chrpos_T
Substring_chrpos_high (T this) {
  /* Need to check for alts, because then alignstart_trim and alignend_trim are set to 0 */
  assert(this->alts_ncoords == 0);

  if (this->plusp == true) {
    return this->alignend_trim - this->chroffset;
  } else {
    return this->alignstart_trim - this->chroffset;
  }
}


Chrpos_T
Substring_alignstart_trim_chr (T this) {
  /* Need to check for alts, because then alignstart_trim and alignend_trim are set to 0 */
  assert(this->alts_ncoords == 0);

  /* Previously had provision for alts_coords */
  return this->alignstart_trim - this->chroffset;
}

Chrpos_T
Substring_alignend_trim_chr (T this) {
  /* Need to check for alts, because then alignstart_trim and alignend_trim are set to 0 */
  assert(this->alts_ncoords == 0);

  /* Previously had provision for alts_coords */
  return this->alignend_trim - this->chroffset;
}

Univcoord_T
Substring_alignstart_trim (T this) {
  return this->alignstart_trim;
}

Univcoord_T
Substring_alignend_trim (T this) {
  return this->alignend_trim;
}


Univcoord_T
Substring_left_genomicseg (T this) {
  return this->left;
}

Univcoord_T
Substring_genomicstart (T this) {
  return this->genomicstart;
}

Univcoord_T
Substring_genomicend (T this) {
  return this->genomicend;
}

double
Substring_amb_prob (T this) {
  if (this->amb_type == DON) {
    return this->siteD_prob;
  } else if (this->amb_type == ACC) {
    return this->siteA_prob;
  } else {
    return 0.0;
  }
}


double
Substring_amb_donor_prob (T this) {
  double max;
  int i;

  if (this->alts_ncoords == 0) {
    return this->siteD_prob;
  } else if (this->amb_type == DON) {
    max = this->alts_probs[0];
    for (i = 1; i < this->alts_ncoords; i++) {
      if (this->alts_probs[i] > max) {
	max = this->alts_probs[i];
      }
    }
    return max;
  } else {
    return this->siteD_prob;
  }
}

double
Substring_amb_acceptor_prob (T this) {
  double max;
  int i;

  if (this->alts_ncoords == 0) {
    return this->siteA_prob;
  } else if (this->amb_type == ACC) {
    max = this->alts_probs[0];
    for (i = 1; i < this->alts_ncoords; i++) {
      if (this->alts_probs[i] > max) {
	max = this->alts_probs[i];
      }
    }
    return max;
  } else {
    return this->siteA_prob;
  }
}



double
Substring_siteD_prob (T this) {
  return this->siteD_prob;
}  

double
Substring_siteA_prob (T this) {
  return this->siteA_prob;
}  

int
Substring_siteD_pos (T this) {
  return this->siteD_pos;
}

int
Substring_siteA_pos (T this) {
  return this->siteA_pos;
}

int
Substring_siteN_pos (T this) {
  return this->siteN_pos;
}


bool
Substring_ambiguous_p (T this) {
  if (this->alts_ncoords > 0) {
    return false;
  } else if (this->amb_type == END) {
    return false;
  } else {
    return true;
  }
}

bool
Substring_has_alts_p (T this) {
  if (this->alts_ncoords > 0) {
    return true;
  } else {
    return false;
  }
}

int
Substring_alts_ncoords (T this) {
  return this->alts_ncoords;
}

Univcoord_T *
Substring_alts_coords (T this) {
  return this->alts_coords;
}

double
Substring_alts_common_prob (T this) {
  if (this->amb_type == DON) {
    return this->siteA_prob;
  } else if (this->amb_type == ACC) {
    return this->siteD_prob;
  } else {
    return 0.0;
  }
}

void
Substring_print_alts_coords (T this) {
  int i;

  for (i = 0; i < this->alts_ncoords; i++) {
    printf(" %u",this->alts_coords[i] - this->chroffset);
  }

  return;
}

int *
Substring_alts_nmismatches (T this) {
  return this->alts_nmismatches;
}

/* circularpos measures query distance from SAM chrlow to origin */
int
Substring_circularpos (T this) {
  if (this == NULL) {
    return -1;

  } else if (this->plusp == true) {
    debug12(printf("Substring circularpos plus: %u..%u vs %u+%u\n",
		   this->alignstart_trim,this->alignend_trim,this->chroffset,this->chrlength));
    if (this->alignstart_trim > this->chroffset + this->chrlength) {
      debug12(printf("all of substring is in the upper part, so previous substring was in lower part\n"));
      debug12(printf("returning querystart %d\n",this->querystart));
      return this->querystart;

    } else if (this->alignend_trim > this->chroffset + this->chrlength) {
      debug12(printf("substring straddles the circular origin\n"));
      debug12(printf("returning querystart %d + chroffset %u + chrlength %u - alignstart_trim %d\n",
		     this->querystart,this->chroffset,this->chrlength,this->alignstart_trim));
      /* return (this->querystart - this->trim_querystart) + (this->chroffset + this->chrlength) - this->alignstart; */
      return this->querystart + (this->chroffset + this->chrlength) - this->alignstart_trim;

    } else {
      return -1;
    }

  } else {
    debug12(printf("Substring circularpos minus: %u..%u vs %u+%u\n",
		   this->alignstart_trim,this->alignend_trim,this->chroffset,this->chrlength));
    if (this->alignend_trim > this->chroffset + this->chrlength) {
      debug12(printf("all of substring is in the upper part, so previous substring was in lower part\n"));
      debug12(printf("returning querylength %d - queryend %d\n",this->querylength,this->queryend));
      return this->querylength - this->queryend;

    } else if (this->alignstart_trim > this->chroffset + this->chrlength) {
      debug12(printf("returning querylength %d - queryend %d + chroffset %u + chrlength %u - alignend_trim %d\n",
		     this->querylength,this->queryend,this->chroffset,this->chrlength,this->alignend_trim));
      /* return ((this->querylength - this->trim_queryend) - this->queryend) + (this->chroffset + this->chrlength) - this->alignend; */
      return (this->querylength - this->queryend) + (this->chroffset + this->chrlength) - this->alignend_trim;

    } else {
      return -1;
    }
  }
}


T
Substring_copy (T old) {
  T new;

  if (old == NULL) {
    return (T) NULL;
  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug2(printf("substring %p is a copy of %p\n",new,old));

    new->nmismatches_bothdiff = old->nmismatches_bothdiff;
    new->nmismatches_refdiff = old->nmismatches_refdiff;

    new->nmatches_plus_spliced_trims = old->nmatches_plus_spliced_trims;
    new->nmatches_to_trims = old->nmatches_to_trims;

    new->mandatory_trim_querystart = old->mandatory_trim_querystart;
    new->mandatory_trim_queryend = old->mandatory_trim_queryend;
    new->start_amb_length = old->start_amb_length;
    new->end_amb_length = old->end_amb_length;

    new->trim_querystart_splicep = old->trim_querystart_splicep;
    new->trim_queryend_splicep = old->trim_queryend_splicep;

    new->chrnum = old->chrnum;
    new->chroffset = old->chroffset;
    new->chrhigh = old->chrhigh;
    new->chrlength = old->chrlength;

    new->left = old->left;
    new->genomicstart = old->genomicstart;
    new->genomicend = old->genomicend;

    new->start_endtype = old->start_endtype;
    new->end_endtype = old->end_endtype;

    new->querystart_pretrim = old->querystart_pretrim;
    new->queryend_pretrim = old->queryend_pretrim;
    new->querystart = old->querystart;
    new->queryend = old->queryend;
    new->amb_splice_pos = old->amb_splice_pos;
    new->querylength = old->querylength;

    new->alignstart_trim = old->alignstart_trim;
    new->alignend_trim = old->alignend_trim;

    new->plusp = old->plusp;
    new->genestrand = old->genestrand;
    new->query_compress = old->query_compress;

    if (old->genomic_bothdiff == NULL) {
      new->genomic_bothdiff = (char *) NULL;
      new->genomic_refdiff = (char *) NULL;

    } else if (old->genomic_refdiff == old->genomic_bothdiff) {
      assert((int) strlen(old->genomic_bothdiff) == old->querylength);
      new->genomic_bothdiff = (char *) CALLOC_OUT(strlen(old->genomic_bothdiff)+1,sizeof(char));
      strcpy(new->genomic_bothdiff,old->genomic_bothdiff);
      new->genomic_refdiff = new->genomic_bothdiff;

    } else {
      assert((int) strlen(old->genomic_bothdiff) == old->querylength);
      assert((int) strlen(old->genomic_refdiff) == old->querylength);
      new->genomic_bothdiff = (char *) CALLOC_OUT(strlen(old->genomic_bothdiff)+1,sizeof(char));
      strcpy(new->genomic_bothdiff,old->genomic_bothdiff);
      new->genomic_refdiff = (char *) CALLOC_OUT(strlen(old->genomic_refdiff)+1,sizeof(char));
      strcpy(new->genomic_refdiff,old->genomic_refdiff);
    }

    new->mapq_loglik = old->mapq_loglik;

    new->sensedir = old->sensedir;

    new->splicecoord_D = old->splicecoord_D;
    new->splicesitesD_knowni = old->splicesitesD_knowni;
    new->siteD_pos = old->siteD_pos;
    new->siteD_prob = old->siteD_prob;

    new->splicecoord_A = old->splicecoord_A;
    new->splicesitesA_knowni = old->splicesitesA_knowni;
    new->siteA_pos = old->siteA_pos;
    new->siteA_prob = old->siteA_prob;

    new->siteN_pos = old->siteN_pos;

    if (old->alts_ncoords == 0) {
      new->alts_ncoords = 0;
      new->alts_coords = (Univcoord_T *) NULL;
      new->alts_knowni = (int *) NULL;
      new->alts_nmismatches = (int *) NULL;
      new->alts_probs = (double *) NULL;
    } else {
      new->alts_ncoords = old->alts_ncoords;
      new->alts_coords = (Univcoord_T *) MALLOC_OUT(old->alts_ncoords * sizeof(Univcoord_T));
      new->alts_knowni = (int *) MALLOC_OUT(old->alts_ncoords * sizeof(int));
      new->alts_nmismatches = (int *) MALLOC_OUT(old->alts_ncoords * sizeof(int));
      new->alts_probs = (double *) MALLOC_OUT(old->alts_ncoords * sizeof(double));

      memcpy(new->alts_coords,old->alts_coords,old->alts_ncoords * sizeof(Univcoord_T));
      memcpy(new->alts_knowni,old->alts_knowni,old->alts_ncoords * sizeof(int));
      memcpy(new->alts_nmismatches,old->alts_nmismatches,old->alts_ncoords * sizeof(int));
      memcpy(new->alts_probs,old->alts_probs,old->alts_ncoords * sizeof(double));
    }
    new->amb_type = old->amb_type;

    return new;
  }
}


/* nmismatches is between donor_pos and endpos */
T
Substring_new_donor (int nmismatches, Univcoord_T donor_coord, int donor_knowni,
		     int querystart, int queryend, int sitepos, double donor_prob,
		     bool splice_querystart_p, Splicetype_T splicetype_querystart, double ambig_prob_querystart,
		     bool splice_queryend_p, Splicetype_T splicetype_queryend, double ambig_prob_queryend,

		     Univcoord_T left, Compress_T query_compress, int querylength, bool plusp, int genestrand, int sensedir,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength) {
  T new;

  /* Previously checked if left >= chroffset + chrlength to exclude
     the duplicate length, but now excluding all translocations to
     circular chromosomes */

  if (chroffset + chrlength < chrhigh) {
    /* Don't splice to circular chromosomes */
    return (T) NULL;

  } else if ((new = Substring_new(nmismatches,left,querystart,queryend,querylength,
				  plusp,genestrand,query_compress,chrnum,chroffset,chrhigh,chrlength,
				  splice_querystart_p,splicetype_querystart,ambig_prob_querystart,
				  splice_queryend_p,splicetype_queryend,ambig_prob_queryend,sensedir)) == NULL) {
    return (T) NULL;
  }

  debug2(printf("Making new donor with splicesites_i %d, coord %u and left %u, plusp %d, sensedir %d, query %d..%d, trim %d..%d\n",
		donor_knowni,donor_coord,left,plusp,sensedir,new->querystart,new->queryend,
		new->querystart,querylength - new->queryend));
  debug2(printf("Setting siteD_prob to be %f\n",donor_prob));

  new->splicecoord_D = donor_coord;
  new->splicesitesD_knowni = donor_knowni;

  new->sensedir = sensedir;

  new->siteD_pos = sitepos;	/* not a qpos */
  new->siteD_prob = donor_prob;

  return new;
}


/* nmismatches is between acceptor_pos and endpos */
T
Substring_new_acceptor (int nmismatches, Univcoord_T acceptor_coord, int acceptor_knowni,
			int querystart, int queryend, int sitepos, double acceptor_prob,
			bool splice_querystart_p, Splicetype_T splicetype_querystart, double ambig_prob_querystart,
			bool splice_queryend_p, Splicetype_T splicetype_queryend, double ambig_prob_queryend,

			Univcoord_T left, Compress_T query_compress, int querylength, bool plusp, int genestrand, int sensedir,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength) {
  T new;

  /* Previously checked if left >= chroffset + chrlength to exclude
     the duplicate length, but now excluding all translocations to
     circular chromosomes */

  if (chroffset + chrlength < chrhigh) {
    /* Don't splice to circular chromosomes */
    return (T) NULL;

  } else if ((new = Substring_new(nmismatches,left,querystart,queryend,querylength,
				  plusp,genestrand,query_compress,chrnum,chroffset,chrhigh,chrlength,
				  splice_querystart_p,splicetype_querystart,ambig_prob_querystart,
				  splice_queryend_p,splicetype_queryend,ambig_prob_queryend,sensedir)) == NULL) {
    return (T) NULL;
  }
 
  debug2(printf("Making new acceptor with splicesites_i %d, coord %u and left %u, plusp %d, sensedir %d, query %d..%d, trim %d..%d\n",
		acceptor_knowni,acceptor_coord,left,plusp,sensedir,new->querystart,new->queryend,
		new->querystart,querylength - new->queryend));
  debug2(printf("Setting siteA_prob to be %f\n",acceptor_prob));

  new->splicecoord_A = acceptor_coord;
  new->splicesitesA_knowni = acceptor_knowni;

  new->sensedir = sensedir;

  new->siteA_pos = sitepos;	/* not a qpos */
  new->siteA_prob = acceptor_prob;

  return new;
}


void
Substring_label_donor (T this, int splice_pos, double donor_prob, int sensedir_distant_guess) {
  if (this->plusp == true) {
    this->splicecoord_D = this->left + splice_pos;
  } else {
    this->splicecoord_D = this->left + this->querylength - splice_pos;
  }
  /* new->splicesitesD_knowni = donor_knowni; */
  /* new->sensedir = sensedir; */
  this->sensedir = sensedir_distant_guess;

  this->siteD_pos = splice_pos;	/* not a qpos */
  this->siteD_prob = donor_prob;

  /* printf("Labeling donor with splicecoord %llu and prob %f\n",this->splicecoord_D,donor_prob); */

  return;
}


void
Substring_label_acceptor (T this, int splice_pos, double acceptor_prob, int sensedir_distant_guess) {
  if (this->plusp == true) {
    this->splicecoord_A = this->left + splice_pos;
  } else {
    this->splicecoord_A = this->left + this->querylength - splice_pos;
  }
  /* new->splicesitesA_knowni = acceptor_knowni; */
  /* new->sensedir = sensedir; */
  this->sensedir = sensedir_distant_guess;

  this->siteA_pos = splice_pos;	/* not a qpos */
  this->siteA_prob = acceptor_prob;

  /* printf("Labeling acceptor with splicecoord %llu and prob %f\n",this->splicecoord_A,acceptor_prob); */

  return;
}



/* Used for distant DNA splicing.  Uses querystart and queryend, not qstart and qend.  middlepos is queryend */
T
Substring_new_startfrag (int nmismatches, int querystart, int queryend,
			 Univcoord_T left, Compress_T query_compress, int querylength, bool plusp, int genestrand,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength) {
  T new;

  /* Previously checked if left >= chroffset + chrlength to exclude
     the duplicate length, but now excluding all translocations to
     circular chromosomes */

  if (chroffset + chrlength < chrhigh) {
    /* Don't splice to circular chromosomes */
    return (T) NULL;

  } else if ((new = Substring_new(nmismatches,left,querystart,queryend,querylength,
				  plusp,genestrand,query_compress,chrnum,chroffset,chrhigh,chrlength,
				  /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				  /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				  /*sensedir*/SENSE_NULL)) == NULL) {
    return (T) NULL;

  } else {
    debug2(printf("Making new startfrag with left %u, plusp %d, query %d..%d, genome %u..%u\n",
		  left,plusp,querystart,queryend,
		  Substring_alignstart_trim(new),Substring_alignend_trim(new)));
    new->siteN_pos = queryend;
    return new;
  }
}


/* Used for distant DNA splicing.  Uses querystart and queryend, not qstart and qend.  middlepos is querystart */
T
Substring_new_endfrag (int nmismatches, int querystart, int queryend,
		       Univcoord_T left, Compress_T query_compress, int querylength, bool plusp, int genestrand,
		       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength) {
  T new;

  /* Previously checked if left >= chroffset + chrlength to exclude
     the duplicate length, but now excluding all translocations to
     circular chromosomes */

  if (chroffset + chrlength < chrhigh) {
    /* Don't splice to circular chromosomes */
    return (T) NULL;

  } else if ((new = Substring_new(nmismatches,left,querystart,queryend,querylength,
				  plusp,genestrand,query_compress,chrnum,chroffset,chrhigh,chrlength,
				  /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				  /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				  /*sensedir*/SENSE_NULL)) == NULL) {
    return (T) NULL;

  } else {
    debug2(printf("Making new endfrag with left %u, plusp %d, %d..%d, genome %u..%u\n",
		  left,plusp,querystart,queryend,
		  Substring_alignstart_trim(new),Substring_alignend_trim(new)));
    new->siteN_pos = querystart;
    return new;
  }
}


T
Substring_trim_startfrag (int nmismatches, T old, int new_queryend) {
  T new;

  debug2(printf("Entered Substring_trim_startfrag with %d..%d -> %d..%d\n",
		old->querystart,old->queryend,old->querystart,new_queryend));

  if ((new = Substring_new(nmismatches,old->left,old->querystart,new_queryend,old->querylength,
			   old->plusp,old->genestrand,old->query_compress,
			   old->chrnum,old->chroffset,old->chrhigh,old->chrlength,
			   /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
			   /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
			   /*sensedir*/SENSE_NULL)) == NULL) {
    debug2(printf("Returning NULL\n"));
    return (T) NULL;

  } else {
    debug2(printf("Making new startfrag with left %u, plusp %d, query %d..%d, genome %u..%u\n",
		  old->left,old->plusp,old->querystart,new_queryend,
		  Substring_alignstart_trim(new),Substring_alignend_trim(new)));
    new->siteN_pos = new_queryend;
    return new;
  }
}

T
Substring_trim_endfrag (int nmismatches, T old, int new_querystart) {
  T new;

  debug2(printf("Entered Substring_trim_endfrag with %d..%d -> %d..%d\n",
		old->querystart,old->queryend,new_querystart,old->queryend));

  if ((new = Substring_new(nmismatches,old->left,new_querystart,old->queryend,old->querylength,
			   old->plusp,old->genestrand,old->query_compress,
			   old->chrnum,old->chroffset,old->chrhigh,old->chrlength,
			   /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
			   /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
			   /*sensedir*/SENSE_NULL)) == NULL) {
    debug2(printf("Returning NULL\n"));
    return (T) NULL;

  } else {
    debug2(printf("Making new endfrag with left %u, plusp %d, %d..%d, genome %u..%u\n",
		  old->left,old->plusp,new_querystart,old->queryend,
		  Substring_alignstart_trim(new),Substring_alignend_trim(new)));
    new->siteN_pos = new_querystart;
    return new;
  }
}



static int
ascending_siteD_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->siteD_pos < y->siteD_pos) {
    return -1;
  } else if (x->siteD_pos > y->siteD_pos) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}

static int
ascending_siteA_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->siteA_pos < y->siteA_pos) {
    return -1;
  } else if (x->siteA_pos > y->siteA_pos) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}

static int
ascending_siteN_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->siteN_pos < y->siteN_pos) {
    return -1;
  } else if (x->siteN_pos > y->siteN_pos) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}

static int
descending_siteD_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->siteD_pos < y->siteD_pos) {
    return -1;
  } else if (x->siteD_pos > y->siteD_pos) {
    return +1;
  } else if (x->genomicstart > y->genomicstart) {
    return -1;
  } else if (x->genomicstart < y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}

static int
descending_siteA_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->siteA_pos < y->siteA_pos) {
    return -1;
  } else if (x->siteA_pos > y->siteA_pos) {
    return +1;
  } else if (x->genomicstart > y->genomicstart) {
    return -1;
  } else if (x->genomicstart < y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}

static int
descending_siteN_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->siteN_pos < y->siteN_pos) {
    return -1;
  } else if (x->siteN_pos > y->siteN_pos) {
    return +1;
  } else if (x->genomicstart > y->genomicstart) {
    return -1;
  } else if (x->genomicstart < y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Substring_sort_siteD_halves (List_T hitlist, Listpool_T listpool, bool ascendingp) {
  List_T sorted = NULL;
  T x, *hits;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitlist);
  debug(printf("Checking %d spliceends for duplicates...",n));
  if (n == 0) {
    debug(printf("\n"));
    return NULL;
  }

  hits = (T *) MALLOCA(n * sizeof(T));
  List_fill_array((void **) hits,hitlist);
  /* List_free(&hitlist); -- allocated by Listpool_push */

  if (ascendingp == true) {
    qsort(hits,n,sizeof(T),ascending_siteD_pos_cmp);
  } else {
    qsort(hits,n,sizeof(T),descending_siteD_pos_cmp);
  }

  /* Check for duplicates */
  eliminate = (bool *) CALLOCA(n,sizeof(bool));
  for (i = 0; i < n; i++) {
    x = hits[i];
    j = i+1;
    while (j < n && hits[j]->siteD_pos == x->siteD_pos && hits[j]->genomicstart == x->genomicstart) {
      eliminate[j] = true;
      j++;
    }
  }

  debug(j = 0);
  for (i = n-1; i >= 0; i--) {
    x = hits[i];
    if (eliminate[i] == false) {
      sorted = Listpool_push(sorted,listpool,(void *) x);
    } else {
      Substring_free(&x);
      debug(j++);
    }
  }
  debug(printf("%d eliminated\n",j));

  FREEA(hits);
  FREEA(eliminate);

  return sorted;
}

List_T
Substring_sort_siteA_halves (List_T hitlist, Listpool_T listpool, bool ascendingp) {
  List_T sorted = NULL;
  T x, *hits;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitlist);
  debug(printf("Checking %d spliceends for duplicates...",n));
  if (n == 0) {
    debug(printf("\n"));
    return NULL;
  }

  hits = (T *) MALLOCA(n * sizeof(T));
  List_fill_array((void **) hits,hitlist);
  /* List_free(&hitlist); -- allocated by Listpool_push */

  if (ascendingp == true) {
    qsort(hits,n,sizeof(T),ascending_siteA_pos_cmp);
  } else {
    qsort(hits,n,sizeof(T),descending_siteA_pos_cmp);
  }

  /* Check for duplicates */
  eliminate = (bool *) CALLOCA(n,sizeof(bool));
  for (i = 0; i < n; i++) {
    x = hits[i];
    j = i+1;
    while (j < n && hits[j]->siteA_pos == x->siteA_pos && hits[j]->genomicstart == x->genomicstart) {
      eliminate[j] = true;
      j++;
    }
  }

  debug(j = 0);
  for (i = n-1; i >= 0; i--) {
    x = hits[i];
    if (eliminate[i] == false) {
      sorted = Listpool_push(sorted,listpool,(void *) x);
    } else {
      Substring_free(&x);
      debug(j++);
    }
  }
  debug(printf("%d eliminated\n",j));

  FREEA(hits);
  FREEA(eliminate);

  return sorted;
}

List_T
Substring_sort_siteN_halves (List_T hitlist, Listpool_T listpool, bool ascendingp) {
  List_T sorted = NULL;
  T x, *hits;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitlist);
  debug(printf("Checking %d spliceends for duplicates...",n));
  if (n == 0) {
    debug(printf("\n"));
    return NULL;
  }

  hits = (T *) MALLOCA(n * sizeof(T));
  List_fill_array((void **) hits,hitlist);
  /* List_free(&hitlist); -- allocated by Listpool_push */

  if (ascendingp == true) {
    qsort(hits,n,sizeof(T),ascending_siteN_pos_cmp);
  } else {
    qsort(hits,n,sizeof(T),descending_siteN_pos_cmp);
  }

  /* Check for duplicates */
  eliminate = (bool *) CALLOCA(n,sizeof(bool));
  for (i = 0; i < n; i++) {
    x = hits[i];
    j = i+1;
    while (j < n && hits[j]->siteN_pos == x->siteN_pos && hits[j]->genomicstart == x->genomicstart) {
      eliminate[j] = true;
      j++;
    }
  }

  debug(j = 0);
  for (i = n-1; i >= 0; i--) {
    x = hits[i];
    if (eliminate[i] == false) {
      sorted = Listpool_push(sorted,listpool,(void *) x);
    } else {
      Substring_free(&x);
      debug(j++);
    }
  }
  debug(printf("%d eliminated\n",j));

  FREEA(hits);
  FREEA(eliminate);

  return sorted;
}


static void
print_snp_labels (Filestring_T fp, T this, Shortread_T queryseq) {
  int *snps, nsnps, querypos, i;
  char *label, *seq1, *seq2;
  bool allocp, printp = false;
  Interval_T interval;
  Chrpos_T position;

  if (this->plusp == true) {
    snps = IIT_get_with_divno(&nsnps,snps_iit,
			      snps_divint_crosstable[this->chrnum],
			      this->alignstart_trim - this->chroffset + 1,this->alignend_trim - this->chroffset,
			      /*sortp*/false);
  } else {
    snps = IIT_get_with_divno(&nsnps,snps_iit,
			      snps_divint_crosstable[this->chrnum],
			      this->alignend_trim - this->chroffset + 1,this->alignstart_trim - this->chroffset,
			      /*sortp*/false);
  }

  FPRINTF(fp,",snps:");

  seq1 = Shortread_fullpointer_uc(queryseq);
  if (this->genomic_bothdiff == NULL) {
    seq2 = Shortread_fullpointer(queryseq);
  } else {
    seq2 = this->genomic_bothdiff;
  }

  if (this->plusp) {

#if 0
    for (i = 0; i < nsnps; i++) {
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = position - (this->genomicstart-this->chroffset) - 1;
      printf("%d ",querypos);
    }
    printf("\n");
#endif

    for (i = 0; i < nsnps; i++) {
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = position - (this->genomicstart-this->chroffset) - 1;
      /* assert(querypos >= 0 && querypos < this->genomiclength); */

#if 0
      alleles = IIT_typestring(snps_iit,Interval_type(interval));
      if (c == alleles[0] || c == alleles[1]) ;
#endif

      if (isupper(seq2[querypos]) && seq1[querypos] != seq2[querypos]) {
	label = IIT_label(snps_iit,snps[i],&allocp);
	if (printp) {
	  FPRINTF(fp,"|");
	}
	FPRINTF(fp,"%d@",querypos+1);
	FPRINTF(fp,"%s",label);
	printp = true;
	if (allocp) FREE(label);
      }
    }

  } else {

    for (i = nsnps-1; i >= 0; i--) {
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = (this->genomicstart-this->chroffset) - position;
      /* assert(querypos >= 0 && querypos < this->genomiclength); */

#if 0
      /* printf("\n%d%c\n",querypos,c); */
      alleles = IIT_typestring(snps_iit,Interval_type(interval));
      if (c == alleles[0] || c == alleles[1]) ;
#endif

      if (isupper(seq2[querypos]) && seq1[querypos] != seq2[querypos]) {
	label = IIT_label(snps_iit,snps[i],&allocp);
	if (printp) {
	  FPRINTF(fp,"|");
	}
	FPRINTF(fp,"%d@",querypos+1);
	FPRINTF(fp,"%s",label);
	printp = true;
	if (allocp) FREE(label);
      }
    }

  }

  FREE(snps);

  return;
}


Chrpos_T
Substring_compute_chrpos (T this, int hardclip_low, bool hide_soft_clips_p) {
  Chrpos_T chrpos_low;

  if (hide_soft_clips_p == true) {
    if (this->plusp == true) {
      /* Add 1 to report in 1-based coordinates */
      chrpos_low = this->genomicstart - this->chroffset + 1;
      chrpos_low += hardclip_low;
      /* *chrpos_high = this->genomicend - this->chroffset + 1; */
      /* *chrpos_high -= hardclip_high; */

    } else {
      /* Add 1 to report in 1-based coordinates */
      chrpos_low = this->genomicend - this->chroffset + 1;
      chrpos_low += hardclip_low;
      /* *chrpos_high = this->genomicstart - this->chroffset + 1; */
      /* *chrpos_high -= hardclip_high; */
    }

  } else {
    if (this->plusp == true) {
      chrpos_low = this->genomicstart - this->chroffset + 1;
      if (this->querystart > hardclip_low) {
	chrpos_low += this->querystart; /* not querystart_orig */
      } else {
	chrpos_low += hardclip_low;
      }

#if 0
      *chrpos_high = this->genomicend - this->chroffset + 1;
      if (this->querylength - this->queryend > hardclip_high) {
	*chrpos_high -= this->querylength - this->queryend;
      } else {
	*chrpos_high -= hardclip_high;
      }
#endif
	
    } else {
      chrpos_low = this->genomicend - this->chroffset + 1;
      if (this->querylength - this->queryend > hardclip_low) {
	chrpos_low += this->querylength - this->queryend; /* not queryend_orig */
      } else {
	chrpos_low += hardclip_low;
      }

#if 0
      *chrpos_high = this->genomicstart - this->chroffset + 1;
      if (this->querystart > hardclip_high) {
	*chrpos_high -= this->querystart;
      } else {
	*chrpos_high -= hardclip_high;
      }
#endif
    }
  }
    
  return chrpos_low;
}



/* Taken from NCBI Blast 2.2.29, algo/blast/core/blast_stat.c */
/* Karlin-Altschul formula: m n exp(-lambda * S + log k) = k m n exp(-lambda * S) */
/* Also in pair.c */

static double
blast_evalue (int alignlength, int nmismatches) {
  double k = 0.1;
  double lambda = 1.58;		/* For a +1, -1 scoring scheme */
  double score;
  
  score = (double) ((alignlength - nmismatches) /* scored as +1 */ - nmismatches /* scored as -1 */);

  return k * (double) alignlength * genomelength * exp(-lambda * score);
}

static double
blast_bitscore (int alignlength, int nmismatches) {
  double k = 0.1;
  double lambda = 1.58;		/* For a +1, -1 scoring scheme */
  double score;
  
  score = (double) ((alignlength - nmismatches) /* scored as +1 */ - nmismatches /* scored as -1 */);
  return (score * lambda - log(k)) / log(2.0);
}


double
Substring_evalue (T substring) {
  int alignlength_trim;

#if 0
  /* Doesn't work for ambiguous substrings */
  if (substring->plusp == true) {
    alignlength_trim = (int) (substring->alignend_trim - substring->alignstart_trim);
  } else {
    alignlength_trim = (int) (substring->alignstart_trim - substring->alignend_trim);
  }
  assert(alignlength_trim == substring->queryend - substring->querystart);
#else
  alignlength_trim = substring->queryend - substring->querystart;
#endif

  return blast_evalue(alignlength_trim,substring->nmismatches_bothdiff);
}


void
Substring_print_m8 (Filestring_T fp, T substring, Shortread_T headerseq, char *acc_suffix,
		    char *chr, bool invertp) {
  double identity;
  int alignlength_trim;

  FPRINTF(fp,"%s%s",Shortread_accession(headerseq),acc_suffix); /* field 0: accession */

  FPRINTF(fp,"\t%s",chr);	/* field 1: chr */

  /* field 2: identity */
  if (substring->plusp == true) {
    alignlength_trim = (int) (substring->alignend_trim - substring->alignstart_trim);
  } else {
    alignlength_trim = (int) (substring->alignstart_trim - substring->alignend_trim);
  }

  identity = (double) (alignlength_trim - substring->nmismatches_bothdiff)/(double) alignlength_trim;
  FPRINTF(fp,"\t%.1f",100.0*identity);


  FPRINTF(fp,"\t%d",alignlength_trim); /* field 3: query length */

  FPRINTF(fp,"\t%d",substring->nmismatches_bothdiff); /* field 4: nmismatches */

  FPRINTF(fp,"\t0");		/* field 5: gap openings */

  FPRINTF(fp,"\t%d",substring->querystart + 1); /* field 6: query start */

  FPRINTF(fp,"\t%d",substring->queryend); /* field 7: query end */

  /* fields 8 and 9: chr start and end */
  if (substring->plusp == true) {
    if (invertp == false) {
      FPRINTF(fp,"\t%u\t%u",substring->alignstart_trim - substring->chroffset + 1,
	      substring->alignend_trim - substring->chroffset);
    } else {
      FPRINTF(fp,"\t%u\t%u",substring->alignend_trim - substring->chroffset,
	      substring->alignstart_trim - substring->chroffset + 1);
    }
  } else {
    if (invertp == false) {
      FPRINTF(fp,"\t%u\t%u",substring->alignstart_trim - substring->chroffset,
	      substring->alignend_trim - substring->chroffset + 1);
    } else {
      FPRINTF(fp,"\t%u\t%u",substring->alignend_trim - substring->chroffset + 1,
	      substring->alignstart_trim - substring->chroffset);
    }
  }

  /* field 10: E value */
  FPRINTF(fp,"\t%.2g",blast_evalue(alignlength_trim,substring->nmismatches_bothdiff));

 /* field 11: bit score */
  FPRINTF(fp,"\t%.1f",blast_bitscore(alignlength_trim,substring->nmismatches_bothdiff));
  
  FPRINTF(fp,"\n");

  return;
}



static void
print_forward (Filestring_T fp, char *string, int n) {
  
  FPRINTF(fp,"%.*s",n,string);
  return;
}

static void
print_lc (Filestring_T fp, char *string, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    FPRINTF(fp,"%c",(char) tolower(string[i]));
  }
  return;
}

static void
print_revcomp (Filestring_T fp, char *nt, int len) {

  FPRINTF(fp,"%.*R",len,nt);
  return;
}

static void
print_revcomp_lc (Filestring_T fp, char *nt, int len) {
  int i;

  for (i = len-1; i >= 0; --i) {
    FPRINTF(fp,"%c",(char) tolower(complCode[(int) nt[i]]));
  }
  return;
}


static void
print_genomic (Filestring_T fp, T substring, char *deletion, int deletionlength, bool invertp,
	       Shortread_T queryseq) {
  int i;

  /* Now handled by calculation of queryleft/queryright by Stage3end_display_prep */
  deletion = NULL;
  deletionlength = 0;

  if (invertp == false) {
    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      Shortread_print_subseq_uc(fp,queryseq,substring->querystart,substring->queryend);

    } else if (show_refdiff_p == true) {
      print_forward(fp,substring->genomic_refdiff,substring->queryend);
      if (deletion != NULL) {
	print_lc(fp,deletion,deletionlength);
      }
      print_forward(fp,&(substring->genomic_refdiff[substring->queryend]),substring->querylength - deletionlength - substring->queryend);

    } else {
      print_forward(fp,substring->genomic_bothdiff,substring->queryend);
      if (deletion != NULL) {
	print_lc(fp,deletion,deletionlength);
      }
      print_forward(fp,&(substring->genomic_bothdiff[substring->queryend]),substring->querylength - deletionlength - substring->queryend);
    }

    for (i = 0; i < Shortread_choplength(queryseq); i++) {
      FPRINTF(fp,"*");
    }
    FPRINTF(fp,"\t");
    FPRINTF(fp,"%d..%d",1 + substring->querystart,substring->queryend);

  } else {
    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      Shortread_print_subseq_revcomp_uc(fp,queryseq,substring->querystart,substring->queryend);

    } else if (show_refdiff_p == true) {
      print_revcomp(fp,&(substring->genomic_refdiff[substring->querystart]),substring->querylength - substring->querystart);
      if (deletion != NULL) {
	print_revcomp_lc(fp,deletion,deletionlength);
      }
      print_revcomp(fp,&(substring->genomic_refdiff[deletionlength]),substring->querystart - deletionlength);

    } else {
      print_revcomp(fp,&(substring->genomic_bothdiff[substring->querystart]),substring->querylength - substring->querystart);
      if (deletion != NULL) {
	print_revcomp_lc(fp,deletion,deletionlength);
      }
      print_revcomp(fp,&(substring->genomic_bothdiff[deletionlength]),substring->querystart - deletionlength);
    }
    for (i = 0; i < Shortread_choplength(queryseq); i++) {
      FPRINTF(fp,"*");
    }
    FPRINTF(fp,"\t");
    FPRINTF(fp,"%d..%d",1 + substring->querylength - substring->queryend,
	   substring->querylength - substring->querystart);
  }
  return;
}


static void
print_coordinates (Filestring_T fp, T substring, char *chr, bool invertp) {

  if (substring->plusp == true) {
    if (invertp == false) {
      FPRINTF(fp,"+%s:%u..%u",chr,substring->alignstart_trim - substring->chroffset + 1,
	     substring->alignend_trim - substring->chroffset);
    } else {
      FPRINTF(fp,"-%s:%u..%u",chr,substring->alignend_trim - substring->chroffset,
	     substring->alignstart_trim - substring->chroffset + 1);
    }
  } else {
    if (invertp == false) {
      FPRINTF(fp,"-%s:%u..%u",chr,substring->alignstart_trim - substring->chroffset,
	     substring->alignend_trim - substring->chroffset + 1);
    } else {
      FPRINTF(fp,"+%s:%u..%u",chr,substring->alignend_trim - substring->chroffset + 1,
	     substring->alignstart_trim - substring->chroffset);
    }
  }

  return;
}



void
Substring_print_alignment (Filestring_T fp, Junction_T pre_junction, T substring, Junction_T post_junction,
			   Shortread_T queryseq, Genome_T genome, char *chr, bool invertp) {
  char *deletion_string;
  int deletion_length;
  Junctiontype_T type1, type2;
  Chrpos_T splice_distance_1, splice_distance_2;

  if (post_junction == NULL) {
    deletion_string = (char *) NULL;
    deletion_length = 0;
  } else if (Junction_type(post_junction) != DEL_JUNCTION) {
    deletion_string = (char *) NULL;
    deletion_length = 0;
  } else {
    deletion_string = Junction_deletion_string(post_junction,genome,substring->plusp);
    deletion_length = Junction_nindels(post_junction);
  }

  print_genomic(fp,substring,deletion_string,deletion_length,invertp,queryseq);
  FREE(deletion_string);
  FPRINTF(fp,"\t");
  print_coordinates(fp,substring,chr,invertp);

  FPRINTF(fp,"\t");
  if (pre_junction == NULL) {
    type1 = NO_JUNCTION;
    /* Handle result of substring_trim_novel_spliceends */
    if (invertp == false) {
      if (substring->start_endtype == DON) {
#ifdef NO_COMPARE
	FPRINTF(fp,"donor:NA");
#else
	FPRINTF(fp,"donor:%.2f",substring->siteD_prob);
#endif
      } else if (substring->start_endtype == ACC) {
#ifdef NO_COMPARE
	FPRINTF(fp,"acceptor:NA");
#else
	FPRINTF(fp,"acceptor:%.2f",substring->siteA_prob);
#endif
      } else {
	FPRINTF(fp,"start:%d",substring->querystart);
      }
    } else {
      if (substring->end_endtype == DON) {
#ifdef NO_COMPARE
	FPRINTF(fp,"donor:NA");
#else
	FPRINTF(fp,"donor:%.2f",substring->siteD_prob);
#endif
      } else if (substring->end_endtype == ACC) {
#ifdef NO_COMPARE
	FPRINTF(fp,"acceptor:NA");
#else
	FPRINTF(fp,"acceptor:%.2f",substring->siteA_prob);
#endif
      } else {
	FPRINTF(fp,"start:%d",substring->querylength - substring->queryend);
      }
    }
  } else if ((type1 = Junction_type(pre_junction)) == INS_JUNCTION) {
    FPRINTF(fp,"ins:%d",Junction_nindels(pre_junction));
  } else if (type1 == DEL_JUNCTION) {
    FPRINTF(fp,"del:%d",Junction_nindels(pre_junction));
  } else if (type1 == SPLICE_JUNCTION) {
    if (invertp == false) {
      if (Junction_sensedir(pre_junction) == SENSE_ANTI) {
#ifdef NO_COMPARE
	FPRINTF(fp,"donor:NA",Junction_donor_prob(pre_junction));
#else
	FPRINTF(fp,"donor:%.2f",Junction_donor_prob(pre_junction));
#endif
      } else {
#ifdef NO_COMPARE
	FPRINTF(fp,"acceptor:NA");
#else
	FPRINTF(fp,"acceptor:%.2f",Junction_acceptor_prob(pre_junction));
#endif
      }
    } else {
      if (Junction_sensedir(pre_junction) == SENSE_ANTI) {
#ifdef NO_COMPARE
	FPRINTF(fp,"acceptor:NA");
#else
	FPRINTF(fp,"acceptor:%.2f",Junction_acceptor_prob(pre_junction));
#endif
      } else {
#ifdef NO_COMPARE
	FPRINTF(fp,"donor:NA");
#else
	FPRINTF(fp,"donor:%.2f",Junction_donor_prob(pre_junction));
#endif
      }
    }
  } else if (type1 == CHIMERA_JUNCTION) {
    FPRINTF(fp,"distant:%u",Junction_splice_distance(pre_junction));
  } else {
    abort();
  }

  FPRINTF(fp,"..");

  if (post_junction == NULL) {
    type2 = NO_JUNCTION;
    /* Handle result of substring_trim_novel_spliceends */
    if (invertp == false) {
      if (substring->end_endtype == DON) {
#ifdef NO_COMPARE
	FPRINTF(fp,"donor:NA");
#else
	FPRINTF(fp,"donor:%.2f",substring->siteD_prob);
#endif
      } else if (substring->end_endtype == ACC) {
#ifdef NO_COMPARE
	FPRINTF(fp,"acceptor:NA");
#else
	FPRINTF(fp,"acceptor:%.2f",substring->siteA_prob);
#endif
      } else {
	FPRINTF(fp,"end:%d",substring->querylength - substring->queryend);
      }
    } else {
      if (substring->start_endtype == DON) {
#ifdef NO_COMPARE
	FPRINTF(fp,"donor:NA");
#else
	FPRINTF(fp,"donor:%.2f",substring->siteD_prob);
#endif
      } else if (substring->start_endtype == ACC) {
#ifdef NO_COMPARE
	FPRINTF(fp,"acceptor:NA");
#else
	FPRINTF(fp,"acceptor:%.2f",substring->siteA_prob);
#endif
      } else {
	FPRINTF(fp,"end:%d",substring->querystart);
      }
    }
  } else if ((type2 = Junction_type(post_junction)) == INS_JUNCTION) {
    FPRINTF(fp,"ins:%d",Junction_nindels(post_junction));
  } else if (type2 == DEL_JUNCTION) {
    FPRINTF(fp,"del:%d",Junction_nindels(post_junction));
  } else if (type2 == SPLICE_JUNCTION) {
    if (invertp == false) {
      if (Junction_sensedir(post_junction) == SENSE_ANTI) {
#ifdef NO_COMPARE
	FPRINTF(fp,"acceptor:NA");
#else
	FPRINTF(fp,"acceptor:%.2f",Junction_acceptor_prob(post_junction));
#endif
      } else {
#ifdef NO_COMPARE
	FPRINTF(fp,"donor:NA");
#else
	FPRINTF(fp,"donor:%.2f",Junction_donor_prob(post_junction));
#endif
      }
    } else {
      if (Junction_sensedir(post_junction) == SENSE_ANTI) {
#ifdef NO_COMPARE
	FPRINTF(fp,"donor:NA");
#else
	FPRINTF(fp,"donor:%.2f",Junction_donor_prob(post_junction));
#endif
      } else {
#ifdef NO_COMPARE
	FPRINTF(fp,"acceptor:NA");
#else
	FPRINTF(fp,"acceptor:%.2f",Junction_acceptor_prob(post_junction));
#endif
      }
    }
  } else if (type2 == CHIMERA_JUNCTION) {
    FPRINTF(fp,"distant:%u",Junction_splice_distance(post_junction));
  } else {
    abort();
  }


  FPRINTF(fp,",matches:%d,sub:%d",Substring_nmatches_to_trims(substring),substring->nmismatches_bothdiff);
  if (print_nsnpdiffs_p) {
    FPRINTF(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }

  if (type1 == SPLICE_JUNCTION && type2 == SPLICE_JUNCTION) {
    if (invertp == false) {
      if (Junction_sensedir(pre_junction) == SENSE_FORWARD) {
	FPRINTF(fp,",dir:sense");
      } else if (Junction_sensedir(pre_junction) == SENSE_ANTI) {
	FPRINTF(fp,",dir:antisense");
      } else {
	FPRINTF(fp,",dir:unknown");
      }
    } else {
      if (Junction_sensedir(pre_junction) == SENSE_FORWARD) {
	FPRINTF(fp,",dir:antisense");
      } else if (Junction_sensedir(pre_junction) == SENSE_ANTI) {
	FPRINTF(fp,",dir:sense");
      } else {
	FPRINTF(fp,",dir:unknown");
      }
    }
    splice_distance_1 = Junction_splice_distance(pre_junction);
    splice_distance_2 = Junction_splice_distance(post_junction);
    if (splice_distance_1 == 0 && splice_distance_2 == 0) {
      /* Skip */
    } else if (splice_distance_1 == 0) {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_2:%u",splice_distance_2);
    } else if (splice_distance_2 == 0) {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_1:%u",splice_distance_1);
    } else {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_1:%u",splice_distance_1);
      FPRINTF(fp,",splice_dist_2:%u",splice_distance_2);
    }

  } else if (type1 == SPLICE_JUNCTION) {
    if (invertp == false) {
      if (Junction_sensedir(pre_junction) == SENSE_FORWARD) {
	FPRINTF(fp,",dir:sense");
      } else if (Junction_sensedir(pre_junction) == SENSE_ANTI) {
	FPRINTF(fp,",dir:antisense");
      } else {
	FPRINTF(fp,",dir:unknown");
      }
    } else {
      if (Junction_sensedir(pre_junction) == SENSE_FORWARD) {
	FPRINTF(fp,",dir:antisense");
      } else if (Junction_sensedir(pre_junction) == SENSE_ANTI) {
	FPRINTF(fp,",dir:sense");
      } else {
	FPRINTF(fp,",dir:unknown");
      }
    }
    if ((splice_distance_1 = Junction_splice_distance(pre_junction)) > 0) {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_1:%u",splice_distance_1);
    }

  } else if (type2 == SPLICE_JUNCTION) {
    if (invertp == false) {
      if (Junction_sensedir(post_junction) == SENSE_FORWARD) {
	FPRINTF(fp,",dir:sense");
      } else if (Junction_sensedir(post_junction) == SENSE_ANTI) {
	FPRINTF(fp,",dir:antisense");
      } else {
	FPRINTF(fp,",dir:unknown");
      }
    } else {
      if (Junction_sensedir(post_junction) == SENSE_FORWARD) {
	FPRINTF(fp,",dir:antisense");
      } else if (Junction_sensedir(post_junction) == SENSE_ANTI) {
	FPRINTF(fp,",dir:sense");
      } else {
	FPRINTF(fp,",dir:unknown");
      }
    }
    if ((splice_distance_2 = Junction_splice_distance(post_junction)) > 0) {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_2:%u",splice_distance_2);
    }
  }

  return;
}


/************************************************************************
 *   Multiclean
 ************************************************************************/

static char *
get_total_tally (long int *tally, char *ptr) {
  int n;
  char *end;

  if ((end = index(ptr,'\n')) == NULL) {
    fprintf(stderr,"Premature end of line %s\n",ptr);
    return 0;
  }
  /* fprintf(stderr,"Getting tally for %.*s\n",end-ptr,ptr); */

  while (ptr < end) {
    while (ptr < end && !isdigit((int) *ptr)) {
      ptr++;
    }
    if (ptr < end) {
      sscanf(ptr,"%d",&n);
#if 0
      debug(if (n > 0) printf(" %d",n));
#endif
      (*tally) += n;
      while (ptr < end && !isspace(*ptr)) {
	ptr++;
      }
      while (ptr < end && isspace(*ptr)) {
	ptr++;
      }
    }
  }

  return ptr;
}

long int
Substring_tally (T this, IIT_T tally_iit, int *tally_divint_crosstable) {
  long int total = 0U;
  Interval_T interval;
  char *annotation, *restofheader, *ptr;
#if 0
  bool alloc_chr_p;
#endif
  bool allocp;
  Chrpos_T chrpos, intervalend;

  Chrpos_T coordstart, coordend, pos5, pos3;
  int *matches;
  int nmatches, i;
  
  pos5 = this->alignstart_trim - this->chroffset;
  pos3 = this->alignend_trim - this->chroffset;

  if (pos5 < pos3) {
    coordstart = pos5;
    coordend = pos3;
  } else {
    coordstart = pos3;
    coordend = pos5;
  }
  coordstart += 1;		/* Because tally IIT is 1-based */
  debug(printf("coordstart = %u, coordend = %u\n",coordstart,coordend));

#if 0
  chr = Univ_IIT_label(chromosome_iit,this->chrnum,&alloc_chr_p);
#endif
  matches = IIT_get_with_divno(&nmatches,tally_iit,tally_divint_crosstable[this->chrnum],
			       coordstart,coordend,/*sortp*/false);

  for (i = 0; i < nmatches; i++) {
    annotation = IIT_annotation(&restofheader,tally_iit,matches[i],&allocp);

    interval = IIT_interval(tally_iit,matches[i]);
    chrpos = Interval_low(interval);
    intervalend = Interval_high(interval);

    ptr = annotation;

    while (chrpos < coordstart) {
      if ((ptr = index(ptr,'\n')) == NULL) {
	fprintf(stderr,"Premature end of tally from %u to %u\n",
		Interval_low(interval),Interval_high(interval));
	return total;
      } else {
	ptr++;
      }
      chrpos++;
    }

    while (chrpos <= intervalend && chrpos <= coordend) {
      ptr = get_total_tally(&total,ptr);
      ptr++;
      chrpos++;
    }

    if (allocp == true) {
      FREE(restofheader);
    }
  }

  FREE(matches);

#if 0
  if (alloc_chr_p) {
    FREE(chr);
  }
#endif

  debug(printf("Subtotal = %ld\n",total));
  return total;
}


bool
Substring_runlength_p (T this, IIT_T runlength_iit, int *runlength_divint_crosstable) {
  Chrpos_T coordstart, coordend, pos5, pos3;
  
  pos5 = this->alignstart_trim - this->chroffset;
  pos3 = this->alignend_trim - this->chroffset;

  if (pos5 < pos3) {
    coordstart = pos5;
    coordend = pos3;
  } else {
    coordstart = pos3;
    coordend = pos5;
  }
  coordstart += 1;		/* Because runlength IIT is 1-based */
  debug(printf("coordstart = %u, coordend = %u\n",coordstart,coordend));

  /* chr = Univ_IIT_label(chromosome_iit,this->chrnum,&alloc_chr_p); */
  return IIT_exists_with_divno(runlength_iit,runlength_divint_crosstable[this->chrnum],
				coordstart,coordend);
}


int
Substring_count_mismatches_region (T this, int trim_querystart, int trim_queryend) {
  int left_bound, right_bound;

  if (this == NULL) {
    return 0;
  } else {
    left_bound = trim_querystart;
    right_bound = this->querylength - trim_queryend;

    if (this->queryend_pretrim < left_bound) {
      return 0;
    } else if (this->querystart_pretrim > right_bound) {
      return 0;
    } else {
      if (this->querystart_pretrim > left_bound) {
	left_bound = this->querystart_pretrim;
      }
      if (this->queryend_pretrim < right_bound) {
	right_bound = this->queryend_pretrim;
      }
      
#if 0
      /* Not sure why we have this assertion.  Values set by eventrim procedure in Stage3end_optimal_score_prefinal.  If assertion fails, then just return 0 */
      assert(left_bound <= right_bound);
#endif

      if (left_bound > right_bound) {
	return 0;
      } else if (this->plusp) {
	return Genome_count_mismatches_substring(genomebits,genomebits_alt,this->query_compress,this->left,
						 /*pos5*/left_bound,/*pos3*/right_bound,/*plusp*/true,this->genestrand);
      } else {
	return Genome_count_mismatches_substring(genomebits,genomebits_alt,this->query_compress,this->left,
						 /*pos5*/this->querylength - right_bound,
						 /*pos3*/this->querylength - left_bound,
						 /*plusp*/false,this->genestrand);
      }
    }
  }
}


/************************************************************************
 *   Conversion to Pair_T format
 ************************************************************************/

#ifdef RESOLVE_INSIDE_GENERAL
List_T
Substring_convert_to_pairs (List_T pairs, T substring, char *queryuc_ptr,
			    Chrpos_T chrlength, Pairpool_T pairpool) {
  int querystart, queryend, querypos, i;
  Chrpos_T chrpos;
  char genome, cdna;

  if (substring == NULL) {
    return pairs;
  } else if (substring->alts_ncoords > 0) {
    return pairs;
  }

  debug6(printf("*** Entered Substring_convert_to_pairs with querylength %d\n",querylength));

  if (substring->plusp == true) {
    querystart = substring->querystart;
    queryend = substring->queryend;

    /* Pairs are all zero-based, so do not add 1 */
#if 0
    chrpos = substring->genomicstart_adj + querystart - substring->chroffset /*+ 1*/;
#else
    chrpos = substring->genomicstart + querystart - substring->chroffset /*+ 1*/;
#endif

    debug6(printf("plus conversion\n"));
    debug6(printf("querystart %d, queryend %d, plusp %d\n",querystart,queryend,substring->plusp));
    debug6(printf("alignstart %u, alignend %u\n",substring->alignstart_trim - substring->chroffset,
		  substring->alignend_trim - substring->chroffset));
    debug6(printf("chrpos %u\n",chrpos));

    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      for (i = querystart, querypos = /*queryseq_offset +*/ querystart; i < queryend; i++, querypos++) {
	pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
			      queryuc_ptr[i],/*comp*/MATCH_COMP,queryuc_ptr[i],/*g_alt*/queryuc_ptr[i],/*dynprogindex*/0);
      }
    } else if (show_refdiff_p == true) {
      for (i = querystart, querypos = /*queryseq_offset +*/ querystart; i < queryend; i++, querypos++) {
	if (isupper(genome = substring->genomic_refdiff[i])) {
	  /* assert(queryuc_ptr[i] == genome || queryuc_ptr[i] == 'N'); -- Doesn't hold for SNPs */
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				queryuc_ptr[i],/*comp*/MATCH_COMP,genome,/*g_alt*/genome,/*dynprogindex*/0);
	} else if ((cdna = queryuc_ptr[i]) == 'N' && genome == 'n') {
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				cdna,/*comp*/MATCH_COMP,toupper(genome),/*g_alt*/toupper(genome),/*dynprogindex*/0);
#if 0
	} else if (cdna == 'N') {
	  /* Use query_unk_mismatch_p ? */
	} else if (genome == 'n') {
	  /* Use genome_unk_mismatch_p ? */
#endif
	} else {
	  assert(queryuc_ptr[i] != toupper(genome));
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				cdna,/*comp*/MISMATCH_COMP,toupper(genome),/*g_alt*/toupper(genome),
				/*dynprogindex*/0);
	}
      }
    } else {
      /* printf("querystart %d, queryend %d\n",querystart,queryend); */
      /* printf("seq1   %s\n",queryuc_ptr); */
      /* printf("genome %s\n",substring->genomic_bothdiff); */
      for (i = querystart, querypos = /*queryseq_offset +*/ querystart; i < queryend; i++, querypos++) {
	if (isupper(genome = substring->genomic_bothdiff[i])) {
	  /* assert(queryuc_ptr[i] == genome || queryuc_ptr[i] == 'N'); -- Doesn't hold for SNPs */
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				queryuc_ptr[i],/*comp*/MATCH_COMP,genome,/*g_alt*/genome,/*dynprogindex*/0);
	} else if ((cdna = queryuc_ptr[i]) == 'N' && genome == 'n') {
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				cdna,/*comp*/MATCH_COMP,toupper(genome),/*g_alt*/toupper(genome),/*dynprogindex*/0);
#if 0
	} else if (cdna == 'N') {
	  /* Use query_unk_mismatch_p ? */
	} else if (genome == 'n') {
	  /* Use genome_unk_mismatch_p ? */
#endif
	} else {
	  assert(queryuc_ptr[i] != toupper(genome));
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				cdna,/*comp*/MISMATCH_COMP,toupper(genome),/*g_alt*/toupper(genome),
				/*dynprogindex*/0);
	}
      }
    }

  } else {
    querystart = substring->querystart;
    queryend = substring->queryend;

    /* For minus, to get 0-based coordinates, subtract 1 */
#if 0
    chrpos = substring->genomicstart_adj - querystart - substring->chroffset - 1;
#else
    chrpos = substring->genomicstart - querystart - substring->chroffset - 1;
    chrpos = chrlength - chrpos;
#endif

    debug6(printf("minus conversion\n"));
    debug6(printf("querystart %d, queryend %d, plusp %d\n",querystart,queryend,substring->plusp));
    debug6(printf("alignstart %u, alignend %u\n",substring->alignstart_trim - substring->chroffset,
		  substring->alignend_trim - substring->chroffset));
    debug6(printf("chrpos %u\n",chrpos));

    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      for (i = querystart, querypos = /*queryseq_offset +*/ querystart; i < queryend; i++, querypos++) {
	pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
			      queryuc_ptr[i],/*comp*/MATCH_COMP,queryuc_ptr[i],/*g_alt*/queryuc_ptr[i],/*dynprogindex*/0);
      }
    } else if (show_refdiff_p == true) {
      for (i = querystart, querypos = /*queryseq_offset +*/ querystart; i < queryend; i++, querypos++) {
	if (isupper(genome = substring->genomic_refdiff[i])) {
	  assert(queryuc_ptr[i] == genome || queryuc_ptr[i] == 'N');
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				queryuc_ptr[i],/*comp*/MATCH_COMP,genome,/*g_alt*/genome,/*dynprogindex*/0);
	} else if ((cdna = queryuc_ptr[i]) == 'N' && genome == 'n') {
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				cdna,/*comp*/MATCH_COMP,toupper(genome),/*g_alt*/toupper(genome),/*dynprogindex*/0);
#if 0
	} else if (cdna == 'N') {
	  /* Use query_unk_mismatch_p ? */
	} else if (genome == 'n') {
	  /* Use genome_unk_mismatch_p ? */
#endif
	} else {
	  assert(queryuc_ptr[i] != toupper(genome));
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				cdna,/*comp*/MISMATCH_COMP,toupper(genome),/*g_alt*/toupper(genome),
				/*dynprogindex*/0);
	}
      }
    } else {
      for (i = querystart, querypos = /*queryseq_offset +*/ querystart; i < queryend; i++, querypos++) {
	if (isupper(genome = substring->genomic_bothdiff[i])) {
	  /* assert(queryuc_ptr[i] == genome || queryuc_ptr[i] == 'N'); */
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				queryuc_ptr[i],/*comp*/MATCH_COMP,genome,/*g_alt*/genome,/*dynprogindex*/0);
	} else if ((cdna = queryuc_ptr[i]) == 'N' && genome == 'n') {
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				cdna,/*comp*/MATCH_COMP,toupper(genome),/*g_alt*/toupper(genome),/*dynprogindex*/0);
#if 0
	} else if (cdna == 'N') {
	  /* Use query_unk_mismatch_p ? */
	} else if (genome == 'n') {
	  /* Use genome_unk_mismatch_p ? */
#endif
	} else {
	  assert(queryuc_ptr[i] != toupper(genome));
	  pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrpos++,
				cdna,/*comp*/MISMATCH_COMP,toupper(genome),/*g_alt*/toupper(genome),
				/*dynprogindex*/0);
	}
      }
    }
  }

  debug6(Simplepair_dump_list(pairs,true));
  return pairs;
}
#endif



List_T
Substring_convert_to_pairs_out (List_T pairs, T substring, int querylength, Shortread_T queryseq,
				int hardclip_low, int hardclip_high, int queryseq_offset) {
  int querystart, queryend, querypos, i;
  Chrpos_T chrpos;
  char *seq1;
  char genome;

  if (substring == NULL) {
    return pairs;
  } else if (substring->alts_ncoords > 0) {
    return pairs;
  }

  debug6(printf("*** Entered Substring_convert_to_pairs with querylength %d, hardclip_low %d, hardclip_high %d, queryseq_offset %d\n",
		querylength,hardclip_low,hardclip_high,queryseq_offset));

  seq1 = Shortread_fullpointer_uc(queryseq);
  if (substring->plusp == true) {
    if (hardclip_low > substring->querystart) {
      querystart = hardclip_low;
    } else {
      querystart = substring->querystart;
    }

    if (querylength - hardclip_high < substring->queryend) {
      queryend = querylength - hardclip_high;
    } else {
      queryend = substring->queryend;
    }
    /* Pairs are all zero-based, so do not add 1 */
#if 0
    chrpos = substring->genomicstart_adj + querystart - substring->chroffset /*+ 1*/;
#else
    chrpos = substring->genomicstart + querystart - substring->chroffset /*+ 1*/;
#endif

    debug6(printf("plus conversion\n"));
    debug6(printf("querystart %d, queryend %d, plusp %d\n",querystart,queryend,substring->plusp));
    debug6(printf("alignstart %u, alignend %u\n",substring->alignstart_trim - substring->chroffset,
		  substring->alignend_trim - substring->chroffset));
    debug6(printf("chrpos %u\n",chrpos));

    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++, querypos++) {
	pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos++,
								seq1[i],/*comp*/MATCH_COMP,seq1[i]));
      }
    } else if (show_refdiff_p == true) {
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++, querypos++) {
	if (isupper(genome = substring->genomic_refdiff[i])) {
	  assert(seq1[i] == genome || seq1[i] == 'N');
	  pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos++,
								  seq1[i],/*comp*/MATCH_COMP,genome));
	} else {
	  assert(seq1[i] != toupper(genome));
	  pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos++,
								  seq1[i],/*comp*/MISMATCH_COMP,toupper(genome)));
	}
      }
    } else {
      /* printf("querystart %d, queryend %d\n",querystart,queryend); */
      /* printf("seq1   %s\n",seq1); */
      /* printf("genome %s\n",substring->genomic_bothdiff); */
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++, querypos++) {
	if (isupper(genome = substring->genomic_bothdiff[i])) {
	  assert(seq1[i] == genome || seq1[i] == 'N');
	  pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos++,
								  seq1[i],/*comp*/MATCH_COMP,genome));
	} else {
	  assert(seq1[i] != toupper(genome));
	  pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos++,
								  seq1[i],/*comp*/MISMATCH_COMP,toupper(genome)));
	}
      }
    }

  } else {
    if (hardclip_high > substring->querystart) {
      querystart = hardclip_high;
    } else {
      querystart = substring->querystart;
    }

    if (querylength - hardclip_low < substring->queryend) {
      queryend = querylength - hardclip_low;
    } else {
      queryend = substring->queryend;
    }
    /* For minus, to get 0-based coordinates, subtract 1 */
#if 0
    chrpos = substring->genomicstart_adj - querystart - substring->chroffset - 1;
#else
    chrpos = substring->genomicstart - querystart - substring->chroffset - 1;
#endif

    debug6(printf("minus conversion\n"));
    debug6(printf("querystart %d, queryend %d, plusp %d\n",querystart,queryend,substring->plusp));
    debug6(printf("alignstart %u, alignend %u\n",substring->alignstart_trim - substring->chroffset,
		  substring->alignend_trim - substring->chroffset));
    debug6(printf("chrpos %u\n",chrpos));

    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++, querypos++) {
	pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos--,
								seq1[i],/*comp*/MATCH_COMP,seq1[i]));
      }
    } else if (show_refdiff_p == true) {
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++, querypos++) {
	if (isupper(genome = substring->genomic_refdiff[i])) {
	  assert(seq1[i] == genome || seq1[i] == 'N');
	  pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos--,
								  seq1[i],/*comp*/MATCH_COMP,genome));
	} else {
	  assert(seq1[i] != toupper(genome));
	  pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos--,
								  seq1[i],/*comp*/MISMATCH_COMP,toupper(genome)));
	}
      }
    } else {
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++, querypos++) {
	if (isupper(genome = substring->genomic_bothdiff[i])) {
	  /* assert(seq1[i] == genome || seq1[i] == 'N'); */
	  pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos--,
								  seq1[i],/*comp*/MATCH_COMP,genome));
	} else {
	  /* assert(seq1[i] != toupper(genome)); */
	  pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrpos--,
								  seq1[i],/*comp*/MISMATCH_COMP,toupper(genome)));
	}
      }
    }
  }

  debug6(Simplepair_dump_list(pairs,true));
  return pairs;
}


#ifdef RESOLVE_INSIDE_GENERAL
List_T
Substring_add_insertion (List_T pairs, T substringA, T substringB,
			 int insertionlength, char *queryuc_ptr,
			 Pairpool_T pairpool) {
  int querystartA, queryendA, querystartB, queryendB, querypos, i;
  Chrpos_T chrendA;


  if (substringA->plusp == true) {
    querystartA = substringA->querystart;
    queryendA = substringA->queryend;
    querystartB = substringB->querystart;
    queryendB = substringB->queryend;

    /* Pairs are all zero-based, so do not add 1 */
#if 0
    chrendA = substringA->genomicstart_adj + queryendA - substringA->chroffset /*+ 1*/;
#else
    chrendA = substringA->genomicstart + queryendA - substringA->chroffset /*+ 1*/;
#endif

  } else {
    querystartA = substringA->querystart;
    queryendA = substringA->queryend;
    querystartB = substringB->querystart;
    queryendB = substringB->queryend;

    /* Pairs are all zero-based, so subtract 1 */
#if 0
    chrendA = substringA->genomicstart_adj - queryendA - substringA->chroffset - 1;
#else
    chrendA = substringA->genomicstart - queryendA - substringA->chroffset - 1;
#endif
  }

  if (querystartA <= queryendA && querystartB <= queryendB) {
    querypos = queryendA /*+ queryseq_offset*/;
    i = queryendA;
    while (--insertionlength >= 0) {
      pairs = Pairpool_push(pairs,pairpool,querypos++,/*genomepos*/chrendA,
			    queryuc_ptr[i++],/*comp*/INDEL_COMP,' ',/*g_alt*/' ',/*dynprogindex*/0);
    }
  }

  return pairs;
}
#endif


#ifdef RESOLVE_INSIDE_GENERAL
List_T
Substring_add_deletion (List_T pairs, T substringA, T substringB,
			char *deletion, int deletionlength,
			Pairpool_T pairpool) {
  int querystartA, queryendA, querystartB, queryendB, querypos, k;
  Chrpos_T chrendA;

  if (substringA->plusp == true) {
    querystartA = substringA->querystart;
    queryendA = substringA->queryend;
    querystartB = substringB->querystart;
    queryendB = substringB->queryend;

    /* Pairs are all zero-based, so do not add 1 */
#if 0
    chrendA = substringA->genomicstart_adj + queryendA - substringA->chroffset /*+ 1*/;
#else
    chrendA = substringA->genomicstart + queryendA - substringA->chroffset /*+ 1*/;
#endif

    if (querystartA < queryendA && querystartB < queryendB) {
      querypos = queryendA /*+ queryseq_offset*/;
      for (k = 0; k < deletionlength; k++) {
	pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrendA++,
			      ' ',/*comp*/INDEL_COMP,deletion[k],/*g_alt*/deletion[k],
			      /*dynprogindex*/0);
      }
    }

  } else {
    querystartA = substringA->querystart;
    queryendA = substringA->queryend;
    querystartB = substringB->querystart;
    queryendB = substringB->queryend;

    /* Pairs are all zero-based, so subtract 1 */
#if 0
    chrendA = substringA->genomicstart_adj - queryendA - substringA->chroffset - 1;
#else
    chrendA = substringA->genomicstart - queryendA - substringA->chroffset - 1;
#endif

    if (querystartA <= queryendA && querystartB <= queryendB) {
      querypos = queryendA /*+ queryseq_offset*/;
      for (k = 0; k < deletionlength; k++) {
	pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrendA++,
			      ' ',/*comp*/INDEL_COMP,deletion[k],/*g_alt*/deletion[k],
			      /*dynprogindex*/0);
      }
    }
  }

  return pairs;
}
#endif


#ifdef RESOLVE_INSIDE_GENERAL
List_T
Substring_add_intron (List_T pairs, T substringA, T substringB, Pairpool_T pairpool) {
  int querystartA, queryendA, querystartB, queryendB, querypos;
  Chrpos_T chrendA;

  if (substringA->plusp == true) {
    querystartA = substringA->querystart;
    queryendA = substringA->queryend;
    querystartB = substringB->querystart;
    queryendB = substringB->queryend;

    /* Pairs are all zero-based, so do not add 1 */
#if 0
    chrendA = substringA->genomicstart_adj + queryendA - substringA->chroffset /*+ 1*/;
#else
    chrendA = substringA->genomicstart + queryendA - substringA->chroffset /*+ 1*/;
#endif

  } else {
    querystartA = substringA->querystart;
    queryendA = substringA->queryend;
    querystartB = substringB->querystart;
    queryendB = substringB->queryend;


    /* Pairs are all zero-based, so subtract 1 */
#if 0
    chrendA = substringA->genomicstart_adj - queryendA - substringA->chroffset - 1;
#else
    chrendA = substringA->genomicstart - queryendA - substringA->chroffset - 1;
#endif
  }

  if (querystartA <= queryendA && querystartB <= queryendB) {
    /* Add gapholder */
    /* All we really need for Pair_print_sam is to set gapp to be true */
    querypos = queryendA /*+ queryseq_offset*/;
    pairs = Pairpool_push(pairs,pairpool,querypos,/*genomepos*/chrendA,
			  ' ',/*comp*/FWD_CANONICAL_INTRON_COMP,' ',/*g_alt*/' ',
			  /*dynprogindex*/0);
  }

  return pairs;
}
#endif



List_T
Substring_add_insertion_out (List_T pairs, T substringA, T substringB, int querylength,
			     int insertionlength, Shortread_T queryseq,
			     int hardclip_low, int hardclip_high, int queryseq_offset) {
  int querystartA, queryendA, querystartB, queryendB, querypos, i;
  Chrpos_T chrendA;
  char *seq1;


  if (substringA->plusp == true) {
    if (hardclip_low > substringA->querystart) {
      querystartA = hardclip_low;
    } else {
      querystartA = substringA->querystart;
    }

    if (querylength - hardclip_high < substringA->queryend) {
      queryendA = querylength - hardclip_high;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_low > substringB->querystart) {
      querystartB = hardclip_low;
    } else {
      querystartB = substringB->querystart;
    }

    if (querylength - hardclip_high < substringB->queryend) {
      queryendB = querylength - hardclip_high;
    } else {
      queryendB = substringB->queryend;
    }

    /* Pairs are all zero-based, so do not add 1 */
#if 0
    chrendA = substringA->genomicstart_adj + queryendA - substringA->chroffset /*+ 1*/;
#else
    chrendA = substringA->genomicstart + queryendA - substringA->chroffset /*+ 1*/;
#endif

  } else {
    if (hardclip_high > substringA->querystart) {
      querystartA = hardclip_high;
    } else {
      querystartA = substringA->querystart;
    }

    if (querylength - hardclip_low < substringA->queryend) {
      queryendA = querylength - hardclip_low;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_high > substringB->querystart) {
      querystartB = hardclip_high;
    } else {
      querystartB = substringB->querystart;
    }

    if (querylength - hardclip_low < substringB->queryend) {
      queryendB = querylength - hardclip_low;
    } else {
      queryendB = substringB->queryend;
    }

    /* Pairs are all zero-based, so subtract 1 */
#if 0
    chrendA = substringA->genomicstart_adj - queryendA - substringA->chroffset - 1;
#else
    chrendA = substringA->genomicstart - queryendA - substringA->chroffset - 1;
#endif
  }

  if (querystartA <= queryendA && querystartB <= queryendB) {
    seq1 = Shortread_fullpointer_uc(queryseq);
    querypos = queryendA + queryseq_offset;
    i = queryendA;
    while (--insertionlength >= 0) {
      pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos++,/*genomepos*/chrendA,
							      seq1[i++],/*comp*/INDEL_COMP,' '));
    }
  }

  return pairs;
}


List_T
Substring_add_deletion_out (List_T pairs, T substringA, T substringB, int querylength,
			    char *deletion, int deletionlength,
			    int hardclip_low, int hardclip_high, int queryseq_offset) {
  int querystartA, queryendA, querystartB, queryendB, querypos, k;
  Chrpos_T chrendA;

  if (substringA->plusp == true) {
    if (hardclip_low > substringA->querystart) {
      querystartA = hardclip_low;
    } else {
      querystartA = substringA->querystart;
    }

    if (querylength - hardclip_high < substringA->queryend) {
      queryendA = querylength - hardclip_high;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_low > substringB->querystart) {
      querystartB = hardclip_low;
    } else {
      querystartB = substringB->querystart;
    }

    if (querylength - hardclip_high < substringB->queryend) {
      queryendB = querylength - hardclip_high;
    } else {
      queryendB = substringB->queryend;
    }

    /* Pairs are all zero-based, so do not add 1 */
#if 0
    chrendA = substringA->genomicstart_adj + queryendA - substringA->chroffset /*+ 1*/;
#else
    chrendA = substringA->genomicstart + queryendA - substringA->chroffset /*+ 1*/;
#endif

    if (querystartA < queryendA && querystartB < queryendB) {
      querypos = queryendA + queryseq_offset;
      for (k = 0; k < deletionlength; k++) {
	pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrendA++,
								' ',/*comp*/INDEL_COMP,deletion[k]));
      }
    }

  } else {
    if (hardclip_high > substringA->querystart) {
      querystartA = hardclip_high;
    } else {
      querystartA = substringA->querystart;
    }

    if (querylength - hardclip_low < substringA->queryend) {
      queryendA = querylength - hardclip_low;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_high > substringB->querystart) {
      querystartB = hardclip_high;
    } else {
      querystartB = substringB->querystart;
    }

    if (querylength - hardclip_low < substringB->queryend) {
      queryendB = querylength - hardclip_low;
    } else {
      queryendB = substringB->queryend;
    }

    /* Pairs are all zero-based, so subtract 1 */
#if 0
    chrendA = substringA->genomicstart_adj - queryendA - substringA->chroffset - 1;
#else
    chrendA = substringA->genomicstart - queryendA - substringA->chroffset - 1;
#endif

    if (querystartA <= queryendA && querystartB <= queryendB) {
      querypos = queryendA + queryseq_offset;
      for (k = 0; k < deletionlength; k++) {
	pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrendA--,
								' ',/*comp*/INDEL_COMP,deletion[k]));
      }
    }
  }

  return pairs;
}



List_T
Substring_add_intron_out (List_T pairs, T substringA, T substringB, int querylength,
			  int hardclip_low, int hardclip_high, int queryseq_offset) {
  int querystartA, queryendA, querystartB, queryendB, querypos;
  Chrpos_T chrendA;

  if (substringA->plusp == true) {
    if (hardclip_low > substringA->querystart) {
      querystartA = hardclip_low;
    } else {
      querystartA = substringA->querystart;
    }

    if (querylength - hardclip_high < substringA->queryend) {
      queryendA = querylength - hardclip_high;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_low > substringB->querystart) {
      querystartB = hardclip_low;
    } else {
      querystartB = substringB->querystart;
    }

    if (querylength - hardclip_high < substringB->queryend) {
      queryendB = querylength - hardclip_high;
    } else {
      queryendB = substringB->queryend;
    }

    /* Pairs are all zero-based, so do not add 1 */
#if 0
    chrendA = substringA->genomicstart_adj + queryendA - substringA->chroffset /*+ 1*/;
#else
    chrendA = substringA->genomicstart + queryendA - substringA->chroffset /*+ 1*/;
#endif

  } else {
    if (hardclip_high > substringA->querystart) {
      querystartA = hardclip_high;
    } else {
      querystartA = substringA->querystart;
    }

    if (querylength - hardclip_low < substringA->queryend) {
      queryendA = querylength - hardclip_low;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_high > substringB->querystart) {
      querystartB = hardclip_high;
    } else {
      querystartB = substringB->querystart;
    }

    if (querylength - hardclip_low < substringB->queryend) {
      queryendB = querylength - hardclip_low;
    } else {
      queryendB = substringB->queryend;
    }

    /* Pairs are all zero-based, so subtract 1 */
#if 0
    chrendA = substringA->genomicstart_adj - queryendA - substringA->chroffset - 1;
#else
    chrendA = substringA->genomicstart - queryendA - substringA->chroffset - 1;
#endif
  }

  if (querystartA <= queryendA && querystartB <= queryendB) {
    /* Add gapholder */
    /* All we really need for Pair_print_sam is to set gapp to be true */
    querypos = queryendA + queryseq_offset;
    pairs = List_push_out(pairs,(void *) Simplepair_new_out(querypos,/*genomepos*/chrendA,
							    ' ',/*comp*/FWD_CANONICAL_INTRON_COMP,' '));
  }

  return pairs;
}



void
Substring_setup (bool print_nsnpdiffs_p_in, bool print_snplabels_p_in,
		 bool show_refdiff_p_in, IIT_T snps_iit_in, int *snps_divint_crosstable_in,
		 Genome_T genomebits_in, Genome_T genomebits_alt_in,
		 Univ_IIT_T chromosome_iit_in, int nchromosomes_in, int circular_typeint_in,
		 IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		 IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
		 int donor_typeint_in, int acceptor_typeint_in,
		 bool novelsplicingp_in, bool knownsplicingp_in,
		 Outputtype_T output_type_in, Mode_T mode_in, Univcoord_T genomelength_in) {
  print_nsnpdiffs_p = print_nsnpdiffs_p_in;
  print_snplabels_p = print_snplabels_p_in;
  show_refdiff_p = show_refdiff_p_in;

  snps_iit = snps_iit_in;
  snps_divint_crosstable = snps_divint_crosstable_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

  chromosome_iit = chromosome_iit_in;
  nchromosomes = nchromosomes_in;
  circular_typeint = circular_typeint_in;

  genes_iit = genes_iit_in;
  genes_divint_crosstable = genes_divint_crosstable_in;

  splicesites_iit = splicesites_iit_in;
  splicesites_divint_crosstable = splicesites_divint_crosstable_in;

  donor_typeint = donor_typeint_in;
  acceptor_typeint = acceptor_typeint_in;

  novelsplicingp = novelsplicingp_in;
  knownsplicingp = knownsplicingp_in;

  output_type = output_type_in;
  mode = mode_in;

  genomelength = (double) genomelength_in;

  return;
}



