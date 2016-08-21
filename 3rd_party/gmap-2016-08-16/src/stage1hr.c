static char rcsid[] = "$Id: stage1hr.c 196433 2016-08-16 20:20:51Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "stage1hr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset() */
#include <math.h>
#include <ctype.h>		/* for tolower() */
#include "assert.h"
#include "mem.h"
#include "reader.h"
#include "oligo.h"
#include "indexdb.h"
#include "indexdb_hr.h"
#include "list.h"
#include "listdef.h"
#include "intlist.h"
#include "intlistdef.h"
#include "splice.h"
#include "indel.h"
#include "stage3hr.h"
#include "substring.h"
#include "complement.h"
#include "compress.h"
#include "genome128_hr.h"
#include "genome_sites.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "iitdef.h"
#include "univinterval.h"
#ifdef LARGE_GENOMES
#include "uint8list.h"
#else
#include "uintlist.h"
#include "sarray-read.h"
#endif

#include "spanningelt.h"
#include "cmet.h"
#include "atoi.h"

#include "stage2.h"
#include "stage3.h"
#include "comp.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#endif


#ifdef HAVE_64_BIT
#ifdef LARGE_GENOMES
#else
#define DIAGONAL_ADD_QUERYPOS 1
#endif
#endif

#define USE_HEAP 1

#define SPEED 1

/* Note: Heapsort still turns out to be a bit faster than a global
   qsort, because it takes advantage of the fact that the positions
   within each batch are already sorted.  Also, heapsort can handle
   8-byte positions. */

/* Not useful to have LONG_ENDSPLICES, which allow for two pairmax
   distances, because long insert lengths will get pruned anyway by
   Stage3pair_new.  In future, could have LONG_ENDSPLICES, but require
   inside alignment to start within one pairmax */
/* #define LONG_ENDSPLICES 1 */	/* Necessary to get outside splices correctly */

#define NO_EXTENSIONS_BEFORE_ZERO 1

#define ALLOW_MIDDLE_ALIGNMENTS 1

/* #define EXTRACT_GENOMICSEG 1 */
#ifdef EXTRACT_GENOMICSEG
#define MAX_INDEXSIZE 8
#endif


/* MAX_NALIGNMENTS of 2 vs 1 gets 1600 improvements in 275,000 reads */
/* MAX_NALIGNMENTS of 3 vs 2 gets 96 improvements in 275,000 reads */
#define MAX_NALIGNMENTS 3

#define MAX_NTERMINALS 100
#define MAX_ALLOCATION 200

#define PAIRMAX_ADDITIONAL 10000 /* Allows for finding of unpaired GMAP alignments beyond pairmax */

static bool use_sarray_p = true;
static bool use_only_sarray_p = true;


/* Mode */
static Mode_T mode;
static bool snpp;
static int maxpaths_search;

/* For spliceable (really "joinable", if we consider indels) */
static Chrpos_T overall_max_distance;
static Chrpos_T overall_max_distance_novelend;

/* Other distances */
static int min_indel_end_matches;
static int max_middle_insertions_default; /* If negative, then compute querylength - 2*min_indel_end_matches */
static int max_middle_deletions;
static Chrpos_T shortsplicedist;
static Chrpos_T shortsplicedist_known;


/* Penalties */
static int subopt_levels;

static bool novelsplicingp;
static bool knownsplicingp;
static bool find_dna_chimeras_p;
static bool distances_observed_p;

static Chrpos_T min_intronlength;

/* Splicing */
static Univcoord_T *splicesites;
static Splicetype_T *splicetypes;
static Chrpos_T *splicedists;
static int nsplicesites;

static int min_distantsplicing_end_matches;
static int min_distantsplicing_identity;


/* GMAP parameters */
static bool gmap_segments_p;	/* previously called gmap_terminal_p.  Should move earlier (1). */
static bool gmap_pairsearch_p;	/* controls halfmapping.  Should move later (2). */
static bool gmap_improvement_p;	/* Should be at end (3). */
static bool gmap_indel_knownsplice_p;
static bool gmap_rerun_p = true;

static int antistranded_penalty;

static int nullgap;
static int maxpeelback;
static int maxpeelback_distalmedial;
static int extramaterial_end;
static int extramaterial_paired;
static int trigger_score_for_gmap;
static int gmap_allowance;
static int max_gmap_pairsearch;
static int max_gmap_segments;	/* Not used */
static int max_gmap_improvement;

static int max_floors_readlength;


#define A_CHAR 0x0
#define C_CHAR 0x1
#define G_CHAR 0x2
#define T_CHAR 0x3


/* Originally allowed only 1, to print only unique translocations.
   But need to allow enough to avoid missing some translocations. */
#define MAXCHIMERAPATHS 100

#define NREQUIRED_FAST 2	/* For candidate generation using
				   multimiss.  A value of 2 implies 
				   specificity of a 24-mer, which
				   should be low for a human-sized
				   genome */

#define MAX_INDEX1INTERVAL 3
#define STAGE2_MIN_OLIGO 3	/* Actually 6, but we are adding index1interval to this */
#define GOOD_GMAP_END 6
#define GMAP_TERMINAL_TRIM 6

#define GREEDY_SHORTSPLICEDIST 30000

static int index1part;
static int index1interval;
static int spansize;
static int max_anchors;

static int two_index1intervals;
static int min_kmer_readlength;
static Univ_IIT_T chromosome_iit;
static int circular_typeint;

static Univcoord_T *chroffsets;
static Univcoord_T *chrhighs;
static Chrpos_T *chrlengths; /* May differ from chrhigh - chroffset in circular chromosomes */
static int nchromosomes;
static Genome_T genome;

static int leftreadshift;
static Oligospace_T oligobase_mask; /* same as kmer_mask */
static int one_miss_querylength;

static int end_miss_one;	/* Used for computing max_terminal_length */
static int end_miss_two;	/* Used for computing max_terminal_length */


/* On 5' end, x = querypos.  On 3' end, x = (query_lastpos - querypos). */
#define FLOOR_END(x) ((x < index1interval) ? 0 : (x + spansize - index1interval)/spansize)

/* Here, x = (querypos - last_querypos).  Was (x-3)/12, but the new formula handles indels. */
#define FLOOR_MIDDLE(x) ((x < two_index1intervals) ? 0 : (x + spansize - two_index1intervals)/spansize)


#define MAX_LOCALSPLICING_POTENTIAL 1000
#if 0
/* Creates issues with ambiguous substrings */
#define LOCALSPLICING_NMATCHES_SLOP 1
#else
#define LOCALSPLICING_NMATCHES_SLOP 0
#endif
#define LOCALSPLICING_PROB_SLOP 0.05


/* Overall flow */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* identify_segments */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Indels */ 
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Indels, end */ 
#ifdef DEBUG2E
#define debug2e(x) x
#else
#define debug2e(x)
#endif

/* Floors */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* find_singlesplices */ 
#ifdef DEBUG4P
#define debug4p(x) x
#else
#define debug4p(x)
#endif

/* find_doublesplices */
#ifdef DEBUG4D
#define debug4d(x) x
#else
#define debug4d(x)
#endif

/* find_singlesplices */
#ifdef DEBUG4S
#define debug4s(x) x
#else
#define debug4s(x)
#endif

/* find_known_doublesplices */
#ifdef DEBUG4K
#define debug4k(x) x
#else
#define debug4k(x)
#endif


/* find_splicepairs_distant */
#ifdef DEBUG4L
#define debug4l(x) x
#else
#define debug4l(x)
#endif

/* find_splicepairs_distant (details) */
#ifdef DEBUG4LD
#define debug4ld(x) x
#else
#define debug4ld(x)
#endif

/* find_spliceends_shortend and find_spliceends_distant */
#ifdef DEBUG4E
#define debug4e(x) x
#else
#define debug4e(x)
#endif

/* Terminals */
#ifdef DEBUG4T
#define debug4t(x) x
#else
#define debug4t(x)
#endif

/* Short overlaps */ 
#ifdef DEBUG4H
#define debug4h(x) x
#else
#define debug4h(x)
#endif

/* Pairing up segments */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Heapify */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* Identify exact/onemiss/multimiss matches */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Identify onemiss matches, list contents */
#ifdef DEBUG7A
#define debug7a(x) x
#else
#define debug7a(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* straddling at beginning of genome.  May want to turn on DEBUG11 in indexdb_hr.c */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif

/* dual_search for known splice sites */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* GMAP */ 
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif

/* GMAP dump pairs */ 
#ifdef DEBUG13A
#define debug13a(x) x
#else
#define debug13a(x)
#endif

/* GMAP details of diagonals used */ 
#ifdef DEBUG13B
#define debug13b(x) x
#else
#define debug13b(x)
#endif

/* identify_all_segments */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

/* binary search method for updating chrnum */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif

/* consolidate_paired_results and choose_among_paired */ 
#ifdef DEBUG16
#define debug16(x) x
#else
#define debug16(x)
#endif


typedef struct Segment_T *Segment_T;
struct Segment_T {
  int splicesites_i;		/* if no splicesites_iit, then splicesites_i is -1 */
  Univcoord_T diagonal;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;
  Chrnum_T chrnum;

  int querypos5;
  int querypos3;
  Univcoord_T lowpos;		/* Needed for dynamic programming in converting segment to GMAP */
  Univcoord_T highpos;		/* Needed for dynamic programming in converting segment to GMAP */

  int floor;
  int floor_xfirst;
  int floor_xlast;

  int floor_left;
  int floor_right;

  int leftmost;			/* For segmenti of local splice */
  int rightmost;		/* For segmentj of local splice */
  bool left_splice_p; /* Set by find_singlesplices, used by find_doublesplices for speed */
  bool right_splice_p; /* Set by find_singlesplices, used by find_doublesplices for speed */
  bool spliceable_low_p;
  bool spliceable_high_p;

  bool usedp;
  bool pairablep;

#if 0
  int leftspan;			/* For segmentm of double splice */
  int rightspan;
#endif
};


static int
Segment_mod_length_cmp (const void *a, const void *b) {
  Segment_T x = * (Segment_T *) a;
  Segment_T y = * (Segment_T *) b;

  int xlength, ylength;
  int xmod, ymod;
  
  xmod = x->querypos5 % index1interval;
  ymod = y->querypos5 % index1interval;

  if (xmod < ymod) {
    return -1;
  } else if (ymod < xmod) {
    return +1;
  } else {
    xlength = x->querypos3 - x->querypos5;
    ylength = y->querypos3 - y->querypos5;
    if (xlength > ylength) {
      return -1;
    } else if (ylength > xlength) {
      return +1;
    } else {
      return 0;
    }
  }
}

static int
Segment_length_cmp (const void *a, const void *b) {
  Segment_T x = * (Segment_T *) a;
  Segment_T y = * (Segment_T *) b;

  int xlength, ylength;
  
  xlength = x->querypos3 - x->querypos5;
  ylength = y->querypos3 - y->querypos5;

  if (xlength > ylength) {
    return -1;
  } else if (ylength > xlength) {
    return +1;
  } else {
    return 0;
  }
}

static int
Segment_diagonal_cmp (const void *a, const void *b) {
  Segment_T x = * (Segment_T *) a;
  Segment_T y = * (Segment_T *) b;

  if (x->diagonal < y->diagonal) {
    return -1;
  } else if (y->diagonal < x->diagonal) {
    return +1;
  } else {
    return 0;
  }
}

static int
Segment_querypos5_ascending_cmp (const void *a, const void *b) {
  Segment_T x = * (Segment_T *) a;
  Segment_T y = * (Segment_T *) b;

  if (x->querypos5 < y->querypos5) {
    return -1;
  } else if (y->querypos5 < x->querypos5) {
    return +1;
  } else {
    return 0;
  }
}

static int
Segment_querypos3_ascending_cmp (const void *a, const void *b) {
  Segment_T x = * (Segment_T *) a;
  Segment_T y = * (Segment_T *) b;

  if (x->querypos3 < y->querypos3) {
    return -1;
  } else if (y->querypos3 < x->querypos3) {
    return +1;
  } else {
    return 0;
  }
}

static int
Segment_querypos5_descending_cmp (const void *a, const void *b) {
  Segment_T x = * (Segment_T *) a;
  Segment_T y = * (Segment_T *) b;

  if (x->querypos5 > y->querypos5) {
    return -1;
  } else if (y->querypos5 > x->querypos5) {
    return +1;
  } else {
    return 0;
  }
}

static int
Segment_querypos3_descending_cmp (const void *a, const void *b) {
  Segment_T x = * (Segment_T *) a;
  Segment_T y = * (Segment_T *) b;

  if (x->querypos3 > y->querypos3) {
    return -1;
  } else if (y->querypos3 > x->querypos3) {
    return +1;
  } else {
    return 0;
  }
}



struct Floors_T {
  int *allocated0;
  int *prev_omitted;

  int **allocated2;
  int *allocated1;
  int **scorefrom;		/* [from][to] */

  int **allocated4;
  int *allocated3;
  int **scoreto;		/* [to][from] */
};


void
Floors_free (Floors_T *old) {
  FREE((*old)->allocated1);
  FREE((*old)->allocated2);

  FREE((*old)->allocated3);
  FREE((*old)->allocated4);

  if ((*old)->allocated0) {
    FREE((*old)->allocated0);
  }

  FREE(*old);

  return;
}

void
Floors_free_keep (Floors_T *old) {

  FREE_KEEP((*old)->allocated1);
  FREE_KEEP((*old)->allocated2);

  FREE_KEEP((*old)->allocated3);
  FREE_KEEP((*old)->allocated4);

  if ((*old)->allocated0) {
    FREE_KEEP((*old)->allocated0);
  }

  FREE_KEEP(*old);

  return;
}

#ifdef DEBUG3
static void
Floors_print (Floors_T floors, int query_lastpos) {
  int from, to;

  if (floors->prev_omitted) {
    for (to = -index1interval; to <= query_lastpos+index1interval; to++) {
      printf("querypos %d, prev_omitted %d\n",to,floors->prev_omitted[to]);
    }
  }

  for (from = -index1interval; from <= query_lastpos+index1interval; from++) {
    for (to = from+1; to <= query_lastpos+index1interval; to++) {
      printf("from %d to %d, floor_score %d or %d",
	     from,to,floors->scorefrom[from][to],floors->scoreto[to][from]);
      if (floors->prev_omitted) {
	printf(" (prev %d)",floors->prev_omitted[to]);
      }
      printf("\n");
    }
  }

  return;
}
#endif


static Floors_T
Floors_new_standard (int querylength, int max_end_insertions, bool keep_floors_p) {
  Floors_T new;
  int query_lastpos, pos, from, to;
  int halfextra, extra;

  if (max_end_insertions < index1interval) {
    halfextra = index1interval;
  } else {
    halfextra = max_end_insertions;
  }
  extra = 1 + halfextra + halfextra;


  query_lastpos = querylength - index1part;

  if (keep_floors_p == true) {
    new = (Floors_T) MALLOC_KEEP(sizeof(*new));
  } else {
    new = (Floors_T) MALLOC(sizeof(*new));
  }
  new->allocated0 = (int *) NULL;
  new->prev_omitted = (int *) NULL;

  if (keep_floors_p == true) {
    new->allocated2 = (int **) CALLOC_KEEP(query_lastpos+extra,sizeof(int *));
    new->allocated1 = (int *) CALLOC_KEEP((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  } else {
    new->allocated2 = (int **) CALLOC(query_lastpos+extra,sizeof(int *));
    new->allocated1 = (int *) CALLOC((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  }
  new->allocated2[0] = &(new->allocated1[halfextra]);
  for (pos = 1; pos < query_lastpos+extra; pos++) {
    new->allocated2[pos] = &(new->allocated2[pos-1][query_lastpos+extra]);
  }
  new->scorefrom = &(new->allocated2[halfextra]);


  if (keep_floors_p == true) {
    new->allocated4 = (int **) CALLOC_KEEP(query_lastpos+extra,sizeof(int *));
    new->allocated3 = (int *) CALLOC_KEEP((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  } else {
    new->allocated4 = (int **) CALLOC(query_lastpos+extra,sizeof(int *));
    new->allocated3 = (int *) CALLOC((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  }
  new->allocated4[0] = &(new->allocated3[halfextra]);
  for (pos = 1; pos < query_lastpos+extra; pos++) {
    new->allocated4[pos] = &(new->allocated4[pos-1][query_lastpos+extra]);
  }
  new->scoreto = &(new->allocated4[halfextra]);

  for (to = -halfextra; to <= query_lastpos+halfextra; to++) {
    for (from = -halfextra; from < to; from++) {
      new->scorefrom[from][to] = new->scoreto[to][from] = FLOOR_MIDDLE(to - from);
    }
  }

  debug3(printf("Floors standard:\n"));
  debug3(Floors_print(new,query_lastpos));
  return new;
}


static Floors_T
Floors_new_omitted (int querylength, int max_end_insertions, bool *omitted) {
  Floors_T new;
  int query_lastpos, querypos, pos, from, to;
  int prev;
  int halfextra, extra;

  if (max_end_insertions < index1interval) {
    halfextra = index1interval;
  } else {
    halfextra = max_end_insertions;
  }
  extra = 1 + halfextra + halfextra;

  query_lastpos = querylength - index1part;
  new = (Floors_T) MALLOC(sizeof(*new));
  new->allocated0 = (int *) CALLOC(query_lastpos+extra,sizeof(int));
  new->prev_omitted = &(new->allocated0[halfextra]);

  new->allocated2 = (int **) CALLOC(query_lastpos+extra,sizeof(int *));
  new->allocated1 = (int *) CALLOC((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  new->allocated2[0] = &(new->allocated1[halfextra]);
  for (pos = 1; pos < query_lastpos+extra; pos++) {
    new->allocated2[pos] = &(new->allocated2[pos-1][query_lastpos+extra]);
  }
  new->scorefrom = &(new->allocated2[halfextra]);


  new->allocated4 = (int **) CALLOC(query_lastpos+extra,sizeof(int *));
  new->allocated3 = (int *) CALLOC((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  new->allocated4[0] = &(new->allocated3[halfextra]);
  for (pos = 1; pos < query_lastpos+extra; pos++) {
    new->allocated4[pos] = &(new->allocated4[pos-1][query_lastpos+extra]);
  }
  new->scoreto = &(new->allocated4[halfextra]);


  /* Set up omitted.  Save for middle_indels computation. */
  prev = -1;
  for (querypos = -halfextra; querypos < 0; querypos++) {
    new->prev_omitted[querypos] = -1;
  }
  for ( ; querypos <= query_lastpos; querypos++) {
    new->prev_omitted[querypos] = prev;
    if (omitted[querypos] == true) {
      prev = querypos;
    }
  }
  for ( ; querypos <= query_lastpos+halfextra; querypos++) {
    new->prev_omitted[querypos] = prev;
  }


  for (to = -halfextra; to <= query_lastpos+halfextra; to++) {
    prev = new->prev_omitted[to];
    for (from = -halfextra; from < prev; from++) {
      new->scorefrom[from][to] = new->scorefrom[from][prev] + FLOOR_MIDDLE(to - prev);
    }
    for ( ; from < to; from++) {
      new->scorefrom[from][to] = FLOOR_MIDDLE(to - from);
    }
  }

  for (to = -halfextra; to <= query_lastpos+halfextra; to++) {
    for (from = -halfextra; from < to; from++) {
      new->scoreto[to][from] = new->scorefrom[from][to];
    }
  }

  debug3(
	 printf("Floors omitted:");
	 for (pos = 0; pos <= query_lastpos; pos++) {
	   if (omitted[pos] == true) {
	     printf(" %d",pos);
	   }
	 }
	 printf("\n");
	 )
  debug3(Floors_print(new,query_lastpos));
  return new;
}


/************************************************************************/


#define T Stage1_T
struct T {
  Spanningelt_T *plus_spanningset[MAX_INDEX1INTERVAL];
  Spanningelt_T *minus_spanningset[MAX_INDEX1INTERVAL];

  struct Spanningelt_T *plus_spanningset_allocated[MAX_INDEX1INTERVAL];
  struct Spanningelt_T *minus_spanningset_allocated[MAX_INDEX1INTERVAL];

  int plus_spanningset_nelts[MAX_INDEX1INTERVAL];
  int minus_spanningset_nelts[MAX_INDEX1INTERVAL];

  bool read_oligos_p;

#ifdef LARGE_GENOMES
  unsigned char **positions_high_allocated;
  unsigned char **plus_positions_high; /* points to above[index1interval-1] */
  unsigned char **minus_positions_high; /* points to above[index1interval-1] */

  UINT4 **positions_low_allocated;
  UINT4 **plus_positions_low; /* points to above[index1interval-1] */
  UINT4 **minus_positions_low; /* points to above[index1interval-1] */
#else
  Univcoord_T **positions_allocated;
  Univcoord_T **plus_positions; /* points to above[index1interval-1] */
  Univcoord_T **minus_positions; /* points to above[index1interval-1] */
#endif

  int *npositions_allocated;
  int *plus_npositions;		/* points to above[index1interval-1] */
  int *minus_npositions;	/* points to above[index1interval-1] */

  bool *retrievedp_allocated;
  bool *plus_retrievedp;	/* points to above[index1interval-1] */
  bool *minus_retrievedp;	/* points to above[index1interval-1] */

#ifdef USE_ALLOCP
  bool *allocp_allocated;
  bool *plus_allocp;		/* points to above[index1interval-1] */
  bool *minus_allocp;		/* points to above[index1interval-1] */
#endif

#ifdef USE_VALIDP
  bool *validp;
#endif
  bool *omitted;

  Oligospace_T *oligos_allocated;
  Oligospace_T *forward_oligos; /* points to above[index1interval-1] */
  Oligospace_T *revcomp_oligos; /* points to above[index1interval-1] */

  struct Segment_T *plus_segments;
  struct Segment_T *minus_segments;
  int plus_nsegments;
  int minus_nsegments;

  Segment_T *plus_spliceable; /* plus_segments with a following diagonal within shortsplicedist or splicedists[j] */
  Segment_T *minus_spliceable; /* minus_segments with a following diagonal within shortsplicedist or splicedists[j] */
  int plus_nspliceable;
  int minus_nspliceable;

  bool all_positions_fetched_p;
};


static void
stage3list_gc (List_T *old) {
  List_T p;
  Stage3end_T hit;

  for (p = *old; p != NULL; p = p->rest) {
    hit = (Stage3end_T) p->first;
    Stage3end_free(&hit);
  }
  List_free(&(*old));
  return;
}

static void
substringlist_gc (List_T *old) {
  List_T p;
  Substring_T hit;

  for (p = *old; p != NULL; p = p->rest) {
    hit = (Substring_T) p->first;
    Substring_free(&hit);
  }
  List_free(&(*old));
  return;
}


static bool free_positions_p;

void
Stage1_init_positions_free (bool positions_fileio_p) {
  if (positions_fileio_p == true) {
    free_positions_p = true;
  } else {
    free_positions_p = false;
  }
  return;
}


void
Stage1_free (T *old, int querylength) {
  Spanningelt_T spanningelt;
  int mod, i;

  /* Stage1hr_check(*old); */

  if (*old) {
    FREE((*old)->plus_spliceable);
    FREE((*old)->minus_spliceable);

    FREE((*old)->plus_segments);
    FREE((*old)->minus_segments);

    for (mod = 0; mod < index1interval; mod++) {
      for (i = 0; i < (*old)->plus_spanningset_nelts[mod]; i++) {
	spanningelt = (Spanningelt_T) (*old)->plus_spanningset[mod][i];
	Spanningelt_gc(spanningelt);
      }
      FREE((*old)->plus_spanningset_allocated[mod]);
      FREE((*old)->plus_spanningset[mod]);

      for (i = 0; i < (*old)->minus_spanningset_nelts[mod]; i++) {
	spanningelt = (Spanningelt_T) (*old)->minus_spanningset[mod][i];
	Spanningelt_gc(spanningelt);
      }
      FREE((*old)->minus_spanningset_allocated[mod]);
      FREE((*old)->minus_spanningset[mod]);
    }

    if (free_positions_p == true) {
      for (i = -index1interval+1; i < querylength; i++) {
	if ((*old)->plus_retrievedp[i] == true) {
#ifdef LARGE_GENOMES
	  FREE((*old)->plus_positions_high[i]);
	  FREE((*old)->plus_positions_low[i]);
#else
	  FREE((*old)->plus_positions[i]);
#endif
	}
	if ((*old)->minus_retrievedp[i] == true) {
#ifdef LARGE_GENOMES
	  FREE((*old)->minus_positions_high[i]);
	  FREE((*old)->minus_positions_low[i]);
#else
	  FREE((*old)->minus_positions[i]);
#endif
	}
      }
#ifdef USE_ALLOCP
    } else {
      for (i = -index1interval+1; i < querylength; i++) {
	if ((*old)->plus_allocp[i] == true) {
#ifdef LARGE_GENOMES
	  FREE((*old)->plus_positions_high[i]);
	  FREE((*old)->plus_positions_low[i]);
#else
	  FREE((*old)->plus_positions[i]);
#endif
	}
	if ((*old)->minus_allocp[i] == true) {
#ifdef LARGE_GENOMES
	  FREE((*old)->minus_positions_high[i]);
	  FREE((*old)->minus_positions_low[i]);
#else
	  FREE((*old)->minus_positions[i]);
#endif
	}
      }
#endif
    }

    FREE((*old)->oligos_allocated);
    FREE((*old)->omitted);
#ifdef USE_VALIDP
    FREE((*old)->validp);
#endif
#ifdef LARGE_GENOMES
    FREE((*old)->positions_high_allocated);
    FREE((*old)->positions_low_allocated);
#else
    FREE((*old)->positions_allocated);
#endif
    FREE((*old)->npositions_allocated);
#ifdef USE_ALLOCP
    FREE((*old)->allocp_allocated);
#endif
    FREE((*old)->retrievedp_allocated);

    FREE(*old);
  }

  return;
}


/************************************************************************/

static bool
check_dinucleotides (char *sequence, int querylength) {
  int index, firsti;
  int n = 0, i;
  int c1, c2;
  bool validp;
  int dinucl_counts[16], first_count, second_count;

  for (index = 0; index < 16; index++) {
    dinucl_counts[index] = 0;
  }

  for (i = 0; i < querylength - 1; i++) {
    c1 = sequence[i];
    c2 = sequence[i+1];
    validp = true;

    switch (c1) {
    case 'A': index = 0; break;
    case 'C': index = 4; break;
    case 'G': index = 8; break;
    case 'T': index = 12; break;
    default: validp = false;
    }

    switch (c2) {
    case 'A': break;
    case 'C': index += 1; break;
    case 'G': index += 2; break;
    case 'T': index += 3; break;
    default: validp = false;
    }
      
    if (validp == true) {
      dinucl_counts[index] += 1;
      n++;
    }
  }

  if (n == 0) {
    return false;
  } else {
    first_count = 0;
    for (index = 0; index < 16; index++) {
      debug(printf("%d: %d\n",index,dinucl_counts[index]));
      if (dinucl_counts[index] > first_count) {
	first_count = dinucl_counts[index];
	firsti = index;
      }
    }

    second_count = 0;
    for (index = 0; index < 16; index++) {
      if (index == firsti) {
	/* Skip */
      } else if (dinucl_counts[index] > second_count) {
	second_count = dinucl_counts[index];
      }
    }

    debug(printf("first count: %d, second count: %d",first_count,second_count));
    if (first_count + second_count > 0.80 * (double) n) {
      debug(printf(" > 0.80*%d = %.2f\n",n,0.80*n));
      return false;
    } else {
      debug(printf("\n"));
      return true;
    }
  }
}



static int
read_oligos (bool *allvalidp, T this, char *queryuc_ptr, int querylength,
	     int query_lastpos, int genestrand) {
  Reader_T reader;
  int querypos, noligos = 0;
  Oligostate_T last_state = INIT;
  Oligospace_T forward = 0, revcomp = 0;

  /* This estimate may be too high */
  /* this->maxfloor = 1 + querylength/oligobase * 2; */

  if (use_only_sarray_p == true) {
    *allvalidp = false;
    return 1;
  } else if (use_sarray_p == true && querylength < min_kmer_readlength) {
    *allvalidp = false;
    return 1;
  } else {
    reader = Reader_new(queryuc_ptr,/*querystart*/0,/*queryend*/querylength);
  }

  /* Prevents us from processing invalid query 12-mers */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    this->plus_retrievedp[querypos] = true;
#ifdef LARGE_GENOMES
    this->plus_positions_high[querypos] = (unsigned char *) NULL;
    this->plus_positions_low[querypos] = (UINT4 *) NULL;
#else
    this->plus_positions[querypos] = (Univcoord_T *) NULL;
#endif
    this->plus_npositions[querypos] = 0;

    this->minus_retrievedp[querypos] = true;
#ifdef LARGE_GENOMES
    this->minus_positions_high[querypos] = (unsigned char *) NULL;
    this->minus_positions_low[querypos] = (UINT4 *) NULL;
#else
    this->minus_positions[querypos] = (Univcoord_T *) NULL;
#endif
    this->minus_npositions[querypos] = 0;
  }

  /* Note: leftshifting is done here, rather than in Oligo_lookup */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  debug(printf("oligobase_mask: %010llX\n",oligobase_mask));
#if 0
  *any_omitted_p = false;
  *all_omitted_p = true;
#endif
  if (mode == STANDARD) {
    while ((last_state = Oligo_next(last_state,&querypos,&forward,&revcomp,
				    reader,/*cdnaend*/FIVE)) != DONE) {
#ifdef LARGE_GENOMES
      this->plus_positions_high[querypos] = (unsigned char *) NULL;
      this->plus_positions_low[querypos] = (UINT4 *) NULL;
      this->minus_positions_high[querypos] = (unsigned char *) NULL;
      this->minus_positions_low[querypos] = (UINT4 *) NULL;
#else
      this->plus_positions[querypos] = (Univcoord_T *) NULL;
      this->minus_positions[querypos] = (Univcoord_T *) NULL;
#endif
      this->plus_npositions[querypos] = 0;
      this->minus_npositions[querypos] = 0;

      if (last_state == VALID) {
#ifdef USE_VALIDP
	this->validp[querypos] = true;
#endif
	this->plus_retrievedp[querypos] = false;
	this->minus_retrievedp[querypos] = false;

	this->forward_oligos[querypos] = forward & oligobase_mask;
	this->revcomp_oligos[querypos] = (revcomp >> leftreadshift) & oligobase_mask;

	debug(printf("At querypos %d, read oligo = %010llX\n",querypos,this->forward_oligos[querypos]));
	noligos++;
      }
    }

  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    if (genestrand == +2) {
      while ((last_state = Oligo_next(last_state,&querypos,&forward,&revcomp,
				      reader,/*cdnaend*/FIVE)) != DONE) {
#ifdef LARGE_GENOMES
	this->plus_positions_high[querypos] = (unsigned char *) NULL;
	this->plus_positions_low[querypos] = (UINT4 *) NULL;
	this->minus_positions_high[querypos] = (unsigned char *) NULL;
	this->minus_positions_low[querypos] = (UINT4 *) NULL;
#else
	this->plus_positions[querypos] = (Univcoord_T *) NULL;
	this->minus_positions[querypos] = (Univcoord_T *) NULL;
#endif
	this->plus_npositions[querypos] = 0;
	this->minus_npositions[querypos] = 0;

	if (last_state == VALID) {
#ifdef USE_VALIDP
	  this->validp[querypos] = true;
#endif
	  this->plus_retrievedp[querypos] = false;
	  this->minus_retrievedp[querypos] = false;

	  this->forward_oligos[querypos] = Cmet_reduce_ga(forward) & oligobase_mask;
	  this->revcomp_oligos[querypos] = Cmet_reduce_ct(revcomp >> leftreadshift) & oligobase_mask;

	  debug(printf("At querypos %d, read oligo = %010llX\n",querypos,this->forward_oligos[querypos]));
	  noligos++;
	}
      }

    } else {
      while ((last_state = Oligo_next(last_state,&querypos,&forward,&revcomp,
				      reader,/*cdnaend*/FIVE)) != DONE) {
#ifdef LARGE_GENOMES
	this->plus_positions_high[querypos] = (unsigned char *) NULL;
	this->plus_positions_low[querypos] = (UINT4 *) NULL;
	this->minus_positions_high[querypos] = (unsigned char *) NULL;
	this->minus_positions_low[querypos] = (UINT4 *) NULL;
#else
	this->plus_positions[querypos] = (Univcoord_T *) NULL;
	this->minus_positions[querypos] = (Univcoord_T *) NULL;
#endif
	this->plus_npositions[querypos] = 0;
	this->minus_npositions[querypos] = 0;

	if (last_state == VALID) {
#ifdef USE_VALIDP
	  this->validp[querypos] = true;
#endif
	  this->plus_retrievedp[querypos] = false;
	  this->minus_retrievedp[querypos] = false;

	  this->forward_oligos[querypos] = Cmet_reduce_ct(forward) & oligobase_mask;
	  this->revcomp_oligos[querypos] = Cmet_reduce_ga(revcomp >> leftreadshift) & oligobase_mask;

	  debug(printf("At querypos %d, read oligo = %010llX\n",querypos,this->forward_oligos[querypos]));
	  noligos++;
	}
      }
    }

  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    if (genestrand == +2) {
      while ((last_state = Oligo_next(last_state,&querypos,&forward,&revcomp,
				      reader,/*cdnaend*/FIVE)) != DONE) {
#ifdef LARGE_GENOMES
	this->plus_positions_high[querypos] = (unsigned char *) NULL;
	this->plus_positions_low[querypos] = (UINT4 *) NULL;
	this->minus_positions_high[querypos] = (unsigned char *) NULL;
	this->minus_positions_low[querypos] = (UINT4 *) NULL;
#else
	this->plus_positions[querypos] = (Univcoord_T *) NULL;
	this->minus_positions[querypos] = (Univcoord_T *) NULL;
#endif
	this->plus_npositions[querypos] = 0;
	this->minus_npositions[querypos] = 0;

	if (last_state == VALID) {
#ifdef USE_VALIDP
	  this->validp[querypos] = true;
#endif
	  this->plus_retrievedp[querypos] = false;
	  this->minus_retrievedp[querypos] = false;
	  
	  this->forward_oligos[querypos] = Atoi_reduce_tc(forward) & oligobase_mask;
	  this->revcomp_oligos[querypos] = Atoi_reduce_ag(revcomp >> leftreadshift) & oligobase_mask;

	  debug(printf("At querypos %d, read oligo = %010llX\n",querypos,this->forward_oligos[querypos]));
	  noligos++;
	}
      }
    } else {
      while ((last_state = Oligo_next(last_state,&querypos,&forward,&revcomp,
				      reader,/*cdnaend*/FIVE)) != DONE) {
#ifdef LARGE_GENOMES
	this->plus_positions_high[querypos] = (unsigned char *) NULL;
	this->plus_positions_low[querypos] = (UINT4 *) NULL;
	this->minus_positions_high[querypos] = (unsigned char *) NULL;
	this->minus_positions_low[querypos] = (UINT4 *) NULL;
#else
	this->plus_positions[querypos] = (Univcoord_T *) NULL;
	this->minus_positions[querypos] = (Univcoord_T *) NULL;
#endif
	this->plus_npositions[querypos] = 0;
	this->minus_npositions[querypos] = 0;

	if (last_state == VALID) {
#ifdef USE_VALIDP
	  this->validp[querypos] = true;
#endif
	  this->plus_retrievedp[querypos] = false;
	  this->minus_retrievedp[querypos] = false;
	  
	  this->forward_oligos[querypos] = Atoi_reduce_ag(forward) & oligobase_mask;
	  this->revcomp_oligos[querypos] = Atoi_reduce_tc(revcomp >> leftreadshift) & oligobase_mask;

	  debug(printf("At querypos %d, read oligo = %010llX\n",querypos,this->forward_oligos[querypos]));
	  noligos++;
	}
      }
    }

  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    if (genestrand == +2) {
      while ((last_state = Oligo_next(last_state,&querypos,&forward,&revcomp,
				      reader,/*cdnaend*/FIVE)) != DONE) {
#ifdef LARGE_GENOMES
	this->plus_positions_high[querypos] = (unsigned char *) NULL;
	this->plus_positions_low[querypos] = (UINT4 *) NULL;
	this->minus_positions_high[querypos] = (unsigned char *) NULL;
	this->minus_positions_low[querypos] = (UINT4 *) NULL;
#else
	this->plus_positions[querypos] = (Univcoord_T *) NULL;
	this->minus_positions[querypos] = (Univcoord_T *) NULL;
#endif
	this->plus_npositions[querypos] = 0;
	this->minus_npositions[querypos] = 0;

	if (last_state == VALID) {
#ifdef USE_VALIDP
	  this->validp[querypos] = true;
#endif
	  this->plus_retrievedp[querypos] = false;
	  this->minus_retrievedp[querypos] = false;
	  
	  this->forward_oligos[querypos] = Atoi_reduce_ag(forward) & oligobase_mask;
	  this->revcomp_oligos[querypos] = Atoi_reduce_tc(revcomp >> leftreadshift) & oligobase_mask;

	  debug(printf("At querypos %d, read oligo = %010llX\n",querypos,this->forward_oligos[querypos]));
	  noligos++;
	}
      }
    } else {
      while ((last_state = Oligo_next(last_state,&querypos,&forward,&revcomp,
				      reader,/*cdnaend*/FIVE)) != DONE) {
#ifdef LARGE_GENOMES
	this->plus_positions_high[querypos] = (unsigned char *) NULL;
	this->plus_positions_low[querypos] = (UINT4 *) NULL;
	this->minus_positions_high[querypos] = (unsigned char *) NULL;
	this->minus_positions_low[querypos] = (UINT4 *) NULL;
#else
	this->plus_positions[querypos] = (Univcoord_T *) NULL;
	this->minus_positions[querypos] = (Univcoord_T *) NULL;
#endif
	this->plus_npositions[querypos] = 0;
	this->minus_npositions[querypos] = 0;

	if (last_state == VALID) {
#ifdef USE_VALIDP
	  this->validp[querypos] = true;
#endif
	  this->plus_retrievedp[querypos] = false;
	  this->minus_retrievedp[querypos] = false;
	  
	  this->forward_oligos[querypos] = Atoi_reduce_tc(forward) & oligobase_mask;
	  this->revcomp_oligos[querypos] = Atoi_reduce_ag(revcomp >> leftreadshift) & oligobase_mask;

	  debug(printf("At querypos %d, read oligo = %010llX\n",querypos,this->forward_oligos[querypos]));
	  noligos++;
	}
      }
    }
  }

  if (noligos < query_lastpos + 1) {
    debug(printf("Read only %d oligos due to non-ACGT; expected %d\n",noligos,query_lastpos + 1));
    *allvalidp = false;
  } else {
    *allvalidp = true;
  }

  Reader_free(&reader);

  this->read_oligos_p = true;
  return noligos;
}




/************************************************************************
 *   Omitted:
 *   In all cases, want to omit poly-AT.
 *   For purposes of finding mismatches, may want to omit frequent oligomers also
 *   For purposes of finding indels, may want to omit repetitive oligomers at ends also
 ************************************************************************/

static void
omit_oligos_clear (T this, int query_lastpos) {
  int querypos;

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    this->omitted[querypos] = false;
  }
  return;
}


#if 0
static void
omit_oligos_polyat (bool *all_omitted_p, bool *any_omitted_p, T this, int query_lastpos) {
  int querypos;

  *all_omitted_p = true;
  *any_omitted_p = false;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->forward_oligos[querypos] == 0U || this->revcomp_oligos[querypos] == 0U) {
      this->omitted[querypos] = true;
      *any_omitted_p = true;
    } else {
      this->omitted[querypos] = false;
      *all_omitted_p = false;
    }
  }

  return;
}
#endif


static void
omit_oligos (bool *all_omitted_p, bool *any_omitted_p, T this, int query_lastpos,
	     int indexdb_size_threshold, bool frequentp, bool repetitivep) {
  int querypos;
  bool still_repetitive_p;

  *any_omitted_p = false;

  /* Always omit poly-AT */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->forward_oligos[querypos] == 0U || this->revcomp_oligos[querypos] == 0U) {
      debug(printf("Querypos %d is poly-A or poly-T\n",querypos));
      this->omitted[querypos] = true;
      *any_omitted_p = true;
    } else {
      this->omitted[querypos] = false;
    }
  }

  if (frequentp == true) {
    /* Omit frequent oligos, but only in the middle */
    for (querypos = index1interval; querypos <= query_lastpos-index1interval; querypos++) {
      if (this->plus_npositions[querypos] > indexdb_size_threshold &&
	  this->minus_npositions[querypos] > indexdb_size_threshold) {
	debug(printf("Querypos %d is frequent with %d plus positions > %d and %d minus positions > %d\n",
		     querypos,this->plus_npositions[querypos],indexdb_size_threshold,
		     this->minus_npositions[querypos],indexdb_size_threshold));
	this->omitted[querypos] = true;
	*any_omitted_p = true;
      }
    }

#if 0
    /* This avoids too many consecutive omitted, but slows down finding double splicing significantly */
    if (*any_omitted_p == true) {
      nconsecutive = 0;
      for (querypos = index1interval; querypos <= query_lastpos-index1interval; querypos += index1interval) {
	if (this->omitted[querypos] == false) {
	  nconsecutive = 0;
	} else {
	  nconsecutive++;
	  if (nconsecutive == 4) { /* corresponds to 4*3 = 12 positions */
	    debug(printf("Consecutive frequent from %d to %d.  ",querypos-11,querypos));
	    this->omitted[querypos] = false;
	    nconsecutive = 0;
	  }
	}
      }

      nconsecutive = 0;
      for (querypos = index1interval+1; querypos <= query_lastpos-index1interval; querypos += index1interval) {
	if (this->omitted[querypos] == false) {
	  nconsecutive = 0;
	} else {
	  nconsecutive++;
	  if (nconsecutive == 4) { /* corresponds to 4*3 = 12 positions */
	    debug(printf("Consecutive frequent from %d to %d.  ",querypos-11,querypos));
	    this->omitted[querypos] = false;
	    nconsecutive = 0;
	  }
	}
      }

      nconsecutive = 0;
      for (querypos = index1interval+index1interval-1; querypos <= query_lastpos-index1interval; querypos += index1interval) {
	if (this->omitted[querypos] == false) {
	  nconsecutive = 0;
	} else {
	  nconsecutive++;
	  if (nconsecutive == 4) { /* corresponds to 4*3 = 12 positions */
	    debug(printf("Consecutive frequent from %d to %d.  ",querypos-11,querypos));
	    this->omitted[querypos] = false;
	    nconsecutive = 0;
	  }
	}
      }
    }
#endif
  }

  if (repetitivep == true) {
    /* Omit repetitive oligos at the ends */
    still_repetitive_p = true;
    querypos = 0;
    while (querypos <= query_lastpos && still_repetitive_p == true) {
      if (Oligo_repetitive_p(this->forward_oligos[querypos])) {
	debug(printf("Querypos %d is repetitive\n",querypos));
	this->omitted[querypos] = true;
	*any_omitted_p = true;
      } else {
	still_repetitive_p = false;
      }
      querypos++;
    }

    still_repetitive_p = true;
    querypos = query_lastpos;
    while (querypos >= 0 && still_repetitive_p == true) {
      if (Oligo_repetitive_p(this->forward_oligos[querypos])) {
	debug(printf("Querypos %d is repetitive\n",querypos));
	this->omitted[querypos] = true;
	*any_omitted_p = true;
      } else {
	still_repetitive_p = false;
      }
      querypos--;
    }
  }

  if (*any_omitted_p == false) {
    debug(printf("No oligos are omitted\n"));
    *all_omitted_p = false;
  } else {
    debug(
	  printf("Omitted oligos:");
	  for (querypos = 0; querypos <= query_lastpos; querypos++) {
	    if (this->omitted[querypos] == true) {
	      printf(" %d",querypos);
	    }
	  }
	  printf("\n"));

    *all_omitted_p = true;
    for (querypos = 0; querypos <= query_lastpos; querypos++) {
      if (this->omitted[querypos] == false) {
	*all_omitted_p = false;
      }
    }
  }

  return;
}


#if 0
static void
omit_oligos_repetitive (bool *all_omitted_p, bool *any_omitted_p, T this, int query_lastpos) {
  int querypos;

  *all_omitted_p = true;
  *any_omitted_p = false;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (Oligo_repetitive_p(this->forward_oligos[querypos])) {
      this->omitted[querypos] = true;
      *any_omitted_p = true;
    } else {
      this->omitted[querypos] = false;
      *all_omitted_p = false;
    }
  }

  return;
}
#endif


static T
Stage1_new (int querylength) {
  T new = (T) MALLOC(sizeof(*new));
  int mod;
  int overhang = index1interval-1;

  for (mod = 0; mod < index1interval; mod++) {
    new->plus_spanningset[mod] = (Spanningelt_T *) NULL;
    new->minus_spanningset[mod] = (Spanningelt_T *) NULL;
    new->plus_spanningset_allocated[mod] = (struct Spanningelt_T *) NULL;
    new->minus_spanningset_allocated[mod] = (struct Spanningelt_T *) NULL;
    new->plus_spanningset_nelts[mod] = 0;
    new->minus_spanningset_nelts[mod] = 0;
  }

  new->read_oligos_p = false;

#ifdef LARGE_GENOMES
  new->positions_high_allocated = (unsigned char **) MALLOC(2 * (querylength+overhang) * sizeof(unsigned char *));
  new->plus_positions_high = &(new->positions_high_allocated[overhang]);
  new->minus_positions_high = &(new->positions_high_allocated[(querylength+overhang)+overhang]);

  new->positions_low_allocated = (UINT4 **) MALLOC(2 * (querylength+overhang) * sizeof(UINT4 *));
  new->plus_positions_low = &(new->positions_low_allocated[overhang]);
  new->minus_positions_low = &(new->positions_low_allocated[(querylength+overhang)+overhang]);
#else
  new->positions_allocated = (Univcoord_T **) MALLOC(2 * (querylength+overhang) * sizeof(Univcoord_T *));
  new->plus_positions = &(new->positions_allocated[overhang]);
  new->minus_positions = &(new->positions_allocated[(querylength+overhang)+overhang]);
#endif

  new->npositions_allocated = (int *) MALLOC(2 * (querylength+overhang) * sizeof(int));
  new->plus_npositions = &(new->npositions_allocated[overhang]);
  new->minus_npositions = &(new->npositions_allocated[(querylength+overhang)+overhang]);

#if 0
  /* No need to initialize, since we assign all values below */
  for (querypos = -index1interval+1; querypos < querylength; querypos++) {
    new->plus_positions[querypos] = (Univcoord_T *) NULL;
    new->plus_npositions[querypos] = 0;
    new->minus_positions[querypos] = (Univcoord_T *) NULL;
    new->minus_npositions[querypos] = 0;
  }
#endif

  /* Can be MALLOC, since we initialize in read_oligos() */
  new->retrievedp_allocated = (bool *) MALLOC(2 * (querylength+overhang) * sizeof(bool));
  new->plus_retrievedp = &(new->retrievedp_allocated[overhang]);
  new->minus_retrievedp = &(new->retrievedp_allocated[(querylength+overhang)+overhang]);

#ifdef USE_ALLOCP
  /* Never set to true, so never used */
  new->allocp_allocated = (bool *) CALLOC(2 * (querylength+overhang),sizeof(bool));
  new->plus_allocp = &(new->allocp_allocated[overhang]);
  new->minus_allocp = &(new->allocp_allocated[(querylength+overhang)+overhang]);
#endif

#ifdef USE_VALIDP
  new->validp = (bool *) CALLOC(querylength,sizeof(bool));
#endif
  new->omitted = (bool *) CALLOC(querylength,sizeof(bool));

  new->oligos_allocated = (Oligospace_T *) CALLOC(2 * (querylength+overhang),sizeof(Oligospace_T));
  new->forward_oligos = &(new->oligos_allocated[overhang]);
  new->revcomp_oligos = &(new->oligos_allocated[(querylength+overhang)+overhang]);

  new->plus_segments = (struct Segment_T *) NULL;
  new->minus_segments = (struct Segment_T *) NULL;
  new->plus_nsegments = 0;
  new->minus_nsegments = 0;

  new->plus_spliceable = (Segment_T *) NULL;
  new->minus_spliceable = (Segment_T *) NULL;
  new->plus_nspliceable = 0;
  new->minus_nspliceable = 0;

  new->all_positions_fetched_p = false;

  return new;
}


/************************************************************************/

static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}


#if 0
static void
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

  return;
}
#endif

/************************************************************************/

#define PARENT(i) ((i) >> 1)
#define LEFT(i) ((i) << 1)
#define RIGHT(i) (((i) << 1) | 1)


typedef struct Batch_T *Batch_T;

struct Batch_T {
#ifndef USE_HEAP
  int nodei;			/* Node in loser tree.  Also used for debugging */
#endif
  int querypos;
  int diagterm;
  int npositions;

  Univcoord_T diagonal;
#ifdef LARGE_GENOMES
  unsigned char *positions_high;
  UINT4 *positions_low;
#endif
  Univcoord_T *positions;
#ifdef DIAGONAL_ADD_QUERYPOS
  UINT8 diagonal_add_querypos;
#endif
};


static void
Batch_init (Batch_T batch, int querypos, int diagterm,
#ifdef LARGE_GENOMES
	    unsigned char *positions_high, UINT4 *positions_low,
#else
	    Univcoord_T *positions,
#endif
	    int npositions, int querylength
#ifndef USE_HEAP
	    , int nodei
#endif
	    ) {

#ifndef USE_HEAP
  batch->nodei = nodei;
#endif
  batch->querypos = querypos;
  batch->diagterm = diagterm;
  batch->npositions = npositions;

#ifdef LARGE_GENOMES
  batch->positions_high = positions_high;
  batch->positions_low = positions_low;
  batch->diagonal = (((Univcoord_T) *positions_high) << 32) + (*positions_low) + diagterm;
#elif defined(WORDS_BIGENDIAN)
  batch->positions = positions;
  batch->diagonal = Bigendian_convert_univcoord(*positions) + diagterm;
#else
  batch->positions = positions;
  batch->diagonal = *positions + diagterm;
#endif


#ifdef NO_EXTENSIONS_BEFORE_ZERO
  /* This prevents us from finding insertions at the beginning of the genome */
  while (batch->npositions > 0 && batch->diagonal < (unsigned int) querylength) {
    debug11(printf("Eliminating diagonal %llu as straddling beginning of genome (Batch_init)\n",
		   (unsigned long long) batch->diagonal));
    batch->npositions--;
    if (batch->npositions > 0) {
#ifdef LARGE_GENOMES
      batch->diagonal = ((Univcoord_T) *(++batch->positions_high) << 32) + *(++batch->positions_low) + diagterm;
#elif defined(WORDS_BIGENDIAN)
      batch->diagonal = Bigendian_convert_univcoord(*(++batch->positions)) + diagterm;
#else
      batch->diagonal = *(++batch->positions) + diagterm;
#endif
    }
  }
#endif


#ifdef DIAGONAL_ADD_QUERYPOS
  batch->diagonal_add_querypos = (UINT8) batch->diagonal;
  batch->diagonal_add_querypos <<= 32;
  batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif

  return;
}


static void
Batch_init_simple (Batch_T batch, Univcoord_T *diagonals, int ndiagonals, int querylength,
#ifdef USE_HEAP
		   int querypos
#else
		   int nodei
#endif
		   ) {

#ifdef USE_HEAP
  batch->querypos = querypos;
#else
  batch->nodei = nodei;
#endif
  batch->positions = diagonals;
  batch->npositions = ndiagonals;
#ifdef USE_HEAP
  batch->diagonal = *diagonals;	/* Already in correct endianness */
#else
  if (batch->npositions == 0) {
    batch->diagonal = -1U;
  } else {
    batch->diagonal = *diagonals;	/* Already in correct endianness */
  }
#endif

  while (batch->npositions > 0 && batch->diagonal < (unsigned int) querylength) {
    debug11(printf("Eliminating diagonal %llu as straddling beginning of genome (Batch_init)\n",
		   (unsigned long long) batch->diagonal));
    batch->npositions--;
    if (batch->npositions > 0) {
      /* positions are really diagonals, already in correct endianness */
      batch->diagonal = *(++batch->positions);
    }
  }

  return;
}


static void
min_heap_insert (Batch_T *heap, int *heapsize, Batch_T batch) {
  int i;
#ifdef DIAGONAL_ADD_QUERYPOS
  UINT8 diagonal_add_querypos;
#else
  int querypos;
  Univcoord_T diagonal;
#endif

  i = ++(*heapsize);
#ifdef DIAGONAL_ADD_QUERYPOS
  diagonal_add_querypos = batch->diagonal_add_querypos;
  while (i > 1 && (heap[PARENT(i)]->diagonal_add_querypos > diagonal_add_querypos)) {
    heap[i] = heap[PARENT(i)];
    i = PARENT(i);
  }
#else
  querypos = batch->querypos;
  diagonal = batch->diagonal;
  /* sort primarily by diagonal, then by querypos */
  while (i > 1 && (heap[PARENT(i)]->diagonal > diagonal ||
		   (heap[PARENT(i)]->diagonal == diagonal && heap[PARENT(i)]->querypos > querypos))) {
    heap[i] = heap[PARENT(i)];
    i = PARENT(i);
  }
#endif
  heap[i] = batch;

  return;
}


static void
min_heap_insert_simple (Batch_T *heap, int *heapsize, Batch_T batch) {
  int i;
  Univcoord_T diagonal;

  i = ++(*heapsize);
  diagonal = batch->diagonal;
  while (i > 1 && (heap[PARENT(i)]->diagonal > diagonal)) {
    heap[i] = heap[PARENT(i)];
    i = PARENT(i);
  }
  heap[i] = batch;

  return;
}



/* Note FORMULA: formulas for querypos <-> diagonal (diagterm in call to Indexdb_read) are:

plus: diagonal = position + querylength - querypos
minus: diagonal = position + querypos + index1part

For minus, the index1part is needed in call to Indexdb_read because
position is stored at beginning of plus oligomer, which corresponds to
end of minus oligomer.  As a result, we have the following formulas:

high genomic position = diagonal (corresponds to querypos =
querylength for plus, and querypos = 0 for minus)

low genomic position = diagonal - querylength (corresponds to querypos
= 0 for plus, and querypos = querylength for minus)

*/


static List_T
report_perfect_segment (int *found_score, int *nhits, List_T hits, Univcoord_T left,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			Chrpos_T chrlength, int querylength, Compress_T query_compress,
			int nmisses_allowed, bool plusp, int genestrand) {
  Stage3end_T hit;
  int nmismatches;

  if (snpp == true) {
    if ((hit = Stage3end_new_substitution(&(*found_score),/*nmismatches*/0,
					  left,/*genomiclength*/querylength,query_compress,
					  plusp,genestrand,chrnum,chroffset,chrhigh,chrlength,
					  /*sarrayp*/false)) == NULL) {
      return hits;
    } else {
      *nhits += 1;
      return List_push(hits,(void *) hit);
    }

  } else if (mode != STANDARD || spansize != index1part) {
    /* Count actual number of mismatches.  May not be a perfect segment. */
    nmismatches = Genome_count_mismatches_limit(query_compress,left,/*pos5*/0,/*pos3*/querylength,
						/*max_mismatches_allowed*/nmisses_allowed,
						plusp,genestrand);
    debug(printf("Got %d mismatches\n",nmismatches));
    if (nmismatches > nmisses_allowed) {
      return hits;
    } else {
      /* Don't use Stage3end_new_exact, because need to mark mismatches */
      if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					    left,/*genomiclength*/querylength,
					    query_compress,plusp,genestrand,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*sarrayp*/false)) == NULL) {
	return hits;
      } else {
	*nhits += 1;
	return List_push(hits,(void *) hit);
      }
    }

  } else {
    /* mode == STANDARD && spansize == index1part */
    if ((hit = Stage3end_new_exact(&(*found_score),left,/*genomiclength*/querylength,
				   query_compress,plusp,genestrand,
				   chrnum,chroffset,chrhigh,chrlength,/*sarrayp*/false)) == NULL) {
      return hits;
    } else {
      *nhits += 1;
      return List_push(hits,(void *) hit);
    }
  }
}


#if 0
static List_T
report_perfect_segment_dibase (int *found_score, int *nhits, List_T hits, Univcoord_T left, Univcoord_T diagonal,
			       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			       char *queryptr, int querylength, Compress_T query_compress,
			       int nmisses_allowed, bool plusp) {
  Stage3end_T hit;

#if 0
  int ncolordiffs;
  Dibase_count_mismatches_substring(&ncolordiffs,query,pos5,pos3,blocks,
				    /*startpos*/left+pos5,/*endpos*/left+pos3);
#endif

  /* Need to fill buffer with nucleotide genome anyway */
  if ((hit = Stage3end_new_substitution(&(*found_score),/*nmismatches*/0,
					left,/*genomiclength*/querylength,
					query_compress,plusp,genestrand,
					chrnum,chroffset,chrhigh,chrlength,
					/*sarrayp*/false)) == NULL) {
    return hits;
  } else {
    *nhits += 1;
    return List_push(hits,(void *) hit);
  }
}
#endif


/* Called only by exact/sub:1 procedures, so need to do Bigendian conversion */
#ifdef WORDS_BIGENDIAN
static int
binary_search_bigendian (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%llu\n",lowi,highi,(unsigned long long) goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%llu %d:%llu %d:%llu   vs. %llu\n",
		   lowi,(unsigned long long) Bigendian_convert_univcoord(positions[lowi]),
		   middlei,(unsigned long long) Bigendian_convert_univcoord(positions[middlei]),
		   highi,(unsigned long long) Bigendian_convert_univcoord(positions[highi]),goal));
    if (goal < Bigendian_convert_univcoord(positions[middlei])) {
      highi = middlei;
    } else if (goal > Bigendian_convert_univcoord(positions[middlei])) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}
#endif


#ifdef LARGE_GENOMES
static int
binary_search_large (int lowi, int highi, unsigned char *positions_high, UINT4 *positions_low, Univcoord_T goal) {
  int middlei;
  Univcoord_T position;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%llu\n",
		 lowi,highi,(unsigned long long) goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    position = ((Univcoord_T) positions_high[middlei] << 32) + positions_low[middlei];
    debug10(printf("  binary: %d:%llu %d:%llu %d:%llu   vs. %llu\n",
		   lowi,(unsigned long long) ((positions_high[lowi] << 32) + positions_low[lowi]),
		   middlei,(unsigned long long) position,
		   highi,(unsigned long long) ((positions_high[highi] << 32) + positions_low[highi]),
		   (unsigned long long) goal));
    if (goal < position) {
      highi = middlei;
    } else if (goal > position) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}
#endif


static int
binary_search (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%llu\n",
		 lowi,highi,(unsigned long long) goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%llu %d:%llu %d:%llu   vs. %llu\n",
		   lowi,(unsigned long long) positions[lowi],
		   middlei,(unsigned long long) positions[middlei],
		   highi,(unsigned long long) positions[highi],
		   (unsigned long long) goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


static int
binary_search_segments (int lowi, int highi, struct Segment_T *segments, Univcoord_T goal) {
  int middlei, middlei_up, middlei_down;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%llu\n",
		 lowi,highi,(unsigned long long) goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    if (segments[middlei].diagonal == (Univcoord_T) -1) {
      middlei_up = middlei + 1;
      middlei_down = middlei - 1;
    } else {
      middlei_up = middlei_down = middlei;
    }
    debug10(printf("  binary: %d:%llu %d:%llu %d:%llu   vs. %llu\n",
		   lowi,(unsigned long long) segments[lowi].diagonal,
		   middlei,(unsigned long long) segments[middlei].diagonal,
		   highi,(unsigned long long) segments[highi].diagonal,
		   (unsigned long long) goal));
    if (goal < segments[middlei_down].diagonal) {
      highi = middlei_down;
    } else if (goal > segments[middlei_up].diagonal) {
      lowi = middlei_up + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}



/* Generalization of identify_exact_iter and identify_onemiss_iter */
static List_T
identify_multimiss_iter (int *found_score, Chrnum_T *chrnum, Univcoord_T *chroffset, Univcoord_T *chrhigh,
			 Chrpos_T *chrlength, int *nhits, List_T hits, Univcoord_T goal, struct List_T *prev, int *nempty,
			 int *global_miss_querypos5, int *global_miss_querypos3,
			 int querylength, Compress_T query_compress, bool plusp, int genestrand,
			 int nmisses_allowed, int nmisses_seen, int miss_querypos5, int miss_querypos3) {
  List_T spanningset;
  Stage3end_T hit;
  Spanningelt_T elt;
  Compoundpos_T compoundpos;
  Univcoord_T local_goal, left;
  Univcoord_T position;
  int nmismatches, j;


  debug7(printf("identify_multimiss_iter on diagonal %llu with spanningset of length %d and %d misses seen initially\n",
		(unsigned long long) goal,List_length(prev),nmisses_seen));

  if (nmisses_seen > nmisses_allowed) {
    debug7(printf("Result: skipping because %d misses seen > %d allowed\n",nmisses_seen,nmisses_allowed));
    return hits;
  }

  for (spanningset = prev->rest; spanningset /* != NULL */; prev = spanningset, spanningset = spanningset->rest) {
    elt = (Spanningelt_T) spanningset->first;
    debug7(printf("nmisses seen %d, allowed %d, remaining %d, goal %llu: ",
		  nmisses_seen,nmisses_allowed,List_length(prev->rest),(unsigned long long) goal));

    if (elt->intersection_diagonals != NULL) {
      /* Intersection diagonals already computed */
      if (elt->intersection_ndiagonals > 0 && *elt->intersection_diagonals < goal) {
	debug7(printf("  (%d>>",elt->intersection_ndiagonals));
	j = 1;
	while (j < elt->intersection_ndiagonals && elt->intersection_diagonals[j] < goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= elt->intersection_ndiagonals) {
	  j = binary_search(j >> 1,elt->intersection_ndiagonals,elt->intersection_diagonals,goal);
	} else {
	  j = binary_search(j >> 1,j,elt->intersection_diagonals,goal);
	}
	elt->intersection_diagonals += j;
	elt->intersection_ndiagonals -= j;
	debug7(printf("  >>%d)",elt->intersection_ndiagonals));
      }

      if (elt->intersection_ndiagonals <= 0) {
	/* List is empty, so modify spanningset and continue with one more miss seen. */
	prev->rest = spanningset->rest;
	spanningset = prev;
	*nempty += 1;
	if (elt->miss_querypos5 < *global_miss_querypos5) *global_miss_querypos5 = elt->miss_querypos5;
	if (elt->miss_querypos3 > *global_miss_querypos3) *global_miss_querypos3 = elt->miss_querypos3;

	debug7(printf(" intersection empty, counts as one miss --"));
	if (++nmisses_seen > nmisses_allowed) {
	  debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	  return hits;
	} else {
	  debug7(printf("  nmisses seen %d <= allowed %d, so continuing (1)\n",nmisses_seen,nmisses_allowed));
	  if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	  /* continue; -- naturally falls to end of loop */
	}
      } else if (*elt->intersection_diagonals > local_goal) {
	/* Already advanced past goal, so continue with one more miss seen. */
	debug7(printf(" one miss --"));
	if (++nmisses_seen > nmisses_allowed) {
	  debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	  return hits;
	} else {
	  debug7(printf("  nmisses seen %d <= allowed %d, so continuing (2)\n",nmisses_seen,nmisses_allowed));
	  if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	  /* continue; -- naturally falls to end of loop */
	}
      } else {
	/* Found goal.  Advance past goal and continue with loop. */
	debug7(printf(" advancing\n"));
	++elt->intersection_diagonals;
	--elt->intersection_ndiagonals;
	/* continue; -- naturally falls to end of loop */
      }

    } else {
      if (elt->partnerp == true) {
	/* Partner is guaranteed to be atomic */
	local_goal = goal - elt->partner_diagterm;

#ifdef LARGE_GENOMES
	if (elt->partner_npositions > 0 &&
	    (((Univcoord_T) *elt->partner_positions_high) << 32) + (*elt->partner_positions_low) < local_goal) {
	  debug7(printf("  (%d>>",elt->partner_npositions));
	  j = 1;
	  while (j < elt->partner_npositions &&
		 ((Univcoord_T) elt->partner_positions_high[j] << 32) + elt->partner_positions_low[j] < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->partner_npositions) {
	    j = binary_search_large(j >> 1,elt->partner_npositions,elt->partner_positions_high,elt->partner_positions_low,local_goal);
	  } else {
	    j = binary_search_large(j >> 1,j,elt->partner_positions_high,elt->partner_positions_low,local_goal);
	  }
	  elt->partner_positions_high += j;
	  elt->partner_positions_low += j;
	  elt->partner_npositions -= j;
	  debug7(printf("  >>%d)",elt->partner_npositions));
	}
#elif defined(WORDS_BIGENDIAN)
	if (elt->partner_npositions > 0 && Bigendian_convert_univcoord(*elt->partner_positions) < local_goal) {
	  debug7(printf("  (%d>>",elt->partner_npositions));
	  j = 1;
	  while (j < elt->partner_npositions && Bigendian_convert_univcoord(elt->partner_positions[j]) < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->partner_npositions) {
	    j = binary_search_bigendian(j >> 1,elt->partner_npositions,elt->partner_positions,local_goal);
	  } else {
	    j = binary_search_bigendian(j >> 1,j,elt->partner_positions,local_goal);
	  }
	  elt->partner_positions += j;
	  elt->partner_npositions -= j;
	  debug7(printf("  >>%d)",elt->partner_npositions));
	}
#else
	if (elt->partner_npositions > 0 && *elt->partner_positions < local_goal) {
	  debug7(printf("  (%d>>",elt->partner_npositions));
	  j = 1;
	  while (j < elt->partner_npositions && elt->partner_positions[j] < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->partner_npositions) {
	    j = binary_search(j >> 1,elt->partner_npositions,elt->partner_positions,local_goal);
	  } else {
	    j = binary_search(j >> 1,j,elt->partner_positions,local_goal);
	  }
	  elt->partner_positions += j;
	  elt->partner_npositions -= j;
	  debug7(printf("  >>%d)",elt->partner_npositions));
	}
#endif

	if (elt->partner_npositions <= 0) {
	  /* Empty, so modify spanningset and continue with one more miss seen. */
	  prev->rest = spanningset->rest;
	  spanningset = prev;
	  *nempty += 1;
	  if (elt->miss_querypos5 < *global_miss_querypos5) *global_miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > *global_miss_querypos3) *global_miss_querypos3 = elt->miss_querypos3;

	  debug7(printf(" partner empty --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (3)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    continue;		/* Don't need to check main list below */
	  }
#ifdef LARGE_GENOMES
	} else if ((((Univcoord_T) *elt->partner_positions_high) << 32) + (*elt->partner_positions_low) > local_goal) {
	  /* Advanced past local_goal, so continue with one more miss seen. */
	  debug7(printf(" not in partner --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (4)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    continue;		/* Don't need to check main list below */
	  }

#elif defined(WORDS_BIGENDIAN)
	} else if (Bigendian_convert_univcoord(*elt->partner_positions) > local_goal) {
	  /* Advanced past local_goal, so continue with one more miss seen. */
	  debug7(printf(" not in partner --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (4)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    continue;		/* Don't need to check main list below */
	  }
#else
	} else if (*elt->partner_positions > local_goal) {
	  /* Advanced past local_goal, so continue with one more miss seen. */
	  debug7(printf(" not in partner --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (4)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    continue;		/* Don't need to check main list below */
	  }
#endif
	} else {
	  /* Found local_goal.  Advance past local_goal and continue with rest of compound querypos */
	  debug7(printf(" found in partner, so continue with rest of compound querypos\n"));
#ifdef LARGE_GENOMES
	  ++elt->partner_positions_high;
	  ++elt->partner_positions_low;
#else
	  ++elt->partner_positions;
#endif
	  --elt->partner_npositions;
	  /* Continue below with main list */
	}
      }

      if ((compoundpos = elt->compoundpos) != NULL) {
	local_goal = goal - elt->compoundpos_diagterm;
	if (Compoundpos_search(&position,compoundpos,local_goal) <= 0) {
	  /* Empty, so modify spanningset and continue with one more miss seen. */
	  prev->rest = spanningset->rest;
	  spanningset = prev;
	  *nempty += 1;
	  if (elt->miss_querypos5 < *global_miss_querypos5) *global_miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > *global_miss_querypos3) *global_miss_querypos3 = elt->miss_querypos3;

	  debug7(printf("  compoundpos empty --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (5)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
	} else if (position > local_goal) {
	  /* Advanced past goal.  Continue with one more miss seen. */
	  debug7(printf("  compoundpos failed %llu > %llu --",(unsigned long long) position,(unsigned long long) local_goal));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (6)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
	} else {
	  /* Found goal.  Advance past goal and continue with loop.  */
	  debug7(printf("  found %llu, advancing...",(unsigned long long) local_goal));
	  /* continue; -- Naturally falls to end of loop */
	}

      } else {
	/* Ordinary querypos */
	local_goal = goal - elt->diagterm;

#ifdef LARGE_GENOMES
	if (elt->npositions > 0 && (((Univcoord_T) *elt->positions_high) << 32) + (*elt->positions_low) < local_goal) {
	  debug7(printf("  (%d>>",elt->npositions));
	  j = 1;
	  while (j < elt->npositions &&
		 ((Univcoord_T) elt->positions_high[j] << 32) + elt->positions_low[j] < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->npositions) {
	    j = binary_search_large(j >> 1,elt->npositions,elt->positions_high,elt->positions_low,local_goal);
	  } else {
	    j = binary_search_large(j >> 1,j,elt->positions_high,elt->positions_low,local_goal);
	  }
	  elt->positions_high += j;
	  elt->positions_low += j;
	  elt->npositions -= j;
	  debug7(printf("  >>%d)",elt->npositions));
	}

#elif defined(WORDS_BIGENDIAN)
	if (elt->npositions > 0 && Bigendian_convert_univcoord(*elt->positions) < local_goal) {
	  debug7(printf("  (%d>>",elt->npositions));
	  j = 1;
	  while (j < elt->npositions && Bigendian_convert_univcoord(elt->positions[j]) < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->npositions) {
	    j = binary_search_bigendian(j >> 1,elt->npositions,elt->positions,local_goal);
	  } else {
	    j = binary_search_bigendian(j >> 1,j,elt->positions,local_goal);
	  }
	  elt->positions += j;
	  elt->npositions -= j;
	  debug7(printf("  >>%d)",elt->npositions));
	}
#else
	if (elt->npositions > 0 && *elt->positions < local_goal) {
	  debug7(printf("  (%d>>",elt->npositions));
	  j = 1;
	  while (j < elt->npositions && elt->positions[j] < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->npositions) {
	    j = binary_search(j >> 1,elt->npositions,elt->positions,local_goal);
	  } else {
	    j = binary_search(j >> 1,j,elt->positions,local_goal);
	  }
	  elt->positions += j;
	  elt->npositions -= j;
	  debug7(printf("  >>%d)",elt->npositions));
	}
#endif

	if (elt->npositions <= 0) {
	  /* List is empty, so continue with one more miss seen. */
	  prev->rest = spanningset->rest;
	  spanningset = prev;
	  *nempty += 1;
	  if (elt->miss_querypos5 < *global_miss_querypos5) *global_miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > *global_miss_querypos3) *global_miss_querypos3 = elt->miss_querypos3;

	  debug7(printf(" positions empty, counts as one miss --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (7)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
#ifdef LARGE_GENOMES
	} else if ((((Univcoord_T) *elt->positions_high) << 32) + (*elt->positions_low) > local_goal) {
	  /* Already advanced past goal, so continue with one more miss seen. */
	  debug7(printf(" one miss %llu > %llu --",
			(unsigned long long) ((((Univcoord_T) *elt->positions_high) << 32) + (*elt->positions_low)),
			(unsigned long long) local_goal));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (8)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
#elif defined(WORDS_BIGENDIAN)
	} else if (Bigendian_convert_univcoord(*elt->positions) > local_goal) {
	  /* Already advanced past goal, so continue with one more miss seen. */
	  debug7(printf(" one miss %llu > %llu --",
			(unsigned long long) Bigendian_convert_univcoord(*elt->positions),(unsigned long long) local_goal));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (8)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
#else
	} else if (*elt->positions > local_goal) {
	  /* Already advanced past goal, so continue with one more miss seen. */
	  debug7(printf(" one miss %llu > %llu --",
			(unsigned long long) *elt->positions,(unsigned long long) local_goal));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing (8)\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
#endif
	} else {
	  /* Found goal.  Advance past goal and continue with loop. */
	  debug7(printf(" advancing\n"));
#ifdef LARGE_GENOMES
	  ++elt->positions_high;
	  ++elt->positions_low;
#else
	  ++elt->positions;
#endif
	  --elt->npositions;
	  /* continue; -- Naturally falls to end of loop */
	}
      }
    }
    /* End of loop */
  }

  /* success */
  debug7(printf("  successful candidate found, with >= %d misses, %d allowed\n",nmisses_seen,nmisses_allowed));
  if (nmisses_seen == 0) {
    left = goal - querylength;
    if (goal > *chrhigh) {
      *chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
      Univ_IIT_interval_bounds(&(*chroffset),&(*chrhigh),&(*chrlength),chromosome_iit,*chrnum,circular_typeint);
      /* *chrhigh += 1; */
    }
    debug(printf("Reporting perfect segment at left %llu and diagonal %llu, with chroffset %llu and chrhigh %llu\n",
		 (unsigned long long) left,(unsigned long long) goal,(unsigned long long) *chroffset,(unsigned long long) *chrhigh));
    if (goal > *chrhigh) {
      /* Query goes over end of chromosome */
      debug(printf("  Ignore: goes over end of chromosome\n"));
      return hits;
    } else {
      return report_perfect_segment(&(*found_score),&(*nhits),hits,left,*chrnum,*chroffset,*chrhigh,*chrlength,
				    querylength,query_compress,nmisses_allowed,plusp,genestrand);
    }
  } else {
    if (goal < (unsigned int) querylength) {
      debug7(printf("  Goes over beginning of chromosome\n"));
      return hits;
    } else {
      left = goal - querylength;	/* goal here is diagonal */
    }

    if (goal > *chrhigh) {
      *chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
      Univ_IIT_interval_bounds(&(*chroffset),&(*chrhigh),&(*chrlength),chromosome_iit,*chrnum,circular_typeint);
      /* *chrhigh += 1; */
    }
    if (goal > *chrhigh) {
      debug7(printf("  Goes over end of chromosome\n"));
      return hits;

    } else {
      if (snpp || mode != STANDARD || spansize != index1part) {
	debug7(printf("  Testing in entire query\n"));
	nmismatches = Genome_count_mismatches_substring(query_compress,left,/*pos5*/0,/*pos3*/querylength,
							plusp,genestrand);
      } else {
	debug7(printf("  Testing in query bounds %d..%d\n",miss_querypos5,miss_querypos3));
	nmismatches = Genome_count_mismatches_substring(query_compress,left,/*pos5*/miss_querypos5,/*pos3*/miss_querypos3,
							plusp,genestrand);

      }
      debug7(printf("nmismatches = %d (vs %d misses allowed)\n",nmismatches,nmisses_allowed));

      if (nmismatches > nmisses_allowed) {
	debug7(printf("Result: too many mismatches\n"));
	return hits;
      } else {
	debug7(printf("Result: successful hit saved\n"));
	debug(printf("Reporting hit with %d mismatches\n",nmismatches));
	if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					      left,/*genomiclength*/querylength,
					      query_compress,plusp,genestrand,
					      *chrnum,*chroffset,*chrhigh,*chrlength,
					      /*sarrayp*/false)) == NULL) {
	  return hits;
	} else {
	  *nhits += 1;
	  return List_push(hits,(void *) hit);
	}
      }
    }
  }
}


/* Since querypos -1 and query_lastpos+1 are now
   stored as compoundpos, we no longer want to use them for boosting */
static void
most_specific_oligomer_1 (int *best_plus_querypos, int *best_minus_querypos, T this,
			  int query_lastpos, Indexdb_T plus_indexdb, Indexdb_T minus_indexdb) {
  int querypos;
  int best_plus_count, best_minus_count;

  /* Not needed, since this is the first procedure called */
  /* Block_restore(this->block5); */

  best_plus_querypos[0] = -1;
  best_minus_querypos[0] = -1;

  best_plus_count = 0;
  best_minus_count = 0;

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
#if 0
    if (this->validp[querypos] == true) {
#endif
      this->plus_npositions[querypos] = Indexdb_count_no_subst(plus_indexdb,this->forward_oligos[querypos]);
      this->minus_npositions[querypos] = Indexdb_count_no_subst(minus_indexdb,this->revcomp_oligos[querypos]);
      debug(printf("Counting at querypos %d, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		   querypos,this->plus_npositions[querypos],this->forward_oligos[querypos],
		   this->minus_npositions[querypos],this->revcomp_oligos[querypos]));

      if (best_plus_querypos[0] < 0 || this->plus_npositions[querypos] < best_plus_count) {
	best_plus_querypos[0] = querypos;
	best_plus_count = this->plus_npositions[querypos];
      }
      if (best_minus_querypos[0] < 0 || this->minus_npositions[querypos] < best_minus_count) {
	best_minus_querypos[0] = querypos;
	best_minus_count = this->minus_npositions[querypos];
      }
#if 0
    }
#endif
  }

  return;
}


/* Since querypos -1 and query_lastpos+1 are now
   stored as compoundpos, we no longer want to use them for boosting */
static void
most_specific_oligomer_2 (int *best_plus_querypos, int *best_minus_querypos, T this,
			  int query_lastpos, Indexdb_T plus_indexdb, Indexdb_T minus_indexdb) {
  int querypos, mod;
  int best_plus_count[2], best_minus_count[2];

  /* Not needed, since this is the first procedure called */
  /* Block_restore(this->block5); */

  best_plus_querypos[0] = -2;
  best_plus_querypos[1] = -2;
  best_minus_querypos[0] = -2;
  best_minus_querypos[1] = -2;

  best_plus_count[0] = best_plus_count[1] = 0;
  best_minus_count[0] = best_minus_count[1] = 0;

#if 0
  if (this->validp[0] == false) {
    debug(printf("Not counting at querypos 0, neg 2 or neg 1 because validp is false at querypos 0\n"));
    this->plus_npositions[-1] = 0;
    this->minus_npositions[-1] = 0;
  } else {
#endif
    this->plus_npositions[-1] = Indexdb_count_left_subst_1(plus_indexdb,this->forward_oligos[0]);
    this->minus_npositions[-1] = Indexdb_count_right_subst_1(minus_indexdb,this->revcomp_oligos[0]);
    debug(printf("Counting at querypos 0, neg 1, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		 this->plus_npositions[-1],this->forward_oligos[0],this->minus_npositions[-1],this->revcomp_oligos[0]));
    best_plus_count[1] = this->plus_npositions[-1];
    best_minus_count[1] = this->minus_npositions[-1];
#if 0
  }
#endif

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
#if 0
    if (this->validp[querypos] == true) {
#endif
      mod = querypos % 2;
      this->plus_npositions[querypos] = Indexdb_count_no_subst(plus_indexdb,this->forward_oligos[querypos]);
      this->minus_npositions[querypos] = Indexdb_count_no_subst(minus_indexdb,this->revcomp_oligos[querypos]);
      debug(printf("Counting at querypos %d, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		   querypos,this->plus_npositions[querypos],this->forward_oligos[querypos],
		   this->minus_npositions[querypos],this->revcomp_oligos[querypos]));

      if (best_plus_querypos[mod] < 0 || this->plus_npositions[querypos] < best_plus_count[mod]) {
	best_plus_querypos[mod] = querypos;
	best_plus_count[mod] = this->plus_npositions[querypos];
      }
      if (best_minus_querypos[mod] < 0 || this->minus_npositions[querypos] < best_minus_count[mod]) {
	best_minus_querypos[mod] = querypos;
	best_minus_count[mod] = this->minus_npositions[querypos];
      }
#if 0
    }
#endif
  }

  querypos = query_lastpos;
#if 0
  if (this->validp[querypos] == false) {
    debug(printf("Not counting at querypos %d (pos 1) or %d (pos 2) because validp is false at querypos %d\n",
		 querypos+1,querypos+2,querypos));
    this->plus_npositions[querypos+1] = 0;
    this->minus_npositions[querypos+1] = 0;
  } else {
#endif
    mod = (querypos+1) % 2;
    this->plus_npositions[querypos+1] = Indexdb_count_right_subst_1(plus_indexdb,this->forward_oligos[querypos]);
    this->minus_npositions[querypos+1] = Indexdb_count_left_subst_1(minus_indexdb,this->revcomp_oligos[querypos]);
    debug(printf("Counting at querypos %d, pos 1, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		 querypos,this->plus_npositions[querypos+1],this->forward_oligos[querypos],
		 this->minus_npositions[querypos+1],this->revcomp_oligos[querypos]));

#if 0
    /* Don't want boostpos to be a compoundpos */
    if (best_plus_querypos[mod] < 0 || this->plus_npositions[querypos+1] < best_plus_count[mod]) {
      best_plus_querypos[mod] = querypos+1;
      best_plus_count[mod] = this->plus_npositions[querypos+1];
    }
    if (best_minus_querypos[mod] < 0 || this->minus_npositions[querypos+1] < best_minus_count[mod]) {
      best_minus_querypos[mod] = querypos+1;
      best_minus_count[mod] = this->minus_npositions[querypos+1];
    }
#endif

#if 0
  }
#endif

  return;
}

/* Since querypos -2, -1, query_lastpos+1, and query_lastpos+2 are now
   stored as compoundpos, we no longer want to use them for boosting */
static void
most_specific_oligomer_3 (int *best_plus_querypos, int *best_minus_querypos, T this,
			  int query_lastpos, Indexdb_T plus_indexdb, Indexdb_T minus_indexdb) {
  int querypos, mod;
  int best_plus_count[3], best_minus_count[3];

  /* Not needed, since this is the first procedure called */
  /* Block_restore(this->block5); */

  best_plus_querypos[0] = -3;
  best_plus_querypos[1] = -3;
  best_plus_querypos[2] = -3;
  best_minus_querypos[0] = -3;
  best_minus_querypos[1] = -3;
  best_minus_querypos[2] = -3;

  best_plus_count[0] = best_plus_count[1] = best_plus_count[2] = 0;
  best_minus_count[0] = best_minus_count[1] = best_minus_count[2] = 0;

#if 0
  if (this->validp[0] == false) {
    debug(printf("Not counting at querypos 0, neg 2 or neg 1 because validp is false at querypos 0\n"));
    this->plus_npositions[-2] = 0;
    this->plus_npositions[-1] = 0;
    this->minus_npositions[-2] = 0;
    this->minus_npositions[-1] = 0;
  } else {
#endif
    this->plus_npositions[-2] = Indexdb_count_left_subst_2(plus_indexdb,this->forward_oligos[0]);
    this->minus_npositions[-2] = Indexdb_count_right_subst_2(minus_indexdb,this->revcomp_oligos[0]);
    debug(printf("Counting at querypos 0, neg 2, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		 this->plus_npositions[-2],this->forward_oligos[0],this->minus_npositions[-2],this->revcomp_oligos[0]));
    best_plus_count[1] = this->plus_npositions[-2];
    best_minus_count[1] = this->minus_npositions[-2];

    this->plus_npositions[-1] = Indexdb_count_left_subst_1(plus_indexdb,this->forward_oligos[0]);
    this->minus_npositions[-1] = Indexdb_count_right_subst_1(minus_indexdb,this->revcomp_oligos[0]);
    debug(printf("Counting at querypos 0, neg 1, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		 this->plus_npositions[-1],this->forward_oligos[0],this->minus_npositions[-1],this->revcomp_oligos[0]));
    best_plus_count[2] = this->plus_npositions[-1];
    best_minus_count[2] = this->minus_npositions[-1];
#if 0
  }
#endif

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
#if 0
    if (this->validp[querypos] == true) {
#endif
      mod = querypos % 3;
      this->plus_npositions[querypos] = Indexdb_count_no_subst(plus_indexdb,this->forward_oligos[querypos]);
      this->minus_npositions[querypos] = Indexdb_count_no_subst(minus_indexdb,this->revcomp_oligos[querypos]);
      debug(printf("Counting at querypos %d, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		   querypos,this->plus_npositions[querypos],this->forward_oligos[querypos],
		   this->minus_npositions[querypos],this->revcomp_oligos[querypos]));

      if (best_plus_querypos[mod] < 0 || this->plus_npositions[querypos] < best_plus_count[mod]) {
	best_plus_querypos[mod] = querypos;
	best_plus_count[mod] = this->plus_npositions[querypos];
      }
      if (best_minus_querypos[mod] < 0 || this->minus_npositions[querypos] < best_minus_count[mod]) {
	best_minus_querypos[mod] = querypos;
	best_minus_count[mod] = this->minus_npositions[querypos];
      }
#if 0
    }
#endif
  }

  querypos = query_lastpos;
#if 0
  if (this->validp[querypos] == false) {
    debug(printf("Not counting at querypos %d (pos 1) or %d (pos 2) because validp is false at querypos %d\n",
		 querypos+1,querypos+2,querypos));
    this->plus_npositions[querypos+1] = 0;
    this->plus_npositions[querypos+2] = 0;
    this->minus_npositions[querypos+1] = 0;
    this->minus_npositions[querypos+2] = 0;
  } else {
#endif
    mod = (querypos+1) % 3;
    this->plus_npositions[querypos+1] = Indexdb_count_right_subst_1(plus_indexdb,this->forward_oligos[querypos]);
    this->minus_npositions[querypos+1] = Indexdb_count_left_subst_1(minus_indexdb,this->revcomp_oligos[querypos]);
    debug(printf("Counting at querypos %d, pos 1, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		 querypos,this->plus_npositions[querypos+1],this->forward_oligos[querypos],
		 this->minus_npositions[querypos+1],this->revcomp_oligos[querypos]));

#if 0
    /* Don't want boostpos to be a compoundpos */
    if (best_plus_querypos[mod] < 0 || this->plus_npositions[querypos+1] < best_plus_count[mod]) {
      best_plus_querypos[mod] = querypos+1;
      best_plus_count[mod] = this->plus_npositions[querypos+1];
    }
    if (best_minus_querypos[mod] < 0 || this->minus_npositions[querypos+1] < best_minus_count[mod]) {
      best_minus_querypos[mod] = querypos+1;
      best_minus_count[mod] = this->minus_npositions[querypos+1];
    }
#endif

    mod = (querypos+2) % 3;
    this->plus_npositions[querypos+2] = Indexdb_count_right_subst_2(plus_indexdb,this->forward_oligos[querypos]);
    this->minus_npositions[querypos+2] = Indexdb_count_left_subst_2(minus_indexdb,this->revcomp_oligos[querypos]);
    debug(printf("Counting at querypos %d, pos 2, plus_npositions = %d (oligo %010llX), minus_npositions = %d (oligo %010llX)\n",
		 querypos,this->plus_npositions[querypos+2],this->forward_oligos[querypos],
		 this->minus_npositions[querypos+2],this->revcomp_oligos[querypos]));

#if 0
    /* Don't want boostpos to be a compoundpos */
    if (best_plus_querypos[mod] < 0 || this->plus_npositions[querypos+2] < best_plus_count[mod]) {
      best_plus_querypos[mod] = querypos+2;
      best_plus_count[mod] = this->plus_npositions[querypos+2];
    }
    if (best_minus_querypos[mod] < 0 || this->minus_npositions[querypos+2] < best_minus_count[mod]) {
      best_minus_querypos[mod] = querypos+2;
      best_minus_count[mod] = this->minus_npositions[querypos+2];
    }
#endif

#if 0
  }
#endif

  return;
}


static List_T
find_spanning_exact_matches (int *found_score, int *nhits, List_T hits, T this, int genestrand,
			     int querylength, int query_lastpos, Indexdb_T plus_indexdb, Indexdb_T minus_indexdb,
			     Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  Spanningelt_T *array;
  List_T prev;
  int best_plus_querypos[MAX_INDEX1INTERVAL], best_minus_querypos[MAX_INDEX1INTERVAL];
  Univcoord_T *diagonals0, diagonal0;
#ifdef LARGE_GENOMES
  unsigned char *positions0_high;
  UINT4 *positions0_low;
#else
  Univcoord_T *positions0;
#endif
  int diagterm0, ndiagonals0, npositions0;
  int boostpos, mod, nelts, max_nelts, elti, minscore;
  int global_miss_querypos5, global_miss_querypos3, elt_miss_querypos5, elt_miss_querypos3;
  int nempty;
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;

  debug(printf("Starting find_spanning_exact_matches\n"));

  /* Use shortest list for candidate generation */
  if (index1interval == 3) {
    most_specific_oligomer_3(best_plus_querypos,best_minus_querypos,this,query_lastpos,plus_indexdb,minus_indexdb);
  } else if (index1interval == 2) {
    most_specific_oligomer_2(best_plus_querypos,best_minus_querypos,this,query_lastpos,plus_indexdb,minus_indexdb);
  } else {
    most_specific_oligomer_1(best_plus_querypos,best_minus_querypos,this,query_lastpos,plus_indexdb,minus_indexdb);
  }
  max_nelts = (querylength + spansize - 1)/spansize;

  /* Plus */
  for (mod = 0; mod < index1interval; mod++) {
    chrhigh = 0U;

    this->plus_spanningset[mod] = array = (Spanningelt_T *) MALLOC(max_nelts * sizeof(Spanningelt_T));
    this->plus_spanningset_allocated[mod] = (struct Spanningelt_T *) MALLOC(max_nelts * sizeof(struct Spanningelt_T));
    for (elti = 0; elti < max_nelts; elti++) {
      array[elti] = &(this->plus_spanningset_allocated[mod][elti]);
    }

    Spanningelt_set(&minscore,&nelts,array,this->forward_oligos,&this->plus_retrievedp,
#ifdef LARGE_GENOMES
		    &this->plus_positions_high,&this->plus_positions_low,
#else
		    &this->plus_positions,
#endif
		    &this->plus_npositions,plus_indexdb,query_lastpos,querylength,mod,/*plusp*/true);
    assert(nelts <= max_nelts);
    this->plus_spanningset_nelts[mod] = nelts;

    boostpos = best_plus_querypos[mod];
    debug(printf("exact_matches, plus mod %d: proposed boostpos is %d\n",mod,boostpos));
    if (this->plus_npositions[boostpos] < minscore &&
	this->plus_retrievedp[boostpos] == false) {
      /* Boost */
      qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);

      /* Get boost positions */
#ifdef LARGE_GENOMES
      this->plus_positions_low[boostpos] =
	Indexdb_read_inplace(&(this->plus_npositions[boostpos]),&(this->plus_positions_high[boostpos]),
			     plus_indexdb,this->forward_oligos[boostpos]);
#else
      this->plus_positions[boostpos] =
	Indexdb_read_inplace(&(this->plus_npositions[boostpos]),plus_indexdb,this->forward_oligos[boostpos]);
#endif
      this->plus_retrievedp[boostpos] = true;
#ifdef LARGE_GENOMES
      positions0_high = this->plus_positions_high[boostpos];
      positions0_low = this->plus_positions_low[boostpos];
#else
      positions0 = this->plus_positions[boostpos];
#endif
      npositions0 = this->plus_npositions[boostpos];
      diagterm0 = querylength - boostpos; /* FORMULA */

      debug(printf("*** find_spanning_exact_matches, plus mod %d, with boost @ %d (%d positions)\n",
		   mod,boostpos,npositions0);
	    Spanningelt_print_array(array,nelts));

#if 0
      prev = List_push(List_copy(sorted),(void *) NULL); /* Add a dummy list elt to front */
#else
      prev = (struct List_T *) MALLOCA((nelts + 1) * sizeof(struct List_T));
      List_fill_array_with_handle(prev,(void *) array,nelts);
#endif

      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--npositions0 >= 0 && nempty == 0 && *nhits <= maxpaths_search) {
#ifdef LARGE_GENOMES
	debug7(printf("diag0 %d:%llu+%d advancing\n",
		      npositions0,(unsigned long long) ((((Univcoord_T) *positions0_high++) << 32) + (*positions0_low++)),diagterm0));
	diagonal0 = (((Univcoord_T) *positions0_high++) << 32) + (*positions0_low++) + diagterm0;
#elif defined(WORDS_BIGENDIAN)
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,Bigendian_convert_univcoord(*positions0),diagterm0));
	diagonal0 = Bigendian_convert_univcoord(*positions0++) + diagterm0;
#else
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,(*positions0),diagterm0));
	diagonal0 = (*positions0++) + diagterm0;
#endif
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal0,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_fwd,
				       /*plusp*/true,genestrand,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3);
      }

      FREEA(prev);

    } else {
      qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
      if (nelts > 1) {
	qsort(&(array[1]),nelts-1,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
      }

      debug(printf("*** find_spanning_exact_matches, plus mod %d, no boosting\n",mod));
      debug(Spanningelt_print_array(this->plus_spanningset[mod],this->plus_spanningset_nelts[mod]));

      /* diagonals0 is now in correct endianness */
      diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) array[0],&elt_miss_querypos5,&elt_miss_querypos3);

#if 0
      prev = List_push(List_copy(sorted->rest),(void *) NULL); /* Add a dummy list elt to front */
#else
      prev = (struct List_T *) MALLOCA((nelts - 1 + 1) * sizeof(struct List_T));
      List_fill_array_with_handle(prev,(void *) &(array[1]),nelts-1);
#endif

      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--ndiagonals0 >= 0 && nempty == 0 && *nhits <= maxpaths_search) {
	debug7(printf("diag0 %d:%llu advancing\n",ndiagonals0,(unsigned long long) (*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,*diagonals0++,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_fwd,
				       /*plusp*/true,genestrand,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3);
      }

      FREEA(prev);
    }
  }

  /* Minus */
  for (mod = 0; mod < index1interval; mod++) {
    chrhigh = 0U;

    this->minus_spanningset[mod] = array = (Spanningelt_T *) MALLOC(max_nelts * sizeof(Spanningelt_T));
    this->minus_spanningset_allocated[mod] = (struct Spanningelt_T *) MALLOC(max_nelts * sizeof(struct Spanningelt_T));
    for (elti = 0; elti < max_nelts; elti++) {
      array[elti] = &(this->minus_spanningset_allocated[mod][elti]);
    }
    Spanningelt_set(&minscore,&nelts,array,this->revcomp_oligos,&this->minus_retrievedp,
#ifdef LARGE_GENOMES
		    &this->minus_positions_high,&this->minus_positions_low,
#else
		    &this->minus_positions,
#endif
		    &this->minus_npositions,minus_indexdb,query_lastpos,querylength,mod,/*plusp*/false);
    assert(nelts <= max_nelts);
    this->minus_spanningset_nelts[mod] = nelts;

    boostpos = best_minus_querypos[mod];
    debug(printf("exact_matches, minus mod %d: proposed boostpos is %d\n",mod,boostpos));
    if (this->minus_npositions[boostpos] < minscore &&
	this->minus_retrievedp[boostpos] == false) {
      /* Boost */
      qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);

      /* Get boost positions */
#ifdef LARGE_GENOMES
      this->minus_positions_low[boostpos] =
	Indexdb_read_inplace(&(this->minus_npositions[boostpos]),&(this->minus_positions_high[boostpos]),
			     minus_indexdb,this->revcomp_oligos[boostpos]);
#else
      this->minus_positions[boostpos] =
	Indexdb_read_inplace(&(this->minus_npositions[boostpos]),minus_indexdb,this->revcomp_oligos[boostpos]);
#endif
      this->minus_retrievedp[boostpos] = true;
#ifdef LARGE_GENOMES
      positions0_high = this->minus_positions_high[boostpos];
      positions0_low = this->minus_positions_low[boostpos];
#else
      positions0 = this->minus_positions[boostpos];
#endif
      npositions0 = this->minus_npositions[boostpos];
      diagterm0 = boostpos + index1part; /* FORMULA */

      debug(printf("*** find_spanning_exact_matches, minus mod %d, with boost @ %d (%d positions)\n",
		   mod,boostpos,npositions0);
	    Spanningelt_print_array(array,nelts));

#if 0
      prev = List_push(List_copy(sorted),(void *) NULL);/* Add a dummy list elt to front */
#else
      prev = (struct List_T *) MALLOCA((nelts + 1) * sizeof(struct List_T));
      List_fill_array_with_handle(prev,(void *) array,nelts);
#endif

      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--npositions0 >= 0 && nempty == 0 && *nhits <= maxpaths_search) {
#ifdef LARGE_GENOMES
	debug7(printf("diag0 %d:%llu+%d advancing\n",
		      npositions0,(unsigned long long) ((((Univcoord_T) *positions0_high++) << 32) + (*positions0_low++)),diagterm0));
	diagonal0 = (((Univcoord_T) *positions0_high++) << 32) + (*positions0_low++) + diagterm0;
#elif defined(WORDS_BIGENDIAN)
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,Bigendian_convert_univcoord(*positions0),diagterm0));
	diagonal0 = Bigendian_convert_univcoord(*positions0++) + diagterm0;
#else
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,(*positions0),diagterm0));
	diagonal0 = (*positions0++) + diagterm0;
#endif
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal0,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_rev,
				       /*plusp*/false,genestrand,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3);
      }

      FREEA(prev);

    } else {
      qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
      if (nelts > 1) {
	qsort(&(array[1]),nelts-1,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
      }

      debug(printf("*** find_spanning_exact_matches, minus mod %d, no boosting\n",mod));
      debug(Spanningelt_print_array(this->minus_spanningset[mod],this->minus_spanningset_nelts[mod]));

      /* diagonals0 is now in correct endianness */
      diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) array[0],&elt_miss_querypos5,&elt_miss_querypos3);

#if 0
      prev = List_push(List_copy(sorted->rest),(void *) NULL); /* Add a dummy list elt to front */
#else
      prev = (struct List_T *) MALLOCA((nelts - 1 + 1) * sizeof(struct List_T));
      List_fill_array_with_handle(prev,(void *) &(array[1]),nelts-1);
#endif

      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--ndiagonals0 >= 0 && nempty == 0 && *nhits <= maxpaths_search) {
	debug7(printf("diag0 %d:%llu advancing\n",ndiagonals0,(unsigned long long) (*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,*diagonals0++,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_rev,
				       /*plusp*/false,genestrand,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3);
      }

      FREEA(prev);
    }
  }

  return hits;
}


static List_T
find_spanning_onemiss_matches (int *found_score, int *nhits, List_T hits, T this, int genestrand,
			       int querylength, Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  Spanningelt_T *array;
  List_T prev;
  Univcoord_T *diagonals0, *diagonals1, diagonal0, diagonal1;
  int global_miss_querypos5, global_miss_querypos3;
  int miss0_querypos5, miss0_querypos3, miss1_querypos5, miss1_querypos3;
  int mod, nelts, elti;
  int ndiagonals0, ndiagonals1;
  int nempty;
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;

  debug(printf("Starting find_spanning_onemiss_matches\n"));

  /* Plus */
  for (mod = 0; mod < index1interval; mod++) {
    debug(printf("Onemiss plus mod %d\n",mod));

    array = this->plus_spanningset[mod];
    nelts = this->plus_spanningset_nelts[mod];

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > 2) {
      qsort(&(array[2]),nelts-2,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    for (elti = 0; elti < nelts; elti++) {
      Spanningelt_reset(array[elti]);
    }

    debug(printf("*** find_spanning_onemiss_matches, plus mod %d\n",mod));
    debug(Spanningelt_print_array(array,nelts));

    /* diagonals0 and diagonals1 are now in correct endianness */
    diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) array[0],&miss0_querypos5,&miss0_querypos3);
    diagonals1 = Spanningelt_diagonals(&ndiagonals1,(Spanningelt_T) array[1],&miss1_querypos5,&miss1_querypos3);

#if 0
    prev = List_push(List_copy(sorted->rest->rest),(void *) NULL); /* Add a dummy list elt to front */
#else
    prev = (struct List_T *) MALLOCA((nelts - 2 + 1) * sizeof(struct List_T));
    List_fill_array_with_handle(prev,(void *) &(array[2]),nelts-2);
#endif

    nempty = 0;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    chrhigh = 0U;

    while (ndiagonals0 > 0 && ndiagonals1 > 0 && nempty <= 1 && *nhits <= maxpaths_search) {
      if ((diagonal0 = (*diagonals0)) < (diagonal1 = (*diagonals1))) {
	debug7(printf("diag0 %d:%llu advancing\n",ndiagonals0,(unsigned long long) diagonal0));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal0,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_fwd,
				       /*plusp*/true,genestrand,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3);
	++diagonals0;
	--ndiagonals0;

      } else if (diagonal1 < diagonal0) {
	debug7(printf("diag1 %d:%llu advancing\n",ndiagonals1,(unsigned long long) diagonal1));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal1,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_fwd,
				       /*plusp*/true,genestrand,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3);
	++diagonals1;
	--ndiagonals1;

      } else {
	debug7(printf("diag0&1 %d:%llu == %d:%llu advancing\n",
		      ndiagonals0,(unsigned long long) diagonal0,ndiagonals1,(unsigned long long) diagonal1));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal0,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_fwd,
				       /*plusp*/true,genestrand,/*nmisses_allowed*/1,
				       /*nmisses_seen*/nempty,global_miss_querypos5,global_miss_querypos3);
	++diagonals0;
	--ndiagonals0;
	++diagonals1;
	--ndiagonals1;
      }
    }

    while (--ndiagonals0 >= 0 && nempty == 0 && *nhits <= maxpaths_search) {
      debug7(printf("diag0 %d:%llu advancing\n",ndiagonals0+1,(unsigned long long) (*diagonals0)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,*diagonals0++,
				     prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     querylength,/*query_compress*/query_compress_fwd,
				     /*plusp*/true,genestrand,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3);
    }

    while (--ndiagonals1 >= 0 && nempty == 0 && *nhits <= maxpaths_search) {
      debug7(printf("diag1 %d:%llu advancing\n",ndiagonals1+1,(unsigned long long) (*diagonals1)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,*diagonals1++,
				     prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     querylength,/*query_compress*/query_compress_fwd,
				     /*plusp*/true,genestrand,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3);
    }

    FREEA(prev);
  }

  /* Minus */
  for (mod = 0; mod < index1interval; mod++) {
    debug(printf("Onemiss minus mod %d\n",mod));

    array = this->minus_spanningset[mod];
    nelts = this->minus_spanningset_nelts[mod];

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > 2) {
      qsort(&(array[2]),nelts-2,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    for (elti = 0; elti < nelts; elti++) {
      Spanningelt_reset(array[elti]);
    }

    debug(printf("*** find_spanning_onemiss_matches, minus mod %d\n",mod));
    debug(Spanningelt_print_array(array,nelts));

    /* diagonals0 and diagonals1 are now in correct endianness */
    diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) array[0],&miss0_querypos5,&miss0_querypos3);
    diagonals1 = Spanningelt_diagonals(&ndiagonals1,(Spanningelt_T) array[1],&miss1_querypos5,&miss1_querypos3);

#if 0
    prev = List_push(List_copy(sorted->rest->rest),(void *) NULL); /* Add a dummy list to front */
#else
    prev = (struct List_T *) MALLOCA((nelts - 2 + 1) * sizeof(struct List_T));
    List_fill_array_with_handle(prev,(void *) &(array[2]),nelts-2);
#endif

    nempty = 0;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    chrhigh = 0U;

    while (ndiagonals0 > 0 && ndiagonals1 > 0 && nempty <= 1 && *nhits <= maxpaths_search) {
      if ((diagonal0 = (*diagonals0)) < (diagonal1 = (*diagonals1))) {
	debug7(printf("diag0 %d:%llu advancing\n",ndiagonals0,(unsigned long long) (*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal0,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_rev,
				       /*plusp*/false,genestrand,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3);
	++diagonals0;
	--ndiagonals0;

      } else if (diagonal1 < diagonal0) {
	debug7(printf("diag1 %d:%llu advancing\n",ndiagonals1,(unsigned long long) (*diagonals1)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal1,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_rev,
				       /*plusp*/false,genestrand,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3);
	++diagonals1;
	--ndiagonals1;

      } else {
	debug7(printf("diag0&1 %d:%llu == %d:%llu advancing\n",
		      ndiagonals0,(unsigned long long) diagonal0,ndiagonals1,(unsigned long long) diagonal1));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal0,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_rev,
				       /*plusp*/false,genestrand,/*nmisses_allowed*/1,
				       /*nmisses_seen*/nempty,global_miss_querypos5,global_miss_querypos3);
	++diagonals0;
	--ndiagonals0;
	++diagonals1;
	--ndiagonals1;
      }
    }

    while (--ndiagonals0 >= 0 && nempty == 0 && *nhits <= maxpaths_search) {
      debug7(printf("diag0 %d:%llu advancing\n",ndiagonals0+1,(unsigned long long) (*diagonals0)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,*diagonals0++,
				     prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     querylength,/*query_compress*/query_compress_rev,
				     /*plusp*/false,genestrand,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3);
    }

    while (--ndiagonals1 >= 0 && nempty == 0 && *nhits <= maxpaths_search) {
      debug7(printf("diag1 %d:%llu advancing\n",ndiagonals1+1,(unsigned long long) (*diagonals1)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,*diagonals1++,
				     prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     querylength,/*query_compress*/query_compress_rev,
				     /*plusp*/false,genestrand,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3);
    }

    FREEA(prev);
  }

  return hits;
}


static void
init_tree (Batch_T *losers, int heapsize) {
  Batch_T *winners;
  int nodei, lefti, righti;

  winners = (Batch_T *) MALLOCA((2*heapsize + 1) * sizeof(Batch_T));
  for (nodei = 2*heapsize - 1; nodei >= heapsize; nodei--) {
    winners[nodei] = losers[nodei];
  }

  for (nodei = heapsize - 1; nodei >= 1; --nodei) {
    lefti = LEFT(nodei);
    righti = lefti + 1;

    if (winners[lefti]->diagonal < winners[righti]->diagonal) {
      winners[nodei] = winners[lefti];
      losers[nodei] = winners[righti];
      
    } else {
      winners[nodei] = winners[righti];
      losers[nodei] = winners[lefti];
    }
  }

  losers[0] = winners[1];
  FREEA(winners);

  return;
}


#ifdef USE_HEAP
static List_T
find_spanning_multimiss_matches (int *found_score, int *nhits, List_T hits, T this, int genestrand,
				 int nrequired, int querylength, Compress_T query_compress_fwd, Compress_T query_compress_rev,
				 int nmisses_allowed) {
  Univcoord_T *diagonals, diagonal;
  Spanningelt_T *array;
  List_T prev;
  int nunion = nmisses_allowed + nrequired, nelts, elti;
  int heapsize, count, mod, i;
  int ndiagonals, nempty;
  int parenti;
  int global_miss_querypos5, global_miss_querypos3;
  int elt_miss_querypos5, elt_miss_querypos3;
  struct Batch_T *batchpool;

  struct Batch_T sentinel_struct;
  Batch_T *heap, sentinel;
  int smallesti, righti;

  Batch_T batch;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;
  Chrnum_T chrnum;

  debug(printf("Starting find_spanning_multimiss_matches with %d misses allowed\n",nmisses_allowed));

  sentinel_struct.diagonal = (Univcoord_T) -1; /* infinity */
  sentinel = &sentinel_struct;

  batchpool = (struct Batch_T *) MALLOCA(nunion * sizeof(struct Batch_T));
  heap = (Batch_T *) MALLOCA((2*(nunion+1)+1+1) * sizeof(Batch_T)); /* being liberal with allocation */


  /* Plus */
  for (mod = 0; mod < index1interval; mod++) {
    array = this->plus_spanningset[mod];
    nelts = this->plus_spanningset_nelts[mod];
    debug(printf("Multimiss plus mod %d, nelts %d\n",mod,nelts));

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > nunion) {
      qsort(&(array[nunion]),nelts-nunion,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    for (elti = 0; elti < nelts; elti++) {
      Spanningelt_reset(array[elti]);
    }

    debug(printf("*** find_spanning_multimiss_matches, %d misses allowed, plus mod %d\n",nmisses_allowed,mod));
    debug(Spanningelt_print_array(array,nelts));

    /* Put first few pointers into heap */
    heapsize = 0;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    for (elti = 0; elti < nelts && elti < nunion; elti++) {
      /* Get list as a special one, and perform conversion if necessary */
      diagonals = Spanningelt_diagonals(&ndiagonals,(Spanningelt_T) array[elti],&elt_miss_querypos5,&elt_miss_querypos3);
      if (elt_miss_querypos5 < global_miss_querypos5) global_miss_querypos5 = elt_miss_querypos5;
      if (elt_miss_querypos3 > global_miss_querypos3) global_miss_querypos3 = elt_miss_querypos3;

      batch = &(batchpool[elti]);
      debug(printf("Adding batch %d of size %d...",elti,ndiagonals));
      if (ndiagonals > 0) {
	Batch_init_simple(batch,diagonals,ndiagonals,querylength,/*querypos*/elti);
	if (batch->npositions > 0) {
	  debug(printf("inserting into heap"));
	  min_heap_insert_simple(heap,&heapsize,batch);
	}
      }
      debug(printf("\n"));
    }
    debug(printf("heapsize is %d\n",heapsize));
    if (heapsize > 0) {
#if 0
      prev = List_push(List_copy(spanningset),(void *) NULL); /* Add a dummy list elt to front */
#else
      prev = (struct List_T *) MALLOCA((nelts - elti + 1) * sizeof(struct List_T));
      List_fill_array_with_handle(prev,(void *) &(array[elti]),nelts - elti);
#endif
      nempty = 0;

      /* Set up rest of heap */
      for (i = heapsize+1; i <= 2*heapsize+1; i++) {
	heap[i] = sentinel;
      }

      debug7(printf("*** multimiss mod %d plus:\n",mod));

      /* Initialize loop */
      batch = heap[1];
      diagonal = batch->diagonal;
      count = 1;
      debug7(printf("at #%d, initial diagonal is %llu\n",batch->nodei,(unsigned long long) diagonal));

      /* Update batch */
      if (--batch->npositions <= 0) {
	/* Use last entry in heap for heapify */
	batch = heap[heapsize];
	heap[heapsize--] = sentinel;
      } else {
	/* Use this batch for heapify */
	/* These positions are diagonals, and already in correct endianness */
	batch->diagonal = *(++batch->positions);
      }

      /* Heapify down */
      debug6(printf("Starting heapify with %llu\n",(unsigned long long) diagonal));
      parenti = 1;
      smallesti = (heap[3]->diagonal < heap[2]->diagonal) ? 3 : 2;
      debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		    2,3,(unsigned long long) heap[2]->diagonal,(unsigned long long)heap[3]->diagonal));
      while (batch->diagonal > heap[smallesti]->diagonal) {
	heap[parenti] = heap[smallesti];
	parenti = smallesti;
	smallesti = LEFT(parenti);
	righti = smallesti+1;
	debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		      smallesti,righti,(unsigned long long) heap[smallesti]->diagonal,
		      (unsigned long long) heap[righti]->diagonal));
	if (heap[righti]->diagonal < heap[smallesti]->diagonal) {
	  smallesti = righti;
	}
      }
      heap[parenti] = batch;
      debug6(printf("Inserting at %d\n\n",parenti));

      /* Iterate through heap */
      chrhigh = 0U;
      while (heapsize > 0 && *nhits <= maxpaths_search) {
	batch = heap[1];

	if (batch->diagonal == diagonal) {
	  count++;
	  debug7(printf("at #%d, incrementing diagonal %llu to count %d\n",
			batch->nodei,(unsigned long long) diagonal,count));
	} else {
	  /* End of diagonal */
	  if (count >= nrequired) {
	    /* printf("Testing %d..%d\n",miss_querypos5,miss_querypos3); */
	    hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal,
					   prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
					   querylength,/*query_compress*/query_compress_fwd,
					   /*plusp*/true,genestrand,nmisses_allowed,
					   /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3);
	  }
	  diagonal = batch->diagonal;
	  count = 1;
	  debug7(printf("at #%d, next diagonal is %llu\n",batch->nodei,(unsigned long long) diagonal));
	}

	/* Update batch */
	if (--batch->npositions <= 0) {
	  /* Use last entry in heap for heapify */
	  batch = heap[heapsize];
	  heap[heapsize--] = sentinel;
	} else {
	  /* Use this batch for heapify */
	  /* These positions are diagonals, and already in correct endianness */
	  batch->diagonal = *(++batch->positions);
	}

	/* Heapify down */
	debug6(printf("Starting heapify with %llu\n",(unsigned long long) diagonal));
	parenti = 1;
	smallesti = (heap[3]->diagonal < heap[2]->diagonal) ? 3 : 2;
	debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		      2,3,(unsigned long long) heap[2]->diagonal,(unsigned long long) heap[3]->diagonal));
	while (batch->diagonal > heap[smallesti]->diagonal) {
	  heap[parenti] = heap[smallesti];
	  parenti = smallesti;
	  smallesti = LEFT(parenti);
	  righti = smallesti+1;
	  debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
			smallesti,righti,(unsigned long long) heap[smallesti]->diagonal,
			(unsigned long long) heap[righti]->diagonal));
	  if (heap[righti]->diagonal < heap[smallesti]->diagonal) {
	    smallesti = righti;
	  }
	}
	heap[parenti] = batch;
	debug6(printf("Inserting at %d\n\n",parenti));
      }

      /* Terminate loop */
      if (count >= nrequired && *nhits <= maxpaths_search) {
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_fwd,
				       /*plusp*/true,genestrand,nmisses_allowed,
				       /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3);
      }

      FREEA(prev);
    }
  }

  /* Minus */
  for (mod = 0; mod < index1interval; mod++) {
    array = this->minus_spanningset[mod];
    nelts = this->minus_spanningset_nelts[mod];
    debug(printf("Multimiss minus mod %d, nelts %d\n",mod,nelts));

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > nunion) {
      qsort(&(array[nunion]),nelts-nunion,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    for (elti = 0; elti < nelts; elti++) {
      Spanningelt_reset(array[elti]);
    }

    debug(printf("*** find_spanning_multimiss_matches, %d misses_allowed, minus mod %d\n",nmisses_allowed,mod));
    debug(Spanningelt_print_array(array,nelts));

    /* Put first few pointers into heap */
    heapsize = 0;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    for (elti = 0; elti < nelts && elti < nunion; elti++) {
      /* Get list as a special one, and perform conversion if necessary */
      diagonals = Spanningelt_diagonals(&ndiagonals,(Spanningelt_T) array[elti],&elt_miss_querypos5,&elt_miss_querypos3);
      if (elt_miss_querypos5 < global_miss_querypos5) global_miss_querypos5 = elt_miss_querypos5;
      if (elt_miss_querypos3 > global_miss_querypos3) global_miss_querypos3 = elt_miss_querypos3;

      batch = &(batchpool[elti]);
      debug(printf("Adding batch %d of size %d...",elti,ndiagonals));
      if (ndiagonals > 0) {
	Batch_init_simple(batch,diagonals,ndiagonals,querylength,/*querypos*/elti);
	if (batch->npositions > 0) {
	  debug(printf("inserting into heap"));
	  min_heap_insert_simple(heap,&heapsize,batch);
	}
      }
      debug(printf("\n"));
    }
    debug(printf("heapsize is %d\n",heapsize));

    if (heapsize > 0) {
#if 0
      prev = List_push(List_copy(spanningset),(void **) NULL); /* Add a dummy list elt to front */
#else
      prev = (struct List_T *) MALLOCA((nelts - elti + 1) * sizeof(struct List_T));
      List_fill_array_with_handle(prev,(void *) &(array[elti]),nelts - elti);
#endif
      nempty = 0;

      /* Set up rest of heap */
      for (i = heapsize+1; i <= 2*heapsize+1; i++) {
	heap[i] = sentinel;
      }

      debug7(printf("*** multimiss mod %d minus:\n",mod));

      /* Initialize loop */
      batch = heap[1];
      diagonal = batch->diagonal;
      count = 1;
      debug7(printf("at #%d, initial diagonal is %llu\n",batch->nodei,(unsigned long long) diagonal));

      /* Update batch */
      if (--batch->npositions <= 0) {
	/* Use last entry in heap for heapify */
	batch = heap[heapsize];
	heap[heapsize--] = sentinel;
      } else {
	/* Use this batch for heapify */
	/* These positions are diagonals, and already in correct endianness */
	batch->diagonal = *(++batch->positions);
      }

      /* Heapify down */
      debug6(printf("Starting heapify with %llu\n",(unsigned long long) diagonal));
      parenti = 1;
      smallesti = (heap[3]->diagonal < heap[2]->diagonal) ? 3 : 2;
      debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		    2,3,(unsigned long long) heap[2]->diagonal,(unsigned long long) heap[3]->diagonal));
      while (batch->diagonal > heap[smallesti]->diagonal) {
	heap[parenti] = heap[smallesti];
	parenti = smallesti;
	smallesti = LEFT(parenti);
	righti = smallesti+1;
	debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		      smallesti,righti,(unsigned long long) heap[smallesti]->diagonal,
		      (unsigned long long) heap[righti]->diagonal));
	if (heap[righti]->diagonal < heap[smallesti]->diagonal) {
	  smallesti = righti;
	}
      }
      heap[parenti] = batch;
      debug6(printf("Inserting at %d\n\n",parenti));

      /* Iterate through heap */
      chrhigh = 0U;
      while (heapsize > 0 && *nhits <= maxpaths_search) {
	batch = heap[1];

	if (batch->diagonal == diagonal) {
	  count++;
	  debug7(printf("at #%d, incrementing diagonal %llu to count %d\n",
			batch->nodei,(unsigned long long) diagonal,count));
	} else {
	  /* End of diagonal */
	  if (count >= nrequired) {
	    hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal,
					   prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
					   querylength,/*query_compress*/query_compress_rev,
					   /*plusp*/false,genestrand,nmisses_allowed,
					   /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3);
	  }
	  diagonal = batch->diagonal;
	  count = 1;
	  debug7(printf("at #%d, next diagonal is %llu\n",batch->nodei,(unsigned long long) diagonal));
	}

	/* Update batch */
	if (--batch->npositions <= 0) {
	  /* Use last entry in heap for heapify */
	  batch = heap[heapsize];
	  heap[heapsize--] = sentinel;
	} else {
	  /* Use this batch for heapify */
	  /* These positions are diagonals, and already in correct endianness */
	  batch->diagonal = *(++batch->positions);
	}

	/* Heapify down */
	debug6(printf("Starting heapify with %llu\n",(unsigned long long) diagonal));
	parenti = 1;
	smallesti = (heap[3]->diagonal < heap[2]->diagonal) ? 3 : 2;
	debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
		      2,3,(unsigned long long) heap[2]->diagonal,(unsigned long long) heap[3]->diagonal));
	while (batch->diagonal > heap[smallesti]->diagonal) {
	  heap[parenti] = heap[smallesti];
	  parenti = smallesti;
	  smallesti = LEFT(parenti);
	  righti = smallesti+1;
	  debug6(printf("Comparing left %d/right %d: %llu and %llu\n",
			smallesti,righti,(unsigned long long) heap[smallesti]->diagonal,
			(unsigned long long) heap[righti]->diagonal));
	  if (heap[righti]->diagonal < heap[smallesti]->diagonal) {
	    smallesti = righti;
	  }
	}
	heap[parenti] = batch;
	debug6(printf("Inserting at %d\n\n",parenti));
      }

      /* Terminate loop */
      if (count >= nrequired && *nhits <= maxpaths_search) {
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_rev,
				       /*plusp*/false,genestrand,nmisses_allowed,
				       /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3);
      }

      FREEA(prev);
    }
  }

  FREEA(heap);
  FREEA(batchpool);

  return hits;
}
#endif

#ifndef USE_HEAP
/* Uses a loser tree */
static List_T
find_spanning_multimiss_matches (int *found_score, int *nhits, List_T hits, T this, int genestrand,
				 int nrequired, int querylength, Compress_T query_compress_fwd, Compress_T query_compress_rev,
				 int nmisses_allowed) {
  Univcoord_T *diagonals, diagonal;
  Spanningelt_T *array;
  List_T prev;
  int nunion = nmisses_allowed + nrequired, nelts, elti;
  int heapsize, count, mod, i;
  int ndiagonals, nempty;
  int parenti;
  int global_miss_querypos5, global_miss_querypos3;
  int elt_miss_querypos5, elt_miss_querypos3;
  struct Batch_T *batchpool;

  Batch_T *losers, current;
  Univcoord_T current_diagonal, winner_diagonal;

  Batch_T batch;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;
  Chrnum_T chrnum;

  debug(printf("Starting find_spanning_multimiss_matches with %d misses allowed\n",nmisses_allowed));

  batchpool = (struct Batch_T *) MALLOCA(nunion * sizeof(struct Batch_T));
  losers = (Batch_T *) MALLOCA((2*(nunion+1)+1+1) * sizeof(Batch_T)); /* being liberal with allocation */


  /* Plus */
  for (mod = 0; mod < index1interval; mod++) {
    array = this->plus_spanningset[mod];
    nelts = this->plus_spanningset_nelts[mod];
    debug(printf("Multimiss plus mod %d, nelts %d\n",mod,nelts));

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > nunion) {
      qsort(&(array[nunion]),nelts-nunion,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    for (elti = 0; elti < nelts; elti++) {
      Spanningelt_reset(array[elti]);
    }

    debug(printf("*** find_spanning_multimiss_matches, %d misses allowed, plus mod %d\n",nmisses_allowed,mod));
    debug(Spanningelt_print_array(array,nelts));

    /* Put first few pointers into heap */
    heapsize = (nelts < nunion) ? nelts : nunion;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    for (elti = 0; elti < nelts && elti < nunion; elti++) {
      /* Get list as a special one, and perform conversion if necessary */
      diagonals = Spanningelt_diagonals(&ndiagonals,(Spanningelt_T) array[elti],&elt_miss_querypos5,&elt_miss_querypos3);
      if (elt_miss_querypos5 < global_miss_querypos5) global_miss_querypos5 = elt_miss_querypos5;
      if (elt_miss_querypos3 > global_miss_querypos3) global_miss_querypos3 = elt_miss_querypos3;

      batch = &(batchpool[elti]);
      debug(printf("Adding batch %d of size %d...",elti,ndiagonals));
      Batch_init_simple(batch,diagonals,ndiagonals,querylength,/*nodei*/heapsize + elti);
      losers[heapsize + elti] = batch;
      debug(printf("\n"));
    }
    debug(printf("heapsize is %d\n",heapsize));
    if (heapsize > 0) {
#if 0
      prev = List_push(List_copy(spanningset),(void *) NULL); /* Add a dummy list elt to front */
#else
      prev = (struct List_T *) MALLOCA((nelts - elti + 1) * sizeof(struct List_T));
      List_fill_array_with_handle(prev,(void *) &(array[elti]),nelts - elti);
#endif
      nempty = 0;

      init_tree(losers,heapsize);

      debug7(printf("*** multimiss mod %d plus:\n",mod));

      /* Initialize loop */
      batch = losers[0];
      diagonal = batch->diagonal;
      count = 1;
      debug7(printf("at #%d, initial diagonal is %llu\n",batch->nodei,(unsigned long long) diagonal));

      /* Update batch */
      if (--batch->npositions <= 0) {
	batch->diagonal = -1U;
      } else {
	/* Use this batch for heapify */
	/* These positions are diagonals, and already in correct endianness */
	batch->diagonal = *(++batch->positions);
      }

      /* Update tree */
      winner_diagonal = batch->diagonal;
      for (parenti = PARENT(batch->nodei); parenti > 0; parenti = PARENT(parenti)) {
	current = losers[parenti];
	current_diagonal = current->diagonal;
	if (current_diagonal < winner_diagonal) {
	  losers[parenti] = batch;
	  batch = current;
	  winner_diagonal = current_diagonal;
	}
      }
      losers[0] = batch;

      /* Iterate through heap */
      chrhigh = 0U;
      while ((batch = losers[0])->diagonal != -1U && *nhits <= maxpaths_search) {

	if (batch->diagonal == diagonal) {
	  count++;
	  debug7(printf("at #%d, incrementing diagonal %llu to count %d\n",
			batch->nodei,(unsigned long long) diagonal,count));
	} else {
	  /* End of diagonal */
	  if (count >= nrequired) {
	    /* printf("Testing %d..%d\n",miss_querypos5,miss_querypos3); */
	    hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal,
					   prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
					   querylength,/*query_compress*/query_compress_fwd,
					   /*plusp*/true,genestrand,nmisses_allowed,
					   /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3);
	  }
	  diagonal = batch->diagonal;
	  count = 1;
	  debug7(printf("at #%d, next diagonal is %llu\n",batch->nodei,(unsigned long long) diagonal));
	}

	/* Update batch */
	if (--batch->npositions <= 0) {
	  batch->diagonal = -1U;
	} else {
	  /* Use this batch for heapify */
	  /* These positions are diagonals, and already in correct endianness */
	  batch->diagonal = *(++batch->positions);
	}

	/* Update tree */
	winner_diagonal = batch->diagonal;
	for (parenti = PARENT(batch->nodei); parenti > 0; parenti = PARENT(parenti)) {
	  current = losers[parenti];
	  current_diagonal = current->diagonal;
	  if (current_diagonal < winner_diagonal) {
	    losers[parenti] = batch;
	    batch = current;
	    winner_diagonal = current_diagonal;
	  }
	}
	losers[0] = batch;
      }

      /* Terminate loop */
      if (count >= nrequired && *nhits <= maxpaths_search) {
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_fwd,
				       /*plusp*/true,genestrand,nmisses_allowed,
				       /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3);
      }

      FREEA(prev);
    }
  }

  /* Minus */
  for (mod = 0; mod < index1interval; mod++) {
    array = this->minus_spanningset[mod];
    nelts = this->minus_spanningset_nelts[mod];
    debug(printf("Multimiss minus mod %d, nelts %d\n",mod,nelts));

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > nunion) {
      qsort(&(array[nunion]),nelts-nunion,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    for (elti = 0; elti < nelts; elti++) {
      Spanningelt_reset(array[elti]);
    }

    debug(printf("*** find_spanning_multimiss_matches, %d misses_allowed, minus mod %d\n",nmisses_allowed,mod));
    debug(Spanningelt_print_array(array,nelts));

    /* Put first few pointers into heap */
    heapsize = (nelts < nunion) ? nelts : nunion;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    for (elti = 0; elti < nelts && elti < nunion; elti++) {
      /* Get list as a special one, and perform conversion if necessary */
      diagonals = Spanningelt_diagonals(&ndiagonals,(Spanningelt_T) array[elti],&elt_miss_querypos5,&elt_miss_querypos3);
      if (elt_miss_querypos5 < global_miss_querypos5) global_miss_querypos5 = elt_miss_querypos5;
      if (elt_miss_querypos3 > global_miss_querypos3) global_miss_querypos3 = elt_miss_querypos3;

      batch = &(batchpool[elti]);
      debug(printf("Adding batch %d of size %d...",elti,ndiagonals));
      Batch_init_simple(batch,diagonals,ndiagonals,querylength,/*querypos*/heapsize + elti);
      losers[heapsize + elti] = batch;
      debug(printf("\n"));
    }
    debug(printf("heapsize is %d\n",heapsize));

    if (heapsize > 0) {
#if 0
      prev = List_push(List_copy(spanningset),(void **) NULL); /* Add a dummy list elt to front */
#else
      prev = (struct List_T *) MALLOCA((nelts - elti + 1) * sizeof(struct List_T));
      List_fill_array_with_handle(prev,(void *) &(array[elti]),nelts - elti);
#endif
      nempty = 0;

      init_tree(losers,heapsize);

      debug7(printf("*** multimiss mod %d minus:\n",mod));

      /* Initialize loop */
      batch = losers[0];
      diagonal = batch->diagonal;
      count = 1;
      debug7(printf("at #%d, initial diagonal is %llu\n",batch->nodei,(unsigned long long) diagonal));

      /* Update batch */
      if (--batch->npositions <= 0) {
	batch->diagonal = -1U;
      } else {
	/* Use this batch for heapify */
	/* These positions are diagonals, and already in correct endianness */
	batch->diagonal = *(++batch->positions);
      }

      /* Update tree */
      winner_diagonal = batch->diagonal;
      for (parenti = PARENT(batch->nodei); parenti > 0; parenti = PARENT(parenti)) {
	current = losers[parenti];
	current_diagonal = current->diagonal;
	if (current_diagonal < winner_diagonal) {
	  losers[parenti] = batch;
	  batch = current;
	  winner_diagonal = current_diagonal;
	}
      }
      losers[0] = batch;

      /* Iterate through heap */
      chrhigh = 0U;
      while ((batch = losers[0])->diagonal != -1U && *nhits <= maxpaths_search) {

	if (batch->diagonal == diagonal) {
	  count++;
	  debug7(printf("at #%d, incrementing diagonal %llu to count %d\n",
			batch->nodei,(unsigned long long) diagonal,count));
	} else {
	  /* End of diagonal */
	  if (count >= nrequired) {
	    hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal,
					   prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
					   querylength,/*query_compress*/query_compress_rev,
					   /*plusp*/false,genestrand,nmisses_allowed,
					   /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3);
	  }
	  diagonal = batch->diagonal;
	  count = 1;
	  debug7(printf("at #%d, next diagonal is %llu\n",batch->nodei,(unsigned long long) diagonal));
	}

	/* Update batch */
	if (--batch->npositions <= 0) {
	  batch->diagonal = -1U;
	} else {
	  /* Use this batch for heapify */
	  /* These positions are diagonals, and already in correct endianness */
	  batch->diagonal = *(++batch->positions);
	}

	/* Update tree */
	winner_diagonal = batch->diagonal;
	for (parenti = PARENT(batch->nodei); parenti > 0; parenti = PARENT(parenti)) {
	  current = losers[parenti];
	  current_diagonal = current->diagonal;
	  if (current_diagonal < winner_diagonal) {
	    losers[parenti] = batch;
	    batch = current;
	    winner_diagonal = current_diagonal;
	  }
	}
	losers[0] = batch;
      }

      /* Terminate loop */
      if (count >= nrequired && *nhits <= maxpaths_search) {
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,&chrlength,&(*nhits),hits,diagonal,
				       prev,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       querylength,/*query_compress*/query_compress_rev,
				       /*plusp*/false,genestrand,nmisses_allowed,
				       /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3);
      }

      FREEA(prev);
    }
  }

  FREEA(losers);
  FREEA(batchpool);

  return hits;
}
#endif


/************************************************************************/


#if 0
static void
trim_ends_unknowns_only (int *trim5, int *trim3, char *sequence1, char *sequence2, int length) {
  int pos;

  pos = 0;
  while (pos < length && sequence2[pos] == OUTOFBOUNDS) {
    pos++;
  }
  debug8(printf("outofbounds: trim 5': at %d: %c != %c\n",pos,sequence2[pos],OUTOFBOUNDS));
  *trim5 = pos;

  pos = length-1;
  debug8(printf("outofbounds: trim 3': %d:%c\n",pos,sequence2[pos]));
  while (pos >= 0 && sequence2[pos] == OUTOFBOUNDS) {
    pos--;
  }
  *trim3 = pos+1;
  debug8(printf("outofbounds: trim 3': %d - %d\n",length,*trim3));
  *trim3 = length - (*trim3);

  debug8(
	 printf("At query ->: %.*s\n",length,sequence1);
	 printf("At genome->: %.*s\n",length,sequence2);
	 printf("%02d %02d    ->: ",*trim5,*trim3);
	 for (pos = 0; pos < *trim5; pos++) {
	   printf(" ");
	 }
	 for ( ; pos < length - (*trim3); pos++) {
	   printf("*");
	 }
	 for ( ; pos < length; pos++) {
	   printf(" ");
	 }
	 printf("\n");
	 );

  return;
}
#endif


/************************************************************************/


/* Returns a master pointer (segments) to the block of segments */
/* If end_indel_mismatches_allowed set to 0, won't save any segments for end indels. */
static List_T
find_complete_mm (int *found_score, int *nhits, List_T hits, Segment_T *anchor_segments, int nanchors,
		  int querylength, Compress_T query_compress,
		  int max_mismatches_allowed, bool plusp, int genestrand) {
  Stage3end_T hit;
  int nmismatches;
  Univcoord_T left;
  Segment_T segmenti, *p;

  for (p = &(anchor_segments[0]); p < &(anchor_segments[nanchors]); p++) {
    segmenti = *p;
    assert(segmenti->diagonal != (Univcoord_T) -1);
    if (segmenti->floor <= max_mismatches_allowed) {
      left = segmenti->diagonal - querylength;
      nmismatches = Genome_count_mismatches_limit(query_compress,left,/*pos5*/0,/*pos3*/querylength,
						  max_mismatches_allowed,plusp,genestrand);
      if (nmismatches <= max_mismatches_allowed) {
	if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					      left,/*genomiclength*/querylength,
					      query_compress,plusp,genestrand,segmenti->chrnum,
					      segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
					      /*sarrayp*/false)) != NULL) {
	  segmenti->usedp = true;
	  *nhits += 1;
	  hits = List_push(hits,(void *) hit);
	}
      }
    }
  }

  return hits;
}


#ifdef USE_HEAP
/* TODO: Change spliceable to be an attribute of the segment.  Then we
   can loop over anchor_segments only */
static struct Segment_T *
identify_all_segments (int *nsegments, Segment_T **anchor_segments, int *nanchors,
		       Segment_T **spliceable, int *nspliceable,
#ifdef LARGE_GENOMES
		       unsigned char **positions_high, UINT4 **positions_low,
#else
		       Univcoord_T **positions,
#endif
		       int *npositions, bool *omitted, int querylength, int query_lastpos,
		       Floors_T floors, bool plusp) {
  struct Segment_T *segments = NULL;
  Segment_T *all_segments, *ptr_all, *ptr_anchor, *dest, *src;
  int length_threshold;
  int n_all_segments, n;
  int nanchors_bymod[MAX_INDEX1INTERVAL], naccept_bymod[MAX_INDEX1INTERVAL];
  int mod;
  int k;

  struct Batch_T *batchpool;
  struct Batch_T sentinel_struct;
  Batch_T *heap, sentinel;
  int smallesti, righti;
  Batch_T batch;
  int heapsize = 0;
  int parenti, i;
  int querypos, first_querypos, last_querypos;
  int floor_left, floor_right, floor_incr;
  int floor, floor_xfirst, floor_xlast, *floors_from_xfirst, *floors_to_xlast;
  int *floors_from_neg3, *floors_to_pos3;
  /* int exclude_xfirst, exclude_xlast; */
  Univcoord_T diagonal, segment_left, last_diagonal, chroffset = 0U, chrhigh = 0U;
  Chrpos_T chrlength, max_distance;
  Chrnum_T chrnum = 1;
#ifdef OLD_FLOOR_ENDS
  int halfquerylength, halfquery_lastpos;
#endif

#ifdef DIAGONAL_ADD_QUERYPOS
  UINT8 diagonal_add_querypos;
#endif
  int total_npositions = 0;
  int joffset = 0, j;

#ifdef DEBUG
  Segment_T segment, *p;
#endif

  Segment_T ptr, ptr_chrstart;
  Segment_T *ptr_spliceable;
  bool last_spliceable_p = false;
  /* bool next_spliceable_p; */
#ifdef DEBUG19
  Segment_T ptr0;
#endif
#ifndef SLOW_CHR_UPDATE
  Univcoord_T goal;
  int nchromosomes_local = nchromosomes;
  Univcoord_T *chrhighs_local = chrhighs;
#endif

  Univcoord_T *splicesites_local, splicesites_static[1];
  int nsplicesites_local;

  debug(printf("*** Starting identify_all_segments on %s ***\n",plusp ? "plus" : "minus"));

  if (floors == NULL) {
    *nsegments = 0;
    *anchor_segments = (Segment_T *) NULL;
    *nanchors = 0;
    *spliceable = (Segment_T *) NULL;
    *nspliceable = 0;
    return (struct Segment_T *) NULL;
  }

  if (splicesites == NULL) {
    splicesites_local = splicesites_static;
    splicesites_local[0] = (Univcoord_T) -1;
    nsplicesites_local = 0;
  } else {
    splicesites_local = splicesites;
    nsplicesites_local = nsplicesites;
  }

#ifdef OLD_FLOOR_ENDS
  halfquerylength = querylength / 2;
  halfquery_lastpos = halfquerylength - index1part;
#endif

  /* Create sentinel */
#ifdef DIAGONAL_ADD_QUERYPOS
  sentinel_struct.diagonal_add_querypos = (UINT8) -1; /* infinity */
  sentinel_struct.diagonal_add_querypos <<= 32;
#else
  sentinel_struct.querypos = querylength; /* essentially infinity */
  sentinel_struct.diagonal = (Univcoord_T) -1; /* infinity */
#endif
  sentinel = &sentinel_struct;

  /* Set up batches */
  batchpool = (struct Batch_T *) MALLOCA((query_lastpos+1) * sizeof(struct Batch_T));
  heap = (Batch_T *) MALLOCA((2*(query_lastpos+1)+1+1) * sizeof(Batch_T));

  /* Don't add entries for compoundpos positions (skip querypos -2, -1, lastpos+1, lastpos+2) */
  if (plusp) {
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      } else if (npositions[querypos] > 0) {
	debug1(printf("Adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
	batch = &(batchpool[i]);
#ifdef LARGE_GENOMES
	Batch_init(batch,querypos,/*diagterm*/querylength - querypos,positions_high[querypos],positions_low[querypos],
		   npositions[querypos],querylength);
#else
	Batch_init(batch,querypos,/*diagterm*/querylength - querypos,positions[querypos],
		   npositions[querypos],querylength);
#endif
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      }
    }
  } else {
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      } else if (npositions[querypos] > 0) {
	debug1(printf("Adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
	batch = &(batchpool[i]);
#ifdef LARGE_GENOMES
	Batch_init(batch,querypos,/*diagterm*/querypos + index1part,positions_high[querypos],positions_low[querypos],
		   npositions[querypos],querylength);
#else
	Batch_init(batch,querypos,/*diagterm*/querypos + index1part,positions[querypos],
		   npositions[querypos],querylength);
#endif
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      }
    }
  }
  debug14(printf("Initial total_npositions = %d\n",total_npositions));


  if (i == 0) {
    FREEA(heap);
    FREEA(batchpool);
    *nsegments = 0;
    return (struct Segment_T *) NULL;
  }

  /* Set up rest of heap */
  for (i = heapsize+1; i <= 2*heapsize+1; i++) {
    heap[i] = sentinel;
  }

  /* Putting chr marker "segments" after each chromosome */
  segments = (struct Segment_T *) MALLOC((total_npositions + nchromosomes) * sizeof(struct Segment_T));
  ptr_chrstart = ptr = &(segments[0]);
  all_segments = (Segment_T *) MALLOC(total_npositions * sizeof(Segment_T));
  ptr_all = &(all_segments[0]);
  *anchor_segments = (Segment_T *) MALLOC(total_npositions * sizeof(Segment_T));
  ptr_anchor = &((*anchor_segments)[0]);
  if (overall_max_distance == 0) {
    ptr_spliceable = *spliceable = (Segment_T *) NULL;
  } else {
    ptr_spliceable = *spliceable = (Segment_T *) CALLOC(total_npositions,sizeof(Segment_T));
  }

  /*
  if ((exclude_xfirst = firstbound-2-index1part-max_end_insertions) < 3) {
    exclude_xfirst = 3;
  }
  if ((exclude_xlast = lastbound+1+max_end_insertions) > query_lastpos-3) {
    exclude_xlast = query_lastpos-3;
  }
  */

#if 0
  /* Should account for firstbound and lastbound */
  floors_from_xfirst = floors->scorefrom[/* xfirst_from = */ firstbound-index1interval+max_end_insertions];
  floors_to_xlast = floors->scoreto[/* xlast_to = */ lastbound+1+index1interval-index1part-max_end_insertions];
#else
  /* This was previously run in identify_all_segments and not in identify_all_segments_for_terminals */
  if (spansize /* +max_end_insertions */ > query_lastpos + index1interval) {
    floors_from_xfirst = floors->scorefrom[query_lastpos+index1interval];
  } else {
    floors_from_xfirst = floors->scorefrom[spansize /* +max_end_insertions */];
  }
  if (query_lastpos-spansize /* -max_end_insertions */ < -index1interval) {
    floors_to_xlast = floors->scoreto[-index1interval];
  } else {
    floors_to_xlast = floors->scoreto[query_lastpos-spansize /* -max_end_insertions */];
  }
#endif
  floors_from_neg3 = floors->scorefrom[-index1interval];
  floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];


  /* Initialize loop */
  batch = heap[1];
  first_querypos = last_querypos = querypos = batch->querypos;
  last_diagonal = diagonal = batch->diagonal;

  floor_incr = floors_from_neg3[first_querypos];
  floor = floor_incr;
  floor_xlast = floor_incr;
  floor_xfirst = floors_from_xfirst[first_querypos] /* floors->scorefrom[xfirst_from][first_querypos] */;

#ifdef OLD_FLOOR_ENDS
  if (querypos < halfquery_lastpos) {
    floor_left = floor_incr;
  } else {
    floor_left = floors->scorefrom[-index1interval][halfquery_lastpos];
  }
  if (querypos < halfquerylength) {
    floor_right = floors->scorefrom[halfquerylength-index1interval][query_lastpos];
  } else {
    floor_right = floors->scorefrom[halfquerylength-index1interval][first_querypos];
  }
#else
  floor_left = floor_incr;
#ifdef DEBUG1
  floor_right = -99;
#endif
#endif


  debug1(printf("multiple_mm_%s, diagonal %llu, querypos %d\n",
		plusp ? "plus" : "minus",(unsigned long long) diagonal,querypos));
  debug1(printf("first_querypos = %d => initial values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
	        first_querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));

  if (--batch->npositions <= 0) {
    /* Use last entry in heap for insertion */
    batch = heap[heapsize];
    querypos = batch->querypos;
    heap[heapsize--] = sentinel;

  } else {
    /* Use this batch for insertion (same querypos) */
#ifdef LARGE_GENOMES
    batch->diagonal = ((Univcoord_T) *(++batch->positions_high) << 32) + *(++batch->positions_low) + batch->diagterm;
#elif defined(WORDS_BIGENDIAN)
    batch->diagonal = Bigendian_convert_univcoord(*(++batch->positions)) + batch->diagterm;
#else
    batch->diagonal = *(++batch->positions) + batch->diagterm;
#endif
#ifdef DIAGONAL_ADD_QUERYPOS
    batch->diagonal_add_querypos = (UINT8) batch->diagonal;
    batch->diagonal_add_querypos <<= 32;
    batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif
  }

  /* heapify */
  parenti = 1;
#ifdef DIAGONAL_ADD_QUERYPOS
  diagonal_add_querypos = batch->diagonal_add_querypos;
  smallesti = (heap[3]->diagonal_add_querypos < heap[2]->diagonal_add_querypos) ? 3 : 2;
  while (diagonal_add_querypos > heap[smallesti]->diagonal_add_querypos) {
    heap[parenti] = heap[smallesti];
    parenti = smallesti;
    smallesti = LEFT(parenti);
    righti = smallesti+1;
    if (heap[righti]->diagonal_add_querypos < heap[smallesti]->diagonal_add_querypos) {
      smallesti = righti;
    }
  }
#else
  diagonal = batch->diagonal;
  smallesti = ((heap[3]->diagonal < heap[2]->diagonal) ||
	       ((heap[3]->diagonal == heap[2]->diagonal) &&
		(heap[3]->querypos < heap[2]->querypos))) ? 3 : 2;
  /* Note that diagonal/querypos will never exceed a sentinel diagonal/querypos */
  while (diagonal > heap[smallesti]->diagonal ||
	 (diagonal == heap[smallesti]->diagonal &&
	  querypos > heap[smallesti]->querypos)) {
    heap[parenti] = heap[smallesti];
    parenti = smallesti;
    smallesti = LEFT(parenti);
    righti = smallesti+1;
    if ((heap[righti]->diagonal < heap[smallesti]->diagonal) ||
		  ((heap[righti]->diagonal == heap[smallesti]->diagonal) &&
		   (heap[righti]->querypos < heap[smallesti]->querypos))) {
      smallesti = righti;
    }
  }
#endif
  heap[parenti] = batch;


  /* Continue after initialization */
  while (heapsize > 0) {
    batch = heap[1];
    querypos = batch->querypos;
    diagonal = batch->diagonal;
    debug14(printf("diagonal = %u, querypos = %d\n",last_diagonal,last_querypos));

    if (diagonal == last_diagonal) {
      /* Continuing exact match or substitution */
      floor_incr = floors->scorefrom[last_querypos][querypos];
      floor += floor_incr;
      floor_xfirst += floor_incr;
      floor_xlast += floor_incr;

#ifdef OLD_FLOOR_ENDS
      /* Why is this here?  Just set floor_left at start and floor_right at end. */
      if (querypos < halfquery_lastpos) {
	floor_left += floor_incr;
      } else if (last_querypos < halfquery_lastpos) {
	/* Finish floor_left */
	floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+index1interval];
      }
      if (querypos >= halfquerylength) {
	if (last_querypos < halfquerylength) {
	  /* Start floor_right */
	  floor_right = floors->scorefrom[halfquerylength-index1interval][querypos];
	} else {
	  floor_right += floor_incr;
	}
      }
#endif

      debug1(printf("diagonal %llu unchanged: last_querypos = %d, querypos = %d => floor increments by %d\n",
		    (unsigned long long) diagonal,last_querypos,querypos,floor_incr));
      debug1(printf("*multiple_mm_%s, diagonal %llu, querypos %d, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    plusp ? "plus" : "minus",(unsigned long long) diagonal,querypos,
		    floor,floor_xfirst,floor_xlast,floor_left,floor_right));

    } else {
      /* End of diagonal */
      floor_incr = floors_to_pos3[last_querypos]  /* floors->score[last_querypos][query_lastpos+index1interval] */;
      floor += floor_incr;
      floor_xfirst += floor_incr;
      floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */

#ifdef OLD_FLOOR_ENDS
      if (last_querypos < halfquery_lastpos) {
	floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+index1interval];
	floor_right = floors->scorefrom[halfquerylength-index1interval][query_lastpos+index1interval];
      }
      if (last_querypos >= halfquerylength) {
	floor_right += floor_incr;
      }
#else
      floor_right = floor_incr;
#endif

      debug1(printf("new diagonal %llu > last diagonal %llu: last_querypos = %d => final values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    (unsigned long long) diagonal,(unsigned long long) last_diagonal,last_querypos,
		    floor,floor_xfirst,floor_xlast,floor_left,floor_right));

      if (last_diagonal > chrhigh) {
	if (ptr > ptr_chrstart) {
	  /* Add chr marker segment */
	  debug14(printf("=== ptr %p > ptr_chrstart %p, so adding chr marker segment\n",ptr,ptr_chrstart));
	  ptr->diagonal = (Univcoord_T) -1;
	  ptr_chrstart = ++ptr;
	}

	/* update chromosome bounds, based on low end */
#ifdef SLOW_CHR_UPDATE
	chrnum = Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* chrhigh += 1; */
#else
	j = 1;
#ifdef NO_EXTENSIONS_BEFORE_ZERO
	goal = last_diagonal - querylength + 1;
#else
	goal = last_diagonal + 1;
#endif
	while (j < nchromosomes_local && chrhighs_local[j] < goal) {
	  j <<= 1;			/* gallop by 2 */
	}
	if (j >= nchromosomes_local) {
	  j = binary_search(j >> 1,nchromosomes_local,chrhighs_local,goal);
	} else {
	  j = binary_search(j >> 1,j,chrhighs_local,goal);
	}
	chrnum += j;
#ifdef DEBUG15
	if (chrnum != Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength)) {
	  fprintf(stderr,"Got chrnum %d, but wanted %d\n",
		  chrnum,Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength));
	  abort();
	}
#endif
	chroffset = chroffsets[chrnum-1];
	chrhigh = chrhighs[chrnum-1];
	chrlength = chrlengths[chrnum-1];
	chrhighs_local += j;
	nchromosomes_local -= j;
#endif
      }
      if (last_diagonal <= chrhigh) { /* FORMULA for high position */
	/* position of high end is within current chromosome */
	debug1(printf("  => multiple_mm, diagonal %llu, query %d..%d, chrbounds %llu..%llu, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		      (unsigned long long) last_diagonal,first_querypos,last_querypos,
		      (unsigned long long) chroffset,(unsigned long long) chrhigh,
		      floor,floor_xfirst,floor_xlast,floor_left,floor_right));

	/* Save segment, but first advance splicesites past segment_left */
	segment_left = last_diagonal - querylength;
	max_distance = overall_max_distance;
	if (splicesites_local[0] >= last_diagonal) {
	  ptr->splicesites_i = -1;
	} else if (Splicetrie_splicesite_p(segment_left,/*pos5*/1,/*pos3*/querylength) == false) {
	  ptr->splicesites_i = -1;
	} else {
	  if (splicesites_local[0] < segment_left) {
	    j = 1;
	    while (j < nsplicesites_local && splicesites_local[j] < segment_left) {
	      j <<= 1;		/* gallop by 2 */
	    }
	    if (j >= nsplicesites_local) {
	      j = binary_search(j >> 1,nsplicesites_local,splicesites_local,segment_left);
	    } else {
	      j = binary_search(j >> 1,j,splicesites_local,segment_left);
	    }
	    joffset += j;
	    splicesites_local += j;
	    nsplicesites_local -= j;
	  }
	    
	  if (splicesites_local[0] >= last_diagonal) {
	    ptr->splicesites_i = -1;
	  } else {
	    ptr->splicesites_i = joffset;
	    j = joffset;
	    while (j < nsplicesites && splicesites[j] < last_diagonal) {
	      if (splicedists[j] > max_distance) {
		max_distance = splicedists[j];
	      }
	      j++;
	    }
	  }
	}

	/* Save segment */
	ptr->diagonal = last_diagonal;
	ptr->chrnum = chrnum;
	ptr->chroffset = chroffset;
	ptr->chrhigh = chrhigh;
	ptr->chrlength = chrlength;
	ptr->querypos5 = first_querypos;
	ptr->querypos3 = last_querypos;

	/* FORMULA */
	if (plusp) {
	  ptr->lowpos = ptr->diagonal - querylength + ptr->querypos5;
	  ptr->highpos = ptr->diagonal - querylength + ptr->querypos3 + index1part;
	} else {
	  ptr->lowpos = ptr->diagonal - ptr->querypos3 - index1part - index1part;
	  ptr->highpos = ptr->diagonal - ptr->querypos5 - index1part;
	}

	ptr->floor = floor;
	ptr->floor_xfirst = floor_xfirst;
	ptr->floor_xlast = floor_xlast;
	ptr->floor_left = floor_left;
	ptr->floor_right = floor_right;
	ptr->leftmost = ptr->rightmost = -1;
	ptr->left_splice_p = ptr->right_splice_p = false;
	ptr->spliceable_low_p = last_spliceable_p;
	/* ptr->spliceable_high_p = false; */
#if 0
	ptr->leftspan = ptr->rightspan = -1;
#endif
	ptr->usedp = false;
	ptr->pairablep = false;

#if 0
	/* Not doing this, because the max_distance test is already good enough */
	if (plusp) {
	  /* For plus-strand splicing, require segmenti->querypos3 < segmentj->querypos5,
	     so if segmenti->querypos3 is too high, then it is not spliceable */
	  if (last_querypos > query_lastpos) {
	    /* Not spliceable */
	    last_spliceable_p = false;
	  } else if (diagonal <= last_diagonal + max_distance) {
	    *ptr_spliceable++ = ptr;
	    ptr->spliceable_high_p = last_spliceable_p = true;
	  }
	} else {
	  /* For minus-strand splicing, require segmenti->querypos5 > segmentj->querypos3,
	     so if segmenti->querypos5 is too low, then it is not spliceable */
	  if (first_querypos < index1part) {
	    /* Not spliceable */
	    last_spliceable_p = false;
	  } else if (diagonal <= last_diagonal + max_distance) {
	    *ptr_spliceable++ = ptr;
	    ptr->spliceable_high_p = last_spliceable_p = true;
	  }
	}
#endif
	if (diagonal <= last_diagonal + max_distance) {
	  *ptr_spliceable++ = ptr;
	  ptr->spliceable_high_p = last_spliceable_p = true;
	  debug4s(printf("%s diagonal %u is spliceable because next one is at %u\n",
			 plusp ? "plus" : "minus",last_diagonal,diagonal));
	} else {
	  ptr->spliceable_high_p = last_spliceable_p = false;
	  debug4s(printf("%s diagonal %u is not spliceable because next one is at %u\n",
			 plusp ? "plus" : "minus",last_diagonal,diagonal));
	}
	debug14(printf("Saving segment at %u (%u), query %d..%d",last_diagonal,last_diagonal-chroffset,ptr->querypos5,ptr->querypos3));
	*ptr_all++ = ptr;
	if (last_querypos >= first_querypos + /*min_segment_length*/1) {
	  debug14(printf(" ANCHOR"));
	  *ptr_anchor++ = ptr;
	}
	debug14(printf("\n"));
	ptr++;
      }

      /* Prepare next diagonal */
      first_querypos = querypos;
      last_diagonal = diagonal;
      floor_incr = floors_from_neg3[first_querypos] /* floors->score[-index1interval][first_querypos] */;
      floor = floor_incr;
      floor_xlast = floor_incr;
      floor_xfirst = floors_from_xfirst[first_querypos];  /* floors->score[xfirst_from][first_querypos]; */

#ifdef OLD_FLOOR_ENDS
      if (querypos < halfquery_lastpos) {
	floor_left = floor_incr;
      } else {
	floor_left = floors->scorefrom[-index1interval][halfquery_lastpos];
      }
      if (querypos < halfquerylength) {
	floor_right = floors->scorefrom[halfquerylength-index1interval][query_lastpos];
      } else {
	floor_right = floors->scorefrom[halfquerylength-index1interval][first_querypos];
      }
#else
      floor_left = floor_incr;
#ifdef DEBUG1
      floor_right = -99;	/* For debugging output */
#endif
#endif

      debug1(printf("*multiple_mm_%s, diagonal %llu, querypos %d\n",
		    plusp ? "plus" : "minus",(unsigned long long) diagonal,querypos));
      debug1(printf("start of diagonal %llu, first_querypos = %d => initial values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    (unsigned long long) diagonal,first_querypos,
		    floor,floor_xfirst,floor_xlast,floor_left,floor_right));

    }
    last_querypos = querypos;


    if (--batch->npositions <= 0) {
      /* Use last entry in heap for insertion */
      batch = heap[heapsize];
      querypos = batch->querypos;
      heap[heapsize--] = sentinel;

    } else {
      /* Use this batch for insertion (same querypos) */
#ifdef LARGE_GENOMES
    batch->diagonal = ((Univcoord_T) *(++batch->positions_high) << 32) + *(++batch->positions_low) + batch->diagterm;
#elif defined(WORDS_BIGENDIAN)
      batch->diagonal = Bigendian_convert_univcoord(*(++batch->positions)) + batch->diagterm;
#else
      batch->diagonal = *(++batch->positions) + batch->diagterm;
#endif
#ifdef DIAGONAL_ADD_QUERYPOS
      batch->diagonal_add_querypos = (UINT8) batch->diagonal;
      batch->diagonal_add_querypos <<= 32;
      batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif
    }

    /* heapify */
    parenti = 1;
#ifdef DIAGONAL_ADD_QUERYPOS
    diagonal_add_querypos = batch->diagonal_add_querypos;
    smallesti = (heap[3]->diagonal_add_querypos < heap[2]->diagonal_add_querypos) ? 3 : 2;
    while (diagonal_add_querypos > heap[smallesti]->diagonal_add_querypos) {
      heap[parenti] = heap[smallesti];
      parenti = smallesti;
      smallesti = LEFT(parenti);
      righti = smallesti+1;
      if (heap[righti]->diagonal_add_querypos < heap[smallesti]->diagonal_add_querypos) {
	smallesti = righti;
      }
    }
#else
    diagonal = batch->diagonal;
    smallesti = ((heap[3]->diagonal < heap[2]->diagonal) ||
		 ((heap[3]->diagonal == heap[2]->diagonal) &&
		  (heap[3]->querypos < heap[2]->querypos))) ? 3 : 2;
    /* Note that diagonal/querypos will never exceed a sentinel diagonal/querypos */
    while (diagonal > heap[smallesti]->diagonal ||
	   (diagonal == heap[smallesti]->diagonal &&
	    querypos > heap[smallesti]->querypos)) {
      heap[parenti] = heap[smallesti];
      parenti = smallesti;
      smallesti = LEFT(parenti);
      righti = smallesti+1;
      if ((heap[righti]->diagonal < heap[smallesti]->diagonal) ||
	  ((heap[righti]->diagonal == heap[smallesti]->diagonal) &&
	   (heap[righti]->querypos < heap[smallesti]->querypos))) {
	smallesti = righti;
      }
    }
#endif
    heap[parenti] = batch;
  }
  debug14(printf("diagonal = %u, querypos = %d\n",last_diagonal,last_querypos));
  debug14(printf("\n"));

  /* Terminate loop. */
  floor_incr = floors_to_pos3[last_querypos];   /* floors->score[last_querypos][query_lastpos+index1interval]; */
  floor += floor_incr;
  floor_xfirst += floor_incr;
  floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */

#ifdef OLD_FLOOR_ENDS
  if (last_querypos < halfquery_lastpos) {
    floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+index1interval];
    floor_right = floors->scorefrom[halfquerylength-index1interval][query_lastpos+index1interval];
  }
  if (last_querypos >= halfquerylength) {
    floor_right += floor_incr;
  }
#else
  floor_right = floor_incr;
#endif
  
  debug1(printf("no more diagonals: last_querypos = %d => terminal values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		last_querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));

  debug1(printf("last_diagonal %u vs chrhigh %u (looking for >)\n",last_diagonal,chrhigh));
  if (last_diagonal > chrhigh) {
    if (ptr > ptr_chrstart) {
      /* Add chr marker segment */
      debug14(printf("=== ptr %p > ptr_chrstart %p, so adding chr marker segment\n",ptr,ptr_chrstart));
      ptr->diagonal = (Univcoord_T) -1;
      ptr_chrstart = ++ptr;
    }

    /* update chromosome bounds, based on low end */
#ifdef SLOW_CHR_UPDATE
    chrnum = Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength);
    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
    /* chrhigh += 1; */
#else
    j = 1;
#ifdef NO_EXTENSIONS_BEFORE_ZERO
    goal = last_diagonal - querylength + 1;
#else
    goal = last_diagonal + 1;
#endif
    while (j < nchromosomes_local && chrhighs_local[j] < goal) {
      j <<= 1;			/* gallop by 2 */
    }
    if (j >= nchromosomes_local) {
      j = binary_search(j >> 1,nchromosomes_local,chrhighs_local,goal);
    } else {
      j = binary_search(j >> 1,j,chrhighs_local,goal);
    }
    chrnum += j;
#ifdef DEBUG15
    if (chrnum != Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength)) {
      fprintf(stderr,"Got chrnum %d, but wanted %d\n",
	      chrnum,Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength));
      abort();
    }
#endif
    chroffset = chroffsets[chrnum-1];
    chrhigh = chrhighs[chrnum-1];
    chrlength = chrlengths[chrnum-1];
    chrhighs_local += j;
    nchromosomes_local -= j;
#endif
  }

  debug1(printf("last_diagonal %u vs chrhigh %u (looking for <=)\n",last_diagonal,chrhigh));
  if (last_diagonal <= chrhigh) { /* FORMULA for high position */
    /* position of high end is within current chromosome */
    debug1(printf("  => multiple_mm, diagonal %llu, query %d..%d, chrbounds %llu..%llu, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		  (unsigned long long) last_diagonal,first_querypos,last_querypos,
		  (unsigned long long) chroffset,(unsigned long long) chrhigh,
		  floor,floor_xfirst,floor_xlast,floor_left,floor_right));

    /* Save segment, but first advance splicesites past segment_left */
    segment_left = last_diagonal - querylength;
#if 0
    /* Last segment is not spliceable */
    max_distance = overall_max_distance;
#endif
    if (splicesites_local[0] >= last_diagonal) {
      ptr->splicesites_i = -1;
    } else if (Splicetrie_splicesite_p(segment_left,/*pos5*/1,/*pos3*/querylength) == false) {
      ptr->splicesites_i = -1;
    } else {
      if (splicesites_local[0] < segment_left) {
	j = 1;
	while (j < nsplicesites_local && splicesites_local[j] < segment_left) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= nsplicesites_local) {
	  j = binary_search(j >> 1,nsplicesites_local,splicesites_local,segment_left);
	} else {
	  j = binary_search(j >> 1,j,splicesites_local,segment_left);
	}
	joffset += j;
	splicesites_local += j;
	nsplicesites_local -= j;
      }

      if (splicesites_local[0] >= last_diagonal) {
	ptr->splicesites_i = -1;
      } else {
	ptr->splicesites_i = joffset;
#if 0
	/* Last segment is not spliceable */
	if (splicedists[joffset] > overall_max_distance) {
	  max_distance = splicedists[joffset];
	}
#endif
      }
    }

    /* Save segment */
    ptr->diagonal = last_diagonal;
    ptr->chrnum = chrnum;
    ptr->chroffset = chroffset;
    ptr->chrhigh = chrhigh;
    ptr->chrlength = chrlength;
    ptr->querypos5 = first_querypos;
    ptr->querypos3 = last_querypos;

    /* FORMULA */
    if (plusp) {
      ptr->lowpos = ptr->diagonal - querylength + ptr->querypos5;
      ptr->highpos = ptr->diagonal - querylength + ptr->querypos3 + index1part;
    } else {
      ptr->lowpos = ptr->diagonal - ptr->querypos3 - index1part - index1part;
      ptr->highpos = ptr->diagonal - ptr->querypos5 - index1part;
    }

    ptr->floor = floor;
    ptr->floor_xfirst = floor_xfirst;
    ptr->floor_xlast = floor_xlast;
    ptr->floor_left = floor_left;
    ptr->floor_right = floor_right;
    ptr->leftmost = ptr->rightmost = -1;
    ptr->left_splice_p = ptr->right_splice_p = false;
    ptr->spliceable_low_p = last_spliceable_p;
    ptr->spliceable_high_p = false;
#if 0
    ptr->leftspan = ptr->rightspan = -1;
#endif
    ptr->usedp = false;
    ptr->pairablep = false;

    /* Last segment is not spliceable */
    debug14(printf("Saving segment at %u (%u), query %d..%d",last_diagonal,last_diagonal - chroffset,ptr->querypos5,ptr->querypos3));
    *ptr_all++ = ptr;
    if (last_querypos >= first_querypos + /*min_segment_length*/1) {
      debug14(printf(" ANCHOR"));
      *ptr_anchor++ = ptr;
    }
    debug14(printf("\n"));
    ptr++;
  }


  if (ptr > ptr_chrstart) {
    /* Final chr marker segment */
    debug14(printf("=== ptr %p > ptr_chrstart %p, so adding final chr marker segment\n",ptr,ptr_chrstart));
    ptr->diagonal = (Univcoord_T) -1;
    /* ptr_chrstart = */ ++ptr;
  }

#ifdef DEBUG19
  for (k = 0, ptr0 = segments; ptr0 < ptr; k++, ptr0++) {
    printf("%d %llu\n",k,(unsigned long long) ptr0->diagonal);
  }
  printf("total_npositions = %d, nchromosomes = %d\n",total_npositions,nchromosomes);
#endif

  FREEA(heap);
  FREEA(batchpool);

  /* Note: segments is in descending diagonal order.  Will need to
     reverse before solving middle deletions */

  *nsegments = ptr - segments;
  *nanchors = ptr_anchor - *anchor_segments;
  *nspliceable = ptr_spliceable - *spliceable;
  debug(printf("nsegments = %d, of which %d are spliceable (total_npositions = %d, nchromosomes = %d)\n",
	       *nsegments,*nspliceable,total_npositions,nchromosomes));
  debug1(printf("nsegments = %d, of which %d are spliceable (total_npositions = %d, nchromosomes = %d)\n",
		*nsegments,*nspliceable,total_npositions,nchromosomes));

  assert(*nsegments <= total_npositions + nchromosomes);
  assert(*nanchors <= total_npositions);
  assert(*nspliceable <= total_npositions);

  n_all_segments = ptr_all - all_segments;
  debug(printf("%d all segments\n",n_all_segments));
  debug(printf("%d anchor segments\n",*nanchors));

  if (n_all_segments <= max_anchors) {
    /* Might as well use all segments */
    FREE(*anchor_segments);
    *anchor_segments = all_segments;
    *nanchors = n_all_segments;

  } else if (*nanchors <= max_anchors) {
    /* Use only the good anchor segments */
    FREE(all_segments);

  } else {
    /* Need to limit anchor segments */
    FREE(all_segments);

    /* Treat each mod separately */
    qsort(*anchor_segments,*nanchors,sizeof(Segment_T),Segment_mod_length_cmp);

    mod = 0;
    i = 0;
    while (mod < index1interval) {
      j = i;
      while (j < *nanchors && (*anchor_segments)[j]->querypos5 % index1interval == mod) {
	j++;
      }
      nanchors_bymod[mod] = j - i;

      if (j - i <= max_anchors) {
	naccept_bymod[mod] = j - i;
      } else {
	k = i + max_anchors;
	length_threshold = (*anchor_segments)[k]->querypos3 - (*anchor_segments)[k]->querypos5;
	while (k < j && k < i + max_anchors + /*ties*/100 &&
	       (*anchor_segments)[k]->querypos3 - (*anchor_segments)[k]->querypos5 == length_threshold) {
	  k++;
	}
	naccept_bymod[mod] = k - i;
      }

      debug(printf("For mod %d, accepting %d out of %d anchor segments with length threshold %d\n",
		   mod,naccept_bymod[mod],nanchors_bymod[mod],length_threshold));
      i = j;
      mod++;
    }

    /* Move good anchors to start of array */
    dest = src = &((*anchor_segments)[0]);
    *nanchors = 0;
    for (mod = 0; mod < index1interval; mod++) {
      memmove((void *) dest,(void *) src,naccept_bymod[mod] * sizeof(Segment_T));
      dest += naccept_bymod[mod];
      src += nanchors_bymod[mod];
      *nanchors += naccept_bymod[mod];
    }

    /* Re-sort in diagonal order */
    qsort(*anchor_segments,*nanchors,sizeof(Segment_T),Segment_diagonal_cmp);
  }


#ifdef DEBUG19
  printf("%d total segments\n",*nsegments);
  for (ptr0 = segments; ptr0 < ptr; ptr0++) {
    printf("%u %d..%d\n",ptr0->diagonal,ptr0->querypos5,ptr0->querypos3);
  }
#endif

#ifdef DEBUG
  printf("%d selected anchor segments\n",*nanchors);
  for (p = &(*anchor_segments)[0]; p< &((*anchor_segments)[*nanchors]); p++) {
    segment = (Segment_T) *p;
    printf("%u %d..%d spliceable_low:%d spliceable_high:%d\n",
	   segment->diagonal,segment->querypos5,segment->querypos3,segment->spliceable_low_p,segment->spliceable_high_p);
  }
#endif

  return segments;
}
#endif



#ifndef USE_HEAP
/* Uses a loser tree */
/* TODO: Change spliceable to be an attribute of the segment.  Then we
   can loop over anchor_segments only */
static struct Segment_T *
identify_all_segments (int *nsegments, Segment_T **anchor_segments, int *nanchors,
		       Segment_T **spliceable, int *nspliceable,
#ifdef LARGE_GENOMES
		       unsigned char **positions_high, UINT4 **positions_low,
#else
		       Univcoord_T **positions,
#endif
		       int *npositions, bool *omitted, int querylength, int query_lastpos,
		       Floors_T floors, bool plusp) {
  struct Segment_T *segments = NULL;
  Segment_T *all_segments, *ptr_all, *ptr_anchor, *dest, *src;
  int length_threshold;
  int n_all_segments, n;
  int nanchors_bymod[MAX_INDEX1INTERVAL], naccept_bymod[MAX_INDEX1INTERVAL];
  int mod;
  int k;

  struct Batch_T *batchpool;
  Batch_T *losers, current;
  Batch_T batch;
  int heapsize;
  int parenti, elti, i;
  int querypos, first_querypos, last_querypos;
  int floor_left, floor_right, floor_incr;
  int floor, floor_xfirst, floor_xlast, *floors_from_xfirst, *floors_to_xlast;
  int *floors_from_neg3, *floors_to_pos3;
  /* int exclude_xfirst, exclude_xlast; */
  Univcoord_T diagonal, segment_left, last_diagonal, chroffset = 0U, chrhigh = 0U;
  Chrpos_T chrlength, max_distance;
  Chrnum_T chrnum = 1;
#ifdef OLD_FLOOR_ENDS
  int halfquerylength, halfquery_lastpos;
#endif

#ifdef DIAGONAL_ADD_QUERYPOS
  UINT8 winner_diagonal_add_querypos, current_diagonal_add_querypos;
#else
  Univcoord_T winner_diagonal, current_diagonal;
  int winner_querypos, current_querypos;
#endif
  int total_npositions = 0;
  int joffset = 0, j;

#ifdef DEBUG
  Segment_T segment, *p;
#endif

  Segment_T ptr, ptr_chrstart;
  Segment_T *ptr_spliceable;
  /* bool next_spliceable_p; */
#ifdef DEBUG19
  Segment_T ptr0;
#endif
#ifndef SLOW_CHR_UPDATE
  Univcoord_T goal;
  int nchromosomes_local = nchromosomes;
  Univcoord_T *chrhighs_local = chrhighs;
#endif

  Univcoord_T *splicesites_local, splicesites_static[1];
  int nsplicesites_local;

  debug(printf("*** Starting identify_all_segments on %s ***\n",plusp ? "plus" : "minus"));

  if (floors == NULL) {
    *nsegments = 0;
    *anchor_segments = (Segment_T *) NULL;
    *nanchors = 0;
    *spliceable = (Segment_T *) NULL;
    *nspliceable = 0;
    return (struct Segment_T *) NULL;
  }

  if (splicesites == NULL) {
    splicesites_local = splicesites_static;
    splicesites_local[0] = (Univcoord_T) -1;
    nsplicesites_local = 0;
  } else {
    splicesites_local = splicesites;
    nsplicesites_local = nsplicesites;
  }

#ifdef OLD_FLOOR_ENDS
  halfquerylength = querylength / 2;
  halfquery_lastpos = halfquerylength - index1part;
#endif

  heapsize = 0;
  for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
    if (omitted[querypos] == true) {
      /* Skip */
    } else if (npositions[querypos] > 0) {
      heapsize++;
    } else {
      /* Skip */
    }
  }
  if (heapsize == 0) {
    *nsegments = 0;
    return (struct Segment_T *) NULL;
  } else {
    /* Set up batches */
    batchpool = (struct Batch_T *) MALLOCA(heapsize * sizeof(struct Batch_T));
    losers = (Batch_T *) MALLOCA((2*heapsize+1) * sizeof(Batch_T));
  }

  /* Don't add entries for compoundpos positions (skip querypos -2, -1, lastpos+1, lastpos+2) */
  if (plusp) {
    elti = 0;
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      } else if (npositions[querypos] > 0) {
	debug1(printf("Adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
	batch = &(batchpool[elti]);
#ifdef LARGE_GENOMES
	Batch_init(batch,querypos,/*diagterm*/querylength - querypos,positions_high[querypos],positions_low[querypos],
		   npositions[querypos],querylength,/*nodei*/heapsize + elti);
#else
	Batch_init(batch,querypos,/*diagterm*/querylength - querypos,positions[querypos],
		   npositions[querypos],querylength,/*nodei*/heapsize + elti);
#endif
	losers[heapsize + elti++] = batch;
	total_npositions += npositions[querypos];
      } else {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      }
    }
  } else {
    elti = 0;
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      } else if (npositions[querypos] > 0) {
	debug1(printf("Adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
	batch = &(batchpool[elti]);
#ifdef LARGE_GENOMES
	Batch_init(batch,querypos,/*diagterm*/querypos + index1part,positions_high[querypos],positions_low[querypos],
		   npositions[querypos],querylength,/*nodei*/heapsize + elti);
#else
	Batch_init(batch,querypos,/*diagterm*/querypos + index1part,positions[querypos],
		   npositions[querypos],querylength,/*nodei*/heapsize + elti);
#endif
	losers[heapsize + elti++] = batch;
	total_npositions += npositions[querypos];
      } else {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      }
    }
  }
  debug14(printf("Initial total_npositions = %d\n",total_npositions));

  init_tree(losers,heapsize);


  /* Putting chr marker "segments" after each chromosome */
  segments = (struct Segment_T *) MALLOC((total_npositions + nchromosomes) * sizeof(struct Segment_T));
  ptr_chrstart = ptr = &(segments[0]);
  all_segments = (Segment_T *) MALLOC(total_npositions * sizeof(Segment_T));
  ptr_all = &(all_segments[0]);
  *anchor_segments = (Segment_T *) MALLOC(total_npositions * sizeof(Segment_T));
  ptr_anchor = &((*anchor_segments)[0]);
  if (overall_max_distance == 0) {
    ptr_spliceable = *spliceable = (Segment_T *) NULL;
  } else {
    ptr_spliceable = *spliceable = (Segment_T *) CALLOC(total_npositions,sizeof(Segment_T));
  }

  /*
  if ((exclude_xfirst = firstbound-2-index1part-max_end_insertions) < 3) {
    exclude_xfirst = 3;
  }
  if ((exclude_xlast = lastbound+1+max_end_insertions) > query_lastpos-3) {
    exclude_xlast = query_lastpos-3;
  }
  */

#if 0
  /* Should account for firstbound and lastbound */
  floors_from_xfirst = floors->scorefrom[/* xfirst_from = */ firstbound-index1interval+max_end_insertions];
  floors_to_xlast = floors->scoreto[/* xlast_to = */ lastbound+1+index1interval-index1part-max_end_insertions];
#else
  /* This was previously run in identify_all_segments and not in identify_all_segments_for_terminals */
  if (spansize /* +max_end_insertions */ > query_lastpos + index1interval) {
    floors_from_xfirst = floors->scorefrom[query_lastpos+index1interval];
  } else {
    floors_from_xfirst = floors->scorefrom[spansize /* +max_end_insertions */];
  }
  if (query_lastpos-spansize /* -max_end_insertions */ < -index1interval) {
    floors_to_xlast = floors->scoreto[-index1interval];
  } else {
    floors_to_xlast = floors->scoreto[query_lastpos-spansize /* -max_end_insertions */];
  }
#endif
  floors_from_neg3 = floors->scorefrom[-index1interval];
  floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];


  /* Initialize loop */
  batch = losers[0];
  first_querypos = last_querypos = querypos = batch->querypos;
  last_diagonal = diagonal = batch->diagonal;

  floor_incr = floors_from_neg3[first_querypos];
  floor = floor_incr;
  floor_xlast = floor_incr;
  floor_xfirst = floors_from_xfirst[first_querypos] /* floors->scorefrom[xfirst_from][first_querypos] */;

#ifdef OLD_FLOOR_ENDS
  if (querypos < halfquery_lastpos) {
    floor_left = floor_incr;
  } else {
    floor_left = floors->scorefrom[-index1interval][halfquery_lastpos];
  }
  if (querypos < halfquerylength) {
    floor_right = floors->scorefrom[halfquerylength-index1interval][query_lastpos];
  } else {
    floor_right = floors->scorefrom[halfquerylength-index1interval][first_querypos];
  }
#else
  floor_left = floor_incr;
#ifdef DEBUG1
  floor_right = -99;
#endif
#endif


  debug1(printf("multiple_mm_%s, diagonal %llu, querypos %d\n",
		plusp ? "plus" : "minus",(unsigned long long) diagonal,querypos));
  debug1(printf("first_querypos = %d => initial values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
	        first_querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));

  if (--batch->npositions <= 0) {
#ifdef DIAGONAL_ADD_QUERYPOS
    batch->diagonal_add_querypos = (UINT8) -1; /* infinity */
    /* batch->diagonal_add_querypos <<= 32; */
#else
    batch->querypos = querylength; /* essentially infinity */
    batch->diagonal = (Univcoord_T) -1; /* infinity */
#endif

  } else {
    /* Use this batch for insertion (same querypos) */
#ifdef LARGE_GENOMES
    batch->diagonal = ((Univcoord_T) *(++batch->positions_high) << 32) + *(++batch->positions_low) + batch->diagterm;
#elif defined(WORDS_BIGENDIAN)
    batch->diagonal = Bigendian_convert_univcoord(*(++batch->positions)) + batch->diagterm;
#else
    batch->diagonal = *(++batch->positions) + batch->diagterm;
#endif
#ifdef DIAGONAL_ADD_QUERYPOS
    batch->diagonal_add_querypos = (UINT8) batch->diagonal;
    batch->diagonal_add_querypos <<= 32;
    batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif
  }

  /* Update tree */
#ifdef DIAGONAL_ADD_QUERYPOS
  winner_diagonal_add_querypos = batch->diagonal_add_querypos;
  for (parenti = PARENT(batch->nodei); parenti > 0; parenti = PARENT(parenti)) {
    current = losers[parenti];
    current_diagonal_add_querypos = current->diagonal_add_querypos;
    if (current_diagonal_add_querypos < winner_diagonal_add_querypos) {
      losers[parenti] = batch;
      batch = current;
      winner_diagonal_add_querypos = current_diagonal_add_querypos;
    }
  }
#else
  winner_diagonal = batch->diagonal;
  winner_querypos = batch->querypos;
  for (parenti = PARENT(batch->nodei); parenti > 0; parenti = PARENT(parenti)) {
    current = losers[parenti];
    current_diagonal = current->diagonal;
    current_querypos = current->querypos;
    if (current_diagonal < winner_diagonal ||
	(current_diagonal == winner_diagonal && current_querypos < winner_querypos)) {
      losers[parenti] = batch;
      batch = current;
      winner_diagonal = current_diagonal;
      winner_querypos = current_querypos;
    }
  }
#endif
  losers[0] = batch;


  /* Continue after initialization */
  while (
#ifdef DIAGONAL_ADD_QUERYPOS
	 (batch = losers[0])->diagonal_add_querypos != (UINT8) -1
#else
	 (batch = losers[0])->diagonal != (Univcoord_T) -1
#endif
	 ) {
   /* batch = losers[0]; */
    querypos = batch->querypos;
    diagonal = batch->diagonal;
    debug14(printf("diagonal = %u, querypos = %d\n",last_diagonal,last_querypos));

    if (diagonal == last_diagonal) {
      /* Continuing exact match or substitution */
      floor_incr = floors->scorefrom[last_querypos][querypos];
      floor += floor_incr;
      floor_xfirst += floor_incr;
      floor_xlast += floor_incr;

#ifdef OLD_FLOOR_ENDS
      /* Why is this here?  Just set floor_left at start and floor_right at end. */
      if (querypos < halfquery_lastpos) {
	floor_left += floor_incr;
      } else if (last_querypos < halfquery_lastpos) {
	/* Finish floor_left */
	floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+index1interval];
      }
      if (querypos >= halfquerylength) {
	if (last_querypos < halfquerylength) {
	  /* Start floor_right */
	  floor_right = floors->scorefrom[halfquerylength-index1interval][querypos];
	} else {
	  floor_right += floor_incr;
	}
      }
#endif

      debug1(printf("diagonal %llu unchanged: last_querypos = %d, querypos = %d => floor increments by %d\n",
		    (unsigned long long) diagonal,last_querypos,querypos,floor_incr));
      debug1(printf("*multiple_mm_%s, diagonal %llu, querypos %d, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    plusp ? "plus" : "minus",(unsigned long long) diagonal,querypos,
		    floor,floor_xfirst,floor_xlast,floor_left,floor_right));

    } else {
      /* End of diagonal */
      floor_incr = floors_to_pos3[last_querypos]  /* floors->score[last_querypos][query_lastpos+index1interval] */;
      floor += floor_incr;
      floor_xfirst += floor_incr;
      floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */

#ifdef OLD_FLOOR_ENDS
      if (last_querypos < halfquery_lastpos) {
	floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+index1interval];
	floor_right = floors->scorefrom[halfquerylength-index1interval][query_lastpos+index1interval];
      }
      if (last_querypos >= halfquerylength) {
	floor_right += floor_incr;
      }
#else
      floor_right = floor_incr;
#endif

      debug1(printf("new diagonal %llu > last diagonal %llu: last_querypos = %d => final values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    (unsigned long long) diagonal,(unsigned long long) last_diagonal,last_querypos,
		    floor,floor_xfirst,floor_xlast,floor_left,floor_right));

      if (last_diagonal > chrhigh) {
	if (ptr > ptr_chrstart) {
	  /* Add chr marker segment */
	  debug14(printf("=== ptr %p > ptr_chrstart %p, so adding chr marker segment\n",ptr,ptr_chrstart));
	  ptr->diagonal = (Univcoord_T) -1;
	  ptr_chrstart = ++ptr;
	}

	/* update chromosome bounds, based on low end */
#ifdef SLOW_CHR_UPDATE
	chrnum = Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* chrhigh += 1; */
#else
	j = 1;
#ifdef NO_EXTENSIONS_BEFORE_ZERO
	goal = last_diagonal - querylength + 1;
#else
	goal = last_diagonal + 1;
#endif
	while (j < nchromosomes_local && chrhighs_local[j] < goal) {
	  j <<= 1;			/* gallop by 2 */
	}
	if (j >= nchromosomes_local) {
	  j = binary_search(j >> 1,nchromosomes_local,chrhighs_local,goal);
	} else {
	  j = binary_search(j >> 1,j,chrhighs_local,goal);
	}
	chrnum += j;
#ifdef DEBUG15
	if (chrnum != Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength)) {
	  fprintf(stderr,"Got chrnum %d, but wanted %d\n",
		  chrnum,Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength));
	  abort();
	}
#endif
	chroffset = chroffsets[chrnum-1];
	chrhigh = chrhighs[chrnum-1];
	chrlength = chrlengths[chrnum-1];
	chrhighs_local += j;
	nchromosomes_local -= j;
#endif
      }
      if (last_diagonal <= chrhigh) { /* FORMULA for high position */
	/* position of high end is within current chromosome */
	debug1(printf("  => multiple_mm, diagonal %llu, query %d..%d, chrbounds %llu..%llu, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		      (unsigned long long) last_diagonal,first_querypos,last_querypos,
		      (unsigned long long) chroffset,(unsigned long long) chrhigh,
		      floor,floor_xfirst,floor_xlast,floor_left,floor_right));

	/* Save segment, but first advance splicesites past segment_left */
	segment_left = last_diagonal - querylength;
	max_distance = overall_max_distance;
	if (splicesites_local[0] >= last_diagonal) {
	  ptr->splicesites_i = -1;
	} else if (Splicetrie_splicesite_p(segment_left,/*pos5*/1,/*pos3*/querylength) == false) {
	  ptr->splicesites_i = -1;
	} else {
	  if (splicesites_local[0] < segment_left) {
	    j = 1;
	    while (j < nsplicesites_local && splicesites_local[j] < segment_left) {
	      j <<= 1;		/* gallop by 2 */
	    }
	    if (j >= nsplicesites_local) {
	      j = binary_search(j >> 1,nsplicesites_local,splicesites_local,segment_left);
	    } else {
	      j = binary_search(j >> 1,j,splicesites_local,segment_left);
	    }
	    joffset += j;
	    splicesites_local += j;
	    nsplicesites_local -= j;
	  }
	    
	  if (splicesites_local[0] >= last_diagonal) {
	    ptr->splicesites_i = -1;
	  } else {
	    ptr->splicesites_i = joffset;
	    j = joffset;
	    while (j < nsplicesites && splicesites[j] < last_diagonal) {
	      if (splicedists[j] > max_distance) {
		max_distance = splicedists[j];
	      }
	      j++;
	    }
	  }
	}

	/* Save segment */
	ptr->diagonal = last_diagonal;
	ptr->chrnum = chrnum;
	ptr->chroffset = chroffset;
	ptr->chrhigh = chrhigh;
	ptr->chrlength = chrlength;
	ptr->querypos5 = first_querypos;
	ptr->querypos3 = last_querypos;

	/* FORMULA */
	if (plusp) {
	  ptr->lowpos = ptr->diagonal - querylength + ptr->querypos5;
	  ptr->highpos = ptr->diagonal - querylength + ptr->querypos3 + index1part;
	} else {
	  ptr->lowpos = ptr->diagonal - ptr->querypos3 - index1part - index1part;
	  ptr->highpos = ptr->diagonal - ptr->querypos5 - index1part;
	}

	ptr->floor = floor;
	ptr->floor_xfirst = floor_xfirst;
	ptr->floor_xlast = floor_xlast;
	ptr->floor_left = floor_left;
	ptr->floor_right = floor_right;
	ptr->leftmost = ptr->rightmost = -1;
	ptr->left_splice_p = ptr->right_splice_p = false;
#if 0
	ptr->leftspan = ptr->rightspan = -1;
#endif
	ptr->usedp = false;
	ptr->pairablep = false;

#if 0
	/* Not doing this, because the max_distance test is already good enough */
	if (plusp) {
	  /* For plus-strand splicing, require segmenti->querypos3 < segmentj->querypos5,
	     so if segmenti->querypos3 is too high, then it is not spliceable */
	  if (last_querypos > query_lastpos) {
	    /* Not spliceable */
	  } else if (diagonal <= last_diagonal + max_distance) {
	    *ptr_spliceable++ = ptr;
	  }
	} else {
	  /* For minus-strand splicing, require segmenti->querypos5 > segmentj->querypos3,
	     so if segmenti->querypos5 is too low, then it is not spliceable */
	  if (first_querypos < index1part) {
	    /* Not spliceable */
	  } else if (diagonal <= last_diagonal + max_distance) {
	    *ptr_spliceable++ = ptr;
	  }
	}
#endif
	if (diagonal <= last_diagonal + max_distance) {
	  *ptr_spliceable++ = ptr;
	  debug4s(printf("%s diagonal %u is spliceable because next one is at %u\n",
			 plusp ? "plus" : "minus",last_diagonal,diagonal));
	} else {
	  debug4s(printf("%s diagonal %u is not spliceable because next one is at %u\n",
			 plusp ? "plus" : "minus",last_diagonal,diagonal));
	}
	debug14(printf("Saving segment at %u, query %d..%d",last_diagonal,ptr->querypos5,ptr->querypos3));
	*ptr_all++ = ptr;
	if (last_querypos >= first_querypos + /*min_segment_length*/1) {

	  *ptr_anchor++ = ptr;
	  debug14(printf(" ANCHOR"));
	}
	debug14(printf("\n"));
	ptr++;
      }

      /* Prepare next diagonal */
      first_querypos = querypos;
      last_diagonal = diagonal;
      floor_incr = floors_from_neg3[first_querypos] /* floors->score[-index1interval][first_querypos] */;
      floor = floor_incr;
      floor_xlast = floor_incr;
      floor_xfirst = floors_from_xfirst[first_querypos];  /* floors->score[xfirst_from][first_querypos]; */

#ifdef OLD_FLOOR_ENDS
      if (querypos < halfquery_lastpos) {
	floor_left = floor_incr;
      } else {
	floor_left = floors->scorefrom[-index1interval][halfquery_lastpos];
      }
      if (querypos < halfquerylength) {
	floor_right = floors->scorefrom[halfquerylength-index1interval][query_lastpos];
      } else {
	floor_right = floors->scorefrom[halfquerylength-index1interval][first_querypos];
      }
#else
      floor_left = floor_incr;
#ifdef DEBUG1
      floor_right = -99;	/* For debugging output */
#endif
#endif

      debug1(printf("*multiple_mm_%s, diagonal %llu, querypos %d\n",
		    plusp ? "plus" : "minus",(unsigned long long) diagonal,querypos));
      debug1(printf("start of diagonal %llu, first_querypos = %d => initial values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    (unsigned long long) diagonal,first_querypos,
		    floor,floor_xfirst,floor_xlast,floor_left,floor_right));

    }
    last_querypos = querypos;


    if (--batch->npositions <= 0) {
#ifdef DIAGONAL_ADD_QUERYPOS
      batch->diagonal_add_querypos = (UINT8) -1; /* infinity */
      /* batch->diagonal_add_querypos <<= 32; */
#else
      batch->querypos = querylength; /* essentially infinity */
      batch->diagonal = (Univcoord_T) -1; /* infinity */
#endif

    } else {
      /* Use this batch for insertion (same querypos) */
#ifdef LARGE_GENOMES
      batch->diagonal = ((Univcoord_T) *(++batch->positions_high) << 32) + *(++batch->positions_low) + batch->diagterm;
#elif defined(WORDS_BIGENDIAN)
      batch->diagonal = Bigendian_convert_univcoord(*(++batch->positions)) + batch->diagterm;
#else
      batch->diagonal = *(++batch->positions) + batch->diagterm;
#endif
#ifdef DIAGONAL_ADD_QUERYPOS
      batch->diagonal_add_querypos = (UINT8) batch->diagonal;
      batch->diagonal_add_querypos <<= 32;
      batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif
    }

    /* Update tree */
#ifdef DIAGONAL_ADD_QUERYPOS
    winner_diagonal_add_querypos = batch->diagonal_add_querypos;
    for (parenti = PARENT(batch->nodei); parenti > 0; parenti = PARENT(parenti)) {
      current = losers[parenti];
      current_diagonal_add_querypos = current->diagonal_add_querypos;
      if (current_diagonal_add_querypos < winner_diagonal_add_querypos) {
	losers[parenti] = batch;
	batch = current;
	winner_diagonal_add_querypos = current_diagonal_add_querypos;
      }
    }
#else
    winner_diagonal = batch->diagonal;
    winner_querypos = batch->querypos;
    for (parenti = PARENT(batch->nodei); parenti > 0; parenti = PARENT(parenti)) {
      current = losers[parenti];
      current_diagonal = current->diagonal;
      current_querypos = current->querypos;
      if (current_diagonal < winner_diagonal ||
	  (current_diagonal == winner_diagonal && current_querypos < winner_querypos)) {
	losers[parenti] = batch;
	batch = current;
	winner_diagonal = current_diagonal;
	winner_querypos = current_querypos;
      }
    }
#endif
  losers[0] = batch;
  }
  debug14(printf("diagonal = %u, querypos = %d\n",last_diagonal,last_querypos));
  debug14(printf("\n"));

  /* Terminate loop. */
  floor_incr = floors_to_pos3[last_querypos];   /* floors->score[last_querypos][query_lastpos+index1interval]; */
  floor += floor_incr;
  floor_xfirst += floor_incr;
  floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */

#ifdef OLD_FLOOR_ENDS
  if (last_querypos < halfquery_lastpos) {
    floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+index1interval];
    floor_right = floors->scorefrom[halfquerylength-index1interval][query_lastpos+index1interval];
  }
  if (last_querypos >= halfquerylength) {
    floor_right += floor_incr;
  }
#else
  floor_right = floor_incr;
#endif
  
  debug1(printf("no more diagonals: last_querypos = %d => terminal values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		last_querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));

  debug1(printf("last_diagonal %u vs chrhigh %u (looking for >)\n",last_diagonal,chrhigh));
  if (last_diagonal > chrhigh) {
    if (ptr > ptr_chrstart) {
      /* Add chr marker segment */
      debug14(printf("=== ptr %p > ptr_chrstart %p, so adding chr marker segment\n",ptr,ptr_chrstart));
      ptr->diagonal = (Univcoord_T) -1;
      ptr_chrstart = ++ptr;
    }

    /* update chromosome bounds, based on low end */
#ifdef SLOW_CHR_UPDATE
    chrnum = Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength);
    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
    /* chrhigh += 1; */
#else
    j = 1;
#ifdef NO_EXTENSIONS_BEFORE_ZERO
    goal = last_diagonal - querylength + 1;
#else
    goal = last_diagonal + 1;
#endif
    while (j < nchromosomes_local && chrhighs_local[j] < goal) {
      j <<= 1;			/* gallop by 2 */
    }
    if (j >= nchromosomes_local) {
      j = binary_search(j >> 1,nchromosomes_local,chrhighs_local,goal);
    } else {
      j = binary_search(j >> 1,j,chrhighs_local,goal);
    }
    chrnum += j;
#ifdef DEBUG15
    if (chrnum != Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength)) {
      fprintf(stderr,"Got chrnum %d, but wanted %d\n",
	      chrnum,Univ_IIT_get_one(chromosome_iit,last_diagonal-querylength,last_diagonal-querylength));
      abort();
    }
#endif
    chroffset = chroffsets[chrnum-1];
    chrhigh = chrhighs[chrnum-1];
    chrlength = chrlengths[chrnum-1];
    chrhighs_local += j;
    nchromosomes_local -= j;
#endif
  }

  debug1(printf("last_diagonal %u vs chrhigh %u (looking for <=)\n",last_diagonal,chrhigh));
  if (last_diagonal <= chrhigh) { /* FORMULA for high position */
    /* position of high end is within current chromosome */
    debug1(printf("  => multiple_mm, diagonal %llu, query %d..%d, chrbounds %llu..%llu, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		  (unsigned long long) last_diagonal,first_querypos,last_querypos,
		  (unsigned long long) chroffset,(unsigned long long) chrhigh,
		  floor,floor_xfirst,floor_xlast,floor_left,floor_right));

    /* Save segment, but first advance splicesites past segment_left */
    segment_left = last_diagonal - querylength;
#if 0
    /* Last segment is not spliceable */
    max_distance = overall_max_distance;
#endif
    if (splicesites_local[0] >= last_diagonal) {
      ptr->splicesites_i = -1;
    } else if (Splicetrie_splicesite_p(segment_left,/*pos5*/1,/*pos3*/querylength) == false) {
      ptr->splicesites_i = -1;
    } else {
      if (splicesites_local[0] < segment_left) {
	j = 1;
	while (j < nsplicesites_local && splicesites_local[j] < segment_left) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= nsplicesites_local) {
	  j = binary_search(j >> 1,nsplicesites_local,splicesites_local,segment_left);
	} else {
	  j = binary_search(j >> 1,j,splicesites_local,segment_left);
	}
	joffset += j;
	splicesites_local += j;
	nsplicesites_local -= j;
      }

      if (splicesites_local[0] >= last_diagonal) {
	ptr->splicesites_i = -1;
      } else {
	ptr->splicesites_i = joffset;
#if 0
	/* Last segment is not spliceable */
	if (splicedists[joffset] > overall_max_distance) {
	  max_distance = splicedists[joffset];
	}
#endif
      }
    }

    /* Save segment */
    ptr->diagonal = last_diagonal;
    ptr->chrnum = chrnum;
    ptr->chroffset = chroffset;
    ptr->chrhigh = chrhigh;
    ptr->chrlength = chrlength;
    ptr->querypos5 = first_querypos;
    ptr->querypos3 = last_querypos;

    /* FORMULA */
    if (plusp) {
      ptr->lowpos = ptr->diagonal - querylength + ptr->querypos5;
      ptr->highpos = ptr->diagonal - querylength + ptr->querypos3 + index1part;
    } else {
      ptr->lowpos = ptr->diagonal - ptr->querypos3 - index1part - index1part;
      ptr->highpos = ptr->diagonal - ptr->querypos5 - index1part;
    }

    ptr->floor = floor;
    ptr->floor_xfirst = floor_xfirst;
    ptr->floor_xlast = floor_xlast;
    ptr->floor_left = floor_left;
    ptr->floor_right = floor_right;
    ptr->leftmost = ptr->rightmost = -1;
    ptr->left_splice_p = ptr->right_splice_p = false;
#if 0
    ptr->leftspan = ptr->rightspan = -1;
#endif
    ptr->usedp = false;
    ptr->pairablep = false;

    /* Last segment is not spliceable */
    debug14(printf("Saving segment at %u, query %d..%d",last_diagonal,ptr->querypos5,ptr->querypos3));
    *ptr_all++ = ptr;
    if (last_querypos >= first_querypos + /*min_segment_length*/1) {
      debug14(printf(" ANCHOR"));
      *ptr_anchor++ = ptr;
    }
    debug14(printf("\n"));
    ptr++;
  }


  if (ptr > ptr_chrstart) {
    /* Final chr marker segment */
    debug14(printf("=== ptr %p > ptr_chrstart %p, so adding final chr marker segment\n",ptr,ptr_chrstart));
    ptr->diagonal = (Univcoord_T) -1;
    /* ptr_chrstart = */ ++ptr;
  }

#ifdef DEBUG19
  for (k = 0, ptr0 = segments; ptr0 < ptr; k++, ptr0++) {
    printf("%d %llu\n",k,(unsigned long long) ptr0->diagonal);
  }
  printf("total_npositions = %d, nchromosomes = %d\n",total_npositions,nchromosomes);
#endif

  FREEA(losers);
  FREEA(batchpool);

  /* Note: segments is in descending diagonal order.  Will need to
     reverse before solving middle deletions */

  *nsegments = ptr - segments;
  *nanchors = ptr_anchor - *anchor_segments;
  *nspliceable = ptr_spliceable - *spliceable;
  debug(printf("nsegments = %d, of which %d are spliceable (total_npositions = %d, nchromosomes = %d)\n",
	       *nsegments,*nspliceable,total_npositions,nchromosomes));
  debug1(printf("nsegments = %d, of which %d are spliceable (total_npositions = %d, nchromosomes = %d)\n",
		*nsegments,*nspliceable,total_npositions,nchromosomes));

  assert(*nsegments <= total_npositions + nchromosomes);
  assert(*nanchors <= total_npositions);
  assert(*nspliceable <= total_npositions);

  n_all_segments = ptr_all - all_segments;
  debug(printf("%d all segments\n",n_all_segments));
  debug(printf("%d anchor segments\n",*nanchors));

  if (n_all_segments <= max_anchors) {
    /* Might as well use all segments */
    FREE(*anchor_segments);
    *anchor_segments = all_segments;
    *nanchors = n_all_segments;

  } else if (*nanchors <= max_anchors) {
    /* Use only the good anchor segments */
    FREE(all_segments);

  } else {
    /* Need to limit anchor segments */
    FREE(all_segments);

    /* Treat each mod separately */
    qsort(*anchor_segments,*nanchors,sizeof(Segment_T),Segment_mod_length_cmp);

    mod = 0;
    i = 0;
    while (mod < index1interval) {
      j = i;
      while (j < *nanchors && (*anchor_segments)[j]->querypos5 % index1interval == mod) {
	j++;
      }
      nanchors_bymod[mod] = j - i;

      if (j - i <= max_anchors) {
	naccept_bymod[mod] = j - i;
      } else {
	k = i + max_anchors;
	length_threshold = (*anchor_segments)[k]->querypos3 - (*anchor_segments)[k]->querypos5;
	while (k < j && k < i + max_anchors + /*ties*/100 &&
	       (*anchor_segments)[k]->querypos3 - (*anchor_segments)[k]->querypos5 == length_threshold) {
	  k++;
	}
	naccept_bymod[mod] = k - i;
      }

      debug(printf("For mod %d, accepting %d out of %d anchor segments with length threshold %d\n",
		   mod,naccept_bymod[mod],nanchors_bymod[mod],length_threshold));
      i = j;
      mod++;
    }

    /* Move good anchors to start of array */
    dest = src = &((*anchor_segments)[0]);
    *nanchors = 0;
    for (mod = 0; mod < index1interval; mod++) {
      memmove((void *) dest,(void *) src,naccept_bymod[mod] * sizeof(Segment_T));
      dest += naccept_bymod[mod];
      src += nanchors_bymod[mod];
      *nanchors += naccept_bymod[mod];
    }

    /* Re-sort in diagonal order */
    qsort(*anchor_segments,*nanchors,sizeof(Segment_T),Segment_diagonal_cmp);
  }


#ifdef DEBUG19
  printf("%d total segments\n",*nsegments);
  for (ptr0 = segments; ptr0 < ptr; ptr0++) {
    printf("%u %d..%d\n",ptr0->diagonal,ptr0->querypos5,ptr0->querypos3);
  }
#endif

#ifdef DEBUG
  printf("%d selected anchor segments\n",*nanchors);
  for (p = &(*anchor_segments)[0]; p< &((*anchor_segments)[*nanchors]); p++) {
    segment = (Segment_T) *p;
    printf("%u %d..%d\n",segment->diagonal,segment->querypos5,segment->querypos3);
  }
#endif

  return segments;
}
#endif



#if 0
/* Modified from pair_up_concordant_aux in stage3hr.c */
static void
pair_up_segments (struct Segment_T *plus_segments_5, int plus_nsegments_5,
		  struct Segment_T *minus_segments_5, int minus_nsegments_5,
		  struct Segment_T *plus_segments_3, int plus_nsegments_3,
		  struct Segment_T *minus_segments_3, int minus_nsegments_3,
		  int querylength5, int querylength3, Chrpos_T pairmax) {
  int i, j;
  Univcoord_T insert_start;
  Segment_T segment5, segment3;	/* Need pointers, because we are changing the pairable value */

  debug(printf("Entered pair_up_segments\n"));

  /* plus/plus */
  j = 0;
  for (i = 0; i < plus_nsegments_5; i++) {
    segment5 = &(plus_segments_5[i]);
    if ((insert_start = segment5->diagonal) == (Univcoord_T) -1) {
      /* Skip chromosomal end marker */
    } else {
#ifdef DEBUG5
      printf("plus/plus: i=%d/%d %u %d..%d\n",
	     i,plus_nsegments_5,segment5->diagonal,segment5->querypos5,segment5->querypos3);
      if (j >= plus_nsegments_3) {
	printf("  current: j=%d/%d\n",j,plus_nsegments_3);
      } else if (plus_segments_3[j].diagonal == (Univcoord_T) -1) {
	printf("  current: j=%d/%d %u\n",j,plus_nsegments_3,plus_segments_3[j].diagonal);
      } else {
	printf("  current: j=%d/%d %u %d..%d\n",
	       j,plus_nsegments_3,plus_segments_3[j].diagonal,plus_segments_3[j].querypos5,plus_segments_3[j].querypos3);
      }
#endif

      /* Get to correct chrnum */
      while (j < plus_nsegments_3 && (plus_segments_3[j].diagonal == (Univcoord_T) -1 || plus_segments_3[j].diagonal < segment5->diagonal)) {
#ifdef DEBUG5
	if (plus_segments_3[j].diagonal == (Univcoord_T) -1) {
	  printf("  advancing: j=%d/%d %u\n",j,plus_nsegments_3,plus_segments_3[j].diagonal);
	} else {
	  printf("  advancing: j=%d/%d %u %d..%d\n",
		 j,plus_nsegments_3,plus_segments_3[j].diagonal,plus_segments_3[j].querypos5,plus_segments_3[j].querypos3);
	}
#endif
	j++;
      }

      if (j < plus_nsegments_3) {
	while (j >= 0 && plus_segments_3[j].diagonal != (Univcoord_T) -1 && plus_segments_3[j].diagonal > segment5->diagonal) {
	  debug5(printf("  backup: j=%d/%d %u %d..%d\n",
			j,plus_nsegments_3,plus_segments_3[j].diagonal,plus_segments_3[j].querypos5,plus_segments_3[j].querypos3));
	  j--;
	}
	j++;		/* Finish backup */

	/* Cannot perform arithmetic on diagonal, because we want to preserve -1 as being the largest value */
	/* Ignore inclusion of querylength inside pairmax */
	while (j < plus_nsegments_3 && plus_segments_3[j].diagonal <= insert_start + pairmax /*- querylength3*/) {
	  debug5(printf("  overlap: j=%d/%d, %u <= %u + %u, %d..%d\n",
			j,plus_nsegments_3,plus_segments_3[j].diagonal,
			insert_start,pairmax,plus_segments_3[j].querypos5,plus_segments_3[j].querypos3));
	  debug5(printf("Setting plus segments %d and %d to be pairable: %u and %u\n",i,j,segment5->diagonal,plus_segments_3[j].diagonal));
	  segment5->pairablep = true;
	  plus_segments_3[j].pairablep = true;
	  j++;
	}
      }
    }
  }
		
  /* minus/minus */
  j = 0;
  for (i = 0; i < minus_nsegments_3; i++) {
    segment3 = &(minus_segments_3[i]);
    if ((insert_start = segment3->diagonal) == (Univcoord_T) -1) {
      /* Skip chromosomal end marker */
    } else {
#ifdef DEBUG5
      printf("minus/minus: i=%d/%d %u %d..%d\n",
	     i,minus_nsegments_3,segment3->diagonal,segment3->querypos5,segment3->querypos3);
      if (j >= minus_nsegments_5) {
	printf("  current: j=%d/%d\n",j,minus_nsegments_5);
      } else if (minus_segments_5[j].diagonal == (Univcoord_T) -1) {
	printf("  current: j=%d/%d %u\n",j,minus_nsegments_5,minus_segments_5[j].diagonal);
      } else {
	printf("  current: j=%d/%d %u %d..%d\n",
	       j,minus_nsegments_5,minus_segments_5[j].diagonal,minus_segments_5[j].querypos5,minus_segments_5[j].querypos3);
      }
#endif
      
      /* Get to correct chrnum */
      while (j < minus_nsegments_5 && (minus_segments_5[j].diagonal == (Univcoord_T) -1 || minus_segments_5[j].diagonal < segment3->diagonal)) {
#ifdef DEBUG5
	if (minus_segments_5[j].diagonal == (Univcoord_T) -1) {
	  printf("  advancing: j=%d/%d %u\n",j,minus_nsegments_5,minus_segments_5[j].diagonal);
	} else {
	  printf("  advancing: j=%d/%d %u %d..%d\n",
		 j,minus_nsegments_5,minus_segments_5[j].diagonal,minus_segments_5[j].querypos5,minus_segments_5[j].querypos3);
	}
#endif
	j++;
      }

      if (j < minus_nsegments_5) {
	while (j >= 0 && minus_segments_5[j].diagonal != (Univcoord_T) -1 && minus_segments_5[j].diagonal > segment3->diagonal) {
	  debug5(printf("  backup: j=%d/%d %u %d..%d\n",
			j,minus_nsegments_5,minus_segments_5[j].diagonal,minus_segments_5[j].querypos5,minus_segments_5[j].querypos3));
	  j--;
	}
	j++;		/* Finish backup */

	/* Cannot perform arithmetic on diagonal, because we want to preserve -1 as being the largest value */
	/* Ignore inclusion of querylength inside pairmax */
	while (j < minus_nsegments_5 && minus_segments_5[j].diagonal <= insert_start + pairmax /*- querylength5*/) {
	  debug5(printf("  overlap: j=%d/%d %u %d..%d\n",
			j,minus_nsegments_5,minus_segments_5[j].diagonal,minus_segments_5[j].querypos5,minus_segments_5[j].querypos3));
	  debug5(printf("Setting minus segments %d and %d to be pairable: %u and %u\n",i,j,segment3->diagonal,minus_segments_5[j].diagonal));
	  segment3->pairablep = true;
	  minus_segments_5[j].pairablep = true;
	  j++;
	}
      }
    }
  }

  return;
}
#endif


static void
pair_up_anchor_segments (Segment_T *plus_anchor_segments_5, Segment_T *minus_anchor_segments_5,
			 Segment_T *plus_anchor_segments_3, Segment_T *minus_anchor_segments_3,
			 int n_plus_anchors_5, int n_minus_anchors_5,
			 int n_plus_anchors_3, int n_minus_anchors_3, Chrpos_T pairmax) {
  /* Univcoord_T insert_start; */
  Segment_T segment5, segment3;
  Segment_T *q, *pstart, *pend, *p;

  debug(printf("Entering pair_up_anchor_segments\n"));

  /* plus/plus */
  pstart = &(plus_anchor_segments_3[0]);
  for (q = &(plus_anchor_segments_5[0]);
       q < &(plus_anchor_segments_5[n_plus_anchors_5]) && pstart < &(plus_anchor_segments_3[n_plus_anchors_3]);
       q++) {
    segment5 = *q;
    assert(segment5->diagonal != (Univcoord_T) -1);
    /* insert_start = segment5->diagonal; */

    while (pstart < &(plus_anchor_segments_3[n_plus_anchors_3]) && (*pstart)->diagonal < segment5->diagonal) {
      pstart++;
    }

    pend = pstart;
    while (pend < &(plus_anchor_segments_3[n_plus_anchors_3]) && (*pend)->diagonal < segment5->diagonal + pairmax) {
      pend++;
    }
	
    for (p = pstart; p != pend; p++) {
      segment3 = *p;
      assert(segment3->diagonal - segment5->diagonal < pairmax);
      debug5(printf("Setting plus segments to be pairable: %u and %u (distance %u)\n",
		    segment5->diagonal,segment3->diagonal,segment3->diagonal - segment5->diagonal));
      segment5->pairablep = true;
      segment3->pairablep = true;
    }
  }
		
  /* minus/minus */
  pstart = &(minus_anchor_segments_5[0]);
  for (q = &(minus_anchor_segments_3[0]);
       q < &(minus_anchor_segments_3[n_minus_anchors_3]) && pstart < &(minus_anchor_segments_5[n_minus_anchors_5]);
       q++) {
    segment3 = *q;
    assert(segment3->diagonal != (Univcoord_T) -1);
    /* insert_start = segment3->diagonal; */

    while (pstart < &(minus_anchor_segments_5[n_minus_anchors_5]) && (*pstart)->diagonal < segment3->diagonal) {
      pstart++;
    }

    pend = pstart;
    while (pend < &(minus_anchor_segments_5[n_minus_anchors_5]) && (*pend)->diagonal < segment3->diagonal + pairmax) {
      pend++;
    }

    for (p = pstart; p != pend; p++) {
      segment5 = *p;
      assert(segment5->diagonal - segment3->diagonal < pairmax);
      debug5(printf("Setting minus segments to be pairable: %u and %u (distance %u)\n",
		    segment3->diagonal,segment5->diagonal,segment5->diagonal - segment3->diagonal));
      segment3->pairablep = true;
      segment5->pairablep = true;
    }
  }

  debug(printf("Exiting pair_up_anchor_segments\n"));

  return;
}



/*

The pattern below is a middle insertion on plus strand, or middle deletion on minus strand:

diagonal 2354, querypos 18
diagonal 2354, querypos 19
diagonal 2354, querypos 20
diagonal 2354, querypos 21
diagonal 2354, querypos 22
diagonal 2354, querypos 23
diagonal 2354, querypos 24
diagonal 2356, querypos 0
diagonal 2356, querypos 1
diagonal 2356, querypos 2
diagonal 2356, querypos 3
diagonal 2356, querypos 4
diagonal 2356, querypos 5


The pattern below is a middle deletion on plus strand, or middle insertion on minus strand:

diagonal 2354, querypos 0
diagonal 2354, querypos 1
diagonal 2354, querypos 2
diagonal 2354, querypos 3
diagonal 2354, querypos 4
diagonal 2354, querypos 5
diagonal 2354, querypos 6
diagonal 2354, querypos 7
diagonal 2356, querypos 18
diagonal 2356, querypos 19
diagonal 2356, querypos 20
diagonal 2356, querypos 21
diagonal 2356, querypos 22
diagonal 2356, querypos 23
diagonal 2356, querypos 24

*/


static List_T
find_middle_indels (int *found_score, int *nhits, List_T hits,
		    Segment_T *plus_spliceable, int plus_nspliceable,
		    Segment_T *minus_spliceable, int minus_nspliceable,
#ifdef DEBUG2
		    char *queryuc_ptr, char *queryrc, 
#endif
		    Floors_T floors, int querylength, int query_lastpos,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    int max_mismatches_allowed, int genestrand) {
  Segment_T segmenti, segmentj, segmentj_end, *ptr;
  int indels, floor, pos, prev, middle;
  int *floors_from_neg3, *floors_to_pos3;
  bool foundp;
  int max_middle_insertions;

  if (max_middle_insertions_default >= 0) {
    max_middle_insertions = max_middle_insertions_default;
  } else {
    max_middle_insertions = querylength - 2*min_indel_end_matches;
  }

  debug(printf("*** find_middle_indels with querylength %d and max_mismatches_allowed %d ***\n",
	       querylength,max_mismatches_allowed));

  /* Plus segments */
  if (floors != NULL) {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    debug2(printf("plus_nspliceable = %d\n",plus_nspliceable));
    for (ptr = plus_spliceable; ptr < &(plus_spliceable[plus_nspliceable]); ptr++) {
      segmenti = *ptr;
      debug2(printf("\nplus segmenti:  diagonal %llu, querypos %d..%d\n",
		    (unsigned long long) segmenti->diagonal,segmenti->querypos5,segmenti->querypos3));
      if (1 || segmenti->diagonal < (Univcoord_T) -1) { /* No markers were stored in spliceable */
	/* Identify potential segmentj for segmenti */
	segmentj_end = segmenti+1;
	while (
#ifdef NO_MARKER_SEGMENTS
	       segmentj_end < &(plus_segments[plus_nsegments]) && segmentj_end->chrnum == segmenti->chrnum &&
#endif
	       segmentj_end->diagonal <= segmenti->diagonal + max_middle_insertions) {
	  segmentj_end++;
	}
	  
	for (segmentj = segmenti+1; segmentj < segmentj_end; segmentj++) {
	  debug2(printf("plus insertion?  diagonal %llu, querypos %d..%d => diagonal %llu, querypos %d..%d => ",
			(unsigned long long) segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
			(unsigned long long) segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	  /* j5 j3 i5 i3 */
	  if (segmentj->querypos3 < segmenti->querypos5) {
	    indels = segmentj->diagonal - segmenti->diagonal; /* positive */
	    floor = floors_from_neg3[segmentj->querypos5] + floors_to_pos3[segmenti->querypos3]
	      /* floors->score[-index1interval][segmentj->querypos5] + floors->score[segmenti->querypos3][query_lastpos+index1interval] */ ;
	    if (floors->prev_omitted == NULL) {
	      if ((middle = FLOOR_MIDDLE(segmenti->querypos5 - segmentj->querypos3 - indels)) > 0) {
		middle--;	/* for insertion, which looks like a mismatch */
	      }
	      debug2(printf("\nmiddle (no omission): %d\n",middle));
	      floor += middle;
	    } else {
	      pos = segmenti->querypos5;
	      debug2(printf("\nmiddle (omission):"));
	      while (pos > segmentj->querypos3) {
		if ((prev = floors->prev_omitted[pos]) < segmentj->querypos3) {
		  prev = segmentj->querypos3;
		}
		if ((middle = FLOOR_MIDDLE(pos - prev - indels)) > 0) {
		  middle--;	/* for insertion, which looks like a mismatch */
		}
		floor += middle;
		debug2(printf("(%d..%d)+%d,",prev,pos,middle));
		pos = prev;
	      }
	      debug2(printf("\n"));
	    }
	    if (floor <= max_mismatches_allowed) {
	      debug2(printf("successful insertion, floor = %d+middle+%d=%d, indels = %d\n",
			    floors->scorefrom[-index1interval][segmentj->querypos5],
			    floors->scorefrom[segmenti->querypos3][query_lastpos+index1interval],
			    floor,indels));
	      hits = Indel_solve_middle_insertion(&foundp,&(*found_score),&(*nhits),hits,
						  /*left*/segmenti->diagonal - querylength,
						  segmenti->chrnum,segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
						  indels,/*query_compress*/query_compress_fwd,
						  querylength,max_mismatches_allowed,
						  /*plusp*/true,genestrand,/*sarrayp*/false);
	    } else {
	      debug2(printf("too many mismatches, because floor %d+middle+%d=%d > %d\n",
			    floors->scorefrom[-index1interval][segmentj->querypos5],
			    floors->scorefrom[segmenti->querypos3][query_lastpos+index1interval],
			    floor,max_mismatches_allowed));
	    }
	  } else {
	    debug2(printf("garbage, because querypos3 %d >= querypos5 %d\n",
			  segmentj->querypos3,segmenti->querypos5));
	  }
	}

	/* Identify potential segmentj for segmenti */
	segmentj_end = segmenti+1;
	while (
#ifdef NO_MARKER_SEGMENTS
	       segmentj_end < &(plus_segments[plus_nsegments]) && segmentj_end->chrnum == segmenti->chrnum &&
#endif
	       segmentj_end->diagonal <= segmenti->diagonal + max_middle_deletions) {
	  segmentj_end++;
	}

	for (segmentj = segmenti+1; segmentj < segmentj_end; segmentj++) {
	  debug2(printf("plus deletion?  diagonal %llu, querypos %d..%d => diagonal %llu, querypos %d..%d => ",
			(unsigned long long) segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
			(unsigned long long) segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	  /* i5 i3 j5 j3 */
	  if (segmenti->querypos3 < segmentj->querypos5) {
	    indels = segmenti->diagonal - segmentj->diagonal; /* negative */
	    floor = floors_from_neg3[segmenti->querypos5] + floors_to_pos3[segmentj->querypos3]
	      /* floors->score[-index1interval][segmenti->querypos5] + floors->score[segmentj->querypos3][query_lastpos+index1interval] */;
	    if (floors->prev_omitted == NULL) {
	      if ((middle = FLOOR_MIDDLE(segmentj->querypos5 - segmenti->querypos3 /*- indels*/)) > 0) {
		middle--;	/* for deletion, which looks like a mismatch */
	      }
	      debug2(printf("\nmiddle (no omission): %d\n",middle));
	      floor += middle;
	    } else {
	      pos = segmentj->querypos5;
	      debug2(printf("\nmiddle (omission):"));
	      while (pos > segmenti->querypos3) {
		if ((prev = floors->prev_omitted[pos]) < segmenti->querypos3) {
		  prev = segmenti->querypos3;
		}
		if ((middle = FLOOR_MIDDLE(pos - prev /*- indels*/)) > 0) {
		  middle--;	/* for deletion, which looks like a mismatch */
		}
		floor += middle;
		debug2(printf("(%d..%d)+%d,",prev,pos,middle));
		pos = prev;
	      }
	      debug2(printf("\n"));
	    }
	    if (floor <= max_mismatches_allowed) {
	      debug2(printf("successful deletion, floor = %d+middle+%d=%d, indels = %d\n",
			    floors->scorefrom[-index1interval][segmenti->querypos5],
			    floors->scorefrom[segmentj->querypos3][query_lastpos+index1interval],
			    floor,indels));
	      hits = Indel_solve_middle_deletion(&foundp,&(*found_score),&(*nhits),hits,
						 /*left*/segmenti->diagonal - querylength,
						 segmenti->chrnum,segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
						 indels,/*query_compress*/query_compress_fwd,querylength,
						 max_mismatches_allowed,
						 /*plusp*/true,genestrand,/*sarrayp*/false);
	    } else {
	      debug2(printf("too many mismatches, because floor = %d+middle+%d=%d > %d\n",
			    floors->scorefrom[-index1interval][segmenti->querypos5],
			    floors->scorefrom[segmentj->querypos3][query_lastpos+index1interval],
			    floor,max_mismatches_allowed));
	    }
	  } else {
	    debug2(printf("garbage, because querypos3 %d >= querypos5 %d\n",
			  segmenti->querypos3,segmentj->querypos5));
	  }
	}
      }
    }

    /* Minus segments */
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    debug2(printf("minus_nspliceable = %d\n",minus_nspliceable));
    for (ptr = minus_spliceable; ptr < &(minus_spliceable[minus_nspliceable]); ptr++) {
      segmenti = *ptr;
      debug2(printf("\nminus segmenti:  diagonal %llu, querypos %d..%d\n",
		    (unsigned long long) segmenti->diagonal,segmenti->querypos5,segmenti->querypos3));
      if (1 || segmenti->diagonal < (Univcoord_T) -1) { /* No markers were stored in spliceable */
	/* Identify potential segmentj for segmenti */
	segmentj_end = segmenti+1;
	while (
#ifdef NO_MARKER_SEGMENTS
	       segmentj_end < &(minus_segments[minus_nsegments]) && segmentj_end->chrnum == segmenti->chrnum &&
#endif
	       segmentj_end->diagonal <= segmenti->diagonal + max_middle_deletions) {
	  segmentj_end++;
	}

	for (segmentj = segmenti+1; segmentj < segmentj_end; segmentj++) {
	  debug2(printf("minus deletion?  diagonal %llu, querypos %d..%d => diagonal %llu, querypos %d..%d => ",
			(unsigned long long) segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
			(unsigned long long) segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	  /* j5 j3 i5 i3 */
	  if (segmentj->querypos3 < segmenti->querypos5) {
	    indels = segmenti->diagonal - segmentj->diagonal; /* negative */
	    floor = floors_from_neg3[segmentj->querypos5] + floors_to_pos3[segmenti->querypos3]
	      /* floors->score[-index1interval][segmentj->querypos5] + floors->score[segmenti->querypos3][query_lastpos+index1interval] */;
	    if (floors->prev_omitted == NULL) {
	      if ((middle = FLOOR_MIDDLE(segmenti->querypos5 - segmentj->querypos3 /*- indels*/)) > 0) {
		middle--;	/* for deletion, which looks like a mismatch */
	      }
	      debug2(printf("\nmiddle (no omission): %d\n",middle));
	      floor += middle;
	    } else {
	      pos = segmenti->querypos5;
	      debug2(printf("\nmiddle (omission):"));
	      while (pos > segmentj->querypos3) {
		if ((prev = floors->prev_omitted[pos]) < segmentj->querypos3) {
		  prev = segmentj->querypos3;
		}
		if ((middle = FLOOR_MIDDLE(pos - prev /*- indels*/)) > 0) {
		  middle--; /* for deletion, which looks like a mismatch */
		}
		floor += middle;
		debug2(printf("(%d..%d)+%d,",prev,pos,middle));
		pos = prev;
	      }
	      debug2(printf("\n"));
	    }
	    if (floor <= max_mismatches_allowed) {
	      debug2(printf("successful deletion, floor = %d+middle+%d=%d, indels = %d\n",
			    floors->scorefrom[-index1interval][segmentj->querypos5],
			    floors->scorefrom[segmenti->querypos3][query_lastpos+index1interval],
			    floor,indels));
	      hits = Indel_solve_middle_deletion(&foundp,&(*found_score),&(*nhits),hits,
						 /*left*/segmenti->diagonal - querylength,
						 segmenti->chrnum,segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
						 indels,/*query_compress*/query_compress_rev,querylength,
						 max_mismatches_allowed,
						 /*plusp*/false,genestrand,/*sarrayp*/false);
	    } else {
	      debug2(printf("too many mismatches, because floor = %d+middle+%d=%d > %d\n",
			    floors->scorefrom[-index1interval][segmentj->querypos5],
			    floors->scorefrom[segmenti->querypos3][query_lastpos+index1interval],
			    floor,max_mismatches_allowed));
	      debug2(printf("too many mismatches, because floor %d > %d\n",floor,max_mismatches_allowed));
	    }
	  } else {
	    debug2(printf("garbage, because querypos3 %d >= querypos5 %d\n",
			  segmentj->querypos3,segmenti->querypos5));
	  }
	}

	/* Identify potential segmentj for segmenti */
	segmentj_end = segmenti+1;
	while (
#ifdef NO_MARKER_SEGMENTS
	       segmentj_end < &(minus_segments[minus_nsegments]) && segmentj_end->chrnum == segmenti->chrnum &&
#endif
	       segmentj_end->diagonal <= segmenti->diagonal + max_middle_insertions) {
	  segmentj_end++;
	}

	for (segmentj = segmenti+1; segmentj < segmentj_end; segmentj++) {
	  debug2(printf("minus insertion?  diagonal %llu, querypos %d..%d => diagonal %llu, querypos %d..%d => ",
			(unsigned long long) segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
			(unsigned long long) segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	  /* i5 i3 j5 j3 */
	  if (segmenti->querypos3 < segmentj->querypos5) {
	    indels = segmentj->diagonal - segmenti->diagonal; /* positive */
	    floor = floors_from_neg3[segmenti->querypos5] + floors_to_pos3[segmentj->querypos3]
	      /* floors->score[-index1interval][segmenti->querypos5] + floors->score[segmentj->querypos3][query_lastpos+index1interval] */;
	    if (floors->prev_omitted == NULL) {
	      if ((middle = FLOOR_MIDDLE(segmentj->querypos5 - segmenti->querypos3 - indels)) > 0) {
		middle--;	/* for insertion, which looks like a mismatch */
	      }
	      debug2(printf("\nmiddle (no omission): %d\n",middle));
	      floor += middle;
	    } else {
	      pos = segmentj->querypos5;
	      debug2(printf("\nmiddle (omission):"));
	      while (pos > segmenti->querypos3) {
		if ((prev = floors->prev_omitted[pos]) < segmenti->querypos3) {
		  prev = segmenti->querypos3;
		}
		if ((middle = FLOOR_MIDDLE(pos - prev - indels)) > 0) {
		  middle--;	/* for insertion, which looks like a mismatch */
		}
		floor += middle;
		debug2(printf("(%d..%d)+%d,",prev,pos,middle));
		pos = prev;
	      }
	      debug2(printf("\n"));
	    }
	    if (floor <= max_mismatches_allowed) {
	      debug2(printf("successful insertion, floor = %d+middle+%d=%d, indels = %d\n",
			    floors->scorefrom[-index1interval][segmenti->querypos5],
			    floors->scorefrom[segmentj->querypos3][query_lastpos+index1interval],
			    floor,indels));
	      hits = Indel_solve_middle_insertion(&foundp,&(*found_score),&(*nhits),hits,
						  /*left*/segmenti->diagonal - querylength,
						  segmenti->chrnum,segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
						  indels,/*query_compress*/query_compress_rev,
						  querylength,max_mismatches_allowed,
						  /*plusp*/false,genestrand,/*sarrayp*/false);
	    } else {
	      debug2(printf("too many mismatches, because floor %d+middle+%d=%d > %d\n",
			    floors->scorefrom[-index1interval][segmenti->querypos5],
			    floors->scorefrom[segmentj->querypos3][query_lastpos+index1interval],
			    floor,max_mismatches_allowed));
	    }
	  } else {
	    debug2(printf("garbage, because querypos3 %d >= querypos5 %d\n",
			  segmenti->querypos3,segmentj->querypos5));
	  }
	}
      }
    }
  }

  return hits;
}


/************************************************************************/

/************************************************************************
 *   right deletion: use <  / indel_pos = [conti]
 *   right insertion: use - sep <  / indel_pos = [conti]
 *   left deletion: use >  / indel_pos = [conti] + 1
 *   left insertion: use + sep >   / indel_pos = [conti]-sep+1
 ************************************************************************/

static int
compute_end_indels_right (int *indels, int *nmismatches_longcont, int *nmismatches_shift,
			  int *mismatch_positions_long, int nmismatches_avail_long,
			  int breakpoint, int querylength, Univcoord_T left, Compress_T query_compress,
			  int min_indel_end_matches, int max_end_insertions, int max_end_deletions,
			  int max_mismatches_short, bool plusp, int genestrand) {
#ifdef DEBUG2E
  int i;
#endif
  int length1;
  int sep, end;
  int nmismatches_avail_shift, nmatches;
  int sum, best_sum = querylength;
  int conti, shifti;
  int best_indel_pos = -1, endlength;
#ifdef OLD_END_INDELS
  int indel_pos;
#else
  int indel_pos_cont, indel_pos_shift;
#endif
  int *mismatch_positions_shift;


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions_shift = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions_shift = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions_shift = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


  debug2e(printf("Entered compute_end_indels_right with breakpoint = %d, max_mismatches_short %d\n",
		 breakpoint,max_mismatches_short));
  length1 = querylength - breakpoint;
#if 0
  /* Should not need to reset max_end_deletions */
  if (max_end_deletions > length1 - min_indel_end_matches) {
    max_end_deletions = length1 - min_indel_end_matches;
  }
#endif
  if (max_end_insertions > length1 - min_indel_end_matches) {
    max_end_insertions = length1 - min_indel_end_matches;
  }

  if (max_end_deletions > 0) {
    for (sep = 1; sep <= max_end_deletions; sep++) {
      /* *indels = -sep; */
      nmismatches_avail_shift = Genome_mismatches_right(mismatch_positions_shift,
							max_mismatches_short,query_compress,
							left-(-sep),/*pos5*/0,/*pos3*/querylength,
							plusp,genestrand);
      debug2e(
	      printf("A. Trying deletion of %d.  ",-sep);
	      printf("%d mismatches on right at:",nmismatches_avail_shift);
	      for (i = 0; i <= nmismatches_avail_shift; i++) {
		printf(" %d",mismatch_positions_shift[i]);
	      }
	      printf(".  ");
	      );

      if (nmismatches_avail_shift == 0) {
	/* Skip, because we have an exact match on the shift */
      } else {
#ifdef OLD_END_INDELS
	/* Compute over mismatch_positions_shift[n-1] to querylength */
	/* A. Right deletion: Primary loop along Genome_mismatches_right (shifti) to get lowest coordinate */
	shifti = 0;
	conti = nmismatches_avail_long - 1;
	while (shifti < nmismatches_avail_shift) {
	  while (conti >= 0 && mismatch_positions_long[conti] > mismatch_positions_shift[shifti]) {
	    conti--;
	  }
	  sum = shifti + conti + 1;
	  debug2e(printf("sum %d=%d+%d at indel_pos %d.  ",sum,conti+1,shifti,mismatch_positions_shift[shifti]+1));
	  if (sum < best_sum) {
	    indel_pos = mismatch_positions_shift[shifti] + 1;
	    if ((endlength = querylength - indel_pos) >= min_indel_end_matches && endlength >= sep) {
	      /* Don't want to delete more than the amount matched */
	      nmatches = endlength - shifti;
	      if (nmatches - 3*shifti - 4 >= 0) {
		debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		/* Want more matches than mismatches */
		best_indel_pos = indel_pos;
		*indels = -sep;
		*nmismatches_longcont = conti + 1;
		*nmismatches_shift = shifti;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	  shifti++;
	}
	debug2e(printf("\n"));

	/* A. Right deletion: Primary loop along Genome_mismatches_left (conti) to see if we missed anything */
	shifti = nmismatches_avail_shift - 1;
	conti = 0;

	/* Start in region */
	while (conti < nmismatches_avail_long && mismatch_positions_long[conti] < mismatch_positions_shift[shifti]) {
	  conti++;
	}

	while (conti < nmismatches_avail_long) {
	  while (shifti >= 0 && mismatch_positions_shift[shifti] < mismatch_positions_long[conti]) {
	    shifti--;
	  }
	  sum = conti + shifti + 1;
	  debug2e(printf("sum %d=%d+%d at indel_pos %d.  ",sum,conti,shifti+1,mismatch_positions_long[conti]));
	  if (sum < best_sum) {
	    indel_pos = mismatch_positions_long[conti];
	    if ((endlength = querylength - indel_pos) >= min_indel_end_matches && endlength >= sep) {
	      /* Don't want to delete more than the amount matched */
	      nmatches = endlength - (shifti + 1);
	      if (nmatches - 3*(shifti+1) - 4 >= 0) {
		/* Want more matches than mismatches */
		debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti+1,nmatches-3*(shifti+1)-4));
		best_indel_pos = indel_pos;
		*indels = -sep;
		*nmismatches_longcont = conti;
		*nmismatches_shift = shifti + 1;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	  conti++;
	}
	debug2e(printf("\n"));
#else
	shifti = nmismatches_avail_shift - 1;
	conti = 0;
	while (conti < nmismatches_avail_long && mismatch_positions_long[conti] < mismatch_positions_shift[shifti]) {
	  conti++;
	}
	indel_pos_cont = mismatch_positions_long[conti];
	indel_pos_shift = mismatch_positions_shift[shifti] + 1;

	while (conti < nmismatches_avail_long && shifti >= 0) {
	  if (indel_pos_cont < indel_pos_shift) {
	    sum = shifti + conti + 1;
	    debug2e(printf("cont %d=%d+%d at indel_pos %d.  ",sum,conti,shifti+1,mismatch_positions_long[conti]));
	    if (sum < best_sum) {
	      if ((endlength = querylength - indel_pos_cont) >= min_indel_end_matches && endlength >= sep) {
		/* Don't want to delete more than the amount matched */
		nmatches = endlength - (shifti + 1);
		if (nmatches - 3*(shifti+1) - 4 >= 0) {
		  /* Want more matches than mismatches */
		  /* Values -3 and -4 correspond to defaults for trim_mismatch_score and trim_indel_score */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti+1,nmatches-3*(shifti+1)-4));
		  best_indel_pos = indel_pos_cont;
		  *indels = -sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti + 1;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    conti++;
	    indel_pos_cont = mismatch_positions_long[conti];

	  } else if (indel_pos_shift < indel_pos_cont) {
	    sum = shifti + conti;
	    debug2e(printf("shift %d=%d+%d at indel_pos %d.  ",sum,conti,shifti,mismatch_positions_shift[shifti]+1));
	    if (sum < best_sum) {
	      if ((endlength = querylength - indel_pos_shift) >= min_indel_end_matches && endlength >= sep) {
		/* Don't want to delete more than the amount matched */
		nmatches = endlength - shifti;
		if (nmatches - 3*shifti - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		  best_indel_pos = indel_pos_shift;
		  *indels = -sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    if (--shifti >= 0) {
	      indel_pos_shift = mismatch_positions_shift[shifti] + 1;
	    }

	  } else {
	    sum = shifti + conti;
	    debug2e(printf("both %d=%d+%d at indel_pos %d.  ",sum,conti,shifti,indel_pos_cont));
	    if (sum < best_sum) {
	      if ((endlength = querylength - indel_pos_shift) >= min_indel_end_matches && endlength >= sep) {
		/* Don't want to delete more than the amount matched */
		nmatches = endlength - shifti;
		if (nmatches - 3*shifti - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		  best_indel_pos = indel_pos_shift;
		  *indels = -sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    conti++;
	    if (--shifti >= 0) {
	      indel_pos_cont = mismatch_positions_long[conti];
	      indel_pos_shift = mismatch_positions_shift[shifti] + 1;
	    }
	  }
	}
	
	if (shifti < 0) {
	  sum = conti /*+ shifti + 1*/;
	  debug2e(printf("last %d=%d at indel_pos %d.  ",sum,conti,indel_pos_cont));
	  if (sum < best_sum) {
	    if ((endlength = querylength - indel_pos_cont) >= min_indel_end_matches && endlength >= sep) {
	      /* Don't want to delete more than the amount matched */
	      nmatches = endlength /*- (shifti + 1)*/;
	      if (nmatches >= /*shifti + 1*/ + 4) {
		/* Want more matches than mismatches */
		best_indel_pos = indel_pos_cont;
		*indels = -sep;
		*nmismatches_longcont = conti;
		*nmismatches_shift = 0 /*shifti + 1*/;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	}

	debug2e(printf("\n"));
#endif

      }
    }
  }
  
  if (max_end_insertions > 0) {
    if (left < (unsigned int) max_end_insertions) {
      debug2e(printf("left %llu < max_end_insertions %d, so end = left\n",
		     (unsigned long long) left,max_end_insertions));
      end = left;
    } else {
      end = max_end_insertions;
    }

    for (sep = 1; sep <= end; sep++) {
      /* *indels = +sep; */
      nmismatches_avail_shift = Genome_mismatches_right(mismatch_positions_shift,
							max_mismatches_short,query_compress,
							left-(+sep),/*pos5*/0,/*pos3*/querylength,
							plusp,genestrand);

      debug2e(
	      printf("B. Trying insertion of %d.  ",+sep);
	      printf("%d mismatches on right at:",nmismatches_avail_shift);
	      for (i = 0; i <= nmismatches_avail_shift; i++) {
		printf(" %d",mismatch_positions_shift[i]);
	      }
	      printf(".  ");
	      );

      if (nmismatches_avail_shift == 0) {
	/* Skip, because we have an exact match on the shift */
      } else {
#ifdef OLD_END_INDELS
	/* Compute over mismatch_positions_shift[n-1] to querylength */
	/* B. Right insertion: First, try primary loop along Genome_mismatches_right (shifti) to get lowest coordinate */
	shifti = 0;
	conti = nmismatches_avail_long - 1;
	while (shifti < nmismatches_avail_shift) {
	  while (conti >= 0 && mismatch_positions_long[conti] > mismatch_positions_shift[shifti] - sep) {
	    conti--;
	  }
	  sum = shifti + conti + 1;
	  debug2e(printf("sum %d=%d+%d at indel_pos %d.  ",sum,conti+1,shifti,mismatch_positions_shift[shifti]-sep+1));
	  if (sum < best_sum) {
	    indel_pos = mismatch_positions_shift[shifti] - sep + 1;
	    endlength = querylength - (indel_pos + sep);
	    if (endlength >= min_indel_end_matches) {
	      nmatches = endlength - shifti;
	      if (nmatches - 3*shifti - 4 >= 0) {
		/* Want more matches than mismatches */
		debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		best_indel_pos = indel_pos;
		*indels = +sep;
		*nmismatches_longcont = conti + 1;
		*nmismatches_shift = shifti;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	  shifti++;
	}
	debug2e(printf("\n"));


	/* B. Right insertion: Try primary loop along Genome_mismatches_left (conti) to see if we missed anything */
	shifti = nmismatches_avail_shift - 1;
	conti = 0;

	/* Start in region */
	while (conti < nmismatches_avail_long && mismatch_positions_long[conti] < mismatch_positions_shift[shifti]) {
	  conti++;
	}

	while (conti < nmismatches_avail_long) {
	  while (shifti >= 0 && mismatch_positions_shift[shifti] < mismatch_positions_long[conti] + sep) {
	    shifti--;
	  }
	  sum = conti + shifti + 1;
	  debug2e(printf("sum %d=%d+%d at indel_pos %d.  ",sum,conti,shifti+1,mismatch_positions_long[conti]));
	  if (sum < best_sum) {
	    indel_pos = mismatch_positions_long[conti];
	    endlength = querylength - indel_pos - sep;
	    if (endlength >= min_indel_end_matches) {
	      nmatches = endlength - (shifti + 1);
	      if (nmatches - 3*(shifti+1) - 4 >= 0) {
		/* Want more matches than mismatches */
		debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti+1,nmatches-3*(shifti+1)-4));
		best_indel_pos = indel_pos;
		*indels = +sep;
		*nmismatches_longcont = conti;
		*nmismatches_shift = shifti + 1;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	  conti++;
	}
	debug2e(printf("\n"));

#else
	shifti = nmismatches_avail_shift - 1;
	conti = 0;
	while (conti < nmismatches_avail_long && mismatch_positions_long[conti] < mismatch_positions_shift[shifti]) {
	  conti++;
	}
	indel_pos_cont = mismatch_positions_long[conti];
	indel_pos_shift = mismatch_positions_shift[shifti] - sep + 1;

	while (conti < nmismatches_avail_long && shifti >= 0) {
	  if (indel_pos_cont < indel_pos_shift) {
	    sum = conti + shifti + 1;
	    debug2e(printf("cont %d=%d+%d at indel_pos %d.  ",sum,conti,shifti+1,mismatch_positions_long[conti]));
	    if (sum < best_sum) {
	      endlength = querylength - indel_pos_cont - sep;
	      if (endlength >= min_indel_end_matches) {
		nmatches = endlength - (shifti + 1);
		if (nmatches - 3*(shifti+1) - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti+1,nmatches-3*(shifti+1)-4));
		  best_indel_pos = indel_pos_cont;
		  *indels = +sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti + 1;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    conti++;
	    indel_pos_cont = mismatch_positions_long[conti];

	  } else if (indel_pos_shift < indel_pos_cont) {
	    sum = shifti + conti;
	    debug2e(printf("shift %d=%d+%d at indel_pos %d.  ",sum,conti,shifti,mismatch_positions_shift[shifti]-sep+1));
	    if (sum < best_sum) {
	      endlength = querylength - (indel_pos_shift + sep);
	      if (endlength >= min_indel_end_matches) {
		nmatches = endlength - shifti;
		if (nmatches - 3*shifti - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		  best_indel_pos = indel_pos_shift;
		  *indels = +sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    if (--shifti >= 0) {
	      indel_pos_shift = mismatch_positions_shift[shifti] - sep + 1;
	    }

	  } else {
	    sum = shifti + conti;
	    debug2e(printf("both %d=%d+%d at indel_pos %d.  ",sum,conti,shifti,mismatch_positions_shift[shifti]-sep+1));
	    if (sum < best_sum) {
	      endlength = querylength - (indel_pos_shift + sep);
	      if (endlength >= min_indel_end_matches) {
		nmatches = endlength - shifti;
		if (nmatches - 3*shifti - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		  best_indel_pos = indel_pos_shift;
		  *indels = +sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    conti++;
	    if (--shifti >= 0) {
	      indel_pos_cont = mismatch_positions_long[conti];
	      indel_pos_shift = mismatch_positions_shift[shifti] - sep + 1;
	    }
	  }
	}

	if (shifti < 0) {
	  sum = conti /*+ shifti + 1*/;
	  debug2e(printf("last %d=%d at indel_pos %d.  ",sum,conti,mismatch_positions_long[conti]));
	  if (sum < best_sum) {
	    endlength = querylength - indel_pos_cont - sep;
	    if (endlength >= min_indel_end_matches) {
	      nmatches = endlength /*- (shifti + 1)*/;
	      if (nmatches >= /*shifti + 1*/ + 4) {
		/* Want more matches than mismatches */
		best_indel_pos = indel_pos_cont;
		*indels = +sep;
		*nmismatches_longcont = conti;
		*nmismatches_shift = 0 /*shifti + 1*/;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	}

	debug2e(printf("\n"));
#endif
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions_shift);
  } else {
    FREE(mismatch_positions_shift);
  }
#else
  FREE(mismatch_positions_shift);
#endif

  debug2e(printf("compute_end_indels_right returning with nmismatches_longcont %d + nmismatches_shift %d for %d indels at indel_pos %d\n",
		 *nmismatches_longcont,*nmismatches_shift,*indels,best_indel_pos));

  return best_indel_pos;
}


/* Want genomic low position for indel, so check deletions and insertions in reverse order of preference
   and check for sum <= best_sum */
static int
compute_end_indels_left (int *indels, int *nmismatches_longcont, int *nmismatches_shift,
			 int *mismatch_positions_long, int nmismatches_avail_long,
			 int breakpoint, int querylength, Univcoord_T left, Compress_T query_compress,
			 int min_indel_end_matches, int max_end_insertions, int max_end_deletions,
			 int max_mismatches_short, bool plusp, int genestrand) {
#ifdef DEBUG2E
  int i;
#endif
  int length1;
  int sep;
  int nmismatches_avail_shift, nmatches;
  int sum, best_sum = querylength;
  int conti, shifti;
  int best_indel_pos = -1;
#ifdef OLD_END_INDELS
  int indel_pos;
#else
  int indel_pos_cont, indel_pos_shift;
#endif
  int *mismatch_positions_shift;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions_shift = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions_shift = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions_shift = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


  debug2e(printf("Entered compute_end_indels_left with breakpoint = %d, max_mismatches_short %d\n",
		 breakpoint,max_mismatches_short));
  length1 = breakpoint;
#if 0
  /* Should not need to reset max_end_deletions */
  if (max_end_deletions > length1 - min_indel_end_matches) {
    max_end_deletions = length1 - min_indel_end_matches;
    debug2e(printf("Resetting max_end_deletions to be %d - %d = %d\n",length1,min_indel_end_matches,max_end_deletions));
  }
#endif
  if (max_end_insertions > length1 - min_indel_end_matches) {
    max_end_insertions = length1 - min_indel_end_matches;
  }

  if (max_end_insertions > 0) {
    for (sep = max_end_insertions; sep >= 1; sep--) {
      /* *indels = +sep; */
      nmismatches_avail_shift = Genome_mismatches_left(mismatch_positions_shift,
						       max_mismatches_short,query_compress,
						       left+(+sep),/*pos5*/0,/*pos3*/querylength,
						       plusp,genestrand);
      debug2e(
	      printf("D. Trying insertion of %d.  ",+sep);
	      printf("%d mismatches on left at:",nmismatches_avail_shift);
	      for (i = 0; i <= nmismatches_avail_shift; i++) {
		printf(" %d",mismatch_positions_shift[i]);
	      }
	      printf(".  ");
	      );

      if (nmismatches_avail_shift == 0) {
	/* Skip, because we have an exact match on the shift */
      } else {
#ifdef OLD_END_INDELS
	/* Compute over 0 to mismatch_positions_shift[n-1] */
	/* D. Left insertion.  First, try primary loop is on Genome_mismatches_right (conti), to get lowest coordinate */
	shifti = nmismatches_avail_shift - 1;
	conti = 0;

	/* Start in region */
	while (conti < nmismatches_avail_long && mismatch_positions_long[conti] > mismatch_positions_shift[shifti]) {
	  conti++;
	}

	while (conti < nmismatches_avail_long) {
	  while (shifti >= 0 && mismatch_positions_shift[shifti] > mismatch_positions_long[conti] - sep) {
	    shifti--;
	  }
	  sum = conti + shifti + 1;
	  debug2e(printf("sum %d=%d+%d at indel_pos %d.  ",sum,conti,shifti+1,mismatch_positions_long[conti]-sep+1));
	  if (sum < best_sum) {
	    indel_pos = mismatch_positions_long[conti] - sep + 1;
	    if (indel_pos >= min_indel_end_matches) {
	      nmatches = indel_pos - (shifti + 1);
	      if (nmatches - 3*(shifti+1)- 4 >= 0) {
		/* Want more matches than mismatches */
		debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti+1,nmatches-3*(shifti+1)-4));
		best_indel_pos = indel_pos;
		*indels = +sep;
		*nmismatches_longcont = conti;
		*nmismatches_shift = shifti + 1;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	  conti++;
	}
	debug2e(printf("\n"));


	/* D. Left insertion.  Then, try primary loop is on Genome_mismatches_left (shifti), to see if we missed anything */
	shifti = 0;
	conti = nmismatches_avail_long - 1;
	while (shifti < nmismatches_avail_shift) {
	  while (conti >= 0 && mismatch_positions_long[conti] < mismatch_positions_shift[shifti] + sep) {
	    conti--;
	  }
	  sum = shifti + conti + 1;
	  debug2e(printf("sum %d=%d+%d at indel_pos %d.  ",sum,conti+1,shifti,mismatch_positions_shift[shifti]));
	  if (sum < best_sum) {
	    indel_pos = mismatch_positions_shift[shifti];
	    if (indel_pos >= min_indel_end_matches) {
	      nmatches = indel_pos - shifti;
	      if (nmatches - 3*shifti - 4 >= 0) {
		/* Want more matches than mismatches */
		debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		best_indel_pos = indel_pos;
		*indels = +sep;
		*nmismatches_longcont = conti + 1;
		*nmismatches_shift = shifti;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	  shifti++;
	}
	debug2e(printf("\n"));

#else
	shifti = nmismatches_avail_shift - 1;
	conti = 0;
	while (conti < nmismatches_avail_long && mismatch_positions_long[conti] > mismatch_positions_shift[shifti]) { /* ? >= */
	  conti++;
	}
	indel_pos_cont = mismatch_positions_long[conti] - sep + 1;
	indel_pos_shift = mismatch_positions_shift[shifti];

	while (conti < nmismatches_avail_long && shifti >= 0) {
	  if (indel_pos_cont > indel_pos_shift) {
	    sum = conti + shifti + 1;
	    debug2e(printf("cont %d=%d+%d at indel_pos %d.  ",sum,conti,shifti+1,mismatch_positions_long[conti]-sep+1));
	    if (sum <= best_sum) {
	      if (indel_pos_cont >= min_indel_end_matches) {
		nmatches = indel_pos_cont - (shifti + 1);
		if (nmatches - 3*(shifti+1) - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti+1,nmatches-3*(shifti+1)-4));
		  best_indel_pos = indel_pos_cont;
		  *indels = +sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti + 1;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    conti++;
	    indel_pos_cont = mismatch_positions_long[conti] - sep + 1;

	  } else if (indel_pos_shift > indel_pos_cont) {
	    sum = shifti + conti;
	    debug2e(printf("shift %d=%d+%d at indel_pos %d.  ",sum,conti,shifti,mismatch_positions_shift[shifti]));
	    if (sum <= best_sum) {
	      if (indel_pos_shift >= min_indel_end_matches) {
		nmatches = indel_pos_shift - shifti;
		if (nmatches - 3*shifti - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		  best_indel_pos = indel_pos_shift;
		  *indels = +sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    if (--shifti >= 0) {
	      indel_pos_shift = mismatch_positions_shift[shifti];
	    }

	  } else {
	    sum = shifti + conti;
	    debug2e(printf("both %d=%d+%d at indel_pos %d.  ",sum,conti,shifti,mismatch_positions_shift[shifti]));
	    if (sum <= best_sum) {
	      if (indel_pos_shift >= min_indel_end_matches) {
		nmatches = indel_pos_shift - shifti;
		if (nmatches - 3*shifti - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		  best_indel_pos = indel_pos_shift;
		  *indels = +sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    conti++;
	    if (--shifti >= 0) {
	      indel_pos_cont = mismatch_positions_long[conti] - sep + 1;
	      indel_pos_shift = mismatch_positions_shift[shifti];
	    }
	  }
	}

	if (shifti < 0) {
	  sum = conti /*+ shifti + 1*/;
	  debug2e(printf("last %d=%d at indel_pos %d.  ",sum,conti,mismatch_positions_long[conti]-sep+1));
	  if (sum <= best_sum) {
	    if (indel_pos_cont >= min_indel_end_matches) {
	      nmatches = indel_pos_cont /*- (shifti + 1)*/;
	      if (nmatches >= /*shifti + 1*/ + 4) {
		/* Want more matches than mismatches */
		best_indel_pos = indel_pos_cont;
		*indels = +sep;
		*nmismatches_longcont = conti;
		*nmismatches_shift = 0 /*shifti + 1*/;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	}

	debug2e(printf("\n"));
#endif
      }
    }
  }


  if (max_end_deletions > 0) {
#if 0    
    /* start is never used */
    if (left < (unsigned int) max_end_deletions) {
      debug2e(printf("left %llu < max_end_deletions %d, so start = left\n",
		     (unsigned long long) left,max_end_deletions));
      start = left;
    } else {
      start = 1;
    }
#endif

    for (sep = max_end_deletions; sep >= 1; sep--) {
      /* *indels = -sep; */
      nmismatches_avail_shift = Genome_mismatches_left(mismatch_positions_shift,
						       max_mismatches_short,query_compress,
						       left+(-sep),/*pos5*/0,/*pos3*/querylength,
						       plusp,genestrand);
      debug2e(
	      printf("C. Trying deletion of %d.  ",-sep);
	      printf("%d mismatches on left at:",nmismatches_avail_shift);
	      for (i = 0; i <= nmismatches_avail_shift; i++) {
		printf(" %d",mismatch_positions_shift[i]);
	      }
	      printf(".  ");
	      );

      if (nmismatches_avail_shift == 0) {
	/* Skip, because we have an exact match on the shift */
      } else {
#ifdef OLD_END_INDELS
	/* Compute over 0 to mismatch_positions_shift[n-1] */
	/* C. Left deletion.  First, try primary loop (cont) on Genome_mismatches_right, to get lowest coordinate */
	shifti = nmismatches_avail_shift - 1;
	conti = 0;

	/* Start in region */
	while (conti < nmismatches_avail_long && mismatch_positions_long[conti] > mismatch_positions_shift[shifti]) {
	  conti++;
	}

	while (conti < nmismatches_avail_long) {
	  while (shifti >= 0 && mismatch_positions_shift[shifti] > mismatch_positions_long[conti]) {
	    shifti--;
	  }
	  sum = conti + shifti + 1;
	  debug2e(printf("sum %d=%d+%d at indel_pos %d.  ",sum,conti,shifti+1,mismatch_positions_long[conti]+1));
	  if (sum < best_sum) {
	    indel_pos = mismatch_positions_long[conti] + 1;
	    if (indel_pos >= min_indel_end_matches && indel_pos >= sep) {
	      nmatches = indel_pos - (shifti + 1);
	      if (nmatches - 3*(shifti+1) - 4 >= 0) {
		/* Want more matches than mismatches */
		debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti+1,nmatches-3*(shifti+1)-4));
		best_indel_pos = indel_pos;
		*indels = -sep;
		*nmismatches_longcont = conti;
		*nmismatches_shift = shifti + 1;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	  conti++;
	}
	debug2e(printf("\n"));


	/* C. Left deletion.  Then, try primary loop on Genome_mismatches_left (shifti) to see if we missed anything */
	shifti = 0;
	conti = nmismatches_avail_long - 1;
	while (shifti < nmismatches_avail_shift) {
	  while (conti >= 0 && mismatch_positions_long[conti] < mismatch_positions_shift[shifti]) {
	    conti--;
	  }
	  sum = shifti + conti + 1;
	  debug2e(printf("sum %d=%d+%d at indel_pos %d.  ",sum,conti+1,shifti,mismatch_positions_shift[shifti]));
	  if (sum < best_sum) {
	    indel_pos = mismatch_positions_shift[shifti];
	    if (indel_pos >= min_indel_end_matches && indel_pos >= sep) {
	      nmatches = indel_pos - shifti;
	      if (nmatches - 3*shifti - 4 >= 0) {
		/* Want more matches than mismatches */
		debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		best_indel_pos = indel_pos;
		*indels = -sep;
		*nmismatches_longcont = conti + 1;
		*nmismatches_shift = shifti;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	  shifti++;
	}
	debug2e(printf("\n"));

#else
	shifti = nmismatches_avail_shift - 1;
	conti = 0;
	while (conti < nmismatches_avail_long && mismatch_positions_long[conti] > mismatch_positions_shift[shifti]) { /* ? >= */
	  conti++;
	}
	indel_pos_cont = mismatch_positions_long[conti] + 1;
	indel_pos_shift = mismatch_positions_shift[shifti];

	while (conti < nmismatches_avail_long && shifti >= 0) {
	  if (indel_pos_cont > indel_pos_shift) {
	    sum = conti + shifti + 1;
	    debug2e(printf("cont %d=%d+%d at indel_pos %d.  ",sum,conti,shifti+1,mismatch_positions_long[conti]+1));
	    if (sum <= best_sum) {
	      if (indel_pos_cont >= min_indel_end_matches && indel_pos_cont >= sep) {
		nmatches = indel_pos_cont - (shifti + 1);
		if (nmatches - 3*(shifti+1) - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti+1,nmatches-3*(shifti+1)-4));
		  best_indel_pos = indel_pos_cont;
		  *indels = -sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti + 1;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    conti++;
	    indel_pos_cont = mismatch_positions_long[conti] + 1;

	  } else if (indel_pos_shift > indel_pos_cont) {
	    sum = shifti + conti;
	    debug2e(printf("shift %d=%d+%d at indel_pos %d.  ",sum,conti,shifti,mismatch_positions_shift[shifti]));
	    if (sum <= best_sum) {
	      if (indel_pos_shift >= min_indel_end_matches && indel_pos_shift >= sep) {
		nmatches = indel_pos_shift - shifti;
		if (nmatches - 3*shifti - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		  best_indel_pos = indel_pos_shift;
		  *indels = -sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    if (--shifti >= 0) {
	      indel_pos_shift = mismatch_positions_shift[shifti];
	    }

	  } else {
	    sum = shifti + conti;
	    debug2e(printf("both %d=%d+%d at indel_pos %d.  ",sum,conti,shifti,mismatch_positions_shift[shifti]));
	    if (sum <= best_sum) {
	      if (indel_pos_shift >= min_indel_end_matches && indel_pos_shift >= sep) {
		nmatches = indel_pos_shift - shifti;
		if (nmatches - 3*shifti - 4 >= 0) {
		  /* Want more matches than mismatches */
		  debug2e(printf("nmatches %d - 3*%d - 4 = %d.  ",nmatches,shifti,nmatches-3*shifti-4));
		  best_indel_pos = indel_pos_shift;
		  *indels = -sep;
		  *nmismatches_longcont = conti;
		  *nmismatches_shift = shifti;
		  debug2e(printf("**"));
		  best_sum = sum;
		}
	      }
	    }
	    conti++;
	    if (--shifti >= 0) {
	      indel_pos_cont = mismatch_positions_long[conti] + 1;
	      indel_pos_shift = mismatch_positions_shift[shifti];
	    }
	  }
	}

	if (shifti < 0) {
	  sum = conti /*+ shifti + 1*/;
	  debug2e(printf("last %d=%d at indel_pos %d.  ",sum,conti,mismatch_positions_long[conti]+1));
	  if (sum <= best_sum) {
	    if (indel_pos_cont >= min_indel_end_matches && indel_pos_cont >= sep) {
	      nmatches = indel_pos_cont /*- (shifti + 1)*/;
	      if (nmatches >= /*shifti + 1*/ + 4) {
		/* Want more matches than mismatches */
		best_indel_pos = indel_pos_cont;
		*indels = -sep;
		*nmismatches_longcont = conti;
		*nmismatches_shift = 0 /*shifti + 1*/;
		debug2e(printf("**"));
		best_sum = sum;
	      }
	    }
	  }
	}

	debug2e(printf("\n"));
#endif
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions_shift);
  } else {
    FREE(mismatch_positions_shift);
  }
#else
  FREE(mismatch_positions_shift);
#endif

  debug2e(printf("compute_end_indels_left returning with nmismatches_cont %d + nmismatches_shift %d for %d indels at indel_pos %d\n",
		 *nmismatches_longcont,*nmismatches_shift,*indels,best_indel_pos));

  return best_indel_pos;
}


/************************************************************************/

/* Was solve_first_indel_plus and solve_last_indel_minus */
static List_T
solve_end_indel_low (int *found_score, int *nhits, List_T hits, Segment_T ptr,
		     Univcoord_T diagonal, int firstbound,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
#ifdef DEBUG2E
		     char *queryptr,
#endif
		     int querylength, Compress_T query_compress,
		     int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		     int indel_penalty_end, int max_mismatches, bool plusp, int genestrand) {
#ifdef DEBUG2E
  char *gbuffer;
#endif
  int i;
  Stage3end_T hit;
  Univcoord_T left;
  int indels, query_indel_pos, indel_pos, breakpoint;
  int nmismatches, nmismatches_long, nmismatches_longcont, nmismatches_shift;
  int nmismatches1, nmismatches2;
  int *mismatch_positions;


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


  left = diagonal - querylength;
  if ((unsigned int) max_end_deletions > left - chroffset) {
    max_end_deletions = left - chroffset;
    /* diagonal - querylength guaranteed to be >= chroffset, so max_end_deletions >= 0 */
  }

  debug2e(
	  if (plusp == true) {
	    printf("\nsolve_end_indel_low: Getting genome at diagonal %llu - querylength %d - max_end_deletions %d = %llu.\n",
		   (unsigned long long) diagonal,querylength,max_end_deletions,(unsigned long long) (left-max_end_deletions));
	  } else {
	    printf("\nsolve_end_indel_low: Getting genome at diagonal %llu + 12 - querylength %d = %llu, max_end_deletions = %d.\n",
		   (unsigned long long) diagonal,querylength,(unsigned long long) left,max_end_deletions);
	  });

  debug2e(gbuffer = (char *) CALLOC(querylength+max_end_deletions+1,sizeof(char)));
  debug2e(Genome_fill_buffer_blocks(left-max_end_deletions,querylength+max_end_deletions,gbuffer));
  debug2e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
  debug2e(FREE(gbuffer));

  /* No need to check chromosome bounds */
  nmismatches = Genome_mismatches_right(mismatch_positions,max_mismatches,
					query_compress,left,/*pos5*/0,/*pos3*/querylength,
					plusp,genestrand);

  debug2e(
	  printf("full read: %d (max %d) mismatches from right:",nmismatches,max_mismatches);
	  for (i = 0; i <= nmismatches; i++) {
	    printf(" %d",mismatch_positions[i]);
	  }
	  printf("\n");
	  );

  /* Find first mismatch past firstbound */
  i = 0;
  while (i <= nmismatches && mismatch_positions[i] > firstbound) {
    i++;
  }
  nmismatches_long = i;

#if 0
  nmismatches_short = nmismatches - i;
  /* Previously checked if nmismatches_short <= 0 */
#endif
  
  /* Should be >= */
  if (i >= nmismatches) {
    debug2e(printf("i %d >= nmismatches %d => no indel\n",i,nmismatches));
  } else {
#if 0
    mismatch_positions_short = &(mismatch_positions[i]);
    debug2e(
	    printf("nmismatches_long = %d, short = %d\n",nmismatches_long,nmismatches_short);
	    printf("short end from firstbound %d:",firstbound);
	    for (i = 0; i <= nmismatches_short; i++) {
	      printf(" %d",mismatch_positions_short[i]);
	    }
	    printf("\n");
	    );
    breakpoint = mismatch_positions_short[0] + 1;
#else
    breakpoint = mismatch_positions[i] + 1;
#endif

    if ((indel_pos = compute_end_indels_left(&indels,&nmismatches_longcont,&nmismatches_shift,
					     mismatch_positions,nmismatches,
					     breakpoint,querylength,left,query_compress,
					     min_indel_end_matches,max_end_insertions,max_end_deletions,
					     /*max_mismatches_allowed*/max_mismatches-nmismatches_long+1,
					     plusp,genestrand)) >= 0) {
      debug2e(printf("Got indel_pos %d.\n",indel_pos));

      /* formulas for query_indel_pos have been checked on examples */
      if (indels > 0) {
	if (plusp == true) {
	  query_indel_pos = indel_pos /* + indels */;
	  debug2e(printf("case 1: query_indel_pos = %d = %d\n",indel_pos,query_indel_pos));
	  nmismatches1 = nmismatches_shift;
	  nmismatches2 = nmismatches_longcont;
	  /* end1_indel_p = true; */
	  /* end2_indel_p = false; */
	} else {
	  query_indel_pos = querylength - indel_pos - indels;
	  debug2e(printf("case 2: query_indel_pos = %d - %d - %d = %d\n",querylength,indel_pos,indels,query_indel_pos));
	  nmismatches1 = nmismatches_longcont;
	  nmismatches2 = nmismatches_shift;
	  /* end1_indel_p = false; */
	  /* end2_indel_p = true; */
	}

	/* genomiclength = querylength-indels; */
	if ((hit = Stage3end_new_insertion(&(*found_score),indels,query_indel_pos,
					   nmismatches1,nmismatches2,left+indels,
					   query_compress,querylength,plusp,genestrand,
					   chrnum,chroffset,chrhigh,chrlength,indel_penalty_end,
					   /*sarrayp*/false)) != NULL) {
	  ptr->usedp = true;
	  *nhits += 1;
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful insertion at %d with %d long/cont + %d shift mismatches\n",
			 indel_pos,nmismatches_longcont,nmismatches_shift));
	}

      } else {
	if (plusp == true) {
	  query_indel_pos = indel_pos;
	  nmismatches1 = nmismatches_shift;
	  nmismatches2 = nmismatches_longcont;
	  /* end1_indel_p = true; */
	  /* end2_indel_p = false; */
	} else {
	  query_indel_pos = querylength - indel_pos;
	  nmismatches1 = nmismatches_longcont;
	  nmismatches2 = nmismatches_shift;
	  /* end1_indel_p = false; */
	  /* end2_indel_p = true; */
	}

	/*genomiclength = querylength-indels; */
	if ((hit = Stage3end_new_deletion(&(*found_score),-indels,query_indel_pos,
					  nmismatches1,nmismatches2,left+indels,
					  query_compress,querylength,plusp,genestrand,
					  chrnum,chroffset,chrhigh,chrlength,indel_penalty_end,
					  /*sarrayp*/false)) != NULL) {
	  ptr->usedp = true;
	  *nhits += 1;
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful end deletion at %d with %d long/cont + %d shift mismatches and nindels %d\n",
			 indel_pos,nmismatches_longcont,nmismatches_shift,-indels));
	}
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
  } else {
    FREE(mismatch_positions);
  }
#else
  FREE(mismatch_positions);
#endif

  return hits;
}


/* Was solve_first_indel_minus and solve_last_indel_plus */
static List_T
solve_end_indel_high (int *found_score, int *nhits, List_T hits, Segment_T ptr,
		      Univcoord_T diagonal, int lastbound,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
#ifdef DEBUG2E
		      char *queryptr,
#endif
		      int querylength, Compress_T query_compress,
		      int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		      int indel_penalty_end, int max_mismatches, bool plusp, int genestrand) {
#ifdef DEBUG2E
  char *gbuffer;
#endif
  int i;
  Stage3end_T hit;
  Univcoord_T left;
  int indels, query_indel_pos, indel_pos, breakpoint;
  int nmismatches, nmismatches_long, nmismatches_longcont, nmismatches_shift;
  int nmismatches1, nmismatches2;
  int *mismatch_positions;


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


  left = diagonal - querylength;
  if ((unsigned int) max_end_deletions > chrhigh - diagonal) {
    max_end_deletions = chrhigh - diagonal;
    /* diagonal guaranteed to be <= chrhigh, so max_end_deletions >= 0 */
  }

  debug2e(
	  if (plusp == true) {
	    printf("\nsolve_end_indel_high: Getting genome at diagonal %llu - querylength %d + max_end_deletions %d = %llu.\n",
		   (unsigned long long) diagonal,querylength,max_end_deletions,(unsigned long long) (left+max_end_deletions));
	  } else {
	    printf("\nsolve_end_indel_high: Getting genome at diagonal %llu + 12 - querylength %d = %llu, max_end_deletions = %d.\n",
		   (unsigned long long) diagonal,querylength,(unsigned long long) left,max_end_deletions);
	  });

  debug2e(gbuffer = (char *) CALLOC(querylength+max_end_deletions+1,sizeof(char)));
  debug2e(Genome_fill_buffer_blocks(left,querylength+max_end_deletions,gbuffer));
  debug2e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
  debug2e(FREE(gbuffer));

  /* No need to check chromosome bounds */
  /* Previously checked from 0 to lastbound */
  nmismatches = Genome_mismatches_left(mismatch_positions,max_mismatches,
				       query_compress,left,/*pos5*/0,/*pos3*/querylength,
				       plusp,genestrand);

  debug2e(
	  printf("full read: %d (max %d) mismatches from left:",nmismatches,max_mismatches);
	  for (i = 0; i <= nmismatches; i++) {
	    printf(" %d",mismatch_positions[i]);
	  }
	  printf("\n");
	  );

  /* Find first mismatch past lastbound */
  i = 0;
  while (i <= nmismatches && mismatch_positions[i] < lastbound) {
    i++;
  }
  nmismatches_long = i;
#if 0
  /* Previously checked if nmismatches_short <= 0 */
  nmismatches_short = nmismatches - i;
#endif
  
  /* Should be >= */
  if (i >= nmismatches) {
    debug2e(printf("i %d >= nmismatches %d => no indel\n",i,nmismatches));
  } else {
#if 0
    mismatch_positions_short = &(mismatch_positions[i]);
    debug2e(
	    printf("nmismatches_long = %d, short = %d\n",nmismatches_long,nmismatches_short);
	    printf("short end from lastbound %d:",lastbound);
	    for (i = 0; i <= nmismatches_short; i++) {
	      printf(" %d",mismatch_positions_short[i]);
	    }
	    printf("\n");
	    );
    breakpoint = mismatch_positions_short[0] - 1;
#else
    breakpoint = mismatch_positions[i] - 1;
#endif

    if ((indel_pos = compute_end_indels_right(&indels,&nmismatches_longcont,&nmismatches_shift,
					      mismatch_positions,nmismatches,
					      breakpoint,querylength,left,query_compress,
					      min_indel_end_matches,max_end_insertions,max_end_deletions,
					      /*max_mismatches_allowed*/max_mismatches-nmismatches_long+1,
					      plusp,genestrand)) >= 0) {
      debug2e(printf("Got indel_pos %d\n",indel_pos));

      /* formulas for query_indel_pos have been checked on examples */
      if (indels > 0) {
	if (plusp == true) {
	  query_indel_pos = indel_pos /* + indels */;
	  debug2e(printf("case 3: query_indel_pos = %d = %d\n",indel_pos,query_indel_pos));
	  nmismatches1 = nmismatches_longcont;
	  nmismatches2 = nmismatches_shift;
	  /* end1_indel_p = false; */
	  /* end2_indel_p = true; */
	} else {
	  query_indel_pos = querylength - indel_pos - indels;
	  debug2e(printf("case 4: query_indel_pos = %d - %d - %d = %d\n",querylength,indel_pos,indels,query_indel_pos));
	  nmismatches1 = nmismatches_shift;
	  nmismatches2 = nmismatches_longcont;
	  /* end1_indel_p = true; */
	  /* end2_indel_p = false; */
	}

	/* genomiclength = querylength-indels; */
	if ((hit = Stage3end_new_insertion(&(*found_score),indels,query_indel_pos,
					   nmismatches1,nmismatches2,left,
					   query_compress,querylength,plusp,genestrand,
					   chrnum,chroffset,chrhigh,chrlength,indel_penalty_end,
					   /*sarrayp*/false)) != NULL) {
	  ptr->usedp = true;
	  *nhits += 1;
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful end insertion at %d with %d long/cont + %d shift mismatches\n",
			 indel_pos,nmismatches_longcont,nmismatches_shift));
	}

      } else {
	if (plusp == true) {
	  query_indel_pos = indel_pos;
	  nmismatches1 = nmismatches_longcont;
	  nmismatches2 = nmismatches_shift;
	  /* end1_indel_p = false; */
	  /* end2_indel_p = true; */
	} else {
	  query_indel_pos = querylength - indel_pos;
	  nmismatches1 = nmismatches_shift;
	  nmismatches2 = nmismatches_longcont;
	  /* end1_indel_p = true; */
	  /* end2_indel_p = false; */
	}

	/*genomiclength = querylength-indels; */
	if ((hit = Stage3end_new_deletion(&(*found_score),-indels,query_indel_pos,
					  nmismatches1,nmismatches2,left,
					  query_compress,querylength,plusp,genestrand,
					  chrnum,chroffset,chrhigh,chrlength,indel_penalty_end,
					  /*sarrayp*/false)) != NULL) {
	  ptr->usedp = true;
	  *nhits += 1;
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful end deletion at %d with %d long/cont + %d shift mismatches and nindels %d\n",
			 indel_pos,nmismatches_longcont,nmismatches_shift,-indels));
	}
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
  } else {
    FREE(mismatch_positions);
  }
#else
  FREE(mismatch_positions);
#endif

  return hits;
}


/* Note: plus_anchor_segments and minus_anchor_segments point to anchors,
   but can use smaller segments for the ends because ptr points to all
   of them */
static List_T
find_end_indels (int *found_score, int *nhits, List_T hits,
		 Segment_T *plus_anchor_segments, Segment_T *minus_anchor_segments,
		 int n_plus_anchors, int n_minus_anchors,
#ifdef DEBUG2E
		 char *queryuc_ptr, char *queryrc,
#endif
		 int querylength, int firstbound, int lastbound,
		 Compress_T query_compress_fwd, Compress_T query_compress_rev,
		 int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		 int indel_penalty_end, int max_mismatches_allowed, int genestrand) {
  Segment_T ptr, *p;

  debug(printf("*** find_end_indels with max_mismatches_allowed %d ***\n",
	       max_mismatches_allowed));

  for (p = &(plus_anchor_segments[0]); p < &(plus_anchor_segments[n_plus_anchors]); p++) {
    ptr = *p;
    if (ptr->diagonal < (Univcoord_T) -1) {

      if (ptr->floor_xfirst <= max_mismatches_allowed) {
	/* First indel, plus */
	debug2e(printf("floor_xfirst %d <= mismatches allowed %d\n",ptr->floor_xfirst,max_mismatches_allowed));
	hits = solve_end_indel_low(&(*found_score),&(*nhits),hits,ptr,ptr->diagonal,firstbound,
				   ptr->chrnum,ptr->chroffset,ptr->chrhigh,ptr->chrlength,
#ifdef DEBUG2E
				   /*queryptr*/queryuc_ptr,
#endif
				   querylength,/*query_compress*/query_compress_fwd,
				   max_end_insertions,max_end_deletions,min_indel_end_matches,
				   indel_penalty_end,max_mismatches_allowed,/*plusp*/true,genestrand);
      }

      if (ptr->floor_xlast <= max_mismatches_allowed) {
	/* Last indel, plus */
	debug2e(printf("floor_xlast %d <= mismatches allowed %d\n",ptr->floor_xlast,max_mismatches_allowed));
	hits = solve_end_indel_high(&(*found_score),&(*nhits),hits,ptr,ptr->diagonal,lastbound,
				    ptr->chrnum,ptr->chroffset,ptr->chrhigh,ptr->chrlength,
#ifdef DEBUG2E
				    /*queryptr*/queryuc_ptr,
#endif
				    querylength,/*query_compress*/query_compress_fwd,
				    max_end_insertions,max_end_deletions,min_indel_end_matches,
				    indel_penalty_end,max_mismatches_allowed,/*plusp*/true,genestrand);
      }
    }
  }

  for (p = &(minus_anchor_segments[0]); p < &(minus_anchor_segments[n_minus_anchors]); p++) {
    ptr = *p;
    if (ptr->diagonal < (Univcoord_T) -1) {

      if (ptr->floor_xfirst <= max_mismatches_allowed) {
	/* First indel, minus */
	debug2e(printf("floor_xfirst %d <= mismatches allowed %d\n",ptr->floor_xfirst,max_mismatches_allowed));
	hits = solve_end_indel_high(&(*found_score),&(*nhits),hits,ptr,ptr->diagonal,lastbound,
				    ptr->chrnum,ptr->chroffset,ptr->chrhigh,ptr->chrlength,
#ifdef DEBUG2E
				    /*queryptr*/queryrc,
#endif
				    querylength,/*query_compress*/query_compress_rev,
				    max_end_insertions,max_end_deletions,min_indel_end_matches,
				    indel_penalty_end,max_mismatches_allowed,/*plusp*/false,genestrand);
      }

      if (ptr->floor_xlast <= max_mismatches_allowed) {
	/* Last indel, minus */
	debug2e(printf("floor_xlast %d <= mismatches allowed %d\n",ptr->floor_xlast,max_mismatches_allowed));
	hits = solve_end_indel_low(&(*found_score),&(*nhits),hits,ptr,ptr->diagonal,firstbound,
				   ptr->chrnum,ptr->chroffset,ptr->chrhigh,ptr->chrlength,
#ifdef DEBUG2E
				   /*queryptr*/queryrc,
#endif
				   querylength,/*query_compress*/query_compress_rev,
				   max_end_insertions,max_end_deletions,min_indel_end_matches,
				   indel_penalty_end,max_mismatches_allowed,/*plusp*/false,genestrand);
      }
    }
  }

  return hits;
}



/************************************************************************
 *   Splicing
 ************************************************************************/


/* Do not compare against true or false */
/* Moderate criterion */
static int
sufficient_splice_prob_distant (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < min_distantsplicing_end_matches) {
    return 0;
  } else if (support < 30) {
    return (spliceprob > 0.95);
  } else if (support < 35) {
    return (spliceprob > 0.90);
  } else if (support < 40) {
    return (spliceprob > 0.85);
  } else {
    return (spliceprob > 0.70);
  }
}

/* Do not compare against true or false */

#ifdef HALFINTRON
/* Strictest criterion */
static int
sufficient_splice_prob_halfintron (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < 20) {
    return 0;
  } else if (support < 26) {
    return (spliceprob > 0.95);
  } else if (support < 32) {
    return (spliceprob > 0.90);
  } else if (support < 38) {
    return (spliceprob > 0.85);
  } else if (support < 44) {
    return (spliceprob > 0.80);
  } else if (support < 50) {
    return (spliceprob > 0.50);
  } else {
    return 1;
  }
}
#endif



/* Copied from sarray-read.c */
static int
donor_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length_orig(Stage3end_substring_donor(x));
  int y_length = Substring_match_length_orig(Stage3end_substring_donor(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    return 0;
  }
}

/* Copied from sarray-read.c */
static int
acceptor_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length_orig(Stage3end_substring_acceptor(x));
  int y_length = Substring_match_length_orig(Stage3end_substring_acceptor(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    return 0;
  }
}


static List_T
find_singlesplices_plus (int *found_score, List_T hits, List_T *ambiguous, List_T *lowprob,
			 Segment_T *plus_spliceable, int plus_nspliceable,
			 Floors_T floors, int querylength, int query_lastpos,
			 Compress_T query_compress /* expecting fwd */,
			 int splicing_penalty, int max_mismatches_allowed, bool first_read_p, int genestrand,
			 bool subs_or_indels_p) {
  int k, j, i, n;
  Segment_T segmenti, segmentj, segmentj_end, *ptr;
  Univcoord_T segmenti_left, segmentj_left;
  int nmismatches_left, nmismatches_right;
  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;
  int *mismatch_positions_left, *mismatch_positions_right;
  int *segmenti_donor_knownpos, *segmentj_acceptor_knownpos, *segmentj_antidonor_knownpos, *segmenti_antiacceptor_knownpos,
    *segmenti_donor_knowni, *segmentj_acceptor_knowni,
    *segmentj_antidonor_knowni, *segmenti_antiacceptor_knowni;
  
#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions_left = (int *) ALLOCA(querylength*sizeof(int));
    mismatch_positions_right = (int *) ALLOCA(querylength*sizeof(int));
    segmenti_donor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions_left = (int *) MALLOC(querylength*sizeof(int));
    mismatch_positions_right = (int *) MALLOC(querylength*sizeof(int));
    segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions_left = (int *) MALLOC(querylength*sizeof(int));
  mismatch_positions_right = (int *) MALLOC(querylength*sizeof(int));
  segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  Chrpos_T max_distance;

  int floor_outer_i;
  int *floors_from_neg3, *floors_to_pos3;
  int nhits_local /*= 0*/;

  List_T accepted_hits, rejected_hits;
  List_T spliceends_sense, spliceends_antisense, p;
  List_T donor_hits, acceptor_hits;
  int donor_length, acceptor_length;
  Stage3end_T hit, *hitarray;
  int n_good_spliceends;
  int best_nmismatches, nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob, donor_prob, acceptor_prob;
  Substring_T donor, acceptor;

#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Intlist_T amb_knowni, amb_nmismatches;
  Doublelist_T amb_probs;


  debug4s(printf("*** Starting find_singlesplices_plus on %d spliceable segments ***\n",plus_nspliceable));
  /* debug(printf("Initially have %d hits\n",List_length(hits))); */

  if (floors != NULL) {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    for (ptr = plus_spliceable; ptr < &(plus_spliceable[plus_nspliceable]); ptr++) {
      segmenti = *ptr;
      debug4s(printf("plus_spliceable segmenti at diagonal %u\n",segmenti->diagonal));
      if (1 || segmenti->diagonal < (Univcoord_T) -1) { /* No markers were stored in spliceable */
	segmenti_left = segmenti->diagonal - querylength;
	floor_outer_i = floors_from_neg3[segmenti->querypos5];

	segmenti_donor_nknown = 0;
	segmenti_antiacceptor_nknown = 0;
	max_distance = shortsplicedist;

	if ((j = segmenti->splicesites_i) >= 0) {
	  /* Ends 1 (donor, plus) and 8 (antiacceptor, plus): mark known splice sites in segmenti */
	  while (j < nsplicesites && splicesites[j] < segmenti->diagonal) {
	    if (splicetypes[j] == DONOR) {
	      debug4s(printf("Setting known donor %d for segmenti at %llu\n",j,(unsigned long long) splicesites[j]));
	      segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - segmenti_left;
	      segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4s(printf("Setting known antiacceptor %d for segmenti at %llu\n",j,(unsigned long long) splicesites[j]));
	      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - segmenti_left;
	      segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	    }

	    /* This computation was already made in identify_all_segments */
	    if (splicedists[j] > max_distance) {
	      debug4s(printf("Setting max_distance for known i %d to be %u\n",j,splicedists[j]));
	      max_distance = splicedists[j];
	    }

	    j++;
	  }
	}
	segmenti_donor_knownpos[segmenti_donor_nknown] = querylength;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength;


	/* Identify potential segmentj for segmenti */
	segmentj_end = segmenti+1;
	while (
#ifdef NO_MARKER_SEGMENTS
	       segmentj_end < &(plus_segments[plus_nsegments]) && segmentj_end->chrnum == segmenti->chrnum &&
#endif
	       segmentj_end->diagonal <= segmenti->diagonal + max_distance) {
	  segmentj_end++;
	}

	spliceends_sense = spliceends_antisense = (List_T) NULL;

	if (segmentj_end - segmenti >= MAX_LOCALSPLICING_POTENTIAL) {
	  /* Too many to check */
	  /* segmentj_end = segmenti+1 + MAX_LOCALSPLICING_POTENTIAL; */
	  segmentj = segmentj_end; /* Don't process any */
	} else {
	  segmentj = segmenti+1;
	}
	for ( ; segmentj < segmentj_end; segmentj++) {
	  debug4s(printf("plus local?  diagonal %llu, querypos %d..%d => diagonal %llu, querypos %d..%d => ",
			 (unsigned long long) segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
			 (unsigned long long) segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	  /* i5 i3 j5 j3 */
	  assert(segmenti->diagonal < segmentj->diagonal);
	  if (segmenti->querypos3 >= segmentj->querypos5) {
	    /* Fail querypos test */
	    debug4s(printf("Bad querypos\n"));

	  } else if (segmenti->diagonal + min_intronlength > segmentj->diagonal) {
	    /* Too short to be an intron */
	    debug4s(printf("Too short\n"));

	  } else {
	    segmenti->right_splice_p = true;
	    segmentj->left_splice_p = true;
	    if (floor_outer_i + floors_to_pos3[segmentj->querypos3] > max_mismatches_allowed) {
	      /* Fail outer floor test */
	      /* floors->score[-index1interval][segmenti->querypos5] +floors->score[segmentj->querypos3][query_lastpos+index1interval] */

	      debug4s(printf("too many mismatches, outer floor = %d+%d=%d > %d\n",
			     floors->scorefrom[-index1interval][segmenti->querypos5],
			     floors->scorefrom[segmentj->querypos3][query_lastpos+index1interval],
			     floors->scorefrom[-index1interval][segmenti->querypos5] +
			     floors->scorefrom[segmentj->querypos3][query_lastpos+index1interval],
			     max_mismatches_allowed));

	    } else {
	      /* Apply leftmost/rightmost test */
	      if (segmenti->leftmost < 0) {
		nmismatches_left = Genome_mismatches_left(mismatch_positions_left,max_mismatches_allowed,
							  query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/querylength,
							  /*plusp*/true,genestrand);
		segmenti->leftmost = (nmismatches_left == 0) ? 0 : mismatch_positions_left[nmismatches_left-1];
		debug4s(printf("%d mismatches on left at:",nmismatches_left);
			for (i = 0; i <= nmismatches_left; i++) {
			  printf(" %d",mismatch_positions_left[i]);
			}
			printf("\n"));
	      }
	  
	      segmentj_left = segmentj->diagonal - querylength;
	      if (segmentj->rightmost < 0) {
		debug4s(printf("left %u\n",segmentj_left));
		nmismatches_right = Genome_mismatches_right(mismatch_positions_right,max_mismatches_allowed,
							    query_compress,/*left*/segmentj_left,/*pos5*/0,/*pos3*/querylength,
							    /*plusp*/true,genestrand);
		segmentj->rightmost = (nmismatches_right == 0) ? 0 : mismatch_positions_right[nmismatches_right-1];
		debug4s(printf("%d mismatches on right [plus] at:",nmismatches_right);
			for (i = 0; i <= nmismatches_right; i++) {
			  printf(" %d",mismatch_positions_right[i]);
			}
			printf("\n"));
	      }
	  
	      debug4s(printf("For a single splice, want leftmost %d > rightmost %d\n",segmenti->leftmost,segmentj->rightmost));
	    
	      if (segmenti->leftmost > segmentj->rightmost) {
		/* Single splice is possible */

		segmentj_acceptor_nknown = 0;
		segmentj_antidonor_nknown = 0;
		if ((j = segmentj->splicesites_i) >= 0) {
		  /* Ends 2 (acceptor, plus) and 7 (antidonor, plus): mark known splice sites in segmentj */
		  while (j < nsplicesites && splicesites[j] < segmentj->diagonal) {
		    if (splicetypes[j] == ACCEPTOR) {
		      debug4s(printf("Setting known acceptor %d for segmentj at %llu\n",j,(unsigned long long) splicesites[j]));
		      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - segmentj_left;
		      segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
		    } else if (splicetypes[j] == ANTIDONOR) {
		      debug4s(printf("Setting known antidonor %d for segmentj at %llu\n",j,(unsigned long long) splicesites[j]));
		      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - segmentj_left;
		      segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
		    }
		    j++;
		  }
		}
		segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength;
		segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength;


		debug4s(printf("  => checking for single splice: Splice_solve_single_plus\n"));
		spliceends_sense =
		  Splice_solve_single_sense(&(*found_score),&nhits_local,spliceends_sense,&(*lowprob),
					    &segmenti->usedp,&segmentj->usedp,
					    /*segmenti_left*/segmenti->diagonal - querylength,
					    /*segmentj_left*/segmentj->diagonal - querylength,
					    segmenti->chrnum,segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
					    segmentj->chrnum,segmentj->chroffset,segmentj->chrhigh,segmentj->chrlength,
					    querylength,query_compress,
					    segmenti_donor_knownpos,segmentj_acceptor_knownpos,
					    segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					    segmenti_donor_knowni,segmentj_acceptor_knowni,
					    segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
					    segmenti_donor_nknown,segmentj_acceptor_nknown,
					    segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
					    splicing_penalty,max_mismatches_allowed,
					    /*plusp*/true,genestrand,first_read_p,subs_or_indels_p,
					    /*sarrayp*/false);
		spliceends_antisense =
		  Splice_solve_single_antisense(&(*found_score),&nhits_local,spliceends_antisense,&(*lowprob),
						&segmenti->usedp,&segmentj->usedp,
						/*segmenti_left*/segmenti->diagonal - querylength,
						/*segmentj_left*/segmentj->diagonal - querylength,
						segmenti->chrnum,segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
						segmentj->chrnum,segmentj->chroffset,segmentj->chrhigh,segmentj->chrlength,
						querylength,query_compress,
						segmenti_donor_knownpos,segmentj_acceptor_knownpos,
						segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
						segmenti_donor_knowni,segmentj_acceptor_knowni,
						segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
						segmenti_donor_nknown,segmentj_acceptor_nknown,
						segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
						splicing_penalty,max_mismatches_allowed,
						/*plusp*/true,genestrand,first_read_p,subs_or_indels_p,
						/*sarrayp*/false);
	      }
	    }
	  }
	}

	/* Process results for segmenti, sense.  Modified from collect_elt_matches in sarray-read.c. */
	if (spliceends_sense != NULL) {
	  /* nmismatches here may be different for spliceends from Splice_solve, so pick based on prob and nmismatches */
	  best_nmismatches = querylength;
	  best_prob = 0.0;
	  for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			  Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
			  Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			  Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
			  Stage3end_nmismatches_whole(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			  Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	    if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	      best_nmismatches = nmismatches;
	    }
	    if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	      best_prob = prob;
	    }
	  }

	  n_good_spliceends = 0;
	  accepted_hits = rejected_hits = (List_T) NULL;
	  for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
		Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	      debug7(printf("accepting distance %d, probabilities %f and %f\n",
			    Stage3end_distance(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			    Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	      n_good_spliceends += 1;
	      accepted_hits = List_push(accepted_hits,(void *) hit);
	    } else {
	      rejected_hits = List_push(rejected_hits,(void *) hit);
	    }
	  }

	  if (n_good_spliceends == 0) {
	    /* Conjunction is too strict.  Allow for disjunction instead. */
	    List_free(&rejected_hits);
	    for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	      hit = (Stage3end_T) List_head(p);
	      if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
		  Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
		debug7(printf("accepting distance %d, probabilities %f and %f\n",
			      Stage3end_distance(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
		n_good_spliceends += 1;
		accepted_hits = List_push(accepted_hits,(void *) hit);
	      } else {
		rejected_hits = List_push(rejected_hits,(void *) hit);
	      }
	    }
	  }

	  for (p = rejected_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    Stage3end_free(&hit);
	  }
	  List_free(&rejected_hits);
	  List_free(&spliceends_sense);

	  if (n_good_spliceends == 1) {
	    hits = List_push(hits,List_head(accepted_hits));
	    List_free(&accepted_hits);

	  } else {
	    /* 1.  Multiple hits, sense, left1 (segmenti_left) */
	    debug7(printf("multiple splice hits, sense, plus\n"));
	    donor_hits = acceptor_hits = (List_T) NULL;

	    /* plus branch from collect_elt_matches */
	    for (p = accepted_hits; p != NULL; p = List_next(p)) {
	      hit = (Stage3end_T) List_head(p);
	      donor = Stage3end_substring_donor(hit);
	      acceptor = Stage3end_substring_acceptor(hit);
	      if (Substring_genomicstart(donor) == segmenti_left) {
		donor_hits = List_push(donor_hits,(void *) hit);
	      } else if (Substring_genomicstart(acceptor) == segmenti_left) {
		acceptor_hits = List_push(acceptor_hits,(void *) hit);
	      } else {
		abort();
		Stage3end_free(&hit);
	      }
	    }

	    if (donor_hits != NULL) {
	      hitarray = (Stage3end_T *) List_to_array_n(&n,donor_hits);
	      qsort(hitarray,n,sizeof(Stage3end_T),donor_match_length_cmp);
	      i = 0;
	      while (i < n) {
		hit = hitarray[i];
		donor = Stage3end_substring_donor(hit);
		donor_length = Substring_match_length_orig(donor);
		j = i + 1;
		while (j < n && Substring_match_length_orig(Stage3end_substring_donor(hitarray[j])) == donor_length) {
		  j++;
		}
		if (j == i + 1) {
		  hits = List_push(hits,(void *) hit);
		} else {
#ifdef LARGE_GENOMES
		  ambcoords = (Uint8list_T) NULL;
#else
		  ambcoords = (Uintlist_T) NULL;
#endif
		  amb_knowni = (Intlist_T) NULL;
		  amb_nmismatches = (Intlist_T) NULL;
		  amb_probs = (Doublelist_T) NULL;

		  for (k = i; k < j; k++) {
		    acceptor = Stage3end_substring_acceptor(hitarray[k]);
#ifdef LARGE_GENOMES
		    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_A(acceptor));
#else
		    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_A(acceptor));
#endif
		    amb_knowni = Intlist_push(amb_knowni,-1);
		    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		    amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
		  }

		  nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
		  donor_prob = Junction_donor_prob(Stage3end_junctionA(hit));
		  prob = best_prob - donor_prob;
		  *ambiguous = List_push(*ambiguous,
					 (void *) Stage3end_new_splice(&(*found_score),
								       /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								       donor,/*acceptor*/NULL,donor_prob,/*acceptor_prob*/prob,/*distance*/0U,
								       /*shortdistancep*/false,/*penalty*/0,querylength,
								       /*ambcoords_donor*/NULL,ambcoords,
								       /*amb_knowni_donor*/NULL,amb_knowni,
								       /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								       /*amb_probs_donor*/NULL,amb_probs,
								       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								       Stage3end_sensedir(hit),/*sarrayp*/false));
		  Doublelist_free(&amb_probs);
		  Intlist_free(&amb_nmismatches);
		  Intlist_free(&amb_knowni);
#ifdef LARGE_GENOMES
		  Uint8list_free(&ambcoords);
#else
		  Uintlist_free(&ambcoords);
#endif

		  for (k = i; k < j; k++) {
		    hit = hitarray[k];
		    Stage3end_free(&hit);
		  }
		}

		i = j;
	      }
	      FREE(hitarray);
	      List_free(&donor_hits);
	    }

	    if (acceptor_hits != NULL) {
	      hitarray = (Stage3end_T *) List_to_array_n(&n,acceptor_hits);
	      qsort(hitarray,n,sizeof(Stage3end_T),acceptor_match_length_cmp);
	      i = 0;
	      while (i < n) {
		hit = hitarray[i];
		acceptor = Stage3end_substring_acceptor(hit);
		acceptor_length = Substring_match_length_orig(acceptor);
		j = i + 1;
		while (j < n && Substring_match_length_orig(Stage3end_substring_acceptor(hitarray[j])) == acceptor_length) {
		  j++;
		}
		if (j == i + 1) {
		  hits = List_push(hits,(void *) hit);
		} else {
#ifdef LARGE_GENOMES
		  ambcoords = (Uint8list_T) NULL;
#else
		  ambcoords = (Uintlist_T) NULL;
#endif
		  amb_knowni = (Intlist_T) NULL;
		  amb_nmismatches = (Intlist_T) NULL;
		  amb_probs = (Doublelist_T) NULL;

		  for (k = i; k < j; k++) {
		    donor = Stage3end_substring_donor(hitarray[k]);
#ifdef LARGE_GENOMES
		    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_D(donor));
#else
		    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_D(donor));
#endif
		    amb_knowni = Intlist_push(amb_knowni,-1);
		    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		    amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
		  }
	    
		  nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
		  acceptor_prob = Junction_acceptor_prob(Stage3end_junctionD(hit));
		  prob = best_prob - acceptor_prob;
		  *ambiguous = List_push(*ambiguous,
					 (void *) Stage3end_new_splice(&(*found_score),
								       nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								       /*donor*/NULL,acceptor,/*donor_prob*/prob,acceptor_prob,/*distance*/0U,
								       /*shortdistancep*/false,/*penalty*/0,querylength,
								       ambcoords,/*ambcoords_acceptor*/NULL,
								       amb_knowni,/*amb_knowni_acceptor*/NULL,
								       amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								       amb_probs,/*amb_probs_acceptor*/NULL,
								       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								       Stage3end_sensedir(hit),/*sarrayp*/false));
		  Doublelist_free(&amb_probs);
		  Intlist_free(&amb_nmismatches);
		  Intlist_free(&amb_knowni);
#ifdef LARGE_GENOMES
		  Uint8list_free(&ambcoords);
#else
		  Uintlist_free(&ambcoords);
#endif

		  for (k = i; k < j; k++) {
		    hit = hitarray[k];
		    Stage3end_free(&hit);
		  }
		}

		i = j;
	      }
	      FREE(hitarray);
	      List_free(&acceptor_hits);
	    }

	    List_free(&accepted_hits);
	  }
	}

	/* Process results for segmenti, antisense.  Modified from collect_elt_matches in sarray-read.c. */
	if (spliceends_antisense != NULL) {
	  /* nmismatches here may be different for spliceends from Splice_solve, so pick based on prob and nmismatches */
	  best_nmismatches = querylength;
	  best_prob = 0.0;
	  for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			  Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
			  Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			  Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
			  Stage3end_nmismatches_whole(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			  Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	    if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	      best_nmismatches = nmismatches;
	    }
	    if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	      best_prob = prob;
	    }
	  }

	  n_good_spliceends = 0;
	  accepted_hits = rejected_hits = (List_T) NULL;
	  for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
		Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	      debug7(printf("accepting distance %d, donor length %d and acceptor length %d, probabilities %f and %f\n",
			    Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			    Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			    Substring_siteD_prob(Stage3end_substring_donor(hit)),
			    Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	      n_good_spliceends += 1;
	      accepted_hits = List_push(accepted_hits,(void *) hit);
	    } else {
	      rejected_hits = List_push(rejected_hits,(void *) hit);
	    }
	  }

	  if (n_good_spliceends == 0) {
	    /* Conjunction is too strict.  Allow for disjunction instead. */
	    List_free(&rejected_hits);
	    for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	      hit = (Stage3end_T) List_head(p);
	      if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
		  Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
		debug7(printf("accepting distance %d, donor length %d and acceptor length %d, probabilities %f and %f\n",
			      Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			      Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			      Substring_siteD_prob(Stage3end_substring_donor(hit)),
			      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
		n_good_spliceends += 1;
		accepted_hits = List_push(accepted_hits,(void *) hit);
	      } else {
		rejected_hits = List_push(rejected_hits,(void *) hit);
	      }
	    }
	  }

	  for (p = rejected_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    Stage3end_free(&hit);
	  }
	  List_free(&rejected_hits);
	  List_free(&spliceends_antisense);

	  if (n_good_spliceends == 1) {
	    hits = List_push(hits,List_head(accepted_hits));
	    List_free(&accepted_hits);

	  } else {
	    /* 2.  Multiple hits, antisense, left1 (segmenti_left) */
	    debug7(printf("multiple splice hits, antisense, plus\n"));
	    donor_hits = acceptor_hits = (List_T) NULL;

	    /* plus branch from collect_elt_matches */
	    for (p = accepted_hits; p != NULL; p = List_next(p)) {
	      hit = (Stage3end_T) List_head(p);
	      donor = Stage3end_substring_donor(hit);
	      acceptor = Stage3end_substring_acceptor(hit);
	      if (Substring_genomicstart(donor) == segmenti_left) {
		donor_hits = List_push(donor_hits,(void *) hit);
	      } else if (Substring_genomicstart(acceptor) == segmenti_left) {
		acceptor_hits = List_push(acceptor_hits,(void *) hit);
	      } else {
		abort();
		Stage3end_free(&hit);
	      }
	    }

	    if (donor_hits != NULL) {
	      hitarray = (Stage3end_T *) List_to_array_n(&n,donor_hits);
	      qsort(hitarray,n,sizeof(Stage3end_T),donor_match_length_cmp);
	      i = 0;
	      while (i < n) {
		hit = hitarray[i];
		donor = Stage3end_substring_donor(hit);
		donor_length = Substring_match_length_orig(donor);
		j = i + 1;
		while (j < n && Substring_match_length_orig(Stage3end_substring_donor(hitarray[j])) == donor_length) {
		  j++;
		}
		if (j == i + 1) {
		  hits = List_push(hits,(void *) hit);
		} else {
#ifdef LARGE_GENOMES
		  ambcoords = (Uint8list_T) NULL;
#else
		  ambcoords = (Uintlist_T) NULL;
#endif
		  amb_knowni = (Intlist_T) NULL;
		  amb_nmismatches = (Intlist_T) NULL;
		  amb_probs = (Doublelist_T) NULL;
	      
		  for (k = i; k < j; k++) {
		    acceptor = Stage3end_substring_acceptor(hitarray[k]);
#ifdef LARGE_GENOMES
		    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_A(acceptor));
#else
		    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_A(acceptor));
#endif
		    amb_knowni = Intlist_push(amb_knowni,-1);
		    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		    amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
		  }
	      
		  nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
		  donor_prob = Junction_donor_prob(Stage3end_junctionA(hit));
		  prob = best_prob - donor_prob;
		  *ambiguous = List_push(*ambiguous,
					 (void *) Stage3end_new_splice(&(*found_score),
								       /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								       donor,/*acceptor*/NULL,donor_prob,/*acceptor_prob*/prob,/*distance*/0U,
								       /*shortdistancep*/false,/*penalty*/0,querylength,
								       /*ambcoords_donor*/NULL,ambcoords,
								       /*amb_knowni_donor*/NULL,amb_knowni,
								       /*amb_nmismatches_donort*/NULL,amb_nmismatches,
								       /*amb_probs_donor*/NULL,amb_probs,
								       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								       Stage3end_sensedir(hit),/*sarrayp*/false));
		  Doublelist_free(&amb_probs);
		  Intlist_free(&amb_nmismatches);
		  Intlist_free(&amb_knowni);
#ifdef LARGE_GENOMES
		  Uint8list_free(&ambcoords);
#else
		  Uintlist_free(&ambcoords);
#endif

		  for (k = i; k < j; k++) {
		    hit = hitarray[k];
		    Stage3end_free(&hit);
		  }
		}

		i = j;
	      }
	      FREE(hitarray);
	      List_free(&donor_hits);
	    }

	    if (acceptor_hits != NULL) {
	      hitarray = (Stage3end_T *) List_to_array_n(&n,acceptor_hits);
	      qsort(hitarray,n,sizeof(Stage3end_T),acceptor_match_length_cmp);
	      i = 0;
	      while (i < n) {
		hit = hitarray[i];
		acceptor = Stage3end_substring_acceptor(hit);
		acceptor_length = Substring_match_length_orig(acceptor);
		j = i + 1;
		while (j < n && Substring_match_length_orig(Stage3end_substring_acceptor(hitarray[j])) == acceptor_length) {
		  j++;
		}
		if (j == i + 1) {
		  hits = List_push(hits,(void *) hit);
		} else {
#ifdef LARGE_GENOMES
		  ambcoords = (Uint8list_T) NULL;
#else
		  ambcoords = (Uintlist_T) NULL;
#endif
		  amb_knowni = (Intlist_T) NULL;
		  amb_nmismatches = (Intlist_T) NULL;
		  amb_probs = (Doublelist_T) NULL;

		  for (k = i; k < j; k++) {
		    donor = Stage3end_substring_donor(hitarray[k]);
#ifdef LARGE_GENOMES
		    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_D(donor));
#else
		    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_D(donor));
#endif
		    amb_knowni = Intlist_push(amb_knowni,-1);
		    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		    amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
		  }
	    
		  nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
		  acceptor_prob = Junction_acceptor_prob(Stage3end_junctionD(hit));
		  prob = best_prob - acceptor_prob;
		  *ambiguous = List_push(*ambiguous,
					 (void *) Stage3end_new_splice(&(*found_score),
								       nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								       /*donor*/NULL,acceptor,/*donor_prob*/prob,acceptor_prob,/*distance*/0U,
								       /*shortdistancep*/false,/*penalty*/0,querylength,
								       ambcoords,/*ambcoords_acceptor*/NULL,
								       amb_knowni,/*amb_knowni_acceptor*/NULL,
								       amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								       amb_probs,/*amb_probs_acceptor*/NULL,
								       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								       Stage3end_sensedir(hit),/*sarrayp*/false));
		  Doublelist_free(&amb_probs);
		  Intlist_free(&amb_nmismatches);
		  Intlist_free(&amb_knowni);
#ifdef LARGE_GENOMES
		  Uint8list_free(&ambcoords);
#else
		  Uintlist_free(&ambcoords);
#endif

		  for (k = i; k < j; k++) {
		    hit = hitarray[k];
		    Stage3end_free(&hit);
		  }
		}

		i = j;
	      }
	      FREE(hitarray);
	      List_free(&acceptor_hits);
	    }

	    List_free(&accepted_hits);
	  }
	}

      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions_left);
    FREEA(mismatch_positions_right);
    FREEA(segmenti_donor_knownpos);
    FREEA(segmentj_acceptor_knownpos);
    FREEA(segmentj_antidonor_knownpos);
    FREEA(segmenti_antiacceptor_knownpos);
    FREEA(segmenti_donor_knowni);
    FREEA(segmentj_acceptor_knowni);
    FREEA(segmentj_antidonor_knowni);
    FREEA(segmenti_antiacceptor_knowni);
  } else {
    FREE(mismatch_positions_left);
    FREE(mismatch_positions_right);
    FREE(segmenti_donor_knownpos);
    FREE(segmentj_acceptor_knownpos);
    FREE(segmentj_antidonor_knownpos);
    FREE(segmenti_antiacceptor_knownpos);
    FREE(segmenti_donor_knowni);
    FREE(segmentj_acceptor_knowni);
    FREE(segmentj_antidonor_knowni);
    FREE(segmenti_antiacceptor_knowni);
  }
#else
  FREE(mismatch_positions_left);
  FREE(mismatch_positions_right);
  FREE(segmenti_donor_knownpos);
  FREE(segmentj_acceptor_knownpos);
  FREE(segmentj_antidonor_knownpos);
  FREE(segmenti_antiacceptor_knownpos);
  FREE(segmenti_donor_knowni);
  FREE(segmentj_acceptor_knowni);
  FREE(segmentj_antidonor_knowni);
  FREE(segmenti_antiacceptor_knowni);
#endif

  debug(printf("Finished find_singlesplices_plus with %d hits and %d lowprob\n",
	       List_length(hits),List_length(*lowprob)));

  return hits;
}


static List_T
find_singlesplices_minus (int *found_score, List_T hits, List_T *ambiguous, List_T *lowprob,
			  Segment_T *minus_spliceable, int minus_nspliceable,
			  Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress /* expecting rev */,
			  int splicing_penalty, int max_mismatches_allowed, bool first_read_p, int genestrand,
			  bool subs_or_indels_p) {
  int k, j, i, n;
  Segment_T segmenti, segmentj, segmentj_end, *ptr;
  Univcoord_T segmenti_left, segmentj_left;
  int nmismatches_left, nmismatches_right;
  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;
  int *mismatch_positions_left, *mismatch_positions_right;
  int *segmenti_donor_knownpos, *segmentj_acceptor_knownpos, *segmentj_antidonor_knownpos, *segmenti_antiacceptor_knownpos,
    *segmenti_donor_knowni, *segmentj_acceptor_knowni,
    *segmentj_antidonor_knowni, *segmenti_antiacceptor_knowni;
  
#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions_left = (int *) ALLOCA(querylength*sizeof(int));
    mismatch_positions_right = (int *) ALLOCA(querylength*sizeof(int));
    segmenti_donor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions_left = (int *) MALLOC(querylength*sizeof(int));
    mismatch_positions_right = (int *) MALLOC(querylength*sizeof(int));
    segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions_left = (int *) MALLOC(querylength*sizeof(int));
  mismatch_positions_right = (int *) MALLOC(querylength*sizeof(int));
  segmenti_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmentj_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segmenti_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  Chrpos_T max_distance;

  int floor_outer_i;
  int *floors_from_neg3, *floors_to_pos3;
  int nhits_local /*= 0*/;

  List_T accepted_hits, rejected_hits;
  List_T spliceends_sense, spliceends_antisense, p;
  List_T donor_hits, acceptor_hits;
  int donor_length, acceptor_length;
  Stage3end_T hit, *hitarray;
  int n_good_spliceends;
  int best_nmismatches, nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob, donor_prob, acceptor_prob;
  Substring_T donor, acceptor;

#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Intlist_T amb_knowni, amb_nmismatches;
  Doublelist_T amb_probs;


  debug4s(printf("*** Starting find_singlesplices_minus on %d spliceable segments ***\n",minus_nspliceable));
  /* debug(printf("Initially have %d hits\n",List_length(hits))); */

  if (floors != NULL) {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    for (ptr = minus_spliceable; ptr < &(minus_spliceable[minus_nspliceable]); ptr++) {
      segmenti = *ptr;
      debug4s(printf("minus_spliceable segmenti at diagonal %u\n",segmenti->diagonal));
      if (1 || segmenti->diagonal < (Univcoord_T) -1) { /* No markers were stored in spliceable */
	segmenti_left = segmenti->diagonal - querylength;
	floor_outer_i = floors_to_pos3[segmenti->querypos3];

	segmenti_antiacceptor_nknown = 0;
	segmenti_donor_nknown = 0;
	max_distance = shortsplicedist;

	if ((j = segmenti->splicesites_i) >= 0) {
	  /* Ends 4 and 5: mark known splice sites in segmenti */
	  while (j < nsplicesites && splicesites[j] < segmenti->diagonal) {
	    if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4s(printf("Setting known antiacceptor %d for segmenti at %llu\n",j,(unsigned long long) splicesites[j]));
	      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - segmenti_left;
	      segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	    } else if (splicetypes[j] == DONOR) {
	      debug4s(printf("Setting known donor %d for segmenti at %llu\n",j,(unsigned long long) splicesites[j]));
	      segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - segmenti_left;
	      segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	    }

	    /* This computation was already made in identify_all_segments */
	    if (splicedists[j] > max_distance) {
	      debug4s(printf("Setting max_distance for known %d to be %u\n",j,splicedists[j]));
	      max_distance = splicedists[j];
	    }

	    j++;
	  }
	}
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength;
	segmenti_donor_knownpos[segmenti_donor_nknown] = querylength;

	/* Identify potential segmentj for segmenti */
	segmentj_end = segmenti+1;
	while (
#ifdef NO_MARKER_SEGMENTS
	       segmentj_end < &(minus_segments[minus_nsegments]) && segmentj_end->chrnum == segmenti->chrnum &&
#endif
	       segmentj_end->diagonal <= segmenti->diagonal + max_distance) {
	  segmentj_end++;
	}


	spliceends_sense = spliceends_antisense = (List_T) NULL;

	if (segmentj_end - segmenti >= MAX_LOCALSPLICING_POTENTIAL) {
	  /* Too many to check */
	  /* segmentj_end = segmenti+1 + MAX_LOCALSPLICING_POTENTIAL; */
	  segmentj = segmentj_end; /* Don't process any */
	} else {
	  segmentj = segmenti+1;
	}
	for ( ; segmentj < segmentj_end; segmentj++) {
	  debug4s(printf("minus local?  diagonal %llu, querypos %d..%d => diagonal %llu, querypos %d..%d => ",
			 (unsigned long long) segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
			 (unsigned long long) segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	  /* j5 j3 i5 i3 */
	  assert(segmenti->diagonal < segmentj->diagonal);
	  if (segmentj->querypos3 >= segmenti->querypos5) {
	    /* Fail querypos test */
	    debug4s(printf("Bad querypos\n"));

	  } else if (segmenti->diagonal + min_intronlength > segmentj->diagonal) {
	    /* Too short to be an intron */
	    debug4s(printf("Too short\n"));

	  } else {
	    segmenti->right_splice_p = true;
	    segmentj->left_splice_p = true;
	    if (floors_from_neg3[segmentj->querypos5] + floor_outer_i > max_mismatches_allowed) {
	      /* Fail outer floor test */
	      /* floors->score[-index1interval][segmentj->querypos5] + floors->score[segmenti->querypos3][query_lastpos+index1interval] */;
	  
	      debug4s(printf("too many mismatches, outer floor = %d+%d=%d > %d\n",
			     floors->scorefrom[-index1interval][segmentj->querypos5],
			     floors->scorefrom[segmenti->querypos3][query_lastpos+index1interval],
			     floors->scorefrom[-index1interval][segmentj->querypos5] +
			     floors->scorefrom[segmenti->querypos3][query_lastpos+index1interval],
			     max_mismatches_allowed));

	    } else {
	      /* Apply leftmost/rightmost test */
	      if (segmenti->leftmost < 0) {
		nmismatches_left = Genome_mismatches_left(mismatch_positions_left,max_mismatches_allowed,
							  query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/querylength,
							  /*plusp*/false,genestrand);
		segmenti->leftmost = (nmismatches_left == 0) ? 0 : mismatch_positions_left[nmismatches_left-1];
		debug4s(printf("%d mismatches on left at:",nmismatches_left);
			for (i = 0; i <= nmismatches_left; i++) {
			  printf(" %d",mismatch_positions_left[i]);
			}
			printf("\n"));
	      }

	      segmentj_left = segmentj->diagonal - querylength;
	      if (segmentj->rightmost < 0) {
		debug4s(printf("left %u\n",segmentj_left));
		nmismatches_right = Genome_mismatches_right(mismatch_positions_right,max_mismatches_allowed,
							    query_compress,/*left*/segmentj_left,/*pos5*/0,/*pos3*/querylength,
							    /*plusp*/false,genestrand);
		segmentj->rightmost = (nmismatches_right == 0) ? 0 : mismatch_positions_right[nmismatches_right-1];
		debug4s(printf("%d mismatches on right [minus] at:",nmismatches_right);
			for (i = 0; i <= nmismatches_right; i++) {
			  printf(" %d",mismatch_positions_right[i]);
			}
			printf("\n"));
	      }

	      debug4s(printf("For a single splice, want leftmost %d > rightmost %d\n",segmenti->leftmost,segmentj->rightmost));

	      if (segmenti->leftmost > segmentj->rightmost) {
		/* Single splice is possible */

		segmentj_antidonor_nknown = 0;
		segmentj_acceptor_nknown = 0;
		if ((j = segmentj->splicesites_i) >= 0) {
		  /* Ends 3 and 6: mark known splice sites in segmentj */
		  while (j < nsplicesites && splicesites[j] < segmentj->diagonal) {
		    if (splicetypes[j] == ANTIDONOR) {
		      debug4s(printf("Setting known antidonor %d for segmentj at %llu\n",j,(unsigned long long) splicesites[j]));
		      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - segmentj_left;
		      segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
		    } else if (splicetypes[j] == ACCEPTOR) {
		      debug4s(printf("Setting known acceptor %d for segmentj at %llu\n",j,(unsigned long long) splicesites[j]));
		      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - segmentj_left;
		      segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
		    }
		    j++;
		  }
		}
		segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength;
		segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength;

		debug4s(printf("  => checking for single splice: Splice_solve_single_minus\n"));
		spliceends_sense =
		  Splice_solve_single_sense(&(*found_score),&nhits_local,spliceends_sense,&(*lowprob),
					    &segmenti->usedp,&segmentj->usedp,
					    /*segmenti_left*/segmenti->diagonal - querylength,
					    /*segmentj_left*/segmentj->diagonal - querylength,
					    segmenti->chrnum,segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
					    segmentj->chrnum,segmentj->chroffset,segmentj->chrhigh,segmentj->chrlength,
					    querylength,query_compress,					 
					    segmenti_donor_knownpos,segmentj_acceptor_knownpos,
					    segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					    segmenti_donor_knowni,segmentj_acceptor_knowni,
					    segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
					    segmenti_donor_nknown,segmentj_acceptor_nknown,
					    segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
					    splicing_penalty,max_mismatches_allowed,
					    /*plusp*/false,genestrand,first_read_p,subs_or_indels_p,
					    /*sarrayp*/false);
		spliceends_antisense =
		  Splice_solve_single_antisense(&(*found_score),&nhits_local,spliceends_antisense,&(*lowprob),
						&segmenti->usedp,&segmentj->usedp,
						/*segmenti_left*/segmenti->diagonal - querylength,
						/*segmentj_left*/segmentj->diagonal - querylength,
						segmenti->chrnum,segmenti->chroffset,segmenti->chrhigh,segmenti->chrlength,
						segmentj->chrnum,segmentj->chroffset,segmentj->chrhigh,segmentj->chrlength,
						querylength,query_compress,					 
						segmenti_donor_knownpos,segmentj_acceptor_knownpos,
						segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
						segmenti_donor_knowni,segmentj_acceptor_knowni,
						segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
						segmenti_donor_nknown,segmentj_acceptor_nknown,
						segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
						splicing_penalty,max_mismatches_allowed,
						/*plusp*/false,genestrand,first_read_p,subs_or_indels_p,
						/*sarrayp*/false);
	      }
	    }
	  }
	}

	/* Process results for segmenti, sense.  Modified from collect_elt_matches in sarray-read.c. */
	if (spliceends_sense != NULL) {
	  /* nmismatches here may be different for spliceends from Splice_solve, so pick based on prob and nmismatches */
	  best_nmismatches = querylength;
	  best_prob = 0.0;
	  for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			  Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
			  Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			  Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
			  Stage3end_nmismatches_whole(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			  Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	    if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	      best_nmismatches = nmismatches;
	    }
	    if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	      best_prob = prob;
	    }
	  }

	  n_good_spliceends = 0;
	  accepted_hits = rejected_hits = (List_T) NULL;
	  for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
		Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	      debug7(printf("accepting distance %d, probabilities %f and %f\n",
			    Stage3end_distance(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			    Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	      n_good_spliceends += 1;
	      accepted_hits = List_push(accepted_hits,(void *) hit);
	    } else {
	      rejected_hits = List_push(rejected_hits,(void *) hit);
	    }
	  }

	  if (n_good_spliceends == 0) {
	    /* Conjunction is too strict.  Allow for disjunction instead. */
	    List_free(&rejected_hits);
	    for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	      hit = (Stage3end_T) List_head(p);
	      if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
		  Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
		debug7(printf("accepting distance %d, probabilities %f and %f\n",
			      Stage3end_distance(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
		n_good_spliceends += 1;
		accepted_hits = List_push(accepted_hits,(void *) hit);
	      } else {
		rejected_hits = List_push(rejected_hits,(void *) hit);
	      }
	    }
	  }

	  for (p = rejected_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    Stage3end_free(&hit);
	  }
	  List_free(&rejected_hits);
	  List_free(&spliceends_sense);

	  if (n_good_spliceends == 1) {
	    hits = List_push(hits,List_head(accepted_hits));
	    List_free(&accepted_hits);

	  } else {
	    /* 1.  Multiple hits, sense, left1 (segmenti_left) */
	    debug7(printf("multiple splice hits, sense, minus\n"));
	    donor_hits = acceptor_hits = (List_T) NULL;

	    /* minus branch from collect_elt_matches */
	    for (p = accepted_hits; p != NULL; p = List_next(p)) {
	      hit = (Stage3end_T) List_head(p);
	      donor = Stage3end_substring_donor(hit);
	      acceptor = Stage3end_substring_acceptor(hit);
	      if (Substring_genomicend(donor) == segmenti_left) {
		donor_hits = List_push(donor_hits,(void *) hit);
	      } else if (Substring_genomicend(acceptor) == segmenti_left) {
		acceptor_hits = List_push(acceptor_hits,(void *) hit);
	      } else {
		abort();
		Stage3end_free(&hit);
	      }
	    }

	    if (donor_hits != NULL) {
	      hitarray = (Stage3end_T *) List_to_array_n(&n,donor_hits);
	      qsort(hitarray,n,sizeof(Stage3end_T),donor_match_length_cmp);
	      i = 0;
	      while (i < n) {
		hit = hitarray[i];
		donor = Stage3end_substring_donor(hit);
		donor_length = Substring_match_length_orig(donor);
		j = i + 1;
		while (j < n && Substring_match_length_orig(Stage3end_substring_donor(hitarray[j])) == donor_length) {
		  j++;
		}
		if (j == i + 1) {
		  hits = List_push(hits,(void *) hit);
		} else {
#ifdef LARGE_GENOMES
		  ambcoords = (Uint8list_T) NULL;
#else
		  ambcoords = (Uintlist_T) NULL;
#endif
		  amb_knowni = (Intlist_T) NULL;
		  amb_nmismatches = (Intlist_T) NULL;
		  amb_probs = (Doublelist_T) NULL;

		  for (k = i; k < j; k++) {
		    acceptor = Stage3end_substring_acceptor(hitarray[k]);
#ifdef LARGE_GENOMES
		    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_A(acceptor));
#else
		    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_A(acceptor));
#endif
		    amb_knowni = Intlist_push(amb_knowni,-1);
		    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		    amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
		  }

		  nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
		  donor_prob = Junction_donor_prob(Stage3end_junctionA(hit));
		  prob = best_prob - donor_prob;
		  *ambiguous = List_push(*ambiguous,
					 (void *) Stage3end_new_splice(&(*found_score),
								       /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								       donor,/*acceptor*/NULL,donor_prob,/*acceptor_prob*/prob,/*distance*/0U,
								       /*shortdistancep*/false,/*penalty*/0,querylength,
								       /*ambcoords_donor*/NULL,ambcoords,
								       /*amb_knowni_donor*/NULL,amb_knowni,
								       /*amb_nmismatches_donort*/NULL,amb_nmismatches,
								       /*amb_probs_donor*/NULL,amb_probs,
								       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								       Stage3end_sensedir(hit),/*sarrayp*/false));
		  Doublelist_free(&amb_probs);
		  Intlist_free(&amb_nmismatches);
		  Intlist_free(&amb_knowni);
#ifdef LARGE_GENOMES
		  Uint8list_free(&ambcoords);
#else
		  Uintlist_free(&ambcoords);
#endif

		  for (k = i; k < j; k++) {
		    hit = hitarray[k];
		    Stage3end_free(&hit);
		  }
		}

		i = j;
	      }
	      FREE(hitarray);
	      List_free(&donor_hits);
	    }

	    if (acceptor_hits != NULL) {
	      hitarray = (Stage3end_T *) List_to_array_n(&n,acceptor_hits);
	      qsort(hitarray,n,sizeof(Stage3end_T),acceptor_match_length_cmp);
	      i = 0;
	      while (i < n) {
		hit = hitarray[i];
		acceptor = Stage3end_substring_acceptor(hit);
		acceptor_length = Substring_match_length_orig(acceptor);
		j = i + 1;
		while (j < n && Substring_match_length_orig(Stage3end_substring_acceptor(hitarray[j])) == acceptor_length) {
		  j++;
		}
		if (j == i + 1) {
		  hits = List_push(hits,(void *) hit);
		} else {
#ifdef LARGE_GENOMES
		  ambcoords = (Uint8list_T) NULL;
#else
		  ambcoords = (Uintlist_T) NULL;
#endif
		  amb_knowni = (Intlist_T) NULL;
		  amb_nmismatches = (Intlist_T) NULL;
		  amb_probs = (Doublelist_T) NULL;

		  for (k = i; k < j; k++) {
		    donor = Stage3end_substring_donor(hitarray[k]);
#ifdef LARGE_GENOMES
		    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_D(donor));
#else
		    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_D(donor));
#endif
		    amb_knowni = Intlist_push(amb_knowni,-1);
		    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		    amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
		  }
	    
		  nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
		  acceptor_prob = Junction_acceptor_prob(Stage3end_junctionD(hit));
		  prob = best_prob - acceptor_prob;
		  *ambiguous = List_push(*ambiguous,
					 (void *) Stage3end_new_splice(&(*found_score),
								       nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								       /*donor*/NULL,acceptor,/*donor_prob*/prob,acceptor_prob,/*distance*/0U,
								       /*shortdistancep*/false,/*penalty*/0,querylength,
								       ambcoords,/*ambcoords_acceptor*/NULL,
								       amb_knowni,/*amb_knowni_acceptor*/NULL,
								       amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								       amb_probs,/*amb_probs_acceptor*/NULL,
								       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								       Stage3end_sensedir(hit),/*sarrayp*/false));
		  Doublelist_free(&amb_probs);
		  Intlist_free(&amb_nmismatches);
		  Intlist_free(&amb_knowni);
#ifdef LARGE_GENOMES
		  Uint8list_free(&ambcoords);
#else
		  Uintlist_free(&ambcoords);
#endif

		  for (k = i; k < j; k++) {
		    hit = hitarray[k];
		    Stage3end_free(&hit);
		  }
		}

		i = j;
	      }
	      FREE(hitarray);
	      List_free(&acceptor_hits);
	    }

	    List_free(&accepted_hits);
	  }
	}

	/* Process results for segmenti, antisense.  Modified from collect_elt_matches in sarray-read.c. */
	if (spliceends_antisense != NULL) {
	  /* nmismatches here may be different for spliceends from Splice_solve, so pick based on prob and nmismatches */
	  best_nmismatches = querylength;
	  best_prob = 0.0;
	  for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    debug7(printf("analyzing distance %d, donor length %d (%llu..%llu) and acceptor length %d (%llu..%llu), nmismatches %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			  Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
			  Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			  Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
			  Stage3end_nmismatches_whole(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
			  Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	    if ((nmismatches = Stage3end_nmismatches_whole(hit)) < best_nmismatches) {
	      best_nmismatches = nmismatches;
	    }
	    if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	      best_prob = prob;
	    }
	  }

	  n_good_spliceends = 0;
	  accepted_hits = rejected_hits = (List_T) NULL;
	  for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
		Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
	      debug7(printf("accepting distance %d, donor length %d and acceptor length %d, probabilities %f and %f\n",
			    Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			    Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			    Substring_siteD_prob(Stage3end_substring_donor(hit)),
			    Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	      n_good_spliceends += 1;
	      accepted_hits = List_push(accepted_hits,(void *) hit);
	    } else {
	      rejected_hits = List_push(rejected_hits,(void *) hit);
	    }
	  }

	  if (n_good_spliceends == 0) {
	    /* Conjunction is too strict.  Allow for disjunction instead. */
	    List_free(&rejected_hits);
	    for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	      hit = (Stage3end_T) List_head(p);
	      if (Stage3end_nmismatches_whole(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
		  Stage3end_chimera_prob(hit) >= best_prob - LOCALSPLICING_PROB_SLOP) {
		debug7(printf("accepting distance %d, donor length %d and acceptor length %d, probabilities %f and %f\n",
			      Stage3end_distance(hit),Substring_match_length_orig(Stage3end_substring_donor(hit)),
			      Substring_match_length_orig(Stage3end_substring_acceptor(hit)),
			      Substring_siteD_prob(Stage3end_substring_donor(hit)),
			      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
		n_good_spliceends += 1;
		accepted_hits = List_push(accepted_hits,(void *) hit);
	      } else {
		rejected_hits = List_push(rejected_hits,(void *) hit);
	      }
	    }
	  }

	  for (p = rejected_hits; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    Stage3end_free(&hit);
	  }
	  List_free(&rejected_hits);
	  List_free(&spliceends_antisense);

	  if (n_good_spliceends == 1) {
	    hits = List_push(hits,List_head(accepted_hits));
	    List_free(&accepted_hits);

	  } else {
	    /* 2.  Multiple hits, antisense, left1 (segmenti_left) */
	    debug7(printf("multiple splice hits, antisense, minus\n"));
	    donor_hits = acceptor_hits = (List_T) NULL;

	    /* minus branch from collect_elt_matches */
	    for (p = accepted_hits; p != NULL; p = List_next(p)) {
	      hit = (Stage3end_T) List_head(p);
	      donor = Stage3end_substring_donor(hit);
	      acceptor = Stage3end_substring_acceptor(hit);
	      if (Substring_genomicend(donor) == segmenti_left) {
		donor_hits = List_push(donor_hits,(void *) hit);
	      } else if (Substring_genomicend(acceptor) == segmenti_left) {
		acceptor_hits = List_push(acceptor_hits,(void *) hit);
	      } else {
		abort();
		Stage3end_free(&hit);
	      }
	    }

	    if (donor_hits != NULL) {
	      hitarray = (Stage3end_T *) List_to_array_n(&n,donor_hits);
	      qsort(hitarray,n,sizeof(Stage3end_T),donor_match_length_cmp);
	      i = 0;
	      while (i < n) {
		hit = hitarray[i];
		donor = Stage3end_substring_donor(hit);
		donor_length = Substring_match_length_orig(donor);
		j = i + 1;
		while (j < n && Substring_match_length_orig(Stage3end_substring_donor(hitarray[j])) == donor_length) {
		  j++;
		}
		if (j == i + 1) {
		  hits = List_push(hits,(void *) hit);
		} else {
#ifdef LARGE_GENOMES
		  ambcoords = (Uint8list_T) NULL;
#else
		  ambcoords = (Uintlist_T) NULL;
#endif
		  amb_knowni = (Intlist_T) NULL;
		  amb_nmismatches = (Intlist_T) NULL;
		  amb_probs = (Doublelist_T) NULL;
	      
		  for (k = i; k < j; k++) {
		    acceptor = Stage3end_substring_acceptor(hitarray[k]);
#ifdef LARGE_GENOMES
		    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_A(acceptor));
#else
		    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_A(acceptor));
#endif
		    amb_knowni = Intlist_push(amb_knowni,-1);
		    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
		    amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
		  }
	      
		  nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
		  donor_prob = Junction_donor_prob(Stage3end_junctionA(hit));
		  prob = best_prob - donor_prob;
		  *ambiguous = List_push(*ambiguous,
					 (void *) Stage3end_new_splice(&(*found_score),
								       /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
								       donor,/*acceptor*/NULL,donor_prob,/*acceptor_prob*/prob,/*distance*/0U,
								       /*shortdistancep*/false,/*penalty*/0,querylength,
								       /*ambcoords_donor*/NULL,ambcoords,
								       /*amb_knowni_donor*/NULL,amb_knowni,
								       /*amb_nmismatches_donort*/NULL,amb_nmismatches,
								       /*amb_probs_donor*/NULL,amb_probs,
								       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								       Stage3end_sensedir(hit),/*sarrayp*/false));
		  Doublelist_free(&amb_probs);
		  Intlist_free(&amb_nmismatches);
		  Intlist_free(&amb_knowni);
#ifdef LARGE_GENOMES
		  Uint8list_free(&ambcoords);
#else
		  Uintlist_free(&ambcoords);
#endif
		  for (k = i; k < j; k++) {
		    hit = hitarray[k];
		    Stage3end_free(&hit);
		  }
		}

		i = j;
	      }
	      FREE(hitarray);
	      List_free(&donor_hits);
	    }

	    if (acceptor_hits != NULL) {
	      hitarray = (Stage3end_T *) List_to_array_n(&n,acceptor_hits);
	      qsort(hitarray,n,sizeof(Stage3end_T),acceptor_match_length_cmp);
	      i = 0;
	      while (i < n) {
		hit = hitarray[i];
		acceptor = Stage3end_substring_acceptor(hit);
		acceptor_length = Substring_match_length_orig(acceptor);
		j = i + 1;
		while (j < n && Substring_match_length_orig(Stage3end_substring_acceptor(hitarray[j])) == acceptor_length) {
		  j++;
		}
		if (j == i + 1) {
		  hits = List_push(hits,(void *) hit);
		} else {
#ifdef LARGE_GENOMES
		  ambcoords = (Uint8list_T) NULL;
#else
		  ambcoords = (Uintlist_T) NULL;
#endif
		  amb_knowni = (Intlist_T) NULL;
		  amb_nmismatches = (Intlist_T) NULL;
		  amb_probs = (Doublelist_T) NULL;

		  for (k = i; k < j; k++) {
		    donor = Stage3end_substring_donor(hitarray[k]);
#ifdef LARGE_GENOMES
		    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord_D(donor));
#else
		    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord_D(donor));
#endif
		    amb_knowni = Intlist_push(amb_knowni,-1);
		    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
		    amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
		  }
	    
		  nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
		  acceptor_prob = Junction_acceptor_prob(Stage3end_junctionD(hit));
		  prob = best_prob - acceptor_prob;
		  *ambiguous = List_push(*ambiguous,
					 (void *) Stage3end_new_splice(&(*found_score),
								       nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
								       /*donor*/NULL,acceptor,/*donor_prob*/prob,acceptor_prob,/*distance*/0U,
								       /*shortdistancep*/false,/*penalty*/0,querylength,
								       ambcoords,/*ambcoords_acceptor*/NULL,
								       amb_knowni,/*amb_knowni_acceptor*/NULL,
								       amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								       amb_probs,/*amb_probs_acceptor*/NULL,
								       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								       Stage3end_sensedir(hit),/*sarrayp*/false));
		  Doublelist_free(&amb_probs);
		  Intlist_free(&amb_nmismatches);
		  Intlist_free(&amb_knowni);
#ifdef LARGE_GENOMES
		  Uint8list_free(&ambcoords);
#else
		  Uintlist_free(&ambcoords);
#endif

		  for (k = i; k < j; k++) {
		    hit = hitarray[k];
		    Stage3end_free(&hit);
		  }
		}

		i = j;
	      }
	      FREE(hitarray);
	      List_free(&acceptor_hits);
	    }

	    List_free(&accepted_hits);
	  }
	}

      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions_left);
    FREEA(mismatch_positions_right);
    FREEA(segmenti_donor_knownpos);
    FREEA(segmentj_acceptor_knownpos);
    FREEA(segmentj_antidonor_knownpos);
    FREEA(segmenti_antiacceptor_knownpos);
    FREEA(segmenti_donor_knowni);
    FREEA(segmentj_acceptor_knowni);
    FREEA(segmentj_antidonor_knowni);
    FREEA(segmenti_antiacceptor_knowni);
  } else {
    FREE(mismatch_positions_left);
    FREE(mismatch_positions_right);
    FREE(segmenti_donor_knownpos);
    FREE(segmentj_acceptor_knownpos);
    FREE(segmentj_antidonor_knownpos);
    FREE(segmenti_antiacceptor_knownpos);
    FREE(segmenti_donor_knowni);
    FREE(segmentj_acceptor_knowni);
    FREE(segmentj_antidonor_knowni);
    FREE(segmenti_antiacceptor_knowni);
  }
#else
  FREE(mismatch_positions_left);
  FREE(mismatch_positions_right);
  FREE(segmenti_donor_knownpos);
  FREE(segmentj_acceptor_knownpos);
  FREE(segmentj_antidonor_knownpos);
  FREE(segmenti_antiacceptor_knownpos);
  FREE(segmenti_donor_knowni);
  FREE(segmentj_acceptor_knowni);
  FREE(segmentj_antidonor_knowni);
  FREE(segmenti_antiacceptor_knowni);
#endif

  debug(printf("Finished find_singlesplices_minus with %d hits and %d lowprob\n",
	       List_length(hits),List_length(*lowprob)));

  return hits;
}


#ifdef LARGE_GENOMES
static Uint8list_T
lookup_splicesites (Doublelist_T *probs_list, Intlist_T splicesites_i, Univcoord_T *splicesites) {
  Uint8list_T coords = NULL;
  Intlist_T p;

  *probs_list = (Doublelist_T) NULL;
  for (p = splicesites_i; p != NULL; p = Intlist_next(p)) {
    coords = Uint8list_push(coords,splicesites[Intlist_head(p)]);
    *probs_list = Doublelist_push(*probs_list,2.0);
  }

  return Uint8list_reverse(coords);
}
#else
static Uintlist_T
lookup_splicesites (Doublelist_T *probs_list, Intlist_T splicesites_i, Univcoord_T *splicesites) {
  Uintlist_T coords = NULL;
  Intlist_T p;

  *probs_list = (Doublelist_T) NULL;
  for (p = splicesites_i; p != NULL; p = Intlist_next(p)) {
    coords = Uintlist_push(coords,splicesites[Intlist_head(p)]);
    *probs_list = Doublelist_push(*probs_list,2.0);
  }

  return Uintlist_reverse(coords);
}
#endif


#if 0
static int
substringD_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length_orig(Stage3end_substringD(x));
  int y_length = Substring_match_length_orig(Stage3end_substringD(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    x_length = Substring_match_length_orig(Stage3end_substringA(x));
    y_length = Substring_match_length_orig(Stage3end_substringA(y));
    if (x_length < y_length) {
      return -1;
    } else if (y_length < x_length) {
      return +1;
    } else {
      return 0;
    }
  }
}
#endif

#if 0
static int
substringA_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length_orig(Stage3end_substringA(x));
  int y_length = Substring_match_length_orig(Stage3end_substringA(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    x_length = Substring_match_length_orig(Stage3end_substringD(x));
    y_length = Substring_match_length_orig(Stage3end_substringD(y));
    if (x_length < y_length) {
      return -1;
    } else if (y_length < x_length) {
      return +1;
    } else {
      return 0;
    }
  }
}
#endif



static void
find_spliceends_shortend (List_T **shortend_donors, List_T **shortend_antidonors,
			  List_T **shortend_acceptors, List_T **shortend_antiacceptors,
			  Segment_T *anchor_segments, int nanchors,
#ifdef DEBUG4E
			  char *queryptr,
#endif
			  Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress,
			  int max_mismatches_allowed, bool plusp, int genestrand) {
#ifdef DEBUG4E
  char *gbuffer;
#endif

  Segment_T segment, *p;
  Substring_T hit;
  Univcoord_T segment_left;
  int nmismatches, jstart, jend, j;
  int splice_pos;

  int nmismatches_left, nmismatches_right;
  int *floors_from_neg3, *floors_to_pos3;
  int sensedir;

  int splice_pos_start, splice_pos_end;
#ifdef DEBUG4E
  int i;
#endif

  int *mismatch_positions;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug4e(printf("Entering find_spliceends_shortend with %d anchor segments\n",nanchors));

  if (floors != NULL) {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    for (p = &(anchor_segments[0]); p < &(anchor_segments[nanchors]); p++) {
      segment = *p;
      assert(segment->diagonal != (Univcoord_T) -1);
      if (segment->splicesites_i >= 0) {
	segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
	debug4e(printf("find_spliceends_shortend: Checking up to %d mismatches at diagonal %llu (querypos %d..%d) - querylength %d = %llu, floors %d and %d\n",
		       max_mismatches_allowed,(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,querylength,(unsigned long long) segment_left,
		       floors_from_neg3[segment->querypos5],floors_to_pos3[segment->querypos3]));

	debug4e(
		gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		printf("genome 0..: %s\n",gbuffer);
		printf("query  0..: %s\n",queryptr);
		FREE(gbuffer);
		);

	/* Splice ends from left to splice site */
	if ((plusp == true && floors_from_neg3[segment->querypos5] <= max_mismatches_allowed) ||
	    (plusp == false && floors_to_pos3[segment->querypos3] <= max_mismatches_allowed)) {

	  /* pos3 was trimpos */
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						    plusp,genestrand);

	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_start = 1;  /* not index1part */
	  if (nmismatches_left <= max_mismatches_allowed) {
	    splice_pos_end = querylength - 1;
	  } else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1) {
	    splice_pos_end = querylength - 1;
	  }

	  debug4e(printf("Search for splice sites from %d up (%llu) to %d (%llu)\n",
			 splice_pos_start,(unsigned long long) segment_left+splice_pos_start,
			 splice_pos_end,(unsigned long long) segment_left+splice_pos_end));

	  jstart = segment->splicesites_i;
	  while (jstart < nsplicesites && splicesites[jstart] < segment_left + splice_pos_start) {
	    jstart++;
	  }
	  jend = jstart;
	  while (jend < nsplicesites && splicesites[jend] <= segment_left + splice_pos_end) { /* Needs to be <= */
	    jend++;
	  }

	  nmismatches = 0;
	  for (j = jstart; j < jend; j++) {
	    debug4e(printf("splicesites_i #%d is at %llu\n",j,(unsigned long long) splicesites[j]));
	    splice_pos = splicesites[j] - segment_left;
	    while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos) { /* Changed from <= to < */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }
#if 0
	    assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/0,/*pos3*/splice_pos,
								    plusp,genestrand,first_read_p));
#endif
	    if (nmismatches > max_mismatches_allowed) {
	      debug4e(printf("nmismatches %d > max_mismatches_allowed %d\n",nmismatches,max_mismatches_allowed));
	    } else if (splicetypes[j] == DONOR) {
	      debug4e(printf("Known donor #%d at querypos %d\n",j,splicesites[j] - segment_left));
	      debug4e(printf("Known donor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     segment_left,(unsigned long long) splice_pos,nmismatches,splice_pos_end));
	      sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	      if ((hit = Substring_new_donor(/*donor_coord*/splicesites[j],/*donor_knowni*/j,
					     splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
					     nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
					     querylength,plusp,genestrand,
					     sensedir,segment->chrnum,segment->chroffset,
					     segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s donor: known at %d (%d mismatches)\n",
			       plusp == true ? "plus" : "minus",Substring_siteD_pos(hit),nmismatches));
		(*shortend_donors)[nmismatches] = List_push((*shortend_donors)[nmismatches],(void *) hit);
	      }

	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4e(printf("Known antiacceptor #%d at querypos %d\n",j,splicesites[j] - segment_left));
	      debug4e(printf("Known antiacceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     segment_left,(unsigned long long) splice_pos,nmismatches,splice_pos_end));
	      sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	      if ((hit = Substring_new_acceptor(/*acceptor_coord*/splicesites[j],/*acceptor_knowni*/j,
						splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
						nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
						querylength,plusp,genestrand,
						sensedir,segment->chrnum,segment->chroffset,
						segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s antiacceptor : known at %d (%d mismatches)\n",
			       plusp == true ? "plus" : "minus",Substring_siteA_pos(hit),nmismatches));
		(*shortend_antiacceptors)[nmismatches] = List_push((*shortend_antiacceptors)[nmismatches],(void *) hit);
	      }
	    }
	  }
	}

	/* Splice ends from splice site to right end */
	if ((plusp == true && floors_to_pos3[segment->querypos3] <= max_mismatches_allowed) ||
	    (plusp == false && floors_from_neg3[segment->querypos5] <= max_mismatches_allowed)) {

	  /* pos5 was trimpos+1 */
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      query_compress,/*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						      plusp,genestrand);

	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_end = querylength - 1;  /* not query_lastpos */
	  if (nmismatches_right <= max_mismatches_allowed) {
	    splice_pos_start = 1;
	  } else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    splice_pos_start = 1;
	  }
	  /* splice_pos_start = 1; */

	  debug4e(printf("Search for splice sites from %d (%llu) down to %d (%llu)\n",
			 splice_pos_end,(unsigned long long) segment_left+splice_pos_end,
			 splice_pos_start,(unsigned long long) segment_left+splice_pos_start));

	  jstart = segment->splicesites_i;
	  while (jstart < nsplicesites && splicesites[jstart] < segment_left + splice_pos_start) {
	    jstart++;
	  }
	  jend = jstart;
	  while (jend < nsplicesites && splicesites[jend] <= segment_left + splice_pos_end) { /* Needs to be <= */
	    jend++;
	  }

	  nmismatches = 0;
	  for (j = jend - 1; j >= jstart; j--) {
	    debug4e(printf("splicesites_i #%d is at %llu\n",j,(unsigned long long) splicesites[j]));
	    splice_pos = splicesites[j] - segment_left;
	    while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }
#if 0
	    assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
								    plusp,genestrand,first_read_p));
#endif
	    if (nmismatches > max_mismatches_allowed) {
	      debug4e(printf("nmismatches %d > max_mismatches_allowed %d\n",nmismatches,max_mismatches_allowed));
	    } else if (splicetypes[j] == ACCEPTOR) {
	      debug4e(printf("Known acceptor #%d at querypos %d\n",j,splicesites[j] - segment_left));
	      debug4e(printf("Known acceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
	      sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	      if ((hit = Substring_new_acceptor(/*acceptor_coord*/splicesites[j],/*acceptor_knowni*/j,
						splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
						nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
						querylength,plusp,genestrand,
						sensedir,segment->chrnum,segment->chroffset,
						segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s acceptor: known at %d (%d mismatches)\n",
			       plusp == true ? "plus" : "minus",Substring_siteA_pos(hit),nmismatches));
		(*shortend_acceptors)[nmismatches] = List_push((*shortend_acceptors)[nmismatches],(void *) hit);
	      }

	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4e(printf("Known antidonor #%d at querypos %d\n",j,splicesites[j] - segment_left));
	      debug4e(printf("Known antidonor for segmenti at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
	      sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	      if ((hit = Substring_new_donor(/*donor_coord*/splicesites[j],/*donor_knowni*/j,
					     splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
					     nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
					     querylength,plusp,genestrand,
					     sensedir,segment->chrnum,segment->chroffset,
					     segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s antidonor: known at %d (%d mismatches)\n",
			       plusp == true ? "plus" : "minus",Substring_siteD_pos(hit),nmismatches));
		(*shortend_antidonors)[nmismatches] = List_push((*shortend_antidonors)[nmismatches],(void *) hit);
	      }
	    }
	  }
	}
      }
    }
  }
    
#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
  } else {
    FREE(mismatch_positions);
  }
#else
  FREE(mismatch_positions);
#endif
  
  return;
}


static void
find_spliceends_distant_dna_plus (List_T **distant_startfrags, List_T **distant_endfrags,
				  Segment_T *anchor_segments, int nanchors,
#ifdef DEBUG4E
				  char *queryptr,
#endif
				  Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress,
				  int max_mismatches_allowed, int genestrand) {
#ifdef DEBUG4E
  char *gbuffer;
  int i;
#endif

  Segment_T segment, *p;
  Substring_T hit;
  Univcoord_T segment_left;
  int nmismatches;
  int splice_pos;

  int nmismatches_left, nmismatches_right;
  int *floors_from_neg3, *floors_to_pos3;

  int splice_pos_start, splice_pos_end;
  int *mismatch_positions;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug4e(printf("Entering find_spliceends_distant_dna with %d anchor segments\n",nanchors));

  if (floors != NULL) {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    for (p = &(anchor_segments[0]); p < &(anchor_segments[nanchors]); p++) {
      segment = *p;
      assert(segment->diagonal != (Univcoord_T) -1);

      segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      debug4e(printf("find_spliceends_distant_dna: Checking up to %d mismatches at diagonal %llu (querypos %d..%d) - querylength %d = %llu, floors %d and %d\n",
		     max_mismatches_allowed,(unsigned long long) segment->diagonal,
		     segment->querypos5,segment->querypos3,querylength,(unsigned long long) segment_left,
		     floors_from_neg3[segment->querypos5],floors_to_pos3[segment->querypos3]));
    
      debug4e(
	      gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
	      Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
	      printf("genome 0..: %s\n",gbuffer);
	      printf("query  0..: %s\n",queryptr);
	      FREE(gbuffer);
	      );

      /* Splice ends from left to splice site */
      if (floors_from_neg3[segment->querypos5] <= max_mismatches_allowed) {

	/* pos3 was trimpos */
	nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						  query_compress,/*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						  /*plusp*/true,genestrand);

	debug4e(
		printf("%d mismatches on left (%d allowed) at:",
		       nmismatches_left,max_mismatches_allowed);
		for (i = 0; i <= nmismatches_left; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);

	splice_pos_start = index1part;
	if (nmismatches_left <= max_mismatches_allowed) {
	  splice_pos_end = querylength - 1;
	} else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1) {
	  splice_pos_end = querylength - 1;
	}

	debug4e(printf("Allow all splice points from %d up to %d\n",splice_pos_start,splice_pos_end));

	nmismatches = 0;
	while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos_start) {
	  debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	  nmismatches++;
	}

	splice_pos = splice_pos_start;
	while (splice_pos <= splice_pos_end && nmismatches <= max_mismatches_allowed) {
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/0,/*pos3*/splice_pos,
								  /*plusp*/true,genestrand));
	  if ((hit = Substring_new_startfrag(/*startfrag_coord*/segment_left + splice_pos,
					     splice_pos,nmismatches,/*left*/segment_left,query_compress,
					     querylength,/*plusp*/true,genestrand,
					     segment->chrnum,segment->chroffset,
					     segment->chrhigh,segment->chrlength)) != NULL) {
	    debug4e(printf("=> plus startfrag: at %d (%d mismatches)\n",Substring_siteN_pos(hit),nmismatches));
	    debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	    (*distant_startfrags)[nmismatches] = List_push((*distant_startfrags)[nmismatches],(void *) hit);
	  }

	  /* use splice_pos in the above loop because splice_pos defines a right substring boundary */
	  if (splice_pos++ == mismatch_positions[nmismatches]) {
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	}
      }

      /* Splice ends from splice site to right end */
      if (floors_to_pos3[segment->querypos3] <= max_mismatches_allowed) {

	/* pos5 was trimpos+1 */
	nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						    /*plusp*/true,genestrand);

	debug4e(
		printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		for (i = 0; i <= nmismatches_right; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);

	splice_pos_end = query_lastpos;
	if (nmismatches_right <= max_mismatches_allowed) {
	  splice_pos_start = 1;
	} else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	  splice_pos_start = 1;
	}

	debug4e(printf("Allow all splice sites from %d down to %d\n",splice_pos_end,splice_pos_start));

	nmismatches = 0;
	while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos_end) {
	  debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	  nmismatches++;
	}

	splice_pos = splice_pos_end;
	while (splice_pos >= splice_pos_start && nmismatches <= max_mismatches_allowed) {
	  debug4e(printf(" splice pos %d, nmismatches (quick) %d, nmismatches (goldstd) %d\n",splice_pos,nmismatches,
			 Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
							   /*plusp*/true,genestrand)));
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
								  /*plusp*/true,genestrand));
	  if ((hit = Substring_new_endfrag(/*endfrag_coord*/segment_left + splice_pos,
					   splice_pos,nmismatches,/*left*/segment_left,query_compress,
					   querylength,/*plusp*/true,genestrand,
					   segment->chrnum,segment->chroffset,
					   segment->chrhigh,segment->chrlength)) != NULL) {
	    debug4e(printf("=> plus endfrag: at %d (%d mismatches)\n",Substring_siteN_pos(hit),nmismatches));
	    debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	    (*distant_endfrags)[nmismatches] = List_push((*distant_endfrags)[nmismatches],(void *) hit);
	  }

	  /* use splice_pos for the next loop because splice_pos defines a left substring boundary */
	  if (--splice_pos == mismatch_positions[nmismatches]) {
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }

	}
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
  } else {
    FREE(mismatch_positions);
  }
#else
  FREE(mismatch_positions);
#endif

  return;
}


static void
find_spliceends_distant_dna_minus (List_T **distant_startfrags, List_T **distant_endfrags,
				   Segment_T *anchor_segments, int nanchors,
#ifdef DEBUG4E
				   char *queryptr,
#endif
				   Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress,
				   int max_mismatches_allowed, int genestrand) {
#ifdef DEBUG4E
  char *gbuffer;
  int i;
#endif

  Segment_T segment, *p;
  Substring_T hit;
  Univcoord_T segment_left;
  int nmismatches;
  int splice_pos;

  int nmismatches_left, nmismatches_right;
  int *floors_from_neg3, *floors_to_pos3;

  int splice_pos_start, splice_pos_end;
  int *mismatch_positions;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


  debug4e(printf("Entering find_spliceends_distant_dna with %d anchor segments\n",nanchors));

  if (floors != NULL) {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    for (p = &(anchor_segments[0]); p < &(anchor_segments[nanchors]); p++) {
      segment = *p;
      assert(segment->diagonal != (Univcoord_T) -1);

      segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      debug4e(printf("find_spliceends_distant_dna: Checking up to %d mismatches at diagonal %llu (querypos %d..%d) - querylength %d = %llu, floors %d and %d\n",
		     max_mismatches_allowed,(unsigned long long) segment->diagonal,
		     segment->querypos5,segment->querypos3,querylength,(unsigned long long) segment_left,
		     floors_from_neg3[segment->querypos5],floors_to_pos3[segment->querypos3]));

      debug4e(
	      gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
	      Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
	      printf("genome 0..: %s\n",gbuffer);
	      printf("query  0..: %s\n",queryptr);
	      FREE(gbuffer);
	      );

      /* Splice ends from left to splice site */
      if (floors_to_pos3[segment->querypos3] <= max_mismatches_allowed) {

	/* pos3 was trimpos */
	nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						  query_compress,/*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						  /*plusp*/false,genestrand);

	debug4e(
		printf("%d mismatches on left (%d allowed) at:",
		       nmismatches_left,max_mismatches_allowed);
		for (i = 0; i <= nmismatches_left; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);

	splice_pos_start = index1part;
	if (nmismatches_left <= max_mismatches_allowed) {
	  splice_pos_end = querylength - 1;
	} else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1) {
	  splice_pos_end = querylength - 1;
	}

	debug4e(printf("Allow all splice points from %d up to %d\n",splice_pos_start,splice_pos_end));

	nmismatches = 0;
	while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos_start) {
	  debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	  nmismatches++;
	}

	splice_pos = splice_pos_start;
	while (splice_pos <= splice_pos_end && nmismatches <= max_mismatches_allowed) {
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/0,/*pos3*/splice_pos,
								  /*plusp*/false,genestrand));
	  if ((hit = Substring_new_endfrag(/*endfrag_coord*/segment_left + splice_pos,
					   splice_pos,nmismatches,/*left*/segment_left,query_compress,
					   querylength,/*plusp*/false,genestrand,
					   segment->chrnum,segment->chroffset,
					   segment->chrhigh,segment->chrlength)) != NULL) {
	    debug4e(printf("=> minus endfrag: at %d (%d mismatches)\n",Substring_siteN_pos(hit),nmismatches));
	    debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	    (*distant_endfrags)[nmismatches] = List_push((*distant_endfrags)[nmismatches],(void *) hit);
	  }

	  /* use splice_pos in the above loop because splice_pos defines a right substring boundary */
	  if (splice_pos++ == mismatch_positions[nmismatches]) {
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	}
      }

      /* Splice ends from splice site to right end */
      if (floors_from_neg3[segment->querypos5] <= max_mismatches_allowed) {

	/* pos5 was trimpos+1 */
	nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						    /*plusp*/false,genestrand);

	debug4e(
		printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		for (i = 0; i <= nmismatches_right; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);

	splice_pos_end = query_lastpos;
	if (nmismatches_right <= max_mismatches_allowed) {
	  splice_pos_start = 1;
	} else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	  splice_pos_start = 1;
	}

	debug4e(printf("Allow all splice sites from %d down to %d\n",splice_pos_end,splice_pos_start));

	nmismatches = 0;
	while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos_end) {
	  debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	  nmismatches++;
	}

	splice_pos = splice_pos_end;
	while (splice_pos >= splice_pos_start && nmismatches <= max_mismatches_allowed) {
	  debug4e(printf(" splice pos %d, nmismatches (quick) %d, nmismatches (goldstd) %d\n",splice_pos,nmismatches,
			 Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
							   /*plusp*/false,genestrand)));
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
								  /*plusp*/false,genestrand));
	  if ((hit = Substring_new_startfrag(/*startfrag_coord*/segment_left + splice_pos,
					     splice_pos,nmismatches,/*left*/segment_left,query_compress,
					     querylength,/*plusp*/false,genestrand,
					     segment->chrnum,segment->chroffset,
					     segment->chrhigh,segment->chrlength)) != NULL) {
	    debug4e(printf("=> minus startfrag: at %d (%d mismatches)\n",Substring_siteN_pos(hit),nmismatches));
	    debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	    (*distant_startfrags)[nmismatches] = List_push((*distant_startfrags)[nmismatches],(void *) hit);
	  }

	  /* use splice_pos for the next loop because splice_pos defines a left substring boundary */
	  if (--splice_pos == mismatch_positions[nmismatches]) {
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	}
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
  } else {
    FREE(mismatch_positions);
  }
#else
  FREE(mismatch_positions);
#endif

  return;
}


#if 0
/* Using a streamlined version instead */
/* Produces lists of distant_donors and distant_acceptors that are substrings */
/* TODO: Change to lists of Stage3end_T objects, including GMAP.
   Change definition of a chimera to be two Stage3end_T objects, instead
   of two substrings. */
static void
find_spliceends_distant_rna (List_T **distant_donors, List_T **distant_antidonors,
			     List_T **distant_acceptors, List_T **distant_antiacceptors,
			     Segment_T *anchor_segments, int nanchors,
#ifdef DEBUG4E
			     char *queryptr,
#endif
			     Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress,
			     int max_mismatches_allowed, bool plusp, int genestrand) {
#ifdef DEBUG4E
  char *gbuffer;
#endif

  Segment_T segment, *p;
  Substring_T hit;
  Univcoord_T segment_left;
  int nmismatches, j, i;
  int splice_pos, last_querypos;
  double prob;

  int nmismatches_left, nmismatches_right;
  int *floors_from_neg3, *floors_to_pos3;
  int sensedir;

  int splice_pos_start, splice_pos_end;
  int segment_donor_nknown, segment_acceptor_nknown, segment_antidonor_nknown, segment_antiacceptor_nknown;
  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;

  int *mismatch_positions;
  int *segment_donor_knownpos, *segment_acceptor_knownpos, *segment_antidonor_knownpos, *segment_antiacceptor_knownpos,
    *segment_donor_knowni, *segment_acceptor_knowni, *segment_antidonor_knowni, *segment_antiacceptor_knowni;
  int *positions_alloc, *knowni_alloc;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_donor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_acceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_antidonor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_antiacceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_donor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_acceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_antidonor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_antiacceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug4e(printf("Entering find_spliceends_distant_rna with %d anchor segments\n",nanchors));

  if (floors != NULL) {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    for (p = &(anchor_segments[0]); p < &(anchor_segments[nanchors]); p++) {
      segment = *p;
      assert(segment->diagonal != (Univcoord_T) -1);

      segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      if ((first_querypos = segment->querypos5 - (index1interval - 1)) < 0) {
	first_querypos = 0;
      }
      if ((last_querypos = segment->querypos3 + index1part + (index1interval - 1)) > querylength) {
	last_querypos = querylength;
      }

      debug4e(printf("find_spliceends_distant_rna: Checking up to %d mismatches at diagonal %llu (querypos %d..%d) - querylength %d = %llu, floors %d and %d, plusp %d\n",
		     max_mismatches_allowed,(unsigned long long) segment->diagonal,
		     segment->querypos5,segment->querypos3,querylength,(unsigned long long) segment_left,
		     floors_from_neg3[segment->querypos5],floors_to_pos3[segment->querypos3],
		     plusp));

      debug4e(
	      gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
	      Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
	      printf("genome 0..: %s\n",gbuffer);
	      printf("query  0..: %s\n",queryptr);
	      FREE(gbuffer);
	      );

      /* Find splices on genomic right */
      if (plusp) {
	/* ? require that floors_from_neg3[segment->querypos5] <= max_mismatches_allowed */
	if (first_querypos < index1part && last_querypos < query_lastpos) {
	  /* genomic left anchor */
	  debug4e(printf("Searching genomic right: plus genomic left anchor\n"));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,
#if 0
						    /*pos5*/0,/*pos3*/querylength,
#else
						    /*pos5 (was 0)*/first_querypos,/*pos3*/querylength,
#endif
						    plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

#if 0
	  splice_pos_start = index1part;
#else
	  splice_pos_start = first_querypos;
#endif
	  if (nmismatches_left <= max_mismatches_allowed) {
	    splice_pos_end = querylength - 1;
	  } else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1) {
	    splice_pos_end = querylength - 1;
	  }
	  
	} else if (first_querypos > index1part && last_querypos > query_lastpos) {
	  /* genomic right anchor.  No need to find splices on genomic right */
	  debug4e(printf("Searching genomic right: plus genomic right anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	  
	} else if (first_querypos > index1part && last_querypos < query_lastpos &&
		   segment->spliceable_low_p == true) {
	  /* middle anchor */
	  debug4e(printf("Searching genomic right: plus middle anchor\n"));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,
						    /*pos5*/first_querypos,/*pos3*/querylength,
						    plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_start = first_querypos;
	  if (nmismatches_left <= max_mismatches_allowed) {
	    splice_pos_end = querylength - 1;
	  } else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1) {
	    splice_pos_end = querylength - 1;
	  }

	} else {
	  /* full anchor */
	  debug4e(printf("Searching genomic right: plus full anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	}

      } else {
	/* ? require that floors_to_pos3[segment->querypos3] <= max_mismatches_allowed */
	if (first_querypos < index1part && last_querypos < query_lastpos) {
	  /* genomic right anchor.  No need to find splices on genomic right */
	  debug4e(printf("Searching genomic right: minus genomic right anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;

	} else if (first_querypos > index1part && last_querypos > query_lastpos) {
	  /* genomic left anchor */
	  debug4e(printf("Searching genomic right: minus genomic left anchor\n"));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,
#if 0
						    /*pos5*/0,/*pos3*/querylength,
#else
						    /*pos5 (was 0)*/querylength - last_querypos,
						    /*pos3*/querylength,
#endif
						    plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

#if 0
	  splice_pos_start = index1part;
#else
	  splice_pos_start = querylength - last_querypos;
#endif
	  if (nmismatches_left <= max_mismatches_allowed) {
	    splice_pos_end = querylength - 1;
	  } else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1) {
	    splice_pos_end = querylength - 1;
	  }

	} else if (first_querypos > index1part && last_querypos < query_lastpos &&
		   segment->spliceable_low_p == true) {
	  /* middle anchor */
	  debug4e(printf("Searching genomic right: minus middle anchor\n"));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,
						    /*pos5*/querylength - last_querypos,
						    /*pos3*/querylength,plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_start = querylength - last_querypos;
	  if (nmismatches_left <= max_mismatches_allowed) {
	    splice_pos_end = querylength - 1;
	  } else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1) {
	    splice_pos_end = querylength - 1;
	  }

	} else {
	  /* full anchor */
	  debug4e(printf("Searching genomic right: minus full anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	}
      }

      if (splice_pos_start <= splice_pos_end) {
	debug4e(printf("Search for splice sites from %d up to %d\n",splice_pos_start,splice_pos_end));

	segment_donor_nknown = 0;
	segment_antiacceptor_nknown = 0;
	if ((j = segment->splicesites_i) >= 0) {
	  /* Known splicing */
	  /* Ends 1 (donor, plus) and 8 (antiacceptor, plus): mark known splice sites in segment */
	  while (j < nsplicesites && splicesites[j] <= segment_left + splice_pos_end) { /* Needs to be <= */
	    if (splicetypes[j] == DONOR) {
	      segment_donor_knownpos[segment_donor_nknown] = splicesites[j] - segment_left;
	      segment_donor_knowni[segment_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      segment_antiacceptor_knownpos[segment_antiacceptor_nknown] = splicesites[j] - segment_left;
	      segment_antiacceptor_knowni[segment_antiacceptor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segment_donor_knownpos[segment_donor_nknown] = querylength;
	segment_antiacceptor_knownpos[segment_antiacceptor_nknown] = querylength;

	/* Originally on plus strand.  No complement */
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	if (novelsplicingp && segment_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
	  donori_nsites = Genome_donor_positions(positions_alloc,knowni_alloc,
						 segment_donor_knownpos,segment_donor_knowni,
						 segment_left,splice_pos_start,splice_pos_end);
	  donori_positions = positions_alloc;
	  donori_knowni = knowni_alloc;
	  debug4e(
		  printf("Donor dinucleotides:");
		  for (i = 0; i < donori_nsites; i++) {
		    printf(" %d",donori_positions[i]);
		  }
		  printf("\n");
		  );
	} else {
	  donori_nsites = segment_donor_nknown;
	  donori_positions = segment_donor_knownpos;
	  donori_knowni = segment_donor_knowni;
	}

	i = 0;
	nmismatches = 0;
	while (i < donori_nsites && nmismatches <= max_mismatches_allowed) {
	  splice_pos = donori_positions[i];
	  while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos) { /* Changed from <= to < */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
#if 0
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/0,/*pos3*/splice_pos,
								  plusp,genestrand,first_read_p));
#endif
	  if (nmismatches <= max_mismatches_allowed) {
	    if (donori_knowni[i] >= 0) {
	      debug4e(printf("Known donor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_end));
		
	      if ((hit = Substring_new_donor(/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/donori_knowni[i],
					     splice_pos,/*substring_querystart*/first_querypos,
					     /*substring_queryend*/last_querypos,
					     nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
					     querylength,plusp,genestrand,
					     sensedir,segment->chrnum,segment->chroffset,
					     segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s donor: %f at %d (%d mismatches) %d..%d\n",
			       plusp == true ? "plus" : "minus",Maxent_hr_donor_prob(segment_left + splice_pos,segment->chroffset),
			       Substring_siteD_pos(hit),nmismatches,Substring_querystart(hit),Substring_queryend(hit)));
		debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		(*distant_donors)[nmismatches] = List_push((*distant_donors)[nmismatches],(void *) hit);
	      }

	    } else {
	      prob = Maxent_hr_donor_prob(segment_left + splice_pos,segment->chroffset);
	      debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
			     splice_pos,nmismatches,prob,sufficient_splice_prob_distant(splice_pos,nmismatches,prob)));
	      if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		debug4e(printf("Novel donor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			       (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_end));
		if ((hit = Substring_new_donor(/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/-1,
					       splice_pos,/*substring_querystart*/first_querypos,
					       /*substring_queryend*/last_querypos,
					       nmismatches,prob,/*left*/segment_left,query_compress,
					       querylength,plusp,genestrand,
					       sensedir,segment->chrnum,segment->chroffset,
					       segment->chrhigh,segment->chrlength)) != NULL) {
		  debug4e(printf("=> %s donor: %f at %d (%d mismatches) %d..%d\n",
				 plusp == true ? "plus" : "minus",prob,Substring_siteD_pos(hit),nmismatches,
				 Substring_querystart(hit),Substring_queryend(hit)));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*distant_donors)[nmismatches] = List_push((*distant_donors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  }

	  i++;
	}


	/* Splicing originally on minus strand.  Complement */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	if (novelsplicingp && segment_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	  antiacceptori_nsites = Genome_antiacceptor_positions(positions_alloc,knowni_alloc,
							       segment_antiacceptor_knownpos,segment_antiacceptor_knowni,
							       segment_left,splice_pos_start,splice_pos_end);
	  antiacceptori_positions = positions_alloc;
	  antiacceptori_knowni = knowni_alloc;
	  debug4e(
		  printf("Antiacceptor dinucleotides:");
		  for (i = 0; i < antiacceptori_nsites; i++) {
		    printf(" %d",antiacceptori_positions[i]);
		  }
		  printf("\n");
		  );
	} else {
	  antiacceptori_nsites = segment_antiacceptor_nknown;
	  antiacceptori_positions = segment_antiacceptor_knownpos;
	  antiacceptori_knowni = segment_antiacceptor_knowni;
	}

	i = 0;
	nmismatches = 0;
	while (i < antiacceptori_nsites && nmismatches <= max_mismatches_allowed) {
	  splice_pos = antiacceptori_positions[i];
	  while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos) { /* Changed from <= to < */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
#if 0
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/0,/*pos3*/splice_pos,
								  plusp,genestrand,first_read_p));
#endif
	  if (nmismatches <= max_mismatches_allowed) {
	    if (antiacceptori_knowni[i] >= 0) {
	      debug4e(printf("Known antiacceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_end));
	      if ((hit = Substring_new_acceptor(/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/antiacceptori_knowni[i],
						splice_pos,/*substring_querystart*/first_querypos,
						/*substring_queryend*/last_querypos,
						nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
						querylength,plusp,genestrand,
						sensedir,segment->chrnum,segment->chroffset,
						segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s antiacceptor : %f at %d (%d mismatches) %d..%d\n",
			       plusp == true ? "plus" : "minus",Maxent_hr_antiacceptor_prob(segment_left + splice_pos,segment->chroffset),
			       Substring_siteA_pos(hit),nmismatches,Substring_querystart(hit),Substring_queryend(hit)));
		debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		(*distant_antiacceptors)[nmismatches] = List_push((*distant_antiacceptors)[nmismatches],(void *) hit);
	      }

	    } else {
	      prob = Maxent_hr_antiacceptor_prob(segment_left + splice_pos,segment->chroffset);
	      debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
			     splice_pos,nmismatches,prob,sufficient_splice_prob_distant(splice_pos,nmismatches,prob)));
	      if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		debug4e(printf("Novel antiacceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			       (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_end));
		if ((hit = Substring_new_acceptor(/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/-1,
						  splice_pos,/*substring_querystart*/first_querypos,
						  /*substring_queryend*/last_querypos,
						  nmismatches,prob,/*left*/segment_left,query_compress,
						  querylength,plusp,genestrand,
						  sensedir,segment->chrnum,segment->chroffset,
						  segment->chrhigh,segment->chrlength)) != NULL) {
		  debug4e(printf("=> %s antiacceptor : %f at %d (%d mismatches) %d..%d\n",
				 plusp == true ? "plus" : "minus",prob,Substring_siteA_pos(hit),nmismatches,
				 Substring_querystart(hit),Substring_queryend(hit)));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*distant_antiacceptors)[nmismatches] = List_push((*distant_antiacceptors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  }

	  i++;
	}
      }


      /* Find splices on genomic left */
      if (plusp) {
	/* ? require that floors_to_pos3[segment->querypos3] <= max_mismatches_allowed */
	if (first_querypos < index1part && last_querypos < query_lastpos) {
	  /* genomic left anchor.  No need to find splices on genomic left. */
	  debug4e(printf("Searching genomic left: plus genomic left anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;

	} else if (first_querypos > index1part && last_querypos > query_lastpos) {
	  /* genomic right anchor */
	  debug4e(printf("Searching genomic left: plus genomic right anchor\n"));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      query_compress,/*left*/segment_left,
#if 0
						      /*pos5*/0,/*pos3*/querylength,
#else
						      /*pos5*/0,/*pos3 (was querylength)*/last_querypos,
#endif
						      plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

#if 0
	  splice_pos_end = query_lastpos;
#else
	  splice_pos_end = last_querypos;
#endif
	  if (nmismatches_right <= max_mismatches_allowed) {
	    splice_pos_start = 1;
	  } else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    splice_pos_start = 1;
	  }

	} else if (first_querypos > index1part && last_querypos < query_lastpos &&
		   segment->spliceable_high_p == true) {
	  /* middle anchor */
	  debug4e(printf("Searching genomic left: plus middle anchor\n"));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      query_compress,/*left*/segment_left,
						      /*pos5*/0,/*pos3*/last_querypos,
						      plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_end = last_querypos;
	  if (nmismatches_right <= max_mismatches_allowed) {
	    splice_pos_start = 1;
	  } else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    splice_pos_start = 1;
	  }

	} else {
	  /* full anchor */
	  debug4e(printf("Searching genomic left: plus full anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	}

      } else {
	/* ? require that floors_from_neg3[segment->querypos5] <= max_mismatches_allowed */
	if (first_querypos < index1part && last_querypos < query_lastpos) {
	  /* genomic right anchor */
	  debug4e(printf("Searching genomic left: minus genomic right anchor\n"));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      query_compress,/*left*/segment_left,
#if 0
						      /*pos5*/0,/*pos3*/querylength,
#else						      
						      /*pos5*/0,/*pos3 (was querylength)*/querylength - first_querypos,
#endif
						      plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

#if 0
	  splice_pos_end = query_lastpos;
#else
	  splice_pos_end = querylength - first_querypos;
#endif
	  if (nmismatches_right <= max_mismatches_allowed) {
	    splice_pos_start = 1;
	  } else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    splice_pos_start = 1;
	  }

	} else if (first_querypos > index1part && last_querypos > query_lastpos) {
	  /* genomic left anchor.  No need to find splices on genomic left. */
	  debug4e(printf("Searching genomic left: minus genomic left anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;

	} else if (first_querypos > index1part && last_querypos < query_lastpos &&
		   segment->spliceable_high_p == true) {
	  /* middle anchor */
	  debug4e(printf("Searching genomic left: minus middle anchor\n"));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      query_compress,/*left*/segment_left,
						      /*pos5*/0,/*pos3*/querylength - first_querypos,
						      plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_end = querylength - first_querypos;
	  if (nmismatches_right <= max_mismatches_allowed) {
	    splice_pos_start = 1;
	  } else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    splice_pos_start = 1;
	  }

	} else {
	  /* full anchor */
	  debug4e(printf("Searching genomic left: minus full anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	}
      }

      if (splice_pos_start <= splice_pos_end) {
	debug4e(printf("Search for splice sites from %d down to %d\n",splice_pos_end,splice_pos_start));

	segment_acceptor_nknown = 0;
	segment_antidonor_nknown = 0;
	if ((j = segment->splicesites_i) >= 0) {
	  /* Known splicing */
	  while (j < nsplicesites && splicesites[j] <= segment_left + splice_pos_end) { /* Needs to be <= */
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4k(printf("Setting known acceptor %d for segment at %llu\n",j,(unsigned long long) splicesites[j]));
	      segment_acceptor_knownpos[segment_acceptor_nknown] = splicesites[j] - segment_left;
	      segment_acceptor_knowni[segment_acceptor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4k(printf("Setting known antidonor %d for segment at %llu\n",j,(unsigned long long) splicesites[j]));
	      segment_antidonor_knownpos[segment_antidonor_nknown] = splicesites[j] - segment_left;
	      segment_antidonor_knowni[segment_antidonor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segment_acceptor_knownpos[segment_acceptor_nknown] = querylength;
	segment_antidonor_knownpos[segment_antidonor_nknown] = querylength;

	/* Splicing originally on plus strand.  No complement. */
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	if (novelsplicingp && segment_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	  acceptorj_nsites = Genome_acceptor_positions(positions_alloc,knowni_alloc,
						       segment_acceptor_knownpos,segment_acceptor_knowni,
						       segment_left,splice_pos_start,splice_pos_end);
	  acceptorj_positions = positions_alloc;
	  acceptorj_knowni = knowni_alloc;
	  debug4e(
		  printf("Acceptor dinucleotides:");
		  for (i = 0; i < acceptorj_nsites; i++) {
		    printf(" %d",acceptorj_positions[i]);
		  }
		  printf("\n");
		  );
	} else {
	  acceptorj_nsites = segment_acceptor_nknown;
	  acceptorj_positions = segment_acceptor_knownpos;
	  acceptorj_knowni = segment_acceptor_knowni;
	}

	i = acceptorj_nsites - 1;
	nmismatches = 0;
	while (i >= 0 && nmismatches <= max_mismatches_allowed) {
	  splice_pos = acceptorj_positions[i];
	  while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
#if 0
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
								  plusp,genestrand,first_read_p));
#endif
	  if (nmismatches <= max_mismatches_allowed) {
	    if (acceptorj_knowni[i] >= 0) {
	      debug4e(printf("Known acceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
	      if ((hit = Substring_new_acceptor(/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/acceptorj_knowni[i],
						splice_pos,/*substring_querystart*/segment->querypos5,
						/*substring_queryend*/last_querypos,
						nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
						querylength,plusp,genestrand,
						sensedir,segment->chrnum,segment->chroffset,
						segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s acceptor: %f at %d (%d mismatches) %d..%d\n",
			       plusp == true ? "plus" : "minus",Maxent_hr_acceptor_prob(segment_left + splice_pos,segment->chroffset),
			       Substring_siteA_pos(hit),nmismatches,Substring_querystart(hit),Substring_queryend(hit)));
		debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		(*distant_acceptors)[nmismatches] = List_push((*distant_acceptors)[nmismatches],(void *) hit);
	      }

	    } else {
	      prob = Maxent_hr_acceptor_prob(segment_left + splice_pos,segment->chroffset);
	      debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
			     splice_pos,nmismatches,prob,sufficient_splice_prob_distant(querylength - splice_pos,nmismatches,prob)));
	      if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		debug4e(printf("Novel acceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			       (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
		if ((hit = Substring_new_acceptor(/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/-1,
						  splice_pos,/*substring_querystart*/segment->querypos5,
						  /*substring_queryend*/last_querypos,
						  nmismatches,prob,/*left*/segment_left,query_compress,
						  querylength,plusp,genestrand,
						  sensedir,segment->chrnum,segment->chroffset,
						  segment->chrhigh,segment->chrlength)) != NULL) {
		  debug4e(printf("=> %s acceptor: %f at %d (%d mismatches) %d..%d\n",
				 plusp == true ? "plus" : "minus",prob,Substring_siteA_pos(hit),nmismatches,
				 Substring_querystart(hit),Substring_queryend(hit)));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*distant_acceptors)[nmismatches] = List_push((*distant_acceptors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  }
	  
	  i--;
	}


	/* Splicing originally on minus strand.  Complement.  */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	if (novelsplicingp && segment_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
	  antidonorj_nsites = Genome_antidonor_positions(positions_alloc,knowni_alloc,
							 segment_antidonor_knownpos,segment_antidonor_knowni,
							 segment_left,splice_pos_start,splice_pos_end);
	  antidonorj_positions = positions_alloc;
	  antidonorj_knowni = knowni_alloc;
	  debug4e(
		  printf("Antidonor dinucleotides:");
		  for (i = 0; i < antidonorj_nsites; i++) {
		    printf(" %d",antidonorj_positions[i]);
		  }
		  printf("\n");
		  );
	} else {
	  antidonorj_nsites = segment_antidonor_nknown;
	  antidonorj_positions = segment_antidonor_knownpos;
	  antidonorj_knowni = segment_antidonor_knowni;
	}

	i = antidonorj_nsites - 1;
	nmismatches = 0;
	while (i >= 0 && nmismatches <= max_mismatches_allowed) {
	  splice_pos = antidonorj_positions[i];
	  while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
#if 0
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
								  plusp,genestrand,first_read_p));
#endif
	  if (nmismatches <= max_mismatches_allowed) {
	    if (antidonorj_knowni[i] >= 0) {
	      debug4e(printf("Known antidonor for segmenti at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
	      if ((hit = Substring_new_donor(/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/antidonorj_knowni[i],
					     splice_pos,/*substring_querystart*/segment->querypos5,
					     /*substring_queryend*/last_querypos,
					     nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
					     querylength,plusp,genestrand,
					     sensedir,segment->chrnum,segment->chroffset,
					     segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s antidonor: %f at %d (%d mismatches) %d..%d\n",
			       plusp == true ? "plus" : "minus",Maxent_hr_antidonor_prob(segment_left + splice_pos,segment->chroffset),
			       Substring_siteD_pos(hit),nmismatches,Substring_querystart(hit),Substring_queryend(hit)));
		debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		(*distant_antidonors)[nmismatches] = List_push((*distant_antidonors)[nmismatches],(void *) hit);
	      }
		  
	    } else {
	      prob = Maxent_hr_antidonor_prob(segment_left + splice_pos,segment->chroffset);
	      debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
			     splice_pos,nmismatches,prob,sufficient_splice_prob_distant(querylength - splice_pos,nmismatches,prob)));
	      if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		debug4e(printf("Novel antidonor for segmenti at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			       (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
		if ((hit = Substring_new_donor(/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/-1,
					       splice_pos,/*substring_querystart*/segment->querypos5,
					       /*substring_queryend*/last_querypos,
					       nmismatches,prob,/*left*/segment_left,query_compress,
					       querylength,plusp,genestrand,
					       sensedir,segment->chrnum,segment->chroffset,
					       segment->chrhigh,segment->chrlength)) != NULL) {
		  debug4e(printf("=> %s antidonor: %f at %d (%d mismatches) %d..%d\n",
				 plusp == true ? "plus" : "minus",prob,Substring_siteD_pos(hit),nmismatches,
				 Substring_querystart(hit),Substring_queryend(hit)));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*distant_antidonors)[nmismatches] = List_push((*distant_antidonors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  }

	  i--;
	}
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
    FREEA(segment_donor_knownpos);
    FREEA(segment_acceptor_knownpos);
    FREEA(segment_antidonor_knownpos);
    FREEA(segment_antiacceptor_knownpos);
    FREEA(segment_donor_knowni);
    FREEA(segment_acceptor_knowni);
    FREEA(segment_antidonor_knowni);
    FREEA(segment_antiacceptor_knowni);
    FREEA(positions_alloc);
    FREEA(knowni_alloc);
  } else {
    FREE(mismatch_positions);
    FREE(segment_donor_knownpos);
    FREE(segment_acceptor_knownpos);
    FREE(segment_antidonor_knownpos);
    FREE(segment_antiacceptor_knownpos);
    FREE(segment_donor_knowni);
    FREE(segment_acceptor_knowni);
    FREE(segment_antidonor_knowni);
    FREE(segment_antiacceptor_knowni);
    FREE(positions_alloc);
    FREE(knowni_alloc);
#else
    FREE(mismatch_positions);
    FREE(segment_donor_knownpos);
    FREE(segment_acceptor_knownpos);
    FREE(segment_antidonor_knownpos);
    FREE(segment_antiacceptor_knownpos);
    FREE(segment_donor_knowni);
    FREE(segment_acceptor_knowni);
    FREE(segment_antidonor_knowni);
    FREE(segment_antiacceptor_knowni);
    FREE(positions_alloc);
    FREE(knowni_alloc);
#endif

  return;
}
#endif



/* Streamlined version */
/* Produces lists of distant_donors and distant_acceptors that are substrings */
/* Note: Call to Genome_donor_positions and similar functions can
   yield positions at (splice_pos_start - 1) through (splice_pos_end +
   1).  Therefore, need to have splice_pos_start > 0 and
   splice_pos_end < querylength. */
/* TODO: Change to lists of Stage3end_T objects, including GMAP.
   Change definition of a chimera to be two Stage3end_T objects, instead
   of two substrings. */
static void
find_spliceends_distant_rna (List_T **distant_donors, List_T **distant_antidonors,
			     List_T **distant_acceptors, List_T **distant_antiacceptors,
			     Segment_T *anchor_segments, int nanchors,
#ifdef DEBUG4E
			     char *queryptr,
#endif
			     Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress,
			     int max_mismatches_allowed, bool plusp, int genestrand) {
#ifdef DEBUG4E
  char *gbuffer;
#endif

  Segment_T segment, *p;
  Substring_T hit;
  Univcoord_T segment_left;
  int nmismatches, j, i;
  int splice_pos, first_querypos, last_querypos;
  double prob;

  int nmismatches_left, nmismatches_right;
  int *floors_from_neg3, *floors_to_pos3;
  int sensedir;

  int splice_pos_start, splice_pos_end;
  int segment_donor_nknown, segment_acceptor_nknown, segment_antidonor_nknown, segment_antiacceptor_nknown;
  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;

  int *mismatch_positions;
  int *segment_donor_knownpos, *segment_acceptor_knownpos, *segment_antidonor_knownpos, *segment_antiacceptor_knownpos,
    *segment_donor_knowni, *segment_acceptor_knowni, *segment_antidonor_knowni, *segment_antiacceptor_knowni;
  int *positions_alloc, *knowni_alloc;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_donor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_acceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_antidonor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_antiacceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_donor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_acceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_antidonor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    segment_antiacceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
    positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    segment_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
    positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_donor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_acceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_antidonor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_antiacceptor_knownpos = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_donor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_acceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_antidonor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  segment_antiacceptor_knowni = (int *) MALLOC((querylength+1)*sizeof(int));
  positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug4e(printf("Entering find_spliceends_distant_rna with %d anchor segments\n",nanchors));

  if (floors != NULL) {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    for (p = &(anchor_segments[0]); p < &(anchor_segments[nanchors]); p++) {
      segment = *p;
      assert(segment->diagonal != (Univcoord_T) -1);

      segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      if ((first_querypos = segment->querypos5 - (index1interval - 1)) < 0) {
	first_querypos = 0;
      }
      if ((last_querypos = segment->querypos3 + index1part + (index1interval - 1)) > querylength) {
	last_querypos = querylength;
      }

      debug4e(printf("find_spliceends_distant_rna: Checking up to %d mismatches at diagonal %llu (querypos %d..%d) - querylength %d = %llu, floors %d and %d, plusp %d\n",
		     max_mismatches_allowed,(unsigned long long) segment->diagonal,
		     segment->querypos5,segment->querypos3,querylength,(unsigned long long) segment_left,
		     floors_from_neg3[segment->querypos5],floors_to_pos3[segment->querypos3],
		     plusp));

      debug4e(
	      gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
	      Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
	      printf("genome 0..: %s\n",gbuffer);
	      printf("query  0..: %s\n",queryptr);
	      FREE(gbuffer);
	      );

      /* Find splices on genomic right */
      if (plusp) {
	/* ? require that floors_from_neg3[segment->querypos5] <= max_mismatches_allowed */
	if (last_querypos < query_lastpos &&
	    (first_querypos < index1part || segment->spliceable_low_p == true)) {
	  /* genomic left anchor or middle anchor */
	  debug4e(printf("Searching genomic right: plus genomic left anchor or middle anchor: %d..%d\n",
			 segment->querypos5,querylength));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,
						    /*pos5*/segment->querypos5,/*pos3*/querylength,
						    plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_start = first_querypos + 1;
	  if (nmismatches_left <= max_mismatches_allowed) {
	    splice_pos_end = querylength - 1 - 1;
	  } else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1 - 1) {
	    splice_pos_end = querylength - 1 - 1;
	  }
	  
	} else {
	  /* genomic right anchor or full anchor */
	  debug4e(printf("Searching genomic right: plus genomic right anchor or full anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	}

      } else {
	/* ? require that floors_to_pos3[segment->querypos3] <= max_mismatches_allowed */
	if (first_querypos > index1part &&
	    (last_querypos > query_lastpos || segment->spliceable_low_p == true)) {
	  /* genomic left anchor or middle anchor */
	  debug4e(printf("Searching genomic right: minus genomic left anchor or middle anchor: %d..%d\n",
			 querylength - segment->querypos3 - index1part,querylength));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    query_compress,/*left*/segment_left,
						    /*pos5*/querylength - segment->querypos3 - index1part,
						    /*pos3*/querylength,plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_start = querylength - last_querypos + 1;
	  if (nmismatches_left <= max_mismatches_allowed) {
	    splice_pos_end = querylength - 1 - 1;
	  } else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1 - 1) {
	    splice_pos_end = querylength - 1 - 1;
	  }

	} else {
	  /* genomic right anchor or full anchor */
	  debug4e(printf("Searching genomic right: minus genomic right anchor or full anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	}
      }

      if (splice_pos_start <= splice_pos_end) {
	debug4e(printf("Search for splice sites from %d up to %d\n",splice_pos_start,splice_pos_end));

	segment_donor_nknown = 0;
	segment_antiacceptor_nknown = 0;
	if ((j = segment->splicesites_i) >= 0) {
	  /* Known splicing */
	  /* Ends 1 (donor, plus) and 8 (antiacceptor, plus): mark known splice sites in segment */
	  while (j < nsplicesites && splicesites[j] <= segment_left + splice_pos_end) { /* Needs to be <= */
	    if (splicetypes[j] == DONOR) {
	      segment_donor_knownpos[segment_donor_nknown] = splicesites[j] - segment_left;
	      segment_donor_knowni[segment_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      segment_antiacceptor_knownpos[segment_antiacceptor_nknown] = splicesites[j] - segment_left;
	      segment_antiacceptor_knowni[segment_antiacceptor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segment_donor_knownpos[segment_donor_nknown] = querylength;
	segment_antiacceptor_knownpos[segment_antiacceptor_nknown] = querylength;

	/* Originally on plus strand.  No complement */
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	if (novelsplicingp && segment_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
	  assert(splice_pos_start > 0);
	  assert(splice_pos_end < querylength - 1);
	  donori_nsites = Genome_donor_positions(positions_alloc,knowni_alloc,
						 segment_donor_knownpos,segment_donor_knowni,
						 segment_left,splice_pos_start,splice_pos_end);
	  donori_positions = positions_alloc;
	  donori_knowni = knowni_alloc;
	  debug4e(
		  printf("Donor dinucleotides:");
		  for (i = 0; i < donori_nsites; i++) {
		    printf(" %d",donori_positions[i]);
		  }
		  printf("\n");
		  );
	} else {
	  donori_nsites = segment_donor_nknown;
	  donori_positions = segment_donor_knownpos;
	  donori_knowni = segment_donor_knowni;
	}

	i = 0;
	nmismatches = 0;
	while (i < donori_nsites && nmismatches <= max_mismatches_allowed) {
	  splice_pos = donori_positions[i];
	  while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos) { /* Changed from <= to < */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
#if 0
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/0,/*pos3*/splice_pos,
								  plusp,genestrand,first_read_p));
#endif
	  if (nmismatches <= max_mismatches_allowed) {
	    if (donori_knowni[i] >= 0) {
	      debug4e(printf("Known donor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_end));
		
	      if ((hit = Substring_new_donor(/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/donori_knowni[i],
					     splice_pos,/*substring_querystart*/first_querypos,
					     /*substring_queryend*/last_querypos,
					     nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
					     querylength,plusp,genestrand,
					     sensedir,segment->chrnum,segment->chroffset,
					     segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s donor: %f at %d (%d mismatches) %d..%d\n",
			       plusp == true ? "plus" : "minus",Maxent_hr_donor_prob(segment_left + splice_pos,segment->chroffset),
			       Substring_siteD_pos(hit),nmismatches,Substring_querystart(hit),Substring_queryend(hit)));
		debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		(*distant_donors)[nmismatches] = List_push((*distant_donors)[nmismatches],(void *) hit);
	      }

	    } else {
	      prob = Maxent_hr_donor_prob(segment_left + splice_pos,segment->chroffset);
	      debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
			     splice_pos,nmismatches,prob,sufficient_splice_prob_distant(splice_pos,nmismatches,prob)));
	      if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		debug4e(printf("Novel donor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			       (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_end));
		if ((hit = Substring_new_donor(/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/-1,
					       splice_pos,/*substring_querystart*/first_querypos,
					       /*substring_queryend, as last_querypos*/querylength,
					       nmismatches,prob,/*left*/segment_left,query_compress,
					       querylength,plusp,genestrand,
					       sensedir,segment->chrnum,segment->chroffset,
					       segment->chrhigh,segment->chrlength)) != NULL) {
		  debug4e(printf("=> %s donor: %f at %d (%d mismatches) %d..%d\n",
				 plusp == true ? "plus" : "minus",prob,Substring_siteD_pos(hit),nmismatches,
				 Substring_querystart(hit),Substring_queryend(hit)));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*distant_donors)[nmismatches] = List_push((*distant_donors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  }

	  i++;
	}


	/* Splicing originally on minus strand.  Complement */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	if (novelsplicingp && segment_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	  assert(splice_pos_start > 0);
	  assert(splice_pos_end < querylength - 1);
	  antiacceptori_nsites = Genome_antiacceptor_positions(positions_alloc,knowni_alloc,
							       segment_antiacceptor_knownpos,segment_antiacceptor_knowni,
							       segment_left,splice_pos_start,splice_pos_end);
	  antiacceptori_positions = positions_alloc;
	  antiacceptori_knowni = knowni_alloc;
	  debug4e(
		  printf("Antiacceptor dinucleotides:");
		  for (i = 0; i < antiacceptori_nsites; i++) {
		    printf(" %d",antiacceptori_positions[i]);
		  }
		  printf("\n");
		  );
	} else {
	  antiacceptori_nsites = segment_antiacceptor_nknown;
	  antiacceptori_positions = segment_antiacceptor_knownpos;
	  antiacceptori_knowni = segment_antiacceptor_knowni;
	}

	i = 0;
	nmismatches = 0;
	while (i < antiacceptori_nsites && nmismatches <= max_mismatches_allowed) {
	  splice_pos = antiacceptori_positions[i];
	  while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos) { /* Changed from <= to < */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
#if 0
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/0,/*pos3*/splice_pos,
								  plusp,genestrand,first_read_p));
#endif
	  if (nmismatches <= max_mismatches_allowed) {
	    if (antiacceptori_knowni[i] >= 0) {
	      debug4e(printf("Known antiacceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_end));
	      if ((hit = Substring_new_acceptor(/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/antiacceptori_knowni[i],
						splice_pos,/*substring_querystart*/first_querypos,
						/*substring_queryend*/last_querypos,
						nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
						querylength,plusp,genestrand,
						sensedir,segment->chrnum,segment->chroffset,
						segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s antiacceptor : %f at %d (%d mismatches) %d..%d\n",
			       plusp == true ? "plus" : "minus",Maxent_hr_antiacceptor_prob(segment_left + splice_pos,segment->chroffset),
			       Substring_siteA_pos(hit),nmismatches,Substring_querystart(hit),Substring_queryend(hit)));
		debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		(*distant_antiacceptors)[nmismatches] = List_push((*distant_antiacceptors)[nmismatches],(void *) hit);
	      }

	    } else {
	      prob = Maxent_hr_antiacceptor_prob(segment_left + splice_pos,segment->chroffset);
	      debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
			     splice_pos,nmismatches,prob,sufficient_splice_prob_distant(splice_pos,nmismatches,prob)));
	      if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		debug4e(printf("Novel antiacceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			       (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_end));
		if ((hit = Substring_new_acceptor(/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/-1,
						  splice_pos,/*substring_querystart*/first_querypos,
						  /*substring_queryend*/last_querypos,
						  nmismatches,prob,/*left*/segment_left,query_compress,
						  querylength,plusp,genestrand,
						  sensedir,segment->chrnum,segment->chroffset,
						  segment->chrhigh,segment->chrlength)) != NULL) {
		  debug4e(printf("=> %s antiacceptor : %f at %d (%d mismatches) %d..%d\n",
				 plusp == true ? "plus" : "minus",prob,Substring_siteA_pos(hit),nmismatches,
				 Substring_querystart(hit),Substring_queryend(hit)));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*distant_antiacceptors)[nmismatches] = List_push((*distant_antiacceptors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  }

	  i++;
	}
      }


      /* Find splices on genomic left */
      if (plusp) {
	/* ? require that floors_to_pos3[segment->querypos3] <= max_mismatches_allowed */
	if (first_querypos > index1part && 
	    (last_querypos > query_lastpos || segment->spliceable_high_p == true)) {
	  /* genomic right anchor or middle anchor */
	  debug4e(printf("Searching genomic left: plus genomic right anchor or middle anchor: %d..%d\n",
			 0,segment->querypos3 + index1part));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      query_compress,/*left*/segment_left,
						      /*pos5*/0,/*pos3*/segment->querypos3 + index1part,
						      plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );
	  
	  splice_pos_end = last_querypos - 1 - 1;
	  if (nmismatches_right <= max_mismatches_allowed) {
	    splice_pos_start = 1;
	  } else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    splice_pos_start = 1;
	  }
	  
	} else {
	  /* genomic left anchor or full anchor */
	  debug4e(printf("Searching genomic left: plus genomic left anchor or full anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	}

      } else {
	/* ? require that floors_from_neg3[segment->querypos5] <= max_mismatches_allowed */
	if (last_querypos < query_lastpos &&
	    (first_querypos < index1part || segment->spliceable_high_p == true)) {
	  /* genomic right anchor or middle anchor*/
	  debug4e(printf("Searching genomic left: minus genomic right anchor or middle anchor: %d..%d\n",
			 0,querylength - segment->querypos5));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      query_compress,/*left*/segment_left,
						      /*pos5*/0,/*pos3*/querylength - segment->querypos5,
						      plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  splice_pos_end = querylength - first_querypos - 1 - 1;
	  if (nmismatches_right <= max_mismatches_allowed) {
	    splice_pos_start = 1;
	  } else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    splice_pos_start = 1;
	  }

	} else {
	  /* genomic left anchor or full anchor */
	  debug4e(printf("Searching genomic left: minus genomic left anchor or full anchor\n"));
	  splice_pos_start = querylength;
	  splice_pos_end = 0;
	}
      }

      if (splice_pos_start <= splice_pos_end) {
	debug4e(printf("Search for splice sites from %d down to %d\n",splice_pos_end,splice_pos_start));

	segment_acceptor_nknown = 0;
	segment_antidonor_nknown = 0;
	if ((j = segment->splicesites_i) >= 0) {
	  /* Known splicing */
	  while (j < nsplicesites && splicesites[j] <= segment_left + splice_pos_end) { /* Needs to be <= */
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4k(printf("Setting known acceptor %d for segment at %llu\n",j,(unsigned long long) splicesites[j]));
	      segment_acceptor_knownpos[segment_acceptor_nknown] = splicesites[j] - segment_left;
	      segment_acceptor_knowni[segment_acceptor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4k(printf("Setting known antidonor %d for segment at %llu\n",j,(unsigned long long) splicesites[j]));
	      segment_antidonor_knownpos[segment_antidonor_nknown] = splicesites[j] - segment_left;
	      segment_antidonor_knowni[segment_antidonor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segment_acceptor_knownpos[segment_acceptor_nknown] = querylength;
	segment_antidonor_knownpos[segment_antidonor_nknown] = querylength;

	/* Splicing originally on plus strand.  No complement. */
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	if (novelsplicingp && segment_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	  assert(splice_pos_start > 0);
	  assert(splice_pos_end < querylength - 1);
	  acceptorj_nsites = Genome_acceptor_positions(positions_alloc,knowni_alloc,
						       segment_acceptor_knownpos,segment_acceptor_knowni,
						       segment_left,splice_pos_start,splice_pos_end);
	  acceptorj_positions = positions_alloc;
	  acceptorj_knowni = knowni_alloc;
	  debug4e(
		  printf("Acceptor dinucleotides:");
		  for (i = 0; i < acceptorj_nsites; i++) {
		    printf(" %d",acceptorj_positions[i]);
		  }
		  printf("\n");
		  );
	} else {
	  acceptorj_nsites = segment_acceptor_nknown;
	  acceptorj_positions = segment_acceptor_knownpos;
	  acceptorj_knowni = segment_acceptor_knowni;
	}

	i = acceptorj_nsites - 1;
	nmismatches = 0;
	while (i >= 0 && nmismatches <= max_mismatches_allowed) {
	  splice_pos = acceptorj_positions[i];
	  while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
#if 0
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
								  plusp,genestrand,first_read_p));
#endif
	  if (nmismatches <= max_mismatches_allowed) {
	    if (acceptorj_knowni[i] >= 0) {
	      debug4e(printf("Known acceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
	      if ((hit = Substring_new_acceptor(/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/acceptorj_knowni[i],
						splice_pos,/*substring_querystart*/first_querypos,
						/*substring_queryend*/last_querypos,
						nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
						querylength,plusp,genestrand,
						sensedir,segment->chrnum,segment->chroffset,
						segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s acceptor: %f at %d (%d mismatches) %d..%d\n",
			       plusp == true ? "plus" : "minus",Maxent_hr_acceptor_prob(segment_left + splice_pos,segment->chroffset),
			       Substring_siteA_pos(hit),nmismatches,Substring_querystart(hit),Substring_queryend(hit)));
		debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		(*distant_acceptors)[nmismatches] = List_push((*distant_acceptors)[nmismatches],(void *) hit);
	      }

	    } else {
	      prob = Maxent_hr_acceptor_prob(segment_left + splice_pos,segment->chroffset);
	      debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
			     splice_pos,nmismatches,prob,sufficient_splice_prob_distant(querylength - splice_pos,nmismatches,prob)));
	      if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		debug4e(printf("Novel acceptor for segment at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			       (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
		if ((hit = Substring_new_acceptor(/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/-1,
						  splice_pos,/*substring_querystart*/first_querypos,
						  /*substring_queryend*/last_querypos,
						  nmismatches,prob,/*left*/segment_left,query_compress,
						  querylength,plusp,genestrand,
						  sensedir,segment->chrnum,segment->chroffset,
						  segment->chrhigh,segment->chrlength)) != NULL) {
		  debug4e(printf("=> %s acceptor: %f at %d (%d mismatches) %d..%d\n",
				 plusp == true ? "plus" : "minus",prob,Substring_siteA_pos(hit),nmismatches,
				 Substring_querystart(hit),Substring_queryend(hit)));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*distant_acceptors)[nmismatches] = List_push((*distant_acceptors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  }
	  
	  i--;
	}


	/* Splicing originally on minus strand.  Complement.  */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	if (novelsplicingp && segment_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
	  assert(splice_pos_start > 0);
	  assert(splice_pos_end < querylength - 1);
	  antidonorj_nsites = Genome_antidonor_positions(positions_alloc,knowni_alloc,
							 segment_antidonor_knownpos,segment_antidonor_knowni,
							 segment_left,splice_pos_start,splice_pos_end);
	  antidonorj_positions = positions_alloc;
	  antidonorj_knowni = knowni_alloc;
	  debug4e(
		  printf("Antidonor dinucleotides:");
		  for (i = 0; i < antidonorj_nsites; i++) {
		    printf(" %d",antidonorj_positions[i]);
		  }
		  printf("\n");
		  );
	} else {
	  antidonorj_nsites = segment_antidonor_nknown;
	  antidonorj_positions = segment_antidonor_knownpos;
	  antidonorj_knowni = segment_antidonor_knowni;
	}

	i = antidonorj_nsites - 1;
	nmismatches = 0;
	while (i >= 0 && nmismatches <= max_mismatches_allowed) {
	  splice_pos = antidonorj_positions[i];
	  while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }
	  debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
#if 0
	  assert(nmismatches == Genome_count_mismatches_substring(query_compress,segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
								  plusp,genestrand,first_read_p));
#endif
	  if (nmismatches <= max_mismatches_allowed) {
	    if (antidonorj_knowni[i] >= 0) {
	      debug4e(printf("Known antidonor for segmenti at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			     (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
	      if ((hit = Substring_new_donor(/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/antidonorj_knowni[i],
					     splice_pos,/*substring_querystart*/first_querypos,
					     /*substring_queryend*/last_querypos,
					     nmismatches,/*prob*/2.0,/*left*/segment_left,query_compress,
					     querylength,plusp,genestrand,
					     sensedir,segment->chrnum,segment->chroffset,
					     segment->chrhigh,segment->chrlength)) != NULL) {
		debug4e(printf("=> %s antidonor: %f at %d (%d mismatches) %d..%d\n",
			       plusp == true ? "plus" : "minus",Maxent_hr_antidonor_prob(segment_left + splice_pos,segment->chroffset),
			       Substring_siteD_pos(hit),nmismatches,Substring_querystart(hit),Substring_queryend(hit)));
		debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		(*distant_antidonors)[nmismatches] = List_push((*distant_antidonors)[nmismatches],(void *) hit);
	      }
		  
	    } else {
	      prob = Maxent_hr_antidonor_prob(segment_left + splice_pos,segment->chroffset);
	      debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
			     splice_pos,nmismatches,prob,sufficient_splice_prob_distant(querylength - splice_pos,nmismatches,prob)));
	      if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		debug4e(printf("Novel antidonor for segmenti at %llu, splice_pos %d (%d mismatches), stopi = %d\n",
			       (unsigned long long) segment_left,splice_pos,nmismatches,splice_pos_start));
		if ((hit = Substring_new_donor(/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/-1,
					       splice_pos,/*substring_querystart*/first_querypos,
					       /*substring_queryend*/last_querypos,
					       nmismatches,prob,/*left*/segment_left,query_compress,
					       querylength,plusp,genestrand,
					       sensedir,segment->chrnum,segment->chroffset,
					       segment->chrhigh,segment->chrlength)) != NULL) {
		  debug4e(printf("=> %s antidonor: %f at %d (%d mismatches) %d..%d\n",
				 plusp == true ? "plus" : "minus",prob,Substring_siteD_pos(hit),nmismatches,
				 Substring_querystart(hit),Substring_queryend(hit)));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*distant_antidonors)[nmismatches] = List_push((*distant_antidonors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  }

	  i--;
	}
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
    FREEA(segment_donor_knownpos);
    FREEA(segment_acceptor_knownpos);
    FREEA(segment_antidonor_knownpos);
    FREEA(segment_antiacceptor_knownpos);
    FREEA(segment_donor_knowni);
    FREEA(segment_acceptor_knowni);
    FREEA(segment_antidonor_knowni);
    FREEA(segment_antiacceptor_knowni);
    FREEA(positions_alloc);
    FREEA(knowni_alloc);
  } else {
    FREE(mismatch_positions);
    FREE(segment_donor_knownpos);
    FREE(segment_acceptor_knownpos);
    FREE(segment_antidonor_knownpos);
    FREE(segment_antiacceptor_knownpos);
    FREE(segment_donor_knowni);
    FREE(segment_acceptor_knowni);
    FREE(segment_antidonor_knowni);
    FREE(segment_antiacceptor_knowni);
    FREE(positions_alloc);
    FREE(knowni_alloc);
  }
#else
  FREE(mismatch_positions);
  FREE(segment_donor_knownpos);
  FREE(segment_acceptor_knownpos);
  FREE(segment_antidonor_knownpos);
  FREE(segment_antiacceptor_knownpos);
  FREE(segment_donor_knowni);
  FREE(segment_acceptor_knowni);
  FREE(segment_antidonor_knowni);
  FREE(segment_antiacceptor_knowni);
  FREE(positions_alloc);
  FREE(knowni_alloc);
#endif

  return;
}



/* Integrates terminals found from ends by counting mismatches, and
   those where querypos3 - querypos5 is long enough */
static List_T
find_terminals (Segment_T *plus_anchor_segments, Segment_T *minus_anchor_segments,
		int n_plus_anchors, int n_minus_anchors,
		int querylength, int query_lastpos,
		Compress_T query_compress_fwd, Compress_T query_compress_rev,
		int max_mismatches_allowed, int genestrand) {
#ifdef DEBUG4T
  char *gbuffer;
#endif
  List_T plus_terminals_middle = NULL, plus_terminals_left = NULL, plus_terminals_right = NULL,
    minus_terminals_middle = NULL, minus_terminals_left = NULL, minus_terminals_right = NULL, q;
  Segment_T segment, *p;
  Stage3end_T hit;
  Univcoord_T segment_left;
  int nmismatches_left, nmismatches_right;
  Endtype_T start_endtype, end_endtype;
  /* int *floors_from_neg3, *floors_to_pos3; */
  int max_terminal_length;
  int nterminals_left, nterminals_right, nterminals_middle;

#ifdef DEBUG4T
  int i;
#endif

  int *mismatch_positions;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


  debug(printf("identify_terminals: Checking up to %d mismatches\n",max_mismatches_allowed));

  /* floors_from_neg3 = floors->scorefrom[-index1interval]; */
  /* floors_to_pos3 = floors->scoreto[query_lastpos+index1interval]; */

  /* Needs to be /3 for long_terminals and short_terminals to work */
  max_terminal_length = querylength/3;
  if (max_terminal_length < index1part) {
    max_terminal_length = index1part;
  }

  nterminals_left = nterminals_right = nterminals_middle = 0;
  for (p = &(plus_anchor_segments[0]);
       p < &(plus_anchor_segments[n_plus_anchors]) && (/*nterminals_middle < MAX_NTERMINALS ||*/ nterminals_left < MAX_NTERMINALS || nterminals_right < MAX_NTERMINALS);
       p++) {
    segment = *p;
    if (0 && segment->usedp == true) {
      /* Previously skipped, but looks like a bad idea */
    } else if (segment->diagonal < (Univcoord_T) -1) {
      debug4t(printf("plus: %llu, %d..%d\n",(unsigned long long) segment->diagonal,segment->querypos5,segment->querypos3));
      segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      debug4t(printf("identify_terminals_plus: Checking up to %d mismatches at diagonal %llu (querypos %d..%d) - querylength %d = %llu\n",
		     max_mismatches_allowed,(unsigned long long) segment->diagonal,
		     segment->querypos5,segment->querypos3,querylength,(unsigned long long) segment_left));
      debug4t(
	      gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
	      Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
	      printf("genome 0..: %s\n",gbuffer);
	      /* printf("query  0..: %s\n",queryuc_ptr); */
	      FREE(gbuffer);
	      );

#ifdef ALLOW_MIDDLE_ALIGNMENTS
      if (segment->querypos3 - segment->querypos5 > max_terminal_length /* was index1part */) {
	/* Check for middle section */
	debug4t(printf(" => ? Middle alignment based on querypos3 %d - querypos5 %d > max_terminal_length %d",
		       segment->querypos3,segment->querypos5,max_terminal_length));
	if (nterminals_middle >= MAX_NTERMINALS) {
	  /* Skip */
	  debug4t(printf(" => Skipping because too many nterminals_middle"));
	} else {
	  start_endtype = (segment->querypos5 < index1interval) ? END : TERM;
	  end_endtype = (segment->querypos3 >= query_lastpos - index1interval) ? END : TERM;
	  debug4t(printf("  querypos3 %d vs index1interval %d => start_endtype %s\n",
			 segment->querypos3,index1interval,Endtype_string(start_endtype)));
	  debug4t(printf("  querypos5 %d vs query_lastpos %d - index1interval %d => end_endtype %s\n",
			 segment->querypos5,query_lastpos,index1interval,Endtype_string(end_endtype)));
	  
	  if ((hit = Stage3end_new_terminal(/*querystart*/0,/*queryend*//*truncate_pos_left*/querylength,
					    /*left*/segment_left,query_compress_fwd,
					    querylength,/*plusp*/true,genestrand,
					    start_endtype,end_endtype,segment->chrnum,segment->chroffset,
					    segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					    /*sarrayp*/false)) != NULL) {
	    debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	    plus_terminals_middle = List_push(plus_terminals_middle,(void *) hit);
	    nterminals_middle += 1;
	  } else {
	    debug4t(printf(" => no"));
	  }
	}
	debug4t(printf("\n"));
	
      } else {
#endif
	
	if (nterminals_left >= MAX_NTERMINALS) {
	  /* Skip */
	} else if (segment->floor_left > max_mismatches_allowed) {
	  debug4t(printf("Not checking left because floor_left %d > max_mismatches_allowed %d\n",
			 segment->floor_left,max_mismatches_allowed));
	} else {
	  /* Check from left */
	  debug4t(printf("Checking left because floor_left %d <= max_mismatches_allowed %d\n",
			 segment->floor_left,max_mismatches_allowed));

	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    query_compress_fwd,/*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						    /*plusp*/true,genestrand);
	    
	  debug4t(
		  printf("%d mismatches on left at:",nmismatches_left);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  if (nmismatches_left == 0 || nmismatches_left <= max_mismatches_allowed || 
	      mismatch_positions[nmismatches_left-1] > querylength - max_terminal_length) {
	    debug4t(printf(" => Long terminal at left: nmismatches_left %d vs max_mismatches_allowed %d, last mismatch %d vs terminal pos %d",
			   nmismatches_left,max_mismatches_allowed,mismatch_positions[nmismatches_left-1],querylength - max_terminal_length));
	    if ((hit = Stage3end_new_terminal(/*querystart*/0,/*queryend*//*truncate_pos_left*/querylength,
					      /*left*/segment_left,query_compress_fwd,
					      querylength,/*plusp*/true,genestrand,
					      /*start_endtype*/END,/*end_endtype*/TERM,
					      segment->chrnum,segment->chroffset,
					      segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					      /*sarrayp*/false)) != NULL) {
	      debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	      plus_terminals_left = List_push(plus_terminals_left,(void *) hit);
	      nterminals_left += 1;
	    } else {
	      debug4t(printf(" => no"));
	    }
	    debug4t(printf("\n"));

	  } else if (mismatch_positions[(nmismatches_left-1)/2] > max_terminal_length) {
	    debug4t(printf(" => Short terminal at left: nmismatches_left %d vs max_mismatches_allowed %d, last mismatch %d vs terminal pos %d",
			   nmismatches_left,max_mismatches_allowed,mismatch_positions[(nmismatches_left-1)/2],max_terminal_length));

	    if ((hit = Stage3end_new_terminal(/*querystart*/0,/*queryend*//*truncate_pos_left*/querylength,
					      /*left*/segment_left,query_compress_fwd,
					      querylength,/*plusp*/true,genestrand,
					      /*start_endtype*/END,/*end_endtype*/TERM,
					      segment->chrnum,segment->chroffset,
					      segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					      /*sarrayp*/false)) != NULL) {
	      debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	      plus_terminals_left = List_push(plus_terminals_left,(void *) hit);
	      nterminals_left += 1;
	    } else {
	      debug4t(printf(" => no"));
	    }
	    debug4t(printf("\n"));

	  }
	}

	if (nterminals_right >= MAX_NTERMINALS) {
	  /* Skip */
	} else if (segment->floor_right > max_mismatches_allowed) {
	  debug4t(printf("Not checking right because floor_right %d > max_mismatches_allowed %d\n",
			 segment->floor_right,max_mismatches_allowed));
	} else {
	  /* Check from right */
	  debug4t(printf("Checking right because floor_right %d <= max_mismatches_allowed %d\n",
			 segment->floor_right,max_mismatches_allowed));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      /*query_compress*/query_compress_fwd,
						      /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						      /*plusp*/true,genestrand);
	    
	  debug4t(
		  printf("%d mismatches on right at:",nmismatches_right);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );
	    
	  debug4t(printf("last mismatch %d, half mismatch %d, long terminalpos %d, short terminalpos %d\n",
			 mismatch_positions[nmismatches_right-1],mismatch_positions[(nmismatches_right-1)/2],
			 max_terminal_length,querylength - max_terminal_length));

	  if (nmismatches_right == 0 || nmismatches_right <= max_mismatches_allowed ||
	      mismatch_positions[nmismatches_right-1] < max_terminal_length) {
	    debug4t(printf(" => Long terminal at right: nmismatches_right %d vs max_mismatches_allowed %d, last mismatch %d vs terminal pos %d",
			   nmismatches_right,max_mismatches_allowed,mismatch_positions[nmismatches_right-1],max_terminal_length));
	    if ((hit = Stage3end_new_terminal(/*querystart*//*truncate_pos_right*/0,/*queryend*/querylength,
					      /*left*/segment_left,query_compress_fwd,
					      querylength,/*plusp*/true,genestrand,
					      /*start_endtype*/TERM,/*end_endtype*/END,
					      segment->chrnum,segment->chroffset,
					      segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					      /*sarrayp*/false)) != NULL) {
	      debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	      plus_terminals_right = List_push(plus_terminals_right,(void *) hit);
	      nterminals_right += 1;
	    } else {
	      debug4t(printf(" => no"));
	    }
	    debug4t(printf("\n"));

	  } else if (mismatch_positions[(nmismatches_right-1)/2] < querylength - max_terminal_length) {
	    debug4t(printf(" => Short terminal at right: nmismatches_right %d vs max_mismatches_allowed %d, last mismatch %d vs terminal pos %d",
			   nmismatches_right,max_mismatches_allowed,mismatch_positions[(nmismatches_right-1)/2],querylength-max_terminal_length));
	    if ((hit = Stage3end_new_terminal(/*querystart*//*truncate_pos_right*/0,/*queryend*/querylength,
					      /*left*/segment_left,query_compress_fwd,
					      querylength,/*plusp*/true,genestrand,
					      /*start_endtype*/TERM,/*end_endtype*/END,
					      segment->chrnum,segment->chroffset,
					      segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					      /*sarrayp*/false)) != NULL) {
	      debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	      plus_terminals_right = List_push(plus_terminals_right,(void *) hit);
	      nterminals_right += 1;
	    } else {
	      debug4t(printf(" => no"));
	    }
	    debug4t(printf("\n"));


	  }
	}
#ifdef ALLOW_MIDDLE_ALIGNMENTS
      }
#endif
    }
  }


  if (nterminals_middle >= MAX_NTERMINALS) {
    for (q = plus_terminals_middle; q != NULL; q = q->rest) {
      hit = (Stage3end_T) q->first;
      Stage3end_free(&hit);
    }
    List_free(&plus_terminals_middle);
    plus_terminals_middle = (List_T) NULL;
  }

  if (nterminals_left >= MAX_NTERMINALS) {
    for (q = plus_terminals_left; q != NULL; q = q->rest) {
      hit = (Stage3end_T) q->first;
      Stage3end_free(&hit);
    }
    List_free(&plus_terminals_left);
    plus_terminals_left = (List_T) NULL;
  }

  if (nterminals_right >= MAX_NTERMINALS) {
    for (q = plus_terminals_right; q != NULL; q = q->rest) {
      hit = (Stage3end_T) q->first;
      Stage3end_free(&hit);
    }
    List_free(&plus_terminals_right);
    plus_terminals_right = (List_T) NULL;
  }

  nterminals_left = nterminals_right = nterminals_middle = 0;
  for (p = &(minus_anchor_segments[0]);
       p < &(minus_anchor_segments[n_minus_anchors]) && (/*nterminals_middle < MAX_NTERMINALS ||*/ nterminals_left < MAX_NTERMINALS || nterminals_right < MAX_NTERMINALS);
       p++) {
    segment = *p;
    if (0 && segment->usedp == true) {
      /* Previously skipped, but looks like a bad idea */
      debug4t(printf("segment used\n"));
    } else if (segment->diagonal < (Univcoord_T) -1) {
      debug4t(printf("minus: %llu, %d..%d\n",(unsigned long long) segment->diagonal,segment->querypos5,segment->querypos3));
      segment_left = segment->diagonal - querylength;
      debug4t(printf("identify_terminals_minus: Getting genome at diagonal %llu (querypos %d..%d) + 12 - querylength %d = %llu\n",
		     (unsigned long long) segment->diagonal,segment->querypos5,segment->querypos3,querylength,
		     (unsigned long long) segment_left));
      debug4t(
	      gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
	      Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
	      printf("genome   0..: %s\n",gbuffer);
	      /* printf("query.rc 0..: %s\n",queryrc); */
	      FREE(gbuffer);
	      );


#ifdef ALLOW_MIDDLE_ALIGNMENTS
      if (segment->querypos3 - segment->querypos5 > max_terminal_length /* was index1part */) {
	/* Check for a middle section */
	debug4t(printf(" => ? Middle alignment based on querypos3 %d - querypos5 %d > max_terminal_length %d",
		       segment->querypos3,segment->querypos5,max_terminal_length));
	if (nterminals_middle >= MAX_NTERMINALS) {
	  /* Skip */
	  debug4t(printf(" => Skipping because too many nterminals_middle"));
	} else {
	  start_endtype = (segment->querypos5 < index1interval) ? END : TERM;
	  end_endtype = (segment->querypos3 >= query_lastpos - index1interval) ? END : TERM;
	  debug4t(printf("  querypos3 %d vs index1interval %d => start_endtype %s\n",
			 segment->querypos3,index1interval,Endtype_string(start_endtype)));
	  debug4t(printf("  querypos5 %d vs query_lastpos %d - index1interval %d => end_endtype %s\n",
			 segment->querypos5,query_lastpos,index1interval,Endtype_string(end_endtype)));
	    
	  if ((hit = Stage3end_new_terminal(/*querystart*/0,/*queryend*/querylength,
					    /*left*/segment_left,query_compress_rev,
					    querylength,/*plusp*/false,genestrand,
					    start_endtype,end_endtype,segment->chrnum,segment->chroffset,
					    segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					    /*sarrayp*/false)) != NULL) {
	    debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	    minus_terminals_middle = List_push(minus_terminals_middle,(void *) hit);
	    nterminals_middle += 1;
	  } else {
	    debug4t(printf(" => no"));
	  }
	}
	debug4t(printf("\n"));

      } else {
#endif

	/* Need to reverse floor_left and floor_right */
	if (nterminals_left >= MAX_NTERMINALS) {
	  /* Skip */
	} else if (segment->floor_right > max_mismatches_allowed) {
	  debug4t(printf("Not checking left because floor_right %d > max_mismatches_allowed %d\n",
			 segment->floor_right,max_mismatches_allowed));
	} else {
	  /* Check from left */
	  debug4t(printf("Checking left because floor_right %d <= max_mismatches_allowed %d\n",
			 segment->floor_right,max_mismatches_allowed));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    /*query_compress*/query_compress_rev,
						    /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						    /*plusp*/false,genestrand);

	  debug4t(
		  printf("%d mismatches on left at:",nmismatches_left);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  if (nmismatches_left == 0 || nmismatches_left <= max_mismatches_allowed || 
	      mismatch_positions[nmismatches_left-1] > querylength - max_terminal_length) {
	    debug4t(printf(" => Long terminal at left: nmismatches_left %d vs max_mismatches_allowed %d, last mismatch %d vs terminal pos %d",
			   nmismatches_left,max_mismatches_allowed,mismatch_positions[nmismatches_left-1],querylength - max_terminal_length));
	    if ((hit = Stage3end_new_terminal(/*querystart*//*querylength-truncate_pos_left*/0,/*queryend*/querylength,
					      /*left*/segment_left,query_compress_rev,
					      querylength,/*plusp*/false,genestrand,
					      /*start_endtype*/TERM,/*end_endtype*/END,
					      segment->chrnum,segment->chroffset,
					      segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					      /*sarrayp*/false)) != NULL) {
	      debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	      minus_terminals_left = List_push(minus_terminals_left,(void *) hit);
	      nterminals_left += 1;
	    } else {
	      debug4t(printf(" => no"));
	    }
	    debug4t(printf("\n"));

	  } else if (mismatch_positions[(nmismatches_left-1)/2] > max_terminal_length) {
	    debug4t(printf(" => Short terminal at left: nmismatches_left %d vs max_mismatches_allowed %d, last mismatch %d vs terminal pos %d",
			   nmismatches_left,max_mismatches_allowed,mismatch_positions[(nmismatches_left-1)/2],max_terminal_length));
	    if ((hit = Stage3end_new_terminal(/*querystart*//*querylength-truncate_pos_left*/0,/*queryend*/querylength,
					      /*left*/segment_left,query_compress_rev,
					      querylength,/*plusp*/false,genestrand,
					      /*start_endtype*/TERM,/*end_endtype*/END,
					      segment->chrnum,segment->chroffset,
					      segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					      /*sarrayp*/false)) != NULL) {
	      debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	      minus_terminals_left = List_push(minus_terminals_left,(void *) hit);
	      nterminals_left += 1;
	    } else {
	      debug4t(printf(" => no"));
	    }
	    debug4t(printf("\n"));
	  }
	}

	if (nterminals_right >= MAX_NTERMINALS) {
	  /* Skip */
	} else if (segment->floor_left > max_mismatches_allowed) {
	  debug4t(printf("Not checking right because floor_left %d > max_mismatches_allowed %d\n",
			 segment->floor_left,max_mismatches_allowed));
	} else {
	  /* Check from right */
	  debug4t(printf("Checking right because floor_left %d <= max_mismatches_allowed %d\n",
			 segment->floor_left,max_mismatches_allowed));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      /*query_compress*/query_compress_rev,
						      /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						      /*plusp*/false,genestrand);

	  debug4t(
		  printf("%d mismatches on right at:",nmismatches_right);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  if (nmismatches_right == 0 || nmismatches_right <= max_mismatches_allowed ||
	      mismatch_positions[nmismatches_right-1] < max_terminal_length) {
	    debug4t(printf(" => Long terminal at right: nmismatches_right %d vs max_mismatches_allowed %d, last mismatch %d vs terminal pos %d",
			   nmismatches_right,max_mismatches_allowed,mismatch_positions[nmismatches_right-1],max_terminal_length));
	    if ((hit = Stage3end_new_terminal(/*querystart*/0,/*queryend*//*querylength-truncate_pos_right*/querylength,
					      /*left*/segment_left,query_compress_rev,
					      querylength,/*plusp*/false,genestrand,
					      /*start_endtype*/END,/*end_endtype*/TERM,
					      segment->chrnum,segment->chroffset,
					      segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					      /*sarrayp*/false)) != NULL) {
	      debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	      minus_terminals_right = List_push(minus_terminals_right,(void *) hit);
	      nterminals_right += 1;
	    } else {
	      debug4t(printf(" => no"));
	    }
	    debug4t(printf("\n"));

	  } else if (mismatch_positions[(nmismatches_right-1)/2] < querylength - max_terminal_length) {
	    debug4t(printf(" => Short terminal at right: nmismatches_right %d vs max_mismatches_allowed %d, last mismatch %d vs terminal pos %d",
			   nmismatches_right,max_mismatches_allowed,mismatch_positions[(nmismatches_right-1)/2],querylength-max_terminal_length));
	    if ((hit = Stage3end_new_terminal(/*querystart*/0,/*queryend*//*querylength-truncate_pos_right*/querylength,
					      /*left*/segment_left,query_compress_rev,
					      querylength,/*plusp*/false,genestrand,
					      /*start_endtype*/END,/*end_endtype*/TERM,
					      segment->chrnum,segment->chroffset,
					      segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					      /*sarrayp*/false)) != NULL) {
	      debug4t(printf(" => yes, with %d matches",Stage3end_nmatches_posttrim(hit)));
	      minus_terminals_right = List_push(minus_terminals_right,(void *) hit);
	      nterminals_right += 1;
	    } else {
	      debug4t(printf(" => no"));
	    }
	    debug4t(printf("\n"));
	  }
	}
#ifdef ALLOW_MIDDLE_ALIGNMENTS
      }
#endif
    }
  }

  if (nterminals_middle >= MAX_NTERMINALS) {
    for (q = minus_terminals_middle; q != NULL; q = q->rest) {
      hit = (Stage3end_T) q->first;
      Stage3end_free(&hit);
    }
    List_free(&minus_terminals_middle);
    minus_terminals_middle = (List_T) NULL;
  }

  if (nterminals_left >= MAX_NTERMINALS) {
    for (q = minus_terminals_left; q != NULL; q = q->rest) {
      hit = (Stage3end_T) q->first;
      Stage3end_free(&hit);
    }
    List_free(&minus_terminals_left);
    minus_terminals_left = (List_T) NULL;
  }

  if (nterminals_right >= MAX_NTERMINALS) {
    for (q = minus_terminals_right; q != NULL; q = q->rest) {
      hit = (Stage3end_T) q->first;
      Stage3end_free(&hit);
    }
    List_free(&minus_terminals_right);
    minus_terminals_right = (List_T) NULL;
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
  } else {
    FREE(mismatch_positions);
  }
#else
  FREE(mismatch_positions);
#endif

  return List_append(plus_terminals_middle,
		     List_append(plus_terminals_left,
				 List_append(plus_terminals_right,
					     List_append(minus_terminals_middle,
							 List_append(minus_terminals_left,minus_terminals_right)))));
}



#if 0
/* Finds terminals just by looking at querypos5 and querypos3 */
static List_T
find_terminals_by_width_only (struct Segment_T *plus_segments, int plus_nsegments,
			      struct Segment_T *minus_segments, int minus_nsegments,
#ifdef DEBUG4T
			      char *queryuc_ptr, /* for debugging */ char *queryrc,
#endif
			      Floors_T floors, int querylength, int query_lastpos,
			      Compress_T query_compress_fwd, Compress_T query_compress_rev,
			      int max_mismatches_allowed, int max_terminal_length,
			      int genestrand) {
#ifdef DEBUG4T
  char *gbuffer;
#endif
  List_T terminals = (List_T) NULL;
  Segment_T segment;
  Stage3end_T hit;
  Univcoord_T segment_left;

#ifdef DEBUG4T
  int i;
#endif

  debug(printf("identify_terminals: Checking up to %d mismatches\n",max_mismatches_allowed));

#if 0
  if (floors == NULL) {
    return (List_T) NULL;

  } else {
    floors_from_neg3 = floors->scorefrom[-index1interval];
    floors_to_pos3 = floors->scoreto[query_lastpos+index1interval];

    if (max_terminal_length > querylength/3) {
      max_terminal_length = querylength/3;
    }
  }
#endif

  if (plus_nsegments > 0) {
    for (segment = plus_segments; segment < &(plus_segments[plus_nsegments]); segment++) {
      if (0 && segment->usedp == true) {
	/* Previously skipped, but looks like a bad idea */
      } else if (segment->diagonal < (Univcoord_T) -1) {
	debug4t(printf("identify_terminals_plus: Checking up to %d mismatches at diagonal %llu (querypos %d..%d)\n",
		       max_mismatches_allowed,(unsigned long long) segment->diagonal,segment->querypos5,segment->querypos3));

	if (segment->querypos3 - segment->querypos5 > index1part) {
	  segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
	  debug4t(
		  gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		  Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		  printf("genome 0..: %s\n",gbuffer);
		  printf("query  0..: %s\n",queryuc_ptr);
		  FREE(gbuffer);
		  );

	  if ((hit = Stage3end_new_terminal(/*querystart*/0,/*queryend*/querylength,
					    /*left*/segment_left,query_compress_fwd,
					    querylength,/*plusp*/true,genestrand,
					    segment->chrnum,segment->chroffset,
					    segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					    /*sarrayp*/false)) != NULL) {
	    debug4t(printf(" => Terminal\n"));
	    terminals = List_push(terminals,(void *) hit);
	  }
	}
      }
    }
  }

  if (minus_nsegments > 0) {
    for (segment = minus_segments; segment < &(minus_segments[minus_nsegments]); segment++) {
      if (0 && segment->usedp == true) {
	/* Previously skipped, but looks like a bad idea */
      } else if (segment->diagonal < (Univcoord_T) -1) {
	debug4t(printf("identify_terminals_minus: Getting genome at diagonal %llu (querypos %d..%d)\n",
		       (unsigned long long) segment->diagonal,segment->querypos5,segment->querypos3));

	if (segment->querypos3 - segment->querypos5 > index1part) {
	  segment_left = segment->diagonal - querylength;
	  debug4t(
		  gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		  Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		  printf("genome   0..: %s\n",gbuffer);
		  printf("query.rc 0..: %s\n",queryrc);
		  FREE(gbuffer);
		  );

	  if ((hit = Stage3end_new_terminal(/*querystart*/0,/*queryend*/querylength,
					    /*left*/segment_left,query_compress_rev,
					    querylength,/*plusp*/false,genestrand,
					    segment->chrnum,segment->chroffset,
					    segment->chrhigh,segment->chrlength,max_mismatches_allowed,
					    /*sarrayp*/false)) != NULL) {
	    debug4t(printf(" => Terminal\n"));
	    terminals = List_push(terminals,(void *) hit);
	  }
	}
      }
    }
  }

  debug4t(printf("Total number of terminals: %d\n",List_length(terminals)));
  fprintf(stderr,"Total number of terminals: %d\n",List_length(terminals));
  return terminals;
}
#endif



static void
fetch_positions_for_all_12mers (T this, Indexdb_T plus_indexdb, Indexdb_T minus_indexdb, 
				char *queryuc_ptr, int querylength, int query_lastpos) {
  int querypos;
  bool allvalidp;

  if (this->read_oligos_p == false) {
    read_oligos(&allvalidp,this,queryuc_ptr,querylength,query_lastpos,/*genestrand*/0);
  }

  /* querypos -2, -1, query_lastpos+1, and query_lastpos+2 are special cases */
  /* if allvalidp is true, then 0 and query_lastpos should have been done already */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->plus_retrievedp[querypos] == false) {
      /* FORMULA */
#ifdef LARGE_GENOMES
      this->plus_positions_low[querypos] =
	Indexdb_read_inplace(&(this->plus_npositions[querypos]),&(this->plus_positions_high[querypos]),
			     plus_indexdb,this->forward_oligos[querypos]);
#else
      this->plus_positions[querypos] =
	Indexdb_read_inplace(&(this->plus_npositions[querypos]),plus_indexdb,this->forward_oligos[querypos]);
#endif
      debug(printf("Retrieving at querypos %d, plus_npositions = %d\n",
		   querypos,this->plus_npositions[querypos]));
      this->plus_retrievedp[querypos] = true;
#ifdef USE_ALLOCP
      this->plus_allocp[querypos] = false;
#endif
    }
    if (this->minus_retrievedp[querypos] == false) {
      /* FORMULA */
#ifdef LARGE_GENOMES
      this->minus_positions_low[querypos] =
	Indexdb_read_inplace(&(this->minus_npositions[querypos]),&(this->minus_positions_high[querypos]),
			     minus_indexdb,this->revcomp_oligos[querypos]);
#else
      this->minus_positions[querypos] =
	Indexdb_read_inplace(&(this->minus_npositions[querypos]),minus_indexdb,this->revcomp_oligos[querypos]);
#endif
      debug(printf("Retrieving at querypos %d, minus_npositions = %d\n",
		   querypos,this->minus_npositions[querypos]));
      this->minus_retrievedp[querypos] = true;
#ifdef USE_ALLOCP
      this->minus_allocp[querypos] = false;
#endif
    }
  }

  this->all_positions_fetched_p = true;

  return;
}



static bool
intragenic_splice_p (Chrpos_T splicedistance, Substring_T donor, Substring_T acceptor) {
  int knowni;

  if ((knowni = Substring_splicesitesD_knowni(donor)) >= 0) {
    if (splicedists[knowni] >= splicedistance) {
      return true;
    }
  }

  if ((knowni = Substring_splicesitesA_knowni(acceptor)) >= 0) {
    if (splicedists[knowni] >= splicedistance) {
      return true;
    }
  }

  return false;
}



static List_T
find_splicepairs_distant_dna (int *found_score, int *ndistantsplicepairs,
			      List_T *localsplicing, List_T distantsplicing_orig,
			      List_T *startfrags_plus, List_T *endfrags_plus,
			      List_T *startfrags_minus, List_T *endfrags_minus,
			      int localsplicing_penalty, int distantsplicing_penalty,
			      int querylength, int nmismatches_allowed, bool first_read_p) {
  List_T distantsplicing = NULL, p, q, qsave;
  Substring_T startfrag, endfrag;
  int min_endlength_1, min_endlength_2, nmismatches1, nmismatches2, pos;
  Chrpos_T distance;
  Univcoord_T startfrag_genomicstart, endfrag_genomicstart;
  bool shortdistancep;
  double nonidentity = 1.0 - min_distantsplicing_identity;
  Chrnum_T chrnum;

  debug(printf("Starting find_splicepairs_distant_dna with nonidentity %f\n",nonidentity));
  debug4l(printf("Starting find_splicepairs_distant_dna with nonidentity %f\n",nonidentity));

  if (nonidentity == 0.0) {
    nmismatches_allowed = 0;
  }

  for (nmismatches1 = 0; nmismatches1 <= nmismatches_allowed; nmismatches1++) {
    nmismatches2 = nmismatches_allowed - nmismatches1;

    if (nonidentity == 0.0) {
      min_endlength_1 = min_endlength_2 = min_distantsplicing_end_matches;
    } else {
      min_endlength_1 = rint((double) nmismatches1/nonidentity);
      if (min_endlength_1 < min_distantsplicing_end_matches) {
	min_endlength_1 = min_distantsplicing_end_matches;
      }
      min_endlength_2 = rint((double) nmismatches2/nonidentity);
      if (min_endlength_2 < min_distantsplicing_end_matches) {
	min_endlength_2 = min_distantsplicing_end_matches;
      }
    }

    debug4l(printf("  nmismatches1 = %d, nmismatches2 = %d, min_endlength_1 = %d, min_endlength_2 = %d\n",
		   nmismatches1,nmismatches2,min_endlength_1,min_endlength_2));

    /************************************************************************
     *   Same strands
     ************************************************************************/

    /* 1.  End 1 to End 2.  Same strands. */
    p = startfrags_plus[nmismatches1];
    q = endfrags_plus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_dna (%d+%d mismatches): startfrags+ (%d) to endfrags+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      startfrag = (Substring_T) p->first;
      endfrag = (Substring_T) q->first;
      debug4ld(printf("end1-end2: startfrag at %llu %d..%d and endfrag at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(startfrag),
		      Substring_querystart(startfrag),Substring_queryend(startfrag),
		      (unsigned long long) Substring_genomicstart(endfrag),
		      Substring_querystart(endfrag),Substring_queryend(endfrag)));

      if ((pos = Substring_siteN_pos(startfrag)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of startfrag < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of startfrag > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteN_pos(endfrag)) {
	debug4ld(printf("chimera_pos of startfrag %d < chimera_pos of endfrag %d\n",pos,Substring_siteN_pos(endfrag)));
	p = p->rest;
      } else if (pos > Substring_siteN_pos(endfrag)) {
	debug4ld(printf("chimera_pos of startfrag %d > chimera_pos of endfrag %d\n",pos,Substring_siteN_pos(endfrag)));
	q = q->rest;
      } else {
	/* Generate all pairs at this splice_pos */
	qsave = q;
	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteN_pos(((Substring_T) p->first)) == pos) {
	  startfrag = (Substring_T) p->first;
	  debug4ld(printf("startfrag at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(startfrag),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteN_pos(((Substring_T) q->first)) == pos) {
	    endfrag = (Substring_T) q->first;
	    debug4ld(printf("endfrag at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(endfrag),pos));
	    if (Substring_genomicstart(endfrag) == Substring_genomicstart(startfrag)) {
	      /* Skip.  Really a continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(startfrag)) != Substring_chrnum(endfrag)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((endfrag_genomicstart = Substring_genomicstart(endfrag)) > (startfrag_genomicstart = Substring_genomicstart(startfrag))) {
		distance = endfrag_genomicstart - startfrag_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
		} else if (distances_observed_p == true &&
			   intragenic_splice_p(distance,startfrag,endfrag) == true) {
		  shortdistancep = true;
		} else {
		  shortdistancep = false;
		}
	      } else {
		distance = startfrag_genomicstart - endfrag_genomicstart;
		shortdistancep = false; /* scramble */
	      }
	      debug4ld(printf("1-2. Pushing a candidate at splice_pos %d (%d..%d), startfrag %llu to endfrag %llu.  shortdistancep = %d\n",
			      pos,min_endlength_1,querylength-min_endlength_2,
			      (unsigned long long) Substring_genomicstart(startfrag),
			      (unsigned long long) Substring_genomicstart(endfrag),shortdistancep));

	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									 startfrag,endfrag,/*donor_prob*/0.0,/*acceptor_prob*/0.0,distance,
									 /*shortdistancep*/true,localsplicing_penalty,querylength,
									 /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									 /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									 /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									 /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									 /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									 /*sensedir*/SENSE_NULL,/*sarrayp*/false));
	      } else if (*ndistantsplicepairs <= MAXCHIMERAPATHS) {
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									  startfrag,endfrag,/*donor_prob*/0.0,/*acceptor_prob*/0.0,distance,
									  /*shortdistancep*/false,distantsplicing_penalty,querylength,
									  /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									  /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									  /*sensedir*/SENSE_NULL,/*sarrayp*/false));
		(*ndistantsplicepairs)++;
	      }

	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }

    /* 4. End 3 to End 4.  Same strands. */
    p = startfrags_minus[nmismatches1];
    q = endfrags_minus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_dna (%d+%d mismatches): startfrags- (%d) to endfrags- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      startfrag = (Substring_T) p->first;
      endfrag = (Substring_T) q->first;
      debug4ld(printf("end3-end4: startfrag at %llu %d..%d and endfrag at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(startfrag),
		      Substring_querystart(startfrag),Substring_queryend(startfrag),
		      (unsigned long long) Substring_genomicstart(endfrag),
		      Substring_querystart(endfrag),Substring_queryend(endfrag)));

      if ((pos = Substring_siteN_pos(startfrag)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of startfrag < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of startfrag > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteN_pos(endfrag)) {
	debug4ld(printf("chimera_pos of startfrag %d < chimera_pos of endfrag %d\n",pos,Substring_siteN_pos(endfrag)));
	p = p->rest;
      } else if (pos > Substring_siteN_pos(endfrag)) {
	debug4ld(printf("chimera_pos of startfrag %d > chimera_pos of endfrag %d\n",pos,Substring_siteN_pos(endfrag)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteN_pos(((Substring_T) p->first)) == pos) {
	  startfrag = (Substring_T) p->first;
	  debug4ld(printf("startfrag at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(startfrag),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteN_pos(((Substring_T) q->first)) == pos) {
	    endfrag = (Substring_T) q->first;
	    debug4ld(printf("endfrag at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(endfrag),pos));
	    if (Substring_genomicstart(endfrag) == Substring_genomicstart(startfrag)) {
	      /* Skip.  Really a continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(startfrag)) != Substring_chrnum(endfrag)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((endfrag_genomicstart = Substring_genomicstart(endfrag)) > (startfrag_genomicstart = Substring_genomicstart(startfrag))) {
		distance = endfrag_genomicstart - startfrag_genomicstart;
		shortdistancep = false; /* scramble */
	      } else {
		distance = startfrag_genomicstart - endfrag_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
		} else if (distances_observed_p == true &&
			   intragenic_splice_p(distance,startfrag,endfrag) == true) {
		  shortdistancep = true;
		} else {
		  shortdistancep = false;
		}
	      }
	      debug4ld(printf("3-4. Pushing a candidate at splice_pos %d (%d..%d), startfrag %llu to endfrag %llu.  shortdistancep = %d.\n",
			      pos,min_endlength_1,querylength-min_endlength_2,
			      (unsigned long long) Substring_genomicstart(startfrag),
			      (unsigned long long) Substring_genomicstart(endfrag),shortdistancep));
	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									 startfrag,endfrag,/*donor_prob*/0.0,/*acceptor_prob*/0.0,distance,
									 /*shortdistancep*/true,localsplicing_penalty,querylength,
									 /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									 /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									 /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									 /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									 /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									 /*sensedir*/SENSE_NULL,/*sarrayp*/false));
	      } else if (*ndistantsplicepairs <= MAXCHIMERAPATHS) {
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									  startfrag,endfrag,/*donor_prob*/0.0,/*acceptor_prob*/0.0,distance,
									  /*shortdistancep*/false,distantsplicing_penalty,querylength,
									  /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									  /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									  /*sensedir*/SENSE_NULL,/*sarrayp*/false));
		(*ndistantsplicepairs)++;
	      }
	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }

    /* 5. End 5 to End 6.  Same strands. */
    /* 8. End 7 to End 8.  Same strands. */


    /************************************************************************
     *   Different strands
     ************************************************************************/

    /* 2. End 1 to End 4.  Different strands. */
    p = startfrags_plus[nmismatches1];
    q = endfrags_minus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_dna (%d+%d mismatches): startfrags+ (%d) to endfrags- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS) {
      startfrag = (Substring_T) p->first;
      endfrag = (Substring_T) q->first;
      debug4ld(printf("end1-end4: startfrag at %llu %d..%d and endfrag at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(startfrag),
		      Substring_querystart(startfrag),Substring_queryend(startfrag),
		      (unsigned long long) Substring_genomicstart(endfrag),
		      Substring_querystart(endfrag),Substring_queryend(endfrag)));

      if ((pos = Substring_siteN_pos(startfrag)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of startfrag < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of startfrag > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteN_pos(endfrag)) {
	debug4ld(printf("chimera_pos of startfrag %d < chimera_pos of endfrag %d\n",pos,Substring_siteN_pos(endfrag)));
	p = p->rest;
      } else if (pos > Substring_siteN_pos(endfrag)) {
	debug4ld(printf("chimera_pos of startfrag %d > chimera_pos of endfrag %d\n",pos,Substring_siteN_pos(endfrag)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteN_pos(((Substring_T) p->first)) == pos) {
	  startfrag = (Substring_T) p->first;
	  debug4ld(printf("startfrag at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(startfrag),pos));
	  q = qsave;
	  while (q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteN_pos(((Substring_T) q->first)) == pos) {
	    endfrag = (Substring_T) q->first;
	    debug4ld(printf("endfrag at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(endfrag),pos));
	    if (Substring_chrnum(startfrag) != Substring_chrnum(endfrag)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(endfrag) - pos) > (Substring_genomicstart(startfrag) + pos)) {
	      distance = (Substring_genomicstart(endfrag) - pos) - (Substring_genomicstart(startfrag) + pos);
	    } else {
	      distance = (Substring_genomicstart(startfrag) + pos) - (Substring_genomicstart(endfrag) - pos);
	    }
	    debug4ld(printf("1-4. Pushing a candidate at splice_pos %d (%d..%d), startfrag %llu to endfrag %llu.  Different strands, so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    (unsigned long long) Substring_genomicstart(startfrag),
			    (unsigned long long) Substring_genomicstart(endfrag)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      startfrag,endfrag,/*donor_prob*/0.0,/*acceptor_prob*/0.0,distance,
								      /*shortdistancep*/false,distantsplicing_penalty,querylength,
								      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								      /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_NULL,/*sarrayp*/false));
	    (*ndistantsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }

    /* 3. End 3 to End 2.  Different strands. */
    p = startfrags_minus[nmismatches1];
    q = endfrags_plus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_dna (%d+%d mismatches): startfrags- (%d) to endfrags+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS) {
      startfrag = (Substring_T) p->first;
      endfrag = (Substring_T) q->first;
      debug4ld(printf("end3-end2: startfrag at %llu %d..%d and endfrag at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(startfrag),
		      Substring_querystart(startfrag),Substring_queryend(startfrag),
		      (unsigned long long) Substring_genomicstart(endfrag),
		      Substring_querystart(endfrag),Substring_queryend(endfrag)));
		      
      if ((pos = Substring_siteN_pos(startfrag)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of startfrag < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of startfrag > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteN_pos(endfrag)) {
	debug4ld(printf("chimera_pos of startfrag %d < chimera_pos of endfrag %d\n",pos,Substring_siteN_pos(endfrag)));
	p = p->rest;
      } else if (pos > Substring_siteN_pos(endfrag)) {
	debug4ld(printf("chimera_pos of startfrag %d > chimera_pos of endfrag %d\n",pos,Substring_siteN_pos(endfrag)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteN_pos(((Substring_T) p->first)) == pos) {
	  startfrag = (Substring_T) p->first;
	  debug4ld(printf("startfrag at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(startfrag),pos));
	  q = qsave;
	  while (q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteN_pos(((Substring_T) q->first)) == pos) {
	    endfrag = (Substring_T) q->first;
	    debug4ld(printf("endfrag at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(endfrag),pos));
	    if (Substring_chrnum(startfrag) != Substring_chrnum(endfrag)) {
	      distance = 0U;
	    } else if (Substring_genomicstart(endfrag) > Substring_genomicstart(startfrag)) {
	      distance = (Substring_genomicstart(endfrag) + pos) - (Substring_genomicstart(startfrag) - pos);
	    } else {
	      distance = (Substring_genomicstart(startfrag) - pos) - (Substring_genomicstart(endfrag) + pos);
	    }
	    debug4ld(printf("3-2. Pushing a candidate at splice_pos %d (%d..%d), startfrag %llu to endfrag %llu.  Different strands so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    (unsigned long long) Substring_genomicstart(startfrag),
			    (unsigned long long) Substring_genomicstart(endfrag)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      startfrag,endfrag,/*donor_prob*/0.0,/*acceptor_prob*/0.0,distance,
								      /*shortdistancep*/false,distantsplicing_penalty,querylength,
								      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								      /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_NULL,/*sarrayp*/false));
	    (*ndistantsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }

    /* 6. End 5 to End 8.  Different strands. */
    /* 7. End 7 to End 6.  Different strands. */
  }

  debug4l(printf("ndistantsplicepairs %d, maxchimerapaths %d\n",*ndistantsplicepairs,MAXCHIMERAPATHS));
  if (*ndistantsplicepairs > MAXCHIMERAPATHS) {
    /* Can afford to ignore these if MAXCHIMERAPATHS is set high enough */
    stage3list_gc(&distantsplicing);
    return distantsplicing_orig;
  } else {
    return List_append(distantsplicing_orig,distantsplicing);
  }
}



static List_T
find_splicepairs_distant_rna (int *found_score, int *ndistantsplicepairs,
			      List_T *localsplicing, List_T distantsplicing_orig,
			      List_T *donors_plus, List_T *antidonors_plus,
			      List_T *acceptors_plus, List_T *antiacceptors_plus,
			      List_T *donors_minus, List_T *antidonors_minus,
			      List_T *acceptors_minus, List_T *antiacceptors_minus,
			      int localsplicing_penalty, int distantsplicing_penalty,
			      int querylength, int nmismatches_allowed, bool first_read_p) {
  List_T distantsplicing = NULL, p, q, qsave;
  Substring_T donor, acceptor;
  int min_endlength_1, min_endlength_2, nmismatches1, nmismatches2, pos;
  Chrpos_T distance;
  Univcoord_T donor_genomicstart, acceptor_genomicstart;
  bool shortdistancep;
  double nonidentity = 1.0 - min_distantsplicing_identity;
  Chrnum_T chrnum;

  debug(printf("Starting find_splicepairs_distant_rna with nonidentity %f\n",nonidentity));
  debug4l(printf("Starting find_splicepairs_distant_rna with nonidentity %f\n",nonidentity));

  if (nonidentity == 0.0) {
    nmismatches_allowed = 0;
  }

  for (nmismatches1 = 0; nmismatches1 <= nmismatches_allowed; nmismatches1++) {
    nmismatches2 = nmismatches_allowed - nmismatches1;

    if (nonidentity == 0.0) {
      min_endlength_1 = min_endlength_2 = min_distantsplicing_end_matches;
    } else {
      min_endlength_1 = rint((double) nmismatches1/nonidentity);
      if (min_endlength_1 < min_distantsplicing_end_matches) {
	min_endlength_1 = min_distantsplicing_end_matches;
      }
      min_endlength_2 = rint((double) nmismatches2/nonidentity);
      if (min_endlength_2 < min_distantsplicing_end_matches) {
	min_endlength_2 = min_distantsplicing_end_matches;
      }
    }

    debug4l(printf("  nmismatches1 = %d, nmismatches2 = %d, min_endlength_1 = %d, min_endlength_2 = %d\n",
		   nmismatches1,nmismatches2,min_endlength_1,min_endlength_2));

    /************************************************************************
     *   Same strands
     ************************************************************************/

    /* 1.  End 1 to End 2.  Same strands. */
    p = donors_plus[nmismatches1];
    q = acceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_rna (%d+%d mismatches): donors+ (%d) to acceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end1-end2: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of donor < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	/* Generate all pairs at this splice_pos */
	qsave = q;
	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(donor)) != Substring_chrnum(acceptor)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((acceptor_genomicstart = Substring_genomicstart(acceptor)) > (donor_genomicstart = Substring_genomicstart(donor))) {
		distance = acceptor_genomicstart - donor_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
		} else if (distances_observed_p == true &&
			   intragenic_splice_p(distance,donor,acceptor) == true) {
		  shortdistancep = true;
		} else {
		  shortdistancep = false;
		}
	      } else {
		distance = donor_genomicstart - acceptor_genomicstart;
		shortdistancep = false; /* scramble */
	      }
	      debug4ld(printf("1-2. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  shortdistancep = %d\n",
			      pos,min_endlength_1,querylength-min_endlength_2,
			      (unsigned long long) Substring_genomicstart(donor),
			      (unsigned long long) Substring_genomicstart(acceptor),shortdistancep));

	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									 donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
									 /*shortdistancep*/true,localsplicing_penalty,querylength,
									 /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									 /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									 /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									 /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									 /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									 /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      } else if (*ndistantsplicepairs <= MAXCHIMERAPATHS) {
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									  donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
									  /*shortdistancep*/false,distantsplicing_penalty,querylength,
									  /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									  /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									  /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
		(*ndistantsplicepairs)++;
	      }

	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }

    /* 4. End 3 to End 4.  Same strands. */
    p = donors_minus[nmismatches1];
    q = acceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_rna (%d+%d mismatches): donors- (%d) to acceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end3-end4: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of donor < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(donor)) != Substring_chrnum(acceptor)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((acceptor_genomicstart = Substring_genomicstart(acceptor)) > (donor_genomicstart = Substring_genomicstart(donor))) {
		distance = acceptor_genomicstart - donor_genomicstart;
		shortdistancep = false; /* scramble */
	      } else {
		distance = donor_genomicstart - acceptor_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
		} else if (distances_observed_p == true &&
			   intragenic_splice_p(distance,donor,acceptor) == true) {
		  shortdistancep = true;
		} else {
		  shortdistancep = false;
		}
	      }
	      debug4ld(printf("3-4. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  shortdistancep = %d.\n",
			      pos,min_endlength_1,querylength-min_endlength_2,
			      (unsigned long long) Substring_genomicstart(donor),
			      (unsigned long long) Substring_genomicstart(acceptor),shortdistancep));
	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									 donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
									 /*shortdistancep*/true,localsplicing_penalty,querylength,
									 /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									 /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									 /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									 /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									 /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									 /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      } else if (*ndistantsplicepairs <= MAXCHIMERAPATHS) {
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									  donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
									  /*shortdistancep*/false,distantsplicing_penalty,querylength,
									  /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									  /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									  /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
		(*ndistantsplicepairs)++;
	      }
	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }

    /* 5. End 5 to End 6.  Same strands. */
    p = antidonors_plus[nmismatches1];
    q = antiacceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_rna (%d+%d mismatches): antidonors+ (%d) to antiacceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end5-end6: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_2) {
	debug4ld(printf("chimera_pos of donor < min_endlength_2\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_1\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really an continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(donor)) != Substring_chrnum(acceptor)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((acceptor_genomicstart = Substring_genomicstart(acceptor)) > (donor_genomicstart = Substring_genomicstart(donor))) {
		distance = acceptor_genomicstart - donor_genomicstart;
		shortdistancep = false; /* scramble */
	      } else {
		distance = donor_genomicstart - acceptor_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
		} else if (distances_observed_p == true &&
			   intragenic_splice_p(distance,donor,acceptor) == true) {
		  shortdistancep = true;
		} else {
		  shortdistancep = false;
		}
	      }

	      debug4ld(printf("5-6. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  shortdistancep = %d\n",
			      pos,min_endlength_2,querylength-min_endlength_1,
			      (unsigned long long) Substring_genomicstart(donor),
			      (unsigned long long) Substring_genomicstart(acceptor),shortdistancep));
	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									 donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
									 /*shortdistancep*/true,localsplicing_penalty,querylength,
									 /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									 /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									 /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									 /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									 /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									 /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      } else if (*ndistantsplicepairs <= MAXCHIMERAPATHS) {
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									  donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
									  /*shortdistancep*/false,distantsplicing_penalty,querylength,
									  /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									  /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									  /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
		(*ndistantsplicepairs)++;
	      }
	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }

    /* 8. End 7 to End 8.  Same strands. */
    p = antidonors_minus[nmismatches1];
    q = antiacceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_rna (%d+%d mismatches): antidonors- (%d) to antiacceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end7-end8: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_2) {
	debug4ld(printf("chimera_pos of donor < min_endlength_2\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_1\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;

	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(donor)) != Substring_chrnum(acceptor)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((acceptor_genomicstart = Substring_genomicstart(acceptor)) > (donor_genomicstart = Substring_genomicstart(donor))) {
		distance = acceptor_genomicstart - donor_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
		} else if (distances_observed_p == true &&
			   intragenic_splice_p(distance,donor,acceptor) == true) {
		  shortdistancep = true;
		} else {
		  shortdistancep = false;
		}
	      } else {
		distance = donor_genomicstart - acceptor_genomicstart;
		shortdistancep = false; /* scramble */
	      }
	      debug4ld(printf("7-8. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  shortdistancep = %d.\n",
			      pos,min_endlength_2,querylength-min_endlength_1,
			      (unsigned long long) Substring_genomicstart(donor),
			      (unsigned long long) Substring_genomicstart(acceptor),shortdistancep));
	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									 donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
									 /*shortdistancep*/true,localsplicing_penalty,querylength,
									 /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									 /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									 /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									 /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									 /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									 /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      } else if (*ndistantsplicepairs <= MAXCHIMERAPATHS) {
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
									  donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
									  /*shortdistancep*/false,distantsplicing_penalty,querylength,
									  /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
									  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
									  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
									  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
									  /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
									  /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
		(*ndistantsplicepairs)++;
	      }
	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }


    /************************************************************************
     *   Different strands
     ************************************************************************/

    /* 2. End 1 to End 4.  Different strands. */
    p = donors_plus[nmismatches1];
    q = acceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_rna (%d+%d mismatches): donors+ (%d) to acceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end1-end4: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of donor < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) - pos) > (Substring_genomicstart(donor) + pos)) {
	      distance = (Substring_genomicstart(acceptor) - pos) - (Substring_genomicstart(donor) + pos);
	    } else {
	      distance = (Substring_genomicstart(donor) + pos) - (Substring_genomicstart(acceptor) - pos);
	    }
	    debug4ld(printf("1-4. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  Different strands, so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    (unsigned long long) Substring_genomicstart(donor),
			    (unsigned long long) Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
								      /*shortdistancep*/false,distantsplicing_penalty,querylength,
								      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								      /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	    (*ndistantsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }

    /* 3. End 3 to End 2.  Different strands. */
    p = donors_minus[nmismatches1];
    q = acceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_rna (%d+%d mismatches): donors- (%d) to acceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end3-end2: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of donor < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if (Substring_genomicstart(acceptor) > Substring_genomicstart(donor)) {
	      distance = (Substring_genomicstart(acceptor) + pos) - (Substring_genomicstart(donor) - pos);
	    } else {
	      distance = (Substring_genomicstart(donor) - pos) - (Substring_genomicstart(acceptor) + pos);
	    }
	    debug4ld(printf("3-2. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  Different strands so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    (unsigned long long) Substring_genomicstart(donor),
			    (unsigned long long) Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
								      /*shortdistancep*/false,distantsplicing_penalty,querylength,
								      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								      /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	    (*ndistantsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }


    /* 6. End 5 to End 8.  Different strands. */
    p = antidonors_plus[nmismatches1];
    q = antiacceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_rna (%d+%d mismatches): antidonors+ (%d) to antiacceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end5-end8: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_2) {
	debug4ld(printf("chimera_pos of donor < min_endlength_2\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_1\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) - pos) > (Substring_genomicstart(donor) + pos)) {
	      distance = (Substring_genomicstart(acceptor) - pos) - (Substring_genomicstart(donor) + pos);
	    } else {
	      distance = (Substring_genomicstart(donor) + pos) - (Substring_genomicstart(acceptor) - pos);
	    }
	    debug4ld(printf("5-8. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  Different strands so not shortdistance.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    (unsigned long long) Substring_genomicstart(donor),
			    (unsigned long long) Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
								      /*shortdistancep*/false,distantsplicing_penalty,querylength,
								      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								      /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	    (*ndistantsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }

    /* 7. End 7 to End 6.  Different strands. */
    p = antidonors_minus[nmismatches1];
    q = antiacceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_distant_rna (%d+%d mismatches): antidonors- (%d) to antiacceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end7-end6: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_2) {
	debug4ld(printf("chimera_pos of donor < min_endlength_2\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_1\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) + pos) > (Substring_genomicstart(donor) - pos)) {
	      distance = (Substring_genomicstart(acceptor) + pos) - (Substring_genomicstart(donor) - pos);
	    } else {
	      distance = (Substring_genomicstart(donor) - pos) - (Substring_genomicstart(acceptor) + pos);
	    }
	    debug4ld(printf("7-6. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  Different strands so not shortdistance.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    (unsigned long long) Substring_genomicstart(donor),
			    (unsigned long long) Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3end_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      donor,acceptor,Substring_siteD_prob(donor),Substring_siteA_prob(acceptor),distance,
								      /*shortdistancep*/false,distantsplicing_penalty,querylength,
								      /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								      /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								      /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								      /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	    (*ndistantsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  }

  debug4l(printf("ndistantsplicepairs %d, maxchimerapaths %d\n",*ndistantsplicepairs,MAXCHIMERAPATHS));
  if (*ndistantsplicepairs > MAXCHIMERAPATHS) {
    /* Can afford to ignore these if MAXCHIMERAPATHS is set high enough */
    stage3list_gc(&distantsplicing);
    return distantsplicing_orig;
  } else {
    return List_append(distantsplicing_orig,distantsplicing);
  }
}


static List_T
find_splicepairs_shortend (int *found_score, List_T hits,
			   List_T *donors_plus, List_T *antidonors_plus,
			   List_T *acceptors_plus, List_T *antiacceptors_plus,
			   List_T *donors_minus, List_T *antidonors_minus,
			   List_T *acceptors_minus, List_T *antiacceptors_minus,
			   Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   char *queryuc_ptr, char *queryrc, int min_shortend, int localsplicing_penalty,
			   int max_mismatches_allowed, int querylength, bool pairedp, bool first_read_p,
			   int genestrand) {
  List_T p;
  Substring_T donor, acceptor;
#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Intlist_T splicesites_i;
  Intlist_T nmismatches_list;
  Doublelist_T probs_list;
  int nmismatches, nmismatches_shortend, nmisses_allowed, support, endlength;
#ifdef DEBUG4H
  int amb_length;
  Univcoord_T leftbound, rightbound;
#endif
  Univcoord_T bestleft, origleft, chroffset, chrhigh;
  int i;
  int bestj = 0;


  debug(printf("Starting find_splicepairs_shortend\n"));
  debug(
	for (nmismatches = 0; nmismatches <= max_mismatches_allowed; nmismatches++) {
	  printf("At %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		 nmismatches,
		 List_length(donors_plus[nmismatches]),
		 List_length(acceptors_plus[nmismatches]),
		 List_length(antidonors_plus[nmismatches]),
		 List_length(antiacceptors_plus[nmismatches]),
		 List_length(donors_minus[nmismatches]),
		 List_length(acceptors_minus[nmismatches]),
		 List_length(antidonors_minus[nmismatches]),
		 List_length(antiacceptors_minus[nmismatches]));
	});

  /* Donors and antiacceptors => Want chimera_pos to be at end */
  /* Acceptors and antidonors => Want chimera_pos to be at beginning */

  /* Note: Previously checked endlength <=
     min_localsplicing_end_matches.  But this missed some ends.  Now
     just checking endlength <= support, to see if we are more than
     halfway. */

  /* Don't want to end when first set of hits found */
  for (nmismatches = 0; /* hits == NULL && */ nmismatches <= max_mismatches_allowed;
       nmismatches++) {
    nmisses_allowed = max_mismatches_allowed - nmismatches;

    /* End 1 */
    for (p = donors_plus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      support = Substring_siteD_pos(donor);
      endlength = querylength - support;
      chrhigh = Substring_chrhigh(donor);
      
#ifdef DEBUG4H
      chroffset = Substring_chroffset(donor);
      leftbound = Substring_alignend_trim(donor) + 1;
#endif
      debug4h(printf("End 1: short-overlap donor_plus: #%d:%u, endlength %d\n",
		     Substring_chrnum(donor),(Chrpos_T) (leftbound-1-chroffset),endlength));

      if (endlength <= support) {
	debug4h(printf("End 1: short-overlap donor_plus: #%d:%u (%d mismatches) => searching right\n",
		       Substring_chrnum(donor),(Chrpos_T) (leftbound-1-chroffset),Substring_nmismatches_whole(donor)));

	if ((i = Substring_splicesitesD_knowni(donor)) >= 0) {
	  origleft = Substring_genomicstart(donor);
	  if ((splicesites_i = 
	       Splicetrie_find_right(&nmismatches_shortend,&nmismatches_list,i,
				     origleft,/*pos5*/support,/*pos3*/querylength,chrhigh,
				     query_compress_fwd,/*queryptr*/queryuc_ptr,
				     nmisses_allowed,/*plusp*/true,genestrand,first_read_p,
				     /*collect_all_p*/pairedp == true && first_read_p == true)) != NULL) {
	    
	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      ambcoords = lookup_splicesites(&probs_list,splicesites_i,splicesites);
	      debug4h(amb_length = endlength /*- nmismatches_shortend*/);
	      debug4h(printf("End 1: short-overlap donor_plus: Successful ambiguous from donor #%d with amb_length %d\n",
			     Substring_splicesitesD_knowni(donor),amb_length));
	      hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								  donor,/*acceptor*/NULL,Substring_siteD_prob(donor),Doublelist_max(probs_list),/*distance*/0U,
								  /*shortdistancep*/false,/*penalty*/0,querylength,
								  /*ambcoords_donor*/NULL,ambcoords,
								  /*ambi_donor*/NULL,/*ambi_acceptor*/splicesites_i,
								  /*amb_nmismatches_donor*/NULL,/*nmismatches_acceptor*/nmismatches_list,
								  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/probs_list,
								  /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								  /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      Doublelist_free(&probs_list);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - support;
	      if ((acceptor = Substring_new_acceptor(/*acceptor_coord*/splicesites[bestj],/*acceptor_knowni*/bestj,
						     Substring_siteD_pos(donor),/*substring_querystart*/0,/*substring_queryend*/querylength,
						     nmismatches_shortend,/*prob*/2.0,/*left*/bestleft,query_compress_fwd,
						     querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
						     Substring_chrnum(donor),Substring_chroffset(donor),
						     Substring_chrhigh(donor),Substring_chrlength(donor))) != NULL) {
		debug4h(printf("End 1: short-overlap donor_plus: Successful splice from donor #%d to acceptor #%d\n",
			       Substring_splicesitesD_knowni(donor),Substring_splicesitesA_knowni(acceptor)));
		hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								    donor,acceptor,Substring_siteD_prob(donor),/*acceptor_prob*/2.0,/*distance*/bestleft-origleft,
								    /*shortdistancep*/true,localsplicing_penalty,querylength,
								    /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								    /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								    /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								    /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								    /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								    /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 2 */
    for (p = acceptors_plus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      endlength = Substring_siteA_pos(acceptor);
      support = querylength - endlength;
      chroffset = Substring_chroffset(acceptor);

#ifdef DEBUG4H
      rightbound = Substring_alignstart_trim(acceptor);
#endif

      debug4h(printf("End 2: short-overlap acceptor_plus: #%d:%u, endlength %d\n",
		     Substring_chrnum(acceptor),(Chrpos_T) (rightbound+1-chroffset),endlength));

      if (endlength <= support) {
	debug4h(printf("End 2: short-overlap acceptor_plus: #%d:%u (%d mismatches) => searching left\n",
		       Substring_chrnum(acceptor),(Chrpos_T) (rightbound+1-chroffset),Substring_nmismatches_whole(acceptor)));

	if ((i = Substring_splicesitesA_knowni(acceptor)) >= 0) {
	  origleft = Substring_genomicstart(acceptor);
	  if ((splicesites_i =
	       Splicetrie_find_left(&nmismatches_shortend,&nmismatches_list,i,
				    origleft,/*pos5*/0,/*pos3*/endlength,chroffset,
				    query_compress_fwd,/*queryptr*/queryuc_ptr,querylength,
				    nmisses_allowed,/*plusp*/true,genestrand,first_read_p,
				    /*collect_all_p*/pairedp == true && first_read_p == false)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      ambcoords = lookup_splicesites(&probs_list,splicesites_i,splicesites);
	      debug4h(amb_length = endlength /*- nmismatches_shortend*/);
	      debug4h(printf("End 2: short-overlap acceptor_plus: Successful ambiguous from acceptor #%d with amb_length %d\n",
			     Substring_splicesitesA_knowni(acceptor),amb_length));
	      hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								  /*donor*/NULL,acceptor,Doublelist_max(probs_list),Substring_siteA_prob(acceptor),/*distance*/0U,
								  /*shortdistancep*/false,/*penalty*/0,querylength,
								  ambcoords,/*ambcoords_acceptor*/NULL,
								  /*amb_knowni_donor*/splicesites_i,/*amb_knowni_acceptor*/NULL,
								  /*amb_nmismatches_donor*/nmismatches_list,/*amb_nmismatches_acceptor*/NULL,
								  /*amb_probs_donor*/probs_list,/*amb_probs_acceptor*/NULL,
								  /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								  /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      Doublelist_free(&probs_list);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - endlength;
	      if ((donor = Substring_new_donor(/*donor_coord*/splicesites[bestj],/*donor_knowni*/bestj,
					       Substring_siteA_pos(acceptor),/*substring_querystart*/0,/*substring_queryend*/querylength,
					       nmismatches_shortend,/*prob*/2.0,/*left*/bestleft,query_compress_fwd,
					       querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
					       Substring_chrnum(acceptor),Substring_chroffset(acceptor),
					       Substring_chrhigh(acceptor),Substring_chrlength(acceptor))) != NULL) {
		debug4h(printf("End 2: short-overlap acceptor_plus: Successful splice from acceptor #%d to donor #%d\n",
			       Substring_splicesitesA_knowni(acceptor),Substring_splicesitesD_knowni(donor)));
		hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								    donor,acceptor,/*donor_prob*/2.0,Substring_siteA_prob(acceptor),/*distance*/origleft-bestleft,
								    /*shortdistancep*/true,localsplicing_penalty,querylength,
								    /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								    /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								    /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								    /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								    /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								    /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 3 */
    for (p = donors_minus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      support = Substring_siteD_pos(donor);
      endlength = querylength - support;
      chroffset = Substring_chroffset(donor);

#ifdef DEBUG4H
      rightbound = Substring_alignend_trim(donor);
#endif

      debug4h(printf("End 3: short-overlap donor_minus: #%d:%u, endlength %d\n",
		     Substring_chrnum(donor),(Chrpos_T) (rightbound+1-chroffset),endlength));

      if (endlength <= support) {
	debug4h(printf("End 3: short-overlap donor_minus: #%d:%u (%d mismatches) => searching left\n",
		       Substring_chrnum(donor),(Chrpos_T) (rightbound+1-chroffset),Substring_nmismatches_whole(donor)));

	if ((i = Substring_splicesitesD_knowni(donor)) >= 0) {
	  origleft = Substring_genomicend(donor);
	  if ((splicesites_i =
	       Splicetrie_find_left(&nmismatches_shortend,&nmismatches_list,i,
				    origleft,/*pos5*/0,/*pos3*/endlength,chroffset,
				    query_compress_rev,/*queryptr*/queryrc,querylength,
				    nmisses_allowed,/*plusp*/false,genestrand,first_read_p,
				    /*collect_all_p*/pairedp == true && first_read_p == true)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      ambcoords = lookup_splicesites(&probs_list,splicesites_i,splicesites);
	      debug4h(amb_length = endlength /*- nmismatches_shortend*/);
	      debug4h(printf("End 3: short-overlap donor_minus: Successful ambiguous from donor #%d with amb_length %d\n",
			     Substring_splicesitesD_knowni(donor),amb_length));
	      hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								  donor,/*acceptor*/NULL,Substring_siteD_prob(donor),Doublelist_max(probs_list),/*distance*/0U,
								  /*shortdistancep*/false,/*penalty*/0,querylength,
								  /*ambcoords_donor*/NULL,ambcoords,
								  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/splicesites_i,
								  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/nmismatches_list,
								  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/probs_list,
								  /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								  /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      Doublelist_free(&probs_list);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - endlength;
	      if ((acceptor = Substring_new_acceptor(/*acceptor_coord*/splicesites[bestj],/*acceptor_knowni*/bestj,
						     querylength-Substring_siteD_pos(donor),
						     /*substring_querystart*/0,/*substring_queryend*/querylength,
						     nmismatches_shortend,/*prob*/2.0,/*left*/bestleft,query_compress_rev,
						     querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
						     Substring_chrnum(donor),Substring_chroffset(donor),
						     Substring_chrhigh(donor),Substring_chrlength(donor))) != NULL) {
		debug4h(printf("End 3: short-overlap donor_minus: Successful splice from donor #%d to acceptor #%d\n",
			       Substring_splicesitesD_knowni(donor),Substring_splicesitesA_knowni(acceptor)));
		hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								    donor,acceptor,Substring_siteD_prob(donor),/*acceptor_prob*/2.0,/*distance*/origleft-bestleft,
								    /*shortdistancep*/true,localsplicing_penalty,querylength,
								    /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								    /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								    /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								    /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								    /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								    /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 4 */
    for (p = acceptors_minus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      endlength = Substring_siteA_pos(acceptor);
      support = querylength - endlength;
      chrhigh = Substring_chrhigh(acceptor);

#ifdef DEBUG4H
      chroffset = Substring_chroffset(acceptor);
      leftbound = Substring_alignstart_trim(acceptor) + 1;
#endif

      debug4h(printf("End 4: short-overlap acceptor_minus: #%d:%u, endlength %d\n",
		     Substring_chrnum(acceptor),(Chrpos_T) (leftbound-1-chroffset),endlength));

      if (endlength <= support) {
	debug4h(printf("End 4: short-overlap acceptor_minus: #%d:%u (%d mismatches) => searching right\n",
		       Substring_chrnum(acceptor),(Chrpos_T) (leftbound-1-chroffset),Substring_nmismatches_whole(acceptor)));

	if ((i = Substring_splicesitesA_knowni(acceptor)) >= 0) {
	  origleft = Substring_genomicend(acceptor);
	  if ((splicesites_i =
	       Splicetrie_find_right(&nmismatches_shortend,&nmismatches_list,i,
				     origleft,/*pos5*/support,/*pos3*/querylength,chrhigh,
				     query_compress_rev,/*queryptr*/queryrc,
				     nmisses_allowed,/*plusp*/false,genestrand,first_read_p,
				     /*collect_all_p*/pairedp == true && first_read_p == false)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      ambcoords = lookup_splicesites(&probs_list,splicesites_i,splicesites);
	      debug4h(amb_length = endlength /*- nmismatches_shortend*/);
	      debug4h(printf("End 4: short-overlap acceptor_minus: Successful ambiguous from acceptor #%d with amb_length %d\n",
			     Substring_splicesitesA_knowni(acceptor),amb_length));
	      hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								  /*donor*/NULL,acceptor,Doublelist_max(probs_list),Substring_siteA_prob(acceptor),/*distance*/0U,
								  /*shortdistancep*/false,/*penalty*/0,querylength,
								  ambcoords,/*ambcoords_acceptor*/NULL,
								  /*amb_knowni_donor*/splicesites_i,/*amb_knowni_acceptor*/NULL,
								  /*amb_nmismatches_donor*/nmismatches_list,/*amb_nmismatches_acceptor*/NULL,
								  /*amb_probs_donor*/probs_list,/*amb_probs_acceptor*/NULL,
								  /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								  /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      Doublelist_free(&probs_list);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - support;
	      if ((donor = Substring_new_donor(/*donor_coord*/splicesites[bestj],/*donor_knowni*/bestj,
					       querylength-Substring_siteA_pos(acceptor),
					       /*substring_querystart*/0,/*substring_queryend*/querylength,
					       nmismatches_shortend,/*prob*/2.0,/*left*/bestleft,query_compress_rev,
					       querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
					       Substring_chrnum(acceptor),Substring_chroffset(acceptor),
					       Substring_chrhigh(acceptor),Substring_chrlength(acceptor))) != NULL) {
		debug4h(printf("End 4: short-overlap acceptor_minus: Successful splice from acceptor #%d to #%d\n",
			       Substring_splicesitesA_knowni(acceptor),Substring_splicesitesD_knowni(donor)));
		hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								    donor,acceptor,/*donor_prob*/2.0,Substring_siteA_prob(acceptor),/*distance*/bestleft-origleft,
								    /*shortdistancep*/true,localsplicing_penalty,querylength,
								    /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								    /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								    /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								    /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								    /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								    /*sensedir*/SENSE_FORWARD,/*sarrayp*/false));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 5 */
    for (p = antidonors_plus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      endlength = Substring_siteD_pos(donor);
      support = querylength - endlength;
      chroffset = Substring_chroffset(donor);

#ifdef DEBUG4H
      rightbound = Substring_alignstart_trim(donor);
#endif

      debug4h(printf("End 5: short-overlap antidonor_plus: #%d:%u, endlength %d\n",
		     Substring_chrnum(donor),(Chrpos_T) (rightbound+1-chroffset),endlength));

      if (endlength <= support) {
	debug4h(printf("End 5: short-overlap antidonor_plus: #%d:%u (%d mismatches) => searching left\n",
		       Substring_chrnum(donor),(Chrpos_T) (rightbound+1-chroffset),Substring_nmismatches_whole(donor)));

	if ((i = Substring_splicesitesD_knowni(donor)) >= 0) {
	  origleft = Substring_genomicstart(donor);
	  if ((splicesites_i =
	       Splicetrie_find_left(&nmismatches_shortend,&nmismatches_list,i,
				    origleft,/*pos5*/0,/*pos3*/endlength,chroffset,
				    query_compress_fwd,/*queryptr*/queryuc_ptr,querylength,
				    nmisses_allowed,/*plusp*/true,genestrand,first_read_p,
				    /*collect_all_p*/pairedp == true && first_read_p == false)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      ambcoords = lookup_splicesites(&probs_list,splicesites_i,splicesites);
	      debug4h(amb_length = endlength /*- nmismatches_shortend*/);
	      debug4h(printf("End 5: short-overlap antidonor_plus: Successful ambiguous from antidonor #%d with amb_length %d\n",
			     Substring_splicesitesD_knowni(donor),amb_length));
	      hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								  donor,/*acceptor*/NULL,Substring_siteD_prob(donor),Doublelist_max(probs_list),/*distance*/0U,
								  /*shortdistancep*/false,/*penalty*/0,querylength,
								  /*ambcoords_donor*/NULL,ambcoords,
								  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/splicesites_i,
								  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/nmismatches_list,
								  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/probs_list,
								  /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								  /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      Doublelist_free(&probs_list);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - endlength;
	      if ((acceptor = Substring_new_acceptor(/*acceptor_coord*/splicesites[bestj],/*acceptor_knowni*/bestj,
						     Substring_siteD_pos(donor),
						     /*substring_querystart*/0,/*substring_queryend*/querylength,
						     nmismatches_shortend,/*prob*/2.0,/*left*/bestleft,query_compress_fwd,
						     querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
						     Substring_chrnum(donor),Substring_chroffset(donor),
						     Substring_chrhigh(donor),Substring_chrlength(donor))) != NULL) {
		debug4h(printf("End 5: short-overlap antidonor_plus: Successful splice from antidonor #%d to antiacceptor #%d\n",
			       Substring_splicesitesD_knowni(donor),Substring_splicesitesA_knowni(acceptor)));
		hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								    donor,acceptor,Substring_siteD_prob(donor),/*acceptor_prob*/2.0,/*distance*/origleft-bestleft,
								    /*shortdistancep*/true,localsplicing_penalty,querylength,
								    /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								    /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								    /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								    /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								    /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								    /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 6 */
    for (p = antiacceptors_plus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      support = Substring_siteA_pos(acceptor);
      endlength = querylength - support;
      chrhigh = Substring_chrhigh(acceptor);

#ifdef DEBUG4H
      chroffset = Substring_chroffset(acceptor);
      leftbound = Substring_alignend_trim(acceptor) + 1;
#endif

      debug4h(printf("End 6: short-overlap antiacceptor_plus: #%d:%u, endlength %d\n",
		     Substring_chrnum(acceptor),(Chrpos_T) (leftbound-1-chroffset),endlength));

      if (endlength <= support) {
	debug4h(printf("End 6: short-overlap antiacceptor_plus: #%d:%u (%d mismatches) => searching right\n",
		       Substring_chrnum(acceptor),(Chrpos_T) (leftbound-1-chroffset),Substring_nmismatches_whole(acceptor)));

	if ((i = Substring_splicesitesA_knowni(acceptor)) >= 0) {
	  origleft = Substring_genomicstart(acceptor);
	  if ((splicesites_i =
	       Splicetrie_find_right(&nmismatches_shortend,&nmismatches_list,i,
				     origleft,/*pos5*/support,/*pos3*/querylength,chrhigh,
				     query_compress_fwd,/*queryptr*/queryuc_ptr,
				     nmisses_allowed,/*plusp*/true,genestrand,first_read_p,
				     /*collect_all_p*/pairedp == true && first_read_p == true)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      ambcoords = lookup_splicesites(&probs_list,splicesites_i,splicesites);
	      debug4h(amb_length = endlength /*- nmismatches_shortend*/);
	      debug4h(printf("End 6: short-overlap antiacceptor_plus: Successful ambiguous from antiacceptor #%d with amb_length %d\n",
			     Substring_splicesitesA_knowni(acceptor),amb_length));
	      hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								  /*donor*/NULL,acceptor,Doublelist_max(probs_list),Substring_siteA_prob(acceptor),/*distance*/0U,
								  /*shortdistancep*/false,/*penalty*/0,querylength,
								  ambcoords,/*ambcoords_acceptor*/NULL,
								  /*amb_knowni_donor*/splicesites_i,/*amb_knowni_acceptor*/NULL,
								  /*amb_nmismatches_donor*/nmismatches_list,/*amb_nmismatches_acceptor*/NULL,
								  /*amb_probs_donor*/probs_list,/*amb_probs_acceptor*/NULL,
								  /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								  /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      Doublelist_free(&probs_list);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - support;
	      if ((donor = Substring_new_donor(/*donor_coord*/splicesites[bestj],/*donor_knowni*/bestj,
					       Substring_siteA_pos(acceptor),
					       /*substring_querystart*/0,/*substring_queryend*/querylength,
					       nmismatches_shortend,/*prob*/2.0,/*left*/bestleft,query_compress_fwd,
					       querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
					       Substring_chrnum(acceptor),Substring_chroffset(acceptor),
					       Substring_chrhigh(acceptor),Substring_chrlength(acceptor))) != NULL) {
		debug4h(printf("End 6: short-overlap antiacceptor_plus: Successful splice from antiacceptor #%d to antidonor #%d\n",
			       Substring_splicesitesA_knowni(acceptor),Substring_splicesitesD_knowni(donor)));
		hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								    donor,acceptor,/*donor_prob*/2.0,Substring_siteA_prob(acceptor),/*distance*/bestleft-origleft,
								    /*shortdistancep*/true,localsplicing_penalty,querylength,
								    /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								    /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								    /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								    /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								    /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								    /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 7 */
    for (p = antidonors_minus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      endlength = Substring_siteD_pos(donor);
      support = querylength - endlength;
      chrhigh = Substring_chrhigh(donor);

#ifdef DEBUG4H
      chroffset = Substring_chroffset(donor);
      leftbound = Substring_alignstart_trim(donor) + 1;
#endif

      debug4h(printf("End 7: short-overlap antidonor_minus: #%d:%u, endlength %d\n",
		     Substring_chrnum(donor),(Chrpos_T) (leftbound-1-chroffset),endlength));

      if (endlength <= support) {
	debug4h(printf("End 7: short-overlap antidonor_minus: #%d:%u (%d mismatches) => searching right\n",
		       Substring_chrnum(donor),(Chrpos_T) (leftbound-1-chroffset),Substring_nmismatches_whole(donor)));

	if ((i = Substring_splicesitesD_knowni(donor)) >= 0) {
	  origleft = Substring_genomicend(donor);
	  if ((splicesites_i =
	       Splicetrie_find_right(&nmismatches_shortend,&nmismatches_list,i,
				     origleft,/*pos5*/support,/*pos3*/querylength,chrhigh,
				     query_compress_rev,/*queryptr*/queryrc,
				     nmisses_allowed,/*plusp*/false,genestrand,first_read_p,
				     /*collect_all_p*/pairedp == true && first_read_p == false)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      ambcoords = lookup_splicesites(&probs_list,splicesites_i,splicesites);
	      debug4h(amb_length = endlength /*- nmismatches_shortend*/);
	      debug4h(printf("End 7: short-overlap antidonor_minus: Successful ambiguous from antidonor #%d with amb_length %d\n",
			     Substring_splicesitesD_knowni(donor),amb_length));
	      hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								  donor,/*acceptor*/NULL,Substring_siteD_prob(donor),Doublelist_max(probs_list),/*distance*/0U,
								  /*shortdistancep*/false,/*penalty*/0,querylength,
								  /*ambcoords_donor*/NULL,ambcoords,
								  /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/splicesites_i,
								  /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/nmismatches_list,
								  /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/probs_list,
								  /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								  /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      Doublelist_free(&probs_list);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - support;
	      if ((acceptor = Substring_new_acceptor(/*acceptor_coord*/splicesites[bestj],/*acceptor_knowni*/bestj,
						     querylength-Substring_siteD_pos(donor),
						     /*substring_querystart*/0,/*substring_queryend*/querylength,
						     nmismatches_shortend,/*prob*/2.0,/*left*/bestleft,query_compress_rev,
						     querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
						     Substring_chrnum(donor),Substring_chroffset(donor),
						     Substring_chrhigh(donor),Substring_chrlength(donor))) != NULL) {
		debug4h(printf("End 7: short-overlap antidonor_minus: Successful splice from antidonor #%d to antiacceptor #%d\n",
			       Substring_splicesitesD_knowni(donor),Substring_splicesitesA_knowni(acceptor)));
		hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								    donor,acceptor,Substring_siteD_prob(donor),/*acceptor_prob*/2.0,/*distance*/bestleft-origleft,
								    /*shortdistancep*/true,localsplicing_penalty,querylength,
								    /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								    /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								    /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								    /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								    /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								    /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 8 */
    for (p = antiacceptors_minus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      support = Substring_siteA_pos(acceptor);
      endlength = querylength - support;
      chroffset = Substring_chroffset(acceptor);

#ifdef DEBUG4H
      rightbound = Substring_alignend_trim(acceptor);
#endif

      debug4h(printf("End 8: short-overlap antiacceptor_minus: #%d:%u, endlength %d\n",
		     Substring_chrnum(acceptor),(Chrpos_T) (rightbound+1-chroffset),endlength));

      if (endlength <= support) {
	debug4h(printf("End 8: short-overlap antiacceptor_minus: #%d:%u (%d mismatches) => searching left\n",
		       Substring_chrnum(acceptor),(Chrpos_T) (rightbound+1-chroffset),Substring_nmismatches_whole(acceptor)));

	if ((i = Substring_splicesitesA_knowni(acceptor)) >= 0) {
	  origleft = Substring_genomicend(acceptor);
	  if ((splicesites_i =
	       Splicetrie_find_left(&nmismatches_shortend,&nmismatches_list,i,
				    origleft,/*pos5*/0,/*pos3*/endlength,chroffset,
				    query_compress_rev,/*queryptr*/queryrc,querylength,
				    nmisses_allowed,/*plusp*/false,genestrand,first_read_p,
				    /*collect_all_p*/pairedp == true && first_read_p == true)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      ambcoords = lookup_splicesites(&probs_list,splicesites_i,splicesites);
	      debug4h(amb_length = endlength /*- nmismatches_shortend*/);
	      debug4h(printf("End 8: short-overlap antiacceptor_minus: Successful ambiguous from antiacceptor #%d with amb_length %d\n",
			     Substring_splicesitesA_knowni(acceptor),amb_length));
	      hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								  /*donor*/NULL,acceptor,Doublelist_max(probs_list),Substring_siteA_prob(acceptor),/*distance*/0U,
								  /*shortdistancep*/false,/*penalty*/0,querylength,
								  ambcoords,/*ambcoords_acceptor*/NULL,
								  /*amb_knowni_donor*/splicesites_i,/*amb_knowni_acceptor*/NULL,
								  /*amb_nmismatches_donor*/nmismatches_list,/*amb_nmismatches_acceptor*/NULL,
								  /*amb_probs_donor*/probs_list,/*amb_probs_acceptor*/NULL,
								  /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								  /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      Doublelist_free(&probs_list);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - endlength;
	      if ((donor = Substring_new_donor(/*donor_coord*/splicesites[bestj],/*donor_knowni*/bestj,
					       querylength-Substring_siteA_pos(acceptor),
					       /*substring_querystart*/0,/*substring_queryend*/querylength,
					       nmismatches_shortend,/*prob*/2.0,/*left*/bestleft,query_compress_rev,
					       querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
					       Substring_chrnum(acceptor),Substring_chroffset(acceptor),
					       Substring_chrhigh(acceptor),Substring_chrlength(acceptor))) != NULL) {
		debug4h(printf("End 8: short-overlap antiacceptor_minus: Successful splice from antiacceptor #%d to antidonor #%d\n",
			       Substring_splicesitesA_knowni(acceptor),Substring_splicesitesD_knowni(donor)));
		hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								    donor,acceptor,/*donor_prob*/2.0,Substring_siteA_prob(acceptor),/*distance*/origleft-bestleft,
								    /*shortdistancep*/true,localsplicing_penalty,querylength,
								    /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								    /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								    /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								    /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								    /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								    /*sensedir*/SENSE_ANTI,/*sarrayp*/false));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

  }
  debug(printf("Ending find_splicepairs_shortend\n"));

  return hits;
}



static void
find_12mer_bounds (int *firstbound, int *lastbound, bool *omitted, int query_lastpos) {
  int nconsecutive;
  int querypos;

  querypos = 0;
  nconsecutive = 0;
  while (nconsecutive < index1interval && querypos <= query_lastpos) {
    if (omitted[querypos]) {
      nconsecutive = 0;
    } else {
      nconsecutive++;
    }
    querypos++;
  }
  if (nconsecutive < index1interval) {
    *firstbound = 0;
    /* *lastbound = query_lastpos; */
    debug(printf("nconsecutive from left < 3, so setting firstbound 0\n"));
    /* return false; */
  } else {
    *firstbound = querypos + (index1interval-1) + index1part; /* From trial-and-error, this is the correct value */
    debug(printf("Assigning firstbound to be querypos %d + (%d - 1) + index1part = %d\n",
		 querypos,index1interval,*firstbound));
  }

  querypos = query_lastpos;
  nconsecutive = 0;
  while (nconsecutive < index1interval && querypos >= 0) {
    if (omitted[querypos]) {
      nconsecutive = 0;
    } else {
      nconsecutive++;
    }
    querypos--;
  }

  if (nconsecutive < index1interval) {
    /* *firstbound = 0 + index1part; */
    *lastbound = query_lastpos;
    debug(printf("nconsecutive from right < %d, so setting lastbound %d\n",
		 index1interval,query_lastpos));
    /* return false; */
  } else {
    *lastbound = querypos-1; /* From trial-and-error, this is the correct value */
    debug(printf("Assigning lastbound to be querypos %d - 1 = %d\n",
		 querypos,*lastbound));
  }

  return;

  /* The following extensions were causing misses on test set */
#if 0
  if (*firstbound > query_lastpos) {
    debug(printf("firstbound %d > query_lastpos %d, so setting firstbound 0 and lastbound %d\n",
		 *firstbound,query_lastpos,query_lastpos));
    *firstbound = 0;
    *lastbound = query_lastpos;
    return false;
#if 0
  } else if (*firstbound - index1part > *lastbound) {
    *firstbound = 0;
    *lastbound = query_lastpos;
    return false;
#endif
  } else if (*lastbound <= index1part) {
    debug(printf("lastbound %d <= %d, so setting firstbound 0 and lastbound %d\n",
		 *lastbound,index1part,query_lastpos));
    *firstbound = 0;
    *lastbound = query_lastpos;
    return false;
  } else {
    return true;
  }
#endif
}



static Floors_T
compute_floors (bool *any_omitted_p, bool *alloc_floors_p, Floors_T *floors_array,
		T this, char *queryuc_ptr, int querylength, int query_lastpos,
		Indexdb_T plus_indexdb, Indexdb_T minus_indexdb,
		int indexdb_size_threshold, int max_end_insertions,
		bool omit_frequent_p, bool omit_repetitive_p, bool keep_floors_p) {
  Floors_T floors;
  bool all_omitted_p;


  if (this->all_positions_fetched_p == true) {
    omit_oligos_clear(this,query_lastpos);
  } else {
    fetch_positions_for_all_12mers(this,plus_indexdb,minus_indexdb,
				   queryuc_ptr,querylength,query_lastpos);
  }

  debug(printf("Omitting frequent/repetitive oligos\n"));
  omit_oligos(&all_omitted_p,&(*any_omitted_p),this,query_lastpos,indexdb_size_threshold,
	      omit_frequent_p,omit_repetitive_p);

  if (all_omitted_p == true) {
    debug(printf("Aborting because all oligos are omitted\n"));
    *alloc_floors_p = false;
    return (Floors_T) NULL;
  } else if (*any_omitted_p) {
    floors = Floors_new_omitted(querylength,max_end_insertions,this->omitted);
    *alloc_floors_p = true;
  } else if (querylength > max_floors_readlength) {
    floors = Floors_new_standard(querylength,max_end_insertions,/*keep_floors_p*/false);
    *alloc_floors_p = true;
  } else if (keep_floors_p == false) {
    floors = Floors_new_standard(querylength,max_end_insertions,/*keep_floors_p*/false);
    *alloc_floors_p = true;
  } else {
    if (floors_array[querylength] == NULL) {
      floors_array[querylength] = Floors_new_standard(querylength,max_end_insertions,/*keep_floors_p*/true);
    }
    floors = floors_array[querylength];
    *alloc_floors_p = false;
  }

  return floors;
}


static void
complete_set_mm_indels (int *found_score, bool *segments_computed_p,
			Segment_T **plus_anchor_segments, Segment_T **minus_anchor_segments,
			int *n_plus_anchors, int *n_minus_anchors,
			int *opt_level, int *done_level, int user_maxlevel,
			bool revise_levels_p, int *nhits, List_T *subs, List_T *indels, T this,
			Compress_T query_compress_fwd, Compress_T query_compress_rev,
#if defined(DEBUG2) || defined(DEBUG2E)
			char *queryuc_ptr, char *queryrc,
#endif
			int querylength, int query_lastpos, Floors_T floors,
			int indel_penalty_middle, int indel_penalty_end,
			bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
			int fast_level, int genestrand) {
  int firstbound, lastbound;
  int max_mismatches_allowed;
#if 0
  int indel_level;
#endif

  debug(printf("Starting complete_set_mm_indels with found_score %d\n",*found_score));

  this->plus_segments = NULL;
  this->minus_segments = NULL;

  /* 4 and 5. Mismatches and indels via complete set.  Requires compress and
     all positions fetched.  Omits oligos and creates segments for some diagonals. */

#if 0
  if (find_12mer_bounds(&firstbound,&lastbound,this->omitted,query_lastpos) == false) {
    debug(printf("Cannot find 12_mer bounds\n"));
    /* This was allowing end indels to be missed */
    /* allow_end_indels_p = false; */
  } else {
    debug(printf("Found firstbound %d and lastbound %d\n",firstbound,lastbound));
  }
#else
  find_12mer_bounds(&firstbound,&lastbound,this->omitted,query_lastpos);
#endif

  if (*done_level < indel_penalty_end) {
    /* Prevents accumulation of segments for indels */
    allow_end_indels_p = false;
  }

  /* 4. Complete set mismatches */
  /* Done as a single batch */
  max_mismatches_allowed = (*done_level <= fast_level) ? -1 : *done_level;
  debug(printf("*** Stage 5.  Complete set mismatches up to %d (done_level %d, fast_level %d) ***\n",
	       max_mismatches_allowed,*done_level,fast_level));

  if (1 || max_mismatches_allowed >= 0) {
    this->plus_segments = identify_all_segments(&this->plus_nsegments,&(*plus_anchor_segments),&(*n_plus_anchors),
						&this->plus_spliceable,&this->plus_nspliceable,
#ifdef LARGE_GENOMES
						this->plus_positions_high,this->plus_positions_low,
#else
						this->plus_positions,
#endif
						this->plus_npositions,this->omitted,querylength,query_lastpos,
						floors,/*plusp*/true);
    this->minus_segments = identify_all_segments(&this->minus_nsegments,&(*minus_anchor_segments),&(*n_minus_anchors),
						 &this->minus_spliceable,&this->minus_nspliceable,
#ifdef LARGE_GENOMES
						 this->minus_positions_high,this->minus_positions_low,
#else
						 this->minus_positions,
#endif
						 this->minus_npositions,this->omitted,querylength,query_lastpos,
						 floors,/*plusp*/false);

    *subs = find_complete_mm(&(*found_score),&(*nhits),*subs,*plus_anchor_segments,*n_plus_anchors,
			     querylength,/*queryptr:queryuc_ptr,*/
			     /*query_compress*/query_compress_fwd,
			     /*max_mismatches_allowed*/*done_level,/*plusp*/true,genestrand);

    *subs = find_complete_mm(&(*found_score),&(*nhits),*subs,*minus_anchor_segments,*n_minus_anchors,
			     querylength,/*queryptr:queryrc,*/
			     /*query_compress*/query_compress_rev,
			     /*max_mismatches_allowed*/*done_level,/*plusp*/false,genestrand);

    *segments_computed_p = true;

    debug(printf("5> found_score = %d, opt_level %d, done_level %d\n",*found_score,*opt_level,*done_level));
    debug(printf("plus_nsegments = %d, minus_nsegments = %d\n",this->plus_nsegments,this->minus_nsegments));
  }

  if (revise_levels_p == true) {
    *opt_level = (*found_score < *opt_level) ? *found_score : *opt_level;
    if ((*done_level = *opt_level + subopt_levels) > user_maxlevel) {
      *done_level = user_maxlevel;
    }
  }

#if 0
    opt_level = (found_score < opt_level) ? found_score : opt_level;
    if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
      done_level = user_maxlevel;
    }
    debug(printf("10> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
#endif

  if (*done_level >= indel_penalty_middle || *done_level >= indel_penalty_end) {
    /* 6. Indels */
    /* Need to reverse, because middle indelsplicing procedure depends on ascending diagonal order */

    if (*segments_computed_p == false) {
      this->plus_segments = identify_all_segments(&this->plus_nsegments,&(*plus_anchor_segments),&(*n_plus_anchors),
						  &this->plus_spliceable,&this->plus_nspliceable,
#ifdef LARGE_GENOMES
						  this->plus_positions_high,this->plus_positions_low,
#else
						  this->plus_positions,
#endif
						  this->plus_npositions,this->omitted,querylength,query_lastpos,
						  floors,/*plusp*/true);
      this->minus_segments = identify_all_segments(&this->minus_nsegments,&(*minus_anchor_segments),&(*n_minus_anchors),
						   &this->minus_spliceable,&this->minus_nspliceable,
#ifdef LARGE_GENOMES
						   this->minus_positions_high,this->minus_positions_low,
#else
						   this->minus_positions,
#endif
						   this->minus_npositions,this->omitted,querylength,query_lastpos,
						   floors,/*plusp*/false);
      *segments_computed_p = true;
    }

#if 0
    /* Done iteratively */
    indel_level = indel_penalty;
    while (indel_level <= *done_level) {
      debug(printf("*** Stage 6.  Middle indels with %d-%d mismatches allowed\n",indel_level,indel_penalty));
      *indels = find_middle_indels(&(*found_score),&(*nhits),*indels,
				   this->plus_spliceable,this->plus_nspliceable,
				   this->minus_spliceable,this->minus_nspliceable,
#ifdef DEBUG2
				   queryuc_ptr,queryrc,
#endif
				   floors,querylength,query_lastpos,query_compress_fwd,query_compress_rev,genome,
				   /*indel_mismatches_allowed*/indel_level - indel_penalty,genestrand);

      if (allow_end_indels_p == true) {
	debug(printf("*** Stage 6.  End indels with %d-%d mismatches allowed\n",indel_level,indel_penalty));
	*indels = find_end_indels(&(*found_score),&(*nhits),*indels,*plus_anchor_segments,*minus_anchor_segments,
				  *n_plus_anchors,*n_minus_anchors,
#ifdef DEBUG2E
				  queryuc_ptr,queryrc,
#endif
				  querylength,firstbound,lastbound,query_compress_fwd,query_compress_rev,genome,
				  max_end_insertions,max_end_deletions,min_indel_end_matches,indel_penalty,
				  /*indel_mismatches_allowed*/indel_level - indel_penalty,genestrand);
      }
      if (revise_levels_p == true) {
	*opt_level = (*found_score < *opt_level) ? *found_score : *opt_level;
	if ((*done_level = *opt_level + subopt_levels) > user_maxlevel) {
	  *done_level = user_maxlevel;
	}
      }
      indel_level++;
      debug(printf("6> found_score = %d, opt_level %d, done_level %d\n",*found_score,*opt_level,*done_level));
    }
#else
    /* Do all in one sweep */
    debug(printf("*** Stage 6 (middle).  Middle indels with %d-%d mismatches allowed, found_score = %d\n",
		 *done_level,indel_penalty_middle,*found_score));
    *indels = find_middle_indels(&(*found_score),&(*nhits),*indels,
				 this->plus_spliceable,this->plus_nspliceable,
				 this->minus_spliceable,this->minus_nspliceable,
#ifdef DEBUG2
				 queryuc_ptr,queryrc,
#endif
				 floors,querylength,query_lastpos,query_compress_fwd,query_compress_rev,
				 /*indel_mismatches_allowed*/(*done_level) - indel_penalty_middle,genestrand);
    if (revise_levels_p == true) {
      *opt_level = (*found_score < *opt_level) ? *found_score : *opt_level;
      if ((*done_level = *opt_level + subopt_levels) > user_maxlevel) {
	*done_level = user_maxlevel;
      }
    }
    debug(printf("6 (middle)> found_score = %d, opt_level %d, done_level %d\n",*found_score,*opt_level,*done_level));

    if (allow_end_indels_p == true) {
      debug(printf("*** Stage 6 (end).  End indels with %d-%d mismatches allowed, found_score = %d\n",
		   *done_level,indel_penalty_end,*found_score));
      *indels = find_end_indels(&(*found_score),&(*nhits),*indels,*plus_anchor_segments,*minus_anchor_segments,
				*n_plus_anchors,*n_minus_anchors,
#ifdef DEBUG2E
				queryuc_ptr,queryrc,
#endif
				querylength,firstbound,lastbound,query_compress_fwd,query_compress_rev,
				max_end_insertions,max_end_deletions,min_indel_end_matches,indel_penalty_end,
				/*indel_mismatches_allowed*/(*done_level) - indel_penalty_end,genestrand);
      if (revise_levels_p == true) {
	*opt_level = (*found_score < *opt_level) ? *found_score : *opt_level;
	if ((*done_level = *opt_level + subopt_levels) > user_maxlevel) {
	  *done_level = user_maxlevel;
	}
      }
      debug(printf("6 (end)> found_score = %d, opt_level %d, done_level %d\n",*found_score,*opt_level,*done_level));
    }
    /* Calling procedure will invoke Stage3_remove_duplicates */
#endif

  }

  debug(printf("Finished with complete_set_mm_indels\n"));

  return;
}



static List_T
complete_set_singlesplicing (int *found_score, List_T localsplicing_orig, Floors_T floors, T this,
			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     int querylength, int query_lastpos, int localsplicing_penalty,
			     int max_mismatches_allowed, int genestrand, bool first_read_p,
			     bool subs_or_indels_p) {
  List_T localsplicing, localsplicing_plus, localsplicing_minus;
  List_T ambiguous = NULL, lowprob = NULL, p;
  Stage3end_T hit;
  int worst_nmatches;

  debug(printf("Starting complete_set_singlesplicing with %d mismatches allowed\n",max_mismatches_allowed));

  if (floors == NULL) {
    /* Returning NULL results in a memory leak */
    return localsplicing_orig;
  }

  localsplicing_plus = find_singlesplices_plus(&(*found_score),/*localsplicing*/NULL,&ambiguous,&lowprob,
					       this->plus_spliceable,this->plus_nspliceable,
					       floors,querylength,query_lastpos,
					       /*query_compress*/query_compress_fwd,
					       /*splicing_penalty*/localsplicing_penalty,
					       max_mismatches_allowed,first_read_p,genestrand,
					       subs_or_indels_p);
  localsplicing_plus = Splice_group_by_segmentj(&(*found_score),localsplicing_plus,&ambiguous,querylength,
						first_read_p,/*sarrayp*/false);

  localsplicing_minus = find_singlesplices_minus(&(*found_score),/*localsplicing*/NULL,&ambiguous,&lowprob,
						 this->minus_spliceable,this->minus_nspliceable,
						 floors,querylength,query_lastpos,
						 /*query_compress*/query_compress_rev,
						 /*splicing_penalty*/localsplicing_penalty,
						 max_mismatches_allowed,first_read_p,genestrand,
						 subs_or_indels_p);
  localsplicing_minus = Splice_group_by_segmentj(&(*found_score),localsplicing_minus,&ambiguous,querylength,
						 first_read_p,/*sarrayp*/false);

  debug(printf("Finished with complete_set_singlesplicing\n"));

  if (localsplicing_plus == NULL && localsplicing_minus == NULL) {
    return List_append(localsplicing_orig,List_append(ambiguous,lowprob));

  } else {
    localsplicing = List_append(localsplicing_orig,List_append(localsplicing_plus,localsplicing_minus));

    worst_nmatches = querylength;
    for (p = localsplicing; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      if (Stage3end_nmatches_posttrim(hit) < worst_nmatches) {
	worst_nmatches = Stage3end_nmatches_posttrim(hit);
      }
    }

    for (p = lowprob; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      if (Stage3end_nmatches_posttrim(hit) < worst_nmatches) {
	/* Dominated by both nmatches and probability */
	Stage3end_free(&hit);
      } else {
	/* Has worse probability but more matches, so keep */
	localsplicing = List_push(localsplicing,(void *) hit);
      }
    }
    List_free(&lowprob);

    return List_append(localsplicing,ambiguous);
  }
}


#if 0
static List_T
complete_set_doublesplicing (int *found_score, List_T localsplicing_orig, Floors_T floors, T this,
			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     char *queryuc_ptr, char *queryrc, int querylength, int query_lastpos,
			     int localsplicing_penalty,
			     int min_shortend, int max_mismatches_allowed, bool pairedp,
			     int genestrand, bool first_read_p, bool subs_or_indels_p) {
  List_T localsplicing, localsplicing_plus, localsplicing_minus;
  List_T lowprob = NULL, p;
  Stage3end_T hit;
  int worst_nmatches;
  
  debug(printf("Starting complete_set_doublesplicing with %d mismatches allowed\n",
	       max_mismatches_allowed));

  if (floors == NULL) {
    /* Returning NULL results in a memory leak */
    return localsplicing_orig;
  }

  localsplicing_plus = find_doublesplices(&(*found_score),/*localsplicing*/NULL,&lowprob,
					  this->plus_spliceable,this->plus_nspliceable,this->plus_segments,
					  /*queryptr*/queryuc_ptr,querylength,query_lastpos,
					  /*query_compress*/query_compress_fwd,
					  /*max_distance*/shortsplicedist,/*splicing_penalty*/localsplicing_penalty,
					  min_shortend,max_mismatches_allowed,pairedp,first_read_p,
					  /*plusp*/true,genestrand,subs_or_indels_p);

  localsplicing_minus = find_doublesplices(&(*found_score),/*localsplicing*/NULL,&lowprob,
					   this->minus_spliceable,this->minus_nspliceable,this->minus_segments,
					   /*queryptr*/queryrc,querylength,query_lastpos,
					   /*query_compress*/query_compress_rev,
					   /*max_distance*/shortsplicedist,/*splicing_penalty*/localsplicing_penalty,
					   min_shortend,max_mismatches_allowed,pairedp,first_read_p,
					   /*plusp*/false,genestrand,subs_or_indels_p);

  debug(printf("Finished with complete_set_doublesplicing\n"));

  if (localsplicing_plus == NULL && localsplicing_minus == NULL) {
    return List_append(localsplicing_orig,lowprob);

  } else {
    localsplicing = List_append(localsplicing_orig,List_append(localsplicing_plus,localsplicing_minus));

    worst_nmatches = querylength;
    for (p = localsplicing; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      if (Stage3end_nmatches_posttrim(hit) < worst_nmatches) {
	worst_nmatches = Stage3end_nmatches_posttrim(hit);
      }
    }

    for (p = lowprob; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      if (Stage3end_nmatches_posttrim(hit) < worst_nmatches) {
	/* Dominated by both nmatches and probability */
	Stage3end_free(&hit);
      } else {
	/* Has worse probability but more matches, so keep */
	localsplicing = List_push(localsplicing,(void *) hit);
      }
    }
    List_free(&lowprob);

    return localsplicing;
  }
}
#endif


/* Simple table */
typedef struct History_T *History_T;
struct History_T {
  List_T keys;			/* List of Univinterval_T */
  List_T values;		/* List of List_T of Stage3end_T */
};

static void
History_free (History_T *old) {
  Univinterval_T interval;
  List_T hitlist, p;

  for (p = (*old)->keys; p != NULL; p = p->rest) {
    interval = (Univinterval_T) p->first;
    Univinterval_free(&interval);
  }
  List_free(&(*old)->keys);
  
  for (p = (*old)->values; p != NULL; p = p->rest) {
    hitlist = (List_T) p->first;
    Stage3end_list_free(&hitlist);
  }
  List_free(&(*old)->values);

  FREE(*old);
  return;
}

static History_T
History_new () {
  History_T new = (History_T) MALLOC(sizeof(*new));

  new->keys = (List_T) NULL;
  new->values = (List_T) NULL;
  return new;
}

static List_T
History_get (History_T this, Univinterval_T interval) {
  List_T p, q;

  for (p = this->keys, q = this->values; p != NULL; p = p->rest, q = q->rest) {
    if (Univinterval_equal(interval,(Univinterval_T) p->first) == true) {
      return (List_T) q->first;
    }
  }
  return (List_T) NULL;
}

static void
History_put (History_T this, Univinterval_T interval, List_T gmap_hits) {
  this->keys = List_push(this->keys,(void *) interval);
  this->values = List_push(this->values,(void *) gmap_hits);
  return;
}


/* Also defined in sarray-read.c and stage3hr.c */
#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


static List_T
run_gmap_for_region (bool *good_start_p, bool *good_end_p, History_T gmap_history,
		     List_T hits, char *accession, char *queryuc_ptr, int querylength,
		     int sense_try, bool favor_right_p, 
		     
		     Univcoord_T mappingstart, Univcoord_T mappingend,
		     Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		     bool watsonp, int genestrand,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		     
		     Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		     Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		     Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  Stage3end_T hit;
#ifdef EXTRACT_GENOMICSEG
  char *genomicseg, *genomicseg_alloc;
#endif

  List_T stored_hits;
  List_T all_stage2results, p;
  Stage2_T stage2;

  Univinterval_T interval;
  List_T pairs1, pairs2;
  struct Pair_T *pairarray1, *pairarray2;
  Univcoord_T start, end;
  double min_splice_prob_1, min_splice_prob_2;
  int goodness1, goodness2;
  int npairs1, npairs2, nsegments, nmismatches_whole, nindels, nintrons, nindelbreaks;
  int cdna_direction, sensedir;
  int matches1, unknowns1, mismatches1, qopens1, qindels1, topens1, tindels1,
    ncanonical1, nsemicanonical1, nnoncanonical1;
  int matches2, unknowns2, mismatches2, qopens2, qindels2, topens2, tindels2,
    ncanonical2, nsemicanonical2, nnoncanonical2;
  int nmatches_posttrim_1, max_match_length_1, ambig_end_length_5_1, ambig_end_length_3_1;
  int nmatches_posttrim_2, max_match_length_2, ambig_end_length_5_2, ambig_end_length_3_2;
  Splicetype_T ambig_splicetype_5_1, ambig_splicetype_3_1;
  Splicetype_T ambig_splicetype_5_2, ambig_splicetype_3_2;
  double ambig_prob_5_1, ambig_prob_3_1;
  double ambig_prob_5_2, ambig_prob_3_2;
  /* int maxintronlen_bound; */


  debug13(printf("Running GMAP at mappingstart %u + %d = mappingend %u, watsonp %d, sense_try %d, querylength %d, limits %u..%u\n",
		 (Chrpos_T) (mappingstart-chroffset),mappingend-mappingstart,
		 (Chrpos_T) (mappingend-chroffset),watsonp,sense_try,querylength,
		 (Chrpos_T) (knownsplice_limit_low-chroffset),(Chrpos_T) (knownsplice_limit_high-chroffset)));

  *good_start_p = *good_end_p = false;

  /* It is possible for mappingend to equal mappingstart if the read
     is forced to the beginning or end of a chromosome */
  if (mappingend > mappingstart) {
    interval = Univinterval_new(mappingstart,mappingend,sense_try);
    debug13(printf("Checking history for interval at %u..%u (sense_try %d)\n",
		   mappingstart,mappingend,sense_try));
    if ((stored_hits = History_get(gmap_history,interval)) != NULL) {
      debug13(printf("Already ran these coordinates, and have results\n"));
      for (p = stored_hits; p != NULL; p = List_next(p)) {
	if ((hit = (Stage3end_T) List_head(p)) != NULL) {
	  if (Stage3end_trim_left(hit) < GOOD_GMAP_END) {
	    *good_start_p = true;
	  }
	  if (Stage3end_trim_right(hit) < GOOD_GMAP_END) {
	    *good_end_p = true;
	  }
	  hits = List_push(hits,(void *) Stage3end_copy(hit));
	}
      }
      Univinterval_free(&interval);
      return hits;
    } else {
      debug13(printf("New coordinates\n"));
      /* stored_hits = (List_T) NULL; -- Already NULL */
    }

#ifdef EXTRACT_GENOMICSEG
    if (watsonp == true) {
      printf("Allocating %u bytes\n",genomiclength);
      genomicseg_alloc = (char *) CALLOC(genomiclength+MAX_INDEXSIZE+1,sizeof(char));
      genomicseg = &(genomicseg_alloc[MAX_INDEXSIZE]);
      Genome_fill_buffer_blocks(genomicstart-MAX_INDEXSIZE,genomiclength+MAX_INDEXSIZE,genomicseg_alloc);
    } else {
      printf("Allocating %u bytes\n",genomiclength);
      genomicseg_alloc = (char *) CALLOC(genomiclength+MAX_INDEXSIZE+1,sizeof(char));
      genomicseg = &(genomicseg_alloc[MAX_INDEXSIZE]);
      Genome_fill_buffer_blocks(genomicstart,genomiclength+MAX_INDEXSIZE,genomicseg_alloc);
      make_complement_inplace(genomicseg_alloc,genomiclength+MAX_INDEXSIZE);
    }
#endif

#if 0
    /* Should be able to have splicing on a circular chromosome */
    if (chroffset + chrlength < chrhigh) {
      debug13(printf("Chromosome is circular because chroffset %u + chrlength %u < chrhigh %u\n",
		     chroffset,chrlength,chrhigh));
      maxintronlen_bound = 0;
    } else {
      maxintronlen_bound = shortsplicedist;
    }
#else
    /* maxintronlen_bound = shortsplicedist; */
    /* maxintronlen_bound = overall_max_distance;  */ /* was shortsplicedist, but misses microinversions */
#endif



    /* Note: Use nmatches post-trim to decide if the alignment is high
       quality or worth keeping.  But if so, then use nmatches_pretrim
       for ranking and scoring purposes. */

    /* use_shifted_canonical_p == true can be slow and can give wrong answers */
    all_stage2results = Stage2_compute(/*queryseq_ptr*/queryuc_ptr,queryuc_ptr,querylength,/*query_offset*/0,
				       /*chrstart*/mappingstart-chroffset,/*chrend*/mappingend-chroffset,
				       chroffset,chrhigh,/*plusp*/watsonp,genestrand,
				       oligoindices_major,pairpool,diagpool,cellpool,/*localp*/true,
				       /*skip_repetitive_p*/true,favor_right_p,/*max_nalignments*/MAX_NALIGNMENTS,
				       /*debug_graphic_p*/false,/*worker_stopwatch*/NULL,/*diag_debug*/false);

    debug13(printf("Got %d stage2 results\n",List_length(all_stage2results)));

    if (all_stage2results == NULL) {
      stored_hits = List_push(stored_hits,(void *) NULL);
    }

    for (p = all_stage2results; p != NULL; p = List_next(p)) {
      stage2 = (Stage2_T) List_head(p);
      if ((pairarray1 = Stage3_compute(&cdna_direction,&sensedir,&pairs1,&npairs1,&goodness1,
				       &matches1,&nmatches_posttrim_1,&max_match_length_1,
				       &ambig_end_length_5_1,&ambig_end_length_3_1,
				       &ambig_splicetype_5_1,&ambig_splicetype_3_1,
				       &ambig_prob_5_1,&ambig_prob_3_1,
				       &unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
				       &ncanonical1,&nsemicanonical1,&nnoncanonical1,&min_splice_prob_1,

				       &pairarray2,&pairs2,&npairs2,&goodness2,
				       &matches2,&nmatches_posttrim_2,&max_match_length_2,
				       &ambig_end_length_5_2,&ambig_end_length_3_2,
				       &ambig_splicetype_5_2,&ambig_splicetype_3_2,
				       &ambig_prob_5_2,&ambig_prob_3_2,
				       &unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
				       &ncanonical2,&nsemicanonical2,&nnoncanonical2,&min_splice_prob_2,

				       Stage2_middle(stage2),Stage2_all_starts(stage2),Stage2_all_ends(stage2),
#ifdef END_KNOWNSPLICING_SHORTCUT
				       cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
#endif
				       /*queryseq_ptr*/queryuc_ptr,queryuc_ptr,querylength,/*skiplength*/0,
#ifdef EXTRACT_GENOMICSEG
				       /*query_subseq_offset*/0,
#else
				       /*query_subseq_offset*/0,
#endif
				       chrnum,chroffset,chrhigh,
				       knownsplice_limit_low,knownsplice_limit_high,watsonp,genestrand,
				       /*jump_late_p*/watsonp ? false : true,
				       
				       maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
				       sense_try,/*sense_filter*/0,
				       oligoindices_minor,diagpool,cellpool)) == NULL) {
	debug13(printf("stage3 is NULL\n"));
	stored_hits = List_push(stored_hits,(void *) NULL);

      } else if (cdna_direction == 0) {
	debug13(printf("stage3 is not NULL, and cdna direction not determined\n"));
	debug13a(Pair_dump_array(pairarray1,npairs1,true));

	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray1,npairs1);
	if (watsonp == true) {
	  start = subtract_bounded(chroffset + Pair_genomepos(&(pairarray1[0])),
				   /*minusterm*/Pair_querypos(&(pairarray1[0])),chroffset);
	  end = add_bounded(chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			    /*plusterm*/querylength - 1 - Pair_querypos(&(pairarray1[npairs1-1])),chrhigh);
	  if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_1,max_match_length_1,
					ambig_end_length_5_1,ambig_end_length_3_1,
					ambig_splicetype_5_1,ambig_splicetype_3_1,
					min_splice_prob_1,
					pairarray1,npairs1,nsegments,nintrons,nindelbreaks,
					/*left*/start,/*genomiclength*/end - start + 1,
					/*plusp*/watsonp,genestrand,
					accession,querylength,chrnum,chroffset,chrhigh,chrlength,
					/*cdna_direction*/+1,/*sensedir*/SENSE_FORWARD,
					/*gmap_source*/GMAP_VIA_REGION)) == NULL) {
	    debug13(printf("Stage3end_new_gmap returns NULL\n"));
	    stored_hits = List_push(stored_hits,(void *) NULL);
	    FREE_OUT(pairarray1);

	  } else {
	    if (Stage3end_trim_left(hit) < GOOD_GMAP_END) {
	      *good_start_p = true;
	    }
	    if (Stage3end_trim_right(hit) < GOOD_GMAP_END) {
	      *good_end_p = true;
	    }
	    debug13(printf("Trim at start: %d, trim at end: %d\n",
			   Stage3end_trim_left(hit),Stage3end_trim_right(hit)));
	    stored_hits = List_push(stored_hits,(void *) Stage3end_copy(hit));
	    hits = List_push(hits,(void *) hit);
	  }

	} else {
	  start = add_bounded(chroffset + Pair_genomepos(&(pairarray1[0])),
			      /*plusterm*/Pair_querypos(&(pairarray1[0])),chrhigh);
	  end = subtract_bounded(chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
				 /*minusterm*/querylength - 1 - Pair_querypos(&(pairarray1[npairs1-1])),chroffset);
	  if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_1,max_match_length_1,
					ambig_end_length_5_1,ambig_end_length_3_1,
					ambig_splicetype_5_1,ambig_splicetype_3_1,
					min_splice_prob_1,
					pairarray1,npairs1,nsegments,nintrons,nindelbreaks,
					/*left*/end,/*genomiclength*/start - end + 1,
					/*plusp*/watsonp,genestrand,
					accession,querylength,chrnum,chroffset,chrhigh,chrlength,
					/*cdna_direction*/+1,/*sensedir*/SENSE_FORWARD,
					/*gmap_source*/GMAP_VIA_REGION)) == NULL) {
	    debug13(printf("Stage3end_new_gmap returns NULL\n"));
	    stored_hits = List_push(stored_hits,(void *) NULL);
	    FREE_OUT(pairarray1);

	  } else {
	    if (Stage3end_trim_right(hit) < GOOD_GMAP_END) {
	      *good_start_p = true;
	    }
	    if (Stage3end_trim_left(hit) < GOOD_GMAP_END) {
	      *good_end_p = true;
	    }
	    debug13(printf("Trim at start: %d, trim at end: %d (raw %d and %d)\n",
			   Stage3end_trim_right(hit),Stage3end_trim_left(hit),
			   Stage3end_trim_right_raw(hit),Stage3end_trim_left_raw(hit)));
	    /* Don't throw away GMAP hits */
	    stored_hits = List_push(stored_hits,(void *) Stage3end_copy(hit));
	    hits = List_push(hits,(void *) hit);
	  }
	}
	/* Don't free pairarray1 */

	debug13a(Pair_dump_array(pairarray2,npairs2,true));

	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray2,npairs2);
	if (watsonp == true) {
	  start = subtract_bounded(chroffset + Pair_genomepos(&(pairarray2[0])),
				   /*minusterm*/Pair_querypos(&(pairarray2[0])),chroffset);
	  end = add_bounded(chroffset + Pair_genomepos(&(pairarray2[npairs2-1])),
			    /*plusterm*/querylength - 1 - Pair_querypos(&(pairarray2[npairs2-1])),chrhigh);
	  if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_2,max_match_length_2,
					ambig_end_length_5_2,ambig_end_length_3_2,
					ambig_splicetype_5_2,ambig_splicetype_3_2,
					min_splice_prob_2,
					pairarray2,npairs2,nsegments,nintrons,nindelbreaks,
					/*left*/start,/*genomiclength*/end - start + 1,
					/*plusp*/watsonp,genestrand,
					accession,querylength,chrnum,chroffset,chrhigh,chrlength,
					/*cdna_direction*/-1,/*sensedir*/SENSE_ANTI,
					/*gmap_source*/GMAP_VIA_REGION)) == NULL) {
	    debug13(printf("Stage3end_new_gmap returns NULL\n"));
	    stored_hits = List_push(stored_hits,(void *) NULL);
	    FREE_OUT(pairarray2);

	  } else {
	    if (Stage3end_trim_left(hit) < GOOD_GMAP_END) {
	      *good_start_p = true;
	    }
	    if (Stage3end_trim_right(hit) < GOOD_GMAP_END) {
	      *good_end_p = true;
	    }
	    debug13(printf("Trim at start: %d, trim at end: %d\n",
			   Stage3end_trim_left(hit),Stage3end_trim_right(hit)));
	    stored_hits = List_push(stored_hits,(void *) Stage3end_copy(hit));
	    hits = List_push(hits,(void *) hit);
	  }

	} else {
	  start = add_bounded(chroffset + Pair_genomepos(&(pairarray2[0])),
			      /*plusterm*/Pair_querypos(&(pairarray2[0])),chrhigh);
	  end = subtract_bounded(chroffset + Pair_genomepos(&(pairarray2[npairs2-1])),
				 /*minusterm*/querylength - 1 - Pair_querypos(&(pairarray2[npairs2-1])),chroffset);
	  if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_2,max_match_length_2,
					ambig_end_length_5_2,ambig_end_length_3_2,
					ambig_splicetype_5_2,ambig_splicetype_3_2,
					min_splice_prob_2,
					pairarray2,npairs2,nsegments,nintrons,nindelbreaks,
					/*left*/end,/*genomiclength*/start - end + 1,
					/*plusp*/watsonp,genestrand,
					accession,querylength,chrnum,chroffset,chrhigh,chrlength,
					/*cdna_direction*/-1,/*sensedir*/SENSE_ANTI,
					/*gmap_source*/GMAP_VIA_REGION)) == NULL) {
	    debug13(printf("Stage3end_new_gmap returns NULL\n"));
	    stored_hits = List_push(stored_hits,(void *) NULL);
	    FREE_OUT(pairarray2);

	  } else {
	    if (Stage3end_trim_right(hit) < GOOD_GMAP_END) {
	      *good_start_p = true;
	    }
	    if (Stage3end_trim_left(hit) < GOOD_GMAP_END) {
	      *good_end_p = true;
	    }
	    debug13(printf("Trim at start: %d, trim at end: %d (raw %d and %d)\n",
			   Stage3end_trim_right(hit),Stage3end_trim_left(hit),
			   Stage3end_trim_right_raw(hit),Stage3end_trim_left_raw(hit)));
	    /* Don't throw away GMAP hits */
	    stored_hits = List_push(stored_hits,(void *) Stage3end_copy(hit));
	    hits = List_push(hits,(void *) hit);
	  }
	}
	/* Don't free pairarray2 */

      } else {
	debug13(printf("stage3 is not NULL, and cdna direction is determined to be %d\n",cdna_direction));
	debug13a(Pair_dump_array(pairarray1,npairs1,true));

	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray1,npairs1);
	if (watsonp == true) {
	  start = subtract_bounded(chroffset + Pair_genomepos(&(pairarray1[0])),
				   /*minusterm*/Pair_querypos(&(pairarray1[0])),chroffset);
	  end = add_bounded(chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			    /*plusterm*/querylength - 1 - Pair_querypos(&(pairarray1[npairs1-1])),chrhigh);
	  if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_1,max_match_length_1,
					ambig_end_length_5_1,ambig_end_length_3_1,
					ambig_splicetype_5_1,ambig_splicetype_3_1,
					min_splice_prob_1,
					pairarray1,npairs1,nsegments,nintrons,nindelbreaks,
					/*left*/start,/*genomiclength*/end - start + 1,
					/*plusp*/watsonp,genestrand,
					accession,querylength,chrnum,chroffset,chrhigh,chrlength,
					cdna_direction,sensedir,/*gmap_source*/GMAP_VIA_REGION)) == NULL) {
	    debug13(printf("Stage3end_new_gmap returns NULL\n"));
	    stored_hits = List_push(stored_hits,(void *) NULL);
	    FREE_OUT(pairarray1);

	  } else {
	    if (Stage3end_trim_left(hit) < GOOD_GMAP_END) {
	      *good_start_p = true;
	    }
	    if (Stage3end_trim_right(hit) < GOOD_GMAP_END) {
	      *good_end_p = true;
	    }
	    debug13(printf("Trim at start: %d, trim at end: %d\n",
			   Stage3end_trim_left(hit),Stage3end_trim_right(hit)));
	    stored_hits = List_push(stored_hits,(void *) Stage3end_copy(hit));
	    hits = List_push(hits,(void *) hit);
	  }

	} else {
	  start = add_bounded(chroffset + Pair_genomepos(&(pairarray1[0])),
			      /*plusterm*/Pair_querypos(&(pairarray1[0])),chrhigh);
	  end = subtract_bounded(chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
				 /*minusterm*/querylength - 1 - Pair_querypos(&(pairarray1[npairs1-1])),chroffset);
	  if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_1,max_match_length_1,
					ambig_end_length_5_1,ambig_end_length_3_1,
					ambig_splicetype_5_1,ambig_splicetype_3_1,
					min_splice_prob_1,
					pairarray1,npairs1,nsegments,nintrons,nindelbreaks,
					/*left*/end,/*genomiclength*/start - end + 1,
					/*plusp*/watsonp,genestrand,
					accession,querylength,chrnum,chroffset,chrhigh,chrlength,
					cdna_direction,sensedir,/*gmap_source*/GMAP_VIA_REGION)) == NULL) {
	    debug13(printf("Stage3end_new_gmap returns NULL\n"));
	    stored_hits = List_push(stored_hits,(void *) NULL);
	    FREE_OUT(pairarray1);

	  } else {
	    if (Stage3end_trim_right(hit) < GOOD_GMAP_END) {
	      *good_start_p = true;
	    }
	    if (Stage3end_trim_left(hit) < GOOD_GMAP_END) {
	      *good_end_p = true;
	    }
	    debug13(printf("Trim at start: %d, trim at end: %d (raw %d and %d)\n",
			   Stage3end_trim_right(hit),Stage3end_trim_left(hit),
			   Stage3end_trim_right_raw(hit),Stage3end_trim_left_raw(hit)));
	    /* Don't throw away GMAP hits */
	    stored_hits = List_push(stored_hits,(void *) Stage3end_copy(hit));
	    hits = List_push(hits,(void *) hit);
	  }
	}
	/* Don't free pairarray1 */

      }

      Stage2_free(&stage2);
    }
    List_free(&all_stage2results);


#ifdef EXTRACT_GENOMICSEG
    FREE(genomicseg_alloc);
#endif

    debug13(printf(" => Got good_start_p %d, good_end_p %d\n",*good_start_p,*good_end_p));
    debug13(printf("Storing history for interval at %u..%u (sense_try %d)\n",
		   mappingstart,mappingend,sense_try));
    History_put(gmap_history,interval,stored_hits);
  }


  return hits;
}


static void
align_single_hit_with_gmap (Stage3end_T *result1, Stage3end_T *result2, Stage3end_T hit,
			    char *queryuc_ptr, int querylength,
#ifdef END_KNOWNSPLICING_SHORTCUT
			    char *queryrc, bool invertedp,
#endif
			    Oligoindex_array_T oligoindices_minor,
			    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			    int genestrand, bool extend_ends_p) {
  /* Both events are tested by Stage3end_anomalous_splice_p */
  if (Stage3end_chrnum(hit) == 0) {
    /* Translocation */
    *result1 = *result2 = (Stage3end_T) NULL;

  } else if (Stage3end_hittype(hit) == SAMECHR_SPLICE) {
    /* A genomic event that doesn't get reflected in chrnum */
    *result1 = *result2 = (Stage3end_T) NULL;

  } else if (Stage3end_hittype(hit) == GMAP) {
    if (Stage3end_plusp(hit) == true) {
      *result1 = Stage3end_gmap_run_gmap_plus(&(*result2),hit,queryuc_ptr,querylength,genestrand,
					      maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
					      oligoindices_minor,diagpool,cellpool);
    } else {
      *result1 = Stage3end_gmap_run_gmap_minus(&(*result2),hit,queryuc_ptr,querylength,genestrand,
					       maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
					       oligoindices_minor,diagpool,cellpool);
    }

  } else if (Stage3end_plusp(hit) == true) {
    debug13(printf("Running Stage3end_substrings_run_gmap_plus\n"));
    *result1 = Stage3end_substrings_run_gmap_plus(&(*result2),hit,queryuc_ptr,querylength,genestrand,
						  maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
						  oligoindices_minor,diagpool,cellpool,extend_ends_p);
  } else {
    debug13(printf("Running Stage3end_substrings_run_gmap_minus\n"));
    *result1 = Stage3end_substrings_run_gmap_minus(&(*result2),hit,queryuc_ptr,querylength,genestrand,
						   maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
						   oligoindices_minor,diagpool,cellpool,extend_ends_p);
  }
  return;
}


#if 0
static List_T
convert_plus_segments_to_gmap_via_region (History_T gmap_history, List_T hits,
					  char *accession, char *queryuc_ptr, int querylength, int query_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
					  char *queryrc, bool invertedp,
#endif
					  Compress_T query_compress_fwd, Compress_T query_compress_rev,
					  List_T anchor_segments, struct Segment_T *plus_segments, int plus_nsegments,
					  Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
					  Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
					  Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
					  int user_maxlevel, int genestrand, bool first_read_p,
					  bool require_pairing_p) {
  Univcoord_T segmentstart, segmentend, left;
  Univcoord_T mappingstart, mappingend, chroffset, chrhigh, mappingpos;
  Univcoord_T origlow, orighigh;
  Univcoord_T close_mappingstart_last, close_mappingend_last,
    middle_mappingstart_last, middle_mappingend_last;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T close_knownsplice_limit_low, close_knownsplice_limit_high;
  Chrpos_T chrlength;
  Chrnum_T chrnum;
  bool close_mappingstart_p = false, close_mappingend_p = false;
  bool middle_mappingstart_p = false, middle_mappingend_p = false;
  bool fallback_mappingstart_p, fallback_mappingend_p;
  bool good_start_p, good_end_p, favor_right_p = false;
  bool novelp;		 /* Want any of the segments in (startk+1)..(endk-1) to not be used */
  bool pairablep;		/* Want any of the segments in (startk+1)..(endk-1) to be pairable */

  List_T p;
  Segment_T anchor_segment;
  int anchork, startk, endk, k;

  anchork = 0;
  for (p = anchor_segments; p != NULL; p = List_next(p)) {
    anchor_segment = (Segment_T) List_head(p);
    assert(anchor_segment->diagonal != (Univcoord_T) -1);
    while (plus_segments[anchork].diagonal != anchor_segment->diagonal) {
      anchork++;
    }

    novelp = (anchor_segment->usedp == true) ? false : true;
    pairablep = anchor_segment->pairablep;
    anchor_segment->usedp = true;

    startk = anchork - 1;
    while (startk >= 0 && plus_segments[startk].diagonal != (Univcoord_T) -1 &&
	   plus_segments[startk].diagonal + shortsplicedist > anchor_segment->diagonal) {
      if (plus_segments[startk].usedp == false) {
	novelp = true;
      }
      plus_segments[startk].usedp = true;
      if (plus_segments[startk].pairablep == true) {
	pairablep = true;
      }
      startk--;
    }

    endk = anchork + 1;
    while (endk < plus_nsegments && plus_segments[endk].diagonal < anchor_segment->diagonal + shortsplicedist) {
      if (plus_segments[endk].usedp == false) {
	novelp = true;
      }
      plus_segments[endk].usedp = true;
      if (plus_segments[endk].pairablep == true) {
	pairablep = true;
      }
      endk++;
    }

    if (novelp == true && (pairablep == true || require_pairing_p == false)) {
      debug13(printf("Processing segments %d to %d inclusive\n",startk+1,endk-1));
      chrnum = anchor_segment->chrnum;
      chroffset = anchor_segment->chroffset;
      chrhigh = anchor_segment->chrhigh;
      chrlength = anchor_segment->chrlength;

      left = anchor_segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      origlow = left - anchor_segment->querypos5;
      orighigh = left + (querylength - anchor_segment->querypos3);

      /* extend left */
      knownsplice_limit_low = subtract_bounded(origlow,overall_max_distance,chroffset);
      mappingstart = segmentstart = subtract_bounded(origlow,overall_max_distance,chroffset);
      debug13(printf("Original bounds A: knownsplice_limit_low %u, mappingstart %u\n",
		     knownsplice_limit_low - chroffset,mappingstart - chroffset));

      /* extend right */
      knownsplice_limit_high = add_bounded(orighigh,shortsplicedist,chrhigh);
      mappingend = segmentend = add_bounded(orighigh,shortsplicedist,chrhigh);
      debug13(printf("Original bounds B: knownsplice_limit_high %u, mappingend %u\n",
		     knownsplice_limit_high - chroffset,mappingend - chroffset));
      
      close_mappingstart_last = middle_mappingstart_last = origlow;
      close_mappingend_last = middle_mappingend_last = orighigh;
      close_mappingstart_p = close_mappingend_p = false;
      middle_mappingstart_p = middle_mappingend_p = false;

      /* 1 */
      for (k = startk + 1; k < endk; k++) {
	debug13(printf("1. plus diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (plus_segments[k].diagonal - chroffset),(unsigned long long) plus_segments[k].diagonal,
		       plus_segments[k].querypos5,plus_segments[k].querypos3,plus_segments[k].usedp,plus_segments[k].pairablep));
	if (plus_segments[k].querypos5 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 3. Missing start of query, so there could be a middle splice */
	  debug13b(printf("  querypos5 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  plus_segments[k].querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(plus_segments[k].diagonal,querylength + shortsplicedist,chroffset)) < middle_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    middle_mappingstart_last = mappingpos;
	    middle_mappingstart_p = true;
	    debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	  }

	} else {
	  debug13b(printf("  querypos5 %d < %d + %d, so using this diagonal\n",
			  plus_segments[k].querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(plus_segments[k].diagonal,querylength,chroffset)) < close_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    close_mappingstart_last = mappingpos;
	    close_mappingstart_p = true;
	    debug13(printf("  Redefining close mappingstart last to %u\n",close_mappingstart_last - chroffset));
	  }
	}


	if (query_lastpos - plus_segments[k].querypos3 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 1. Missing end of query, so there could be a middle splice */
	  debug13b(printf("  query_lastpos %d - querypos3 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  query_lastpos,plus_segments[k].querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = add_bounded(plus_segments[k].diagonal,shortsplicedist,chrhigh)) > middle_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    middle_mappingend_last = mappingpos;
	    middle_mappingend_p = true;
	    debug13(printf("  Redefining middle mappingend last to %u\n",middle_mappingend_last - chroffset));
	  }

	} else {
	  debug13b(printf("  query_lastpos %d - querypos3 %d < %d + %d, so using this diagonal\n",
			  query_lastpos,plus_segments[k].querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = plus_segments[k].diagonal) > close_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    close_mappingend_last = mappingpos;
	    close_mappingend_p = true;
	    debug13(printf("  Redefining close mappingend last to %u\n",close_mappingend_last - chroffset));
	  }
	}
      }

      /* 2 */
      if (close_mappingstart_p == true) {
	close_knownsplice_limit_low = subtract_bounded(close_mappingstart_last,shortsplicedist,chroffset);
      } else if (middle_mappingstart_p == true) {
	debug13(printf("Using middle mappingstart\n"));
	close_knownsplice_limit_low = middle_mappingstart_last;
	close_mappingstart_last = middle_mappingstart_last;
	close_mappingstart_p = true;
      }

      if (middle_mappingstart_p == true && middle_mappingstart_last < close_mappingstart_last) {
	knownsplice_limit_low = middle_mappingstart_last;
	mappingstart = middle_mappingstart_last;
      }
      if (close_mappingstart_p == false) {
	fallback_mappingstart_p = false;
      } else {
	debug13(printf("Fallback mappingstart = %u\n",mappingstart - chroffset));
	fallback_mappingstart_p = true;
      }

      /* 3 */
      if (close_mappingend_p == true) {
	close_knownsplice_limit_high = add_bounded(close_mappingend_last,shortsplicedist,chrhigh);
      } else if (middle_mappingend_p == true) {
	close_knownsplice_limit_high = middle_mappingend_last;
	close_mappingend_last = middle_mappingend_last;
	close_mappingend_p = true;
	debug13(printf("Using middle mappingend => close_mappingend %u\n",close_mappingend_last));
      }
      if (middle_mappingend_p == true && middle_mappingend_last > close_mappingend_last) {
	knownsplice_limit_high = middle_mappingend_last;
	mappingend = middle_mappingend_last;
      }
      if (close_mappingend_p == false) {
	fallback_mappingend_p = false;
      } else {
	debug13(printf("Fallback mappingend = %u\n",mappingend - chroffset));
	fallback_mappingend_p = true;
      }

      /* 4 */
      if (close_mappingstart_p == true && close_mappingend_p == true) {
	debug13(printf("Single hit: Running gmap with close mappingstart and close mappingend\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				   /*sense_try*/0,favor_right_p,close_mappingstart_last,close_mappingend_last,
				   close_knownsplice_limit_low,close_knownsplice_limit_high,
				   /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
      
	if (good_start_p == true && good_end_p == true) {
	  /* Success */
	} else if (gmap_rerun_p == false) {
	  debug13(printf("Skipping re-run of gmap\n"));
	} else if (good_start_p == true) {
	  if (fallback_mappingend_p == true) {
	    debug13(printf("Single hit: Re-running gmap with close mappingstart only\n"));
	    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				       /*sense_try*/0,favor_right_p,close_mappingstart_last,mappingend,
				       close_knownsplice_limit_low,knownsplice_limit_high,
				       /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				       oligoindices_major,oligoindices_minor,
				       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	  }
	} else if (good_end_p == true) {
	  if (fallback_mappingstart_p == true) {
	    debug13(printf("Single hit: Re-running gmap with close mappingend only\n"));
	    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				       /*sense_try*/0,favor_right_p,mappingstart,close_mappingend_last,
				       knownsplice_limit_low,close_knownsplice_limit_high,
				       /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				       oligoindices_major,oligoindices_minor,
				       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	  }
	} else {
	  if (fallback_mappingstart_p == true && fallback_mappingend_p == true) {
	    debug13(printf("Single hit: Re-running gmap with far mappingstart and mappingend\n"));
	    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				       /*sense_try*/0,favor_right_p,mappingstart,mappingend,
				       knownsplice_limit_low,close_knownsplice_limit_high,
				       /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				       oligoindices_major,oligoindices_minor,
				       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	  }
	}

      } else if (close_mappingstart_p == true) {
	debug13(printf("Single hit: Running gmap with close mappingstart\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				   /*sense_try*/0,favor_right_p,close_mappingstart_last,mappingend,
				   close_knownsplice_limit_low,knownsplice_limit_high,
				   /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	if (good_start_p == true) {
	  /* Success */
	} else if (gmap_rerun_p == false) {
	  debug13(printf("Skipping re-run of gmap\n"));
	} else if (fallback_mappingstart_p == true) {
	  debug13(printf("Single hit: Re-running gmap with far mappingstart\n"));
	  hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				     /*sense_try*/0,favor_right_p,mappingstart,mappingend,
				     knownsplice_limit_low,knownsplice_limit_high,
				     /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				     oligoindices_major,oligoindices_minor,
				     pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	}

      } else if (close_mappingend_p == true) {
	debug13(printf("Single hit: Running gmap with close mappingend\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				   /*sense_try*/0,favor_right_p,mappingstart,close_mappingend_last,
				   knownsplice_limit_low,close_knownsplice_limit_high,
				   /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	if (good_end_p == true) {
	  /* Success */
	} else if (gmap_rerun_p == false) {
	  debug13(printf("Skipping re-run of gmap\n"));
	} else if (fallback_mappingend_p == true) {
	  debug13(printf("Single hit: Re-running gmap with far mappingend\n"));
	  hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				     /*sense_try*/0,favor_right_p,mappingstart,mappingend,
				     knownsplice_limit_low,knownsplice_limit_high,
				     /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				     oligoindices_major,oligoindices_minor,
				     pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	}

      } else {
	debug13(printf("Single hit: Running gmap with far mappingstart and mappingend\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				   /*sense_try*/0,favor_right_p,mappingstart,mappingend,
				   knownsplice_limit_low,knownsplice_limit_high,
				   /*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
      }
    }
  }

  return hits;
}
#endif


#if 0
static List_T
convert_minus_segments_to_gmap_via_region (History_T gmap_history, List_T hits,
					   char *accession, char *queryuc_ptr, int querylength, int query_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
					   char *queryrc, bool invertedp,
#endif
					   Compress_T query_compress_fwd, Compress_T query_compress_rev,
					   List_T anchor_segments, struct Segment_T *minus_segments, int minus_nsegments,
					   Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
					   Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
					   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
					   int user_maxlevel, int genestrand, bool first_read_p,
					   bool require_pairing_p) {
  Univcoord_T segmentstart, segmentend, left;
  Univcoord_T mappingstart, mappingend, chroffset, chrhigh, mappingpos;
  Univcoord_T origlow, orighigh;
  Univcoord_T close_mappingstart_last, close_mappingend_last,
    middle_mappingstart_last, middle_mappingend_last;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T close_knownsplice_limit_low, close_knownsplice_limit_high;
  Chrpos_T chrlength;
  Chrnum_T chrnum;
  bool close_mappingstart_p, close_mappingend_p;
  bool middle_mappingstart_p, middle_mappingend_p;
  bool fallback_mappingstart_p, fallback_mappingend_p;
  bool good_start_p, good_end_p, favor_right_p = true;
  bool novelp;		 /* Want any of the segments in starti..(endi-1) to not be used */
  bool pairablep;		/* Want any of the segments in starti..(endi-1) to be pairable */

  List_T p;
  Segment_T anchor_segment;
  int anchork, startk, endk, k;


  anchork = 0;
  for (p = anchor_segments; p != NULL; p = List_next(p)) {
    anchor_segment = (Segment_T) List_head(p);
    assert(anchor_segment->diagonal != (Univcoord_T) -1);
    while (minus_segments[anchork].diagonal != anchor_segment->diagonal) {
      anchork++;
    }

    novelp = (anchor_segment->usedp == true) ? false : true;
    pairablep = anchor_segment->pairablep;
    anchor_segment->usedp = true;

    startk = anchork - 1;
    while (startk >= 0 && minus_segments[startk].diagonal != (Univcoord_T) -1 &&
	   minus_segments[startk].diagonal + shortsplicedist > anchor_segment->diagonal) {
      if (minus_segments[startk].usedp == false) {
	novelp = true;
      }
      minus_segments[startk].usedp = true;
      if (minus_segments[startk].pairablep == true) {
	pairablep = true;
      }
      startk--;
    }

    endk = anchork + 1;
    while (endk < minus_nsegments && minus_segments[endk].diagonal < anchor_segment->diagonal + shortsplicedist) {
      if (minus_segments[endk].usedp == false) {
	novelp = true;
      }
      minus_segments[endk].usedp = true;
      if (minus_segments[endk].pairablep == true) {
	pairablep = true;
      }
      endk++;
    }
      
    if (novelp == true && (pairablep == true || require_pairing_p == false)) {
      debug13(printf("Processing segments %d to %d inclusive\n",startk+1,endk-1));
      chrnum = anchor_segment->chrnum;
      chroffset = anchor_segment->chroffset;
      chrhigh = anchor_segment->chrhigh;
      chrlength = anchor_segment->chrlength;

      left = anchor_segment->diagonal - querylength; /* FORMULA */
      origlow = left - (querylength - anchor_segment->querypos3);
      orighigh = left + anchor_segment->querypos5;
      
      /* extend right */
      knownsplice_limit_low = subtract_bounded(origlow,shortsplicedist,chroffset);
      mappingstart = segmentstart = subtract_bounded(origlow,shortsplicedist,chroffset);
      debug13(printf("Original bounds C: knownsplice_limit_low %u, mappingstart %u\n",
		     knownsplice_limit_low - chroffset,mappingstart - chroffset));
      
      /* extend left */
      knownsplice_limit_high = add_bounded(orighigh,overall_max_distance,chrhigh);
      mappingend = segmentend =	add_bounded(orighigh,overall_max_distance,chrhigh);
      debug13(printf("Original bounds D: knownsplice_limit_high %u, mappingend %u\n",
		     knownsplice_limit_high - chroffset,mappingend - chroffset));
      
      close_mappingstart_last = middle_mappingstart_last = origlow;
      close_mappingend_last = middle_mappingend_last = orighigh;
      close_mappingstart_p = close_mappingend_p = false;
      middle_mappingstart_p = middle_mappingend_p = false;
      
      /* 1 */
      for (k = startk + 1; k < endk; k++) {
	debug13(printf("1. minus diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (minus_segments[k].diagonal - chroffset),(unsigned long long) minus_segments[k].diagonal,
		       minus_segments[k].querypos5,minus_segments[k].querypos3,minus_segments[k].usedp,minus_segments[k].pairablep));
	if (query_lastpos - minus_segments[k].querypos3 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 2. Missing end of query, so there could be a middle splice */
	  debug13b(printf("  query_lastpos %d - querypos3 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  query_lastpos,minus_segments[k].querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(minus_segments[k].diagonal,querylength + shortsplicedist,chroffset)) < middle_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    middle_mappingstart_last = mappingpos;
	    middle_mappingstart_p = true;
	    debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	  }
	  
	} else {
	  debug13b(printf("  query_lastpos %d - querypos3 %d < %d + %d, so using this diagonal\n",
			  query_lastpos,minus_segments[k].querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(minus_segments[k].diagonal,querylength,chroffset)) < close_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    close_mappingstart_last = mappingpos;
	    close_mappingstart_p = true;
	    debug13(printf("  Redefining close mappingstart last to %u\n",close_mappingstart_last - chroffset));
	  }
	}
      

	if (minus_segments[k].querypos5 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 4. Missing start of query, so there could be a middle splice */
	  debug13b(printf("  querypos5 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  minus_segments[k].querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = add_bounded(minus_segments[k].diagonal,shortsplicedist,chrhigh)) > middle_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    middle_mappingend_last = mappingpos;
	    middle_mappingend_p = true;
	    debug13(printf("  Redefining middle mappingend last to %u\n",middle_mappingend_last - chroffset));
	  }
	  
	} else {
	  debug13b(printf("  querypos5 %d < %d + %d, so using this diagonal\n",
			  minus_segments[k].querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = minus_segments[k].diagonal) > close_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    close_mappingend_last = mappingpos;
	    close_mappingend_p = true;
	    debug13(printf("  Redefining close mappingend last to %u\n",close_mappingend_last - chroffset));
	  }
	}
      }
      
      /* 2 */
      if (close_mappingstart_p == true) {
	close_knownsplice_limit_low = subtract_bounded(close_mappingstart_last,shortsplicedist,chroffset);
      } else if (middle_mappingstart_p == true) {
	debug13(printf("Using middle mappingstart\n"));
	close_knownsplice_limit_low = middle_mappingstart_last;
	close_mappingstart_last = middle_mappingstart_last;
	close_mappingstart_p = true;
      }
      if (middle_mappingstart_p == true && middle_mappingstart_last < close_mappingstart_last) {
	knownsplice_limit_low = middle_mappingstart_last;
	mappingstart = middle_mappingstart_last;
      }
      if (close_mappingstart_p == false) {
	fallback_mappingstart_p = false;
      } else {
	debug13(printf("Fallback mappingstart = %u\n",mappingstart - chroffset));
	fallback_mappingstart_p = true;
      }
      
      /* 3 */
      if (close_mappingend_p == true) {
	close_knownsplice_limit_high = add_bounded(close_mappingend_last,shortsplicedist,chrhigh);
      } else if (middle_mappingend_p == true) {
	debug13(printf("Using middle mappingend\n"));
	close_knownsplice_limit_high = middle_mappingend_last;
	close_mappingend_last = middle_mappingend_last;
	close_mappingend_p = true;
      }
      
      if (middle_mappingend_p == true && middle_mappingend_last > close_mappingend_last) {
	knownsplice_limit_high = middle_mappingend_last;
	mappingend = middle_mappingend_last;
      }
      if (close_mappingend_p == false) {
	fallback_mappingend_p = false;
      } else {
	debug13(printf("Fallback mappingend = %u\n",mappingend - chroffset));
	fallback_mappingend_p = true;
      }
      
      /* 4 */
      if (close_mappingstart_p == true && close_mappingend_p == true) {
	debug13(printf("Single hit: Running gmap with close mappingstart and close mappingend\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				   /*sense_try*/0,favor_right_p,query_compress_fwd,query_compress_rev,
				   close_mappingstart_last,close_mappingend_last,
				   close_knownsplice_limit_low,close_knownsplice_limit_high,
				   /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	
	if (good_start_p == true && good_end_p == true) {
	  /* Success */
	} else if (gmap_rerun_p == false) {
	  debug13(printf("Skipping re-run of gmap\n"));
	} else if (good_start_p == true) {
	  if (fallback_mappingend_p == true) {
	    debug13(printf("Single hit: Re-running gmap with close mappingstart only\n"));
	    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				       /*sense_try*/0,favor_right_p,
				       query_compress_fwd,query_compress_rev,close_mappingstart_last,mappingend,
				       close_knownsplice_limit_low,knownsplice_limit_high,
				       /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				       oligoindices_major,oligoindices_minor,
				       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	  }
	} else if (good_end_p == true) {
	  if (fallback_mappingstart_p == true) {
	    debug13(printf("Single hit: Re-running gmap with close mappingend only\n"));
	    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				       /*sense_try*/0,favor_right_p,
				       query_compress_fwd,query_compress_rev,mappingstart,close_mappingend_last,
				       knownsplice_limit_low,close_knownsplice_limit_high,
				       /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				       oligoindices_major,oligoindices_minor,
				       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	  }
	} else {
	  if (fallback_mappingstart_p == true && fallback_mappingend_p == true) {
	    debug13(printf("Single hit: Re-running gmap with far mappingstart and mappingend\n"));
	    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				       /*sense_try*/0,favor_right_p,
				       query_compress_fwd,query_compress_rev,mappingstart,mappingend,
				       knownsplice_limit_low,close_knownsplice_limit_high,
				       /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				       oligoindices_major,oligoindices_minor,
				       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	  }
	}
	
      } else if (close_mappingstart_p == true) {
	debug13(printf("Single hit: Running gmap with close mappingstart\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				   /*sense_try*/0,favor_right_p,
				   query_compress_fwd,query_compress_rev,close_mappingstart_last,mappingend,
				   close_knownsplice_limit_low,knownsplice_limit_high,
				   /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	if (good_start_p == true) {
	  /* Success */
	} else if (gmap_rerun_p == false) {
	  debug13(printf("Skipping re-run of gmap\n"));
	} else if (fallback_mappingstart_p == true) {
	  debug13(printf("Single hit: Re-running gmap with far mappingstart\n"));
	  hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				     /*sense_try*/0,favor_right_p,
				     query_compress_fwd,query_compress_rev,mappingstart,mappingend,
				     knownsplice_limit_low,knownsplice_limit_high,
				     /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				     oligoindices_major,oligoindices_minor,
				     pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	}
	
      } else if (close_mappingend_p == true) {
	debug13(printf("Single hit: Running gmap with close mappingend\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				   /*sense_try*/0,favor_right_p,
				   query_compress_fwd,query_compress_rev,mappingstart,close_mappingend_last,
				   knownsplice_limit_low,close_knownsplice_limit_high,
				   /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	if (good_end_p == true) {
	  /* Success */
	} else if (gmap_rerun_p == false) {
	  debug13(printf("Skipping re-run of gmap\n"));
	} else if (fallback_mappingend_p == true) {
	  debug13(printf("Single hit: Re-running gmap with far mappingend\n"));
	  hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				     /*sense_try*/0,favor_right_p,
				     query_compress_fwd,query_compress_rev,mappingstart,mappingend,
				     knownsplice_limit_low,knownsplice_limit_high,
				     /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				     oligoindices_major,oligoindices_minor,
				     pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
	}
	
      } else {
	debug13(printf("Single hit: Running gmap with far mappingstart and mappingend\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,hits,accession,queryuc_ptr,querylength,
				   /*sense_try*/0,favor_right_p,
				   query_compress_fwd,query_compress_rev,mappingstart,mappingend,
				   knownsplice_limit_low,knownsplice_limit_high,
				   /*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
      }
    }
  }
    
  return hits;
}
#endif
  
  
/* Segment chaining */
static List_T
convert_plus_segments_to_gmap (List_T hits, char *queryuc_ptr, int querylength, int query_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
			       char *queryrc, bool invertedp,
#endif
			       Segment_T *anchor_segments, int nanchors,
			       struct Segment_T *plus_segments, int plus_nsegments,
			       Oligoindex_array_T oligoindices_minor,
			       Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			       Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			       int genestrand, bool require_pairing_p) {
  Univcoord_T chroffset, chrhigh, mappingpos;
  Univcoord_T origlow, orighigh;
  Univcoord_T close_mappingstart_last, close_mappingend_last,
    middle_mappingstart_last, middle_mappingend_last;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T close_knownsplice_limit_low, close_knownsplice_limit_high;
  Chrpos_T chrlength;
  Chrnum_T chrnum;
  bool close_mappingstart_p = false, close_mappingend_p = false;
  bool middle_mappingstart_p = false, middle_mappingend_p = false;
  bool novelp;		 /* Want any of the segments in startk..(endk-1) to not be used */
  bool pairablep;		/* Want any of the segments in startk..(endk-1) to be pairable */

  Segment_T anchor_segment, segment, *p;
  int anchork, startk, endk, n, i, j, firstj, lastj, k, best_starti, best_endi;

  Stage3end_T hit;  
  Pair_T *array;
  struct Pair_T *pairarray1, *pairarray2;
  List_T pairs1, pairs2, stage2pairs, unsorted_pairs;
  int querypos, boundpos, seglength;
  Chrpos_T genomepos, min_genomepos, max_genomepos;
  char c, g, g_alt;
  char *gsequence_orig, *gsequence_alt;

  Segment_T *sorted, *sorted_allocated;
  int *scores, *scores_allocated, best_score, score;
  int *prev_left, *prev_right, *prev_allocated, besti;

  int npairs;
  int cdna_direction, sensedir;
  int npairs1, goodness1, matches1, nmatches_posttrim_1,
    max_match_length_1, ambig_end_length_5_1, ambig_end_length_3_1,
    unknowns1, mismatches1, qopens1, qindels1, topens1, tindels1,
    ncanonical1, nsemicanonical1, nnoncanonical1;
  int npairs2, goodness2, matches2, nmatches_posttrim_2,
    max_match_length_2, ambig_end_length_5_2, ambig_end_length_3_2,
    unknowns2, mismatches2, qopens2, qindels2, topens2, tindels2,
    ncanonical2, nsemicanonical2, nnoncanonical2;
  double ambig_prob_5_1, ambig_prob_3_1, min_splice_prob_1;
  double ambig_prob_5_2, ambig_prob_3_2, min_splice_prob_2;
  Splicetype_T ambig_splicetype_5_1, ambig_splicetype_3_1;
  Splicetype_T ambig_splicetype_5_2, ambig_splicetype_3_2;
  Univcoord_T start, end, left;
  int nsegments, nmismatches_whole, nindels, nintrons, nindelbreaks;


  if (plus_nsegments > 0) {
#ifdef HAVE_ALLOCA
    if (plus_nsegments < MAX_ALLOCATION) {
      prev_allocated = (int *) ALLOCA(plus_nsegments*sizeof(int));
      scores_allocated = (int *) ALLOCA(plus_nsegments*sizeof(int));
      sorted_allocated = (Segment_T *) ALLOCA(plus_nsegments*sizeof(Segment_T));
    } else {
      prev_allocated = (int *) MALLOC(plus_nsegments*sizeof(int));
      scores_allocated = (int *) MALLOC(plus_nsegments*sizeof(int));
      sorted_allocated = (Segment_T *) MALLOC(plus_nsegments*sizeof(Segment_T));
    }
    gsequence_orig = (char *) MALLOCA((querylength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOCA((querylength+1) * sizeof(char));
#else
    prev_allocated = (int *) MALLOC(plus_nsegments*sizeof(int));
    scores_allocated = (int *) MALLOC(plus_nsegments*sizeof(int));
    sorted_allocated = (Segment_T *) MALLOC(plus_nsegments*sizeof(Segment_T));
    gsequence_orig = (char *) MALLOC((querylength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOC((querylength+1) * sizeof(char));
#endif
  }

  anchork = 0;
  for (p = &(anchor_segments[0]); p < &(anchor_segments[nanchors]); p++) {
    anchor_segment = *p;
    assert(anchor_segment->diagonal != (Univcoord_T) -1);
    while (plus_segments[anchork].diagonal != anchor_segment->diagonal) {
      anchork++;
    }

    startk = anchork - 1;
    while (startk >= 0 && plus_segments[startk].diagonal != (Univcoord_T) -1 &&
	   plus_segments[startk].diagonal + shortsplicedist > anchor_segment->diagonal) {
      startk--;
    }

    endk = anchork + 1;
    while (endk < plus_nsegments && plus_segments[endk].diagonal < anchor_segment->diagonal + shortsplicedist) {
      endk++;
    }
    debug13(printf("Found plus segments %d to %d inclusive for anchor %d\n",startk+1,endk-1,anchork));


    /* Dynamic programming on left (low) side (querypos5) */
    if ((n = (anchork - 1) - (startk + 1) + 1) == 0) {
      debug13(printf("On querypos5 side, n == 0\n"));
      best_starti = -1;
    } else {
      prev_left = &(prev_allocated[startk+1]);
      scores = &(scores_allocated[startk+1]);
      sorted = &(sorted_allocated[startk+1]);

      for (k = startk + 1, i = 0; k < anchork; k++) {
	sorted[i++] = &(plus_segments[k]);
      }
      qsort(sorted,n,sizeof(Segment_T),Segment_querypos5_ascending_cmp);

      lastj = 0;
      while (lastj < n && sorted[lastj]->querypos5 < anchor_segment->querypos5) {
	lastj++;
      }

      /* Go leftward to favor shorter splices */
      for (j = lastj - 1; j >= 0; --j) {
	debug13(printf("querypos5 side: %d..%d %u..%u\n",
		       sorted[j]->querypos5,sorted[j]->querypos3,sorted[j]->lowpos,sorted[j]->highpos));

	best_score = 0;
	besti = -1;
	for (i = j - 1; i >= 0; --i) {
	  if (sorted[i]->lowpos >= sorted[j]->lowpos) {
	    /* Skip, since doesn't add nucleotides to left */
	  } else if (sorted[i]->highpos < sorted[j]->lowpos) {
	    if ((score = (sorted[i]->highpos - sorted[i]->lowpos)) > best_score) {
	      best_score = score;
	      besti = i;
	    }
	  } else if ((score = (sorted[j]->lowpos - sorted[i]->lowpos)) > best_score) {
	    best_score = score;
	    besti = i;
	  }
	}
	scores[j] = sorted[j]->highpos - sorted[j]->lowpos;
	debug13(printf("Best prev is %d with score %d\n",besti,best_score));
	if ((prev_left[j] = besti) >= 0) {
	  scores[j] += best_score;
	}
      }
      
      /* Anchor segment */
      best_score = 0;
      best_starti = -1;
      for (i = lastj - 1; i >= 0; --i) {
	if (sorted[i]->lowpos >= anchor_segment->lowpos) {
	  /* Skip, since doesn't add nucleotides to left */
	} else if (sorted[i]->highpos < anchor_segment->lowpos) {
	  if ((score = (sorted[i]->highpos - sorted[i]->lowpos)) > best_score) {
	    best_score = score;
	    best_starti = i;
	  }
	} else if ((score = (anchor_segment->lowpos - sorted[i]->lowpos)) > best_score) {
	  best_score = score;
	  best_starti = i;
	}
      }
    }


    /* Dynamic programming on right (high) side (querypos3) */
    if ((n = (endk - 1) - (anchork + 1) + 1) == 0) {
      debug13(printf("On querypos3 side, n == 0\n"));
      best_endi = -1;
    } else {
      prev_right = &(prev_allocated[anchork+1]);
      scores = &(scores_allocated[anchork+1]);
      sorted = &(sorted_allocated[anchork+1]);

      for (k = anchork + 1, i = 0; k < endk; k++) {
	sorted[i++] = &(plus_segments[k]);
      }
      qsort(sorted,n,sizeof(Segment_T),Segment_querypos3_ascending_cmp);

      firstj = n - 1;
      while (firstj >= 0 && sorted[firstj]->querypos3 > anchor_segment->querypos3) {
	firstj--;
      }

      /* Go rightward to favor shorter splices */
      for (j = firstj + 1; j < n; j++) {
	debug13(printf("querypos3 side: %d..%d %u..%u\n",
		       sorted[j]->querypos5,sorted[j]->querypos3,sorted[j]->lowpos,sorted[j]->highpos));

	best_score = 0;
	besti = -1;
	for (i = j + 1; i < n; i++) {
	  if (sorted[i]->highpos <= sorted[i]->highpos) {
	    /* Skip, since doesn't add nucleotides to right */
	  } else if (sorted[i]->lowpos > sorted[j]->highpos) {
	    if ((score = (sorted[i]->highpos - sorted[i]->lowpos)) > best_score) {
	      best_score = score;
	      besti = i;
	    }
	  } else if ((score = (sorted[i]->highpos - sorted[j]->highpos)) > best_score) {
	    best_score = score;
	    besti = i;
	  }
	}
	scores[j] = sorted[j]->highpos - sorted[j]->lowpos;
	debug13(printf("Best prev is %d with score %d\n",besti,best_score));
	if ((prev_right[j] = besti) >= 0) {
	  scores[j] += best_score;
	}
      }

      /* Anchor segment */
      best_score = 0;
      best_endi = -1;
      for (i = firstj + 1; i < n; i++) {
	if (sorted[i]->highpos <= anchor_segment->highpos) {
	  /* Skip, since doesn't add nucleotides to right */
	} else if (sorted[i]->lowpos > anchor_segment->highpos) {
	  if ((score = (sorted[i]->highpos - sorted[i]->lowpos)) > best_score) {
	    best_score = score;
	    best_endi = i;
	  }
	} else if ((score = (sorted[i]->highpos - anchor_segment->highpos)) > best_score) {
	  best_score = score;
	  best_endi = i;
	}
      }
    }


    /* Evaluate set of segments */
    novelp = pairablep = false;
    if (anchor_segment->usedp == false) {
      novelp = true;
    }
    if (anchor_segment->pairablep == true) {
      pairablep = true;
    }

    sorted = &(sorted_allocated[startk+1]);
    for (k = best_starti; k >= 0; k = prev_left[k]) {
      if (sorted[k]->usedp == false) {
	novelp = true;
      }
      if (sorted[k]->pairablep == true) {
	pairablep = true;
      }
    }

    sorted = &(sorted_allocated[anchork+1]);
    for (k = best_endi; k >= 0; k = prev_right[k]) {
      if (sorted[k]->usedp == false) {
	novelp = true;
      }
      if (sorted[k]->pairablep == true) {
	pairablep = true;
      }
    }

    debug13(printf("Processing plus segments %d to %d inclusive: novelp %d, pairablep %d\n",
		   startk+1,endk-1,novelp,pairablep));
    if (novelp == true && (pairablep == true || require_pairing_p == false)) {
      anchor_segment->usedp = true;
      chrnum = anchor_segment->chrnum;
      chroffset = anchor_segment->chroffset;
      chrhigh = anchor_segment->chrhigh;
      chrlength = anchor_segment->chrlength;
      
      left = anchor_segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      origlow = left - anchor_segment->querypos5;
      orighigh = left + (querylength - anchor_segment->querypos3);

      /* extend left */
      knownsplice_limit_low = subtract_bounded(origlow,shortsplicedist,chroffset);
      debug13(printf("Original bounds A: knownsplice_limit_low %u\n",knownsplice_limit_low - chroffset));
      
      /* extend right */
      knownsplice_limit_high = add_bounded(orighigh,shortsplicedist,chrhigh);
      debug13(printf("Original bounds B: knownsplice_limit_high %u\n",knownsplice_limit_high - chroffset));
      
      close_mappingstart_last = middle_mappingstart_last = origlow;
      close_mappingend_last = middle_mappingend_last = orighigh;
      close_mappingstart_p = close_mappingend_p = false;
      middle_mappingstart_p = middle_mappingend_p = false;

      /* 1 */
      sorted = &(sorted_allocated[startk+1]);
      for (k = best_starti; k >= 0; k = prev_left[k]) {
	segment = sorted[k];
	segment->usedp = true;
	debug13(printf("1. plus diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (segment->diagonal - chroffset),(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,segment->usedp,segment->pairablep));
	if (segment->querypos5 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 3. Missing start of query, so there could be a middle splice */
	  debug13b(printf("  querypos5 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  segment->querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(segment->diagonal,querylength + shortsplicedist,chroffset)) < middle_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    middle_mappingstart_last = mappingpos;
	    middle_mappingstart_p = true;
	    debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	  }

	} else {
	  debug13b(printf("  querypos5 %d < %d + %d, so using this diagonal\n",
			  segment->querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(segment->diagonal,querylength,chroffset)) < close_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    close_mappingstart_last = mappingpos;
	    close_mappingstart_p = true;
	    debug13(printf("  Redefining close mappingstart last to %u\n",close_mappingstart_last - chroffset));
	  }
	}


	if (query_lastpos - segment->querypos3 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 1. Missing end of query, so there could be a middle splice */
	  debug13b(printf("  query_lastpos %d - querypos3 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  query_lastpos,segment->querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = add_bounded(segment->diagonal,shortsplicedist,chrhigh)) > middle_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    middle_mappingend_last = mappingpos;
	    middle_mappingend_p = true;
	    debug13(printf("  Redefining middle mappingend last to %u\n",middle_mappingend_last - chroffset));
	  }

	} else {
	  debug13b(printf("  query_lastpos %d - querypos3 %d < %d + %d, so using this diagonal\n",
			  query_lastpos,segment->querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = segment->diagonal) > close_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    close_mappingend_last = mappingpos;
	    close_mappingend_p = true;
	    debug13(printf("  Redefining close mappingend last to %u\n",close_mappingend_last - chroffset));
	  }
	}
      }

      sorted = &(sorted_allocated[anchork+1]);
      for (k = best_endi; k >= 0; k = prev_right[k]) {
	segment = sorted[k];
	segment->usedp = true;
	debug13(printf("1. plus diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (segment->diagonal - chroffset),(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,segment->usedp,segment->pairablep));
	if (segment->querypos5 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 3. Missing start of query, so there could be a middle splice */
	  debug13b(printf("  querypos5 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  segment->querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(segment->diagonal,querylength + shortsplicedist,chroffset)) < middle_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    middle_mappingstart_last = mappingpos;
	    middle_mappingstart_p = true;
	    debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	  }

	} else {
	  debug13b(printf("  querypos5 %d < %d + %d, so using this diagonal\n",
			  segment->querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(segment->diagonal,querylength,chroffset)) < close_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    close_mappingstart_last = mappingpos;
	    close_mappingstart_p = true;
	    debug13(printf("  Redefining close mappingstart last to %u\n",close_mappingstart_last - chroffset));
	  }
	}


	if (query_lastpos - segment->querypos3 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 1. Missing end of query, so there could be a middle splice */
	  debug13b(printf("  query_lastpos %d - querypos3 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  query_lastpos,segment->querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = add_bounded(segment->diagonal,shortsplicedist,chrhigh)) > middle_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    middle_mappingend_last = mappingpos;
	    middle_mappingend_p = true;
	    debug13(printf("  Redefining middle mappingend last to %u\n",middle_mappingend_last - chroffset));
	  }

	} else {
	  debug13b(printf("  query_lastpos %d - querypos3 %d < %d + %d, so using this diagonal\n",
			  query_lastpos,segment->querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = segment->diagonal) > close_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    close_mappingend_last = mappingpos;
	    close_mappingend_p = true;
	    debug13(printf("  Redefining close mappingend last to %u\n",close_mappingend_last - chroffset));
	  }
	}
      }


      /* 2 */
      if (close_mappingstart_p == true) {
	close_knownsplice_limit_low = subtract_bounded(close_mappingstart_last,shortsplicedist,chroffset);
      } else if (middle_mappingstart_p == true) {
	debug13(printf("Using middle mappingstart\n"));
	close_knownsplice_limit_low = middle_mappingstart_last;
	close_mappingstart_last = middle_mappingstart_last;
	close_mappingstart_p = true;
      }

      if (middle_mappingstart_p == true && middle_mappingstart_last < close_mappingstart_last) {
	knownsplice_limit_low = middle_mappingstart_last;
      }

      /* 3 */
      if (close_mappingend_p == true) {
	close_knownsplice_limit_high = add_bounded(close_mappingend_last,shortsplicedist,chrhigh);
      } else if (middle_mappingend_p == true) {
	close_knownsplice_limit_high = middle_mappingend_last;
	close_mappingend_last = middle_mappingend_last;
	close_mappingend_p = true;
	debug13(printf("Using middle mappingend => close_mappingend %u\n",close_mappingend_last));
      }
      if (middle_mappingend_p == true && middle_mappingend_last > close_mappingend_last) {
	knownsplice_limit_high = middle_mappingend_last;
      }

      /* 4 */
      if (close_mappingstart_p == true) {
	knownsplice_limit_low = close_knownsplice_limit_low;
      }
      if (close_mappingend_p == true) {
	knownsplice_limit_high = close_knownsplice_limit_high;
      }


      /* F.  Make stage2pairs (anchor) */
      unsorted_pairs = (List_T) NULL;

      debug13(printf("plus anchor diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		     (Chrpos_T) (anchor_segment->diagonal - chroffset),(unsigned long long) anchor_segment->diagonal,
		     anchor_segment->querypos5,anchor_segment->querypos3,anchor_segment->usedp,anchor_segment->pairablep));
      querypos = anchor_segment->querypos5;
      seglength = (anchor_segment->querypos3 + index1part) - querypos;
      
      left = anchor_segment->diagonal - querylength; /* FORMULA */
      min_genomepos = genomepos = (left - chroffset) + querypos;
      Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/chroffset+genomepos,
				      seglength,chrhigh,/*revcomp*/false);
      
      for (i = 0; i < seglength; i++) {
	c = queryuc_ptr[querypos];
	g = gsequence_orig[i];
	g_alt = gsequence_alt[i];
	if (g == c || g_alt == c) {
	  debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MATCH_COMP,g,querypos,genomepos));
	  unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					 /*cdna*/c,/*comp*/MATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					 /*dynprogindex*/0);
#if 0
	} else {
	  /* Let stage 3 handle mismatches */
	  debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MISMATCH_COMP,g,querypos,genomepos));
	  unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					 /*cdna*/c,/*comp*/MISMATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					 /*dynprogindex*/0);
#endif
	}
	querypos++;
	genomepos++;
      }
      max_genomepos = genomepos - 1;

      /* F.  Make stage2pairs (left) */
      sorted = &(sorted_allocated[startk+1]);
      boundpos = anchor_segment->querypos5;
      for (k = best_starti; k >= 0; k = prev_left[k]) {
	segment = sorted[k];
	debug13(printf("plus left diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (segment->diagonal - chroffset),(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,segment->usedp,segment->pairablep));

	querypos = segment->querypos5;
	if (querypos < boundpos) {
	  left = segment->diagonal - querylength; /* FORMULA */
	  genomepos = (left - chroffset) + querypos;
	  if (genomepos < min_genomepos) {
	    seglength = (segment->querypos3 + index1part) - querypos;
	    Genome_get_segment_blocks_left(gsequence_orig,gsequence_alt,/*right*/chroffset+genomepos+seglength,
					   seglength,chroffset,/*revcomp*/false);
	    
	    i = 0;
	    while (i < seglength && querypos < boundpos && genomepos < min_genomepos) {
	      c = queryuc_ptr[querypos];
	      g = gsequence_orig[i];
	      g_alt = gsequence_alt[i];
	      if (g == c || g_alt == c) {
		debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MATCH_COMP,g,querypos,genomepos));
		unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					       /*cdna*/c,/*comp*/MATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					       /*dynprogindex*/0);
#if 0
	      } else {
		/* Let stage 3 handle mismatches */
		debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MISMATCH_COMP,g,querypos,genomepos));
		unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					       /*cdna*/c,/*comp*/MISMATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					       /*dynprogindex*/0);
#endif
	      }
	      querypos++;
	      genomepos++;
	      i++;
	    }
	    boundpos = segment->querypos5;
	    min_genomepos = (left - chroffset) + segment->querypos5;
	  }
	}
      }

      /* F.  Make stage2pairs (right) */
      sorted = &(sorted_allocated[anchork+1]);
      boundpos = anchor_segment->querypos3 + index1part;
      for (k = best_endi; k >= 0; k = prev_right[k]) {
	segment = sorted[k];
	debug13(printf("plus right diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (segment->diagonal - chroffset),(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,segment->usedp,segment->pairablep));
	querypos = segment->querypos5;
	seglength = (segment->querypos3 + index1part) - querypos;
      
	left = segment->diagonal - querylength; /* FORMULA */
	genomepos = left - chroffset + querypos;
	Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/chroffset+genomepos,
					seglength,chrhigh,/*revcomp*/false);
      
	i = 0;
	while (i < seglength && (querypos <= boundpos || genomepos <= max_genomepos)) {
	  querypos++;
	  genomepos++;
	  i++;
	}

	while (i < seglength) {
	  c = queryuc_ptr[querypos];
	  g = gsequence_orig[i];
	  g_alt = gsequence_alt[i];
	  if (g == c || g_alt == c) {
	    debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MATCH_COMP,g,querypos,genomepos));
	    unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					   /*cdna*/c,/*comp*/MATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					   /*dynprogindex*/0);
#if 0
	  } else {
	    /* Let stage 3 handle mismatches */
	    debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MISMATCH_COMP,g,querypos,genomepos));
	    unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					   /*cdna*/c,/*comp*/MISMATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					   /*dynprogindex*/0);
#endif
	  }
	  querypos++;
	  genomepos++;
	  i++;
	}
	boundpos = segment->querypos3 + index1part;
	max_genomepos = genomepos - 1;
      }


      /* Sort pairs and get unique ones */
      array = (Pair_T *) List_to_array_n(&npairs,unsorted_pairs);
      qsort(array,npairs,sizeof(Pair_T),Pair_cmp);
	
      stage2pairs = (List_T) NULL;
      i = 0;
      while (i < npairs) {
	j = i + 1;
	while (j < npairs && array[j]->querypos == array[i]->querypos) {
	  j++;
	}
	if (j == i + 1) {
	  /* Only a single pair at this querypos */
	  debug13(Pair_dump_one(array[i],true));
	  debug13(printf("\n"));
	  stage2pairs = Pairpool_push_existing(stage2pairs,pairpool,array[i]);
	}
	i = j;
      }
      stage2pairs = List_reverse(stage2pairs);
      FREE(array);


      /* Run GMAP */
      if (stage2pairs == NULL) {
	/* hit = (T) NULL; */
      } else if ((pairarray1 = Stage3_compute(&cdna_direction,&sensedir,&pairs1,&npairs1,&goodness1,
					      &matches1,&nmatches_posttrim_1,&max_match_length_1,
					      &ambig_end_length_5_1,&ambig_end_length_3_1,
					      &ambig_splicetype_5_1,&ambig_splicetype_3_1,
					      &ambig_prob_5_1,&ambig_prob_3_1,
					      &unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
					      &ncanonical1,&nsemicanonical1,&nnoncanonical1,&min_splice_prob_1,

					      &pairarray2,&pairs2,&npairs2,&goodness2,
					      &matches2,&nmatches_posttrim_2,&max_match_length_2,
					      &ambig_end_length_5_2,&ambig_end_length_3_2,
					      &ambig_splicetype_5_2,&ambig_splicetype_3_2,
					      &ambig_prob_5_2,&ambig_prob_3_2,
					      &unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
					      &ncanonical2,&nsemicanonical2,&nnoncanonical2,&min_splice_prob_2,

					      stage2pairs,/*all_stage2_starts*/NULL,/*all_stage2_ends*/NULL,
#ifdef END_KNOWNSPLICING_SHORTCUT
					      cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
					      watsonp ? query_compress_fwd : query_compress_rev,
#endif
					      /*queryseq_ptr*/queryuc_ptr,queryuc_ptr,querylength,/*skiplength*/0,
#ifdef EXTRACT_GENOMICSEG
					      /*query_subseq_offset*/0,
#else
					      /*query_subseq_offset*/0,
#endif
					      chrnum,chroffset,chrhigh,
					      knownsplice_limit_low,knownsplice_limit_high,/*plusp*/true,genestrand,
					      /*jump_late_p*/false,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
					      /*sense_try*/0,/*sense_filter*/0,
					      oligoindices_minor,diagpool,cellpool)) == NULL) {
	/* hit = (T) NULL; */

      } else if (cdna_direction == 0) {
	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray1,npairs1);
	start = subtract_bounded(chroffset + Pair_genomepos(&(pairarray1[0])),
				 /*minusterm*/Pair_querypos(&(pairarray1[0])),chroffset);
	end = add_bounded(chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			  /*plusterm*/querylength - 1 - Pair_querypos(&(pairarray1[npairs1-1])),chrhigh);

	if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_1,max_match_length_1,
				      ambig_end_length_5_1,ambig_end_length_3_1,
				      ambig_splicetype_5_1,ambig_splicetype_3_1,
				      min_splice_prob_1,
				      pairarray1,npairs1,nsegments,nintrons,nindelbreaks,
				      /*left*/start,/*genomiclength*/end - start + 1,
				      /*plusp*/true,genestrand,
				      /*accession*/NULL,querylength,chrnum,chroffset,chrhigh,chrlength,
				      /*cdna_direction*/+1,/*sensedir*/SENSE_FORWARD,
				      /*gmap_source*/GMAP_VIA_SEGMENTS)) == NULL) {

	  FREE_OUT(pairarray1);
	} else {
	  hits = List_push(hits,(void *) hit);
	}

	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray2,npairs2);
	start = subtract_bounded(chroffset + Pair_genomepos(&(pairarray2[0])),
				 /*minusterm*/Pair_querypos(&(pairarray2[0])),chroffset);
	end = add_bounded(chroffset + Pair_genomepos(&(pairarray2[npairs2-1])),
			  /*plusterm*/querylength - 1 - Pair_querypos(&(pairarray2[npairs2-1])),chrhigh);

	if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_2,max_match_length_2,
				      ambig_end_length_5_2,ambig_end_length_3_2,
				      ambig_splicetype_5_2,ambig_splicetype_3_2,
				      min_splice_prob_2,
				      pairarray2,npairs2,nsegments,nintrons,nindelbreaks,
				      /*left*/start,/*genomiclength*/end - start + 1,
				      /*plusp*/true,genestrand,
				      /*accession*/NULL,querylength,chrnum,chroffset,chrhigh,chrlength,
				      /*cdna_direction*/-1,/*sensedir*/SENSE_ANTI,
				      /*gmap_source*/GMAP_VIA_SEGMENTS)) == NULL) {

	  FREE_OUT(pairarray2);
	} else {
	  hits = List_push(hits,(void *) hit);
	}

      } else {
	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray1,npairs1);
	start = subtract_bounded(chroffset + Pair_genomepos(&(pairarray1[0])),
				 /*minusterm*/Pair_querypos(&(pairarray1[0])),chroffset);
	end = add_bounded(chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			  /*plusterm*/querylength - 1 - Pair_querypos(&(pairarray1[npairs1-1])),chrhigh);

	if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_1,max_match_length_1,
				      ambig_end_length_5_1,ambig_end_length_3_1,
				      ambig_splicetype_5_1,ambig_splicetype_3_1,
				      min_splice_prob_1,
				      pairarray1,npairs1,nsegments,nintrons,nindelbreaks,
				      /*left*/start,/*genomiclength*/end - start + 1,
				      /*plusp*/true,genestrand,
				      /*accession*/NULL,querylength,chrnum,chroffset,chrhigh,chrlength,
				      cdna_direction,sensedir,/*gmap_source*/GMAP_VIA_SEGMENTS)) == NULL) {

	  FREE_OUT(pairarray1);
	} else {
	  hits = List_push(hits,(void *) hit);
	}
      }
    }
  }

  if (plus_nsegments > 0) {
#ifdef HAVE_ALLOCA
    FREEA(gsequence_alt);
    FREEA(gsequence_orig);
    if (plus_nsegments < MAX_ALLOCATION) {
      FREEA(sorted_allocated);
      FREEA(scores_allocated);
      FREEA(prev_allocated);
    } else {
      FREE(sorted_allocated);
      FREE(scores_allocated);
      FREE(prev_allocated);
    }
#else
    FREE(gsequence_alt);
    FREE(gsequence_orig);
    FREE(sorted_allocated);
    FREE(scores_allocated);
    FREE(prev_allocated);
#endif
  }

  return hits;
}


/* Segment chaining */
static List_T
convert_minus_segments_to_gmap (List_T hits, char *queryuc_ptr, int querylength, int query_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
				char *queryrc, bool invertedp,
#endif
				Segment_T *anchor_segments, int nanchors,
				struct Segment_T *minus_segments, int minus_nsegments,
				Oligoindex_array_T oligoindices_minor,
				Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
				Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				int genestrand, bool require_pairing_p) {
  Univcoord_T chroffset, chrhigh, mappingpos;
  Univcoord_T origlow, orighigh;
  Univcoord_T close_mappingstart_last, close_mappingend_last,
    middle_mappingstart_last, middle_mappingend_last;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T close_knownsplice_limit_low, close_knownsplice_limit_high;
  Chrpos_T chrlength;
  Chrnum_T chrnum;
  bool close_mappingstart_p, close_mappingend_p;
  bool middle_mappingstart_p, middle_mappingend_p;
  bool novelp;		 /* Want any of the segments in startk..(endk-1) to not be used */
  bool pairablep;		/* Want any of the segments in startk..(endk-1) to be pairable */

  Segment_T anchor_segment, segment, *p;
  int anchork, startk, endk, n, i, j, firstj, lastj, k, best_starti, best_endi;

  Stage3end_T hit;  
  Pair_T *array;
  struct Pair_T *pairarray1, *pairarray2;
  List_T pairs1, pairs2, stage2pairs, unsorted_pairs;
  int querypos, boundpos, seglength;
  Chrpos_T genomepos, min_genomepos, max_genomepos;
  char c, g, g_alt;
  char *gsequence_orig, *gsequence_alt;

  Segment_T *sorted, *sorted_allocated;
  int *scores, *scores_allocated, best_score, score;
  int *prev_left, *prev_right, *prev_allocated, besti;

  int npairs;
  int cdna_direction, sensedir;
  int npairs1, goodness1, matches1, nmatches_posttrim_1,
    max_match_length_1, ambig_end_length_5_1, ambig_end_length_3_1,
    unknowns1, mismatches1, qopens1, qindels1, topens1, tindels1,
    ncanonical1, nsemicanonical1, nnoncanonical1;
  int npairs2, goodness2, matches2, nmatches_posttrim_2,
    max_match_length_2, ambig_end_length_5_2, ambig_end_length_3_2,
    unknowns2, mismatches2, qopens2, qindels2, topens2, tindels2,
    ncanonical2, nsemicanonical2, nnoncanonical2;
  double ambig_prob_5_1, ambig_prob_3_1, min_splice_prob_1;
  double ambig_prob_5_2, ambig_prob_3_2, min_splice_prob_2;
  Splicetype_T ambig_splicetype_5_1, ambig_splicetype_3_1;
  Splicetype_T ambig_splicetype_5_2, ambig_splicetype_3_2;
  Univcoord_T start, end, left;
  int nsegments, nmismatches_whole, nindels, nintrons, nindelbreaks;


  if (minus_nsegments > 0) {
#ifdef HAVE_ALLOCA
    if (minus_nsegments < MAX_ALLOCATION) {
      prev_allocated = (int *) ALLOCA(minus_nsegments*sizeof(int));
      scores_allocated = (int *) ALLOCA(minus_nsegments*sizeof(int));
      sorted_allocated = (Segment_T *) ALLOCA(minus_nsegments*sizeof(Segment_T));
    } else {
      prev_allocated = (int *) MALLOC(minus_nsegments*sizeof(int));
      scores_allocated = (int *) MALLOC(minus_nsegments*sizeof(int));
      sorted_allocated = (Segment_T *) MALLOC(minus_nsegments*sizeof(Segment_T));
    }
    gsequence_orig = (char *) MALLOCA((querylength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOCA((querylength+1) * sizeof(char));
#else
    prev_allocated = (int *) MALLOC(minus_nsegments*sizeof(int));
    scores_allocated = (int *) MALLOC(minus_nsegments*sizeof(int));
    sorted_allocated = (Segment_T *) MALLOC(minus_nsegments*sizeof(Segment_T));
    gsequence_orig = (char *) MALLOC((querylength+1) * sizeof(char));
    gsequence_alt = (char *) MALLOC((querylength+1) * sizeof(char));
#endif
  }

  anchork = 0;
  for (p = &(anchor_segments[0]); p < &(anchor_segments[nanchors]); p++) {
    anchor_segment = *p;
    assert(anchor_segment->diagonal != (Univcoord_T) -1);
    while (minus_segments[anchork].diagonal != anchor_segment->diagonal) {
      anchork++;
    }

    startk = anchork - 1;
    while (startk >= 0 && minus_segments[startk].diagonal != (Univcoord_T) -1 &&
	   minus_segments[startk].diagonal + shortsplicedist > anchor_segment->diagonal) {
      startk--;
    }

    endk = anchork + 1;
    while (endk < minus_nsegments && minus_segments[endk].diagonal < anchor_segment->diagonal + shortsplicedist) {
      endk++;
    }
    debug13(printf("Found minus segments %d to %d inclusive for anchor %d\n",startk+1,endk-1,anchork));


    /* Dynamic programming on left (low) side (querypos3) */
    if ((n = (anchork - 1) - (startk + 1) + 1) == 0) {
      debug13(printf("On querypos5 side, n == 0\n"));
      best_starti = -1;
    } else {
      prev_left = &(prev_allocated[startk+1]);
      scores = &(scores_allocated[startk+1]);
      sorted = &(sorted_allocated[startk+1]);

      for (k = startk + 1, i = 0; k < anchork; k++) {
	sorted[i++] = &(minus_segments[k]);
      }
      qsort(sorted,n,sizeof(Segment_T),Segment_querypos3_descending_cmp);

      lastj = 0;
      while (lastj < n && sorted[lastj]->querypos3 > anchor_segment->querypos3) {
	lastj++;
      }

      /* Go rightward to favor shorter splices */
      for (j = lastj - 1; j >= 0; --j) {
	debug13(printf("querypos3 side: %d..%d %u..%u\n",
		       sorted[j]->querypos5,sorted[j]->querypos3,sorted[j]->lowpos,sorted[j]->highpos));

	best_score = 0;
	besti = -1;
	for (i = j - 1; i >= 0; --i) {
	  if (sorted[i]->lowpos >= sorted[j]->lowpos) {
	    /* Skip, since doesn't add nucleotides to left */
	  } else if (sorted[i]->highpos < sorted[j]->lowpos) {
	    if ((score = (sorted[i]->highpos - sorted[i]->lowpos)) > best_score) {
	      best_score = score;
	      besti = i;
	    }
	  } else if ((score = (sorted[j]->lowpos - sorted[i]->lowpos)) > best_score) {
	    best_score = score;
	    besti = i;
	  }
	}
	scores[j] = sorted[j]->highpos - sorted[j]->lowpos;
	debug13(printf("Best prev is %d with score %d\n",besti,best_score));
	if ((prev_left[j] = besti) >= 0) {
	  scores[j] += best_score;
	}
      }

      /* Anchor segment */
      best_score = 0;
      best_starti = -1;
      for (i = lastj - 1; i >= 0; --i) {
	if (sorted[i]->lowpos >= anchor_segment->lowpos) {
	  /* Skip, since doesn't add nucleotides to left */
	} else if (sorted[i]->highpos < anchor_segment->lowpos) {
	  if ((score = (sorted[i]->highpos - sorted[i]->lowpos)) > best_score) {
	    best_score = score;
	    best_starti = i;
	  }
	} else if ((score = (anchor_segment->lowpos - sorted[i]->lowpos)) > best_score) {
	  best_score = score;
	  best_starti = i;
	}
      }
    }


    /* Dynamic programming on right (high) side (querypos5) */
    if ((n = (endk - 1) - (anchork + 1) + 1) == 0) {
      debug13(printf("On querypos5 side, n == 0\n"));
      best_endi = -1;
    } else {
      prev_right = &(prev_allocated[anchork+1]);
      scores = &(scores_allocated[anchork+1]);
      sorted = &(sorted_allocated[anchork+1]);

      for (k = anchork + 1, i = 0; k < endk; k++) {
	sorted[i++] = &(minus_segments[k]);
      }
      qsort(sorted,n,sizeof(Segment_T),Segment_querypos5_descending_cmp);

      firstj = n - 1;
      while (firstj >= 0 && sorted[firstj]->querypos5 < anchor_segment->querypos5) {
	firstj--;
      }

      /* Go leftward to favor shorter splices */
      for (j = firstj + 1; j < n; j++) {
	debug13(printf("querypos5 side: %d..%d %u..%u\n",
		       sorted[j]->querypos5,sorted[j]->querypos3,sorted[j]->lowpos,sorted[j]->highpos));

	best_score = 0;
	besti = -1;
	for (i = j + 1; i < n; i++) {
	  if (sorted[i]->highpos <= sorted[i]->highpos) {
	    /* Skip, since doesn't add nucleotides to right */
	  } else if (sorted[i]->lowpos > sorted[j]->highpos) {
	    if ((score = (sorted[i]->highpos - sorted[i]->lowpos)) > best_score) {
	      best_score = score;
	      besti = i;
	    }
	  } else if ((score = (sorted[i]->highpos - sorted[j]->highpos)) > best_score) {
	    best_score = score;
	    besti = i;
	  }
	}
	scores[j] = sorted[j]->highpos - sorted[j]->lowpos;
	debug13(printf("Best prev is %d with score %d\n",besti,best_score));
	if ((prev_right[j] = besti) >= 0) {
	  scores[j] += best_score;
	}
      }

      /* Anchor segment */
      best_score = 0;
      best_endi = -1;
      for (i = firstj + 1; i < n; i++) {
	if (sorted[i]->highpos <= anchor_segment->highpos) {
	  /* Skip, since doesn't add nucleotides to right */
	} else if (sorted[i]->lowpos > anchor_segment->highpos) {
	  if ((score = (sorted[i]->highpos - sorted[i]->lowpos)) > best_score) {
	    best_score = score;
	    best_endi = i;
	  }
	} else if ((score = (sorted[i]->highpos - anchor_segment->highpos)) > best_score) {
	  best_score = score;
	  best_endi = i;
	}
      }
    }


    /* Evaluate set of segments */
    novelp = pairablep = false;
    if (anchor_segment->usedp == false) {
      novelp = true;
    }
    if (anchor_segment->pairablep == true) {
      pairablep = true;
    }

    sorted = &(sorted_allocated[startk+1]);
    for (k = best_starti; k >= 0; k = prev_left[k]) {
      if (sorted[k]->usedp == false) {
	novelp = true;
      }
      if (sorted[k]->pairablep == true) {
	pairablep = true;
      }
    }

    sorted = &(sorted_allocated[anchork+1]);
    for (k = best_endi; k >= 0; k = prev_right[k]) {
      if (sorted[k]->usedp == false) {
	novelp = true;
      }
      if (sorted[k]->pairablep == true) {
	pairablep = true;
      }
    }


    debug13(printf("Processing minus segments %d to %d inclusive: novelp %d, pairablep %d\n",
		   startk+1,endk-1,novelp,pairablep));
    if (novelp == true && (pairablep == true || require_pairing_p == false)) {
      anchor_segment->usedp = true;
      chrnum = anchor_segment->chrnum;
      chroffset = anchor_segment->chroffset;
      chrhigh = anchor_segment->chrhigh;
      chrlength = anchor_segment->chrlength;

      left = anchor_segment->diagonal - querylength; /* FORMULA */
      origlow = left - (querylength - anchor_segment->querypos3);
      orighigh = left + anchor_segment->querypos5;
      
      /* extend right */
      knownsplice_limit_low = subtract_bounded(origlow,shortsplicedist,chroffset);
      debug13(printf("Original bounds C: knownsplice_limit_low %u\n",knownsplice_limit_low - chroffset));
      
      /* extend left */
      knownsplice_limit_high = add_bounded(orighigh,shortsplicedist,chrhigh);
      debug13(printf("Original bounds D: knownsplice_limit_high %u\n",knownsplice_limit_high - chroffset));
      
      close_mappingstart_last = middle_mappingstart_last = origlow;
      close_mappingend_last = middle_mappingend_last = orighigh;
      close_mappingstart_p = close_mappingend_p = false;
      middle_mappingstart_p = middle_mappingend_p = false;
      
      /* 1 */
      sorted = &(sorted_allocated[startk+1]);
      for (k = best_starti; k >= 0; k = prev_left[k]) {
	segment = sorted[k];
	segment->usedp = true;
	debug13(printf("1. minus diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (segment->diagonal - chroffset),(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,segment->usedp,segment->pairablep));
	if (query_lastpos - segment->querypos3 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 2. Missing end of query, so there could be a middle splice */
	  debug13b(printf("  query_lastpos %d - querypos3 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  query_lastpos,segment->querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(segment->diagonal,querylength + shortsplicedist,chroffset)) < middle_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    middle_mappingstart_last = mappingpos;
	    middle_mappingstart_p = true;
	    debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	  }
	  
	} else {
	  debug13b(printf("  query_lastpos %d - querypos3 %d < %d + %d, so using this diagonal\n",
			  query_lastpos,segment->querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(segment->diagonal,querylength,chroffset)) < close_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    close_mappingstart_last = mappingpos;
	    close_mappingstart_p = true;
	    debug13(printf("  Redefining close mappingstart last to %u\n",close_mappingstart_last - chroffset));
	  }
	}
      

	if (segment->querypos5 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 4. Missing start of query, so there could be a middle splice */
	  debug13b(printf("  querypos5 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  segment->querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = add_bounded(segment->diagonal,shortsplicedist,chrhigh)) > middle_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    middle_mappingend_last = mappingpos;
	    middle_mappingend_p = true;
	    debug13(printf("  Redefining middle mappingend last to %u\n",middle_mappingend_last - chroffset));
	  }
	  
	} else {
	  debug13b(printf("  querypos5 %d < %d + %d, so using this diagonal\n",
			  segment->querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = segment->diagonal) > close_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    close_mappingend_last = mappingpos;
	    close_mappingend_p = true;
	    debug13(printf("  Redefining close mappingend last to %u\n",close_mappingend_last - chroffset));
	  }
	}
      }

      sorted = &(sorted_allocated[anchork+1]);
      for (k = best_endi; k >= 0; k = prev_right[k]) {
	segment = sorted[k];
	segment->usedp = true;
	debug13(printf("1. minus diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (segment->diagonal - chroffset),(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,segment->usedp,segment->pairablep));
	if (query_lastpos - segment->querypos3 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 2. Missing end of query, so there could be a middle splice */
	  debug13b(printf("  query_lastpos %d - querypos3 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  query_lastpos,segment->querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(segment->diagonal,querylength + shortsplicedist,chroffset)) < middle_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    middle_mappingstart_last = mappingpos;
	    middle_mappingstart_p = true;
	    debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	  }
	  
	} else {
	  debug13b(printf("  query_lastpos %d - querypos3 %d < %d + %d, so using this diagonal\n",
			  query_lastpos,segment->querypos3,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = subtract_bounded(segment->diagonal,querylength,chroffset)) < close_mappingstart_last) {
	    /* Use < for NOT_GREEDY */
	    close_mappingstart_last = mappingpos;
	    close_mappingstart_p = true;
	    debug13(printf("  Redefining close mappingstart last to %u\n",close_mappingstart_last - chroffset));
	  }
	}
      

	if (segment->querypos5 >= STAGE2_MIN_OLIGO + index1interval) {
	  /* Case 4. Missing start of query, so there could be a middle splice */
	  debug13b(printf("  querypos5 %d >= %d + %d, so using this diagonal plus shortsplicedist\n",
			  segment->querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = add_bounded(segment->diagonal,shortsplicedist,chrhigh)) > middle_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    middle_mappingend_last = mappingpos;
	    middle_mappingend_p = true;
	    debug13(printf("  Redefining middle mappingend last to %u\n",middle_mappingend_last - chroffset));
	  }
	  
	} else {
	  debug13b(printf("  querypos5 %d < %d + %d, so using this diagonal\n",
			  segment->querypos5,STAGE2_MIN_OLIGO,index1interval));
	  if ((mappingpos = segment->diagonal) > close_mappingend_last) {
	    /* Use > for NOT_GREEDY */
	    close_mappingend_last = mappingpos;
	    close_mappingend_p = true;
	    debug13(printf("  Redefining close mappingend last to %u\n",close_mappingend_last - chroffset));
	  }
	}
      }
      
      /* 2 */
      if (close_mappingstart_p == true) {
	close_knownsplice_limit_low = subtract_bounded(close_mappingstart_last,shortsplicedist,chroffset);
      } else if (middle_mappingstart_p == true) {
	debug13(printf("Using middle mappingstart\n"));
	close_knownsplice_limit_low = middle_mappingstart_last;
	close_mappingstart_last = middle_mappingstart_last;
	close_mappingstart_p = true;
      }
      if (middle_mappingstart_p == true && middle_mappingstart_last < close_mappingstart_last) {
	knownsplice_limit_low = middle_mappingstart_last;
      }
      
      /* 3 */
      if (close_mappingend_p == true) {
	close_knownsplice_limit_high = add_bounded(close_mappingend_last,shortsplicedist,chrhigh);
      } else if (middle_mappingend_p == true) {
	debug13(printf("Using middle mappingend\n"));
	close_knownsplice_limit_high = middle_mappingend_last;
	close_mappingend_last = middle_mappingend_last;
	close_mappingend_p = true;
      }
      if (middle_mappingend_p == true && middle_mappingend_last > close_mappingend_last) {
	knownsplice_limit_high = middle_mappingend_last;
      }

      /* 4 */
      if (close_mappingstart_p == true) {
	knownsplice_limit_low = close_knownsplice_limit_low;
      }
      if (close_mappingend_p == true) {
	knownsplice_limit_high = close_knownsplice_limit_high;
      }


      /* F.  Make stage2pairs (anchor) */
      unsorted_pairs = (List_T) NULL;

      debug13(printf("minus anchor diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		     (Chrpos_T) (anchor_segment->diagonal - chroffset),(unsigned long long) anchor_segment->diagonal,
		     anchor_segment->querypos5,anchor_segment->querypos3,anchor_segment->usedp,anchor_segment->pairablep));
      querypos = anchor_segment->querypos5;
      seglength = (anchor_segment->querypos3 + index1part) - querypos;

      /* left = anchor_segment->diagonal - querylength; -- FORMULA */
      min_genomepos = genomepos = chrhigh - (anchor_segment->diagonal - 1) + querypos;
      Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/anchor_segment->diagonal - querypos - seglength,
                                      seglength,chrhigh,/*revcomp*/true);

      for (i = 0; i < seglength; i++) {
	c = queryuc_ptr[querypos];
	g = gsequence_orig[i];
	g_alt = gsequence_alt[i];
	if (g == c || g_alt == c) {
	  debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MATCH_COMP,g,querypos,genomepos));
	  unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					 /*cdna*/c,/*comp*/MATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					 /*dynprogindex*/0);
#if 0
	} else {
	  /* Let stage 3 handle mismatches */
	  debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MISMATCH_COMP,g,querypos,genomepos));
	  unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					 /*cdna*/c,/*comp*/MISMATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					 /*dynprogindex*/0);
#endif
	}
	querypos++;
	genomepos++;
      }
      max_genomepos = genomepos - 1;

      /* F.  Make stage2pairs (left) */
      sorted = &(sorted_allocated[startk+1]);
      boundpos = anchor_segment->querypos3 + index1part;
      for (k = best_starti; k >= 0; k = prev_left[k]) {
	segment = sorted[k];
	debug13(printf("minus left diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (segment->diagonal - chroffset),(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,segment->usedp,segment->pairablep));
	querypos = segment->querypos5;
	seglength = (segment->querypos3 + index1part) - querypos;

	/* left = segment->diagonal - querylength; -- FORMULA */
	genomepos = chrhigh - (segment->diagonal - 1) + querypos;
	Genome_get_segment_blocks_left(gsequence_orig,gsequence_alt,/*right*/segment->diagonal - querypos /*- seglength*/,
				       seglength,chroffset,/*revcomp*/true);

	i = 0;
	while (i < seglength && (querypos <= boundpos || genomepos <= max_genomepos)) {
	  querypos++;
	  genomepos++;
	  i++;
	}

	while (i < seglength) {
	  c = queryuc_ptr[querypos];
	  g = gsequence_orig[i];
	  g_alt = gsequence_alt[i];
	  if (g == c || g_alt == c) {
	    debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MATCH_COMP,g,querypos,genomepos));
	    unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					   /*cdna*/c,/*comp*/MATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					   /*dynprogindex*/0);
#if 0
	  } else {
	    /* Let stage 3 handle mismatches */
	    debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MISMATCH_COMP,g,querypos,genomepos));
	    unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					   /*cdna*/c,/*comp*/MISMATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					   /*dynprogindex*/0);
#endif
	  }
	  querypos++;
	  genomepos++;
	  i++;
	}
	boundpos = segment->querypos3 + index1part;
	max_genomepos = genomepos - 1;
      }

      /* F.  Make stage2pairs (right) */
      sorted = &(sorted_allocated[anchork+1]);
      boundpos = anchor_segment->querypos5;
      for (k = best_endi; k >= 0; k = prev_right[k]) {
	segment = sorted[k];
	debug13(printf("minus right diagonal %u (%llu), querypos %d..%d, usedp %d, pairablep %d\n",
		       (Chrpos_T) (segment->diagonal - chroffset),(unsigned long long) segment->diagonal,
		       segment->querypos5,segment->querypos3,segment->usedp,segment->pairablep));
	querypos = segment->querypos5;
	if (querypos < boundpos) {
	  /* left = segment->diagonal - querylength; -- FORMULA */
	  genomepos = chrhigh - (segment->diagonal - 1) + querypos;
	  if (genomepos < min_genomepos) {
	    seglength = (segment->querypos3 + index1part) - querypos;
	    Genome_get_segment_blocks_right(gsequence_orig,gsequence_alt,/*left*/segment->diagonal - querypos - seglength,
					    seglength,chrhigh,/*revcomp*/true);

	    i = 0;
	    while (i < seglength && querypos < boundpos && genomepos < min_genomepos) {
	      c = queryuc_ptr[querypos];
	      g = gsequence_orig[i];
	      g_alt = gsequence_alt[i];
	      if (g == c || g_alt == c) {
		debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MATCH_COMP,g,querypos,genomepos));
		unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					       /*cdna*/c,/*comp*/MATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					       /*dynprogindex*/0);
#if 0
	      } else {
		/* Let stage 3 handle mismatches */
		debug13(printf("Pushing %c %c %c at %d,%u\n",c,/*comp*/MISMATCH_COMP,g,querypos,genomepos));
		unsorted_pairs = Pairpool_push(unsorted_pairs,pairpool,querypos,genomepos,
					       /*cdna*/c,/*comp*/MISMATCH_COMP,/*genome*/g,/*genomealt*/g_alt,
					       /*dynprogindex*/0);
#endif
	      }
	      querypos++;
	      genomepos++;
	      i++;
	    }
	    boundpos = segment->querypos5;
	    min_genomepos = chrhigh - (segment->diagonal - 1) + segment->querypos5;
	  }
	}
      }

      /* Sort pairs and get unique ones */
      array = (Pair_T *) List_to_array_n(&npairs,unsorted_pairs);
      qsort(array,npairs,sizeof(Pair_T),Pair_cmp);

      stage2pairs = (List_T) NULL;
      i = 0;
      while (i < npairs) {
	j = i + 1;
	while (j < npairs && array[j]->querypos == array[i]->querypos) {
	  j++;
	}
	if (j == i + 1) {
	  /* Only a single pair at this querypos */
	  debug13(Pair_dump_one(array[i],true));
	  debug13(printf("\n"));
	  stage2pairs = Pairpool_push_existing(stage2pairs,pairpool,array[i]);
	}
	i = j;
      }
      stage2pairs = List_reverse(stage2pairs);
      FREE(array);
	

      /* Run GMAP */
      if (stage2pairs == NULL) {
	/* hit = (T) NULL; */
	  
      } else if ((pairarray1 = Stage3_compute(&cdna_direction,&sensedir,&pairs1,&npairs1,&goodness1,
					      &matches1,&nmatches_posttrim_1,&max_match_length_1,
					      &ambig_end_length_5_1,&ambig_end_length_3_1,
					      &ambig_splicetype_5_1,&ambig_splicetype_3_1,
					      &ambig_prob_5_1,&ambig_prob_3_1,
					      &unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
					      &ncanonical1,&nsemicanonical1,&nnoncanonical1,&min_splice_prob_1,

					      &pairarray2,&pairs2,&npairs2,&goodness2,
					      &matches2,&nmatches_posttrim_2,&max_match_length_2,
					      &ambig_end_length_5_2,&ambig_end_length_3_2,
					      &ambig_splicetype_5_2,&ambig_splicetype_3_2,
					      &ambig_prob_5_2,&ambig_prob_3_2,
					      &unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
					      &ncanonical2,&nsemicanonical2,&nnoncanonical2,&min_splice_prob_2,

					      stage2pairs,/*all_stage2_starts*/NULL,/*all_stage2_ends*/NULL,
#ifdef END_KNOWNSPLICING_SHORTCUT
					      cutoff_level,/*queryptr*/watsonp ? queryuc_ptr : queryrc,
					      watsonp ? query_compress_fwd : query_compress_rev,
#endif
					      /*queryseq_ptr*/queryuc_ptr,queryuc_ptr,querylength,/*skiplength*/0,
#ifdef EXTRACT_GENOMICSEG
					      /*query_subseq_offset*/0,
#else
					      /*query_subseq_offset*/0,
#endif
					      chrnum,chroffset,chrhigh,
					      knownsplice_limit_low,knownsplice_limit_high,/*plusp*/false,genestrand,
					      /*jump_late_p*/true,maxpeelback,pairpool,dynprogL,dynprogM,dynprogR,
					      /*sense_try*/0,/*sense_filter*/0,
					      oligoindices_minor,diagpool,cellpool)) == NULL) {
	/* hit = (T) NULL; */
	  
      } else if (cdna_direction == 0) {
	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray1,npairs1);
	start = add_bounded(chroffset + Pair_genomepos(&(pairarray1[0])),
			    /*plusterm*/Pair_querypos(&(pairarray1[0])),chrhigh);
	end = subtract_bounded(chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			       /*minusterm*/querylength - 1 - Pair_querypos(&(pairarray1[npairs1-1])),chroffset);

	if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_1,max_match_length_1,
				      ambig_end_length_5_1,ambig_end_length_3_1,
				      ambig_splicetype_5_1,ambig_splicetype_3_1,
				      min_splice_prob_1,
				      pairarray1,npairs1,nsegments,nintrons,nindelbreaks,
				      /*left*/end,/*genomiclength*/start - end + 1,
				      /*plusp*/false,genestrand,
				      /*accession*/NULL,querylength,chrnum,chroffset,chrhigh,chrlength,
				      /*cdna_direction*/+1,/*sensedir*/SENSE_FORWARD,
				      /*gmap_source*/GMAP_VIA_SEGMENTS)) == NULL) {
	  FREE_OUT(pairarray1);
	} else {
	  hits = List_push(hits,(void *) hit);
	}

	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray2,npairs2);
	start = add_bounded(chroffset + Pair_genomepos(&(pairarray2[0])),
			    /*plusterm*/Pair_querypos(&(pairarray2[0])),chrhigh);
	end = subtract_bounded(chroffset + Pair_genomepos(&(pairarray2[npairs2-1])),
			       /*minusterm*/querylength - 1 - Pair_querypos(&(pairarray2[npairs2-1])),chroffset);

	if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_2,max_match_length_2,
				      ambig_end_length_5_2,ambig_end_length_3_2,
				      ambig_splicetype_5_2,ambig_splicetype_3_2,
				      min_splice_prob_2,
				      pairarray2,npairs2,nsegments,nintrons,nindelbreaks,
				      /*left*/end,/*genomiclength*/start - end + 1,
				      /*plusp*/false,genestrand,
				      /*accession*/NULL,querylength,chrnum,chroffset,chrhigh,chrlength,
				      /*cdna_direction*/-1,/*sensedir*/SENSE_ANTI,
				      /*gmap_source*/GMAP_VIA_SEGMENTS)) == NULL) {
	  FREE_OUT(pairarray2);
	} else {
	  hits = List_push(hits,(void *) hit);
	}

      } else {
	nsegments = Pair_gsnap_nsegments(&nmismatches_whole,&nindels,&nintrons,&nindelbreaks,
					 pairarray1,npairs1);
	start = add_bounded(chroffset + Pair_genomepos(&(pairarray1[0])),
			    /*plusterm*/Pair_querypos(&(pairarray1[0])),chrhigh);
	end = subtract_bounded(chroffset + Pair_genomepos(&(pairarray1[npairs1-1])),
			       /*minusterm*/querylength - 1 - Pair_querypos(&(pairarray1[npairs1-1])),chroffset);

	if ((hit = Stage3end_new_gmap(nmismatches_whole,nmatches_posttrim_1,max_match_length_1,
				      ambig_end_length_5_1,ambig_end_length_3_1,
				      ambig_splicetype_5_1,ambig_splicetype_3_1,
				      min_splice_prob_1,
				      pairarray1,npairs1,nsegments,nintrons,nindelbreaks,
				      /*left*/end,/*genomiclength*/start - end + 1,
				      /*plusp*/false,genestrand,
				      /*accession*/NULL,querylength,chrnum,chroffset,chrhigh,chrlength,
				      cdna_direction,sensedir,/*gmap_source*/GMAP_VIA_SEGMENTS)) == NULL) {
	  FREE_OUT(pairarray1);
	} else {
	  hits = List_push(hits,(void *) hit);
	}
      }
    }
  }

  if (minus_nsegments > 0) {
#ifdef HAVE_ALLOCA
    FREEA(gsequence_alt);
    FREEA(gsequence_orig);
    if (minus_nsegments < MAX_ALLOCATION) {
      FREEA(sorted_allocated);
      FREEA(scores_allocated);
      FREEA(prev_allocated);
    } else {
      FREE(sorted_allocated);
      FREE(scores_allocated);
      FREE(prev_allocated);
    }
#else
    FREE(gsequence_alt);
    FREE(gsequence_orig);
    FREE(sorted_allocated);
    FREE(scores_allocated);
    FREE(prev_allocated);
#endif
  }

  return hits;
}


static List_T
align_singleend_with_gmap (List_T result, char *queryuc_ptr, int querylength,
			   Oligoindex_array_T oligoindices_minor,
			   Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			   int user_maxlevel) {
  List_T new_result = NULL;
  Stage3end_T hit, gmap1, gmap2;
  List_T p;
  int genestrand;
  int i;
#ifdef DEBUG13
  int missing_hit, missing_gmap;
#endif
  
  
  debug13(printf("Sorting hits by nmatches\n"));
  result = Stage3end_sort_bymatches(result);
  
  for (p = result, i = 0; p != NULL && i < max_gmap_improvement; p = p->rest, i++) {
    hit = (Stage3end_T) List_head(p);
    genestrand = Stage3end_genestrand(hit);
    
    debug13(printf("GMAP improvement: Entering align_singleend_with_gmap with hittype %s\n",
		   Stage3end_hittype_string(hit)));
    
    /* Was querylength5 - Stage3end_matches(hit5) > 5 */
    if (Stage3end_hittype(hit) == GMAP) {
      /* Skip */
      debug13(printf("Skipping hit of type GMAP\n"));
      new_result = List_push(new_result,(void *) hit);
      
    } else if (Stage3end_improved_by_gmap_p(hit) == true) {
      /* Skip */
      debug13(printf("Skipping hit already improved by GMAP\n"));
      new_result = List_push(new_result,(void *) hit);
      
#if 0
      /* Don't skip on final align_singleend_with_gmap */
    } else if (Stage3end_hittype(hit) == TERMINAL) {
      /* Skip */
      debug13(printf("Skipping hit of type TERMINAL\n"));
      new_result = List_push(new_result,(void *) hit);
#endif
      
    } else if (querylength - Stage3end_nmatches_posttrim(hit) <= user_maxlevel) {
      /* Skip */
      debug13(printf("Skipping hit with nmismatches %d - %d <= user_maxlevel %d\n",
		     querylength,Stage3end_nmatches_posttrim(hit),user_maxlevel));
      new_result = List_push(new_result,(void *) hit);
      
    } else if (Stage3end_terminal_trim(hit) <= GMAP_TERMINAL_TRIM
	       && Stage3end_contains_known_splicesite(hit) == false
	       ) {
      debug13(printf("Skipping good hit\n"));
      new_result = List_push(new_result,(void *) hit);
      
    } else {
      debug13(printf("To correct hit terminalp %d or known_splicesite %d, running GMAP on 5' to match with 3' end\n",
		     Stage3end_hittype(hit) == TERMINAL,
		     Stage3end_contains_known_splicesite(hit)));
      
      /* Want high quality because we already have a pretty good answer */
      align_single_hit_with_gmap(&gmap1,&gmap2,hit,queryuc_ptr,querylength,
#ifdef END_KNOWNSPLICING_SHORTCUT
				 queryrc,Shortread_invertedp(queryseq),
#endif
				 oligoindices_minor,
				 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				 genestrand,/*extend_ends_p*/false);
      if (gmap1 == NULL) {
	new_result = List_push(new_result,(void *) hit);
      } else {
	debug13(missing_hit = querylength - Stage3end_nmatches_posttrim(hit));
	debug13(missing_gmap = querylength - Stage3end_nmatches_posttrim(gmap1));
	debug13(printf("GMAP %p with %d matches, %d missing compared with original hit with %d matches, %d missing\n",
		       gmap1,Stage3end_nmatches_posttrim(gmap1),missing_gmap,Stage3end_nmatches_posttrim(hit),missing_hit));
	new_result = List_push(new_result,(void *) gmap1);

	if (gmap2 != NULL) {
	  new_result = List_push(new_result,(void *) gmap2);
	}
	Stage3end_set_improved_by_gmap(hit);
      }
    }
  }
  
  for ( ; p != NULL; p = p->rest) {
    hit = (Stage3end_T) List_head(p);
    new_result = List_push(new_result,(void *) hit);
  }
  
  List_free(&result);
  return new_result;
}


/* Search order for single-end reads:

   1.  suffix array
   2.  exact/subs, via spanning set algorithm
   3.  subs/indels, via complete set algorithm
   4.  segments -> single splicing
   5.  segments -> double splicing (currently disabled)

   (6).  paired segments -> GMAP via segments (not applicable for single-end reads)
   7.  distant splicing (needs to be before terminals, or we won't find them)
   8.  terminals
   
   (9).  if still no concordance: GMAP pairsearch (not applicable for single-end reads)
   9.  if found score is low: anchor segments -> GMAP via segments

   10.  GMAP improvement
*/


/* done_level should probably be renamed final_level.  opt_level
   should probably be renamed found_level or opt_level. */
static List_T
align_end (int *cutoff_level, T this,
	   Compress_T query_compress_fwd, Compress_T query_compress_rev,
	   char *queryuc_ptr, char *queryrc, int querylength, int query_lastpos,
	   Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev, int indexdb_size_threshold, Floors_T *floors_array,
	   
	   Oligoindex_array_T oligoindices_minor,
	   Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
	   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	   
	   int user_maxlevel, int indel_penalty_middle, int indel_penalty_end,
	   int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
	   bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
	   bool keep_floors_p, int genestrand, bool first_read_p) {
  List_T hits, greedy = NULL, subs = NULL, terminals = NULL, indels = NULL,
    singlesplicing = NULL, doublesplicing = NULL, shortendsplicing = NULL,
    longsinglesplicing = NULL, distantsplicing = NULL, gmap_hits = NULL, q;
  Segment_T *plus_anchor_segments = NULL, *minus_anchor_segments = NULL;
  int n_plus_anchors = 0, n_minus_anchors = 0;
  Stage3end_T hit, gmap1, gmap2;
  int found_score, done_level, opt_level, fast_level, mismatch_level, nmismatches;
  int max_splice_mismatches, i;
  int nhits = 0, nsplicepairs = 0;
  List_T *startfrags_plus, *endfrags_plus, *startfrags_minus, *endfrags_minus;
  List_T *donors_plus, *antidonors_plus, *acceptors_plus, *antiacceptors_plus,
    *donors_minus, *antidonors_minus, *acceptors_minus, *antiacceptors_minus;
  bool any_omitted_p, ambiguousp, alloc_floors_p = false, floors_computed_p = false;
  Floors_T floors;
  bool spanningsetp, completesetp, gmapp;
  bool segments_computed_p = false;
  Indexdb_T plus_indexdb, minus_indexdb;
  bool allvalidp;
  int min_trim, trim;
#ifndef LARGE_GENOMES
  int nmisses_allowed_sarray;
#endif
#ifdef DEBUG13
  int missing_hit, missing_gmap;
#endif

  
  if (genestrand == +2) {
    plus_indexdb = indexdb_rev;
    minus_indexdb = indexdb_fwd;
  } else {
    plus_indexdb = indexdb_fwd;
    minus_indexdb = indexdb_rev;
  }
  
  found_score = querylength;
  if (querylength < min_kmer_readlength) {
    fast_level = querylength - 1 - NREQUIRED_FAST;
    debug(printf("fast_level %d = querylength %d - 1 - nrequired_fast %d\n",
		 fast_level,querylength,NREQUIRED_FAST));
  } else {
    fast_level = (querylength + index1interval - 1)/spansize - NREQUIRED_FAST;
    debug(printf("fast_level %d = (querylength %d + index1interval %d - 1)/spansize %d - nrequired_fast %d\n",
		 fast_level,querylength,index1interval,spansize,NREQUIRED_FAST));
  }
  
#if 0
  /* This prevents complete_mm procedure, needed for short reads */
  if (fast_level < 1 && user_maxlevel < 0) {
    debug(printf("Changing fast_level to 0\n"));
    fast_level = 1;		/* Do at least 1 mismatch */
  }
#endif
  
  if (user_maxlevel >= 0) {
    *cutoff_level = user_maxlevel;
  } else if (fast_level >= 0) {
    *cutoff_level = fast_level;
  } else {
    *cutoff_level = 0;
  }
  debug(printf("cutoff_level = %d\n",*cutoff_level));
  
  if (user_maxlevel < 0) {
    if (fast_level >= 0) {
      user_maxlevel = fast_level;
    } else {
      user_maxlevel = 0;
    }
  }
  debug(printf("user_maxlevel = %d\n",user_maxlevel));
  
#if 0
  if (dibasep) {
    opt_level = querylength;	/* Allow extra because color errors may exceed nt errors */
  }
#endif
  opt_level = user_maxlevel;
  done_level = user_maxlevel /* + subopt_levels.  -- Initially the same */;
  debug(printf("0> opt_level %d, done_level %d\n",opt_level,done_level));
  
  nhits = 0;
  
#ifndef LARGE_GENOMES
  nmisses_allowed_sarray = *cutoff_level;
  
  if (use_only_sarray_p == true || (use_sarray_p == true && querylength < min_kmer_readlength)) {
    hits = Sarray_search_greedy(&(*cutoff_level),
				queryuc_ptr,queryrc,querylength,query_compress_fwd,query_compress_rev,
				nmisses_allowed_sarray,genestrand);
    
    hits = Stage3end_remove_overlaps(hits,/*finalp*/true);
    hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				   querylength,/*keep_gmap_p*/false,/*finalp*/true);
    hits = Stage3end_resolve_multimapping(hits);
    
    hits = Stage3end_remove_circular_alias(hits);
    hits = Stage3end_remove_duplicates(hits); /* Aliases can cause duplicates */
    
    return hits;
  }
#endif
  

  /* Search 1: Suffix array */
#ifdef LARGE_GENOMES
  spanningsetp = true;
#else
  if (use_sarray_p == false) {
    spanningsetp = true;
  } else {
    spanningsetp = false;	/* Suffix array search replaces spanning set */
    
    debug(printf("Trying suffix array\n"));
    greedy = Sarray_search_greedy(&found_score,queryuc_ptr,queryrc,querylength,query_compress_fwd,query_compress_rev,
				  nmisses_allowed_sarray,genestrand);
    
    opt_level = (found_score < opt_level) ? found_score : opt_level;
    if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
      done_level = user_maxlevel;
    }
    debug(printf("SA> found_score %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
    
    for (q = greedy; q != NULL; q = List_next(q)) {
      hit = (Stage3end_T) List_head(q);
      if (Stage3end_total_trim(hit) > 15) {
	spanningsetp = true;
      }
    }
    debug(printf("spanningsetp %d\n",spanningsetp));
  }
#endif
  
  if (querylength < min_kmer_readlength) {
    spanningsetp = false;
  }

  /* Search 2: Exact/subs via spanning set */

  if (spanningsetp == true) {
    read_oligos(&allvalidp,this,queryuc_ptr,querylength,query_lastpos,/*genestrand*/0);

    /* 1. Exact.  Requires compress if cmet or genomealt.  Creates and uses spanning set. */
    mismatch_level = 0;
    if (done_level == 0 && snpp == false) {
      debug(printf("Suffix array already found exact matches and no SNPs, so spanning set can't do any better\n"));
    } else if (allvalidp == false) {
      debug(printf("Not all oligos are valid, so cannot perform spanning set\n"));
      fast_level = -1;
      spanningsetp = false;
    } else {
      debug(printf("fast_level = %d\n",fast_level));
      debug(printf("*** Stage 1.  Exact ***\n"));
      subs = find_spanning_exact_matches(&found_score,&nhits,/*hits*/NULL,this,genestrand,
					 querylength,query_lastpos,plus_indexdb,minus_indexdb,
					 query_compress_fwd,query_compress_rev);
      opt_level = (found_score < opt_level) ? found_score : opt_level;
      if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	done_level = user_maxlevel;
      }
      mismatch_level = 1;
      debug(printf("1> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
    }

    /* 2. One mismatch.  Requires spanning set and compress. */
    if (allvalidp && querylength >= one_miss_querylength && done_level >= 1) {
      debug(printf("*** Stage 2.  One miss ***\n"));
      subs = find_spanning_onemiss_matches(&found_score,&nhits,subs,this,genestrand,
					   querylength,query_compress_fwd,query_compress_rev);
      opt_level = (found_score < opt_level) ? found_score : opt_level;
      if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	done_level = user_maxlevel;
      }
      mismatch_level = 2;
      debug(printf("2> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
    }

    /* 3. Mismatches via spanning set.  Requires spanning set and compress. */
    if (allvalidp && done_level >= 2) {
      while (mismatch_level <= fast_level && mismatch_level <= done_level) {
	debug(printf("*** Stage 3 (level %d).  Spanning set mismatches ***\n",mismatch_level));
	subs = find_spanning_multimiss_matches(&found_score,&nhits,subs,this,genestrand,
					       NREQUIRED_FAST,querylength,query_compress_fwd,query_compress_rev,
					       /*nmisses_allowed*/mismatch_level);
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
	mismatch_level++;
	debug(printf("3> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
      }
    }
  }


  /* Search 3: Subs/indels via complete set */

  /* 4, 5.  Complete set mismatches and indels, omitting frequent oligos */
  if (subs == NULL) {
    completesetp = true;
  } else {
    completesetp = false;
    for (q = subs; q != NULL; q = List_next(q)) {
      hit = (Stage3end_T) List_head(q);
      debug(printf("Hit has total score of %d\n",Stage3end_score(hit)));
      if (Stage3end_score(hit) > done_level) {
	completesetp = true;
      }
    }
  }
  debug(printf("completesetp %d\n",completesetp));

#if 0
  if (found_score <= done_level) {
    debug(printf("Test for completeset: false because found_score %d >done_level %d\n",found_score,done_level));
    completesetp = false;
  }
#endif

  if (querylength < min_kmer_readlength) {
    completesetp = false;
  }

  if (completesetp == true) {
    if (this->read_oligos_p == false) {
      read_oligos(&allvalidp,this,queryuc_ptr,querylength,query_lastpos,/*genestrand*/0);
    }

    floors = compute_floors(&any_omitted_p,&alloc_floors_p,floors_array,this,queryuc_ptr,querylength,query_lastpos,
			    plus_indexdb,minus_indexdb,indexdb_size_threshold,max_end_insertions,
			    /*omit_frequent_p*/true,/*omit_repetitive_p*/true,keep_floors_p);
    floors_computed_p = true;
    complete_set_mm_indels(&found_score,&segments_computed_p,
			   &plus_anchor_segments,&minus_anchor_segments,
			   &n_plus_anchors,&n_minus_anchors,
			   &opt_level,&done_level,user_maxlevel,/*revise_levels_p*/true,
			   &nhits,&subs,&indels,this,query_compress_fwd,query_compress_rev,
#if defined(DEBUG2) || defined(DEBUG2E)
			   queryuc_ptr,queryrc,
#endif
			   querylength,query_lastpos,floors,indel_penalty_middle,indel_penalty_end,
			   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			   fast_level,genestrand);
    if (found_score <= done_level) {
      debug(printf("Test for completeset: false because found_score %d >done_level %d\n",found_score,done_level));
      completesetp = false;
    }
  }


  /* Search 4: Segments -> single splicing */

  /* 6/7/8/9.  Splicing.  Requires compress and all positions fetched */
  /* SPEED: For more hits, turn off first branch */
  if (use_sarray_p == true && completesetp == false) {
      /* Skip.  Suffix array already found something.  Also, get memory errors if run both algorithms.  */

  } else if (knownsplicingp || novelsplicingp || find_dna_chimeras_p) {
    /* 6.  Single splicing */
    debug(printf("Deciding whether to do singlesplicing: done_level %d >=? localsplicing_penalty %d\n",
		 done_level,localsplicing_penalty));
    if (done_level >= localsplicing_penalty) {
      debug(printf("*** Stage 6.  Single splicing masking frequent oligos with done_level %d ***\n",done_level));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      if (floors_computed_p == false) {
	floors = compute_floors(&any_omitted_p,&alloc_floors_p,floors_array,this,queryuc_ptr,querylength,query_lastpos,
				plus_indexdb,minus_indexdb,indexdb_size_threshold,max_end_insertions,
				/*omit_frequent_p*/true,/*omit_repetitive_p*/true,keep_floors_p);
	floors_computed_p = true;
      }

      if (segments_computed_p == false) {
	this->plus_segments = identify_all_segments(&this->plus_nsegments,&plus_anchor_segments,&n_plus_anchors,
						    &this->plus_spliceable,&this->plus_nspliceable,
#ifdef LARGE_GENOMES
						    this->plus_positions_high,this->plus_positions_low,
#else
						    this->plus_positions,
#endif
						    this->plus_npositions,this->omitted,querylength,query_lastpos,
						    floors,/*plusp*/true);
	this->minus_segments = identify_all_segments(&this->minus_nsegments,&minus_anchor_segments,&n_minus_anchors,
						     &this->minus_spliceable,&this->minus_nspliceable,
#ifdef LARGE_GENOMES
						     this->minus_positions_high,this->minus_positions_low,
#else
						     this->minus_positions,
#endif
						     this->minus_npositions,this->omitted,querylength,query_lastpos,
						     floors,/*plusp*/false);
	segments_computed_p = true;
      }

      singlesplicing = complete_set_singlesplicing(&found_score,singlesplicing,floors,this,
						   query_compress_fwd,query_compress_rev,
						   querylength,query_lastpos,
						   localsplicing_penalty,
						   /*max_mismatches_allowed*/done_level - localsplicing_penalty,
						   genestrand,first_read_p,
						   /*subs_or_indels_p*/(subs != NULL || indels != NULL) ? true : false);

#if 0
      /* Mark ambiguous splices only for single-end reads */
      singlesplicing = Stage3end_mark_ambiguous_splices(&ambiguousp,singlesplicing);
#endif
      singlesplicing = Stage3end_optimal_score(singlesplicing,query_compress_fwd,query_compress_rev,querylength,
					       /*keep_gmap_p*/true,/*finalp*/false);

      if (singlesplicing) {
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
      }
    }


    /* Search 5: Segments -> double splicing (currently disabled) */
#ifdef PERFORM_DOUBLESPLICING
    /* 7.  Double splicing */
    debug(printf("Deciding whether to do doublesplicing: done_level %d >=? localsplicing_penalty %d\n",
		 done_level,localsplicing_penalty));
    if (done_level >= localsplicing_penalty) {
      debug(printf("*** Stage 7.  Double splicing masking frequent oligos with done_level %d ***\n",done_level));
      if (floors_computed_p == false) {
	floors = compute_floors(&any_omitted_p,&alloc_floors_p,floors_array,this,queryuc_ptr,querylength,query_lastpos,
				plus_indexdb,minus_indexdb,indexdb_size_threshold,max_end_insertions,
				/*omit_frequent_p*/true,/*omit_repetitive_p*/true,keep_floors_p);
	floors_computed_p = true;
      }
      doublesplicing = complete_set_doublesplicing(&found_score,doublesplicing,floors,this,
						   query_compress_fwd,query_compress_rev,
						   queryuc_ptr,queryrc,querylength,query_lastpos,
						   localsplicing_penalty,min_shortend,
						   /*max_mismatches_allowed*/done_level - localsplicing_penalty,
						   /*pairedp*/false,genestrand,first_read_p,
						   /*subs_or_indels_p*/(subs != NULL || indels != NULL) ? true : false);
      
#if 0
      /* Mark ambiguous splices only for single-end reads */
      doublesplicing = Stage3end_mark_ambiguous_splices(&ambiguousp,doublesplicing);
#endif
      doublesplicing = Stage3end_optimal_score(doublesplicing,query_compress_fwd,query_compress_rev,querylength,
					       /*keep_gmap_p*/true,/*finalp*/false);

      if (doublesplicing) {
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
      }
    }
#endif

    if (knownsplicingp == true && done_level >= localsplicing_penalty &&
	(max_splice_mismatches = done_level - localsplicing_penalty) >= 0) {
      /* Want >= and not > to give better results.  Negligible effect on speed. */
      /* 8.  Shortend splicing */
      debug(printf("*** Stage 8.  Short-end splicing, allowing %d mismatches ***\n",max_splice_mismatches));

      donors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      antidonors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      acceptors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      antiacceptors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      donors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      antidonors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      acceptors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      antiacceptors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));

      if (floors_computed_p == false) {
	floors = compute_floors(&any_omitted_p,&alloc_floors_p,floors_array,this,queryuc_ptr,querylength,query_lastpos,
				plus_indexdb,minus_indexdb,indexdb_size_threshold,max_end_insertions,
				/*omit_frequent_p*/true,/*omit_repetitive_p*/true,keep_floors_p);
	floors_computed_p = true;
      }

      debug(printf("Starting find_spliceends (plus) with %d anchor segments\n",n_plus_anchors));
      find_spliceends_shortend(&donors_plus,&antidonors_plus,&acceptors_plus,&antiacceptors_plus,
			       plus_anchor_segments,n_plus_anchors,
#ifdef DEBUG4E
			       /*queryptr*/queryuc_ptr,
#endif
			       floors,querylength,query_lastpos,/*query_compress*/query_compress_fwd,
			       max_splice_mismatches,/*plusp*/true,genestrand);
      debug(printf("Finished find_spliceends (plus)\n"));

      debug(printf("Starting find_spliceends (minus) with %d anchor segments\n",n_minus_anchors));
      find_spliceends_shortend(&antidonors_minus,&donors_minus,&antiacceptors_minus,&acceptors_minus,
			       minus_anchor_segments,n_minus_anchors,
#ifdef DEBUG4E
			       /*queryptr*/queryrc,
#endif
			       floors,querylength,query_lastpos,/*query_compress*/query_compress_rev,
			       max_splice_mismatches,/*plusp*/false,genestrand);
      debug(printf("Finished find_spliceends (minus)\n"));

      shortendsplicing = find_splicepairs_shortend(&found_score,shortendsplicing,
						   donors_plus,antidonors_plus,
						   acceptors_plus,antiacceptors_plus,
						   donors_minus,antidonors_minus,
						   acceptors_minus,antiacceptors_minus,
						   query_compress_fwd,query_compress_rev,
						   queryuc_ptr,queryrc,min_shortend,
						   localsplicing_penalty,
						   /*max_mismatches_allowed*/max_splice_mismatches,querylength,
						   /*pairedp*/false,genestrand,first_read_p);
      opt_level = (found_score < opt_level) ? found_score : opt_level;
      if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	done_level = user_maxlevel;
      }
      debug(printf("8> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));

      for (i = 0; i <= max_splice_mismatches; i++) {
	substringlist_gc(&(donors_plus[i]));
	substringlist_gc(&(antidonors_plus[i]));
	substringlist_gc(&(acceptors_plus[i]));
	substringlist_gc(&(antiacceptors_plus[i]));
	substringlist_gc(&(donors_minus[i]));
	substringlist_gc(&(antidonors_minus[i]));
	substringlist_gc(&(acceptors_minus[i]));
	substringlist_gc(&(antiacceptors_minus[i]));
      }
      FREEA(donors_plus);
      FREEA(antidonors_plus);
      FREEA(acceptors_plus);
      FREEA(antiacceptors_plus);
      FREEA(donors_minus);
      FREEA(antidonors_minus);
      FREEA(acceptors_minus);
      FREEA(antiacceptors_minus);
    }
  }

  debug(printf("Before distant splicing:\n"));
  debug(printf("  greedy: %d\n",List_length(greedy)));
  debug(printf("  subs: %d\n",List_length(subs)));
  debug(printf("  indels: %d\n",List_length(indels)));
  debug(printf("  singlesplicing %d\n",List_length(singlesplicing)));
  debug(printf("  doublesplicing %d\n",List_length(doublesplicing)));
  debug(printf("  shortendsplicing: %d\n",List_length(shortendsplicing)));
  debug(printf("  done_level: %d\n",done_level));

  hits = List_append(greedy,
		     List_append(subs,
				 List_append(terminals,
					     List_append(indels,
							 List_append(singlesplicing,shortendsplicing)))));

  if (knownsplicingp || novelsplicingp || find_dna_chimeras_p) {
    /* Search 7: Distant splicing */
    min_trim = querylength;
    for (q = hits; q != NULL; q = q->rest) {
      hit = (Stage3end_T) q->first;
      if ((trim = Stage3end_total_trim(hit)) < min_trim) {
	min_trim = trim;
      }
    }

    if (done_level < distantsplicing_penalty && min_trim < min_distantsplicing_end_matches) {
      /* Want < and not <=, because otherwise distant splicing does not work on 50-bp reads */
      /* Want <= and not <, because distant splicing needs to be better than other alternatives */
      /* Don't find distant splicing */
      debug(printf("Skipping distant splicing because done_level %d < distantsplicing_penalty %d and min_trim %d < %d\n",
		   done_level,distantsplicing_penalty,min_trim,min_distantsplicing_end_matches));

    } else if (find_dna_chimeras_p == true &&
	       (max_splice_mismatches = done_level - distantsplicing_penalty) >= 0) {
      /* 9 (DNA).  Find distant splicing for DNA */
      debug(printf("*** Stage 9 (DNA).  Distant splice ends, allowing %d mismatches ***\n",max_splice_mismatches));

      startfrags_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      endfrags_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      startfrags_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      endfrags_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));

      if (floors_computed_p == false) {
	floors = compute_floors(&any_omitted_p,&alloc_floors_p,floors_array,this,queryuc_ptr,querylength,query_lastpos,
				plus_indexdb,minus_indexdb,indexdb_size_threshold,max_end_insertions,
				/*omit_frequent_p*/true,/*omit_repetitive_p*/true,keep_floors_p);
	floors_computed_p = true;
      }

      debug(printf("Starting find_spliceends_distant_dna_plus\n"));
      find_spliceends_distant_dna_plus(&startfrags_plus,&endfrags_plus,plus_anchor_segments,n_plus_anchors,
#ifdef DEBUG4E
				       /*queryptr*/queryuc_ptr,
#endif
				       floors,querylength,query_lastpos,/*query_compress*/query_compress_fwd,
				       max_splice_mismatches,genestrand);
      debug(printf("Finished find_spliceends_distant_dna_plus\n"));

      debug(printf("Starting find_spliceends_distant_dna_minus\n"));
      find_spliceends_distant_dna_minus(&startfrags_minus,&endfrags_minus,minus_anchor_segments,n_minus_anchors,
#ifdef DEBUG4E
					/*queryptr*/queryrc,
#endif
					floors,querylength,query_lastpos,/*query_compress*/query_compress_rev,
					max_splice_mismatches,genestrand);
      debug(printf("Finished find_spliceends_distant_dna_minus\n"));

      nmismatches = 0;
      ambiguousp = false;
      while (longsinglesplicing == NULL &&
	     nmismatches <= done_level - distantsplicing_penalty &&
	     nsplicepairs < MAXCHIMERAPATHS && ambiguousp == false) {
	debug(printf("*** Stage 9 (DNA).  Distant splicing, allowing %d mismatches ***\n",nmismatches));

	debug4e(printf("Sorting splice ends\n"));
	startfrags_plus[nmismatches] = Substring_sort_siteN_halves(startfrags_plus[nmismatches],/*ascendingp*/true);
	endfrags_plus[nmismatches] = Substring_sort_siteN_halves(endfrags_plus[nmismatches],/*ascendingp*/true);

	startfrags_minus[nmismatches] = Substring_sort_siteN_halves(startfrags_minus[nmismatches],/*ascendingp*/false);
	endfrags_minus[nmismatches] = Substring_sort_siteN_halves(endfrags_minus[nmismatches],/*ascendingp*/false);

	debug4e(printf("Splice ends at %d nmismatches: +startfrags/endfrags %d/%d, -startfrags/endfrags %d/%d\n",
		       nmismatches,
		       List_length(startfrags_plus[nmismatches]),List_length(endfrags_plus[nmismatches]),
		       List_length(startfrags_minus[nmismatches]),List_length(endfrags_minus[nmismatches])));

	distantsplicing = find_splicepairs_distant_dna(&found_score,&nsplicepairs,&longsinglesplicing,distantsplicing,
						       startfrags_plus,endfrags_plus,startfrags_minus,endfrags_minus,
						       localsplicing_penalty,distantsplicing_penalty,
						       querylength,nmismatches,first_read_p);
#if 0
	assert(List_length(distantsplicing) <= 1);
#endif

#if 0
	/* Mark ambiguous splices only for single-end reads */
	distantsplicing = Stage3end_mark_ambiguous_splices(&ambiguousp,distantsplicing);
#endif

	/* Excess distant splicing should be freed already in find_splicepairs_distant_rna */
	debug(printf("Entering Stage3end_optimal_score with %d hits\n",List_length(distantsplicing)));
	distantsplicing = Stage3end_optimal_score(distantsplicing,query_compress_fwd,query_compress_rev,querylength,
						  /*keep_gmap_p*/true,/*finalp*/false);
	debug(printf("Exiting Stage3end_optimal_score with %d hits\n",List_length(distantsplicing)));

	if (distantsplicing) {
	  opt_level = (found_score < opt_level) ? found_score : opt_level;
	  if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	    done_level = user_maxlevel;
	  }
	  debug(printf("9 (DNA)> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
	}
	nmismatches++;

      }

      if (longsinglesplicing != NULL) {
	debug(printf("Entering Stage3end_optimal_score with %d longsinglesplicing hits\n",List_length(longsinglesplicing)));
	longsinglesplicing = Stage3end_optimal_score(longsinglesplicing,query_compress_fwd,query_compress_rev,querylength,
						     /*keep_gmap_p*/true,/*finalp*/false);
	debug(printf("Exiting Stage3end_optimal_score with %d hits\n",List_length(longsinglesplicing)));

	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
	debug(printf("9 (DNA)> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
      }

      for (i = 0; i <= max_splice_mismatches; i++) {
	substringlist_gc(&(startfrags_plus[i]));
	substringlist_gc(&(endfrags_plus[i]));
	substringlist_gc(&(startfrags_minus[i]));
	substringlist_gc(&(endfrags_minus[i]));
      }
      FREEA(startfrags_plus);
      FREEA(endfrags_plus);
      FREEA(startfrags_minus);
      FREEA(endfrags_minus);

    } else if ((knownsplicingp || novelsplicingp) &&
	       (max_splice_mismatches = done_level - distantsplicing_penalty) >= 0) {
      /* 9 (RNA).  Find distant splicing for RNA iteratively using both known and novel splice sites */
      debug(printf("*** Stage 9 (RNA).  Distant splice ends, allowing %d mismatches ***\n",max_splice_mismatches));

      donors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      antidonors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      acceptors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      antiacceptors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      donors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      antidonors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      acceptors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
      antiacceptors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));

      if (floors_computed_p == false) {
	floors = compute_floors(&any_omitted_p,&alloc_floors_p,floors_array,this,queryuc_ptr,querylength,query_lastpos,
				plus_indexdb,minus_indexdb,indexdb_size_threshold,max_end_insertions,
				/*omit_frequent_p*/true,/*omit_repetitive_p*/true,keep_floors_p);
	floors_computed_p = true;
      }

      debug(printf("Starting find_spliceends_distant_rna (plus)\n"));
      find_spliceends_distant_rna(&donors_plus,&antidonors_plus,&acceptors_plus,&antiacceptors_plus,
				  plus_anchor_segments,n_plus_anchors,
#ifdef DEBUG4E
				  /*queryptr*/queryuc_ptr,
#endif
				  floors,querylength,query_lastpos,/*query_compress*/query_compress_fwd,
				  max_splice_mismatches,/*plusp*/true,genestrand);
      debug(printf("Finished find_spliceends_distant_rna (plus)\n"));


      debug(printf("Starting find_spliceends_distant_rna (minus)\n"));
      find_spliceends_distant_rna(&antidonors_minus,&donors_minus,&antiacceptors_minus,&acceptors_minus,
				  minus_anchor_segments,n_minus_anchors,
#ifdef DEBUG4E
				  /*queryptr*/queryrc,
#endif
				  floors,querylength,query_lastpos,/*query_compress*/query_compress_rev,
				  max_splice_mismatches,/*plusp*/false,genestrand);
      debug(printf("Finished find_spliceends_distant_rna (minus)\n"));


      nmismatches = 0;
      ambiguousp = false;
      while (longsinglesplicing == NULL &&
	     nmismatches <= done_level - distantsplicing_penalty &&
	     nsplicepairs < MAXCHIMERAPATHS && ambiguousp == false) {
	debug(printf("*** Stage 9 (RNA).  Distant splicing, allowing %d mismatches ***\n",nmismatches));

	debug4e(printf("Sorting splice ends\n"));
	donors_plus[nmismatches] = Substring_sort_siteD_halves(donors_plus[nmismatches],/*ascendingp*/true);
	acceptors_plus[nmismatches] = Substring_sort_siteA_halves(acceptors_plus[nmismatches],/*ascendingp*/true);

	antidonors_plus[nmismatches] = Substring_sort_siteD_halves(antidonors_plus[nmismatches],/*ascendingp*/false);
	antiacceptors_plus[nmismatches] = Substring_sort_siteA_halves(antiacceptors_plus[nmismatches],/*ascendingp*/false);

	donors_minus[nmismatches] = Substring_sort_siteD_halves(donors_minus[nmismatches],/*ascendingp*/false);
	acceptors_minus[nmismatches] = Substring_sort_siteA_halves(acceptors_minus[nmismatches],/*ascendingp*/false);

	antidonors_minus[nmismatches] = Substring_sort_siteD_halves(antidonors_minus[nmismatches],/*ascendingp*/true);
	antiacceptors_minus[nmismatches] = Substring_sort_siteA_halves(antiacceptors_minus[nmismatches],/*ascendingp*/true);

	debug4e(printf("Splice ends at %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		       nmismatches,
		       List_length(donors_plus[nmismatches]),List_length(acceptors_plus[nmismatches]),
		       List_length(antidonors_plus[nmismatches]),List_length(antiacceptors_plus[nmismatches]),
		       List_length(donors_minus[nmismatches]),List_length(acceptors_minus[nmismatches]),
		       List_length(antidonors_minus[nmismatches]),List_length(antiacceptors_minus[nmismatches])));

	distantsplicing = find_splicepairs_distant_rna(&found_score,&nsplicepairs,&longsinglesplicing,distantsplicing,
						       donors_plus,antidonors_plus,acceptors_plus,antiacceptors_plus,
						       donors_minus,antidonors_minus,acceptors_minus,antiacceptors_minus,
						       localsplicing_penalty,distantsplicing_penalty,
						       querylength,nmismatches,first_read_p);
#if 0
	assert(List_length(distantsplicing) <= 1);
#endif

#if 0
	/* Mark ambiguous splices only for single-end reads */
	distantsplicing = Stage3end_mark_ambiguous_splices(&ambiguousp,distantsplicing);
#endif


	/* Excess distant splicing should be freed already in find_splicepairs_distant_rna */
	debug(printf("Entering Stage3end_optimal_score with %d hits\n",List_length(distantsplicing)));
	distantsplicing = Stage3end_optimal_score(distantsplicing,query_compress_fwd,query_compress_rev,querylength,
						  /*keep_gmap_p*/true,/*finalp*/false);
	debug(printf("Exiting Stage3end_optimal_score with %d hits\n",List_length(distantsplicing)));

	if (distantsplicing) {
	  opt_level = (found_score < opt_level) ? found_score : opt_level;
	  if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	    done_level = user_maxlevel;
	  }
	  debug(printf("9 (RNA)> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
	}
	nmismatches++;

      }

      if (longsinglesplicing != NULL) {
	debug(printf("Entering Stage3end_optimal_score with %d longsinglesplicing hits\n",List_length(longsinglesplicing)));
	longsinglesplicing = Stage3end_optimal_score(longsinglesplicing,query_compress_fwd,query_compress_rev,querylength,
						     /*keep_gmap_p*/true,/*finalp*/false);
	debug(printf("Exiting Stage3end_optimal_score with %d hits\n",List_length(longsinglesplicing)));

	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
	debug(printf("9 (RNA)> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
      }

      for (i = 0; i <= max_splice_mismatches; i++) {
	substringlist_gc(&(donors_plus[i]));
	substringlist_gc(&(antidonors_plus[i]));
	substringlist_gc(&(acceptors_plus[i]));
	substringlist_gc(&(antiacceptors_plus[i]));
	substringlist_gc(&(donors_minus[i]));
	substringlist_gc(&(antidonors_minus[i]));
	substringlist_gc(&(acceptors_minus[i]));
	substringlist_gc(&(antiacceptors_minus[i]));
      }
      FREEA(donors_plus);
      FREEA(antidonors_plus);
      FREEA(acceptors_plus);
      FREEA(antiacceptors_plus);
      FREEA(donors_minus);
      FREEA(antidonors_minus);
      FREEA(acceptors_minus);
      FREEA(antiacceptors_minus);
    }

    debug(printf("%d single splices, %d double splices, %d long single splices, %d distant splices\n",
		 List_length(singlesplicing),List_length(doublesplicing),
		 List_length(longsinglesplicing),List_length(distantsplicing)));
  }


  hits = List_append(hits,
		     List_append(longsinglesplicing,distantsplicing));

  /* Search 8: Terminals */

  /* Previously criterion for skipping find_terminals was (greedy ||
     subs || indels || singlesplicing || doublesplicing ||
     shortendsplicing || longsinglesplicing || distantsplicing) */
  if (0 && found_score > opt_level) { 
    terminals = find_terminals(plus_anchor_segments,minus_anchor_segments,
			       n_plus_anchors,n_minus_anchors,querylength,query_lastpos,
			       query_compress_fwd,query_compress_rev,
			       /*max_mismatches_allowed*/done_level,genestrand);

  }

  /* Search 9: GMAP via segments */
  gmapp = true;
  if (gmap_segments_p == false) {
    gmapp = false;
  } else if (found_score < trigger_score_for_gmap) {
    debug(printf("Test for stage 9: false because found_score %d < trigger_score_for_gmap %d\n",found_score,trigger_score_for_gmap));
    gmapp = false;
  } else if (min_trim < min_distantsplicing_end_matches) {
    gmapp = false;
  } else if (distantsplicing != NULL) {
    gmapp = false;
  }


  gmap_hits = (List_T) NULL;
  if (gmapp == true) {
    gmap_hits = convert_plus_segments_to_gmap(/*hits*/NULL,queryuc_ptr,querylength,query_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
					      queryrc,Shortread_invertedp(queryseq),
#endif
					      plus_anchor_segments,n_plus_anchors,this->plus_segments,this->plus_nsegments,
					      oligoindices_minor,
					      pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
					      genestrand,/*require_pairing_p*/false);
    gmap_hits = convert_minus_segments_to_gmap(/*hits*/gmap_hits,queryuc_ptr,querylength,query_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
					       queryrc,Shortread_invertedp(queryseq),
#endif
					       minus_anchor_segments,n_minus_anchors,this->minus_segments,this->minus_nsegments,
					       oligoindices_minor,
					       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
					       genestrand,/*require_pairing_p*/false);

    opt_level = (found_score < opt_level) ? found_score : opt_level;
    if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
      done_level = user_maxlevel;
    }
    debug(printf("10> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
  }


  /* Search 10: GMAP improvement */
  debug13(printf("%d hits (vs max_gmap_improvement %d)\n",List_length(gmap_hits),max_gmap_improvement));
  if (hits != NULL && gmap_improvement_p == true) {
    /* 11.  GMAP terminal */
    
    /* This is done for paired-ends, but should not be necessary for single-end */
    debug13(printf("Before remove overlaps at cutoff level %d: %d hits\n",opt_level,List_length(hits)));
    hits = Stage3end_sort_bymatches(Stage3end_remove_overlaps(hits,/*finalp*/false));
    debug13(printf("After remove overlaps: %d\n",List_length(hits)));

    i = 0;
    debug13(printf("%d hits\n",List_length(hits)));
    debug13(printf("For each hit, running GMAP on single end to match with hit\n"));

    for (q = hits; q != NULL && i < max_gmap_improvement; q = q->rest) {
      hit = (Stage3end_T) q->first;
      align_single_hit_with_gmap(&gmap1,&gmap2,hit,queryuc_ptr,querylength,
#ifdef END_KNOWNSPLICING_SHORTCUT
				 queryrc,Shortread_invertedp(queryseq),
#endif
				 oligoindices_minor,
				 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				 genestrand,/*extend_ends_p*/true);
      if (gmap1 != NULL) {
	debug13(missing_hit = querylength - Stage3end_nmatches_posttrim(hit));
	debug13(missing_gmap = querylength - Stage3end_nmatches_posttrim(gmap1));
	debug13(printf("GMAP %p with %d matches, %d missing compared with original terminal with %d matches, %d missing\n",
		       gmap1,Stage3end_nmatches_posttrim(gmap1),missing_gmap,Stage3end_nmatches_posttrim(hit),missing_hit));
	gmap_hits = List_push(gmap_hits,(void *) gmap1);
	if (gmap2 != NULL) {
	  gmap_hits = List_push(gmap_hits,(void *) gmap2);
	}
	Stage3end_set_improved_by_gmap(hit);
      }
    }
  }
  debug13(printf("Have %d GMAP hits\n",List_length(gmap_hits)));

  if (alloc_floors_p == true) {
    Floors_free(&floors);
  }

  /* Keep gmap_hits found in search 9 and 10 */
  if (gmap_hits != NULL) {
    hits = List_append(hits,gmap_hits);
  }

  if (gmap_improvement_p == false) {
    debug(printf("No GMAP improvement: Before remove_overlaps at cutoff level %d: %d\n",*cutoff_level,List_length(hits)));
    hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				   querylength,/*keep_gmap_p*/true,/*finalp*/true);
    /* hits = Stage3end_reject_trimlengths(hits); */
    hits = Stage3end_remove_overlaps(hits,/*finalp*/true);
    hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				   querylength,/*keep_gmap_p*/false,/*finalp*/true);
    hits = Stage3end_resolve_multimapping(hits);
    debug(printf("After remove_overlaps (1): %d\n",List_length(hits)));

  } else {
    debug(printf("GMAP improvement: Before remove_overlaps at cutoff level %d: %d\n",*cutoff_level,List_length(hits)));
    hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				   querylength,/*keep_gmap_p*/true,/*finalp*/false);
    /* Don't reject based on trimlength until after GMAP improvements */
    hits = Stage3end_remove_overlaps(hits,/*finalp*/false);
    hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				   querylength,/*keep_gmap_p*/false,/*finalp*/false);
    hits = Stage3end_resolve_multimapping(hits);
    debug(printf("After remove_overlaps (2): %d\n",List_length(hits)));

    hits = align_singleend_with_gmap(hits,queryuc_ptr,querylength,oligoindices_minor,
				     pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,user_maxlevel);
    hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				   querylength,/*keep_gmap_p*/true,/*finalp*/true);
    /* hits = Stage3end_reject_trimlengths(hits); */
    hits = Stage3end_remove_overlaps(hits,/*finalp*/true);
    hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				   querylength,/*keep_gmap_p*/false,/*finalp*/true);
    hits = Stage3end_resolve_multimapping(hits);
  }

  hits = Stage3end_remove_circular_alias(hits);
  hits = Stage3end_remove_duplicates(hits); /* Aliases can cause duplicates */

  FREE(plus_anchor_segments);
  FREE(minus_anchor_segments);

  return hits;
}


static Stage3end_T *
single_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
	     Shortread_T queryseq, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev,
	     int indexdb_size_threshold, Floors_T *floors_array,
	     double user_maxlevel_float, double user_mincoverage_float,
	     int indel_penalty_middle, int indel_penalty_end,
	     bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
	     int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
	     Oligoindex_array_T oligoindices_minor,
	     Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
	     Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	     bool keep_floors_p) {
  Stage3end_T *stage3array;
  History_T gmap_history;
  List_T hits = NULL;
  T this = NULL;
  int user_maxlevel, min_coverage;
  int querylength, query_lastpos, cutoff_level;
  char *queryuc_ptr, *quality_string;
  Compress_T query_compress_fwd = NULL, query_compress_rev = NULL;
  char *queryrc;

  querylength = Shortread_fulllength(queryseq);

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    queryrc = (char *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    queryrc = (char *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  queryrc = (char *) MALLOC((querylength+1)*sizeof(int));
#endif

  if (user_maxlevel_float < 0.0) {
    user_maxlevel = -1;
  } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
    user_maxlevel = (int) rint(user_maxlevel_float * (double) querylength);
  } else {
    user_maxlevel = (int) user_maxlevel_float;
  }

  if (user_mincoverage_float < 0.0) {
    min_coverage = 0;
  } else if (user_mincoverage_float > 0.0 && user_mincoverage_float < 1.0) {
    min_coverage = (int) rint(user_mincoverage_float * (double) querylength);
  } else {
    min_coverage = (int) user_mincoverage_float;
  }

  /* Limit search on repetitive sequences */
  queryuc_ptr = Shortread_fullpointer_uc(queryseq);
  quality_string = Shortread_quality_string(queryseq);
  if (check_dinucleotides(queryuc_ptr,querylength) == false) {
    user_maxlevel = 0;
  }

  query_compress_fwd = Compress_new_fwd(queryuc_ptr,querylength);
  query_compress_rev = Compress_new_rev(queryuc_ptr,querylength);
  make_complement_buffered(queryrc,queryuc_ptr,querylength);

  this = Stage1_new(querylength);
  query_lastpos = querylength - index1part;

  gmap_history = History_new();
  hits = align_end(&cutoff_level,this,query_compress_fwd,query_compress_rev,
		   queryuc_ptr,queryrc,querylength,query_lastpos,
		   indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
		   oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
		   user_maxlevel,indel_penalty_middle,indel_penalty_end,
		   localsplicing_penalty,distantsplicing_penalty,min_shortend,
		   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
		   keep_floors_p,/*genestrand*/0,/*first_read_p*/true);

  hits = Stage3end_filter_coverage(hits,min_coverage);
  if (hits == NULL) {
    *npaths_primary = *npaths_altloc = 0;
    stage3array = (Stage3end_T *) NULL;
  } else {
    Stage3end_count_hits(&(*npaths_primary),&(*npaths_altloc),hits);
    stage3array = (Stage3end_T *) List_to_array_out(hits,NULL); List_free(&hits); /* Return value */
    stage3array = Stage3end_eval_and_sort(/*npaths*/(*npaths_primary) + (*npaths_altloc),
					  &(*first_absmq),&(*second_absmq),
					  stage3array,maxpaths_search,queryseq,queryuc_ptr,
					  query_compress_fwd,query_compress_rev,
					  quality_string,/*displayp*/true);
  }
     
  History_free(&gmap_history);
  Compress_free(&query_compress_fwd);
  Compress_free(&query_compress_rev);
  Stage1_free(&this,querylength); 

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(queryrc);
  } else {
    FREE(queryrc);
  }
#else
  FREE(queryrc);
#endif

  return stage3array;
}


static Stage3end_T *
single_read_tolerant_nonstranded (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
				  Shortread_T queryseq, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev,
				  int indexdb_size_threshold, Floors_T *floors_array,
				  double user_maxlevel_float, double user_mincoverage_float,
				  int indel_penalty_middle, int indel_penalty_end,
				  bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
				  int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
				  Oligoindex_array_T oligoindices_minor, 
				  Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
				  Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				  bool keep_floors_p) {
  Stage3end_T *stage3array;
  History_T gmap_history;
  List_T hits, hits_geneplus = NULL, hits_geneminus = NULL;
  T this_geneplus = NULL, this_geneminus = NULL;
  int user_maxlevel, min_coverage;
  int querylength, query_lastpos, cutoff_level;
  char *queryuc_ptr, *quality_string;
  Compress_T query_compress_fwd = NULL, query_compress_rev = NULL;
  bool allvalidp;
  char *queryrc;

  querylength = Shortread_fulllength(queryseq);

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    queryrc = (char *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    queryrc = (char *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  queryrc = (char *) MALLOC((querylength+1)*sizeof(int));
#endif

  if (user_maxlevel_float < 0.0) {
    user_maxlevel = -1;
  } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
    user_maxlevel = (int) rint(user_maxlevel_float * (double) querylength);
  } else {
    user_maxlevel = (int) user_maxlevel_float;
  }

  if (user_mincoverage_float < 0.0) {
    min_coverage = 0;
  } else if (user_mincoverage_float > 0.0 && user_mincoverage_float < 1.0) {
    min_coverage = (int) rint(user_mincoverage_float * (double) querylength);
  } else {
    min_coverage = (int) user_mincoverage_float;
  }

  this_geneplus = Stage1_new(querylength);
  this_geneminus = Stage1_new(querylength);

  queryuc_ptr = Shortread_fullpointer_uc(queryseq);
  quality_string = Shortread_quality_string(queryseq);
  query_lastpos = querylength - index1part;

  /* Limit search on repetitive sequences */
  if (check_dinucleotides(queryuc_ptr,querylength) == false) {
    user_maxlevel = 0;
  }

  query_compress_fwd = Compress_new_fwd(queryuc_ptr,querylength);
  query_compress_rev = Compress_new_rev(queryuc_ptr,querylength);
  gmap_history = History_new();
  make_complement_buffered(queryrc,queryuc_ptr,querylength);

  if (read_oligos(&allvalidp,this_geneplus,queryuc_ptr,querylength,query_lastpos,/*genestrand*/+1) > 0) {
    hits_geneplus = align_end(&cutoff_level,this_geneplus,query_compress_fwd,query_compress_rev,
			      queryuc_ptr,queryrc,querylength,query_lastpos,
			      indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
			      oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
			      user_maxlevel,indel_penalty_middle,indel_penalty_end,
			      localsplicing_penalty,distantsplicing_penalty,min_shortend,
			      allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			      keep_floors_p,/*genestrand*/+1,/*first_read_p*/true);
  }

  if (read_oligos(&allvalidp,this_geneminus,queryuc_ptr,querylength,query_lastpos,/*genestrand*/+2) > 0) {
    hits_geneminus = align_end(&cutoff_level,this_geneminus,query_compress_fwd,query_compress_rev,
			       queryuc_ptr,queryrc,querylength,query_lastpos,
			       indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
			       oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
			       user_maxlevel,indel_penalty_middle,indel_penalty_end,
			       localsplicing_penalty,distantsplicing_penalty,min_shortend,
			       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			       keep_floors_p,/*genestrand*/+2,/*first_read_p*/true);
  }

  hits = List_append(hits_geneplus,hits_geneminus);
  hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				 querylength,/*keep_gmap_p*/true,/*finalp*/true);
  /* hits = Stage3end_reject_trimlengths(hits); */
  hits = Stage3end_remove_overlaps(hits,/*finalp*/true);
  hits = Stage3end_optimal_score(hits,query_compress_fwd,query_compress_rev,
				 querylength,/*keep_gmap_p*/false,/*finalp*/true);
  hits = Stage3end_resolve_multimapping(hits);

  hits = Stage3end_filter_coverage(hits,min_coverage);
  if (hits == NULL) {
    *npaths_primary = *npaths_altloc = 0;
    stage3array = (Stage3end_T *) NULL;
  } else {
    Stage3end_count_hits(&(*npaths_primary),&(*npaths_altloc),hits);
    stage3array = (Stage3end_T *) List_to_array_out(hits,NULL); List_free(&hits); /* Return value */
    stage3array = Stage3end_eval_and_sort(/*npaths*/(*npaths_primary) + (*npaths_altloc),
					  &(*first_absmq),&(*second_absmq),
					  stage3array,maxpaths_search,queryseq,queryuc_ptr,
					  query_compress_fwd,query_compress_rev,
					  quality_string,/*displayp*/true);
  }

  History_free(&gmap_history);
  Compress_free(&query_compress_fwd);
  Compress_free(&query_compress_rev);
  Stage1_free(&this_geneminus,querylength); 
  Stage1_free(&this_geneplus,querylength); 

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(queryrc);
  } else {
    FREE(queryrc);
  }
#else
  FREE(queryrc);
#endif

  return stage3array;
}


Stage3end_T *
Stage1_single_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    Shortread_T queryseq, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev,
		    int indexdb_size_threshold, Floors_T *floors_array,
		    double user_maxlevel_float, double user_mincoverage_float,
		    int indel_penalty_middle, int indel_penalty_end,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
		    Oligoindex_array_T oligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    bool keep_floors_p) {

  if (mode == STANDARD || mode == CMET_STRANDED || mode == ATOI_STRANDED || mode == TTOC_STRANDED) {
    return single_read(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
		       queryseq,indexdb_fwd,indexdb_rev,indexdb_size_threshold,
		       floors_array,user_maxlevel_float,user_mincoverage_float,
		       indel_penalty_middle,indel_penalty_end,
		       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
		       localsplicing_penalty,distantsplicing_penalty,min_shortend,
		       oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,keep_floors_p);
  } else if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
    return single_read_tolerant_nonstranded(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),queryseq,
					    indexdb_fwd,indexdb_rev,indexdb_size_threshold,
					    floors_array,user_maxlevel_float,user_mincoverage_float,
					    indel_penalty_middle,indel_penalty_end,
					    allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
					    localsplicing_penalty,distantsplicing_penalty,min_shortend,
					    oligoindices_minor,
					    pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,keep_floors_p);
  } else {
    fprintf(stderr,"Do not recognize mode %d\n",mode);
    abort();
  }
}



/* #define HITARRAY_SHORTENDSPLICING 4 */
/* #define HITARRAY_DISTANTSPLICING 4 */


static List_T
align_halfmapping_with_gmap (History_T gmap_history, Stage3end_T hit5, Stage3end_T hit3, 
			     Shortread_T queryseq5, Shortread_T queryseq3,
			     char *queryuc_ptr, int querylength, int query_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
			     char *queryrc, bool invertedp,
#endif
			     struct Segment_T *plus_segments, int plus_nsegments,
			     struct Segment_T *minus_segments, int minus_nsegments,
			     Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
			     Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			     Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			     Chrpos_T pairmax, Chrpos_T shortsplicedist, int genestrand) {
  List_T hits = NULL;
  int sensedir, sense_try;

  /* int zero_offset = 0; */
  Univcoord_T segmentstart, segmentend;
  Univcoord_T genomicbound, mappingstart, mappingend,
    chroffset, chrhigh, mappingpos;
#ifdef USE_GREEDY
  Univcoord_T close_mappingstart_greedy, close_mappingend_greedy,
    middle_mappingstart_greedy, middle_mappingend_greedy;
#endif
  Univcoord_T close_mappingstart_last, close_mappingend_last,
    middle_mappingstart_last, middle_mappingend_last;
  Univcoord_T knownsplice_limit_low, knownsplice_limit_high;
  Univcoord_T close_knownsplice_limit_low, close_knownsplice_limit_high;
  Chrpos_T chrlength;
  Chrnum_T chrnum;
  bool close_mappingstart_p = false, close_mappingend_p = false;
  bool middle_mappingstart_p = false, middle_mappingend_p = false;
  bool fallback_mappingstart_p, fallback_mappingend_p;
  bool good_start_p, good_end_p, watsonp, favor_right_p;

  int starti, endi, i;

  if (hit3 == NULL) {
    /* Both events are tested by Stage3end_anomalous_splice_p */
    if ((chrnum = Stage3end_chrnum(hit5)) == 0) {
      /* Translocation */
      return (List_T) NULL;

    } else if (Stage3end_hittype(hit5) == SAMECHR_SPLICE) {
      /* A genomic event that doesn't get reflected in chrnum */
      return (List_T) NULL;

    } else if ((watsonp = Stage3end_plusp(hit5)) == true) {
      chroffset = Stage3end_chroffset(hit5);
      chrhigh = Stage3end_chrhigh(hit5);
      chrlength = Stage3end_chrlength(hit5);

      if (Shortread_find_primers(queryseq5,queryseq3) == true) {
	/* Go from genomicstart */
	debug13(printf("Found primers\n"));
	genomicbound = Stage3end_genomicstart(hit5);

      } else if (Stage3end_anomalous_splice_p(hit5) == true) {
	/* Go from genomicstart */
	debug13(printf("Anomalous splice\n"));
	genomicbound = Stage3end_genomicstart(hit5);

      } else {
	genomicbound = Stage3end_genomicend(hit5);

#if 0
	/* TODO: Previously called Shortread_find_overlap.  Now with Shortread_max_overlap, can optimize this code */
	if ((overlap = Shortread_max_overlap(queryseq5,queryseq3)) > 0 &&
	    Stage3end_genomicbound_from_end(&genomicbound2,hit5,overlap,chroffset) == true) {
	  debug13(printf("Found overlap of %d\n",overlap));
	  if (genomicbound2 < genomicbound) {
	    zero_offset = genomicbound - genomicbound2;
	    genomicbound = genomicbound2;
	  }
	}
#endif
      }

      debug13(printf("Case 1: hit5 plus %s %u..%u (sensedir %d) => genomicbound %u\n",
		     Stage3end_hittype_string(hit5),
		     Stage3end_genomicstart(hit5) - chroffset,Stage3end_genomicend(hit5) - chroffset,
		     Stage3end_sensedir(hit5),genomicbound - chroffset));

      knownsplice_limit_low = mappingstart = segmentstart = genomicbound;
      knownsplice_limit_high =  add_bounded(Stage3end_genomicend(hit5),pairmax + overall_max_distance,chrhigh);
      segmentend = add_bounded(Stage3end_genomicend(hit5),pairmax + PAIRMAX_ADDITIONAL,chrhigh);
#ifdef LONG_ENDSPLICES
      mappingend = add_bounded(Stage3end_genomicend(hit5),pairmax + overall_max_distance,chrhigh);
#else
      mappingend = add_bounded(Stage3end_genomicend(hit5),pairmax + overall_max_distance_novelend,chrhigh);
      debug13(printf("Original bounds E: knownsplice_limit_low %u, knownsplice_limit_high %u, mappingend %u\n",
		     knownsplice_limit_low - chroffset,knownsplice_limit_high - chroffset,mappingend - chroffset));
#endif

      close_mappingend_last = middle_mappingend_last = Stage3end_genomicend(hit5);
#ifdef USE_GREEDY
      close_mappingend_greedy = middle_mappingend_greedy = segmentend;
#endif

      if (plus_nsegments > 0) {
	/* Use segments to bound */
	debug13(printf("Finding segments from segmentstart %u to segmentend %u (plus_nsegments %d)\n",
		       segmentstart - chroffset,segmentend - chroffset,plus_nsegments));
	starti = endi = -1;
	i = binary_search_segments(0,plus_nsegments-1,plus_segments,segmentstart);
	while (i < plus_nsegments - 1 && plus_segments[i].diagonal == (Univcoord_T) -1) {
	  i++;
	}
	starti = i;
	while (plus_segments[i].diagonal < segmentend) {
	  endi = i;
	  i++;
	}
	if (starti >= 0 && endi >= 0) {
	  debug13(printf("starti = %d, endi = %d\n",starti,endi));
	  assert(starti <= endi);
	  for (i = starti; i <= endi; i++) {
	    debug13(printf("diagonal %u (%llu), querypos %d..%d\n",
			   (Chrpos_T) (plus_segments[i].diagonal - chroffset),(unsigned long long) plus_segments[i].diagonal,
			   plus_segments[i].querypos5,plus_segments[i].querypos3));
	    if (query_lastpos - plus_segments[i].querypos3 >= STAGE2_MIN_OLIGO + index1interval) {
	      /* Case 1. Missing end of query, so there could be a middle splice */
	      debug13b(printf("  query_lastpos %d - querypos3 %d >= %d + %d, so using this diagonal plus overall_max_distance\n",
			      query_lastpos,plus_segments[i].querypos3,STAGE2_MIN_OLIGO,index1interval));
#ifdef USE_GREEDY
	      if ((mappingpos = add_bounded(plus_segments[i].diagonal,shortsplicedist_novelend,chrhigh)) < middle_mappingend_greedy &&
		  mappingpos > genomicbound) {
		middle_mappingend_greedy = mappingpos;
		middle_mappingend_p = true;
		debug13(printf("  Redefining middle mappingend greedy to %u\n",middle_mappingend_greedy - chroffset));
	      }
#endif

#ifdef LONG_ENDSPLICES
	      if ((mappingpos = add_bounded(plus_segments[i].diagonal,overall_max_distance,chrhigh)) > middle_mappingend_last) {
		/* Use > for NOT_GREEDY */
		middle_mappingend_last = mappingpos;
		middle_mappingend_p = true;
		debug13(printf("  Redefining middle mappingend last to %u\n",middle_mappingend_last - chroffset));
	      }
#else
	      if ((mappingpos = plus_segments[i].diagonal) > middle_mappingend_last) {
		/* Use > for NOT_GREEDY */
		middle_mappingend_last = mappingpos;
		middle_mappingend_p = true;
		debug13(printf("  Redefining middle mappingend last to %u\n",middle_mappingend_last - chroffset));
	      }
#endif

	    } else {
	      debug13b(printf("  query_lastpos %d - querypos3 %d < %d + %d, so using this diagonal\n",
			      query_lastpos,plus_segments[i].querypos3,STAGE2_MIN_OLIGO,index1interval));
#ifdef USE_GREEDY
	      if ((mappingpos = plus_segments[i].diagonal) < close_mappingend_greedy &&
		  mappingpos > genomicbound) {
		close_mappingend_greedy = mappingpos;
		close_mappingend_p = true;
		debug13(printf("  Redefining close mappingend greedy to %u\n",close_mappingend_greedy - chroffset));
	      }
#endif
	      if ((mappingpos = plus_segments[i].diagonal) > close_mappingend_last) {
		/* Use > for NOT_GREEDY */
		close_mappingend_last = mappingpos;
		close_mappingend_p = true;
		debug13(printf("  Redefining close mappingend last to %u\n",close_mappingend_last - chroffset));
	      }
	    }
	  }

#ifdef USE_GREEDY
	  if (close_mappingend_p == true) {
	    close_knownsplice_limit_high = add_bounded(close_mappingend_greedy,overall_max_distance,chrhigh);
	  } else if (middle_mappingend_p == true) {
	    debug13(printf("Using middle mappingend\n"));
	    close_knownsplice_limit_high = middle_mappingend_greedy;
	    close_mappingend_greedy = middle_mappingend_greedy;
	    close_mappingend_p = true;
	  }
#else
	  if (close_mappingend_p == true) {
	    close_knownsplice_limit_high = add_bounded(close_mappingend_last,shortsplicedist,chrhigh);
	  } else if (middle_mappingend_p == true) {
	    debug13(printf("Using middle mappingend\n"));
	    close_knownsplice_limit_high = middle_mappingend_last;
	    close_mappingend_last = middle_mappingend_last;
	    close_mappingend_p = true;
	  }
#endif
#ifdef USE_GREEDY
	  if (middle_mappingend_p == true && middle_mappingend_last > close_mappingend_greedy) {
	    knownsplice_limit_high = middle_mappingend_last;
	    mappingend = middle_mappingend_last;
	  } else if (close_mappingend_p == true && close_mappingend_last != close_mappingend_greedy) {
	    knownsplice_limit_high = add_bounded(close_mappingend_last,overall_max_distance,chrhigh);
	    mappingend = close_mappingend_last;
	  }
#else
	  if (middle_mappingend_p == true && middle_mappingend_last > close_mappingend_last) {
	    knownsplice_limit_high = middle_mappingend_last;
	    mappingend = middle_mappingend_last;
	  }
#endif

	  if (close_mappingend_p == false) {
	    fallback_mappingend_p = false;
#ifdef USE_GREEDY
	  } else if (mappingend <= close_mappingend_greedy) {
	    fallback_mappingend_p = false;
#endif
	  } else {
	    debug13(printf("Fallback mappingend = %u\n",mappingend - chroffset));
	    fallback_mappingend_p = true;
	  }
	}
      }

      favor_right_p = false;

    } else {
      chroffset = Stage3end_chroffset(hit5);
      chrhigh = Stage3end_chrhigh(hit5);
      chrlength = Stage3end_chrlength(hit5);

      if (Shortread_find_primers(queryseq5,queryseq3) == true) {
	/* Go from genomicstart */
	debug13(printf("Found primers\n"));
	genomicbound = Stage3end_genomicstart(hit5);

      } else if (Stage3end_anomalous_splice_p(hit5) == true) {
	/* Go from genomicstart */
	debug13(printf("Anomalous splice\n"));
	genomicbound = Stage3end_genomicstart(hit5);

      } else {
	genomicbound = Stage3end_genomicend(hit5);

#if 0
	/* TODO: Previously called Shortread_find_overlap.  Now with Shortread_max_overlap, can optimize this code */
	if ((overlap = Shortread_max_overlap(queryseq5,queryseq3)) > 0 &&
	    Stage3end_genomicbound_from_end(&genomicbound2,hit5,overlap,chroffset) == true) {
	  debug13(printf("Found overlap of %d\n",overlap));
	  if (genomicbound2 > genomicbound) {
	    zero_offset = genomicbound2 - genomicbound;
	    genomicbound = genomicbound2;
	  }
	}
#endif
      }

      debug13(printf("Case 2: hit5 minus %s %u..%u (sensedir %d) => genomicbound %u\n",
		     Stage3end_hittype_string(hit5),
		     Stage3end_genomicstart(hit5) - chroffset,Stage3end_genomicend(hit5) - chroffset,
		     Stage3end_sensedir(hit5),genomicbound - chroffset));

      knownsplice_limit_high = mappingend = segmentend = genomicbound;
      knownsplice_limit_low = subtract_bounded(Stage3end_genomicend(hit5),pairmax + overall_max_distance,chroffset);
      segmentstart = subtract_bounded(Stage3end_genomicend(hit5),pairmax + PAIRMAX_ADDITIONAL,chroffset);
#ifdef LONG_ENDSPLICES
      mappingstart = subtract_bounded(Stage3end_genomicend(hit5),pairmax + overall_max_distance,chroffset);
#else
      mappingstart = subtract_bounded(Stage3end_genomicend(hit5),pairmax + overall_max_distance_novelend,chroffset);
#endif
      debug13(printf("Original bounds F: knownsplice_limit_low %u, knownsplice_limit_high %u, mappingstart %u\n",
		     knownsplice_limit_low - chroffset,knownsplice_limit_high - chroffset,mappingstart - chroffset));

      close_mappingstart_last = middle_mappingstart_last = Stage3end_genomicend(hit5);
#ifdef USE_GREEDY
      close_mappingstart_greedy = middle_mappingstart_greedy = segmentstart;
#endif

      if (minus_nsegments > 0) {
	/* Use segments to bound */
	debug13(printf("Finding segments from segmentstart %u to segmentend %u (minus_nsegments %d)\n",
		       segmentstart - chroffset,segmentend - chroffset,minus_nsegments));
	starti = endi = -1;
	i = binary_search_segments(0,minus_nsegments-1,minus_segments,segmentend);
	while (i >= 0 && minus_segments[i].diagonal >= segmentend) {
	  i--;
	}
	starti = i;
	while (i >= 0 && minus_segments[i].diagonal > segmentstart) {
	  if (minus_segments[i].diagonal < (Univcoord_T) -1) {
	    endi = i;
	  }
	  i--;
	}
	if (starti >= 0 && endi >= 0) {
	  debug13(printf("starti = %d, endi = %d\n",starti,endi));
	  assert(starti >= endi);
	  for (i = starti; i >= endi; i--) {
	    debug13(printf("diagonal %u (%llu), querypos %d..%d\n",
			   (Chrpos_T) (minus_segments[i].diagonal - chroffset),(unsigned long long) minus_segments[i].diagonal,
			   minus_segments[i].querypos5,minus_segments[i].querypos3));
	    if (query_lastpos - minus_segments[i].querypos3 >= STAGE2_MIN_OLIGO + index1interval) {
	      /* Case 2. Missing end of query, so there could be a middle splice */
	      debug13b(printf("  query_lastpos %d - querypos3 %d >= %d + %d, so using this diagonal plus overall_max_distance\n",
			      query_lastpos,minus_segments[i].querypos3,STAGE2_MIN_OLIGO,index1interval));
#ifdef USE_GREEDY
	      if ((mappingpos = subtract_bounded(minus_segments[i].diagonal,querylength + shortsplicedist_novelend,chroffset)) > middle_mappingstart_greedy &&
		  mappingpos < genomicbound) {
		middle_mappingstart_greedy = mappingpos;
		middle_mappingstart_p = true;
		debug13(printf("  Redefining middle mappingstart greedy to %u\n",middle_mappingstart_greedy - chroffset));
	      }
#endif
#ifdef LONG_ENDSPLICES
	      if ((mappingpos = subtract_bounded(minus_segments[i].diagonal,querylength + overall_max_distance,chroffset)) < middle_mappingstart_last) {
		/* Use < for NOT_GREEDY */
		middle_mappingstart_last = mappingpos;
		middle_mappingstart_p = true;
		debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	      }
#else
	      if ((mappingpos = subtract_bounded(minus_segments[i].diagonal,querylength,chroffset)) < middle_mappingstart_last) {
		/* Use < for NOT_GREEDY */
		middle_mappingstart_last = mappingpos;
		middle_mappingstart_p = true;
		debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	      }
#endif

	    } else {
	      debug13b(printf("  query_lastpos %d - querypos3 %d < %d + %d, so using this diagonal\n",
			      query_lastpos,minus_segments[i].querypos3,STAGE2_MIN_OLIGO,index1interval));
#ifdef USE_GREEDY
	      if ((mappingpos = subtract_bounded(minus_segments[i].diagonal,querylength,chroffset)) > close_mappingstart_greedy &&
		  mappingpos < genomicbound) {
		close_mappingstart_greedy = mappingpos;
		close_mappingstart_p = true;
		debug13(printf("  Redefining close mappingstart greedy to %u\n",close_mappingstart_greedy - chroffset));
	      }
#endif
	      if ((mappingpos = subtract_bounded(minus_segments[i].diagonal,querylength,chroffset)) < close_mappingstart_last) {
		/* Use < for NOT_GREEDY */
		close_mappingstart_last = mappingpos;
		close_mappingstart_p = true;
		debug13(printf("  Redefining close mappingstart last to %u\n",close_mappingstart_last - chroffset));
	      }
	    }
	  }

#ifdef USE_GREEDY
	  if (close_mappingstart_p == true) {
	    close_knownsplice_limit_low = subtract_bounded(close_mappingstart_greedy,overall_max_distance,chroffset);
	  } else if (middle_mappingstart_p == true) {
	    debug13(printf("Using middle mappingstart\n"));
	    close_knownsplice_limit_low = middle_mappingstart_greedy;
	    close_mappingstart_greedy = middle_mappingstart_greedy;
	    close_mappingstart_p = true;
	  }
#else
	  if (close_mappingstart_p == true) {
	    close_knownsplice_limit_low = subtract_bounded(close_mappingstart_last,shortsplicedist,chroffset);
	  } else if (middle_mappingstart_p == true) {
	    debug13(printf("Using middle mappingstart\n"));
	    close_knownsplice_limit_low = middle_mappingstart_last;
	    close_mappingstart_last = middle_mappingstart_last;
	    close_mappingstart_p = true;
	  }
#endif
#ifdef USE_GREEDY
	  if (middle_mappingstart_p == true && middle_mappingstart_last < close_mappingstart_greedy) {
	    knownsplice_limit_low = middle_mappingstart_last;
	    mappingstart = middle_mappingstart_last;
	  } else if (close_mappingstart_p == true && close_mappingstart_last != close_mappingstart_greedy) {
	    knownsplice_limit_low = subtract_bounded(close_mappingstart_last,overall_max_distance,chroffset);
	    mappingstart = close_mappingstart_last;
	  }
#else
	  if (middle_mappingstart_p == true && middle_mappingstart_last < close_mappingstart_last) {
	    knownsplice_limit_low = middle_mappingstart_last;
	    mappingstart = middle_mappingstart_last;
	  }
#endif
	  if (close_mappingstart_p == false) {
	    fallback_mappingstart_p = false;
#ifdef USE_GREEDY
	  } else if (mappingstart >= close_mappingstart_greedy) {
	    fallback_mappingstart_p = false;
#endif
	  } else {
	    debug13(printf("Fallback mappingstart = %u\n",mappingstart - chroffset));
	    fallback_mappingstart_p = true;
	  }
	}
      }

      favor_right_p = false;
    }

    if ((sensedir = Stage3end_sensedir(hit5)) == SENSE_FORWARD) {
      sense_try = +1;
    } else if (sensedir == SENSE_ANTI) {
      sense_try = -1;
    } else {
      sense_try = 0;
    }

  } else if (hit5 == NULL) {
    /* Both events are tested by Stage3end_anomalous_splice_p */
    if ((chrnum = Stage3end_chrnum(hit3)) == 0) {
      /* Translocation */
      return (List_T) NULL;

    } else if (Stage3end_hittype(hit3) == SAMECHR_SPLICE) {
      /* A genomic event that doesn't get reflected in chrnum */
      return (List_T) NULL;

    } else if ((watsonp = Stage3end_plusp(hit3)) == true) {
      chroffset = Stage3end_chroffset(hit3);
      chrhigh = Stage3end_chrhigh(hit3);
      chrlength = Stage3end_chrlength(hit3);

      if (Shortread_find_primers(queryseq5,queryseq3) == true) {
	/* Go from genomicend */
	debug13(printf("Found primers\n"));
	genomicbound = Stage3end_genomicend(hit3);

      } else if (Stage3end_anomalous_splice_p(hit3) == true) {
	/* Go from genomicend */
	debug13(printf("Anomalous splice\n"));
	genomicbound = Stage3end_genomicend(hit3);

      } else {
	genomicbound = Stage3end_genomicstart(hit3);

#if 0
	/* TODO: Previously called Shortread_find_overlap.  Now with Shortread_max_overlap, can optimize this code */
	if ((overlap = Shortread_max_overlap(queryseq5,queryseq3)) > 0 &&
	    Stage3end_genomicbound_from_start(&genomicbound2,hit3,overlap,chroffset) == true) {
	  debug13(printf("Found overlap of %d\n",overlap));
	  if (genomicbound2 > genomicbound) {
	    zero_offset = genomicbound2 - genomicbound;
	    genomicbound = genomicbound2;
	  }
	}
#endif
      }

      debug13(printf("Case 3: hit3 plus %s %u..%u (sensedir %d) => genomicbound %u\n",
		     Stage3end_hittype_string(hit3),
		     Stage3end_genomicstart(hit3) - chroffset,Stage3end_genomicend(hit3) - chroffset,
		     Stage3end_sensedir(hit3),genomicbound - chroffset));

      knownsplice_limit_high = mappingend = segmentend = genomicbound;
      knownsplice_limit_low = subtract_bounded(Stage3end_genomicstart(hit3),pairmax + overall_max_distance,chroffset);
      segmentstart = subtract_bounded(Stage3end_genomicstart(hit3),pairmax + PAIRMAX_ADDITIONAL,chroffset);
#ifdef LONG_ENDSPLICES
      mappingstart = subtract_bounded(Stage3end_genomicstart(hit3),pairmax + overall_max_distance,chroffset);
#else
      mappingstart = subtract_bounded(Stage3end_genomicstart(hit3),pairmax + overall_max_distance_novelend,chroffset);
#endif

      close_mappingstart_last = middle_mappingstart_last = Stage3end_genomicstart(hit3);
#ifdef USE_GREEDY
      close_mappingstart_greedy = middle_mappingstart_greedy = segmentstart;
#endif

      if (plus_nsegments > 0) {
	/* Use segments to bound */
	debug13(printf("Finding segments from segmentstart %u to segmentend %u (plus_nsegments %d)\n",
		       segmentstart - chroffset,segmentend - chroffset,plus_nsegments));
	starti = endi = -1;
	i = binary_search_segments(0,plus_nsegments-1,plus_segments,segmentend);
	while (i >= 0 && plus_segments[i].diagonal >= segmentend) {
	  i--;
	}
	starti = i;
	while (i >= 0 && plus_segments[i].diagonal > segmentstart) {
	  if (plus_segments[i].diagonal < (Univcoord_T) -1) {
	    endi = i;
	  }
	  i--;
	}
	if (starti >= 0 && endi >= 0) {
	  debug13(printf("starti = %d, endi = %d\n",starti,endi));
	  assert(starti >= endi);
	  for (i = starti; i >= endi; i--) {
	    debug13(printf("diagonal %u (%llu), querypos %d..%d\n",
			   (Chrpos_T) (plus_segments[i].diagonal - chroffset),(unsigned long long) plus_segments[i].diagonal,
			   plus_segments[i].querypos5,plus_segments[i].querypos3));
	    if (plus_segments[i].querypos5 >= STAGE2_MIN_OLIGO + index1interval) {
	      /* Case 3. Missing start of query, so there could be a middle splice */
	      debug13b(printf("  querypos5 %d >= %d + %d, so using this diagonal plus overall_max_distance\n",
			      plus_segments[i].querypos5,STAGE2_MIN_OLIGO,index1interval));
#ifdef USE_GREEDY
	      if ((mappingpos = subtract_bounded(plus_segments[i].diagonal,querylength + shortsplicedist_novelend,chroffset)) > middle_mappingstart_greedy &&
		  mappingpos < genomicbound) {
		middle_mappingstart_greedy = mappingpos;
		middle_mappingstart_p = true;
		debug13(printf("  Redefining middle mappingstart greedy to %u\n",middle_mappingstart_greedy - chroffset));
	      }
#endif
#ifdef LONG_ENDSPLICES
	      if ((mappingpos = subtract_bounded(plus_segments[i].diagonal,querylength + overall_max_distance,chroffset)) < middle_mappingstart_last) {
		/* Use < for NOT_GREEDY */
		middle_mappingstart_last = mappingpos;
		middle_mappingstart_p = true;
		debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	      }
#else
	      if ((mappingpos = subtract_bounded(plus_segments[i].diagonal,querylength,chroffset)) < middle_mappingstart_last) {
		/* Use < for NOT_GREEDY */
		middle_mappingstart_last = mappingpos;
		middle_mappingstart_p = true;
		debug13(printf("  Redefining middle mappingstart last to %u\n",middle_mappingstart_last - chroffset));
	      }
#endif

	    } else {
	      debug13b(printf("  querypos5 %d < %d + %d, so using this diagonal\n",
			      plus_segments[i].querypos5,STAGE2_MIN_OLIGO,index1interval));
#ifdef USE_GREEDY
	      if ((mappingpos = subtract_bounded(plus_segments[i].diagonal,querylength,chroffset)) > close_mappingstart_greedy &&
		  mappingpos < genomicbound) {
		close_mappingstart_greedy = mappingpos;
		close_mappingstart_p = true;
		debug13(printf("  Redefining close mappingstart greedy to %u\n",close_mappingstart_greedy - chroffset));
	      }
#endif
	      if ((mappingpos = subtract_bounded(plus_segments[i].diagonal,querylength,chroffset)) < close_mappingstart_last) {
		/* Use < for NOT_GREEDY */
		close_mappingstart_last = mappingpos;
		close_mappingstart_p = true;
		debug13(printf("  Redefining close mappingstart last to %u\n",close_mappingstart_last - chroffset));
	      }
	    }
	  }

#ifdef USE_GREEDY
	  if (close_mappingstart_p == true) {
	    close_knownsplice_limit_low = subtract_bounded(close_mappingstart_greedy,overall_max_distance,chroffset);
	  } else if (middle_mappingstart_p == true) {
	    debug13(printf("Using middle mappingstart\n"));
	    close_knownsplice_limit_low = middle_mappingstart_greedy;
	    close_mappingstart_greedy = middle_mappingstart_greedy;
	    close_mappingstart_p = true;
	  }
#else
	  if (close_mappingstart_p == true) {
	    close_knownsplice_limit_low = subtract_bounded(close_mappingstart_last,shortsplicedist,chroffset);
	  } else if (middle_mappingstart_p == true) {
	    debug13(printf("Using middle mappingstart\n"));
	    close_knownsplice_limit_low = middle_mappingstart_last;
	    close_mappingstart_last = middle_mappingstart_last;
	    close_mappingstart_p = true;
	  }
#endif
#ifdef USE_GREEDY
	  if (middle_mappingstart_p == true && middle_mappingstart_last < close_mappingstart_greedy) {
	    knownsplice_limit_low = middle_mappingstart_last;
	    mappingstart = middle_mappingstart_last;
	  } else if (close_mappingstart_p == true && close_mappingstart_last != close_mappingstart_greedy) {
	    knownsplice_limit_low = subtract_bounded(close_mappingstart_last,overall_max_distance,chroffset);
	    mappingstart = close_mappingstart_last;
	  }
#else
	  if (middle_mappingstart_p == true && middle_mappingstart_last < close_mappingstart_last) {
	    knownsplice_limit_low = middle_mappingstart_last;
	    mappingstart = middle_mappingstart_last;
	  }
#endif
	  if (close_mappingstart_p == false) {
	    fallback_mappingstart_p = false;
#ifdef USE_GREEDY
	  } else if (mappingstart >= close_mappingstart_greedy) {
	    fallback_mappingstart_p = false;
#endif
	  } else {
	    debug13(printf("Fallback mappingstart = %u\n",mappingstart - chroffset));
	    fallback_mappingstart_p = true;
	  }
	}
      }

      favor_right_p = true;

    } else {
      chroffset = Stage3end_chroffset(hit3);
      chrhigh = Stage3end_chrhigh(hit3);
      chrlength = Stage3end_chrlength(hit3);

      if (Shortread_find_primers(queryseq5,queryseq3) == true) {
	/* Go from genomicend */
	debug13(printf("Found primers\n"));
	genomicbound = Stage3end_genomicend(hit3);

      } else if (Stage3end_anomalous_splice_p(hit3) == true) {
	/* Go from genomicend */
	debug13(printf("Anomalous splice\n"));
	genomicbound = Stage3end_genomicend(hit3);

      } else {
	genomicbound = Stage3end_genomicstart(hit3);

#if 0
	/* TODO: Previously called Shortread_find_overlap.  Now with Shortread_max_overlap, can optimize this code */
	if ((overlap = Shortread_max_overlap(queryseq5,queryseq3)) > 0 &&
	    Stage3end_genomicbound_from_start(&genomicbound2,hit3,overlap,chroffset) == true) {
	  debug13(printf("Found overlap of %d\n",overlap));
	  if (genomicbound2 < genomicbound) {
	    zero_offset = genomicbound - genomicbound2;
	    genomicbound = genomicbound2;
	  }
	}
#endif
      }

      debug13(printf("Case 4: hit3 minus %s %u..%u (sensedir %d) => genomicbound %u\n",
		     Stage3end_hittype_string(hit3),
		     Stage3end_genomicstart(hit3) - chroffset,Stage3end_genomicend(hit3) - chroffset,
		     Stage3end_sensedir(hit3),genomicbound - chroffset));

      knownsplice_limit_low = mappingstart = segmentstart = genomicbound;
      knownsplice_limit_high = add_bounded(Stage3end_genomicstart(hit3),pairmax + overall_max_distance,chrhigh);
      segmentend = add_bounded(Stage3end_genomicstart(hit3),pairmax + PAIRMAX_ADDITIONAL,chrhigh);
#ifdef LONG_ENDSPLICES
      mappingend = add_bounded(Stage3end_genomicstart(hit3),pairmax + overall_max_distance,chrhigh);
#else
      mappingend = add_bounded(Stage3end_genomicstart(hit3),pairmax + overall_max_distance_novelend,chrhigh);
#endif

      close_mappingend_last = middle_mappingend_last = Stage3end_genomicstart(hit3);
#ifdef USE_GREEDY
      close_mappingend_greedy = middle_mappingend_greedy = segmentend;
#endif

      if (minus_nsegments > 0) {
	/* Use segments to bound */
	debug13(printf("Finding segments from segmentstart %u to segmentend %u (minus_nsegments %d)\n",
		       segmentstart - chroffset,segmentend - chroffset,minus_nsegments));
	starti = endi = -1;
	i = binary_search_segments(0,minus_nsegments-1,minus_segments,segmentstart);
	while (i < minus_nsegments - 1 && minus_segments[i].diagonal == (Univcoord_T) -1) {
	  i++;
	}
	starti = i;
	while (minus_segments[i].diagonal < segmentend) {
	  endi = i;
	  i++;
	}
	if (starti >= 0 && endi >= 0) {
	  debug13(printf("starti = %d, endi = %d\n",starti,endi));
	  assert(starti <= endi);
	  for (i = starti; i <= endi; i++) {
	    debug13(printf("diagonal %u (%llu), querypos %d..%d\n",
			   (Chrpos_T) (minus_segments[i].diagonal - chroffset),(unsigned long long) minus_segments[i].diagonal,
			   minus_segments[i].querypos5,minus_segments[i].querypos3));
	    if (minus_segments[i].querypos5 >= STAGE2_MIN_OLIGO + index1interval) {
	      /* Case 4. Missing start of query, so there could be a middle splice */
	      debug13b(printf("  querypos5 %d >= %d + %d, so using this diagonal plus overall_max_distance\n",
			      minus_segments[i].querypos5,STAGE2_MIN_OLIGO,index1interval));
#ifdef USE_GREEDY
	      if ((mappingpos = add_bounded(minus_segments[i].diagonal,shortsplicedist_novelend,chrhigh)) < middle_mappingend_greedy &&
		  mappingpos > genomicbound) {
		middle_mappingend_greedy = mappingpos;
		middle_mappingend_p = true;
		debug13(printf("  Redefining middle mappingend greedy to %u\n",middle_mappingend_greedy - chroffset));
	      }
#endif
#ifdef LONG_ENDSPLICES
	      if ((mappingpos = add_bounded(minus_segments[i].diagonal,overall_max_distance,chrhigh)) > middle_mappingend_last) {
		/* Use > for NOT_GREEDY */
		middle_mappingend_last = mappingpos;
		middle_mappingend_p = true;
		debug13(printf("  Redefining middle mappingend to %u\n",middle_mappingend_last - chroffset));
	      }
#else
	      if ((mappingpos = minus_segments[i].diagonal) > middle_mappingend_last) {
		/* Use > for NOT_GREEDY */
		middle_mappingend_last = mappingpos;
		middle_mappingend_p = true;
		debug13(printf("  Redefining middle mappingend to %u\n",middle_mappingend_last - chroffset));
	      }
#endif

	    } else {
	      debug13b(printf("  querypos5 %d < %d + %d, so using this diagonal\n",
			      minus_segments[i].querypos5,STAGE2_MIN_OLIGO,index1interval));
#ifdef USE_GREEDY
	      if ((mappingpos = minus_segments[i].diagonal) < close_mappingend_greedy &&
		  mappingpos > genomicbound) {
		close_mappingend_greedy = mappingpos;
		close_mappingend_p = true;
		debug13(printf("  Redefining close mappingend greedy to %u\n",close_mappingend_greedy - chroffset));
	      }
#endif
	      if ((mappingpos = minus_segments[i].diagonal) > close_mappingend_last) {
		/* Use > for NOT_GREEDY */
		close_mappingend_last = mappingpos;
		close_mappingend_p = true;
		debug13(printf("  Redefining close mappingend last to %u\n",close_mappingend_last - chroffset));
	      }
	    }
	  }

#ifdef USE_GREEDY
	  if (close_mappingend_p == true) {
	    close_knownsplice_limit_high = add_bounded(close_mappingend_greedy,overall_max_distance,chrhigh);
	  } else if (middle_mappingend_p == true) {
	    debug13(printf("Using middle mappingend\n"));
	    close_knownsplice_limit_high = middle_mappingend_greedy;
	    close_mappingend_greedy = middle_mappingend_greedy;
	    close_mappingend_p = true;
	  }
#else
	  if (close_mappingend_p == true) {
	    close_knownsplice_limit_high = add_bounded(close_mappingend_last,shortsplicedist,chrhigh);
	  } else if (middle_mappingend_p == true) {
	    debug13(printf("Using middle mappingend\n"));
	    close_knownsplice_limit_high = middle_mappingend_last;
	    close_mappingend_last = middle_mappingend_last;
	    close_mappingend_p = true;
	  }
#endif
#ifdef USE_GREEDY
	  if (middle_mappingend_p == true && middle_mappingend_last > close_mappingend_greedy) {
	    knownsplice_limit_high = middle_mappingend_last;
	    mappingend = middle_mappingend_last;
	  } else if (close_mappingend_p == true && close_mappingend_last != close_mappingend_greedy) {
	    knownsplice_limit_high = add_bounded(close_mappingend_last,overall_max_distance,chrhigh);
	    mappingend = close_mappingend_last;
	  }
#else
	  if (middle_mappingend_p == true && middle_mappingend_last > close_mappingend_last) {
	    knownsplice_limit_high = middle_mappingend_last;
	    mappingend = middle_mappingend_last;
	  }
#endif
	  if (close_mappingend_p == false) {
	    fallback_mappingend_p = false;
#ifdef USE_GREEDY
	  } else if (mappingend <= close_mappingend_greedy) {
	    fallback_mappingend_p = false;
#endif
	  } else {
	    debug13(printf("Fallback mappingend = %u\n",mappingend - chroffset));
	    fallback_mappingend_p = true;
	  }
	}
      }

      favor_right_p = true;
    }

    if ((sensedir = Stage3end_sensedir(hit3)) == SENSE_FORWARD) {
      sense_try = +1;
    } else if (sensedir == SENSE_ANTI) {
      sense_try = -1;
    } else {
      sense_try = 0;
    }

  } else {
    abort();
  }

#ifdef OLD_GENOMICBOUND
  knownsplice_limit_low = genomicstart + querylength;
  knownsplice_limit_high = genomicend - querylength;
#endif

  if (close_mappingstart_p == true && close_mappingend_p == true) {
    debug13(printf("Halfmapping: Running gmap with close mappingstart and close mappingend\n"));
    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
			       hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
			       close_mappingstart_last,close_mappingend_last,
			       close_knownsplice_limit_low,close_knownsplice_limit_high,
			       watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
			       oligoindices_major,oligoindices_minor,
			       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);

    if (good_start_p == true && good_end_p == true) {
      /* Success */
    } else if (gmap_rerun_p == false) {
      debug13(printf("Skipping re-run of gmap\n"));
    } else if (/* require both ends to be good */ 0 && good_start_p == true) {
      if (fallback_mappingend_p == true) {
	debug13(printf("Halfmapping: Re-running gmap with close mappingstart only\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
				   hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
				   close_mappingstart_last,mappingend,
				   close_knownsplice_limit_low,knownsplice_limit_high,
				   watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
      }

    } else if (/* require both ends to be good */ 0 && good_end_p == true) {
      if (fallback_mappingstart_p == true) {
	debug13(printf("Halfmapping: Re-running gmap with close mappingend only\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
				   hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
				   mappingstart,close_mappingend_last,
				   knownsplice_limit_low,close_knownsplice_limit_high,
				   watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
      }
    } else {
      if (fallback_mappingstart_p == true && fallback_mappingend_p == true) {
	debug13(printf("Halfmapping: Re-running gmap with far mappingstart and mappingend\n"));
	hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
				   hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
				   mappingstart,mappingend,
				   knownsplice_limit_low,knownsplice_limit_high,
				   watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
				   oligoindices_major,oligoindices_minor,
				   pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
      }
    }

  } else if (close_mappingstart_p == true) {
    debug13(printf("Halfmapping: Running gmap with close mappingstart\n"));
    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
			       hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
			       close_mappingstart_last,mappingend,
			       close_knownsplice_limit_low,knownsplice_limit_high,
			       watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
			       oligoindices_major,oligoindices_minor,
			       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);

    if (good_start_p == true && /* require both ends to be good */ good_end_p == true) {
      /* Success */
    } else if (gmap_rerun_p == false) {
      debug13(printf("Skipping re-run of gmap\n"));
    } else if (fallback_mappingstart_p == true) {
      debug13(printf("Halfmapping: Re-running gmap with far mappingstart\n"));
      hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
				 hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
				 mappingstart,mappingend,
				 knownsplice_limit_low,knownsplice_limit_high,
				 watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
				 oligoindices_major,oligoindices_minor,
				 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
    }

  } else if (close_mappingend_p == true) {
    debug13(printf("Halfmapping: Running gmap with close mappingend\n"));
    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
			       hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
			       mappingstart,close_mappingend_last,
			       knownsplice_limit_low,close_knownsplice_limit_high,
			       watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
			       oligoindices_major,oligoindices_minor,
			       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);

    if (good_end_p == true && /* require both ends to be good */ good_start_p == true) {
      /* Success */
    } else if (gmap_rerun_p == false) {
      debug13(printf("Skipping re-run of gmap\n"));
    } else if (fallback_mappingend_p == true) {
      debug13(printf("Halfmapping: Re-running gmap with far mappingend\n"));
      hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
				 hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
				 mappingstart,mappingend,
				 knownsplice_limit_low,knownsplice_limit_high,
				 watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
				 oligoindices_major,oligoindices_minor,
				 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
    }

  } else {
    debug13(printf("Halfmapping: Running gmap with far mappingstart and mappingend\n"));
    hits = run_gmap_for_region(&good_start_p,&good_end_p,gmap_history,
			       hits,Shortread_accession(queryseq5),queryuc_ptr,querylength,sense_try,favor_right_p,
			       mappingstart,mappingend,
			       knownsplice_limit_low,knownsplice_limit_high,
			       watsonp,genestrand,chrnum,chroffset,chrhigh,chrlength,
			       oligoindices_major,oligoindices_minor,
			       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
  }

  return hits;
}


static List_T
align_pair_with_gmap (Pairtype_T *final_pairtype, List_T result,
		      char *queryuc_ptr_5, char *queryuc_ptr_3,
		      int querylength5, int querylength3,
		      Oligoindex_array_T oligoindices_minor,
		      Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		      int cutoff_level_5, int cutoff_level_3, bool expect_concordant_p) {
  Stage3pair_T newpair, stage3pair;
  List_T gmap5_hits = NULL, gmap3_hits = NULL;
  Stage3end_T hit5, hit3, gmap5, gmap3, gmap5_1, gmap5_2, gmap3_1, gmap3_2;
  List_T p, a, b, rest;
  int genestrand;
  int i;
  bool replacedp;
#ifdef DEBUG13
  int missing_hit, missing_gmap;
#endif


  debug16(printf("Running align_pair_with GMAP on %d results\n",List_length(result)));
  debug13(printf("Sorting %d hitpairs by nmatches\n",List_length(result)));
  result = Stage3pair_sort_bymatches(result);

  for (p = result, i = 0; p != NULL && i < max_gmap_improvement; p = p->rest, i++) {
    stage3pair = (Stage3pair_T) List_head(p);
    genestrand = Stage3pair_genestrand(stage3pair);
    hit5 = Stage3pair_hit5(stage3pair);
    hit3 = Stage3pair_hit3(stage3pair);
    gmap5_1 = gmap5_2 = gmap3_1 = gmap3_2 = (Stage3end_T) NULL;

    debug13(printf("GMAP improvement #%d: Entering align_pair_with_gmap with hittypes %s and %s\n",
		   i,Stage3end_hittype_string(hit5),Stage3end_hittype_string(hit3)));

    /* Was querylength5 - Stage3end_matches(hit5) > 5 */
    debug13(printf("Looking at hit5 %p with nmismatches %d - %d ?<= cutoff_level %d\n",
		   hit5,querylength5,Stage3end_nmatches_posttrim(hit5),cutoff_level_5));

#if 0
    if (Stage3end_sarrayp(hit5) == true && redo_for_sense_p == false) {
      /* Skip */
      debug13(printf("Skipping hit5 from sarray search\n"));

    } else if (Stage3end_hittype(hit5) == GMAP && redo_for_sense_p == false) {
      /* Skip */
      debug13(printf("Skipping hit5 of type GMAP\n"));

      /* Don't skip on final align_concordant_with_gmap */
    } else if (Stage3end_hittype(hit5) == TERMINAL) {
      /* Skip */
      debug13(printf("Skipping hit5 of type TERMINAL\n"));

    } /* else  */
#endif

    if (Stage3end_hittype(hit5) == GMAP) {
      /* Skip */
      debug13(printf("Skipping hit5 already of type GMAP\n"));

    } else if (querylength5 - Stage3end_nmatches_posttrim(hit5) <= cutoff_level_5) {
      /* Skip, because already good enough */
      debug13(printf("Skipping hit5 with nmismatches %d - %d <= cutoff_level %d\n",
		     querylength5,Stage3end_nmatches_posttrim(hit5),cutoff_level_5));

    } else {
      debug13(printf("Aligning hit5 with GMAP\n"));
      align_single_hit_with_gmap(&gmap5_1,&gmap5_2,hit5,queryuc_ptr_5,querylength5,
#ifdef END_KNOWNSPLICING_SHORTCUT
				 queryrc5,Shortread_invertedp(queryseq5),
#endif
				 oligoindices_minor,
				 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				 genestrand,/*extend_ends_p*/true);  /* if false, then causes major slowdown */
      if (gmap5_1 != NULL) {
	debug13(missing_hit = querylength5 - Stage3end_nmatches_posttrim(hit5));
	debug13(missing_gmap = querylength5 - Stage3end_nmatches_posttrim(gmap5_1));
	debug13(printf("GMAP %p with %d matches, %d missing compared with original 5' hit with %d matches, %d missing\n",
		       gmap5_1,Stage3end_nmatches_posttrim(gmap5_1),missing_gmap,Stage3end_nmatches_posttrim(hit5),missing_hit));
	gmap5_hits = List_push(gmap5_hits,(void *) gmap5_1);
	if (gmap5_2 != NULL) {
	  gmap5_hits = List_push(gmap5_hits,(void *) gmap5_2);
	}
	Stage3end_set_improved_by_gmap(hit5);
      }
    }

    debug13(printf("Looking at hit3 %p with nmismatches %d - %d ?<= cutoff_level %d\n",
		   hit3,querylength3,Stage3end_nmatches_posttrim(hit3),cutoff_level_3));

    if (Stage3end_hittype(hit3) == GMAP) {
      /* Skip */
      debug13(printf("Skipping hit3 already of type GMAP\n"));

    } else if (querylength3 - Stage3end_nmatches_posttrim(hit3) <= cutoff_level_3) {
      /* Skip, because already good enough */
      debug13(printf("Skipping hit3 with nmismatches %d - %d <= cutoff_level %d\n",
		     querylength3,Stage3end_nmatches_posttrim(hit3),cutoff_level_3));

    } else {
      debug13(printf("Aligning hit3 with GMAP\n"));
      align_single_hit_with_gmap(&gmap3_1,&gmap3_2,hit3,queryuc_ptr_3,querylength3,
#ifdef END_KNOWNSPLICING_SHORTCUT
				 queryrc3,Shortread_invertedp(queryseq3),
#endif
				 oligoindices_minor,
				 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				 genestrand,/*extend_ends_p*/true);  /* if false, then causes major slowdown */
      if (gmap3_1 != NULL) {
	debug13(missing_hit = querylength3 - Stage3end_nmatches_posttrim(hit3));
	debug13(missing_gmap = querylength3 - Stage3end_nmatches_posttrim(gmap3_1));
	debug13(printf("GMAP %p with %d matches, %d missing compared with original 3' hit with %d matches, %d missing\n",
		       gmap3_1,Stage3end_nmatches_posttrim(gmap3_1),missing_gmap,Stage3end_nmatches_posttrim(hit3),missing_hit));
	gmap3_hits = List_push(gmap3_hits,(void *) gmap3_1);
	if (gmap3_2 != NULL) {
	  gmap3_hits = List_push(gmap3_hits,(void *) gmap3_2);
	}
	Stage3end_set_improved_by_gmap(hit3);
      }
    }

    if (gmap5_hits != NULL && gmap3_hits != NULL) {
      replacedp = false;
      for (a = gmap5_hits; a != NULL; a = List_next(a)) {
	gmap5 = (Stage3end_T) List_head(a);

	for (b = gmap3_hits; b != NULL; b = List_next(b)) {
	  gmap3 = (Stage3end_T) List_head(b);

	  debug13(printf("Imperfect concordant uniq: Double GMAP on hit5 and hit3"));
	  if ((newpair = Stage3pair_new(Stage3end_copy(gmap5),Stage3end_copy(gmap3),genestrand,/*pairtype*/UNSPECIFIED,
					/*private5p*/true,/*private3p*/true,expect_concordant_p)) == NULL) {
	    /* Stage3end_free(&gmap3); -- done by Stage3pair_new */
	    /* Stage3end_free(&gmap5); -- done by Stage3pair_new */
	    debug13(printf("  => NULL, so eliminating\n"));

	  } else if (Stage3pair_determine_pairtype(newpair) != CONCORDANT) {
	    debug13(printf("  => not concordant, so eliminating\n"));

	  } else if (replacedp == false) {
	    /* Convert to gmap-gmap */
	    debug13(printf("  => replacement\n"));
	    *final_pairtype = CONCORDANT;
	    List_head_set(p,(void *) newpair);
	    replacedp = true;

	  } else {
	    debug13(printf("  => addition\n"));
	    *final_pairtype = CONCORDANT;
	    rest = List_push(List_next(p),(void *) newpair);
	    List_tail_set(p,rest);
	    p = rest;
	  }
	}
      }

      if (replacedp == true) {
	Stage3pair_free(&stage3pair); /* Also frees hit5 and hit3 */
      }
      for (a = gmap5_hits; a != NULL; a = List_next(a)) {
	gmap5 = (Stage3end_T) List_head(a);
	Stage3end_free(&gmap5);
      }
      for (b = gmap3_hits; b != NULL; b = List_next(b)) {
	gmap3 = (Stage3end_T) List_head(b);
	Stage3end_free(&gmap3);
      }

      List_free(&gmap3_hits);
      List_free(&gmap5_hits);

    } else {
      debug13(printf("Have %d GMAP 5' hits and %d GMAP 3' hits\n",
		     List_length(gmap5_hits),List_length(gmap3_hits)));

      /* Handle gmap5 hits */
      replacedp = false;
      for (a = gmap5_hits; a != NULL; a = List_next(a)) {
	gmap5 = (Stage3end_T) List_head(a);

	debug13(printf("Imperfect concordant uniq: Single GMAP on hit5"));
	if ((newpair = Stage3pair_new(gmap5,Stage3end_copy(hit3),genestrand,/*pairtype*/UNSPECIFIED,
				      /*private5p*/true,/*private3p*/true,expect_concordant_p)) == NULL) {
	  /* Stage3end_free(&gmap5); -- done by Stage3pair_new */
	  debug13(printf(" => NULL\n"));

	} else if (Stage3pair_determine_pairtype(newpair) != CONCORDANT) {
	  debug13(printf("  => not concordant, so eliminating\n"));

	} else if (replacedp == false) {
	  /* Convert to gmap-xx */
	  debug13(printf(" => replacement\n"));
	  *final_pairtype = CONCORDANT;
	  List_head_set(p,(void *) newpair);
	  replacedp = true;

	} else {
	  debug13(printf(" => addition\n"));
	  *final_pairtype = CONCORDANT;
	  rest = List_push(List_next(p),(void *) newpair);
	  List_tail_set(p,rest);
	  p = rest;
	}
      }

      if (replacedp == true) {
	Stage3pair_free(&stage3pair);
      }
      /* Do not free gmap5 objects, since not copied */
      List_free(&gmap5_hits);


      /* Handle gmap3 hits */
      replacedp = false;
      for (b = gmap3_hits; b != NULL; b = List_next(b)) {
	gmap3 = (Stage3end_T) List_head(b);

	debug13(printf("Imperfect concordant uniq: Single GMAP on hit3"));
	if ((newpair = Stage3pair_new(Stage3end_copy(hit5),gmap3,genestrand,/*pairtype*/UNSPECIFIED,
				      /*private5p*/true,/*private3p*/true,expect_concordant_p)) == NULL) {
	  /* Stage3end_free(&gmap3); -- done by Stage3pair_new */
	  debug13(printf(" => NULL\n"));

	} else if (Stage3pair_determine_pairtype(newpair) != CONCORDANT) {
	  debug13(printf("  => not concordant, so eliminating\n"));

	} else if (replacedp == false) {
	  /* Convert to xx-gmap */
	  debug13(printf(" => replacement\n"));
	  *final_pairtype = CONCORDANT;
	  List_head_set(p,(void *) newpair);
	  replacedp = true;

	} else {
	  debug13(printf(" => addition\n"));
	  *final_pairtype = CONCORDANT;
	  rest = List_push(List_next(p),(void *) newpair);
	  List_tail_set(p,rest);
	  p = rest;
	}
      }

      if (replacedp == true) {
	Stage3pair_free(&stage3pair);
      }
      /* Do not free gmap3 objects, since not copied */
      List_free(&gmap3_hits);

    }
  }

  debug13(printf("End of align_pair_with_gmap with %d results\n",List_length(result)));

  return result;
}


/* Need to have this to resolve asymmetry between plus and minus
   searches for suffix array.  This will invoke deeper methods when
   necessary. */
static bool
better_free_end_exists_p (List_T greedy, List_T subs, List_T terminals,
			  List_T indels, List_T singlesplicing, List_T doublesplicing,
			  int querylength) {
  int best_concordant_score = querylength, score;

  /* SPEED */
  return false;

  if ((score = Stage3end_best_score_paired(greedy)) < best_concordant_score) {
    best_concordant_score = score;
  }
  if ((score = Stage3end_best_score_paired(subs)) < best_concordant_score) {
    best_concordant_score = score;
  }
  if ((score = Stage3end_best_score_paired(terminals)) < best_concordant_score) {
    best_concordant_score = score;
  }
  if ((score = Stage3end_best_score_paired(indels)) < best_concordant_score) {
    best_concordant_score = score;
  }
  if ((score = Stage3end_best_score_paired(singlesplicing)) < best_concordant_score) {
    best_concordant_score = score;
  }
  if ((score = Stage3end_best_score_paired(doublesplicing)) < best_concordant_score) {
    best_concordant_score = score;
  }
  debug(printf("Best concordant score = %d\n",best_concordant_score));

  if (Stage3end_equiv_score_unpaired_p(greedy,best_concordant_score) == true) {
    debug(printf("Better or equivalent score found in greedy\n"));
    return true;
  } else if (Stage3end_equiv_score_unpaired_p(subs,best_concordant_score) == true) {
    debug(printf("Better or equivalent score found in subs\n"));
    return true;
  } else if (Stage3end_equiv_score_unpaired_p(terminals,best_concordant_score) == true) {
    debug(printf("Better or equivalent score found in terminals\n"));
    return true;
  } else if (Stage3end_equiv_score_unpaired_p(indels,best_concordant_score) == true) {
    debug(printf("Better or equivalent score found in indels\n"));
    return true;
  } else if (Stage3end_equiv_score_unpaired_p(singlesplicing,best_concordant_score) == true) {
    debug(printf("Better or equivalent score found in singlesplicing\n"));
    return true;
  } else if (Stage3end_equiv_score_unpaired_p(doublesplicing,best_concordant_score) == true) {
    debug(printf("Better or equivalent score found in doublesplicing\n"));
    return true;
  } else {
    return false;
  }
}



/* Search order for paired-end reads:

   1.  suffix array
   2.  exact/subs, via spanning set algorithm
   3.  subs/indels, via complete set algorithm
   4.  segments -> single splicing
   5.  segments -> double splicing (currently disabled)

   6.  paired segments -> GMAP via segments
   7.  distant splicing (needs to be before terminals, or we won't find them)
   8.  terminals
   
   9.  if still no concordance: GMAP pairsearch

   in caller: consolidate: does GMAP via substrings
*/


#define HITARRAY_GREEDY 0
#define HITARRAY_SUBS 1
#define HITARRAY_INDELS 2
#define HITARRAY_SINGLESPLICING 3
#define HITARRAY_DOUBLESPLICING 4

#if 0
#define HITARRAY_LONGSINGLESPLICING 7
#define HITARRAY_DISTANTSPLICING 8
#define HITARRAY_SEGMENTS_GMAP 9
#define HITARRAY_N 10
#else
#define HITARRAY_N 5
#endif

static List_T
align_pair (bool *abort_pairing_p, int *found_score, int *cutoff_level_5, int *cutoff_level_3,
	    List_T *samechr, List_T *conc_transloc, List_T *terminals,
	    History_T gmap_history_5, History_T gmap_history_3, List_T *hits5, List_T *hits3, T this5, T this3,
	    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
	    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
	    char *queryuc_ptr_5, char *queryuc_ptr_3,
#ifndef LARGE_GENOMES
	    char *queryrc5, char *queryrc3,
#endif
	    int querylength5, int querylength3, int query5_lastpos, int query3_lastpos,
	    Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev, int indexdb_size_threshold, Floors_T *floors_array,

	    Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
	    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
	    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,

	    int user_maxlevel_5, int user_maxlevel_3, int indel_penalty_middle, int indel_penalty_end,
	    int localsplicing_penalty, int distantsplicing_penalty,
	    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
	    Chrpos_T pairmax, int maxpairedpaths, bool keep_floors_p, Shortread_T queryseq5, Shortread_T queryseq3,
	    int genestrand) {

  List_T hitpairs = NULL, p;
  Stage3pair_T newpair;
  List_T halfmapping5, halfmapping3, a;
  Stage3end_T hit5, hit3, gmap5, gmap3;
  List_T hitarray5[HITARRAY_N], hitarray3[HITARRAY_N];
  Segment_T *plus_anchor_segments_5 = NULL, *minus_anchor_segments_5 = NULL,
    *plus_anchor_segments_3 = NULL, *minus_anchor_segments_3 = NULL;
  int n_plus_anchors_5 = 0, n_minus_anchors_5 = 0, n_plus_anchors_3 = 0, n_minus_anchors_3 = 0;
  List_T greedy5 = NULL, subs5 = NULL, terminals5 = NULL,
    indels5 = NULL, singlesplicing5 = NULL, doublesplicing5 = NULL,
    distantsplicing5 = NULL, gmap5_hits = NULL;
  List_T greedy3 = NULL, subs3 = NULL, terminals3 = NULL,
    indels3 = NULL, singlesplicing3 = NULL, doublesplicing3 = NULL,
    distantsplicing3 = NULL, gmap3_hits = NULL;
  List_T longsinglesplicing5 = NULL, longsinglesplicing3 = NULL;
  int ignore_found_score, done_level_5, done_level_3, opt_level, fast_level_5, fast_level_3,
    mismatch_level_5, mismatch_level_3, nmismatches;
  int max_splice_mismatches_5 = -1, max_splice_mismatches_3 = -1, i;
  int nhits5 = 0, nhits3 = 0, nsplicepairs5 = 0, nsplicepairs3 = 0;
  List_T *donors_plus_5, *antidonors_plus_5, *acceptors_plus_5, *antiacceptors_plus_5,
    *donors_minus_5, *antidonors_minus_5, *acceptors_minus_5, *antiacceptors_minus_5;
  List_T *donors_plus_3, *antidonors_plus_3, *acceptors_plus_3, *antiacceptors_plus_3,
    *donors_minus_3, *antidonors_minus_3, *acceptors_minus_3, *antiacceptors_minus_3;

  bool spanningset5p, spanningset3p, completeset5p, completeset3p, gmap5p, gmap3p;
  bool did_alignment_p, did_singlesplicing5_p, did_singlesplicing3_p;
  bool any_omitted_p_5, any_omitted_p_3;
  Floors_T floors5, floors3;
  bool alloc_floors_p_5 = false, alloc_floors_p_3 = false, floors5_computed_p = false, floors3_computed_p = false,
    segments5_computed_p = false, segments3_computed_p = false;
  int best_score_paired;
  /* bool found_terminals_p = false; */
  int nconcordant = 0, nsamechr = 0;
  Indexdb_T plus_indexdb_5, plus_indexdb_3, minus_indexdb_5, minus_indexdb_3;
  bool allvalidp5, allvalidp3;
  int min_trim, trim;
#ifndef LARGE_GENOMES
  int nmisses_allowed_sarray_5, nmisses_allowed_sarray_3;
#endif


  if (genestrand == +2) {
    plus_indexdb_5 = indexdb_rev;
    plus_indexdb_3 = indexdb_rev;
    minus_indexdb_5 = indexdb_fwd;
    minus_indexdb_3 = indexdb_fwd;
  } else {
    plus_indexdb_5 = indexdb_fwd;
    plus_indexdb_3 = indexdb_fwd;
    minus_indexdb_5 = indexdb_rev;
    minus_indexdb_3 = indexdb_rev;
  }

  *samechr = (List_T) NULL;
  *conc_transloc = (List_T) NULL;
  *terminals = (List_T) NULL;
  *abort_pairing_p = false;

  /* For paired-end alignment, ignore found_scores from single-end
     alignments.  Use only the found_score from
     Stage3_pair_up_concordant. */
  *found_score = querylength5 + querylength3;
  ignore_found_score = querylength5 + querylength3;

  if (querylength5 < min_kmer_readlength) {
    fast_level_5 = querylength5 - 1 - NREQUIRED_FAST;
    debug(printf("fast_level_5 %d = querylength %d - 1 - nrequired_fast %d\n",
		 fast_level_5,querylength5,NREQUIRED_FAST));
  } else {
    fast_level_5 = (querylength5 + index1interval - 1)/spansize - NREQUIRED_FAST;
    debug(printf("fast_level_5 %d = (querylength %d + index1interval %d - 1)/spansize %d - nrequired_fast %d\n",
		 fast_level_5,querylength5,index1interval,spansize,NREQUIRED_FAST));
  }

  if (querylength3 < min_kmer_readlength) {
    fast_level_3 = querylength3 - 1 - NREQUIRED_FAST;
    debug(printf("fast_level_3 %d = querylength %d - 1 - nrequired_fast %d\n",
		 fast_level_3,querylength3,NREQUIRED_FAST));
  } else {
    fast_level_3 = (querylength3 + index1interval - 1)/spansize - NREQUIRED_FAST;
    debug(printf("fast_level_3 %d = (querylength %d + index1interval %d - 1)/spansize %d - nrequired_fast %d\n",
		 fast_level_3,querylength3,index1interval,spansize,NREQUIRED_FAST));
  }

#if 0
  /* This prevents complete_mm procedure, needed for short reads */
  if (fast_level_5 < 1 && user_maxlevel_5 < 0) {
    fast_level_5 = 1;		/* Do at least 1 mismatch */
  }
  if (fast_level_3 < 1 && user_maxlevel_3 < 0) {
    fast_level_3 = 1;		/* Do at least 1 mismatch */
  }
#endif

  if (user_maxlevel_5 >= 0) {
    *cutoff_level_5 = user_maxlevel_5;
  } else if (fast_level_5 >= 0) {
    *cutoff_level_5 = fast_level_5;
  } else {
    *cutoff_level_5 = 0;
  }

  if (user_maxlevel_3 >= 0) {
    *cutoff_level_3 = user_maxlevel_3;
  } else if (fast_level_3 >= 0) {
    *cutoff_level_3 = fast_level_3;
  } else {
    *cutoff_level_3 = 0;
  }

  if (user_maxlevel_5 < 0) {
    if (fast_level_5 >= 0) {
      user_maxlevel_5 = fast_level_5;
    } else {
      user_maxlevel_5 = 0;
    }
  }

  if (user_maxlevel_3 < 0) {
    if (fast_level_3 >= 0) {
      user_maxlevel_3 = fast_level_3;
    } else {
      user_maxlevel_3 = 0;
    }
  }

#if 0
  if (dibasep) {
    opt_level = querylength5 + querylength3;
    done_level_5 = querylength5;
    done_level_3 = querylength3;
  }
#endif
  opt_level = user_maxlevel_5 + user_maxlevel_3;
  done_level_5 = user_maxlevel_5 /* + subopt_levels */;
  done_level_3 = user_maxlevel_3 /* + subopt_levels */;
  debug(printf("0> opt_level %d, done_level %d,%d\n",opt_level,done_level_5,done_level_3));

  for (i = 0; i < HITARRAY_N; i++) {
    hitarray5[i] = hitarray3[i] = (List_T) NULL;
  }

  nhits5 = nhits3 = 0;

#ifndef LARGE_GENOMES
  nmisses_allowed_sarray_5 = *cutoff_level_5;
  nmisses_allowed_sarray_3 = *cutoff_level_3;

  if (use_only_sarray_p == true) {
    *hits5 = Sarray_search_greedy(&(*cutoff_level_5),
				  queryuc_ptr_5,queryrc5,querylength5,query5_compress_fwd,query5_compress_rev,
				  nmisses_allowed_sarray_5,genestrand);

    *hits3 = Sarray_search_greedy(&(*cutoff_level_3),
				  queryuc_ptr_3,queryrc3,querylength3,query3_compress_fwd,query3_compress_rev,
				  nmisses_allowed_sarray_3,genestrand);

    /* Need to run Stage3end_remove_duplicates before we append the results together */
    hitarray5[HITARRAY_GREEDY] = *hits5;
    hitarray3[HITARRAY_GREEDY] = *hits3;
    debug(printf("sarray only: 5' end has %d greedy\n",List_length(*hits5)));
    debug(printf("sarray only: 3' end has %d greedy\n",List_length(*hits3)));

    if (*hits5 == NULL || *hits3 == NULL) {
      return (List_T) NULL;

    } else {
      hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					   &(*samechr),&(*conc_transloc),&(*terminals),
					   hitpairs,hitarray5,/*narray5*/HITARRAY_GREEDY+1,
					   hitarray3,/*narray3*/HITARRAY_GREEDY+1,
					   *cutoff_level_5,*cutoff_level_3,
					   querylength5,querylength3,maxpairedpaths,genestrand);
      debug(printf("After pairing sarray, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		   nconcordant,nsamechr,List_length(*terminals),*found_score));
      debug(printf("SA> found_score = %d, done_level %d,%d\n",*found_score,done_level_5,done_level_3));

      return Stage3pair_remove_circular_alias(hitpairs);
    }
  }
#endif

  /* Search 1: Suffix array */
#ifdef LARGE_GENOMES
  spanningset5p = spanningset3p = true;
#else
  if (use_sarray_p == false) {
    spanningset5p = spanningset3p = true;
  } else {
    spanningset5p = spanningset3p = false;	/* By default, suffix array search replaces spanning set */
    
    debug(printf("Trying suffix array on 5' end\n"));
    greedy5 = Sarray_search_greedy(&ignore_found_score,
				   queryuc_ptr_5,queryrc5,querylength5,query5_compress_fwd,query5_compress_rev,
				   nmisses_allowed_sarray_5,genestrand);

    debug(printf("Trying suffix array on 3' end\n"));
    greedy3 = Sarray_search_greedy(&ignore_found_score,
				   queryuc_ptr_3,queryrc3,querylength3,query3_compress_fwd,query3_compress_rev,
				   nmisses_allowed_sarray_3,genestrand);

    hitarray5[HITARRAY_GREEDY] = greedy5;
    hitarray3[HITARRAY_GREEDY] = greedy3;
    debug(printf("sarray initial: 5' end has %d greedy\n",List_length(greedy5)));
    debug(printf("sarray initial: 3' end has %d greedy\n",List_length(greedy3)));
    
    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					 &(*samechr),&(*conc_transloc),&(*terminals),
					 hitpairs,hitarray5,/*narray5*/HITARRAY_GREEDY+1,
					 hitarray3,/*narray3*/HITARRAY_GREEDY+1,
					 *cutoff_level_5,*cutoff_level_3,
					 querylength5,querylength3,maxpairedpaths,genestrand);
      
    debug(printf("After pairing sarray, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		 nconcordant,nsamechr,List_length(*terminals),*found_score));
    if (*abort_pairing_p == true) {
      *hits5 = greedy5;
      *hits3 = greedy3;
      hitpairs = Stage3pair_remove_circular_alias(hitpairs);
#if 0
      hitpairs = Stage3pair_remove_overlaps(hitpairs,/*translocp*/false,/*finalp*/true);
#endif
      return hitpairs;

    } else {
      opt_level = (*found_score < opt_level) ? *found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("SA> found_score = %d, opt_level %d, done_level %d,%d\n",*found_score,opt_level,done_level_5,done_level_3));
    }
    nhits5 = List_length(greedy5);
    nhits3 = List_length(greedy3);

    debug(printf("nconcordant %d\n",nconcordant));
    if (nconcordant == 0) {
      if (*terminals == NULL) {
	/* Need to have this to compensate for greediness of suffix array algorithm */
	debug(printf("nconcordant is 0, so we are doing spanningset on both sides\n"));
	spanningset5p = spanningset3p = true;
      } else {
	/* Have concordant terminals, so may need to refine just one end */
	for (p = *terminals; p != NULL; p = List_next(p)) {
	  newpair = (Stage3pair_T) List_head(p);
	  debug(printf("Terminal has total trims of %d and %d\n",
		       Stage3end_total_trim(Stage3pair_hit5(newpair)),
		       Stage3end_total_trim(Stage3pair_hit3(newpair))));
	  if (Stage3end_total_trim(Stage3pair_hit5(newpair)) > 15) {
	    spanningset5p = true;
	  }
	  if (Stage3end_total_trim(Stage3pair_hit3(newpair)) > 15) {
	    spanningset3p = true;
	  }
	  debug(printf("spanningset5p %d, spanningset3p %d\n",spanningset5p,spanningset3p));
	}
      }
#if 0
    } else if (*found_score >= done_level_5 + done_level_3) {
      /* Can eliminate this, because nconcordant no longer includes terminals, which are our main concern */
      debug(printf("found_score %d >= done_level_5 %d + done_level_3 %d,, so we are doing spanningset\n",
		   *found_score,done_level_5,done_level_3));
      spanningset5p = spanningset3p = true;
#endif
    }
  }
#endif

  if (querylength5 < min_kmer_readlength) {
    spanningset5p = false;
  }
  if (querylength3 < min_kmer_readlength) {
    spanningset3p = false;
  }

  /* Search 2: Exact/subs via spanning set algorithm */
  if (spanningset5p == true || spanningset3p == true) {
    /* 1A. Exact.  Requires compress if cmet or genomealt.  Creates and uses spanning set. */
    debug(printf("Performing spanning set with found_score %d\n",*found_score));

    mismatch_level_5 = 0;
    if (done_level_5 == 0 && snpp == false) {
      debug(printf("Suffix array already found exact matches for 5' end and no SNPs, so spanning set can't do any better\n"));
    } else {
      read_oligos(&allvalidp5,this5,queryuc_ptr_5,querylength5,query5_lastpos,genestrand);
      if (allvalidp5 == false) {
	debug(printf("Not all oligos in 5' end are valid, so cannot perform spanning set\n"));
	fast_level_5 = -1;
	spanningset5p = false;
      } else if (spanningset5p == true) {
	debug(printf("fast_level_5 = %d\n",fast_level_5));
	debug(printf("*** Stage 1.  Exact ***\n"));
	ignore_found_score = *found_score;
	subs5 = find_spanning_exact_matches(&ignore_found_score,&nhits5,subs5,this5,genestrand,
					    querylength5,query5_lastpos,plus_indexdb_5,minus_indexdb_5,
					    query5_compress_fwd,query5_compress_rev);
	mismatch_level_5 = 1;
      }
    }

    /* 1B. Exact.  Requires compress if cmet or genomealt.  Creates and uses spanning set. */
    mismatch_level_3 = 0;
    if (done_level_3 == 0 && snpp == false) {
      debug(printf("Suffix array already found exact matches for 3' end and no SNPs, so spanning set can't do any better\n"));
    } else {
      read_oligos(&allvalidp3,this3,queryuc_ptr_3,querylength3,query3_lastpos,genestrand);
      if (allvalidp3 == false) {
	debug(printf("Not all oligos in 3' end are valid, so cannot perform spanning set\n"));
	fast_level_3 = -1;
	spanningset3p = false;
      } else if (spanningset3p == true) {
	debug(printf("fast_level_3 = %d\n",fast_level_3));
	debug(printf("*** Stage 1.  Exact ***\n"));
	ignore_found_score = *found_score;
	subs3 = find_spanning_exact_matches(&ignore_found_score,&nhits3,subs3,this3,genestrand,
					    querylength3,query3_lastpos,plus_indexdb_3,minus_indexdb_3,
					    query3_compress_fwd,query3_compress_rev);
	mismatch_level_3 = 1;
      }
    }

    /* 1. Pairing after exact */
    /* Should not have duplicates from the spanning set procedure */
    hitarray5[HITARRAY_SUBS] = subs5; /* = Stage3end_remove_duplicates(subs5) */;
    hitarray3[HITARRAY_SUBS] = subs3; /* = Stage3end_remove_duplicates(subs3) */;
    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					 &(*samechr),&(*conc_transloc),&(*terminals),
					 hitpairs,hitarray5,/*narray5*/HITARRAY_SUBS+1,
					 hitarray3,/*narray3*/HITARRAY_SUBS+1,
					 *cutoff_level_5,*cutoff_level_3,
					 querylength5,querylength3,maxpairedpaths,genestrand);
    debug(printf("After pairing exact, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		 nconcordant,nsamechr,List_length(*terminals),*found_score));
    if (*abort_pairing_p == true) {
      *hits5 = List_append(greedy5,subs5);
      *hits3 = List_append(greedy3,subs3);
      return hitpairs;
    } else {
      opt_level = (*found_score < opt_level) ? *found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("1> found_score = %d, opt_level %d, done_level %d,%d\n",*found_score,opt_level,done_level_5,done_level_3));
    }

    did_alignment_p = false;

    /* 2A. One mismatch.  Requires spanning set and compress. */
    if (spanningset5p /*&& allvalidp5*/ && querylength5 >= one_miss_querylength && done_level_5 >= 1) {
      debug(printf("*** Stage 2A.  One miss ***\n"));
      did_alignment_p = true;
      ignore_found_score = *found_score;
      subs5 = find_spanning_onemiss_matches(&ignore_found_score,&nhits5,subs5,this5,genestrand,
					    querylength5,query5_compress_fwd,query5_compress_rev);
      mismatch_level_5 = 2;
    }

    /* 2B. One mismatch.  Requires spanning set and compress. */
    if (spanningset3p /*&& allvalidp3*/ && querylength3 >= one_miss_querylength && done_level_3 >= 1) {
      debug(printf("*** Stage 2B.  One miss ***\n"));
      did_alignment_p = true;
      ignore_found_score = *found_score;
      subs3 = find_spanning_onemiss_matches(&ignore_found_score,&nhits3,subs3,this3,genestrand,
					    querylength3,query3_compress_fwd,query3_compress_rev);
      mismatch_level_3 = 2;
    }

    if (did_alignment_p == true) {
      /* 2. Pairing after one mismatch */
      hitarray5[HITARRAY_SUBS] = subs5 /* = Stage3end_remove_duplicates(subs5,queryseq5,queryseq3) */;
      hitarray3[HITARRAY_SUBS] = subs3 /* = Stage3end_remove_duplicates(subs3,queryseq5,queryseq3) */;
      hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					   &(*samechr),&(*conc_transloc),&(*terminals),
					   hitpairs,hitarray5,/*narray5*/HITARRAY_SUBS+1,
					   hitarray3,/*narray3*/HITARRAY_SUBS+1,
					   *cutoff_level_5,*cutoff_level_3,
					   querylength5,querylength3,maxpairedpaths,genestrand);
      debug(printf("After pairing one mismatch, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		   nconcordant,nsamechr,List_length(*terminals),*found_score));
      if (*abort_pairing_p == true) {
	*hits5 = subs5;
	*hits3 = subs3;
	hitpairs = Stage3pair_remove_circular_alias(hitpairs);
#if 0
	hitpairs = Stage3pair_remove_overlaps(hitpairs,/*translocp*/false,/*finalp*/true);
#endif
	return hitpairs;
      } else {
	opt_level = (*found_score < opt_level) ? *found_score : opt_level;
	if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	  done_level_5 = user_maxlevel_5;
	}
	if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	  done_level_3 = user_maxlevel_3;
	}
	debug(printf("2> found_score = %d, opt_level %d, done_level %d,%d\n",*found_score,opt_level,done_level_5,done_level_3));
      }
    }


    did_alignment_p = false;

    /* 3A. Mismatches via spanning set.  Requires spanning set and compress. */
    if (spanningset5p /*&& allvalidp5*/ && done_level_5 >= 2) {
      /* NOTE: Since done_level isn't updated, can do in one batch instead of iteratively */
      while (mismatch_level_5 <= fast_level_5 && mismatch_level_5 <= done_level_5) {
	debug(printf("*** Stage 3A (level %d).  Spanning set mismatches ***\n",mismatch_level_5));
	did_alignment_p = true;
	ignore_found_score = *found_score;
	subs5 = find_spanning_multimiss_matches(&ignore_found_score,&nhits5,subs5,this5,genestrand,
						NREQUIRED_FAST,querylength5,query5_compress_fwd,query5_compress_rev,
						/*nmisses_allowed*/mismatch_level_5);
	mismatch_level_5++;
      }
    }

    /* 3B. Mismatches via spanning set.  Requires spanning set and compress. */
    if (spanningset3p /*&& allvalidp3*/ && done_level_3 >= 2) {
      /* NOTE: Since done_level isn't updated, can do in one batch instead of iteratively */
      while (mismatch_level_3 <= fast_level_3 && mismatch_level_3 <= done_level_3) {
	debug(printf("*** Stage 3B (level %d).  Spanning set mismatches ***\n",mismatch_level_3));
	did_alignment_p = true;
	ignore_found_score = *found_score;
	subs3 = find_spanning_multimiss_matches(&ignore_found_score,&nhits3,subs3,this3,genestrand,
						NREQUIRED_FAST,querylength3,query3_compress_fwd,query3_compress_rev,
						/*nmisses_allowed*/mismatch_level_3);
	mismatch_level_3++;
      }
    }

    if (did_alignment_p == true) {
      /* 3. Pairing after spanning set subs */
      hitarray5[HITARRAY_SUBS] = subs5 /* = Stage3end_remove_duplicates(subs5,queryseq5,queryseq3) */;
      hitarray3[HITARRAY_SUBS] = subs3 /* = Stage3end_remove_duplicates(subs3,queryseq5,queryseq3) */;
      hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					   &(*samechr),&(*conc_transloc),&(*terminals),
					   hitpairs,hitarray5,/*narray5*/HITARRAY_SUBS+1,
					   hitarray3,/*narray3*/HITARRAY_SUBS+1,
					   *cutoff_level_5,*cutoff_level_3,
					   querylength5,querylength3,maxpairedpaths,genestrand);
      debug(printf("After pairing spanning set, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		   nconcordant,nsamechr,List_length(*terminals),*found_score));
      if (*abort_pairing_p == true) {
	*hits5 = List_append(greedy5,subs5);
	*hits3 = List_append(greedy3,subs3);
	hitpairs = Stage3pair_remove_circular_alias(hitpairs);
#if 0
	hitpairs = Stage3pair_remove_overlaps(hitpairs,/*translocp*/false,/*finalp*/true);
#endif
	return hitpairs;
      } else {
	opt_level = (*found_score < opt_level) ? *found_score : opt_level;
	if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	  done_level_5 = user_maxlevel_5;
	}
	if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	  done_level_3 = user_maxlevel_3;
	}
	debug(printf("3> found_score = %d, opt_level %d, done_level %d,%d\n",*found_score,opt_level,done_level_5,done_level_3));
      }
    }
  }


  /* Search 3: Subs/indels via complete set algorithm */

  /* 4/5A.  Complete set mismatches and indels, omitting frequent oligos */
  completeset5p = completeset3p = false;
  if (nconcordant == 0) {
    if (*terminals == NULL) {
      debug(printf("nconcordant is 0, so we are doing completeset on both sides\n"));
      completeset5p = completeset3p = true;
    } else {
      /* Have concordant terminals, so may need to refine just one
	 end.  A single bad alignment is enough to refine that end;
	 otherwise, the algorithm misses alignments. */
      for (p = *terminals; p != NULL; p = List_next(p)) {
	newpair = (Stage3pair_T) List_head(p);
	debug(printf("Terminal has total scores of %d and %d\n",
		     Stage3end_score(Stage3pair_hit5(newpair)),
		     Stage3end_score(Stage3pair_hit3(newpair))));
	if (Stage3end_score(Stage3pair_hit5(newpair)) > done_level_5) {
	  completeset5p = true;
	}
	if (Stage3end_score(Stage3pair_hit3(newpair)) > done_level_3) {
	  completeset3p = true;
	}
	debug(printf("completeset5p %d, completeset3p %d\n",completeset5p,completeset3p));
      }
    }

#if 0
  } else if (*found_score <= done_level_5 + done_level_3) {
    debug(printf("Test for completeset: false because *found_score %d < done_level_5 %d + done_level_3 %d\n",
		 *found_score,done_level_5,done_level_3));
    completeset5p = completeset3p = false;
#endif

#if 0
  } else {
    if (better_free_end_exists_p(greedy5,subs5,terminals5,indels5,singlesplicing5,doublesplicing5,querylength5) == true) {
      completeset3p = true;	/* Do search on other end using complete set algorithm */
    }
    if (better_free_end_exists_p(greedy3,subs3,terminals3,indels3,singlesplicing3,doublesplicing3,querylength3) == true) {
      completeset5p = true;	/* Do search on other end using complete set algorithm */
    }
    debug(printf("Test for completeset using better_free_end_exists_p: completeset5p %d, completeset3p %d\n",completeset5p,completeset3p));
#endif

#if 0
  } else {
    /* This causes very slow running time */
    if (subs5 == NULL) {
      completeset5p = true;
    }
    if (subs3 == NULL) {
      completeset3p = true;
    }
#endif
  }

  if (querylength5 < min_kmer_readlength) {
    completeset5p = false;
  }
  if (querylength3 < min_kmer_readlength) {
    completeset3p = false;
  }

  if (completeset5p == true) {
    debug(printf("Performing complete set analysis on 5' end\n"));
    if (this5->read_oligos_p == false) {
      read_oligos(&allvalidp5,this5,queryuc_ptr_5,querylength5,query5_lastpos,genestrand);
    }

    floors5 = compute_floors(&any_omitted_p_5,&alloc_floors_p_5,floors_array,this5,queryuc_ptr_5,
			     querylength5,query5_lastpos,plus_indexdb_5,minus_indexdb_5,
			     indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
			     keep_floors_p);
    floors5_computed_p = true;
    ignore_found_score = *found_score;
    complete_set_mm_indels(&ignore_found_score,&segments5_computed_p,
			   &plus_anchor_segments_5,&minus_anchor_segments_5,
			   &n_plus_anchors_5,&n_minus_anchors_5,
			   &opt_level,&done_level_5,user_maxlevel_5,/*revise_levels_p*/false,
			   &nhits5,&subs5,&indels5,this5,query5_compress_fwd,query5_compress_rev,
#if defined(DEBUG2) || defined(DEBUG2E)
			   queryuc_ptr_5,queryrc5,
#endif
			   querylength5,query5_lastpos,floors5,indel_penalty_middle,indel_penalty_end,
			   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			   fast_level_5,genestrand);
  }

  /* 4/5B.  Complete set mismatches and indels, omitting frequent oligos */
  if (completeset3p == true) {
    debug(printf("Performing complete set analysis on 3' end\n"));

    if (this3->read_oligos_p == false) {
      read_oligos(&allvalidp3,this3,queryuc_ptr_3,querylength3,query3_lastpos,genestrand);
    }

    floors3 = compute_floors(&any_omitted_p_3,&alloc_floors_p_3,floors_array,this3,queryuc_ptr_3,
			     querylength3,query3_lastpos,plus_indexdb_3,minus_indexdb_3,
			     indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
			     keep_floors_p);
    floors3_computed_p = true;
    ignore_found_score = *found_score;
    complete_set_mm_indels(&ignore_found_score,&segments3_computed_p,
			   &plus_anchor_segments_3,&minus_anchor_segments_3,
			   &n_plus_anchors_3,&n_minus_anchors_3,
			   &opt_level,&done_level_3,user_maxlevel_3,/*revise_levels_p*/false,
			   &nhits3,&subs3,&indels3,this3,query3_compress_fwd,query3_compress_rev,
#if defined(DEBUG2) || defined(DEBUG2E)
			   queryuc_ptr_3,queryrc3,
#endif
			   querylength3,query3_lastpos,floors3,indel_penalty_middle,indel_penalty_end,
			   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			   fast_level_3,genestrand);
  }

  debug(printf("complete: 5' end has %d subs, %d indels, %d single splices\n",
	       List_length(subs5),List_length(indels5),List_length(singlesplicing5)));
  debug(printf("complete: 3' end has %d subs, %d indels, %d single splices\n",
	       List_length(subs3),List_length(indels3),List_length(singlesplicing3)));

  if (completeset5p == true || completeset3p == true) {
    /* 4/5. Pairing after complete set subs and indels */
    debug(printf("Starting pairing of 4 and 5\n"));
    hitarray5[HITARRAY_SUBS] = subs5;
    hitarray5[HITARRAY_INDELS] = indels5;
    hitarray3[HITARRAY_SUBS] = subs3;
    hitarray3[HITARRAY_INDELS] = indels3;
    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					 &(*samechr),&(*conc_transloc),&(*terminals),
					 hitpairs,hitarray5,/*narray5*/HITARRAY_INDELS+1,
					 hitarray3,/*narray3*/HITARRAY_INDELS+1,
					 *cutoff_level_5,*cutoff_level_3,
					 querylength5,querylength3,maxpairedpaths,genestrand);
    debug(printf("After pairing complete set mismatches and indels, found %d concordant, %d nsamechr, %d terminals, found_score %d\n",
		 nconcordant,nsamechr,List_length(*terminals),*found_score));
    if (*abort_pairing_p == true) {
      *hits5 = List_append(greedy5,List_append(subs5,indels5));
      *hits3 = List_append(greedy3,List_append(subs3,indels3));
#if 0
      hitpairs = Stage3pair_remove_circular_alias(hitpairs);
#endif
      hitpairs = Stage3pair_remove_overlaps(hitpairs,/*translocp*/false,/*finalp*/true);
      return hitpairs;
    } else {
      opt_level = (*found_score < opt_level) ? *found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("4/5> found_score = %d, opt_level %d, done_level %d,%d\n",*found_score,opt_level,done_level_5,done_level_3));
    }

    if (*found_score <= done_level_5 + done_level_3) {
      debug(printf("Test for completeset: false because *found_score %d <done_level_5 %d + done_level_3 %d\n",
		   *found_score,done_level_5,done_level_3));
      completeset5p = completeset3p = false;
    }
  }


  /* Search 4: Segments -> single splicing */

  /* 6/7/8. Local splicing.  Requires compress and all positions fetched. */
  /* Subtract 1 from done_level for previous hits */
  did_singlesplicing5_p = false;
  /* SPEED: For more hits, turn off first branch */
  if (use_sarray_p == true && completeset5p == false) {
    /* Skip.  Suffix array already found something */
    debug(printf("Skipping complete set on 5', because sarray found a hitpair\n"));

  } else if (knownsplicingp || novelsplicingp) {
    debug(printf("Deciding whether to do singlesplicing: done_level_5 %d >=? localsplicing_penalty %d\n",
		 done_level_5,localsplicing_penalty));

    if (done_level_5 >= localsplicing_penalty) {
      debug(printf("*** Stage 6A.  Single splicing masking frequent oligos with done_level %d ***\n",done_level_5));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      if (floors5_computed_p == false) {
	floors5 = compute_floors(&any_omitted_p_5,&alloc_floors_p_5,floors_array,this5,queryuc_ptr_5,
				 querylength5,query5_lastpos,plus_indexdb_5,minus_indexdb_5,
				 indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
				 keep_floors_p);
	floors5_computed_p = true;
      }

      if (segments5_computed_p == false) {
	this5->plus_segments = identify_all_segments(&this5->plus_nsegments,&plus_anchor_segments_5,&n_plus_anchors_5,
						     &this5->plus_spliceable,&this5->plus_nspliceable,
#ifdef LARGE_GENOMES
						     this5->plus_positions_high,this5->plus_positions_low,
#else
						     this5->plus_positions,
#endif
						     this5->plus_npositions,this5->omitted,querylength5,query5_lastpos,
						     floors5,/*plusp*/true);
	this5->minus_segments = identify_all_segments(&this5->minus_nsegments,&minus_anchor_segments_5,&n_minus_anchors_5,
						      &this5->minus_spliceable,&this5->minus_nspliceable,
#ifdef LARGE_GENOMES
						      this5->minus_positions_high,this5->minus_positions_low,
#else
						      this5->minus_positions,
#endif
						      this5->minus_npositions,this5->omitted,querylength5,query5_lastpos,
						      floors5,/*plusp*/false);
	segments5_computed_p = true;
      }

      did_singlesplicing5_p = true;
      ignore_found_score = *found_score;
      singlesplicing5 = complete_set_singlesplicing(&ignore_found_score,singlesplicing5,floors5,this5,
						    query5_compress_fwd,query5_compress_rev,
						    querylength5,query5_lastpos,
						    localsplicing_penalty,
						    /*max_mismatches_allowed*/done_level_5 - localsplicing_penalty,
						    genestrand,/*first_read_p*/true,
						    /*subs_or_indels_p*/(subs5 != NULL || indels5 != NULL) ? true : false);
    }
  }

  did_singlesplicing3_p = false;
  /* SPEED: For more hits, turn off first branch */
  if (use_sarray_p == true && completeset3p == false) {
    /* Skip.  Suffix array already found something */
    debug(printf("Skipping complete set on 3', because sarray found a hitpair\n"));

  } else if (knownsplicingp || novelsplicingp) {
    debug(printf("Deciding whether to do singlesplicing: done_level_3 %d >=? localsplicing_penalty %d\n",
		 done_level_3,localsplicing_penalty));
    if (done_level_3 >= localsplicing_penalty) {
      debug(printf("*** Stage 6B.  Single splicing masking frequent oligos with done_level %d ***\n",done_level_3));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      if (floors3_computed_p == false) {
	floors3 = compute_floors(&any_omitted_p_3,&alloc_floors_p_3,floors_array,this3,queryuc_ptr_3,
				 querylength3,query3_lastpos,plus_indexdb_3,minus_indexdb_3,
				 indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
				 keep_floors_p);
	floors3_computed_p = true;
      }

      if (segments3_computed_p == false) {
	this3->plus_segments = identify_all_segments(&this3->plus_nsegments,&plus_anchor_segments_3,&n_plus_anchors_3,
						     &this3->plus_spliceable,&this3->plus_nspliceable,
#ifdef LARGE_GENOMES
						     this3->plus_positions_high,this3->plus_positions_low,
#else
						     this3->plus_positions,
#endif
						     this3->plus_npositions,this3->omitted,querylength3,query3_lastpos,
						     floors3,/*plusp*/true);
	this3->minus_segments = identify_all_segments(&this3->minus_nsegments,&minus_anchor_segments_3,&n_minus_anchors_3,
						      &this3->minus_spliceable,&this3->minus_nspliceable,
#ifdef LARGE_GENOMES
						      this3->minus_positions_high,this3->minus_positions_low,
#else
						      this3->minus_positions,
#endif
						      this3->minus_npositions,this3->omitted,querylength3,query3_lastpos,
						      floors3,/*plusp*/false);
	segments3_computed_p = true;
      }

      did_singlesplicing3_p = true;
      ignore_found_score = *found_score;
      singlesplicing3 = complete_set_singlesplicing(&ignore_found_score,singlesplicing3,floors3,this3,
						    query3_compress_fwd,query3_compress_rev,
						    querylength3,query3_lastpos,
						    localsplicing_penalty,
						    /*max_mismatches_allowed*/done_level_3 - localsplicing_penalty,
						    genestrand,/*first_read_p*/false,
						    /*subs_or_indels_p*/(subs3 != NULL || indels3 != NULL) ? true : false);
    }
  }

  if (did_singlesplicing5_p == true || did_singlesplicing3_p == true) {
    /* 6.  Pairing after single splicing */
    /* Mark ambiguous splices only for single-end reads */
    hitarray5[HITARRAY_SINGLESPLICING] = singlesplicing5;
    hitarray3[HITARRAY_SINGLESPLICING] = singlesplicing3;

    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					 &(*samechr),&(*conc_transloc),&(*terminals),
					 hitpairs,hitarray5,/*narray5*/HITARRAY_SINGLESPLICING+1,
					 hitarray3,/*narray3*/HITARRAY_SINGLESPLICING+1,
					 *cutoff_level_5,*cutoff_level_3,
					 querylength5,querylength3,maxpairedpaths,genestrand);
    debug(printf("After pairing single splicing, found %d concordant, %d nsamechr, %d terminals, found_score %d\n",
		 nconcordant,nsamechr,List_length(*terminals),*found_score));
    if (*abort_pairing_p == true) {
      if (alloc_floors_p_5 == true) {
	Floors_free(&floors5);
      }
      if (alloc_floors_p_3 == true) {
	Floors_free(&floors3);
      }
      *hits5 = List_append(greedy5,List_append(subs5,List_append(indels5,singlesplicing5)));
      *hits3 = List_append(greedy3,List_append(subs3,List_append(indels3,singlesplicing3)));
#if 0
      hitpairs = Stage3pair_remove_circular_alias(hitpairs);
#endif
      hitpairs = Stage3pair_remove_overlaps(hitpairs,/*translocp*/false,/*finalp*/true);
      return hitpairs;

    } else {
      opt_level = (*found_score < opt_level) ? *found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("Pairing after 6A and 6B> found_score = %d, opt_level %d, done_level %d,%d\n",
		   *found_score,opt_level,done_level_5,done_level_3));
    }


    /* Search 5: Segments -> single splicing */
#ifdef PERFORM_DOUBLESPLICING

    /* 7.  Double splicing.  Probably found instead by segment-to-GMAP algorithm */
    if (done_level_5 >= localsplicing_penalty) {
      debug(printf("*** Stage 7A.  Double splicing masking frequent oligos with done_level %d ***\n",done_level_5));
      if (floors5_computed_p == false) {
	floors5 = compute_floors(&any_omitted_p_5,&alloc_floors_p_5,floors_array,this5,queryuc_ptr_5,
				 querylength5,query5_lastpos,plus_indexdb_5,minus_indexdb_5,
				 indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
				 keep_floors_p);
	floors5_computed_p = true;
      }
      ignore_found_score = *found_score;
      doublesplicing5 = complete_set_doublesplicing(&ignore_found_score,doublesplicing5,floors5,this5,
						    query5_compress_fwd,query5_compress_rev,
						    queryuc_ptr_5,queryrc5,querylength5,query5_lastpos,
						    localsplicing_penalty,min_shortend,
						    /*max_mismatches_allowed*/done_level_5 - localsplicing_penalty,
						    /*pairedp*/true,genestrand,/*first_read_p*/true,
						    /*subs_or_indels_p*/(subs5 != NULL || indels5 != NULL) ? true : false);
    }

    if (done_level_3 >= localsplicing_penalty) {
      debug(printf("*** Stage 7B.  Double splicing masking frequent oligos with done_level %d ***\n",done_level_3));
      if (floors3_computed_p == false) {
	floors3 = compute_floors(&any_omitted_p_3,&alloc_floors_p_3,floors_array,this3,queryuc_ptr_3,
				 querylength3,query3_lastpos,plus_indexdb_3,minus_indexdb_3,
				 indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
				 keep_floors_p);
	floors3_computed_p = true;
      }
      ignore_found_score = *found_score;
      doublesplicing3 = complete_set_doublesplicing(&ignore_found_score,doublesplicing3,floors3,this3,
						    query3_compress_fwd,query3_compress_rev,
						    queryuc_ptr_3,queryrc3,querylength3,query3_lastpos,
						    localsplicing_penalty,min_shortend,
						    /*max_mismatches_allowed*/done_level_3 - localsplicing_penalty,
						    /*pairedp*/true,genestrand,/*first_read_p*/false,
						    /*subs_or_indels_p*/(subs3 != NULL || indels3 != NULL) ? true : false);
    }
    
    /* 7.  Pairing after double splicing */
    /* Mark ambiguous splices only for single-end reads */
    hitarray5[HITARRAY_DOUBLESPLICING] = doublesplicing5;
    hitarray3[HITARRAY_DOUBLESPLICING] = doublesplicing3;

    debug(printf("Starting Stage3_pair_up_concordant\n"));
    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					 &(*samechr),&(*conc_transloc),&(*terminals),
					 hitpairs,hitarray5,/*narray5*/HITARRAY_DOUBLESPLICING+1,
					 hitarray3,/*narray3*/HITARRAY_DOUBLESPLICING+1,
					 *cutoff_level_5,*cutoff_level_3,
					 querylength5,querylength3,maxpairedpaths,genestrand);
    debug(printf("After pairing double splicing, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		 nconcordant,nsamechr,List_length(*terminals),*found_score));
    if (*abort_pairing_p == true) {
      if (alloc_floors_p_5 == true) {
	Floors_free(&floors5);
      }
      if (alloc_floors_p_3 == true) {
	Floors_free(&floors3);
      }
      *hits5 = List_append(greedy5,List_append(subs5,List_append(indels5,List_append(singlesplicing5,doublesplicing5))));
      *hits3 = List_append(greedy3,List_append(subs3,List_append(indels3,List_append(singlesplicing3,doublesplicing3))));
      hitpairs = Stage3pair_remove_circular_alias(hitpairs);
#if 0
      hitpairs = Stage3pair_remove_overlaps(hitpairs,/*translocp*/false,/*finalp*/true);
#endif
      return hitpairs;

    } else {
      opt_level = (*found_score < opt_level) ? *found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("Pairing after 7A and 7B> found_score = %d, opt_level %d, done_level %d,%d\n",
		   *found_score,opt_level,done_level_5,done_level_3));
    }
#endif
  }


  debug(printf("nconcordant = %d.  found_score %d, trigger_score %d, done level %d + %d\n",
	       nconcordant,*found_score,trigger_score_for_gmap,done_level_5,done_level_3));
  
  *hits5 = List_append(greedy5,List_append(subs5,List_append(indels5,List_append(singlesplicing5,doublesplicing5))));
  *hits3 = List_append(greedy3,List_append(subs3,List_append(indels3,List_append(singlesplicing3,doublesplicing3))));


  /* Search 6: Paired segments -> GMAP via segments */

  gmap5p = gmap3p = true;
  if (gmap_segments_p == false) {
    debug(printf("gmap_segments_p is false, so setting gmap5p and gmap3p false\n"));
    gmap5p = gmap3p = false;
  } else if (*abort_pairing_p == true) {
    debug(printf("abort_pairing_p is true, so setting gmap5p and gmap3p false\n"));
    gmap5p = gmap3p = false;
  } else if (nconcordant > 0) {
    /* Rely upon GMAP improvement instead */
    debug(printf("nconcordant > 0, so setting gmap5p and gmap3p false\n"));
    gmap5p = gmap3p = false;
  } else if (*found_score < trigger_score_for_gmap) {
    debug(printf("found_score %d < trigger_score_for_gmap %d, so setting gmap5p and gmap3p false\n",
		 *found_score,trigger_score_for_gmap));
    gmap5p = gmap3p = false;
  } else if (*found_score < done_level_5 + done_level_3) {
    debug(printf("found_score %d < done_level_5 %d + done_level_3 %d, so setting gmap5p and gmap3p false\n",
		 *found_score,done_level_5,done_level_3));
    gmap5p = gmap3p = false;
  }

  if (gmap5p == true || gmap3p == true) {
    debug(printf("***Trying to pair up segments***\n"));
    pair_up_anchor_segments(plus_anchor_segments_5,minus_anchor_segments_5,
			    plus_anchor_segments_3,minus_anchor_segments_3,
			    n_plus_anchors_5,n_minus_anchors_5,
			    n_plus_anchors_3,n_minus_anchors_3,pairmax);

    if (gmap5p == true) {
      gmap5_hits = convert_plus_segments_to_gmap(/*hits*/NULL,queryuc_ptr_5,querylength5,query5_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
						 queryrc5,Shortread_invertedp(queryseq5),
#endif
						 plus_anchor_segments_5,n_plus_anchors_5,this5->plus_segments,this5->plus_nsegments,
						 oligoindices_minor,
						 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						 genestrand,/*require_pairing_p*/true);
      gmap5_hits = convert_minus_segments_to_gmap(/*hits*/gmap5_hits,queryuc_ptr_5,querylength5,query5_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
						  queryrc5,Shortread_invertedp(queryseq5),
#endif
						  minus_anchor_segments_5,n_minus_anchors_5,this5->minus_segments,this5->minus_nsegments,
						  oligoindices_minor,
						  pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						  genestrand,/*require_pairing_p*/true);
#if 0
      /* Note: cannot use hitarray after we have removed overlapping alignments.  Have to point to hits5 and hits3 and set narray5 = narray3 = 1 */
      hitarray5[HITARRAY_SEGMENTS_GMAP] = gmap5_hits;
#else
      *hits5 = List_append(*hits5,gmap5_hits);
#endif
    }

    if (gmap3p == true) {
      gmap3_hits = convert_plus_segments_to_gmap(/*hits*/NULL,queryuc_ptr_3,querylength3,query3_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
						 queryrc3,Shortread_invertedp(queryseq3),
#endif
						 plus_anchor_segments_3,n_plus_anchors_3,this3->plus_segments,this3->plus_nsegments,
						 oligoindices_minor,
						 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						 genestrand,/*require_pairing_p*/true);
      gmap3_hits = convert_minus_segments_to_gmap(/*hits*/gmap3_hits,queryuc_ptr_3,querylength3,query3_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
						  queryrc3,Shortread_invertedp(queryseq3),
#endif
						  minus_anchor_segments_3,n_minus_anchors_3,this3->minus_segments,this3->minus_nsegments,
						  oligoindices_minor,
						  pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						  genestrand,/*require_pairing_p*/true);
#if 0
      /* Note: cannot use hitarray after we have removed overlapping alignments.  Have to point to hits5 and hits3 and set narray5 = narray3 = 1 */
      hitarray3[HITARRAY_SEGMENTS_GMAP] = gmap3_hits;
#else
      *hits3 = List_append(*hits3,gmap3_hits);
#endif
    }
  }

  if (gmap5_hits != NULL || gmap3_hits != NULL) {
    /* found_terminals_p = true; */
    debug4t(printf("Running Stage3_pair_up_concordant\n"));
    /* Note: cannot use hitarray after we have removed overlapping alignments */
    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					 &(*samechr),&(*conc_transloc),&(*terminals),
                                         hitpairs,/*hitarray5*/&(*hits5),/*narray5*/1,
					 /*hitarray3*/&(*hits3),/*narray3*/1,
					 *cutoff_level_5,*cutoff_level_3,
					 querylength5,querylength3,maxpairedpaths,genestrand);
    debug(printf("11> After pairing GMAP, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		 nconcordant,nsamechr,List_length(*terminals),*found_score));
    if (*abort_pairing_p == false) {
      opt_level = (*found_score < opt_level) ? *found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("Pairing after 11A and 11B> found_score = %d, opt_level %d, done_level %d,%d\n",
		   *found_score,opt_level,done_level_5,done_level_3));
    }
  }


  /* Search 7: Distant splicing */
  min_trim = querylength5 + querylength3;
  for (p = hitpairs; p != NULL; p = List_next(p)) {
    newpair = (Stage3pair_T) List_head(p);
    if ((trim = Stage3pair_total_trim(newpair)) < min_trim) {
      min_trim = trim;
    }
  }

  debug(printf("? distant splicing.  nconcordant is %d.  min_trim is %d.  done_level is %d + %d\n",
	       nconcordant,min_trim,done_level_5,done_level_3));
  if (nconcordant > 0 && min_trim < min_distantsplicing_end_matches) {
    /* Skip search for distant splicing */

  } else if (*abort_pairing_p == true) {
    /* Skip further searching */

  } else if (knownsplicingp == false && novelsplicingp == false) {
    /* TODO: Find distant splicing for DNA */

  } else {
    /* Find distant splicing for RNA */
    if (done_level_5 >= distantsplicing_penalty) {
      /* Want >= and not >, because otherwise distant splicing does not work on 50-bp reads */
      /* Want > and not >=, because distant splicing needs to be better than other alternatives */
      max_splice_mismatches_5 = done_level_5 - distantsplicing_penalty;

      donors_plus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
      antidonors_plus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
      acceptors_plus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
      antiacceptors_plus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
      donors_minus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
      antidonors_minus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
      acceptors_minus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
      antiacceptors_minus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));

      if (floors5_computed_p == false) {
	floors5 = compute_floors(&any_omitted_p_5,&alloc_floors_p_5,floors_array,this5,queryuc_ptr_5,
				 querylength5,query5_lastpos,plus_indexdb_5,minus_indexdb_5,
				 indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
				 keep_floors_p);
	floors5_computed_p = true;
      }

      /* 11A.  Distant splicing */
      debug(printf("Starting find_spliceends (plus) on 5' end with %d anchor segments\n",n_plus_anchors_5));
      find_spliceends_distant_rna(&donors_plus_5,&antidonors_plus_5,&acceptors_plus_5,&antiacceptors_plus_5,
				  plus_anchor_segments_5,n_plus_anchors_5,
#ifdef DEBUG4E
				  /*queryptr*/queryuc_ptr_5,
#endif
				  floors5,querylength5,query5_lastpos,/*query_compress*/query5_compress_fwd,
				  max_splice_mismatches_5,/*plusp*/true,genestrand);
      debug(printf("Finished find_spliceends (plus)\n"));

      debug(printf("Starting find_spliceends (minus) on 5' end with %d anchor segments\n",n_minus_anchors_5));
      find_spliceends_distant_rna(&antidonors_minus_5,&donors_minus_5,&antiacceptors_minus_5,&acceptors_minus_5,
				  minus_anchor_segments_5,n_minus_anchors_5,
#ifdef DEBUG4E
				  /*queryptr*/queryrc5,
#endif
				  floors5,querylength5,query5_lastpos,/*query_compress*/query5_compress_rev,
				  max_splice_mismatches_5,/*plusp*/false,genestrand);
      debug(printf("Finished find_spliceends (minus)\n"));


      /* 11A.  Distant splicing */
      nmismatches = 0;
      while (longsinglesplicing5 == NULL &&
	     nmismatches <= max_splice_mismatches_5 /* && nsplicepairs5 < MAXCHIMERAPATHS */) {
	debug(printf("*** Stage 10A.  Distant splicing, allowing %d mismatches out of %d***\n",
		     nmismatches,max_splice_mismatches_5));
	
	debug4e(printf("Sorting splice ends\n"));
	donors_plus_5[nmismatches] = Substring_sort_siteD_halves(donors_plus_5[nmismatches],/*ascendingp*/true);
	acceptors_plus_5[nmismatches] = Substring_sort_siteA_halves(acceptors_plus_5[nmismatches],/*ascendingp*/true);

	antidonors_plus_5[nmismatches] = Substring_sort_siteD_halves(antidonors_plus_5[nmismatches],/*ascendingp*/false);
	antiacceptors_plus_5[nmismatches] = Substring_sort_siteA_halves(antiacceptors_plus_5[nmismatches],/*ascendingp*/false);

	donors_minus_5[nmismatches] = Substring_sort_siteD_halves(donors_minus_5[nmismatches],/*ascendingp*/false);
	acceptors_minus_5[nmismatches] = Substring_sort_siteA_halves(acceptors_minus_5[nmismatches],/*ascendingp*/false);

	antidonors_minus_5[nmismatches] = Substring_sort_siteD_halves(antidonors_minus_5[nmismatches],/*ascendingp*/true);
	antiacceptors_minus_5[nmismatches] = Substring_sort_siteA_halves(antiacceptors_minus_5[nmismatches],/*ascendingp*/true);

	debug4e(printf("Splice ends at %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		       nmismatches,
		       List_length(donors_plus_5[nmismatches]),List_length(acceptors_plus_5[nmismatches]),
		       List_length(antidonors_plus_5[nmismatches]),List_length(antiacceptors_plus_5[nmismatches]),
		       List_length(donors_minus_5[nmismatches]),List_length(acceptors_minus_5[nmismatches]),
		       List_length(antidonors_minus_5[nmismatches]),List_length(antiacceptors_minus_5[nmismatches])));

	ignore_found_score = *found_score;
	distantsplicing5 = find_splicepairs_distant_rna(&ignore_found_score,&nsplicepairs5,&longsinglesplicing5,/*hits*/distantsplicing5,
							donors_plus_5,antidonors_plus_5,acceptors_plus_5,antiacceptors_plus_5,
							donors_minus_5,antidonors_minus_5,acceptors_minus_5,antiacceptors_minus_5,
							localsplicing_penalty,distantsplicing_penalty,
							querylength5,nmismatches,/*first_read_p*/true);
	debug(printf("Found %d distant splices on 5' end\n",List_length(distantsplicing5)));
	nmismatches++;
      }

      /* Clean up 5 */
      for (i = 0; i <= max_splice_mismatches_5; i++) {
	substringlist_gc(&(donors_plus_5[i]));
	substringlist_gc(&(antidonors_plus_5[i]));
	substringlist_gc(&(acceptors_plus_5[i]));
	substringlist_gc(&(antiacceptors_plus_5[i]));
	substringlist_gc(&(donors_minus_5[i]));
	substringlist_gc(&(antidonors_minus_5[i]));
	substringlist_gc(&(acceptors_minus_5[i]));
	substringlist_gc(&(antiacceptors_minus_5[i]));
      }
      FREEA(donors_plus_5);
      FREEA(antidonors_plus_5);
      FREEA(acceptors_plus_5);
      FREEA(antiacceptors_plus_5);
      FREEA(donors_minus_5);
      FREEA(antidonors_minus_5);
      FREEA(acceptors_minus_5);
      FREEA(antiacceptors_minus_5);
    }

    if (done_level_3 >= distantsplicing_penalty) {
      /* Want >= and not >, because otherwise distant splicing does not work on 50-bp reads */
      /* Want > and not >=, because distant splicing needs to be better than other alternatives */
      max_splice_mismatches_3 = done_level_3 - distantsplicing_penalty;

      donors_plus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
      antidonors_plus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
      acceptors_plus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
      antiacceptors_plus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
      donors_minus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
      antidonors_minus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
      acceptors_minus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
      antiacceptors_minus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));

      if (floors3_computed_p == false) {
	floors3 = compute_floors(&any_omitted_p_3,&alloc_floors_p_3,floors_array,this3,queryuc_ptr_3,
				 querylength3,query3_lastpos,plus_indexdb_3,minus_indexdb_3,
				 indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
				 keep_floors_p);
	floors3_computed_p = true;
      }

      /* 11B.  Distant splicing */
      debug(printf("Starting find_spliceends (plus) on 3' end with %d anchor segments\n",n_plus_anchors_3));
      find_spliceends_distant_rna(&donors_plus_3,&antidonors_plus_3,&acceptors_plus_3,&antiacceptors_plus_3,
				  plus_anchor_segments_3,n_plus_anchors_3,
#ifdef DEBUG4E
				  /*queryptr*/queryuc_ptr_3,
#endif
				  floors3,querylength3,query3_lastpos,/*query_compress*/query3_compress_fwd,
				  max_splice_mismatches_3,/*plusp*/true,genestrand);
      debug(printf("Finished find_spliceends (plus)\n"));

      debug(printf("Starting find_spliceends (minus) on 3' end with %d anchor segments\n",n_minus_anchors_3));
      find_spliceends_distant_rna(&antidonors_minus_3,&donors_minus_3,&antiacceptors_minus_3,&acceptors_minus_3,
				  minus_anchor_segments_3,n_minus_anchors_3,
#ifdef DEBUG4E
				  /*queryptr*/queryrc3,
#endif
				  floors3,querylength3,query3_lastpos,/*query_compress*/query3_compress_rev,
				  max_splice_mismatches_3,/*plusp*/false,genestrand);
      debug(printf("Finished find_spliceends (minus)\n"));

      /* 11B.  Distant splicing */
      nmismatches = 0;
      while (longsinglesplicing3 == NULL &&
	     nmismatches <= max_splice_mismatches_3 /* && nsplicepairs3 < MAXCHIMERAPATHS */) {
	debug(printf("*** Stage 10B.  Distant splicing, allowing %d mismatches out of %d***\n",
		     nmismatches,max_splice_mismatches_3));

	debug4e(printf("Sorting splice ends\n"));
	donors_plus_3[nmismatches] = Substring_sort_siteD_halves(donors_plus_3[nmismatches],/*ascendingp*/true);
	acceptors_plus_3[nmismatches] = Substring_sort_siteA_halves(acceptors_plus_3[nmismatches],/*ascendingp*/true);

	antidonors_plus_3[nmismatches] = Substring_sort_siteD_halves(antidonors_plus_3[nmismatches],/*ascendingp*/false);
	antiacceptors_plus_3[nmismatches] = Substring_sort_siteA_halves(antiacceptors_plus_3[nmismatches],/*ascendingp*/false);

	donors_minus_3[nmismatches] = Substring_sort_siteD_halves(donors_minus_3[nmismatches],/*ascendingp*/false);
	acceptors_minus_3[nmismatches] = Substring_sort_siteA_halves(acceptors_minus_3[nmismatches],/*ascendingp*/false);

	antidonors_minus_3[nmismatches] = Substring_sort_siteD_halves(antidonors_minus_3[nmismatches],/*ascendingp*/true);
	antiacceptors_minus_3[nmismatches] = Substring_sort_siteA_halves(antiacceptors_minus_3[nmismatches],/*ascendingp*/true);

	debug4e(printf("Splice ends at %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		       nmismatches,
		       List_length(donors_plus_3[nmismatches]),List_length(acceptors_plus_3[nmismatches]),
		       List_length(antidonors_plus_3[nmismatches]),List_length(antiacceptors_plus_3[nmismatches]),
		       List_length(donors_minus_3[nmismatches]),List_length(acceptors_minus_3[nmismatches]),
		       List_length(antidonors_minus_3[nmismatches]),List_length(antiacceptors_minus_3[nmismatches])));

	ignore_found_score = *found_score;
	distantsplicing3 = find_splicepairs_distant_rna(&ignore_found_score,&nsplicepairs3,&longsinglesplicing3,/*hits*/distantsplicing3,
							donors_plus_3,antidonors_plus_3,acceptors_plus_3,antiacceptors_plus_3,
							donors_minus_3,antidonors_minus_3,acceptors_minus_3,antiacceptors_minus_3,
							localsplicing_penalty,distantsplicing_penalty,
							querylength3,nmismatches,/*first_read_p*/false);
	debug(printf("Found %d distant splices on 5' end\n",List_length(distantsplicing3)));
	nmismatches++;
      }

      /* Clean up 3 */
      for (i = 0; i <= max_splice_mismatches_3; i++) {
	substringlist_gc(&(donors_plus_3[i]));
	substringlist_gc(&(antidonors_plus_3[i]));
	substringlist_gc(&(acceptors_plus_3[i]));
	substringlist_gc(&(antiacceptors_plus_3[i]));
	substringlist_gc(&(donors_minus_3[i]));
	substringlist_gc(&(antidonors_minus_3[i]));
	substringlist_gc(&(acceptors_minus_3[i]));
	substringlist_gc(&(antiacceptors_minus_3[i]));
      }
      FREEA(donors_plus_3);
      FREEA(antidonors_plus_3);
      FREEA(acceptors_plus_3);
      FREEA(antiacceptors_plus_3);
      FREEA(donors_minus_3);
      FREEA(antidonors_minus_3);
      FREEA(acceptors_minus_3);
      FREEA(antiacceptors_minus_3);
    }

    /* 11.  Pairing after distant splicing using longsinglesplicing */

    if (longsinglesplicing5 != NULL || longsinglesplicing3 != NULL) {
#if 0
      /* Note: cannot use hitarray after we have removed overlapping alignments.  Have to point to hits5 and hits3 and set narray5 = narray3 = 1 */
      hitarray5[HITARRAY_LONGSINGLESPLICING] = longsinglesplicing5;
      hitarray3[HITARRAY_LONGSINGLESPLICING] = longsinglesplicing3;
#else
      if (longsinglesplicing5 != NULL) {
	*hits5 = List_append(*hits5,longsinglesplicing5);
      }
      if (longsinglesplicing3 != NULL) {
	*hits3 = List_append(*hits3,longsinglesplicing3);
      }
#endif
      /* Note: cannot use hitarray after we have removed overlapping alignments.  Have to point to hits5 and hits3 and set narray5 = narray3 = 1 */
      hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					   &(*samechr),&(*conc_transloc),&(*terminals),
					   hitpairs,/*hitarray5*/&(*hits5),/*narray5*/1,
					   /*hitarray3*/&(*hits3),/*narray3*/1,
					   *cutoff_level_5,*cutoff_level_3,
					   querylength5,querylength3,maxpairedpaths,genestrand);
      debug(printf("10> After pairing long single splicing, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		   nconcordant,nsamechr,List_length(*terminals),*found_score));

      if (*abort_pairing_p == false) {
	opt_level = (*found_score < opt_level) ? *found_score : opt_level;
	if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	  done_level_5 = user_maxlevel_5;
	}
	if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	  done_level_3 = user_maxlevel_3;
	}
	debug(printf("10> found_score = %d, opt_level %d, done_level %d,%d\n",*found_score,opt_level,done_level_5,done_level_3));
      }

    }

    /* 11.  Pairing after distant splicing using distantsplicing */
#if 0
    /* Note: cannot use hitarray after we have removed overlapping alignments.  Have to point to hits5 and hits3 and set narray5 = narray3 = 1 */
    hitarray5[HITARRAY_DISTANTSPLICING] = distantsplicing5;
    hitarray3[HITARRAY_DISTANTSPLICING] = distantsplicing3;
#else
    if (distantsplicing5 != NULL) {
      *hits5 = List_append(*hits5,distantsplicing5);
    }
    if (distantsplicing3 != NULL) {
      *hits3 = List_append(*hits3,distantsplicing3);
    }
#endif

    if (nconcordant == 0 && (distantsplicing5 != NULL || distantsplicing3 != NULL)) {
      /* Note: cannot use hitarray after we have removed overlapping alignments.  Have to point to hits5 and hits3 and set narray5 = narray3 = 1 */
      hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					   &(*samechr),&(*conc_transloc),&(*terminals),
					   hitpairs,/*hitarray5*/&(*hits5),/*narray5*/1,
					   /*hitarray3*/&(*hits3),/*narray3*/1,
					   *cutoff_level_5,*cutoff_level_3,
					   querylength5,querylength3,maxpairedpaths,genestrand);
      debug(printf("11> After pairing distant splicing, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		   nconcordant,nsamechr,List_length(*terminals),*found_score));

      if (*abort_pairing_p == false) {
	opt_level = (*found_score < opt_level) ? *found_score : opt_level;
	if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	  done_level_5 = user_maxlevel_5;
	}
	if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	  done_level_3 = user_maxlevel_3;
	}
	debug(printf("10> found_score = %d, opt_level %d, done_level %d,%d\n",*found_score,opt_level,done_level_5,done_level_3));
      }
    }
  }

  /* Search 8: Terminals */
  if (0) {
  if (nconcordant == 0 || *found_score > opt_level) {
    terminals5 = find_terminals(plus_anchor_segments_5,minus_anchor_segments_5,
				n_plus_anchors_5,n_minus_anchors_5,
				querylength5,query5_lastpos,
				query5_compress_fwd,query5_compress_rev,
				/*max_mismatches_allowed*/done_level_5,genestrand);
    *hits5 = List_append(*hits5,terminals5);

    terminals3 = find_terminals(plus_anchor_segments_3,minus_anchor_segments_3,
				n_plus_anchors_3,n_minus_anchors_3,
				querylength3,query3_lastpos,
				query3_compress_fwd,query3_compress_rev,
				/*max_mismatches_allowed*/done_level_3,genestrand);
    *hits3 = List_append(*hits3,terminals3);

    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					 &(*samechr),&(*conc_transloc),&(*terminals),
					 hitpairs,/*hitarray5*/&(*hits5),/*narray5*/1,
					 /*hitarray3*/&(*hits3),/*narray3*/1,
					 *cutoff_level_5,*cutoff_level_3,
					 querylength5,querylength3,maxpairedpaths,genestrand);
    debug(printf("After pairing terminals, found %d concordant, %d nsamechr, %d terminals, found_score %d\n",
		 nconcordant,nsamechr,List_length(*terminals),*found_score));
  }
  }


  /* Search 9: GMAP pairsearch/halfmapping */
  if (gmap_pairsearch_p == true) {
    if (*abort_pairing_p == true) {
      /* Don't do GMAP */
      gmap5p = gmap3p = false;
    } else if (nconcordant > 0) {
      gmap5p = gmap3p = false;
    } else if (*found_score >= trigger_score_for_gmap) {
      debug(printf("Test for stage 9: true because found_score %d >= trigger_score_for_gmap %d\n",*found_score,trigger_score_for_gmap));
      gmap5p = gmap3p = true;
    } else {
      gmap5p = gmap3p = false;
      if (better_free_end_exists_p(greedy5,subs5,terminals5,indels5,singlesplicing5,doublesplicing5,querylength5) == true) {
	gmap3p = true;	/* Do GMAP on other end */
      }
      if (better_free_end_exists_p(greedy3,subs3,terminals3,indels3,singlesplicing3,doublesplicing3,querylength3) == true) {
	gmap5p = true;	/* Do GMAP on other end */
      }
      debug(printf("Test for stage 9 using better_free_end_exists_p: gmap5p %d, gmap3p %d\n",gmap5p,gmap3p));
    }


    /* 9A,B.  GMAP pairsearch/halfmapping/unpaired */
    /* Our previous test for doing GMAP was if nconcordant == 0, but
       could lead to a false positive concordant match. */
    /* Do not update nconcordant, because poor GMAP alignments can stop search for a distant splice */

    /* Relying upon trigger_score_for_gmap can occasionally lead to faulty concordant alignments.  However, running it on everything
       speed by half */

    /* Go ahead and resolve overlaps on each end by Stage3end, since
       we cannot do it by Stage3pair, but do not apply optimal
       score */

    /* Previously did pairsearch only if hits were limited, but affected by poor terminals */
    if (gmap3p == true) {
      debug(printf("Before remove_overlaps of 5' at cutoff level %d: %d hits\n",*cutoff_level_5,List_length(*hits5)));
      *hits5 = Stage3end_sort_bymatches(Stage3end_remove_overlaps(*hits5,/*finalp*/false));
      debug(printf("After remove_overlaps (3): %d\n",List_length(*hits5)));

      i = 0;
      best_score_paired = Stage3end_best_score_paired(*hits5);
      debug13(printf("%d hits on 5' end\n",List_length(*hits5)));
      debug13(printf("For pairsearch, running GMAP on 3' end to match with 5' ends with score <= score %d\n",
		     best_score_paired));
      for (p = *hits5; p != NULL && i < max_gmap_pairsearch; p = List_next(p)) {
	hit5 = (Stage3end_T) List_head(p);
	if (Stage3end_hittype(hit5) == TRANSLOC_SPLICE) {
	  debug13(printf("No GMAP on transloc splice\n"));
	} else if (Stage3end_paired_usedp(hit5) == false && Stage3end_score(hit5) <= best_score_paired) {
	  halfmapping3 = align_halfmapping_with_gmap(gmap_history_3,hit5,/*hit3*/NULL,queryseq5,queryseq3,
						     queryuc_ptr_3,/*querylength*/querylength3,query3_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
						     queryrc3,Shortread_invertedp(queryseq3),
#endif
						     this3->plus_segments,this3->plus_nsegments,this3->minus_segments,this3->minus_nsegments,
						     oligoindices_major,oligoindices_minor,
						     pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						     pairmax,shortsplicedist,genestrand);
	  for (a = halfmapping3; a != NULL; a = List_next(a)) {
	    gmap3 = (Stage3end_T) List_head(a);
	    debug13(printf("=> Successful pairsearch GMAP on hit3 with score %d and nmatches %d\n",
			   Stage3end_score(gmap3),Stage3end_nmatches_posttrim(gmap3)));
#if 0
	    if (Stage3end_score(gmap3) > *cutoff_level_3 + gmap_allowance) {
	      /* nsalvage += 1; */
	      debug13(printf("Score is only %d vs cutoff level %d\n",Stage3end_score(gmap3),*cutoff_level_3));
	      Stage3end_free(&gmap3);

	    } else if ((newpair = Stage3pair_new(Stage3end_copy(hit5),gmap3,genestrand,/*pairtype*/CONCORDANT,
						 /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/true)) == NULL) {
	      debug13(printf(  "newpair is NULL\n"));
	      /* Stage3end_free(&gmap3); -- done by Stage3pair_new */

	    } else if (Stage3end_hittype(hit5) != TERMINAL) {
	      if (Stage3end_nmatches_posttrim(gmap3) >= querylength3 - (*cutoff_level_3) &&
		  Stage3end_gmap_max_match_length(gmap3) >= querylength3/2) {
		/* Want high standard for nconcordant, since this precludes finding terminals */
		nconcordant += 1;
		debug13(printf("High quality (nmatches %d >= querylength %d - cutoff level %d) => nconcordant %d\n",
			       Stage3end_nmatches_posttrim(gmap3),querylength3,*cutoff_level_3,nconcordant));
	      }
	      hitpairs = List_push(hitpairs,(void *) newpair);

	    } else if (Stage3end_trimlength(hit5) < reject_trimlength) {
	      if (Stage3end_nmatches_posttrim(gmap3) >= querylength3 - (*cutoff_level_3) &&
		  Stage3end_gmap_max_match_length(gmap3) >= querylength3/2) {
		/* Want high standard for nconcordant, since this precludes finding terminals */
		nconcordant += 1;
		debug13(printf("High quality (nmatches %d >= querylength %d - cutoff level %d) => nconcordant %d\n",
			       Stage3end_nmatches_posttrim(gmap3),querylength3,*cutoff_level_3,nconcordant));
	      }
	      hitpairs = List_push(hitpairs,(void *) newpair);
	    } else {
	      /* Stage3end_free(&gmap3); */
	      Stage3pair_free(&newpair);
	    }
#else
	    if ((newpair = Stage3pair_new(Stage3end_copy(hit5),gmap3,genestrand,/*pairtype*/CONCORDANT,
					  /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/true)) == NULL) {
	      debug13(printf("newpair is NULL\n"));
	      /* Stage3end_free(&gmap3); -- done by Stage3pair_new */

	    } else if (Stage3pair_determine_pairtype(newpair) != CONCORDANT) {
	      debug13(printf("  => not concordant, so eliminating\n"));

	    } else {
	      nconcordant += 1;
	      debug13(printf("New pair => nconcordant %d\n",nconcordant));
	      hitpairs = List_push(hitpairs,(void *) newpair);
	    }
#endif
	  }
	  List_free(&halfmapping3);
	  i++;
	}
      }
    }

    if (gmap5p == true) {
      debug(printf("Before remove_overlaps of 3' at cutoff level %d: %d hits\n",*cutoff_level_3,List_length(*hits3)));
      *hits3 = Stage3end_sort_bymatches(Stage3end_remove_overlaps(*hits3,/*finalp*/false));
      debug(printf("After remove_overlaps (4): %d\n",List_length(*hits3)));

      i = 0;
      best_score_paired = Stage3end_best_score_paired(*hits3);
      debug13(printf("%d hits on 3' end\n",List_length(*hits3)));
      debug13(printf("For pairsearch, running GMAP on 5' end to match with 3' ends with score <= score %d\n",
		     best_score_paired));
      for (p = *hits3; p != NULL && i < max_gmap_pairsearch; p = List_next(p)) {
	hit3 = (Stage3end_T) List_head(p);
	if (Stage3end_hittype(hit3) == TRANSLOC_SPLICE) {
	  debug13(printf("Not GMAP on transloc splice\n"));
	} else if (Stage3end_paired_usedp(hit3) == false && Stage3end_score(hit3) <= best_score_paired) {
	  halfmapping5 = align_halfmapping_with_gmap(gmap_history_5,/*hit5*/NULL,hit3,queryseq5,queryseq3,
						     queryuc_ptr_5,/*querylength*/querylength5,query5_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
						     queryrc5,Shortread_invertedp(queryseq5),
#endif
						     this5->plus_segments,this5->plus_nsegments,this5->minus_segments,this5->minus_nsegments,
						     oligoindices_major,oligoindices_minor,
						     pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						     pairmax,shortsplicedist,genestrand);
	  for (a = halfmapping5; a != NULL; a = List_next(a)) {
	    gmap5 = (Stage3end_T) List_head(a);
	    debug13(printf("=> Successful pairsearch GMAP on hit5 with score %d and nmatches %d\n",
			   Stage3end_score(gmap5),Stage3end_nmatches_posttrim(gmap5)));
#if 0
	    /* Stage3end_nmatches_posttrim(gmap5) >= querylength5 - (*cutoff_level_5); */
	    if (Stage3end_score(gmap5) > *cutoff_level_5 + gmap_allowance) {
	      /* nsalvage += 1; */
	      debug13(printf("Score is only %d vs cutoff level %d\n",Stage3end_score(gmap5),*cutoff_level_5));
	      Stage3end_free(&gmap5);

	    } else if ((newpair = Stage3pair_new(gmap5,Stage3end_copy(hit3),genestrand,/*pairtype*/CONCORDANT,
						 /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/true)) == NULL) {
	      debug13(printf(  "newpair is NULL\n"));
	      /* Stage3end_free(&gmap5); -- done by Stage3pair_new */

	    } else if (Stage3end_hittype(hit3) != TERMINAL) {
	      if (Stage3end_nmatches_posttrim(gmap5) >= querylength5 - (*cutoff_level_5) &&
		  Stage3end_gmap_max_match_length(gmap5) >= querylength5/2) {
		/* Want high standard for nconcordant, since this precludes finding terminals */
		nconcordant += 1;
		debug13(printf("High quality (nmatches %d >= querylength %d - cutoff level %d) => nconcordant %d\n",
			       Stage3end_nmatches_posttrim(gmap5),querylength5,*cutoff_level_5,nconcordant));
	      }
	      hitpairs = List_push(hitpairs,(void *) newpair);
	    } else if (Stage3end_trimlength(hit3) < reject_trimlength) {
	      if (Stage3end_nmatches_posttrim(gmap5) >= querylength5 - (*cutoff_level_5) &&
		  Stage3end_gmap_max_match_length(gmap5) >= querylength5/2) {
		/* Want high standard for nconcordant, since this precludes finding terminals */
		nconcordant += 1;
		debug13(printf("High quality (nmatches %d >= querylength %d - cutoff level %d) => nconcordant %d\n",
			       Stage3end_nmatches_posttrim(gmap5),querylength5,*cutoff_level_5,nconcordant));
	      }
	      hitpairs = List_push(hitpairs,(void *) newpair);
	    } else {
	      /* Stage3end_free(&gmap5); */
	      Stage3pair_free(&newpair);
	    }
#else
	    if ((newpair = Stage3pair_new(gmap5,Stage3end_copy(hit3),genestrand,/*pairtype*/CONCORDANT,
					  /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/true)) == NULL) {
	      debug13(printf(  "newpair is NULL\n"));
	      /* Stage3end_free(&gmap5); -- done by Stage3pair_new */

	    } else if (Stage3pair_determine_pairtype(newpair) != CONCORDANT) {
	      debug13(printf("  => not concordant, so eliminating\n"));

	    } else {
	      nconcordant += 1;
	      debug13(printf("new pair => nconcordant %d\n",nconcordant));
	      hitpairs = List_push(hitpairs,(void *) newpair);
	    }
#endif
	  }
	  List_free(&halfmapping5);
	  i++;
	}
      }
    }

    debug(printf("9> After GMAP pairsearch, found %d concordant\n",nconcordant));
  }


#if 0
  /* Unused code */
  alloc5p = false;
  if (knownsplicingp == true && done_level_5 >= localsplicing_penalty) {
    /* Want >= and not > to give better results.  Negligible effect on speed. */
    /* 8A.  Shortend splicing */
    max_splice_mismatches_5 = done_level_5 - localsplicing_penalty;
    
    alloc5p = true;
    donors_plus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
    antidonors_plus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
    acceptors_plus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
    antiacceptors_plus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
    donors_minus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
    antidonors_minus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
    acceptors_minus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));
    antiacceptors_minus_5 = (List_T *) CALLOCA(max_splice_mismatches_5+1,sizeof(List_T));

    if (floors5_computed_p == false) {
      floors5 = compute_floors(&any_omitted_p_5,&alloc_floors_p_5,floors_array,this5,queryuc_ptr_5,
			       querylength5,query5_lastpos,plus_indexdb_5,minus_indexdb_5,
			       indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
			       keep_floors_p);
      floors5_computed_p = true;
    }

    find_spliceends_shortend(&donors_plus_5,&antidonors_plus_5,&acceptors_plus_5,&antiacceptors_plus_5,
			     plus_anchor_segments_5,n_plus_anchors_5,
#ifdef DEBUG4E
			     queryuc_ptr_5,
#endif
			     floors5,querylength5,query5_lastpos,/*query_compress*/query5_compress_fwd,
			     /*max_mismatches_allowed*/max_splice_mismatches_5,/*plusp*/true,genestrand);

    find_spliceends_shortend(&antidonors_minus_5,&donors_minus_5,&antiacceptors_minus_5,&acceptors_minus_5,
			     minus_anchor_segments_5,n_minus_anchors_5,
#ifdef DEBUG4E
			     /*queryptr*/queryrc5,
#endif
			     floors5,querylength5,query5_lastpos,/*query_compress*/query5_compress_rev,
			     /*max_mismatches_allowed*/max_splice_mismatches_5,/*plusp*/false,genestrand);

    ignore_found_score = *found_score;
    singlesplicing5 = find_splicepairs_shortend(&ignore_found_score,/*hits*/singlesplicing5,
						donors_plus_5,antidonors_plus_5,acceptors_plus_5,antiacceptors_plus_5,
						donors_minus_5,antidonors_minus_5,acceptors_minus_5,antiacceptors_minus_5,
						query5_compress_fwd,query5_compress_rev,
						queryuc_ptr_5,queryrc5,min_shortend,localsplicing_penalty,
						/*max_mismatches_allowed*/max_splice_mismatches_5,querylength5,
						/*pairedp*/true,/*first_read_p*/true,genestrand);
  }


  alloc3p = false;
  if (knownsplicingp == true && done_level_3 >= localsplicing_penalty) {
    /* Want >= and not > to give better results.  Negligible effect on speed. */
    /* 8B.  Short-Overlap splicing */
    max_splice_mismatches_3 = done_level_3 - localsplicing_penalty;

    alloc3p = true;
    donors_plus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
    antidonors_plus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
    acceptors_plus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
    antiacceptors_plus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
    donors_minus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
    antidonors_minus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
    acceptors_minus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));
    antiacceptors_minus_3 = (List_T *) CALLOCA(max_splice_mismatches_3+1,sizeof(List_T));

    if (floors3_computed_p == false) {
      floors3 = compute_floors(&any_omitted_p_3,&alloc_floors_p_3,floors_array,this3,queryuc_ptr_3,
			       querylength3,query3_lastpos,plus_indexdb_3,minus_indexdb_3,
			       indexdb_size_threshold,max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true,
			       keep_floors_p);
      floors3_computed_p = true;
    }

    find_spliceends_shortend(&donors_plus_3,&antidonors_plus_3,&acceptors_plus_3,&antiacceptors_plus_3,
			     plus_anchor_segments_3,n_plus_anchors_3,
#ifdef DEBUG4E
			     queryuc_ptr_3,
#endif
			     floors3,querylength3,query3_lastpos,/*query_compress*/query3_compress_fwd,
			     /*max_mismatches_allowed*/max_splice_mismatches_3,/*plusp*/true,genestrand);

    find_spliceends_shortend(&antidonors_minus_3,&donors_minus_3,&antiacceptors_minus_3,&acceptors_minus_3,
			     minus_anchor_segments_3,n_minus_anchors_3,
#ifdef DEBUG4E
			     /*queryptr*/queryrc3,
#endif
			     floors3,querylength3,query3_lastpos,/*query_compress*/query3_compress_rev,
			     /*max_mismatches_allowed*/max_splice_mismatches_3,/*plusp*/false,genestrand);
      
    ignore_found_score = *found_score;
    singlesplicing3 = find_splicepairs_shortend(&ignore_found_score,/*hits*/singlesplicing3,
						donors_plus_3,antidonors_plus_3,acceptors_plus_3,antiacceptors_plus_3,
						donors_minus_3,antidonors_minus_3,acceptors_minus_3,antiacceptors_minus_3,
						query3_compress_fwd,query3_compress_rev,
						queryuc_ptr_3,queryrc3,min_shortend,localsplicing_penalty,
						/*max_mismatches_allowed*/max_splice_mismatches_3,querylength3,
						/*pairedp*/true,/*first_read_p*/false,genestrand);
  }

  if (singlesplicing5 != NULL || singlesplicing3 != NULL) {
    /* 8.  Pairing after short-overlaps */
    hitarray5[HITARRAY_SINGLESPLICING] = singlesplicing5 /* = Stage3end_remove_duplicates(singlesplicing5,queryseq5,queryseq3) */;
    hitarray3[HITARRAY_SINGLESPLICING] = singlesplicing3 /* = Stage3end_remove_duplicates(singlesplicing3,queryseq5,queryseq3) */;
    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&(*found_score),&nconcordant,&nsamechr,
					 &(*samechr),&(*conc_transloc),&(*terminals),
					 hitpairs,hitarray5,/*narray5*/HITARRAY_DOUBLESPLICING+1,
					 hitarray3,/*narray3*/HITARRAY_DOUBLESPLICING+1,
					 *cutoff_level_5,*cutoff_level_3,
					 querylength5,querylength3,maxpairedpaths,genestrand);
    debug(printf("After pairing short-overlap splicing, found %d concordant, %d samechr, %d terminals, found_score %d\n",
		 nconcordant,nsamechr,List_length(*terminals),*found_score));
    if (*abort_pairing_p == false) {
      opt_level = (*found_score < opt_level) ? *found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("Pairing after 8A and 8B> found_score = %d, opt_level %d, done_level %d,%d\n",
		   *found_score,opt_level,done_level_5,done_level_3));
    }
  }

  if (alloc5p == true) {
    /* Clean up 5 */
    for (i = 0; i <= max_splice_mismatches_5; i++) {
      substringlist_gc(&(donors_plus_5[i]));
      substringlist_gc(&(antidonors_plus_5[i]));
      substringlist_gc(&(acceptors_plus_5[i]));
      substringlist_gc(&(antiacceptors_plus_5[i]));
      substringlist_gc(&(donors_minus_5[i]));
      substringlist_gc(&(antidonors_minus_5[i]));
      substringlist_gc(&(acceptors_minus_5[i]));
      substringlist_gc(&(antiacceptors_minus_5[i]));
    }
    FREEA(donors_plus_5);
    FREEA(antidonors_plus_5);
    FREEA(acceptors_plus_5);
    FREEA(antiacceptors_plus_5);
    FREEA(donors_minus_5);
    FREEA(antidonors_minus_5);
    FREEA(acceptors_minus_5);
    FREEA(antiacceptors_minus_5);
  }

  if (alloc3p == true) {
    /* Clean up 3 */
    for (i = 0; i <= max_splice_mismatches_3; i++) {
      substringlist_gc(&(donors_plus_3[i]));
      substringlist_gc(&(antidonors_plus_3[i]));
      substringlist_gc(&(acceptors_plus_3[i]));
      substringlist_gc(&(antiacceptors_plus_3[i]));
      substringlist_gc(&(donors_minus_3[i]));
      substringlist_gc(&(antidonors_minus_3[i]));
      substringlist_gc(&(acceptors_minus_3[i]));
      substringlist_gc(&(antiacceptors_minus_3[i]));
    }
    FREEA(donors_plus_3);
    FREEA(antidonors_plus_3);
    FREEA(acceptors_plus_3);
    FREEA(antiacceptors_plus_3);
    FREEA(donors_minus_3);
    FREEA(antidonors_minus_3);
    FREEA(acceptors_minus_3);
    FREEA(antiacceptors_minus_3);
  }
#endif


#if 0			 
  /* Unused code */
  /* This halfmapping appears to be a duplicate of the previous halfmapping */

  debug13(printf("found_terminals_p = %d\n",found_terminals_p));
  /* nconcordant might include a concordant pair of terminals */
  if (/* nconcordant == 0 && */ found_terminals_p == true && gmap_pairsearch_p == true) {
    /* 13.  GMAP terminal */
    /* Go ahead and resolve overlaps on each end by Stage3end, since
       we cannot do it by Stage3pair, but do not apply optimal
       score */
#if 1
    debug13(printf("Before remove_overlaps of 5' at cutoff level %d: %d hits\n",*cutoff_level_5,List_length(*hits5)));
    *hits5 = Stage3end_sort_bymatches(Stage3end_remove_overlaps(*hits5,/*finalp*/false));
    debug13(printf("After remove_overlaps (5): %d\n",List_length(*hits5)));
    
    debug13(printf("Before remove_overlaps of 3' at cutoff level %d: %d hits\n",*cutoff_level_3,List_length(*hits3)));
    *hits3 = Stage3end_sort_bymatches(Stage3end_remove_overlaps(*hits3,/*finalp*/false));
    debug13(printf("After remove_overlaps (6): %d\n",List_length(*hits3)));
#else
    /* Focus on those terminals not yet processed */
    terminals5 = Stage3end_sort_by_paired_seenp(terminals5);
    terminals3 = Stage3end_sort_by_paired_seenp(terminals3);
#endif

    i = 0;
    debug13(printf("%d hits on 5' end (vs max_gmap_pairsearch %d)\n",List_length(*hits5),max_gmap_pairsearch));
    debug13(printf("For each hit, running GMAP on 3' end to match with 5' ends\n"));
    *hits5 = Stage3end_sort_bymatches(*hits5);
    for (p = *hits5; p != NULL && i < max_gmap_pairsearch; p = List_next(p)) {
      hit5 = (Stage3end_T) List_head(p);
      debug13(printf("#%d/%d with nmatches %d\n",i,max_gmap_pairsearch,Stage3end_nmatches_posttrim(hit5)));
      halfmapping3 = align_halfmapping_with_gmap(gmap_history_3,hit5,/*hit3*/NULL,queryseq5,queryseq3,
						 queryuc_ptr_3,/*querylength*/querylength3,query3_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
						 queryrc3,Shortread_invertedp(queryseq3),
#endif
						 this3->plus_segments,this3->plus_nsegments,this3->minus_segments,this3->minus_nsegments,
						 oligoindices_major,oligoindices_minor,
						 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						 pairmax,shortsplicedist,genestrand);
      
      for (a = halfmapping3; a != NULL; a = List_next(a)) {
	gmap3 = (Stage3end_T) List_head(a);
	debug13(printf("=> Successful terminal GMAP on hit3 %p with score %d and nmatches %d.  Copying hit5 %p\n",
		       gmap3,Stage3end_score(gmap3),Stage3end_nmatches_posttrim(gmap3),hit5));
	if (Stage3end_score(gmap3) > *cutoff_level_3 + gmap_allowance) {
	  debug13(printf("Score is only %d vs cutoff level %d\n",Stage3end_score(gmap3),*cutoff_level_3));
	  Stage3end_free(&gmap3);
	} else if ((newpair = Stage3pair_new(Stage3end_copy(hit5),gmap3,genestrand,/*pairtype*/CONCORDANT,
					     /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/true)) == NULL) {
	  /* Stage3end_free(&gmap3); -- done by Stage3pair_new */
#if 0
	} else if (Stage3end_trimlength(hit5) < reject_trimlength) {
	  /* Save hit5-gmap3 */
	  hitpairs = List_push(hitpairs,(void *) newpair);
#endif
	} else {
	  /* Stage3end_free(&gmap3); */
	  Stage3pair_free(&newpair);
	}
      }
      List_free(&halfmapping3);
      i++;
    }

    i = 0;
    debug13(printf("%d hits on 3' end (vs max_gmap_pairsearch %d)\n",List_length(*hits3),max_gmap_pairsearch));
    debug13(printf("For each hit, running GMAP on 5' end to match with 3' ends\n"));
    *hits3 = Stage3end_sort_bymatches(*hits3);
    for (p = *hits3; p != NULL && i < max_gmap_pairsearch; p = List_next(p)) {
      hit3 = (Stage3end_T) List_head(p);
      debug13(printf("#%d/%d with nmatches %d\n",i,max_gmap_pairsearch,Stage3end_nmatches_posttrim(hit3)));
      halfmapping5 = align_halfmapping_with_gmap(gmap_history_5,/*hit5*/NULL,hit3,queryseq5,queryseq3,
						 queryuc_ptr_5,/*querylength*/querylength5,query5_lastpos,
#ifdef END_KNOWNSPLICING_SHORTCUT
						 queryrc5,Shortread_invertedp(queryseq5),
#endif
						 this5->plus_segments,this5->plus_nsegments,this5->minus_segments,this5->minus_nsegments,
						 oligoindices_major,oligoindices_minor,
						 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
						 pairmax,shortsplicedist,genestrand);
      for (a = halfmapping5; a != NULL; a = List_next(a)) {
	gmap5 = (Stage3end_T) List_head(a);
	debug13(printf("=> Successful terminal GMAP on hit5 %p with score %d and nmatches %d.  Copying hit3 %p\n",
		       hit5,Stage3end_score(gmap5),Stage3end_nmatches_posttrim(gmap5),hit3));
	if (Stage3end_score(gmap5) > *cutoff_level_5 + gmap_allowance) {
	  debug13(printf("Score is only %d vs cutoff level %d\n",Stage3end_score(gmap5),*cutoff_level_5));
	  Stage3end_free(&gmap5);

	} else if ((newpair = Stage3pair_new(gmap5,Stage3end_copy(hit3),genestrand,/*pairtype*/CONCORDANT,
					     /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/true)) == NULL) {
	  /* Stage3end_free(&gmap5); -- done by Stage3pair_new */
#if 0
	} else if (Stage3end_trimlength(hit3) < reject_trimlength) {
	  /* Save gmap5-hit3 */
	  hitpairs = List_push(hitpairs,(void *) newpair);
#endif
	} else {
	  /* Stage3end_free(&gmap5); */
	  Stage3pair_free(&newpair);
	}
      }
      List_free(&halfmapping5);
      i++;
    }

    debug(printf("13> After GMAP terminals, found %d concordant\n",nconcordant));
  }
#endif



  if (alloc_floors_p_5 == true) {
    Floors_free(&floors5);
  }
  if (alloc_floors_p_3 == true) {
    Floors_free(&floors3);
  }

  debug(printf("Ending with %d hitpairs, %d samechr, %d conc_transloc, %d terminals\n",
	       List_length(hitpairs),List_length(*samechr),List_length(*conc_transloc),List_length(*terminals)));

  hitpairs = Stage3pair_remove_circular_alias(hitpairs);
#if 0
  hitpairs = Stage3pair_remove_overlaps(hitpairs,/*translocp*/false,/*finalp*/true);
#endif

  FREE(plus_anchor_segments_5);
  FREE(minus_anchor_segments_5);
  FREE(plus_anchor_segments_3);
  FREE(minus_anchor_segments_3);

  return hitpairs;
}


static Pairtype_T
choose_among_paired (int *best_nmatches_paired, int *best_nmatches_5, int *best_nmatches_3,
		     List_T hitpairs, List_T samechr, List_T conc_transloc, List_T terminals) {
  Pairtype_T final_pairtype = UNPAIRED;
  List_T p;
  Stage3pair_T hitpair;
  int nmatches, nmatches5, nmatches3;

  debug16(printf("choose: %d hitpairs, %d conc_transloc, %d samechr, %d terminals\n",
		 List_length(hitpairs),List_length(conc_transloc),List_length(samechr),List_length(terminals)));

  *best_nmatches_paired = 0;
  for (p = hitpairs; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if ((nmatches = Stage3pair_nmatches_posttrim(&nmatches5,&nmatches3,hitpair)) > *best_nmatches_paired) {
      final_pairtype = CONCORDANT;
      *best_nmatches_paired = nmatches;
      *best_nmatches_5 = nmatches5;
      *best_nmatches_3 = nmatches3;
    }
  }

  if (hitpairs == NULL) {
    for (p = terminals; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((nmatches = Stage3pair_nmatches_posttrim(&nmatches5,&nmatches3,hitpair)) > *best_nmatches_paired) {
	final_pairtype = PAIRED_TERMINALS;
	*best_nmatches_paired = nmatches;
	*best_nmatches_5 = nmatches5;
	*best_nmatches_3 = nmatches3;
      }
    }

    *best_nmatches_paired += 1; /* penalty for choosing translocation over others */

    for (p = conc_transloc; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((nmatches = Stage3pair_nmatches_posttrim(&nmatches5,&nmatches3,hitpair)) > *best_nmatches_paired) {
	final_pairtype = CONCORDANT_TRANSLOCATIONS;
	*best_nmatches_paired = nmatches;
	*best_nmatches_5 = nmatches5;
	*best_nmatches_3 = nmatches3;
      }
    }

    for (p = samechr; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((nmatches = Stage3pair_nmatches_posttrim(&nmatches5,&nmatches3,hitpair)) > *best_nmatches_paired) {
	final_pairtype = PAIRED_UNSPECIFIED;
	*best_nmatches_paired = nmatches;
	*best_nmatches_5 = nmatches5;
	*best_nmatches_3 = nmatches3;
      }
    }
  }

  debug16(printf("best_nmatches_paired among paired = %d = %d + %d\n",
		 *best_nmatches_paired,*best_nmatches_5,*best_nmatches_3));
  debug16(printf("final pairtype: %s\n",Pairtype_string(final_pairtype)));

  return final_pairtype;
}


static int
best_nmatches_singleend (List_T hits) {
  int best_nmatches = 0;
  List_T p;
  Stage3end_T hit;
  int nmatches;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (Stage3end_T) p->first;
    if ((nmatches = Stage3end_nmatches_posttrim(hit)) > best_nmatches) {
      best_nmatches = nmatches;
    }
  }

  return best_nmatches;
}



/* Clean up all previous calculations */
static void
paired_results_free (T this5, T this3, List_T hitpairs, List_T samechr, List_T conc_transloc, List_T terminals,
		     List_T hits5, List_T hits3, int querylength5, int querylength3) {
  List_T p;
  Stage3pair_T stage3pair;

  for (p = hitpairs; p != NULL; p = List_next(p)) {
    stage3pair = (Stage3pair_T) List_head(p);
    Stage3pair_free(&stage3pair);
  }
  List_free(&hitpairs);

  for (p = samechr; p != NULL; p = List_next(p)) {
    stage3pair = (Stage3pair_T) List_head(p);
    Stage3pair_free(&stage3pair);
  }
  List_free(&samechr);

  for (p = conc_transloc; p != NULL; p = List_next(p)) {
    stage3pair = (Stage3pair_T) List_head(p);
    Stage3pair_free(&stage3pair);
  }
  List_free(&conc_transloc);

  for (p = terminals; p != NULL; p = List_next(p)) {
    stage3pair = (Stage3pair_T) List_head(p);
    Stage3pair_free(&stage3pair);
  }
  List_free(&conc_transloc);

  stage3list_gc(&hits3);
  stage3list_gc(&hits5);
  Stage1_free(&this3,querylength3);
  Stage1_free(&this5,querylength5);

  return;
}

static void
realign_separately (Stage3end_T **stage3array5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
		    Stage3end_T **stage3array3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
		    T this5, T this3, Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		    Shortread_T queryseq5, char *queryuc_ptr_5, char *queryrc5, char *quality_string_5, int querylength5, int query5_lastpos,
		    Shortread_T queryseq3, char *queryuc_ptr_3, char *queryrc3, char *quality_string_3, int querylength3, int query3_lastpos,
		    Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev, int indexdb_size_threshold,
		    Floors_T *floors_array,
		    int user_maxlevel_5, int user_maxlevel_3, int min_coverage_5, int min_coverage_3,
		    int indel_penalty_middle, int indel_penalty_end,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
		    Oligoindex_array_T oligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    bool keep_floors_p, int genestrand) {
  List_T singlehits5, singlehits3;
  int cutoff_level_5, cutoff_level_3;
  bool allvalidp5, allvalidp3;

  /* Re-align 5' end as a single end */
  if (read_oligos(&allvalidp5,this5,queryuc_ptr_5,querylength5,query5_lastpos,genestrand) == 0) {
    debug(printf("Aborting because no hits found anywhere\n"));
    singlehits5 = (List_T) NULL;
  } else {
    singlehits5 = align_end(&cutoff_level_5,this5,query5_compress_fwd,query5_compress_rev,
			    queryuc_ptr_5,queryrc5,querylength5,query5_lastpos,
			    indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
			    oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
			    user_maxlevel_5,indel_penalty_middle,indel_penalty_end,
			    localsplicing_penalty,distantsplicing_penalty,min_shortend,
			    allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			    keep_floors_p,genestrand,/*first_read_p*/true);
  }

  singlehits5 = Stage3end_filter_coverage(singlehits5,min_coverage_5);
  if (singlehits5 == NULL) {
    *nhits5_primary = *nhits5_altloc = 0;
    *stage3array5 = (Stage3end_T *) NULL;
  } else {
    Stage3end_count_hits(&(*nhits5_primary),&(*nhits5_altloc),singlehits5);
    *stage3array5 = (Stage3end_T *) List_to_array_out(singlehits5,NULL); List_free(&singlehits5); /* Return value */
    *stage3array5 = Stage3end_eval_and_sort((*nhits5_primary) + (*nhits5_altloc),
					    &(*first_absmq5),&(*second_absmq5),*stage3array5,maxpaths_search,
					    queryseq5,queryuc_ptr_5,query5_compress_fwd,query5_compress_rev,
					    quality_string_5,/*displayp*/true);
  }

  /* Re-align 3' end as a single end */
  if (read_oligos(&allvalidp3,this3,queryuc_ptr_3,querylength3,query3_lastpos,genestrand) == 0) {
    debug(printf("Aborting because no hits found anywhere\n"));
    singlehits3 = (List_T) NULL;
  } else {
    singlehits3 = align_end(&cutoff_level_3,this3,query3_compress_fwd,query3_compress_rev,
			    queryuc_ptr_3,queryrc3,querylength3,query3_lastpos,
			    indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
			    oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
			    user_maxlevel_3,indel_penalty_middle,indel_penalty_end,
			    localsplicing_penalty,distantsplicing_penalty,min_shortend,
			    allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			    keep_floors_p,genestrand,/*first_read_p*/false);
  }

  singlehits3 = Stage3end_filter_coverage(singlehits3,min_coverage_3);
  if (singlehits3 == NULL) {
    *nhits3_primary = *nhits3_altloc = 0;
    *stage3array3 = (Stage3end_T *) NULL;
  } else {
    Stage3end_count_hits(&(*nhits3_primary),&(*nhits3_altloc),singlehits3);
    *stage3array3 = (Stage3end_T *) List_to_array_out(singlehits3,NULL); List_free(&singlehits3); /* Return value */
    *stage3array3 = Stage3end_eval_and_sort((*nhits3_primary) + (*nhits3_altloc),
					    &(*first_absmq3),&(*second_absmq3),*stage3array3,maxpaths_search,
					    queryseq3,queryuc_ptr_3,query3_compress_fwd,query3_compress_rev,
					    quality_string_3,/*displayp*/true);
  }

  return;
}


/* Have four lists: hitpairs + terminals, samechr, and conc_transloc => result */
static Stage3pair_T *
consolidate_paired_results (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
			    Stage3end_T **stage3array5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
			    Stage3end_T **stage3array3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
			    List_T hitpairs, List_T samechr, List_T conc_transloc, List_T terminals,
			    List_T hits5, List_T hits3,
			    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

			    Shortread_T queryseq5, char *queryuc_ptr_5,
			    char *quality_string_5, int querylength5,
			    Shortread_T queryseq3, char *queryuc_ptr_3,
			    char *quality_string_3, int querylength3,

			    int cutoff_level_5, int cutoff_level_3, int min_coverage_5, int min_coverage_3,
			    Oligoindex_array_T oligoindices_minor,
			    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  Stage3pair_T *stage3pairarray, stage3pair, newpair;
  Stage3end_T hit5, hit3;
  List_T result, singlehits5, singlehits3, p;
  Pairtype_T pairtype;
  double percentage_paired, percentage_single;
  int best_nmatches_paired, best_nmatches_paired_5, best_nmatches_paired_3, best_nmatches_5, best_nmatches_3;

  
  debug16(printf("Entered consolidate_paired_results.  Passing pointer %p\n",&best_nmatches_paired));
  *final_pairtype = choose_among_paired(&best_nmatches_paired,&best_nmatches_paired_5,&best_nmatches_paired_3,
					hitpairs,samechr,conc_transloc,terminals);

  if (*final_pairtype == CONCORDANT) {
    /* Have concordant results */
    debug16(printf("Have %d concordant results\n",List_length(hitpairs)));
    for (p = samechr; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&samechr);

    for (p = conc_transloc; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&conc_transloc);
	  
    hitpairs = List_append(hitpairs,terminals);
    if (novelsplicingp || knownsplicingp) {
      hitpairs = Stage3pair_remove_excess_terminals(hitpairs);
    }

    if (gmap_improvement_p == false) {
      debug16(printf("No GMAP improvement: Before removing overlaps, %d results\n",List_length(hitpairs)));
      result = Stage3pair_optimal_score(hitpairs,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);
      result = Stage3pair_resolve_multimapping(result);
      /* result = Stage3pair_sort_distance(result); */
      debug16(printf("After removing overlaps, %d results\n",List_length(result)));

    } else {
      debug16(printf("GMAP improvement: Before removing overlaps, %d results\n",List_length(hitpairs)));
      result = Stage3pair_optimal_score(hitpairs,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/false);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/false);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/false);
      result = Stage3pair_resolve_multimapping(result);
      /* result = Stage3pair_sort_distance(result); */
      debug16(printf("After removing overlaps, %d results\n",List_length(result)));

      result = align_pair_with_gmap(&(*final_pairtype),result,
				    queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
				    oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				    cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/true);
      if (Stage3pair_sense_consistent_p(result) == false) {
	debug16(printf("sense is not consistent\n"));
	result = align_pair_with_gmap(&(*final_pairtype),result,
				      queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
				      oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				      cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/true);
      }

      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);

      /* TODO: Resolve terminals by doing full GMAP, and then redo optimal_score */

      result = Stage3pair_resolve_multimapping(result);
    }

  } else if (*final_pairtype == PAIRED_UNSPECIFIED) {
    /* Have paired results */
    debug16(printf("Have paired unspecified\n"));
    for (p = hitpairs; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&hitpairs);

    for (p = terminals; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&terminals);

    for (p = conc_transloc; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&conc_transloc);

    if (gmap_improvement_p == false) {
      debug16(printf("No GMAP improvement: Before removing overlaps, %d results\n",List_length(samechr)));
      result = Stage3pair_optimal_score(samechr,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);
      result = Stage3pair_resolve_multimapping(result);
    } else {
      debug16(printf("GMAP improvement: Before removing overlaps, %d results\n",List_length(samechr)));
      result = Stage3pair_optimal_score(samechr,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/false);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/false);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/false);
      result = Stage3pair_resolve_multimapping(result);

      result = align_pair_with_gmap(&(*final_pairtype),result,
				    queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
				    oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				    cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/false);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);
      result = Stage3pair_resolve_multimapping(result);

      if (Stage3pair_concordantp(result) == true) {
	debug16(printf("Found remaining concordant solution, so removing non-concordant ones\n"));
	*final_pairtype = CONCORDANT;
	result = Stage3pair_filter_nonconcordant(result);
	debug16(printf("Concordant results: %d\n",List_length(result)));
      } else {
	*final_pairtype = PAIRED_UNSPECIFIED;
      }
    }

  } else if (*final_pairtype == PAIRED_TERMINALS) {
    /* Have paired results */
    debug16(printf("Have paired terminals: %d\n",List_length(terminals)));
    for (p = hitpairs; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&hitpairs);

    for (p = samechr; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&samechr);

    for (p = conc_transloc; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&conc_transloc);

    if (gmap_improvement_p == false) {
      debug16(printf("No GMAP improvement: Before removing overlaps, %d results\n",List_length(terminals)));
      result = Stage3pair_optimal_score(terminals,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);
      result = Stage3pair_resolve_multimapping(result);

    } else {
      debug16(printf("GMAP improvement: Before removing overlaps, %d results\n",List_length(terminals)));
      result = Stage3pair_optimal_score(terminals,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/false);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/false);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/false);
      result = Stage3pair_resolve_multimapping(result);

      result = align_pair_with_gmap(&(*final_pairtype),result,
				    queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
				    oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				    cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/false);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);
      result = Stage3pair_resolve_multimapping(result);

#if 0
      debug16(printf("Assigning pairtypes\n"));
      *final_pairtype = PAIRED_TERMINALS;
      for (p = result; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	if (Stage3pair_assign_pairtype(stage3pair) == CONCORDANT) {
	  *final_pairtype = CONCORDANT;
	}
      }
      if (*final_pairtype == CONCORDANT) {
	debug16(printf("Found remaining concordant solution, so removing non-concordant ones\n"));
	result = Stage3pair_filter_nonconcordant(result);
	debug16(printf("Concordant results: %d\n",List_length(result)));
      }
#endif
    }

    *final_pairtype = CONCORDANT; /* Because paired terminals are each concordant */


  } else if (*final_pairtype == CONCORDANT_TRANSLOCATIONS) {
    debug16(printf("Have %d concordant translocation results\n",List_length(conc_transloc)));
    for (p = hitpairs; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&hitpairs);

    for (p = terminals; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&terminals);

    for (p = samechr; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&samechr);

    result = Stage3pair_optimal_score(conc_transloc,
				      query5_compress_fwd,query5_compress_rev,
				      query3_compress_fwd,query3_compress_rev,
				      querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/true);
    result = Stage3pair_remove_overlaps(result,/*translocp*/true,/*finalp*/true);
    result = Stage3pair_optimal_score(result,
				      query5_compress_fwd,query5_compress_rev,
				      query3_compress_fwd,query3_compress_rev,
				      querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);
    result = Stage3pair_resolve_multimapping(result);
    debug16(printf("Finally, have %d concordant translocation results\n",List_length(result)));

  } else if (*final_pairtype == CONCORDANT_TERMINAL) {
    debug16(printf("Have concordant terminal results.  Need to check against halfmapping\n"));

    /* Combine evidence */
    percentage_paired = 1.0 - (1.0 - (double) best_nmatches_paired_5/(double) querylength5) * (1.0 - (double) best_nmatches_paired_3/(double) querylength3);

    best_nmatches_5 = best_nmatches_singleend(hits5);
    best_nmatches_3 = best_nmatches_singleend(hits3);
    percentage_single = (double) best_nmatches_5/(double) querylength5;
    if ((double) best_nmatches_3/(double) querylength3 > percentage_single) {
      percentage_single = (double) best_nmatches_3/(double) querylength3;
    }

    debug16(printf("best_nmatches_5 = %d, best_nmatches_3 = %d\n",best_nmatches_5,best_nmatches_3));
    debug16(printf("=> misses_paired = %d, misses_5 = %d, misses_3 = %d\n",
		   querylength5 + querylength3 - best_nmatches_paired,querylength5 - best_nmatches_5,querylength3 - best_nmatches_3));
    debug16(printf("=> match_percentage paired = %f + %f => %f, percentage_5 = %f, percentage_3 = %f\n",
		   (double) best_nmatches_paired_5/(double) querylength5,
		   (double) best_nmatches_paired_3/(double) querylength3,
		   percentage_paired,(double) best_nmatches_5/(double) querylength5,(double) best_nmatches_3/(double) querylength3));

    if (percentage_single > percentage_paired) {
      debug16(printf("Unpaired is better than concordant_terminal\n"));
      for (p = hitpairs; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	Stage3pair_free(&stage3pair);
      }
      List_free(&hitpairs);

      for (p = samechr; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	Stage3pair_free(&stage3pair);
      }
      List_free(&samechr);

      for (p = conc_transloc; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	Stage3pair_free(&stage3pair);
      }
      List_free(&conc_transloc);

      for (p = terminals; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	Stage3pair_free(&stage3pair);
      }
      List_free(&terminals);

      result = (List_T) NULL;

    } else {
      debug16(printf("Concordant_terminal is better than unpaired\n"));
      *final_pairtype = CONCORDANT; /* CONCORDANT_TERMINAL used just to rank results behind translocations */

      for (p = hitpairs; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	Stage3pair_free(&stage3pair);
      }
      List_free(&hitpairs);

      for (p = samechr; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	Stage3pair_free(&stage3pair);
      }
      List_free(&samechr);

      for (p = conc_transloc; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	Stage3pair_free(&stage3pair);
      }
      List_free(&conc_transloc);

      for (p = terminals; p != NULL; p = List_next(p)) {
	stage3pair = (Stage3pair_T) List_head(p);
	Stage3pair_free(&stage3pair);
      }
      List_free(&terminals);

      if (gmap_improvement_p == false) {
	debug16(printf("No GMAP improvement: Before removing overlaps, %d results\n",List_length(result)));
	result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/true);
	result = Stage3pair_optimal_score(result,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);
	result = Stage3pair_resolve_multimapping(result);
	/* result = Stage3pair_sort_distance(result); */
	debug16(printf("After removing overlaps, %d results\n",List_length(result)));

      } else {
	debug16(printf("GMAP improvement: Before removing overlaps, %d results\n",List_length(result)));
	result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/false);
	result = Stage3pair_optimal_score(result,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/false);
	result = Stage3pair_resolve_multimapping(result);
	/* result = Stage3pair_sort_distance(result); */
	debug16(printf("After removing overlaps, %d results\n",List_length(result)));

	result = align_pair_with_gmap(&(*final_pairtype),result,
				      queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
				      oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				      cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/true);
	if (Stage3pair_sense_consistent_p(result) == false) {
	  debug16(printf("sense is not consistent\n"));
	  result = align_pair_with_gmap(&(*final_pairtype),result,
					queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
					oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
					cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/true);
	}

	result = Stage3pair_optimal_score(result,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  querylength5,querylength3,/*keep_gmap_p*/true,/*finalp*/true);
	result = Stage3pair_remove_overlaps(result,/*translocp*/false,/*finalp*/true);
	result = Stage3pair_optimal_score(result,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  querylength5,querylength3,/*keep_gmap_p*/false,/*finalp*/true);
	result = Stage3pair_resolve_multimapping(result);
      }
    }

  } else {
    debug16(printf("Have unpaired results\n"));
    /* Need to free conc_transloc, since we can get here with multiple results */
    for (p = conc_transloc; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    List_free(&conc_transloc);

    result = (List_T) NULL;
  }


  if (result == NULL) {
    singlehits5 = Stage3end_optimal_score(hits5,query5_compress_fwd,query5_compress_rev,
					  querylength5,/*keep_gmap_p*/true,/*finalp*/true);
    /* singlehits5 = Stage3end_reject_trimlengths(singlehits5); */
    singlehits5 = Stage3end_linearize_5(singlehits5);
    singlehits5 = Stage3end_remove_overlaps(singlehits5,/*finalp*/true);
    singlehits5 = Stage3end_optimal_score(singlehits5,query5_compress_fwd,query5_compress_rev,
					  querylength5,/*keep_gmap_p*/false,/*finalp*/true);
    singlehits5 = Stage3end_resolve_multimapping(singlehits5);

    singlehits3 = Stage3end_optimal_score(hits3,query3_compress_fwd,query3_compress_rev,
					  querylength3,/*keep_gmap_p*/true,/*finalp*/true);
    /* singlehits3 = Stage3end_reject_trimlengths(singlehits3); */
    singlehits3 = Stage3end_linearize_3(singlehits3);
    singlehits3 = Stage3end_remove_overlaps(singlehits3,/*finalp*/true);
    singlehits3 = Stage3end_optimal_score(singlehits3,query3_compress_fwd,query3_compress_rev,
					  querylength3,/*keep_gmap_p*/false,/*finalp*/true);
    singlehits3 = Stage3end_resolve_multimapping(singlehits3);

    debug16(printf("5' end has %d hits and 3' end has %d hits\n",
		   List_length(singlehits5),List_length(singlehits3)));

    if (List_length(singlehits5) == 1 && List_length(singlehits3) == 1 &&
	(pairtype = Stage3_determine_pairtype(hit5=(Stage3end_T) List_head(singlehits5),hit3=(Stage3end_T) List_head(singlehits3))) != UNPAIRED) {
      /* Convert unpaired uniq to a paired uniq */
      debug16(printf("Converting unpaired uniq to paired uniq, with initial pairtype %s\n",Pairtype_string(pairtype)));
      if ((newpair = Stage3pair_new(hit5,hit3,/*genestrand*/0,pairtype,
				    /*private5p*/false,/*private3p*/false,
				    /*expect_concordant_p*/pairtype == CONCORDANT ? true : false)) != NULL) {
	stage3pairarray = (Stage3pair_T *) CALLOC_OUT(1,sizeof(Stage3pair_T));
	stage3pairarray[0] = newpair;
	    
	*nhits5_primary = *nhits5_altloc = 0;
	*nhits3_primary = *nhits3_altloc = 0;
	*stage3array5 = *stage3array3 = (Stage3end_T *) NULL;
	    
	if (Stage3pair_altlocp(newpair) == true) {
	  *npaths_primary = 0;
	  *npaths_altloc = 1;
	} else {
	  *npaths_primary = 1;
	  *npaths_altloc = 0;
	}
	if (pairtype == CONCORDANT) {
	  debug16(printf("final pairtype is CONCORDANT\n"));
	  *final_pairtype = CONCORDANT;
	} else {
	  debug16(printf("final pairtype is PAIRED_UNSPECIFIED\n"));
	  *final_pairtype = PAIRED_UNSPECIFIED;
	}
	Stage3pair_privatize(stage3pairarray,/*npairs*/1);
	Stage3pair_eval_and_sort(/*npaths*/(*npaths_primary) + (*npaths_altloc),
				 &(*first_absmq),&(*second_absmq),stage3pairarray,maxpaths_search,queryseq5,queryseq3,
				 queryuc_ptr_5,queryuc_ptr_3,
				 query5_compress_fwd,query5_compress_rev,
				 query3_compress_fwd,query3_compress_rev,
				 quality_string_5,quality_string_3);
	    
	stage3list_gc(&singlehits3);
	stage3list_gc(&singlehits5);
	debug16(printf("1 Exiting consolidate_paired_results\n"));
	return stage3pairarray;
      }
    }

    /* Fall through: halfmapping or unpaired */
    *npaths_primary = *npaths_altloc = 0;
    *final_pairtype = UNPAIRED;
	  
    singlehits5 = Stage3end_filter_coverage(singlehits5,min_coverage_5);
    if (singlehits5 == NULL) {
      *nhits5_primary = *nhits5_altloc = 0;
      *stage3array5 = (Stage3end_T *) NULL;
    } else {
#if 0
      singlehits5 = Stage3end_unalias_circular(singlehits5);
#else
      singlehits5 = Stage3end_remove_circular_alias(singlehits5); /* Contains a call to unalias_circular */
      singlehits5 = Stage3end_remove_duplicates(singlehits5); /* Aliases can cause duplicates */
      Stage3end_count_hits(&(*nhits5_primary),&(*nhits5_altloc),singlehits5);
#endif
      *stage3array5 = (Stage3end_T *) List_to_array_out(singlehits5,NULL); List_free(&singlehits5); /* Return value */
    }

    singlehits3 = Stage3end_filter_coverage(singlehits3,min_coverage_3);
    if (singlehits3 == NULL) {
      *nhits3_primary = *nhits3_altloc = 0;
      *stage3array3 = (Stage3end_T *) NULL;
    } else {
#if 0
      singlehits3 = Stage3end_unalias_circular(singlehits3);
#else
      singlehits3 = Stage3end_remove_circular_alias(singlehits3); /* Contains a call to unalias_circular */
      singlehits3 = Stage3end_remove_duplicates(singlehits3); /* Aliases can cause duplicates */
      Stage3end_count_hits(&(*nhits3_primary),&(*nhits3_altloc),singlehits3);
#endif
      *stage3array3 = (Stage3end_T *) List_to_array_out(singlehits3,NULL); List_free(&singlehits3); /* Return value */
    }

    if ((*nhits5_primary) + (*nhits5_altloc) > 0) {
      if ((*nhits3_primary) + (*nhits3_altloc) == 1) {
      /* Use single 3' hit to guide sorting of multiple 5' hits */
        *stage3array5 = Stage3end_eval_and_sort_guided((*nhits5_primary) + (*nhits5_altloc),
   	                                               &(*first_absmq5),&(*second_absmq5),/*guide*/(*stage3array3)[0],
						       *stage3array5,maxpaths_search,queryseq5,
						       queryuc_ptr_5,query5_compress_fwd,query5_compress_rev,
						       quality_string_5,/*displayp*/true);
      } else {
        *stage3array5 = Stage3end_eval_and_sort((*nhits5_primary) + (*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
 						*stage3array5,maxpaths_search,queryseq5,
						queryuc_ptr_5,query5_compress_fwd,query5_compress_rev,
						quality_string_5,/*displayp*/true);
      }
    }

    if ((*nhits3_primary) + (*nhits3_altloc) > 0) {
      if ((*nhits5_primary) + (*nhits5_altloc) == 1) {
	/* Use single 5' hit to guide sorting of multiple 3' hits */
        *stage3array3 = Stage3end_eval_and_sort_guided((*nhits3_primary) + (*nhits3_altloc),
                                                       &(*first_absmq3),&(*second_absmq3),/*guide*/(*stage3array5)[0],
						       *stage3array3,maxpaths_search,queryseq3,
						       queryuc_ptr_3,query3_compress_fwd,query3_compress_rev,
						       quality_string_3,/*displayp*/true);
      } else {
        *stage3array3 = Stage3end_eval_and_sort((*nhits3_primary) + (*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
						*stage3array3,maxpaths_search,queryseq3,
						queryuc_ptr_3,query3_compress_fwd,query3_compress_rev,
						quality_string_3,/*displayp*/true);
      }
    }
    debug16(printf("Result is NULL, and we have %d hits on 5' end and %d hits on 3' end\n",
		   (*nhits5_primary) + (*nhits5_altloc),(*nhits3_primary) + (*nhits3_altloc)));
    debug16(printf("2 Exiting consolidate_paired_results\n"));
    return (Stage3pair_T *) NULL;

  } else {
    debug16(printf("Result is not NULL (%d paths), and we fall through to concordant, paired, or transloc pairs\n",
		   List_length(result)));

    result = Stage3pair_filter_coverage(result,min_coverage_5,min_coverage_3);
    if (result == NULL) {
      *npaths_primary = *npaths_altloc = 0;
      stage3list_gc(&hits3);
      stage3list_gc(&hits5);
      *nhits5_primary = *nhits5_altloc = 0;
      *nhits3_primary = *nhits3_altloc = 0;
      *stage3array5 = *stage3array3 = (Stage3end_T *) NULL;
      debug16(printf("3 Exiting consolidate_paired_results\n"));
      return (Stage3pair_T *) NULL;

    } else {
      /* result != NULL */
      /* Concordant, paired, or transloc pairs found.  Remove single hits. */
      Stage3pair_count_hits(&(*npaths_primary),&(*npaths_altloc),result);
      stage3pairarray = (Stage3pair_T *) List_to_array_out(result,NULL); List_free(&result); /* Return value */
      Stage3pair_privatize(stage3pairarray,/*npaths*/(*npaths_primary) + (*npaths_altloc));
      Stage3pair_eval_and_sort(/*npaths*/(*npaths_primary) + (*npaths_altloc),
			       &(*first_absmq),&(*second_absmq),stage3pairarray,maxpaths_search,queryseq5,queryseq3,
			       queryuc_ptr_5,queryuc_ptr_3,
			       query5_compress_fwd,query5_compress_rev,
			       query3_compress_fwd,query3_compress_rev,
			       quality_string_5,quality_string_3);
      stage3list_gc(&hits3);
      stage3list_gc(&hits5);

      *nhits5_primary = *nhits5_altloc = 0;
      *nhits3_primary = *nhits3_altloc = 0;
      *stage3array5 = *stage3array3 = (Stage3end_T *) NULL;
      debug16(printf("4 Exiting consolidate_paired_results\n"));

      return stage3pairarray;
    }
  }
}


static Stage3pair_T *
paired_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
	     Stage3end_T **stage3array5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
	     Stage3end_T **stage3array3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
	     Shortread_T queryseq5, Shortread_T queryseq3,
	     Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev, int indexdb_size_threshold,
	     Floors_T *floors_array,
	     double user_maxlevel_float, double user_mincoverage_float, int indel_penalty_middle, int indel_penalty_end,
	     bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
	     int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
	     Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
	     Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
	     Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	     Chrpos_T pairmax, bool keep_floors_p) {
  Stage3pair_T *stage3pairarray;
  List_T hitpairs = NULL, samechr = NULL, conc_transloc = NULL, terminals = NULL, hits5 = NULL, hits3 = NULL;
  T this5, this3;
  char *queryuc_ptr_5, *queryuc_ptr_3, *quality_string_5, *quality_string_3;
  Compress_T query5_compress_fwd = NULL, query5_compress_rev = NULL, query3_compress_fwd = NULL, query3_compress_rev = NULL;
  History_T gmap_history_5, gmap_history_3;
  int user_maxlevel_5, user_maxlevel_3, min_coverage_5, min_coverage_3;
  int found_score, cutoff_level_5, cutoff_level_3;
  int querylength5, querylength3, query5_lastpos, query3_lastpos;
#if 0
  int maxpairedpaths = 10*maxpaths; /* For computation, not for printing. */
#else
  int maxpairedpaths = maxpaths_search;  /* 100000 */
#endif
  bool abort_pairing_p;
  char *queryrc5, *queryrc3;

  querylength5 = Shortread_fulllength(queryseq5);
  querylength3 = Shortread_fulllength(queryseq3);

#ifdef HAVE_ALLOCA
  if (querylength5 <= MAX_STACK_READLENGTH) {
    queryrc5 = (char *) ALLOCA((querylength5+1)*sizeof(char));
  } else {
    queryrc5 = (char *) MALLOC((querylength5+1)*sizeof(char));
  }
  if (querylength3 <= MAX_STACK_READLENGTH) {
    queryrc3 = (char *) ALLOCA((querylength3+1)*sizeof(char));
  } else {
    queryrc3 = (char *) MALLOC((querylength3+1)*sizeof(char));
  }
#else
  queryrc5 = (char *) MALLOC((querylength5+1)*sizeof(char));
  queryrc3 = (char *) MALLOC((querylength3+1)*sizeof(char));
#endif

  if (user_maxlevel_float < 0.0) {
    user_maxlevel_5 = user_maxlevel_3 = -1;
  } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
    user_maxlevel_5 = (int) rint(user_maxlevel_float * (double) querylength5);
    user_maxlevel_3 = (int) rint(user_maxlevel_float * (double) querylength3);
  } else {
    user_maxlevel_5 = user_maxlevel_3 = (int) user_maxlevel_float;
  }

  if (user_mincoverage_float < 0.0) {
    min_coverage_5 = min_coverage_3 = 0;
  } else if (user_mincoverage_float > 0.0 && user_mincoverage_float < 1.0) {
    min_coverage_5 = (int) rint(user_mincoverage_float * (double) querylength5);
    min_coverage_3 = (int) rint(user_mincoverage_float * (double) querylength3);
  } else {
    min_coverage_5 = min_coverage_3 = (int) user_mincoverage_float;
  }


  this5 = Stage1_new(querylength5);
  this3 = Stage1_new(querylength3);
  queryuc_ptr_5 = Shortread_fullpointer_uc(queryseq5);
  queryuc_ptr_3 = Shortread_fullpointer_uc(queryseq3);
  quality_string_5 = Shortread_quality_string(queryseq5);
  quality_string_3 = Shortread_quality_string(queryseq3);
  query5_lastpos = querylength5 - index1part;
  query3_lastpos = querylength3 - index1part;

  /* Limit search on repetitive sequences */
  if (check_dinucleotides(queryuc_ptr_5,querylength5) == false) {
    user_maxlevel_5 = 0;
  }
  if (check_dinucleotides(queryuc_ptr_3,querylength3) == false) {
    user_maxlevel_3 = 0;
  }

  query5_compress_fwd = Compress_new_fwd(queryuc_ptr_5,querylength5);
  query5_compress_rev = Compress_new_rev(queryuc_ptr_5,querylength5);
  query3_compress_fwd = Compress_new_fwd(queryuc_ptr_3,querylength3);
  query3_compress_rev = Compress_new_rev(queryuc_ptr_3,querylength3);
  gmap_history_5 = History_new();
  gmap_history_3 = History_new();
  make_complement_buffered(queryrc5,queryuc_ptr_5,querylength5);
  make_complement_buffered(queryrc3,queryuc_ptr_3,querylength3);

  hitpairs = align_pair(&abort_pairing_p,&found_score,&cutoff_level_5,&cutoff_level_3,
			&samechr,&conc_transloc,&terminals,gmap_history_5,gmap_history_3,
			&hits5,&hits3,this5,this3,query5_compress_fwd,query5_compress_rev,
			query3_compress_fwd,query3_compress_rev,
			queryuc_ptr_5,queryuc_ptr_3,
#ifndef LARGE_GENOMES
			queryrc5,queryrc3,
#endif
			querylength5,querylength3,query5_lastpos,query3_lastpos,
			indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,

			oligoindices_major,oligoindices_minor,
			pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,

			user_maxlevel_5,user_maxlevel_3,indel_penalty_middle,indel_penalty_end,
			localsplicing_penalty,distantsplicing_penalty,
			allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			pairmax,maxpairedpaths,keep_floors_p,queryseq5,queryseq3,/*genestrand*/0);

  if (abort_pairing_p == true) {
    debug16(printf("abort_pairing_p is true\n"));
    paired_results_free(this5,this3,hitpairs,samechr,conc_transloc,terminals,
			hits5,hits3,querylength5,querylength3);

    this5 = Stage1_new(querylength5);
    this3 = Stage1_new(querylength3);
    realign_separately(stage3array5,&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
		       stage3array3,&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
		       this5,this3,
		       query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
		       queryseq5,queryuc_ptr_5,queryrc5,quality_string_5,querylength5,query5_lastpos,
		       queryseq3,queryuc_ptr_3,queryrc3,quality_string_3,querylength3,query3_lastpos,
		       indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
		       user_maxlevel_5,user_maxlevel_3,min_coverage_5,min_coverage_3,
		       indel_penalty_middle,indel_penalty_end,
		       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
		       localsplicing_penalty,distantsplicing_penalty,min_shortend,
		       oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
		       keep_floors_p,/*genestrand*/0);

    *npaths_primary = *npaths_altloc = 0;
    *final_pairtype = UNPAIRED;
    History_free(&gmap_history_3);
    History_free(&gmap_history_5);
    Compress_free(&query5_compress_fwd);
    Compress_free(&query5_compress_rev);
    Compress_free(&query3_compress_fwd);
    Compress_free(&query3_compress_rev);
    Stage1_free(&this5,querylength5);
    Stage1_free(&this3,querylength3);

    stage3pairarray = (Stage3pair_T *) NULL;

  } else {
    stage3pairarray =
      consolidate_paired_results(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),&(*final_pairtype),
				 &(*stage3array5),&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
				 &(*stage3array3),&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
				 hitpairs,samechr,conc_transloc,terminals,hits5,hits3,
				 query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				 queryseq5,queryuc_ptr_5,quality_string_5,querylength5,
				 queryseq3,queryuc_ptr_3,quality_string_3,querylength3,
				 cutoff_level_5,cutoff_level_3,min_coverage_5,min_coverage_3,
				 oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);

    History_free(&gmap_history_3);
    History_free(&gmap_history_5);
    Compress_free(&query5_compress_fwd);
    Compress_free(&query5_compress_rev);
    Compress_free(&query3_compress_fwd);
    Compress_free(&query3_compress_rev);
    Stage1_free(&this5,querylength5);
    Stage1_free(&this3,querylength3);
  }

#ifdef HAVE_ALLOCA
  if (querylength5 <= MAX_STACK_READLENGTH) {
    FREEA(queryrc5);
  } else {
    FREE(queryrc5);
  }
  if (querylength3 <= MAX_STACK_READLENGTH) {
    FREEA(queryrc3);
  } else {
    FREE(queryrc3);
  }
#else
  FREE(queryrc5);
  FREE(queryrc3);
#endif

  return stage3pairarray;
}


static Stage3pair_T *
paired_read_tolerant_nonstranded (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
				  Stage3end_T **stage3array5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
				  Stage3end_T **stage3array3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
				  Shortread_T queryseq5, Shortread_T queryseq3,
				  Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev, int indexdb_size_threshold,
				  Floors_T *floors_array,
				  double user_maxlevel_float, double user_mincoverage_float, int indel_penalty_middle, int indel_penalty_end,
				  bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
				  int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
				  Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
				  Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
				  Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				  Chrpos_T pairmax, bool keep_floors_p) {
  Stage3pair_T *stage3pairarray;
  List_T hitpairs, hitpairs_geneplus = NULL, hitpairs_geneminus = NULL;
  List_T samechr, samechr_geneplus = NULL, samechr_geneminus = NULL;
  List_T conc_transloc, conc_transloc_geneplus = NULL, conc_transloc_geneminus = NULL;
  List_T terminals, terminals_geneplus = NULL, terminals_geneminus = NULL;
  List_T hits5, hits3, hits_geneplus_5 = NULL, hits_geneplus_3 = NULL, hits_geneminus_5 = NULL, hits_geneminus_3 = NULL;
  T this_geneplus_5, this_geneplus_3, this_geneminus_5, this_geneminus_3;
  char *queryuc_ptr_5, *queryuc_ptr_3, *quality_string_5, *quality_string_3;
  Compress_T query5_compress_fwd = NULL, query5_compress_rev = NULL, query3_compress_fwd = NULL, query3_compress_rev = NULL;
  History_T gmap_history_5, gmap_history_3;
  int user_maxlevel_5, user_maxlevel_3, min_coverage_5, min_coverage_3;
  int found_score_geneplus, found_score_geneminus;
  int cutoff_level_5, cutoff_level_3;
  int querylength5, querylength3, query5_lastpos, query3_lastpos;
#if 0
  int maxpairedpaths = 10*maxpaths; /* For computation, not for printing. */
#else
  int maxpairedpaths = maxpaths_search;  /* 100000 */
#endif
  bool abort_pairing_p_geneplus, abort_pairing_p_geneminus;
  char *queryrc5, *queryrc3;


  querylength5 = Shortread_fulllength(queryseq5);
  querylength3 = Shortread_fulllength(queryseq3);

#ifdef HAVE_ALLOCA
  if (querylength5 <= MAX_STACK_READLENGTH) {
    queryrc5 = (char *) ALLOCA((querylength5+1)*sizeof(char));
  } else {
    queryrc5 = (char *) MALLOC((querylength5+1)*sizeof(char));
  }
  if (querylength3 <= MAX_STACK_READLENGTH) {
    queryrc3 = (char *) ALLOCA((querylength3+1)*sizeof(char));
  } else {
    queryrc3 = (char *) MALLOC((querylength3+1)*sizeof(char));
  }
#else
  queryrc5 = (char *) MALLOC((querylength5+1)*sizeof(char));
  queryrc3 = (char *) MALLOC((querylength3+1)*sizeof(char));
#endif

  if (user_maxlevel_float < 0.0) {
    user_maxlevel_5 = user_maxlevel_3 = -1;
  } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
    user_maxlevel_5 = (int) rint(user_maxlevel_float * (double) querylength5);
    user_maxlevel_3 = (int) rint(user_maxlevel_float * (double) querylength3);
  } else {
    user_maxlevel_5 = user_maxlevel_3 = (int) user_maxlevel_float;
  }

  if (user_mincoverage_float < 0.0) {
    min_coverage_5 = min_coverage_3 = 0;
  } else if (user_mincoverage_float > 0.0 && user_mincoverage_float < 1.0) {
    min_coverage_5 = (int) rint(user_mincoverage_float * (double) querylength5);
    min_coverage_3 = (int) rint(user_mincoverage_float * (double) querylength3);
  } else {
    min_coverage_5 = min_coverage_3 = (int) user_mincoverage_float;
  }

  this_geneplus_5 = Stage1_new(querylength5);
  this_geneplus_3 = Stage1_new(querylength3);
  this_geneminus_5 = Stage1_new(querylength5);
  this_geneminus_3 = Stage1_new(querylength3);

  queryuc_ptr_5 = Shortread_fullpointer_uc(queryseq5);
  queryuc_ptr_3 = Shortread_fullpointer_uc(queryseq3);
  quality_string_5 = Shortread_quality_string(queryseq5);
  quality_string_3 = Shortread_quality_string(queryseq3);
  query5_lastpos = querylength5 - index1part;
  query3_lastpos = querylength3 - index1part;

  /* Limit search on repetitive sequences */
  if (check_dinucleotides(queryuc_ptr_5,querylength5) == false) {
    user_maxlevel_5 = 0;
  }
  if (check_dinucleotides(queryuc_ptr_3,querylength3) == false) {
    user_maxlevel_3 = 0;
  }

  query5_compress_fwd = Compress_new_fwd(queryuc_ptr_5,querylength5);
  query5_compress_rev = Compress_new_rev(queryuc_ptr_5,querylength5);
  query3_compress_fwd = Compress_new_fwd(queryuc_ptr_3,querylength3);
  query3_compress_rev = Compress_new_rev(queryuc_ptr_3,querylength3);
  gmap_history_5 = History_new();
  gmap_history_3 = History_new();
  make_complement_buffered(queryrc5,queryuc_ptr_5,querylength5);
  make_complement_buffered(queryrc3,queryuc_ptr_3,querylength3);

  abort_pairing_p_geneplus = false;
  hitpairs_geneplus = align_pair(&abort_pairing_p_geneplus,&found_score_geneplus,
				 &cutoff_level_5,&cutoff_level_3,
				 &samechr_geneplus,&conc_transloc_geneplus,&terminals_geneplus,
				 gmap_history_5,gmap_history_3,
				 &hits_geneplus_5,&hits_geneplus_3,this_geneplus_5,this_geneplus_3,
				 query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				 queryuc_ptr_5,queryuc_ptr_3,
#ifndef LARGE_GENOMES
				 queryrc5,queryrc3,
#endif
				 querylength5,querylength3,query5_lastpos,query3_lastpos,
				 indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
				 
				 oligoindices_major,oligoindices_minor,
				 pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				 
				 user_maxlevel_5,user_maxlevel_3,indel_penalty_middle,indel_penalty_end,
				 localsplicing_penalty,distantsplicing_penalty,
				 allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				 pairmax,maxpairedpaths,keep_floors_p,
				 queryseq5,queryseq3,/*genestrand*/+1);

  abort_pairing_p_geneminus = false;
  hitpairs_geneminus = align_pair(&abort_pairing_p_geneminus,&found_score_geneminus,
				  &cutoff_level_5,&cutoff_level_3,
				  &samechr_geneminus,&conc_transloc_geneminus,&terminals_geneminus,
				  gmap_history_5,gmap_history_3,
				  &hits_geneminus_5,&hits_geneminus_3,this_geneminus_5,this_geneminus_3,
				  query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				  queryuc_ptr_5,queryuc_ptr_3,
#ifndef LARGE_GENOMES
				  queryrc5,queryrc3,
#endif
				  querylength5,querylength3,query5_lastpos,query3_lastpos,
				  indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
				  
				  oligoindices_major,oligoindices_minor,
				  pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				  
				  user_maxlevel_5,user_maxlevel_3,indel_penalty_middle,indel_penalty_end,
				  localsplicing_penalty,distantsplicing_penalty,
				  allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				  pairmax,maxpairedpaths,keep_floors_p,queryseq5,queryseq3,/*genestrand*/+2);

  if (found_score_geneplus < found_score_geneminus) {
    paired_results_free(this_geneminus_5,this_geneminus_3,hitpairs_geneminus,samechr_geneminus,conc_transloc_geneminus,
			terminals_geneminus,hits_geneminus_5,hits_geneminus_3,querylength5,querylength3);

    if (abort_pairing_p_geneplus == true) {
      debug16(printf("abort_pairing_p_geneplus is true\n"));
      paired_results_free(this_geneplus_5,this_geneplus_3,hitpairs_geneplus,samechr_geneplus,conc_transloc_geneplus,
			  terminals_geneplus,hits_geneplus_5,hits_geneplus_3,querylength5,querylength3);

      this_geneplus_5 = Stage1_new(querylength5);
      this_geneplus_3 = Stage1_new(querylength3);
      realign_separately(stage3array5,&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
			 stage3array3,&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
			 this_geneplus_5,this_geneplus_3,
			 query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
			 queryseq5,queryuc_ptr_5,queryrc5,quality_string_5,querylength5,query5_lastpos,
			 queryseq3,queryuc_ptr_3,queryrc3,quality_string_3,querylength3,query3_lastpos,
			 indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
			 user_maxlevel_5,user_maxlevel_3,min_coverage_5,min_coverage_3,
			 indel_penalty_middle,indel_penalty_end,
			 allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			 localsplicing_penalty,distantsplicing_penalty,min_shortend,
			 oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
			 keep_floors_p,/*genestrand*/+1);

      *npaths_primary = *npaths_altloc = 0;
      *final_pairtype = UNPAIRED;
      History_free(&gmap_history_3);
      History_free(&gmap_history_5);
      Compress_free(&query5_compress_fwd);
      Compress_free(&query5_compress_rev);
      Compress_free(&query3_compress_fwd);
      Compress_free(&query3_compress_rev);
      Stage1_free(&this_geneplus_5,querylength5);
      Stage1_free(&this_geneplus_3,querylength3);
      
      stage3pairarray = (Stage3pair_T *) NULL;

    } else {
      stage3pairarray =
	consolidate_paired_results(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),&(*final_pairtype),
				   &(*stage3array5),&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
				   &(*stage3array3),&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
				   hitpairs_geneplus,samechr_geneplus,conc_transloc_geneplus,terminals_geneplus,
				   hits_geneplus_5,hits_geneplus_3,
				   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				   queryseq5,queryuc_ptr_5,quality_string_5,querylength5,
				   queryseq3,queryuc_ptr_3,quality_string_3,querylength3,
				   cutoff_level_5,cutoff_level_3,min_coverage_5,min_coverage_3,
				   oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
      History_free(&gmap_history_3);
      History_free(&gmap_history_5);
      Compress_free(&query5_compress_fwd);
      Compress_free(&query5_compress_rev);
      Compress_free(&query3_compress_fwd);
      Compress_free(&query3_compress_rev);
      Stage1_free(&this_geneplus_5,querylength5);
      Stage1_free(&this_geneplus_3,querylength3);
      /* return stage3pairarray; */
    }

  } else if (found_score_geneminus < found_score_geneplus) {
    paired_results_free(this_geneplus_5,this_geneplus_3,hitpairs_geneplus,samechr_geneplus,conc_transloc_geneplus,
			terminals_geneplus,hits_geneplus_5,hits_geneplus_3,querylength5,querylength3);

    if (abort_pairing_p_geneminus == true) {
      debug16(printf("abort_pairing_p_geneminus is true\n"));
      paired_results_free(this_geneminus_5,this_geneminus_3,hitpairs_geneminus,samechr_geneminus,conc_transloc_geneminus,
			  terminals_geneminus,hits_geneminus_5,hits_geneminus_3,querylength5,querylength3);

      this_geneminus_5 = Stage1_new(querylength5);
      this_geneminus_3 = Stage1_new(querylength3);
      realign_separately(stage3array5,&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
			 stage3array3,&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
			 this_geneminus_5,this_geneminus_3,
			 query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
			 queryseq5,queryuc_ptr_5,queryrc5,quality_string_5,querylength5,query5_lastpos,
			 queryseq3,queryuc_ptr_3,queryrc3,quality_string_3,querylength3,query3_lastpos,
			 indexdb_fwd,indexdb_rev,indexdb_size_threshold,floors_array,
			 user_maxlevel_5,user_maxlevel_3,min_coverage_5,min_coverage_3,
			 indel_penalty_middle,indel_penalty_end,
			 allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			 localsplicing_penalty,distantsplicing_penalty,min_shortend,
			 oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
			 keep_floors_p,/*genestrand*/+2);

      *npaths_primary = *npaths_altloc = 0;
      *final_pairtype = UNPAIRED;
      History_free(&gmap_history_3);
      History_free(&gmap_history_5);
      Compress_free(&query5_compress_fwd);
      Compress_free(&query5_compress_rev);
      Compress_free(&query3_compress_fwd);
      Compress_free(&query3_compress_rev);
      Stage1_free(&this_geneminus_5,querylength5);
      Stage1_free(&this_geneminus_3,querylength3);

      stage3pairarray = (Stage3pair_T *) NULL;

    } else {
      stage3pairarray =
	consolidate_paired_results(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),&(*final_pairtype),
				   &(*stage3array5),&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
				   &(*stage3array3),&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
				   hitpairs_geneminus,samechr_geneminus,conc_transloc_geneminus,terminals_geneminus,
				   hits_geneminus_5,hits_geneminus_3,
				   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				   queryseq5,queryuc_ptr_5,quality_string_5,querylength5,
				   queryseq3,queryuc_ptr_3,quality_string_3,querylength3,
				   cutoff_level_5,cutoff_level_3,min_coverage_5,min_coverage_3,
				   oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
      History_free(&gmap_history_3);
      History_free(&gmap_history_5);
      Compress_free(&query5_compress_fwd);
      Compress_free(&query5_compress_rev);
      Compress_free(&query3_compress_fwd);
      Compress_free(&query3_compress_rev);
      Stage1_free(&this_geneminus_5,querylength5);
      Stage1_free(&this_geneminus_3,querylength3);
      /* return stage3pairarray; */
    }

  } else {
    hitpairs = List_append(hitpairs_geneplus,hitpairs_geneminus);
    samechr = List_append(samechr_geneplus,samechr_geneminus);
    conc_transloc = List_append(conc_transloc_geneplus,conc_transloc_geneminus);
    terminals = List_append(terminals_geneplus,terminals_geneminus);
    hits5 = List_append(hits_geneplus_5,hits_geneminus_5);
    hits3 = List_append(hits_geneplus_3,hits_geneminus_3);

    stage3pairarray =
      consolidate_paired_results(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),&(*final_pairtype),
				 &(*stage3array5),&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
				 &(*stage3array3),&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
				 hitpairs,samechr,conc_transloc,terminals,hits5,hits3,
				 query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				 queryseq5,queryuc_ptr_5,quality_string_5,querylength5,
				 queryseq3,queryuc_ptr_3,quality_string_3,querylength3,
				 cutoff_level_5,cutoff_level_3,min_coverage_5,min_coverage_3,
				 oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR);
    History_free(&gmap_history_3);
    History_free(&gmap_history_5);
    Compress_free(&query5_compress_fwd);
    Compress_free(&query5_compress_rev);
    Compress_free(&query3_compress_fwd);
    Compress_free(&query3_compress_rev);
    Stage1_free(&this_geneminus_5,querylength5);
    Stage1_free(&this_geneminus_3,querylength3);
    Stage1_free(&this_geneplus_5,querylength5);
    Stage1_free(&this_geneplus_3,querylength3);
    /* return stage3pairarray */
  }

#ifdef HAVE_ALLOCA
  if (querylength5 <= MAX_STACK_READLENGTH) {
    FREEA(queryrc5);
  } else {
    FREE(queryrc5);
  }
  if (querylength3 <= MAX_STACK_READLENGTH) {
    FREEA(queryrc3);
  } else {
    FREE(queryrc3);
  }
#else
  FREE(queryrc5);
  FREE(queryrc3);
#endif

  return stage3pairarray;
}


Stage3pair_T *
Stage1_paired_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
		    Stage3end_T **stage3array5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
		    Stage3end_T **stage3array3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
		    Shortread_T queryseq5, Shortread_T queryseq3,
		    Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev, int indexdb_size_threshold,
		    Floors_T *floors_array,
		    double user_maxlevel_float, double user_mincoverage_float, int indel_penalty_middle, int indel_penalty_end,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_shortend,
		    Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    Chrpos_T pairmax, bool keep_floors_p) {

  if (mode == STANDARD || mode == CMET_STRANDED || mode == ATOI_STRANDED || mode == TTOC_STRANDED) {
    return paired_read(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),&(*final_pairtype),
		       &(*stage3array5),&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
		       &(*stage3array3),&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
		       queryseq5,queryseq3,indexdb_fwd,indexdb_rev,indexdb_size_threshold,
		       floors_array,user_maxlevel_float,user_mincoverage_float,indel_penalty_middle,indel_penalty_end,
		       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
		       localsplicing_penalty,distantsplicing_penalty,min_shortend,
		       oligoindices_major,oligoindices_minor,
		       pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,pairmax,keep_floors_p);

  } else if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
    return paired_read_tolerant_nonstranded(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),&(*final_pairtype),
					    &(*stage3array5),&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
					    &(*stage3array3),&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
					    queryseq5,queryseq3,indexdb_fwd,indexdb_rev,indexdb_size_threshold,
					    floors_array,user_maxlevel_float,user_mincoverage_float,indel_penalty_middle,indel_penalty_end,
					    allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
					    localsplicing_penalty,distantsplicing_penalty,min_shortend,
					    oligoindices_major,oligoindices_minor,
					    pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,pairmax,keep_floors_p);
  } else {
    fprintf(stderr,"Do not recognize mode %d\n",mode);
    abort();
  }
}


void
Stage1hr_cleanup () {
  FREE(chroffsets);
  FREE(chrhighs);
  FREE(chrlengths);
  return;
}


void
Stage1hr_setup (bool use_sarray_p_in, bool use_only_sarray_p_in, int index1part_in, int index1interval_in,
		int spansize_in, int max_anchors_in, Univ_IIT_T chromosome_iit_in, int nchromosomes_in,
		Genome_T genome_in, Genome_T genomealt, Mode_T mode_in, int maxpaths_search_in,

		Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		Chrpos_T *splicedists_in, int nsplicesites_in,
		
		bool novelsplicingp_in, bool knownsplicingp_in, bool find_dna_chimeras_p_in,
		bool distances_observed_p_in, int subopt_levels_in, int min_indel_end_matches_in,
		int max_middle_insertions_in, int max_middle_deletions_in,
		Chrpos_T shortsplicedist_in, Chrpos_T shortsplicedist_known_in, Chrpos_T shortsplicedist_novelend_in,
		Chrpos_T min_intronlength_in,

		int min_distantsplicing_end_matches_in, int min_distantsplicing_identity_in,

		int nullgap_in, int maxpeelback_in, int maxpeelback_distalmedial_in,
		int extramaterial_end_in, int extramaterial_paired_in,
		int gmap_mode, int trigger_score_for_gmap_in, int gmap_allowance_in,
		int max_gmap_pairsearch_in, int max_gmap_segments_in,
		int max_gmap_improvement_in, int antistranded_penalty_in,
		int max_floors_readlength_in) {
  bool gmapp = false;

  use_sarray_p = use_sarray_p_in;
  use_only_sarray_p = use_only_sarray_p_in;

  index1part = index1part_in;
  index1interval = index1interval_in;
  two_index1intervals = index1interval_in + index1interval_in;
  spansize = spansize_in;
  max_anchors = max_anchors_in;

  min_kmer_readlength = index1part_in + index1interval_in - 1;
  chromosome_iit = chromosome_iit_in;
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  nchromosomes = nchromosomes_in;
  genome = genome_in;

  if (use_only_sarray_p == false) {
    Univ_IIT_intervals_setup(&chroffsets,&chrhighs,&chrlengths,chromosome_iit,nchromosomes,circular_typeint);
  }

#ifdef HAVE_64_BIT
  leftreadshift = 64 - index1part - index1part;
  oligobase_mask = ~(~0ULL << 2*index1part);
#else
  leftreadshift = 32 - index1part - index1part;
  oligobase_mask = ~(~0U << 2*index1part);
#endif
  one_miss_querylength = spansize + spansize - (index1interval - 1); /* For 12-mers, 22 */

#if 0
  /* Should be 2 index1parts to handle a mismatch, plus a shift of 2, minus 3 if second one aligns */
  /* But this leads to many false positives, and GMAP can handle these other cases. */
  end_miss_two = spansize + (index1interval - 1) + spansize - index1interval;
#else
  end_miss_one = spansize + (index1interval - 1);
  end_miss_two = spansize + (index1interval - 1) + spansize - index1interval;
#endif
  
  mode = mode_in;
  maxpaths_search = maxpaths_search_in;

  splicesites = splicesites_in;
  splicetypes = splicetypes_in;
  splicedists = splicedists_in;
  nsplicesites = nsplicesites_in;

  novelsplicingp = novelsplicingp_in;
  knownsplicingp = knownsplicingp_in;
  find_dna_chimeras_p = find_dna_chimeras_p_in;
  distances_observed_p = distances_observed_p_in;

  subopt_levels = subopt_levels_in;
  min_indel_end_matches = min_indel_end_matches_in;
  max_middle_insertions_default = max_middle_insertions_in;
  max_middle_deletions = max_middle_deletions_in;

  shortsplicedist = shortsplicedist_in;

  overall_max_distance = shortsplicedist;
  if (max_middle_deletions > (int) overall_max_distance) {
    overall_max_distance = max_middle_deletions;
  }
  if (max_middle_insertions_default > (int) overall_max_distance) {
    overall_max_distance = max_middle_insertions_default;
  }

  min_intronlength = min_intronlength_in;
  min_distantsplicing_end_matches = min_distantsplicing_end_matches_in;
  min_distantsplicing_identity = min_distantsplicing_identity_in;

  nullgap = nullgap_in;
  maxpeelback = maxpeelback_in;
  maxpeelback_distalmedial = maxpeelback_distalmedial_in;
  extramaterial_end = extramaterial_end_in;
  extramaterial_paired = extramaterial_paired_in;

  gmap_segments_p = false;
  gmap_pairsearch_p = false;
  gmap_indel_knownsplice_p = false;
  gmap_improvement_p = false;

  fprintf(stderr,"GMAP modes:");
  if ((gmap_mode & GMAP_PAIRSEARCH) != 0) {
    if (gmapp == true) {
      fprintf(stderr,",");
    } else {
      gmapp = true;
    }
    fprintf(stderr," pairsearch");
    gmap_pairsearch_p = true;
  }
  if ((gmap_mode & GMAP_INDEL_KNOWNSPLICE) != 0) {
    if (gmapp == true) {
      fprintf(stderr,",");
    } else {
      gmapp = true;
    }
    fprintf(stderr," indel_knownsplice");
    gmap_indel_knownsplice_p = true;
  }
  if ((gmap_mode & GMAP_TERMINAL) != 0) {
    if (gmapp == true) {
      fprintf(stderr,",");
    } else {
      gmapp = true;
    }
    fprintf(stderr," segments");
    gmap_segments_p = true;
  }
  if ((gmap_mode & GMAP_IMPROVEMENT) != 0) {
    if (gmapp == true) {
      fprintf(stderr,",");
    } else {
      gmapp = true;
    }
    fprintf(stderr," improvement");
    gmap_improvement_p = true;
  }
  if (gmapp == false) {
    fprintf(stderr," none");
  }
  fprintf(stderr,"\n");


  trigger_score_for_gmap = trigger_score_for_gmap_in;
  gmap_allowance = gmap_allowance_in;

  max_gmap_pairsearch = max_gmap_pairsearch_in;
  max_gmap_segments = max_gmap_segments_in;
  max_gmap_improvement = max_gmap_improvement_in;

  antistranded_penalty = antistranded_penalty_in;

  if (genomealt != NULL) {
    snpp = true;
  } else {
    snpp = false;
  }

  max_floors_readlength = max_floors_readlength_in;

  return;
}
