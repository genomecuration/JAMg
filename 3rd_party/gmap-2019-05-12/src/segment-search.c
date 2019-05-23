static char rcsid[] = "$Id: segment-search.c 218694 2019-03-19 17:40:26Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "segment-search.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset() */
#include <math.h>
#include <ctype.h>		/* for tolower() */
#include "assert.h"
#include "mem.h"
#include "oligo.h"
#include "comp.h"

#include "list.h"
#include "intlist.h"
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
#endif
#include "univdiag.h"
#include "univdiagdef.h"

#ifndef LARGE_GENOMES
#include "merge-diagonals-simd-uint4.h"
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
#include "merge-diagonals-simd-uint8.h"
#else
#include "merge-diagonals-heap.h"
#endif

#ifdef LARGE_GENOMES
#include "intersect-large.h"
#endif
#include "intersect.h"
#include "sedgesort.h"
#include "path-solve.h"


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

/* Records */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Record_overlap_p */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif


#define MAX_INDEX1INTERVAL 3
#define STAGE2_MIN_OLIGO 3	/* Actually 6, but we are adding index1interval to this */
#define MAX_ALLOCATION 200

#define ABSOLUTE_LIMIT 5000


static Mode_T mode;
static int leftreadshift;
static Oligospace_T oligobase_mask; /* same as kmer_mask */

static int index1part;
static int index1interval;
static int max_anchors;

static int nchromosomes;
static Univ_IIT_T chromosome_iit;
static int circular_typeint;

static Univcoord_T *chroffsets;
static Univcoord_T *chrhighs;
static Chrpos_T *chrlengths; /* May differ from chrhigh - chroffset in circular chromosomes */

/* For spliceable */
static int max_deletionlen;
static Chrpos_T overall_max_distance;
static Chrpos_T shortsplicedist;

/* Splicing */
static Univcoord_T *splicesites;
static Splicetype_T *splicetypes;
static Chrpos_T *splicedists;
static int nsplicesites;


#if 0
static int
Record_diagonal_cmp (const void *a, const void *b) {
  Record_T x = * (Record_T *) a;
  Record_T y = * (Record_T *) b;

  if (x->diagonal < y->diagonal) {
    return -1;
  } else if (y->diagonal < x->diagonal) {
    return +1;
  } else {
    return 0;
  }
}
#endif

static int
Record_querypos5_ascending_cmp (const void *a, const void *b) {
  Record_T x = * (Record_T *) a;
  Record_T y = * (Record_T *) b;

  if (x->querypos < y->querypos) {
    return -1;
  } else if (y->querypos < x->querypos) {
    return +1;
  } else {
    return 0;
  }
}

static int
Record_querypos3_ascending_cmp (const void *a, const void *b) {
  Record_T x = * (Record_T *) a;
  Record_T y = * (Record_T *) b;

  if (x->queryend < y->queryend) {
    return -1;
  } else if (y->queryend < x->queryend) {
    return +1;
  } else {
    return 0;
  }
}


static bool
Record_overlap_p (Record_T segment1, Record_T segment2) {
  int querypos5, querypos3;

  debug2(printf("Entered Record_overlap_p with %d..%d and %d..%d\n",
		segment1->querypos,segment1->queryend,
		segment2->querypos,segment2->queryend));

  if (segment1->queryend < segment2->querypos) {
    debug2(printf("No overlap => false\n"));
    return false;
  } else if (segment2->queryend < segment1->querypos) {
    debug2(printf("No overlap => false\n"));
    return false;
  } else {
    if (segment1->queryend < segment2->queryend) {
      querypos3 = segment1->queryend;
    } else {
      querypos3 = segment2->queryend;
    }

    if (segment1->querypos > segment2->querypos) {
      querypos5 = segment1->querypos;
    } else {
      querypos5 = segment2->querypos;
    }

    debug2(printf("querypos5 %d, querypos3 %d\n",querypos5,querypos3));

    if (3 * (querypos3 - querypos5) > (segment1->queryend - segment1->querypos)) {
      debug2(printf("Amount of overlap is significant compared to segment1\n"));
      return true;
    } else if (3 * (querypos3 - querypos5) > (segment2->queryend - segment2->querypos)) {
      debug2(printf("Amount of overlap is significant compared to segment2\n"));
      return true;
    } else {
      debug2(printf("Amount of overlap is insignificant\n"));
      return false;
    }
  }
}


#if 0
/* Previously used for LARGE_GENOMES */
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



#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))

#if 0
/* Tries to do too much.  Instead, just push diagonals onto left and right, and let Path_solve take care of it */
/* Segment chaining */
/* Modified from previous convert_plus_segments_to_gmap */
/* Generic, for either plus or minus genomic strand */
static List_T
solve_segments_old (int *found_score, List_T hits, List_T *complete_paths, int querylength,
		    Record_T *segments, int nanchors, int nsegments, int *order,
		    char *queryptr, Compress_T query_compress,
		    bool plusp, int genestrand, bool require_pairing_p,
		    bool paired_end_p, bool first_read_p, int level) {
  int nmisses_allowed = 5;

  Record_T anchor_segment, base_segment, trial_segment, segment;
  int anchori, anchork, startk, endk, n, i, j, firstj, lastj, k, l, best_starti, best_endi;

  Univdiag_T middle_diagonal;
  List_T complete_path, middle_path, right_diagonals, left_diagonals;

  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;
  Chrnum_T chrnum;
  bool novelp;		 /* Want any of the segments in startk..(endk-1) to not be used */
  bool pairablep;		/* Want any of the segments in startk..(endk-1) to be pairable */

  int querypos, max_querypos5, max_querypos3, min_querypos5, min_querypos3, boundpos;
  Univcoord_T genomepos, min_genomepos, max_genomepos;

  Record_T *sorted5, *sorted5_allocated, *sorted3, *sorted3_allocated;
  int *scores, *scores_allocated, best_score, score;
  int *prev_left, *prev_right, *prev5_allocated, *prev3_allocated, besti;

  Univcoord_T left;
  bool foundp;


  if (nsegments > 0) {
#ifdef HAVE_ALLOCA
    if (nsegments < MAX_ALLOCATION) {
      prev5_allocated = (int *) ALLOCA(nsegments*sizeof(int));
      prev3_allocated = (int *) ALLOCA(nsegments*sizeof(int));
      scores_allocated = (int *) ALLOCA(nsegments*sizeof(int));
      sorted5_allocated = (Record_T *) ALLOCA(nsegments*sizeof(Record_T));
      sorted3_allocated = (Record_T *) ALLOCA(nsegments*sizeof(Record_T));
    } else {
      prev5_allocated = (int *) MALLOC(nsegments*sizeof(int));
      prev3_allocated = (int *) MALLOC(nsegments*sizeof(int));
      scores_allocated = (int *) MALLOC(nsegments*sizeof(int));
      sorted5_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
      sorted3_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
    }
#else
    prev5_allocated = (int *) MALLOC(nsegments*sizeof(int));
    prev3_allocated = (int *) MALLOC(nsegments*sizeof(int));
    scores_allocated = (int *) MALLOC(nsegments*sizeof(int));
    sorted5_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
    sorted3_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
#endif
  }

  for (anchork = nsegments - 1; anchork >= nsegments - nanchors; anchork--) {
    anchori = order[anchork];
    anchor_segment = segments[anchori];

    startk = anchor_segment->starti;
    endk = anchor_segment->endi;
    debug13(printf("Found segments %d to %d inclusive for anchor #%d, %d..%d %u, chrpos %u, (range %u..%u)\n",
		   startk,endk,anchori,anchor_segment->querypos,anchor_segment->queryend,
		   anchor_segment->diagonal,anchor_segment->diagonal - chroffset,
		   anchor_segment->lowpos,anchor_segment->highpos));
#ifdef DEBUG13
    for (k = startk; k <= endk; k++) {
      printf("  #%d: %d..%d %u\n",k,segments[k]->querypos,segments[k]->queryend,segments[k]->diagonal);
    }
    printf("\n");
#endif

    novelp = pairablep = false;
    if (anchor_segment->usedp == false) {
      novelp = true;
    }
    if (anchor_segment->pairablep == true) {
      pairablep = true;
    }

    n = endk - startk + 1;
    debug13(printf("n = %d\n",n));
    prev_left = &(prev5_allocated[startk]);
    prev_right = &(prev3_allocated[startk]);
    scores = &(scores_allocated[startk]);
    sorted5 = &(sorted5_allocated[startk]);
    sorted3 = &(sorted3_allocated[startk]);


    /* Dynamic programming on left (low) side (querypos5) */
    for (k = startk, i = 0; k <= endk; k++) {
      sorted5[i++] = segments[k];
    }
    qsort(sorted5,n,sizeof(Record_T),Record_querypos5_ascending_cmp);
#ifdef DEBUG13
    printf("Sorted by querypos5\n");
    for (i = 0; i < n; i++) {
      printf("  #%d: %d..%d %u\n",i,sorted5[i]->querypos,sorted5[i]->queryend,sorted5[i]->diagonal);
    }
#endif

    lastj = 0;
    while (lastj < n && sorted5[lastj]->querypos < anchor_segment->querypos) {
      lastj++;
    }
#ifdef DEBUG13
    printf("On the querypos5 side, considering:\n");
    for (i = 0; i < lastj; i++) {
      printf("  #%d: %d..%d left %u, range %u..%u\n",
	     i,sorted5[i]->querypos,sorted5[i]->queryend,sorted5[i]->diagonal,sorted5[i]->lowpos,sorted5[i]->highpos);
    }
#endif
    /* TODO: Re-sort by coordinates within querypos clusters */

    /* Go leftward to favor shorter splices */
    for (j = lastj - 1; j >= 0; --j) {
      base_segment = sorted5[j];
      debug13(printf("Base on the querypos5 side: %d..%d %u (range %u..%u)\n",
		     base_segment->querypos,base_segment->queryend,base_segment->diagonal,base_segment->lowpos,base_segment->highpos));
      
      /* Appending sorted5[i] to left of base_segment */
      best_score = 0;
      besti = -1;
      for (i = j - 1; i >= 0; --i) {
	trial_segment = sorted5[i];
	debug13(printf("BASE SEGMENT %d..%d left %u (range %u..%u)\n",
		       base_segment->querypos,base_segment->queryend,base_segment->diagonal,base_segment->lowpos,base_segment->highpos));
	debug13(printf("TRIAL SEGMENT ON 5' %d..%d left %u (range %u..%u)\n",
		       trial_segment->querypos,trial_segment->queryend,trial_segment->diagonal,trial_segment->lowpos,trial_segment->highpos));

	if (trial_segment->lowpos >= base_segment->lowpos + max_deletionlen) {
	  debug13(printf("(1) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE LEFT OF BASE\n"));
	  
	} else if (trial_segment->lowpos >= base_segment->lowpos && trial_segment->highpos <= base_segment->highpos) {
	  debug13(printf("(1) SKIPPING, SINCE BASE SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));
	  
	} else if (base_segment->lowpos >= trial_segment->lowpos && base_segment->highpos <= trial_segment->highpos) {
	  debug13(printf("(1) SKIPPING, SINCE TRIAL SUBSUMES BASE (PROBABLE GENOMIC REPEAT)\n"));

	} else if (trial_segment->highpos < base_segment->lowpos) {
	  debug13(printf("Query skip, so score is the segment itself %d\n",trial_segment->highpos - trial_segment->lowpos));
	  if ((score = (trial_segment->highpos - trial_segment->lowpos)) > best_score) {
	    best_score = score;
	    besti = i;
	  }
	} else if ((score = (base_segment->lowpos - trial_segment->lowpos)) > best_score) {
	  debug13(printf("Overlap, so score is the new part %d\n",base_segment->lowpos - trial_segment->lowpos));
	  best_score = score;
	  besti = i;
	} else {
	  debug13(printf("Overlap, so score is the new part %d\n",base_segment->lowpos - trial_segment->lowpos));
	}
      }
      scores[j] = base_segment->highpos - base_segment->lowpos;
      debug13(printf("Best prev is %d with score %d\n",besti,best_score));
      if ((prev_left[j] = besti) >= 0) {
	scores[j] += best_score;
      }
    }

    /* Go leftward to favor shorter splices */
    debug13(printf("Now trying to append to anchor segment\n"));
    best_score = 0;
    best_starti = -1;
    for (j = lastj - 1; j >= 0; --j) {
      trial_segment = sorted5[j];
      debug13(printf("ANCHOR SEGMENT %d..%d left %u (range %u..%u)\n",
		     anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal,anchor_segment->lowpos,anchor_segment->highpos));
      debug13(printf("TRIAL SEGMENT ON 5' %d..%d left %u (range %u..%u)\n",
		     trial_segment->querypos,trial_segment->queryend,trial_segment->diagonal,trial_segment->lowpos,trial_segment->highpos));

      if (trial_segment->lowpos >= anchor_segment->lowpos + max_deletionlen) {
	debug13(printf("(2) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE LEFT OF ANCHOR\n"));
	
	} else if (trial_segment->lowpos >= anchor_segment->lowpos && trial_segment->highpos <= anchor_segment->highpos) {
	  debug13(printf("(2) SKIPPING, SINCE ANCHOR SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));
	  
	} else if (anchor_segment->lowpos >= trial_segment->lowpos && anchor_segment->highpos <= trial_segment->highpos) {
	  debug13(printf("(2) SKIPPING, SINCE TRIAL SUBSUMES ANCHOR (PROBABLE GENOMIC REPEAT)\n"));


      } else if (trial_segment->highpos < anchor_segment->lowpos) {
	debug13(printf("Query skip, so score is from the path %d\n",scores[j]));
	if ((score = scores[j]) > best_score) {
	  best_score = score;
	  best_starti = j;
	}
      } else if ((score = scores[j] - (trial_segment->highpos - anchor_segment->lowpos)) > best_score) {
	debug13(printf("Overlap, so score is from the path %d minus (%d - %d)\n",
		       scores[j],trial_segment->highpos,anchor_segment->lowpos));
	best_score = score;
	best_starti = j;
      } else {
	debug13(printf("Score %d does not exceed the best score %d\n",score,best_score));
      }
    }
    debug13(printf("Found best_starti to be %d\n",best_starti));

    /* Evaluate set of segments */
    for (k = best_starti; k >= 0; k = prev_left[k]) {
      if (sorted5[k]->usedp == false) {
	novelp = true;
      }
      if (sorted5[k]->pairablep == true) {
	pairablep = true;
      }
    }


    /* Dynamic programming on right (high) side (querypos3) */
    for (k = startk, i = 0; k <= endk; k++) {
      sorted3[i++] = segments[k];
    }
    qsort(sorted3,n,sizeof(Record_T),Record_querypos3_ascending_cmp);
#ifdef DEBUG13
    printf("Sorted by querypos3\n");
    for (i = 0; i < n; i++) {
      printf("  #%d: %d..%d %u\n",i,sorted3[i]->querypos,sorted3[i]->queryend,sorted3[i]->diagonal);
    }
#endif

    firstj = n - 1;
    while (firstj >= 0 && sorted3[firstj]->queryend > anchor_segment->queryend) {
      firstj--;
    }
#ifdef DEBUG13
    printf("On the querypos3 side, considering:\n");
    for (j = firstj + 1; j < n; j++) {
      printf("  #%d: %d..%d left %u, range %u..%u\n",
	     i,sorted3[j]->querypos,sorted3[j]->queryend,sorted3[j]->diagonal,sorted3[j]->lowpos,sorted3[j]->highpos);
    }
#endif
    /* TODO: Re-sort by coordinates within querypos clusters */


    /* Go rightward to favor shorter splices */
    for (j = firstj + 1; j < n; j++) {
      base_segment = sorted3[j];
      debug13(printf("Base on the querypos3 side: %d..%d %u (range %u..%u)\n",
		     base_segment->querypos,base_segment->queryend,base_segment->diagonal,base_segment->lowpos,base_segment->highpos));

      /* Appending sorted3[i] to right of base_segment */
      best_score = 0;
      besti = -1;
      for (i = j + 1; i < n; i++) {
	trial_segment = sorted3[i];
	debug13(printf("BASE SEGMENT %d..%d left %u (range %u..%u)\n",
		       base_segment->querypos,base_segment->queryend,base_segment->diagonal,base_segment->lowpos,base_segment->highpos));
	debug13(printf("TRIAL SEGMENT ON 3' %d..%d left %u (range %u..%u)\n",
		       trial_segment->querypos,trial_segment->queryend,trial_segment->diagonal,trial_segment->lowpos,trial_segment->highpos));

	if (trial_segment->highpos + max_deletionlen <= base_segment->highpos) {
	  debug13(printf("(3) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE RIGHT OF BASE\n"));
	  
	} else if (trial_segment->lowpos >= base_segment->lowpos && trial_segment->highpos <= base_segment->highpos) {
	  debug13(printf("(3) SKIPPING, SINCE BASE SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));
	  
	} else if (base_segment->lowpos >= trial_segment->lowpos && base_segment->highpos <= trial_segment->highpos) {
	  debug13(printf("(3) SKIPPING, SINCE TRIAL SUBSUMES BASE (PROBABLE GENOMIC REPEAT)\n"));

	} else if (trial_segment->lowpos > base_segment->highpos) {
	  debug13(printf("Query skip, so score is the segment itself %d\n",trial_segment->highpos - trial_segment->lowpos));
	  if ((score = (trial_segment->highpos - trial_segment->lowpos)) > best_score) {
	    best_score = score;
	    besti = i;
	  }
	} else if ((score = (trial_segment->highpos - base_segment->highpos)) > best_score) {
	  debug13(printf("Overlap, so score is the new part %d\n",trial_segment->highpos - base_segment->highpos));
	  best_score = score;
	  besti = i;
	} else {
	  debug13(printf("Overlap, so score is the new part %d\n",trial_segment->highpos - base_segment->highpos));
	}
      }
      scores[j] = base_segment->highpos - base_segment->lowpos;
      debug13(printf("Best prev is %d with score %d\n",besti,best_score));
      if ((prev_right[j] = besti) >= 0) {
	scores[j] += best_score;
      }
    }

#ifdef DEBUG13
    for (j = firstj + 1; j < n; j++) {
      printf("scores[%d] = %d.  prev_right is %d\n",j,scores[j],prev_right[j]);
    }
#endif

    /* Go rightward to favor shorter splices */
    debug13(printf("Now trying to append to anchor segment\n"));
    best_score = 0;
    best_endi = -1;
    for (j = firstj + 1; j < n; j++) {
      trial_segment = sorted3[j];
      debug13(printf("ANCHOR SEGMENT %d..%d left %u (range %u..%u)\n",
		     anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal,anchor_segment->lowpos,anchor_segment->highpos));
      debug13(printf("TRIAL SEGMENT ON 5' %d..%d left %u (range %u..%u)\n",
		     trial_segment->querypos,trial_segment->queryend,trial_segment->diagonal,trial_segment->lowpos,trial_segment->highpos));

      if (trial_segment->highpos + max_deletionlen <= anchor_segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE RIGHT OF ANCHOR\n"));
	
      } else if (trial_segment->lowpos >= anchor_segment->lowpos && trial_segment->highpos <= anchor_segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE ANCHOR SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));

      } else if (anchor_segment->lowpos >= trial_segment->lowpos && anchor_segment->highpos <= trial_segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE TRIAL SUBSUMES ANCHOR (PROBABLE GENOMIC REPEAT)\n"));

      } else if (trial_segment->lowpos > anchor_segment->highpos) {
	debug13(printf("Query skip, so score is from the path %d\n",scores[j]));
	if ((score = scores[j]) > best_score) {
	  best_score = score;
	  best_endi = j;
	}
	
      } else if ((score = scores[j] - (anchor_segment->highpos - trial_segment->lowpos)) > best_score) {
	debug13(printf("Overlap, so score is from the path %d minus (%d - %d)\n",
		       scores[j],anchor_segment->highpos,trial_segment->lowpos));
	best_score = score;
	best_endi = j;
      } else {
	debug13(printf("Score %d does not exceed the best score %d\n",score,best_score));
      }
    }
    debug13(printf("Found best_endi to be %d\n",best_endi));

    /* Evaluate set of segments */
    for (k = best_endi; k >= 0; k = prev_right[k]) {
      if (sorted3[k]->usedp == false) {
	novelp = true;
      }
      if (sorted3[k]->pairablep == true) {
	pairablep = true;
      }
    }

    debug13(printf("Processing segments %d to %d inclusive: novelp %d, pairablep %d\n",
		   startk,endk,novelp,pairablep));
    if (novelp == true && (pairablep == true || require_pairing_p == false)) {

      left = anchor_segment->diagonal /*- querylength*/; /* NEW FORMULA: Corresponds to querypos 0 */

      /* left is safer than anchor_segment->lowpos, in case trim goes further left than lowpos */
      chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

      middle_path = (List_T) NULL;
      right_diagonals = left_diagonals = (List_T) NULL;

      /* F.  Compute middle diagonal first */
      debug13(printf("Anchor diagonal %llu, querypos %d..%d, usedp %d, pairablep %d\n",
		     (unsigned long long) anchor_segment->diagonal,
		     anchor_segment->querypos,anchor_segment->queryend,anchor_segment->usedp,anchor_segment->pairablep));
      middle_diagonal = Univdiag_new(univdiagpool,anchor_segment->querypos,
				     anchor_segment->queryend,anchor_segment->diagonal);
      min_genomepos = left + anchor_segment->querypos;
      max_genomepos = left + anchor_segment->queryend;

      /* F.  Process right diagonals */
      boundpos = anchor_segment->queryend;
      for (k = best_endi; k >= 0; k = prev_right[k]) {
	segment = sorted3[k];
	segment->allowablep = false;
	segment->ambiguousp = false;
	debug13(printf("Right diagonal %llu, querypos %d..%d, usedp %d, pairablep %d\n",
		       (unsigned long long) segment->diagonal,
		       segment->querypos,segment->queryend,segment->usedp,segment->pairablep));
#if 0
	if (anchor_segment->diagonal < chroffset + chrlength && segment->diagonal > chroffset + chrlength) {
	  debug13(printf("Cannot cross circular origin\n"));
	} else {
#endif
	  querypos = segment->queryend;
	  if (querypos > boundpos) {
	    left = segment->diagonal /*- querylength*/; /* NEW FORMULA */
	    genomepos = left + querypos;
	    if (genomepos <= max_genomepos) {
	      debug13(printf("Not allowing right diagonal that doesn't advance max_genomepos\n"));
	    } else if ((anchor_segment->queryend - segment->querypos) > (segment->queryend - segment->querypos) / 2) {
	      debug13(printf("Not allowing right diagonal that mainly overlaps anchor diagonal\n"));
	    } else if (Record_overlap_p(segment,anchor_segment) == true) {
	      debug13(printf("Not allowing right diagonal that overlaps anchor diagonal\n"));
	    } else {
	      debug13(printf("Marking right diagonal as allowable\n"));
	      segment->allowablep = true;
	      boundpos = segment->queryend;
	      max_genomepos = left + boundpos;
	    }
	  }
#if 0
	}
#endif
      }
	      
      /* Put overlapping segments into right_diagonals */
      min_querypos5 = min_querypos3 = querylength;
      for (k = best_endi; k >= 0; k = prev_right[k]) {
	if (sorted3[k]->allowablep == true) {
	  for (l = prev_right[k]; l >= 0; l = prev_right[l]) {
	    if (sorted3[l]->allowablep == true) {
	      if (Record_overlap_p(sorted3[k],sorted3[l]) == true) {
		if (sorted3[k]->ambiguousp == false) {
		  right_diagonals = Univdiagpool_push(right_diagonals,univdiagpool,
						      sorted3[k]->querypos,sorted3[k]->queryend,sorted3[k]->diagonal);
		  if (sorted3[k]->querypos < min_querypos5) {
		    min_querypos5 = sorted3[k]->querypos;
		  }
		  if (sorted3[k]->queryend < min_querypos3) {
		    min_querypos3 = sorted3[k]->queryend;
		  }
		  sorted3[k]->ambiguousp = true;
		}
		if (sorted3[l]->ambiguousp == false) {
		  right_diagonals = Univdiagpool_push(right_diagonals,univdiagpool,
						      sorted3[l]->querypos,sorted3[l]->queryend,sorted3[l]->diagonal);
		  if (sorted3[l]->querypos < min_querypos5) {
		    min_querypos5 = sorted3[l]->querypos;
		  }
		  if (sorted3[l]->queryend < min_querypos3) {
		    min_querypos3 = sorted3[l]->queryend;
		  }
		  sorted3[l]->ambiguousp = true;
		}
	      }
	    }
	  }
	}
      }

      /* Put non-overlapping segments into middle_path (if they are to the left of all right diagonals) */
      for (k = best_endi; k >= 0; k = prev_right[k]) {
	segment = sorted3[k];
	if (segment->allowablep == true && segment->ambiguousp == false) {
	  if (segment->querypos < min_querypos5 && segment->queryend < min_querypos3) {
	    debug13(printf("Putting right diagonal onto the middle path\n"));
	    middle_path = Univdiagpool_push(middle_path,univdiagpool,
					    segment->querypos,segment->queryend,segment->diagonal);
	  } else {
	    debug13(printf("Putting right diagonal into right diagonals\n"));
	    right_diagonals = Univdiagpool_push(right_diagonals,univdiagpool,
						segment->querypos,segment->queryend,segment->diagonal);
	  }
	}
      }

      middle_path = List_reverse(middle_path);

      /* F.  Push on middle diagonal */
      middle_path = List_push(middle_path,(void *) middle_diagonal);

      /* F.  Process left diagonals */
      boundpos = anchor_segment->querypos;
      for (k = best_starti; k >= 0; k = prev_left[k]) {
	segment = sorted5[k];
	segment->allowablep = false;
	segment->ambiguousp = false;
	debug13(printf("Left diagonal %llu, querypos %d..%d, usedp %d, pairablep %d\n",
		       (unsigned long long) segment->diagonal,
		       segment->querypos,segment->queryend,segment->usedp,segment->pairablep));
#if 0
	if (anchor_segment->diagonal > chroffset + chrlength && segment->diagonal < chroffset + chrlength) {
	  debug13(printf("Cannot cross circular origin\n"));
	} else {
#endif
	  querypos = segment->querypos;
	  if (querypos < boundpos) {
	    left = segment->diagonal /*- querylength*/; /* NEW FORMULA */
	    genomepos = left + querypos;
	    if (genomepos >= min_genomepos) {
	      debug13(printf("Not allowing left diagonal that doesn't advance min_genomepos\n"));
	    } else if ((segment->queryend - anchor_segment->querypos) > (segment->queryend - segment->querypos) / 2) {
	      debug13(printf("Not allowing left diagonal that mainly overlaps middle diagonal\n"));
	    } else if (Record_overlap_p(segment,anchor_segment) == true) {
	      debug13(printf("Not allowing left diagonal that overlaps anchor diagonal\n"));
	    } else {
	      debug13(printf("Marking left diagonal as allowable\n"));
	      segment->allowablep = true;
	      boundpos = segment->querypos;
	      min_genomepos = left + boundpos;
	    }
	  }
#if 0
	}
#endif
      }

      /* Put overlapping segments into left_diagonals */
      max_querypos5 = max_querypos3 = 0;
      for (k = best_starti; k >= 0; k = prev_left[k]) {
	if (sorted5[k]->allowablep == true) {
	  for (l = prev_left[k]; l >= 0; l = prev_left[l]) {
	    if (sorted5[l]->allowablep == true) {
	      if (Record_overlap_p(sorted5[k],sorted5[l]) == true) {
		if (sorted5[k]->ambiguousp == false) {
		  left_diagonals = Univdiagpool_push(left_diagonals,univdiagpool,
						     sorted5[k]->querypos,sorted5[k]->queryend,sorted5[k]->diagonal);
		  if (sorted5[k]->querypos > max_querypos5) {
		    max_querypos5 = sorted5[k]->querypos;
		  }
		  if (sorted5[k]->queryend > max_querypos3) {
		    max_querypos3 = sorted5[k]->queryend;
		  }
		  sorted5[k]->ambiguousp = true;
		}
		if (sorted5[l]->ambiguousp == false) {
		  left_diagonals = Univdiagpool_push(left_diagonals,univdiagpool,
						     sorted5[l]->querypos,sorted5[l]->queryend,sorted5[l]->diagonal);
		  if (sorted5[l]->querypos > max_querypos5) {
		    max_querypos5 = sorted5[l]->querypos;
		  }
		  if (sorted5[l]->queryend > max_querypos3) {
		    max_querypos3 = sorted5[l]->queryend;
		  }
		  sorted5[l]->ambiguousp = true;
		}
	      }
	    }
	  }
	}
      }

      /* Put non-overlapping segments into middle_path (if they are to the right of all left diagonals) */
      for (k = best_starti; k >= 0; k = prev_left[k]) {
	segment = sorted5[k];
	if (segment->allowablep == true && segment->ambiguousp == false) {
	  if (segment->querypos > max_querypos5 && segment->queryend > max_querypos3) {
	    debug13(printf("Putting left diagonal onto the middle path\n"));
	    middle_path = Univdiagpool_push(middle_path,univdiagpool,
					    segment->querypos,segment->queryend,segment->diagonal);
	  } else {
	    debug13(printf("Putting left diagonal into left diagonals\n"));
	    left_diagonals = Univdiagpool_push(left_diagonals,univdiagpool,
					       segment->querypos,segment->queryend,segment->diagonal);
	  }
	}
      }

      /* No need to reverse middle_path */

      /* TODO: Use localdb to add left and right diagonals.  Then modify
	 Path_solve_from_diagonals to handle a path, and call that */
      debug13(printf("SOLVING A PATH WITH %d LEFT, ONE MIDDLE, and %d RIGHT DIAGONALS\n\n",
		     List_length(left_diagonals),List_length(right_diagonals)));

      left = anchor_segment->diagonal;

      hits = Path_solve_from_diagonals(&foundp,&(*found_score_overall),&(*found_score_within_trims),hits,
				       middle_path->univdiagonal,middle_path->qstart,middle_path->qend,
				       right_diagonals,left_diagonals,queryptr,querylength,
				       mismatch_positions_alloc,spliceinfo,stream_alloc,streamsize_alloc,
				       query_compress,chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
				       /*nmismatches_allowed*/nmisses_allowed,paired_end_p,first_read_p,
				       intlistpool,univcoordlistpool,listpool,univdiagpool,
				       hitlistpool,/*method*/SEGMENT,level);
      /* Univdiag_gc(&right_diagonals); -- allocated by Univdiagpool_push */
      /* Univdiag_gc(&left_diagonals); -- allocated by Univdiagpool_push */
      /* Univdiag_gc(&middle_path); -- allocated by Univdiagpool_push */

      if (foundp == false) {
	/* *complete_paths = List_push(*complete_paths,(void *) complete_path); */

      } else {
	/* Mark segments from middle_path as used */
	anchor_segment->usedp = true;

	for (k = best_endi; k >= 0; k = prev_right[k]) {
	  segment = sorted3[k];
	  if (segment->allowablep == true && segment->ambiguousp == false) {
	    segment->usedp = true;
	  }
	}

	for (k = best_starti; k >= 0; k = prev_left[k]) {
	  segment = sorted5[k];
	  if (segment->allowablep == true && segment->ambiguousp == false) {
	    segment->usedp = true;
	  }
	}
      }
    }
  }

  if (nsegments > 0) {
#ifdef HAVE_ALLOCA
    if (nsegments < MAX_ALLOCATION) {
      FREEA(sorted3_allocated);
      FREEA(sorted5_allocated);
      FREEA(scores_allocated);
      FREEA(prev3_allocated);
      FREEA(prev5_allocated);
    } else {
      FREE(sorted3_allocated);
      FREE(sorted5_allocated);
      FREE(scores_allocated);
      FREE(prev3_allocated);
      FREE(prev5_allocated);
    }
#else
    FREE(sorted3_allocated);
    FREE(sorted5_allocated);
    FREE(scores_allocated);
    FREE(prev3_allocated);
    FREE(prev5_allocated);
#endif
  }

  return hits;
}
#endif


#if 0
/* Segment chaining */
/* Modified from previous convert_plus_segments_to_gmap */
/* Generic, for either plus or minus genomic strand */
static List_T
solve_segments (int *found_score, List_T hits, List_T *complete_paths, int querylength,
		Record_T *segments, int nanchors, int nsegments, int *order,
		char *queryptr, Compress_T query_compress,
		bool plusp, int genestrand, bool paired_end_p, bool first_read_p, int level) {
  int nmisses_allowed = 5;

  Record_T anchor_segment, segment;
  int anchori, anchork, startk, endk, n, i, firstj, lastj, j, k;

  Univdiag_T middle_diagonal;
  List_T complete_path, middle_path, right_diagonals, left_diagonals;

  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;
  Chrnum_T chrnum;

  Record_T *sorted5, *sorted5_allocated, *sorted3, *sorted3_allocated;

  Univcoord_T left;
  bool foundp;


  if (nsegments > 0) {
#ifdef HAVE_ALLOCA
    if (nsegments < MAX_ALLOCATION) {
      sorted5_allocated = (Record_T *) ALLOCA(nsegments*sizeof(Record_T));
      sorted3_allocated = (Record_T *) ALLOCA(nsegments*sizeof(Record_T));
    } else {
      sorted5_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
      sorted3_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
    }
#else
    sorted5_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
    sorted3_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
#endif
  }

  for (anchork = nsegments - 1; anchork >= nsegments - nanchors; anchork--) {
    anchori = order[anchork];
    anchor_segment = segments[anchori];
    left = anchor_segment->diagonal /*- querylength*/; /* NEW FORMULA: Corresponds to querypos 0 */

    /* left is safer than anchor_segment->lowpos, in case trim goes further left than lowpos */
    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
    
    startk = anchor_segment->starti;
    endk = anchor_segment->endi;
    debug13(printf("Found segments %d to %d inclusive for anchor #%d, %d..%d %u, chrpos %u (range %u..%u)\n",
		   startk,endk,anchori,anchor_segment->querypos,anchor_segment->queryend,
		   anchor_segment->diagonal,anchor_segment->diagonal - chroffset,
		   anchor_segment->lowpos,anchor_segment->highpos));
#ifdef DEBUG13
    for (k = startk; k <= endk; k++) {
      printf("  #%d: %d..%d %u\n",k,segments[k]->querypos,segments[k]->queryend,segments[k]->diagonal);
    }
    printf("\n");
#endif

    n = endk - startk + 1;
    debug13(printf("n = %d\n",n));
    sorted5 = &(sorted5_allocated[startk]);
    sorted3 = &(sorted3_allocated[startk]);


    /* middle_path = (List_T) NULL; */
    right_diagonals = left_diagonals = (List_T) NULL;

    /* Add diagonals to left (low) side (querypos5) */
    for (k = startk, i = 0; k <= endk; k++) {
      sorted5[i++] = segments[k];
    }
    qsort(sorted5,n,sizeof(Record_T),Record_querypos5_ascending_cmp);
#ifdef DEBUG13
    printf("Sorted by querypos5\n");
    for (i = 0; i < n; i++) {
      printf("  #%d: %d..%d %u\n",i,sorted5[i]->querypos,sorted5[i]->queryend,sorted5[i]->diagonal);
    }
#endif

    /* Skip obvious ones */
    lastj = 0;
    while (lastj < n && sorted5[lastj]->querypos < anchor_segment->querypos) {
      lastj++;
    }
#ifdef DEBUG13
    printf("On the querypos5 side, considering:\n");
    for (i = 0; i < lastj; i++) {
      printf("  #%d: %d..%d left %u, range %u..%u\n",
	     i,sorted5[i]->querypos,sorted5[i]->queryend,sorted5[i]->diagonal,sorted5[i]->lowpos,sorted5[i]->highpos);
    }
#endif

    /* Go leftward to favor shorter splices */
    debug13(printf("Trying to append to start of anchor segment\n"));
    for (j = lastj - 1; j >= 0; --j) {
      segment = sorted5[j];
      debug13(printf("ANCHOR SEGMENT %d..%d left %u (range %u..%u)\n",
		     anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal,anchor_segment->lowpos,anchor_segment->highpos));
      debug13(printf("TRIAL SEGMENT ON 5' %d..%d left %u (range %u..%u)\n",
		     segment->querypos,segment->queryend,segment->diagonal,segment->lowpos,segment->highpos));

      if (segment->lowpos >= anchor_segment->lowpos + max_deletionlen) {
	debug13(printf("(2) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE LEFT OF ANCHOR\n"));
	
      } else if (segment->lowpos >= anchor_segment->lowpos && segment->highpos <= anchor_segment->highpos) {
	debug13(printf("(2) SKIPPING, SINCE ANCHOR SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));
	  
      } else if (anchor_segment->lowpos >= segment->lowpos && anchor_segment->highpos <= segment->highpos) {
	debug13(printf("(2) SKIPPING, SINCE TRIAL SUBSUMES ANCHOR (PROBABLE GENOMIC REPEAT)\n"));

#if 0
      } else if (segment->highpos < chroffset) {
	/* Now handled by path-solve, which accounts for the best chromosome for the anchor segment */
	debug13(printf("Not allowing left diagonal that goes into next chromosome\n"));
#endif

      } else {
	debug13(printf("Candidate for start diagonal\n"));
	if (segment->diagonal + segment->querypos >= left + anchor_segment->querypos) {
	  debug13(printf("Not allowing left diagonal that doesn't advance min_genomepos\n"));
	} else if ((segment->queryend - anchor_segment->querypos) > (segment->queryend - segment->querypos) / 2) {
	  debug13(printf("Not allowing left diagonal that mainly overlaps middle diagonal\n"));
	} else if (Record_overlap_p(segment,anchor_segment) == true) {
	  debug13(printf("Not allowing left diagonal that overlaps anchor diagonal\n"));
	} else {
	  debug13(printf("Left diagonal is allowable\n"));
	  left_diagonals = List_push(left_diagonals,
				     (void *) Univdiag_new(segment->querypos,segment->queryend,segment->diagonal));
	}
      }
    }


    /* Add diagonals to right (high) side (querypos3) */
    for (k = startk, i = 0; k <= endk; k++) {
      sorted3[i++] = segments[k];
    }
    qsort(sorted3,n,sizeof(Record_T),Record_querypos3_ascending_cmp);
#ifdef DEBUG13
    printf("Sorted by querypos3\n");
    for (i = 0; i < n; i++) {
      printf("  #%d: %d..%d %u\n",i,sorted3[i]->querypos,sorted3[i]->queryend,sorted3[i]->diagonal);
    }
#endif

    /* Skip obvious ones */
    firstj = n - 1;
    while (firstj >= 0 && sorted3[firstj]->queryend > anchor_segment->queryend) {
      firstj--;
    }
#ifdef DEBUG13
    printf("On the querypos3 side, considering:\n");
    for (j = firstj + 1; j < n; j++) {
      printf("  #%d: %d..%d left %u, range %u..%u\n",
	     i,sorted3[j]->querypos,sorted3[j]->queryend,sorted3[j]->diagonal,sorted3[j]->lowpos,sorted3[j]->highpos);
    }
#endif

    /* Go rightward to favor shorter splices */
    debug13(printf("Trying to append to end of anchor segment\n"));
    for (j = firstj + 1; j < n; j++) {
      segment = sorted3[j];
      debug13(printf("ANCHOR SEGMENT %d..%d left %u (range %u..%u)\n",
		     anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal,anchor_segment->lowpos,anchor_segment->highpos));
      debug13(printf("TRIAL SEGMENT ON 5' %d..%d left %u (range %u..%u)\n",
		     segment->querypos,segment->queryend,segment->diagonal,segment->lowpos,segment->highpos));

      if (segment->highpos + max_deletionlen <= anchor_segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE RIGHT OF ANCHOR\n"));
	
      } else if (segment->lowpos >= anchor_segment->lowpos && segment->highpos <= anchor_segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE ANCHOR SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));

      } else if (anchor_segment->lowpos >= segment->lowpos && anchor_segment->highpos <= segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE TRIAL SUBSUMES ANCHOR (PROBABLE GENOMIC REPEAT)\n"));

#if 0
      } else if (segment->lowpos >= chrhigh) {
	/* Now handled by path-solve, which accounts for the best chromosome for the anchor segment */
	debug13(printf("Not allowing right diagonal that goes into next chromosome\n"));
#endif

      } else {
	debug13(printf("Candidate for end diagonal\n"));
	if (segment->diagonal + segment->queryend <= left + anchor_segment->queryend) {
	  debug13(printf("Not allowing right diagonal that doesn't advance max_genomepos\n"));
	} else if ((anchor_segment->queryend - segment->querypos) > (segment->queryend - segment->querypos) / 2) {
	  debug13(printf("Not allowing right diagonal that mainly overlaps anchor diagonal\n"));
	} else if (Record_overlap_p(segment,anchor_segment) == true) {
	  debug13(printf("Not allowing right diagonal that overlaps anchor diagonal\n"));
	} else {
	  debug13(printf("Right diagonal is allowable\n"));
	  right_diagonals = List_push(right_diagonals,
				      (void *) Univdiag_new(segment->querypos,segment->queryend,segment->diagonal));
	}
      }
    }

    /* Compute middle diagonal */
    debug13(printf("Anchor diagonal %llu, querypos %d..%d\n",
		   (unsigned long long) anchor_segment->diagonal,anchor_segment->querypos,anchor_segment->queryend));
    middle_diagonal = Univdiag_new(anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal);
    
    /* TODO: Use localdb to add left and right diagonals.  Then modify
       Path_solve_from_diagonals to handle a path, and call that */
    debug13(printf("SOLVING A PATH WITH %d LEFT, ONE MIDDLE, and %d RIGHT DIAGONALS\n\n",
		   List_length(left_diagonals),List_length(right_diagonals)));
    
    hits = Path_solve_from_diagonals(&foundp,&(*found_score_overall),&(*found_score_within_trims),hits,
				     middle_diagonal->univdiagonal,middle_diagonal->qstart,middle_diagonal->qend,
				     right_diagonals,left_diagonals,queryptr,querylength,
				     mismatch_positions_alloc,spliceinfo,stream_alloc,streamsize_alloc,
				     query_compress,chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
				     /*nmismatches_allowed*/nmisses_allowed,paired_end_p,first_read_p,
				     intlistpool,univcoordlistpool,listpool,univdiagpool,
				     hitlistpool,/*method*/SEGMENT,level);
    Univdiag_free(&middle_diagonal);
    Univdiag_gc(&right_diagonals);
    Univdiag_gc(&left_diagonals);
    
    if (foundp == false) {
      /* *complete_paths = List_push(*complete_paths,(void *) complete_path); */
    }
  }

  if (nsegments > 0) {
#ifdef HAVE_ALLOCA
    if (nsegments < MAX_ALLOCATION) {
      FREEA(sorted3_allocated);
      FREEA(sorted5_allocated);
    } else {
      FREE(sorted3_allocated);
      FREE(sorted5_allocated);
    }
#else
    FREE(sorted3_allocated);
    FREE(sorted5_allocated);
#endif
  }

  return hits;
}
#endif


#define SINGLETON 0
#define DUPLICATE -2
#define ADJACENT -1		/* Within overall_max_distance */


/* Returns status as count (if positive) or one of the flags SINGLETON, DUPLICATE, or ADJACENT */
static int *
compute_status_greedy (int *nanchors, int *nadjacent, Univcoord_T *diagonals, int ndiagonals) {
  int *status;
  int i, j, k, newi;
  int count, max_count = 0;

  status = (int *) CALLOC(ndiagonals,sizeof(int));
  *nanchors = *nadjacent = 0;

  /* Compute max_count */
  i = 0;
  while (i < ndiagonals) {
    j = i;
    while (j + 1 < ndiagonals && diagonals[j+1] == diagonals[i]) {
      j++;
    }
    newi = j+1;

    if (j > i) {  /* (count j - i + 1 > 1) */
      if ((count = j - i + 1) > max_count) {
	max_count = count;
      }
    }

    i = newi;
  }


  i = 0;
  while (i < ndiagonals) {
    j = i;
    while (j + 1 < ndiagonals && diagonals[j+1] == diagonals[i]) {
      j++;
    }
    newi = j+1;

    if (j > i) {  /* (count j - i + 1 > 1) */
      /* Mark duplicates */
      for (k = i+1; k <= j; k++) {
	debug1(printf("Changing status of %u from %d to %d => ",diagonals[k],status[k],DUPLICATE));
	if (status[k] == ADJACENT) {
	  (*nadjacent)--;
	}
	debug1(printf("nadjacent = %d\n",*nadjacent));
	status[k] = DUPLICATE;
      }

      if ((count = j - i + 1) == max_count) {
	/* Anchor */
	debug1(printf("Changing status of %u from %d to %d => ",diagonals[i],status[i],count));
	if (status[i] == ADJACENT) {
	  (*nadjacent)--;
	}
	debug1(printf("nadjacent = %d\n",*nadjacent));
	status[i] = count;
	(*nanchors)++;
	
	/* Mark segments ahead (adjacent is lowest status and might be overwritten) */
	while (j + 1 < ndiagonals && diagonals[j+1] < diagonals[i] + overall_max_distance) {
	  debug1(printf("Changing status of %u from %d to %d => ",diagonals[j+1],status[j+1],ADJACENT));
	  if (status[++j] != ADJACENT) {
	    status[j] = ADJACENT;
	    (*nadjacent)++;
	  }
	  debug1(printf("nadjacent = %d\n",*nadjacent));
	}
	
	/* Mark segments behind */
	j = i;
	while (j - 1 >= 0 && diagonals[j-1] + overall_max_distance > diagonals[i]) {
	  debug1(printf("Changing status of %u from %d to %d => ",diagonals[j-1],status[j-1],ADJACENT));
	  if (status[--j] == SINGLETON) {
	    status[j] = ADJACENT;
	    (*nadjacent)++;
	  }
	  debug1(printf("nadjacent = %d\n",*nadjacent));
	}
      }
    }

    i = newi;
  }

#ifdef DEBUG1
  printf("Status for %d diagonals:\n",ndiagonals);
  for (i = 0; i < ndiagonals; i++) {
    printf("i=%d %d %u\n",i,status[i],diagonals[i]);
  }
  printf("nanchors %d, nadjacent %d\n",*nanchors,*nadjacent);
#endif

#if 0
  nanchors_check = nadjacent_check = 0;
  for (i = 0; i < ndiagonals; i++) {
    if (status[i] > 0) {
      nanchors_check += 1;
    } else if (status[i] == ADJACENT) {
      nadjacent_check += 1;
    }
  }
  printf("nanchors_check %d, nadjacent_check %d\n",nanchors_check,nadjacent_check);

  if (*nanchors != nanchors_check) {
    abort();
  }
  if (*nadjacent != nadjacent_check) {
    abort();
  }
#endif

  return status;
}


/* Returns status as a positive count (if exceeds count_threshold) or
   one of the non-positive flags SINGLETON, DUPLICATE, or ADJACENT */
static int *
compute_status_threshold (int *nanchors, int *nadjacent, Univcoord_T *diagonals, int ndiagonals,
			  int count_threshold) {
  int *status;
  int i, j, k, newi;
  int count;

  status = (int *) CALLOC(ndiagonals,sizeof(int));
  *nanchors = *nadjacent = 0;

  i = 0;
  while (i < ndiagonals) {
    j = i;
    while (j + 1 < ndiagonals && diagonals[j+1] == diagonals[i]) {
      j++;
    }
    newi = j+1;

    if (j > i) {  /* (count j - i + 1 > 1) */
      /* Mark duplicates */
      for (k = i+1; k <= j; k++) {
	debug1(printf("Changing status of %u from %d to %d => ",diagonals[k],status[k],DUPLICATE));
	if (status[k] == ADJACENT) {
	  (*nadjacent)--;
	}
	debug1(printf("nadjacent = %d\n",*nadjacent));
	status[k] = DUPLICATE;
      }

      if ((count = j - i + 1) >= count_threshold) {
	/* Anchor */
	debug1(printf("Changing status of %u from %d to %d => ",diagonals[i],status[i],count));
	if (status[i] == ADJACENT) {
	  (*nadjacent)--;
	}
	debug1(printf("nadjacent = %d\n",*nadjacent));
	status[i] = count;
	(*nanchors)++;
	
	/* Mark segments ahead (adjacent is lowest status and might be overwritten) */
	while (j + 1 < ndiagonals && diagonals[j+1] < diagonals[i] + overall_max_distance) {
	  debug1(printf("Changing status of %u from %d to %d => ",diagonals[j+1],status[j+1],ADJACENT));
	  if (status[++j] != ADJACENT) {
	    status[j] = ADJACENT;
	    (*nadjacent)++;
	  }
	  debug1(printf("nadjacent = %d\n",*nadjacent));
	}
	
	/* Mark segments behind */
	j = i;
	while (j - 1 >= 0 && diagonals[j-1] + overall_max_distance > diagonals[i]) {
	  debug1(printf("Changing status of %u from %d to %d => ",diagonals[j-1],status[j-1],ADJACENT));
	  if (status[--j] == SINGLETON) {
	    status[j] = ADJACENT;
	    (*nadjacent)++;
	  }
	  debug1(printf("nadjacent = %d\n",*nadjacent));
	}
      }
    }

    i = newi;
  }

#ifdef DEBUG1
  printf("Status for %d diagonals:\n",ndiagonals);
  for (i = 0; i < ndiagonals; i++) {
    printf("%d %u\n",status[i],diagonals[i]);
  }
  printf("nanchors %d, nadjacent %d\n",*nanchors,*nadjacent);
#endif

#if 0
  nanchors_check = nadjacent_check = 0;
  for (i = 0; i < ndiagonals; i++) {
    if (status[i] > 0) {
      nanchors_check += 1;
    } else if (status[i] == ADJACENT) {
      nadjacent_check += 1;
    }
  }
  printf("nanchors_check %d, nadjacent_check %d\n",nanchors_check,nadjacent_check);

  if (*nanchors != nanchors_check) {
    abort();
  }
  if (*nadjacent != nadjacent_check) {
    abort();
  }
#endif

  return status;
}



/* Create targets for intersection, and all_records */
static struct Record_T *
make_records (int *status, Univcoord_T *diagonals, int ndiagonals, int nrecords,
#ifdef LARGE_GENOMES
	      unsigned char **stream_high_alloc, UINT4 **stream_low_alloc,
#endif
	      Univcoord_T **stream_alloc, int *streamsize_alloc, int *diagterm_alloc,
	      int nstreams, bool streams_are_diagonals_p, int querylength) {
  struct Record_T *records;
  Record_T record;
  Univcoord_T *targets, diagonal;
  int querypos, diagterm;
  int *indices;
  int streami, nindices, i, j, testj, endj, k;

  Chrnum_T chrnum = 1;
  Univcoord_T chrhigh = 0;
  Chrpos_T chrlength;
  Univcoord_T chroffset;
  Chrpos_T most_inbounds, inbounds;

#ifndef SLOW_CHR_UPDATE
  Univcoord_T goal;
  int nchromosomes_local = nchromosomes;
  Univcoord_T *chrhighs_local = chrhighs;
#endif


  records = (struct Record_T *) MALLOC(nrecords*sizeof(struct Record_T));
  targets = (Univcoord_T *) MALLOC(nrecords*sizeof(Univcoord_T));
  indices = MALLOC(nrecords*sizeof(int));

  k = 0;
  for (i = 0; i < ndiagonals; i++) {
    if (status[i] > 0) {
      /* status is a count */
      targets[k] = diagonals[i];
      records[k].diagonal = diagonals[i];
      records[k].querypos = -1; /* to indicate that it is not yet set */
      records[k].anchorp = true;
      k++;
    } else if (status[i] == ADJACENT) {
      targets[k] = diagonals[i];
      records[k].diagonal = diagonals[i];
      records[k].querypos = -1; /* to indicate that it is not yet set */
      records[k].anchorp = false;
      k++;
    }
  }
  
  for (streami = 0; streami < nstreams; streami++) {
    debug1(printf("Computing intersection on stream %d\n",streami));
    if (streams_are_diagonals_p == true) {
      nindices = Intersect_exact_indices_univcoord(indices,stream_alloc[streami],streamsize_alloc[streami],
						   targets,nrecords);
    } else {
      diagterm = diagterm_alloc[streami];
#ifdef LARGE_GENOMES
      nindices = Intersect_exact_indices_large(indices,stream_high_alloc[streami],stream_low_alloc[streami],
					       streamsize_alloc[streami],diagterm,targets,nrecords);
#else
      nindices = Intersect_exact_indices_small(indices,stream_alloc[streami],streamsize_alloc[streami],
					       diagterm,targets,nrecords);
#endif
    }

    diagterm = diagterm_alloc[streami];
    querypos = -diagterm;
    for (i = 0; i < nindices; i++) {
      k = indices[i];
      if (records[k].querypos < 0) {
	records[k].querypos = querypos;
      }
      records[k].queryend = querypos + index1part;
      debug1(printf("i=%d k=%d %u %d..%d\n",i,k,records[k].diagonal,records[k].querypos,records[k].queryend));
    }
    debug1(printf("\n"));
  }
  
#ifdef DEBUG1
  for (k = 0; k < nrecords; k++) {
    printf("k=%d %u %d..%d\n",k,records[k].diagonal,records[k].querypos,records[k].queryend);
    assert(records[k].querypos >= 0);
  }
#endif

  FREE(indices);
  FREE(targets);


  for (k = 0; k < nrecords; k++) {
    record = &(records[k]);
    diagonal = record->diagonal;

#ifdef SLOW_CHR_UPDATE
    chrhigh = Univ_IIT_update_chrnum(&chrnum,&chroffset,chrhigh,&chrlength,chromosome_iit,
				     diagonal,querylength,circular_typeint);
#else
    /* Code modeled after Univ_IIT_update_chrnum */
    if (diagonal >= chrhigh) {
      /* update chromosome bounds, based on diagonal */
      debug1(printf("\nUpdating chrhigh because diagonal %u >= chrhigh %u\n",diagonal,chrhigh));

      j = 1;
      goal = diagonal + 1;
      while (j < nchromosomes_local && chrhighs_local[j] < goal) {
	j <<= 1;			/* gallop by 2 */
      }
      if (j >= nchromosomes_local) {
	j = binary_search(j >> 1,nchromosomes_local,chrhighs_local,goal);
      } else {
	j = binary_search(j >> 1,j,chrhighs_local,goal);
      }
      chrnum += j;

      chrhigh = chrhighs[chrnum-1];
      chroffset = chroffsets[chrnum-1];
      chrlength = chrlengths[chrnum-1];
      chrhighs_local += j;
      nchromosomes_local -= j;
      debug1(printf("Got chrnum %d, chroffset %u, chrhigh %u\n",chrnum,chroffset,chrhigh));
    }

    if ((Univcoord_T) (diagonal + querylength) >= chrhigh) {
      /* Straddles two or more chromosomes */
      debug1(printf("Diagonal %u, query %d..%d straddles two or more chromosomes\n",diagonal,0,querylength));
      endj = 0;
      while (endj+1 < nchromosomes_local && diagonal + querylength >= chrhighs_local[endj+1]) {
	debug1(printf("For diagonal %u + querylength %d, advancing to chrhigh %u\n",
		      diagonal,querylength,chrhighs_local[endj+1]));
	endj++;
      }
      if (endj+1 < nchromosomes_local) {
	endj++;
      }
      debug1(printf("endj is %d\n",endj));

      /* Test first chromosome */
      j = 0;
      most_inbounds = chrhigh - diagonal;
      debug1(printf("  First chromosome inbounds: %u.  ",most_inbounds));

      /* Test middle chromosomes */
      for (testj = 1; testj < endj; testj++) {
	debug1(printf("Next inbounds: %u.  ",chrhighs_local[testj] - chroffsets[(chrnum+testj)-1]));
	if ((inbounds = chrhighs_local[testj] - chroffsets[(chrnum+testj)-1]) > most_inbounds) {
	  j = testj;
	  most_inbounds = inbounds;
	}
      }

      /* Test last chromosome */
      debug1(printf("Last inbounds: %u = %u + %d - chroffset %u\n",
		    (diagonal + querylength) - chroffsets[(chrnum+endj)-1],
		    diagonal,querylength,chroffsets[(chrnum+endj)-1]));
      if ((/*inbounds = */(diagonal + querylength) - chroffsets[(chrnum+endj)-1]) > most_inbounds) {
	j = endj;
      }

      chrnum += j;
      chrhigh = chrhighs[chrnum-1];
      chroffset = chroffsets[chrnum-1];
      chrlength = chrlengths[chrnum-1];
      chrhighs_local += j;
      nchromosomes_local -= j;

      debug1(printf("Modifying query %d..%d to ",record->querypos,record->queryend));
      if ((record->querypos = (int) (chroffset - diagonal)) < 0) {
	record->querypos = 0;
      }
      if ((record->queryend = (int) (chrhigh - diagonal)) > querylength) {
	record->queryend = querylength;
      }
      debug1(printf("%d..%d\n",record->querypos,record->queryend));
      debug1(printf("For diagonal %u + querylength %d, advancing to chrnum %d, chrhigh %u\n",
		    diagonal,querylength,chrnum,chrhigh));
    }
#endif
	  
#if 0
    /* Necessary only when a straddle has occurred */
    if (chroffset > diagonal + record->querypos) {
      record->querypos = (int) (chroffset - diagonal);
    }
    if (diagonal + record->queryend >= chrhigh) {
      record->queryend = (int) (chrhigh - diagonal);
    }
#endif

    debug1(printf("For diagonal %u, making record %d..%d at chrnum %d, chroffset %u, chrhigh %u\n",
		  diagonal,record->querypos,record->queryend,chrnum,chroffset,chrhigh));
    assert(record->querypos < record->queryend);

    record->chrnum = chrnum;
    record->chroffset = chroffset;
    record->chrhigh = chrhigh;
    record->chrlength = chrlength;
      
    record->lowpos = diagonal + record->querypos;
    record->highpos = diagonal + record->queryend;

#if 0
    if (plusp == true) {
      /* record->genomicstart = diagonal; */
      /* record->genomicend = diagonal + querylength; */
    } else {
      /* record->genomicstart = diagonal + querylength; */
      /* record->genomicend = diagonal; */
    }
#endif
  }

  return records;
}



/* gplus_5 and gminus_3 */
struct Record_T *
Segment_identify_lower (int *nrecords,
#ifdef LARGE_GENOMES
			unsigned char **positions_high, UINT4 **positions_low,
#else
			Univcoord_T **positions,
#endif		  
			int *npositions, bool *validp,

			Univcoord_T **stream_alloc, int *streamsize_alloc, int *diagterm_alloc,
			Chrpos_T max_pairlength, int querylength,
			Univcoord_T *ref_diagonals, int ref_ndiagonals) {
  struct Record_T *records;
  int query_lastpos = querylength - index1part;

  int total_npositions, used_npositions, nstreams, streami;

  int querypos;
  Univcoord_T *diagonals_alloc, *diagonals_used, *diagonals, *ptr;
  int *status;
  int ndiagonals, nanchors, nadjacent;


  debug1(printf("*** Starting Segment_identify_lower ***\n"));

  total_npositions = 0;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (validp[querypos] == false) {
      debug1(printf("Skipping batch for querypos %d with %d positions, but not valid\n",
		    querypos,npositions[querypos]));
    } else {
      assert(npositions[querypos] >= 0);
      total_npositions += npositions[querypos];
    }
  }

  if (total_npositions == 0) {
    *nrecords = 0;
    records = (struct Record_T *) NULL;

  } else {
    /* Large allocation */
    debug1(printf("Total npositions = %d\n",total_npositions));
    ptr = diagonals_alloc = (Univcoord_T *) MALLOC(total_npositions * sizeof(Univcoord_T));
    
    used_npositions = 0;
    nstreams = 0;
    for (querypos = 0; querypos <= query_lastpos; querypos++) {
      if (validp[querypos] == false) {
	debug1(printf("Skipping batch for querypos %d with %d positions, but not valid\n",
		      querypos,npositions[querypos]));
      } else if (npositions[querypos] <= 0) {
	debug1(printf("Skipping batch for querypos %d with %d positions, but not valid\n",
		      querypos,npositions[querypos]));
      } else {
	debug1(printf("Adding batch for querypos %d with %d positions\n",querypos,npositions[querypos]));
	ndiagonals = Intersect_approx_lower(ptr,
#ifdef LARGE_GENOMES
					    positions_high[querypos],positions_low[querypos],
#else
					    positions[querypos],
#endif
					    npositions[querypos],/*diagterm*/-querypos,
					    ref_diagonals,ref_ndiagonals,max_pairlength);
	if (ndiagonals > 0) {
	  /* stream_alloc[nstreams] = ptr; -- assigned in positions_used */
	  used_npositions += streamsize_alloc[nstreams] = ndiagonals;
	  ptr += ndiagonals;
	  diagterm_alloc[nstreams++] = -querypos;
	}
      }
    }
    debug1(printf("Used npositions = %d\n",used_npositions));
    
    
    if (nstreams == 0) {
      FREE(diagonals_alloc);
      *nrecords = 0;
      records = (struct Record_T *) NULL;
      
    } else {
      /* Move memory to a smaller allocation */
      ptr = diagonals_used = (Univcoord_T *) MALLOC(used_npositions * sizeof(Univcoord_T));
      memcpy(diagonals_used,diagonals_alloc,used_npositions*sizeof(Univcoord_T));
      FREE(diagonals_alloc);
      
      for (streami = 0; streami < nstreams; streami++) {
	stream_alloc[streami] = ptr;
	ptr += streamsize_alloc[streami];
      }
      
#ifdef LARGE_GENOMES
      diagonals = Merge_diagonals_uint8(&ndiagonals,stream_alloc,streamsize_alloc,nstreams);
#else
      diagonals = Merge_diagonals_uint4(&ndiagonals,stream_alloc,streamsize_alloc,nstreams);
#endif
      
      status = compute_status_threshold(&nanchors,&nadjacent,diagonals,ndiagonals,/*count_threshold*/2);
      if (nanchors == 0) {
	FREE(diagonals_used);
	*nrecords = 0;
	records = (struct Record_T *) NULL;
      } else {
	*nrecords = nanchors + nadjacent;
#ifdef LARGE_GENOMES
	records = make_records(status,diagonals,ndiagonals,*nrecords,
			       /*stream_high_alloc*/NULL,/*stream_low_alloc*/NULL,
			       stream_alloc,streamsize_alloc,diagterm_alloc,nstreams,
			       /*streams_are_diagonals_p*/true,querylength);
#else
	records = make_records(status,diagonals,ndiagonals,*nrecords,
			       stream_alloc,streamsize_alloc,diagterm_alloc,nstreams,
			       /*streams_are_diagonals_p*/true,querylength);
#endif
	FREE(diagonals_used);
      }
      
      FREE(status);
#ifndef LARGE_GENOMES
      FREE_ALIGN(diagonals);
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
      FREE_ALIGN(diagonals);
#else
      FREE(diagonals);
#endif
    }
  }

  return records;
}


/* gplus_3 and gminus_5 */
struct Record_T *
Segment_identify_higher (int *nrecords,
#ifdef LARGE_GENOMES
			 unsigned char **positions_high, UINT4 **positions_low,
#else
			 Univcoord_T **positions,
#endif		  
			 int *npositions, bool *validp,

			 Univcoord_T **stream_alloc, int *streamsize_alloc, int *diagterm_alloc,
			 Chrpos_T max_pairlength, int querylength,
			 Univcoord_T *ref_diagonals, int ref_ndiagonals) {
  struct Record_T *records;
  int query_lastpos = querylength - index1part;

  int total_npositions, used_npositions, nstreams, streami;
  int querypos;
  Univcoord_T *diagonals_alloc, *diagonals_used, *diagonals, *ptr;
  int *status;
  int ndiagonals, nanchors, nadjacent;


  debug(printf("*** Starting Segment_identify_higher ***\n"));

  total_npositions = 0;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (validp[querypos] == false) {
      debug1(printf("Skipping batch for querypos %d with %d positions, but not valid\n",
		    querypos,npositions[querypos]));
    } else {
      assert(npositions[querypos] >= 0);
      total_npositions += npositions[querypos];
    }
  }

  if (total_npositions == 0) {
    *nrecords = 0;
    records = (struct Record_T *) NULL;
  } else {
    /* Large allocation */
    ptr = diagonals_alloc = (Univcoord_T *) MALLOC(total_npositions * sizeof(Univcoord_T));
    
    used_npositions = 0;
    nstreams = 0;
    for (querypos = 0; querypos <= query_lastpos; querypos++) {
      if (validp[querypos] == false) {
	debug1(printf("Skipping batch for querypos %d with %d positions, but not valid\n",
		      querypos,npositions[querypos]));
      } else if (npositions[querypos] <= 0) {
	debug1(printf("Skipping batch for querypos %d with %d positions, but not valid\n",
		      querypos,npositions[querypos]));
      } else {
	debug1(printf("Adding batch for querypos %d with %d positions\n",querypos,npositions[querypos]));
	ndiagonals = Intersect_approx_higher(ptr,
#ifdef LARGE_GENOMES
					     positions_high[querypos],positions_low[querypos],
#else
					     positions[querypos],
#endif
					     npositions[querypos],/*diagterm*/-querypos,
					     ref_diagonals,ref_ndiagonals,max_pairlength);
	if (ndiagonals > 0) {
	  /* stream_alloc[nstreams] = ptr; -- assigned in positions_used */
	  used_npositions += streamsize_alloc[nstreams] = ndiagonals;
	  ptr += ndiagonals;
	  diagterm_alloc[nstreams++] = -querypos;
	}
      }
    }
    debug1(printf("Used npositions = %d\n",used_npositions));
    
    
    if (nstreams == 0) {
      FREE(diagonals_alloc);
      *nrecords = 0;
      records = (struct Record_T *) NULL;
      
    } else {
      /* Move memory to a smaller allocation */
      ptr = diagonals_used = (Univcoord_T *) MALLOC(used_npositions * sizeof(Univcoord_T));
      memcpy(diagonals_used,diagonals_alloc,used_npositions*sizeof(Univcoord_T));
      FREE(diagonals_alloc);
      
      for (streami = 0; streami < nstreams; streami++) {
	stream_alloc[streami] = ptr;
	ptr += streamsize_alloc[streami];
      }
      
#ifdef LARGE_GENOMES
      diagonals = Merge_diagonals_uint8(&ndiagonals,stream_alloc,streamsize_alloc,nstreams);
#else
      diagonals = Merge_diagonals_uint4(&ndiagonals,stream_alloc,streamsize_alloc,nstreams);
#endif
      
      status = compute_status_threshold(&nanchors,&nadjacent,diagonals,ndiagonals,/*count_threshold*/2);
      if (nanchors == 0) {
	FREE(diagonals_used);
	*nrecords = 0;
	records = (struct Record_T *) NULL;
      } else {
	*nrecords = nanchors + nadjacent;
#ifdef LARGE_GENOMES
	records = make_records(status,diagonals,ndiagonals,*nrecords,
			       /*stream_high_alloc*/NULL,/*stream_low_alloc*/NULL,
			       stream_alloc,streamsize_alloc,diagterm_alloc,nstreams,
			       /*streams_are_diagonals_p*/true,querylength);
#else
	records = make_records(status,diagonals,ndiagonals,*nrecords,
			       stream_alloc,streamsize_alloc,diagterm_alloc,nstreams,
			       /*streams_are_diagonals_p*/true,querylength);
#endif
	FREE(diagonals_used);
      }
      
      FREE(status);
#ifndef LARGE_GENOMES
      FREE_ALIGN(diagonals);
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
      FREE_ALIGN(diagonals);
#else
      FREE(diagonals);
#endif
    }
  }

  return records;
}


struct Record_T *
Segment_identify (int *nrecords,
#ifdef LARGE_GENOMES
		  unsigned char **positions_high,
#endif		  
		  UINT4 **positions, int *npositions, bool *validp,

#ifdef LARGE_GENOMES
		  unsigned char **stream_high_alloc, UINT4 **stream_low_alloc,
#else
		  Univcoord_T **stream_alloc,
#endif
		  int *streamsize_alloc, int *diagterm_alloc,
		  int querylength, int sizelimit) {
  struct Record_T *records;
  int query_lastpos = querylength - index1part;

  int total_npositions = 0, nstreams;

  int querypos;
  Univcoord_T *diagonals;
  int ndiagonals;
  int *status;
  int nanchors, nadjacent;


  debug(printf("*** Starting Segment_identify ***\n"));

  nstreams = 0;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (sizelimit > 0 && npositions[querypos] > sizelimit) {
      debug1(printf("Skipping batch for querypos %d with %d positions, over sizelimit %d\n",
		    querypos,npositions[querypos],sizelimit));
    } else if (npositions[querypos] > ABSOLUTE_LIMIT) {
      debug1(printf("Skipping batch for querypos %d with %d positions, over absolute_limit %d\n",
		    querypos,npositions[querypos],ABSOLUTE_LIMIT));
    } else if (validp[querypos] == false) {
      debug1(printf("Skipping batch for querypos %d with %d positions, but not valid\n",
		    querypos,npositions[querypos]));
    } else if (npositions[querypos] > 0) {
      debug1(printf("Adding batch for querypos %d with %d positions\n",querypos,npositions[querypos]));
      /* Previously had a diagterm_list here filled with
	 (querylength - querypos) for plus or (querylength -
	 querypos) for minus when we called Indexdb_fill_inplace,
	 but now we use Indexdb_read_with_diagterm */
      
#ifdef LARGE_GENOMES
      stream_high_alloc[nstreams] = positions_high[querypos];
      stream_low_alloc[nstreams] = positions[querypos];
#else
      stream_alloc[nstreams] = positions[querypos];
#endif      
      streamsize_alloc[nstreams] = npositions[querypos];
      diagterm_alloc[nstreams] = -querypos;

      total_npositions += npositions[querypos];
      nstreams++;
    } else {
      debug1(printf("Not adding batch for querypos %d with %d positions\n",querypos,npositions[querypos]));
    }
  }
  debug1(printf("Initial total_npositions = %d\n",total_npositions));


  if (nstreams == 0) {
    *nrecords = 0;
    records = (struct Record_T *) NULL;
  } else {
#ifdef LARGE_GENOMES
    diagonals = Merge_diagonals_large(&ndiagonals,stream_high_alloc,stream_low_alloc,
				      streamsize_alloc,diagterm_alloc,nstreams);
#else
    diagonals = Merge_diagonals(&ndiagonals,stream_alloc,streamsize_alloc,diagterm_alloc,nstreams);
#endif

    status = compute_status_greedy(&nanchors,&nadjacent,diagonals,ndiagonals);
    if (nanchors == 0) {
      *nrecords = 0;
      records = (struct Record_T *) NULL;
    } else {
      *nrecords = nanchors + nadjacent;
#ifdef LARGE_GENOMES
      records = make_records(status,diagonals,ndiagonals,*nrecords,
			     stream_high_alloc,stream_low_alloc,
			     /*stream_alloc*/NULL,streamsize_alloc,diagterm_alloc,nstreams,
			     /*streams_are_diagonals_p*/false,querylength);
#else
      records = make_records(status,diagonals,ndiagonals,*nrecords,
			     stream_alloc,streamsize_alloc,diagterm_alloc,nstreams,
			     /*streams_are_diagonals_p*/false,querylength);
#endif
    }
    FREE(status);

#ifndef LARGE_GENOMES
    FREE_ALIGN(diagonals);
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
    FREE_ALIGN(diagonals);
#else
    FREE(diagonals);
#endif
  }

  return records;
}


#if 0
/* Segment chaining */
/* Modified from previous convert_plus_segments_to_gmap */
/* Generic, for either plus or minus genomic strand */
static List_T
solve_filtered (int *found_score, List_T hits, List_T *complete_paths, int querylength,
		Record_T *segments, int nsegments, Intlist_T anchors,
		char *queryptr, Compress_T query_compress,
		bool plusp, int genestrand, bool paired_end_p, bool first_read_p, int level) {
  int nmisses_allowed = 5;
  Intlist_T p;

  Record_T anchor_segment, segment;
  int anchori, startk, endk, n, i, firstj, lastj, j, k;

  Univdiag_T middle_diagonal;
  List_T complete_path, right_diagonals, left_diagonals;

  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;
  Chrnum_T chrnum;

  Record_T *sorted5, *sorted5_allocated, *sorted3, *sorted3_allocated;

  Univcoord_T left;
  bool foundp;


  debug(printf("Entered solve_filtered\n"));

  if (nsegments > 0) {
#ifdef HAVE_ALLOCA
    if (nsegments < MAX_ALLOCATION) {
      sorted5_allocated = (Record_T *) ALLOCA(nsegments*sizeof(Record_T));
      sorted3_allocated = (Record_T *) ALLOCA(nsegments*sizeof(Record_T));
    } else {
      sorted5_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
      sorted3_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
    }
#else
    sorted5_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
    sorted3_allocated = (Record_T *) MALLOC(nsegments*sizeof(Record_T));
#endif
  }

  for (p = anchors; p != NULL; p = Intlist_next(p)) {
    anchori = Intlist_head(p);
    anchor_segment = segments[anchori];
    left = anchor_segment->diagonal /*- querylength*/; /* NEW FORMULA: Corresponds to querypos 0 */

    /* left is safer than anchor_segment->lowpos, in case trim goes further left than lowpos */
    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
    
    j = anchori - 1;
    while (j >= 0 && left <= segments[j]->diagonal + overall_max_distance) {
      segments[j]->lowpos = segments[j]->diagonal + segments[j]->querypos;
      segments[j]->highpos = segments[j]->diagonal + segments[j]->queryend;
      j--;
    }
    startk = j + 1;

    j = anchori + 1;
    while (j < nsegments && segments[j]->diagonal <= left + overall_max_distance) {
      segments[j]->lowpos = segments[j]->diagonal + segments[j]->querypos;
      segments[j]->highpos = segments[j]->diagonal + segments[j]->queryend;
      j++;
    }
    endk = j - 1;

    debug13(printf("Found segments %d to %d inclusive for anchor #%d, %d..%d %u, chrpos %u (range %u..%u)\n",
		   startk,endk,anchori,anchor_segment->querypos,anchor_segment->queryend,
		   anchor_segment->diagonal,anchor_segment->diagonal - chroffset,
		   anchor_segment->lowpos,anchor_segment->highpos));
#ifdef DEBUG13
    for (k = startk; k <= endk; k++) {
      printf("  #%d: %d..%d %u\n",k,segments[k]->querypos,segments[k]->queryend,segments[k]->diagonal);
    }
    printf("\n");
#endif

    n = endk - startk + 1;
    debug13(printf("n = %d\n",n));
    sorted5 = &(sorted5_allocated[startk]);
    sorted3 = &(sorted3_allocated[startk]);


    /* middle_path = (List_T) NULL; */
    right_diagonals = left_diagonals = (List_T) NULL;

    /* Add diagonals to left (low) side (querypos5) */
    for (k = startk, i = 0; k <= endk; k++) {
      sorted5[i++] = segments[k];
    }
    qsort(sorted5,n,sizeof(Record_T),Record_querypos5_ascending_cmp);
#ifdef DEBUG13
    printf("Sorted by querypos5\n");
    for (i = 0; i < n; i++) {
      printf("  #%d: %d..%d %u\n",i,sorted5[i]->querypos,sorted5[i]->queryend,sorted5[i]->diagonal);
    }
#endif

    /* Skip obvious ones */
    lastj = 0;
    while (lastj < n && sorted5[lastj]->querypos < anchor_segment->querypos) {
      lastj++;
    }
#ifdef DEBUG13
    printf("On the querypos5 side, considering:\n");
    for (i = 0; i < lastj; i++) {
      printf("  #%d: %d..%d left %u, range %u..%u\n",
	     i,sorted5[i]->querypos,sorted5[i]->queryend,sorted5[i]->diagonal,sorted5[i]->lowpos,sorted5[i]->highpos);
    }
#endif

    /* Go leftward to favor shorter splices */
    debug13(printf("Trying to append to start of anchor segment\n"));
    for (j = lastj - 1; j >= 0; --j) {
      segment = sorted5[j];
      debug13(printf("ANCHOR SEGMENT %d..%d left %u (range %u..%u)\n",
		     anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal,anchor_segment->lowpos,anchor_segment->highpos));
      debug13(printf("(1) TRIAL SEGMENT ON 5' %d..%d left %u (range %u..%u)\n",
		     segment->querypos,segment->queryend,segment->diagonal,segment->lowpos,segment->highpos));

      if (segment->lowpos >= anchor_segment->lowpos + max_deletionlen) {
	debug13(printf("(2) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE LEFT OF ANCHOR\n"));
	
      } else if (segment->lowpos >= anchor_segment->lowpos && segment->highpos <= anchor_segment->highpos) {
	debug13(printf("(2) SKIPPING, SINCE ANCHOR SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));
	  
      } else if (anchor_segment->lowpos >= segment->lowpos && anchor_segment->highpos <= segment->highpos) {
	debug13(printf("(2) SKIPPING, SINCE TRIAL SUBSUMES ANCHOR (PROBABLE GENOMIC REPEAT)\n"));

#if 0
      } else if (segment->highpos < chroffset) {
	/* Now handled by path-solve, which accounts for the best chromosome for the anchor segment */
	debug13(printf("Not allowing left diagonal that goes into next chromosome\n"));
#endif

      } else {
	debug13(printf("Candidate for start diagonal\n"));
	if (segment->diagonal + segment->querypos >= left + anchor_segment->querypos) {
	  debug13(printf("Not allowing left diagonal that doesn't advance min_genomepos\n"));
	} else if ((segment->queryend - anchor_segment->querypos) > (segment->queryend - segment->querypos) / 2) {
	  debug13(printf("Not allowing left diagonal that mainly overlaps middle diagonal\n"));
	} else if (Record_overlap_p(segment,anchor_segment) == true) {
	  debug13(printf("Not allowing left diagonal that overlaps anchor diagonal\n"));
	} else {
	  debug13(printf("Left diagonal is allowable: %u..%u vs chromosome %u..%u\n",
			 segment->lowpos,segment->highpos,chroffset,chrhigh));
	  left_diagonals = List_push(left_diagonals,
				     (void *) Univdiag_new(segment->querypos,segment->queryend,segment->diagonal));
	}
      }
    }


    /* Add diagonals to right (high) side (querypos3) */
    for (k = startk, i = 0; k <= endk; k++) {
      sorted3[i++] = segments[k];
    }
    qsort(sorted3,n,sizeof(Record_T),Record_querypos3_ascending_cmp);
#ifdef DEBUG13
    printf("Sorted by querypos3\n");
    for (i = 0; i < n; i++) {
      printf("  #%d: %d..%d %u\n",i,sorted3[i]->querypos,sorted3[i]->queryend,sorted3[i]->diagonal);
    }
#endif

    /* Skip obvious ones */
    firstj = n - 1;
    while (firstj >= 0 && sorted3[firstj]->queryend > anchor_segment->queryend) {
      firstj--;
    }
#ifdef DEBUG13
    printf("On the querypos3 side, considering:\n");
    for (j = firstj + 1; j < n; j++) {
      printf("  #%d: %d..%d left %u, range %u..%u\n",
	     i,sorted3[j]->querypos,sorted3[j]->queryend,sorted3[j]->diagonal,sorted3[j]->lowpos,sorted3[j]->highpos);
    }
#endif

    /* Go rightward to favor shorter splices */
    debug13(printf("Trying to append to end of anchor segment\n"));
    for (j = firstj + 1; j < n; j++) {
      segment = sorted3[j];
      debug13(printf("ANCHOR SEGMENT %d..%d left %u (range %u..%u)\n",
		     anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal,anchor_segment->lowpos,anchor_segment->highpos));
      debug13(printf("(2) TRIAL SEGMENT ON 3' %d..%d left %u (range %u..%u)\n",
		     segment->querypos,segment->queryend,segment->diagonal,segment->lowpos,segment->highpos));

      if (segment->highpos + max_deletionlen <= anchor_segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE RIGHT OF ANCHOR\n"));
	
      } else if (segment->lowpos >= anchor_segment->lowpos && segment->highpos <= anchor_segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE ANCHOR SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));

      } else if (anchor_segment->lowpos >= segment->lowpos && anchor_segment->highpos <= segment->highpos) {
	debug13(printf("(4) SKIPPING, SINCE TRIAL SUBSUMES ANCHOR (PROBABLE GENOMIC REPEAT)\n"));

#if 0
      } else if (segment->lowpos >= chrhigh) {
	/* Now handled by path-solve, which accounts for the best chromosome for the anchor segment */
	debug13(printf("Not allowing right diagonal that goes into next chromosome\n"));
#endif

      } else {
	debug13(printf("Candidate for end diagonal\n"));
	if (segment->diagonal + segment->queryend <= left + anchor_segment->queryend) {
	  debug13(printf("Not allowing right diagonal that doesn't advance max_genomepos\n"));
	} else if ((anchor_segment->queryend - segment->querypos) > (segment->queryend - segment->querypos) / 2) {
	  debug13(printf("Not allowing right diagonal that mainly overlaps anchor diagonal\n"));
	} else if (Record_overlap_p(segment,anchor_segment) == true) {
	  debug13(printf("Not allowing right diagonal that overlaps anchor diagonal\n"));
	} else {
	  debug13(printf("Right diagonal is allowable: %u..%u vs chromosome %u..%u\n",
			 segment->lowpos,segment->highpos,chroffset,chrhigh));
	  right_diagonals = List_push(right_diagonals,
				      (void *) Univdiag_new(segment->querypos,segment->queryend,segment->diagonal));
	}
      }
    }

    /* Compute middle diagonal */
    debug13(printf("Anchor diagonal %llu, querypos %d..%d\n",
		   (unsigned long long) anchor_segment->diagonal,anchor_segment->querypos,anchor_segment->queryend));
    middle_diagonal = Univdiag_new(anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal);
    
    /* TODO: Use localdb to add left and right diagonals.  Then modify
       Path_solve_from_diagonals to handle a path, and call that */
    debug13(printf("SOLVING A PATH WITH %d LEFT, ONE MIDDLE, and %d RIGHT DIAGONALS\n\n",
		   List_length(left_diagonals),List_length(right_diagonals)));
    
    hits = Path_solve_from_diagonals(&foundp,&(*found_score_overall),&(*found_score_within_trims),hits,
				     middle_diagonal->univdiagonal,middle_diagonal->qstart,middle_diagonal->qend,
				     right_diagonals,left_diagonals,queryptr,querylength,
				     mismatch_positions_alloc,spliceinfo,stream_alloc,streamsize_alloc,
				     query_compress,chrnum,chroffset,chrhigh,chrlength,plusp,genestrand,
				     /*nmismatches_allowed*/nmisses_allowed,paired_end_p,first_read_p,
				     intlistpool,univcoordlistpool,listpool,univdiagpool,
				     hitlistpool,/*method*/SEGMENT,level);
    Univdiag_free(&middle_diagonal);
    Univdiag_gc(&right_diagonals);
    Univdiag_gc(&left_diagonals);
    
    if (foundp == false) {
      /* *complete_paths = List_push(*complete_paths,(void *) complete_path); */
    }
  }

  if (nsegments > 0) {
#ifdef HAVE_ALLOCA
    if (nsegments < MAX_ALLOCATION) {
      FREEA(sorted3_allocated);
      FREEA(sorted5_allocated);
    } else {
      FREE(sorted3_allocated);
      FREE(sorted5_allocated);
    }
#else
    FREE(sorted3_allocated);
    FREE(sorted5_allocated);
#endif
  }

  debug(printf("Returning from solve_filtered\n"));

  return hits;
}
#endif


#if 0
/* Assumes that stage1 has been filled with all positions */
void
Segment_search_filtered (int *found_score, List_T *plus_hits, List_T *minus_hits,
			 List_T *plus_complete_paths, List_T *minus_complete_paths,

			 Record_T *plus_records, int plus_nrecords, Intlist_T plus_anchors, 
			 Record_T *minus_records, int minus_nrecords, Intlist_T minus_anchors,

			 char *queryuc_ptr, char *queryrc, int querylength,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 int genestrand, bool paired_end_p, int level) {

  if (plus_records != NULL) {
    *plus_hits = solve_filtered(&(*found_score),*plus_hits,&(*plus_complete_paths),querylength,
				plus_records,plus_nrecords,plus_anchors,
				/*queryptr*/queryuc_ptr,query_compress_fwd,
				/*plusp*/true,genestrand,paired_end_p,level);
  }

  if (minus_records != NULL) {
    *minus_hits = solve_filtered(&(*found_score),*minus_hits,&(*minus_complete_paths),querylength,
				 minus_records,minus_nrecords,minus_anchors,
				 /*queryptr*/queryrc,query_compress_rev,
				 /*plusp*/false,genestrand,paired_end_p,level);
  }

  return;
}
#endif


/* Generic, for either plus or minus genomic strand */
static List_T
solve_all (int *found_score_overall, int *found_score_within_trims,
	   List_T hits, int querylength, struct Record_T *records, int nrecords,
	   char *queryptr, int *mismatch_positions_alloc, Spliceinfo_T spliceinfo,
	   Univcoord_T **stream_alloc, int *streamsize_alloc, Compress_T query_compress,
	   bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
	   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	   Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
	   Method_T method, int level) {
  int nmisses_allowed = 5;

  Record_T anchor_segment, segment;
  int anchori, startk, endk, n, i, firstj, lastj, j, k;

  Univdiag_T middle_diagonal;
  List_T right_diagonals, left_diagonals;

  Chrnum_T chrnum;

  Record_T *sorted5, *sorted5_allocated, *sorted3, *sorted3_allocated;

  Univcoord_T left;
  bool foundp;


  debug(printf("Entered solve_all\n"));
  assert(nrecords > 0);

  sorted5_allocated = (Record_T *) MALLOC(nrecords*sizeof(Record_T));
  sorted3_allocated = (Record_T *) MALLOC(nrecords*sizeof(Record_T));


  for (anchori = 0; anchori < nrecords; anchori++) {
    if (records[anchori].anchorp == true) {
      anchor_segment = &(records[anchori]);
      left = anchor_segment->diagonal /*- querylength*/; /* NEW FORMULA: Corresponds to querypos 0 */

      /* left is safer than anchor_segment->lowpos, in case trim goes further left than lowpos */
      chrnum = anchor_segment->chrnum;
    
      j = anchori - 1;
      while (j >= 0 && left <= records[j].diagonal + overall_max_distance && records[j].chrnum == chrnum) {
	j--;
      }
      startk = j + 1;

      j = anchori + 1;
      while (j < nrecords && records[j].diagonal <= left + overall_max_distance && records[j].chrnum == chrnum) {
	j++;
      }
      endk = j - 1;

      debug13(printf("Found records %d to %d inclusive for anchor #%d, %d..%d %u, chrpos %u (range %u..%u)\n",
		     startk,endk,anchori,anchor_segment->querypos,anchor_segment->queryend,
		     anchor_segment->diagonal,anchor_segment->diagonal - anchor_segment->chroffset,
		     anchor_segment->lowpos,anchor_segment->highpos));
#ifdef DEBUG13
      for (k = startk; k <= endk; k++) {
	printf("  #%d: %d..%d %u\n",k,records[k].querypos,records[k].queryend,records[k].diagonal);
      }
      printf("\n");
#endif

      n = endk - startk + 1;
      debug13(printf("n = %d\n",n));
      sorted5 = &(sorted5_allocated[startk]);
      sorted3 = &(sorted3_allocated[startk]);

      /* middle_path = (List_T) NULL; */
      right_diagonals = left_diagonals = (List_T) NULL;

      /* Add diagonals to left (low) side (querypos5) */
      for (k = startk, i = 0; k <= endk; k++) {
	sorted5[i++] = &(records[k]);
      }
      qsort(sorted5,n,sizeof(Record_T),Record_querypos5_ascending_cmp);
#ifdef DEBUG13
      printf("Sorted by querypos5\n");
      for (i = 0; i < n; i++) {
	printf("  #%d: %d..%d %u\n",i,sorted5[i]->querypos,sorted5[i]->queryend,sorted5[i]->diagonal);
      }
#endif

      /* Skip obvious ones */
      lastj = 0;
      while (lastj < n && sorted5[lastj]->querypos < anchor_segment->querypos) {
	lastj++;
      }
#ifdef DEBUG13
      printf("On the querypos5 side, considering:\n");
      for (i = 0; i < lastj; i++) {
	printf("  #%d: %d..%d left %u, range %u..%u\n",
	       i,sorted5[i]->querypos,sorted5[i]->queryend,sorted5[i]->diagonal,sorted5[i]->lowpos,sorted5[i]->highpos);
      }
#endif

      /* Go leftward to favor shorter splices */
      debug13(printf("Trying to append to start of anchor segment\n"));
      for (j = lastj - 1; j >= 0; --j) {
	segment = sorted5[j];
	debug13(printf("ANCHOR SEGMENT %d..%d left %u (range %u..%u)\n",
		       anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal,anchor_segment->lowpos,anchor_segment->highpos));
	debug13(printf("(3) TRIAL SEGMENT ON 5' %d..%d left %u (range %u..%u)\n",
		       segment->querypos,segment->queryend,segment->diagonal,segment->lowpos,segment->highpos));
	
	if (segment->lowpos >= anchor_segment->lowpos + max_deletionlen) {
	  debug13(printf("(2) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE LEFT OF ANCHOR\n"));
	  
	} else if (segment->lowpos >= anchor_segment->lowpos && segment->highpos <= anchor_segment->highpos) {
	  debug13(printf("(2) SKIPPING, SINCE ANCHOR SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));
	  
	} else if (anchor_segment->lowpos >= segment->lowpos && anchor_segment->highpos <= segment->highpos) {
	  debug13(printf("(2) SKIPPING, SINCE TRIAL SUBSUMES ANCHOR (PROBABLE GENOMIC REPEAT)\n"));
	  
#if 0
	} else if (segment->highpos < anchor_segment->chroffset) {
	  /* Now handled by path-solve, which accounts for the best chromosome for the anchor segment */
	  debug13(printf("Not allowing left diagonal that goes into next chromosome\n"));
#endif
	  
	} else {
	  debug13(printf("Candidate for start diagonal\n"));
	  if (segment->diagonal + segment->querypos >= left + anchor_segment->querypos) {
	    debug13(printf("Not allowing left diagonal that doesn't advance min_genomepos\n"));
	  } else if ((segment->queryend - anchor_segment->querypos) > (segment->queryend - segment->querypos) / 2) {
	    debug13(printf("Not allowing left diagonal that mainly overlaps middle diagonal\n"));
	  } else if (Record_overlap_p(segment,anchor_segment) == true) {
	    debug13(printf("Not allowing left diagonal that overlaps anchor diagonal\n"));
	  } else {
	    debug13(printf("Left diagonal is allowable: %u..%u vs chromosome %u..%u\n",
			   segment->lowpos,segment->highpos,segment->chroffset,segment->chrhigh));
	    left_diagonals = Univdiagpool_push(left_diagonals,univdiagpool,
					       segment->querypos,segment->queryend,segment->diagonal);
	  }
	}
      }


      /* Add diagonals to right (high) side (querypos3) */
      for (k = startk, i = 0; k <= endk; k++) {
	sorted3[i++] = &(records[k]);
      }
      qsort(sorted3,n,sizeof(Record_T),Record_querypos3_ascending_cmp);
#ifdef DEBUG13
      printf("Sorted by querypos3\n");
      for (i = 0; i < n; i++) {
	printf("  #%d: %d..%d %u\n",i,sorted3[i]->querypos,sorted3[i]->queryend,sorted3[i]->diagonal);
      }
#endif
      
      /* Skip obvious ones */
      firstj = n - 1;
      while (firstj >= 0 && sorted3[firstj]->queryend > anchor_segment->queryend) {
	firstj--;
      }
#ifdef DEBUG13
      printf("On the querypos3 side, considering:\n");
      for (j = firstj + 1; j < n; j++) {
	printf("  #%d: %d..%d left %u, range %u..%u\n",
	       i,sorted3[j]->querypos,sorted3[j]->queryend,sorted3[j]->diagonal,sorted3[j]->lowpos,sorted3[j]->highpos);
      }
#endif
      
      /* Go rightward to favor shorter splices */
      debug13(printf("Trying to append to end of anchor segment\n"));
      for (j = firstj + 1; j < n; j++) {
	segment = sorted3[j];
	debug13(printf("ANCHOR SEGMENT %d..%d left %u (range %u..%u)\n",
		       anchor_segment->querypos,anchor_segment->queryend,anchor_segment->diagonal,anchor_segment->lowpos,anchor_segment->highpos));
	debug13(printf("(4) TRIAL SEGMENT ON 3' %d..%d left %u (range %u..%u)\n",
		       segment->querypos,segment->queryend,segment->diagonal,segment->lowpos,segment->highpos));
	
	if (segment->highpos + max_deletionlen <= anchor_segment->highpos) {
	  debug13(printf("(4) SKIPPING, SINCE TRIAL DOESN'T ADD NUCLEOTIDES TO THE RIGHT OF ANCHOR\n"));
	  
	} else if (segment->lowpos >= anchor_segment->lowpos && segment->highpos <= anchor_segment->highpos) {
	  debug13(printf("(4) SKIPPING, SINCE ANCHOR SUBSUMES TRIAL (PROBABLE GENOMIC REPEAT)\n"));
	  
	} else if (anchor_segment->lowpos >= segment->lowpos && anchor_segment->highpos <= segment->highpos) {
	  debug13(printf("(4) SKIPPING, SINCE TRIAL SUBSUMES ANCHOR (PROBABLE GENOMIC REPEAT)\n"));
	  
#if 0
	} else if (segment->lowpos >= chrhigh) {
	  /* Now handled by path-solve, which accounts for the best chromosome for the anchor segment */
	  debug13(printf("Not allowing right diagonal that goes into next chromosome\n"));
#endif
	  
	} else {
	  debug13(printf("Candidate for end diagonal\n"));
	  if (segment->diagonal + segment->queryend <= left + anchor_segment->queryend) {
	    debug13(printf("Not allowing right diagonal that doesn't advance max_genomepos\n"));
	  } else if ((anchor_segment->queryend - segment->querypos) > (segment->queryend - segment->querypos) / 2) {
	    debug13(printf("Not allowing right diagonal that mainly overlaps anchor diagonal\n"));
	  } else if (Record_overlap_p(segment,anchor_segment) == true) {
	    debug13(printf("Not allowing right diagonal that overlaps anchor diagonal\n"));
	  } else {
	    debug13(printf("Right diagonal is allowable: %u..%u vs chromosome %u..%u\n",
			   segment->lowpos,segment->highpos,segment->chroffset,segment->chrhigh));
	    right_diagonals = Univdiagpool_push(right_diagonals,univdiagpool,
						segment->querypos,segment->queryend,segment->diagonal);
	  }
	}
      }
      
      /* Compute middle diagonal */
      debug13(printf("Anchor diagonal %llu, querypos %d..%d\n",
		     (unsigned long long) anchor_segment->diagonal,anchor_segment->querypos,anchor_segment->queryend));
      middle_diagonal = Univdiag_new(univdiagpool,anchor_segment->querypos,anchor_segment->queryend,
				     anchor_segment->diagonal);
      
      /* TODO: Use localdb to add left and right diagonals.  Then modify
	 Path_solve_from_diagonals to handle a path, and call that */
      debug13(printf("SOLVING A PATH WITH %d LEFT, ONE MIDDLE, and %d RIGHT DIAGONALS\n\n",
		     List_length(left_diagonals),List_length(right_diagonals)));
      
      hits = Path_solve_from_diagonals(&foundp,&(*found_score_overall),&(*found_score_within_trims),hits,
				       middle_diagonal->univdiagonal,middle_diagonal->qstart,middle_diagonal->qend,
				       right_diagonals,left_diagonals,queryptr,querylength,
				       mismatch_positions_alloc,spliceinfo,stream_alloc,streamsize_alloc,
				       query_compress,chrnum,anchor_segment->chroffset,anchor_segment->chrhigh,anchor_segment->chrlength,
				       plusp,genestrand,/*nmismatches_allowed*/nmisses_allowed,paired_end_p,first_read_p,
				       intlistpool,univcoordlistpool,listpool,univdiagpool,
				       hitlistpool,method,level);
      /* Univdiag_free(&middle_diagonal); -- allocated by Univdiagpool_push */
      /* Univdiag_gc(&right_diagonals); -- allocated by Univdiagpool_push */
      /* Univdiag_gc(&left_diagonals); -- allocated by Univdiagpool_push */
    }
  }

  FREE(sorted3_allocated);
  FREE(sorted5_allocated);

  debug(printf("Returning from solve_all\n"));

  return hits;
}


void
Segment_search_all (int *found_score_overall, int *found_score_within_trims,
		    List_T *plus_hits, List_T *minus_hits,

		    struct Record_T *plus_records, int plus_nrecords, 
		    struct Record_T *minus_records, int minus_nrecords,

		    char *queryuc_ptr, char *queryrc, int querylength,
		    int *mismatch_positions_alloc, Spliceinfo_T spliceinfo,
		    Univcoord_T **stream_alloc, int *streamsize_alloc,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    int genestrand, bool paired_end_p, bool first_read_p,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
		    Method_T method, int level) {

  if (plus_records != NULL) {
    *plus_hits = solve_all(&(*found_score_overall),&(*found_score_within_trims),
			   *plus_hits,querylength,plus_records,plus_nrecords,
			   /*queryptr*/queryuc_ptr,mismatch_positions_alloc,spliceinfo,
			   stream_alloc,streamsize_alloc,query_compress_fwd,
			   /*plusp*/true,genestrand,paired_end_p,first_read_p,
			   intlistpool,univcoordlistpool,listpool,univdiagpool,
			   hitlistpool,method,level);
  }

  if (minus_records != NULL) {
    *minus_hits = solve_all(&(*found_score_overall),&(*found_score_within_trims),
			    *minus_hits,querylength,minus_records,minus_nrecords,
			    /*queryptr*/queryrc,mismatch_positions_alloc,spliceinfo,
			    stream_alloc,streamsize_alloc,query_compress_rev,
			    /*plusp*/false,genestrand,paired_end_p,first_read_p,
			    intlistpool,univcoordlistpool,listpool,univdiagpool,
			    hitlistpool,method,level);
  }

  return;
}




#if 0
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
#endif




void
Segment_search_setup (int index1part_in, int index1interval_in,
		      int max_anchors_in, Univ_IIT_T chromosome_iit_in, int nchromosomes_in,
		      int circular_typeint_in, Mode_T mode_in,
		      Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		      Chrpos_T *splicedists_in, int nsplicesites_in,
		      int max_middle_deletions, Chrpos_T shortsplicedist_in) {

  index1part = index1part_in;
  index1interval = index1interval_in;
  max_anchors = max_anchors_in;

  nchromosomes = nchromosomes_in;
  chromosome_iit = chromosome_iit_in;
  circular_typeint = circular_typeint_in;

  Univ_IIT_intervals_setup(&chroffsets,&chrhighs,&chrlengths,chromosome_iit,nchromosomes,circular_typeint);

#ifdef HAVE_64_BIT
  leftreadshift = 64 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#else
  leftreadshift = 32 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#endif

  mode = mode_in;

  splicesites = splicesites_in;
  splicetypes = splicetypes_in;
  splicedists = splicedists_in;
  nsplicesites = nsplicesites_in;

  shortsplicedist = shortsplicedist_in;

  overall_max_distance = shortsplicedist;
  max_deletionlen = max_middle_deletions;
  if (max_middle_deletions > (int) overall_max_distance) {
    overall_max_distance = max_middle_deletions;
  }
#if 0
  if (max_middle_insertions_default > (int) overall_max_distance) {
    overall_max_distance = max_middle_insertions_default;
  }
#endif

  return;
}

void
Segment_search_cleanup () {
  FREE(chroffsets);
  FREE(chrhighs);
  FREE(chrlengths);
  return;
}
