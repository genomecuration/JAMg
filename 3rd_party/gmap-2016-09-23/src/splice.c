static char rcsid[] = "$Id: splice.c 197917 2016-09-16 13:39:50Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "splice.h"

#include <stdio.h>
#include "mem.h"
#include "assert.h"
#include "sense.h"
#include "genome128_hr.h"
#include "genome_sites.h"
#include "substring.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "stage3hr.h"


#define LOWPROB_SUPPORT 20

#if 0
/* Creates issues with ambiguous substrings */
#define LOCALSPLICING_NMATCHES_SLOP 1
#else
#define LOCALSPLICING_NMATCHES_SLOP 0
#endif
#define LOCALSPLICING_PROB_SLOP 0.05


/* Splice_solve_single */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Splice_solve_double */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Group by segmentj */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Trim novel splice ends */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif




static bool novelsplicingp = true;
static int min_shortend;

void
Splice_setup (int min_shortend_in) {
  min_shortend = min_shortend_in;
  return;
}



#if 0
/* Copied from stage3.c  Used also in splice.c */
static int
sufficient_splice_prob_local (int support, int nmatches, int nmismatches, double distal_spliceprob,
			      double medial_spliceprob) {
  debug3(printf("Checking for sufficient splice prob, based on %d matches, %d mismatches, and support %d\n",
		nmatches,nmismatches,support));
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
#endif

/* Do not compare against true or false */
/* Want loose criterion, otherwise, we incur slowdown from having to
   run GSNAP algorithm */
static int
sufficient_splice_prob_local (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support <= 9) {
    return (spliceprob > 0.80);
  } else if (support <= 12) {
    return (spliceprob > 0.70);
  } else if (support <= 15) {
    return (spliceprob > 0.60);
  } else if (support <= 25) {
    return (spliceprob > 0.50);
  } else {
    return (spliceprob > 0.40);
  }
}





/* Note: knowni holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

/* Called only by sarray-read.c, where plusp is always true */
int
Splice_resolve_sense (int *best_knowni_i, int *best_knowni_j,
		      int *best_nmismatches_i, int *best_nmismatches_j,
		      double *best_prob_i, double *best_prob_j,

		      Univcoord_T segmenti_left, Univcoord_T segmentj_left,
		      Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
		     
		      int querystart, int queryend, int querylength, Compress_T query_compress,
		      int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
		      int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		      int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
		      int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		      int segmenti_donor_nknown, int segmentj_acceptor_nknown,
		      int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
		      int max_mismatches_allowed, bool plusp, int genestrand) {
  int best_splice_pos = -1, splice_pos_start, splice_pos_end, splice_pos, i, j;

  int best_nmismatches, nmismatches;
  int segmenti_nmismatches, segmentj_nmismatches;
  double best_prob, probi, probj;
  bool sufficient1p, sufficient2p;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;
  int *donor_positions_alloc, *acceptor_positions_alloc, *donor_knowni_alloc, *acceptor_knowni_alloc;


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    donor_positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    acceptor_positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    donor_knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    acceptor_knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    donor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    acceptor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    donor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    acceptor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  donor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  acceptor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  donor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  acceptor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
#endif


  debug1(printf("Splice_resolve_sense: Getting genome at lefti %u and leftj %u (diff: %d), range %d..%d\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left,querystart,queryend));

  *best_knowni_i = *best_knowni_j = -1;
  *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
  *best_prob_i = *best_prob_j = 0.0;

  splice_pos_start = querystart;
  splice_pos_end = queryend;
  if (splice_pos_start < min_shortend) {
    splice_pos_start = min_shortend;
  }
  if (splice_pos_end > querylength - min_shortend) {
    splice_pos_end = querylength - min_shortend;
  }

  if (plusp == true) {
    /* Originally from plus strand.  No complement.  */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     segmenti_left,splice_pos_start,splice_pos_end);
      donori_positions = donor_positions_alloc;
      donori_knowni = donor_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   segmentj_left,splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor_positions_alloc;
      acceptorj_knowni = acceptor_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < donori_nsites && j < acceptorj_nsites) {
      if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	i++;
      } else if (splice_pos > acceptorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/querystart,/*pos3*/splice_pos,
								 plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/queryend,
								 plusp,genestrand);
	if (segmenti_nmismatches > (splice_pos - querystart)/10) {
	  /* Skip, because too many mismatches in segmenti */
	} else if (segmentj_nmismatches > (queryend - splice_pos)/10) {
	  /* Skip, because too many mismatches in segmentj */
	} else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  debug1(printf("nmismatches %d + %d <= best_nmismatches %d\n",segmenti_nmismatches,segmentj_nmismatches,best_nmismatches));
	  if (donori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	    sufficient1p = true;
	  } else {
	    probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    sufficient1p = sufficient_splice_prob_local(/*support*/splice_pos,segmenti_nmismatches,probi);
	  }

	  if (acceptorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	    sufficient2p = true;
	  } else {
	    probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    sufficient2p = sufficient_splice_prob_local(/*support*/querylength - splice_pos,segmentj_nmismatches,probj);
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		 } else {
		   printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		 });

	  if (sufficient1p && sufficient2p) {
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      /* best_donor_splicecoord = segmenti_left + splice_pos; */
	      /* best_acceptor_splicecoord = segmentj_left + splice_pos; */
	      *best_knowni_i = donori_knowni[i];
	      *best_knowni_j = acceptorj_knowni[j];
	      *best_prob_i = probi; /* donor_prob */
	      *best_prob_j = probj; /* acceptor_prob */
	      best_splice_pos = splice_pos;
	      *best_nmismatches_i = segmenti_nmismatches;
	      *best_nmismatches_j = segmentj_nmismatches;
	    }
	  }
	}
	i++;
	j++;
      }
    }

  } else {
    /* minus */
    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   segmenti_left,splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor_positions_alloc;
      antiacceptori_knowni = acceptor_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
      antidonorj_positions = donor_positions_alloc;
      antidonorj_knowni = donor_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < antiacceptori_nsites && j < antidonorj_nsites) {
      if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	i++;
      } else if (splice_pos > antidonorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/querystart,/*pos3*/splice_pos,
								 plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/queryend,
								 plusp,genestrand);
	if (segmenti_nmismatches > (splice_pos - querystart)/10) {
	  /* Skip, because too many mismatches in segmenti */
	} else if (segmentj_nmismatches > (queryend - splice_pos)/10) {
	  /* Skip, because too many mismatches in segmentj */
	} else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  debug1(printf("nmismatches %d + %d <= best_nmismatches %d\n",segmenti_nmismatches,segmentj_nmismatches,best_nmismatches));
	  if (antiacceptori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	    sufficient1p = true;
	  } else {
	    probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    sufficient1p = sufficient_splice_prob_local(/*support*/splice_pos,segmenti_nmismatches,probi);
	  }

	  if (antidonorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	    sufficient2p = true;
	  } else {
	    probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    sufficient2p = sufficient_splice_prob_local(/*support*/querylength - splice_pos,segmentj_nmismatches,probj);
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		 } else {
		   printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		 });
	  
	  if (sufficient1p && sufficient2p) {
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;
	      
	      /* best_donor_splicecoord = segmentj_left + splice_pos; */
	      /* best_acceptor_splicecoord = segmenti_left + splice_pos; */
	      *best_knowni_j = antidonorj_knowni[j];
	      *best_knowni_i = antiacceptori_knowni[i];
	      *best_prob_j = probj; /* donor_prob */
	      *best_prob_i = probi;
	      best_splice_pos = splice_pos;
	      *best_nmismatches_j = segmentj_nmismatches;
	      *best_nmismatches_i = segmenti_nmismatches;
	    }
	  }
	}
	i++;
	j++;
      }
    }
  }

  debug1(printf("best_knowni_i is %d and best_knowni_j is %d\n",*best_knowni_i,*best_knowni_j));

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(donor_positions_alloc);
    FREEA(acceptor_positions_alloc);
    FREEA(donor_knowni_alloc);
    FREEA(acceptor_knowni_alloc);
  } else {
    FREE(donor_positions_alloc);
    FREE(acceptor_positions_alloc);
    FREE(donor_knowni_alloc);
    FREE(acceptor_knowni_alloc);
  }
#else
  FREE(donor_positions_alloc);
  FREE(acceptor_positions_alloc);
  FREE(donor_knowni_alloc);
  FREE(acceptor_knowni_alloc);
#endif


  if (*best_prob_i > 0.95 && *best_prob_j > 0.70) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    debug1(printf("nmismatches %d and %d\n",*best_nmismatches_i,*best_nmismatches_j));
    return best_splice_pos;
  } else if (*best_prob_i > 0.70 && *best_prob_j > 0.95) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    debug1(printf("nmismatches %d and %d\n",*best_nmismatches_i,*best_nmismatches_j));
    return best_splice_pos;
  } else if (*best_prob_i > 0.40 && *best_prob_j > 0.40) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    debug1(printf("nmismatches %d and %d\n",*best_nmismatches_i,*best_nmismatches_j));
    return best_splice_pos;
  } else {
    debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
    return -1;
  }
}


/* Called only by sarray-read.c, where plusp is always true */
int
Splice_resolve_antisense (int *best_knowni_i, int *best_knowni_j,
			  int *best_nmismatches_i, int *best_nmismatches_j,
			  double *best_prob_i, double *best_prob_j,
			  
			  Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			  Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
			  
			  int querystart, int queryend, int querylength, Compress_T query_compress,
			  int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			  int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			  int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			  int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			  int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			  int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			  int max_mismatches_allowed, bool plusp, int genestrand) {
  int best_splice_pos = -1, splice_pos_start, splice_pos_end, splice_pos, i, j;

  int best_nmismatches, nmismatches;
  int segmenti_nmismatches, segmentj_nmismatches;
  double best_prob, probi, probj;
  bool sufficient1p, sufficient2p;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;
  int *donor_positions_alloc, *acceptor_positions_alloc, *donor_knowni_alloc, *acceptor_knowni_alloc;


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    donor_positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    acceptor_positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    donor_knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    acceptor_knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    donor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    acceptor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    donor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    acceptor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  donor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  acceptor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  donor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  acceptor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug1(printf("Splice_resolve_antisense: Getting genome at lefti %u and leftj %u (diff: %d), range %d..%d\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left,querystart,queryend));

  *best_knowni_i = *best_knowni_j = -1;
  *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
  *best_prob_i = *best_prob_j = 0.0;

  splice_pos_start = querystart;
  splice_pos_end = queryend;
  if (splice_pos_start < min_shortend) {
    splice_pos_start = min_shortend;
  }
  if (splice_pos_end > querylength - min_shortend) {
    splice_pos_end = querylength - min_shortend;
  }

  if (plusp == false) {
    /* minus */
    /* Originally from plus strand.  No complement.  */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     segmenti_left,splice_pos_start,splice_pos_end);
      donori_positions = donor_positions_alloc;
      donori_knowni = donor_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   segmentj_left,splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor_positions_alloc;
      acceptorj_knowni = acceptor_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < donori_nsites && j < acceptorj_nsites) {
      if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	i++;
      } else if (splice_pos > acceptorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/querystart,/*pos3*/splice_pos,
								 plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/queryend,
								 plusp,genestrand);
	if (segmenti_nmismatches > (splice_pos - querystart)/10) {
	  /* Skip, because too many mismatches in segmenti */
	} else if (segmentj_nmismatches > (queryend - splice_pos)/10) {
	  /* Skip, because too many mismatches in segmentj */
	} else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  if (donori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	    sufficient1p = true;
	  } else {
	    probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    sufficient1p = sufficient_splice_prob_local(/*support*/splice_pos,segmenti_nmismatches,probi);
	  }
	  
	  if (acceptorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	    sufficient2p = true;
	  } else {
	    probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    sufficient2p = sufficient_splice_prob_local(/*support*/querylength - splice_pos,segmentj_nmismatches,probj);
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		 } else {
		   printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		 });

	  if (sufficient1p && sufficient2p) {
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;
	      
	      /* best_donor_splicecoord = segmenti_left + splice_pos; */
	      /* best_acceptor_splicecoord = segmentj_left + splice_pos; */
	      *best_knowni_i = donori_knowni[i];
	      *best_knowni_j = acceptorj_knowni[j];
	      *best_prob_i = probi; /* donor_prob */
	      *best_prob_j = probj; /* acceptor_prob */
	      best_splice_pos = splice_pos;
	      *best_nmismatches_i = segmenti_nmismatches;
	      *best_nmismatches_j = segmentj_nmismatches;
	    }
	  }
	}
	i++;
	j++;
      }
    }

  } else {
    /* plus */
    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   segmenti_left,splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor_positions_alloc;
      antiacceptori_knowni = acceptor_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
      antidonorj_positions = donor_positions_alloc;
      antidonorj_knowni = donor_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < antiacceptori_nsites && j < antidonorj_nsites) {
      if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	i++;
      } else if (splice_pos > antidonorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/querystart,/*pos3*/splice_pos,
								 plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/queryend,
								 plusp,genestrand);
	if (segmenti_nmismatches > (splice_pos - querystart)/10) {
	  /* Skip, because too many mismatches in segmenti */
	} else if (segmentj_nmismatches > (queryend - splice_pos)/10) {
	  /* Skip, because too many mismatches in segmentj */
	} else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  if (antiacceptori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	    sufficient1p = true;
	  } else {
	    probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    sufficient1p = sufficient_splice_prob_local(/*support*/splice_pos,segmenti_nmismatches,probi);
	  }

	  if (antidonorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	    sufficient2p = true;
	  } else {
	    probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    sufficient2p = sufficient_splice_prob_local(/*support*/querylength - splice_pos,segmentj_nmismatches,probj);
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		 } else {
		   printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		 });
	  
	  if (sufficient1p && sufficient2p) {
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;
	      
	      /* best_donor_splicecoord = segmentj_left + splice_pos; */
	      /* best_acceptor_splicecoord = segmenti_left + splice_pos; */
	      *best_knowni_j = antidonorj_knowni[j];
	      *best_knowni_i = antiacceptori_knowni[i];
	      *best_prob_j = probj; /* donor_prob */
	      *best_prob_i = probi; /* acceptor_prob */
	      best_splice_pos = splice_pos;
	      *best_nmismatches_j = segmentj_nmismatches;
	      *best_nmismatches_i = segmenti_nmismatches;
	    }
	  }
	}
	i++;
	j++;
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(donor_positions_alloc);
    FREEA(acceptor_positions_alloc);
    FREEA(donor_knowni_alloc);
    FREEA(acceptor_knowni_alloc);
  } else {
    FREE(donor_positions_alloc);
    FREE(acceptor_positions_alloc);
    FREE(donor_knowni_alloc);
    FREE(acceptor_knowni_alloc);
  }
#else
  FREE(donor_positions_alloc);
  FREE(acceptor_positions_alloc);
  FREE(donor_knowni_alloc);
  FREE(acceptor_knowni_alloc);
#endif


  debug1(printf("best_knowni_i is %d and best_knowni_j is %d\n",*best_knowni_i,*best_knowni_j));

  if (*best_prob_i > 0.95 && *best_prob_j > 0.70) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    debug1(printf("nmismatches %d and %d\n",*best_nmismatches_i,*best_nmismatches_j));
    return best_splice_pos;
  } else if (*best_prob_i > 0.70 && *best_prob_j > 0.95) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    debug1(printf("nmismatches %d and %d\n",*best_nmismatches_i,*best_nmismatches_j));
    return best_splice_pos;
  } else if (*best_prob_i > 0.40 && *best_prob_j > 0.40) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    debug1(printf("nmismatches %d and %d\n",*best_nmismatches_i,*best_nmismatches_j));
    return best_splice_pos;
  } else {
    debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
    *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
    return -1;
  }
}



/* Note: knowni holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

List_T
Splice_solve_single_sense (int *found_score, int *nhits, List_T hits, List_T *lowprob,

			   bool *segmenti_usedp, bool *segmentj_usedp,
			   Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			   Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
			   Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
			   Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
			   Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
			   int querylength, Compress_T query_compress,
			   int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			   int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			   int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			   int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			   int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			   int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			   int splicing_penalty, int max_mismatches_allowed,
			   bool plusp, int genestrand, bool first_read_p,
			   bool subs_or_indels_p, bool sarrayp) {
  Substring_T donor, acceptor;
  int best_splice_pos, splice_pos_start, splice_pos_end, splice_pos, i, j;

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentj_nmismatches, segmenti_nmismatches, segmentj_nmismatches;
  int donor_support, acceptor_support;
  Univcoord_T best_donor_splicecoord, best_acceptor_splicecoord;
  int best_donor_knowni, best_acceptor_knowni;
  double best_prob, best_donor_prob, best_acceptor_prob, probi, probj;
  bool sufficient1p, sufficient2p, orig_plusp;
  int sensedir;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;
  int *donor_positions_alloc, *acceptor_positions_alloc, *donor_knowni_alloc, *acceptor_knowni_alloc;


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    donor_positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    acceptor_positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    donor_knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    acceptor_knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    donor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    acceptor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    donor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    acceptor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  donor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  acceptor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  donor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  acceptor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug1(printf("Splice_solve_single: Getting genome at lefti %u and leftj %u (diff: %d)\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left));
  *nhits = 0;

#if 0
  int sum, lefti, righti;
  splice_pos_start = querylength;
  splice_pos_end = 0;
  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right &&
	  mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	debug1(printf("At %d+%d mismatches, splice_pos using right: %d\n",lefti,righti,mismatch_positions_right[righti]+1));
	debug1(printf("At %d+%d mismatches, splice_pos using left: %d\n",lefti,righti,mismatch_positions_left[lefti]));
	if (mismatch_positions_right[righti] + 1 < splice_pos_start) {
	  splice_pos_start = mismatch_positions_right[righti] + 1;	/* This is leftmost position in righti+1 .. lefti */
	}
	if (mismatch_positions_left[lefti] > splice_pos_end) {
	  splice_pos_end = mismatch_positions_left[lefti];	/* This is rightmost position in righti+1 .. lefti */
	}
      }
    }
  }

  /* Exclude ends */
  if (splice_pos_start < min_localsplicing_end_matches) {
    splice_pos_start = min_localsplicing_end_matches;
  }
  if (splice_pos_end > querylength - min_localsplicing_end_matches) {
    splice_pos_end = querylength - min_localsplicing_end_matches;
  }
#else
  /* splice_pos_start = min_localsplicing_end_matches; */
  /* splice_pos_end = querylength - min_localsplicing_end_matches; */
  splice_pos_start = min_shortend;
  splice_pos_end = querylength - min_shortend; /* ? off by 1, so -l 3 allows only ends of up to 2 */
#endif


  if (splice_pos_start <= splice_pos_end) {
    if (plusp == true) {
      /* Originally from plus strand.  No complement.  */
      /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
      if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
	donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					       segmenti_donor_knownpos,segmenti_donor_knowni,
					       segmenti_left,splice_pos_start,splice_pos_end);
	donori_positions = donor_positions_alloc;
	donori_knowni = donor_knowni_alloc;
      } else {
	donori_nsites = segmenti_donor_nknown;
	donori_positions = segmenti_donor_knownpos;
	donori_knowni = segmenti_donor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d donori sites:",donori_nsites);
      for (i = 0; i < donori_nsites; i++) {
	printf(" %d",donori_positions[i]);
	if (donori_knowni[i] >= 0) {
	  printf(" (%d)",donori_knowni[i]);
	}
      }
      printf("\n");
#endif

      if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						     segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
	acceptorj_positions = acceptor_positions_alloc;
	acceptorj_knowni = acceptor_knowni_alloc;
      } else {
	acceptorj_nsites = segmentj_acceptor_nknown;
	acceptorj_positions = segmentj_acceptor_knownpos;
	acceptorj_knowni = segmentj_acceptor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d acceptorj sites:",acceptorj_nsites);
      for (i = 0; i < acceptorj_nsites; i++) {
	printf(" %d",acceptorj_positions[i]);
	if (acceptorj_knowni[i] >= 0) {
	  printf(" (%d)",acceptorj_knowni[i]);
	}
      }
      printf("\n");
#endif

      best_nmismatches = max_mismatches_allowed;
      best_prob = 0.0;

      i = j = 0;
      while (i < donori_nsites && j < acceptorj_nsites) {
	if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	  i++;
	} else if (splice_pos > acceptorj_positions[j]) {
	  j++;
	} else {
	  segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   plusp,genestrand);
	  segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   plusp,genestrand);
	  if (segmenti_nmismatches > (splice_pos - 0)/10) {
	    /* Skip, because too many mismatches in segmenti */
	  } else if (segmentj_nmismatches > (querylength - splice_pos)/10) {
	    /* Skip, because too many mismatches in segmentj */
	  } else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	    if (donori_knowni[i] >= 0) {
	      probi = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    }

	    if (acceptorj_knowni[j] >= 0) {
	      probj = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    }

	    debug1(
		   if (plusp == true) {
		     printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		   } else {
		     printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		   });

	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      best_donor_splicecoord = segmenti_left + splice_pos;
	      best_acceptor_splicecoord = segmentj_left + splice_pos;
	      best_donor_knowni = donori_knowni[i];
	      best_acceptor_knowni = acceptorj_knowni[j];
	      best_donor_prob = probi;
	      best_acceptor_prob = probj;
	      best_splice_pos = splice_pos;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      orig_plusp = true;	/* for sense, require plusp to be true */
	    }
	  }
	  i++;
	  j++;
	}
      }

    } else {
      /* minus */
      /* Originally from minus strand.  Complement. */
      /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
      if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							     segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							     segmenti_left,splice_pos_start,splice_pos_end);
	antiacceptori_positions = acceptor_positions_alloc;
	antiacceptori_knowni = acceptor_knowni_alloc;
      } else {
	antiacceptori_nsites = segmenti_antiacceptor_nknown;
	antiacceptori_positions = segmenti_antiacceptor_knownpos;
	antiacceptori_knowni = segmenti_antiacceptor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d antiacceptori sites:",antiacceptori_nsites);
      for (i = 0; i < antiacceptori_nsites; i++) {
	printf(" %d",antiacceptori_positions[i]);
	if (antiacceptori_knowni[i] >= 0) {
	  printf(" (%d)",antiacceptori_knowni[i]);
	}
      }
      printf("\n");
#endif

      if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
	antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						       segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						       segmentj_left,splice_pos_start,splice_pos_end);
	antidonorj_positions = donor_positions_alloc;
	antidonorj_knowni = donor_knowni_alloc;
      } else {
	antidonorj_nsites = segmentj_antidonor_nknown;
	antidonorj_positions = segmentj_antidonor_knownpos;
	antidonorj_knowni = segmentj_antidonor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d antidonorj sites:",antidonorj_nsites);
      for (i = 0; i < antidonorj_nsites; i++) {
	printf(" %d",antidonorj_positions[i]);
	if (antidonorj_knowni[i] >= 0) {
	  printf(" (%d)",antidonorj_knowni[i]);
	}
      }
      printf("\n");
#endif

      best_nmismatches = max_mismatches_allowed;
      best_prob = 0.0;

      i = j = 0;
      while (i < antiacceptori_nsites && j < antidonorj_nsites) {
	if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	  i++;
	} else if (splice_pos > antidonorj_positions[j]) {
	  j++;
	} else {
	  segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   plusp,genestrand);
	  segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   plusp,genestrand);
	  if (segmenti_nmismatches > (splice_pos - 0)/10) {
	    /* Skip, because too many mismatches in segmenti */
	  } else if (segmentj_nmismatches > (querylength - splice_pos)/10) {
	    /* Skip, because too many mismatches in segmentj */
	  } else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	    if (antiacceptori_knowni[i] >= 0) {
	      probi = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    }

	    if (antidonorj_knowni[j] >= 0) {
	      probj = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    }

	    debug1(
		   if (plusp == true) {
		     printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		   } else {
		     printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		   });
	  
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      best_donor_splicecoord = segmentj_left + splice_pos;
	      best_acceptor_splicecoord = segmenti_left + splice_pos;
	      best_donor_knowni = antidonorj_knowni[j];
	      best_acceptor_knowni = antiacceptori_knowni[i];
	      best_donor_prob = probj;
	      best_acceptor_prob = probi;
	      best_splice_pos = splice_pos;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      orig_plusp = false;	/* for sense, require plusp to be false */
	    }
	  }
	  i++;
	  j++;
	}
      }
    }

    if (best_prob > 0.0) {
      debug1(printf("best_prob = %f at splice_pos %d (%u,%u)\n",
		    best_prob,best_splice_pos,best_donor_splicecoord,best_acceptor_splicecoord));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */

	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;
	assert(plusp == true);
	assert(sensedir == SENSE_FORWARD);

	donor = Substring_new_donor(best_donor_splicecoord,best_donor_knowni,
				    best_splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
				    best_segmenti_nmismatches,
				    best_donor_prob,/*left*/segmenti_left,query_compress,
				    querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
				    segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_splicecoord,best_acceptor_knowni,
					  best_splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
					  best_segmentj_nmismatches,
					  best_acceptor_prob,/*left*/segmentj_left,query_compress,
					  querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
					  segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);
	  
	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  debug1(printf("Splice_solve_single_sense success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  donor_support = best_splice_pos;
	  acceptor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(donor_support,best_segmenti_nmismatches,best_donor_prob);
	  sufficient2p = sufficient_splice_prob_local(acceptor_support,best_segmentj_nmismatches,best_acceptor_prob);

	  if (sufficient1p && sufficient2p) {
	    *nhits += 1;
	    hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
								donor,acceptor,best_donor_prob,best_acceptor_prob,
								/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,
								/*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								/*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								/*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								/*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	    /* return hits; */
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* return hits; */
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* return hits; */
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
							       donor,acceptor,best_donor_prob,best_acceptor_prob,
							       /*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,
							       /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
							       /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
							       /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
							       /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    /* return hits; */
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* ? return hits; */
	  }
	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;
	assert(plusp == false);
	assert(sensedir == SENSE_FORWARD);

	donor = Substring_new_donor(best_donor_splicecoord,best_donor_knowni,
				    best_splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
				    best_segmentj_nmismatches,
				    best_donor_prob,/*left*/segmentj_left,query_compress,
				    querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
				    segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_splicecoord,best_acceptor_knowni,
					  best_splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
					  best_segmenti_nmismatches,
					  best_acceptor_prob,/*left*/segmenti_left,query_compress,
					  querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
					  segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  debug1(printf("Splice_solve_single_sense success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  acceptor_support = best_splice_pos;
	  donor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(acceptor_support,best_segmenti_nmismatches,best_acceptor_prob);
	  sufficient2p = sufficient_splice_prob_local(donor_support,best_segmentj_nmismatches,best_donor_prob);
	  if (sufficient1p && sufficient2p) {
	    *nhits += 1;
	    hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
								donor,acceptor,best_donor_prob,best_acceptor_prob,
								/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,
								/*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								/*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								/*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								/*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	    /* return hits; */
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* return hits; */
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* return hits; */
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
							       donor,acceptor,best_donor_prob,best_acceptor_prob,
							       /*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,
							       /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
							       /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
							       /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
							       /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    /* return hits; */
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* ? return hits; */
	  }
	}
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(donor_positions_alloc);
    FREEA(acceptor_positions_alloc);
    FREEA(donor_knowni_alloc);
    FREEA(acceptor_knowni_alloc);
  } else {
    FREE(donor_positions_alloc);
    FREE(acceptor_positions_alloc);
    FREE(donor_knowni_alloc);
    FREE(acceptor_knowni_alloc);
  }
#else
  FREE(donor_positions_alloc);
  FREE(acceptor_positions_alloc);
  FREE(donor_knowni_alloc);
  FREE(acceptor_knowni_alloc);
#endif

  return hits;
}


List_T
Splice_solve_single_antisense (int *found_score, int *nhits, List_T hits, List_T *lowprob,

			       bool *segmenti_usedp, bool *segmentj_usedp,
			       Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			       Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
			       Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
			       Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
			       Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
			       int querylength, Compress_T query_compress,
			       int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			       int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			       int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			       int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			       int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			       int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			       int splicing_penalty, int max_mismatches_allowed,
			       bool plusp, int genestrand, bool first_read_p,
			       bool subs_or_indels_p, bool sarrayp) {
  Substring_T donor, acceptor;
  int best_splice_pos, splice_pos_start, splice_pos_end, splice_pos, i, j;

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentj_nmismatches, segmenti_nmismatches, segmentj_nmismatches;
  int donor_support, acceptor_support;
  Univcoord_T best_donor_splicecoord, best_acceptor_splicecoord;
  int best_donor_knowni, best_acceptor_knowni;
  double best_prob, best_donor_prob, best_acceptor_prob, probi, probj;
  bool sufficient1p, sufficient2p, orig_plusp;
  int sensedir;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;
  int *donor_positions_alloc, *acceptor_positions_alloc, *donor_knowni_alloc, *acceptor_knowni_alloc;


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    donor_positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    acceptor_positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    donor_knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
    acceptor_knowni_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    donor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    acceptor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    donor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
    acceptor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  donor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  acceptor_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  donor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  acceptor_knowni_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug1(printf("Splice_solve_single: Getting genome at lefti %u and leftj %u (diff: %d)\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left));
  *nhits = 0;

#if 0
  int sum, lefti, righti;
  splice_pos_start = querylength;
  splice_pos_end = 0;
  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right &&
	  mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	debug1(printf("At %d+%d mismatches, splice_pos using right: %d\n",lefti,righti,mismatch_positions_right[righti]+1));
	debug1(printf("At %d+%d mismatches, splice_pos using left: %d\n",lefti,righti,mismatch_positions_left[lefti]));
	if (mismatch_positions_right[righti] + 1 < splice_pos_start) {
	  splice_pos_start = mismatch_positions_right[righti] + 1;	/* This is leftmost position in righti+1 .. lefti */
	}
	if (mismatch_positions_left[lefti] > splice_pos_end) {
	  splice_pos_end = mismatch_positions_left[lefti];	/* This is rightmost position in righti+1 .. lefti */
	}
      }
    }
  }

  /* Exclude ends */
  if (splice_pos_start < min_localsplicing_end_matches) {
    splice_pos_start = min_localsplicing_end_matches;
  }
  if (splice_pos_end > querylength - min_localsplicing_end_matches) {
    splice_pos_end = querylength - min_localsplicing_end_matches;
  }
#else
  /* splice_pos_start = min_localsplicing_end_matches; */
  /* splice_pos_end = querylength - min_localsplicing_end_matches; */
  splice_pos_start = min_shortend;
  splice_pos_end = querylength - min_shortend; /* ? off by 1, so -l 3 allows only ends of up to 2 */
#endif


  if (splice_pos_start <= splice_pos_end) {
    if (plusp == false) {
      /* minus */
      /* Originally from plus strand.  No complement.  */
      /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
      if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
	donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					       segmenti_donor_knownpos,segmenti_donor_knowni,
					       segmenti_left,splice_pos_start,splice_pos_end);
	donori_positions = donor_positions_alloc;
	donori_knowni = donor_knowni_alloc;
      } else {
	donori_nsites = segmenti_donor_nknown;
	donori_positions = segmenti_donor_knownpos;
	donori_knowni = segmenti_donor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d donori sites:",donori_nsites);
      for (i = 0; i < donori_nsites; i++) {
	printf(" %d",donori_positions[i]);
	if (donori_knowni[i] >= 0) {
	  printf(" (%d)",donori_knowni[i]);
	}
      }
      printf("\n");
#endif

      if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						     segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
	acceptorj_positions = acceptor_positions_alloc;
	acceptorj_knowni = acceptor_knowni_alloc;
      } else {
	acceptorj_nsites = segmentj_acceptor_nknown;
	acceptorj_positions = segmentj_acceptor_knownpos;
	acceptorj_knowni = segmentj_acceptor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d acceptorj sites:",acceptorj_nsites);
      for (i = 0; i < acceptorj_nsites; i++) {
	printf(" %d",acceptorj_positions[i]);
	if (acceptorj_knowni[i] >= 0) {
	  printf(" (%d)",acceptorj_knowni[i]);
	}
      }
      printf("\n");
#endif

      best_nmismatches = max_mismatches_allowed;
      best_prob = 0.0;

      i = j = 0;
      while (i < donori_nsites && j < acceptorj_nsites) {
	if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	  i++;
	} else if (splice_pos > acceptorj_positions[j]) {
	  j++;
	} else {
	  segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   plusp,genestrand);
	  segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   plusp,genestrand);
	  if (segmenti_nmismatches > (splice_pos - 0)/10) {
	    /* Skip, because too many mismatches in segmenti */
	  } else if (segmentj_nmismatches > (querylength - splice_pos)/10) {
	    /* Skip, because too many mismatches in segmentj */
	  } else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	    if (donori_knowni[i] >= 0) {
	      probi = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    }

	    if (acceptorj_knowni[j] >= 0) {
	      probj = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    }

	    debug1(
		   if (plusp == true) {
		     printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		   } else {
		     printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		   });

	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      best_donor_splicecoord = segmenti_left + splice_pos;
	      best_acceptor_splicecoord = segmentj_left + splice_pos;
	      best_donor_knowni = donori_knowni[i];
	      best_acceptor_knowni = acceptorj_knowni[j];
	      best_donor_prob = probi;
	      best_acceptor_prob = probj;
	      best_splice_pos = splice_pos;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      orig_plusp = true;	/* for antisense, require plusp to be false */
	    }
	  }
	  i++;
	  j++;
	}
      }

    } else {
      /* plus */
      /* Originally from minus strand.  Complement. */
      /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
      if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							     segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							     segmenti_left,splice_pos_start,splice_pos_end);
	antiacceptori_positions = acceptor_positions_alloc;
	antiacceptori_knowni = acceptor_knowni_alloc;
      } else {
	antiacceptori_nsites = segmenti_antiacceptor_nknown;
	antiacceptori_positions = segmenti_antiacceptor_knownpos;
	antiacceptori_knowni = segmenti_antiacceptor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d antiacceptori sites:",antiacceptori_nsites);
      for (i = 0; i < antiacceptori_nsites; i++) {
	printf(" %d",antiacceptori_positions[i]);
	if (antiacceptori_knowni[i] >= 0) {
	  printf(" (%d)",antiacceptori_knowni[i]);
	}
      }
      printf("\n");
#endif

      if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
	antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						       segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						       segmentj_left,splice_pos_start,splice_pos_end);
	antidonorj_positions = donor_positions_alloc;
	antidonorj_knowni = donor_knowni_alloc;
      } else {
	antidonorj_nsites = segmentj_antidonor_nknown;
	antidonorj_positions = segmentj_antidonor_knownpos;
	antidonorj_knowni = segmentj_antidonor_knowni;
      }

#ifdef DEBUG1
      printf("Found %d antidonorj sites:",antidonorj_nsites);
      for (i = 0; i < antidonorj_nsites; i++) {
	printf(" %d",antidonorj_positions[i]);
	if (antidonorj_knowni[i] >= 0) {
	  printf(" (%d)",antidonorj_knowni[i]);
	}
      }
      printf("\n");
#endif

      best_nmismatches = max_mismatches_allowed;
      best_prob = 0.0;

      i = j = 0;
      while (i < antiacceptori_nsites && j < antidonorj_nsites) {
	if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	  i++;
	} else if (splice_pos > antidonorj_positions[j]) {
	  j++;
	} else {
	  segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   plusp,genestrand);
	  segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   plusp,genestrand);
	  if (segmenti_nmismatches > (splice_pos - 0)/10) {
	    /* Skip, because too many mismatches in segmenti */
	  } else if (segmentj_nmismatches > (querylength - splice_pos)/10) {
	    /* Skip, because too many mismatches in segmentj */
	  } else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	    if (antiacceptori_knowni[i] >= 0) {
	      probi = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	    }

	    if (antidonorj_knowni[j] >= 0) {
	      probj = 1.0; /* Needs to be 1.0 for output */
	    } else {
	      probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	    }

	    debug1(
		   if (plusp == true) {
		     printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		   } else {
		     printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		   });
	  
	    if (nmismatches < best_nmismatches ||
		(nmismatches == best_nmismatches && probi + probj > best_prob)) {
	      /* Success */
	      best_nmismatches = nmismatches;
	      best_prob = probi + probj;

	      best_donor_splicecoord = segmentj_left + splice_pos;
	      best_acceptor_splicecoord = segmenti_left + splice_pos;
	      best_donor_knowni = antidonorj_knowni[j];
	      best_acceptor_knowni = antiacceptori_knowni[i];
	      best_donor_prob = probj;
	      best_acceptor_prob = probi;
	      best_splice_pos = splice_pos;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      orig_plusp = false;	/* for antisense, require plusp to be true */
	    }
	  }
	  i++;
	  j++;
	}
      }
    }

    if (best_prob > 0.0) {
      debug1(printf("best_prob = %f at splice_pos %d (%u,%u)\n",
		    best_prob,best_splice_pos,best_donor_splicecoord,best_acceptor_splicecoord));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;
	assert(plusp == false);
	assert(sensedir == SENSE_ANTI);

	donor = Substring_new_donor(best_donor_splicecoord,best_donor_knowni,
				    best_splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
				    best_segmenti_nmismatches,
				    best_donor_prob,/*left*/segmenti_left,query_compress,
				    querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
				    segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_splicecoord,best_acceptor_knowni,
					  best_splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
					  best_segmentj_nmismatches,
					  best_acceptor_prob,/*left*/segmentj_left,query_compress,
					  querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
					  segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  debug1(printf("Splice_solve_single_antisense success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  donor_support = best_splice_pos;
	  acceptor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(donor_support,best_segmenti_nmismatches,best_donor_prob);
	  sufficient2p = sufficient_splice_prob_local(acceptor_support,best_segmentj_nmismatches,best_acceptor_prob);

	  if (sufficient1p && sufficient2p) {
	    *nhits += 1;
	    hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
								donor,acceptor,best_donor_prob,best_acceptor_prob,
								/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,
								/*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								/*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								/*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								/*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	    /* return hits; */
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* return hits; */
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* return hits; */
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
							       donor,acceptor,best_donor_prob,best_acceptor_prob,
							       /*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,
							       /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
							       /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
							       /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
							       /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    /* return hits; */
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* ? return hits; */
	  }
	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;
	assert(plusp == true);
	assert(sensedir == SENSE_ANTI);

	donor = Substring_new_donor(best_donor_splicecoord,best_donor_knowni,
				    best_splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
				    best_segmentj_nmismatches,
				    best_donor_prob,/*left*/segmentj_left,query_compress,
				    querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
				    segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_splicecoord,best_acceptor_knowni,
					  best_splice_pos,/*substring_querystart*/0,/*substring_queryend*/querylength,
					  best_segmenti_nmismatches,
					  best_acceptor_prob,/*left*/segmenti_left,query_compress,
					  querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
					  segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  debug1(printf("Splice_solve_single_antisense success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  acceptor_support = best_splice_pos;
	  donor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(acceptor_support,best_segmenti_nmismatches,best_acceptor_prob);
	  sufficient2p = sufficient_splice_prob_local(donor_support,best_segmentj_nmismatches,best_donor_prob);
	  if (sufficient1p && sufficient2p) {
	    *nhits += 1;
	    hits = List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
								donor,acceptor,best_donor_prob,best_acceptor_prob,
								/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,
								/*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
								/*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
								/*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
								/*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	    /* return hits; */
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* return hits; */
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* return hits; */
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
							       donor,acceptor,best_donor_prob,best_acceptor_prob,
							       /*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,
							       /*ambcoords_donor*/NULL,/*ambcoords_acceptor*/NULL,
							       /*amb_knowni_donor*/NULL,/*amb_knowni_acceptor*/NULL,
							       /*amb_nmismatches_donor*/NULL,/*amb_nmismatches_acceptor*/NULL,
							       /*amb_probs_donor*/NULL,/*amb_probs_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    /* return hits; */
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    /* ? return hits; */
	  }
	}
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(donor_positions_alloc);
    FREEA(acceptor_positions_alloc);
    FREEA(donor_knowni_alloc);
    FREEA(acceptor_knowni_alloc);
  } else {
    FREE(donor_positions_alloc);
    FREE(acceptor_positions_alloc);
    FREE(donor_knowni_alloc);
    FREE(acceptor_knowni_alloc);
  }
#else
  FREE(donor_positions_alloc);
  FREE(acceptor_positions_alloc);
  FREE(donor_knowni_alloc);
  FREE(acceptor_knowni_alloc);
#endif

  return hits;
}



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
group_by_segmenti_aux (int *found_score, List_T winners, List_T *ambiguous,
		       Stage3end_T *hitarray, int n, int querylength, bool first_read_p, bool sarrayp) {
  Stage3end_T hit, *subarray;
  int i, j, k, ii, jj, kk, nn;
  int n_good_spliceends;
  Univcoord_T segmenti_left;
  Substring_T donor, acceptor;
  int best_nmismatches, nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob, donor_prob, acceptor_prob;
  List_T accepted_hits, rejected_hits, donor_hits, acceptor_hits, p;

  int sensedir;
#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Intlist_T amb_knowni, amb_nmismatches;
  Doublelist_T amb_probs;
  int donor_length, acceptor_length;


  i = 0;
  while (i < n) {
    hit = hitarray[i];
    segmenti_left = Stage3end_chimera_segmenti_left(hit);
    j = i + 1;
    while (j < n && Stage3end_chimera_segmenti_left(hitarray[j]) == segmenti_left) {
      j++;
    }
    if (j == i + 1) {
      /* Singleton */
      debug7(printf("Saving hit %d\n",i));
      winners = List_push(winners,(void *) hit);

    } else {
      best_nmismatches = querylength;
      best_prob = 0.0;
      for (k = i; k < j; k++) {
	hit = hitarray[k];
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
      for (k = i; k < j; k++) {
	hit = hitarray[k];
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
	for (k = i; k < j; k++) {
	  hit = hitarray[k];
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

      if (n_good_spliceends == 1) {
	winners = List_push(winners,List_head(accepted_hits));
	List_free(&accepted_hits);

      } else {
	/* Multiple hits */
	donor_hits = acceptor_hits = (List_T) NULL;
	for (p = accepted_hits; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  donor = Stage3end_substring_donor(hit);
	  acceptor = Stage3end_substring_acceptor(hit);
	  if (Stage3end_plusp(hit) == true) {
	    if (Substring_genomicstart(donor) == segmenti_left) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicstart(acceptor) == segmenti_left) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  } else {
	    if (Substring_genomicend(donor) == segmenti_left) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicend(acceptor) == segmenti_left) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      Stage3end_free(&hit);
	    }
	  }
	}

	if (donor_hits != NULL) {
	  subarray = (Stage3end_T *) List_to_array_n(&nn,donor_hits);
	  qsort(subarray,nn,sizeof(Stage3end_T),donor_match_length_cmp);
	  ii = 0;
	  while (ii < nn) {
	    hit = subarray[ii];
	    donor = Stage3end_substring_donor(hit);
	    donor_length = Substring_match_length_orig(donor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length_orig(Stage3end_substring_donor(subarray[jj])) == donor_length) {
	      jj++;
	    }
	    if (jj == ii + 1) {
	      winners = List_push(winners,(void *) hit);
	    } else {
	      sensedir = Stage3end_sensedir(hit);

	      ambcoords = NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (kk = ii; kk < jj; kk++) {
		acceptor = Stage3end_substring_acceptor(subarray[kk]);
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
								   sensedir,sarrayp));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	      for (kk = ii; kk < jj; kk++) {
		hit = subarray[kk];
		Stage3end_free(&hit);
	      }
	    }

	    ii = jj;
	  }
	  FREE(subarray);
	  List_free(&donor_hits);
	}

	if (acceptor_hits != NULL) {
	  subarray = (Stage3end_T *) List_to_array_n(&nn,acceptor_hits);
	  qsort(subarray,nn,sizeof(Stage3end_T),acceptor_match_length_cmp);
	  ii = 0;
	  while (ii < nn) {
	    hit = subarray[ii];
	    acceptor = Stage3end_substring_acceptor(hit);
	    acceptor_length = Substring_match_length_orig(acceptor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length_orig(Stage3end_substring_acceptor(subarray[jj])) == acceptor_length) {
	      jj++;
	    }
	    if (jj == ii + 1) {
	      winners = List_push(winners,(void *) hit);
	    } else {
	      sensedir = Stage3end_sensedir(hit);

	      ambcoords = NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (kk = ii; kk < jj; kk++) {
		donor = Stage3end_substring_donor(subarray[kk]);
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
								   sensedir,sarrayp));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	      for (kk = ii; kk < jj; kk++) {
		hit = subarray[kk];
		Stage3end_free(&hit);
	      }
	    }

	    ii = jj;
	  }
	  FREE(subarray);
	  List_free(&acceptor_hits);
	}

	List_free(&accepted_hits);
      }
    }

    i = j;
  }

  return winners;
}

List_T
Splice_group_by_segmenti (int *found_score, List_T localsplicing, List_T *ambiguous,
			  int querylength, bool first_read_p, bool sarrayp) {
  List_T winners = NULL, p;
  Stage3end_T *array_forward, *array_anti, hit;
  int n_sense_forward = 0, n_sense_anti = 0, k_forward, k_anti;

  for (p = localsplicing; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_sensedir(hit) == SENSE_FORWARD) {
      n_sense_forward++;
    } else {
      assert(Stage3end_sensedir(hit) == SENSE_ANTI);
      n_sense_anti++;
    }
  }

  if (n_sense_forward > 0) {
    array_forward = (Stage3end_T *) MALLOCA(n_sense_forward * sizeof(Stage3end_T));
    k_forward = 0;
  }
  if (n_sense_anti > 0) {
    array_anti = (Stage3end_T *) MALLOCA(n_sense_anti * sizeof(Stage3end_T));
    k_anti = 0;
  }

  for (p = localsplicing; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_sensedir(hit) == SENSE_FORWARD) {
      array_forward[k_forward++] = (Stage3end_T) List_head(p);
    } else {
      array_anti[k_anti++] = (Stage3end_T) List_head(p);
    }
  }

  if (n_sense_forward > 0) {
    qsort(array_forward,n_sense_forward,sizeof(Stage3end_T),Stage3end_chimera_segmenti_cmp);
    winners = group_by_segmenti_aux(&(*found_score),winners,&(*ambiguous),array_forward,n_sense_forward,
				    querylength,first_read_p,sarrayp);
    FREEA(array_forward);
  }

  if (n_sense_anti > 0) {
    qsort(array_anti,n_sense_anti,sizeof(Stage3end_T),Stage3end_chimera_segmenti_cmp);
    winners = group_by_segmenti_aux(&(*found_score),winners,&(*ambiguous),array_anti,n_sense_anti,
				    querylength,first_read_p,sarrayp);
    FREEA(array_anti);
  }

  List_free(&localsplicing);

  return winners;
}




static List_T
group_by_segmentj_aux (int *found_score, List_T winners, List_T *ambiguous, 
		       Stage3end_T *hitarray, int n, int querylength, bool first_read_p, bool sarrayp) {
  Stage3end_T hit, *subarray;
  int i, j, k, ii, jj, kk, nn;
  int n_good_spliceends;
  Univcoord_T segmentj_left;
  Substring_T donor, acceptor;
  int best_nmismatches, nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob, donor_prob, acceptor_prob;
  List_T accepted_hits, rejected_hits, donor_hits, acceptor_hits, p;
  int donor_length, acceptor_length;

  int sensedir;
#ifdef LARGE_GENOMES
  Uint8list_T ambcoords;
#else
  Uintlist_T ambcoords;
#endif
  Intlist_T amb_knowni, amb_nmismatches;
  Doublelist_T amb_probs;

  i = 0;
  while (i < n) {
    hit = hitarray[i];
    segmentj_left = Stage3end_chimera_segmentj_left(hit);
    j = i + 1;
    while (j < n && Stage3end_chimera_segmentj_left(hitarray[j]) == segmentj_left) {
      j++;
    }
    if (j == i + 1) {
      /* Singleton */
      debug7(printf("Saving hit %d\n",i));
      winners = List_push(winners,(void *) hit);

    } else {
      best_nmismatches = querylength;
      best_prob = 0.0;
      for (k = i; k < j; k++) {
	hit = hitarray[k];
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
      for (k = i; k < j; k++) {
	hit = hitarray[k];
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
	for (k = i; k < j; k++) {
	  hit = hitarray[k];
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

      if (n_good_spliceends == 1) {
	assert(List_length(accepted_hits) == 1);
	winners = List_push(winners,List_head(accepted_hits));
	List_free(&accepted_hits);

      } else {
	/* Multiple hits */
	donor_hits = acceptor_hits = (List_T) NULL;
	for (p = accepted_hits; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  donor = Stage3end_substring_donor(hit);
	  acceptor = Stage3end_substring_acceptor(hit);
	  if (Stage3end_plusp(hit) == true) {
	    if (Substring_genomicstart(donor) == segmentj_left) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicstart(acceptor) == segmentj_left) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      abort();
	      Stage3end_free(&hit);
	    }
	  } else {
	    if (Substring_genomicend(donor) == segmentj_left) {
	      donor_hits = List_push(donor_hits,(void *) hit);
	    } else if (Substring_genomicend(acceptor) == segmentj_left) {
	      acceptor_hits = List_push(acceptor_hits,(void *) hit);
	    } else {
	      abort();
	      Stage3end_free(&hit);
	    }
	  }
	}

	if (donor_hits != NULL) {
	  subarray = (Stage3end_T *) List_to_array_n(&nn,donor_hits);
	  qsort(subarray,nn,sizeof(Stage3end_T),donor_match_length_cmp);
	  ii = 0;
	  while (ii < nn) {
	    hit = subarray[ii];
	    donor = Stage3end_substring_donor(hit);
	    donor_length = Substring_match_length_orig(donor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length_orig(Stage3end_substring_donor(subarray[jj])) == donor_length) {
	      jj++;
	    }
	    if (jj == ii + 1) {
	      winners = List_push(winners,(void *) hit);
	    } else {
	      sensedir = Stage3end_sensedir(hit);

	      ambcoords = NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (kk = ii; kk < jj; kk++) {
		acceptor = Stage3end_substring_acceptor(subarray[kk]);
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
								   sensedir,sarrayp));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	      for (kk = ii; kk < jj; kk++) {
		hit = subarray[kk];
		Stage3end_free(&hit);
	      }
	    }

	    ii = jj;
	  }
	  FREE(subarray);
	  List_free(&donor_hits);
	}

	if (acceptor_hits != NULL) {
	  subarray = (Stage3end_T *) List_to_array_n(&nn,acceptor_hits);
	  qsort(subarray,nn,sizeof(Stage3end_T),acceptor_match_length_cmp);
	  ii = 0;
	  while (ii < nn) {
	    hit = subarray[ii];
	    acceptor = Stage3end_substring_acceptor(hit);
	    acceptor_length = Substring_match_length_orig(acceptor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length_orig(Stage3end_substring_acceptor(subarray[jj])) == acceptor_length) {
	      jj++;
	    }
	    if (jj == ii + 1) {
	      winners = List_push(winners,(void *) hit);
	    } else {
	      sensedir = Stage3end_sensedir(hit);

	      ambcoords = NULL;
	      amb_knowni = (Intlist_T) NULL;
	      amb_nmismatches = (Intlist_T) NULL;
	      amb_probs = (Doublelist_T) NULL;

	      for (kk = ii; kk < jj; kk++) {
		donor = Stage3end_substring_donor(subarray[kk]);
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
								   sensedir,sarrayp));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
#ifdef LARGE_GENOMES
	      Uint8list_free(&ambcoords);
#else
	      Uintlist_free(&ambcoords);
#endif
	      for (kk = ii; kk < jj; kk++) {
		hit = subarray[kk];
		Stage3end_free(&hit);
	      }
	    }

	    ii = jj;
	  }
	  FREE(subarray);
	  List_free(&acceptor_hits);
	}

	List_free(&accepted_hits);
      }
    }

    i = j;
  }

  return winners;
}

List_T
Splice_group_by_segmentj (int *found_score, List_T localsplicing, List_T *ambiguous,
			  int querylength, bool first_read_p, bool sarrayp) {
  List_T winners = NULL, p;
  Stage3end_T *array_forward, *array_anti, hit;
  int n_sense_forward = 0, n_sense_anti = 0, k_forward, k_anti;

  for (p = localsplicing; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_sensedir(hit) == SENSE_FORWARD) {
      n_sense_forward++;
    } else {
      assert(Stage3end_sensedir(hit) == SENSE_ANTI);
      n_sense_anti++;
    }
  }

  if (n_sense_forward > 0) {
    array_forward = (Stage3end_T *) MALLOCA(n_sense_forward * sizeof(Stage3end_T));
    k_forward = 0;
  }
  if (n_sense_anti > 0) {
    array_anti = (Stage3end_T *) MALLOCA(n_sense_anti * sizeof(Stage3end_T));
    k_anti = 0;
  }

  for (p = localsplicing; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_sensedir(hit) == SENSE_FORWARD) {
      array_forward[k_forward++] = (Stage3end_T) List_head(p);
    } else {
      array_anti[k_anti++] = (Stage3end_T) List_head(p);
    }
  }

  if (n_sense_forward > 0) {
    qsort(array_forward,n_sense_forward,sizeof(Stage3end_T),Stage3end_chimera_segmentj_cmp);
    winners = group_by_segmentj_aux(&(*found_score),winners,&(*ambiguous),array_forward,n_sense_forward,
				    querylength,first_read_p,sarrayp);
    FREEA(array_forward);
  }

  if (n_sense_anti > 0) {
    qsort(array_anti,n_sense_anti,sizeof(Stage3end_T),Stage3end_chimera_segmentj_cmp);
    winners = group_by_segmentj_aux(&(*found_score),winners,&(*ambiguous),array_anti,n_sense_anti,
				    querylength,first_read_p,sarrayp);
    FREEA(array_anti);
  }

  List_free(&localsplicing);

  return winners;
}



#define END_SPLICESITE_PROB_MATCH 0.90
#define END_SPLICESITE_PROB_MISMATCH 0.95

/* Derived from substring_trim_novel_spliceends in substring.c, which
   was modified from trim_novel_spliceends in stage3.c */
/* Note: If substring does not extend to ends of query, then region
   beyond querystart and queryend might actually be matching, and not
   mismatches.  Could fix in the future. */
int
Splice_trim_novel_spliceends (int *ambig_end_length_5, int *ambig_end_length_3,
			      Splicetype_T *ambig_splicetype_5, Splicetype_T *ambig_splicetype_3,
			      double *ambig_prob_5, double *ambig_prob_3, int orig_sensedir,
			      Univcoord_T start5, Univcoord_T middle5, Univcoord_T end5, bool solve5p,
			      Univcoord_T start3, Univcoord_T middle3, Univcoord_T end3, bool solve3p,
			      Univcoord_T genomicstart5, Univcoord_T genomicend3,
			      Univcoord_T chroffset, bool plusp) {

  int new_sensedir;
  Univcoord_T genomicpos, start_genomicpos, middle_genomicpos, end_genomicpos;
  Univcoord_T splice_genomepos_5, splice_genomepos_3, splice_genomepos_5_mm, splice_genomepos_3_mm;
  double donor_prob, acceptor_prob;
  double max_prob_5 = 0.0, max_prob_3 = 0.0,
    max_prob_sense_forward_5 = 0.0, max_prob_sense_anti_5 = 0.0,
    max_prob_sense_forward_3 = 0.0, max_prob_sense_anti_3 = 0.0;
  double max_prob_5_mm = 0.0, max_prob_3_mm = 0.0,
    max_prob_sense_forward_5_mm = 0.0, max_prob_sense_anti_5_mm = 0.0,
    max_prob_sense_forward_3_mm = 0.0, max_prob_sense_anti_3_mm = 0.0;
  Splicetype_T splicetype5, splicetype3, splicetype5_mm, splicetype3_mm;
  int splice_sensedir_5, splice_sensedir_3, splice_sensedir_5_mm, splice_sensedir_3_mm;


  debug13(printf("\nEntered Splice_trim_novel_spliceends with orig_sensedir %d\n",orig_sensedir));
  *ambig_end_length_5 = 0;
  *ambig_end_length_3 = 0;
  *ambig_prob_5 = 0.0;
  *ambig_prob_3 = 0.0;

#if 0
  /* Responsibility of caller */
  /* start is distal, end is medial */
  if (solve3p == false) {
    /* Skip 3' end*/
  } else if (plusp == true) {
    middle = substringN->alignend_trim + 1;
    if ((start = middle + END_SPLICESITE_SEARCH) > substringN->genomicend) {
      start = substringN->genomicend;
    }
    if ((end = middle - END_SPLICESITE_SEARCH) < substringN->alignstart_trim + MIN_EXON_LENGTH) {
      end = substringN->alignstart_trim + MIN_EXON_LENGTH;
    }
    debug13(printf("\n1 Set end points for 3' trim to be %u..%u..%u\n",start,middle,end));

  } else {
    middle = substringN->alignend_trim - 1;
    if ((start = middle - END_SPLICESITE_SEARCH) < substringN->genomicend) {
      start = substringN->genomicend;
    }
    if ((end = middle + END_SPLICESITE_SEARCH) > substringN->alignstart_trim - MIN_EXON_LENGTH) {
      end = substringN->alignstart_trim - MIN_EXON_LENGTH;
    }
    debug13(printf("\n2 Set end points for 3' trim to be %u..%u..%u\n",start,middle,end));
  }
#endif

  new_sensedir = SENSE_NULL;

  if (solve3p == false) {
    /* Skip 3' end */
  } else if (orig_sensedir == SENSE_FORWARD) {
    if (plusp) {
      splicetype3 = splicetype3_mm = DONOR;
      
      start_genomicpos = start3;
      middle_genomicpos = middle3;
      end_genomicpos = end3;
      
      /* assert(start_genomicpos >= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos >= middle_genomicpos) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	debug13(printf("3', watson, sense anti %u %u %f mm\n",genomicpos,genomicpos-chroffset,donor_prob));
	if (donor_prob > max_prob_3_mm) {
	  max_prob_3_mm = donor_prob;
	  splice_genomepos_3_mm = genomicpos;
	}
	genomicpos--;
      }
      while (genomicpos >= end_genomicpos) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	debug13(printf("3', watson, sense anti %u %u %f\n",genomicpos,genomicpos-chroffset,donor_prob));
	if (donor_prob > max_prob_3) {
	  max_prob_3 = donor_prob;
	  splice_genomepos_3 = genomicpos;
	}
	genomicpos--;
      }
      debug13(printf("\n"));

    } else {
      splicetype3 = splicetype3_mm = ANTIDONOR;

      start_genomicpos = start3;
      middle_genomicpos = middle3;
      end_genomicpos = end3;

      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos <= middle_genomicpos) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 3 */
	debug13(printf("3', crick, sense forward %u %u %f mm\n",genomicpos,genomicpos-chroffset,donor_prob));
	if (donor_prob > max_prob_3_mm) {
	  max_prob_3_mm = donor_prob;
	  splice_genomepos_3_mm = genomicpos;
	}
	genomicpos++;
      }
      while (genomicpos <= end_genomicpos) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 3 */
	debug13(printf("3', crick, sense forward %u %u %f\n",genomicpos,genomicpos-chroffset,donor_prob));
	if (donor_prob > max_prob_3) {
	  max_prob_3 = donor_prob;
	  splice_genomepos_3 = genomicpos;
	}
	genomicpos++;
      }
      debug13(printf("\n"));
    }

  } else if (orig_sensedir == SENSE_ANTI) {
    if (plusp) {
      splicetype3 = splicetype3_mm = ANTIACCEPTOR;

      start_genomicpos = start3;
      middle_genomicpos = middle3;
      end_genomicpos = end3;

      /* assert(start_genomicpos >= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos >= middle_genomicpos) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	debug13(printf("3', watson, sense forward %u %u %f mm\n",genomicpos,genomicpos-chroffset,acceptor_prob));
	if (acceptor_prob > max_prob_3_mm) {
	  max_prob_3_mm = acceptor_prob;
	  splice_genomepos_3_mm = genomicpos;
	}
	genomicpos--;
      }
      while (genomicpos >= end_genomicpos) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	debug13(printf("3', watson, sense forward %u %u %f\n",genomicpos,genomicpos-chroffset,acceptor_prob));
	if (acceptor_prob > max_prob_3) {
	  max_prob_3 = acceptor_prob;
	  splice_genomepos_3 = genomicpos;
	}
	genomicpos--;
      }
      debug13(printf("\n"));

    } else {
      splicetype3 = splicetype3_mm = ACCEPTOR;

      start_genomicpos = start3;
      middle_genomicpos = middle3;
      end_genomicpos = end3;

      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos <= middle_genomicpos) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 7 */
	debug13(printf("3', crick, sense anti %u %u %f mm\n",genomicpos,genomicpos-chroffset,acceptor_prob));
	if (acceptor_prob > max_prob_3_mm) {
	  max_prob_3_mm = acceptor_prob;
	  splice_genomepos_3_mm = genomicpos;
	}
	genomicpos++;
      }
      while (genomicpos <= end_genomicpos) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 7 */
	debug13(printf("3', crick, sense anti %u %u %f\n",genomicpos,genomicpos-chroffset,acceptor_prob));
	if (acceptor_prob > max_prob_3) {
	  max_prob_3 = acceptor_prob;
	  splice_genomepos_3 = genomicpos;
	}
	genomicpos++;
      }
      debug13(printf("\n"));
    }
      
  } else {
    if (plusp) {
      start_genomicpos = start3;
      middle_genomicpos = middle3;
      end_genomicpos = end3;

      /* assert(start_genomicpos >= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos >= middle_genomicpos) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	debug13(printf("3', watson, sense null %u %u %f %f mm\n",genomicpos,genomicpos-chroffset,donor_prob,acceptor_prob));
	if (donor_prob > max_prob_sense_forward_3_mm) {
	  max_prob_sense_forward_3_mm = donor_prob;
	  if (donor_prob > max_prob_3_mm) {
	    max_prob_3_mm = donor_prob;
	    splice_genomepos_3_mm = genomicpos;
	    /* splice_cdna_direction_3_mm = +1; */
	    splice_sensedir_3_mm = SENSE_FORWARD;
	    splicetype3_mm = DONOR;
	  }
	}
	if (acceptor_prob > max_prob_sense_anti_3_mm) {
	  max_prob_sense_anti_3_mm = acceptor_prob;
	  if (acceptor_prob > max_prob_3_mm) {
	    max_prob_3_mm = acceptor_prob;
	    splice_genomepos_3_mm = genomicpos;
	    /* splice_cdna_direction_3_mm = -1; */
	    splice_sensedir_3_mm = SENSE_ANTI;
	    splicetype3_mm = ANTIACCEPTOR;
	  }
	}
	genomicpos--;
      }
      while (genomicpos >= end_genomicpos) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	debug13(printf("3', watson, sense null %u %u %f %f\n",genomicpos,genomicpos-chroffset,donor_prob,acceptor_prob));
	if (donor_prob > max_prob_sense_forward_3) {
	  max_prob_sense_forward_3 = donor_prob;
	  if (donor_prob > max_prob_3) {
	    max_prob_3 = donor_prob;
	    splice_genomepos_3 = genomicpos;
	    /* splice_cdna_direction_3 = +1; */
	    splice_sensedir_3 = SENSE_FORWARD;
	    splicetype3 = DONOR;
	  }
	}
	if (acceptor_prob > max_prob_sense_anti_3) {
	  max_prob_sense_anti_3 = acceptor_prob;
	  if (acceptor_prob > max_prob_3) {
	    max_prob_3 = acceptor_prob;
	    splice_genomepos_3 = genomicpos;
	    /* splice_cdna_direction_3 = -1; */
	    splice_sensedir_3 = SENSE_ANTI;
	    splicetype3 = ANTIACCEPTOR;
	  }
	}
	genomicpos--;
      }
      debug13(printf("\n"));

    } else {
      start_genomicpos = start3;
      middle_genomicpos = middle3;
      end_genomicpos = end3;

      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos <= middle_genomicpos) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 3 */
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 7 */
	debug13(printf("3', crick, sense null %u %u %f %f mm\n",genomicpos,genomicpos-chroffset,donor_prob,acceptor_prob));
	if (donor_prob > max_prob_sense_forward_3_mm) {
	  max_prob_sense_forward_3_mm = donor_prob;
	  if (donor_prob > max_prob_3_mm) {
	    max_prob_3_mm = donor_prob;
	    splice_genomepos_3_mm = genomicpos;
	    /* splice_cdna_direction_3_mm = +1; */
	    splice_sensedir_3_mm = SENSE_FORWARD;
	    splicetype3_mm = ANTIDONOR;
	  }
	}
	if (acceptor_prob > max_prob_sense_anti_3_mm) {
	  max_prob_sense_anti_3_mm = acceptor_prob;
	  if (acceptor_prob > max_prob_3_mm) {
	    max_prob_3_mm = acceptor_prob;
	    splice_genomepos_3_mm = genomicpos;
	    /* splice_cdna_direction_3_mm = -1; */
	    splice_sensedir_3_mm = SENSE_ANTI;
	    splicetype3_mm = ACCEPTOR;
	  }
	}
	genomicpos++;
      }
      while (genomicpos <= end_genomicpos) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 3 */
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 7 */
	debug13(printf("3', crick, sense null %u %u %f %f\n",genomicpos,genomicpos-chroffset,donor_prob,acceptor_prob));
	if (donor_prob > max_prob_sense_forward_3) {
	  max_prob_sense_forward_3 = donor_prob;
	  if (donor_prob > max_prob_3) {
	    max_prob_3 = donor_prob;
	    splice_genomepos_3 = genomicpos;
	    /* splice_cdna_direction_3 = +1; */
	    splice_sensedir_3 = SENSE_FORWARD;
	    splicetype3 = ANTIDONOR;
	  }
	}
	if (acceptor_prob > max_prob_sense_anti_3) {
	  max_prob_sense_anti_3 = acceptor_prob;
	  if (acceptor_prob > max_prob_3) {
	    max_prob_3 = acceptor_prob;
	    splice_genomepos_3 = genomicpos;
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

  if (solve3p == false) {
    /* Skip 3' end */
  } else if (orig_sensedir != SENSE_NULL) {
    if (max_prob_3 > END_SPLICESITE_PROB_MATCH) {
      debug13(printf("Found good splice %s on 3' end at %u with probability %f\n",
		     Splicetype_string(splicetype3),splice_genomepos_3-chroffset,max_prob_3));
      if (plusp) {
	*ambig_end_length_3 = genomicend3 - splice_genomepos_3;
	debug13(printf("Set ambig_end_length_3 to be %d = %u - %u\n",*ambig_end_length_3,genomicend3,splice_genomepos_3));
      } else {
	*ambig_end_length_3 = splice_genomepos_3 - genomicend3;
	debug13(printf("Set ambig_end_length_3 to be %d = %u - %u\n",*ambig_end_length_3,splice_genomepos_3,genomicend3));
      }
      *ambig_splicetype_3 = splicetype3;
      *ambig_prob_3 = max_prob_3;

    } else if (max_prob_3_mm > END_SPLICESITE_PROB_MISMATCH) {
      debug13(printf("Found good mismatch splice %s on 3' end at %u with probability %f\n",
		     Splicetype_string(splicetype3_mm),splice_genomepos_3_mm-chroffset,max_prob_3_mm));
      if (plusp) {
	*ambig_end_length_3 = genomicend3 - splice_genomepos_3_mm;
	debug13(printf("Set ambig_end_length_3 to be %d = %u - %u\n",*ambig_end_length_3,genomicend3,splice_genomepos_3_mm));
      } else {
	*ambig_end_length_3 = splice_genomepos_3_mm - genomicend3;
	debug13(printf("Set ambig_end_length_3 to be %d = %u - %u\n",*ambig_end_length_3,splice_genomepos_3_mm,genomicend3));
      }
      *ambig_splicetype_3 = splicetype3_mm;
      *ambig_prob_3 = max_prob_3_mm;
    }
  }


#if 0
  /* Responsibility of caller */
  /* start is distal, end is medial */
  if (solve5p == false) {
    /* Skip 5' end */
  } else if (plusp == true) {
    middle = substring1->alignstart_trim - 1;
    if ((start = middle - END_SPLICESITE_SEARCH) < substring1->genomicstart) {
      start = substring1->genomicstart;
    }
    if ((end = middle + END_SPLICESITE_SEARCH) > substring1->alignend_trim - MIN_EXON_LENGTH) {
      end = substring1->alignend_trim - MIN_EXON_LENGTH;
    }
    debug13(printf("\n1 Set end points for 5' trim to be %u..%u..%u\n",start,middle,end));

  } else {
    middle = substring1->alignstart_trim + 1;
    if ((start = middle + END_SPLICESITE_SEARCH) > substring1->genomicstart) {
      start = substring1->genomicstart;
    }
    if ((end = middle - END_SPLICESITE_SEARCH) < substring1->alignend_trim + MIN_EXON_LENGTH) {
      end = substring1->alignend_trim + MIN_EXON_LENGTH;
    }
    debug13(printf("\n2 Set end points for 5' trim to be %u..%u..%u\n",start,middle,end));
  }
#endif

  if (solve5p == false) {
    /* Skip 5' end */
  } else if (orig_sensedir == SENSE_FORWARD) {
    if (plusp) {
      splicetype5 = splicetype5_mm = ACCEPTOR;

      start_genomicpos = start5;
      middle_genomicpos = middle5;
      end_genomicpos = end5;

      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos <= middle_genomicpos) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	debug13(printf("5', watson, sense forward %u %u %f mm\n",genomicpos,genomicpos-chroffset,acceptor_prob));
	if (acceptor_prob > max_prob_5_mm) {
	  max_prob_5_mm = acceptor_prob;
	  splice_genomepos_5_mm = genomicpos;
	}
	genomicpos++;
      }
      while (genomicpos <= end_genomicpos) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	debug13(printf("5', watson, sense forward %u %u %f\n",genomicpos,genomicpos-chroffset,acceptor_prob));
	if (acceptor_prob > max_prob_5) {
	  max_prob_5 = acceptor_prob;
	  splice_genomepos_5 = genomicpos;
	}
	genomicpos++;
      }
      debug13(printf("\n"));

    } else {
      splicetype5 = splicetype5_mm = ANTIACCEPTOR;

      start_genomicpos = start5;
      middle_genomicpos = middle5;
      end_genomicpos = end5;

      /* assert(start_genomicpos >= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos >= middle_genomicpos) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 4 */
	debug13(printf("5', crick, sense anti %u %u %f mm\n",genomicpos,genomicpos-chroffset,acceptor_prob));
	if (acceptor_prob > max_prob_5_mm) {
	  max_prob_5_mm = acceptor_prob;
	  splice_genomepos_5_mm = genomicpos;
	}
	genomicpos--;
      }
      while (genomicpos >= end_genomicpos) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 4 */
	debug13(printf("5', crick, sense anti %u %u %f\n",genomicpos,genomicpos-chroffset,acceptor_prob));
	if (acceptor_prob > max_prob_5) {
	  max_prob_5 = acceptor_prob;
	  splice_genomepos_5 = genomicpos;
	}
	genomicpos--;
      }
      debug13(printf("\n"));
    }

  } else if (orig_sensedir == SENSE_ANTI) {
    if (plusp) {
      splicetype5 = splicetype5_mm = ANTIDONOR;

      start_genomicpos = start5;
      middle_genomicpos = middle5;
      end_genomicpos = end5;
	
      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos <= middle_genomicpos) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	debug13(printf("5', watson, sense anti %u %u %f mm\n",genomicpos,genomicpos-chroffset,donor_prob));
	if (donor_prob > max_prob_5_mm) {
	  max_prob_5_mm = donor_prob;
	  splice_genomepos_5_mm = genomicpos;
	}
	genomicpos++;
      }
      while (genomicpos <= end_genomicpos) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	debug13(printf("5', watson, sense anti %u %u %f\n",genomicpos,genomicpos-chroffset,donor_prob));
	if (donor_prob > max_prob_5) {
	  max_prob_5 = donor_prob;
	  splice_genomepos_5 = genomicpos;
	}
	genomicpos++;
      }
      debug13(printf("\n"));

    } else {
      splicetype5 = splicetype5_mm = DONOR;

      start_genomicpos = start5;
      middle_genomicpos = middle5;
      end_genomicpos = end5;

      /* assert(start_genomicpos >= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos >= middle_genomicpos) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 8 */
	debug13(printf("5', crick, sense forward %u %u %f mm\n",genomicpos,genomicpos-chroffset,donor_prob));
	if (donor_prob > max_prob_5_mm) {
	  max_prob_5_mm = donor_prob;
	  splice_genomepos_5_mm = genomicpos;
	}
	genomicpos--;
      }
      while (genomicpos >= end_genomicpos) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 8 */
	debug13(printf("5', crick, sense forward %u %u %f\n",genomicpos,genomicpos-chroffset,donor_prob));
	if (donor_prob > max_prob_5) {
	  max_prob_5 = donor_prob;
	  splice_genomepos_5 = genomicpos;
	}
	genomicpos--;
      }
      debug13(printf("\n"));
    }
      
  } else {
    if (plusp) {
      start_genomicpos = start5;
      middle_genomicpos = middle5;
      end_genomicpos = end5;

      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos <= middle_genomicpos) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	debug13(printf("5', watson, sense null %u %u %f %f mm\n",genomicpos,genomicpos-chroffset,donor_prob,acceptor_prob));
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
      }
      while (genomicpos <= end_genomicpos) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	debug13(printf("5', watson, sense null %u %u %f %f\n",genomicpos,genomicpos-chroffset,donor_prob,acceptor_prob));
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
      start_genomicpos = start5;
      middle_genomicpos = middle5;
      end_genomicpos = end5;

      /* assert(start_genomicpos >= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos >= middle_genomicpos) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 4 */
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 8 */
	debug13(printf("5', crick, sense null %u %u %f %f mm\n",genomicpos,genomicpos-chroffset,donor_prob,acceptor_prob));
	if (acceptor_prob > max_prob_sense_forward_5_mm) {
	  max_prob_sense_forward_5_mm = acceptor_prob;
	  if (acceptor_prob > max_prob_5_mm) {
	    max_prob_5_mm = acceptor_prob;
	    splice_genomepos_5_mm = genomicpos;
	    /* splice_cdna_direction_5_mm = +1; */
	    splice_sensedir_5_mm = SENSE_FORWARD;
	    splicetype5_mm = ANTIACCEPTOR;
	  }
	}
	if (donor_prob > max_prob_sense_anti_5_mm) {
	  max_prob_sense_anti_5_mm = donor_prob;
	  if (donor_prob > max_prob_5_mm) {
	    max_prob_5_mm = donor_prob;
	    splice_genomepos_5_mm = genomicpos;
	    /* splice_cdna_direction_5_mm = -1; */
	    splice_sensedir_5_mm = SENSE_ANTI;
	    splicetype5_mm = DONOR;
	  }
	}
	genomicpos--;
      }
      while (genomicpos >= end_genomicpos) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 4 */
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 8 */
	debug13(printf("5', crick, sense null %u %u %f %f\n",genomicpos,genomicpos-chroffset,donor_prob,acceptor_prob));
	if (acceptor_prob > max_prob_sense_forward_5) {
	  max_prob_sense_forward_5 = acceptor_prob;
	  if (acceptor_prob > max_prob_5) {
	    max_prob_5 = acceptor_prob;
	    splice_genomepos_5 = genomicpos;
	    /* splice_cdna_direction_5 = +1; */
	    splice_sensedir_5 = SENSE_FORWARD;
	    splicetype5 = ANTIACCEPTOR;
	  }
	}
	if (donor_prob > max_prob_sense_anti_5) {
	  max_prob_sense_anti_5 = donor_prob;
	  if (donor_prob > max_prob_5) {
	    max_prob_5 = donor_prob;
	    splice_genomepos_5 = genomicpos;
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

  if (solve5p == false) {
    /* Skip 5' end */
  } else if (orig_sensedir != SENSE_NULL) {
    if (max_prob_5 > END_SPLICESITE_PROB_MATCH) {
      debug13(printf("Found good splice %s on 5' end at %u with probability %f\n",
		     Splicetype_string(splicetype5),splice_genomepos_5-chroffset,max_prob_5));
      if (plusp) {
	*ambig_end_length_5 = splice_genomepos_5 - genomicstart5;
	debug13(printf("1 Set ambig_end_length_5 to be %d = %u - %u\n",*ambig_end_length_5,splice_genomepos_5,genomicstart5));
      } else {
	*ambig_end_length_5 = genomicstart5 - splice_genomepos_5;
	debug13(printf("2 Set ambig_end_length_5 to be %d = %u - %u\n",*ambig_end_length_5,genomicstart5,splice_genomepos_5));
      }
      *ambig_splicetype_5 = splicetype5;
      *ambig_prob_5 = max_prob_5;
    } else if (max_prob_5_mm > END_SPLICESITE_PROB_MISMATCH) {
      debug13(printf("Found good mismatch splice %s on 5' end at %u with probability %f\n",
		     Splicetype_string(splicetype5_mm),splice_genomepos_5_mm-chroffset,max_prob_5_mm));
      if (plusp) {
	*ambig_end_length_5 = splice_genomepos_5_mm - genomicstart5;
	debug13(printf("3 Set ambig_end_length_5 to be %d = %u - %u\n",*ambig_end_length_5,splice_genomepos_5_mm,genomicstart5));
      } else {
	*ambig_end_length_5 = genomicstart5 - splice_genomepos_5_mm;
	debug13(printf("4 Set ambig_end_length_5 to be %d = %u - %u\n",*ambig_end_length_5,genomicstart5,splice_genomepos_5_mm));
      }
      *ambig_splicetype_5 = splicetype5_mm;
      *ambig_prob_5 = max_prob_5_mm;
    }
  }

  if (orig_sensedir == SENSE_NULL) {
    if (max_prob_3 >= END_SPLICESITE_PROB_MATCH || max_prob_5 >= END_SPLICESITE_PROB_MATCH) {
      if (max_prob_3 >= END_SPLICESITE_PROB_MATCH && max_prob_5 >= END_SPLICESITE_PROB_MATCH
	  && max_prob_sense_forward_3 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_anti_3 < END_SPLICESITE_PROB_MATCH
	  && max_prob_sense_forward_5 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_anti_5 < END_SPLICESITE_PROB_MATCH) {
	/* Forward sense wins on both sides */
	if (plusp) {
	  *ambig_end_length_3 = genomicend3 - splice_genomepos_3;
	  *ambig_end_length_5 = splice_genomepos_5 - genomicstart5;
	} else {
	  *ambig_end_length_3 = splice_genomepos_3 - genomicend3;
	  *ambig_end_length_5 = genomicstart5 - splice_genomepos_5;
	}
	*ambig_splicetype_3 = splicetype3;
	*ambig_prob_3 = max_prob_3;
	debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	*ambig_splicetype_5 = splicetype5;
	*ambig_prob_5 = max_prob_5;
	debug13(printf("Set ambig_end_length_5 to be %d\n",*ambig_end_length_5));
	new_sensedir = SENSE_FORWARD; /* = splice_sensedir_3 */

      } else if (max_prob_3 >= END_SPLICESITE_PROB_MATCH && max_prob_5 >= END_SPLICESITE_PROB_MATCH
		 && max_prob_sense_anti_3 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_forward_3 < END_SPLICESITE_PROB_MATCH
		 && max_prob_sense_anti_5 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_forward_5 < END_SPLICESITE_PROB_MATCH) {
	/* Anti sense wins on both sides */
	if (plusp) {
	  *ambig_end_length_3 = genomicend3 - splice_genomepos_3;
	  *ambig_end_length_5 = splice_genomepos_5 - genomicstart5;
	} else {
	  *ambig_end_length_3 = splice_genomepos_3 - genomicend3;
	  *ambig_end_length_5 = genomicstart5 - splice_genomepos_5;
	}
	*ambig_splicetype_3 = splicetype3;
	*ambig_prob_3 = max_prob_3;
	debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	*ambig_splicetype_5 = splicetype5;
	*ambig_prob_5 = max_prob_5;
	debug13(printf("Set ambig_end_length_5 to be %d\n",*ambig_end_length_5));
	new_sensedir = SENSE_ANTI; /* = splice_sensedir_3 */

      } else if (max_prob_3 > max_prob_5) {
	/* Consider just 3' end */
	debug13(printf("Found good splice %s on 3' end at %u with probability %f\n",
		       Splicetype_string(splicetype3),splice_genomepos_3-chroffset,max_prob_3));
	if (plusp) {
	  *ambig_end_length_3 = genomicend3 - splice_genomepos_3;
	} else {
	  *ambig_end_length_3 = splice_genomepos_3 - genomicend3;
	}
	*ambig_splicetype_3 = splicetype3;
	*ambig_prob_3 = max_prob_3;
	/* *cdna_direction = splice_cdna_direction_3; */
	debug13(printf("Set ambig_end_length_3 to be %d\n",*ambig_end_length_3));
	if (max_prob_sense_forward_3 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_anti_3 < END_SPLICESITE_PROB_MATCH
	    && max_prob_sense_anti_5 < END_SPLICESITE_PROB_MATCH) {
	  new_sensedir = splice_sensedir_3;
	} else if (max_prob_sense_anti_3 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_forward_3 < END_SPLICESITE_PROB_MATCH
		   && max_prob_sense_forward_5 < END_SPLICESITE_PROB_MATCH) {
	  new_sensedir = splice_sensedir_3;
	} else {
	  /* Not enough evidence to set sensedir */
	}

      } else {
	/* Consider just 5' end */
	debug13(printf("Found good splice %s on 5' end at %u with probability %f\n",
		       Splicetype_string(splicetype5),splice_genomepos_5-chroffset,max_prob_5));
	if (plusp) {
	  *ambig_end_length_5 = splice_genomepos_5 - genomicstart5;
	  debug13(printf("5 Set ambig_end_length_5 to be %d = %u - %u\n",*ambig_end_length_5,splice_genomepos_5,genomicstart5));
	} else {
	  *ambig_end_length_5 = genomicstart5 - splice_genomepos_5;
	  debug13(printf("6 Set ambig_end_length_5 to be %d = %u - %u\n",*ambig_end_length_5,genomicstart5,splice_genomepos_5));
	}
	*ambig_splicetype_5 = splicetype5;
	*ambig_prob_5 = max_prob_5;
	/* *cdna_direction = splice_cdna_direction_5; */
	if (max_prob_sense_forward_5 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_anti_5 < END_SPLICESITE_PROB_MATCH
	    && max_prob_sense_anti_3 < END_SPLICESITE_PROB_MATCH) {
	  new_sensedir = splice_sensedir_5;
	} else if (max_prob_sense_anti_5 >= END_SPLICESITE_PROB_MATCH && max_prob_sense_forward_5 < END_SPLICESITE_PROB_MATCH
		   && max_prob_sense_forward_3 < END_SPLICESITE_PROB_MATCH) {
	  new_sensedir = splice_sensedir_5;
	} else {
	  /* Not enough evidence to set sensedir */
	}
      }

    } else if (max_prob_3_mm >= END_SPLICESITE_PROB_MISMATCH || max_prob_5_mm >= END_SPLICESITE_PROB_MISMATCH) {
      if (max_prob_3_mm > max_prob_5_mm) {
	debug13(printf("Found good mismatch splice %s on 3' end at %u with probability %f\n",
		       Splicetype_string(splicetype3_mm),splice_genomepos_3_mm-chroffset,max_prob_3_mm));
	if (plusp) {
	  *ambig_end_length_3 = genomicend3 - splice_genomepos_3_mm;
	  debug13(printf("Set ambig_end_length_3 to be %d = %u - %u\n",*ambig_end_length_3,genomicend3,splice_genomepos_3_mm));
	} else {
	  *ambig_end_length_3 = splice_genomepos_3_mm - genomicend3;
	  debug13(printf("Set ambig_end_length_3 to be %d = %u - %u\n",*ambig_end_length_3,splice_genomepos_3_mm,genomicend3));
	}
	*ambig_splicetype_3 = splicetype3_mm;
	*ambig_prob_3 = max_prob_3_mm;
	/* *cdna_direction = splice_cdna_direction_3_mm; */
	if (max_prob_sense_forward_3_mm >= END_SPLICESITE_PROB_MISMATCH && max_prob_sense_anti_3_mm < END_SPLICESITE_PROB_MISMATCH
	    && max_prob_sense_anti_5_mm < END_SPLICESITE_PROB_MISMATCH) {
	  new_sensedir = splice_sensedir_3_mm;
	} else if (max_prob_sense_anti_3_mm >= END_SPLICESITE_PROB_MISMATCH && max_prob_sense_forward_3_mm < END_SPLICESITE_PROB_MISMATCH
		   && max_prob_sense_forward_5_mm < END_SPLICESITE_PROB_MISMATCH) {
	  new_sensedir = splice_sensedir_3_mm;
	} else {
	  /* Not enough evidence to set sensedir */
	}
      } else {
	debug13(printf("Found good mismatch splice %s on 5' end at %u with probability %f\n",
		       Splicetype_string(splicetype5_mm),splice_genomepos_5_mm-chroffset,max_prob_5_mm));
	if (plusp) {
	  *ambig_end_length_5 = splice_genomepos_5_mm - genomicstart5;
	  debug13(printf("7 Set ambig_end_length_5 to be %d = %u - %u\n",*ambig_end_length_5,splice_genomepos_5_mm,genomicstart5));
	} else {
	  *ambig_end_length_5 = genomicstart5 - splice_genomepos_5_mm;
	  debug13(printf("8 Set ambig_end_length_5 to be %d = %u - %u\n",*ambig_end_length_5,genomicstart5,splice_genomepos_5_mm));
	}
	*ambig_splicetype_5 = splicetype5_mm;
	*ambig_prob_5 = max_prob_5_mm;
	/* *cdna_direction = splice_cdna_direction_5_mm; */
	if (max_prob_sense_forward_5_mm >= END_SPLICESITE_PROB_MISMATCH && max_prob_sense_anti_5_mm < END_SPLICESITE_PROB_MISMATCH
	    && max_prob_sense_anti_3_mm < END_SPLICESITE_PROB_MISMATCH) {
	  new_sensedir = splice_sensedir_5_mm;
	} else if (max_prob_sense_anti_5_mm >= END_SPLICESITE_PROB_MISMATCH && max_prob_sense_forward_5_mm < END_SPLICESITE_PROB_MISMATCH
		   && max_prob_sense_forward_3_mm < END_SPLICESITE_PROB_MISMATCH) {
	  new_sensedir = splice_sensedir_5_mm;
	} else {
	  /* Not enough evidence to set sensedir */
	}
      }
    }
  }

  debug13(printf("Returning ambig_end_length_5 %d and ambig_end_length_3 %d, probs %f and %f, new_sensedir %d\n",
		 *ambig_end_length_5,*ambig_end_length_3,*ambig_prob_5,*ambig_prob_3,new_sensedir));
  return new_sensedir;
}

