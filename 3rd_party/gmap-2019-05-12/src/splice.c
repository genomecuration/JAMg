static char rcsid[] = "$Id: splice.c 218675 2019-03-16 01:25:48Z twu $";
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
#include "indel.h"
#include "univcoord.h"


#define LOWPROB_SUPPORT 20
#define MIN_SUPPORT_SPLICE 8
#define MIN_SUPPORT_SPLICE_PLUS_INDEL 12

#define END_SPLICESITE_SEARCH_MM 1 /* Amount to search in the trimmed area */
#define END_SPLICESITE_SEARCH 10
#define MIN_EXON_LENGTH 9

#define MAX_NCONSECUTIVE 20	/* Needs to be generous to find splices */

#if 0
/* Creates issues with ambiguous substrings */
#define LOCALSPLICING_NMATCHES_SLOP 1
#else
#define LOCALSPLICING_NMATCHES_SLOP 0
#endif
#define LOCALSPLICING_PROB_SLOP 0.05

#define SLOP 1


/* Splice_resolve_sense and Splice_resolve_antisense */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Splice_resolve_distant */
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

static Genome_T genomebits;
static Genome_T genomebits_alt;

static int min_shortend;
static int min_shortend_distant = 20;


void
Splice_setup (Genome_T genomebits_in, Genome_T genomebits_alt_in, int min_shortend_in) {
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  min_shortend = min_shortend_in;
  /* max_end_length = 2*(index1part_in + index1interval_in - 1); */
  return;
}


void
Spliceinfo_free (Spliceinfo_T *old) {
  FREE((*old)->double_memory);
  FREE((*old)->int_memory);
  FREE(*old);
  return;
}


Spliceinfo_T
Spliceinfo_new (int querylength) {
  Spliceinfo_T new = (Spliceinfo_T) MALLOC(sizeof(*new));
  int a, b = (querylength + 1);

  new->int_memory = (int *) MALLOC(12 * b * sizeof(int));

  /* 12 arrays of (int *) */
  new->donor_positions_alloc = &(new->int_memory[0]); a = b;
  new->donor_knowni_alloc = &(new->int_memory[a]); a += b;
  new->acceptor_positions_alloc = &(new->int_memory[a]); a += b;
  new->acceptor_knowni_alloc = &(new->int_memory[a]); a += b;

  new->segmenti_donor_knownpos = &(new->int_memory[a]); a += b;
  new->segmentj_acceptor_knownpos = &(new->int_memory[a]); a += b;
  new->segmentj_antidonor_knownpos = &(new->int_memory[a]); a += b;
  new->segmenti_antiacceptor_knownpos = &(new->int_memory[a]); a += b;

  new->segmenti_donor_knowni = &(new->int_memory[a]); a += b;
  new->segmentj_acceptor_knowni = &(new->int_memory[a]); a += b;
  new->segmentj_antidonor_knowni = &(new->int_memory[a]); a += b;
  new->segmenti_antiacceptor_knowni = &(new->int_memory[a]); 

  new->double_memory = (double *) MALLOC(2 * b * sizeof(double));
  new->donor_probs_alloc = &(new->double_memory[0]); a = b;
  new->acceptor_probs_alloc = &(new->double_memory[a]);

  return new;
}

/* Note: contents of spliceinfo are filled in by kmer-search and path-solve procedures */


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

#if 0
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
#endif




/* Note: knowni holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

int
Splice_resolve_sense (int *best_nindels, int *best_indel_pos, int *best_knowni_i, int *best_knowni_j,
		      int *best_nmismatches_i, int *best_nmismatches_j, int *best_nmismatches_indel,
		      double *best_prob_i, double *best_prob_j,

		      Univcoord_T segmenti_left, Univcoord_T segmentj_left,
		      Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
		     
		      int querystart, int queryend, int querylength, Compress_T query_compress,
		      Spliceinfo_T spliceinfo, int max_mismatches_allowed,
		      bool plusp, int genestrand, int max_deletionlen, int max_insertionlen,
		      bool allow_indel_p) {
  int best_splice_qpos = -1, splice_qpos_start, splice_qpos_end, splice_qpos, i, j;

  int best_nmismatches, nmismatches, nmismatches1, nmismatches2;
  int segmenti_nmismatches, segmentj_nmismatches;
  double best_prob, probi, probj;
  /* bool sufficient1p, sufficient2p; */

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;
  double *donori_probs, *acceptorj_probs, *antiacceptori_probs, *antidonorj_probs;

  int supporti, supportj;
  int nindels, indel_pos;
  int nmismatches_indel, nmismatches_i, nmismatches_j;


  debug1(printf("Splice_resolve_sense: Getting genome at lefti %u and leftj %u (diff: %d), range %d..%d, max_mismatches_allowed %d\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left,querystart,queryend,max_mismatches_allowed));

  *best_nindels = 0;
  *best_indel_pos = -1;
  *best_knowni_i = *best_knowni_j = -1;
  *best_nmismatches_i = *best_nmismatches_j = *best_nmismatches_indel = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
  *best_prob_i = *best_prob_j = 0.0;

  /* Require separation from endpoints */
  splice_qpos_start = querystart + 1;
  splice_qpos_end = queryend - 1;
  if (splice_qpos_start < min_shortend) {
    splice_qpos_start = min_shortend;
  }
  if (splice_qpos_end > querylength - min_shortend) {
    splice_qpos_end = querylength - min_shortend;
  }

  if (plusp == true) {
    /* Originally from plus strand.  No complement.  */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
    if (novelsplicingp && segmenti_left + splice_qpos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(spliceinfo->donor_positions_alloc,spliceinfo->donor_knowni_alloc,
					     spliceinfo->segmenti_donor_knownpos,
					     spliceinfo->segmenti_donor_knowni,
					     segmenti_left,splice_qpos_start,splice_qpos_end);
      donori_positions = spliceinfo->donor_positions_alloc;
      donori_knowni = spliceinfo->donor_knowni_alloc;
    } else {
      donori_nsites = spliceinfo->segmenti_donor_nknown;
      donori_positions = spliceinfo->segmenti_donor_knownpos;
      donori_knowni = spliceinfo->segmenti_donor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      } else {
	probi = Maxent_hr_donor_prob(segmenti_left + donori_positions[i],segmenti_chroffset);
	printf(" (%.6f)",probi);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_qpos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(spliceinfo->acceptor_positions_alloc,spliceinfo->acceptor_knowni_alloc,
						   spliceinfo->segmentj_acceptor_knownpos,
						   spliceinfo->segmentj_acceptor_knowni,
						   segmentj_left,splice_qpos_start,splice_qpos_end);
      acceptorj_positions = spliceinfo->acceptor_positions_alloc;
      acceptorj_knowni = spliceinfo->acceptor_knowni_alloc;
    } else {
      acceptorj_nsites = spliceinfo->segmentj_acceptor_nknown;
      acceptorj_positions = spliceinfo->segmentj_acceptor_knownpos;
      acceptorj_knowni = spliceinfo->segmentj_acceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      } else {
	probj = Maxent_hr_acceptor_prob(segmentj_left + acceptorj_positions[i],segmentj_chroffset);
	printf(" (%.6f)",probj);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < donori_nsites && j < acceptorj_nsites) {
      supporti = donori_positions[i] - querystart;
      supportj = queryend - acceptorj_positions[j];

      if ((splice_qpos = donori_positions[i]) < acceptorj_positions[j]) {
	i++;
      } else if (splice_qpos > acceptorj_positions[j]) {
	j++;
      } else if (supporti < MIN_SUPPORT_SPLICE || supportj < MIN_SUPPORT_SPLICE) {
	/* Skip */
	i++; j++;
      } else {
	debug1(printf("splice matches at %d\n",splice_qpos));
	segmenti_nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								 /*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								 /*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	debug1(printf("%d mismatches on segmenti (%d..%d)\n",segmenti_nmismatches,querystart,splice_qpos));
	debug1(printf("%d mismatches on segmentj (%d..%d)\n",segmentj_nmismatches,splice_qpos,queryend));

	if (supporti - 3*segmenti_nmismatches < MIN_SUPPORT_SPLICE) {
	  /* Skip, because too many mismatches in segmenti */
	  debug1(printf("Skipping, because too many mismatches %d in segmenti\n",
			segmenti_nmismatches));
	} else if (supportj - 3*segmentj_nmismatches < MIN_SUPPORT_SPLICE) {
	  /* Skip, because too many mismatches in segmentj */
	  debug1(printf("Skipping, because too many mismatches %d in segmentj\n",
			segmentj_nmismatches));
	} else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) > best_nmismatches) {
	  debug1(printf("Skipping, because too many nmismatches %d > best_nmismatches %d\n",
			nmismatches,best_nmismatches));
	} else {
	  debug1(printf("nmismatches %d + %d <= best_nmismatches %d\n",segmenti_nmismatches,segmentj_nmismatches,best_nmismatches));
	  if (donori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	    /* sufficient1p = true; */
	  } else {
	    probi = Maxent_hr_donor_prob(segmenti_left + splice_qpos,segmenti_chroffset);
	    /* sufficient1p = sufficient_splice_prob_local(splice_qpos,segmenti_nmismatches,probi); */
	  }

	  if (acceptorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	    /* sufficient2p = true; */
	  } else {
	    probj = Maxent_hr_acceptor_prob(segmentj_left + splice_qpos,segmentj_chroffset);
	    /* sufficient2p = sufficient_splice_prob_local(querylength - splice_qpos,segmentj_nmismatches,probj); */
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus sense splice_qpos  %d, i.donor %f, j.acceptor %f\n",splice_qpos,probi,probj);
		 } else {
		   printf("minus antisense splice_qpos  %d, i.donor %f, j.acceptor %f\n",splice_qpos,probi,probj);
		 });

	  if (nmismatches < best_nmismatches ||
	      (nmismatches == best_nmismatches && probi + probj > best_prob)) {
	    /* Success */
	    best_nmismatches = nmismatches;
	    best_prob = probi + probj;
	    
	    /* best_donor_splicecoord = segmenti_left + splice_qpos; */
	    /* best_acceptor_splicecoord = segmentj_left + splice_qpos; */
	    *best_knowni_i = donori_knowni[i];
	    *best_knowni_j = acceptorj_knowni[j];
	    *best_prob_i = probi; /* donor_prob */
	    *best_prob_j = probj; /* acceptor_prob */
	    best_splice_qpos = splice_qpos;
	    *best_nmismatches_i = segmenti_nmismatches;
	    *best_nmismatches_j = segmentj_nmismatches;
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
    if (novelsplicingp && segmenti_left + splice_qpos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(spliceinfo->acceptor_positions_alloc,spliceinfo->acceptor_knowni_alloc,
							   spliceinfo->segmenti_antiacceptor_knownpos,
							   spliceinfo->segmenti_antiacceptor_knowni,
							   segmenti_left,splice_qpos_start,splice_qpos_end);
      antiacceptori_positions = spliceinfo->acceptor_positions_alloc;
      antiacceptori_knowni = spliceinfo->acceptor_knowni_alloc;
    } else {
      antiacceptori_nsites = spliceinfo->segmenti_antiacceptor_nknown;
      antiacceptori_positions = spliceinfo->segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = spliceinfo->segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      } else {
	probi = Maxent_hr_antiacceptor_prob(segmenti_left + antiacceptori_positions[i],segmenti_chroffset);
	printf(" (%.6f)",probi);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_qpos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(spliceinfo->donor_positions_alloc,spliceinfo->donor_knowni_alloc,
						     spliceinfo->segmentj_antidonor_knownpos,
						     spliceinfo->segmentj_antidonor_knowni,
						     segmentj_left,splice_qpos_start,splice_qpos_end);
      antidonorj_positions = spliceinfo->donor_positions_alloc;
      antidonorj_knowni = spliceinfo->donor_knowni_alloc;
    } else {
      antidonorj_nsites = spliceinfo->segmentj_antidonor_nknown;
      antidonorj_positions = spliceinfo->segmentj_antidonor_knownpos;
      antidonorj_knowni = spliceinfo->segmentj_antidonor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      } else {
	probj = Maxent_hr_antidonor_prob(segmentj_left + antidonorj_positions[i],segmentj_chroffset);
	printf(" (%.6f)",probj);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < antiacceptori_nsites && j < antidonorj_nsites) {
      supporti = antiacceptori_positions[i] - querystart;
      supportj = queryend - antidonorj_positions[j];

      if ((splice_qpos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	i++;
      } else if (splice_qpos > antidonorj_positions[j]) {
	j++;
      } else if (supporti < MIN_SUPPORT_SPLICE || supportj < MIN_SUPPORT_SPLICE) {
	/* Skip */
	i++; j++;
      } else {
	debug1(printf("splice matches at %d\n",splice_qpos));
	segmenti_nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								 /*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								 /*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	debug1(printf("%d mismatches on segmenti (%d..%d)\n",segmenti_nmismatches,querystart,splice_qpos));
	debug1(printf("%d mismatches on segmentj (%d..%d)\n",segmentj_nmismatches,splice_qpos,queryend));

	if (supporti - 3*segmenti_nmismatches < MIN_SUPPORT_SPLICE) {
	  /* Skip, because too many mismatches in segmenti */
	  debug1(printf("Skipping, because too many mismatches %d in segmenti\n",
			segmenti_nmismatches));
	} else if (supportj - 3*segmentj_nmismatches < MIN_SUPPORT_SPLICE) {
	  /* Skip, because too many mismatches in segmentj */
	  debug1(printf("Skipping, because too many mismatches %d in segmentj\n",
			segmentj_nmismatches));
	} else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) > best_nmismatches) {
	  debug1(printf("Skipping, because too many nmismatches %d > best_nmismatches %d\n",
			nmismatches,best_nmismatches));
	} else {
	  debug1(printf("nmismatches %d + %d <= best_nmismatches %d\n",segmenti_nmismatches,segmentj_nmismatches,best_nmismatches));
	  if (antiacceptori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	    /* sufficient1p = true; */
	  } else {
	    probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_qpos,segmenti_chroffset);
	    /* sufficient1p = sufficient_splice_prob_local(splice_qpos,segmenti_nmismatches,probi); */
	  }

	  if (antidonorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	    /* sufficient2p = true; */
	  } else {
	    probj = Maxent_hr_antidonor_prob(segmentj_left + splice_qpos,segmentj_chroffset);
	    /* sufficient2p = sufficient_splice_prob_local(querylength - splice_qpos,segmentj_nmismatches,probj); */
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus antisense splice_qpos  %d, j.donor %f, i.acceptor %f\n",splice_qpos,probj,probi);
		 } else {
		   printf("minus sense splice_qpos  %d, j.donor %f, i.acceptor %f\n",splice_qpos,probj,probi);
		 });
	  
	  if (nmismatches < best_nmismatches ||
	      (nmismatches == best_nmismatches && probi + probj > best_prob)) {
	    /* Success */
	    best_nmismatches = nmismatches;
	    best_prob = probi + probj;
	    
	    /* best_donor_splicecoord = segmentj_left + splice_qpos; */
	    /* best_acceptor_splicecoord = segmenti_left + splice_qpos; */
	    *best_knowni_j = antidonorj_knowni[j];
	    *best_knowni_i = antiacceptori_knowni[i];
	    *best_prob_j = probj; /* donor_prob */
	    *best_prob_i = probi;
	    best_splice_qpos = splice_qpos;
	    *best_nmismatches_j = segmentj_nmismatches;
	    *best_nmismatches_i = segmenti_nmismatches;
	  }
	}
	i++;
	j++;
      }
    }
  }

  debug1(printf("best_knowni_i is %d and best_knowni_j is %d\n",*best_knowni_i,*best_knowni_j));

  if (*best_prob_i >= 0.9 || *best_prob_j >= 0.9) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
    debug1(printf("nmismatches %d and %d\n",*best_nmismatches_i,*best_nmismatches_j));
    return best_splice_qpos;
  } else if (allow_indel_p == false) {
    debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
    *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
    return -1;
  } else {
    /* Find an indel below */
  }

  if (plusp == true) {
    /* Genomic sense (inferred from antisense code) */
    if (donori_nsites == 0 || acceptorj_nsites == 0) {
      debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
      *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
      return -1;
    } else {
      donori_probs = spliceinfo->donor_probs_alloc;
      for (i = 0; i < donori_nsites; i++) {
	donori_probs[i] = Maxent_hr_donor_prob(segmenti_left + donori_positions[i],segmenti_chroffset);
      }

      acceptorj_probs = spliceinfo->acceptor_probs_alloc;
      for (i = 0; i < acceptorj_nsites; i++) {
	acceptorj_probs[i] = Maxent_hr_acceptor_prob(segmentj_left + acceptorj_positions[i],segmentj_chroffset);
      }

      i = j = 0;
      best_nmismatches = max_mismatches_allowed + 1;
      while (i < donori_nsites && donori_positions[i] - querystart < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	i++;
      }

      while (i < donori_nsites) {
	supporti = donori_positions[i] - querystart;

	/* Backup */
	while (j >= 0 && acceptorj_positions[j] + max_deletionlen > donori_positions[i]) {
	  j--;
	}
	j++;			/* Finish backup */

	/* Advance */
	while (j < acceptorj_nsites && acceptorj_positions[j] + max_deletionlen <= donori_positions[i]) {
	  j++;
	}

	/* Deletions */
	while (j < acceptorj_nsites && acceptorj_positions[j] < donori_positions[i] &&
	       (supportj = queryend - acceptorj_positions[j]) >= MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	  debug1(printf("Deletion: %d (%.6f) -- %d (%.6f).  Support: %d and %d\n",
			donori_positions[i],donori_probs[i],acceptorj_positions[j],acceptorj_probs[j],
			donori_positions[i] - querystart,queryend - acceptorj_positions[j]));
	  probi = donori_probs[i];
	  probj = acceptorj_probs[j];
	  if (probi >= 0.85 && probj >= 0.85) {
	    nindels = acceptorj_positions[j] - donori_positions[i]; /* Should be negative */

	    /* Try deletion on segmenti */
	    splice_qpos = acceptorj_positions[j];
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,segmenti_left,nindels,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							   querystart,/*queryend*/splice_qpos,querylength,
							   max_mismatches_allowed,/*plusp:true*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > querystart);
	      assert(indel_pos < splice_qpos);
	      nmismatches_indel = nmismatches1;
	      nmismatches_i = nmismatches2;
	      nmismatches_j = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
									 /*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	      debug1(printf("  => deletion indel_pos %d on segmenti with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_indel,nmismatches_i,nmismatches_j));

	      if (supporti - 3*(nmismatches_indel + nmismatches_i) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supportj - 3*nmismatches_j < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }

	    /* Try deletion on segmentj */
	    splice_qpos = donori_positions[i];
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,segmentj_left + nindels,nindels,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							   /*querystart*/splice_qpos,queryend,querylength,
							   max_mismatches_allowed,/*plusp:true*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos < queryend);
	      nmismatches_i = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								/*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	      nmismatches_j = nmismatches1;
	      nmismatches_indel = nmismatches2;
	      debug1(printf("  => deletion indel_pos %d on segmentj with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_i,nmismatches_j,nmismatches_indel));

	      if (supportj - 3*(nmismatches_indel + nmismatches_j) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supporti - 3*nmismatches_i < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_i + nmismatches_j + nmismatches_indel) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }
	  }
	  j++;
	}

	if (j < acceptorj_nsites && acceptorj_positions[j] == donori_positions[i]) {
	  /* Not an indel */
	  j++;
	}

	/* Insertions */
	while (j < acceptorj_nsites && acceptorj_positions[j] <= donori_positions[i] + max_insertionlen &&
	       (supportj = queryend - acceptorj_positions[j]) >= MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	  debug1(printf("Insertion: %d (%.6f) -- %d (%.6f).  Support: %d and %d\n",
			donori_positions[i],donori_probs[i],acceptorj_positions[j],acceptorj_probs[j],
			donori_positions[i] - querystart,queryend - acceptorj_positions[j]));
	  probi = donori_probs[i];
	  probj = acceptorj_probs[j];
	  if (probi >= 0.85 && probj >= 0.85) {
	    nindels = acceptorj_positions[j] - donori_positions[i]; /* Should be positive */

	    /* Try insertion on segmenti */
	    splice_qpos = acceptorj_positions[j];
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,segmenti_left,nindels,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    querystart,/*queryend*/splice_qpos,querylength,
							    max_mismatches_allowed,/*plusp:true*/true,genestrand,
							    /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > querystart);
	      assert(indel_pos + nindels < splice_qpos);
	      nmismatches_indel = nmismatches1;
	      nmismatches_i = nmismatches2;
	      nmismatches_j = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								/*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	      debug1(printf("  => insertion indel_pos %d on segmenti with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_indel,nmismatches_i,nmismatches_j));

	      if (supporti - 3*(nmismatches_indel + nmismatches_i) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supportj - 3*nmismatches_j < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }

	    /* Try insertion on segmentj */
	    splice_qpos = donori_positions[i];
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,segmentj_left + nindels,nindels,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    /*querystart*/splice_qpos,queryend,querylength,
							    max_mismatches_allowed,/*plusp:true*/true,genestrand,
							    /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos + nindels < queryend);
	      nmismatches_i = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								/*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	      nmismatches_j = nmismatches1;
	      nmismatches_indel = nmismatches2;
	      debug1(printf("  => insertion indel_pos %d on segmentj with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_i,nmismatches_j,nmismatches_indel));
	      
	      if (supportj - 3*(nmismatches_indel + nmismatches_j) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supporti - 3*nmismatches_i < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_i + nmismatches_j + nmismatches_indel) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }
	  }
	  j++;
	}
	j--;			/* Finish advance */

	i++;
      }

      debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
      *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
      return -1;
    }

  } else {
    /* Genomic antisense (verified) */
    if (antiacceptori_nsites == 0 || antidonorj_nsites == 0) {
      debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
      *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
      return -1;
    } else {
      antiacceptori_probs = spliceinfo->acceptor_probs_alloc;
      for (i = 0; i < antiacceptori_nsites; i++) {
	antiacceptori_probs[i] = Maxent_hr_antiacceptor_prob(segmenti_left + antiacceptori_positions[i],segmenti_chroffset);
      }

      antidonorj_probs = spliceinfo->donor_probs_alloc;
      for (i = 0; i < antidonorj_nsites; i++) {
	antidonorj_probs[i] = Maxent_hr_antidonor_prob(segmentj_left + antidonorj_positions[i],segmentj_chroffset);
      }

      i = j = 0;
      best_nmismatches = max_mismatches_allowed + 1;
      while (i < antiacceptori_nsites && antiacceptori_positions[i] - querystart < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	i++;
      }

      while (i < antiacceptori_nsites) {
	supporti = antiacceptori_positions[i] - querystart;

	/* Backup */
	while (j >= 0 && antidonorj_positions[j] + max_deletionlen > antiacceptori_positions[i]) {
	  j--;
	}
	j++;			/* Finish backup */

	/* Advance */
	while (j < antidonorj_nsites && antidonorj_positions[j] + max_deletionlen <= antiacceptori_positions[i]) {
	  j++;
	}

	/* Deletions */
	while (j < antidonorj_nsites && antidonorj_positions[j] < antiacceptori_positions[i] &&
	       (supportj = queryend - antidonorj_positions[j]) >= MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	  debug1(printf("Deletion: %d (%.6f) -- %d (%.6f).  Support: %d and %d\n",
			antiacceptori_positions[i],antiacceptori_probs[i],antidonorj_positions[j],antidonorj_probs[j],
			antiacceptori_positions[i] - querystart,queryend - antidonorj_positions[j]));
	  probi = antiacceptori_probs[i];
	  probj = antidonorj_probs[j];
	  if (probi >= 0.85 && probj >= 0.85) {
	    nindels = antidonorj_positions[j] - antiacceptori_positions[i]; /* Should be negative */

	    /* Try deletion on segmenti */
	    splice_qpos = antidonorj_positions[j];
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,segmenti_left,nindels,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							   querystart,/*queryend*/splice_qpos,querylength,
							   max_mismatches_allowed,/*plusp:true*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > querystart);
	      assert(indel_pos < splice_qpos);
	      nmismatches_indel = nmismatches1;
	      nmismatches_i = nmismatches2;
	      nmismatches_j = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								/*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	      debug1(printf("  => deletion indel_pos %d on segmenti with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_indel,nmismatches_i,nmismatches_j));

	      if (supporti - 3*(nmismatches_indel + nmismatches_i) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supportj - 3*nmismatches_j < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }

	    /* Try deletion on segmentj */
	    splice_qpos = antiacceptori_positions[i];
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,segmentj_left + nindels,nindels,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							   /*querystart*/splice_qpos,queryend,querylength,
							   max_mismatches_allowed,/*plusp:true*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos < queryend);
	      nmismatches_i = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								/*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	      nmismatches_j = nmismatches1;
	      nmismatches_indel = nmismatches2;
	      debug1(printf("  => deletion indel_pos %d on segmentj with mismatches %d+%d+%d,  ",
			    indel_pos,nmismatches_i,nmismatches_j,nmismatches_indel));
	      
	      if (supportj - 3*(nmismatches_indel + nmismatches_j) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supporti - 3*nmismatches_i < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_i + nmismatches_j + nmismatches_indel) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }
	  }
	  j++;
	}

	if (j < antidonorj_nsites && antidonorj_positions[j] == antiacceptori_positions[i]) {
	  /* Not an indel */
	  j++;
	}

	/* Insertions */
	while (j < antidonorj_nsites && antidonorj_positions[j] <= antiacceptori_positions[i] + max_insertionlen &&
	       (supportj = queryend - antidonorj_positions[j]) >= MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	  debug1(printf("Insertion: %d (%.6f) -- %d (%.6f).  Support: %d and %d\n",
			antiacceptori_positions[i],antiacceptori_probs[i],antidonorj_positions[j],antidonorj_probs[j],
			antiacceptori_positions[i] - querystart,queryend - antidonorj_positions[j]));
	  probi = antiacceptori_probs[i];
	  probj = antidonorj_probs[j];
	  if (probi >= 0.85 && probj >= 0.85) {
	    nindels = antidonorj_positions[j] - antiacceptori_positions[i]; /* Should be positive */

	    /* Try insertion on segmenti */
	    splice_qpos = antidonorj_positions[j];
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,segmenti_left,nindels,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    querystart,/*queryend*/splice_qpos,querylength,
							    max_mismatches_allowed,/*plusp:true*/true,genestrand,
							    /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > querystart);
	      assert(indel_pos + nindels < splice_qpos);
	      nmismatches_indel = nmismatches1;
	      nmismatches_i = nmismatches2;
	      nmismatches_j = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								/*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	      debug1(printf("  => insertion indel_pos %d on segmenti with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_indel,nmismatches_i,nmismatches_j));

	      if (supporti - 3*(nmismatches_indel + nmismatches_i) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supportj - 3*nmismatches_j < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }

	    /* Try insertion on segmentj */
	    splice_qpos = antiacceptori_positions[i];
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,segmentj_left + nindels,nindels,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    /*querystart*/splice_qpos,queryend,querylength,
							    max_mismatches_allowed,/*plusp:true*/true,genestrand,
							    /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos + nindels < queryend);
	      nmismatches_i = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								/*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	      nmismatches_j = nmismatches1;
	      nmismatches_indel = nmismatches2;
	      debug1(printf("  => insertion indel_pos %d on segmentj with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_i,nmismatches_j,nmismatches_indel));
	      
	      if (supportj - 3*(nmismatches_indel + nmismatches_j) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supporti - 3*nmismatches_i < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_i + nmismatches_j + nmismatches_indel) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }
	  }
	  j++;
	}
	j--;			/* Finish advance */

	i++;
      }
    }
  }

  if (*best_nindels == 0) {
    debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
    *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
    return -1;					    /* Could not find a splice plus indel */
  } else {
    return best_splice_qpos;
  }
}


int
Splice_resolve_antisense (int *best_nindels, int *best_indel_pos, int *best_knowni_i, int *best_knowni_j,
			  int *best_nmismatches_i, int *best_nmismatches_j, int *best_nmismatches_indel,
			  double *best_prob_i, double *best_prob_j,
			  
			  Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			  Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
			  
			  int querystart, int queryend, int querylength, Compress_T query_compress,
			  Spliceinfo_T spliceinfo, int max_mismatches_allowed,
			  bool plusp, int genestrand, int max_deletionlen, int max_insertionlen,
			  bool allow_indel_p) {
  int best_splice_qpos = -1, splice_qpos_start, splice_qpos_end, splice_qpos, i, j;

  int best_nmismatches, nmismatches, nmismatches1, nmismatches2;
  int segmenti_nmismatches, segmentj_nmismatches;
  double best_prob, probi, probj;
  /* bool sufficient1p, sufficient2p; */

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;
  double *donori_probs, *acceptorj_probs, *antiacceptori_probs, *antidonorj_probs;

  int supporti, supportj;
  int nindels, indel_pos;
  int nmismatches_indel, nmismatches_i, nmismatches_j;


  debug1(printf("Splice_resolve_antisense: Getting genome at lefti %u and leftj %u (diff: %d), range %d..%d, max_mismatches_allowed %d\n",
		segmenti_left,segmentj_left,segmentj_left-segmenti_left,querystart,queryend,max_mismatches_allowed));

  *best_nindels = 0;
  *best_indel_pos = -1;
  *best_knowni_i = *best_knowni_j = -1;
  *best_nmismatches_i = *best_nmismatches_j = *best_nmismatches_indel = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
  *best_prob_i = *best_prob_j = 0.0;

  /* Require separation from endpoints */
  splice_qpos_start = querystart + 1;
  splice_qpos_end = queryend - 1;
  if (splice_qpos_start < min_shortend) {
    splice_qpos_start = min_shortend;
  }
  if (splice_qpos_end > querylength - min_shortend) {
    splice_qpos_end = querylength - min_shortend;
  }

  if (plusp == false) {
    /* minus */
    /* Originally from plus strand.  No complement.  */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
    if (novelsplicingp && segmenti_left + splice_qpos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(spliceinfo->donor_positions_alloc,spliceinfo->donor_knowni_alloc,
					     spliceinfo->segmenti_donor_knownpos,
					     spliceinfo->segmenti_donor_knowni,
					     segmenti_left,splice_qpos_start,splice_qpos_end);
      donori_positions = spliceinfo->donor_positions_alloc;
      donori_knowni = spliceinfo->donor_knowni_alloc;
    } else {
      donori_nsites = spliceinfo->segmenti_donor_nknown;
      donori_positions = spliceinfo->segmenti_donor_knownpos;
      donori_knowni = spliceinfo->segmenti_donor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      } else {
	probi = Maxent_hr_donor_prob(segmenti_left + donori_positions[i],segmenti_chroffset);
	printf(" (%.6f)",probi);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_qpos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(spliceinfo->acceptor_positions_alloc,spliceinfo->acceptor_knowni_alloc,
						   spliceinfo->segmentj_acceptor_knownpos,
						   spliceinfo->segmentj_acceptor_knowni,
						   segmentj_left,splice_qpos_start,splice_qpos_end);
      acceptorj_positions = spliceinfo->acceptor_positions_alloc;
      acceptorj_knowni = spliceinfo->acceptor_knowni_alloc;
    } else {
      acceptorj_nsites = spliceinfo->segmentj_acceptor_nknown;
      acceptorj_positions = spliceinfo->segmentj_acceptor_knownpos;
      acceptorj_knowni = spliceinfo->segmentj_acceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      } else {
	probj = Maxent_hr_acceptor_prob(segmentj_left + acceptorj_positions[i],segmentj_chroffset);
	printf(" (%.6f)",probj);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < donori_nsites && j < acceptorj_nsites) {
      supporti = donori_positions[i] - querystart;
      supportj = queryend - acceptorj_positions[j];

      if ((splice_qpos = donori_positions[i]) < acceptorj_positions[j]) {
	i++;
      } else if (splice_qpos > acceptorj_positions[j]) {
	j++;
      } else if (supporti < MIN_SUPPORT_SPLICE || supportj < MIN_SUPPORT_SPLICE) {
	/* Skip */
	i++; j++;
      } else {
	debug1(printf("splice matches at %d\n",splice_qpos));
	segmenti_nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								 /*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								 /*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	debug1(printf("%d mismatches on segmenti (%d..%d)\n",segmenti_nmismatches,querystart,splice_qpos));
	debug1(printf("%d mismatches on segmentj (%d..%d)\n",segmentj_nmismatches,splice_qpos,queryend));

	if (supporti - 3*segmenti_nmismatches < MIN_SUPPORT_SPLICE) {
	  /* Skip, because too many mismatches in segmenti */
	  debug1(printf("Skipping, because too many mismatches %d in segmenti\n",
			segmenti_nmismatches));
	} else if (supportj - 3*segmentj_nmismatches < MIN_SUPPORT_SPLICE) {
	  /* Skip, because too many mismatches in segmentj */
	  debug1(printf("Skipping, because too many mismatches %d in segmentj\n",
			segmentj_nmismatches));
	} else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) > best_nmismatches) {
	  debug1(printf("Skipping, because too many nmismatches %d > best_nmismatches %d\n",
			nmismatches,best_nmismatches));
	} else {
	  if (donori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	    /* sufficient1p = true; */
	  } else {
	    probi = Maxent_hr_donor_prob(segmenti_left + splice_qpos,segmenti_chroffset);
	    /* sufficient1p = sufficient_splice_prob_local(splice_qpos,segmenti_nmismatches,probi); */
	  }
	  
	  if (acceptorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	    /* sufficient2p = true; */
	  } else {
	    probj = Maxent_hr_acceptor_prob(segmentj_left + splice_qpos,segmentj_chroffset);
	    /* sufficient2p = sufficient_splice_prob_local(querylength - splice_qpos,segmentj_nmismatches,probj); */
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus sense splice_qpos  %d, i.donor %f, j.acceptor %f\n",splice_qpos,probi,probj);
		 } else {
		   printf("minus antisense splice_qpos  %d, i.donor %f, j.acceptor %f\n",splice_qpos,probi,probj);
		 });

	  if (nmismatches < best_nmismatches ||
	      (nmismatches == best_nmismatches && probi + probj > best_prob)) {
	    /* Success */
	    best_nmismatches = nmismatches;
	    best_prob = probi + probj;
	    
	    /* best_donor_splicecoord = segmenti_left + splice_qpos; */
	    /* best_acceptor_splicecoord = segmentj_left + splice_qpos; */
	    *best_knowni_i = donori_knowni[i];
	    *best_knowni_j = acceptorj_knowni[j];
	    *best_prob_i = probi; /* donor_prob */
	    *best_prob_j = probj; /* acceptor_prob */
	    best_splice_qpos = splice_qpos;
	    *best_nmismatches_i = segmenti_nmismatches;
	    *best_nmismatches_j = segmentj_nmismatches;
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
    if (novelsplicingp && segmenti_left + splice_qpos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(spliceinfo->acceptor_positions_alloc,spliceinfo->acceptor_knowni_alloc,
							   spliceinfo->segmenti_antiacceptor_knownpos,
							   spliceinfo->segmenti_antiacceptor_knowni,
							   segmenti_left,splice_qpos_start,splice_qpos_end);
      antiacceptori_positions = spliceinfo->acceptor_positions_alloc;
      antiacceptori_knowni = spliceinfo->acceptor_knowni_alloc;
    } else {
      antiacceptori_nsites = spliceinfo->segmenti_antiacceptor_nknown;
      antiacceptori_positions = spliceinfo->segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = spliceinfo->segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      } else {
	probi = Maxent_hr_antiacceptor_prob(segmenti_left + antiacceptori_positions[i],segmenti_chroffset);
	printf(" (%.6f)",probi);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_qpos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(spliceinfo->donor_positions_alloc,spliceinfo->donor_knowni_alloc,
						     spliceinfo->segmentj_antidonor_knownpos,
						     spliceinfo->segmentj_antidonor_knowni,
						     segmentj_left,splice_qpos_start,splice_qpos_end);
      antidonorj_positions = spliceinfo->donor_positions_alloc;
      antidonorj_knowni = spliceinfo->donor_knowni_alloc;
    } else {
      antidonorj_nsites = spliceinfo->segmentj_antidonor_nknown;
      antidonorj_positions = spliceinfo->segmentj_antidonor_knownpos;
      antidonorj_knowni = spliceinfo->segmentj_antidonor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      } else {
	probj = Maxent_hr_antidonor_prob(segmentj_left + antidonorj_positions[i],segmentj_chroffset);
	printf(" (%.6f)",probj);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;

    i = j = 0;
    while (i < antiacceptori_nsites && j < antidonorj_nsites) {
      supporti = antiacceptori_positions[i] - querystart;
      supportj = queryend - antidonorj_positions[j];

      if ((splice_qpos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	i++;
      } else if (splice_qpos > antidonorj_positions[j]) {
	j++;
      } else if (supporti < MIN_SUPPORT_SPLICE || supportj < MIN_SUPPORT_SPLICE) {
	/* Skip */
	i++; j++;
      } else {
	debug1(printf("splice matches at %d\n",splice_qpos));
	segmenti_nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								 /*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								 /*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	debug1(printf("%d mismatches on segmenti (%d..%d)\n",segmenti_nmismatches,querystart,splice_qpos));
	debug1(printf("%d mismatches on segmentj (%d..%d)\n",segmentj_nmismatches,splice_qpos,queryend));

	if (supporti - 3*segmenti_nmismatches < MIN_SUPPORT_SPLICE) {
	  /* Skip, because too many mismatches in segmenti */
	  debug1(printf("Skipping, because too many mismatches %d in segmenti\n",
			segmenti_nmismatches));
	} else if (supportj - 3*segmentj_nmismatches < MIN_SUPPORT_SPLICE) {
	  /* Skip, because too many mismatches in segmentj */
	  debug1(printf("Skipping, because too many mismatches %d in segmentj\n",
			segmentj_nmismatches));
	} else if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) > best_nmismatches) {
	  debug1(printf("Skipping, because too many nmismatches %d > best_nmismatches %d\n",
			nmismatches,best_nmismatches));
	} else {
	  if (antiacceptori_knowni[i] >= 0) {
	    probi = 1.0; /* Needs to be 1.0 for output */
	    /* sufficient1p = true; */
	  } else {
	    probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_qpos,segmenti_chroffset);
	    /* sufficient1p = sufficient_splice_prob_local(splice_qpos,segmenti_nmismatches,probi); */
	  }

	  if (antidonorj_knowni[j] >= 0) {
	    probj = 1.0; /* Needs to be 1.0 for output */
	    /* sufficient2p = true; */
	  } else {
	    probj = Maxent_hr_antidonor_prob(segmentj_left + splice_qpos,segmentj_chroffset);
	    /* sufficient2p = sufficient_splice_prob_local(querylength - splice_qpos,segmentj_nmismatches,probj); */
	  }

	  debug1(
		 if (plusp == true) {
		   printf("plus antisense splice_qpos  %d, j.donor %f, i.acceptor %f\n",splice_qpos,probj,probi);
		 } else {
		   printf("minus sense splice_qpos  %d, j.donor %f, i.acceptor %f\n",splice_qpos,probj,probi);
		 });
	  
	  if (nmismatches < best_nmismatches ||
	      (nmismatches == best_nmismatches && probi + probj > best_prob)) {
	    /* Success */
	    best_nmismatches = nmismatches;
	    best_prob = probi + probj;
	    
	    /* best_donor_splicecoord = segmentj_left + splice_qpos; */
	    /* best_acceptor_splicecoord = segmenti_left + splice_qpos; */
	    *best_knowni_j = antidonorj_knowni[j];
	    *best_knowni_i = antiacceptori_knowni[i];
	    *best_prob_j = probj; /* donor_prob */
	    *best_prob_i = probi; /* acceptor_prob */
	    best_splice_qpos = splice_qpos;
	    *best_nmismatches_j = segmentj_nmismatches;
	    *best_nmismatches_i = segmenti_nmismatches;
	  }
	}
	i++;
	j++;
      }
    }
  }

  debug1(printf("best_knowni_i is %d and best_knowni_j is %d\n",*best_knowni_i,*best_knowni_j));

  if (*best_prob_i >= 0.9 || *best_prob_j >= 0.9) {
    debug1(printf("Returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
    debug1(printf("nmismatches %d and %d\n",*best_nmismatches_i,*best_nmismatches_j));
    return best_splice_qpos;
  } else if (allow_indel_p == false) {
    debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
    *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
    return -1;
  } else {
    /* Find an indel below */
  }

  if (plusp == false) {
    /* Genomic sense (inferred from antisense code) */
    if (donori_nsites == 0 || acceptorj_nsites == 0) {
      debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
      *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
      return -1;
    } else {
      donori_probs = spliceinfo->donor_probs_alloc;
      for (i = 0; i < donori_nsites; i++) {
	donori_probs[i] = Maxent_hr_donor_prob(segmenti_left + donori_positions[i],segmenti_chroffset);
      }

      acceptorj_probs = spliceinfo->acceptor_probs_alloc;
      for (i = 0; i < acceptorj_nsites; i++) {
	acceptorj_probs[i] = Maxent_hr_acceptor_prob(segmentj_left + acceptorj_positions[i],segmentj_chroffset);
      }

      i = j = 0;
      best_nmismatches = max_mismatches_allowed + 1;
      while (i < donori_nsites && donori_positions[i] - querystart < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	i++;
      }

      while (i < donori_nsites) {
	supporti = donori_positions[i] - querystart;

	/* Backup */
	while (j >= 0 && acceptorj_positions[j] + max_deletionlen > donori_positions[i]) {
	  j--;
	}
	j++;			/* Finish backup */

	/* Advance */
	while (j < acceptorj_nsites && acceptorj_positions[j] + max_deletionlen <= donori_positions[i]) {
	  j++;
	}

	/* Deletions */
	while (j < acceptorj_nsites && acceptorj_positions[j] < donori_positions[i] &&
	       (supportj = queryend - acceptorj_positions[j]) >= MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	  debug1(printf("Deletion: %d (%.6f) -- %d (%.6f).  Support: %d and %d\n",
			donori_positions[i],donori_probs[i],acceptorj_positions[j],acceptorj_probs[j],
			donori_positions[i] - querystart,queryend - acceptorj_positions[j]));
	  probi = donori_probs[i];
	  probj = acceptorj_probs[j];
	  if (probi >= 0.85 && probj >= 0.85) {
	    nindels = acceptorj_positions[j] - donori_positions[i]; /* Should be negative */

	    /* Try deletion on segmenti */
	    splice_qpos = acceptorj_positions[j];
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,segmenti_left,nindels,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							   querystart,/*queryend*/splice_qpos,querylength,
							   max_mismatches_allowed,/*plusp:true*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > querystart);
	      assert(indel_pos < splice_qpos);
	      nmismatches_indel = nmismatches1;
	      nmismatches_i = nmismatches2;
	      nmismatches_j = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								      /*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	      debug1(printf("  => deletion indel_pos %d on segmenti with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_indel,nmismatches_i,nmismatches_j));

	      if (supporti - 3*(nmismatches_indel + nmismatches_i) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supportj - 3*nmismatches_j < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }

	    /* Try deletion on segmentj */
	    splice_qpos = donori_positions[i];
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,segmentj_left + nindels,nindels,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							   /*querystart*/splice_qpos,queryend,querylength,
							   max_mismatches_allowed,/*plusp:true*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos < queryend);
	      nmismatches_i = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								/*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	      nmismatches_j = nmismatches1;
	      nmismatches_indel = nmismatches2;
	      debug1(printf("  => deletion indel_pos %d on segmentj with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_i,nmismatches_j,nmismatches_indel));
	      
	      if (supportj - 3*(nmismatches_indel + nmismatches_j) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supporti - 3*nmismatches_i < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_i + nmismatches_j + nmismatches_indel) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }
	  }
	  j++;
	}

	if (j < acceptorj_nsites && acceptorj_positions[j] == donori_positions[i]) {
	  /* Not an indel */
	  j++;
	}

	/* Insertions */
	while (j < acceptorj_nsites && acceptorj_positions[j] <= donori_positions[i] + max_insertionlen &&
	       (supportj = queryend - acceptorj_positions[j]) >= MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	  debug1(printf("Insertion: %d (%.6f) -- %d (%.6f).  Support: %d and %d\n",
			donori_positions[i],donori_probs[i],acceptorj_positions[j],acceptorj_probs[j],
			donori_positions[i] - querystart,queryend - acceptorj_positions[j]));
	  probi = donori_probs[i];
	  probj = acceptorj_probs[j];
	  if (probi >= 0.85 && probj >= 0.85) {
	    nindels = acceptorj_positions[j] - donori_positions[i]; /* Should be positive */

	    /* Try insertion on segmenti */
	    splice_qpos = acceptorj_positions[j];
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,segmenti_left,nindels,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    querystart,/*queryend*/splice_qpos,querylength,
							    max_mismatches_allowed,/*plusp:true*/true,genestrand,
							    /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > querystart);
	      assert(indel_pos + nindels < splice_qpos);
	      nmismatches_indel = nmismatches1;
	      nmismatches_i = nmismatches2;
	      nmismatches_j = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								/*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	      debug1(printf("  => insertion indel_pos %d on segmenti with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_indel,nmismatches_i,nmismatches_j));

	      if (supporti - 3*(nmismatches_indel + nmismatches_i) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supportj - 3*nmismatches_j < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }

	    /* Try insertion on segmentj */
	    splice_qpos = donori_positions[i];
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,segmentj_left + nindels,nindels,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    /*querystart*/splice_qpos,queryend,querylength,
							    max_mismatches_allowed,/*plusp:true*/true,genestrand,
							    /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos + nindels < queryend);
	      nmismatches_i = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								/*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	      nmismatches_j = nmismatches1;
	      nmismatches_indel = nmismatches2;
	      debug1(printf("  => insertion indel_pos %d on segmentj with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_i,nmismatches_j,nmismatches_indel));
	      
	      if (supportj - 3*(nmismatches_indel + nmismatches_j) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supporti - 3*nmismatches_i < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_i + nmismatches_j + nmismatches_indel) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }
	  }
	  j++;
	}
	j--;			/* Finish advance */

	i++;
      }

      debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
      *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
      return -1;
    }

  } else {
    /* Genomic antisense (verified) */
    if (antiacceptori_nsites == 0 || antidonorj_nsites == 0) {
      debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
      *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
      return -1;
    } else {
      antiacceptori_probs = spliceinfo->acceptor_probs_alloc;
      for (i = 0; i < antiacceptori_nsites; i++) {
	antiacceptori_probs[i] = Maxent_hr_antiacceptor_prob(segmenti_left + antiacceptori_positions[i],segmenti_chroffset);
      }

      antidonorj_probs = spliceinfo->donor_probs_alloc;
      for (i = 0; i < antidonorj_nsites; i++) {
	antidonorj_probs[i] = Maxent_hr_antidonor_prob(segmentj_left + antidonorj_positions[i],segmentj_chroffset);
      }

      i = j = 0;
      best_nmismatches = max_mismatches_allowed + 1;
      while (i < antiacceptori_nsites && antiacceptori_positions[i] - querystart < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	i++;
      }

      while (i < antiacceptori_nsites) {
	supporti = antiacceptori_positions[i] - querystart;

	/* Backup */
	while (j >= 0 && antidonorj_positions[j] + max_deletionlen > antiacceptori_positions[i]) {
	  j--;
	}
	j++;			/* Finish backup */

	/* Advance */
	while (j < antidonorj_nsites && antidonorj_positions[j] + max_deletionlen <= antiacceptori_positions[i]) {
	  j++;
	}

	/* Deletions */
	while (j < antidonorj_nsites && antidonorj_positions[j] < antiacceptori_positions[i] &&
	       (supportj = queryend - antidonorj_positions[j]) >= MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	  debug1(printf("Deletion: %d (%.6f) -- %d (%.6f).  Support: %d and %d\n",
			antiacceptori_positions[i],antiacceptori_probs[i],antidonorj_positions[j],antidonorj_probs[j],
			antiacceptori_positions[i] - querystart,queryend - antidonorj_positions[j]));
	  probi = antiacceptori_probs[i];
	  probj = antidonorj_probs[j];
	  if (probi >= 0.85 && probj >= 0.85) {
	    nindels = antidonorj_positions[j] - antiacceptori_positions[i]; /* Should be negative */

	    /* Try deletion on segmenti */
	    splice_qpos = antidonorj_positions[j];
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,segmenti_left,nindels,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							   querystart,/*queryend*/splice_qpos,querylength,
							   max_mismatches_allowed,/*plusp:true*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > querystart);
	      assert(indel_pos < splice_qpos);
	      nmismatches_indel = nmismatches1;
	      nmismatches_i = nmismatches2;
	      nmismatches_j = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								/*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	      debug1(printf("  => deletion indel_pos %d on segmenti with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_indel,nmismatches_i,nmismatches_j));

	      if (supporti - 3*(nmismatches_indel + nmismatches_i) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supportj - 3*nmismatches_j < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }

	    /* Try deletion on segmentj */
	    splice_qpos = antiacceptori_positions[i];
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,segmentj_left + nindels,nindels,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							   /*querystart*/splice_qpos,queryend,querylength,
							   max_mismatches_allowed,/*plusp:true*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos < queryend);
	      nmismatches_i = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								/*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	      nmismatches_j = nmismatches1;
	      nmismatches_indel = nmismatches2;
	      debug1(printf("  => deletion indel_pos %d on segmentj with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_i,nmismatches_j,nmismatches_indel));
	      
	      if (supportj - 3*(nmismatches_indel + nmismatches_j) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supporti - 3*nmismatches_i < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_i + nmismatches_j + nmismatches_indel) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }
	  }
	  j++;
	}

	if (j < antidonorj_nsites && antidonorj_positions[j] == antiacceptori_positions[i]) {
	  /* Not an indel */
	  j++;
	}

	/* Insertions */
	while (j < antidonorj_nsites && antidonorj_positions[j] <= antiacceptori_positions[i] + max_insertionlen &&
	       (supportj = queryend - antidonorj_positions[j]) >= MIN_SUPPORT_SPLICE_PLUS_INDEL) {
	  debug1(printf("Insertion: %d (%.6f) -- %d (%.6f).  Support: %d and %d\n",
			antiacceptori_positions[i],antiacceptori_probs[i],antidonorj_positions[j],antidonorj_probs[j],
			antiacceptori_positions[i] - querystart,queryend - antidonorj_positions[j]));
	  probi = antiacceptori_probs[i];
	  probj = antidonorj_probs[j];
	  if (probi >= 0.85 && probj >= 0.85) {
	    nindels = antidonorj_positions[j] - antiacceptori_positions[i]; /* Should be positive */

	    /* Try insertion on segmenti */
	    splice_qpos = antidonorj_positions[j];
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,segmenti_left,nindels,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    querystart,/*queryend*/splice_qpos,querylength,
							    max_mismatches_allowed,/*plusp:true*/true,genestrand,
							    /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > querystart);
	      assert(indel_pos + nindels < splice_qpos);
	      nmismatches_indel = nmismatches1;
	      nmismatches_i = nmismatches2;
	      nmismatches_j = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmentj_left,
								/*pos5*/splice_qpos,/*pos3*/queryend,plusp,genestrand);
	      debug1(printf("  => insertion indel_pos %d on segmenti with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_indel,nmismatches_i,nmismatches_j));

	      if (supporti - 3*(nmismatches_indel + nmismatches_i) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supportj - 3*nmismatches_j < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }

	    /* Try insertion on segmentj */
	    splice_qpos = antiacceptori_positions[i];
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,segmentj_left + nindels,nindels,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    /*querystart*/splice_qpos,queryend,querylength,
							    max_mismatches_allowed,/*plusp:true*/true,genestrand,
							    /*want_lowest_coordinate_p*/true)) >= 0) {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos + nindels < queryend);
	      nmismatches_i = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,/*left*/segmenti_left,
								/*pos5*/querystart,/*pos3*/splice_qpos,plusp,genestrand);
	      nmismatches_j = nmismatches1;
	      nmismatches_indel = nmismatches2;
	      debug1(printf("  => insertion indel_pos %d on segmentj with mismatches %d+%d+%d.  ",
			    indel_pos,nmismatches_i,nmismatches_j,nmismatches_indel));
	      
	      if (supportj - 3*(nmismatches_indel + nmismatches_j) < MIN_SUPPORT_SPLICE_PLUS_INDEL || supporti - 3*nmismatches_i < MIN_SUPPORT_SPLICE_PLUS_INDEL) {
		debug1(printf("Not enough support after mismatches\n"));
	      } else if ((nmismatches = nmismatches_i + nmismatches_j + nmismatches_indel) < best_nmismatches) {
		debug1(printf("Enough support after mismatches\n"));
		*best_nindels = nindels;
		*best_indel_pos = indel_pos;
		*best_prob_i = probi;
		*best_prob_j = probj;
		*best_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_i = nmismatches_i;
		*best_nmismatches_j = nmismatches_j;
		best_splice_qpos = splice_qpos;
		best_nmismatches = nmismatches;
	      }
	    }
	  }
	  j++;
	}
	j--;			/* Finish advance */

	i++;
      }
    }
  }

  if (*best_nindels == 0) {
    debug1(printf("Not returning %d with probi %f and probj %f\n",best_splice_qpos,*best_prob_i,*best_prob_j));
    *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
    return -1;					    /* Could not find a splice plus indel */
  } else {
    return best_splice_qpos;
  }
}


static double
splice_prob_eval (int *sensedir, double *best_prob_i, double *best_prob_j, int splice_pos,
		  Univcoord_T segmenti_left, Univcoord_T segmentj_left,
		  Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
		  int querylength, bool plusp_i, bool plusp_j) {
  double prob, probi, probj;

  debug2(printf("\nplusp %d and %d\n",plusp_i,plusp_j));

  /* SENSE_FORWARD */
  if (plusp_i == true) {
    probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
  } else {
    probi = Maxent_hr_antidonor_prob(segmenti_left + querylength - splice_pos,segmenti_chroffset); /* correct */
  }
  if (plusp_j == true) {
    probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
  } else {
    probj = Maxent_hr_antiacceptor_prob(segmentj_left + querylength - splice_pos,segmentj_chroffset); /* correct */
  }
#if 0
  debug2(printf("sense forward: %f + %f at %u..%u\n",
		probi,probj,segmenti_left + querylength - splice_pos,segmentj_left + querylength - splice_pos));
#endif
  *best_prob_i = probi;
  *best_prob_j = probj;
  prob = probi + probj;
  
  /* SENSE_ANTI */
  if (plusp_i == true) {
    probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset); /* correct */
  } else {
    probi = Maxent_hr_acceptor_prob(segmenti_left + querylength - splice_pos,segmenti_chroffset);
  }
  if (plusp_j == true) {
    probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset); /* correct */
  } else {
    probj = Maxent_hr_donor_prob(segmentj_left + querylength - splice_pos,segmentj_chroffset);
  }
#if 0
  debug2(printf("sense anti: %f + %f at %u..%u\n",
		probi,probj,segmenti_left + splice_pos,segmentj_left + splice_pos));
#endif
  if (probi + probj > prob) {
    *best_prob_i = probi;
    *best_prob_j = probj;
    *sensedir = SENSE_ANTI;
    return probi + probj;
  } else {
    *sensedir = SENSE_FORWARD;
    return prob;
  }
}


int
Splice_resolve_distant (int *best_nmismatches_i, int *best_nmismatches_j,
			double *best_prob_i, double *best_prob_j, int *sensedir_distant_guess,

			int *mismatch_positions_i, int segmenti_nmismatches,
			int *mismatch_positions_j, int segmentj_nmismatches,

			Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,

			int querystart, int queryend,
			int splice_pos_start, int splice_pos_end, int querylength,
			bool plusp_i, bool plusp_j) {

  int best_splice_pos = -1, prev_splice_pos, next_splice_pos, low_splice_pos, high_splice_pos, splice_pos;
  int lowi, highi, lowj, highj, i, j;
  int supporti, supportj;
  int best_nmismatches, nmismatches, nmismatches_i, nmismatches_j;
  double best_prob, prob, prob_i, prob_j;
  int sensedir;
  bool donep;
#ifdef DEBUG2
  int k;
#endif


#ifdef DEBUG2
  printf("Entering Splice_resolve_distant with splice_pos_start %d and splice_pos_end %d\n",
	 splice_pos_start,splice_pos_end);

  printf("startfrag has %d mismatches at:",segmenti_nmismatches);
  for (k = 0; k <= segmenti_nmismatches; k++) {
    printf(" %d",mismatch_positions_i[k]);
  }
  printf("\n");

  printf("endfrag has %d mismatches at:",segmentj_nmismatches);
  for (k = 0; k <= segmentj_nmismatches; k++) {
    printf(" %d",mismatch_positions_j[k]);
  }
  printf("\n");
  printf("\n");
#endif


  *best_nmismatches_i = *best_nmismatches_j = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
  *best_prob_i = *best_prob_j = 0.0;

  /* Require separation from endpoints */
  /* May want a different min_shortend for distant */
  debug2(printf("splice_pos boundaries: %d..%d => ",splice_pos_start,splice_pos_end));
  if (splice_pos_start < querystart + min_shortend_distant) {
    splice_pos_start = querystart + min_shortend_distant;
  }
  if (splice_pos_end > queryend - min_shortend_distant) {
    splice_pos_end = queryend - min_shortend_distant;
  }
  debug2(printf("%d..%d\n",splice_pos_start,splice_pos_end));

  if (splice_pos_start >= splice_pos_end) {
    debug2(printf("Endpoints do not work\n"));
    return -1;
  }

  /* Initialize to splice_pos_start */
  nmismatches_i = 0;
  nmismatches_j = segmentj_nmismatches + 1;

  i = j = 0;
  lowi = highi = 0;
  lowj = highj = 0;
  while (i <= segmenti_nmismatches && mismatch_positions_i[i] < splice_pos_start) {
    lowi = highi = i;
    nmismatches_i++; i++;
  }
  while (j <= segmentj_nmismatches && mismatch_positions_j[j] < splice_pos_start) {
    lowj = highj = j;
    nmismatches_j--; j++;
  }
  best_nmismatches = querylength;
    

  donep = false;
  while (donep == false && i <= segmenti_nmismatches && j <= segmentj_nmismatches) {
    if ((splice_pos = mismatch_positions_i[i]) < mismatch_positions_j[j]) {
      debug2(printf("(1 i) Considering break pos prev through %d, which would give %d+%d mismatches\n",
		    splice_pos,nmismatches_i,nmismatches_j));
      if (splice_pos > splice_pos_end) {
	donep = true;
      } else if (nmismatches_i + nmismatches_j < best_nmismatches) {
	lowi = highi = i;
	lowj = highj = j;
	best_nmismatches = nmismatches_i + nmismatches_j;
      } else if (nmismatches_i + nmismatches_j == best_nmismatches) {
	highi = i;
	highj = j;
      }
      nmismatches_i++; i++;
	
    } else if ((splice_pos = mismatch_positions_j[j]) < mismatch_positions_i[i]) {
      debug2(printf("(1 j) Considering break pos prev through %d, which would give %d+%d mismatches\n",
		    splice_pos,nmismatches_i,nmismatches_j));
      if (splice_pos > splice_pos_end) {
	donep = true;
      } else if (nmismatches_i + nmismatches_j < best_nmismatches) {
	lowi = highi = i;
	lowj = highj = j;
	best_nmismatches = nmismatches_i + nmismatches_j;
      } else if (nmismatches_i + nmismatches_j == best_nmismatches) {
	highi = i;
	highj = j;
      }
      nmismatches_j--; j++;
	
    } else {
      debug2(printf("(1 equal) Considering break pos prev through %d, which would give %d+%d mismatches\n",
		    splice_pos,nmismatches_i,nmismatches_j));
      if (splice_pos > splice_pos_end) {
	donep = true;
      } else if (nmismatches_i + nmismatches_j < best_nmismatches) {
	lowi = highi = i;
	lowj = highj = j;
	best_nmismatches = nmismatches_i + nmismatches_j;
      } else if (nmismatches_i + nmismatches_j == best_nmismatches) {
	highi = i;
	highj = j;
      }
      nmismatches_i++; i++;
      nmismatches_j--; j++;
    }
  }


  if (i > segmenti_nmismatches && j > segmentj_nmismatches) {
    debug2(printf("Reached both ends\n"));
    
  } else if (j > segmentj_nmismatches) {
    while (donep == false && i <= segmenti_nmismatches) {
      debug2(printf("(2 i) Considering break pos prev through %d, which would give %d+%d mismatches\n",
		    mismatch_positions_i[i],nmismatches_i,nmismatches_j));
      if ((splice_pos = mismatch_positions_i[i]) > splice_pos_end) {
	donep = true;
      } else if (nmismatches_i + nmismatches_j < best_nmismatches) {
	lowi = highi = i;
	best_nmismatches = nmismatches_i + nmismatches_j;	
      } else if (nmismatches_i + nmismatches_j == best_nmismatches) {
	highi = i;
      }
      nmismatches_i++; i++;
    }
    
  } else if (i > segmenti_nmismatches) {
    while (donep == false && j <= segmentj_nmismatches) {
      debug2(printf("(2 j) Considering break pos prev through %d, which would give %d+%d mismatches\n",
		    mismatch_positions_j[j],nmismatches_i,nmismatches_j));
      if ((splice_pos = mismatch_positions_j[j]) > splice_pos_end) {
	donep = true;
      } else if (nmismatches_i + nmismatches_j < best_nmismatches) {
	lowj = highj = j;
	best_nmismatches = nmismatches_i + nmismatches_j;
      } else if (nmismatches_i + nmismatches_j == best_nmismatches) {
	highj = j;
      }
      nmismatches_j--; j++;
    }
  }
  debug2(printf("i range: %d..%d.  j range: %d..%d\n",lowi,highi,lowj,highj));
  
  if (mismatch_positions_i[lowi] < mismatch_positions_j[lowj]) {
    low_splice_pos = mismatch_positions_i[lowi];
  } else {
    low_splice_pos = mismatch_positions_j[lowj];
  }

  if (mismatch_positions_i[highi] > mismatch_positions_j[highj]) {
    high_splice_pos = mismatch_positions_i[highi];
  } else {
    high_splice_pos = mismatch_positions_j[highj];
  }
  debug2(printf("(low_splice_pos %d, high_splice_pos %d]\n",low_splice_pos,high_splice_pos));

  if (novelsplicingp == false) {
    *best_prob_i = *best_prob_j = 0.0;
    *sensedir_distant_guess = SENSE_NULL;
    return (low_splice_pos + high_splice_pos)/2;
  }


  /* Find choices outside of this range */
  nmismatches = best_nmismatches;
  prev_splice_pos = low_splice_pos;
  while ((lowi > 0 || lowj > 0) && nmismatches < best_nmismatches + SLOP) {
    if (lowi == 0) {
      prev_splice_pos = mismatch_positions_j[--lowj];
    } else if (lowj == 0) {
      prev_splice_pos = mismatch_positions_i[--lowi];
    } else if (mismatch_positions_i[lowi] > mismatch_positions_j[lowj]) {
      prev_splice_pos = mismatch_positions_i[--lowi];
    } else if (mismatch_positions_j[lowj] > mismatch_positions_i[lowi]) {
      prev_splice_pos = mismatch_positions_j[--lowj];
    } else {
      prev_splice_pos = mismatch_positions_i[--lowi];
      --lowj;
    }

    debug2(printf("Going back to prev_splice_pos %d, with %d mismatches vs best %d\n",
		  prev_splice_pos,nmismatches,best_nmismatches));

    if (prev_splice_pos == mismatch_positions_i[lowi]) {
      nmismatches--;
    }
    if (prev_splice_pos == mismatch_positions_j[lowj]) {
      nmismatches++;
    }
  }

  assert(highi <= segmenti_nmismatches);
  assert(highj <= segmentj_nmismatches);
  nmismatches = best_nmismatches;
  next_splice_pos = high_splice_pos;
  while ((highi < segmenti_nmismatches || highj < segmentj_nmismatches) && nmismatches < best_nmismatches + SLOP) {
    if (highi == segmenti_nmismatches) {
      next_splice_pos = mismatch_positions_j[++highj];
    } else if (highj == segmentj_nmismatches) {
      next_splice_pos = mismatch_positions_i[++highi];
    } else if (mismatch_positions_i[highi] < mismatch_positions_j[highj]) {
      next_splice_pos = mismatch_positions_i[++highi];
    } else if (mismatch_positions_j[highj] < mismatch_positions_i[highi]) {
      next_splice_pos = mismatch_positions_j[++highj];
    } else {
      next_splice_pos = mismatch_positions_i[++highi];
      ++highj;
    }

    debug2(printf("Going forward to next_splice_pos %d, with %d mismatches vs best %d\n",
		  next_splice_pos,nmismatches,best_nmismatches));
    if (next_splice_pos == mismatch_positions_i[highi]) {
      nmismatches++;
    }
    if (next_splice_pos == mismatch_positions_j[highj]) {
      nmismatches--;
    }
  }

  debug2(printf("(prev_splice_pos %d, next_splice_pos %d]\n",
		prev_splice_pos,next_splice_pos));

  if (prev_splice_pos < splice_pos_start) {
    prev_splice_pos = splice_pos_start;
  }
  if (next_splice_pos > splice_pos_end) {
    next_splice_pos = splice_pos_end;
  }
  debug2(printf("After bounds, (prev_splice_pos %d, next_splice_pos %d]\n",
		prev_splice_pos,next_splice_pos));


  /* Check for best splice probabilities from prev_best to next_best */
  nmismatches_i = 0;
  nmismatches_j = segmentj_nmismatches + 1;

  i = j = 0;
  lowi = highi = 0;
  lowj = highj = 0;
  while (i <= segmenti_nmismatches && mismatch_positions_i[i] <= prev_splice_pos) {
    lowi = highi = i;
    nmismatches_i++; i++;
  }
  while (j <= segmentj_nmismatches && mismatch_positions_j[j] <= prev_splice_pos) {
    lowj = highj = j;
    nmismatches_j--; j++;
  }

  best_prob = 0.0;
  for (splice_pos = prev_splice_pos + 1; splice_pos <= next_splice_pos; splice_pos++) {
    supporti = splice_pos - querystart;
    supportj = queryend - splice_pos;

    if (nmismatches_i >= supporti/15) {
      /* Skip, because too many mismatches in segmenti */
      debug2(printf("Skipping, because too many mismatches %d in segmenti\n",nmismatches_i));

    } else if (nmismatches_j >= supportj/15) {
      /* Skip, because too many mismatches in segmentj */
      debug2(printf("Skipping, because too many mismatches %d in segmentj\n",nmismatches_j));

    } else if ((prob = splice_prob_eval(&sensedir,&prob_i,&prob_j,splice_pos,
					segmenti_left,segmentj_left,segmenti_chroffset,segmentj_chroffset,
					querylength,plusp_i,plusp_j)) > best_prob) {
      *best_prob_i = prob_i;
      *best_prob_j = prob_j;
      *best_nmismatches_i = nmismatches_i;
      *best_nmismatches_j = nmismatches_j;
      *sensedir_distant_guess = sensedir;
      best_splice_pos = splice_pos;
      best_prob = prob;
    }

    debug2(printf("splice_pos %d, supporti %d and %d, mismatches %d+%d, sensedir %d, prob %f + %f\n",
		  splice_pos,supporti,supportj,nmismatches_i,nmismatches_j,sensedir,prob_i,prob_j));
    
    if (i <= segmenti_nmismatches && mismatch_positions_i[i] == splice_pos) {
      nmismatches_i++; i++;
    }
    if (j <= segmentj_nmismatches && mismatch_positions_j[j] == splice_pos) {
      nmismatches_j--; j++;
    }
  }

  debug2(printf("Returning best break pos at %d, probs %f + %f\n",best_splice_pos,*best_prob_i,*best_prob_j));
  return best_splice_pos;
}



#if 0
static int
donor_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length(Stage3end_substring_donor(x));
  int y_length = Substring_match_length(Stage3end_substring_donor(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    return 0;
  }
}
#endif


#if 0
static int
acceptor_match_length_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  int x_length = Substring_match_length(Stage3end_substring_acceptor(x));
  int y_length = Substring_match_length(Stage3end_substring_acceptor(y));

  if (x_length < y_length) {
    return -1;
  } else if (y_length < x_length) {
    return +1;
  } else {
    return 0;
  }
}
#endif


#if 0
static List_T
group_by_segmenti_aux (int *found_score, List_T winners, List_T *ambiguous,
		       Stage3end_T *hitarray, int n, int querylength, bool first_read_p,
		       Method_T method, int level) {
  Stage3end_T hit, *subarray;
  int i, j, k, ii, jj, kk, nn;
  int n_good_spliceends;
  Univcoord_T segmenti_left;
  Substring_T donor, acceptor;
  int best_nmismatches, nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob, donor_prob, acceptor_prob;
  List_T accepted_hits, rejected_hits, donor_hits, acceptor_hits, p;

  int sensedir;
  Univcoordlist_T ambcoords;
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
		      Stage3end_distance(hit),Substring_match_length(Stage3end_substring_donor(hit)),
		      Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
		      Substring_match_length(Stage3end_substring_acceptor(hit)),
		      Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
		      Stage3end_nmismatches_bothdiff(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
		      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	if ((nmismatches = Stage3end_nmismatches_bothdiff(hit)) < best_nmismatches) {
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
	if (Stage3end_nmismatches_bothdiff(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
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
	  if (Stage3end_nmismatches_bothdiff(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
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
	    donor_length = Substring_match_length(donor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length(Stage3end_substring_donor(subarray[jj])) == donor_length) {
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
		ambcoords = Univcoordlist_push(ambcoords,Substring_splicecoord_A(acceptor));
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_bothdiff(acceptor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
	      }

	      nmismatches_acceptor = best_nmismatches - Substring_nmismatches_bothdiff(donor);
	      donor_prob = Junction_donor_prob(Stage3end_junctionA(hit));
	      prob = best_prob - donor_prob;
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   /*nmismatches_donor*/Substring_nmismatches_bothdiff(donor),nmismatches_acceptor,
								   donor,/*acceptor*/NULL,donor_prob,/*acceptor_prob*/prob,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*ambcoords_donor*/NULL,ambcoords,
								   /*amb_knowni_donor*/NULL,amb_knowni,
								   /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								   /*amb_probs_donor*/NULL,amb_probs,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								   sensedir,method,level));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
	      Univcoordlist_free(&ambcoords);
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
	    acceptor_length = Substring_match_length(acceptor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length(Stage3end_substring_acceptor(subarray[jj])) == acceptor_length) {
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
		ambcoords = Univcoordlist_push(ambcoords,Substring_splicecoord_D(donor));
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_bothdiff(donor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
	      }

	      nmismatches_donor = best_nmismatches - Substring_nmismatches_bothdiff(acceptor);
	      acceptor_prob = Junction_acceptor_prob(Stage3end_junctionD(hit));
	      prob = best_prob - acceptor_prob;
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_bothdiff(acceptor),
								   /*donor*/NULL,acceptor,/*donor_prob*/prob,acceptor_prob,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   ambcoords,/*ambcoords_acceptor*/NULL,
								   amb_knowni,/*amb_knowni_acceptor*/NULL,
								   amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								   amb_probs,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								   sensedir,method,level));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
	      Univcoordlist_free(&ambcoords);
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
#endif


#if 0
List_T
Splice_group_by_segmenti (int *found_score, List_T localsplicing, List_T *ambiguous,
			  int querylength, bool first_read_p, Method_T method, int level) {
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
				    querylength,first_read_p,method,level);
    FREEA(array_forward);
  }

  if (n_sense_anti > 0) {
    qsort(array_anti,n_sense_anti,sizeof(Stage3end_T),Stage3end_chimera_segmenti_cmp);
    winners = group_by_segmenti_aux(&(*found_score),winners,&(*ambiguous),array_anti,n_sense_anti,
				    querylength,first_read_p,method,level);
    FREEA(array_anti);
  }

  List_free(&localsplicing);

  return winners;
}
#endif


#if 0
static List_T
group_by_segmentj_aux (int *found_score, List_T winners, List_T *ambiguous, 
		       Stage3end_T *hitarray, int n, int querylength, bool first_read_p,
		       Method_T method, int level) {
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
  Univcoordlist_T ambcoords;
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
		      Stage3end_distance(hit),Substring_match_length(Stage3end_substring_donor(hit)),
		      Substring_genomicstart(Stage3end_substring_donor(hit)),Substring_genomicend(Stage3end_substring_donor(hit)),
		      Substring_match_length(Stage3end_substring_acceptor(hit)),
		      Substring_genomicstart(Stage3end_substring_acceptor(hit)),Substring_genomicend(Stage3end_substring_acceptor(hit)),
		      Stage3end_nmismatches_bothdiff(hit),Substring_siteD_prob(Stage3end_substring_donor(hit)),
		      Substring_siteA_prob(Stage3end_substring_acceptor(hit))));
	if ((nmismatches = Stage3end_nmismatches_bothdiff(hit)) < best_nmismatches) {
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
	if (Stage3end_nmismatches_bothdiff(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP &&
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
	  if (Stage3end_nmismatches_bothdiff(hit) <= best_nmismatches + LOCALSPLICING_NMATCHES_SLOP ||
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
	    donor_length = Substring_match_length(donor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length(Stage3end_substring_donor(subarray[jj])) == donor_length) {
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
		ambcoords = Univcoordlist_push(ambcoords,Substring_splicecoord_A(acceptor));
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_bothdiff(acceptor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteA_prob(acceptor));
	      }

	      nmismatches_acceptor = best_nmismatches - Substring_nmismatches_bothdiff(donor);
	      donor_prob = Junction_donor_prob(Stage3end_junctionA(hit));
	      prob = best_prob - donor_prob;
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   /*nmismatches_donor*/Substring_nmismatches_bothdiff(donor),nmismatches_acceptor,
								   donor,/*acceptor*/NULL,donor_prob,/*acceptor_prob*/prob,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   /*ambcoords_donor*/NULL,ambcoords,
								   /*amb_knowni_donor*/NULL,amb_knowni,
								   /*amb_nmismatches_donor*/NULL,amb_nmismatches,
								   /*amb_probs_donor*/NULL,amb_probs,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								   sensedir,method,level));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
	      Univcoordlist_free(&ambcoords);
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
	    acceptor_length = Substring_match_length(acceptor);
	    jj = ii + 1;
	    while (jj < nn && Substring_match_length(Stage3end_substring_acceptor(subarray[jj])) == acceptor_length) {
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
		ambcoords = Univcoordlist_push(ambcoords,Substring_splicecoord_D(donor));
		amb_knowni = Intlist_push(amb_knowni,-1);
		amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_bothdiff(donor));
		amb_probs = Doublelist_push(amb_probs,Substring_siteD_prob(donor));
	      }

	      nmismatches_donor = best_nmismatches - Substring_nmismatches_bothdiff(acceptor);
	      acceptor_prob = Junction_acceptor_prob(Stage3end_junctionD(hit));
	      prob = best_prob - acceptor_prob;
	      *ambiguous = List_push(*ambiguous,
				     (void *) Stage3end_new_splice(&(*found_score),
								   nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_bothdiff(acceptor),
								   /*donor*/NULL,acceptor,/*donor_prob*/prob,acceptor_prob,/*distance*/0U,
								   /*shortdistancep*/false,/*penalty*/0,querylength,
								   ambcoords,/*ambcoords_acceptor*/NULL,
								   amb_knowni,/*amb_knowni_acceptor*/NULL,
								   amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
								   amb_probs,/*amb_probs_acceptor*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								   sensedir,method,level));
	      Doublelist_free(&amb_probs);
	      Intlist_free(&amb_knowni);
	      Intlist_free(&amb_nmismatches);
	      Univcoordlist_free(&ambcoords);
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
#endif


#if 0
List_T
Splice_group_by_segmentj (int *found_score, List_T localsplicing, List_T *ambiguous,
			  int querylength, bool first_read_p, Method_T method, int level) {
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
				    querylength,first_read_p,method,level);
    FREEA(array_forward);
  }

  if (n_sense_anti > 0) {
    qsort(array_anti,n_sense_anti,sizeof(Stage3end_T),Stage3end_chimera_segmentj_cmp);
    winners = group_by_segmentj_aux(&(*found_score),winners,&(*ambiguous),array_anti,n_sense_anti,
				    querylength,first_read_p,method,level);
    FREEA(array_anti);
  }

  List_free(&localsplicing);

  return winners;
}
#endif



#define END_SPLICESITE_PROB_MATCH 0.90

/* Derived from substring_trim_novel_spliceends in substring.c, which
   was modified from trim_novel_spliceends in stage3.c */
/* Note: If substring does not extend to ends of query, then region
   beyond querystart and queryend might actually be matching, and not
   mismatches.  Could fix in the future. */

/* TODO: Consider whether there is a conflict between the regular prob
   and mismatch prob, and if so, set splicedir to be SENSE_NULL */
int
Splice_trim_novel_spliceends_5 (Splicetype_T *splicetype, int **ambig_qstarts, double **ambig_probs_5,
				Univcoord_T left, int qstart, int qend,
				int *mismatch_positions, int nmismatches,
				Univcoord_T chroffset, bool plusp, int sensedir) {
  int nspliceends = 0;
  Univcoord_T start_genomicpos, middle_genomicpos, end_genomicpos, genomicpos;
  double donor_prob, acceptor_prob;
  int querypos;
  int nconsecutive, mismatchi;


  debug13(printf("\nEntered Splice_trim_novel_spliceends_5 with sensedir %d\n",sensedir));
  assert(sensedir != SENSE_NULL);

  middle_genomicpos = left + qstart;

  if (middle_genomicpos < left + END_SPLICESITE_SEARCH_MM) {
    start_genomicpos = left;
  } else {
    start_genomicpos = middle_genomicpos - END_SPLICESITE_SEARCH_MM;
  }

  if ((end_genomicpos = middle_genomicpos + END_SPLICESITE_SEARCH) > left + qend) {
    end_genomicpos = left + qend;
  }

  if (left + qend < MIN_EXON_LENGTH) {
    /* Skip */
  } else if (end_genomicpos < left + qend - MIN_EXON_LENGTH) {
    end_genomicpos = left + qend - MIN_EXON_LENGTH;
  }

  debug13(printf("\n1 Set end points for 5' trim to be %u..%u..%u\n",
		 start_genomicpos - chroffset,middle_genomicpos - chroffset,end_genomicpos - chroffset));

  if (sensedir == SENSE_FORWARD) {
    if (plusp) {
      debug13(printf("Case 2\n"));
      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos < middle_genomicpos) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	debug13(printf("5', watson, sense forward %u %u %d %f mm\n",genomicpos,genomicpos-chroffset,genomicpos-left,acceptor_prob));
	if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	genomicpos++;
      }
      
      querypos = genomicpos - left;
      mismatchi = 0;
      while (mismatchi < nmismatches && mismatch_positions[mismatchi] < querypos) {
	mismatchi++;
      }
      nconsecutive = 0;
      while (genomicpos < end_genomicpos && (nconsecutive < MAX_NCONSECUTIVE || nspliceends == 0)) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	debug13(printf("5', watson, sense forward %u %u %d %f.  nconsecutive %d\n",
		       genomicpos,genomicpos-chroffset,querypos,acceptor_prob,nconsecutive));
	if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	if (mismatchi < nmismatches && querypos == mismatch_positions[mismatchi]) {
	  nconsecutive = 0;
	  mismatchi++;
	} else {
	  nconsecutive++;
	}
	querypos++;
	genomicpos++;
      }
      end_genomicpos = genomicpos;
      debug13(printf("\n"));

    } else {
      debug13(printf("Case 6\n"));
      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos < middle_genomicpos) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	debug13(printf("5', watson, sense anti %u %u %d %f mm\n",genomicpos,genomicpos-chroffset,genomicpos-left,donor_prob));
	if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	genomicpos++;
      }
      
      querypos = genomicpos - left;
      mismatchi = 0;
      while (mismatchi < nmismatches && mismatch_positions[mismatchi] < querypos) {
	mismatchi++;
      }
      nconsecutive = 0;
      while (genomicpos < end_genomicpos && (nconsecutive < MAX_NCONSECUTIVE || nspliceends == 0)) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	debug13(printf("5', watson, sense anti %u %u %d %f.  nconsecutive %d\n",
		       genomicpos,genomicpos-chroffset,querypos,donor_prob,nconsecutive));
	if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	if (mismatchi < nmismatches && querypos == mismatch_positions[mismatchi]) {
	  nconsecutive = 0;
	  mismatchi++;
	} else {
	  nconsecutive++;
	}
	querypos++;
	genomicpos++;
      }
      end_genomicpos = genomicpos;
      debug13(printf("\n"));
    }

  } else {
    /* SENSE_ANTI */
    if (plusp) {
      debug13(printf("Case 6\n"));
      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos < middle_genomicpos) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	debug13(printf("5', watson, sense anti %u %u %d %f mm\n",genomicpos,genomicpos-chroffset,genomicpos-left,donor_prob));
	if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	genomicpos++;
      }
      
      querypos = genomicpos - left;
      mismatchi = 0;
      while (mismatchi < nmismatches && mismatch_positions[mismatchi] < querypos) {
	mismatchi++;
      }
      nconsecutive = 0;
      while (genomicpos < end_genomicpos && (nconsecutive < MAX_NCONSECUTIVE || nspliceends == 0)) {
	donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	debug13(printf("5', watson, sense anti %u %u %d %f.  nconsecutive %d\n",
		       genomicpos,genomicpos-chroffset,genomicpos-left,donor_prob,nconsecutive));
	if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	if (mismatchi < nmismatches && querypos == mismatch_positions[mismatchi]) {
	  nconsecutive = 0;
	  mismatchi++;
	} else {
	  nconsecutive++;
	}
	querypos++;
	genomicpos++;
      }
      end_genomicpos = genomicpos;
      debug13(printf("\n"));

    } else {
      debug13(printf("Case 2\n"));
      /* assert(start_genomicpos <= end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos < middle_genomicpos) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	debug13(printf("5', watson, sense forward %u %u %d %f mm\n",genomicpos,genomicpos-chroffset,genomicpos-left,acceptor_prob));
	if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	genomicpos++;
      }
      
      querypos = genomicpos - left;
      mismatchi = 0;
      while (mismatchi < nmismatches && mismatch_positions[mismatchi] < querypos) {
	mismatchi++;
      }
      nconsecutive = 0;
      while (genomicpos < end_genomicpos && (nconsecutive < MAX_NCONSECUTIVE || nspliceends == 0)) {
	acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	debug13(printf("5', watson, sense forward %u %u %d %f.  nconsecutive %d\n",
		       genomicpos,genomicpos-chroffset,querypos,acceptor_prob,nconsecutive));
	if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	if (mismatchi < nmismatches && querypos == mismatch_positions[mismatchi]) {
	  nconsecutive = 0;
	  mismatchi++;
	} else {
	  nconsecutive++;
	}
	querypos++;
	genomicpos++;
      }
      end_genomicpos = genomicpos;
      debug13(printf("\n"));
    }
  }

  if (nspliceends == 0) {
    debug13(printf("Got no spliceends\n"));
    return 0;
  } else {
    debug13(printf("Going from genomicpos %u up to %u\n",start_genomicpos,end_genomicpos));
    *ambig_qstarts = (int *) MALLOC(nspliceends*sizeof(int));
    *ambig_probs_5 = (double *) MALLOC(nspliceends*sizeof(double));
    nspliceends = 0;

    if (sensedir == SENSE_FORWARD) {
      if (plusp) {
	*splicetype = ACCEPTOR;
	genomicpos = start_genomicpos;
	while (genomicpos < end_genomicpos) {
	  acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	  if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	    (*ambig_qstarts)[nspliceends] = (int) (genomicpos - left);
	    (*ambig_probs_5)[nspliceends] = acceptor_prob;
	    nspliceends += 1;
	  }
	  genomicpos++;
	}
      } else {
	*splicetype = ANTIDONOR;

	genomicpos = start_genomicpos;
	while (genomicpos < end_genomicpos) {
	  donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	  if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	    (*ambig_qstarts)[nspliceends] = (int) (genomicpos - left);
	    (*ambig_probs_5)[nspliceends] = donor_prob;
	    nspliceends += 1;
	  }
	  genomicpos++;
	}
      }

    } else {
      /* SENSE_ANTI */
      if (plusp) {
	*splicetype = ANTIDONOR;

	genomicpos = start_genomicpos;
	while (genomicpos < end_genomicpos) {
	  donor_prob = Maxent_hr_antidonor_prob(genomicpos,chroffset); /* Case 6 */
	  if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	    (*ambig_qstarts)[nspliceends] = (int) (genomicpos - left);
	    (*ambig_probs_5)[nspliceends] = donor_prob;
	    nspliceends += 1;
	  }
	  genomicpos++;
	}
      } else {
	*splicetype = ACCEPTOR;
	
	genomicpos = start_genomicpos;
	while (genomicpos < end_genomicpos) {
	  acceptor_prob = Maxent_hr_acceptor_prob(genomicpos,chroffset); /* Case 2 */
	  if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	    (*ambig_qstarts)[nspliceends] = (int) (genomicpos - left);
	    (*ambig_probs_5)[nspliceends] = acceptor_prob;
	    nspliceends += 1;
	  }
	  genomicpos++;
	}
      }
    }

    debug13(printf("Got %d spliceends\n",nspliceends));
    return nspliceends;
  }
}



int
Splice_trim_novel_spliceends_3 (Splicetype_T *splicetype, int **ambig_qends, double **ambig_probs_3,
				Univcoord_T left, int qstart, int qend, int querylength,
				int *mismatch_positions, int nmismatches,
				Univcoord_T chroffset, bool plusp, int sensedir) {

  int nspliceends = 0;
  Univcoord_T start_genomicpos, middle_genomicpos, end_genomicpos, genomicpos;
  double donor_prob, acceptor_prob;
  int querypos;
  int nconsecutive, mismatchi;


  debug13(printf("\nEntered Splice_trim_novel_spliceends_3 with sensedir %d\n",sensedir));
  assert(sensedir != SENSE_NULL);

  middle_genomicpos = left + qend;

  if ((start_genomicpos = middle_genomicpos + END_SPLICESITE_SEARCH_MM) > left + querylength) {
    start_genomicpos = left + querylength;
  }

  if (middle_genomicpos < left + qstart + END_SPLICESITE_SEARCH) {
    end_genomicpos = left + qstart;
  } else {
    end_genomicpos = middle_genomicpos - END_SPLICESITE_SEARCH;
  }

  if (end_genomicpos < left + qstart + MIN_EXON_LENGTH) {
    end_genomicpos = left + qstart + MIN_EXON_LENGTH;
  }

  debug13(printf("\n1 Set end points for 3' trim to be %u..%u..%u\n",
		 start_genomicpos - chroffset,middle_genomicpos - chroffset,end_genomicpos - chroffset));

  if (sensedir == SENSE_FORWARD) {
    if (plusp) {
      debug13(printf("Case 1\n"));
      /* assert(start_genomicpos > end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos > middle_genomicpos) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	debug13(printf("3', watson, sense anti %u %u %d %f mm\n",genomicpos,genomicpos-chroffset,genomicpos - left,donor_prob));
	if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	genomicpos--;
      }
      
      querypos = genomicpos - left;
      mismatchi = 0;
      while (mismatchi < nmismatches && mismatch_positions[mismatchi] > querypos) {
	mismatchi++;
      }
      nconsecutive = 0;
      while (genomicpos > end_genomicpos && (nconsecutive < MAX_NCONSECUTIVE || nspliceends == 0)) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	debug13(printf("3', watson, sense anti %u %u %d %f.  nconsecutive %d\n",
		       genomicpos,genomicpos-chroffset,querypos,donor_prob,nconsecutive));
	if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	if (mismatchi < nmismatches && querypos == mismatch_positions[mismatchi]) {
	  nconsecutive = 0;
	  mismatchi++;
	} else {
	  nconsecutive++;
	}
	querypos--;
	genomicpos--;
      }
      end_genomicpos = genomicpos;
      debug13(printf("\n"));

    } else {
      debug13(printf("Case 5\n"));
      /* assert(start_genomicpos > end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos > middle_genomicpos) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	debug13(printf("3', watson, sense forward %u %u %d %f mm\n",genomicpos,genomicpos-chroffset,genomicpos - left,acceptor_prob));
	if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	genomicpos--;
      }
      
      querypos = genomicpos - left;
      mismatchi = 0;
      while (mismatchi < nmismatches && mismatch_positions[mismatchi] > querypos) {
	mismatchi++;
      }
      nconsecutive = 0;
      while (genomicpos > end_genomicpos && (nconsecutive < MAX_NCONSECUTIVE || nspliceends == 0)) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	debug13(printf("3', watson, sense forward %u %u %d %f.  nconsecutive %d\n",
		       genomicpos,genomicpos-chroffset,querypos,acceptor_prob,nconsecutive));
	if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	if (mismatchi < nmismatches && querypos == mismatch_positions[mismatchi]) {
	  nconsecutive = 0;
	  mismatchi++;
	} else {
	  nconsecutive++;
	}
	querypos--;
	genomicpos--;
      }
      end_genomicpos = genomicpos;
      debug13(printf("\n"));
    }

  } else {
    /* SENSE_ANTI */
    if (plusp) {
      debug13(printf("Case 5\n"));
      /* assert(start_genomicpos > end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos > middle_genomicpos) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	debug13(printf("3', watson, sense forward %u %u %d %f mm\n",genomicpos,genomicpos-chroffset,genomicpos - left,acceptor_prob));
	if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	genomicpos--;
      }
      
      querypos = genomicpos - left;
      mismatchi = 0;
      while (mismatchi < nmismatches && mismatch_positions[mismatchi] > querypos) {
	mismatchi++;
      }
      nconsecutive = 0;
      while (genomicpos > end_genomicpos && (nconsecutive < MAX_NCONSECUTIVE || nspliceends == 0)) {
	acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	debug13(printf("3', watson, sense forward %u %u %d %f.  nconsecutive %d\n",
		       genomicpos,genomicpos-chroffset,querypos,acceptor_prob,nconsecutive));
	if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	if (mismatchi < nmismatches && querypos == mismatch_positions[mismatchi]) {
	  nconsecutive = 0;
	  mismatchi++;
	} else {
	  nconsecutive++;
	}
	querypos--;
	genomicpos--;
      }
      end_genomicpos = genomicpos;
      debug13(printf("\n"));

    } else {
      debug13(printf("Case 1\n"));
      /* assert(start_genomicpos > end_genomicpos); */
      genomicpos = start_genomicpos;
      while (genomicpos > middle_genomicpos) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	debug13(printf("3', watson, sense anti %u %u %d %f mm\n",genomicpos,genomicpos-chroffset,genomicpos - left,donor_prob));
	if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	genomicpos--;
      }
      
      querypos = genomicpos - left;
      mismatchi = 0;
      while (mismatchi < nmismatches && mismatch_positions[mismatchi] > querypos) {
	mismatchi++;
      }
      nconsecutive = 0;
      while (genomicpos > end_genomicpos && (nconsecutive < MAX_NCONSECUTIVE || nspliceends == 0)) {
	donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	debug13(printf("3', watson, sense anti %u %u %d %f.  nconsecutive %d\n",
		       genomicpos,genomicpos-chroffset,querypos,donor_prob,nconsecutive));
	if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	  nspliceends += 1;
	}
	if (mismatchi < nmismatches && querypos == mismatch_positions[mismatchi]) {
	  nconsecutive = 0;
	  mismatchi++;
	} else {
	  nconsecutive++;
	}
	querypos--;
	genomicpos--;
      }
      end_genomicpos = genomicpos;
      debug13(printf("\n"));
    }
  }

  if (nspliceends == 0) {
    debug13(printf("Got no spliceends\n"));
    return 0;
  } else {
    debug13(printf("Going from genomicpos %u down to %u\n",start_genomicpos,end_genomicpos));
    *ambig_qends = (int *) MALLOC(nspliceends*sizeof(int));
    *ambig_probs_3 = (double *) MALLOC(nspliceends*sizeof(double));
    nspliceends = 0;

    if (sensedir == SENSE_FORWARD) {
      if (plusp) {
	*splicetype = DONOR;
	
	genomicpos = start_genomicpos;
	while (genomicpos > end_genomicpos) {
	  donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	  if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	    (*ambig_qends)[nspliceends] = (int) (genomicpos - left);
	    (*ambig_probs_3)[nspliceends] = donor_prob;
	    nspliceends += 1;
	  }
	  genomicpos--;
	}

      } else {
	*splicetype = ANTIACCEPTOR;
      
	genomicpos = start_genomicpos;
	while (genomicpos > end_genomicpos) {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	  if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	    (*ambig_qends)[nspliceends] = (int) (genomicpos - left);
	    (*ambig_probs_3)[nspliceends] = acceptor_prob;
	    nspliceends += 1;
	  }
	  genomicpos--;
	}
      }

    } else {
      if (plusp) {
	*splicetype = ANTIACCEPTOR;
      
	genomicpos = start_genomicpos;
	while (genomicpos > end_genomicpos) {
	  acceptor_prob = Maxent_hr_antiacceptor_prob(genomicpos,chroffset); /* Case 5 */
	  if (acceptor_prob >= END_SPLICESITE_PROB_MATCH) {
	    (*ambig_qends)[nspliceends] = (int) (genomicpos - left);
	    (*ambig_probs_3)[nspliceends] = acceptor_prob;
	    nspliceends += 1;
	  }
	  genomicpos--;
	}

      } else {
	*splicetype = DONOR;
	
	genomicpos = start_genomicpos;
	while (genomicpos > end_genomicpos) {
	  donor_prob = Maxent_hr_donor_prob(genomicpos,chroffset); /* Case 1 */
	  if (donor_prob >= END_SPLICESITE_PROB_MATCH) {
	    (*ambig_qends)[nspliceends] = (int) (genomicpos - left);
	    (*ambig_probs_3)[nspliceends] = donor_prob;
	    nspliceends += 1;
	  }
	  genomicpos--;
	}
      }
    }

    debug13(printf("Got %d spliceends\n",nspliceends));
    return nspliceends;
  }
}

