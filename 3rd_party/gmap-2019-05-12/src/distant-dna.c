static char rcsid[] = "$Id: distant-dna.c 218675 2019-03-16 01:25:48Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "distant-dna.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For rint */
#include "mem.h"
#include "assert.h"
#include "types.h"
#include "sense.h"

#include "genome128_hr.h"
#include "substring.h"
#include "splice.h"
#include "stage3hr.h"


/* Originally allowed only 1, to print only unique translocations.
   But need to allow enough to avoid missing some translocations. */
/* For transcript splicing, need to increase MAXCHIMERAPATHS */
/* #define MAXCHIMERAPATHS 100 */
#define MAXCHIMERAPATHS 10000


static Genome_T genomebits;
static Genome_T genomebits_alt;
static Chrpos_T shortsplicedist;


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG4E
#define debug4e(x) x
#else
#define debug4e(x)
#endif

#ifdef DEBUG4L
#define debug4l(x) x
#else
#define debug4l(x)
#endif

#ifdef DEBUG4LD
#define debug4ld(x) x
#else
#define debug4ld(x)
#endif


/* plus startfrag: Genomic_mismatches_left
   minus startfrag: Genomic_mismatches_right
   plus_endfrag: Genomic_mismatches_right
   minus_endfrag: Genomic_mismatches_left */

void
Distant_dna_solve (int *found_score_overall, int *found_score_within_trims,
		   List_T *hits_plus, List_T *hits_minus,

		   List_T startfrags_plus, List_T endfrags_plus,
		   List_T startfrags_minus, List_T endfrags_minus,

		   int *mismatch_positions_alloc,
		   Compress_T query_compress_fwd, Compress_T query_compress_rev,
		   int querylength, bool first_read_p,
		   Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  List_T localhits_plus = NULL, localhits_minus = NULL, distanthits_plus = NULL, distanthits_minus = NULL, p, q;
  double best_prob_i, best_prob_j;
  int sensedir_distant_guess;
  int ndistantsplicepairs = 0;
  Substring_T startfrag, endfrag;
  int **mismatch_positions_startfrags_plus, **mismatch_positions_endfrags_plus,
    **mismatch_positions_startfrags_minus, **mismatch_positions_endfrags_minus;
  int *nmismatches_startfrags_plus, *nmismatches_endfrags_plus,
    *nmismatches_startfrags_minus, *nmismatches_endfrags_minus;
  int best_nmismatches_i, best_nmismatches_j, total_nmismatches_i, total_nmismatches_j;
  int splice_pos;
  Chrpos_T distance;
  Univcoord_T startfrag_genomicstart, endfrag_genomicstart;
  bool shortdistancep;
  Chrnum_T chrnum;
  Stage3end_T splice;
  int n1, n2, n3, n4, n, i, j, k;

  debug(printf("Starting find_splicepairs_dna\n"));
  debug(printf("Splice ends: +startfrags/endfrags %d/%d, -startfrags/endfrags %d/%d\n",
	       List_length(startfrags_plus),List_length(endfrags_plus),
	       List_length(startfrags_minus),List_length(endfrags_minus)));
  
  n1 = List_length(startfrags_plus);
  n2 = List_length(endfrags_plus);
  n3 = List_length(startfrags_minus);
  n4 = List_length(endfrags_minus);

  debug(printf("Possibilities %f = (%d+%d) * (%d+%d)\n",((float) (n1 + n3) * (float) (n2 + n4)),n1,n3,n2,n4));

  /* Have to convert to float, because value could overflow into negative integers */
  if (((float) (n1 + n3) * (float) (n2 + n4)) > 100.0) {
    debug(printf("Too many possibilities (%f), so skipping\n",((float) (n1 + n3) * (float) (n2 + n4))));
    return;

  } else {
    if (n1 == 0) {
      mismatch_positions_startfrags_plus = (int **) NULL;
      nmismatches_startfrags_plus = (int *) NULL;
    } else {
      mismatch_positions_startfrags_plus = (int **) CALLOC(n1,sizeof(int *));
      nmismatches_startfrags_plus = (int *) MALLOC(n1*sizeof(int));
    }

    if (n2 == 0) {
      mismatch_positions_endfrags_plus = (int **) NULL;
      nmismatches_endfrags_plus = (int *) NULL;
    } else {
      mismatch_positions_endfrags_plus = (int **) CALLOC(n2,sizeof(int *));
      nmismatches_endfrags_plus = (int *) MALLOC(n2*sizeof(int));
    }
    
    if (n3 == 0) {
      mismatch_positions_startfrags_minus = (int **) NULL;
      nmismatches_startfrags_minus = (int *) NULL;
    } else {
      mismatch_positions_startfrags_minus = (int **) CALLOC(n3,sizeof(int *));
      nmismatches_startfrags_minus = (int *) MALLOC(n3*sizeof(int));
    }
    
    if (n4 == 0) {
      mismatch_positions_endfrags_minus = (int **) NULL;
      nmismatches_endfrags_minus = (int *) NULL;
    } else {
      mismatch_positions_endfrags_minus = (int **) CALLOC(n4,sizeof(int *));
      nmismatches_endfrags_minus = (int *) MALLOC(n4*sizeof(int));
    }
  }


  /************************************************************************
   *   Same strands
   ************************************************************************/

  /* 1.  End 1 to End 2.  Same strands. */
  debug4l(printf("find_splicepairs_distant_dna: startfrags+ (%d) to endfrags+ (%d)\n",
		 List_length(startfrags_plus),List_length(endfrags_plus)));
  for (p = startfrags_plus, i = 0; p != NULL; p = p->rest, i++) {
    startfrag = (Substring_T) p->first;
    for (q = endfrags_plus, j = 0;
	 q != NULL && Substring_siteN_pos((Substring_T) q->first) <= Substring_siteN_pos(startfrag);
	 q = q->rest, j++) {
      endfrag = (Substring_T) q->first;

      debug4ld(printf("end1-end2: startfrag at %llu %d..%d and endfrag at %llu %d..%d\n",
		      (unsigned long long) Substring_left(startfrag),
		      Substring_querystart(startfrag),Substring_queryend(startfrag),
		      (unsigned long long) Substring_left(endfrag),
		      Substring_querystart(endfrag),Substring_queryend(endfrag)));

      if (mismatch_positions_startfrags_plus[i] == NULL) {
	total_nmismatches_i = Substring_nmismatches_bothdiff(startfrag);
	mismatch_positions_startfrags_plus[i] = MALLOC((total_nmismatches_i+2)*sizeof(int));
	/* Values are ascending */
	n = Genome_mismatches_left(mismatch_positions_alloc,/*max_mismatches_allowed*/total_nmismatches_i,
				   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_fwd,
				   /*left*/Substring_left(startfrag),
				   /*pos5*/Substring_querystart(startfrag),/*pos3*/Substring_queryend(startfrag),
				   /*plusp*/true,Substring_genestrand(startfrag));
#ifdef DEBUG4E
	printf("(1a) For qpos %d..%d:",/*pos5*/Substring_querystart(startfrag),/*pos3*/Substring_queryend(startfrag));
	for (k = 0; k <= n; k++) {
	  printf(" %d",mismatch_positions_alloc[k]);
	}
	printf(" => ");
#endif

	nmismatches_startfrags_plus[i] = n;

	for (k = 0; k <= n; k++) {
	  mismatch_positions_startfrags_plus[i][k] = mismatch_positions_alloc[k];
	}

#ifdef DEBUG4E
	for (k = 0; k <= n; k++) {
	  printf(" %d",mismatch_positions_startfrags_plus[i][k]);
	}
	printf("\n");
#endif
      }

      if (mismatch_positions_endfrags_plus[j] == NULL) {
	total_nmismatches_j = Substring_nmismatches_bothdiff(endfrag);
	mismatch_positions_endfrags_plus[j] = MALLOC((total_nmismatches_j+2)*sizeof(int));
	/* Values are descending */
	n = Genome_mismatches_right(mismatch_positions_alloc,/*max_mismatches_allowed*/total_nmismatches_j,
				    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_fwd,
				    /*left*/Substring_left(endfrag),
				    /*pos5*/Substring_querystart(endfrag),/*pos3*/Substring_queryend(endfrag),
				    /*plusp*/true,Substring_genestrand(endfrag));
#ifdef DEBUG4E
	printf("(2a) For qpos %d..%d:",/*pos5*/Substring_querystart(endfrag),/*pos3*/Substring_queryend(endfrag));
	for (k = 0; k <= n; k++) {
	  printf(" %d",mismatch_positions_alloc[k]);
	}
	printf(" => ");
#endif

	nmismatches_endfrags_plus[j] = n;

	for (k = 0; k <= n; k++) {
	  mismatch_positions_endfrags_plus[j][k] = mismatch_positions_alloc[n - k];
	}

#ifdef DEBUG4E
	for (k = 0; k <= n; k++) {
	  printf(" %d",mismatch_positions_endfrags_plus[j][k]);
	}
	printf("\n");
#endif
      }

      if ((splice_pos = Splice_resolve_distant(&best_nmismatches_i,&best_nmismatches_j,
					       &best_prob_i,&best_prob_j,&sensedir_distant_guess,
					       mismatch_positions_startfrags_plus[i],nmismatches_startfrags_plus[i],
					       mismatch_positions_endfrags_plus[j],nmismatches_endfrags_plus[j],
					       Substring_left(startfrag),Substring_left(endfrag),
					       Substring_chroffset(startfrag),Substring_chroffset(endfrag),
					       /*querystart*/Substring_querystart(startfrag),/*queryend*/Substring_queryend(endfrag),
					       /*splice_pos_start*/Substring_querystart(endfrag),
					       /*splice_pos_end*/Substring_queryend(startfrag),querylength,
					       /*plusp_i*/true,/*plusp_j*/true)) > 0) {
	if ((chrnum = Substring_chrnum(startfrag)) != Substring_chrnum(endfrag)) {
	  distance = 0U;
	  shortdistancep = false;
	} else if ((endfrag_genomicstart = Substring_genomicstart(endfrag)) > (startfrag_genomicstart = Substring_genomicstart(startfrag))) {
	  distance = endfrag_genomicstart - startfrag_genomicstart;
	  if (distance <= shortsplicedist) {
	    shortdistancep = true;
#if 0
	  } else if (distances_observed_p == true && intragenic_splice_p(distance,startfrag,endfrag) == true) {
	    shortdistancep = true;
#endif
	  } else {
	    shortdistancep = false;
	  }
	} else {
	  distance = startfrag_genomicstart - endfrag_genomicstart;
	  shortdistancep = false; /* scramble */
	}
	debug4ld(printf("1-2. Pushing a candidate at break_pos %d (%d..%d), startfrag %llu to endfrag %llu.  shortdistancep = %d, splicep %d and %d\n",
			splice_pos,Substring_queryend(startfrag),Substring_querystart(endfrag),
			(unsigned long long) Substring_left(startfrag),
			(unsigned long long) Substring_left(endfrag),shortdistancep,
			Substring_trim_querystart_splicep(startfrag),Substring_trim_queryend_splicep(endfrag)));
	
	if (shortdistancep) {
	  if ((splice = Stage3end_new_distant(&(*found_score_overall),&(*found_score_within_trims),
					      startfrag,endfrag,splice_pos,best_nmismatches_i,best_nmismatches_j,
					      best_prob_i,best_prob_j,sensedir_distant_guess,
					      distance,/*shortdistancep*/true,querylength,
					      first_read_p,listpool,level)) != NULL) {
	    localhits_plus = Hitlist_push(localhits_plus,hitlistpool,(void *) splice);
	  }
	} else {
	  if ((splice = Stage3end_new_distant(&(*found_score_overall),&(*found_score_within_trims),
					      startfrag,endfrag,splice_pos,best_nmismatches_i,best_nmismatches_j,
					      best_prob_i,best_prob_j,sensedir_distant_guess,
					      distance,/*shortdistancep*/false,querylength,
					      first_read_p,listpool,level)) != NULL) {
	    distanthits_plus = Hitlist_push(distanthits_plus,hitlistpool,(void *) splice);
	    ndistantsplicepairs++;
	  }
	}
      }
    }
  }

  /* 4. End 3 to End 4.  Same strands. */
  debug4l(printf("find_splicepairs_distant_dna: startfrags- (%d) to endfrags- (%d)\n",
		 List_length(startfrags_minus),List_length(endfrags_minus)));
  for (p = startfrags_minus, i = 0; p != NULL; p = p->rest, i++) {
    startfrag = (Substring_T) p->first;
    for (q = endfrags_minus, j = 0;
	 q != NULL && Substring_siteN_pos((Substring_T) q->first) <= Substring_siteN_pos(startfrag);
	 q = q->rest, j++) {
      endfrag = (Substring_T) q->first;

      debug4ld(printf("end3-end4: startfrag at %llu %d..%d and endfrag at %llu %d..%d\n",
		      (unsigned long long) Substring_left(startfrag),
		      Substring_querystart(startfrag),Substring_queryend(startfrag),
		      (unsigned long long) Substring_left(endfrag),
		      Substring_querystart(endfrag),Substring_queryend(endfrag)));
    
      if (mismatch_positions_startfrags_minus[i] == NULL) {
	total_nmismatches_i = Substring_nmismatches_bothdiff(startfrag);
	mismatch_positions_startfrags_minus[i] = MALLOC((total_nmismatches_i+2)*sizeof(int));
	/* Values are descending */
	n = Genome_mismatches_right(mismatch_positions_alloc,/*max_mismatches_allowed*/total_nmismatches_i,
				    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_rev,
				    /*left*/Substring_left(startfrag),
				    /*pos5*/querylength - Substring_queryend(startfrag),
				    /*pos3*/querylength - Substring_querystart(startfrag),
				    /*plusp*/false,Substring_genestrand(startfrag));
#ifdef DEBUG4E
	printf("(3a) For qpos %d..%d:",
	       /*pos5*/querylength - Substring_queryend(startfrag),
	       /*pos3*/querylength - Substring_querystart(startfrag));
	for (k = 0; k <= n; k++) {
	  printf(" %d",mismatch_positions_alloc[k]);
	}
	printf(" => ");
#endif

	nmismatches_startfrags_minus[i] = n;

	for (k = 0; k <= n; k++) {
	  mismatch_positions_startfrags_minus[i][k] = (querylength - 1) - mismatch_positions_alloc[k];
	}

#ifdef DEBUG4E
	for (k = 0; k <= n; k++) {
	  printf(" %d",mismatch_positions_startfrags_minus[i][k]);
	}
	printf("\n");
#endif
      }

      if (mismatch_positions_endfrags_minus[j] == NULL) {
	total_nmismatches_j = Substring_nmismatches_bothdiff(endfrag);
	mismatch_positions_endfrags_minus[j] = MALLOC((total_nmismatches_j+2)*sizeof(int));
	/* Values are ascending */
	n = Genome_mismatches_left(mismatch_positions_alloc,/*max_mismatches_allowed*/total_nmismatches_j,
				   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_rev,
				   /*left*/Substring_left(endfrag),
				   /*pos5*/querylength - Substring_queryend(endfrag),
				   /*pos3*/querylength - Substring_querystart(endfrag),
				   /*plusp*/false,Substring_genestrand(endfrag));
#ifdef DEBUG4E
	printf("(4a) For qpos %d..%d:",
	       /*pos5*/querylength - Substring_queryend(endfrag),
	       /*pos3*/querylength - Substring_querystart(endfrag));
	for (k = 0; k <= n; k++) {
	  printf(" %d",mismatch_positions_alloc[k]);
	}
	printf(" => ");
#endif

	nmismatches_endfrags_minus[j] = n;

	for (k = 0; k <= n; k++) {
	  mismatch_positions_endfrags_minus[j][k] = (querylength - 1) - mismatch_positions_alloc[n - k];
	}

#ifdef DEBUG4E
	for (k = 0; k <= n; k++) {
	  printf(" %d",mismatch_positions_endfrags_minus[j][k]);
	}
	printf("\n");
#endif
      }

      if ((splice_pos = Splice_resolve_distant(&best_nmismatches_i,&best_nmismatches_j,
					       &best_prob_i,&best_prob_j,&sensedir_distant_guess,
					       mismatch_positions_startfrags_minus[i],nmismatches_startfrags_minus[i],
					       mismatch_positions_endfrags_minus[j],nmismatches_endfrags_minus[j],
					       Substring_left(startfrag),Substring_left(endfrag),
					       Substring_chroffset(startfrag),Substring_chroffset(endfrag),
					       /*querystart*/Substring_querystart(startfrag),/*queryend*/Substring_queryend(endfrag),
					       /*splice_pos_start*/Substring_querystart(endfrag),
					       /*splice_pos_end*/Substring_queryend(startfrag),querylength,
					       /*plusp_i*/false,/*plusp_j*/false)) > 0) {
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
#if 0
	  } else if (distances_observed_p == true && intragenic_splice_p(distance,startfrag,endfrag) == true) {
	    shortdistancep = true;
#endif
	  } else {
	    shortdistancep = false;
	  }
	}
	debug4ld(printf("3-4. Pushing a candidate at break_pos %d (%d..%d), startfrag %llu to endfrag %llu.  shortdistancep = %d, splicep %d and %d\n",
			splice_pos,Substring_queryend(startfrag),Substring_querystart(endfrag),
			(unsigned long long) Substring_left(startfrag),
			(unsigned long long) Substring_left(endfrag),shortdistancep,
			Substring_trim_querystart_splicep(startfrag),Substring_trim_queryend_splicep(endfrag)));

	if (shortdistancep) {
	  if ((splice = Stage3end_new_distant(&(*found_score_overall),&(*found_score_within_trims),
					      startfrag,endfrag,splice_pos,best_nmismatches_i,best_nmismatches_j,
					      best_prob_i,best_prob_j,sensedir_distant_guess,
					      distance,/*shortdistancep*/true,querylength,
					      first_read_p,listpool,level)) != NULL) {
	    localhits_minus = Hitlist_push(localhits_minus,hitlistpool,(void *) splice);
	  }
	} else {
	  if ((splice = Stage3end_new_distant(&(*found_score_overall),&(*found_score_within_trims),
					      startfrag,endfrag,splice_pos,best_nmismatches_i,best_nmismatches_j,
					      best_prob_i,best_prob_j,sensedir_distant_guess,
					      distance,/*shortdistancep*/false,querylength,
					      first_read_p,listpool,level)) != NULL) {
	    distanthits_minus = Hitlist_push(distanthits_minus,hitlistpool,(void *) splice);
	    ndistantsplicepairs++;
	  }
	}
      }
    }
  }
  
  /* 5. End 5 to End 6.  Same strands. */
  /* 8. End 7 to End 8.  Same strands. */


  if (localhits_plus != NULL || localhits_minus != NULL) {
    /* A local splice takes precedence over distant splices */
    Stage3end_gc(distanthits_plus);
    Stage3end_gc(distanthits_minus);

    *hits_plus = List_append(*hits_plus,localhits_plus);
    *hits_minus = List_append(*hits_minus,localhits_minus);
    /* Return at bottom after freeing memory */

  } else {
    *hits_plus = List_append(*hits_plus,distanthits_plus);
    *hits_minus = List_append(*hits_minus,distanthits_minus);

    /************************************************************************
     *   Different strands
     ************************************************************************/
  
    /* 2. End 1 to End 4.  Different strands. */
    debug4l(printf("find_splicepairs_distant_dna: startfrags+ (%d) to endfrags- (%d)\n",
		   List_length(startfrags_plus),List_length(endfrags_minus)));
    for (p = startfrags_plus, i = 0; p != NULL; p = p->rest, i++) {
      startfrag = (Substring_T) p->first;
      for (q = endfrags_minus, j = 0;
	   q != NULL && Substring_siteN_pos((Substring_T) q->first) <= Substring_siteN_pos(startfrag);
	   q = q->rest, j++) {
	endfrag = (Substring_T) q->first;

	debug4ld(printf("end1-end4: startfrag at %llu %d..%d and endfrag at %llu %d..%d\n",
			(unsigned long long) Substring_left(startfrag),
			Substring_querystart(startfrag),Substring_queryend(startfrag),
			(unsigned long long) Substring_left(endfrag),
			Substring_querystart(endfrag),Substring_queryend(endfrag)));
	
	if (mismatch_positions_startfrags_plus[i] == NULL) {
	  total_nmismatches_i = Substring_nmismatches_bothdiff(startfrag);
	  mismatch_positions_startfrags_plus[i] = MALLOC((total_nmismatches_i+2)*sizeof(int));
	  /* Values are ascending */
	  n = Genome_mismatches_left(mismatch_positions_alloc,/*max_mismatches_allowed*/total_nmismatches_i,
				     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_fwd,
				     /*left*/Substring_left(startfrag),
				     /*pos5*/Substring_querystart(startfrag),/*pos3*/Substring_queryend(startfrag),
				     /*plusp*/true,Substring_genestrand(startfrag));
#ifdef DEBUG4E
	  printf("(1b) For qpos %d..%d:",/*pos5*/Substring_querystart(startfrag),/*pos3*/Substring_queryend(startfrag));
	  for (k = 0; k <= n; k++) {
	    printf(" %d",mismatch_positions_alloc[k]);
	  }
	  printf(" => ");
#endif
	  
	  nmismatches_startfrags_plus[i] = n;
	  
	  for (k = 0; k <= n; k++) {
	    mismatch_positions_startfrags_plus[i][k] = mismatch_positions_alloc[k];
	  }
	  
#ifdef DEBUG4E
	  for (k = 0; k <= n; k++) {
	    printf(" %d",mismatch_positions_startfrags_plus[i][k]);
	  }
	  printf("\n");
#endif
	}

	if (mismatch_positions_endfrags_minus[j] == NULL) {
	  total_nmismatches_j = Substring_nmismatches_bothdiff(endfrag);
	  mismatch_positions_endfrags_minus[j] = MALLOC((total_nmismatches_j+2)*sizeof(int));
	  /* Values are ascending */
	  n = Genome_mismatches_left(mismatch_positions_alloc,/*max_mismatches_allowed*/total_nmismatches_j,
				     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_rev,
				     /*left*/Substring_left(endfrag),
				     /*pos5*/querylength - Substring_queryend(endfrag),
				     /*pos3*/querylength - Substring_querystart(endfrag),
				     /*plusp*/false,Substring_genestrand(endfrag));
#ifdef DEBUG4E
	  printf("(4b) For qpos %d..%d:",
		 /*pos5*/querylength - Substring_queryend(endfrag),
		 /*pos3*/querylength - Substring_querystart(endfrag));
	  for (k = 0; k <= n; k++) {
	    printf(" %d",mismatch_positions_alloc[k]);
	  }
	  printf(" => ");
#endif

	  nmismatches_endfrags_minus[j] = n;

	  for (k = 0; k <= n; k++) {
	    mismatch_positions_endfrags_minus[j][k] = (querylength - 1) - mismatch_positions_alloc[n - k];
	  }

#ifdef DEBUG4E
	  for (k = 0; k <= n; k++) {
	    printf(" %d",mismatch_positions_endfrags_minus[j][k]);
	  }
	  printf("\n");
#endif
	}

	if ((splice_pos = Splice_resolve_distant(&best_nmismatches_i,&best_nmismatches_j,
						 &best_prob_i,&best_prob_j,&sensedir_distant_guess,
						 mismatch_positions_startfrags_plus[i],nmismatches_startfrags_plus[i],
						 mismatch_positions_endfrags_minus[j],nmismatches_endfrags_minus[j],
						 Substring_left(startfrag),Substring_left(endfrag),
						 Substring_chroffset(startfrag),Substring_chroffset(endfrag),
						 /*querystart*/Substring_querystart(startfrag),/*queryend*/Substring_queryend(endfrag),
						 /*splice_pos_start*/Substring_querystart(endfrag),
						 /*splice_pos_end*/Substring_queryend(startfrag),querylength,
						 /*plusp_i*/true,/*plusp_j*/false)) > 0) {
	  if (Substring_chrnum(startfrag) != Substring_chrnum(endfrag)) {
	    distance = 0U;
	  } else if ((Substring_genomicstart(endfrag) - splice_pos) > (Substring_genomicstart(startfrag) + splice_pos)) {
	    distance = (Substring_genomicstart(endfrag) - splice_pos) - (Substring_genomicstart(startfrag) + splice_pos);
	  } else {
	    distance = (Substring_genomicstart(startfrag) + splice_pos) - (Substring_genomicstart(endfrag) - splice_pos);
	  }
	  debug4ld(printf("1-4. Pushing a candidate at break_pos %d (%d..%d), startfrag %llu to endfrag %llu.  Different strands, so not shortdistance\n",
			  splice_pos,Substring_queryend(startfrag),Substring_querystart(endfrag),
			  (unsigned long long) Substring_left(startfrag),
			  (unsigned long long) Substring_left(endfrag)));
	  if ((splice = Stage3end_new_distant(&(*found_score_overall),&(*found_score_within_trims),
					      startfrag,endfrag,splice_pos,best_nmismatches_i,best_nmismatches_j,
					      best_prob_i,best_prob_j,sensedir_distant_guess,
					      distance,/*shortdistancep*/false,querylength,
					      first_read_p,listpool,level)) != NULL) {
	    if (Stage3end_plusp(splice) == true) {
	      /* Determined by substring_for_concordance */
	      *hits_plus = Hitlist_push(*hits_plus,hitlistpool,(void *) splice);
	    } else {
	      *hits_minus = Hitlist_push(*hits_minus,hitlistpool,(void *) splice);
	    }
	    ndistantsplicepairs++;
	  }
	}
      }
    }
  
    /* 3. End 3 to End 2.  Different strands. */
    debug4l(printf("find_splicepairs_distant_dna: startfrags- (%d) to endfrags+ (%d)\n",
		   List_length(startfrags_minus),List_length(endfrags_plus)));
    for (p = startfrags_minus, i = 0; p != NULL; p = p->rest, i++) {
      startfrag = (Substring_T) p->first;
      for (q = endfrags_plus, j = 0;
	   q != NULL && Substring_siteN_pos((Substring_T) q->first) <= Substring_siteN_pos(startfrag);
	   q = q->rest, j++) {
	endfrag = (Substring_T) q->first;

	debug4ld(printf("end3-end2: startfrag at %llu %d..%d and endfrag at %llu %d..%d\n",
			(unsigned long long) Substring_left(startfrag),
			Substring_querystart(startfrag),Substring_queryend(startfrag),
			(unsigned long long) Substring_left(endfrag),
			Substring_querystart(endfrag),Substring_queryend(endfrag)));
    
	if (mismatch_positions_startfrags_minus[i] == NULL) {
	  total_nmismatches_i = Substring_nmismatches_bothdiff(startfrag);
	  mismatch_positions_startfrags_minus[i] = MALLOC((total_nmismatches_i+2)*sizeof(int));
	  /* Values are descending */
	  n = Genome_mismatches_right(mismatch_positions_alloc,/*max_mismatches_allowed*/total_nmismatches_i,
				      /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_rev,
				      /*left*/Substring_left(startfrag),
				      /*pos5*/querylength - Substring_queryend(startfrag),
				      /*pos3*/querylength - Substring_querystart(startfrag),
				      /*plusp*/false,Substring_genestrand(startfrag));
#ifdef DEBUG4E
	  printf("(3b) For qpos %d..%d:",
		 /*pos5*/querylength - Substring_queryend(startfrag),
		 /*pos3*/querylength - Substring_querystart(startfrag));
	  for (k = 0; k <= n; k++) {
	    printf(" %d",mismatch_positions_alloc[k]);
	  }
	  printf(" => ");
#endif

	  nmismatches_startfrags_minus[i] = n;

	  for (k = 0; k <= n; k++) {
	    mismatch_positions_startfrags_minus[i][k] = (querylength - 1) - mismatch_positions_alloc[k];
	  }

#ifdef DEBUG4E
	  for (k = 0; k <= n; k++) {
	    printf(" %d",mismatch_positions_startfrags_minus[i][k]);
	  }
	  printf("\n");
#endif
	}

	if (mismatch_positions_endfrags_plus[j] == NULL) {
	  total_nmismatches_j = Substring_nmismatches_bothdiff(endfrag);
	  mismatch_positions_endfrags_plus[j] = MALLOC((total_nmismatches_j+2)*sizeof(int));
	  /* Values are descending */
	  n = Genome_mismatches_right(mismatch_positions_alloc,/*max_mismatches_allowed*/total_nmismatches_j,
				      /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_fwd,
				      /*left*/Substring_left(endfrag),
				      /*pos5*/Substring_querystart(endfrag),/*pos3*/Substring_queryend(endfrag),
				      /*plusp*/true,Substring_genestrand(endfrag));
#ifdef DEBUG4E
	  printf("(2b) For qpos %d..%d:",/*pos5*/Substring_querystart(endfrag),/*pos3*/Substring_queryend(endfrag));
	  for (k = 0; k <= n; k++) {
	    printf(" %d",mismatch_positions_alloc[k]);
	  }
	  printf(" => ");
#endif

	  nmismatches_endfrags_plus[j] = n;

	  for (k = 0; k <= n; k++) {
	    mismatch_positions_endfrags_plus[j][k] = mismatch_positions_alloc[n - k];
	  }

#ifdef DEBUG4E
	  for (k = 0; k <= n; k++) {
	    printf(" %d",mismatch_positions_endfrags_plus[j][k]);
	  }
	  printf("\n");
#endif
	}

	if ((splice_pos = Splice_resolve_distant(&best_nmismatches_i,&best_nmismatches_j,
						 &best_prob_i,&best_prob_j,&sensedir_distant_guess,
						 mismatch_positions_startfrags_minus[i],nmismatches_startfrags_minus[i],
						 mismatch_positions_endfrags_plus[j],nmismatches_endfrags_plus[j],
						 Substring_left(startfrag),Substring_left(endfrag),
						 Substring_chroffset(startfrag),Substring_chroffset(endfrag),
						 /*querystart*/Substring_querystart(startfrag),/*queryend*/Substring_queryend(endfrag),
						 /*splice_pos_start*/Substring_querystart(endfrag),
						 /*splice_pos_end*/Substring_queryend(startfrag),querylength,
						 /*plusp_i*/false,/*plusp_j*/true)) > 0) {
	  if (Substring_chrnum(startfrag) != Substring_chrnum(endfrag)) {
	    distance = 0U;
	  } else if (Substring_genomicstart(endfrag) > Substring_genomicstart(startfrag)) {
	    distance = (Substring_genomicstart(endfrag) + splice_pos) - (Substring_genomicstart(startfrag) - splice_pos);
	  } else {
	    distance = (Substring_genomicstart(startfrag) - splice_pos) - (Substring_genomicstart(endfrag) + splice_pos);
	  }
	  debug4ld(printf("3-2. Pushing a candidate at break_pos %d (%d..%d), startfrag %llu to endfrag %llu.  Different strands so not shortdistance.\n",
			  splice_pos,Substring_queryend(startfrag),Substring_querystart(endfrag),
			  (unsigned long long) Substring_left(startfrag),
			  (unsigned long long) Substring_left(endfrag)));
	  if ((splice = Stage3end_new_distant(&(*found_score_overall),&(*found_score_within_trims),
					      startfrag,endfrag,splice_pos,best_nmismatches_i,best_nmismatches_j,
					      best_prob_i,best_prob_j,sensedir_distant_guess,
					      distance,/*shortdistancep*/false,querylength,
					      first_read_p,listpool,level)) != NULL) {
	    if (Stage3end_plusp(splice) == true) {
	      /* Determined by substring_for_concordance */
	      *hits_plus = Hitlist_push(*hits_plus,hitlistpool,(void *) splice);
	    } else {
	      *hits_minus = Hitlist_push(*hits_minus,hitlistpool,(void *) splice);
	    }
	    ndistantsplicepairs++;
	  }
	}
      }
    }
  
    /* 6. End 5 to End 8.  Different strands. */
    /* 7. End 7 to End 6.  Different strands. */
  }


  for (i = 0; i < List_length(startfrags_plus); i++) {
    FREE(mismatch_positions_startfrags_plus[i]);
  }
  FREE(mismatch_positions_startfrags_plus);
  FREE(nmismatches_startfrags_plus);

  for (j = 0; j < List_length(endfrags_plus); j++) {
    FREE(mismatch_positions_endfrags_plus[j]);
  }
  FREE(mismatch_positions_endfrags_plus);
  FREE(nmismatches_endfrags_plus);


  for (i = 0; i < List_length(startfrags_minus); i++) {
    FREE(mismatch_positions_startfrags_minus[i]);
  }
  FREE(mismatch_positions_startfrags_minus);
  FREE(nmismatches_startfrags_minus);

  for (j = 0; j < List_length(endfrags_minus); j++) {
    FREE(mismatch_positions_endfrags_minus[j]);
  }
  FREE(mismatch_positions_endfrags_minus);
  FREE(nmismatches_endfrags_minus);


#if 0
  debug4l(printf("ndistantsplicepairs %d, maxchimerapaths %d\n",ndistantsplicepairs,MAXCHIMERAPATHS));
  debug4ld(printf("ndistantsplicepairs %d, maxchimerapaths %d\n",ndistantsplicepairs,MAXCHIMERAPATHS));
  if (0 && ndistantsplicepairs > MAXCHIMERAPATHS) {
    /* Can afford to ignore these if MAXCHIMERAPATHS is set high enough */
    stage3list_gc(&distantsplicing);
    return (List_T) NULL;
  } else {
    return distantsplicing;
  }
#endif

  return;
}


#if 0
distantsplicing = find_splicepairs_dna(&found_score,&nsplicepairs,&longsinglesplicing,distantsplicing,
				       startfrags_plus,endfrags_plus,startfrags_minus,endfrags_minus,
				       localsplicing_penalty,distantsplicing_penalty,
				       querylength,nmismatches,first_read_p);

if (longsinglesplicing != NULL) {
  debug(printf("Entering Stage3end_optimal_score with %d longsinglesplicing hits\n",List_length(longsinglesplicing)));
  longsinglesplicing = Stage3end_optimal_score(longsinglesplicing,query_compress_fwd,query_compress_rev,querylength,
					       /*keep_gmap_p*/true,/*finalp*/false);
  debug(printf("Exiting Stage3end_optimal_score with %d hits\n",List_length(longsinglesplicing)));
  hits = List_append(hits,longsinglesplicing);
  
  opt_level = (found_score < opt_level) ? found_score : opt_level;
  if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
    done_level = user_maxlevel;
  }
 }

if (distantsplicing != NULL) {
  /* Excess distant splicing should be freed already in find_splicepairs_distant_rna */
  debug(printf("Entering Stage3end_optimal_score with %d hits\n",List_length(distantsplicing)));
  distantsplicing = Stage3end_optimal_score(distantsplicing,query_compress_fwd,query_compress_rev,querylength,
					    /*keep_gmap_p*/true,/*finalp*/false);
  debug(printf("Exiting Stage3end_optimal_score with %d hits\n",List_length(distantsplicing)));
  
  hits = List_append(hits,distantsplicing);

  opt_level = (found_score < opt_level) ? found_score : opt_level;
  if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
    done_level = user_maxlevel;
  }
 }
debug(printf("9 (DNA)> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
#endif



void
Distant_dna_setup (Genome_T genomebits_in, Genome_T genomebits_alt_in,
		   Chrpos_T shortsplicedist_in) {

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  shortsplicedist = shortsplicedist_in;

  return;
}
