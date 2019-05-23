static char rcsid[] = "$Id: indel.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "indel.h"

#include "assert.h"
#include "mem.h"
#include "genome128_hr.h"
#include "stage3hr.h"
#include "intron.h"


/* Resolve indels */ 
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Solve end indels */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


static int min_indel_end_matches;
static int indel_penalty_middle;


void
Indel_setup (int min_indel_end_matches_in, int indel_penalty_middle_in) {
  min_indel_end_matches = min_indel_end_matches_in;
  indel_penalty_middle = indel_penalty_middle_in;
  return;
}


/* For transcriptome alignments, plusp may be true or false */
/* For alignments via middle_path with ascending univdiags, plusp should be true */
/* indels is positive here */
int
Indel_resolve_middle_insertion (int *best_nmismatches_i, int *best_nmismatches_j,
				Univcoord_T left, int indels,
				int *mismatch_positions_left, int nmismatches_left,
				int *mismatch_positions_right, int nmismatches_right,
				Genome_T ome, Genome_T ome_alt, Compress_T query_compress,
				int querystart, int queryend, int querylength,
				int nmismatches_allowed, bool plusp, int genestrand,
				bool want_lowest_coordinate_p) {
  int best_indel_pos = -1, indel_pos;
#ifdef DEBUG2
  int i;
  char *gbuffer;
#endif
  int best_sum, sum, nmismatches_lefti, nmismatches_righti, lefti, righti;
  bool free_left_p = false, free_right_p = false;

  assert(indels > 0);
  if (nmismatches_allowed > querylength) {
    nmismatches_allowed = querylength;
  }

  if (mismatch_positions_left == NULL) {
    free_left_p = true;
#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      mismatch_positions_left = (int *) ALLOCA(querylength * sizeof(int));
    } else {
      mismatch_positions_left = (int *) MALLOC(querylength * sizeof(int));
    }
#else
    mismatch_positions_left = (int *) MALLOC(querylength * sizeof(int));
#endif
  debug2(printf("nmismatches_allowed is %d.  Calling Genome_mismatches_left over %d..%d\n",
		nmismatches_allowed,querystart,queryend));
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,nmismatches_allowed,
					    ome,ome_alt,query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
					    plusp,genestrand);
  debug2(
	 printf("%d mismatches on left at:",nmismatches_left);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  }

  if (mismatch_positions_right == NULL) {
    free_right_p = true;
#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      mismatch_positions_right = (int *) ALLOCA(querylength * sizeof(int));
    } else {
      mismatch_positions_right = (int *) MALLOC(querylength * sizeof(int));
    }
#else
    mismatch_positions_right = (int *) MALLOC(querylength * sizeof(int));
#endif
    debug2(printf("nmismatches_allowed is %d.  Calling Genome_mismatches_right over %d..%d\n",
		  nmismatches_allowed,querystart,queryend));
    nmismatches_right = Genome_mismatches_right(mismatch_positions_right,nmismatches_allowed,
						ome,ome_alt,query_compress,left-indels,/*pos5*/querystart,/*pos3*/queryend,
						plusp,genestrand);
    debug2(
	   printf("%d mismatches on right at:",nmismatches_right);
	   for (i = 0; i <= nmismatches_right; i++) {
	     printf(" %d",mismatch_positions_right[i]);
	   }
	   printf("\n");
	   );
  }

  /* query has insertion.  Get |indels| less from genome; trim from left. */
  /* left = ptr->diagonal - querylength; */

#ifdef DEBUG2
  gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char));
  Genome_fill_buffer_blocks(left+indels,querylength-indels,gbuffer);
  printf("solve_middle_indel, plusp %d, want low %d, insertion: Getting genome at left %llu + indels %d = %llu\n",
	 plusp,want_lowest_coordinate_p,(unsigned long long) left,indels,(unsigned long long) left+indels);
  printf("g1: %s\n",gbuffer);
  printf("g2: %s\n",&(gbuffer[indels]));
  FREE(gbuffer);
#endif

  best_sum = querylength + querylength;

  if (want_lowest_coordinate_p == true) {
    /* Modeled after end D to get lowest possible coordinate */
    righti = 0;
    lefti = nmismatches_left - 1;
    nmismatches_righti = /*righti*/ 0;
    nmismatches_lefti = /*lefti+1*/ nmismatches_left;

    while (righti < nmismatches_right) {
      while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti] - indels) {
	lefti--;
      }
      sum = righti + lefti + 1;
      debug2(printf("  (Case D) sum %d=%d+%d at indel_pos %d.",
		    sum,righti,lefti+1,mismatch_positions_right[righti]-indels+1));
      if (sum <= best_sum) {
	indel_pos = mismatch_positions_right[righti] - indels + 1;
	if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches &&
	    /* Require separation from endpoints */
	    indel_pos > querystart && indel_pos + indels < queryend) {
	  best_indel_pos = indel_pos;
	  nmismatches_righti = righti;
	  nmismatches_lefti = lefti + 1;
	  debug2(printf("**"));
	  best_sum = sum;
	}
      }
      righti++;
    }
    debug2(printf("\n"));


    /* Try from other side to see if we missed anything */
    lefti = 0;
    righti = nmismatches_right - 1;

    while (lefti < nmismatches_left) {
      while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti] + indels) {
	righti--;
      }
      sum = lefti + righti + 1;
      debug2(printf("  (Case D2) sum %d=%d+%d at indel_pos %d.",
		    sum,lefti,righti+1,mismatch_positions_left[lefti]));
      if (sum < best_sum) {
	indel_pos = mismatch_positions_left[lefti];
	if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches &&
	    /* Require separation from endpoints */
	    indel_pos > querystart && indel_pos + indels < queryend) {
	  best_indel_pos = indel_pos;
	  nmismatches_righti = righti + 1;
	  nmismatches_lefti = lefti;
	  debug2(printf("**"));
	  best_sum = sum;
	}
      } else if (sum == best_sum) {
	indel_pos = mismatch_positions_left[lefti];
	if (indel_pos < best_indel_pos) {
	  if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos + indels < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_righti = righti + 1;
	    nmismatches_lefti = lefti;
	    debug2(printf("**"));
	    /* best_sum = sum; */
	  }
	}
      }
      lefti++;
    }
    debug2(printf("\n"));

  } else {
    /* Want highest possible coordinate.  Modified from code above, but not sure if this is correct */
    lefti = 0;
    righti = nmismatches_right - 1;
    nmismatches_lefti = /*lefti*/ 0;
    nmismatches_righti = /*righti+1*/ nmismatches_right;

    while (lefti < nmismatches_left) {
      while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti] + indels) {
	righti--;
      }
      sum = lefti + righti + 1;
      debug2(printf("  (Case X2) sum %d=%d+%d at indel_pos %d.",
		    sum,lefti,righti+1,mismatch_positions_left[lefti]));
      if (sum <= best_sum) {
	indel_pos = mismatch_positions_left[lefti];
	if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches &&
	    /* Require separation from endpoints */
	    indel_pos > querystart && indel_pos + indels < queryend) {
	  best_indel_pos = indel_pos;
	  nmismatches_righti = righti + 1;
	  nmismatches_lefti = lefti;
	  debug2(printf("**"));
	  best_sum = sum;
	}
      }
      lefti++;
    }
    debug2(printf("\n"));


    /* Try from other side to see if we missed anything */
    righti = 0;
    lefti = nmismatches_left - 1;

    while (righti < nmismatches_right) {
      while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti] - indels) {
	lefti--;
      }
      sum = righti + lefti + 1;
      debug2(printf("  (Case X) sum %d=%d+%d at indel_pos %d.",
		    sum,righti,lefti+1,mismatch_positions_right[righti]-indels+1));
      if (sum < best_sum) {
	indel_pos = mismatch_positions_right[righti] - indels + 1;
	if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches &&
	    /* Require separation from endpoints */
	    indel_pos > querystart && indel_pos + indels < queryend) {
	  best_indel_pos = indel_pos;
	  nmismatches_righti = righti;
	  nmismatches_lefti = lefti + 1;
	  debug2(printf("**"));
	  best_sum = sum;
	}
      } else if (sum == best_sum) {
	indel_pos = mismatch_positions_right[righti] - indels + 1;
	if (indel_pos > best_indel_pos) {
	  if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos + indels < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_righti = righti;
	    nmismatches_lefti = lefti + 1;
	    debug2(printf("**"));
	    /* best_sum = sum; */
	  }
	}
      }
      righti++;
    }
    debug2(printf("\n"));

  }


  if (free_left_p == true) {
#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      FREEA(mismatch_positions_left);
    } else {
      FREE(mismatch_positions_left);
    }
#else
    FREE(mismatch_positions_left);
#endif
  }

  if (free_right_p == true) {
#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      FREEA(mismatch_positions_right);
    } else {
      FREE(mismatch_positions_right);
    }
#else
    FREE(mismatch_positions_right);
#endif
  }

  *best_nmismatches_i = nmismatches_lefti;
  *best_nmismatches_j = nmismatches_righti;

  if (best_sum > nmismatches_allowed) {
    debug2(printf("Returning -1\n"));
    return -1;
#if 0
  } else if (plusp == true) {
    return best_indel_pos;
  } else {
    return querylength - best_indel_pos - indels;
#else
  } else {
    debug2(printf("Returning %d with mismatches %d+%d\n",
		  best_indel_pos,*best_nmismatches_i,*best_nmismatches_j));
    return best_indel_pos;
#endif
  }
}


/* For transcriptome alignments, plusp may be true or false */
/* For alignments via middle_path with ascending univdiags, plusp should be true */
/* indels is negative here */

/* Caller can provide either mismatch_positions_left or
   mismatch_positions_right to save on computation */
int
Indel_resolve_middle_deletion (int *best_nmismatches_i, int *best_nmismatches_j,
			       Univcoord_T left, int indels,
			       int *mismatch_positions_left, int nmismatches_left,
			       int *mismatch_positions_right, int nmismatches_right,
			       Genome_T ome, Genome_T ome_alt, Compress_T query_compress,
			       int querystart, int queryend, int querylength,
			       int nmismatches_allowed, bool plusp, int genestrand,
			       bool want_lowest_coordinate_p) {
  int best_indel_pos = -1, indel_pos;
#ifdef DEBUG2
  char *gbuffer;
  int i;
#endif
  int nmismatches_lefti, nmismatches_righti;
  int best_sum, sum, lefti, righti;
  bool free_left_p = false, free_right_p = false;

  assert(indels < 0);
  if (nmismatches_allowed > querylength) {
    nmismatches_allowed = querylength;
  }

  if (mismatch_positions_left == NULL) {
    free_left_p = true;
#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      mismatch_positions_left = (int *) ALLOCA(querylength * sizeof(int));
    } else {
      mismatch_positions_left = (int *) MALLOC(querylength * sizeof(int));
    }
#else
    mismatch_positions_left = (int *) MALLOC(querylength * sizeof(int));
#endif
    debug2(printf("nmismatches_allowed is %d.  Calling Genome_mismatches_left over %d..%d\n",
		  nmismatches_allowed,querystart,queryend));
    nmismatches_left = Genome_mismatches_left(mismatch_positions_left,nmismatches_allowed,
					      ome,ome_alt,query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
					      plusp,genestrand);
    debug2(
	   printf("%d mismatches on left at:",nmismatches_left);
	   for (i = 0; i <= nmismatches_left; i++) {
	     printf(" %d",mismatch_positions_left[i]);
	   }
	   printf("\n");
	   );
  }

  if (mismatch_positions_right == NULL) {
    free_right_p = true;
#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      mismatch_positions_right = (int *) ALLOCA(querylength * sizeof(int));
    } else {
      mismatch_positions_right = (int *) MALLOC(querylength * sizeof(int));
    }
#else
    mismatch_positions_right = (int *) MALLOC(querylength * sizeof(int));
#endif
    debug2(printf("nmismatches_allowed is %d.  Calling Genome_mismatches_right over %d..%d\n",
		  nmismatches_allowed,querystart,queryend));
    nmismatches_right = Genome_mismatches_right(mismatch_positions_right,nmismatches_allowed,
						ome,ome_alt,query_compress,left-indels,/*pos5*/querystart,/*pos3*/queryend,
						plusp,genestrand);
    debug2(
	   printf("%d mismatches on right at:",nmismatches_right);
	   for (i = 0; i <= nmismatches_right; i++) {
	     printf(" %d",mismatch_positions_right[i]);
	   }
	   printf("\n");
	   );
  }

  /* query has deletion.  Get |indels| more from genome; add to right. */
  /* left = ptr->diagonal - querylength; */

#ifdef DEBUG2
  gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char));
  Genome_fill_buffer_blocks(left,querylength-indels,gbuffer);
  printf("solve_middle_indel, plusp %d, want low %d, deletion (indels %d), nmismatches_allowed %d: Getting genome at left %llu\n",
	 plusp,want_lowest_coordinate_p,indels,nmismatches_allowed,(unsigned long long) left);
  printf("g1: %s\n",gbuffer);
  printf("g2: %s\n",&(gbuffer[-indels]));
  FREE(gbuffer);
#endif


  best_sum = querylength + querylength;

  if (want_lowest_coordinate_p == true) {
    /* Modeled after end C to get lowest possible coordinate */
    righti = 0;
    lefti = nmismatches_left - 1;
    nmismatches_righti = /*righti*/ 0;
    nmismatches_lefti = /*lefti+1*/ nmismatches_left;

    while (righti < nmismatches_right) {
      while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	lefti--;
      }
      sum = righti + lefti + 1;

      debug2(printf("  (Case C1) sum %d=%d+%d at indel_pos %d.",
		    sum,righti,lefti+1,mismatch_positions_right[righti]+1));
      if (sum <= best_sum) {
	indel_pos = mismatch_positions_right[righti] + 1;
	if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	    /* Require separation from endpoints */
	    indel_pos > querystart && indel_pos < queryend) {
	  best_indel_pos = indel_pos;
	  nmismatches_righti = righti;
	  nmismatches_lefti = lefti + 1;
	  debug2(printf("**"));
	  best_sum = sum;
	}
      }
      righti++;
    }
    debug2(printf("\n"));

    /* Try from other side to see if we missed anything */
    lefti = 0;
    righti = nmismatches_right - 1;

    while (lefti < nmismatches_left) {
      while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti]) {
	righti--;
      }
      sum = lefti + righti + 1;

      debug2(printf("  (Case C2) sum %d=%d+%d at indel_pos %d.",
		    sum,lefti,righti+1,mismatch_positions_left[lefti]));
      indel_pos = mismatch_positions_left[lefti];
      if (sum < best_sum ||
	  (sum == best_sum && indel_pos < best_indel_pos)) {
	if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	    /* Require separation from endpoints */
	    indel_pos > querystart && indel_pos < queryend) {
	  best_indel_pos = indel_pos;
	  nmismatches_lefti = lefti;
	  nmismatches_righti = righti + 1;
	  debug2(printf("**"));
	  best_sum = sum;
	}
      }
      lefti++;
    }
    debug2(printf("\n"));

  } else {
    /* Want highest possible coordinate.  Modified from code above, but not sure if this is correct */
    lefti = 0;
    righti = nmismatches_right - 1;
    nmismatches_lefti = /*lefti*/ 0;
    nmismatches_righti = /*righti+1*/ nmismatches_right;

    while (lefti < nmismatches_left) {
      while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti]) {
	righti--;
      }
      sum = lefti + righti + 1;

      debug2(printf("  (Case X2) sum %d=%d+%d at indel_pos %d.",
		    sum,lefti,righti+1,mismatch_positions_left[lefti]));
      if (sum <= best_sum) {
	indel_pos = mismatch_positions_left[lefti];
	if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	    /* Require separation from endpoints */
	    indel_pos > querystart && indel_pos < queryend) {
	  best_indel_pos = indel_pos;
	  nmismatches_lefti = lefti;
	  nmismatches_righti = righti + 1;
	  debug2(printf("**"));
	  best_sum = sum;
	}
      }
      lefti++;
    }
    debug2(printf("\n"));

    /* Try from other side to see if we missed anything */
    righti = 0;
    lefti = nmismatches_left - 1;

    while (righti < nmismatches_right) {
      while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	lefti--;
      }
      sum = righti + lefti + 1;

      debug2(printf("  (Case X1) sum %d=%d+%d at indel_pos %d.",
		    sum,righti,lefti+1,mismatch_positions_right[righti]+1));
      indel_pos = mismatch_positions_right[righti] + 1;
      if (sum < best_sum ||
	  (sum == best_sum && indel_pos > best_indel_pos)) {
	if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	    /* Require separation from endpoints */
	    indel_pos > querystart && indel_pos < queryend) {
	  best_indel_pos = indel_pos;
	  nmismatches_righti = righti;
	  nmismatches_lefti = lefti + 1;
	  debug2(printf("**"));
	  best_sum = sum;
	}
      }
      righti++;
    }
    debug2(printf("\n"));
  }


  if (free_left_p == true) {
#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      FREEA(mismatch_positions_left);
    } else {
      FREE(mismatch_positions_left);
    }
#else
    FREE(mismatch_positions_left);
#endif
  }

  if (free_right_p == true) {
#ifdef HAVE_ALLOCA
    if (querylength <= MAX_STACK_READLENGTH) {
      FREEA(mismatch_positions_right);
    } else {
      FREE(mismatch_positions_right);
    }
#else
    FREE(mismatch_positions_right);
#endif
  }

  *best_nmismatches_i = nmismatches_lefti;
  *best_nmismatches_j = nmismatches_righti;

  if (best_sum > nmismatches_allowed) {
    debug2(printf("Returning -1\n"));
    return -1;
#if 0
  } else if (plusp == true) {
    return best_indel_pos;
  } else {
    return querylength - best_indel_pos;
#else
  } else {
    debug2(printf("Returning %d with nmismatches %d+%d\n",
		  best_indel_pos,*best_nmismatches_i,*best_nmismatches_j));
    return best_indel_pos;
#endif
  }
}


int
Indel_resolve_middle_deletion_or_splice (int *best_introntype,
					 int *best_nmismatches_i, int *best_nmismatches_j,
					 Univcoord_T left, int indels,
					 Genome_T ome, Genome_T ome_alt, Compress_T query_compress,
					 int querystart, int queryend, int querylength,
					 int nmismatches_allowed, bool plusp, int genestrand,
					 int min_intronlength, bool want_lowest_coordinate_p) {
  int best_indel_pos = -1, indel_pos;
  char *gbuffer;
#ifdef DEBUG2
  int i;
#endif
  int nmismatches_left, nmismatches_right, nmismatches_lefti, nmismatches_righti;
  int best_sum, sum, lefti, righti;
  int *mismatch_positions_left, *mismatch_positions_right;
  char left1, left2, right2, right1;
  int introntype, intron_level, best_intron_level;

  *best_introntype = NONINTRON;
  best_intron_level = Intron_level(NONINTRON);

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions_left = (int *) ALLOCA(querylength * sizeof(int));
    mismatch_positions_right = (int *) ALLOCA(querylength * sizeof(int));
  } else {
    mismatch_positions_left = (int *) MALLOC(querylength * sizeof(int));
    mismatch_positions_right = (int *) MALLOC(querylength * sizeof(int));
  }
#else
  mismatch_positions_left = (int *) MALLOC(querylength * sizeof(int));
  mismatch_positions_right = (int *) MALLOC(querylength * sizeof(int));
#endif

  if (nmismatches_allowed > querylength) {
    nmismatches_allowed = querylength;
  }

  /* query has deletion.  Get |indels| more from genome; add to right. */
  /* left = ptr->diagonal - querylength; */

  assert(indels < 0);
#ifdef DEBUG2
  gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char));
  Genome_fill_buffer_blocks(left,querylength-indels,gbuffer);
#else
  if (-indels >= min_intronlength) {
    gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char));
    Genome_fill_buffer_blocks(left,querylength-indels,gbuffer);
  }  
#endif
  debug2(printf("solve_middle_indel, plusp %d, want low %d, deletion (indels %d), nmismatches_allowed %d: Getting genome at diagonal - querylength %d = %llu\n",
		plusp,want_lowest_coordinate_p,indels,nmismatches_allowed,querylength,(unsigned long long) left));
  debug2(printf("g1: %s\n",gbuffer));
  debug2(printf("g2: %s\n",&(gbuffer[-indels])));

  /* No need to check chromosome bounds */
  debug2(printf("nmismatches_allowed is %d.  Calling Genome_mismatches_left over %d..%d\n",
		nmismatches_allowed,querystart,queryend));
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,nmismatches_allowed,
					    ome,ome_alt,query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
					    plusp,genestrand);

  debug2(
	 printf("%d mismatches on left at:",nmismatches_left);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  /* No need to check chromosome bounds */
  debug2(printf("nmismatches_allowed is %d.  Calling Genome_mismatches_right over %d..%d\n",
		nmismatches_allowed,querystart,queryend));
  nmismatches_right = Genome_mismatches_right(mismatch_positions_right,nmismatches_allowed,
					      ome,ome_alt,query_compress,left-indels,/*pos5*/querystart,/*pos3*/queryend,
					      plusp,genestrand);

  debug2(
	 printf("%d mismatches on right at:",nmismatches_right);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

  best_sum = querylength + querylength;

  if (want_lowest_coordinate_p == true) {
    /* Modeled after end C to get lowest possible coordinate */
    righti = 0;
    lefti = nmismatches_left - 1;
    nmismatches_righti = /*righti*/ 0;
    nmismatches_lefti = /*lefti+1*/ nmismatches_left;

    while (righti < nmismatches_right) {
      while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	lefti--;
      }
      sum = righti + lefti + 1;

      if (-indels >= min_intronlength) {
	/* Account for introntype in cases of ties */
	indel_pos = mismatch_positions_right[righti] + 1;
	left1 = gbuffer[indel_pos];
	left2 = gbuffer[indel_pos+1];
	right2 = gbuffer[indel_pos-indels-2];
	right1 = gbuffer[indel_pos-indels-1];
	introntype = Intron_type(left1,left2,right2,right1,left1,left2,right2,right1,/*cdna_direction*/0);
	intron_level = Intron_level(introntype);
	debug2(printf("  (Case C1) sum %d=%d+%d at indel_pos %d (%c%c-%c%c, type %s).",
		      sum,righti,lefti+1,mismatch_positions_right[righti]+1,
		      left1,left2,right2,right1,Intron_type_string(introntype)));
	if (sum < best_sum ||
	    (sum == best_sum && intron_level > best_intron_level)) {
	  if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_righti = righti;
	    nmismatches_lefti = lefti + 1;
	    debug2(printf("**"));
	    best_sum = sum;
	    *best_introntype = introntype;
	    best_intron_level = intron_level;
	  }
	}

      } else {
	debug2(printf("  (Case C1) sum %d=%d+%d at indel_pos %d.",
		      sum,righti,lefti+1,mismatch_positions_right[righti]+1));
	if (sum <= best_sum) {
	  indel_pos = mismatch_positions_right[righti] + 1;
	  if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_righti = righti;
	    nmismatches_lefti = lefti + 1;
	    debug2(printf("**"));
	    best_sum = sum;
	  }
	}
      }
      righti++;
    }
    debug2(printf("\n"));

    /* Try from other side to see if we missed anything */
    lefti = 0;
    righti = nmismatches_right - 1;

    while (lefti < nmismatches_left) {
      while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti]) {
	righti--;
      }
      sum = lefti + righti + 1;

      if (-indels >= min_intronlength) {
	/* Account for introntype in cases of ties */
	indel_pos = mismatch_positions_left[lefti];
	left1 = gbuffer[indel_pos];
	left2 = gbuffer[indel_pos+1];
	right2 = gbuffer[indel_pos-indels-2];
	right1 = gbuffer[indel_pos-indels-1];
	introntype = Intron_type(left1,left2,right2,right1,left1,left2,right2,right1,/*cdna_direction*/0);
	intron_level = Intron_level(introntype);
	debug2(printf("  (Case C2) sum %d=%d+%d at indel_pos %d (%c%c-%c%c).",
		      sum,lefti,righti+1,mismatch_positions_left[lefti],
		      gbuffer[indel_pos],gbuffer[indel_pos+1],gbuffer[indel_pos-indels-2],gbuffer[indel_pos-indels-1]));
	if (sum < best_sum ||
	    (sum == best_sum && intron_level > best_intron_level) ||
	    (sum == best_sum && intron_level == best_intron_level && indel_pos < best_indel_pos)) {
	  if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_lefti = lefti;
	    nmismatches_righti = righti + 1;
	    debug2(printf("**"));
	    best_sum = sum;
	    *best_introntype = introntype;
	    best_intron_level = intron_level;
	  }
	}

      } else {
	debug2(printf("  (Case C2) sum %d=%d+%d at indel_pos %d.",
		      sum,lefti,righti+1,mismatch_positions_left[lefti]));
	indel_pos = mismatch_positions_left[lefti];
	if (sum < best_sum ||
	    (sum == best_sum && indel_pos < best_indel_pos)) {
	  if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_lefti = lefti;
	    nmismatches_righti = righti + 1;
	    debug2(printf("**"));
	    best_sum = sum;
	  }
	}
      }
      lefti++;
    }
    debug2(printf("\n"));

  } else {
    /* Want highest possible coordinate.  Modified from code above, but not sure if this is correct */
    lefti = 0;
    righti = nmismatches_right - 1;
    nmismatches_lefti = /*lefti*/ 0;
    nmismatches_righti = /*righti+1*/ nmismatches_right;

    while (lefti < nmismatches_left) {
      while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti]) {
	righti--;
      }
      sum = lefti + righti + 1;

      if (-indels >= min_intronlength) {
	/* Account for introntype in cases of ties */
	indel_pos = mismatch_positions_left[lefti];
	left1 = gbuffer[indel_pos];
	left2 = gbuffer[indel_pos+1];
	right2 = gbuffer[indel_pos-indels-2];
	right1 = gbuffer[indel_pos-indels-1];
	introntype = Intron_type(left1,left2,right2,right1,left1,left2,right2,right1,/*cdna_direction*/0);
	intron_level = Intron_level(introntype);
	debug2(printf("  (Case X2) sum %d=%d+%d at indel_pos %d (%c%c-%c%c).",
		      sum,lefti,righti+1,mismatch_positions_left[lefti],
		      gbuffer[indel_pos],gbuffer[indel_pos+1],gbuffer[indel_pos-indels-2],gbuffer[indel_pos-indels-1]));
	if (sum < best_sum ||
	    (sum == best_sum && intron_level > best_intron_level)) {
	  if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_lefti = lefti;
	    nmismatches_righti = righti + 1;
	    debug2(printf("**"));
	    best_sum = sum;
	    *best_introntype = introntype;
	    best_intron_level = intron_level;
	  }
	}

      } else {
	debug2(printf("  (Case X2) sum %d=%d+%d at indel_pos %d.",
		      sum,lefti,righti+1,mismatch_positions_left[lefti]));
	if (sum <= best_sum) {
	  indel_pos = mismatch_positions_left[lefti];
	  if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_lefti = lefti;
	    nmismatches_righti = righti + 1;
	    debug2(printf("**"));
	    best_sum = sum;
	  }
	}
      }
      lefti++;
    }
    debug2(printf("\n"));

    /* Try from other side to see if we missed anything */
    righti = 0;
    lefti = nmismatches_left - 1;

    while (righti < nmismatches_right) {
      while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	lefti--;
      }
      sum = righti + lefti + 1;

      if (-indels >= min_intronlength) {
	/* Account for introntype in cases of ties */
	indel_pos = mismatch_positions_right[righti] + 1;
	left1 = gbuffer[indel_pos];
	left2 = gbuffer[indel_pos+1];
	right2 = gbuffer[indel_pos-indels-2];
	right1 = gbuffer[indel_pos-indels-1];
	introntype = Intron_type(left1,left2,right2,right1,left1,left2,right2,right1,/*cdna_direction*/0);
	intron_level = Intron_level(introntype);
	debug2(printf("  (Case X1) sum %d=%d+%d at indel_pos %d (%c%c-%c%c, type %s).",
		      sum,righti,lefti+1,mismatch_positions_right[righti]+1,
		      left1,left2,right2,right1,Intron_type_string(introntype)));
	if (sum < best_sum ||
	    (sum == best_sum && intron_level > best_intron_level) ||
	    (sum == best_sum && intron_level == best_intron_level && indel_pos > best_indel_pos)) {
	  if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_righti = righti;
	    nmismatches_lefti = lefti + 1;
	    debug2(printf("**"));
	    best_sum = sum;
	    *best_introntype = introntype;
	    best_intron_level = intron_level;
	  }
	}

      } else {
	debug2(printf("  (Case X1) sum %d=%d+%d at indel_pos %d.",
		      sum,righti,lefti+1,mismatch_positions_right[righti]+1));
	indel_pos = mismatch_positions_right[righti] + 1;
	if (sum < best_sum ||
	    (sum == best_sum && indel_pos > best_indel_pos)) {
	  if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches &&
	      /* Require separation from endpoints */
	      indel_pos > querystart && indel_pos < queryend) {
	    best_indel_pos = indel_pos;
	    nmismatches_righti = righti;
	    nmismatches_lefti = lefti + 1;
	    debug2(printf("**"));
	    best_sum = sum;
	  }
	}
      }
      righti++;
    }
    debug2(printf("\n"));
  }


#ifdef DEBUG2
  FREE(gbuffer);
#else
  if (-indels >= min_intronlength) {
    FREE(gbuffer);
  }
#endif

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions_left);
    FREEA(mismatch_positions_right);
  } else {
    FREE(mismatch_positions_left);
    FREE(mismatch_positions_right);
  }
#else
  FREE(mismatch_positions_left);
  FREE(mismatch_positions_right);
#endif

  *best_nmismatches_i = nmismatches_lefti;
  *best_nmismatches_j = nmismatches_righti;

  if (best_sum > nmismatches_allowed) {
    debug2(printf("Returning -1\n"));
    return -1;
#if 0
  } else if (plusp == true) {
    return best_indel_pos;
  } else {
    return querylength - best_indel_pos;
#else
  } else {
    debug2(printf("Returning %d with nmismatches %d+%d\n",
		  best_indel_pos,*best_nmismatches_i,*best_nmismatches_j));
    return best_indel_pos;
#endif
  }
}



/* TODO: Implement a faster version that looks only at distal indel_pos */

/* firstbound should be trim5 */
int
Indel_solve_end_low (int *best_adj, int *total_nmismatches, Univcoord_T left, int firstbound,
		     int querylength, Compress_T query_compress,
		     Genome_T omebits, Genome_T omebits_alt,
		     int max_end_insertions, int max_end_deletions,
		     int nmismatches_allowed, bool plusp, int genestrand,
		     bool want_lowest_coordinate_p) {
  int best_indel_pos = -1, indel_pos, breakpoint;
#ifdef DEBUG3
  char *gbuffer;
  int i;
#endif
  int nmismatches_right;
  int *mismatch_positions_right;

  int best_nmismatches_i, best_nmismatches_j;
  int adj;

  *total_nmismatches = querylength;
  *best_adj = 0;
  breakpoint = firstbound;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions_right = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions_right = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions_right = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug3(printf("\nsolve_end_indel_low: Getting genome at left %llu - max_end_deletions %d = %llu.\n",
		(unsigned long long) left,max_end_deletions,(unsigned long long) (left-max_end_deletions)));

#ifdef DEBUG3
  gbuffer = (char *) CALLOC(querylength+max_end_deletions+1,sizeof(char));
  Genome_fill_buffer_blocks(left-max_end_deletions,querylength+max_end_deletions,gbuffer);
  printf("g: %s\n",gbuffer);
  FREE(gbuffer);
#endif

  if (nmismatches_allowed > querylength) {
    nmismatches_allowed = querylength;
  }
  nmismatches_right = Genome_mismatches_right(mismatch_positions_right,nmismatches_allowed,
					      /*ome*/omebits,/*ome_alt*/omebits_alt,query_compress,
					      left,/*pos5*/0,/*pos3*/querylength,
					      plusp,genestrand);

  debug3(
	  printf("full read: %d (max %d) mismatches from right:",nmismatches_right,nmismatches_allowed);
	  for (i = 0; i <= nmismatches_right; i++) {
	    printf(" %d",mismatch_positions_right[i]);
	  }
	  printf("\n");
	  );

  if (max_end_insertions > breakpoint - min_indel_end_matches) {
    max_end_insertions = breakpoint - min_indel_end_matches;
  }

  for (adj = 1; adj <= max_end_deletions; adj++) {
    /* *indels = -adj; */
    if ((indel_pos = Indel_resolve_middle_deletion(&best_nmismatches_i,&best_nmismatches_j,
						   /*left*/left-adj,/*indels*/-adj,
						   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						   mismatch_positions_right,nmismatches_right,
						   omebits,omebits_alt,query_compress,
						   /*querystart*/0,/*queryend*/querylength-adj,querylength,
						   nmismatches_allowed,plusp,genestrand,
						   want_lowest_coordinate_p)) > 0) {
      debug3(printf("Got indel_pos %d.\n",indel_pos));
      if (best_nmismatches_i + best_nmismatches_j < *total_nmismatches) {
	best_indel_pos = indel_pos;
	*best_adj = adj;
	*total_nmismatches = best_nmismatches_i + best_nmismatches_j;
      }
    }
  }

  for (adj = -1; adj >= -max_end_insertions; adj--) {
    /* *indels = -adj; */
    if ((indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
						    /*left*/left-adj,/*indels*/-adj,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    mismatch_positions_right,nmismatches_right,
						    omebits,omebits_alt,query_compress,
						    /*querystart*/0,/*queryend*/querylength+adj,querylength,
						    nmismatches_allowed,plusp,genestrand,
						    want_lowest_coordinate_p)) > 0) {
      debug3(printf("Got indel_pos %d.\n",indel_pos));
      if (best_nmismatches_i + best_nmismatches_j < *total_nmismatches) {
	best_indel_pos = indel_pos;
	*best_adj = adj;
	*total_nmismatches = best_nmismatches_i + best_nmismatches_j;
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions_right);
  } else {
    FREE(mismatch_positions_right);
  }
#else
  FREE(mismatch_positions_right);
#endif

  debug3(printf("Returning best_indel_pos %d with %d mismatches\n\n",best_indel_pos,*total_nmismatches));
  return best_indel_pos;
}


/* lastbound should be querylength - trim3 */
int
Indel_solve_end_high (int *best_adj, int *total_nmismatches, Univcoord_T left, int lastbound,
		      int querylength, Compress_T query_compress,
		      Genome_T omebits, Genome_T omebits_alt,
		      int max_end_insertions, int max_end_deletions,
		      int nmismatches_allowed, bool plusp, int genestrand,
		      bool want_lowest_coordinate_p) {
  int best_indel_pos = -1, indel_pos, breakpoint;
#ifdef DEBUG3
  char *gbuffer;
  int i;
#endif
  int nmismatches_left;
  int *mismatch_positions_left;

  int best_nmismatches_i, best_nmismatches_j;
  int adj;

  *total_nmismatches = querylength;
  *best_adj = 0;
  breakpoint = lastbound - 1;

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions_left = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions_left = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions_left = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug3(printf("\nsolve_end_indel_high: Getting genome at  %llu + max_end_deletions %d = %llu.\n",
		(unsigned long long) left,max_end_deletions,(unsigned long long) (left+max_end_deletions)));

#ifdef DEBUG3
  gbuffer = (char *) CALLOC(querylength+max_end_deletions+1,sizeof(char));
  Genome_fill_buffer_blocks(left,querylength+max_end_deletions,gbuffer);
  printf("g: %s\n",gbuffer);
  FREE(gbuffer);
#endif

  if (nmismatches_allowed > querylength) {
    nmismatches_allowed = querylength;
  }
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,nmismatches_allowed,
					    /*ome*/omebits,/*ome_alt*/omebits_alt,query_compress,
					    left,/*pos5*/0,/*pos3*/querylength,
					    plusp,genestrand);

  debug3(
	  printf("full read: %d (max %d) mismatches from left:",nmismatches_left,nmismatches_allowed);
	  for (i = 0; i <= nmismatches_left; i++) {
	    printf(" %d",mismatch_positions_left[i]);
	  }
	  printf("\n");
	  );

  if (max_end_insertions > breakpoint - min_indel_end_matches) {
    max_end_insertions = breakpoint - min_indel_end_matches;
  }

  for (adj = 1; adj <= max_end_deletions; adj++) {
    if ((indel_pos = Indel_resolve_middle_deletion(&best_nmismatches_i,&best_nmismatches_j,
						   /*left*/left,/*indels*/-adj,
						   mismatch_positions_left,nmismatches_left,
						   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   omebits,omebits_alt,query_compress,
						   /*querystart*/0,/*queryend*/querylength-adj,querylength,
						   nmismatches_allowed,plusp,genestrand,
						   want_lowest_coordinate_p)) > 0) {
      debug3(printf("Got indel_pos %d.\n",indel_pos));
      if (best_nmismatches_i + best_nmismatches_j < *total_nmismatches) {
	best_indel_pos = indel_pos;
	*best_adj = adj;
	*total_nmismatches = best_nmismatches_i + best_nmismatches_j;
      }
    }
  }
  
  for (adj = -1; adj >= -max_end_insertions; adj--) {
    if ((indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
						    /*left*/left,/*indels*/-adj,
						    mismatch_positions_left,nmismatches_left,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    omebits,omebits_alt,query_compress,
						    /*querystart*/0,/*queryend*/querylength+adj,querylength,
						    nmismatches_allowed,plusp,genestrand,
						    want_lowest_coordinate_p)) > 0) {
      debug3(printf("Got indel_pos %d.\n",indel_pos));
      if (best_nmismatches_i + best_nmismatches_j < *total_nmismatches) {
	best_indel_pos = indel_pos;
	*best_adj = adj;
	*total_nmismatches = best_nmismatches_i + best_nmismatches_j;
      }
    }
  }

#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions_left);
  } else {
    FREE(mismatch_positions_left);
  }
#else
  FREE(mismatch_positions_left);
#endif

  debug3(printf("Returning best_indel_pos %d with %d mismatches\n\n",best_indel_pos,*total_nmismatches));
  return best_indel_pos;
}
