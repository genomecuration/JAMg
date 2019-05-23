static char rcsid[] = "$Id: intersect.c 218369 2019-02-13 22:41:06Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "intersect.h"
#include "mem.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Intersect_exact_indices */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Intersect_approx */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Intersect_approx_lower */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Intersect_approx_higher */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif


static int
binary_search_uint4 (int lowi, int highi, UINT4 *positions, UINT4 goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
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

#ifdef LARGE_GENOMES
static int
binary_search_uint8 (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
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
#endif


/* LARGE_GENOMES needs this to handle transcriptomes */
UINT4 *
Intersect_exact (int *ndiagonals,
		 UINT4 *positionsa, int npositionsa, int diagterma,
		 UINT4 *positionsb, int npositionsb, int diagtermb) {
  UINT4 *diagonals, local_goal, last_diagonal, this_diagonal;
  UINT4 diagterm, delta;
  UINT4 *positions0, *positions1;
  int npositions0, npositions1, j;


  *ndiagonals = 0;
  if (npositionsa == 0 || npositionsb == 0) {
    debug(printf("Intersect_exact: intersection is null because npositionsa %d or npositionsb %d is zero\n",npositionsa,npositionsb));
    return (UINT4 *) NULL;

  } else if (npositionsa < npositionsb) {
    positions0 = positionsa;
    positions1 = positionsb;
    npositions0 = npositionsa;
    npositions1 = npositionsb;

    diagterm = (UINT4) diagtermb;	/* local_goal based on larger list */
    delta = (UINT4) (diagterma - diagtermb); /* list0 + (diagterm0 - diagterm1) = list1 */

  } else {
    positions0 = positionsb;
    positions1 = positionsa;
    npositions0 = npositionsb;
    npositions1 = npositionsa;

    diagterm = (UINT4) diagterma;	/* local_goal based on larger list */
    delta = (UINT4) (diagtermb - diagterma); /* list0 + (diagterm0 - diagterm1) = list1 */
  }

  debug(printf("Intersect_exact with %d positions <= %d positions.  diagterm %d\n",
	       npositions0,npositions1,(int) diagterm));

  diagonals = (UINT4 *) MALLOC(npositions0 * sizeof(UINT4));

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && (*positions0) < (Univcoord_T) -(diagtermb - diagterma)) {
      ++positions0;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && (*positions1) < (Univcoord_T) -diagterma) {
      ++positions1;
      --npositions1;
    }
  }
#endif

  while (npositions0 > 0) {
    local_goal = (*positions0) + delta;
    debug(printf("intersection list 0: %d:%llu => local_goal %llu\n",npositions0,*positions0,local_goal));
    if (npositions1 > 0 && *positions1 < local_goal) {
      j = 1;
      while (j < npositions1 && positions1[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search_uint4(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search_uint4(j >> 1,j,positions1,local_goal);
      }
      positions1 += j;
      npositions1 -= j;
    }

    if (npositions1 <= 0) {
      return diagonals;

    } else if ((*positions1) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%llu  found\n",npositions1,*positions1));
      if (*ndiagonals == 0) {
	last_diagonal = diagonals[(*ndiagonals)++] = local_goal + diagterm;
      } else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	last_diagonal = diagonals[(*ndiagonals)++] = this_diagonal;
      }
      ++positions1;
      --npositions1;
    }

    ++positions0;
    --npositions0;
  }
  debug(printf("\n"));

  return diagonals;
}


int
Intersect_exact_indices_univcoord (int *indices,
				   Univcoord_T *positions1, int npositions1,
				   Univcoord_T *positions0, int npositions0) {
  int nindices = 0;
  Univcoord_T local_goal, last_diagonal, this_diagonal;
  /* Univcoord_T diagterm, delta; */
  int j, k0;

  debug1(printf("Intersect_exact_indices_univcoord with %d positions against %d targets\n",
		npositions1,npositions0));

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && (*positions0) < (Univcoord_T) -(diagtermb - diagterma)) {
      ++positions0;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && (*positions1) < (Univcoord_T) -diagterma) {
      ++positions1;
      --npositions1;
    }
  }
#endif

  k0 = 0;
  while (npositions0 > 0) {
    local_goal = (*positions0) /*+ delta*/;
    debug1(printf("intersection list 0: %d:%llu => local_goal %llu\n",npositions0,*positions0,local_goal));
    if (npositions1 > 0 && *positions1 < local_goal) {
      j = 1;
      while (j < npositions1 && positions1[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
#ifdef LARGE_GENOMES
      if (j >= npositions1) {
	j = binary_search_uint8(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search_uint8(j >> 1,j,positions1,local_goal);
      }
#else
      if (j >= npositions1) {
	j = binary_search_uint4(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search_uint4(j >> 1,j,positions1,local_goal);
      }
#endif
      positions1 += j;
      npositions1 -= j;
    }

    if (npositions1 <= 0) {
      return nindices;

    } else if ((*positions1) == local_goal) {
      /* Found local goal.  Save and advance */
      debug1(printf("    intersection list 1: %d:%llu  found => index %d\n",npositions1,*positions1,k0));
      if (nindices == 0) {
	indices[nindices++] = k0;
	last_diagonal = local_goal /*+diagterm*/;
      } else if ((this_diagonal = local_goal /*+ diagterm*/) != last_diagonal) {
	indices[nindices++] = k0;
	last_diagonal = this_diagonal;
      }
      ++positions1;
      --npositions1;
    }

    ++positions0;
    --npositions0;
    k0++;
  }
  debug1(printf("\n"));

  return nindices;
}


#ifndef LARGE_GENOMES
int
Intersect_exact_indices_small (int *indices,
			       UINT4 *positions1, int npositions1, int diagterm1,
			       Univcoord_T *positions0, int npositions0) {
  int nindices = 0;
  Univcoord_T local_goal, last_diagonal, this_diagonal;
  Univcoord_T diagterm, delta;
  int j, k0;

  diagterm = (Univcoord_T) diagterm1;	/* local_goal based on larger list */
  delta = (Univcoord_T) (-diagterm1); /* list0 + (diagterm0 - diagterm1) = list1 */

  debug1(printf("Intersect_exact_indices with %d positions against %d targets.  diagterm %d\n",
	       npositions1,npositions0,(int) diagterm));

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && (*positions0) < (Univcoord_T) -(diagtermb - diagterma)) {
      ++positions0;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && (*positions1) < (Univcoord_T) -diagterma) {
      ++positions1;
      --npositions1;
    }
  }
#endif

  k0 = 0;
  while (npositions0 > 0) {
    local_goal = (*positions0) + delta;
    debug1(printf("intersection list 0: %d:%llu => local_goal %llu\n",npositions0,*positions0,local_goal));
    if (npositions1 > 0 && *positions1 < local_goal) {
      j = 1;
      while (j < npositions1 && positions1[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search_uint4(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search_uint4(j >> 1,j,positions1,local_goal);
      }
      positions1 += j;
      npositions1 -= j;
    }

    if (npositions1 <= 0) {
      return nindices;

    } else if ((*positions1) == local_goal) {
      /* Found local goal.  Save and advance */
      debug1(printf("    intersection list 1: %d:%llu  found => index %d\n",npositions1,*positions1,k0));
      if (nindices == 0) {
	indices[nindices++] = k0;
	last_diagonal = local_goal + diagterm;
      } else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	indices[nindices++] = k0;
	last_diagonal = this_diagonal;
      }
      ++positions1;
      --npositions1;
    }

    ++positions0;
    --npositions0;
    k0++;
  }
  debug1(printf("\n"));

  return nindices;
}
#endif



/* TODO: Write a serial version for small lists */


/* LARGE_GENOMES still needs this to handle transcriptome */
/* Returns results as pairs of coordinates */
UINT4 *
Intersect_approx (bool *exactp, int *ndiagpairs,
		  UINT4 *positionsa, int npositionsa, int diagterma,
		  UINT4 *positionsb, int npositionsb, int diagtermb,
		  Chrpos_T maxdistance) {
  int ndiagonals = 0;
  UINT4 *diagonals, *more_diagonals, local_goal;
  UINT4 diagterm, delta;
  UINT4 *positions0, *positions1;
  UINT4 *ptr1, *start1;
  int npositions0, npositions1, j;
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
  UINT4 last_diagonal = 0U, this_diagonal;
#else
  int guess_allocation;
#endif

  *exactp = false;

  if (npositionsa == 0 || npositionsb == 0) {
    debug(printf("Intersect_approx: intersection is null because npositionsa %d or npositionsb %d is zero\n",npositionsa,npositionsb));
    *ndiagpairs = 0;
    return (UINT4 *) NULL;

  } else if (npositionsa < npositionsb) {
    positions0 = positionsa;
    positions1 = positionsb;
    npositions0 = npositionsa;
    npositions1 = npositionsb;

    diagterm = (UINT4) diagtermb;	/* local_goal based on larger list */
    delta = (UINT4) (diagterma - diagtermb); /* list0 + (diagterm0 - diagterm1) = list1 */

  } else {
    positions0 = positionsb;
    positions1 = positionsa;
    npositions0 = npositionsb;
    npositions1 = npositionsa;

    diagterm = (UINT4) diagterma;	/* local_goal based on larger list */
    delta = (UINT4) (diagtermb - diagterma); /* list0 + (diagterm0 - diagterm1) = list1 */
  }

  start1 = positions1;
  /* end1 = &(positions1[npositions1]); */

  debug(printf("Intersect_approx with %d positions <= %d positions.  diagterm %d\n",
	       npositions0,npositions1,(int) diagterm));
#ifdef DEBUG
  for (ptr1 = positions0; ptr1 < &(positions0[npositions0]); ptr1++) {
    printf("%llu + %d\n",*ptr1,diagterma);
  }
  printf("\n");

  for (ptr1 = positions1; ptr1 < &(positions1[npositions1]); ptr1++) {
    printf("%llu\n",*ptr1);
  }
  printf("\n");
#endif

  ndiagonals = 0;
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
  diagonals = (UINT4 *) MALLOC(2 * npositions0 * sizeof(UINT4));
#else
  guess_allocation = 2 * npositions0;
  diagonals = (UINT4 *) MALLOC(guess_allocation * sizeof(UINT4));
#endif

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && (*positions0) < (Univcoord_T) -(diagtermb - diagterma)) {
      ++positions0;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && (*positions1) < (Univcoord_T) -diagterma) {
      ++positions1;
      --npositions1;
    }
  }
#endif

  while (npositions0 > 0) {
    local_goal = (*positions0) + delta;
    debug(printf("intersection list 0: %d:%llu => local_goal %llu.  Searching from npositions1 %d\n",
		 npositions0,*positions0,local_goal,npositions1));
    if (npositions1 > 0 && *positions1 < local_goal) {
      j = 1;
      while (j < npositions1 && positions1[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search_uint4(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search_uint4(j >> 1,j,positions1,local_goal);
      }
      positions1 += j;
      npositions1 -= j;
    }
#ifdef DEBUG
    if (npositions1 > 0) {
      printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,*positions1);
    }
#endif

    if (npositions1 <= 0) {
      /* Check backwards only */
      debug(printf("    intersection list 1 at end  checking for approximate:"));
      ptr1 = &(positions1[-1]);
      if (ptr1 >= start1) {
	debug(printf(" prev %d:%llu?",npositions1-1,*ptr1));
	if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	  debug(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	  if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	    /* Not a duplicate */
	    diagonals[ndiagonals++] = this_diagonal;
	    diagonals[ndiagonals++] = (*ptr1) + diagterm;
	    last_diagonal = this_diagonal;
	  }
#else
	  /* Want duplicates */
	  if (ndiagonals >= guess_allocation) {
	    more_diagonals = (UINT4 *) MALLOC(2 * guess_allocation * sizeof(UINT4));
	    memcpy(more_diagonals,diagonals,guess_allocation*sizeof(UINT4));
	    FREE(diagonals);
	    diagonals = more_diagonals;
	    guess_allocation = 2 * guess_allocation;
	  }
	  diagonals[ndiagonals++] = local_goal + diagterm;
	  diagonals[ndiagonals++] = (*ptr1) + diagterm;
#endif
	}
      }
      debug(printf("\n"));

      /* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
      *ndiagpairs = ndiagonals / 2; /* Number of pairs */
#ifdef DEBUG
      printf("Returning %d diagpairs\n",*ndiagpairs);
      for (j = 0; j < ndiagonals; j++) {
	printf("%u\n",diagonals[j]);
      }
#endif
      return diagonals;

    } else if ((*positions1) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%llu  exact\n",npositions1,*positions1));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
      if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	/* Not a duplicate */
	diagonals[ndiagonals++] = this_diagonal;
	diagonals[ndiagonals++] = this_diagonal;
	last_diagonal = this_diagonal;
      }
#else
      if (ndiagonals >= guess_allocation) {
	more_diagonals = (UINT4 *) MALLOC(2 * guess_allocation * sizeof(UINT4));
	memcpy(more_diagonals,diagonals,guess_allocation*sizeof(UINT4));
	FREE(diagonals);
	diagonals = more_diagonals;
	guess_allocation = 2 * guess_allocation;
      }
      diagonals[ndiagonals++] = local_goal + diagterm;
      diagonals[ndiagonals++] = local_goal + diagterm;
#endif
      *exactp = true;
      ++positions1;
      --npositions1;

    } else {
      debug(printf("    intersection list 1: %d:%llu  checking for approximate:",npositions1,*positions1));
      ptr1 = &(positions1[-1]); /* closest position < local_goal */
      if (ptr1 >= start1) {
	debug(printf(" prev %d:%llu?",npositions1+1,*ptr1));
	if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	  debug(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	  if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	    /* Not a duplicate */
	    diagonals[ndiagonals++] = this_diagonal;
	    diagonals[ndiagonals++] = (*ptr1) + diagterm;
	    last_diagonal = this_diagonal;
	  }
#else
	  if (ndiagonals >= guess_allocation) {
	    more_diagonals = (UINT4 *) MALLOC(2 * guess_allocation * sizeof(UINT4));
	    memcpy(more_diagonals,diagonals,guess_allocation*sizeof(UINT4));
	    FREE(diagonals);
	    diagonals = more_diagonals;
	    guess_allocation = 2 * guess_allocation;
	  }
	  diagonals[ndiagonals++] = local_goal + diagterm;
	  diagonals[ndiagonals++] = (*ptr1) + diagterm;
#endif
	}
      }

      ptr1 = &(positions1[0]); /* closest position > local_goal */
      debug(printf(" at %d:%llu?",npositions1,*ptr1));
      if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	debug(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	  /* Not a duplicate */
	  diagonals[ndiagonals++] = this_diagonal;
	  diagonals[ndiagonals++] = (*ptr1) + diagterm;
	  last_diagonal = this_diagonal;
	}
#else
	if (ndiagonals >= guess_allocation) {
	  more_diagonals = (UINT4 *) MALLOC(2 * guess_allocation * sizeof(UINT4));
	  memcpy(more_diagonals,diagonals,guess_allocation*sizeof(UINT4));
	  FREE(diagonals);
	  diagonals = more_diagonals;
	  guess_allocation = 2 * guess_allocation;
	}
	diagonals[ndiagonals++] = local_goal + diagterm;
	diagonals[ndiagonals++] = (*ptr1) + diagterm;
#endif
      }
      debug(printf("\n"));
    }

    ++positions0;
    --npositions0;
  }
  debug(printf("\n"));

  *ndiagpairs = ndiagonals / 2; /* Number of pairs */
#ifdef DEBUG
  printf("Returning %d diagpairs\n",*ndiagpairs);
  for (j = 0; j < ndiagonals; j++) {
    printf("%u\n",diagonals[j]);
  }
#endif
  return diagonals;
}



/* Returns results as pairs of coordinates */
/* No diagterms.  Used by both large genomes and regular genomes for
   their diagonals, so type should be Univcoord_T */
Univcoord_T *
Intersect_approx_simple (bool *exactp, int *ndiagpairs,
			 Univcoord_T *positionsa, int npositionsa,
			 Univcoord_T *positionsb, int npositionsb,
			 Chrpos_T maxdistance) {
  int ndiagonals = 0;
  Univcoord_T *diagonals, *more_diagonals, local_goal;
  Univcoord_T *positions0, *positions1;
  Univcoord_T *ptr1, *start1;
  int npositions0, npositions1, j;
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
  Univcoord_T last_diagonal = 0U, this_diagonal;
#else
  int guess_allocation;
#endif

  if (npositionsa == 0 || npositionsb == 0) {
    debug2(printf("Intersect_approx_simple: intersection is null because npositionsa %d or npositionsb %d is zero\n",npositionsa,npositionsb));
    *ndiagpairs = 0;
    return (Univcoord_T *) NULL;

  } else if (npositionsa < npositionsb) {
    positions0 = positionsa;
    positions1 = positionsb;
    npositions0 = npositionsa;
    npositions1 = npositionsb;

  } else {
    positions0 = positionsb;
    positions1 = positionsa;
    npositions0 = npositionsb;
    npositions1 = npositionsa;
  }

  *exactp = false;
  start1 = positions1;
  /* end1 = &(positions1[npositions1]); */

  debug2(printf("Intersect_approx_simple with %d positions <= %d positions, maxdistance %u\n",
	       npositions0,npositions1,maxdistance));
#ifdef DEBUG2
  for (ptr1 = positions0; ptr1 < &(positions0[npositions0]); ptr1++) {
    printf("%llu\n",*ptr1);
  }
  printf("\n");

  for (ptr1 = positions1; ptr1 < &(positions1[npositions1]); ptr1++) {
    printf("%llu\n",*ptr1);
  }
  printf("\n");
#endif

  ndiagonals = 0;
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
  diagonals = (Univcoord_T *) MALLOC(2 * npositions0 * sizeof(Univcoord_T));
#else
  guess_allocation = 2 * npositions0;
  diagonals = (Univcoord_T *) MALLOC(guess_allocation * sizeof(Univcoord_T));
#endif

  while (npositions0 > 0) {
    local_goal = (*positions0);
    debug2(printf("intersection list 0: %d:%llu => local_goal %llu.  Searching from npositions1 %d\n",
		 npositions0,*positions0,local_goal,npositions1));
    if (npositions1 > 0 && *positions1 < local_goal) {
      j = 1;
      while (j < npositions1 && positions1[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
#ifdef LARGE_GENOMES
      if (j >= npositions1) {
	j = binary_search_uint8(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search_uint8(j >> 1,j,positions1,local_goal);
      }
#else
      if (j >= npositions1) {
	j = binary_search_uint4(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search_uint4(j >> 1,j,positions1,local_goal);
      }
#endif
      positions1 += j;
      npositions1 -= j;
    }
#ifdef DEBUG2
    if (npositions1 > 0) {
      printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,*positions1);
    }
#endif

    if (npositions1 <= 0) {
      /* Check backwards only */
      debug2(printf("    intersection list 1 at end  checking for approximate:"));
      ptr1 = &(positions1[-1]);
      if (ptr1 >= start1) {
	debug2(printf(" prev %d:%llu?",npositions1-1,*ptr1));
	if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	  debug2(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	  if ((this_diagonal = local_goal) != last_diagonal) {
	    /* Not a duplicate */
	    diagonals[ndiagonals++] = this_diagonal;
	    diagonals[ndiagonals++] = (*ptr1);
	    last_diagonal = this_diagonal;
	  }
#else
	  /* Want duplicates */
	  if (ndiagonals >= guess_allocation) {
	    more_diagonals = (Univcoord_T *) MALLOC(2 * guess_allocation * sizeof(Univcoord_T));
	    memcpy(more_diagonals,diagonals,guess_allocation*sizeof(Univcoord_T));
	    FREE(diagonals);
	    diagonals = more_diagonals;
	    guess_allocation = 2 * guess_allocation;
	  }
	  diagonals[ndiagonals++] = local_goal;
	  diagonals[ndiagonals++] = (*ptr1);
#endif
	}
      }
      debug2(printf("\n"));

      /* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
      *ndiagpairs = ndiagonals / 2; /* Number of pairs */
      return diagonals;

    } else if ((*positions1) == local_goal) {
      /* Found local goal.  Save and advance */
      debug2(printf("    intersection list 1: %d:%llu  exact\n",npositions1,*positions1));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
      if ((this_diagonal = local_goal) != last_diagonal) {
	/* Not a duplicate */
	diagonals[ndiagonals++] = this_diagonal;
	diagonals[ndiagonals++] = this_diagonal;
	last_diagonal = this_diagonal;
      }
#else
      if (ndiagonals >= guess_allocation) {
	more_diagonals = (Univcoord_T *) MALLOC(2 * guess_allocation * sizeof(Univcoord_T));
	memcpy(more_diagonals,diagonals,guess_allocation*sizeof(Univcoord_T));
	FREE(diagonals);
	diagonals = more_diagonals;
	guess_allocation = 2 * guess_allocation;
      }
      diagonals[ndiagonals++] = local_goal;
      diagonals[ndiagonals++] = local_goal;
#endif
      *exactp = true;
      ++positions1;
      --npositions1;

    } else {
      debug2(printf("    intersection list 1: %d:%llu  checking for approximate:",npositions1,*positions1));
      ptr1 = &(positions1[-1]); /* closest position < local_goal */
      if (ptr1 >= start1) {
	debug2(printf(" prev %d:%llu?",npositions1+1,*ptr1));
	if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	  debug2(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	  if ((this_diagonal = local_goal) != last_diagonal) {
	    /* Not a duplicate */
	    diagonals[ndiagonals++] = this_diagonal;
	    diagonals[ndiagonals++] = (*ptr1);
	    last_diagonal = this_diagonal;
	  }
#else
	  if (ndiagonals >= guess_allocation) {
	    more_diagonals = (Univcoord_T *) MALLOC(2 * guess_allocation * sizeof(Univcoord_T));
	    memcpy(more_diagonals,diagonals,guess_allocation*sizeof(Univcoord_T));
	    FREE(diagonals);
	    diagonals = more_diagonals;
	    guess_allocation = 2 * guess_allocation;
	  }
	  diagonals[ndiagonals++] = local_goal;
	  diagonals[ndiagonals++] = (*ptr1);
#endif
	}
      }

      ptr1 = &(positions1[0]); /* closest position > local_goal */
      debug2(printf(" at %d:%llu?",npositions1,*ptr1));
      if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	debug2(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	if ((this_diagonal = local_goal) != last_diagonal) {
	  /* Not a duplicate */
	  diagonals[ndiagonals++] = this_diagonal;
	  diagonals[ndiagonals++] = (*ptr1);
	  last_diagonal = this_diagonal;
	}
#else
	if (ndiagonals >= guess_allocation) {
	  more_diagonals = (Univcoord_T *) MALLOC(2 * guess_allocation * sizeof(Univcoord_T));
	  memcpy(more_diagonals,diagonals,guess_allocation*sizeof(Univcoord_T));
	  FREE(diagonals);
	  diagonals = more_diagonals;
	  guess_allocation = 2 * guess_allocation;
	}
	diagonals[ndiagonals++] = local_goal;
	diagonals[ndiagonals++] = (*ptr1);
#endif
      }
      debug2(printf("\n"));
    }

    ++positions0;
    --npositions0;
  }
  debug2(printf("\n"));

  *ndiagpairs = ndiagonals / 2; /* Number of pairs */
  return diagonals;
}


#ifndef LARGE_GENOMES
/* diagonals is already allocated by caller */
/* Needs to return diagonals in ascending order.  Need to check if diagonals0 has duplicates. */
int
Intersect_approx_lower (Univcoord_T *diagonals,
			UINT4 *positions1, int npositions1, int diagterm1,
			Univcoord_T *positions0, int npositions0,
			Chrpos_T maxdistance) {
  int ndiagonals = 0;
  Univcoord_T local_goal, last_positions0;
  Univcoord_T diagterm, delta;
  Univcoord_T last_diagonal, this_diagonal;
  UINT4 *ptr1, *start1;
  int j;

  diagterm = (Univcoord_T) diagterm1;	/* local_goal based on larger list */
  delta = (Univcoord_T) (-diagterm1); /* list0 + (diagterm0 - diagterm1) = list1 */

  start1 = positions1;
  /* end1 = &(positions1[npositions1]); */

  debug3(printf("Intersect_approx_lower with %d positions against %d positions.  diagterm %d\n",
	       npositions1,npositions0,(int) diagterm));
#ifdef DEBUG3
  for (ptr1 = positions0; ptr1 < &(positions0[npositions0]); ptr1++) {
    printf("%llu\n",*ptr1);
  }
  printf("\n");

  for (ptr1 = positions1; ptr1 < &(positions1[npositions1]); ptr1++) {
    printf("%llu + %d\n",*ptr1,diagterm1);
  }
  printf("\n");
#endif

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagterm1 < 0) {
    while (npositions1 > 0 && (*positions1) < (Univcoord_T) -diagterm1) {
      ++positions1;
      --npositions1;
    }
  }
#endif

  last_positions0 = (Univcoord_T) -1;
  while (npositions0 > 0) {
    if (*positions0 == last_positions0) {
      /* Skip duplicate in positions0 */
      /* last_positions0 = *positions0 */
      ++positions0;
      --npositions0;

    } else {
      local_goal = (*positions0) + delta;
      debug3(printf("intersection list 0: %d:%llu => local_goal %llu.  Searching from npositions1 %d\n",
		    npositions0,*positions0,local_goal,npositions1));
      if (npositions1 > 0 && *positions1 < local_goal) {
	j = 1;
	while (j < npositions1 && positions1[j] < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= npositions1) {
	  j = binary_search_uint4(j >> 1,npositions1,positions1,local_goal);
	} else {
	  j = binary_search_uint4(j >> 1,j,positions1,local_goal);
	}
	positions1 += j;
	npositions1 -= j;
      }
#ifdef DEBUG3
      if (npositions1 > 0) {
	printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,*positions1);
      }
#endif

      if (npositions1 <= 0) {
	/* Check backwards only */
	debug3(printf("    intersection list 1 at end  checking for approximate:"));
	ptr1 = &(positions1[-1]);
	if (ptr1 >= start1) {
	  debug3(printf(" prev %d:%llu?",npositions1-1,*ptr1));
	  if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal /*(omit for lower) + maxdistance*/) {
	    debug3(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	    } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	debug3(printf("\n"));

	/* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
#ifdef DEBUG3
	printf("Returning %d diagonals\n",ndiagonals);
	for (j = 0; j < ndiagonals; j++) {
	  printf("%u\n",diagonals[j]);
	}
#endif
	return ndiagonals;

      } else if ((*positions1) == local_goal) {
	/* Found local goal.  Save and advance */
	debug3(printf("    intersection list 1: %d:%llu  exact\n",npositions1,*positions1));
	/* diagonals[ndiagonals++] = local_goal + diagterm; */
	if (ndiagonals == 0) {
	  last_diagonal = diagonals[ndiagonals++] = local_goal + diagterm;
	} else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	  last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	}
	++positions1;
	--npositions1;

      } else {
	debug3(printf("    intersection list 1: %d:%llu  checking for approximate:",npositions1,*positions1));
	ptr1 = &(positions1[-1]); /* closest position < local_goal */
	if (ptr1 >= start1) {
	  debug3(printf(" prev %d:%llu?",npositions1+1,*ptr1));
	  if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal /*(omit for lower) + maxdistance*/) {
	    debug3(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	    } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}

	ptr1 = &(positions1[0]); /* closest position > local_goal */
	debug3(printf(" at %d:%llu?",npositions1,*ptr1));
	if ((*ptr1) + maxdistance >= local_goal && (*ptr1) <= local_goal /*(omit for lower) + maxdistance*/) {
	  debug3(printf(" yes."));
	  /* diagonals[ndiagonals++] = local_goal + diagterm; */
	  if (ndiagonals == 0) {
	    last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	  } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
	    last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	  }
	}
	debug3(printf("\n"));
      }

      last_positions0 = *positions0;
      ++positions0;
      --npositions0;
    }
  }
  debug3(printf("\n"));

#ifdef DEBUG3
  printf("Returning %d diagonals\n",ndiagonals);
  for (j = 0; j < ndiagonals; j++) {
    printf("%u\n",diagonals[j]);
  }
#endif
  return ndiagonals;
}
#endif


#ifndef LARGE_GENOMES
/* diagonals is already allocated by caller */
/* Needs to return diagonals in ascending order.  Need to check if diagonals0 has duplicates. */
int
Intersect_approx_higher (Univcoord_T *diagonals,
			 UINT4 *positions1, int npositions1, int diagterm1,
			 Univcoord_T *positions0, int npositions0,
			 Chrpos_T maxdistance) {
  int ndiagonals = 0;
  Univcoord_T local_goal, last_positions0;
  Univcoord_T diagterm, delta;
  Univcoord_T last_diagonal, this_diagonal;
  UINT4 *ptr1, *start1;
  int j;

  diagterm = (Univcoord_T) diagterm1;	/* local_goal based on larger list */
  delta = (Univcoord_T) (-diagterm1); /* list0 + (diagterm0 - diagterm1) = list1 */

  start1 = positions1;
  /* end1 = &(positions1[npositions1]); */

  debug4(printf("Intersect_approx_higher with %d positions against %d positions.  diagterm %d\n",
	       npositions1,npositions0,(int) diagterm));
#ifdef DEBUG4
  for (ptr1 = positions0; ptr1 < &(positions0[npositions0]); ptr1++) {
    printf("%llu\n",*ptr1);
  }
  printf("\n");

  for (ptr1 = positions1; ptr1 < &(positions1[npositions1]); ptr1++) {
    printf("%llu + %d\n",*ptr1,diagterm1);
  }
  printf("\n");
#endif


#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagterm1 < 0) {
    while (npositions1 > 0 && (*positions1) < (Univcoord_T) -diagterm1) {
      ++positions1;
      --npositions1;
    }
  }
#endif

  last_positions0 = (Univcoord_T) -1;
  while (npositions0 > 0) {
    if (*positions0 == last_positions0) {
      /* Skip duplicate in positions0 */
      /* last_positions0 = *positions0 */
      ++positions0;
      --npositions0;

    } else {
      local_goal = (*positions0) + delta;
      debug4(printf("intersection list 0: %d:%llu => local_goal %llu.  Searching from npositions1 %d\n",
		    npositions0,*positions0,local_goal,npositions1));
      if (npositions1 > 0 && *positions1 < local_goal) {
	j = 1;
	while (j < npositions1 && positions1[j] < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= npositions1) {
	  j = binary_search_uint4(j >> 1,npositions1,positions1,local_goal);
	} else {
	  j = binary_search_uint4(j >> 1,j,positions1,local_goal);
	}
	positions1 += j;
	npositions1 -= j;
      }
#ifdef DEBUG4
      if (npositions1 > 0) {
	printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,*positions1);
      }
#endif
      
      if (npositions1 <= 0) {
	/* Check backwards only */
	debug4(printf("    intersection list 1 at end  checking for approximate:"));
	ptr1 = &(positions1[-1]);
	if (ptr1 >= start1) {
	  debug4(printf(" prev %d:%llu?",npositions1-1,*ptr1));
	  if ((*ptr1) /*(omit for higher) + maxdistance*/ >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	    debug4(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; -- For diagpairs */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	    } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	debug4(printf("\n"));
	
	/* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
#ifdef DEBUG4
	printf("Returning %d diagonals\n",ndiagonals);
	for (j = 0; j < ndiagonals; j++) {
	  printf("%u\n",diagonals[j]);
	}
#endif
	return ndiagonals;
	
      } else if ((*positions1) == local_goal) {
	/* Found local goal.  Save and advance */
	debug4(printf("    intersection list 1: %d:%llu  exact\n",npositions1,*positions1));
	/* diagonals[ndiagonals++] = local_goal + diagterm; -- For diagpairs */
	if (ndiagonals == 0) {
	  last_diagonal = diagonals[ndiagonals++] = local_goal + diagterm;
	} else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	  last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	}
	++positions1;
	--npositions1;
	
      } else {
	debug4(printf("    intersection list 1: %d:%llu  checking for approximate:",npositions1,*positions1));
	ptr1 = &(positions1[-1]); /* closest position < local_goal */
	if (ptr1 >= start1) {
	  debug4(printf(" prev %d:%llu?",npositions1+1,*ptr1));
	  if ((*ptr1) /*(omit for higher) + maxdistance*/ >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	    debug4(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; -- For diagpairs */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	    } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	
	ptr1 = &(positions1[0]); /* closest position > local_goal */
	debug4(printf(" at %d:%llu?",npositions1,*ptr1));
	if ((*ptr1) /*(omit for higher) + maxdistance*/ >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	  debug4(printf(" yes."));
	  /* diagonals[ndiagonals++] = local_goal + diagterm; -- For diagpairs */
	  if (ndiagonals == 0) {
	    last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	  } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
	    last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	  }
	}
	debug4(printf("\n"));
      }

      last_positions0 = *positions0;
      ++positions0;
      --npositions0;
    }
  }
  debug4(printf("\n"));

#ifdef DEBUG4
  printf("Returning %d diagonals\n",ndiagonals);
  for (j = 0; j < ndiagonals; j++) {
    printf("%u\n",diagonals[j]);
  }
#endif

  return ndiagonals;
}
#endif


