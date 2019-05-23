static char rcsid[] = "$Id: intersect-large.c 218369 2019-02-13 22:41:06Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "intersect-large.h"
#include "mem.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Approx */
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


#define GETPOS(high,low) (((Univcoord_T) high << 32) + low)


static int
binary_search_large (int lowi, int highi, unsigned char *positions_high, UINT4 *positions_low, Univcoord_T goal) {
  int middlei;
  Univcoord_T position;

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    position = ((Univcoord_T) positions_high[middlei] << 32) + positions_low[middlei];
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,(positions_high[lowi] << 32) + positions_low[lowi],
		   middlei,position,
		   highi,(positions_high[highi] << 32) + positions_low[highi],goal));
    if (goal < position) {
      highi = middlei;
    } else if (goal > position) {
      lowi = middlei + 1;
    } else {
      return middlei;
    }
  }

  return highi;
}


Univcoord_T *
Intersect_exact_large (int *ndiagonals,
		       unsigned char *positionsa_high, UINT4 *positionsa_low,
		       int npositionsa, int diagterma,
		       unsigned char *positionsb_high, UINT4 *positionsb_low,
		       int npositionsb, int diagtermb) {
  Univcoord_T *diagonals, local_goal, last_diagonal, this_diagonal;
  Univcoord_T diagterm, delta;
  unsigned char *positions0_high, *positions1_high;
  UINT4 *positions0_low, *positions1_low;
  int npositions0, npositions1, j;


  *ndiagonals = 0;
  if (npositionsa == 0 || npositionsb == 0) {
    debug(printf("Intersect_exact_large: intersection is null because npositionsa %d or npositionsb %d is zero\n",npositionsa,npositionsb));
    return (Univcoord_T *) NULL;

  } else if (npositionsa < npositionsb) {
    positions0_high = positionsa_high;
    positions0_low = positionsa_low;
    positions1_high = positionsb_high;
    positions1_low = positionsb_low;
    npositions0 = npositionsa;
    npositions1 = npositionsb;

    diagterm = (Univcoord_T) diagtermb;	/* local_goal based on larger list */
    delta = (Univcoord_T) (diagterma - diagtermb); /* list0 + (diagterm0 - diagterm1) = list1 */

  } else {
    positions0_high = positionsb_high;
    positions0_low = positionsb_low;
    positions1_high = positionsa_high;
    positions1_low = positionsa_low;
    npositions0 = npositionsb;
    npositions1 = npositionsa;

    diagterm = (Univcoord_T) diagterma;	/* local_goal based on larger list */
    delta = (Univcoord_T) (diagtermb - diagterma); /* list0 + (diagterm0 - diagterm1) = list1 */
  }

  debug(printf("Intersect_exact_large with %d positions <= %d positions.  diagterm %d\n",
	       npositions0,npositions1,(int) diagterm));

  diagonals = (Univcoord_T *) MALLOC(npositions0 * sizeof(Univcoord_T));

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && GETPOS(*positions0_high,*positions0_low) < (Univcoord_T) -(diagtermb - diagterma)) {
      ++positions0_high;
      ++positions0_low;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < (Univcoord_T) -diagterma) {
      ++positions1_high;
      ++positions1_low;
      --npositions1;
    }
  }
#endif

  while (npositions0 > 0) {
    local_goal = GETPOS(*positions0_high,*positions0_low) + delta;
    debug(printf("intersection list 0: %d:%llu => local_goal %llu\n",
		 npositions0,GETPOS(*positions0_high,*positions0_low),local_goal));
    if (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < local_goal) {
      j = 1;
      while (j < npositions1 && GETPOS(positions1_high[j],positions1_low[j]) < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search_large(j >> 1,npositions1,positions1_high,positions1_low,local_goal);
      } else {
	j = binary_search_large(j >> 1,j,positions1_high,positions1_low,local_goal);
      }
      positions1_high += j;
      positions1_low += j;
      npositions1 -= j;
    }

    if (npositions1 <= 0) {
      return diagonals;

    } else if (GETPOS(*positions1_high,*positions1_low) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%llu  found\n",
		   npositions1,GETPOS(*positions1_high,*positions1_low)));
      if (*ndiagonals == 0) {
	last_diagonal = diagonals[(*ndiagonals)++] = local_goal + diagterm;
      } else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	last_diagonal = diagonals[(*ndiagonals)++] = this_diagonal;
      }
      ++positions1_high;
      ++positions1_low;
      --npositions1;
    }

    ++positions0_high;
    ++positions0_low;
    --npositions0;
  }
  debug(printf("\n"));

  return diagonals;
}


int
Intersect_exact_indices_large (int *indices,
			       unsigned char *positions1_high, UINT4 *positions1_low,
			       int npositions1, int diagterm1,
			       Univcoord_T *positions0, int npositions0) {
  int nindices = 0;
  Univcoord_T local_goal, last_diagonal, this_diagonal;
  Univcoord_T diagterm, delta;
  int j, k0;

  diagterm = (Univcoord_T) diagterm1;	/* local_goal based on larger list */
  delta = (Univcoord_T) (-diagterm1); /* list0 + (diagterm0 - diagterm1) = list1 */

  debug(printf("Intersect_exact_indices_large with %d positions against %d positions.  diagterm %d\n",
	       npositions1,npositions0,(int) diagterm));

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && GETPOS(*positions0_high,*positions0_low) < (Univcoord_T) -(diagtermb - diagterma)) {
      ++positions0_high;
      ++positions0_low;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < (Univcoord_T) -diagterma) {
      ++positions1_high;
      ++positions1_low;
      --npositions1;
    }
  }
#endif

  k0 = 0;
  while (npositions0 > 0) {
    local_goal = (*positions0) + delta;
    debug(printf("intersection list 0: %d:%llu => local_goal %llu\n",
		 npositions0,*positions0,local_goal));
    if (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < local_goal) {
      j = 1;
      while (j < npositions1 && GETPOS(positions1_high[j],positions1_low[j]) < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search_large(j >> 1,npositions1,positions1_high,positions1_low,local_goal);
      } else {
	j = binary_search_large(j >> 1,j,positions1_high,positions1_low,local_goal);
      }
      positions1_high += j;
      positions1_low += j;
      npositions1 -= j;
    }

    if (npositions1 <= 0) {
      return nindices;

    } else if (GETPOS(*positions1_high,*positions1_low) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%llu  found\n",
		   npositions1,GETPOS(*positions1_high,*positions1_low)));
      if (nindices == 0) {
	indices[nindices++] = k0;
	last_diagonal = local_goal + diagterm;
      } else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	indices[nindices++] = k0;
	last_diagonal = this_diagonal;
      }
      ++positions1_high;
      ++positions1_low;
      --npositions1;
    }

    ++positions0;
    --npositions0;
    k0++;
  }
  debug(printf("\n"));

  return nindices;
}


/* Returns results as pairs of coordinates */
Univcoord_T *
Intersect_approx_large (bool *exactp, int *ndiagpairs,
			unsigned char *positionsa_high, UINT4 *positionsa_low,
			int npositionsa, int diagterma,
			unsigned char *positionsb_high, UINT4 *positionsb_low,
			int npositionsb, int diagtermb,
			Chrpos_T maxdistance) {
  int ndiagonals = 0;
  Univcoord_T *diagonals, *more_diagonals, local_goal;
  Univcoord_T diagterm, delta;
  unsigned char *positions0_high, *positions1_high, *ptr1_high, *start1;
  UINT4 *positions0_low, *positions1_low, *ptr1_low;
  int npositions0, npositions1, j;
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
  Univcoord_T last_diagonal = 0U, this_diagonal;
#else
  int guess_allocation;
#endif

  *exactp = false;

  if (npositionsa == 0 || npositionsb == 0) {
    debug(printf("Intersect_approx_large: intersection is null because npositionsa %d or npositionsb %d is zero\n",npositionsa,npositionsb));
    *ndiagpairs = 0;
    return (Univcoord_T *) NULL;

  } else if (npositionsa < npositionsb) {
    positions0_high = positionsa_high;
    positions0_low = positionsa_low;
    positions1_high = positionsb_high;
    positions1_low = positionsb_low;
    npositions0 = npositionsa;
    npositions1 = npositionsb;

    diagterm = (Univcoord_T) diagtermb;	/* local_goal based on larger list */
    delta = (Univcoord_T) (diagterma - diagtermb); /* list0 + (diagterm0 - diagterm1) = list1 */

  } else {
    positions0_high = positionsb_high;
    positions0_low = positionsb_low;
    positions1_high = positionsa_high;
    positions1_low = positionsa_low;
    npositions0 = npositionsb;
    npositions1 = npositionsa;

    diagterm = (Univcoord_T) diagterma;	/* local_goal based on larger list */
    delta = (Univcoord_T) (diagtermb - diagterma); /* list0 + (diagterm0 - diagterm1) = list1 */
  }

  start1 = positions1_high;
  /* end1 = &(positions1_high[npositions1]); */

  debug(printf("Intersect_approx_large with %d positions <= %d positions.  diagterm %d\n",
	       npositions0,npositions1,(int) diagterm));
#ifdef DEBUG
  for (ptr1_high = positions0_high, ptr1_low = positions0_low;
       ptr1_high < &(positions0_high[npositions0]); ptr1_high++, ptr1_low++) {
    printf("%llu + %d\n",GETPOS(*ptr1_high,*ptr1_low),diagterma);
  }
  printf("\n");

  for (ptr1_high = positions1_high, ptr1_low = positions1_low;
       ptr1_high < &(positions1_high[npositions1]); ptr1_high++, ptr1_low++) {
    printf("%llu + %d\n",GETPOS(*ptr1_high,*ptr1_low),diagterma);
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

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && GETPOS(*positions0_high,*positions0_low) < (Univcoord_T) -(diagtermb - diagterma)) {
      ++positions0_high;
      ++positions0_low;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < (Univcoord_T) -diagterma) {
      ++positions1_high;
      ++positions1_low;
      --npositions1;
    }
  }
#endif

  while (npositions0 > 0) {
    local_goal = GETPOS(*positions0_high,*positions0_low) + delta;
    debug(printf("intersection list 0: %d:%llu => local_goal %llu.  Searching from npositions1 %d\n",
		 npositions0,GETPOS(*positions0_high,*positions0_low),local_goal,npositions1));
    if (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < local_goal) {
      j = 1;
      while (j < npositions1 && GETPOS(positions1_high[j],positions1_low[j]) < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search_large(j >> 1,npositions1,positions1_high,positions1_low,local_goal);
      } else {
	j = binary_search_large(j >> 1,j,positions1_high,positions1_low,local_goal);
      }
      positions1_high += j;
      positions1_low += j;
      npositions1 -= j;
    }
#ifdef DEBUG
    if (npositions1 > 0) {
      printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,GETPOS(*positions1_high,*positions1_low));
    }
#endif

    if (npositions1 <= 0) {
      /* Check backwards only */
      debug(printf("    intersection list 1 at end  checking for approximate:"));
      ptr1_high = &(positions1_high[-1]);
      ptr1_low = &(positions1_low[-1]);
      if (ptr1_high >= start1) {
	debug(printf(" prev %d:%llu?",npositions1-1,GETPOS(*ptr1_high,*ptr1_low)));
	if (GETPOS(*ptr1_high,*ptr1_low) + maxdistance >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	  debug(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	  if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	    /* Not a duplicate */
	    diagonals[ndiagonals++] = this_diagonal;
	    diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
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
	  diagonals[ndiagonals++] = local_goal + diagterm;
	  diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
#endif
	}
      }
      debug(printf("\n"));

      /* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
      *ndiagpairs = ndiagonals / 2; /* Number of pairs */
      return diagonals;

    } else if (GETPOS(*positions1_high,*positions1_low) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%llu  exact\n",
		   npositions1,GETPOS(*positions1_high,*positions1_low)));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
      if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
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
      diagonals[ndiagonals++] = local_goal + diagterm;
      diagonals[ndiagonals++] = local_goal + diagterm;
#endif
      *exactp = true;
      ++positions1_high;
      ++positions1_low;
      --npositions1;

    } else {
  debug(printf("    intersection list 1: %d:%llu  checking for approximate:",
	       npositions1,GETPOS(*positions1_high,*positions1_low)));
      ptr1_high = &(positions1_high[-1]); /* closest position < local_goal */
      ptr1_low = &(positions1_low[-1]); /* closest position < local_goal */
      if (ptr1_high >= start1) {
	debug(printf(" prev %d:%llu?",npositions1+1,GETPOS(*ptr1_high,*ptr1_low)));
	if (GETPOS(*ptr1_high,*ptr1_low) + maxdistance >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	  debug(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	  if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	    /* Not a duplicate */
	    diagonals[ndiagonals++] = this_diagonal;
	    diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
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
	  diagonals[ndiagonals++] = local_goal + diagterm;
	  diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
#endif
	}
      }

      ptr1_high = &(positions1_high[0]); /* closest position > local_goal */
      ptr1_low = &(positions1_low[0]); /* closest position > local_goal */
      debug(printf(" at %d:%llu?",npositions1,GETPOS(*ptr1_high,*ptr1_low)));
      if (GETPOS(*ptr1_high,*ptr1_low) + maxdistance >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	debug(printf(" yes."));
#ifdef REMOVE_DUPLICATES_OF_FIRST_DIAGONAL
	if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	  /* Not a duplicate */
	  diagonals[ndiagonals++] = this_diagonal;
	  diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
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
	diagonals[ndiagonals++] = local_goal + diagterm;
	diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
#endif
      }
      debug(printf("\n"));
    }

    ++positions0_high;
    ++positions0_low;
    --npositions0;
  }
  debug(printf("\n"));

  *ndiagpairs = ndiagonals / 2; /* Number of pairs */
  return diagonals;
}


/* diagonals is already allocated by caller */
/* Needs to return diagonals in ascending order.  Need to check if diagonals0 has duplicates. */
int
Intersect_approx_lower (Univcoord_T *diagonals,
			unsigned char *positions1_high, UINT4 *positions1_low,
			int npositions1, int diagterm1,
			Univcoord_T *positions0, int npositions0,
			Chrpos_T maxdistance) {
  int ndiagonals = 0;
  Univcoord_T local_goal, last_positions0;
  Univcoord_T diagterm, delta;
  Univcoord_T last_diagonal, this_diagonal;
  unsigned char *ptr1_high, *start1;
  UINT4 *ptr1_low;
  int j;
#ifdef DEBUG
  Univcoord_T *ptr0;
#endif


  diagterm = (Univcoord_T) diagterm1;	/* local_goal based on larger list */
  delta = (Univcoord_T) (-diagterm1); /* list0 + (diagterm0 - diagterm1) = list1 */

  start1 = positions1_high;
  /* end1 = &(positions1_high[npositions1]); */

  debug(printf("Intersect_approx_lower with %d positions against %d positions.  diagterm %d\n",
	       npositions1,npositions0,(int) diagterm));
#ifdef DEBUG
  for (ptr0 = positions0;  ptr0 < &(positions0[npositions0]); ptr0++) {
    printf("%llu\n",*ptr0);
  }
  printf("\n");

  for (ptr1_high = positions1_high, ptr1_low = positions1_low;
       ptr1_high < &(positions1_high[npositions1]); ptr1_high++, ptr1_low++) {
    printf("%llu + %d\n",GETPOS(*ptr1_high,*ptr1_low),diagterm1);
  }
  printf("\n");
#endif


#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagterm1 < 0) {
    while (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < (Univcoord_T) -diagterma) {
      ++positions1_high;
      ++positions1_low;
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
      debug(printf("intersection list 0: %d:%llu => local_goal %llu.  Searching from npositions1 %d\n",
		   npositions0,*positions0,local_goal,npositions1));
      if (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < local_goal) {
	j = 1;
	while (j < npositions1 && GETPOS(positions1_high[j],positions1_low[j]) < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= npositions1) {
	  j = binary_search_large(j >> 1,npositions1,positions1_high,positions1_low,local_goal);
	} else {
	  j = binary_search_large(j >> 1,j,positions1_high,positions1_low,local_goal);
	}
	positions1_high += j;
	positions1_low += j;
	npositions1 -= j;
      }
#ifdef DEBUG
      if (npositions1 > 0) {
	printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,GETPOS(*positions1_high,*positions1_low));
      }
#endif

      if (npositions1 <= 0) {
	/* Check backwards only */
	debug(printf("    intersection list 1 at end  checking for approximate:"));
	ptr1_high = &(positions1_high[-1]);
	ptr1_low = &(positions1_low[-1]);
	if (ptr1_high >= start1) {
	  debug(printf(" prev %d:%llu?",npositions1-1,GETPOS(*ptr1_high,*ptr1_low)));
	  if (GETPOS(*ptr1_high,*ptr1_low) + maxdistance >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal /*(omit for lower) + maxdistance*/) {
	    debug(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	    } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	debug(printf("\n"));
	
	/* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
	return ndiagonals;
	
      } else if (GETPOS(*positions1_high,*positions1_low) == local_goal) {
	/* Found local goal.  Save and advance */
	debug(printf("    intersection list 1: %d:%llu  exact\n",
		     npositions1,GETPOS(*positions1_high,*positions1_low)));
	/* diagonals[ndiagonals++] = local_goal + diagterm; */
	if (ndiagonals == 0) {
	  last_diagonal = diagonals[ndiagonals++] = local_goal + diagterm;
	} else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	  last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	}
	++positions1_high;
	++positions1_low;
	--npositions1;
	
      } else {
	debug(printf("    intersection list 1: %d:%llu  checking for approximate:",
		     npositions1,GETPOS(*positions1_high,*positions1_low)));
	ptr1_high = &(positions1_high[-1]); /* closest position < local_goal */
	ptr1_low = &(positions1_low[-1]); /* closest position < local_goal */
	if (ptr1_high >= start1) {
	  debug(printf(" prev %d:%llu?",npositions1+1,GETPOS(*ptr1_high,*ptr1_low)));
	  if (GETPOS(*ptr1_high,*ptr1_low) + maxdistance >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal /*(omit for lower) + maxdistance*/) {
	    debug(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	    } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}

	ptr1_high = &(positions1_high[0]); /* closest position > local_goal */
	ptr1_low = &(positions1_low[0]); /* closest position > local_goal */
	debug(printf(" at %d:%llu?",npositions1,GETPOS(*ptr1_high,*ptr1_low)));
	if (GETPOS(*ptr1_high,*ptr1_low) + maxdistance >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal /*(omit for lower) + maxdistance*/) {
	  debug(printf(" yes."));
	  /* diagonals[ndiagonals++] = local_goal + diagterm; */
	  if (ndiagonals == 0) {
	    last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	  } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	    last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	  }
	}
	debug(printf("\n"));
      }

      last_positions0 = *positions0;
      ++positions0;
      --npositions0;
    }
  }
  debug(printf("\n"));

  return ndiagonals;
}


/* diagonals is already allocated by caller */
/* Needs to return diagonals in ascending order.  Need to check if diagonals0 has duplicates. */
int
Intersect_approx_higher (Univcoord_T *diagonals,
			 unsigned char *positions1_high, UINT4 *positions1_low,
			 int npositions1, int diagterm1,
			 Univcoord_T *positions0, int npositions0,
			 Chrpos_T maxdistance) {
  int ndiagonals = 0;
  Univcoord_T local_goal, last_positions0;
  Univcoord_T diagterm, delta;
  Univcoord_T last_diagonal, this_diagonal;
  unsigned char *ptr1_high, *start1;
  UINT4 *ptr1_low;
  int j;
#ifdef DEBUG
  Univcoord_T *ptr0;
#endif

  diagterm = (Univcoord_T) diagterm1;	/* local_goal based on larger list */
  delta = (Univcoord_T) (-diagterm1); /* list0 + (diagterm0 - diagterm1) = list1 */

  start1 = positions1_high;
  /* end1 = &(positions1_high[npositions1]); */

  debug(printf("Intersect_approx_higher with %d positions against %d positions.  diagterm %d\n",
	       npositions1,npositions0,(int) diagterm));
#ifdef DEBUG
  for (ptr0 = positions0;  ptr0 < &(positions0[npositions0]); ptr0++) {
    printf("%llu\n",*ptr0);
  }
  printf("\n");

  for (ptr1_high = positions1_high, ptr1_low = positions1_low;
       ptr1_high < &(positions1_high[npositions1]); ptr1_high++, ptr1_low++) {
    printf("%llu + %d\n",GETPOS(*ptr1_high,*ptr1_low),diagterm1);
  }
  printf("\n");
#endif


#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagterm1 < 0) {
    while (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < (Univcoord_T) -diagterm1) {
      ++positions1_high;
      ++positions1_low;
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
      debug(printf("intersection list 0: %d:%llu => local_goal %llu.  Searching from npositions1 %d\n",
		   npositions0,*positions0,local_goal,npositions1));
      if (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < local_goal) {
	j = 1;
	while (j < npositions1 && GETPOS(positions1_high[j],positions1_low[j]) < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= npositions1) {
	  j = binary_search_large(j >> 1,npositions1,positions1_high,positions1_low,local_goal);
	} else {
	  j = binary_search_large(j >> 1,j,positions1_high,positions1_low,local_goal);
	}
	positions1_high += j;
	positions1_low += j;
	npositions1 -= j;
      }
#ifdef DEBUG
      if (npositions1 > 0) {
	printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,GETPOS(*positions1_high,*positions1_low));
      }
#endif
      
      if (npositions1 <= 0) {
	/* Check backwards only */
	debug(printf("    intersection list 1 at end  checking for approximate:"));
	ptr1_high = &(positions1_high[-1]);
	ptr1_low = &(positions1_low[-1]);
	if (ptr1_high >= start1) {
	  debug(printf(" prev %d:%llu?",npositions1-1,GETPOS(*ptr1_high,*ptr1_low)));
	  if (GETPOS(*ptr1_high,*ptr1_low) /*(omit for higher) + maxdistance*/ >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	    debug(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	    } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	debug(printf("\n"));
	
	/* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
	return ndiagonals;
	
      } else if (GETPOS(*positions1_high,*positions1_low) == local_goal) {
	/* Found local goal.  Save and advance */
	debug(printf("    intersection list 1: %d:%llu  exact\n",
		     npositions1,GETPOS(*positions1_high,*positions1_low)));
	/* diagonals[ndiagonals++] = local_goal + diagterm; */
	if (ndiagonals == 0) {
	  last_diagonal = diagonals[ndiagonals++] = local_goal + diagterm;
	} else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	  last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	}
	++positions1_high;
	++positions1_low;
	--npositions1;
	
      } else {
	debug(printf("    intersection list 1: %d:%llu  checking for approximate:",
		     npositions1,GETPOS(*positions1_high,*positions1_low)));
	ptr1_high = &(positions1_high[-1]); /* closest position < local_goal */
	ptr1_low = &(positions1_low[-1]); /* closest position < local_goal */
	if (ptr1_high >= start1) {
	  debug(printf(" prev %d:%llu?",npositions1+1,GETPOS(*ptr1_high,*ptr1_low)));
	  if (GETPOS(*ptr1_high,*ptr1_low) /*(omit for higher) + maxdistance*/ >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	    debug(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	    } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	
	ptr1_high = &(positions1_high[0]); /* closest position > local_goal */
	ptr1_low = &(positions1_low[0]); /* closest position > local_goal */
	debug(printf(" at %d:%llu?",npositions1,GETPOS(*ptr1_high,*ptr1_low)));
	if (GETPOS(*ptr1_high,*ptr1_low) /*(omit for higher) + maxdistance*/ >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	  debug(printf(" yes."));
	  /* diagonals[ndiagonals++] = local_goal + diagterm; */
	  if (ndiagonals == 0) {
	    last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	  } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	    last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	  }
	}
	debug(printf("\n"));
      }
      
      last_positions0 = *positions0;
      ++positions0;
      --npositions0;
    }
  }
  debug(printf("\n"));

  return ndiagonals;
}

