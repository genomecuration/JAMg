/* $Id: intersect.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef INTERSECT_INCLUDED
#define INTERSECT_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "bool.h"
#include "types.h"
#include "univcoord.h"
#include "genomicpos.h"


extern UINT4 *
Intersect_exact (int *ndiagonals,
		 UINT4 *positionsa, int npositionsa, int diagterma,
		 UINT4 *positionsb, int npositionsb, int diagtermb);

extern int
Intersect_exact_indices_univcoord (int *indices,
				   Univcoord_T *positions1, int npositions1,
				   Univcoord_T *positions0, int npositions0);

#ifndef LARGE_GENOMES
extern int
Intersect_exact_indices_small (int *indices,
			       UINT4 *positions1, int npositions1, int diagterm1,
			       Univcoord_T *positions0, int npositions0);
#endif

extern UINT4 *
Intersect_approx (bool *exactp, int *ndiagpairs,
		  UINT4 *positionsa, int npositionsa, int diagterma,
		  UINT4 *positionsb, int npositionsb, int diagtermb,
		  Chrpos_T maxdistance);

extern Univcoord_T *
Intersect_approx_simple (bool *exactp, int *ndiagpairs,
			 Univcoord_T *positionsa, int npositionsa,
			 Univcoord_T *positionsb, int npositionsb,
			 Chrpos_T maxdistance);

#ifndef LARGE_GENOMES
extern int
Intersect_approx_lower (Univcoord_T *diagonals,
			UINT4 *positions1, int npositions1, int diagterm1,
			Univcoord_T *positions0, int npositions0,
			Chrpos_T maxdistance);
#endif

#ifndef LARGE_GENOMES
extern int
Intersect_approx_higher (Univcoord_T *diagonals,
			 UINT4 *positions1, int npositions1, int diagterm1,
			 Univcoord_T *positions0, int npositions0,
			 Chrpos_T maxdistance);
#endif

#endif

