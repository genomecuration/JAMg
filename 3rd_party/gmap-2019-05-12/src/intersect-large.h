/* $Id: intersect-large.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef INTERSECT_LARGE_INCLUDED
#define INTERSECT_LARGE_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "bool.h"
#include "types.h"
#include "univcoord.h"
#include "genomicpos.h"


extern Univcoord_T *
Intersect_exact_large (int *ndiagonals,
		       unsigned char *positionsa_high, UINT4 *positionsa_low,
		       int npositionsa, int diagterma,
		       unsigned char *positionsb_high, UINT4 *positionsb_low,
		       int npositionsb, int diagtermb);

extern int
Intersect_exact_indices_large (int *indices,
			       unsigned char *positions1_high, UINT4 *positions1_low,
			       int npositions1, int diagterm1,
			       Univcoord_T *positions0, int npositions0);

extern Univcoord_T *
Intersect_approx_large (bool *exactp, int *ndiagpairs,
			unsigned char *positionsa_high, UINT4 *positionsa_low,
			int npositionsa, int diagterma,
			unsigned char *positionsb_high, UINT4 *positionsb_low,
			int npositionsb, int diagtermb,
			Chrpos_T maxdistance);

extern int
Intersect_approx_lower (Univcoord_T *diagonals,
			unsigned char *positions1_high, UINT4 *positions1_low,
			int npositions1, int diagterm1,
			Univcoord_T *positions0, int npositions0,
			Chrpos_T maxdistance);

extern int
Intersect_approx_higher (Univcoord_T *diagonals,
			 unsigned char *positions1_high, UINT4 *positions1_low,
			 int npositions1, int diagterm1,
			 Univcoord_T *positions0, int npositions0,
			 Chrpos_T maxdistance);

#endif

