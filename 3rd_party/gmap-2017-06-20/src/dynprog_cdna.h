/* $Id: dynprog_cdna.h 184440 2016-02-17 21:03:49Z twu $ */
#ifndef DYNPROG_CDNA_INCLUDED
#define DYNPROG_CDNA_INCLUDED

#include "bool.h"
#include "list.h"
#include "pairpool.h"
#include "chrnum.h"
#include "iit-read.h"
#include "types.h"
#include "dynprog.h"

#define T Dynprog_T

/* Sequences rsequenceL and rsequenceR represent the two ends of the cDNA insertion */
extern List_T
Dynprog_cdna_gap (int *dynprogindex, int *finalscore, bool *incompletep,
		  T dynprogL, T dynprogR, char *rsequenceL, char *rsequence_ucL, 
		  char *rev_rsequenceR, char *rev_rsequence_ucR,
#if 0
		  char *gsequence, char *gsequence_uc,
#endif
		  int rlengthL, int rlengthR, int glength,
		  int roffsetL, int rev_roffsetR, int goffset,
		  Univcoord_T chroffset, Univcoord_T chrhigh,
		  bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_paired, double defect_rate);


#undef T
#endif

