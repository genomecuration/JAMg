/* $Id: ladder.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef LADDER_INCLUDED
#define LADDER_INCLUDED

#include "bool.h"
#include "list.h"
#include "stage3hr.h"
#include "types.h"
#include "hitlistpool.h"

#define T Ladder_T
typedef struct T *T;

extern int
Ladder_maxscore (T this);
extern int
Ladder_cutoff (T this);
extern Stage3end_T *
Ladder_hits_for_score (int *nhits, T this, Hitlistpool_T hitlistpool, int score);
extern int *
Ladder_nhits (T this);
extern Univcoord_T *
Ladder_genomicstarts (int *ndiagonals, T this);
extern Univcoord_T *
Ladder_genomicends (int *ndiagonals, T this);

extern void
Ladder_gc_hits (T this);
extern int
Ladder_minimax_trim (T ladder_plus, T ladder_minus, int querylength);

extern void
Ladder_free (T *old);
extern void
Ladder_gc_duplicates (T this);
extern void
Ladder_to_hits (List_T *hits5, List_T *hits3,
		T *ladder5_plus, T *ladder5_minus,
		T *ladder3_plus, T *ladder3_minus,
		Hitlistpool_T hitlistpool);
extern T
Ladder_new (List_T hitlist, Hitlistpool_T hitlistpool, bool end5p);

extern void
Ladder_merge (T dest, T source, Hitlistpool_T hitlistpool);

#undef T
#endif
