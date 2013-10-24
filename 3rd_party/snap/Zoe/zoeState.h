/******************************************************************************\
zoeState.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_STATE_H
#define ZOE_STATE_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zoeFeature.h"
#include "zoeTools.h"

typedef enum {
	INTERNAL,
	EXTERNAL,
	SHUTTLE
} zoeStateType;

struct zoeState  {
	zoeStateType type;      /* enumerated above */
	zoeLabel     label;     /* typical state label */
	float        init;      /* initial probability of state */
	float        term;      /* terminal probability of state */
	int          min;       /* minimum length */
	int          max;       /* maximum length, -1 is unlimited */
	int          geometric; /* true/false for geometric length */
};
typedef struct zoeState * zoeState;

void     zoeDeleteState (zoeState);
zoeState zoeNewState (zoeStateType, zoeLabel, float, float, int, int, int);
zoeState zoeReadState (FILE *);
void     zoeWriteState (FILE *, const zoeState);

#endif
