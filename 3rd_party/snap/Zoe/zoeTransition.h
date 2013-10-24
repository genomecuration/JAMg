/******************************************************************************\
zoeTransition.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_TRANSITION_H
#define ZOE_TRANSITION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeMath.h"
#include "zoeFeature.h"
#include "zoeTools.h"

struct zoeTransition  {
	zoeLabel from;  /* state id */
	zoeLabel to;    /* state id */
	float    prob;  /* probability of transition */
	score_t  score; /* score of transition */
};
typedef struct zoeTransition * zoeTransition;

void          zoeDeleteTransition (zoeTransition);
zoeTransition zoeNewTransition (const char *, const char *, float);
zoeTransition zoeReadTransition (FILE *);
void          zoeWriteTransition (FILE *, const zoeTransition);

#endif
