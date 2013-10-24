/******************************************************************************\
zoeDuration.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_DURATION_H
#define ZOE_DURATION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeDistribution.h"
#include "zoeFeature.h"
#include "zoeTools.h"

struct zoeDuration  {
	zoeLabel          label;
	int               distributions;
	zoeDistribution * distribution;
};
typedef struct zoeDuration * zoeDuration;

void        zoeDeleteDuration (zoeDuration);
zoeDuration zoeNewDuration (zoeLabel, int, const zoeDistribution *);
zoeDuration zoeReadDuration (FILE *);
void        zoeWriteDuration (FILE *, const zoeDuration);
score_t     zoeScoreDuration (const zoeDuration, coor_t);

#endif
