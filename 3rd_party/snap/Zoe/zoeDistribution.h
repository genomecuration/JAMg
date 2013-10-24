/******************************************************************************\
zoeDistribution.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_DISTRIBUTION_H
#define ZOE_DISTRIBUTION_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zoeMath.h"
#include "zoeTools.h"

typedef enum {
	UNKNOWN,
	DEFINED,
	GEOMETRIC,
	POISSON,
	CONSTANT
} zoeDistributionType;

struct zoeDistribution  {
	zoeDistributionType   type;   /* the type of function (see above) */
	coor_t                start;  /* starting coordinate */
	coor_t                end;    /* ending coordinate */
	int                   params; /* number of parameters to function */
	float               * param;  /* parameters of function */
};
typedef struct zoeDistribution * zoeDistribution;

void            zoeDeleteDistribution (zoeDistribution);
zoeDistribution zoeNewDistribution (zoeDistributionType, coor_t, coor_t, int, const float *);
zoeDistribution zoeReadDistribution (FILE *);
void            zoeWriteDistribution (FILE *, const zoeDistribution);
score_t         zoeScoreDistribution (const zoeDistribution, coor_t);

#endif
