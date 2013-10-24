/******************************************************************************\
zoeFeatureFactory.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_FEATUREFACTORY_H
#define ZOE_FEATUREFACTORY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeDNA.h"
#include "zoeScanner.h"
#include "zoeFeature.h"
#include "zoeFeatureTable.h"
#include "zoeTools.h"

struct zoeFeatureFactory  {
	/* used by all factories */
	zoeFeatureVec (* create)(struct zoeFeatureFactory *, coor_t);
	zoeLabel         type;
	zoeDNA           dna;
	
	/* RFactory and OFactory/XFactory */
	int        length;  /* min length */
	zoeScanner scanner; /* for those factories that use a scanner */
	
	/* OFactory/XFactory only */
	zoeFeatureVec orfs;
	zoeHash       hash;
	score_t       score;
	strand_t      strand;
	
	/* EFactory */
	int       offset;
	score_t * cds[3]; /* cds score in each FRAME */
	score_t * start;  /* score at each valid position */
	score_t * stop;   /* MIN_SCORE at invalid positions */
	score_t * acc;
	score_t * don;
	int * fstop;  /* position of previous stop in this frame */
	/* changed to int 2004-04-06 */
	
};
typedef struct zoeFeatureFactory * zoeFeatureFactory;

void              zoeDeleteFeatureFactory (zoeFeatureFactory);
zoeFeatureFactory zoeNewEFactory (zoeScanner, zoeScanner, zoeScanner, zoeScanner, zoeScanner);
zoeFeatureFactory zoeNewOFactory (zoeScanner, coor_t, score_t, strand_t);
zoeFeatureFactory zoeNewXFactory (zoeScanner, zoeScanner, zoeScanner, zoeScanner, zoeScanner, coor_t, score_t);
zoeFeatureFactory zoeNewRFactory (zoeScanner, coor_t);
zoeFeatureFactory zoeNewSFactory (zoeScanner, zoeLabel);

#endif
