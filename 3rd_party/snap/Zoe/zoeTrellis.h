/******************************************************************************\
zoeTrellis.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_TRELLIS_H
#define ZOE_TRELLIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeCDS.h"
#include "zoeDNA.h"
#include "zoeFeatureFactory.h"
#include "zoeHMM.h"
#include "zoeModel.h"
#include "zoeScanner.h"
#include "zoeFeature.h"
#include "zoeFeatureTable.h"
#include "zoeTools.h"

struct zoeTrellis {
	zoeDNA              dna;
	zoeDNA              anti;
	zoeHMM              hmm;
	zoeFeatureVec       xdef;
	score_t             max_score;           /* set at the end */
	score_t             exp_score;           /* expected score of null model */
	zoeFeatureVec       keep;                /* max features */
	zoeIVec             jump;                /* trace-back */
	int                 min_len[zoeLABELS];  /* minimum length (internal & external) */
	int                 max_len[zoeLABELS];  /* maximum explicit length (internal only) */
	zoeScanner          scanner[zoeLABELS];  /* map hmm models to scanners here */
	zoeFeatureFactory   factory[zoeLABELS];  /* external feature factory */
	int                 internal[zoeLABELS]; /* internal states used */	
	int               * trace[zoeLABELS];    /* viterbi trace-back */
	score_t           * score[zoeLABELS];    /* viterbi score */
	zoeFeatureVec       features[zoeLABELS]; /* current features[state_label] */
	score_t          (* ext)(struct zoeTrellis *, coor_t, zoeLabel, zoeFeature);
};
typedef struct zoeTrellis * zoeTrellis;

void       zoeDeleteTrellis (zoeTrellis);
zoeTrellis zoeNewTrellis (zoeDNA, zoeHMM, zoeFeatureVec);
zoeVec     zoePredictGenes (zoeTrellis);
void       zoeScoreCDS(zoeTrellis, zoeCDS, int, int);
void       zoeSetTrellisMeter (int);
void       zoeSetTrellisPadding (int);
char*      zoeGetPartialProtein (zoeTrellis, zoeLabel, zoeFeature);
score_t    zoeScoreExon   (zoeTrellis, zoeFeature, int, int);
score_t    zoeScoreIntron (zoeTrellis, zoeFeature, int);

#endif
