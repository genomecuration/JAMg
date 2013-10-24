/******************************************************************************\
zoeHMM.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_HMM_H
#define ZOE_HMM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeDuration.h"
#include "zoeState.h"
#include "zoeModel.h"
#include "zoeFeature.h"
#include "zoePhasePref.h"
#include "zoeTransition.h"
#include "zoeTools.h"

struct zoeHMM  {
	/* general attributes */
	char * name;        /* completely arbitrary */
	int    states;      /* number of states in the HMM */
	int    transitions; /* number of transtions in the HMM */
	int    durations;   /* number of duration models */
	int    models;      /* number of sequence models */
	
	/* object storage */
	zoeState      * state;      /* array of zoeState */
	zoeTransition * transition; /* array of zoeTransition */
	zoeDuration   * duration;   /* array of zoeDuration */
	zoeModel      * model;      /* array of zoeModel */
	zoePhasePref    phasepref;  /* exon-intron phase preferences */
	
	/* object mapping */
	zoeDuration dmap[zoeLABELS];            /* durations */
	zoeModel    mmap[zoeLABELS];            /* models */
	zoeState    smap[zoeLABELS];            /* states */
	zoeIVec     jmap[zoeLABELS][zoeLABELS]; /* jump list (reverse arrows) */
	score_t     imap[zoeLABELS];            /* initial probability */
	score_t     kmap[zoeLABELS];            /* terminal probability */
	score_t     xmap[zoeLABELS];            /* geometric extension score */
	score_t     tmap[zoeLABELS][zoeLABELS]; /* transition score */
	coor_t      cmap[zoeLABELS];            /* coordinate adjustments */
};
typedef struct zoeHMM * zoeHMM;

void   zoeSetNscore (zoeLabel, score_t);
void   zoeSetAscore (zoeLabel, score_t);
void   zoeDeleteHMM (zoeHMM);
zoeHMM zoeNewHMM (void);
zoeHMM zoeReadHMM (FILE *);
void   zoeWriteHMM (FILE *, const zoeHMM);
zoeHMM zoeGetHMM (const char *);
score_t zoeGetAscore (int);

#endif
