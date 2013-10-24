/******************************************************************************\
zoeScanner.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_SCANNER_H
#define ZOE_SCANNER_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zoeDNA.h"
#include "zoeFeature.h"
#include "zoeMath.h"
#include "zoeModel.h"
#include "zoeTools.h"

struct zoeScanner  {
	coor_t               min_pos;    /* minimum scoring position */
	coor_t               max_pos;    /* maximum scoring position */
	zoeDNA               dna;        /* scanners require a seq */
	zoeDNA               anti;       /* reverse-complement */
	zoeModel             model;      /* scanners require a model */
	struct zoeScanner ** subscanner; /* for higher order models */
	char               * sig;        /* binary signature (SDT) */
	score_t            * uscore;     /* user-defined score */
	score_t            * ascore;     /* user-defined anti-parallel score */
	score_t           (* score) (struct zoeScanner *, coor_t);
	score_t           (* scoref)(struct zoeScanner *, zoeFeature);
};
typedef struct zoeScanner * zoeScanner;

void       zoeDeleteScanner (zoeScanner);
zoeScanner zoeNewScanner (zoeDNA, zoeDNA, zoeModel);
void       zoeSetScannerScore (zoeScanner, coor_t, score_t);

#endif
