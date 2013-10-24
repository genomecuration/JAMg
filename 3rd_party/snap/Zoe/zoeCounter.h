/******************************************************************************\
zoeCounter.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_COUNTER_H
#define ZOE_COUNTER_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zoeDNA.h"
#include "zoeFeature.h"
#include "zoeMath.h"
#include "zoeModel.h"
#include "zoeTools.h"

struct zoeCounter  {
	coor_t               min_pos;                   /* minimum scoring position */
	coor_t               max_pos;                   /* maximum scoring position */
	zoeDNA               dna;                       /* counters require a dna */
	zoeDNA               anti;                      /* reverse-complement */
	zoeModel             model;                     /* counters require a model */
	struct zoeCounter ** subcounter;                /* for higher order models */
	char               * sig;                       /* binary signature (SDT) */
	void              (* count) (struct zoeCounter *, coor_t);  /* ptr to counting function */
	void              (* countf)(struct zoeCounter *, zoeFeature);
};
typedef struct zoeCounter * zoeCounter;

void       zoeDeleteCounter (zoeCounter);
zoeCounter zoeNewCounter (zoeDNA, zoeDNA, zoeModel);


#endif

