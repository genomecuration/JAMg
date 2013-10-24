/******************************************************************************\
zoeIsochore.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_ISOCHORE_H
#define ZOE_ISOCHORE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeHMM.h"
#include "zoeTools.h"

struct zoeIsochore  {
	int     count;
	zoeFVec min_GC;
	zoeFVec max_GC;
	zoeTVec hmm_file;
	zoeVec  hmms;
};
typedef struct zoeIsochore * zoeIsochore;

void        zoeDeleteIsochore (zoeIsochore);
zoeIsochore zoeNewIsochore (void);
zoeIsochore zoeReadIsochore (FILE *);
zoeIsochore zoeGetIsochore (const char *);
zoeHMM      zoeSelectIsochore (const zoeIsochore, float);

#endif
