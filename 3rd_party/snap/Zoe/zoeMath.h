/******************************************************************************\
zoeMath.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_MATH_H
#define ZOE_MATH_H

#include <math.h>
#include "zoeTools.h"

score_t zoeFloat2Score (double);
double  zoeScore2Float (score_t);
double  zoeLog2 (double);
double  zoeLnFactorial (int);
double  zoeDivide (double, double);
score_t zoeScoreGeometric (double, double);
score_t zoeScorePoisson (double, double);
void    zoeDecToBase (int, int, char *);
int     zoeBaseToDec (int, const char *);

extern const int zoePOWER[6][8];

#endif
