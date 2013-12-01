/******************************************************************************\
zoeAlignment.h - part of the ZOE library for genomic analysis
 
Copyright (C) Ian Korf 2002-2013.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The MIT License (MIT) - opensource.org/licenses/MIT

\******************************************************************************/

#ifndef ZOE_ALIGNMENT_H
#define ZOE_ALIGNMENT_H

#include <stdio.h>
#include <string.h>

#include "zoeDNA.h"

typedef enum {
	DIAGONAL,
	HORIZONTAL,
	VERTICAL
} zoeSWdirection;

struct zoeHSP  {
	score_t  score;
	score_t  bits;
	score_t  e_value;
	int      match;
	int      mismatch;
	int      q_gap;
	int      s_gap;
	coor_t   q_start;
	coor_t   q_end;
	coor_t   s_start;
	coor_t   s_end;
	char   * q_aln;
	char   * s_aln;
	char   * a_str;
};
typedef struct zoeHSP * zoeHSP;

struct zoeAlignCell01  {
	int  score;
	char trace; /* 0=blank   1=diagonal 2=horizontal 3=vertical */
	char state; /* 0=gapless 1=gap1     2=gap2                  */
};
typedef struct zoeAlignCell01 zoeAlignCell01;

struct zoeAlignCell02  {
	int  score;
	char trace; /* 0=blank   1=diagonal 2=horizontal 3=vertical */
	char state; /* 0=gapless 1=gap1     2=gap2                  */
	char found; /* 0=no      1=yes                              */
};
typedef struct zoeAlignCell02 zoeAlignCell02;

void   zoeDeleteHSP (zoeHSP);
zoeHSP zoeNewHSP (void);
zoeHSP zoeReadHSP (FILE *);
void   zoeWriteHSP (FILE *, const zoeHSP);
void   zoeDecorateHSP (zoeHSP);
zoeHSP zoeAlign01 (const zoeDNA, const zoeDNA, int, int, int, int);
zoeVec zoeAlign02 (const zoeDNA, const zoeDNA, int, int, int, int, int);
zoeHSP zoeAlign03 (const zoeDNA, const zoeDNA, int, int, int, int, int);
/*zoeVec zoeChainHSPs (const zoeVec);*/

struct zoePCell {
	float score;
	char  trace;
	char  state;
};
typedef struct zoePCell zoePCell;

struct zoeScoreSystem {
	float ** score;
	int      gapO;
	int      gapE;
	int      width; /* banded alignment, 0 is infinite */
};
typedef struct zoeScoreSystem * zoeScoreSystem;

zoeScoreSystem zoeNewScoreSystem(char *, int, int, int);
void   zoeDeleteScoreSystem(zoeScoreSystem);
zoeHSP zoeProtAlignLocal(zoeProtein, zoeProtein, zoeScoreSystem);
zoeHSP zoeProtAlignNterm(zoeProtein, zoeProtein, zoeScoreSystem);

#endif

