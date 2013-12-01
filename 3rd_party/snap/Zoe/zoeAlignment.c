/******************************************************************************\
zoeAlignment.c - part of the ZOE library for genomic analysis

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

#ifndef ZOE_ALIGNMENT_C
#define ZOE_ALIGNMENT_C

#include <stdio.h>
#include <string.h>

#include "zoeAlignment.h"

/* private structures */


struct zoeHotSpot {
	int score;
	int i;
	int j;
};
typedef struct zoeHotSpot * zoeHotSpot;


struct zoeHSPcell  {
	int    score;
	zoeHSP hsp;
	char   trace; /* 0=blank   1=diagonal 2=horizontal 3=vertical */
	int    i;
	int    j;
};
typedef struct zoeHSPcell * zoeHSPcell;



void zoeDeleteHSP (zoeHSP hsp) {	
	if (hsp == NULL) return;
	
	if (hsp->q_aln) {
		zoeFree(hsp->q_aln);
		hsp->q_aln = NULL;
	}
	if (hsp->s_aln) {
		zoeFree(hsp->s_aln);
		hsp->s_aln = NULL;
	}
	if (hsp->a_str) {
		zoeFree(hsp->a_str);
		hsp->a_str = NULL;
	}
	zoeFree(hsp);
	hsp = NULL;
}

zoeHSP zoeNewHSP (void) {
	zoeHSP hsp = zoeMalloc(sizeof(struct zoeHSP));
	
	hsp->score   = -1;
	hsp->bits    = -1;
	hsp->e_value = -1;
	
	hsp->match    = -1;
	hsp->mismatch = -1;
	hsp->q_gap    = -1;
	hsp->s_gap    = -1;
	
	hsp->q_start  = -1;
	hsp->q_end    = -1;
	hsp->s_start  = -1;
	hsp->s_end    = -1;
	
	hsp->q_aln = NULL;
	hsp->s_aln = NULL;
	hsp->a_str = NULL;
	return hsp;
}

/*

>query1

sbjct line
HSP query
HSP sbjct

sbjct line
HSP query
HSP sbjct

>query2

sbjct line
HSP query
HSP sbjct



There are 2 forms for the sbjct line
------------------------------------

<sbjct>: <score> <match> <mismatch> <qgap> <sgap>     // long form
<sbjct>: <score> <%ID>                                // short form

There are 4 forms for the HSP lines
-----------------------------------

HSP <start> <end> <+/-> . . . <frame> <aln>           // long form
HSP <start> <end>   +   . . . <frame> <aln>           // long form

HSP <start> <end> <+/-> . . . <frame>                 // long form, no alignment
HSP <start> <end>   +   . . . <frame>                 // long form, no alignment

HSP <start> <end> <aln>                               // short form
HSP <start> <end> <aln>                               // short form

HSP <start> <end>                                     // short form, no alignment
HSP <start> <end>                                     // short form, no alignment

*/

zoeHSP zoeReadHSP (FILE * stream) {
	zoeHSP hsp = zoeMalloc(sizeof(struct zoeHSP));
	return hsp;
}

void zoeWriteHSP (FILE * stream, const zoeHSP hsp) {
	
	/* stats line */
	zoeS(stream, "Stats:");
	if (hsp->score   >= 0) zoeS(stream, " score=%d", (int)hsp->score);
	if (hsp->bits    >= 0) zoeS(stream, " bits=%g", hsp->bits);
	if (hsp->e_value >= 0) zoeS(stream, " E=%g", hsp->e_value);
	zoeS(stream, "\n");
	
	/* alignment info */
	zoeS(stream, "Alignment:");
	zoeS(stream, " Q:%d..%d", hsp->q_start +1, hsp->q_end +1);
	zoeS(stream, " S:%d..%d", hsp->s_start +1, hsp->s_end +1);
	if (hsp->match >= 0 && hsp->mismatch >= 0) {
		zoeS(stream, " %d/%d", hsp->match, hsp->mismatch);
		if (hsp->q_gap >= 0 && hsp->s_gap >= 0) zoeS(stream, " %d,%d", hsp->q_gap, hsp->s_gap);
	}
	zoeS(stream, "\n");
	
	/* alignments */
	if (hsp->q_aln) zoeS(stream, "Q: %s\n", hsp->q_aln);
	if (hsp->a_str) zoeS(stream, "   %s\n", hsp->a_str);
	if (hsp->s_aln) zoeS(stream, "S: %s\n", hsp->s_aln);
	zoeS(stream, "\n");
	
}

void zoeDecorateHSP (zoeHSP hsp) {
	int i, length;
	int qgap = 0, sgap = 0, match = 0, mismatch = 0;
		
	length = strlen(hsp->q_aln);
	for (i = 0; i < length; i++) {
		if      (hsp->q_aln[i] == hsp->s_aln[i]) match++;
		else if (hsp->q_aln[i] == '-')           qgap++;
		else if (hsp->s_aln[i] == '-')           sgap++;
		else                                     mismatch++;
	}
	
	hsp->match    = match;
	hsp->mismatch = mismatch;
	hsp->q_gap    = qgap;
	hsp->s_gap    = sgap;

}

int zoeHotSpotCmp (const zoeHotSpot h1, const zoeHotSpot h2) {
	return h2->score - h1->score; /* sort in descending order */
}

int zoeHotSpotCmpPtr (const void * h1, const void * h2) {
	return zoeHotSpotCmp( *(zoeHotSpot *)h1, *(zoeHotSpot *)h2 );
}

void zoeViewMatrix (zoeAlignCell02 ** matrix, zoeDNA dna1, zoeDNA dna2) {
	int i, j;
	
	for (i = 0; i < dna1->length; i++) zoeO("\t%c", dna1->seq[i]);
	zoeO("\n");
	for (j = 1; j <= dna2->length; j++) {
		zoeO("%c", dna2->seq[j-1]);
		for (i = 1; i <= dna1->length; i++) {
			zoeO("\t");
			switch (matrix[j][i].trace) {
				case 1: zoeO("+"); break;
				case 2: zoeO("<"); break;
				case 3: zoeO("^"); break;
			}
			
			zoeO("%d", matrix[j][i].score);
		}
		zoeO("\n");
	}
}

void flip_string (char * string) {
	int i, length, half;
	char c;
	
	length = strlen(string);
	half   = length / 2;
		
	for (i = 0; i < half; i++) {
		c = string[i]; /* save it */
		string[i] = string[length -i -1];
		string[length -i -1] = c;
	}
}

zoeHSP zoeTraceBack01 (zoeAlignCell01 ** matrix, int i, int j, zoeDNA dna1, zoeDNA dna2) {
	zoeHSP   hsp = NULL;
	char   * a1;         /* alignment sequence 1 */
	char   * a2;         /* alignment sequence 2 */
	char   * as;         /* alignment string */
	int      index = 0;  /* current index of alignment strings */
	int      abort = 0;  /* set to true if not a novel trace back */
	int      I = i, J = j, S = matrix[j][i].score;
		
	/* allocate maximal storage for alignment strings */
	a1 = zoeMalloc(i + j + 1);
	a2 = zoeMalloc(i + j + 1);
	as = zoeMalloc(i + j + 1);
		
	/* trace back */
	while (matrix[j][i].trace != 0) {
		
		switch (matrix[j][i].trace) {
			case 1: /* match state */
				i--;
				j--;
				a1[index] = dna1->seq[i];
				a2[index] = dna2->seq[j];
				as[index] = (a1[index] == a2[index]) ? '|' : ' ';
				break;
			case 2: /* horizontal gap */
				i--;
				a1[index] = dna1->seq[i];
				a2[index] = '-';
				as[index] = ' ';
				break;
			case 3: /* vertical gap */
				a1[index] = '-';
				a2[index] = dna2->seq[j];
				as[index] = ' ';
				j--;
				break;
			default:
				zoeExit("impossible");
		}
		index++;
		
	}
	
	if (abort) {
		zoeFree(a1);
		zoeFree(a2);
		zoeFree(as);
		return NULL;
	}
	
	a1[index] = '\0';
	a2[index] = '\0';
	as[index] = '\0';
	flip_string(a1);
	flip_string(a2);
	flip_string(as);
	
	hsp = zoeNewHSP();
	hsp->score   = S;
	hsp->q_start = i;
	hsp->q_end   = I -1;
	hsp->s_start = j;
	hsp->s_end   = J -1;
	hsp->q_aln   = a1;
	hsp->s_aln   = a2;
	hsp->a_str   = as;
	zoeDecorateHSP(hsp);
	
	return hsp;
}


zoeHSP zoeAlign01 (
	const zoeDNA dna1, /* query dna */
	const zoeDNA dna2, /* sbjct dna */
	int M, /* match */
	int N, /* mismatch */
	int G, /* gap */
	int E) /* gap extend */
{

	/* variables */
	int               i, j, dscore, hscore, vscore, max_i, max_j, max_s;
	zoeAlignCell01 ** matrix;
	zoeHSP            msp;
	
	/* allocate matrix */
	matrix = zoeMalloc((dna2->length +1) * sizeof(struct zoeAlignCell01));
	for (i = 0; i <= dna2->length; i++) {
		matrix[i] = zoeMalloc((dna1->length +1) * sizeof(struct zoeAlignCell01));
	}
	
	/* initialization */
	for (i = 0; i <= dna1->length; i++) {
		matrix[0][i].score = 0;
		matrix[0][i].trace = 0;
		matrix[0][i].state = 0;
	}
	for (j = 0; j <= dna2->length; j++) {
		matrix[j][0].score = 0;
		matrix[j][0].trace = 0;
		matrix[j][0].state = 0;
	}
	max_i = 0;
	max_j = 0;
	max_s = 0;
	
	/* induction */
	for (j = 1; j <= dna2->length; j++) {
		for (i = 1; i <= dna1->length; i++) {
			dscore  = (dna1->s5[i-1] == dna2->s5[j-1]) ? M : -N;
			dscore += matrix[j-1][i-1].score;
			vscore  = matrix[j-1][i].score - G;
			hscore  = matrix[j][i-1].score - G;
			
			if (dscore >= vscore && dscore >= hscore && dscore > 0) {
				matrix[j][i].score = dscore;
				matrix[j][i].trace = 1;
				matrix[j][i].state = 0;
				if (matrix[j][i].score > max_s) {
					max_s = matrix[j][i].score;
					max_j = j;
					max_i = i;
				}
			} else if (vscore >= hscore && vscore > 0) {
				matrix[j][i].score = vscore;
				matrix[j][i].trace = 3;
				matrix[j][i].state = 0;
			} else if (hscore > 0) {
				matrix[j][i].score = hscore;
				matrix[j][i].trace = 2;
				matrix[j][i].state = 0;
			} else {
				matrix[j][i].score = 0;
				matrix[j][i].trace = 0;
				matrix[j][i].state = 0;
			}
		}
	}
	
	/* trace-back */
	msp = zoeTraceBack01(matrix, max_i, max_j, dna1, dna2);
	
	/* clean up */
	for (i = 0; i <= dna2->length; i++) {
		zoeFree(matrix[i]);
	}
	zoeFree(matrix);
	
	return msp;
}

zoeHSP zoeTraceBack02 (zoeAlignCell02 ** matrix, int i, int j, zoeDNA dna1, zoeDNA dna2) {
	zoeHSP   hsp = NULL;
	char   * a1;         /* alignment sequence 1 */
	char   * a2;         /* alignment sequence 2 */
	char   * as;         /* alignment string */
	int      index = 0;  /* current index of alignment strings */
	int      abort = 0;  /* set to true if not a novel trace back */
	int      I = i, J = j, S = matrix[j][i].score;
		
	/* allocate maximal storage for alignment strings */
	a1 = zoeMalloc(i + j + 1);
	a2 = zoeMalloc(i + j + 1);
	as = zoeMalloc(i + j + 1);
	
	/* trace back */
	while (matrix[j][i].trace != 0) {
		if (matrix[j][i].found) {
			abort = 1;
			break;
		}
		
		matrix[j][i].found = 1;
		switch (matrix[j][i].trace) {
			case 1: /* match state */
				i--;
				j--;
				a1[index] = dna1->seq[i];
				a2[index] = dna2->seq[j];
				as[index] = (a1[index] == a2[index]) ? '|' : ' ';
				break;
			case 2: /* horizontal gap */
				i--;
				a1[index] = dna1->seq[i];
				a2[index] = '-';
				as[index] = ' ';
				break;
			case 3: /* vertical gap */
				a1[index] = '-';
				a2[index] = dna2->seq[j];
				as[index] = ' ';
				j--;
				break;
			default:
				zoeExit("impossible");
		}
		index++;
		
	}
	
	if (abort) {
		zoeFree(a1);
		zoeFree(a2);
		zoeFree(as);
		return NULL;
	}
	
	a1[index] = '\0';
	a2[index] = '\0';
	as[index] = '\0';
	flip_string(a1);
	flip_string(a2);
	flip_string(as);
	
	hsp = zoeNewHSP();
	hsp->score   = S;
	hsp->q_start = i;
	hsp->q_end   = I -1;
	hsp->s_start = j;
	hsp->s_end   = J -1;
	hsp->q_aln   = a1;
	hsp->s_aln   = a2;
	hsp->a_str   = as;
	zoeDecorateHSP(hsp);
	
	return hsp;
}

zoeVec zoeAlign02 (
	const zoeDNA dna1, /* query dna */
	const zoeDNA dna2, /* sbjct dna */
	int M, /* match */
	int N, /* mismatch */
	int G, /* gap */
	int E, /* gap extend */
	int S) /* score */
{

	/* variables */
	int               i, j, k, dscore, hscore, vscore;
	zoeAlignCell02 ** matrix;
	zoeVec            HSPs  = zoeNewVec();
	zoeVec            spots = zoeNewVec();
	zoeHotSpot        spot;
	zoeHSP            hsp;	
	
	/* allocate matrix */
	matrix = zoeMalloc((dna2->length +1) * sizeof(struct zoeAlignCell02));
	for (i = 0; i <= dna2->length; i++) {
		matrix[i] = zoeMalloc((dna1->length +1) * sizeof(struct zoeAlignCell02));
	}
	
	/* initialization */
	for (i = 0; i <= dna1->length; i++) {
		matrix[0][i].score = 0;
		matrix[0][i].trace = 0;
		matrix[0][i].state = 0;
		matrix[0][i].found = 0;
	}
	for (j = 0; j <= dna2->length; j++) {
		matrix[j][0].score = 0;
		matrix[j][0].trace = 0;
		matrix[j][0].state = 0;
		matrix[j][0].found = 0;
	}
	
	/* induction */
	for (j = 1; j <= dna2->length; j++) {
		for (i = 1; i <= dna1->length; i++) {
			dscore  = (dna1->s5[i-1] == dna2->s5[j-1]) ? M : -N;
			dscore += matrix[j-1][i-1].score;
			vscore  = matrix[j-1][i].score - G;
			hscore  = matrix[j][i-1].score - G;
			
			if (dscore >= vscore && dscore >= hscore && dscore > 0) {
				matrix[j][i].score = dscore;
				matrix[j][i].trace = 1;
				matrix[j][i].state = 0;
				matrix[j][i].found = 0;
			} else if (vscore >= hscore && vscore > 0) {
				matrix[j][i].score = vscore;
				matrix[j][i].trace = 3;
				matrix[j][i].state = 0;
				matrix[j][i].found = 0;
			} else if (hscore > 0) {
				matrix[j][i].score = hscore;
				matrix[j][i].trace = 2;
				matrix[j][i].state = 0;
				matrix[j][i].found = 0;
			} else {
				matrix[j][i].score = 0;
				matrix[j][i].trace = 0;
				matrix[j][i].state = 0;
				matrix[j][i].found = 0;
			}
		}
	}
	
	/* show matrix */
	if (zoeOption("-debug")) zoeViewMatrix(matrix, dna1, dna2);
	
	/* gather hot spots */
	for (j = 1; j <= dna2->length; j++) {
		for (i = 1; i <= dna1->length; i++) {
			if (matrix[j][i].score >= S) {
				spot = zoeMalloc(sizeof(struct zoeHotSpot));
				spot->i = i;
				spot->j = j;
				spot->score = matrix[j][i].score;
				zoePushVec(spots, spot);
			}
		}
	}
	qsort(spots->elem, spots->size, sizeof(zoeHotSpot), zoeHotSpotCmpPtr);
	
	/* multiple trace-back */
	for (k = 0; k < spots->size; k++) {
		spot = spots->elem[k];
		if (matrix[spot->j][spot->i].found) continue;
		
		hsp = zoeTraceBack02(matrix, spot->i, spot->j, dna1, dna2);
		if (hsp != NULL) {
		
			/* clear neighboring cells */
			for (j = hsp->s_start; j <= hsp->s_end; j++) {
				for (i = hsp->q_start; i <= hsp->q_end; i++) {
					matrix[j][i].found = 1;
				}
			}
			zoePushVec(HSPs, hsp);
		}
	}
	
	for (i = 0; i < spots->size; i++) zoeFree(spots->elem[i]);
	zoeDeleteVec(spots);
	
	/* clean up */
	for (i = 0; i <= dna2->length; i++) {
		zoeFree(matrix[i]);
	}
	zoeFree(matrix);

	return HSPs;
}


zoeHSP zoeAlign03 (
	const zoeDNA dna1, /* query dna */
	const zoeDNA dna2, /* sbjct dna */
	int M, /* match */
	int N, /* mismatch */
	int G, /* gap */
	int E, /* gap extend */
	int GC /* GC score */
	) 
{

	/* variables */
	int               i, j, dscore, hscore, vscore, max_i, max_j, max_s;
	zoeAlignCell01 ** matrix;
	zoeHSP            msp;
	
	/* allocate matrix */
	matrix = zoeMalloc((dna2->length +1) * sizeof(struct zoeAlignCell01));
	for (i = 0; i <= dna2->length; i++) {
		matrix[i] = zoeMalloc((dna1->length +1) * sizeof(struct zoeAlignCell01));
	}
	
	/* initialization */
	for (i = 0; i <= dna1->length; i++) {
		matrix[0][i].score = 0;
		matrix[0][i].trace = 0;
		matrix[0][i].state = 0;
	}
	for (j = 0; j <= dna2->length; j++) {
		matrix[j][0].score = 0;
		matrix[j][0].trace = 0;
		matrix[j][0].state = 0;
	}
	max_i = 0;
	max_j = 0;
	max_s = 0;
	
	/* induction */
	for (j = 1; j <= dna2->length; j++) {
		for (i = 1; i <= dna1->length; i++) {
			if (dna1->s5[i-1] == 1 && dna2->s5[j-1] == 2) {
				dscore = GC;
			} else if (dna1->s5[i-1] == 2 && dna2->s5[j-1] == 1) {
				dscore = GC;
			} else {
				dscore = (dna1->s5[i-1] == dna2->s5[j-1]) ? M : -N;
			}
			dscore += matrix[j-1][i-1].score;
			vscore  = matrix[j-1][i].score - G;
			hscore  = matrix[j][i-1].score - G;
			
			if (dscore >= vscore && dscore >= hscore && dscore > 0) {
				matrix[j][i].score = dscore;
				matrix[j][i].trace = 1;
				matrix[j][i].state = 0;
				if (matrix[j][i].score > max_s) {
					max_s = matrix[j][i].score;
					max_j = j;
					max_i = i;
				}
			} else if (vscore >= hscore && vscore > 0) {
				matrix[j][i].score = vscore;
				matrix[j][i].trace = 3;
				matrix[j][i].state = 0;
			} else if (hscore > 0) {
				matrix[j][i].score = hscore;
				matrix[j][i].trace = 2;
				matrix[j][i].state = 0;
			} else {
				matrix[j][i].score = 0;
				matrix[j][i].trace = 0;
				matrix[j][i].state = 0;
			}
		}
	}
	
	/* trace-back */
	msp = zoeTraceBack01(matrix, max_i, max_j, dna1, dna2);
	
	/* clean up */
	for (i = 0; i <= dna2->length; i++) {
		zoeFree(matrix[i]);
	}
	zoeFree(matrix);
	
	return msp;
}



int zoeHSPCmpQuery (const zoeHSP h1, const zoeHSP h2) {
	if      (h1->q_start < h2->q_start && h1->q_end < h2->q_end) return -1;
	else if (h1->q_start > h2->q_start && h2->q_end > h2->q_end) return  1;
	else    return  0;
}

int zoeHSPCmpSbjct (const zoeHSP h1, const zoeHSP h2) {
	if      (h1->s_start < h2->s_start && h1->s_end < h2->s_end) return -1;
	else if (h1->s_start > h2->s_start && h2->s_end > h2->s_end) return  1;
	else    return  0;
}

int zoeHSPCmpQ (const void * h1, const void * h2) {
	return zoeHSPCmpQuery( *(zoeHSP *)h1, *(zoeHSP *)h2 );
}

int zoeHSPCmpS (const void * h1, const void * h2) {
	return zoeHSPCmpSbjct( *(zoeHSP *)h1, *(zoeHSP *)h2 );
}

zoeVec zoeGroupHSPs(zoeVec HSPs, char type) {
	zoeVec group = zoeNewVec();
	int i;
	
	/* copy HSPs */
	for (i = 0; i < HSPs->size; i++) zoePushVec(group, HSPs->elem[i]);

	
	/* sort HSPs */
	switch (type) {
		case 'q':
			qsort(group->elem, group->size, sizeof(zoeHSP), zoeHSPCmpQ);
			break;
		case 's':
			qsort(group->elem, group->size, sizeof(zoeHSPcell), zoeHSPCmpS);
			break;
		default:
			zoeExit("unsupported");	
	}
	
	/* group HSPs */
	for (i = 0; i < group->size; i++) {
	
	}
	
	
	return group;
}

/*

zoeVec zoeChainHSPs (const zoeVec HSPs) {
	int           i;
	zoeVec        bestHSPs = NULL;
	zoeVec        qvec;
	zoeVec        svec;

	for (i = 0; i < HSPs->size; i++) {
		zoeWriteHSP(stdout, HSPs->elem[i]);
	}
	
	qvec = zoeGroupHSPs(HSPs, 'q');
	svec = zoeGroupHSPs(HSPs, 's');
	
	return bestHSPs;
}

*/

void zoeDeleteScoreSystem (zoeScoreSystem s) {
	int i;
	for (i = 0; i < 22; i++) zoeFree(s->score[i]);
	zoeFree(s->score);
	zoeFree(s);
}

zoeScoreSystem zoeNewScoreSystem (char* filename, int gapO, int gapE, int width) {
	int       i, j, a, b;
	float  ** matrix;
	FILE    * file;
	char      line[256];
	char    * symbol;
	char    * value;
	float     s;
	zoeTVec   symbols = zoeNewTVec();
	zoeFVec   values = zoeNewFVec();
	zoeScoreSystem ss = zoeMalloc(sizeof(struct zoeScoreSystem));
	
	/* read scoring matrix */	
	if ((file = fopen(filename, "r")) == NULL) zoeExit("file error %s\n", filename);
	
	matrix = malloc(22 * sizeof(float *));
	for (i = 0; i < 22; i++) matrix[i] = malloc(22 * sizeof(float));
	for (i = 0; i < 22; i++) {
		for (j = 0; j < 22; j++) matrix[i][j] = -100;
	}
	
	/* strip comments */
	while (fgets(line, sizeof(line), file)) {
		 if (line[0] == '#') continue;
		 break;
	}
	
	/* get list of symbols */
	symbol = strtok(line, " ");
	zoePushTVec(symbols, symbol);
	while ((symbol = strtok(NULL, " "))) {
		zoePushTVec(symbols, symbol);
	}
	
	/* values */
	while (fgets(line, sizeof(line), file)) {
		symbol = strtok(line, " ");
		while ((value = strtok(NULL, " "))) {
			if (isspace(value[0])) continue; /* last token is a space */
			s = atof(value);
			zoePushFVec(values, s);
		}
	}
	
	/* load matrix */
	for (i = 0; i < symbols->size; i++) {
		for (j = 0; j < symbols->size; j++) {
			a = zoe_char2aa(symbols->elem[i][0]);
			b = zoe_char2aa(symbols->elem[j][0]);
			if (a == -1 || b == -1) continue;   /* skip non-canonical */
			s = values->elem[i * symbols->size + j];
			matrix[a][b] = s;
		}
	}
	
	/* clean up */
	zoeDeleteTVec(symbols);
	zoeDeleteFVec(values);
	
	/* assignment */
	ss->score = matrix;
	ss->gapO  = gapO;
	ss->gapE  = gapE;
	ss->width = width;

	return ss;
}

/*

	Protein Alignment Functions

	. 0 1 2 3 i (p1)
	0 
	1
	2
	3
	j
	(p2)

*/


zoeHSP zoeProtAlignLocal(zoeProtein p1, zoeProtein p2, zoeScoreSystem s) {
	int         i, j;
	float       dscore, hscore, vscore, max_i, max_j, max_s;
	zoePCell ** matrix;
	zoeHSP      hsp = zoeNewHSP();
	char       *a1, *a2, *as;
	int        index;
	
	/* setup */
	matrix = zoeMalloc((p1->length +1) * sizeof(struct zoePCell));
	for (i = 0; i <= p1->length; i++) {
		matrix[i] = zoeMalloc((p2->length +1) * sizeof(struct zoePCell));
	}
	a1 = zoeMalloc(p1->length + p2->length +1);
	a2 = zoeMalloc(p1->length + p2->length +1);
	as = zoeMalloc(p1->length + p2->length +1);
	
	/* initialization */
	for (i = 0; i <= p1->length; i++) {
		matrix[i][0].score = 0;
		matrix[i][0].trace = 0;
		matrix[i][0].state = 0;
	}
	for (j = 0; j <= p2->length; j++) {
		matrix[0][j].score = 0;
		matrix[0][j].trace = 0;
		matrix[0][j].state = 0;
	}
	
	max_i = 0;
	max_j = 0;
	max_s = -s->gapO * (p1->length + p2->length);
	
	/* induction */
	for (i = 1; i <= p1->length; i++) {
		for (j = 1; j <= p2->length; j++) {
		
			if (s->width && abs(i-j) > s->width) continue;
			
			dscore  = s->score[(int)p1->s22[i-1]][(int)p2->s22[j-1]];
			dscore += matrix[i-1][j-1].score;
			
			vscore = matrix[i][j-1].score;
			if (matrix[i][j-1].state == 0) vscore -= s->gapO;
			else                           vscore -= s->gapE;
			
			hscore = matrix[i-1][j].score;
			if (matrix[i-1][j].state == 0) hscore -= s->gapO;
			else                           hscore -= s->gapE;
			
			if (dscore >= vscore && dscore >= hscore && dscore > 0) {
				matrix[i][j].score = dscore;
				matrix[i][j].trace = 1;
				matrix[i][j].state = 0;
				if (matrix[i][j].score > max_s) {
					max_s = matrix[i][j].score;
					max_j = j;
					max_i = i;
				}
			} else if (vscore >= hscore && vscore > 0) {
				matrix[i][j].score = vscore;
				matrix[i][j].trace = 3;
				matrix[i][j].state = 1;
			} else if (hscore > 0) {
				matrix[i][j].score = hscore;
				matrix[i][j].trace = 2;
				matrix[i][j].state = 1;
			} else {
				matrix[i][j].score = 0;
				matrix[i][j].trace = 0;
				matrix[i][j].state = 0;
			}
		}
	}
		
	/* trace back */
	i = max_i;
	j = max_j;
	index = 0;
	while (matrix[i][j].trace != 0) {
		switch (matrix[i][j].trace) {
			case 1: /* match state */
				i--;
				j--;
				a1[index] = p1->seq[i];
				a2[index] = p2->seq[j];
				if (a1[index] == a2[index]) {
					as[index] = a1[index];
				} else if (s->score[(int)p1->s22[i]][(int)p2->s22[j]] > 0) {
					as[index] = '+';
				} else {
					as[index] = ' ';
				}
				break;
			case 2: /* horizontal gap */
				i--;
				a1[index] = p1->seq[i];
				a2[index] = '-';
				as[index] = ' ';
				break;
			case 3: /* vertical gap */
				a1[index] = '-';
				a2[index] = p2->seq[j];
				as[index] = ' ';
				j--;
				break;
			default:
				zoeExit("impossible");
		}
		index++;		
	}
	a1[index] = '\0'; a2[index] = '\0'; as[index] = '\0';
	flip_string(a1);  flip_string(a2);  flip_string(as);
	
	/* decorations */
	hsp->score   = max_s;
	hsp->q_start = i;
	hsp->q_end   = max_i -1;
	hsp->s_start = j;
	hsp->s_end   = max_j -1;
	hsp->q_aln   = a1;
	hsp->s_aln   = a2;
	hsp->a_str   = as;
	
	hsp->match    = 0;
	hsp->mismatch = 0;
	hsp->q_gap    = 0;
	hsp->s_gap    = 0;
	for (i = 0; i < index; i++) {
		if (a1[i] == a2[i])    hsp->match++;
		else if (a1[i] == '-') hsp->q_gap++;
		else if (a2[i] == '-') hsp->s_gap++;
		else                   hsp->mismatch++;
	}
	
	/* clean up */	
	for (i = 0; i <= p1->length; i++) zoeFree(matrix[i]);
	zoeFree(matrix);
	
	return hsp;
}

zoeHSP zoeProtAlignNterm(zoeProtein p1, zoeProtein p2, zoeScoreSystem s) {
	int         i, j;
	float       dscore, hscore, vscore, max_i, max_j, max_s;
	zoePCell ** matrix;
	zoeHSP      hsp = zoeNewHSP();
	char       *a1, *a2, *as;
	int        index;
	
	/* setup */
	matrix = zoeMalloc((p1->length +1) * sizeof(struct zoePCell));
	for (i = 0; i <= p1->length; i++) {
		matrix[i] = zoeMalloc((p2->length +1) * sizeof(struct zoePCell));
	}
	a1 = zoeMalloc(p1->length + p2->length +1);
	a2 = zoeMalloc(p1->length + p2->length +1);
	as = zoeMalloc(p1->length + p2->length +1);
	
	/* initialization */
	matrix[0][0].score = 0;
	matrix[0][0].trace = 0;
	matrix[0][0].state = 0;
	for (i = 1; i <= p1->length; i++) {
		matrix[i][0].score = -s->gapO * i;
		matrix[i][0].trace = 0;
		matrix[i][0].state = 0;
	}
	for (j = 1; j <= p2->length; j++) {
		matrix[0][j].score = -s->gapO * i;
		matrix[0][j].trace = 0;
		matrix[0][j].state = 0;
	}
	
	max_i = 0;
	max_j = 0;
	max_s = -s->gapO * (p1->length + p2->length);
	
	/* induction */
	for (i = 1; i <= p1->length; i++) {
		for (j = 1; j <= p2->length; j++) {
		
			if (s->width && abs(i-j) > s->width) continue;
			
			dscore  = s->score[(int)p1->s22[i-1]][(int)p2->s22[j-1]];
			dscore += matrix[i-1][j-1].score;
			
			vscore = matrix[i][j-1].score;
			if (matrix[i][j-1].state == 0) vscore -= s->gapO;
			else                           vscore -= s->gapE;
			
			hscore = matrix[i-1][j].score;
			if (matrix[i-1][j].state == 0) hscore -= s->gapO;
			else                           hscore -= s->gapE;
			
			if (dscore >= vscore && dscore >= hscore) {
				matrix[i][j].score = dscore;
				matrix[i][j].trace = 1;
				matrix[i][j].state = 0;
				if (matrix[i][j].score > max_s) {
					max_s = matrix[i][j].score;
					max_j = j;
					max_i = i;
				}
			} else if (vscore >= hscore) {
				matrix[i][j].score = vscore;
				matrix[i][j].trace = 3;
				matrix[i][j].state = 1;
			} else {
				matrix[i][j].score = hscore;
				matrix[i][j].trace = 2;
				matrix[i][j].state = 1;
			}
		}
	}
	
	
	/* trace back */
	i = max_i;
	j = max_j;
	index = 0;
	while (matrix[i][j].trace != 0) {
		switch (matrix[i][j].trace) {
			case 1: /* match state */
				i--;
				j--;
				a1[index] = p1->seq[i];
				a2[index] = p2->seq[j];
				if (a1[index] == a2[index]) {
					as[index] = a1[index];
				} else if (s->score[(int)p1->s22[i]][(int)p2->s22[j]] > 0) {
					as[index] = '+';
				} else {
					as[index] = ' ';
				}
				break;
			case 2: /* horizontal gap */
				i--;
				a1[index] = p1->seq[i];
				a2[index] = '-';
				as[index] = ' ';
				break;
			case 3: /* vertical gap */
				a1[index] = '-';
				a2[index] = p2->seq[j];
				as[index] = ' ';
				j--;
				break;
			default:
				zoeExit("impossible");
		}
		index++;		
	}
	a1[index] = '\0'; a2[index] = '\0'; as[index] = '\0';
	flip_string(a1);  flip_string(a2);  flip_string(as);
	
	/* decorations */
	hsp->score   = max_s;
	hsp->q_start = i;
	hsp->q_end   = max_i -1;
	hsp->s_start = j;
	hsp->s_end   = max_j -1;
	hsp->q_aln   = a1;
	hsp->s_aln   = a2;
	hsp->a_str   = as;
	
	hsp->match    = 0;
	hsp->mismatch = 0;
	hsp->q_gap    = 0;
	hsp->s_gap    = 0;
	for (i = 0; i < index; i++) {
		if (a1[i] == a2[i])    hsp->match++;
		else if (a1[i] == '-') hsp->q_gap++;
		else if (a2[i] == '-') hsp->s_gap++;
		else                   hsp->mismatch++;
	}
	
	/* clean up */	
	for (i = 0; i <= p1->length; i++) zoeFree(matrix[i]);
	zoeFree(matrix);
	
	return hsp;
}


#endif
