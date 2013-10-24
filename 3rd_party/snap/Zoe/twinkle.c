/*****************************************************************************\
twinkle.c

Pairwise alignment with multiple tracebacks

Copyright (C) 2002-2013 Ian Korf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zoe.h"

zoeVec zoeAlign02 (const zoeDNA, const zoeDNA, int, int, int, int, int);
zoeVec zoeChainHSPs (const zoeVec);

static int MAX_MEMORY = 100; /* default is 100 meg */

int main (int argc, char *argv[]) {	
	zoeDNA dna1, dna2;
	zoeVec HSPs;
	int i, match = 1, mismatch = 1, gapI = 3, gapE = 1, score = 10;
	size_t mem;
	
	zoeSetProgramName(argv[0]);
	zoeSetOption("-match",      1);
	zoeSetOption("-mismatch",   1);
	zoeSetOption("-gap-init",   1);
	zoeSetOption("-gap-extend", 1);
	zoeSetOption("-score",      1);
	zoeSetOption("-debug",      0);
	zoeSetOption("-memory",     1);
	zoeParseOptions(&argc, argv);
	
	if (argc != 3) {
		zoeE("usage: %s [options] <dna1> <dna2>\n", argv[0]);
		zoeM(stderr, 8,
			"options:                      [default]",
			"  -match      <positive int>  [1]",
			"  -mismatch   <positive int>  [1]",
			"  -gap-init   <positive int>  [3]",
			"  -gap-extend <positive int>  [1]",
			"  -score      <positive int>  [10]",
			"  -memory     <positive int>  [100]  (units are megabytes)",
			"  -debug"
		);
		exit(1);
	}
	
	dna1 = zoeGetDNA(argv[1]);
	dna2 = zoeGetDNA(argv[2]);
	
	if (zoeOption("-match"))      match = atoi(zoeOption("-match"));
	if (zoeOption("-mismatch"))   mismatch = atoi(zoeOption("-mismatch"));
	if (zoeOption("-gap-init"))   gapI = atoi(zoeOption("-gap-init"));
	if (zoeOption("-gap-extend")) gapE = atoi(zoeOption("-gap-extend"));
	if (zoeOption("-score"))      score = atoi(zoeOption("-score"));
	if (zoeOption("-memory"))     MAX_MEMORY = atoi(zoeOption("-memory"));
	
	if (match < 1 || mismatch < 1 || gapI < 1 || gapE < 1 || score < 1) {
		zoeExit("all parameters must be positive integers");
	}
	
	mem = (dna1->length +1) * (dna2->length +1) * sizeof(zoeAlignCell02) / 1000000;
	if (mem > MAX_MEMORY) {
		zoeExit("maximum memory limit (%dM) exceeded (%dM)", MAX_MEMORY, mem);
	}
	
	HSPs = zoeAlign02(dna1, dna2, match, mismatch, gapI, gapE, score);
	
	for (i = 0; i < HSPs->size; i++) zoeWriteHSP(stdout, HSPs->elem[i]);
	
	return 0;
}



