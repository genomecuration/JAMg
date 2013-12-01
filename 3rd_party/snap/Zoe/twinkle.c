/*****************************************************************************\
twinkle.c

Pairwise alignment with multiple tracebacks

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



