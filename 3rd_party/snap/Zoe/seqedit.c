/*****************************************************************************\
seqedit.c

Sequence editor

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

int main (int argc, char *argv[]) {	
	FILE         * stream;
	char         * filename;
	zoeFastaFile   fasta;
	zoeDNA         dna = NULL;
	zoeDNA         tmp = NULL;
	zoeProtein     pro = NULL;
	coor_t         cut, length;
	
	zoeSetProgramName(argv[0]);
	zoeSetOption("-reverse", 0);
	zoeSetOption("-complement", 0);
	zoeSetOption("-cut", 1);
	zoeSetOption("-length", 1);
	zoeSetOption("-translate", 0);
	zoeSetOption("-lcfilter", 0);
	zoeSetOption("-lcsmooth", 0);
	zoeSetOption("-flank", 1);
	zoeSetOption("-island", 1);
	zoeSetOption("-min-len", 1);
	zoeParseOptions(&argc, argv);
	
	if (argc == 1) {
		zoeS(stdout, "usage: %s [options] <filename or - for stdin>\n", argv[0]);
		zoeM(stderr, 7,
			"options:",
			"  -reverse",
			"  -complement",
			"  -cut <pos> -length <pos>",
			"  -translate",
			"  -lcfilter",
			"  -lcsmooth -flank <int> -island <int> -min-len <int>");
		exit(1);
	}
	filename = argv[1];
	
	/* file or stdin? */
	if ((strcmp(filename, "-") == 0)) stream = stdin;
	else if ((stream = fopen(filename, "r") ) == NULL)
			zoeExit("file error (%s)", filename);
	
	fasta = zoeReadFastaFile(stream);
	dna   = zoeNewDNA(fasta->def, fasta->seq);
	
	/* process the sequence */
	
	if (zoeOption("-lcsmooth")) {
		if (!zoeOption("-flank")) zoeExit("-flank required for -lcsmooth");
		if (!zoeOption("-island")) zoeExit("-island required for -lcsmooth");
		if (!zoeOption("-min-len")) zoeExit("-min-len required for -lcsmooth");
		zoeLCsmooth(dna, atoi(zoeOption("-flank")), atoi(zoeOption("-island")),
			atoi(zoeOption("-min-len")));
	}
	
	if (zoeOption("-lcfilter")) zoeLCfilter(dna);
	
	if (zoeOption("-reverse")) {
		tmp = zoeReverseDNA(dna->def, dna);
		zoeDeleteDNA(dna);
		dna = tmp;
		tmp = NULL;
	}
	
	if (zoeOption("-complement")) {
		tmp = zoeComplementDNA(dna->def, dna);
		zoeDeleteDNA(dna);
		dna = tmp;
		tmp = NULL;
	}
	
	if (zoeOption("-cut") || zoeOption("-length")) {
		if (zoeOption("-cut") && zoeOption("-length")) {
			cut    = atoi(zoeOption("-cut"));
			length = atoi(zoeOption("-length"));
			tmp = zoeSubseqDNA(dna->def, dna, cut -1, length);
			zoeDeleteDNA(dna);
			dna = tmp;
			tmp = NULL;
		} else {
			zoeExit("you must specify both -cut and -length");
		}
	}
	
	if (zoeOption("-translate")) {
		pro = zoeTranslateDNA(dna->def, dna, 0);
		zoeWriteProtein(stdout, pro);
		zoeDeleteProtein(pro);
	} else {
		zoeWriteDNA(stdout, dna);
		zoeDeleteDNA(dna);
	}
	
	return 0;
}

