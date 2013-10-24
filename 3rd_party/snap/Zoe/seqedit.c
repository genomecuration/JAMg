/*****************************************************************************\
seqedit.c

Sequence editor

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

