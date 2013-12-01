/*****************************************************************************\
 hmm-info.c

HMM information output utility

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

void modelInfo (const zoeHMM);
void generalInfo (const zoeHMM);
void exportDurations (const zoeHMM);

/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char *argv[]) {
	zoeHMM          hmm;
	
	/* set the program name */
	zoeSetProgramName(argv[0]);
	zoeSetOption("-durations", 1);
	zoeSetOption("-general",   0);
	zoeSetOption("-models", 0);
	zoeParseOptions(&argc, argv);
	
	/* usage */
	if (argc != 2) {
		zoeE("usage: %s <hmm file>\n", argv[0]);
		zoeM(stderr, 4,
			"commands:",
			"  -models",
			"  -general",
			"  -durations <length>"
		);
		exit(1);
	}

	/* get HMM */
	if ((hmm = zoeGetHMM(argv[1])) == NULL) zoeExit("error opening hmm file");
	
	if (zoeOption("-models"))    modelInfo(hmm);
	if (zoeOption("-general"))   generalInfo(hmm);
	if (zoeOption("-durations")) exportDurations(hmm);
	
	zoeDeleteHMM(hmm);
	return 0;
}

void modelInfo (const zoeHMM hmm) {
	int i;
	for (i = 0; i < hmm->models; i++) {
		zoeWriteModel(stdout, hmm->model[i]);
	}	
}


void generalInfo (const zoeHMM hmm) {
	int i;
	
	/* models */
	for (i = 0; i < hmm->models; i++) {
		zoeWriteModelHeader(stdout, hmm->model[i]);
	}

}


void exportDurations (const zoeHMM hmm) {
	int  label, i;
	char name[16];
	int  length;
	
	length = atoi(zoeOption("-durations"));
	
	/* should have done this simpler */
	for (label = 0; label < zoeLABELS; label++) {
		if (hmm->dmap[label] == NULL) continue;
		
		switch (label) {
			case Int0:
				strcpy(name, "Intron");
				break;
			case Int1: case Int1T: case Int2: case Int2TA: case Int2TG:
				continue;
			default:
				zoeLabel2Text(label, name);
		}
		
		zoeO("%s\t", name);
	}
	zoeO("\n");
	
	
	for (i = 1; i < length; i++) {
		for (label = 0; label < zoeLABELS; label++) {
			if (hmm->dmap[label] == NULL) continue;
			if (label == Int1 || label == Int1T || label == Int2
				|| label == Int2TA || label == Int2TG) continue;
			
			if (i < hmm->smap[label]->min) {
				zoeO("0.000000\t");
			} else {
				zoeO("%f\t", zoeScore2Float(zoeScoreDuration(hmm->dmap[label], i)));
			}
		}
		zoeO("\n");
	}
	
}

