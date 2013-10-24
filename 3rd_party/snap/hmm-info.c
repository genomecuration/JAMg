/*****************************************************************************\
 hmm-info.c

HMM information output utility

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

