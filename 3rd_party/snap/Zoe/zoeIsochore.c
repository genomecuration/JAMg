/******************************************************************************\
zoeIsochore.c - part of the ZOE library for genomic analysis
 
Copyright (C) 2002-2013 Ian Korf

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

\******************************************************************************/

#ifndef ZOE_ISOCHORE_C
#define ZOE_ISOCHORE_C

#include "zoeHMM.h"
#include "zoeIsochore.h"

void zoeDeleteIsochore (zoeIsochore iso) {
	int i;
	zoeHMM hmm;
	
	if (iso->min_GC) {
		zoeDeleteFVec(iso->min_GC);
		iso->min_GC = NULL;
	}
	if (iso->max_GC) {
		zoeDeleteFVec(iso->max_GC);
		iso->max_GC = NULL;
	}
	if (iso->hmm_file) {
		zoeDeleteTVec(iso->hmm_file);
		iso->hmm_file = NULL;
	}
	if (iso->hmms) {
		for (i = 0; i < iso->hmms->size; i++) {
			hmm = iso->hmms->elem[i];
			zoeDeleteHMM(hmm);
			hmm = NULL;
		}
		zoeDeleteVec(iso->hmms);
	}
}

zoeIsochore zoeNewIsochore (void) {
	zoeIsochore iso = zoeMalloc(sizeof(struct zoeIsochore));
	
	iso->count    = 0;
	iso->min_GC   = zoeNewFVec();
	iso->max_GC   = zoeNewFVec();
	iso->hmm_file = zoeNewTVec();
	iso->hmms     = zoeNewVec();

	return iso;
}


zoeIsochore zoeReadIsochore (FILE * stream) {
	char        name[256];
	int         i, count;
	float       min, max;
	zoeIsochore iso = zoeNewIsochore();
	zoeHMM      hmm;
		
	if (fscanf(stream, "%s %d", name, &count) != 2)
		zoeExit("zoeReadIsochore failed to read header");
	if (strcmp(name, "zoeIsochore") != 0)
		zoeExit("zoeReadIsochore found an unrecognized file type");
		
	iso->count = count;
	for (i = 0; i < count; i++) {
		if (fscanf(stream, "%f %f %s", &min, &max, name) != 3) zoeExit("zoeReadIsochore format error");
		zoePushFVec(iso->min_GC, min);
		zoePushFVec(iso->max_GC, max);
		zoePushTVec(iso->hmm_file, name);
		hmm = zoeGetHMM(name);
		zoePushVec(iso->hmms, hmm);
	}
	
	return iso;
}

zoeIsochore zoeGetIsochore (const char * file) {
	FILE        * stream = NULL;
	zoeIsochore   iso    = NULL;
	char        * ZOE    = getenv("ZOE");
	char          path[1024];
	
	stream = fopen(file, "r");
	if (stream == NULL) {
		sprintf(path, "%s/HMM/%s", ZOE, file);
		stream = fopen(path, "r");
		if (stream == NULL) {
			zoeExit("error opening isochore file");
		}
	}
	
	iso = zoeReadIsochore(stream);
	if (iso == NULL) zoeExit("error reading isochore file");
	
	(void)fclose(stream);
	return(iso);
}

zoeHMM zoeSelectIsochore (const zoeIsochore iso, float GC_fraction) {
	int    i;
	float  min, max;
	zoeHMM hmm;
		
	for (i = 0; i < iso->hmms->size; i++) {
		min = iso->min_GC->elem[i];
		max = iso->max_GC->elem[i];
		hmm = iso->hmms->elem[i];
		if (GC_fraction > min && GC_fraction < max) return hmm;
	}
	
	zoeExit("zoeSelectIsochore could not find a valid isochore for %f", GC_fraction);
	return NULL;
}

#endif
