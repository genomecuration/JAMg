/******************************************************************************\
zoePhasePref.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_PHASEPREF_C
#define ZOE_PHASEPREF_C

#include "zoePhasePref.h"

void zoeDeletePhasePref (zoePhasePref pp) {	
	if (pp == NULL) return;
	zoeFree(pp);
	pp = NULL;
}

zoePhasePref zoeNewPhasePref (void) {
	int          i;
	zoePhasePref pp = zoeMalloc(sizeof(struct zoePhasePref));
	
	for (i = 0; i < zoePHASECOUNT; i++) {
		pp->score[i] = MIN_SCORE;
		pp->prob[i]  = 0;
	}
	return pp;
}

zoePhasePref zoeReadPhasePref (FILE * stream) {
	int          i;
	zoePhasePref pp = zoeNewPhasePref();
	
	for (i = 0; i < zoePHASECOUNT; i++) {
		if (fscanf(stream, "%f", &pp->prob[i]) != 1) {
			zoeWarn("zoeReadPhasePref fscanf error");
			zoeDeletePhasePref(pp);
			return NULL;
		}
		/* convert to score too */
		pp->score[i] = zoeFloat2Score(pp->prob[i]);
	}
	return pp;
}

void zoeWritePhasePref (FILE *stream, const zoePhasePref pp) {
	int i;
	
	for (i = 0; i < zoePHASECOUNT; i++) {
		zoeS(stream, "%f\n", pp->prob[i]);
	}
}

score_t zoeScorePhase (zoePhasePref pp, zoeLabel from, zoeLabel to, int inc5) {
	
	switch (from) {
		case Einit:
			switch (to) {
				case Int0:                           return pp->score[Ei_I0];
				case Int1: case Int1T:               return pp->score[Ei_I1];
				case Int2: case Int2TA: case Int2TG: return pp->score[Ei_I2];
				default: zoeExit("zoeScorePhase does not allow to %d", to);
			}
		case Exon:
			switch (inc5) {
				case 0: /* corresponds to E0 */
					switch (to) {
						case Int0:                           return pp->score[E0_I0];
						case Int1: case Int1T:               return pp->score[E0_I1];
						case Int2: case Int2TA: case Int2TG: return pp->score[E0_I2];	
						default: zoeExit("zoeScorePhase impossible A %d", to);
					}
				case 1: /* corresponds to E2 */
					switch (to) {
						case Int0:                           return pp->score[E2_I0];
						case Int1: case Int1T:               return pp->score[E2_I1];
						case Int2: case Int2TA: case Int2TG: return pp->score[E2_I2];
						default: zoeExit("zoeScorePhase impossible B %d", to);
					}
				case 2: /* corresponds to E1 */
					switch (to) {
						case Int0:                           return pp->score[E1_I0];
						case Int1: case Int1T:               return pp->score[E1_I1];
						case Int2: case Int2TA: case Int2TG: return pp->score[E1_I2];
						default: zoeExit("zoeScorePhase impossible C %d", to);
					}
				default: zoeExit("zoeScorePhase impossible D");
			}
		case Int0:
			switch (to) {
				case Exon:  return pp->score[I0_E0];
				case Eterm: return pp->score[I0_Et];
				default: zoeExit("zoeScorePhase does not allow to %d", to);
			}
		case Int1: case Int1T:
			switch (to) {
				case Exon:  return pp->score[I1_E1];
				case Eterm: return pp->score[I1_Et];
				default: zoeExit("zoeScorePhase does not allow to %d", to);
			}
		case Int2: case Int2TA: case Int2TG:
			switch (to) {
				case Exon:  return pp->score[I2_E2];
				case Eterm: return pp->score[I2_Et];
				default: zoeExit("zoeScorePhase does not allow to %d", to);
			}
		default: return 0;
	}
}

#endif
