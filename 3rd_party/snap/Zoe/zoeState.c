/******************************************************************************\
 zoeState.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_STATE_C
#define ZOE_STATE_C

#include "zoeState.h"

void zoeDeleteState (zoeState s) {
	if (s == NULL) return;
	zoeFree(s);
	s = NULL;
}

zoeState zoeNewState (
	zoeStateType type, zoeLabel label, float i, float t, int min, int max, int geo)
{
	zoeState s = zoeMalloc(sizeof(struct zoeState));
	
	s->type      = type;
	s->label     = label;
	s->init      = i;
	s->term      = t;
	s->min       = min;
	s->max       = max;
	s->geometric = geo;
	
	return s;
}

zoeState zoeReadState (FILE * stream) {
	char         state_name[16];
	char         duration[16];
	float        init, term;
	int          min, max, geometric = 0;
	zoeStateType type;
	zoeLabel     label;
	zoeState     state = NULL;

	if (fscanf(stream, "%s %f %f %d %d %s", state_name, &init, &term, &min, &max, duration) != 6) {
		zoeWarn("zoeReadState header");
		return NULL;
	}
	
	/* name */
	label = zoeText2Label(state_name);
	if (label == None) {
		zoeWarn("zoeReadState name not allowed (%s)", state_name);
		return NULL;
	}
		
	/* assign state type */
	switch (label) {
		case Inter: case Intron: case UTR5: case UTR3:
			type = INTERNAL;
			break;
		case Esngl: case Einit: case Eterm: case Exon:
			type = EXTERNAL;
			break;
		case TSS: case PolyA: case Prom:
			type = EXTERNAL;
			break;
		case Repeat: case CNS: case ORF:
			type = SHUTTLE;
			break;
		default:
			zoeWarn("zoeReadState no type binding for %s", state_name);
			return NULL;
	}
	
	/* assign geometric */
	if (strcmp(duration, "geometric") == 0) {
		geometric = 1;
		if (min || max) zoeExit("geometric min/max must be == 0");
	} else if (strcmp(duration, "explicit") == 0)  {
		geometric = 0;
		if (min < 1) zoeExit("min state length < 1 makes no sense");
		if (max == -1) max = INT_MAX;
	} else zoeExit("must be geometric or explicit in zoeReadState");
	
	state = zoeNewState(type, label, init, term, min, max, geometric);
	return state;
}

void zoeWriteState (FILE * stream, const zoeState state) {
	char name[16];

	zoeLabel2Text(state->label, name);
	zoeS(stream, "%s\t%f\t%f\t%d\n", name, state->init, state->term, state->min, state->max);
}

#endif
