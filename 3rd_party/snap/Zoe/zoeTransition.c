/******************************************************************************\
 zoeTransition.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_TRANSITION_C
#define ZOE_TRANSITION_C

#include "zoeTransition.h"

void zoeDeleteTransition (zoeTransition t) {
	if (t == NULL) return;
	zoeFree(t);
	t = NULL;
}

zoeTransition zoeNewTransition (const char * from, const char * to, float prob) {
	zoeTransition t = zoeMalloc(sizeof(struct zoeTransition));
	
	t->from  = zoeText2Label(from);
	t->to    = zoeText2Label(to);
	t->prob  = prob;
	t->score = zoeFloat2Score(prob);
	return t;
}

zoeTransition zoeReadTransition (FILE * stream) {
	char  from[16];
	char  to[16];
	float prob;
		
	if (fscanf(stream, "%s %s %f", from, to, &prob) != 3) {
		zoeWarn("fscanf error in zoeReadTransition");
		return NULL;
	}
	return zoeNewTransition(from, to, prob);
}

void zoeWriteTransition (FILE * stream, const zoeTransition t) {
	char from[16];
	char to[16];
		
	zoeLabel2Text(t->from, from);
	zoeLabel2Text(t->to,   to);
	zoeS(stream, "%s\t%s\t%f\n", from, to, t->prob);
}

#endif
