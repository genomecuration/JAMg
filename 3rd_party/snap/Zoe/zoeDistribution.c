/******************************************************************************\
zoeDuration.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_DISTRIBUTION_C
#define ZOE_DISTRIBUTION_C

#include "zoeDistribution.h"

void zoeDeleteDistribution(zoeDistribution dist) {	
	if (dist == NULL) return;
	if (dist->param) {
		zoeFree(dist->param);
		dist->param = NULL;
	}
	zoeFree(dist);
	dist = NULL;
}

zoeDistribution zoeNewDistribution (
	zoeDistributionType   type,
	coor_t                start,
	coor_t                end,
	int                   params,
	const float         * param)
{
	zoeDistribution d = zoeMalloc(sizeof(struct zoeDistribution));
	
	if (end == -1) end = INT_MAX;
	
	d->type   = type;
	d->start  = start;
	d->end    = end;
	d->params = params;
	d->param  = zoeMalloc(params * sizeof(float));
	(void)memcpy(d->param, param, params * sizeof(float));
	return d;
}

zoeDistribution zoeReadDistribution (FILE * stream) {
	char                  func_name[32];
	char                  input[32];
	int                   start, end;
	int                   params, k;
	float               * param;
	zoeDistributionType   type;
	zoeDistribution       dist;

	if (fscanf(stream, "%s %d %d", func_name, &start, &end) != 3) {
		zoeWarn("fscanf zoeReadDistribution func");
		return NULL;
	}
				
	params = -1; /* impossible value */
	 if (strcmp(func_name, "DEFINED") == 0) {
	 	 type = DEFINED;
	 	 params = end - start + 1;
	 } else if (strcmp(func_name, "GEOMETRIC") == 0) {
	 	 type = GEOMETRIC;
	 	 params = 1;
	 } else if (strcmp(func_name, "POISSON") == 0) {
	 	 type = POISSON;
	 	 params = 1;
	 } else if (strcmp(func_name, "CONSTANT") == 0) {
	 	 type = CONSTANT;
	 	 params = 1;
	 } else {
	 	 zoeWarn("zoeReadDistribution unknown func_name (%s)", func_name);
		 return NULL;
	 }
	 					 
	 param = zoeMalloc(params * sizeof(float));
	
	 for (k = 0; k < params; k++) {
	 	if (fscanf(stream, "%s", input) != 1) {
	 		zoeWarn("fscanf zoeReadDistribution params");
			 return NULL;
	 	}
	 	param[k] = zoeText2Score(input);
	 }
	 
	 dist = zoeNewDistribution(type, start, end, params, param);
	 zoeFree(param);
	 return dist;
}

void zoeWriteDistribution (FILE * stream, const zoeDistribution dist) {
	int k;
	char score[16];

	switch (dist->type) {
		case DEFINED:   zoeS(stream, "\tDEFINED");   break;
		case GEOMETRIC: zoeS(stream, "\tGEOMETRIC"); break;
		case POISSON:   zoeS(stream, "\tPOISSON");   break;
		case CONSTANT:  zoeS(stream, "\tCONSTANT");  break;
		default: zoeExit("Unknown distribution type");
	}
	zoeS(stream, " %d", dist->start);
	if (dist->end == INT_MAX) zoeS(stream, " -1\t");
	else                      zoeS(stream, " %d\t", dist->end);
	for (k = 0; k < dist->params; k++) {
		if ((k % 5) == 0) zoeS(stream, "\n\t\t");
		zoeScore2Text(dist->param[k], score);
		zoeS(stream, "%s\t", score);
	}
	zoeS(stream, "\n");
}


score_t zoeScoreDistribution(const zoeDistribution dist, coor_t pos) {
	switch (dist->type) {
		case GEOMETRIC:
			return zoeScoreGeometric((double)dist->param[0], (double)pos);
		case POISSON:
			return zoeScorePoisson((double)dist->param[0], (double)pos);
		case CONSTANT:
			return dist->param[0];
		case DEFINED:
			if (pos >= dist->start && pos <= dist->end)
				return dist->param[pos - dist->start];
			else zoeExit("zoeScoreDistribution out of bounds");
		default: zoeExit("zoeScoreDistribution unknown type");
	}
	return MIN_SCORE; /* shush compiler warning */
}

#endif
