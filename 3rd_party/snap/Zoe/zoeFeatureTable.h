/******************************************************************************\
zoeFeatureTable.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_FEATURETABLE_H
#define ZOE_FEATURETABLE_H

#include <stdio.h>
#include <string.h>

#include "zoeCDS.h"
#include "zoeDNA.h"
#include "zoeFeature.h"
#include "zoeTools.h"

struct zoeFeatureTable  {
	char          * def;
	zoeFeatureVec   vec;     /* storage for the actual features */
	zoeHash         region;  /* feature indexing by large-ish region */
	zoeTVec         regions; /* names of indicies */
};
typedef struct zoeFeatureTable * zoeFeatureTable;

void            zoeDeleteFeatureTable (zoeFeatureTable);
zoeFeatureTable zoeNewFeatureTable (const char *, const zoeFeatureVec);
zoeFeatureTable zoeReadFeatureTable (FILE *);
void            zoeWriteFeatureTable (FILE *, const zoeFeatureTable);
void            zoeWriteTriteFeatureTable (FILE *, const zoeFeatureTable);
void            zoeAddFeature (zoeFeatureTable, const zoeFeature);
void            zoeAddFeatures (zoeFeatureTable, const zoeFeatureVec);
zoeFeatureTable zoeSelectExons (const zoeFeatureTable);
zoeFeatureTable zoeSelectByGroup (const char *, const zoeFeatureTable, const char *);
zoeFeatureTable zoeSelectByLabel (const char *, const zoeFeatureTable, zoeLabel);
void            zoeAntiFeatureTable (zoeFeatureTable, coor_t length);
zoeFeatureTable zoeGetFeatureTable (const char *);
zoeTVec         zoeFeatureTableGroups (const zoeFeatureTable);
zoeVec          zoeGetGenes (const zoeFeatureTable, const zoeDNA); /* vector of zoeCDS */
zoeFeatureVec   zoeGetFeaturesNear(const zoeFeatureTable, coor_t, coor_t);
void            zoePadFeatureTable(zoeFeatureTable, int);

#endif

