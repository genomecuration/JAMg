/******************************************************************************\
zoeCDS.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_CDS_H
#define ZOE_CDS_H

#include <stdio.h>
#include <string.h>

#include "zoeDNA.h"
#include "zoeFeature.h"
#include "zoeProtein.h"
#include "zoeTools.h"

struct zoeCDS  {

	/* parent DNA */
	zoeDNA     dna;
	
	/* all genes */
	char          * name;        /* gene name */
	coor_t          start;
	coor_t          end;
	strand_t        strand;
	score_t         score;
	zoeFeatureVec   source;      /* original features */
	zoeFeatureVec   exons;       /* exons */
	zoeFeatureVec   introns;     /* introns */
	zoeDNA          tx;          /* transcript */
	zoeProtein      aa;          /* protein */
	
	/* partial genes */
	int     start_found;
	int     end_found;
	frame_t inc5;      /* number of incomplete codon nucleotides on 5' end */
	frame_t inc3;      /*                                           3' end */	
	
	/* error checking */
	int OK; /* no errors or warnings */
	zoeTVec errors;
	zoeTVec warnings;
	
};
typedef struct zoeCDS * zoeCDS;

void    zoeDeleteCDS (zoeCDS);
zoeCDS  zoeNewCDS (const char *, const zoeDNA, const zoeFeatureVec);
void    zoeAntiCDS (zoeCDS, coor_t);
void    zoeWriteCDS (FILE *, const zoeCDS);
void    zoeWriteFullCDS (FILE *, const zoeCDS);
void    zoeWriteTriteCDS (FILE *, const zoeCDS);
void    zoeReportCDS (FILE *, const zoeCDS);
int     zoeCDScmp (const zoeCDS, const zoeCDS);
int     zoeCDScmpptr (const void *, const void *);
int     zoeCDSsOverlap (const zoeCDS, const zoeCDS);
int     zoeCDSsShareSequence (const zoeCDS, const zoeCDS);
void    zoeSetMinIntron (int);
void    zoeSetMaxIntron (int);
void    zoeSetMinExon (int);
void    zoeSetMaxExon (int);
void    zoeSetMinGene (int);
void    zoeSetMaxGene (int);
void    zoeSetMinCDS (int);
void    zoeSetMaxCDS (int);

#endif

