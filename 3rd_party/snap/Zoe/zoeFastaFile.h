/******************************************************************************\
zoeFastaFile.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_FASTA_FILE_H
#define ZOE_FASTA_FILE_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeTools.h"

struct zoeFastaFile  {
	coor_t   length;
	char   * def;
	char   * seq;
};
typedef struct zoeFastaFile * zoeFastaFile;

void          zoeDeleteFastaFile (zoeFastaFile);
zoeFastaFile  zoeNewFastaFile (const char *, const char *);
zoeFastaFile  zoeReadFastaFile (FILE *);
zoeFastaFile  zoeGetFastaFile (const char *);
void          zoeSetFastaLineLength (unsigned int);
void          zoeWriteFastaFile (FILE *, const zoeFastaFile);
zoeFastaFile  zoeGetFastaFile(const char *);

#endif
