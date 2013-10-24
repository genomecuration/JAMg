/******************************************************************************\
zoeProtein.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2013 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_PROTEIN_H
#define ZOE_PROTEIN_H

#include <stdlib.h>
#include <string.h>

#include "zoeFastaFile.h"
#include "zoeTools.h"

struct zoeProtein  {
	coor_t   length;
	char   * def;
	char   * seq;
	char   * s22; /* 20 amino acids plus X and stop */
};
typedef struct zoeProtein * zoeProtein;

int        zoe_char2aa (const int);
char       zoe_aa2char (const int);
void       zoeDeleteProtein (zoeProtein);
zoeProtein zoeNewProtein (const char *, const char *);
zoeProtein zoeNewTrustedProtein (const char *, const char *);
void       zoeWriteProtein (FILE *, const zoeProtein);
zoeProtein zoeGetProtein (const char *);

#endif
