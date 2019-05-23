/* $Id: reader.h 218153 2019-01-17 05:38:29Z twu $ */
#ifndef READER_INCLUDED
#define READER_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "types.h"

/* Still used by pair.c and chimera.c */
typedef enum {FIVE, THREE, MIDDLE} cDNAEnd_T;

#define T Reader_T
typedef struct T *T;

extern int
Reader_querystart (T this);
extern int
Reader_queryend (T this);
extern int
Reader_startpos (T this);
extern int
Reader_endpos (T this);
extern void
Reader_reset_start (T this, int querypos);
extern void
Reader_reset_end (T this, int querypos);
extern void
Reader_reset_ends (T this);

extern T
Reader_new (char *sequence, int querystart, int queryend);
extern void
Reader_free (T *old);

#ifndef GSNAP
extern char
Reader_getc (T this, cDNAEnd_T cdnaend, int blocksize);
#endif

extern char
Reader_getc_5 (T this);
extern char
Reader_getc_3 (T this);

extern Oligospace_T
Reader_check (T this, int querypos, int indexsize);


#undef T
#endif
