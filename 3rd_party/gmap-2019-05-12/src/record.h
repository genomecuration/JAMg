/* $Id: record.h 216937 2018-10-09 23:00:01Z twu $ */
#ifndef RECORD_INCLUDED
#define RECORD_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bool.h"
#include "types.h"
#include "chrnum.h"

typedef struct Record_T *Record_T;
struct Record_T {
  Univcoord_T diagonal;		/* Primary sort */
  int querypos;			/* Secondary sort */
  int queryend;
  bool anchorp;

  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;

  Univcoord_T lowpos;
  Univcoord_T highpos;
};

#endif

