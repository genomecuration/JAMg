/* $Id: univdiagpool.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef UNIVDIAGPOOL_INCLUDED
#define UNIVDIAGPOOL_INCLUDED

#include "univdiag.h"
#include "list.h"
#include "univcoord.h"


#define T Univdiagpool_T
typedef struct T *T;

extern void
Univdiagpool_free (T *old);
extern void
Univdiagpool_free_memory (T this);
extern void
Univdiagpool_report_memory (T this);
extern T
Univdiagpool_new (void);
extern void
Univdiagpool_reset (T this);
extern Univdiag_T
Univdiag_new (T this, int qstart, int qend, Univcoord_T univdiagonal);
extern List_T
Univdiagpool_push (List_T list, T this, int qstart, int qend, Univcoord_T univdiagonal);
extern List_T
Univdiagpool_pop (List_T list, Univdiag_T *x);
extern List_T
Univdiagpool_push_existing (List_T list, T this, Univdiag_T univdiag);


#undef T
#endif
