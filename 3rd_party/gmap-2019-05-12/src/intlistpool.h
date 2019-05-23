/* $Id: intlistpool.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef INTLISTPOOL_INCLUDED
#define INTLISTPOOL_INCLUDED

typedef struct Intlistpool_T *Intlistpool_T;

#include "types.h"
#include "intlist.h"

#define T Intlistpool_T

extern void
Intlistpool_free (T *old);
extern void
Intlistpool_free_memory (T this);
extern T
Intlistpool_new (void);
extern void
Intlistpool_reset (T this);
extern Intlist_T
Intlistpool_push (Intlist_T intlist, T this, int integer);
extern Intlist_T
Intlistpool_pop (Intlist_T intlist, int *integer);
extern Intlist_T
Intlistpool_copy (Intlist_T source, T this);
extern Intlist_T
Intlistpool_copy_but_last (Intlist_T source, T this);

#undef T
#endif
