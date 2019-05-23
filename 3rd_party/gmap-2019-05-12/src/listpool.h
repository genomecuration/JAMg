/* $Id: listpool.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef LISTPOOL_INCLUDED
#define LISTPOOL_INCLUDED

typedef struct Listpool_T *Listpool_T;

#include "list.h"

#define T Listpool_T

extern void
Listpool_free (T *old);
extern void
Listpool_free_memory (T this);
extern T
Listpool_new (void);
extern void
Listpool_reset (T this);
extern List_T
Listpool_push (List_T list, T this, void *contents);
extern List_T
Listpool_pop (List_T list, void **contents);
extern List_T
Listpool_copy (List_T source, T this);

#undef T
#endif
