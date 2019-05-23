/* $Id: uintlistpool.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef UINTLISTPOOL_INCLUDED
#define UINTLISTPOOL_INCLUDED

typedef struct Uintlistpool_T *Uintlistpool_T;

#include "types.h"
#include "uintlist.h"

#define T Uintlistpool_T

extern void
Uintlistpool_free (T *old);
extern void
Uintlistpool_free_memory (T this);
extern T
Uintlistpool_new (void);
extern void
Uintlistpool_reset (T this);
extern Uintlist_T
Uintlistpool_push (Uintlist_T list, T this, unsigned int integer);
extern Uintlist_T
Uintlistpool_pop (Uintlist_T list, unsigned int *integer);
extern Uintlist_T
Uintlistpool_copy (Uintlist_T source, T this);
extern Uintlist_T
Uintlistpool_copy_but_last (Uintlist_T source, T this);

#undef T
#endif
