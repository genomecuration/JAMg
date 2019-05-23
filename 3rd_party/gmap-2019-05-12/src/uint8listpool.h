/* $Id: uint8listpool.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef UINT8LISTPOOL_INCLUDED
#define UINT8LISTPOOL_INCLUDED

typedef struct Uint8listpool_T *Uint8listpool_T;

#include "types.h"
#include "uint8list.h"

#define T Uint8listpool_T

extern void
Uint8listpool_free (T *old);
extern void
Uint8listpool_free_memory (T this);
extern T
Uint8listpool_new (void);
extern void
Uint8listpool_reset (T this);
extern Uint8list_T
Uint8listpool_push (Uint8list_T list, T this, UINT8 integer);
extern Uint8list_T
Uint8listpool_pop (Uint8list_T list, UINT8 *integer);
extern Uint8list_T
Uint8listpool_copy (Uint8list_T source, T this);
extern Uint8list_T
Uint8listpool_copy_but_last (Uint8list_T source, T this);

#undef T
#endif
