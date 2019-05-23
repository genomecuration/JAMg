/* $Id: hitlistpool.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef HITLISTPOOL_INCLUDED
#define HITLISTPOOL_INCLUDED

#include "list.h"

#define T Hitlistpool_T
typedef struct T *T;

/* #define DEBUG_HITLISTPOOL 1 */

#ifdef DEBUG_HITLISTPOOL

extern List_T
Hitlist_push_actual (List_T list, T this, void *contents,
		     const char *file, int line);
#define Hitlist_push(list,this,contents) Hitlist_push_actual(list,this,contents,__FILE__,__LINE__)


extern void
Hitlist_free_actual (List_T *old, const char *file, int line);
#define Hitlist_free(old) Hitlist_free_actual(old,__FILE__,__LINE__)


#elif !defined(HAVE_INLINE)

extern List_T
Hitlist_push (List_T list, T this, void *contents);

extern void Hitlist_free (List_T *old);

#else

extern List_T
Hitlist_push (List_T list, T this, void *contents);

static inline void
Hitlist_free (List_T *old) {
  *old = (List_T) NULL;
  return;
}

#endif


extern void
Hitlistpool_free (T *old);
extern void
Hitlistpool_free_memory (T this);
extern T
Hitlistpool_new (void);
extern void
Hitlistpool_reset (T this);
extern List_T
Hitlist_pop (List_T list, void **contents);
extern List_T
Hitlist_copy (List_T source, T this);

#undef T
#endif
