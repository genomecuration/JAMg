/* $Id: uint8list.h 213653 2018-02-23 00:13:01Z twu $ */
#ifndef UINT8LIST_INCLUDED
#define UINT8LIST_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_INLINE */
#endif

typedef struct Uint8list_T *Uint8list_T;

#include <stdlib.h>
#include "bool.h"
#include "mem.h"
#include "types.h"

#define T Uint8list_T
struct T {
  UINT8 first;
  struct T *rest;
};

#if !defined(HAVE_INLINE)

extern T Uint8list_push (T list, UINT8 x);
extern T Uint8list_pop (T list, UINT8 *x);
extern UINT8 Uint8list_head (T list);
extern T Uint8list_next (T list);
extern void Uint8list_free (T *list);
extern T Uint8list_reverse (T list);
extern int Uint8list_length (T list);

#else

static inline T
Uint8list_push (T list, UINT8 x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

static inline T
Uint8list_pop (T list, UINT8 *x) {
  T head;

  if (list) {
    head = list->rest;
    *x = list->first;
    FREE(list);
    return head;
  } else {
    return list;
  }
}

static inline UINT8
Uint8list_head (T list) {
  return list->first;
}

static inline T
Uint8list_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

static inline void
Uint8list_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }
}

static inline T
Uint8list_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

static inline int
Uint8list_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}
#endif


extern void 
Uint8list_head_set (T list, UINT8 x);
extern T
Uint8list_keep_one (T list, int i);
extern UINT8
Uint8list_max (T list);
extern UINT8
Uint8list_min (T list);
extern UINT8 *
Uint8list_to_array (int *n, T list);
extern UINT8 *
Uint8list_to_array_out (int *n, T list);
extern void
Uint8list_fill_array (UINT8 *array, T list);
extern void
Uint8list_fill_array_and_free (UINT8 *array, T *list);
extern T
Uint8list_from_array (UINT8 *array, int n);
extern T 
Uint8list_copy (T list);
extern T 
Uint8list_append (T list, T tail);
extern UINT8 
Uint8list_last_value (T this);
extern UINT8 
Uint8list_index (T this, int index);
extern bool
Uint8list_find (T this, UINT8 value);
extern char *
Uint8list_to_string (T this);
#undef T
#endif
