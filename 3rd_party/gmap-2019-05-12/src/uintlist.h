/* $Id: uintlist.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef UINTLIST_INCLUDED
#define UINTLIST_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_INLINE */
#endif

typedef struct Uintlist_T *Uintlist_T;

#include <stdlib.h>
#include "bool.h"
#include "mem.h"

#define T Uintlist_T
struct T {
  unsigned int first;
  struct T *rest;
};


#if !defined(HAVE_INLINE)

extern T Uintlist_push (T list, unsigned int x);
extern T Uintlist_pop (T list, unsigned int *x);
extern unsigned int Uintlist_head (T list);
extern T Uintlist_next (T list);
extern void Uintlist_free (T *list);
extern T Uintlist_reverse (T list);
extern int Uintlist_length (T list);

#else
/* extern inline T Uintlist_push (T list, unsigned int x); */
static inline T
Uintlist_push (T list, unsigned int x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

/* extern inline T Uintlist_pop (T list, unsigned int *x); */
static inline T
Uintlist_pop (T list, unsigned int *x) {
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
  
/* extern inline unsigned int Uintlist_head (T list); */
static inline unsigned int
Uintlist_head (T list) {
  return list->first;
}

/* extern inline T Uintlist_next (T list); */
static inline T
Uintlist_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

/* extern inline void Uintlist_free (T *list); */
static inline void
Uintlist_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

/* extern inline T Uintlist_reverse (T list); */
static inline T
Uintlist_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

/* extern inline int Uintlist_length (T list); */
static inline int
Uintlist_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}
#endif


extern void 
Uintlist_head_set (T list, unsigned int x);
extern T
Uintlist_keep_one (T list, int i);
extern unsigned int
Uintlist_max (T list);
extern unsigned int
Uintlist_min (T list);
extern unsigned int *
Uintlist_to_array (int *n, T list);
extern void
Uintlist_fill_array (unsigned int *array, T list);
extern void
Uintlist_fill_array_and_free (unsigned int *array, T *list);
extern unsigned int *
Uintlist_to_array_out (int *n, T list);
extern T
Uintlist_from_array (unsigned int *array, int n);
extern T 
Uintlist_copy (T list);
extern T 
Uintlist_append (T list, T tail);
extern unsigned int 
Uintlist_last_value (T this);
extern unsigned int 
Uintlist_index (T this, int index);
extern bool
Uintlist_find (T this, unsigned int value);
extern char *
Uintlist_to_string (T this);
extern char *
Uintlist_to_string_offset (T this, unsigned int chroffset);

#undef T
#endif
