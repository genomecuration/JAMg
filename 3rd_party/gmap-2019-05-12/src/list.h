/* $Id: list.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef LIST_INCLUDED
#define LIST_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_INLINE */
#endif

typedef struct List_T *List_T;

#include <stdlib.h>
#include "mem.h"

#define T List_T
struct T {
  void *first;
  struct T *rest;
};


#if !defined(HAVE_INLINE)

extern T List_push (T list, void *x);
extern T List_pop (T list, void **x);
extern void *List_head (T list);
extern T List_next (T list);
extern void List_free (T *list);
extern int List_length (T list);
extern T List_reverse (T list);
extern T List_transfer_one (T dest, T *source);

#else
/* extern inline T List_push (T list, void *x); */

#if 0
/* Code to trace callers of List_push */
static inline T
List_push_actual (T list, void *x, const char *file, int line) {
  T new = (T) MALLOC(sizeof(*new));
  
  printf("List_push %s:%d\n",file,line);
  new->first = x;
  new->rest = list;
  return new;
}

#define List_push(list,x) List_push_actual(list,x,__FILE__,__LINE__)
#else
static inline T
List_push (T list, void *x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}
#endif


/* extern inline T List_pop (T list, void **x); */
static inline T
List_pop (T list, void **x) {
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

/* extern inline void *List_head (T list); */
static inline void *
List_head (T list) {
  return list->first;
}

/* extern inline T List_next (T list); */
static inline T
List_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

/* extern inline void List_free (T *list); */
static inline void
List_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

/* extern inline int List_length (T list); */
static inline int
List_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }

  return n;
}

/* extern inline T List_reverse (T list);*/
static inline T
List_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

/* extern inline T List_transfer_one (T dest, T *source); */
static inline T
List_transfer_one (T dest, T *source) {
  T next;

  next = (*source)->rest;
  (*source)->rest = dest;
  dest = *source;
  *source = next;
  return dest;
}
#endif



extern T List_push_keep (T list, void *x);
extern T List_push_out (T list, void *x);
extern T List_pop_out (T list, void **x);
extern void List_head_set (T list, void *x);
extern void List_tail_set (T this, T rest);
extern void List_free_keep (T *list);
extern void List_free_out (T *list);
extern T List_truncate (T list, int n);
extern void **List_to_array (T list, void *end);
extern void List_fill_array (void **array, T list);
extern void List_fill_array_and_free (void **array, T *list);
extern T List_fill_array_with_handle (struct T *new, void **array, int nelts);
extern void **List_to_array_out (T list, void *end);
extern void **List_to_array_n (int *n, T list);
extern void **List_to_array_out_n (int *n, T list);
extern T List_copy (T list);
extern void
List_dump (T list);
extern T List_append (T list, T tail);
extern void *
List_last_value (T this);
extern T
List_last_item (T this);
extern void *
List_index (T this, int index);
extern T
List_insert (T this, void *x);
extern void
List_reinsert (T *listptr, T cell);
extern T
List_push_existing (T dest, T source);
extern T
List_from_string (char *string);

#undef T
#endif
