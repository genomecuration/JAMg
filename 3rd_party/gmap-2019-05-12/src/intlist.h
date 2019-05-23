/* $Id: intlist.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef INTLIST_INCLUDED
#define INTLIST_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_INLINE */
#endif

typedef struct Intlist_T *Intlist_T;

#include <stdlib.h>
#include "mem.h"
#include "bool.h"

#define T Intlist_T
struct T {
  int first;
  struct T *rest;
};


#if !defined(HAVE_INLINE)

extern T Intlist_push (T list, int x);
extern T Intlist_pop (T list, int *x);
extern int Intlist_head (T list);
extern T Intlist_next (T list);
extern T Intlist_reverse (T list);
extern int Intlist_length (T list);

#else
/* extern inline T Intlist_push (T list, int x); */

#if 0
/* Code to trace callers of Intlist_push */
static inline T
Intlist_push_actual (T list, int x, const char *file, int line) {
  T new = (T) MALLOC(sizeof(*new));
  
  printf("Intlist_push %s:%d\n",file,line);
  new->first = x;
  new->rest = list;
  return new;
}

#define Intlist_push(list,x) Intlist_push_actual(list,x,__FILE__,__LINE__)
#else
static inline T
Intlist_push (T list, int x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}
#endif


/* extern inline T Intlist_pop (T list, int *x); */
static inline T
Intlist_pop (T list, int *x) {
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
  
/* extern inline int Intlist_head (T list); */
static inline int
Intlist_head (T list) {
  return list->first;
}

/* extern inline T Intlist_next (T list); */
static inline T
Intlist_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

/* extern inline void Intlist_free (T *list); */
static inline void
Intlist_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }
}

/* extern inline T Intlist_reverse (T list); */
static inline T
Intlist_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

/* extern inline int Intlist_length (T list); */
static inline int
Intlist_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}
#endif


extern T 
Intlist_push_in (T list, int x);
extern T
Intlist_insert_second (T list, int x);
extern void
Intlist_delete (T prev, T this);
extern void 
Intlist_head_set (T list, int x);
extern void 
Intlist_free_in (T *list);
extern T
Intlist_keep_one (T list, int i);
extern int
Intlist_max (T list);
extern int
Intlist_min (T list);
extern int
Intlist_sum (T list);
extern bool
Intlist_vary (T list);
extern bool
Intlist_exists_p (T list, int x);
extern int *
Intlist_to_array (int *n, T list);
extern int *
Intlist_to_array_plus_extra (int *n, T list);
extern void
Intlist_fill_array (int *array, T list);
extern void
Intlist_fill_array_and_free (int *array, T *list);
extern int *
Intlist_to_array_out (int *n, T list);
extern char *
Intlist_to_char_array (int *n, T list);
extern char *
Intlist_to_char_array_in (int *n, T list);
extern T
Intlist_from_array (int *array, int n);
extern T 
Intlist_copy (T list);
extern T 
Intlist_append (T list, T tail);
extern int 
Intlist_last_value (T this);
extern int
Intlist_penultimate_value (T this);
extern int 
Intlist_index (T this, int index);
extern T
Intlist_from_string (char *string);
extern char *
Intlist_to_string (T this);
extern int *
Intlist_array_ascending_by_key (int *n, T this, T key);
extern void
Intlist_array_dual_ascending_by_key (int *sorted, int *keyarray, int n, T this, T key);
extern T
Intlist_list_ascending_by_key (T this, T key);
extern T
Intlist_list_descending_by_key (T this, T key);
extern T
Intlist_sort_ascending (T this);
extern bool
Intlist_equal (T x, T y);
extern bool
Intlist_intersect_p (T x, T y);
extern T
Intlist_filter (T list, bool *deletep);

#undef T
#endif
