static char rcsid[] = "$I$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "univdiag.h"
#include "univdiagdef.h"

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "assert.h"


#define T Univdiag_T


#if 0
/* Now in univdiagpool.c */
T
Univdiag_new (int qstart, int qend, Univcoord_T univdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  assert(qend > qstart);

  new->univdiagonal = univdiagonal;
  new->qstart = qstart;
  new->qend = qend;
  /* new->nconsecutive = qend - qstart + 1; */

  new->intscore = -1;		/* Unknown number of matches */

  return new;
}
#endif


#if 0
T
Univdiag_copy (T this) {
  return Univdiag_new(this->qstart,this->qend,this->univdiagonal);
}
#endif


#if 0
List_T
Univdiag_copy_list (List_T list) {
  List_T new = NULL, p;
  T this;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);

    new = List_push(new,(void *) Univdiag_new(this->qstart,this->qend,this->univdiagonal));
  }
  
  return List_reverse(new);
}
#endif
  


#if 0
T
Univdiag_new_fillin (int qstart, int qend, int indexsize, Univcoord_T univdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  new->univdiagonal = univdiagonal;
  new->qstart = qstart;
  new->qend = qend + indexsize - 1;
  /* new->nconsecutive = new->qend - qstart + 1; */

  new->intscore = -1;

  return new;
}
#endif


#if 0
/* All Univdiag_T objects and lists now allocated in univdiagpool.h */
void
Univdiag_free (T *old) {
  FREE(*old);
  return;
}
#endif

#if 0
/* All Univdiag_T objects and lists now allocated in univdiagpool.h */
void
Univdiag_gc (List_T *list) {
  T univdiagonal;
  List_T p;

  for (p = *list; p != NULL; p = List_next(p)) {
    univdiagonal = (T) List_head(p);
    FREE(univdiagonal);
  }
  List_free(&(*list));
  return;
}
#endif

#if 0
/* All Univdiag_T objects and lists now allocated in univdiagpool.h */
void
Univdiag_list_gc (List_T *paths) {
  List_T path, p;

  for (p = *paths; p != NULL; p = List_next(p)) {
    path = (List_T) List_head(p);
    Univdiag_gc(&path);
  }
  List_free(&(*paths));
  return;
}
#endif


int
Univdiag_list_length (List_T path) {
  int length = 0;
  T this;
  List_T p;

  for (p = path; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    length += this->qend - this->qstart + 1;
  }
  return length;
}


int
Univdiag_ascending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->qstart < y->qstart) {
    return -1;
  } else if (y->qstart < x->qstart) {
    return +1;
  } else if (x->qend < y->qend) {
    return -1;
  } else if (y->qend < x->qend) {
    return +1;
  } else if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal < x->univdiagonal) {
    return +1;
  } else {
    return 0;
  }
}


int
Univdiag_descending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->qstart > y->qstart) {
    return -1;
  } else if (y->qstart > x->qstart) {
    return +1;
  } else if (x->qend > y->qend) {
    return -1;
  } else if (y->qend > x->qend) {
    return +1;
  } else if (x->univdiagonal > y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal > x->univdiagonal) {
    return +1;
  } else {
    return 0;
  }
}



int
Univdiag_diagonal_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal < x->univdiagonal) {
    return +1;
  } else if (x->qstart < y->qstart) {
    return -1;
  } else if (y->qstart < x->qstart) {
    return +1;
  } else if (x->qend < y->qend) {
    return -1;
  } else if (y->qend < x->qend) {
    return +1;
  } else {
    return 0;
  }
}
