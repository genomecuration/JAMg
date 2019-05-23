static char rcsid[] = "$Id: uintlist.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "uintlist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"

#define T Uintlist_T

#if !defined(HAVE_INLINE)
T
Uintlist_push (T list, unsigned int x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
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
  
unsigned int
Uintlist_head (T list) {
  return list->first;
}

T
Uintlist_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
Uintlist_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

T
Uintlist_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
Uintlist_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}
#endif


void
Uintlist_head_set (T this, unsigned int x) {
  this->first = x;
  return;
}

T
Uintlist_keep_one (T list, int i) {
  T head;

  while (--i >= 0) {
    /* Pop */
    head = list->rest;
    FREE(list);
    list = head;
  }

  Uintlist_free(&list->rest);
  return list;
}

unsigned int
Uintlist_max (T list) {
  unsigned int m = 0;

  while (list) {
    if (list->first > m) {
      m = list->first;
    }
    list = list->rest;
  }

  return m;
}

unsigned int
Uintlist_min (T list) {
  unsigned int m;

  if (list == NULL) {
    return 0;

  } else {
    m = list->first;
    list = list->rest;
    while (list) {
      if (list->first < m) {
	m = list->first;
      }
      list = list->rest;
    }

    return m;
  }
}

unsigned int *
Uintlist_to_array (int *n, T list) {
  unsigned int *array;
  int i;

  *n = Uintlist_length(list);
  array = (unsigned int *) CALLOC(*n,sizeof(unsigned int));
  for (i = 0; i < *n; i++) {
    array[i] = list->first;
    list = list->rest;
  }
  return array;
}

void
Uintlist_fill_array (unsigned int *array, T list) {
  int i = 0;

  while (list) {
    array[i++] = list->first;
    list = list->rest;
  }

  return;
}

void
Uintlist_fill_array_and_free (unsigned int *array, T *list) {
  T prev;
  int i = 0;

  while ((prev = *list) != NULL) {
    array[i++] = prev->first;
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

unsigned int *
Uintlist_to_array_out (int *n, T list) {
  unsigned int *array;
  int i;

  *n = Uintlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (unsigned int *) CALLOC_OUT(*n,sizeof(unsigned int));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

T
Uintlist_from_array (unsigned int *array, int n) {
  T list = NULL, p;

  while (--n >= 0) {
    p = (T) MALLOC(sizeof(*p));
    p->first = array[n];
    p->rest = list;
    list = p;
  }

  return list;
}

T
Uintlist_copy (T list) {
  T head, *p = &head;

  for ( ; list; list = list->rest) {
    *p = (T) MALLOC(sizeof(**p));
    (*p)->first = list->first;
    p = &(*p)->rest;
  }
  *p = NULL;
  return head;
}

T
Uintlist_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

unsigned int
Uintlist_last_value (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

unsigned int
Uintlist_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}


bool
Uintlist_find (T this, unsigned int value) {
  T r;

  for (r = this; r != NULL; r = r->rest) {
    if (r->first == value) {
      return true;
    }
  }
  return false;
}

char *
Uintlist_to_string (T this) {
  char *string, Buffer[256];
  T p;
  int n, i, strlength;

  if ((n = Uintlist_length(this)) == 0) {
    string = (char *) CALLOC(1,sizeof(char));
    string[0] = '\0';
  } else {
    strlength = 0;
    for (i = 0, p = this; i < n-1; i++, p = Uintlist_next(p)) {
      sprintf(Buffer,"%u,",Uintlist_head(p));
      strlength += strlen(Buffer);
    }
    sprintf(Buffer,"%u",Uintlist_head(p));
    strlength += strlen(Buffer);

    string = (char *) CALLOC(strlength + 1,sizeof(char));
    string[0] = '\0';
    for (i = 0, p = this; i < n-1; i++, p = Uintlist_next(p)) {
      sprintf(Buffer,"%u,",Uintlist_head(p));
      strcat(string,Buffer);
    }
    sprintf(Buffer,"%u",Uintlist_head(p));
    strcat(string,Buffer);
  }

  return string;
}


char *
Uintlist_to_string_offset (T this, unsigned int chroffset) {
  char *string, Buffer[256];
  T p;
  int n, i, strlength;

  if ((n = Uintlist_length(this)) == 0) {
    string = (char *) CALLOC(1,sizeof(char));
    string[0] = '\0';
  } else {
    strlength = 0;
    for (i = 0, p = this; i < n-1; i++, p = Uintlist_next(p)) {
      sprintf(Buffer,"%u,",Uintlist_head(p) - chroffset);
      strlength += strlen(Buffer);
    }
    sprintf(Buffer,"%u",Uintlist_head(p));
    strlength += strlen(Buffer);

    string = (char *) CALLOC(strlength + 1,sizeof(char));
    string[0] = '\0';
    for (i = 0, p = this; i < n-1; i++, p = Uintlist_next(p)) {
      sprintf(Buffer,"%u,",Uintlist_head(p) - chroffset);
      strcat(string,Buffer);
    }
    sprintf(Buffer,"%u",Uintlist_head(p) - chroffset);
    strcat(string,Buffer);
  }

  return string;
}
