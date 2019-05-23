static char rcsid[] = "$Id: uintlistpool.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "uintlistpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "list.h"

#define CHUNKSIZE 16384

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* For mechanics of memory allocation and deallocation */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#define T Uintlistpool_T
struct T {
  int nlistcells;
  int listcellctr;
  struct Uintlist_T *listcellptr;
  List_T chunks;
};

void
Uintlistpool_free_memory (T this) {
  List_T p;
  struct Uintlist_T *listcellptr;

  for (p = this->chunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct Uintlist_T *) List_head(p);
    FREE_KEEP(listcellptr);
  }
  List_free_keep(&this->chunks);

  this->nlistcells = 0;
  this->listcellctr = 0;
  this->chunks = NULL;
  /* this->listcellptr = add_new_listcellchunk(this); */

  return;
}

void
Uintlistpool_free (T *old) {
  Uintlistpool_free_memory(*old);
  FREE_KEEP(*old);
  return;
}




static struct Uintlist_T *
add_new_chunk (T this) {
  struct Uintlist_T *chunk;

  chunk = (struct Uintlist_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Uintlist_T));
  this->chunks = List_push_keep(this->chunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  this->nlistcells += CHUNKSIZE;
  return chunk;
}

T
Uintlistpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->chunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Uintlistpool_reset (T this) {
  this->listcellctr = 0;
  return;
}

Uintlist_T
Uintlistpool_push (Uintlist_T list, T this, unsigned int integer) {
  Uintlist_T listcell;
  List_T p;
  int n;

  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_chunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->chunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct Uintlist_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = integer;
  listcell->rest = list;

  return listcell;
}

Uintlist_T
Uintlistpool_pop (Uintlist_T list, unsigned int *integer) {
  Uintlist_T head;

  if (list != NULL) {
    head = list->rest;
    *integer = list->first;
    return head;
  } else {
    return list;
  }
}


Uintlist_T
Uintlistpool_copy (Uintlist_T source, T this) {
  Uintlist_T dest = NULL;

  while (source != NULL) {
    dest = Uintlistpool_push(dest,this,/*orig*/source->first);
    source = source->rest;
  }
  return Uintlist_reverse(dest);
}

Uintlist_T
Uintlistpool_copy_but_last (Uintlist_T source, T this) {
  Uintlist_T dest = NULL;

  if (source == NULL) {
    return (Uintlist_T) NULL;
  } else {
    while (source->rest != NULL) {
      dest = Uintlistpool_push(dest,this,/*orig*/source->first);
      source = source->rest;
    }
    return Uintlist_reverse(dest);
  }
}

