static char rcsid[] = "$Id: intlistpool.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "intlistpool.h"
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


#define T Intlistpool_T
struct T {
  int nlistcells;
  int listcellctr;
  struct Intlist_T *listcellptr;
  List_T chunks;
};

void
Intlistpool_free_memory (T this) {
  List_T p;
  struct Intlist_T *listcellptr;

  for (p = this->chunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct Intlist_T *) List_head(p);
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
Intlistpool_free (T *old) {
  Intlistpool_free_memory(*old);
  FREE_KEEP(*old);
  return;
}




static struct Intlist_T *
add_new_chunk (T this) {
  struct Intlist_T *chunk;

  chunk = (struct Intlist_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Intlist_T));
  this->chunks = List_push_keep(this->chunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  this->nlistcells += CHUNKSIZE;
  return chunk;
}

T
Intlistpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->chunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Intlistpool_reset (T this) {
  this->listcellctr = 0;
  return;
}

Intlist_T
Intlistpool_push (Intlist_T list, T this, int integer) {
  Intlist_T listcell;
  List_T p;
  int n;

  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_chunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->chunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct Intlist_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = integer;
  listcell->rest = list;

  return listcell;
}

Intlist_T
Intlistpool_pop (Intlist_T list, int *integer) {
  Intlist_T head;

  if (list != NULL) {
    head = list->rest;
    *integer = list->first;
    return head;
  } else {
    return list;
  }
}


Intlist_T
Intlistpool_copy (Intlist_T source, T this) {
  Intlist_T dest = NULL;

  while (source != NULL) {
    dest = Intlistpool_push(dest,this,/*orig*/source->first);
    source = source->rest;
  }
  return Intlist_reverse(dest);
}

Intlist_T
Intlistpool_copy_but_last (Intlist_T source, T this) {
  Intlist_T dest = NULL;

  if (source == NULL) {
    return (Intlist_T) NULL;
  } else {
    while (source->rest != NULL) {
      dest = Intlistpool_push(dest,this,/*orig*/source->first);
      source = source->rest;
    }
    return Intlist_reverse(dest);
  }
}

