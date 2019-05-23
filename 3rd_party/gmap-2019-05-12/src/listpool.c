static char rcsid[] = "$Id: listpool.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "listpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"


/* Used for lists of Substring_T and Junction_T objects in Stage3end_T objects */
/* Used for lists of Path_T, Uinttable_T, and Univdiag_T objects in path-solve.c */
/* Used for lists of Elt_T objects in extension-search.c */


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


#define T Listpool_T
struct T {
  int nlistcells;
  int listcellctr;
  struct List_T *listcellptr;
  List_T chunks;
};

void
Listpool_free_memory (T this) {
  List_T p;
  struct List_T *listcellptr;

  for (p = this->chunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct List_T *) List_head(p);
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
Listpool_free (T *old) {
  Listpool_free_memory(*old);
  FREE_KEEP(*old);
  return;
}



static struct List_T *
add_new_chunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct List_T));
  this->chunks = List_push_keep(this->chunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  this->nlistcells += CHUNKSIZE;
  return chunk;
}

T
Listpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->chunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Listpool_reset (T this) {
  this->listcellctr = 0;
  return;
}

List_T
Listpool_push (List_T list, T this, void *contents) {
  List_T listcell;
  List_T p;
  int n;

  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_chunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->chunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = contents;
  listcell->rest = list;

  return listcell;
}

List_T
Listpool_pop (List_T list, void **contents) {
  List_T head;

  if (list != NULL) {
    head = list->rest;
    *contents = list->first;
    return head;
  } else {
    return list;
  }
}


List_T
Listpool_copy (List_T source, T this) {
  List_T dest = NULL;

  while (source != NULL) {
    dest = Listpool_push(dest,this,/*orig*/source->first);
    source = source->rest;
  }
  return List_reverse(dest);
}

