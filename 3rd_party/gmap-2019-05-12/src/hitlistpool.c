static char rcsid[] = "$Id: hitlistpool.c 218195 2019-01-17 13:53:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "hitlistpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"


/* Same as Listpool_T.  Currently used for lists of Stage3end_T and
   Stage3pair_T objects */


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


#define T Hitlistpool_T
struct T {
  int nlistcells;
  int listcellctr;
  struct List_T *listcellptr;
  List_T chunks;
};

void
Hitlistpool_free_memory (T this) {
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
Hitlistpool_free (T *old) {
  Hitlistpool_free_memory(*old);
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
Hitlistpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->chunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Hitlistpool_reset (T this) {
  this->listcellctr = 0;
  return;
}

#ifdef DEBUG_HITLISTPOOL
List_T
Hitlist_push_actual (List_T list, T this, void *contents,
		     const char *file, int line) {
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

  printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
	 listcell,listcell,file,line);
  return listcell;
}

#else
List_T
Hitlist_push (List_T list, T this, void *contents) {
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
#endif


#ifdef DEBUG_HITLISTPOOL
void
Hitlist_free_actual (List_T *old, const char *file, int line) {
  List_T p;

  for (p = *old; p != NULL; p = List_next(p)) {
    printf("Freeing %p in standard pool at %s:%d\n",p,file,line);
  }
  *old = (List_T) NULL;
  return;
}

#elif !defined(HAVE_INLINE)
void
Hitlist_free (List_T *old) {
  *old = (List_T) NULL;
  return;
}
#endif

List_T
Hitlist_pop (List_T list, void **contents) {
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
Hitlist_copy (List_T source, T this) {
  List_T dest = NULL;

  while (source != NULL) {
    dest = Hitlist_push(dest,this,/*orig*/source->first);
    source = source->rest;
  }
  return List_reverse(dest);
}


