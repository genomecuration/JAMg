static char rcsid[] = "$Id: univdiagpool.c 219221 2019-05-12 22:28:21Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "univdiagpool.h"
#include "univdiagdef.h"

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "comp.h"


#define CHUNKSIZE 20000

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

/* For popping */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#define T Univdiagpool_T
struct T {
  int nunivdiags;
  int univdiagctr;
  struct Univdiag_T *univdiagptr;
  List_T univdiagchunks;

  int nlistcells;
  int listcellctr;
  struct List_T *listcellptr;
  List_T listcellchunks;
};

void
Univdiagpool_free_memory (T this) {
  List_T p;
  struct Univdiag_T *univdiagptr;
  struct List_T *listcellptr;

  for (p = this->univdiagchunks; p != NULL; p = List_next(p)) {
    univdiagptr = (struct Univdiag_T *) List_head(p);
    FREE_KEEP(univdiagptr);
  }
  List_free_keep(&this->univdiagchunks);
  for (p = this->listcellchunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct List_T *) List_head(p);
    FREE_KEEP(listcellptr);
  }
  List_free_keep(&this->listcellchunks);

  this->nunivdiags = 0;
  this->univdiagctr = 0;
  this->univdiagchunks = NULL;
  /* this->univdiagptr = add_new_univdiagchunk(this); */

  this->nlistcells = 0;
  this->listcellctr = 0;
  this->listcellchunks = NULL;
  /* this->listcellptr = add_new_listcellchunk(this); */

  return;
}

void
Univdiagpool_free (T *old) {
  Univdiagpool_free_memory(*old);
  FREE_KEEP(*old);
  return;
}


void
Univdiagpool_report_memory (T this) {
  printf("Univdiagpool has %d univdiagchunks and %d listcellchunks\n",
	 List_length(this->univdiagchunks),List_length(this->listcellchunks));
  return;
}

static struct Univdiag_T *
add_new_univdiagchunk (T this) {
  struct Univdiag_T *chunk;

  chunk = (struct Univdiag_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Univdiag_T));
  this->univdiagchunks = List_push_keep(this->univdiagchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of univdiags.  Ptr for univdiag %d is %p\n",
		this->nunivdiags,chunk));

  this->nunivdiags += CHUNKSIZE;

  return chunk;
}

static struct List_T *
add_new_listcellchunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct List_T));
  this->listcellchunks = List_push_keep(this->listcellchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  this->nlistcells += CHUNKSIZE;

  return chunk;
}

T
Univdiagpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->nunivdiags = 0;
  new->univdiagctr = 0;
  new->univdiagchunks = NULL;
  /* new->univdiagptr = add_new_univdiagchunk(new); */

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->listcellchunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Univdiagpool_reset (T this) {
  this->univdiagctr = 0;
  this->listcellctr = 0;
  return;
}

Univdiag_T
Univdiag_new (T this, int qstart, int qend, Univcoord_T univdiagonal) {
  Univdiag_T univdiag;
  List_T p;
  int n;

  if (this->univdiagctr >= this->nunivdiags) {
    this->univdiagptr = add_new_univdiagchunk(this);
  } else if ((this->univdiagctr % CHUNKSIZE) == 0) {
    for (n = this->nunivdiags - CHUNKSIZE, p = this->univdiagchunks;
	 n > this->univdiagctr; p = p->rest, n -= CHUNKSIZE) ;
    this->univdiagptr = (struct Univdiag_T *) p->first;
    debug1(printf("Located univdiag %d at %p\n",this->univdiagctr,this->univdiagptr));
  }    
  univdiag = this->univdiagptr++;
  this->univdiagctr++;

  univdiag->univdiagonal = univdiagonal;
  univdiag->qstart = qstart;
  univdiag->qend = qend;
  univdiag->intscore = -1;	/* Unknown number of matches */

  debug(printf("Creating %p: %u %d..%d\n",univdiag,univdiag->univdiagonal,univdiag->querystart,univdiag->queryend));
  assert(qstart < qend);

  return univdiag;
}


List_T
Univdiagpool_push (List_T list, T this, int qstart, int qend, Univcoord_T univdiagonal) {
  List_T listcell;
  Univdiag_T univdiag;
  List_T p;
  int n;

  if (this->univdiagctr >= this->nunivdiags) {
    this->univdiagptr = add_new_univdiagchunk(this);
  } else if ((this->univdiagctr % CHUNKSIZE) == 0) {
    for (n = this->nunivdiags - CHUNKSIZE, p = this->univdiagchunks;
	 n > this->univdiagctr; p = p->rest, n -= CHUNKSIZE) ;
    this->univdiagptr = (struct Univdiag_T *) p->first;
    debug1(printf("Located univdiag %d at %p\n",this->univdiagctr,this->univdiagptr));
  }    
  univdiag = this->univdiagptr++;
  this->univdiagctr++;

  univdiag->univdiagonal = univdiagonal;
  univdiag->qstart = qstart;
  univdiag->qend = qend;
  univdiag->intscore = -1;	/* Unknown number of matches */

  debug(printf("Creating %p: %u %d..%d\n",univdiag,univdiag->univdiagonal,univdiag->querystart,univdiag->queryend));
  assert(qstart < qend);

  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) univdiag;
  listcell->rest = list;

  return listcell;
}


List_T
Univdiagpool_pop (List_T list, Univdiag_T *x) {
  List_T head;

  if (list != NULL) {
    head = list->rest;
    *x = (Univdiag_T) list->first;
    return head;
  } else {
    return list;
  }
}


List_T
Univdiagpool_push_existing (List_T list, T this, Univdiag_T univdiag) {
  List_T listcell;
  List_T p;
  int n;

  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) univdiag;
  listcell->rest = list;

  return listcell;
}
