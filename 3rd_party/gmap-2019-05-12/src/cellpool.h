/* $Id: cellpool.h 210504 2017-10-12 18:33:46Z twu $ */
#ifndef CELLPOOL_INCLUDED
#define CELLPOOL_INCLUDED

#include "bool.h"
#include "list.h"


typedef struct Cell_T *Cell_T;
struct Cell_T {
  int rootposition;		/* Want to allow for -1 */
  int endposition;
  int querypos;
  int hit;
  bool fwdp;
  int score;
  int pushedp;
};


#define T Cellpool_T
typedef struct T *T;

extern void
Cellpool_free (T *old);
extern void
Cellpool_free_memory (T this);
extern void
Cellpool_report_memory (T this);
extern T
Cellpool_new (void);
extern void
Cellpool_reset (T this);
extern List_T
Cellpool_push (List_T list, T this, int rootposition, int endposition,
	       int querypos, int hit, bool fwdp, int score);
extern List_T
Cellpool_pop (List_T list, Cell_T *x);

#undef T
#endif


