static char rcsid[] = "$Id: result.c 182440 2016-01-15 22:42:45Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "result.h"
#include <stdlib.h>
#include "mem.h"
#include "diag.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Result_T
struct T {
  int id;

  bool mergedp;			/* true if two parts were merged */
  Chimera_T chimera;		/* NULL indicates not a chimera */
  List_T gregionlist;		/* For debugging of stage 1 */
  List_T diagonals;		/* For debugging of diag */
  Stage3_T *array;
  int npaths_primary;
  int npaths_altloc;
  int first_absmq;
  int second_absmq;
  Diagnostic_T diagnostic;
  Failure_T failuretype;
};


int
Result_id (T this) {
  return this->id;
}

bool
Result_mergedp (T this) {
  return this->mergedp;
}

Chimera_T
Result_chimera (T this) {
  return this->chimera;
}


Stage3_T *
Result_array (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, T this) {
  *npaths_primary = this->npaths_primary;
  *npaths_altloc = this->npaths_altloc;
  *first_absmq = this->first_absmq;
  *second_absmq = this->second_absmq;
  return this->array;
}


List_T
Result_gregionlist (T this) {
  return this->gregionlist;
}

List_T
Result_diagonals (T this) {
  return this->diagonals;
}


Diagnostic_T
Result_diagnostic (T this) {
  return this->diagnostic;
}


Failure_T
Result_failuretype (T this) {
  return this->failuretype;
}


T
Result_new (int id, bool mergedp, Chimera_T chimera, Stage3_T *array,
	    int npaths_primary, int npaths_altloc, int first_absmq, int second_absmq,
	    Diagnostic_T diagnostic, Failure_T failuretype) {
  T new = (T) MALLOC_OUT(sizeof(*new));

  new->id = id;
  new->mergedp = mergedp;
  new->chimera = chimera;
  new->gregionlist = (List_T) NULL;
  new->diagonals = (List_T) NULL;
  new->array = array;
  new->npaths_primary = npaths_primary;
  new->npaths_altloc = npaths_altloc;
  new->first_absmq = first_absmq;
  new->second_absmq = second_absmq;
  new->diagnostic = diagnostic;
  new->failuretype = failuretype;

  return new;
}


void
Result_free (T *old) {
  Chimera_T chimera;
  Stage3_T stage3;
  int i;
  List_T p;
  Gregion_T gregion;
#ifndef USE_DIAGPOOL
  Diag_T diag;
#endif

  if (*old) {
    if ((chimera = (*old)->chimera) != NULL) {
      Chimera_free(&chimera);

      stage3 = (*old)->array[0];
      debug(printf("Freeing 0 stage3 %p and pairarray %p\n",
		   stage3,Stage3_pairarray(stage3)));
      Stage3_free(&stage3);

      stage3 = (*old)->array[1];
      debug(printf("Freeing 1 stage3 %p, but not pairarray %p\n",
		   stage3,Stage3_pairarray(stage3)));
      Stage3_free(&stage3);
      
      for (i = 2; i < (*old)->npaths_primary + (*old)->npaths_altloc; i++) {
	stage3 = (*old)->array[i];
	debug(printf("Freeing %d stage3 %p and pairarray %p\n",
		     i,stage3,Stage3_pairarray(stage3)));
	Stage3_free(&stage3);
      }

      FREE_OUT((*old)->array);

    } else {
      if ((*old)->array) {
	for (i = 0; i < (*old)->npaths_primary + (*old)->npaths_altloc; i++) {
	  stage3 = (*old)->array[i];
	  Stage3_free(&stage3);
	}
	FREE_OUT((*old)->array);
      }
    }

    if ((*old)->diagnostic != NULL) {
      Diagnostic_free(&(*old)->diagnostic);
    }
    if ((*old)->gregionlist) {
      for (p = (*old)->gregionlist; p != NULL; p = List_next(p)) {
	gregion = (Gregion_T) List_head(p);
	Gregion_free(&gregion);
      }
      List_free(&((*old)->gregionlist));
    }

#ifndef USE_DIAGPOOL
    /* No need to free since memory is allocated separately */
    if ((*old)->diagonals) {
      for (p = (*old)->diagonals; p != NULL; p = List_next(p)) {
	diag = (Diag_T) List_head(p);
	Diag_free(&diag);
      }
      List_free(&((*old)->diagonals));
    }
#endif

    FREE_OUT(*old);
  }

  return;
}

