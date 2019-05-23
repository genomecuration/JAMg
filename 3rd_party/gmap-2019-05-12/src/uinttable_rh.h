/* $Id: uinttable_rh.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef UINTTABLE_RH_INCLUDED
#define UINTTABLE_RH_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bool.h"

#define T Uinttable_T
typedef struct T *T;

extern T 
Uinttable_new (int n, bool save_contents_p);

extern int
Uinttable_length (T this);

extern void *
Uinttable_get (T this, const unsigned int key);

extern void
Uinttable_put (T this, unsigned int key, void *contents);

extern void
Uinttable_put_and_save (T this, unsigned int key, void *contents);

extern void **
Uinttable_saved_values (int *nvalues, T this);

unsigned int *
Uinttable_keys (T this, bool sortp);

extern void 
Uinttable_free (T *old);

#undef T
#endif


