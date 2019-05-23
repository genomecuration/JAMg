/* $Id: uint8table_rh.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef UINT8TABLE_RH_INCLUDED
#define UINT8TABLE_RH_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types.h"		/* For HAVE_64_BIT and UINT8 */

#ifdef HAVE_64_BIT
typedef struct Uint8table_T *Uint8table_T;
#endif

#include "bool.h"

#ifdef HAVE_64_BIT
#define T Uint8table_T

extern T 
Uint8table_new (int n, bool save_contents_p);

extern int
Uint8table_length (T this);

extern void *
Uint8table_get (T this, const UINT8 key);

extern void
Uint8table_put (T this, UINT8 key, void *contents);

extern void
Uint8table_put_and_save (T this, UINT8 key, void *contents);

extern void **
Uint8table_saved_values (int *nvalues, T this);

extern UINT8 *
Uint8table_keys (T this, bool sortp);

extern void 
Uint8table_free (T *old);

#undef T

#endif /*HAVE_64_BIT */

#endif


