/* $Id: cmet.h 214305 2018-03-19 23:40:43Z twu $ */
#ifndef CMET_INCLUDED
#define CMET_INCLUDED

#include "types.h"

extern Oligospace_T
Cmet_reduce_ct (Oligospace_T oligo);
extern Oligospace_T
Cmet_reduce_ga (Oligospace_T oligo);
extern UINT2
Cmet_reduce_ct_local (UINT2 oligo);
extern UINT2
Cmet_reduce_ga_local (UINT2 oligo);

#endif

