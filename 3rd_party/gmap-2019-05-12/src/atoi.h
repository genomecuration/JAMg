/* $Id: atoi.h 214305 2018-03-19 23:40:43Z twu $ */
#ifndef ATOI_INCLUDED
#define ATOI_INCLUDED

#include "types.h"

extern Oligospace_T
Atoi_reduce_ag (Oligospace_T oligo);
extern Oligospace_T
Atoi_reduce_tc (Oligospace_T oligo);
extern UINT2
Atoi_reduce_tc_local (UINT2 oligo);
extern UINT2
Atoi_reduce_ag_local (UINT2 oligo);

#endif

