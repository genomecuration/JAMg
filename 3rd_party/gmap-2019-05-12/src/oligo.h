/* $Id: oligo.h 214025 2018-03-05 07:06:43Z twu $ */
#ifndef OLIGO_INCLUDED
#define OLIGO_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "reader.h"

typedef enum {INIT, DONE, INVALID, VALID} Oligostate_T;

extern Oligostate_T
Oligo_next_5 (Oligostate_T last_state, int *querypos, Oligospace_T *forward, 
	      Oligospace_T *revcomp, Reader_T reader, int genestrand);
extern Oligostate_T
Oligo_next_3 (Oligostate_T last_state, int *querypos, Oligospace_T *forward, 
	      Oligospace_T *revcomp, Reader_T reader, int genestrand);

extern Oligostate_T
Oligo_skip_5 (Oligostate_T last_state, int *querypos, Oligospace_T *forward,
	      Oligospace_T *revcomp, Reader_T reader, int genestrand, int nskip);
extern Oligostate_T
Oligo_skip_3 (Oligostate_T last_state, int *querypos, Oligospace_T *forward,
	      Oligospace_T *revcomp, Reader_T reader, int genestrand, int nskip);

extern void
Oligo_setup (int index1part, int mode);

#endif
