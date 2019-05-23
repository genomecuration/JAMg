/* $Id: localdb-write.h 212659 2018-01-20 00:58:14Z twu $ */
#ifndef LOCALDB_WRITE_INCLUDED
#define LOCALDB_WRITE_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include "bool.h"
#include "types.h"		/* For Oligospace_T */
#include "iit-read-univ.h"

#ifdef PMAP
#include "alphabet.h"
#endif


void
Localdb_write (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
	       Alphabet_T alphabet, Width_T local1part_aa, bool watsonp,
#else
	       Width_T local1part,
#endif
	       Width_T local1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p);


#endif
