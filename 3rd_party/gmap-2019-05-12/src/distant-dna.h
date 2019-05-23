/* $Id: distant-dna.h 218675 2019-03-16 01:25:48Z twu $ */
#ifndef DISTANT_DNA_INCLUDED
#define DISTANT_DNA_INCLUDED

#include "bool.h"
#include "list.h"
#include "genomicpos.h"
#include "compress.h"
#include "genome.h"
#include "listpool.h"
#include "hitlistpool.h"

extern void
Distant_dna_solve (int *found_score_overall, int *found_score_within_trims,
		   List_T *hits_plus, List_T *hits_minus,

		   List_T startfrags_plus, List_T endfrags_plus,
		   List_T startfrags_minus, List_T endfrags_minus,

		   int *mismatch_positions_alloc,
		   Compress_T query_compress_fwd, Compress_T query_compress_rev,
		   int querylength, bool first_read_p,
		   Listpool_T listpool, Hitlistpool_T hitlistpool, int level);

extern void
Distant_dna_setup (Genome_T genomebits_in, Genome_T genomebits_alt_in,
		   Chrpos_T shortsplicedist_in);

#endif
