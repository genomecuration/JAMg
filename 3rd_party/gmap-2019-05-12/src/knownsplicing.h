/* $Id: knownsplicing.h 218286 2019-01-23 16:46:55Z twu $ */
#ifndef KNOWNSPLICING_INCLUDED
#define KNOWNSPLICING_INCLUDED


#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "genome.h"
#include "univcoord.h"


extern void
Knownsplicing_setup (Genomecomp_T *splicecomp_in);

extern bool
Knownsplicing_splicesite_p (Univcoord_T left, int pos5, int pos3);

extern Univcoord_T *
Knownsplicing_retrieve_via_splicesites (bool *distances_observed_p, Genomecomp_T **splicecomp,
					Splicetype_T **splicetypes, Chrpos_T **splicedists,
					int *nsplicesites, IIT_T splicing_iit, int *splicing_divint_crosstable,
					int donor_typeint, int acceptor_typeint, Univ_IIT_T chromosome_iit,
					Genome_T genome, Genome_T genomealt, Chrpos_T shortsplicedist);

#endif
