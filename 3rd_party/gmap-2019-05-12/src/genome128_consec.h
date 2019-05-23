/* $Id: genome128_consec.h 209672 2017-09-05 22:37:11Z twu $ */
#ifndef GENOME128_CONSEC_INCLUDED
#define GENOME128_CONSEC_INCLUDED
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "compress.h"
#include "genome.h"

extern void
Genome_consec_setup (bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		     Mode_T mode_in);

extern int
Genome_consecutive_matches_rightward (Genome_T ome, Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				      bool plusp, int genestrand);
extern int
Genome_consecutive_matches_leftward (Genome_T ome, Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				     bool plusp, int genestrand);
extern int
Genome_consecutive_matches_pair (Genome_T ome, UINT4 lefta, UINT4 leftb, UINT4 genomelength);

#endif

