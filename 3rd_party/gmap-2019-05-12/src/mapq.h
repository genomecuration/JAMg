/* $Id: mapq.h 209672 2017-09-05 22:37:11Z twu $ */
#ifndef MAPQ_INCLUDED
#define MAPQ_INCLUDED

#include "types.h"
#include "compress.h"
#include "genomicpos.h"
#include "genome.h"


#define MAX_QUALITY_SCORE_INPUT 96	/* Was 40 */
#define MAX_QUALITY_SCORE 40

extern void
MAPQ_init (int quality_score_adj_in);
extern int
MAPQ_max_quality_score (char *quality_string, int querylength);
extern float
MAPQ_loglik (Genome_T ome, Genome_T ome_alt, Compress_T query_compress,
	     Univcoord_T left, int querystart, int queryend,
	     int querylength, char *quality_string, bool plusp, int genestrand);

#endif

