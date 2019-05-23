/* $Id: indel.h 212659 2018-01-20 00:58:14Z twu $ */
#ifndef INDEL_INCLUDED
#define INDEL_INCLUDED

#include "bool.h"
#include "list.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "compress.h"
#include "genome.h"

extern void
Indel_setup (int min_indel_end_matches_in, int indel_penalty_middle_in);

extern int
Indel_resolve_middle_insertion (int *best_nmismatches_i, int *best_nmismatches_j,
				Univcoord_T left, int indels,
				int *mismatch_positions_left, int nmismatches_left,
				int *mismatch_positions_right, int nmismatches_right,
				Genome_T ome, Genome_T ome_alt, Compress_T query_compress,
				int querystart, int queryend, int querylength,
				int max_mismatches_allowed, bool plusp, int genestrand,
				bool want_lowest_coordinate_p);

int
Indel_resolve_middle_deletion (int *best_nmismatches_i, int *best_nmismatches_j,
			       Univcoord_T left, int indels,
			       int *mismatch_positions_left, int nmismatches_left,
			       int *mismatch_positions_right, int nmismatches_right,
			       Genome_T ome, Genome_T ome_alt, Compress_T query_compress,
			       int querystart, int queryend, int querylength,
			       int max_mismatches_allowed, bool plusp, int genestrand,
			       bool want_lowest_coordinate_p);

extern int
Indel_resolve_middle_deletion_or_splice (int *best_introntype, int *best_nmismatches_i, int *best_nmismatches_j,
					 Univcoord_T left, int indels,
					 Genome_T ome, Genome_T ome_alt, Compress_T query_compress,
					 int querystart, int queryend, int querylength,
					 int max_mismatches_allowed, bool plusp, int genestrand,
					 int min_intronlength, bool want_lowest_coordinate_p);


extern int
Indel_solve_end_low (int *best_adj, int *total_nmismatches, Univcoord_T left, int firstbound,
		     int querylength, Compress_T query_compress,
		     Genome_T omebits, Genome_T omebits_alt,
		     int max_end_insertions, int max_end_deletions,
		     int nmismatches_allowed, bool plusp, int genestrand,
		     bool want_lowest_coordinate_p);

extern int
Indel_solve_end_high (int *best_adj, int *total_nmismatches, Univcoord_T left, int lastbound,
		      int querylength, Compress_T query_compress,
		      Genome_T omebits, Genome_T omebits_alt,
		      int max_end_insertions, int max_end_deletions,
		      int nmismatches_allowed, bool plusp, int genestrand,
		      bool want_lowest_coordinate_p);

#endif

