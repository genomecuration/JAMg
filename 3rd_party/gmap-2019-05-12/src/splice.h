/* $Id: splice.h 218675 2019-03-16 01:25:48Z twu $ */
#ifndef SPLICE_INCLUDED
#define SPLICE_INCLUDED

#include "bool.h"
#include "types.h"	/* For Splicetype_T */
#include "list.h"
#include "method.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "compress.h"
#include "genome.h"


typedef struct Spliceinfo_T *Spliceinfo_T;
struct Spliceinfo_T {
  int *int_memory;

  int *donor_positions_alloc;
  int *donor_knowni_alloc;
  int *acceptor_positions_alloc;
  int *acceptor_knowni_alloc;

  int *segmenti_donor_knownpos;
  int *segmentj_acceptor_knownpos;
  int *segmentj_antidonor_knownpos;
  int *segmenti_antiacceptor_knownpos;

  int *segmenti_donor_knowni;
  int *segmentj_acceptor_knowni;
  int *segmentj_antidonor_knowni;
  int *segmenti_antiacceptor_knowni;

  int segmenti_donor_nknown;
  int segmentj_acceptor_nknown;
  int segmentj_antidonor_nknown;
  int segmenti_antiacceptor_nknown;

  double *double_memory;

  double *donor_probs_alloc;
  double *acceptor_probs_alloc;
};


extern void
Splice_setup (Genome_T genomebits_in, Genome_T genomebits_alt_in, int min_shortend_in);

extern void
Spliceinfo_free (Spliceinfo_T *old);

extern Spliceinfo_T
Spliceinfo_new (int querylength);

extern int
Splice_resolve_sense (int *best_nindels, int *best_indel_pos, int *best_knowni_i, int *best_knowni_j,
		      int *best_nmismatches_i, int *best_nmismatches_j, int *best_nmismatches_indel,
		      double *best_prob_i, double *best_prob_j,

		      Univcoord_T segmenti_left, Univcoord_T segmentj_left,
		      Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
		     
		      int querystart, int queryend, int querylength, Compress_T query_compress,
		      Spliceinfo_T spliceinfo, int max_mismatches_allowed,
		      bool plusp, int genestrand, int max_deletionlen, int max_insertionlen,
		      bool allow_indel_p);

extern int
Splice_resolve_antisense (int *best_nindels, int *best_indel_pos, int *best_knowni_i, int *best_knowni_j,
			  int *best_nmismatches_i, int *best_nmismatches_j, int *best_nmismatches_indel,
			  double *best_prob_i, double *best_prob_j,
			  
			  Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			  Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
			  
			  int querystart, int queryend, int querylength, Compress_T query_compress,
			  Spliceinfo_T spliceinfo, int max_mismatches_allowed,
			  bool plusp, int genestrand, int max_deletionlen, int max_insertionlen,
			  bool allow_indel_p);

extern int
Splice_resolve_distant (int *best_nmismatches_i, int *best_nmismatches_j,
			double *best_prob_i, double *best_prob_j, int *sensedir_distant_guess,

			int *mismatch_positions_i, int segmenti_nmismatches,
			int *mismatch_positions_j, int segmentj_nmismatches,

			Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,

			int querystart, int queryend,
			int splice_pos_start, int splice_pos_end, int querylength,
			bool plusp_i, bool plusp_j);

#if 0
extern List_T
Splice_solve_single_sense (int *found_score, int *nhits, List_T hits, List_T *lowprob,

			   bool *segmenti_usedp, bool *segmentj_usedp,
			   Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			   Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
			   Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
			   Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
			   Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
			   int querylength, Compress_T query_compress,
			   Spliceinfo_T spliceinfo, nt max_mismatches_allowed,
			   bool plusp, int genestrand, bool first_read_p,
			   bool subs_or_indels_p, Method_T method, int level);
#endif

#if 0
extern List_T
Splice_solve_single_antisense (int *found_score, int *nhits, List_T hits, List_T *lowprob,

			       bool *segmenti_usedp, bool *segmentj_usedp,
			       Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			       Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
			       Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
			       Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
			       Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
			       int querylength, Compress_T query_compress,
			       Spliceinfo_T spliceinfo, int max_mismatches_allowed,
			       bool plusp, int genestrand, bool first_read_p,
			       bool subs_or_indels_p, Method_T method, int level);
#endif

extern int
Splice_trim_novel_spliceends_5 (Splicetype_T *splicetype, int **ambig_qends, double **ambig_probs_3,
				Univcoord_T left, int qstart, int qend,
				int *mismatch_positions, int nmismatches,
				Univcoord_T chroffset, bool plusp, int sensedir);

extern int
Splice_trim_novel_spliceends_3 (Splicetype_T *splicetype, int **ambig_qends, double **ambig_probs_3,
				Univcoord_T left, int qstart, int qend, int querylength,
				int *mismatch_positions, int nmismatches,
				Univcoord_T chroffset, bool plusp, int sensedir);

#endif

