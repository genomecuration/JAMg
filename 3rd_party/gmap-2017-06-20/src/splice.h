/* $Id: splice.h 197917 2016-09-16 13:39:50Z twu $ */
#ifndef SPLICE_INCLUDED
#define SPLICE_INCLUDED
#include "bool.h"
#include "list.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "compress.h"
#include "splicetrie_build.h"	/* For Splicetype_T */

extern void
Splice_setup (int min_shortend_in);

extern int
Splice_resolve_sense (int *best_knowni_i, int *best_knowni_j,
		      int *best_nmismatches_i, int *best_nmismatches_j,
		      double *best_prob_i, double *best_prob_j,

		      Univcoord_T segmenti_left, Univcoord_T segmentj_left,
		      Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
		     
		      int querystart, int queryend, int querylength, Compress_T query_compress,
		      int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
		      int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		      int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
		      int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		      int segmenti_donor_nknown, int segmentj_acceptor_nknown,
		      int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
		      int max_mismatches_allowed, bool plusp, int genestrand);

extern int
Splice_resolve_antisense (int *best_knowni_i, int *best_knowni_j,
			  int *best_nmismatches_i, int *best_nmismatches_j,
			  double *best_prob_i, double *best_prob_j,
			  
			  Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			  Univcoord_T segmenti_chroffset, Univcoord_T segmentj_chroffset,
			  
			  int querystart, int queryend, int querylength, Compress_T query_compress,
			  int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			  int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			  int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			  int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			  int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			  int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			  int max_mismatches_allowed, bool plusp, int genestrand);

extern List_T
Splice_solve_single_sense (int *found_score, int *nhits, List_T hits, List_T *lowprob,

			   bool *segmenti_usedp, bool *segmentj_usedp,
			   Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			   Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
			   Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
			   Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
			   Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
			   int querylength, Compress_T query_compress,
			   int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			   int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			   int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			   int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			   int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			   int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			   int splicing_penalty, int max_mismatches_allowed,
			   bool plusp, int genestrand, bool first_read_p,
			   bool subs_or_indels_p, bool sarrayp);

extern List_T
Splice_solve_single_antisense (int *found_score, int *nhits, List_T hits, List_T *lowprob,

			       bool *segmenti_usedp, bool *segmentj_usedp,
			       Univcoord_T segmenti_left, Univcoord_T segmentj_left,
			       Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
			       Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
			       Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
			       Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
			       int querylength, Compress_T query_compress,
			       int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
			       int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
			       int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
			       int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
			       int segmenti_donor_nknown, int segmentj_acceptor_nknown,
			       int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
			       int splicing_penalty, int max_mismatches_allowed,
			       bool plusp, int genestrand, bool first_read_p,
			       bool subs_or_indels_p, bool sarrayp);

extern List_T
Splice_group_by_segmenti (int *found_score, List_T localsplicing, List_T *ambiguous, 
			  int querylength, bool first_read_p, bool sarrayp);

extern List_T
Splice_group_by_segmentj (int *found_score, List_T localsplicing, List_T *ambiguous, 
			  int querylength, bool first_read_p, bool sarrayp);

extern int
Splice_trim_novel_spliceends (int *ambig_end_length_5, int *ambig_end_length_3,
			      Splicetype_T *ambig_splicetype_5, Splicetype_T *ambig_splicetype_3,
			      double *ambig_prob_5, double *ambig_prob_3, int orig_sensedir,
			      Univcoord_T start5, Univcoord_T middle5, Univcoord_T end5, bool solve5p,
			      Univcoord_T start3, Univcoord_T middle3, Univcoord_T end3, bool solve3p,
			      Univcoord_T genomicstart5, Univcoord_T genomicend3,
			      Univcoord_T chroffset, bool plusp);

#endif

