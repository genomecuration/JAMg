/* $Id: splicetrie.h 218286 2019-01-23 16:46:55Z twu $ */
#ifndef SPLICETRIE_INCLUDED
#define SPLICETRIE_INCLUDED

#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "list.h"
#include "intlist.h"
#include "splicetrie_build.h"	/* For Splicetype_T */
#include "compress.h"

#ifndef GSNAP
#include "dynprog.h"
#endif
#include "pairpool.h"


extern void
Splicetrie_setup (Univcoord_T *splicesites_in, Genomecomp_T *splicefrags_ref_in, Genomecomp_T *splicefrags_alt_in,
		  Trieoffset_T *trieoffsets_obs_in, Triecontent_T *triecontents_obs_in,
		  Trieoffset_T *trieoffsets_max_in, Triecontent_T *triecontents_max_in,
		  Genome_T genomebits_in, Genome_T genomebits_alt_in,
		  bool snpp_in, bool amb_closest_p_in, bool amb_clip_p_in, int min_shortend_in);

#ifdef GSNAP
extern bool
Splicetrie_splicesite_p (Univcoord_T left, int pos5, int pos3);
#endif

#ifndef GSNAP
extern List_T
Splicetrie_solve_end5 (List_T best_pairs, Triecontent_T *triecontents, Trieoffset_T *trieoffsets, int j,
		       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		       int splicelength, int contlength, Splicetype_T far_splicetype,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *revsequence1, char *revsequenceuc1,
		       int length1, int length2, int revoffset1, int revoffset2,
		       int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		       Pairpool_T pairpool, int extraband_end, double defect_rate);
#endif

#ifndef GSNAP
extern List_T
Splicetrie_solve_end3 (List_T best_pairs, Triecontent_T *triecontents, Trieoffset_T *trieoffsets, int j,
		       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		       int splicelength, int contlength, Splicetype_T far_splicetype,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *sequence1, char *sequenceuc1,
		       int length1, int length2, int offset1, int offset2,
		       int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		       Pairpool_T pairpool, int extraband_end, double defect_rate);
#endif

#ifdef GSNAP
#ifdef END_KNOWNSPLICING_SHORTCUT

extern Univcoordlist_T
Splicetrie_dump_coords_left (int *best_nmismatches, Triecontent_T *triestart, int pos5, int pos3,
			     Compress_T query_compress, char *queryptr, bool plusp,
			     Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high);

extern Univcoordlist_T
Splicetrie_dump_coords_right (int *best_nmismatches, Triecontent_T *triestart, int pos5, int pos3,
			      Compress_T query_compress, char *queryptr, bool plusp,
			      Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high);
#endif
#endif


#endif
