/* $Id: simplepair.h 218286 2019-01-23 16:46:55Z twu $ */
#ifndef SIMPLEPAIR_INCLUDED
#define SIMPLEPAIR_INCLUDED

typedef struct Simplepair_T *Simplepair_T;
typedef enum {CIGAR_ACTION_IGNORE, CIGAR_ACTION_WARNING, CIGAR_ACTION_NOPRINT, CIGAR_ACTION_ABORT} Cigar_action_T;

#include "bool.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "list.h"
#include "iit-read-univ.h"
#include "filestring.h"
#include "stage3hr.h"


#define T Simplepair_T
struct T {
  int querypos;
  Chrpos_T genomepos;

  char cdna;
  char comp;
  char genome;
  char genomealt;

  bool gapp;			/* True if comp is in a big gap (from genomic perspective):
                                   >])([<#= (but not '-' or '~'). */
};



extern Chrpos_T
Simplepair_head_genomepos (List_T pairs);
extern Chrpos_T
Simplepair_last_genomepos (List_T pairs);

extern T
Simplepair_new_out (int querypos, Chrpos_T genomepos, char cdna, char comp, char genome);
extern void
Simplepair_free_out (T *old);
extern void
Simplepair_dump_list (List_T pairs, bool zerobasedp);

extern List_T
Simplepair_strip_gaps_at_head (List_T pairs);
extern List_T
Simplepair_strip_gaps_at_tail (List_T pairs);

extern Chrpos_T
Simplepair_genomicpos_low (int hardclip_low, int hardclip_high,
			   struct T *pairarray, int npairs, int querylength,
			   bool plusp, bool hide_soft_clips_p);
extern void
Simplepair_overlap_print_sam (Filestring_T fp, char *abbrev, struct T *pairarray, int npairs,
			      char *acc1, char *acc2, Chrnum_T chrnum, Univ_IIT_T chromosome_iit,
			      char *queryseq_ptr, char *quality_string,
			      int hardclip_low, int hardclip_high, int querylength_given,
			      bool watsonp, int sensedir,
			      int quality_shift, bool first_read_p, int pathnum, int npaths_primary, int npaths_altloc,
			      int absmq_score, int second_absmq, Chrpos_T chrpos, Chrpos_T chrlength,
			      
			      Shortread_T queryseq, unsigned int flag,
			      int pair_mapq_score, int end_mapq_score,
			      Stage3pair_T stage3pair, Stage3end_T stage3hit,
			      
			      char *sam_read_group_id);

extern void
Simplepair_setup (bool novelsplicingp_in, IIT_T splicesites_iit_in,
		  Univ_IIT_T transcript_iit_in, bool sam_insert_0M_p_in,
		  bool md_lowercase_variant_p_in, bool snps_p_in,
		  bool cigar_extended_p_in, Cigar_action_T cigar_action_in);

#undef T
#endif

