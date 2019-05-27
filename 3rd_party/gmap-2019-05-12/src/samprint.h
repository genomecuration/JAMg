/* $Id: samprint.h 218675 2019-03-16 01:25:48Z twu $ */
#ifndef SAMPRINT_INCLUDED
#define SAMPRINT_INCLUDED

#include <stdio.h>
#include "iit-read-univ.h"
#include "iit-read.h"
#include "genomicpos.h"
#include "types.h"
#include "substring.h"
#include "bool.h"
#include "intlist.h"
#include "filestring.h"

#include "chrnum.h"
#include "shortread.h"
#include "stage3hr.h"
#include "resulthr.h"
#include "genome.h"


extern void
SAM_setup (bool add_paired_nomappers_p_in, bool paired_flag_means_concordant_p_in,
	   bool omit_concordant_uniq_p_in, bool omit_concordant_mult_p_in, 
	   bool quiet_if_excessive_p_in, int maxpaths_report_in,
	   char *failedinput_root_in, bool fastq_format_p_in, bool hide_soft_clips_p_in, bool method_print_p_in,
	   bool clip_overlap_p_in, bool merge_overlap_p_in, bool merge_samechr_p_in,
	   bool sam_multiple_primaries_p_in, bool sam_sparse_secondaries_p_in,
	   bool force_xs_direction_p_in, bool md_lowercase_variant_p_in, IIT_T snps_iit_in,
	   bool find_dna_chimeras_p_in, bool novelsplicingp_in, IIT_T splicing_iit_in, int donor_typeint_in, int acceptor_typeint_in,
	   Univ_IIT_T chromosome_iit_in, Genome_T genome_in, Univ_IIT_T transcript_iit_in);

extern void
SAM_print_mate_cigar (Filestring_T fp, Stage3end_T mate, int mate_querylength, int mate_hardclip_low, int mate_hardclip_high);

extern void
SAM_print_pairs_nomapping (Filestring_T fp, char *abbrev, char *acc1, char *acc2, char *queryseq_ptr,
			   char *quality_string, int querylength, int quality_shift,
			   bool first_read_p, bool sam_paired_p, char *sam_read_group_id);

extern void
SAM_print_nomapping (Filestring_T fp, char *abbrev, Shortread_T queryseq, char *acc1, char *acc2,
		     Univ_IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int pathnum, int npaths_primary, int npaths_altloc, bool artificial_mate_p, int npaths_mate,
		     Stage3end_T mate, int mate_querylength, int mate_hardclip_low, int mate_hardclip_high,
		     int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p);

extern void
SAM_print (Filestring_T fp, Filestring_T fp_failedinput, char *abbrev, Stage3pair_T stage3pair,
	   Stage3end_T this, int querylength, char *acc1, char *acc2, int pathnum, int npaths_primary, int npaths_altloc,
	   int absmq_score, int first_absmq, int second_absmq, int mapq_score, Univ_IIT_T chromosome_iit,
	   Shortread_T queryseq, Shortread_T queryseq_mate, int pairedlength, int pair_relationship,

	   int hardclip_low, int hardclip_high,
	   Stage3end_T mate, int mate_querylength, int mate_hardclip_low, int mate_hardclip_high,

	   Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate,
	   int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p);

extern void
SAM_print_paired (Filestring_T fp, Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2,
		  Result_T result, Resulttype_T resulttype, Univ_IIT_T chromosome_iit,
		  Shortread_T queryseq1, Shortread_T queryseq2, bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, int quality_shift, char *sam_read_group_id);

#endif
