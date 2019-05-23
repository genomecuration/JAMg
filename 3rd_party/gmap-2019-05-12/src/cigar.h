/* $Id: cigar.h 214292 2018-03-19 23:32:03Z twu $ */
#ifndef CIGAR_INCLUDED
#define CIGAR_INCLUDED

#include "bool.h"
#include "filestring.h"
#include "list.h"
#include "substring.h"
#include "stage3hr.h"


extern void
Cigar_setup (bool cigar_extended_p_in, bool hide_soft_clips_p_in,
	     bool merge_samechr_p_in, bool md_lowercase_variant_p_in,
	     bool sam_hardclip_use_S_p_in, bool sam_insert_0M_p_in);
extern int
Cigar_length (List_T tokens);
extern void
Cigar_print_tokens (Filestring_T fp, List_T tokens);
extern int
Cigar_length_substrings (Stage3end_T stage3end, int querylength, int hardclip_low, int hardclip_high);
extern Filestring_T
Cigar_compute_main (int *part_hardclip_low, int *part_hardclip_high,
		    int *nindels, List_T *startp, List_T *startq,
		    List_T *prevp, List_T *nextp, List_T *finalp, List_T *endp,
		    bool *plusp, Chrnum_T *chrnum, Chrpos_T *chrpos,
		    Stage3end_T stage3end, int querylength, bool first_read_p, Stage3end_T mate,
		    int hardclip_low, int hardclip_high, bool hide_soft_clips_p);

extern Filestring_T
Cigar_compute_supplemental (int *part_hardclip_low, int *part_hardclip_high,
			    int *nindels, List_T *startp, List_T *startq,
			    List_T *prevp, List_T *nextp, List_T *finalp, List_T *endp,
			    bool *plusp, Chrnum_T *chrnum, Chrpos_T *chrpos,
			    Stage3end_T stage3end, int querylength, bool first_read_p, Stage3end_T mate,
			    int hardclip_low, int hardclip_high, bool hide_soft_clips_p);

#endif

