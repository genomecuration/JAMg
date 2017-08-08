/* $Id: cigar.h 206761 2017-05-30 17:39:28Z twu $ */
#ifndef CIGAR_INCLUDED
#define CIGAR_INCLUDED

#include "bool.h"
#include "filestring.h"
#include "list.h"
#include "substring.h"
#include "stage3hr.h"


extern int
Cigar_length (List_T tokens);
extern void
Cigar_print_tokens (Filestring_T fp, List_T tokens);
extern void
Cigar_print_substrings (int *nindels, List_T *startp, List_T *startq, List_T *prevp, List_T *nextp, List_T *finalp, List_T *endp,
			Filestring_T fp, Stage3end_T stage3end,
			int querylength, int hardclip_low, int hardclip_high);
extern void
Cigar_print_halfdonor (Filestring_T fp, Substring_T donor, Stage3end_T this,
		       int querylength, int *hardclip_low, int *hardclip_high,
		       bool use_hardclip_p);
extern void
Cigar_print_halfacceptor (Filestring_T fp, Substring_T acceptor, Stage3end_T this,
			  int querylength, int *hardclip_low, int *hardclip_high,
			  bool use_hardclip_p);
extern void
Cigar_print_mate (Filestring_T fp, Stage3end_T mate, int mate_querylength, int mate_hardclip_low, int mate_hardclip_high);


#endif

