/* $Id: transcript.h 212789 2018-01-26 14:08:00Z twu $ */
#ifndef TRANSCRIPT_INCLUDED
#define TRANSCRIPT_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "filestring.h"
#include "list.h"
#include "bool.h"
#include "genomicpos.h"
#include "iit-read-univ.h"


#define T Transcript_T
typedef struct T *T;
struct T {
  int num;	      /* Index in transcriptome */
  int start;	      /* Starting coordinate on transcript, 0-based */
  int end; /* Ending coordinate.  If transcript_plusp is false, trend > trstart */
};


extern int Transcript_num (T this);


extern void
Transcript_setup (Chrpos_T pairmax_transcriptome_in, Chrpos_T expected_pairlength_in);
extern void
Transcript_free (T *old);
extern void
Transcript_gc (List_T *list);
extern T
Transcript_new (int num, int start, int end);
extern List_T
Transcript_copy_list (List_T old);
extern bool
Transcript_in_list_p (T x, List_T list);
extern void
Transcript_print_nums (List_T list);
extern void
Transcript_print_list (List_T list);

/* Depends on pairmax_transcriptome */
extern bool
Transcript_intersect_p (int *min_insertlength, List_T transcripts5, List_T transcripts3);

/* Depends on pairmax_transcriptome */
extern bool
Transcript_concordant_p (T transcript5, T transcript3);

/* Depends on pairmax_transcriptome, expected_pairlength */
extern void
Transcript_concordance (List_T *newtranscripts5, List_T *newtranscripts3, List_T transcripts5, List_T transcripts3);

extern void
Transcript_print_info (Filestring_T fp, List_T transcripts, Univ_IIT_T transcript_iit,
		       bool invertp);
extern void
Transcript_print_diff (Filestring_T fp, List_T transcripts, List_T common_transcripts,
		       Univ_IIT_T transcript_iit, bool invertp);

#undef T
#endif

