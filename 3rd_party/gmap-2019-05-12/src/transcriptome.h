/* $Id: transcriptome.h 210670 2017-10-18 17:26:20Z twu $ */
#ifndef TRANSCRIPTOME_INCLUDED
#define TRANSCRIPTOME_INCLUDED
#include "bool.h"
#include "types.h"
#include "chrnum.h"
#include "iit-read-univ.h"

#define T Transcriptome_T
typedef struct T *T;

extern void
Transcriptome_free (T *old);
extern T
Transcriptome_new (char *genomesubdir, char *fileroot, Univ_IIT_T chromosome_iit,
		   bool sharedp);
extern Chrnum_T
Transcriptome_chrnum (int *transcript_genestrand, T this, int trnum);
extern int
Transcriptome_exons (int **exonbounds, Chrpos_T **exonstarts, T this, int trnum);
extern bool
Transcriptome_genomic_bounded_p (int trnum, Chrnum_T chrbound, Chrpos_T lowbound, Chrpos_T highbound, T this);

#undef T
#endif


