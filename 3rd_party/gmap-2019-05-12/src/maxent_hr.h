/* $Id: maxent_hr.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef MAXENT_HR_INCLUDED
#define MAXENT_HR_INCLUDED

#include "genomicpos.h"
#include "types.h"
#include "univcoord.h"


extern void
Maxent_hr_setup (Genomecomp_T *ref_blocks_in, Genomecomp_T *snp_blocks_in);

extern double
Maxent_hr_donor_prob (Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_acceptor_prob (Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_antidonor_prob (Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_antiacceptor_prob (Univcoord_T splice_pos, Univcoord_T chroffset);

#endif

