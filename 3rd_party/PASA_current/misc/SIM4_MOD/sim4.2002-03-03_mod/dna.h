#ifndef SIM_DNA_H
#define SIM_DNA_H
/* $Id: dna.h,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $ */

#define DEFAULT_E       1 
#define DEFAULT_I       1
#define DEFAULT_M       0
#define DEFAULT_O       0
#define DEFAULT_V       1

void DNA_scores(argv_scores_t *ds, ss_t ss);
void DNA_scores_dflt(argv_scores_t *ds, ss_t ss, const argv_scores_t *dflt);

#endif
