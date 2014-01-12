#include "libc.h"
#include "types.h"
#include "seq.h"
#include "misc.h"
#include "args.h"
#include "dna.h"

#ifndef __lint
static const char rcsid[] =
"$Id: dna.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif

static const argv_scores_t EIMOV = { 
	DEFAULT_E, 
	DEFAULT_I,
	DEFAULT_M,
	DEFAULT_O,
	DEFAULT_V
};

static void set_argv_scores(argv_scores_t *s, const argv_scores_t *const dflt)
{
        *s = *dflt;
}

/* DNA_scores -----------------------------------  set scoring matrix for DNA */
void DNA_scores_dflt(argv_scores_t *ds, ss_t ss, const argv_scores_t *dflt)
{
	int i, j, bad;

	ck_argc("DNA_scores");
	set_argv_scores(ds, dflt);

	for (i = 0; i < NACHARS; ++i)
		for (j = 0; j < NACHARS; ++j)
			ss[i][j] = ds->V;

	bad = -100*ds->M;
	for (i = 0; i < NACHARS; ++i)
		ss['X'][i] = ss[i]['X'] = bad;

	ss['a']['a'] = ss['c']['c'] = ss['g']['g'] = ss['t']['t'] = ds->M;
	ss['a']['A'] = ss['c']['C'] = ss['g']['G'] = ss['t']['T'] = ds->M;
	ss['A']['a'] = ss['C']['c'] = ss['G']['g'] = ss['T']['t'] = ds->M;
	ss['A']['A'] = ss['C']['C'] = ss['G']['G'] = ss['T']['T'] = ds->M;

	ss['a']['g'] = ss['g']['a'] = ss['c']['t'] = ss['t']['c'] = ds->I;
	ss['a']['G'] = ss['g']['A'] = ss['c']['T'] = ss['t']['C'] = ds->I;
	ss['A']['g'] = ss['G']['a'] = ss['C']['t'] = ss['T']['c'] = ds->I;
	ss['A']['G'] = ss['G']['A'] = ss['C']['T'] = ss['T']['C'] = ds->I;
}

void DNA_scores(argv_scores_t *ds, ss_t ss)
{
	DNA_scores_dflt(ds, ss, &EIMOV);
}
