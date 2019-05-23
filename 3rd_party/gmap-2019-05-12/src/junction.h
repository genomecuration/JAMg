/* $Id: junction.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef JUNCTION_INCLUDED
#define JUNCTION_INCLUDED

typedef enum {NO_JUNCTION, INS_JUNCTION, DEL_JUNCTION, SPLICE_JUNCTION,
	      CHIMERA_JUNCTION, AMB_JUNCTION, END_JUNCTION} Junctiontype_T;

#include "types.h"
#include "genomicpos.h"
#include "bool.h"
#include "genome.h"
#include "list.h"
#include "listpool.h"


#define T Junction_T
typedef struct T *T;

extern void
Junction_print (T this);
extern void
Junction_print_list (List_T list);

extern void
Junction_free (T *old);
extern void
Junction_gc (List_T *list);

extern T
Junction_new_insertion (int nindels);
extern T
Junction_new_deletion (int nindels, Univcoord_T deletionpos);
extern T
Junction_new_splice (Chrpos_T splice_distance, int sensedir, double donor_prob, double acceptor_prob);

extern T
Junction_new_chimera (int sensedir, double donor_prob, double acceptor_prob);

extern T
Junction_new_generic (Univcoord_T left1, Univcoord_T left2, int querypos1, int querypos2,
		      Univcoord_T chroffset, bool plusp, int sensedir);

extern T
Junction_copy (T old);
extern List_T
Junction_copy_list (List_T old, Listpool_T listpool);


extern Junctiontype_T
Junction_type (T this);
extern char *
Junction_typestring (T this);
extern int
Junction_sensedir (T this);
extern double
Junction_prob (T this);
extern double
Junction_donor_prob (T this);
extern double
Junction_acceptor_prob (T this);
extern double
Junction_splice_score (T this);

extern int
Junction_nindels (T this);
extern int
Junction_adj (T this);
extern int
Junction_ninserts (T this);
extern int
Junction_total_ninserts (List_T list);

extern Univcoord_T
Junction_deletionpos (T this);
extern void
Junction_set_deletionpos (T this, Univcoord_T deletionpos);
extern char *
Junction_deletion_string (T this, Genome_T genome, bool plusp);
extern Chrpos_T
Junction_splice_distance (T this);
extern void
Junction_set_unambiguous (T this, Chrpos_T distance, double donor_prob, double acceptor_prob);
extern void
Junction_set_ambiguous (T this);

#undef T
#endif

