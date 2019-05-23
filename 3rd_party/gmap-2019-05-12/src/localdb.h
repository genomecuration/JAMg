/* $Id: localdb.h 218530 2019-03-04 23:52:01Z twu $ */
#ifndef LOCALDB_INCLUDED
#define LOCALDB_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

typedef struct Localdb_T *Localdb_T;

#include <stdio.h>
#include "access.h"
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "bool.h"

#include "compress.h"
#include "genome.h"


#define T Localdb_T

typedef struct Localdb_filenames_T *Localdb_filenames_T;
struct Localdb_filenames_T {
  char *loctable_filename;
  char *locpointers_filename;
  char *locoffsets_filename;
  char *locpositions_filename;

  char *loctable_basename_ptr;
  char *locpointers_basename_ptr;
  char *locoffsets_basename_ptr;
  char *locpositions_basename_ptr;

  char *loctable_local1info_ptr;
  char *locpointers_local1info_ptr;
  char *locoffsets_local1info_ptr;
  char *locpositions_local1info_ptr;
};


extern void
Localdb_filenames_free (Localdb_filenames_T *old);


extern Localdb_filenames_T
Localdb_get_filenames (
#ifdef PMAP
		       Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
		       Width_T *local1part, Width_T *local1interval,
		       char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
		       Width_T required_local1part, Width_T required_interval, bool offsets_only_p);


extern T
Localdb_new_genome (Width_T *local1part, Width_T *local1interval,
		    char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		    Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		    Width_T required_local1part, Width_T required_interval,
		    Access_mode_T locoffsetsstrm_access, Access_mode_T locpositions_access, bool sharedp,
		    bool multiple_sequences_p, bool unload_shared_memory_p);

extern void
Localdb_free (T *old);

extern Univcoord_T *
Localdb_get_diagonals (int *nentries, T this, char *queryptr, int pos5, int pos3,
		       Univcoord_T low, Univcoord_T high, bool plusp, int genestrand,
		       Univcoord_T **stream_alloc, int *streamsize_alloc,
		       bool remove_repetitive_p);

extern void
Localdb_test (T this, Oligospace_T oligo, int regioni);

extern void
Localdb_setup (int local1part, int mode);

#undef T
#endif
