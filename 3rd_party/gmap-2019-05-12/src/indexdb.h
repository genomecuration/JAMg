/* $Id: indexdb.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef INDEXDB_INCLUDED
#define INDEXDB_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_64_BIT */
#endif

#include <stdio.h>
#include "access.h"
#include "types.h"
#include "univcoord.h"
#include "mode.h"
#include "genomicpos.h"
#include "bool.h"
#include "iitdef.h"

#ifdef PMAP
#include "alphabet.h"
#endif


#ifdef PMAP
#define SUFFICIENT_SUPPORT 9
#else
#define SUFFICIENT_SUPPORT 18
#endif


#define T Indexdb_T
typedef struct T *T;

#ifdef PMAP
extern void
Indexdb_setup (Width_T index1part_aa_in);
#else
extern void
Indexdb_setup (Width_T index1part_in);
#endif

extern void
Indexdb_free (T *old);
#ifndef PMAP
extern Width_T
Indexdb_interval (T this);
#endif
extern bool
Indexdb_positions_fileio_p (T this);
extern double
Indexdb_mean_size (T this, Mode_T mode, Width_T index1part);


typedef struct Indexdb_filenames_T *Indexdb_filenames_T;
struct Indexdb_filenames_T {
  char *pages_filename;
  char *pointers_filename;
  char *offsets_filename;
  char *positions_high_filename;
  char *positions_filename;

  char *pointers_basename_ptr;
  char *offsets_basename_ptr;
  char *positions_high_basename_ptr;
  char *positions_basename_ptr;

  char *pointers_index1info_ptr;
  char *offsets_index1info_ptr;
  char *positions_high_index1info_ptr;
  char *positions_index1info_ptr;
};


extern void
Indexdb_filenames_free (Indexdb_filenames_T *old);

extern Indexdb_filenames_T
Indexdb_get_filenames (int *compression_type,
#ifdef PMAP
		       Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
		       Width_T *index1part, Width_T *index1interval, char *genomesubdir,
		       char *fileroot, char *idx_filesuffix, char *snps_root,
		       Width_T required_index1part, Width_T required_interval,
		       bool offsets_only_p);

extern Univcoord_T *
Indexdb_point_one_shift (int *nentries, T this, Oligospace_T subst);
extern int
Indexdb_count_one_shift (T this, Oligospace_T subst, int nadjacent);


extern Positionsptr_T *
Indexdb_offsets_from_bitpack (char *offsetsmetafile, char *offsetsstrmfile, 
#ifdef PMAP
			      int alphabet_size, Width_T index1part_aa
#else
			      Width_T index1part
#endif
			      );

#if defined(HAVE_64_BIT) && defined(UTILITYP)
extern Hugepositionsptr_T *
Indexdb_offsets_from_bitpack_huge (char *bitpackpagesfile, char *offsetsmetafile, char *offsetsstrmfile,
#ifdef PMAP
				   int alphabet_size, Width_T index1part_aa
#else
				   Width_T index1part
#endif
				   );
#endif

#if 0
extern void
Indexdb_shmem_remove (char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		      Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		      Width_T required_index1part, Width_T required_interval, bool expand_offsets_p);
#endif

extern T
Indexdb_new_genome (Width_T *index1part, Width_T *index1interval,
		    char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		    Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		    Width_T required_index1part, Width_T required_interval, bool expand_offsets_p,
		    Access_mode_T offsetsstrm_access, Access_mode_T positions_access, bool sharedp,
		    bool multiple_sequences_p, bool preload_shared_memory_p, bool unload_shared_memory_p);
#ifndef UTILITYP
extern T
Indexdb_new_segment (char *genomicseg,
#ifdef PMAP
		     int alphabet_size, Width_T index1part_aa, bool watsonp,
#else
		     Width_T index1part,
#endif
		     Width_T index1interval);
#endif

#ifdef PMAP
extern Univcoord_T *
Indexdb_read (int *nentries, T this, Oligospace_T aaindex);
#else
extern Univcoord_T *
Indexdb_read (int *nentries, T this, Oligospace_T oligo);
extern UINT4 *
Indexdb_read_inplace (int *nentries,
#ifdef LARGE_GENOMES
		      unsigned char **positions_high,
#endif
		      T this, Oligospace_T oligo);
#endif

#ifdef LARGE_GENOMES
extern int
Indexdb_largeptr_with_diagterm (unsigned char **positions_high, UINT4 **positions_low,
				T this, Oligospace_T oligo, int diagterm);
#endif

extern int
Indexdb_ptr_with_diagterm (UINT4 **positions,
			   T this, Oligospace_T oligo, int diagterm);

extern int
Indexdb_count_with_diagterm (Positionsptr_T *ptr0, T this, Oligospace_T oligo, int diagterm);
extern Univcoord_T *
Indexdb_read_with_diagterm (int *nentries, T this, Oligospace_T oligo, int diagterm);
extern Univcoord_T *
Indexdb_read_with_diagterm_sizelimit (int *nentries, T this, Oligospace_T oligo, int diagterm,
				      int size_threshold);
extern Univcoord_T *
Indexdb_fill_with_diagterm (int nentries, Positionsptr_T ptr0, T this, int diagterm);


#undef T
#endif

