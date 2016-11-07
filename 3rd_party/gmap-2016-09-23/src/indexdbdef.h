/* $Id: indexdbdef.h 183995 2016-02-09 18:55:29Z twu $ */
#ifndef INDEXDBDEF_INCLUDED
#define INDEXDBDEF_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For HAVE_PTHREAD */
#endif

#include "genomicpos.h"
#include "access.h"
#include "types.h"

#ifdef PMAP
#include "alphabet.h"
#endif


#define BADVAL (Univcoord_T) -1

/* Compression types */
#define NO_COMPRESSION 0
#define BITPACK64_COMPRESSION 1


#define T Indexdb_T
struct T {
#ifdef PMAP
  Alphabet_T alphabet;
  int alphabet_size;
#endif

  int compression_type;
  Width_T index1part;
  Width_T index1interval;
  Blocksize_T blocksize;	/* e.g., 64 = 4^(15-12) */

#ifdef LARGE_GENOMES
  Access_T offsetspages_access;
  int offsetspages_shmid;
  key_t offsetspages_key;
  int offsetspages_fd;
  size_t offsetspages_len;
  UINT4 *offsetspages;
#endif

  Access_T offsetsmeta_access;
  int offsetsmeta_shmid;
  key_t offsetsmeta_key;
  int offsetsmeta_fd;
  size_t offsetsmeta_len;
  UINT4 *offsetsmeta;

  Access_T offsetsstrm_access;
  int offsetsstrm_shmid;
  key_t offsetsstrm_key;
  int offsetsstrm_fd;
  size_t offsetsstrm_len;
  UINT4 *offsetsstrm;

#ifdef LARGE_GENOMES
  Access_T positions_high_access;
  int positions_high_shmid;
  key_t positions_high_key;
  int positions_high_fd;
  size_t positions_high_len;

  Access_T positions_low_access;
  int positions_low_shmid;
  key_t positions_low_key;
  int positions_low_fd;
  size_t positions_low_len;

  unsigned char *positions_high;
  UINT4 *positions_low;
#else
  Access_T positions_access;
  int positions_shmid;
  key_t positions_key;
  int positions_fd;
  size_t positions_len;
  UINT4 *positions;		/* For small genomes, same as Univcoord_T */
#endif

#ifdef HAVE_PTHREAD
  pthread_mutex_t positions_read_mutex;
#endif
};

#undef T
#endif

