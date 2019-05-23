/* $Id: localdbdef.h 218147 2019-01-16 21:28:41Z twu $ */
#ifndef LOCALDBDEF_INCLUDED
#define LOCALDBDEF_INCLUDED
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


/* A localdb has a region or loctable file, an offsetsmeta file, and
   offsetsstrm file, and a positions file.  It is designed to allow
   GSNAP to find 8-mers within a given genomic region quickly.

   A region is 2^16 = 65536 bp in the genome, and covers 4^8 = 65536
   distinct oligos.

   The loctable file has two UINT4 entries per region.  (Hence
   LOCTABLE_SIZE equals 2.)  The first entry is the pointer into the
   offsetsstrm file for the region, in terms of the number of 128-bit
   registers from the beginning of the offsetsstrm file.  The second
   entry is the pointer into the positions file for the region, in
   terms of the number of UINT2 positions from the beginning of the
   positions file.  There are a pair of terminating entries marking
   the end of the offsetsstrm file and positions file after the last
   genomic region.

   The offsetsmeta file has 65536 / 64 = 1024 entries per region,
   where 64 is the blocksize.  (Therefore regionsize equals 1024.)
   The value of 64 occurs because we store values in an
   epu16-bitpack64 bitpacking format, with 8 16-bit unsigned shorts
   per 128-bit SIMD register.  These shorts are cyclic differences,
   designed to allow a rapid calculation of the cumulative offset for
   the given block.  Description of this cyclic difference format is
   given in my paper in Algorithms for Molecular Biology on bitpacking
   methods for storing genomes (although that covers the case of 4
   32-bit unsigned ints per 128-bit SIMD register).

   There are therefore 1024 blocks per region.  Each block corresponds
   to a bmer value, or oligo / 64.  The offsetsmeta file contains two
   entries per block.  (Hence DIFFERENTIAL_METAINFO_SIZE equals 2.)
   The first entry is the pointer into the offsetsstrm file for the
   block, in terms of the number of 128-bit registers from the
   loctable pointer, as described above.  The second entry is the
   number of UINT2 positions from the beginning of the loctable
   pointer, as described above.  For a bmer value of 0, these entries
   would be 0 and 0.  We do not store these values, but instead store
   the values of bmer 1 through 1023, and then a terminating pair of
   entries.  These terminating entries are the values of the
   offsetsstrm and pointers locations for the next block.  The
   terminating entries make it possible to compute the number of
   128-bit registers by subtracting two adjacent offsetsstrm pointers,
   and the number of positions by subtracting two adjacent positions
   pointers.

   The positions file is a set of UINT2 values, which measure the
   relative position of some oligo from the beginning of the genomic
   region.  The beginning of the genomic region can be calculated as
   the region number * 65536 bp.  Since each genomic region contains
   65536 bp, the relative positions can be stored in an unsigned
   short.
 */



#define T Localdb_T
struct T {
#ifdef PMAP
  Alphabet_T alphabet;
  int alphabet_size;
#endif

  int compression_type;
  Width_T local1part;
  Width_T local1interval;
  Blocksize_T regionsize;	/* e.g., 1024 = 4^(8-3) for 8-mers */
  Blocksize_T blocksize;	/* 64 = 4^3 */

  Access_T loctable_access;
  int loctable_shmid;
  key_t loctable_key;
  int loctable_fd;
  size_t loctable_len;
  UINT4 *loctable;

  Access_T locoffsetsmeta_access;
  int locoffsetsmeta_shmid;
  key_t locoffsetsmeta_key;
  int locoffsetsmeta_fd;
  size_t locoffsetsmeta_len;
  UINT2 *locoffsetsmeta;	/* 4^3 entries with offsets and prefix sums */

  Access_T locoffsetsstrm_access;
  int locoffsetsstrm_shmid;
  key_t locoffsetsstrm_key;
  int locoffsetsstrm_fd;
  size_t locoffsetsstrm_len;
  UINT2 *locoffsetsstrm;	/* Really stored as 128-bit compressed registers */

  Access_T locpositions_access;
  int locpositions_shmid;
  key_t locpositions_key;
  int locpositions_fd;
  size_t locpositions_len;
  UINT2 *locpositions;		/* Values range from 0..65535 */
};

#undef T
#endif

