static char rcsid[] = "$Id: cmetindex.c 218151 2019-01-17 05:36:54Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include <ctype.h>		/* For toupper */
#include <sys/mman.h>		/* For munmap */
#include <math.h>		/* For qsort */
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For lseek and close */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#include "mem.h"
#include "fopen.h"
#include "access.h"
#include "types.h"		/* For Positionsptr_T and Oligospace_T */
#include "mode.h"
#include "filesuffix.h"

#include "cmet.h"

#include "assert.h"
#include "bool.h"
#include "genomicpos.h"
#include "iit-read-univ.h"
#include "indexdb.h"
#include "indexdbdef.h"
#include "indexdb-write.h"
#include "genome.h"
#include "genome128_hr.h"

#include "bitpack64-write.h"
#include "bitpack64-read.h"
#include "bitpack64-readtwo.h"
#include "bitpack64-access.h"
#include "bitpack64-incr.h"

#include "localdb.h"
#include "localdbdef.h"
#include "epu16-bitpack64-write.h"
#include "epu16-bitpack64-read.h"
#include "epu16-bitpack64-access.h"
#include "epu16-bitpack64-incr.h"

#include "datadir.h"
#include "getopt.h"


#define POSITIONS8_HIGH_SHIFT 32
#define POSITIONS8_LOW_MASK 0xFFFFFFFF

#define LOCTABLE_SIZE 2
#define DIFFERENTIAL_METAINFO_SIZE 2
#define BLOCKSIZE 64
#define BUFFER_SIZE 1000000

#define MONITOR_INTERVAL 10000000
#define MONITOR_INTERVAL_REGION 100000


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif



static char *user_sourcedir = NULL;
static char *user_destdir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static int compression_type;

static int index1part = 15;
static int required_index1part = 0;
static int index1interval;
static int required_index1interval = 0;

static int local1part = 8;
static int required_local1part = 0;
static int local1interval = 1;
static int required_local1interval = 0;

static char *snps_root = NULL;


static struct option long_options[] = {
  /* Input options */
  {"sourcedir", required_argument, 0, 'F'},	/* user_sourcedir */
  {"destdir", required_argument, 0, 'D'},	/* user_destdir */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part */
  {"sampling", required_argument, 0, 'q'}, /* required_interval */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"usesnps", required_argument, 0, 'v'}, /* snps_root */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"CMETINDEX: Builds GMAP index files for methylated genomes\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage ();


static Oligospace_T
power (int base, int exponent) {
  Oligospace_T result = 1UL;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


/*                      87654321 */
#define LOW_TWO_BITS  0x00000003


#if 0
static Oligospace_T
reduce_oligo_old (Oligospace_T oligo, int oligosize) {
  Oligospace_T reduced = 0U, lowbits;
  int i;

  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    reduced >>= 2;
    switch (lowbits) {
    case RIGHT_A: break;
    case RIGHT_C: reduced |= LEFT_T; break;
    case RIGHT_G: reduced |= LEFT_G; break;
    case RIGHT_T: reduced |= LEFT_T; break;
    }
    oligo >>= 2;
  }

  for ( ; i < 16; i++) {
    reduced >>= 2;
  }

  return reduced;
}
#endif


#if 0
static char *
shortoligo_nt (Oligospace_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Oligospace_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return nt;
}
#endif


#if 0
/* Old method */
static Positionsptr_T *
compute_offsets_ct_using_array (Positionsptr_T *oldoffsets, Oligospace_T oligospace, Oligospace_T mask) {
  Positionsptr_T *offsets;
  Oligospace_T oligoi, reduced;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  for (oligoi = 0; oligoi < oligospace; oligoi++) {
#if 0
    if (reduce_oligo(oligoi) != reduce_oligo_old(oligoi,index1part)) {
      abort();
    }
#endif
    reduced = Cmet_reduce_ct(oligoi) & mask;
    debug(
	  nt1 = shortoligo_nt(oligoi,index1part);
	  nt2 = shortoligo_nt(reduced,index1part);
	  printf("For oligo %s, updating sizes for %s from %u",nt1,nt2,offsets[reduced+1]);
	  );
#ifdef WORDS_BIGENDIAN
    /*size*/offsets[reduced+1] += (Bigendian_convert_uint(oldoffsets[oligoi+1]) - Bigendian_convert_uint(oldoffsets[oligoi]));
#else
    /*size*/offsets[reduced+1] += (oldoffsets[oligoi+1] - oldoffsets[oligoi]);
#endif
    debug(
	  printf(" to %u\n",offsets[reduced+1]);
	  FREE(nt2);
	  FREE(nt1);
	  );
  }

  offsets[0] = 0U;
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi-1] + /*size*/offsets[oligoi];
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %X: %u\n",oligoi,offsets[oligoi]);
	  });
  }

  return offsets;
}
#endif


#if 0
/* Old method */
static Positionsptr_T *
compute_offsets_ga_using_array (Positionsptr_T *oldoffsets, Oligospace_T oligospace, Oligospace_T mask) {
  Positionsptr_T *offsets;
  Oligospace_T oligoi, reduced;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    reduced = Cmet_reduce_ga(oligoi) & mask;
    debug(
	  nt1 = shortoligo_nt(oligoi,index1part);
	  nt2 = shortoligo_nt(reduced,index1part);
	  printf("For oligo %s, updating sizes for %s from %u",nt1,nt2,offsets[reduced+1]);
	  );
#ifdef WORDS_BIGENDIAN
    /*size*/offsets[reduced+1] += (Bigendian_convert_uint(oldoffsets[oligoi+1]) - Bigendian_convert_uint(oldoffsets[oligoi]));
#else
    /*size*/offsets[reduced+1] += (oldoffsets[oligoi+1] - oldoffsets[oligoi]);
#endif
    debug(
	  printf(" to %u\n",offsets[reduced+1]);
	  FREE(nt2);
	  FREE(nt1);
	  );
  }

  offsets[0] = 0U;
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi-1] + /*size*/offsets[oligoi];
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %X: %u\n",oligoi,offsets[oligoi]);
	  });
  }

  return offsets;
}
#endif


static UINT4
local_offsets_ct_using_bitpack (UINT4 *npositions, FILE *ptrs_fp, FILE *comp_fp,
				UINT2 *oldoffsetsmeta, UINT2 *oldoffsetsstrm,
				Localspace_T localspace, Localspace_T mask) {
  UINT4 nregisters;
  Localspace_T bmerspace;
  char *packsizes;
  UINT2 **bitpacks;
  UINT2 increment;
  int new_packsize;
  UINT4 offsets[BLOCKSIZE+1];
  Localspace_T oligoi, reduced, bmer;
  int ii;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  bmerspace = localspace/BLOCKSIZE;
  packsizes = (char *) CALLOC(bmerspace,sizeof(char));
  bitpacks = (UINT2 **) CALLOC(bmerspace,sizeof(UINT2 *));

  for (oligoi = 0; oligoi < localspace; oligoi += BLOCKSIZE) {
    Epu16_bitpack64_block_offsets(offsets,oligoi,oldoffsetsmeta,oldoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((increment = offsets[ii+1] - offsets[ii]) > 0) {
	reduced = Cmet_reduce_ct_local(oligoi + ii) & mask;
	bmer = reduced/BLOCKSIZE;
	if ((new_packsize = Epu16_bitpack64_access_new_packsize(reduced,/*old_packsize*/(int) packsizes[bmer],
								bitpacks[bmer],increment)) != (int) packsizes[bmer]) {
	  bitpacks[bmer] = Epu16_bitpack64_realloc_multiple((int) packsizes[bmer],new_packsize,bitpacks[bmer]);
	  packsizes[bmer] = (char) new_packsize;
	}
	Epu16_bitpack64_add_bitpack(reduced,(int) packsizes[bmer],bitpacks[bmer],increment);
      }
    }
  }

  nregisters = Epu16_bitpack64_append_differential(&(*npositions),ptrs_fp,comp_fp,
						   packsizes,bitpacks,localspace);

  for (bmer = 0; bmer < bmerspace; bmer++) {
    FREE(bitpacks[bmer]);
  }
  FREE(bitpacks);
  FREE(packsizes);

  return nregisters;
}


static void
index_offsets_ct_using_bitpack (char *new_pointers_filename, char *new_offsets_filename,
				UINT4 *oldoffsetsmeta, UINT4 *oldoffsetsstrm,
				Oligospace_T oligospace, Oligospace_T mask) {
  Oligospace_T bmerspace;
  char *packsizes;
  UINT4 **bitpacks;
  UINT4 increment;
  int new_packsize;
  UINT4 offsets[BLOCKSIZE+1];
  Oligospace_T oligoi, reduced, bmer;
  int ii;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  bmerspace = oligospace/BLOCKSIZE;
  packsizes = (char *) CALLOC(bmerspace,sizeof(char));
  bitpacks = (UINT4 **) CALLOC(bmerspace,sizeof(UINT4 *));

  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    Bitpack64_block_offsets(offsets,oligoi,oldoffsetsmeta,oldoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((increment = offsets[ii+1] - offsets[ii]) > 0) {
	reduced = Cmet_reduce_ct(oligoi + ii) & mask;
	bmer = reduced/BLOCKSIZE;
	if ((new_packsize = Bitpack64_access_new_packsize(reduced,/*old_packsize*/(int) packsizes[bmer],
							  bitpacks[bmer],increment)) != (int) packsizes[bmer]) {
	  bitpacks[bmer] = Bitpack64_realloc_multiple((int) packsizes[bmer],new_packsize,bitpacks[bmer]);
	  packsizes[bmer] = (char) new_packsize;
	}
	Bitpack64_add_bitpack(reduced,(int) packsizes[bmer],bitpacks[bmer],increment);
      }
    }
  }

  Bitpack64_write_differential_bitpacks(new_pointers_filename,new_offsets_filename,packsizes,bitpacks,oligospace);

  for (bmer = 0; bmer < bmerspace; bmer++) {
    FREE(bitpacks[bmer]);
  }
  FREE(bitpacks);
  FREE(packsizes);

  return;
}


static UINT4
local_offsets_ga_using_bitpack (UINT4 *npositions, FILE *ptrs_fp, FILE *comp_fp,
				UINT2 *oldoffsetsmeta, UINT2 *oldoffsetsstrm,
				Localspace_T localspace, Localspace_T mask) {
  UINT4 nregisters;
  Localspace_T bmerspace;
  char *packsizes;
  UINT2 **bitpacks;
  UINT2 increment;
  int new_packsize;
  UINT4 offsets[BLOCKSIZE+1];
  Localspace_T oligoi, reduced, bmer;
  int ii;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  bmerspace = localspace/BLOCKSIZE;
  packsizes = (char *) CALLOC(bmerspace,sizeof(char));
  bitpacks = (UINT2 **) CALLOC(bmerspace,sizeof(UINT2 *));

  for (oligoi = 0; oligoi < localspace; oligoi += BLOCKSIZE) {
    Epu16_bitpack64_block_offsets(offsets,oligoi,oldoffsetsmeta,oldoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((increment = offsets[ii+1] - offsets[ii]) > 0) {
	reduced = Cmet_reduce_ga_local(oligoi + ii) & mask;
	bmer = reduced/BLOCKSIZE;
	if ((new_packsize = Epu16_bitpack64_access_new_packsize(reduced,/*old_packsize*/(int) packsizes[bmer],
								bitpacks[bmer],increment)) != (int) packsizes[bmer]) {
	  bitpacks[bmer] = Epu16_bitpack64_realloc_multiple((int) packsizes[bmer],new_packsize,bitpacks[bmer]);
	  packsizes[bmer] = (char) new_packsize;
	}
	Epu16_bitpack64_add_bitpack(reduced,(int) packsizes[bmer],bitpacks[bmer],increment);
      }
    }
  }

  nregisters = Epu16_bitpack64_append_differential(&(*npositions),ptrs_fp,comp_fp,
						   packsizes,bitpacks,localspace);

  for (bmer = 0; bmer < bmerspace; bmer++) {
    FREE(bitpacks[bmer]);
  }
  FREE(bitpacks);
  FREE(packsizes);

  return nregisters;
}

static void
index_offsets_ga_using_bitpack (char *new_pointers_filename, char *new_offsets_filename,
				UINT4 *oldoffsetsmeta, UINT4 *oldoffsetsstrm,
				Oligospace_T oligospace, Oligospace_T mask) {
  Oligospace_T bmerspace;
  char *packsizes;
  UINT4 **bitpacks;
  UINT4 increment;
  int new_packsize;
  UINT4 offsets[BLOCKSIZE+1];
  Oligospace_T oligoi, reduced, bmer;
  int ii;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  bmerspace = oligospace/BLOCKSIZE;
  packsizes = (char *) CALLOC(bmerspace,sizeof(char));
  bitpacks = (UINT4 **) CALLOC(bmerspace,sizeof(UINT4 *));

  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    Bitpack64_block_offsets(offsets,oligoi,oldoffsetsmeta,oldoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((increment = offsets[ii+1] - offsets[ii]) > 0) {
	reduced = Cmet_reduce_ga(oligoi + ii) & mask;
	bmer = reduced/BLOCKSIZE;
	if ((new_packsize = Bitpack64_access_new_packsize(reduced,/*old_packsize*/(int) packsizes[bmer],
							  bitpacks[bmer],increment)) != (int) packsizes[bmer]) {
	  bitpacks[bmer] = Bitpack64_realloc_multiple((int) packsizes[bmer],new_packsize,bitpacks[bmer]);
	  packsizes[bmer] = (char) new_packsize;
	}
	Bitpack64_add_bitpack(reduced,(int) packsizes[bmer],bitpacks[bmer],increment);
      }
    }
  }

  Bitpack64_write_differential_bitpacks(new_pointers_filename,new_offsets_filename,packsizes,bitpacks,oligospace);

  for (bmer = 0; bmer < bmerspace; bmer++) {
    FREE(bitpacks[bmer]);
  }
  FREE(bitpacks);
  FREE(packsizes);

  return;
}


#if 0
/* Needs work */
static void
sort_positions_local_snps (FILE *positions_fp, UINT2 *positions2,
			   UINT2 *newoffsetsmeta, UINT2 *newoffsetsstrm, Localspace_T localspace) {
  Localspace_T oligoi;
  UINT2 block_start, block_end, j, npositions;
  int ii;

  UINT4 newoffsets[BLOCKSIZE+1];

#if 0
  /* For snps_root */
  FILE *ptrs_fp, *comp_fp;
  UINT4 *ptrs;
  UINT4 ptri;
  int nregisters;
  UINT4 ascending[BLOCKSIZE+1];
  UINT2 diffs[BLOCKSIZE];
  UINT4 totalcount;
  int packsize;

  /* Buffer is used to avoid frequent writes to the file */
  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;
#endif


  /* Sort positions in each block.  For snps_root, use code from Epu16_bitpack64_append_differential */
  if (snps_root) {
    fprintf(stderr,"Not supported\n");
    abort();
#if 0
    ptrs = (UINT4 *) CALLOC(((oligospace + BLOCKSIZE)/BLOCKSIZE + 1) * DIFFERENTIAL_METAINFO_SIZE,sizeof(UINT4));
    ptri = 0;

    if ((comp_fp = FOPEN_WRITE_BINARY(compfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",compfile);
      exit(9);
    }
    buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
    buffer_i = 0;
    nregisters = 0;
#endif
  }

  /* totalcount = 0; -- For snps_root */
  for (oligoi = 0; oligoi + BLOCKSIZE <= oligospace; oligoi += BLOCKSIZE) {
    if (snps_root) {
#if 0
      ptrs[ptri++] = nregisters;
      ptrs[ptri++] = totalcount;
#endif
    }
    /* ascending[0] = totalcount; -- For snps_root */
    Epu16_bitpack64_block_offsets(newoffsets,oligoi,newoffsetsmeta,newoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      block_start = newoffsets[ii];
      block_end = newoffsets[ii+1];
      if ((npositions = block_end - block_start) > 0) {
	qsort(&(positions2[block_start]),npositions,sizeof(UINT2),UINT2_compare);
	if (snps_root == NULL) {
	  FWRITE_UINTS(&(positions2[block_start]),npositions,positions_fp);
	} else {
#if 0
	  FWRITE_UINT(positions4[block_start],positions_fp);
	  for (j = block_start+1; j < block_end; j++) {
	    if (positions4[j] == positions4[j-1]) {
	      npositions--;
	    } else {
	      FWRITE_UINT(positions4[j],positions_fp);
	    }
	  }
#endif
	}
	/* totalcount += npositions; -- For snps_root */
      }
      /* ascending[ii+1] = totalcount; -- For snps_root */
    }
    if (snps_root) {
#if 0
      packsize = Bitpack64_compute_q4_diffs_bidir(diffs,ascending);
      buffer_i = Bitpack64_write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);
      nregisters += packsize / 2;
#endif
    }
  }

  if (snps_root) {
#if 0    
    /* Write the final pointer, which will point after the end of the file */
    ptrs[ptri++] = nregisters;	/* In 128-bit registers */

    /* Value for end of block */
    ptrs[ptri++] = totalcount;
    
    if ((ptrs_fp = FOPEN_WRITE_BINARY(ptrsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",ptrsfile);
      exit(9);
    } else {
      FWRITE_UINTS(ptrs,ptri,ptrs_fp);
      FREE(ptrs);
      fclose(ptrs_fp);
    }
    
    /* Empty buffer */
    if (buffer_i > 0) {
      FWRITE_UINTS(buffer,buffer_i,comp_fp);	
      buffer_i = 0;
    }
    FREE(buffer);
    fclose(comp_fp);
#endif
  }

  return;
}
#endif


#if 0
/* Does not support SNPs.  Written in-line into compute_ct_local and compute_ag_local */
static void
sort_positions_local (FILE *positions_fp, UINT2 *positions2,
		      UINT2 *newoffsetsmeta, UINT2 *newoffsetsstrm, Localspace_T localspace) {
  Localspace_T oligoi;
  UINT2 block_start, block_end, j, npositions;
  UINT4 newoffsets[BLOCKSIZE+1];
  int ii;

  /* Sort positions in each block */
  for (oligoi = 0; oligoi + BLOCKSIZE <= localspace; oligoi += BLOCKSIZE) {
    Epu16_bitpack64_block_offsets(newoffsets,oligoi,newoffsetsmeta,newoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      block_start = newoffsets[ii];
      block_end = newoffsets[ii+1];
      if ((npositions = block_end - block_start) > 0) {
	qsort(&(positions2[block_start]),npositions,sizeof(UINT2),UINT2_compare);
	FWRITE_UINTS(&(positions2[block_start]),npositions,positions_fp);
      }
    }
  }

  return;
}
#endif


static void
sort_8mers (unsigned char *positions8_high, UINT4 *positions8_low, Positionsptr_T npositions) {
  UINT8 *positions8;
  Positionsptr_T i;

  positions8 = (UINT8 *) MALLOC(npositions*sizeof(UINT8));
  for (i = 0; i < npositions; i++) {
    positions8[i] = ((UINT8) positions8_high[i] << 32) + positions8_low[i];
  }
  qsort(positions8,npositions,sizeof(UINT8),UINT8_compare);
  for (i = 0; i < npositions; i++) {
    positions8_high[i] = positions8[i] >> POSITIONS8_HIGH_SHIFT;
    positions8_low[i] = positions8[i] & POSITIONS8_LOW_MASK;
  }
  return;
}


static void
sort_positions (char *ptrsfile, char *compfile, FILE *positions_high_fp, FILE *positions_fp,
		unsigned char *positions8_high, UINT4 *positions8_low, UINT4 *positions4,
		UINT4 *newoffsetsmeta, UINT4 *newoffsetsstrm, Oligospace_T oligospace, bool coord_values_8p) {
  Oligospace_T oligoi;
  Positionsptr_T block_start, block_end, j, npositions;

  FILE *ptrs_fp, *comp_fp;
  UINT4 *ptrs, nregisters;
  UINT4 ptri;
  int ii;

  /* Buffer is used to avoid frequent writes to the file */
  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT4 newoffsets[BLOCKSIZE+1];
  UINT4 diffs[BLOCKSIZE], ascending[BLOCKSIZE+1];
  UINT4 totalcount;
  int packsize;


  /* Sort positions in each block.  For snps_root, use code from Bitpack64_write_differential */
  if (snps_root) {
    ptrs = (UINT4 *) CALLOC(((oligospace + BLOCKSIZE)/BLOCKSIZE + 1) * DIFFERENTIAL_METAINFO_SIZE,sizeof(UINT4));
    ptri = 0;

    if ((comp_fp = FOPEN_WRITE_BINARY(compfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",compfile);
      exit(9);
    }
    buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
    buffer_i = 0;
    nregisters = 0;
  }

  totalcount = 0;
  for (oligoi = 0; oligoi + BLOCKSIZE <= oligospace; oligoi += BLOCKSIZE) {
    if (snps_root) {
      ptrs[ptri++] = nregisters;
      ptrs[ptri++] = totalcount;
    }
    ascending[0] = totalcount;
    Bitpack64_block_offsets(newoffsets,oligoi,newoffsetsmeta,newoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((oligoi + ii) % MONITOR_INTERVAL == 0) {
	fprintf(stderr,".");
      }
      block_start = newoffsets[ii];
      block_end = newoffsets[ii+1];
      if ((npositions = block_end - block_start) > 0) {
	if (coord_values_8p == true) {
	  sort_8mers(&(positions8_high[block_start]),&(positions8_low[block_start]),npositions);
	  if (snps_root == NULL) {
	    /* FWRITE_UINT8S(&(positions8[block_start]),npositions,positions_fp); */
	    FWRITE_CHARS(&(positions8_high[block_start]),npositions,positions_high_fp);
	    FWRITE_UINTS(&(positions8_low[block_start]),npositions,positions_fp);
	  } else {
	    /* FWRITE_UINT8(positions8[block_start],positions_fp); */
	    FWRITE_CHAR(positions8_high[block_start],positions_high_fp);
	    FWRITE_UINT(positions8_low[block_start],positions_fp);
	    for (j = block_start+1; j < block_end; j++) {
	      if (positions8_high[j] == positions8_high[j-1] && positions8_low[j] == positions8_low[j-1]) {
		npositions--;
	      } else {
		/* FWRITE_UINT8(positions8[j],positions_fp); */
		FWRITE_CHAR(positions8_high[j],positions_high_fp);
		FWRITE_UINT(positions8_low[j],positions_fp);
	      }
	    }
	  }

	} else {
	  qsort(&(positions4[block_start]),npositions,sizeof(UINT4),UINT4_compare);
	  if (snps_root == NULL) {
	    FWRITE_UINTS(&(positions4[block_start]),npositions,positions_fp);
	  } else {
	    FWRITE_UINT(positions4[block_start],positions_fp);
	    for (j = block_start+1; j < block_end; j++) {
	      if (positions4[j] == positions4[j-1]) {
		npositions--;
	      } else {
		FWRITE_UINT(positions4[j],positions_fp);
	      }
	    }
	  }
	}
	totalcount += npositions;
      }
      ascending[ii+1] = totalcount;
    }
    if (snps_root) {
      packsize = Bitpack64_compute_q4_diffs_bidir(diffs,ascending);
      buffer_i = Bitpack64_write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);
      nregisters += packsize / 2;
    }
  }

  if (oligoi <= oligospace) {
    /* Finish last block (containing 1 entry, but expand to 64) */
    if (snps_root) {
      ptrs[ptri++] = nregisters;
      ptrs[ptri++] = totalcount;
    }
    ascending[0] = totalcount;
    for (ii = 0; ii <= (int) (oligospace - oligoi); ii++) {
      block_start = Bitpack64_read_two(&block_end,oligoi+ii,newoffsetsmeta,newoffsetsstrm);
      if ((npositions = block_end - block_start) > 0) {
	if (coord_values_8p == true) {
	  sort_8mers(&(positions8_high[block_start]),&(positions8_low[block_start]),npositions);
	  if (snps_root == NULL) {
	    /* FWRITE_UINT8S(&(positions8[block_start]),npositions,positions_fp); */
	    FWRITE_CHARS(&(positions8_high[block_start]),npositions,positions_high_fp);
	    FWRITE_UINTS(&(positions8_low[block_start]),npositions,positions_fp);
	  } else {
	    /* FWRITE_UINT8(positions8[block_start],positions_fp); */
	    FWRITE_CHAR(positions8_high[block_start],positions_high_fp);
	    FWRITE_UINT(positions8_low[block_start],positions_fp);
	    for (j = block_start+1; j < block_end; j++) {
	      if (positions8_high[j] == positions8_high[j-1] && positions8_low[j] == positions8_low[j-1]) {
		npositions--;
	      } else {
		/* FWRITE_UINT8(positions8[j],positions_fp); */
		FWRITE_CHAR(positions8_high[j],positions_high_fp);
		FWRITE_UINT(positions8_low[j],positions_fp);
	      }
	    }
	  }

	} else {
	  qsort(&(positions4[block_start]),npositions,sizeof(UINT4),UINT4_compare);
	  if (snps_root == NULL) {
	    FWRITE_UINTS(&(positions4[block_start]),npositions,positions_fp);
	  } else {
	    FWRITE_UINT(positions4[block_start],positions_fp);
	    for (j = block_start+1; j < block_end; j++) {
	      if (positions4[j] == positions4[j-1]) {
		npositions--;
	      } else {
		FWRITE_UINT(positions4[j],positions_fp);
	      }
	    }
	  }
	}
	totalcount += npositions;
      }
      ascending[ii+1] = totalcount;
    }
    if (snps_root) {
      for ( ; ii < BLOCKSIZE; ii++) {
	ascending[ii+1] = totalcount;
      }
      packsize = Bitpack64_compute_q4_diffs_bidir(diffs,ascending);
      buffer_i = Bitpack64_write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);
      nregisters += packsize / 2;
    }
  }

  if (snps_root) {
    /* Write the final pointer, which will point after the end of the file */
    ptrs[ptri++] = nregisters;	/* In 128-bit registers */

    /* Value for end of block */
    ptrs[ptri++] = totalcount;
    
    if ((ptrs_fp = FOPEN_WRITE_BINARY(ptrsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",ptrsfile);
      exit(9);
    } else {
      FWRITE_UINTS(ptrs,ptri,ptrs_fp);
      FREE(ptrs);
      fclose(ptrs_fp);
    }
    
    /* Empty buffer */
    if (buffer_i > 0) {
      FWRITE_UINTS(buffer,buffer_i,comp_fp);	
      buffer_i = 0;
    }
    FREE(buffer);
    fclose(comp_fp);
  }

  return;
}


#if 0
/*                                       A  T  G  T */
static unsigned char ct_conversion[4] = {0, 3, 2, 3};
static char CT_CHARTABLE[4] = {'A','T','G','T'};
#endif

static void
compute_ct_local (char *new_region_filename,char *new_pointers_filename,
		  char *new_offsets_filename, char *new_positions_filename,
		  Localdb_T localdb, Localspace_T localspace, Localspace_T mask,
		  Univcoord_T genomelength) {
  FILE *region_fp, *ptrs_fp, *comp_fp, *positions_fp;
  UINT4 ascending = 0, total_npositions = 0;

  Univcoord_T position;
  int regioni;
  UINT4 old_region_strm_start, old_region_positions_start, old_region_positions_end;
  UINT4 new_region_strm_start;

  UINT2 *positions2, *oldpositions2;
  Localspace_T oligoi, reduced;
  UINT4 newoffsets[BLOCKSIZE+1], oldoffsets[BLOCKSIZE+1];

  UINT4 *new_loctable;
  UINT2 *new_locoffsetsmeta, *new_locoffsetsstrm;
  UINT2 *newoffsetsmeta, *newoffsetsstrm, *oldoffsetsmeta, *oldoffsetsstrm;

  UINT4 *counter, preunique_totalcounts, npositions, block_start, block_end, j, offset;
  size_t new_loctable_len, new_locoffsetsmeta_len, new_locoffsetsstrm_len;
  int ii;
#ifdef HAVE_MMAP
  int new_loctable_fd, new_locoffsetsmeta_fd, new_locoffsetsstrm_fd;
#else
  Access_T loctable_access, offsetsmeta_access, offsetsstrm_access;
  double seconds;
#endif
#ifdef DEBUG
  char *nt1, *nt2;
#endif


  if ((region_fp = FOPEN_WRITE_BINARY(new_region_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_region_filename);
    exit(9);
  }
  if ((ptrs_fp = FOPEN_WRITE_BINARY(new_pointers_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_pointers_filename);
    exit(9);
  }
  if ((comp_fp = FOPEN_WRITE_BINARY(new_offsets_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_offsets_filename);
    exit(9);
  }
  if ((positions_fp = FOPEN_WRITE_BINARY(new_positions_filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",new_positions_filename);
    exit(9);
  }


  fprintf(stderr,"Rearranging CT offsets for localdb...");
  for (position = 0, regioni = 0; position < genomelength; position += localspace, regioni++) {
    /* regioni = position/65536; */
    if (regioni % MONITOR_INTERVAL_REGION == 0) {
      fprintf(stderr,".");
    }

    FWRITE_UINT(ascending,region_fp);
    FWRITE_UINT(total_npositions,region_fp);

    oldoffsetsmeta = &(localdb->locoffsetsmeta[regioni * localdb->regionsize * DIFFERENTIAL_METAINFO_SIZE]);
    old_region_strm_start = localdb->loctable[regioni * LOCTABLE_SIZE];
    oldoffsetsstrm = &(localdb->locoffsetsstrm[old_region_strm_start * 8]); /* 8 shorts per 128-bit register */

    ascending += local_offsets_ct_using_bitpack(&npositions,ptrs_fp,comp_fp,
						oldoffsetsmeta,oldoffsetsstrm,localspace,mask);
    total_npositions += npositions;
  }

  /* Write final entries to region_fp */
  FWRITE_UINT(ascending,region_fp);
  FWRITE_UINT(total_npositions,region_fp);
  fclose(region_fp);
  fclose(comp_fp);
  fclose(ptrs_fp);
  fprintf(stderr,"done\n");


  /* Point to offsets and revise */
  fprintf(stderr,"Sorting CT positions for localdb...");
#ifdef HAVE_MMAP
  new_loctable = (UINT4 *) Access_mmap(&new_loctable_fd,&new_loctable_len,new_region_filename,/*randomp*/false);
  new_locoffsetsmeta = (UINT2 *) Access_mmap(&new_locoffsetsmeta_fd,&new_locoffsetsmeta_len,new_pointers_filename,/*randomp*/false);
  new_locoffsetsstrm = (UINT2 *) Access_mmap(&new_locoffsetsstrm_fd,&new_locoffsetsstrm_len,new_offsets_filename,/*randomp*/false);
#else
  new_loctable = (UINT4 *) Access_allocate_private(&loctable_access,&new_loctable_len,&seconds,new_region_filename,sizeof(UINT4));
  new_locoffsetsmeta = (UINT2 *) Access_allocate_private(&offsetsmeta_access,&new_locoffsetsmeta_len,&seconds,new_pointers_filename,sizeof(UINT2));
  new_locoffsetsstrm = (UINT2 *) Access_allocate_private(&offsetsstrm_access,&new_locoffsetsstrm_len,&seconds,new_offsets_filename,sizeof(UINT2));
#endif

  for (position = 0, regioni = 0; position < genomelength; position += localspace, regioni++) {
    old_region_strm_start = localdb->loctable[regioni * LOCTABLE_SIZE]; /* In terms of 128-bit registers */
    old_region_positions_start = localdb->loctable[regioni * LOCTABLE_SIZE + 1];
    old_region_positions_end = localdb->loctable[(regioni + 1) * LOCTABLE_SIZE + 1];
    preunique_totalcounts = old_region_positions_end - old_region_positions_start;

    oldoffsetsmeta = &(localdb->locoffsetsmeta[regioni * localdb->regionsize * DIFFERENTIAL_METAINFO_SIZE]);
    oldoffsetsstrm = &(localdb->locoffsetsstrm[old_region_strm_start * 8]); /* 8 shorts per 128-bit register */

    new_region_strm_start = new_loctable[regioni * LOCTABLE_SIZE]; /* In terms of 128-bit registers */
    newoffsetsmeta = &(new_locoffsetsmeta[regioni * localdb->regionsize * DIFFERENTIAL_METAINFO_SIZE]);
    newoffsetsstrm = &(new_locoffsetsstrm[new_region_strm_start * 8]); /* 8 shorts per 128-bit register */
#if 0
    ascending += newoffsetsmeta[(localdb->regionsize - 1) * DIFFERENTIAL_METAINFO_SIZE]; /* In terms of 128-bit registers */
#endif

    oldpositions2 = &(localdb->locpositions[old_region_positions_start]);

    if (preunique_totalcounts == 0) {
      /* Could be a region with all N's */
      /* fprintf(stderr,"Something is wrong with the offsets for region %d.  Total counts is zero.\n",regioni); */
    } else {
      positions2 = (UINT2 *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(UINT2));

      counter = (UINT4 *) CALLOC(localspace,sizeof(UINT4));
      
      /* Rearrange offsets */
      for (oligoi = 0; oligoi < localspace; oligoi += BLOCKSIZE) {
	Epu16_bitpack64_block_offsets(oldoffsets,oligoi,oldoffsetsmeta,oldoffsetsstrm);
	for (ii = 0; ii < BLOCKSIZE; ii++) {
	  block_start = oldoffsets[ii];
	  block_end = oldoffsets[ii+1];
	  assert(block_end >= block_start);
	  if (block_end > block_start) {
	    reduced = Cmet_reduce_ct_local(oligoi + ii) & mask;
	    offset = Epu16_bitpack64_read_one(reduced,newoffsetsmeta,newoffsetsstrm) + counter[reduced];
	    for (j = block_start; j < block_end; j++) {
	      debug(nt1 = shortoligo_nt(oligoi,local1part);
		    nt2 = shortoligo_nt(reduced,local1part);
		    printf("Oligo %s => %s: copying position %u to location %u\n",
			   nt1,nt2,oldpositions2[j],pointers[oligoi]);
		    FREE(nt2);
		    FREE(nt1);
		    );
	      assert(offset < preunique_totalcounts);
	      positions2[offset] = oldpositions2[j];
	      offset++;
	    }
	    counter[reduced] += block_end - block_start;
	  }
	}
      }
      
      FREE(counter);

      /* Sort positions */
      for (oligoi = 0; oligoi < localspace; oligoi += BLOCKSIZE) {
	Epu16_bitpack64_block_offsets(newoffsets,oligoi,newoffsetsmeta,newoffsetsstrm);
	for (ii = 0; ii < BLOCKSIZE; ii++) {
	  block_start = newoffsets[ii];
	  block_end = newoffsets[ii+1];
	  if ((npositions = block_end - block_start) > 0) {
	    qsort(&(positions2[block_start]),npositions,sizeof(UINT2),UINT2_compare);
	    FWRITE_USHORTS(&(positions2[block_start]),npositions,positions_fp);
	  }
	}
      }
      
      FREE(positions2);
    }
  }
  fclose(positions_fp);


  /* Clean up */
#ifdef HAVE_MMAP
  munmap((void *) new_loctable,new_loctable_len);
  munmap((void *) new_locoffsetsmeta,new_locoffsetsmeta_len);
  munmap((void *) new_locoffsetsstrm,new_locoffsetsstrm_len);
  close(new_loctable_fd);
  close(new_locoffsetsmeta_fd);
  close(new_locoffsetsstrm_fd);
#else
  FREE(new_loctable);
  FREE(new_locoffsetsmeta);
  FREE(new_locoffsetsstrm);
#endif
  fprintf(stderr,"done\n");

  return;
}


static void
compute_ct (char *new_pointers_filename, char *new_offsets_filename,
	    char *new_positions_high_filename, char *new_positions_filename,
	    Indexdb_T indexdb, Oligospace_T oligospace, Oligospace_T mask,
	    bool coord_values_8p) {
  char *ptrsfile, *compfile;
  FILE *positions_high_fp, *positions_fp;

  unsigned char *positions8_high;
  UINT4 *positions8_low;
  UINT4 *positions4;

  Oligospace_T oligoi, reduced;
  UINT4 oldoffsets[BLOCKSIZE+1];

  UINT4 *newoffsetsmeta, *newoffsetsstrm, *countermeta, *counterstrm;
  Positionsptr_T preunique_totalcounts, block_start, block_end, j, offset;
  size_t offsetsmeta_len, offsetsstrm_len;
  int ii;
#ifdef HAVE_MMAP
  int offsetsmeta_fd, offsetsstrm_fd;
#else
  Access_T offsetsmeta_access, offsetsstrm_access;
  double seconds;
#endif
#ifdef DEBUG
  char *nt1, *nt2;
#endif


  /* offsets = compute_offsets_ct_using_array(oldoffsets,oligospace,mask); */
  index_offsets_ct_using_bitpack(new_pointers_filename,new_offsets_filename,
				 indexdb->offsetsmeta,indexdb->offsetsstrm,
				 oligospace,mask);

  preunique_totalcounts = Bitpack64_read_one(oligospace,indexdb->offsetsmeta,indexdb->offsetsstrm);
  if (preunique_totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets.  Total counts is zero.\n");
    exit(9);

  } else if (coord_values_8p == true) {
    fprintf(stderr,"Trying to allocate %u*(%d+%d) bytes of memory...",
	    preunique_totalcounts,(int) sizeof(unsigned char),(int) sizeof(UINT4));
    positions4 = (UINT4 *) NULL;
    positions8_high = (unsigned char *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(unsigned char));
    positions8_low = (UINT4 *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(UINT4));
    if (positions8_high == NULL || positions8_low == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }

  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",preunique_totalcounts,(int) sizeof(UINT4));
    positions8_high = (unsigned char *) NULL;
    positions8_low = (UINT4 *) NULL;
    positions4 = (UINT4 *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(UINT4));
    if (positions4 == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  }

  /* Point to offsets and revise */
  fprintf(stderr,"Rearranging CT offsets for indexdb...");
#ifdef HAVE_MMAP
  newoffsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,new_pointers_filename,/*randomp*/false);
  newoffsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,new_offsets_filename,/*randomp*/false);
#else
  newoffsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_access,&offsetsmeta_len,&seconds,new_pointers_filename,sizeof(UINT4));
  newoffsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_access,&offsetsstrm_len,&seconds,new_offsets_filename,sizeof(UINT4));
#endif
  countermeta = Indexdb_bitpack_counter(&counterstrm,newoffsetsmeta,index1part);

  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    Bitpack64_block_offsets(oldoffsets,oligoi,indexdb->offsetsmeta,indexdb->offsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((oligoi + ii) % MONITOR_INTERVAL == 0) {
	fprintf(stderr,".");
      }
      block_start = oldoffsets[ii];
      block_end = oldoffsets[ii+1];

      if (block_end > block_start) {
	reduced = Cmet_reduce_ct(oligoi + ii) & mask;
	offset = Bitpack64_read_one(reduced,newoffsetsmeta,newoffsetsstrm) + Bitpack64_access(reduced,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  for (j = block_start; j < block_end; j++) {
	    positions8_high[offset] = indexdb->positions_high[j];
	    positions8_low[offset] = indexdb->positions[j];
	    offset++;
	  }
	} else {
	  for (j = block_start; j < block_end; j++) {
	    positions4[offset] = indexdb->positions[j];
	    offset++;
	  }
	}
	Bitpack64_add(reduced,countermeta,counterstrm,block_end - block_start);
      }
    }
  }

  FREE(counterstrm);
  FREE(countermeta);
  fprintf(stderr,"done\n");


  fprintf(stderr,"Sorting CT positions for indexdb...");
  if (snps_root == NULL) {
    ptrsfile = compfile = (char *) NULL;
  } else {
    ptrsfile = (char *) MALLOC((strlen(new_pointers_filename) + strlen(".temp") + 1) * sizeof(char));
    sprintf(ptrsfile,"%s.temp",new_pointers_filename);
    compfile = (char *) MALLOC((strlen(new_offsets_filename) + strlen(".temp") + 1) * sizeof(char));
    sprintf(compfile,"%s.temp",new_offsets_filename);
  }

  if (new_positions_high_filename != NULL) {
    if ((positions_high_fp = FOPEN_WRITE_BINARY(new_positions_high_filename)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",new_positions_high_filename);
      exit(9);
    }
  }

  if ((positions_fp = FOPEN_WRITE_BINARY(new_positions_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_positions_filename);
    exit(9);
  }

  sort_positions(ptrsfile,compfile,positions_high_fp,positions_fp,
		 positions8_high,positions8_low,positions4,
		 newoffsetsmeta,newoffsetsstrm,oligospace,coord_values_8p);

  /* Clean up */
#ifdef HAVE_MMAP
  munmap((void *) newoffsetsmeta,offsetsmeta_len);
  munmap((void *) newoffsetsstrm,offsetsstrm_len);
  close(offsetsmeta_fd);
  close(offsetsstrm_fd);
#else
  FREE(newoffsetsmeta);
  FREE(newoffsetsstrm);
#endif

  if (coord_values_8p == true) {
    FREE(positions8_high);
    FREE(positions8_low);
    fclose(positions_high_fp);
    fclose(positions_fp);
  } else {
    FREE(positions4);
    fclose(positions_fp);
  }

  if (snps_root) {
    rename(ptrsfile,new_pointers_filename);
    rename(compfile,new_offsets_filename);
    FREE(ptrsfile);
    FREE(compfile);
  }

  fprintf(stderr,"done\n");

  return;
}

#if 0
/*                                       A  C  A  T */
static unsigned char ga_conversion[4] = {0, 1, 0, 3};
static char GA_CHARTABLE[4] = {'A','C','A','T'};
#endif

static void
compute_ga_local (char *new_region_filename, char *new_pointers_filename,
		  char *new_offsets_filename, char *new_positions_filename,
		  Localdb_T localdb, Localspace_T localspace, Localspace_T mask,
		  Univcoord_T genomelength) {
  FILE *region_fp, *ptrs_fp, *comp_fp, *positions_fp;
  UINT4 ascending = 0, total_npositions = 0;

  Univcoord_T position;
  int regioni;
  UINT4 old_region_strm_start, old_region_positions_start, old_region_positions_end;
  UINT4 new_region_strm_start;

  UINT2 *positions2, *oldpositions2;
  Localspace_T oligoi, reduced;
  UINT4 newoffsets[BLOCKSIZE+1], oldoffsets[BLOCKSIZE+1];

  UINT4 *new_loctable;
  UINT2 *new_locoffsetsmeta, *new_locoffsetsstrm;
  UINT2 *newoffsetsmeta, *newoffsetsstrm, *oldoffsetsmeta, *oldoffsetsstrm;

  UINT4 *counter, preunique_totalcounts, npositions, block_start, block_end, j, offset;
  size_t new_loctable_len, new_locoffsetsmeta_len, new_locoffsetsstrm_len;
  int ii;
#ifdef HAVE_MMAP
  int new_loctable_fd, new_locoffsetsmeta_fd, new_locoffsetsstrm_fd;
#else
  Access_T loctable_access, offsetsmeta_access, offsetsstrm_access;
  double seconds;
#endif
#ifdef DEBUG
  char *nt1, *nt2;
#endif


  if ((region_fp = FOPEN_WRITE_BINARY(new_region_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_region_filename);
    exit(9);
  }
  if ((ptrs_fp = FOPEN_WRITE_BINARY(new_pointers_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_pointers_filename);
    exit(9);
  }
  if ((comp_fp = FOPEN_WRITE_BINARY(new_offsets_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_offsets_filename);
    exit(9);
  }
  if ((positions_fp = FOPEN_WRITE_BINARY(new_positions_filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",new_positions_filename);
    exit(9);
  }


  fprintf(stderr,"Rearranging GA offsets for localdb...");
  for (position = 0, regioni = 0; position < genomelength; position += localspace, regioni++) {
    /* regioni = position/65536; */
    if (regioni % MONITOR_INTERVAL_REGION == 0) {
      fprintf(stderr,".");
    }

    FWRITE_UINT(ascending,region_fp);
    FWRITE_UINT(total_npositions,region_fp);

    oldoffsetsmeta = &(localdb->locoffsetsmeta[regioni * localdb->regionsize * DIFFERENTIAL_METAINFO_SIZE]);
    old_region_strm_start = localdb->loctable[regioni * LOCTABLE_SIZE];
    oldoffsetsstrm = &(localdb->locoffsetsstrm[old_region_strm_start * 8]); /* 8 shorts per 128-bit register */

    ascending += local_offsets_ga_using_bitpack(&npositions,ptrs_fp,comp_fp,
						oldoffsetsmeta,oldoffsetsstrm,localspace,mask);
    total_npositions += npositions;
  }

  /* Write final entries to region_fp */
  FWRITE_UINT(ascending,region_fp);
  FWRITE_UINT(total_npositions,region_fp);
  fclose(region_fp);
  fclose(comp_fp);
  fclose(ptrs_fp);
  fprintf(stderr,"done\n");


  /* Point to offsets and revise */
  fprintf(stderr,"Sorting GA positions for localdb...");
#ifdef HAVE_MMAP
  new_loctable = (UINT4 *) Access_mmap(&new_loctable_fd,&new_loctable_len,new_region_filename,/*randomp*/false);
  new_locoffsetsmeta = (UINT2 *) Access_mmap(&new_locoffsetsmeta_fd,&new_locoffsetsmeta_len,new_pointers_filename,/*randomp*/false);
  new_locoffsetsstrm = (UINT2 *) Access_mmap(&new_locoffsetsstrm_fd,&new_locoffsetsstrm_len,new_offsets_filename,/*randomp*/false);
#else
  new_loctable = (UINT4 *) Access_allocate_private(&loctable_access,&new_loctable_len,&seconds,new_region_filename,sizeof(UINT4));
  new_locoffsetsmeta = (UINT2 *) Access_allocate_private(&offsetsmeta_access,&new_locoffsetsmeta_len,&seconds,new_pointers_filename,sizeof(UINT2));
  new_locoffsetsstrm = (UINT2 *) Access_allocate_private(&offsetsstrm_access,&new_locoffsetsstrm_len,&seconds,new_offsets_filename,sizeof(UINT2));
#endif

  for (position = 0, regioni = 0; position < genomelength; position += localspace, regioni++) {
    old_region_strm_start = localdb->loctable[regioni * LOCTABLE_SIZE]; /* In terms of 128-bit registers */
    old_region_positions_start = localdb->loctable[regioni * LOCTABLE_SIZE + 1];
    old_region_positions_end = localdb->loctable[(regioni + 1) * LOCTABLE_SIZE + 1];
    preunique_totalcounts = old_region_positions_end - old_region_positions_start;

    oldoffsetsmeta = &(localdb->locoffsetsmeta[regioni * localdb->regionsize * DIFFERENTIAL_METAINFO_SIZE]);
    oldoffsetsstrm = &(localdb->locoffsetsstrm[old_region_strm_start * 8]); /* 8 shorts per 128-bit register */

    new_region_strm_start = new_loctable[regioni * LOCTABLE_SIZE]; /* In terms of 128-bit registers */
    newoffsetsmeta = &(new_locoffsetsmeta[regioni * localdb->regionsize * DIFFERENTIAL_METAINFO_SIZE]);
    newoffsetsstrm = &(new_locoffsetsstrm[new_region_strm_start * 8]); /* 8 shorts per 128-bit register */
#if 0
    ascending += newoffsetsmeta[(localdb->regionsize - 1) * DIFFERENTIAL_METAINFO_SIZE]; /* In terms of 128-bit registers */
#endif

    oldpositions2 = &(localdb->locpositions[old_region_positions_start]);

    if (preunique_totalcounts == 0) {
      /* Could be a region with all N's */
      /* fprintf(stderr,"Something is wrong with the offsets in region %d.  Total counts is zero.\n",regioni); */
    } else {
      positions2 = (UINT2 *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(UINT2));

      counter = (UINT4 *) CALLOC(localspace,sizeof(UINT4));

      /* Rearrange offsets */
      for (oligoi = 0; oligoi < localspace; oligoi += BLOCKSIZE) {
	Epu16_bitpack64_block_offsets(oldoffsets,oligoi,oldoffsetsmeta,oldoffsetsstrm);
	for (ii = 0; ii < BLOCKSIZE; ii++) {
	  block_start = oldoffsets[ii];
	  block_end = oldoffsets[ii+1];
	  assert(block_end >= block_start);
	  if (block_end > block_start) {
	    reduced = Cmet_reduce_ga_local(oligoi + ii) & mask;
	    offset = Epu16_bitpack64_read_one(reduced,newoffsetsmeta,newoffsetsstrm) + counter[reduced];
	    for (j = block_start; j < block_end; j++) {
	      debug(nt1 = shortoligo_nt(oligoi,local1part);
		    nt2 = shortoligo_nt(reduced,local1part);
		    printf("Oligo %s => %s: copying position %u to location %u\n",
			   nt1,nt2,oldpositions2[j],pointers[oligoi]);
		    FREE(nt2);
		    FREE(nt1);
		    );
	      assert(offset < preunique_totalcounts);
	      positions2[offset] = oldpositions2[j];
	      offset++;
	    }
	    counter[reduced] += block_end - block_start;
	  }
	}
      }
      
      FREE(counter);
      
      /* Sort positions */
      for (oligoi = 0; oligoi < localspace; oligoi += BLOCKSIZE) {
	Epu16_bitpack64_block_offsets(newoffsets,oligoi,newoffsetsmeta,newoffsetsstrm);
	for (ii = 0; ii < BLOCKSIZE; ii++) {
	  block_start = newoffsets[ii];
	  block_end = newoffsets[ii+1];
	  if ((npositions = block_end - block_start) > 0) {
	    qsort(&(positions2[block_start]),npositions,sizeof(UINT2),UINT2_compare);
	    FWRITE_USHORTS(&(positions2[block_start]),npositions,positions_fp);
	  }
	}
      }
      
      FREE(positions2);
    }
  }
  fclose(positions_fp);


  /* Clean up */
#ifdef HAVE_MMAP
  munmap((void *) new_loctable,new_loctable_len);
  munmap((void *) new_locoffsetsmeta,new_locoffsetsmeta_len);
  munmap((void *) new_locoffsetsstrm,new_locoffsetsstrm_len);
  close(new_loctable_fd);
  close(new_locoffsetsmeta_fd);
  close(new_locoffsetsstrm_fd);
#else
  FREE(new_loctable);
  FREE(new_locoffsetsmeta);
  FREE(new_locoffsetsstrm);
#endif


  fprintf(stderr,"done\n");

  return;
}



static void
compute_ga (char *new_pointers_filename, char *new_offsets_filename,
	    char *new_positions_high_filename, char *new_positions_filename,
	    Indexdb_T indexdb, Oligospace_T oligospace, Oligospace_T mask,
	    bool coord_values_8p) {
  char *ptrsfile, *compfile;
  FILE *positions_high_fp, *positions_fp;

  unsigned char *positions8_high;
  UINT4 *positions8_low;
  UINT4 *positions4;

  Oligospace_T oligoi, reduced;
  UINT4 oldoffsets[BLOCKSIZE+1];

  UINT4 *newoffsetsmeta, *newoffsetsstrm, *countermeta, *counterstrm;
  Positionsptr_T preunique_totalcounts, block_start, block_end, j, offset;
  size_t offsetsmeta_len, offsetsstrm_len;
  int ii;
#ifdef HAVE_MMAP
  int offsetsmeta_fd, offsetsstrm_fd;
#else
  Access_T offsetsmeta_access, offsetsstrm_access;
  double seconds;
#endif
#ifdef DEBUG
  char *nt1, *nt2;
#endif


  /* offsets = compute_offsets_ga_using_array(oldoffsets,oligospace,mask); */
  index_offsets_ga_using_bitpack(new_pointers_filename,new_offsets_filename,
				 indexdb->offsetsmeta,indexdb->offsetsstrm,
				 oligospace,mask);

  preunique_totalcounts = Bitpack64_read_one(oligospace,indexdb->offsetsmeta,indexdb->offsetsstrm);
  if (preunique_totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets.  Total counts is zero.\n");
    exit(9);

  } else if (coord_values_8p == true) {
    fprintf(stderr,"Trying to allocate %u*(%d+%d) bytes of memory...",
	    preunique_totalcounts,(int) sizeof(unsigned char),(int) sizeof(UINT4));
    positions4 = (UINT4 *) NULL;
    positions8_high = (unsigned char *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(unsigned char));
    positions8_low = (UINT4 *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(UINT4));
    if (positions8_high == NULL || positions8_low == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }

  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",preunique_totalcounts,(int) sizeof(UINT4));
    positions8_high = (unsigned char *) NULL;
    positions8_low = (UINT4 *) NULL;
    positions4 = (UINT4 *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(UINT4));
    if (positions4 == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  }

  /* Point to offsets and revise */
  fprintf(stderr,"Rearranging GA offsets for indexdb...");
#ifdef HAVE_MMAP
  newoffsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,new_pointers_filename,/*randomp*/false);
  newoffsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,new_offsets_filename,/*randomp*/false);
#else
  newoffsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_access,&offsetsmeta_len,&seconds,new_pointers_filename,sizeof(UINT4));
  newoffsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_access,&offsetsstrm_len,&seconds,new_offsets_filename,sizeof(UINT4));
#endif
  countermeta = Indexdb_bitpack_counter(&counterstrm,newoffsetsmeta,index1part);

  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    Bitpack64_block_offsets(oldoffsets,oligoi,indexdb->offsetsmeta,indexdb->offsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((oligoi + ii) % MONITOR_INTERVAL == 0) {
	fprintf(stderr,".");
      }
      block_start = oldoffsets[ii];
      block_end = oldoffsets[ii+1];

      if (block_end > block_start) {
	reduced = Cmet_reduce_ga(oligoi + ii) & mask;
	offset = Bitpack64_read_one(reduced,newoffsetsmeta,newoffsetsstrm) + Bitpack64_access(reduced,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  for (j = block_start; j < block_end; j++) {
	    positions8_high[offset] = indexdb->positions_high[j];
	    positions8_low[offset] = indexdb->positions[j];
	    offset++;
	  }
	} else {
	  for (j = block_start; j < block_end; j++) {
	    positions4[offset] = indexdb->positions[j];
	    offset++;
	  }
	}
	Bitpack64_add(reduced,countermeta,counterstrm,block_end - block_start);
      }
    }
  }

  FREE(counterstrm);
  FREE(countermeta);
  fprintf(stderr,"done\n");


  fprintf(stderr,"Sorting GA positions for indexdb...");
  if (snps_root == NULL) {
    ptrsfile = compfile = (char *) NULL;
  } else {
    ptrsfile = (char *) MALLOC((strlen(new_pointers_filename) + strlen(".temp") + 1) * sizeof(char));
    sprintf(ptrsfile,"%s.temp",new_pointers_filename);
    compfile = (char *) MALLOC((strlen(new_offsets_filename) + strlen(".temp") + 1) * sizeof(char));
    sprintf(compfile,"%s.temp",new_offsets_filename);
  }

  if (new_positions_high_filename != NULL) {
    if ((positions_high_fp = FOPEN_WRITE_BINARY(new_positions_high_filename)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",new_positions_high_filename);
      exit(9);
    }
  }

  if ((positions_fp = FOPEN_WRITE_BINARY(new_positions_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_positions_filename);
    exit(9);
  }

  sort_positions(ptrsfile,compfile,positions_high_fp,positions_fp,
		 positions8_high,positions8_low,positions4,
		 newoffsetsmeta,newoffsetsstrm,oligospace,coord_values_8p);

/* Clean up */
#ifdef HAVE_MMAP
  munmap((void *) newoffsetsmeta,offsetsmeta_len);
  munmap((void *) newoffsetsstrm,offsetsstrm_len);
  close(offsetsmeta_fd);
  close(offsetsstrm_fd);
#else
  FREE(newoffsetsmeta);
  FREE(newoffsetsstrm);
#endif

  if (coord_values_8p == true) {
    FREE(positions8_high);
    FREE(positions8_low);
    fclose(positions_high_fp);
    fclose(positions_fp);
  } else {
    FREE(positions4);
    fclose(positions_fp);
  }

  if (snps_root) {
    rename(ptrsfile,new_pointers_filename);
    rename(compfile,new_offsets_filename);
    FREE(ptrsfile);
    FREE(compfile);
  }

  fprintf(stderr,"done\n");

  return;
}


/* Usage: cmetindex -d <genome> */


/* Note: Depends on having gmapindex sampled on mod 3 bounds */
int
main (int argc, char *argv[]) {
  char *sourcedir = NULL, *destdir = NULL, *filename, *fileroot;
  Indexdb_filenames_T ifilenames;
  Localdb_filenames_T lfilenames;
  char *new_region_filename, *new_pointers_filename, *new_offsets_filename,
    *new_positions_high_filename, *new_positions_filename;

  Univ_IIT_T chromosome_iit;
  Univcoord_T genomelength;
  /* size_t totalcounts, i; */
  Localspace_T localspace, local_mask;
  Oligospace_T oligospace, index_mask;
  bool coord_values_8p;

  Indexdb_T indexdb;
  Localdb_T localdb;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"F:D:d:k:q:v:",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0: 
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'cmetindex --help'",long_name);
	exit(9);
      }
      break;

    case 'F': user_sourcedir = optarg; break;
    case 'D': user_destdir = optarg; break;
    case 'd': dbroot = optarg; break;
    case 'k': required_index1part = atoi(optarg); break;
    case 'q': required_index1interval = atoi(optarg); break;
    case 'v': 
      snps_root = optarg;
      fprintf(stderr,"Combination of cmetindex and snps is not yet supported in 2018 versions\n");
      exit(9);
      break;
    default: fprintf(stderr,"Do not recognize flag %c\n",opt); exit(9);
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (dbroot == NULL) {
    fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
    fprintf(stderr,"Usage: cmetindex -d <genome>\n");
    exit(9);
  } else {
    sourcedir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_sourcedir,dbroot);
    fprintf(stderr,"Reading source files from %s\n",sourcedir);
  }

  if (user_destdir == NULL) {
    destdir = sourcedir;
  } else {
    destdir = user_destdir;
  }
  fprintf(stderr,"Writing cmet index files to %s\n",destdir);


  /* Chromosome IIT file.  Need to determine coord_values_8p and genomelength */
  filename = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(filename,"%s/%s.chromosome.iit",sourcedir,fileroot);
  if ((chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"IIT file %s is not valid\n",filename);
    exit(9);
  } else {
    coord_values_8p = Univ_IIT_coord_values_8p(chromosome_iit);
  }
  genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
  Univ_IIT_free(&chromosome_iit);
  FREE(filename);


  lfilenames = Localdb_get_filenames(&local1part,&local1interval,
				     sourcedir,fileroot,IDX_FILESUFFIX,snps_root,
				     required_local1part,required_local1interval,
				     /*offsets_only_p*/false);
  if ((localdb = Localdb_new_genome(&local1part,&local1interval,
				    sourcedir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
				    required_local1part,required_local1interval,
				    /*locoffsetsstrm_access*/USE_ALLOCATE,
				    /*locpositions_access*/USE_ALLOCATE,/*sharedp*/false,
				    /*multiple_sequences_p*/false,/*unload_shared_memory_p*/false)) != NULL) {

    local_mask = ~(~0U << 2*local1part);
    localspace = (Localspace_T) power(4,local1part);

    new_region_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen("metct")+strlen(lfilenames->loctable_local1info_ptr)+1,sizeof(char));
    sprintf(new_region_filename,"%s/%s.%s%s",destdir,fileroot,"metct",lfilenames->loctable_local1info_ptr);

    new_pointers_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					    strlen(".")+strlen("metct")+strlen(lfilenames->locpointers_local1info_ptr)+1,sizeof(char));
    sprintf(new_pointers_filename,"%s/%s.%s%s",destdir,fileroot,"metct",lfilenames->locpointers_local1info_ptr);
    
    new_offsets_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen("metct")+strlen(lfilenames->locoffsets_local1info_ptr)+1,sizeof(char));
    sprintf(new_offsets_filename,"%s/%s.%s%s",destdir,fileroot,"metct",lfilenames->locoffsets_local1info_ptr);

    new_positions_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					     strlen(".")+strlen("metct")+strlen(lfilenames->locpositions_local1info_ptr)+1,sizeof(char));
    sprintf(new_positions_filename,"%s/%s.%s%s",destdir,fileroot,"metct",lfilenames->locpositions_local1info_ptr);

    compute_ct_local(new_region_filename,new_pointers_filename,new_offsets_filename,new_positions_filename,
		     localdb,localspace,local_mask,genomelength);

    FREE(new_region_filename);
    FREE(new_positions_filename);
    FREE(new_offsets_filename);
    FREE(new_pointers_filename);


    new_region_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen("metga")+strlen(lfilenames->loctable_local1info_ptr)+1,sizeof(char));
    sprintf(new_region_filename,"%s/%s.%s%s",destdir,fileroot,"metga",lfilenames->loctable_local1info_ptr);

    new_pointers_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					    strlen(".")+strlen("metga")+strlen(lfilenames->locpointers_local1info_ptr)+1,sizeof(char));
    sprintf(new_pointers_filename,"%s/%s.%s%s",destdir,fileroot,"metga",lfilenames->locpointers_local1info_ptr);
    
    new_offsets_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen("metga")+strlen(lfilenames->locoffsets_local1info_ptr)+1,sizeof(char));
    sprintf(new_offsets_filename,"%s/%s.%s%s",destdir,fileroot,"metga",lfilenames->locoffsets_local1info_ptr);

    new_positions_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					     strlen(".")+strlen("metga")+strlen(lfilenames->locpositions_local1info_ptr)+1,sizeof(char));
    sprintf(new_positions_filename,"%s/%s.%s%s",destdir,fileroot,"metga",lfilenames->locpositions_local1info_ptr);

    compute_ga_local(new_region_filename,new_pointers_filename,new_offsets_filename,new_positions_filename,
		     localdb,localspace,local_mask,genomelength);

    FREE(new_region_filename);
    FREE(new_positions_filename);
    FREE(new_offsets_filename);
    FREE(new_pointers_filename);

    Localdb_free(&localdb);
    Localdb_filenames_free(&lfilenames);
  }


  ifilenames = Indexdb_get_filenames(&compression_type,&index1part,&index1interval,
				     sourcedir,fileroot,IDX_FILESUFFIX,snps_root,
				     required_index1part,required_index1interval,
				     /*offsets_only_p*/false);
  if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
				    sourcedir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
				    required_index1part,required_index1interval,
				    /*expand_offsets_p*/false,/*offsetsstrm_access*/USE_ALLOCATE,
				    /*positions_access*/USE_ALLOCATE,/*sharedp*/false,
				    /*multiple_sequences_p*/false,/*preload_shared_memory_p*/false,/*unload_shared_memory_p*/false)) != NULL) {
#ifdef HAVE_64_BIT
    index_mask = ~(~0ULL << 2*index1part);
#else
    index_mask = ~(~0U << 2*index1part);
#endif
    oligospace = power(4,index1part);


    /* Compute and write CT files */
    new_pointers_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					    strlen(".")+strlen("metct")+strlen(ifilenames->pointers_index1info_ptr)+1,sizeof(char));
    sprintf(new_pointers_filename,"%s/%s.%s%s",destdir,fileroot,"metct",ifilenames->pointers_index1info_ptr);
    
    new_offsets_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen("metct")+strlen(ifilenames->offsets_index1info_ptr)+1,sizeof(char));
    sprintf(new_offsets_filename,"%s/%s.%s%s",destdir,fileroot,"metct",ifilenames->offsets_index1info_ptr);
    
    new_positions_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					     strlen(".")+strlen("metct")+strlen(ifilenames->positions_index1info_ptr)+1,sizeof(char));
    sprintf(new_positions_filename,"%s/%s.%s%s",destdir,fileroot,"metct",ifilenames->positions_index1info_ptr);
    
    if (coord_values_8p == true) {
      new_positions_high_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
						    strlen(".")+strlen("metct")+strlen(ifilenames->positions_high_index1info_ptr)+1,sizeof(char));
      sprintf(new_positions_high_filename,"%s/%s.%s%s",destdir,fileroot,"metct",ifilenames->positions_high_index1info_ptr);
    }

    compute_ct(new_pointers_filename,new_offsets_filename,
	       new_positions_high_filename,new_positions_filename,
	       indexdb,oligospace,index_mask,coord_values_8p);
    
    if (coord_values_8p == true) {
      FREE(new_positions_high_filename);
    }
    FREE(new_positions_filename);
    FREE(new_offsets_filename);
    FREE(new_pointers_filename);


    /* Compute and write GA files */
    new_pointers_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					    strlen(".")+strlen("metga")+strlen(ifilenames->pointers_index1info_ptr)+1,sizeof(char));
    sprintf(new_pointers_filename,"%s/%s.%s%s",destdir,fileroot,"metga",ifilenames->pointers_index1info_ptr);
    
    new_offsets_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen("metga")+strlen(ifilenames->offsets_index1info_ptr)+1,sizeof(char));
    sprintf(new_offsets_filename,"%s/%s.%s%s",destdir,fileroot,"metga",ifilenames->offsets_index1info_ptr);
    
    new_positions_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					     strlen(".")+strlen("metga")+strlen(ifilenames->positions_index1info_ptr)+1,sizeof(char));
    sprintf(new_positions_filename,"%s/%s.%s%s",destdir,fileroot,"metga",ifilenames->positions_index1info_ptr);

    if (coord_values_8p == true) {
      new_positions_high_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
						    strlen(".")+strlen("metga")+strlen(ifilenames->positions_high_index1info_ptr)+1,sizeof(char));
      sprintf(new_positions_high_filename,"%s/%s.%s%s",destdir,fileroot,"metga",ifilenames->positions_high_index1info_ptr);
    }

    compute_ga(new_pointers_filename,new_offsets_filename,
	       new_positions_high_filename,new_positions_filename,
	       indexdb,oligospace,index_mask,coord_values_8p);
    if (coord_values_8p == true) {
      FREE(new_positions_high_filename);
    }
    FREE(new_positions_filename);
    FREE(new_offsets_filename);
    FREE(new_pointers_filename);

    Indexdb_free(&indexdb);
    Indexdb_filenames_free(&ifilenames);
  }

  FREE(dbversion);
  FREE(fileroot);
  FREE(sourcedir);

  return 0;
}



static void
print_program_usage () {
  fprintf(stdout,"\
Usage: cmetindex [OPTIONS...] -d <genome>\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Options (must include -d)\n");
  fprintf(stdout,"\
  -F, --sourcedir=directory      Directory where to read cmet index files (default is\n\
                                   GMAP genome directory specified at compile time)\n\
  -D, --destdir=directory        Directory where to write cmet index files (default is\n\
                                   value of -F, if provided; otherwise the value of the\n\
                                   GMAP genome directory specified at compile time)\n\
  -d, --db=STRING                Genome database\n\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less).\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  -q, --sampling=INT             Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected basesize and k-mer size\n\
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

