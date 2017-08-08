static char rcsid[] = "$Id: cmetindex.c 184156 2016-02-12 18:30:44Z twu $";
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

#include "cmet.h"

#include "bool.h"
#include "genomicpos.h"
#include "iit-read-univ.h"
#include "indexdb.h"
#include "indexdb-write.h"
#include "genome.h"
#include "genome128_hr.h"
#include "bytecoding.h"
#include "sarray-write.h"
#include "bitpack64-write.h"
#include "bitpack64-read.h"
#include "bitpack64-readtwo.h"
#include "bitpack64-access.h"
#include "bitpack64-incr.h"
#include "datadir.h"
#include "getopt.h"


#define POSITIONS8_HIGH_SHIFT 32
#define POSITIONS8_LOW_MASK 0xFFFFFFFF

#define DIFFERENTIAL_METAINFO_SIZE 2
#define BLOCKSIZE 64
#define BUFFER_SIZE 1000000

#define MONITOR_INTERVAL 10000000


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
static int required_interval = 0;

static bool build_suffix_array_p = true;
static char *snps_root = NULL;


static struct option long_options[] = {
  /* Input options */
  {"sourcedir", required_argument, 0, 'F'},	/* user_sourcedir */
  {"destdir", required_argument, 0, 'D'},	/* user_destdir */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part */
  {"sampling", required_argument, 0, 'q'}, /* required_interval */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"usesnps", required_argument, 0, 'v'}, /* snps_root */
  {"build-sarray", required_argument, 0, 0}, /* build_suffix_array_p */

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


static void
compute_offsets_ct_using_bitpack (char *new_pointers_filename, char *new_offsets_filename,
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


static void
compute_offsets_ga_using_bitpack (char *new_pointers_filename, char *new_offsets_filename,
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
sort_positions (char *ptrsfile, char *compfile, FILE *positions_high_fp, FILE *positions_low_fp,
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
	    FWRITE_UINTS(&(positions8_low[block_start]),npositions,positions_low_fp);
	  } else {
	    /* FWRITE_UINT8(positions8[block_start],positions_fp); */
	    FWRITE_CHAR(positions8_high[block_start],positions_high_fp);
	    FWRITE_UINT(positions8_low[block_start],positions_low_fp);
	    for (j = block_start+1; j < block_end; j++) {
	      if (positions8_high[j] == positions8_high[j-1] && positions8_low[j] == positions8_low[j-1]) {
		npositions--;
	      } else {
		/* FWRITE_UINT8(positions8[j],positions_fp); */
		FWRITE_CHAR(positions8_high[j],positions_high_fp);
		FWRITE_UINT(positions8_low[j],positions_low_fp);
	      }
	    }
	  }

	} else {
	  qsort(&(positions4[block_start]),npositions,sizeof(UINT4),UINT4_compare);
	  if (snps_root == NULL) {
	    FWRITE_UINTS(&(positions4[block_start]),npositions,positions_low_fp);
	  } else {
	    FWRITE_UINT(positions4[block_start],positions_low_fp);
	    for (j = block_start+1; j < block_end; j++) {
	      if (positions4[j] == positions4[j-1]) {
		npositions--;
	      } else {
		FWRITE_UINT(positions4[j],positions_low_fp);
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
	    FWRITE_UINTS(&(positions8_low[block_start]),npositions,positions_low_fp);
	  } else {
	    /* FWRITE_UINT8(positions8[block_start],positions_fp); */
	    FWRITE_CHAR(positions8_high[block_start],positions_high_fp);
	    FWRITE_UINT(positions8_low[block_start],positions_low_fp);
	    for (j = block_start+1; j < block_end; j++) {
	      if (positions8_high[j] == positions8_high[j-1] && positions8_low[j] == positions8_low[j-1]) {
		npositions--;
	      } else {
		/* FWRITE_UINT8(positions8[j],positions_fp); */
		FWRITE_CHAR(positions8_high[j],positions_high_fp);
		FWRITE_UINT(positions8_low[j],positions_low_fp);
	      }
	    }
	  }

	} else {
	  qsort(&(positions4[block_start]),npositions,sizeof(UINT4),UINT4_compare);
	  if (snps_root == NULL) {
	    FWRITE_UINTS(&(positions4[block_start]),npositions,positions_low_fp);
	  } else {
	    FWRITE_UINT(positions4[block_start],positions_low_fp);
	    for (j = block_start+1; j < block_end; j++) {
	      if (positions4[j] == positions4[j-1]) {
		npositions--;
	      } else {
		FWRITE_UINT(positions4[j],positions_low_fp);
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


/*                                       A  T  G  T */
static unsigned char ct_conversion[4] = {0, 3, 2, 3};
static char CT_CHARTABLE[4] = {'A','T','G','T'};


static void
compute_ct (char *new_pointers_filename, char *new_offsets_filename,
	    FILE *positions_high_fp, FILE *positions_low_fp,
	    UINT4 *oldoffsetsmeta, UINT4 *oldoffsetsstrm, UINT8 *oldpositions8, UINT4 *oldpositions4,
	    Oligospace_T oligospace, Oligospace_T mask, bool coord_values_8p) {
  char *ptrsfile, *compfile;
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
  double seconds;
#endif
#ifdef DEBUG
  char *nt1, *nt2;
#endif


  /* offsets = compute_offsets_ct_using_array(oldoffsets,oligospace,mask); */
  compute_offsets_ct_using_bitpack(new_pointers_filename,new_offsets_filename,oldoffsetsmeta,oldoffsetsstrm,
				   oligospace,mask);

  preunique_totalcounts = Bitpack64_read_one(oligospace,oldoffsetsmeta,oldoffsetsstrm);
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
  fprintf(stderr,"Rearranging CT positions...");
#ifdef HAVE_MMAP
  newoffsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,new_pointers_filename,/*randomp*/false);
  newoffsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,new_offsets_filename,/*randomp*/false);
#else
  newoffsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_access,&offsetsmeta_len,&seconds,new_pointers_filename,sizeof(UINT4));
  newoffsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_access,&offsetsstrm_len,&seconds,new_offsets_filename,sizeof(UINT4));
#endif
  countermeta = Indexdb_bitpack_counter(&counterstrm,newoffsetsmeta,newoffsetsstrm,index1part);

  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    Bitpack64_block_offsets(oldoffsets,oligoi,oldoffsetsmeta,oldoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((oligoi + ii) % MONITOR_INTERVAL == 0) {
	fprintf(stderr,".");
      }
      block_start = oldoffsets[ii];
      block_end = oldoffsets[ii+1];
#ifdef WORDS_BIGENDIAN
      reduced = Cmet_reduce_ct(oligoi + ii) & mask;
      offset = Bitpack64_read_one(reduced,newoffsetsmeta,newoffsetsstrm) + Bitpack64_access(reduced,countermeta,counterstrm);
      if (coord_values_8p == true) {
	for (j = Bigendian_convert_uint(block_start); j < Bigendian_convert_uint(block_end); j++) {
	  debug(nt1 = shortoligo_nt(oligoi,index1part);
		nt2 = shortoligo_nt(reduced,index1part);
		printf("Oligo %s => %s: copying position %u to location %u\n",
		       nt1,nt2,oldpositions8[j],pointers[oligoi]);
		FREE(nt2);
		FREE(nt1);
		);
	  positions8_high[pointers[reduced]/*++*/] = Bigendian_convert_uint8(oldpositions8[j]) >> POSITIONS8_HIGH_SHIFT;
	  positions8_low[pointers[reduced]++] = Bigendian_convert_uint8(oldpositions8[j]) & POSITIONS8_LOW_MASK;
	}
      } else {
	for (j = Bigendian_convert_uint(oldoffsets[oligoi]); j < Bigendian_convert_uint(oldoffsets[oligoi+1]); j++) {
	  debug(nt1 = shortoligo_nt(oligoi,index1part);
		nt2 = shortoligo_nt(reduced,index1part);
		printf("Oligo %s => %s: copying position %u to location %u\n",
		       nt1,nt2,oldpositions4[j],pointers[oligoi]);
		FREE(nt2);
		FREE(nt1);
		);
	  positions4[pointers[reduced]++] = Bigendian_convert_uint(oldpositions4[j]);
	}
      }
#else
      if (block_end > block_start) {
	reduced = Cmet_reduce_ct(oligoi + ii) & mask;
	offset = Bitpack64_read_one(reduced,newoffsetsmeta,newoffsetsstrm) + Bitpack64_access(reduced,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  for (j = block_start; j < block_end; j++) {
	    debug(nt1 = shortoligo_nt(oligoi,index1part);
		  nt2 = shortoligo_nt(reduced,index1part);
		  printf("Oligo %s => %s: copying position %u to location %u\n",
			 nt1,nt2,oldpositions8[j],pointers[oligoi]);
		  FREE(nt2);
		  FREE(nt1);
		  );
	    positions8_high[offset] = oldpositions8[j] >> POSITIONS8_HIGH_SHIFT;
	    positions8_low[offset] = oldpositions8[j] & POSITIONS8_LOW_MASK;
	    offset++;
	  }
	} else {
	  for (j = block_start; j < block_end; j++) {
	    debug(nt1 = shortoligo_nt(oligoi,index1part);
		  nt2 = shortoligo_nt(reduced,index1part);
		  printf("Oligo %s => %s: copying position %u to location %u\n",
			 nt1,nt2,oldpositions4[j],pointers[oligoi]);
		  FREE(nt2);
		  FREE(nt1);
		  );
	    positions4[offset] = oldpositions4[j];
	    offset++;
	  }
	}
	Bitpack64_add(reduced,countermeta,counterstrm,block_end - block_start);
      }
    }
#endif
  }

  FREE(counterstrm);
  FREE(countermeta);
  fprintf(stderr,"done\n");


  fprintf(stderr,"Sorting CT positions...");
  if (snps_root == NULL) {
    ptrsfile = compfile = (char *) NULL;
  } else {
    ptrsfile = (char *) MALLOC((strlen(new_pointers_filename) + strlen(".temp") + 1) * sizeof(char));
    sprintf(ptrsfile,"%s.temp",new_pointers_filename);
    compfile = (char *) MALLOC((strlen(new_offsets_filename) + strlen(".temp") + 1) * sizeof(char));
    sprintf(compfile,"%s.temp",new_offsets_filename);
  }

  sort_positions(ptrsfile,compfile,positions_high_fp,positions_low_fp,
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
  } else {
    FREE(positions4);
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

/*                                       A  C  A  T */
static unsigned char ga_conversion[4] = {0, 1, 0, 3};
static char GA_CHARTABLE[4] = {'A','C','A','T'};

static void
compute_ga (char *new_pointers_filename, char *new_offsets_filename,
	    FILE *positions_high_fp, FILE *positions_low_fp,
	    UINT4 *oldoffsetsmeta, UINT4 *oldoffsetsstrm, UINT8 *oldpositions8, UINT4 *oldpositions4,
	    Oligospace_T oligospace, Oligospace_T mask, bool coord_values_8p) {
  char *ptrsfile, *compfile;
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
  double seconds;
#endif
#ifdef DEBUG
  char *nt1, *nt2;
#endif


  /* offsets = compute_offsets_ga_using_array(oldoffsets,oligospace,mask); */
  compute_offsets_ga_using_bitpack(new_pointers_filename,new_offsets_filename,oldoffsetsmeta,oldoffsetsstrm,
				   oligospace,mask);

  preunique_totalcounts = Bitpack64_read_one(oligospace,oldoffsetsmeta,oldoffsetsstrm);
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
  fprintf(stderr,"Rearranging GA positions...");
#ifdef HAVE_MMAP
  newoffsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,new_pointers_filename,/*randomp*/false);
  newoffsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,new_offsets_filename,/*randomp*/false);
#else
  newoffsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_access,&offsetsmeta_len,&seconds,new_pointers_filename,sizeof(UINT4));
  newoffsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_access,&offsetsstrm_len,&seconds,new_offsets_filename,sizeof(UINT4));
#endif
  countermeta = Indexdb_bitpack_counter(&counterstrm,newoffsetsmeta,newoffsetsstrm,index1part);

  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    Bitpack64_block_offsets(oldoffsets,oligoi,oldoffsetsmeta,oldoffsetsstrm);
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      if ((oligoi + ii) % MONITOR_INTERVAL == 0) {
	fprintf(stderr,".");
      }
      block_start = oldoffsets[ii];
      block_end = oldoffsets[ii+1];
#ifdef WORDS_BIGENDIAN
      reduced = Cmet_reduce_ga(oligoi + ii) & mask;
      offset = Bitpack64_read_one(reduced,newoffsetsmeta,newoffsetsstrm) + Bitpack64_access(reduced,countermeta,counterstrm);
      if (coord_values_8p == true) {
	for (j = Bigendian_convert_uint(block_start); j < Bigendian_convert_uint(block_end); j++) {
	  debug(nt1 = shortoligo_nt(oligoi,index1part);
		nt2 = shortoligo_nt(reduced,index1part);
		printf("Oligo %s => %s: copying position %u to location %u\n",
		       nt1,nt2,oldpositions8[j],pointers[oligoi]);
		FREE(nt2);
		FREE(nt1);
		);
	  positions8_high[pointers[reduced]/*++*/] = Bigendian_convert_uint8(oldpositions8[j]) >> POSITIONS8_HIGH_SHIFT;
	  positions8_low[pointers[reduced]++] = Bigendian_convert_uint8(oldpositions8[j]) & POSITIONS8_LOW_MASK;
	}
      } else {
	for (j = Bigendian_convert_uint(oldoffsets[oligoi]); j < Bigendian_convert_uint(oldoffsets[oligoi+1]); j++) {
	  debug(nt1 = shortoligo_nt(oligoi,index1part);
		nt2 = shortoligo_nt(reduced,index1part);
		printf("Oligo %s => %s: copying position %u to location %u\n",
		       nt1,nt2,oldpositions4[j],pointers[oligoi]);
		FREE(nt2);
		FREE(nt1);
		);
	  positions4[pointers[reduced]++] = Bigendian_convert_uint(oldpositions4[j]);
	}
      }
#else
      if (block_end > block_start) {
	reduced = Cmet_reduce_ga(oligoi + ii) & mask;
	offset = Bitpack64_read_one(reduced,newoffsetsmeta,newoffsetsstrm) + Bitpack64_access(reduced,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  for (j = block_start; j < block_end; j++) {
	    debug(nt1 = shortoligo_nt(oligoi,index1part);
		  nt2 = shortoligo_nt(reduced,index1part);
		  printf("Oligo %s => %s: copying position %u to location %u\n",
			 nt1,nt2,oldpositions8[j],pointers[oligoi]);
		  FREE(nt2);
		  FREE(nt1);
		  );
	    positions8_high[offset] = oldpositions8[j] >> POSITIONS8_HIGH_SHIFT;
	    positions8_low[offset] = oldpositions8[j] & POSITIONS8_LOW_MASK;
	    offset++;
	  }
	} else {
	  for (j = block_start; j < block_end; j++) {
	    debug(nt1 = shortoligo_nt(oligoi,index1part);
		  nt2 = shortoligo_nt(reduced,index1part);
		  printf("Oligo %s => %s: copying position %u to location %u\n",
			 nt1,nt2,oldpositions4[j],pointers[oligoi]);
		  FREE(nt2);
		  FREE(nt1);
		  );
	    positions4[offset] = oldpositions4[j];
	    offset++;
	  }
	}
	Bitpack64_add(reduced,countermeta,counterstrm,block_end - block_start);
      }
    }
#endif
  }

  FREE(counterstrm);
  FREE(countermeta);
  fprintf(stderr,"done\n");


  fprintf(stderr,"Sorting GA positions...");
  if (snps_root == NULL) {
    ptrsfile = compfile = (char *) NULL;
  } else {
    ptrsfile = (char *) MALLOC((strlen(new_pointers_filename) + strlen(".temp") + 1) * sizeof(char));
    sprintf(ptrsfile,"%s.temp",new_pointers_filename);
    compfile = (char *) MALLOC((strlen(new_offsets_filename) + strlen(".temp") + 1) * sizeof(char));
    sprintf(compfile,"%s.temp",new_offsets_filename);
  }

  sort_positions(ptrsfile,compfile,positions_high_fp,positions_low_fp,
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
  } else {
    FREE(positions4);
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
  Filenames_T filenames;
  char *new_pointers_filename, *new_offsets_filename;
  Univ_IIT_T chromosome_iit;
  UINT4 *ref_offsetsmeta, *ref_offsetsstrm;
  size_t totalcounts, i;
  Oligospace_T mask;
  unsigned char *ref_positions8_high;
  UINT4 *ref_positions8_low;
  UINT8 *ref_positions8;
  UINT4 *ref_positions4;
  Oligospace_T oligospace;
  bool coord_values_8p;

  /* For suffix array */
  Univcoord_T genomelength;
  char *sarrayfile, *lcpexcfile, *lcpguidefile;
  char *childexcfile, *childguidefile;
  char *lcpchilddcfile;
  char *indexijptrsfile, *indexijcompfile;
  Genome_T genomecomp;
  unsigned char *gbuffer;
  UINT4 *SA;
  UINT4 nbytes;

  unsigned char *discrim_chars;
  unsigned char *lcp_bytes;
  UINT4 *lcp_guide, *lcp_exceptions, n_lcp_exceptions;

  unsigned char *child_bytes;
  UINT4 *child_exceptions, n_child_exceptions;

  int sa_fd;
  Access_T lcpguide_access;
  size_t sa_len, lcpguide_len;
  double seconds;


  int offsetsmeta_fd, offsetsstrm_fd;
  size_t offsetsmeta_len, offsetsstrm_len;

  FILE *positions_high_fp, *positions_low_fp;
  int ref_positions_high_fd, ref_positions_low_fd;
  size_t ref_positions_high_len, ref_positions_low_len;
#ifndef HAVE_MMAP
  double seconds;
#endif

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

      } else if (!strcmp(long_name,"build-sarray")) {
	if (!strcmp(optarg,"0")) {
	  build_suffix_array_p = false;
	} else if (!strcmp(optarg,"1")) {
	  build_suffix_array_p = true;
	} else {
	  fprintf(stderr,"Argument to --build-sarray must be 0 or 1\n");
	  exit(9);
	}

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
    case 'q': required_interval = atoi(optarg); break;
    case 'v': snps_root = optarg; break;
    default: fprintf(stderr,"Do not recognize flag %c\n",opt); exit(9);
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

#ifdef MEMUSAGE
  Mem_usage_init();
#endif

  if (dbroot == NULL) {
    fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
    fprintf(stderr,"Usage: cmetindex -d <genome>\n");
    exit(9);
  } else {
    sourcedir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_sourcedir,dbroot);
    fprintf(stderr,"Reading source files from %s\n",sourcedir);
  }

  /* Chromosome IIT file.  Need to determine coord_values_8p */
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


  filenames = Indexdb_get_filenames(&compression_type,&index1part,&index1interval,
				    sourcedir,fileroot,IDX_FILESUFFIX,snps_root,
				    required_index1part,required_interval,
				    /*offsets_only_p*/false);

#ifdef HAVE_64_BIT
  mask = ~(~0ULL << 2*index1part);
#else
  mask = ~(~0U << 2*index1part);
#endif
  oligospace = power(4,index1part);

  /* Read offsets */
#ifdef MEMUSAGE
  printf("%lu\n",Mem_usage_report_std_heap()); /* 421 */
#endif
#ifdef HAVE_MMAP
  ref_offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,filenames->pointers_filename,/*randomp*/false);
  ref_offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,filenames->offsets_filename,/*randomp*/false);
#else
  ref_offsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_access,&offsetsmeta_len,&seconds,filenames->pointers_filename,sizeof(UINT4));
  ref_offsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_access,&offsetsstrm_len,&seconds,filenames->offsets_filename,sizeof(UINT4));
#endif
#ifdef MEMUSAGE
  printf("%lu\n",Mem_usage_report_std_heap()); /* 4,294,967,721 */
#endif

  /* Read positions */
#if 0
  if ((ref_positions_fp = FOPEN_READ_BINARY(filenames->positions_low_filename)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",filenames->positions_low_filename);
    exit(9);
  } else {
    fclose(ref_positions_fp);
  }
#endif

  if (coord_values_8p == true) {
#ifdef HAVE_MMAP
    ref_positions8_high = (unsigned char *) Access_mmap(&ref_positions_high_fd,&ref_positions_high_len,
							filenames->positions_high_filename,/*randomp*/false);
    ref_positions8_low = (UINT4 *) Access_mmap(&ref_positions_low_fd,&ref_positions_low_len,
					       filenames->positions_low_filename,/*randomp*/false);
#else
    ref_positions8_high = (unsigned char *) Access_allocate_private(&ref_positions_high_access,&ref_positions_high_len,&seconds,
								    filenames->positions_high_filename,sizeof(unsigned char));
    ref_positions8_low = (UINT4 *) Access_allocate_private(&ref_positions_low_access,&ref_positions_low_len,&seconds,
							   filenames->positions_low_filename,sizeof(UINT4));
#endif
    /* Unpack */
    totalcounts = ref_positions_high_len/sizeof(unsigned char);
    if (totalcounts > 4294967295) {
      fprintf(stderr,"Program not yet designed to handle huge genomes\n");
      abort();
    }
    ref_positions8 = (UINT8 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT8));
    for (i = 0; i < totalcounts; i++) {
      ref_positions8[i] = ((UINT8) ref_positions8_high[i] << 32) + ref_positions8_low[i];
    }
#ifdef HAVE_MMAP
    munmap((void *) ref_positions8_high,ref_positions_high_len);
    munmap((void *) ref_positions8_low,ref_positions_low_len);
    close(ref_positions_high_fd);
    close(ref_positions_low_fd);
#else
    FREE(ref_positions8_high);
    FREE(ref_positions8_low);
#endif

  } else {
#ifdef HAVE_MMAP
    ref_positions4 = (UINT4 *) Access_mmap(&ref_positions_low_fd,&ref_positions_low_len,
					   filenames->positions_low_filename,/*randomp*/false);
#else
    ref_positions4 = (UINT4 *) Access_allocate_private(&ref_positions_low_access,&ref_positions_low_len,&seconds,
						       filenames->positions_low_filename,sizeof(UINT4));
#endif
  }


  /* Open CT output files */
  if (user_destdir == NULL) {
    destdir = sourcedir;
  } else {
    destdir = user_destdir;
  }
  fprintf(stderr,"Writing cmet index files to %s\n",destdir);

  new_pointers_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					  strlen(".")+strlen("metct")+strlen(filenames->pointers_index1info_ptr)+1,sizeof(char));
  sprintf(new_pointers_filename,"%s/%s.%s%s",destdir,fileroot,"metct",filenames->pointers_index1info_ptr);

  new_offsets_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("metct")+strlen(filenames->offsets_index1info_ptr)+1,sizeof(char));
  sprintf(new_offsets_filename,"%s/%s.%s%s",destdir,fileroot,"metct",filenames->offsets_index1info_ptr);


  if (coord_values_8p == true) {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metct")+strlen(filenames->positions_high_index1info_ptr)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s",destdir,fileroot,"metct",filenames->positions_high_index1info_ptr);
    
    if ((positions_high_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
      fprintf(stderr,"Can't open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);
  }

  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("metct")+strlen(filenames->positions_low_index1info_ptr)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%s",destdir,fileroot,"metct",filenames->positions_low_index1info_ptr);

  if ((positions_low_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);


  /* Compute and write CT files */
  if (build_suffix_array_p == true) {
#ifdef MEMUSAGE
    printf("Before building suffix array: %lu\n",Mem_usage_report_std_heap()); /* 4,294,967,891 */
#endif
    fprintf(stderr,"Building suffix array for CT\n");
    sarrayfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.sarray")+1,sizeof(char));
    sprintf(sarrayfile,"%s/%s.metct.sarray",destdir,fileroot);
    genomecomp = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
			    /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);

    gbuffer = (unsigned char *) CALLOC(genomelength+1,sizeof(unsigned char)); /* Adds 3 GB */
    Genome_fill_buffer_int_string(genomecomp,/*left*/0,/*length*/genomelength,gbuffer,ct_conversion);
    gbuffer[genomelength] = 0;	/* Tried N/X, but SACA_K fails */
    Sarray_write_array_from_genome(sarrayfile,gbuffer,genomelength);
#ifdef MEMUSAGE
    printf("After writing suffix array: %lu\n",Mem_usage_report_std_heap()); /* 7,390,678,634 */
#endif

    /* Bucket array */
    indexijptrsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.saindex64meta")+1,sizeof(char));
    sprintf(indexijptrsfile,"%s/%s.metct.saindex64meta",destdir,fileroot);
    indexijcompfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.saindex64strm")+1,sizeof(char));
    sprintf(indexijcompfile,"%s/%s.metct.saindex64strm",destdir,fileroot);
    Sarray_write_index_interleaved(indexijptrsfile,indexijcompfile,
				   sarrayfile,genomecomp,genomelength,/*compressp*/true,CT_CHARTABLE);
    FREE(indexijcompfile);
    FREE(indexijptrsfile);

    SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,/*randomp*/false);

#if 0
    /* Not needed if we already have gbuffer */
    /* Required for computing LCP, but uses non-SIMD instructions */
    genomebits = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_BITS,
			    /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);
    Genome_hr_setup(Genome_blocks(genomebits),/*snp_blocks*/NULL,
		    /*query_unk_mismatch_p*/false,/*genome_unk_mismatch_p*/false,
		    /*mode*/CMET_STRANDED);
#endif

#ifdef MEMUSAGE
    printf("Before lcp: %lu\n",Mem_usage_report_std_heap()); /* 7,390,678,634 */
#endif
    lcp_bytes = Sarray_compute_lcp_bytes_from_genome(&lcp_exceptions,&n_lcp_exceptions,SA,gbuffer,/*n*/genomelength);
#ifdef MEMUSAGE
    printf("After lcp: %lu\n",Mem_usage_report_std_heap());    /* lcp: 19,773,520,854.  lcp step requires 12 GB.  lcp_bytes: 18,711,482,949.  Using cells: 12,542,531,581  */
#endif

    munmap((void *) SA,sa_len);
    close(sa_fd);
    FREE(gbuffer);
#if 0
    Genome_free(&genomebits);
#endif

    /* Write lcp exceptions/guide */
    lcpexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.salcpexc")+1,sizeof(char));
    sprintf(lcpexcfile,"%s/%s.metct.salcpexc",destdir,fileroot);
    lcpguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.salcpguide1024")+1,sizeof(char));
    sprintf(lcpguidefile,"%s/%s.metct.salcpguide1024",destdir,fileroot);
    Bytecoding_write_bytes(/*bytesfile*/NULL,lcpexcfile,lcpguidefile,lcp_bytes,lcp_exceptions,n_lcp_exceptions,
			   genomelength,/*guide_interval*/1024);
#ifdef MEMUSAGE
    printf("Creating lcp_bytes: %lu\n",Mem_usage_report_std_heap()); /* 19,773,521,016.  Using cells: 9,446,821,188 */
#endif

    FREE(lcpguidefile);
    FREE(lcpexcfile);

#ifdef MEMUSAGE
    printf("After free of lcp: %lu\n",Mem_usage_report_std_heap()); /* 7,390,678,634 */
#endif

    /* DC array */
    /* Assume we have lcp_bytes and lcp_exceptions already in memory.  Don't need to use guide for speed. */
    lcpguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.salcpguide1024")+1,sizeof(char));
    sprintf(lcpguidefile,"%s/%s.metct.salcpguide1024",destdir,fileroot);
    lcp_guide = (UINT4 *) Access_allocate_private(&lcpguide_access,&lcpguide_len,&seconds,lcpguidefile,sizeof(UINT4));
    FREE(lcpguidefile);

    /* Compute discriminating chars (DC) array */
#ifdef MEMUSAGE
    printf("Before computing discrim_chars: %lu\n",Mem_usage_report_std_heap()); /* 9,458,913,650 */
#endif
    discrim_chars = Sarray_discriminating_chars(&nbytes,sarrayfile,genomecomp,lcp_bytes,lcp_guide,
						lcp_exceptions,/*guide_interval*/1024,/*n*/genomelength,
						CT_CHARTABLE);
#ifdef MEMUSAGE
    printf("After computing discrim_chars: %lu\n",Mem_usage_report_std_heap()); /* 11,006,768,928 */
#endif
    FREE(sarrayfile);

#ifdef MEMUSAGE
    printf("Before computing child array: %lu\n",Mem_usage_report_std_heap()); /* 11,006,768,852 */
#endif
    fprintf(stderr,"Building child array\n");
    /* Compute child array (relative values) */
    child_bytes = Sarray_compute_child_bytes(&child_exceptions,&n_child_exceptions,lcp_bytes,lcp_guide,lcp_exceptions,/*n*/genomelength);
#ifdef MEMUSAGE
    printf("After computing child array: %lu\n",Mem_usage_report_std_heap()); /* 23,389,611,072 *** */
#endif
    FREE(lcp_exceptions);
    FREE(lcp_guide);

    childexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.sachildexc")+1,sizeof(char));
    sprintf(childexcfile,"%s/%s.metct.sachildexc",destdir,fileroot);
    childguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.sachildguide1024")+1,sizeof(char));
    sprintf(childguidefile,"%s/%s.metct.sachildguide1024",destdir,fileroot);
    Bytecoding_write_bytes(/*bytesfile*/NULL,childexcfile,childguidefile,child_bytes,child_exceptions,n_child_exceptions,
			   genomelength,/*guide_interval*/1024);
    FREE(childguidefile);
    FREE(childexcfile);
    FREE(child_exceptions);

    /* Write combined lcpchilddc file */
    lcpchilddcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metct.salcpchilddc")+1,sizeof(char));
    sprintf(lcpchilddcfile,"%s/%s.metct.salcpchilddc",destdir,fileroot);
    Bytecoding_interleave_lcpchilddc(lcpchilddcfile,child_bytes,discrim_chars,lcp_bytes,genomelength);
    FREE(lcpchilddcfile);
    
    FREE(child_bytes);
    FREE(discrim_chars);
    FREE(lcp_bytes);
#ifdef MEMUSAGE
    printf("At end of suffix array: %lu\n",Mem_usage_report_std_heap()); /* 4,294,968,003 */
#endif
  }

  /* Read offsets */
  /* ref_offsets = Indexdb_offsets_from_bitpack(filenames->pointers_filename,filenames->offsets_filename,index1part); */
#ifdef HAVE_MMAP
  ref_offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,filenames->pointers_filename,/*randomp*/false);
  ref_offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,filenames->offsets_filename,/*randomp*/false);
#else
  ref_offsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_access,&offsetsmeta_len,&seconds,filenames->pointers_filename,sizeof(UINT4));
  ref_offsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_access,&offsetsstrm_len,&seconds,filenames->offsets_filename,sizeof(UINT4));
#endif

#ifdef MEMUSAGE
  printf("Before compute_ct: %lu\n",Mem_usage_report_std_heap()); /* 4,294,968,003 */
#endif
  compute_ct(new_pointers_filename,new_offsets_filename,
	     positions_high_fp,positions_low_fp,
	     ref_offsetsmeta,ref_offsetsstrm,ref_positions8,ref_positions4,
	     oligospace,mask,coord_values_8p);
#ifdef MEMUSAGE
  printf("After compute_ct: %lu\n",Mem_usage_report_std_heap()); /* 4,294,968,003 */
#endif
  if (coord_values_8p == true) {
    fclose(positions_high_fp);
  }
  fclose(positions_low_fp);
  FREE(new_offsets_filename);
  FREE(new_pointers_filename);


  /* Open GA output files */
  new_pointers_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					  strlen(".")+strlen("metga")+strlen(filenames->pointers_index1info_ptr)+1,sizeof(char));
  sprintf(new_pointers_filename,"%s/%s.%s%s",destdir,fileroot,"metga",filenames->pointers_index1info_ptr);

  new_offsets_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("metga")+strlen(filenames->offsets_index1info_ptr)+1,sizeof(char));
  sprintf(new_offsets_filename,"%s/%s.%s%s",destdir,fileroot,"metga",filenames->offsets_index1info_ptr);


  if (coord_values_8p == true) {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metga")+strlen(filenames->positions_high_index1info_ptr)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s",destdir,fileroot,"metga",filenames->positions_high_index1info_ptr);
    
    if ((positions_high_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
      fprintf(stderr,"Can't open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);
  }

  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("metga")+strlen(filenames->positions_low_index1info_ptr)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%s",destdir,fileroot,"metga",filenames->positions_low_index1info_ptr);

  if ((positions_low_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  /* Compute and write GA files */
  if (build_suffix_array_p == true) {
    fprintf(stderr,"Building suffix array for GA\n");
    sarrayfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.sarray")+1,sizeof(char));
    sprintf(sarrayfile,"%s/%s.metga.sarray",destdir,fileroot);
    /* Already have genomecomp open */
    gbuffer = (unsigned char *) CALLOC(genomelength+1,sizeof(unsigned char));
    Genome_fill_buffer_int_string(genomecomp,/*left*/0,/*length*/genomelength,gbuffer,ga_conversion);
    gbuffer[genomelength] = 0;	/* Tried N/X, but SACA_K fails */
    Sarray_write_array_from_genome(sarrayfile,gbuffer,genomelength);

    indexijptrsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.saindex64meta")+1,sizeof(char));
    sprintf(indexijptrsfile,"%s/%s.metga.saindex64meta",destdir,fileroot);
    indexijcompfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.saindex64strm")+1,sizeof(char));
    sprintf(indexijcompfile,"%s/%s.metga.saindex64strm",destdir,fileroot);
    Sarray_write_index_interleaved(indexijptrsfile,indexijcompfile,
				   sarrayfile,genomecomp,genomelength,/*compressp*/true,GA_CHARTABLE);
    FREE(indexijcompfile);
    FREE(indexijptrsfile);

    SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,/*randomp*/false);

#if 0
    /* Not needed if we already have gbuffer */
    /* Required for computing LCP, but uses non-SIMD instructions */
    genomebits = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_BITS,
			    /*uncompressedp*/false,/*access*/USE_MMAP_ONLY);
    Genome_hr_setup(Genome_blocks(genomebits),/*snp_blocks*/NULL,
		    /*query_unk_mismatch_p*/false,/*genome_unk_mismatch_p*/false,
		    /*mode*/CMET_STRANDED);
#endif

    lcp_bytes = Sarray_compute_lcp_bytes_from_genome(&lcp_exceptions,&n_lcp_exceptions,SA,gbuffer,/*n*/genomelength);
    munmap((void *) SA,sa_len);
    close(sa_fd);
    FREE(gbuffer);
#if 0
    Genome_free(&genomebits);
#endif

    /* Write lcp exceptions/guide */
    lcpexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.salcpexc")+1,sizeof(char));
    sprintf(lcpexcfile,"%s/%s.metga.salcpexc",destdir,fileroot);
    lcpguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.salcpguide1024")+1,sizeof(char));
    sprintf(lcpguidefile,"%s/%s.metga.salcpguide1024",destdir,fileroot);
    Bytecoding_write_bytes(/*bytesfile*/NULL,lcpexcfile,lcpguidefile,lcp_bytes,lcp_exceptions,n_lcp_exceptions,
			   genomelength,/*guide_interval*/1024);

    FREE(lcpguidefile);
    FREE(lcpexcfile);


    /* DC array */
    /* Assume we have lcp_bytes and lcp_exceptions already in memory.  Don't need to use guide for speed. */
    lcpguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.salcpguide1024")+1,sizeof(char));
    sprintf(lcpguidefile,"%s/%s.metga.salcpguide1024",destdir,fileroot);
    lcp_guide = (UINT4 *) Access_allocate_private(&lcpguide_access,&lcpguide_len,&seconds,lcpguidefile,sizeof(UINT4));
    FREE(lcpguidefile);

    /* Compute discriminating chars (DC) array */
    discrim_chars = Sarray_discriminating_chars(&nbytes,sarrayfile,genomecomp,lcp_bytes,lcp_guide,
						lcp_exceptions,/*guide_interval*/1024,/*n*/genomelength,
						GA_CHARTABLE);
    FREE(sarrayfile);

    fprintf(stderr,"Building child array\n");
    /* Compute child array (relative values) */
    child_bytes = Sarray_compute_child_bytes(&child_exceptions,&n_child_exceptions,lcp_bytes,lcp_guide,lcp_exceptions,/*n*/genomelength);
    FREE(lcp_exceptions);
    FREE(lcp_guide);

    childexcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.sachildexc")+1,sizeof(char));
    sprintf(childexcfile,"%s/%s.metga.sachildexc",destdir,fileroot);
    childguidefile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.sachildguide1024")+1,sizeof(char));
    sprintf(childguidefile,"%s/%s.metga.sachildguide1024",destdir,fileroot);
    Bytecoding_write_bytes(/*bytesfile*/NULL,childexcfile,childguidefile,child_bytes,child_exceptions,n_child_exceptions,
			   genomelength,/*guide_interval*/1024);
    FREE(childguidefile);
    FREE(childexcfile);
    FREE(child_exceptions);

    /* Write combined lcpchilddc file */
    lcpchilddcfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+strlen(".metga.salcpchilddc")+1,sizeof(char));
    sprintf(lcpchilddcfile,"%s/%s.metga.salcpchilddc",destdir,fileroot);
    Bytecoding_interleave_lcpchilddc(lcpchilddcfile,child_bytes,discrim_chars,lcp_bytes,genomelength);
    FREE(lcpchilddcfile);
    
    FREE(child_bytes);
    FREE(discrim_chars);
    FREE(lcp_bytes);

    Genome_free(&genomecomp);
  }

  compute_ga(new_pointers_filename,new_offsets_filename,
	     positions_high_fp,positions_low_fp,
	     ref_offsetsmeta,ref_offsetsstrm,ref_positions8,ref_positions4,
	     oligospace,mask,coord_values_8p);
  if (coord_values_8p == true) {
    fclose(positions_high_fp);
  }
  fclose(positions_low_fp);
  FREE(new_offsets_filename);
  FREE(new_pointers_filename);



  /* Clean up */
#ifdef HAVE_MMAP
  munmap((void *) ref_offsetsmeta,offsetsmeta_len);
  munmap((void *) ref_offsetsstrm,offsetsstrm_len);
  close(offsetsmeta_fd);
  close(offsetsstrm_fd);
#else
  FREE(ref_offsetsmeta);
  FREE(ref_offsetsstrm);
#endif

  if (coord_values_8p == true) {
    FREE(ref_positions8);
  } else {
#ifdef HAVE_MMAP
    munmap((void *) ref_positions4,ref_positions_low_len);
    close(ref_positions_low_fd);
#else
    FREE(ref_positions4);
#endif
  }

  Filenames_free(&filenames);

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

