static char rcsid[] = "$Id: snpindex.c 184165 2016-02-12 19:14:32Z twu $";
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
#include "types.h"

#include "compress-write.h"
#include "compress.h"
#include "interval.h"
#include "complement.h"

#include "bool.h"
#include "genomicpos.h"
#include "iitdef.h"
#include "uintlist.h"
#include "chrnum.h"
#include "genome.h"
#include "datadir.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "indexdb.h"
#include "indexdb-write.h"
#include "bitpack64-write.h"
#include "bitpack64-read.h"
#include "bitpack64-readtwo.h"
#include "bitpack64-access.h"
#include "bitpack64-incr.h"
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

/* SNP blocks and writing of positions */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

#ifdef CHECK
#define check(x) x
#else
#define check(x)
#endif



static char *user_sourcedir = NULL;
static char *user_destdir = NULL;
static char *dbroot = NULL;
static char *fileroot = NULL;
static int index1part = 15;
static int required_index1part = 0;
static int index1interval = 3;
static int required_interval = 0;

static char *dbversion = NULL;
static int circular_typeint = -1;

static char *snps_root = NULL;
static bool show_warnings_p = true;
static int max_warnings = -1;


static struct option long_options[] = {
  /* Input options */
  {"sourcedir", required_argument, 0, 'D'},	/* user_sourcedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"kmer", required_argument, 0, 'k'},	   /* required_index1part */
  {"sampling", required_argument, 0, 'q'},  /* required_interval */
  {"destdir", required_argument, 0, 'V'},	/* user_destdir */
  {"snpsdb", required_argument, 0, 'v'}, /* snps_root */

  /* Output options */
  {"max-warnings", required_argument, 0, 'w'}, /* max_warnings */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"SNPINDEX: Builds GMAP index files for known SNPs\n");
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

static char
check_acgt (char nt, IIT_T snps_iit, char *divstring,
	    int origindex, Interval_T interval) {
  char *label;
  bool allocp;

  if (nt == 'N') {
    return nt;
  } else if (nt != 'A' && nt != 'C' && nt != 'G' && nt != 'T') {
    label = IIT_label(snps_iit,origindex,&allocp);
    fprintf(stderr,"\nFor %s at %s:%u, alternate allele %c is not ACGT, so using 'N' as alternate allele",
	    label,divstring,Interval_low(interval),nt);
    if (allocp) {
      FREE(label);
    }
    return 'N';
  } else {
    return nt;
  }
}


static Oligospace_T
revise_oligomer (Oligospace_T oligo, char nt) {
  switch (nt) {
  case 'A': return (oligo << 2);
  case 'C': return (oligo << 2) | 1U;
  case 'G': return (oligo << 2) | 2U;
  case 'T': return (oligo << 2) | 3U;
  }
  fprintf(stderr,"altstring char is %c\n",nt);
  exit(9);
  return 0;
}


#if 0
static bool
samep (char *string1, char *string2, int length) {
  int i;

  for (i = 0; i < length; i++) {
    if (string1[i] != string2[i]) {
      return false;
    }
  }
  return true;
}
#endif


typedef struct Labeled_interval_T *Labeled_interval_T;
struct Labeled_interval_T {
  int origindex;
  Interval_T interval;
};

static Labeled_interval_T
Labeled_interval_new (int origindex, Interval_T interval) {
  Labeled_interval_T new = (Labeled_interval_T) MALLOC(sizeof(*new));

  new->origindex = origindex;
  new->interval = interval;
  return new;
}

static void
Labeled_interval_free (Labeled_interval_T *old, bool free_interval_p) {
  if (free_interval_p == true) {
    Interval_free(&((*old)->interval));
  }
  FREE(*old);
  return;
}

static int
Labeled_interval_cmp (const void *a, const void *b) {
  Labeled_interval_T x = * (Labeled_interval_T *) a;
  Labeled_interval_T y = * (Labeled_interval_T *) b;

  return Interval_cmp((const void *) &(x->interval),(const void *) &(y->interval));
}



/* Same procedure used for writing offsets and writing positions, depending on whether positions4/positions8 is NULL */
static int
process_snp_block (int *nwarnings, char *packsizes, UINT4 **bitpacks,
		   UINT4 *snponly_offsetsmeta, UINT4 *snponly_offsetsstrm, UINT4 *countermeta, UINT4 *counterstrm,
		   UINT4 *positions4, UINT8 *positions8, Labeled_interval_T *intervals, int nintervals,
		   Univcoord_T chroffset, Genome_T genome,
		   Genomecomp_T *snp_blocks, char *divstring,
		   IIT_T snps_iit, int index1part, bool coord_values_8p) {
  int nerrors = 0;
  bool *snpp;
  char *refstring;
  char *altstring;
#ifdef DEBUG1
  char *nt;
#endif

  int length;
  int nsnps, stringi, starti, shift, i, k;
  char *snptype, *label, refnt, altnt;
  Univcoord_T ptr;
  Univcoord_T snpposition, startposition, endposition, first_snppos, last_snppos, position;
  Chrpos_T chrpos;
  Interval_T interval;
  bool badcharp, allocp;
  
  Uintlist_T oligomers, newoligomers, p;
  Oligospace_T oligo, bmer;
  UINT4 offset;
#ifdef WORDS_BIGENDIAN
  Genomecomp_T high, low, flags;
#endif


  /* Subtract 1 because snps_iit is 1-based */
  first_snppos = Interval_low(intervals[0]->interval) - 1U;
  last_snppos = Interval_low(intervals[nintervals-1]->interval) - 1U;
  debug1(
	 if (nintervals == 1) {
	   printf("Processing snp at chrpos %s:%u\n",divstring,first_snppos+1U);
	 } else {
	   printf("Processing block of %d snps from chrpos %s:%u to %u\n",
		  nintervals,divstring,first_snppos+1U,last_snppos+1U);
	 }
	 );

  if (first_snppos < (index1part - 1)) {
    startposition = chroffset;
  } else {
    startposition = chroffset + first_snppos - (index1part - 1);
  }
  endposition = chroffset + last_snppos + (index1part - 1);
  length = endposition - startposition + 1;

  snpp = (bool *) CALLOC(length,sizeof(bool));
  refstring = (char *) CALLOC(length + 1,sizeof(char));
  altstring = (char *) CALLOC(length + 1,sizeof(char));

#if 0
  /* Doesn't work with circular chromosomes */
  Genome_fill_buffer(&chrnum,&nunknowns,genome,/*left*/startposition,length,
		     refstring,chromosome_iit);
#else
  Genome_fill_buffer_simple(genome,/*left*/startposition,length,refstring);
#endif

  for (i = 0; i < nintervals; i++) {
    interval = intervals[i]->interval;
    snpposition = chroffset + Interval_low(interval) - 1U; /* Subtract 1 because snps_iit is 1-based */

    debug1(label = IIT_label(snps_iit,intervals[i]->origindex,&allocp));
    debug1(printf("  Neighbor %s at %s:%u\n",label,divstring,Interval_low(interval)));

    stringi = snpposition - startposition;

    snptype = IIT_typestring(snps_iit,Interval_type(interval));
    if (strlen(snptype) != 2) {
      fprintf(stderr,"Unrecognized snptype %s\n",snptype);
      abort();
    } else {
      refnt = refstring[stringi];
      if (refnt == snptype[0]) {
	if (altstring[stringi] != '\0' && altstring[stringi] != snptype[1]) {
	  nerrors++;
	  if (show_warnings_p == true) {
	    label = IIT_label(snps_iit,intervals[i]->origindex,&allocp);
	    fprintf(stderr,"\nFor %s at %s:%u, saw two different alternate alleles %c and %c, so using N as alternate allele.",
		    label,divstring,Interval_low(interval),altstring[stringi],snptype[1]);
	    if (allocp == true) {
	      FREE(label);
	    }
	    if (++(*nwarnings) == max_warnings) {
	      fprintf(stderr,"\nMaximum of %d warnings reached.  No more warnings will be shown\n",max_warnings);
	      show_warnings_p = false;
	    }
	  }
	  altnt = 'N';
	} else {
	  altnt = check_acgt(snptype[1],snps_iit,divstring,intervals[i]->origindex,intervals[i]->interval);
	}
      } else if (refnt == snptype[1]) {
	if (altstring[stringi] != '\0' && altstring[stringi] != snptype[0]) {
	  nerrors++;
	  if (show_warnings_p == true) {
	    label = IIT_label(snps_iit,intervals[i]->origindex,&allocp);
	    fprintf(stderr,"\nFor %s at %s:%u, saw two different alternate alleles %c and %c, so using N as alternate allele.",
		    label,divstring,Interval_low(interval),altstring[stringi],snptype[0]);
	    if (allocp == true) {
	      FREE(label);
	    }
	    if (++(*nwarnings) == max_warnings) {
	      fprintf(stderr,"\nMaximum of %d warnings reached.  No more warnings will be shown\n",max_warnings);
	      show_warnings_p = false;
	    }
	  }
	  altnt = 'N';
	} else {
	  altnt = check_acgt(snptype[0],snps_iit,divstring,intervals[i]->origindex,intervals[i]->interval);
	}
      } else {
	nerrors++;
	if (show_warnings_p == true) {
	  label = IIT_label(snps_iit,intervals[i]->origindex,&allocp);
	  fprintf(stderr,"\nFor %s at %s:%u, snptype %s not consistent with reference allele %c, so ignoring.",
		  label,divstring,Interval_low(interval),snptype,refstring[stringi]);
	  if (allocp == true) {
	    FREE(label);
	  }
	  if (++(*nwarnings) == max_warnings) {
	    fprintf(stderr,"\nMaximum of %d warnings reached.  No more warnings will be shown\n",max_warnings);
	    show_warnings_p = false;
	  }
	}
	altnt = '\0';		/* Ignoring */
      }

      if (altnt == '\0') {
	/* Skip */
      } else if (altnt == refnt) {
	fprintf(stderr,"\nAt %s:%u, alternate allele %c is same as reference allele\n",
		divstring,Interval_low(interval),altnt);
      } else {
	altstring[stringi] = altnt;
	snpp[stringi] = true;
	if (snp_blocks != NULL) {
	  /* revising genome */
	  ptr = snpposition/32U*3;
	  shift = snpposition % 32U;

#ifdef WORDS_BIGENDIAN
	  flags = Bigendian_convert_uint(snp_blocks[ptr+2]);
	  flags |= (1 << shift);
	  snp_blocks[ptr+2] = Bigendian_convert_uint(flags);
#else
	  snp_blocks[ptr+2] |= (1 << shift);	/* Flags.  Change even for 'N'. */
#endif
	  if (altnt == 'N') {
	    /* refnt + flag indicates 'N' */
	    altnt = refnt;	/* Change back to refnt, if necessary. */
	  }

	  if (shift >= 16) {
	    /* high */
	    shift -= 16;
#ifdef WORDS_BIGENDIAN
	    high = Bigendian_convert_uint(snp_blocks[ptr]);
	    high &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': high |= (0x1 << (2*shift)); break;
	    case 'G': high |= (0x2 << (2*shift)); break;
	    case 'T': high |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
	    snp_blocks[ptr] = Bigendian_convert_uint(high);
#else
	    snp_blocks[ptr] &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': snp_blocks[ptr] |= (0x1 << (2*shift)); break;
	    case 'G': snp_blocks[ptr] |= (0x2 << (2*shift)); break;
	    case 'T': snp_blocks[ptr] |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
#endif
	  } else {
	    /* low */
#ifdef WORDS_BIGENDIAN
	    low = Bigendian_convert_uint(snp_blocks[ptr+1]);
	    low &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': low |= (0x1 << (2*shift)); break;
	    case 'G': low |= (0x2 << (2*shift)); break;
	    case 'T': low |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
	    snp_blocks[ptr+1] = Bigendian_convert_uint(low);
#else
	    snp_blocks[ptr+1] &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': snp_blocks[ptr+1] |= (0x1 << (2*shift)); break;
	    case 'G': snp_blocks[ptr+1] |= (0x2 << (2*shift)); break;
	    case 'T': snp_blocks[ptr+1] |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
#endif
	  }
	}
      }
    }
  }
  
  for (starti = 0, position = startposition, chrpos = startposition - chroffset; 
       position <= endposition-index1part+1U; starti++, position++, chrpos++) {
    if (chrpos % index1interval == 0) {
      /* chrpos % 3 == 0 is same as the condition in indexdb.c for storing a position */
      nsnps = 0;
      badcharp = false;
      for (k = starti; k < starti + index1part; k++) {
	if (snpp[k] == true) {
	  nsnps++;
	}
	if (refstring[k] != 'A' && refstring[k] != 'C' && refstring[k] != 'G' && refstring[k] != 'T') {
	  badcharp = true;
	}
      }
      if (nsnps == 0) {
	/* no snps */
	/* fprintf(stderr,"\nNo snps at position %llu, %s:%u",(unsigned long long) position,divstring,Interval_low(interval)); */
      } else if (nsnps > 4) {
	/* too many snps */
      } else if (badcharp == true) {
	/* bad reference char */
      } else {
	oligomers = Uintlist_push(NULL,0U);
	for (k = starti, i = 0; k < starti + index1part; k++, i++) {
	  newoligomers = NULL;
	  if (snpp[k] == false) {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	    }
	    
	  } else if (altstring[k] != 'N') {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),altstring[k]));
	    }

	  } else {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	      if (refstring[k] != 'A') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'A'));
	      }
	      if (refstring[k] != 'C') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'C'));
	      }
	      if (refstring[k] != 'G') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'G'));
	      }
	      if (refstring[k] != 'T') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'T'));
	      }
	    }
	  }

	  Uintlist_free(&oligomers);
	  oligomers = Uintlist_reverse(newoligomers);
	}

#ifdef CHECK
	for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	  oligo = Uintlist_head(p);
	  nt = shortoligo_nt(oligo,index1part);
	  if (samep(nt,&(refstring[starti]),index1part) == true) {
	    fprintf(stderr,"Storing oligomer %s that is the same as the reference at %llu (%s:%u)\n",
		    nt,(unsigned long long) position,divstring,chrpos+1U);
	    abort();
	  }
	  FREE(nt);
	}
#endif	

	/* Ignore the first element in oligomers, which is all reference */
	if (positions4 == NULL && positions8 == NULL) {
	  /* writing offsets */
	  for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	    oligo = Uintlist_head(p);
	    /* offsets[oligo + 1U] += 1; */
	    bmer = oligo/BLOCKSIZE;
	    if (Bitpack64_access_filledp(oligo,(int) packsizes[bmer],bitpacks[bmer]) == true) {
	      bitpacks[bmer] = Bitpack64_realloc_one((int) packsizes[bmer],bitpacks[bmer]);
	      packsizes[bmer] += 2;
	    }
	    Bitpack64_incr_bitpack(oligo,(int) packsizes[bmer],bitpacks[bmer]);
	    debug1(nt = shortoligo_nt(oligo,index1part);
		   printf("Storing %s at %llu (%s:%u)\n",nt,(unsigned long long) position,divstring,chrpos+1U);
		   FREE(nt));
	  }

	} else {
	  /* writing positions */
	  if (coord_values_8p == true) {
	    for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	      oligo = Uintlist_head(p);
	      offset = Bitpack64_read_one(oligo,snponly_offsetsmeta,snponly_offsetsstrm) + Bitpack64_incr(oligo,countermeta,counterstrm);
	      positions8[offset] = position;
	      debug1(nt = shortoligo_nt(oligo,index1part);
		     printf("Storing %s at %llu (%s:%u)\n",nt,(unsigned long long) position,divstring,chrpos+1U);
		     FREE(nt));
	    }
	    
	  } else {
	    for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	      oligo = Uintlist_head(p);
	      offset = Bitpack64_read_one(oligo,snponly_offsetsmeta,snponly_offsetsstrm) + Bitpack64_incr(oligo,countermeta,counterstrm);
	      positions4[offset] = (UINT4) position;
	      debug1(nt = shortoligo_nt(oligo,index1part);
		     printf("Storing %s at %llu (%s:%u)\n",nt,(unsigned long long) position,divstring,chrpos+1U);
		     FREE(nt));
	    }
	  }
	}
	Uintlist_free(&oligomers);
      }
    }
  }

  FREE(altstring);
  FREE(refstring);
  FREE(snpp);

  return nerrors;
}


/* Copied from bitpack64-write.c */
static UINT4
compute_ascending (UINT4 *ascending, UINT4 *counts) {
  int i;

  ascending[0] = 0;
  for (i = 1; i <= 64; i++) {
    ascending[i] = ascending[i-1] + counts[i-1];
  }

  return ascending[64];
}


static void
write_offsets (char *snponly_offsetsmetafile, char *snponly_offsetsstrmfile,
	       char *snp_offsetsmetafile, char *snp_offsetsstrmfile, IIT_T snps_iit,
	       UINT4 *ref_offsetsmeta, UINT4 *ref_offsetsstrm,
	       Univ_IIT_T chromosome_iit, Genome_T genome,
	       Genomecomp_T *snp_blocks, Oligospace_T oligospace, int index1part) {
  char *packsizes;
  UINT4 **bitpacks;

  Labeled_interval_T *intervals;
  int origindex;
  Interval_T interval, copy;
  int nintervals, nintervals_alias, nerrors, divno, i, j;
  Oligospace_T bmerspace, bmer, oligoi;
  char *divstring;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Chrpos_T chrlength;
  int nwarnings = 0;

  FILE *ptrs_fp, *comp_fp;
  UINT4 *ptrs, nregisters;
  UINT4 totalcount;
  UINT4 ptri;
  UINT4 ptr, end;

  /* Buffer is used to avoid frequent writes to the file */
  UINT4 *buffer;
  int buffer_size = BUFFER_SIZE;
  int buffer_i;

  UINT4 diffs[BLOCKSIZE], ascending[BLOCKSIZE+1], counts[BLOCKSIZE], offsets[BLOCKSIZE+1];
  int packsize;


  bmerspace = oligospace/BLOCKSIZE;
  packsizes = (char *) CALLOC(bmerspace,sizeof(char));
  bitpacks = (UINT4 **) CALLOC(bmerspace,sizeof(UINT4 *));

  for (divno = 1; divno < snps_iit->ndivs; divno++) {
    divstring = IIT_divstring(snps_iit,divno);
    fprintf(stderr,"Processing snp offsets for chromosome %s...",divstring);
    if ((chrnum = Univ_IIT_find_one(chromosome_iit,divstring)) <= 0) {
      fprintf(stderr,"not found in chromosome iit\n");
    } else {
      nerrors = 0;
      fprintf(stderr,"has %d snps...",IIT_nintervals(snps_iit,divno));
      chroffset = Univ_IIT_interval_low(chromosome_iit,chrnum);
      chrlength = Univ_IIT_interval_length(chromosome_iit,chrnum);
      nintervals = IIT_nintervals(snps_iit,divno);

      if (Univ_IIT_interval_type(chromosome_iit,chrnum) == circular_typeint) {
	fprintf(stderr,"and is circular...");
	nintervals_alias = 2*nintervals;
	intervals = (Labeled_interval_T *) CALLOC(nintervals_alias,sizeof(Labeled_interval_T));
	for (i = 0; i < nintervals; i++) {
	  origindex = IIT_index(snps_iit,divno,i);
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/&(snps_iit->intervals[divno][i]));
	}
	for (j = 0; i < nintervals_alias; i++, j++) {
	  origindex = IIT_index(snps_iit,divno,j);
	  interval = intervals[j]->interval;
	  copy = Interval_new(/*low*/Interval_low(interval) + chrlength,
			      /*high*/Interval_high(interval) + chrlength,
			      Interval_type(interval));
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/copy);
	}

      } else {
	nintervals_alias = nintervals;
	intervals = (Labeled_interval_T *) CALLOC(nintervals,sizeof(Labeled_interval_T));
	for (i = 0; i < nintervals; i++) {
	  origindex = IIT_index(snps_iit,divno,i);
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/&(snps_iit->intervals[divno][i]));
	}
      }

      qsort(intervals,nintervals,sizeof(Labeled_interval_T),Labeled_interval_cmp);

      i = 0;
      while (i < nintervals_alias) {
	j = i + 1;
	while (j < nintervals_alias && Interval_low(intervals[j]->interval) < Interval_low(intervals[j-1]->interval) + index1part) {
	  j++;
	}
	nerrors += process_snp_block(&nwarnings,packsizes,bitpacks,
				     /*snponly_offsetsmeta*/NULL,/*snponly_offsetsstrm*/NULL,/*countermeta*/NULL,/*counterstrm*/NULL,
				     /*positions4*/NULL,/*positions8*/NULL,&(intervals[i]),/*nintervals*/j-i,
				     chroffset,genome,snp_blocks,
				     divstring,snps_iit,index1part,
				     /*coord_values_8p (irrelevant)*/false);
	i = j;
      }

      for (i = nintervals; i < nintervals_alias; i++) {
	Labeled_interval_free(&(intervals[i]),/*free_interval_p*/true);
      }
      for (i = 0; i < nintervals; i++) {
	Labeled_interval_free(&(intervals[i]),/*free_interval_p*/false);
      }

      FREE(intervals);
      fprintf(stderr,"done (%d snps inconsistent with reference genome)\n",nerrors);
    }
  }

  Bitpack64_write_differential_bitpacks(snponly_offsetsmetafile,snponly_offsetsstrmfile,packsizes,bitpacks,oligospace);


  /* Code derived from Bitpack64_write_differential_bitpacks */
  fprintf(stderr,"Computing combined ref+snp offsets");
  if ((comp_fp = FOPEN_WRITE_BINARY(snp_offsetsstrmfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",snp_offsetsstrmfile);
    exit(9);
  }

  ptrs = (UINT4 *) CALLOC(((oligospace + BLOCKSIZE)/BLOCKSIZE + 1) * DIFFERENTIAL_METAINFO_SIZE,sizeof(UINT4));
  ptri = 0;
  buffer = (UINT4 *) CALLOC(buffer_size,sizeof(UINT4));
  buffer_i = 0;
  nregisters = 0;

  totalcount = 0;
  for (oligoi = 0, bmer = 0; oligoi + BLOCKSIZE <= oligospace; oligoi += BLOCKSIZE, bmer++) {
    ptrs[ptri++] = nregisters;	/* In 128-bit registers */
    ptrs[ptri++] = totalcount;
    Bitpack64_extract_bitpack(counts,packsizes[bmer],bitpacks[bmer]);
    Bitpack64_block_offsets(offsets,oligoi,ref_offsetsmeta,ref_offsetsstrm);
    for (i = 0; i < 64; i++) {
      if ((oligoi + i) % MONITOR_INTERVAL == 0) {
	fprintf(stderr,".");
      }
      counts[i] += (offsets[i+1] - offsets[i]);
    }
    totalcount += compute_ascending(ascending,counts);
    packsize = Bitpack64_compute_q4_diffs_bidir(diffs,ascending); /* Note: This packsize may differ from packsizes[bmer] */
    buffer_i = Bitpack64_write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

    nregisters += packsize / 2;
  }
  fprintf(stderr,"done\n");

  /* For nucleotides, expect a single final block where oligoi == oligospace */
  if (oligoi <= oligospace) {
    ptrs[ptri++] = nregisters;	/* In 128-bit registers */
    ptrs[ptri++] = totalcount;

    if (oligoi == oligospace) {
      /* Don't have a bitpack at [bmerspace] */
      Bitpack64_extract_bitpack(counts,/*packsize*/0,/*bitpack*/NULL);
    } else {
      Bitpack64_extract_bitpack(counts,packsizes[bmer],bitpacks[bmer]);
    }

    /* For differential, want <=.  For direct, want < */
    for (i = 0; i <= (int) (oligospace - oligoi); i++) {
      ptr = Bitpack64_read_two(&end,oligoi,ref_offsetsmeta,ref_offsetsstrm);
      counts[i] += (end - ptr);
    }
    for ( ; i < BLOCKSIZE; i++) {
      /* Copy last value for rest of block */
      counts[i] = 0;
    }

    /* Pack block of < 64 diffs */
    totalcount += compute_ascending(ascending,counts);
    packsize = Bitpack64_compute_q4_diffs_bidir(diffs,ascending);
    buffer_i = Bitpack64_write_columnar(comp_fp,buffer,buffer_size,buffer_i,diffs,packsize);

    nregisters += packsize / 2;
  }

  /* Write the final pointers, which will point after the end of the file */
  ptrs[ptri++] = nregisters;	/* In 128-bit registers */
  ptrs[ptri++] = totalcount;

  if ((ptrs_fp = FOPEN_WRITE_BINARY(snp_offsetsmetafile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",snp_offsetsmetafile);
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


  /* Clean up */
  FREE(ptrs);
  for (bmer = 0; bmer < bmerspace; bmer++) {
    FREE(bitpacks[bmer]);
  }
  FREE(bitpacks);
  FREE(packsizes);

  return;
}


static void *
compute_positions (UINT4 *snponly_offsetsmeta, UINT4 *snponly_offsetsstrm, IIT_T snps_iit, Univ_IIT_T chromosome_iit,
		   Genome_T genome, Oligospace_T oligospace, bool coord_values_8p) {
  UINT4 *positions4 = NULL;
  UINT8 *positions8 = NULL;
  UINT4 offsets[BLOCKSIZE+1];

  Labeled_interval_T *intervals;
  int origindex;
  Interval_T interval, copy;
  int nintervals, nintervals_alias, divno, i, j;
  Oligospace_T oligoi;
  char *divstring;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Chrpos_T chrlength;
  Positionsptr_T totalcounts, npositions;
  UINT4 *countermeta, *counterstrm;
  int nwarnings = 0;

  totalcounts = Bitpack64_read_one(oligospace,snponly_offsetsmeta,snponly_offsetsstrm);
  if (totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the SNP-only offsets.  Total counts is zero.\n");
    fprintf(stderr,"Do the chromosomes in the IIT file match those in the genome?\n");
    fprintf(stderr,"Here are known chromosomes in the genome: ");
    Univ_IIT_dump_labels(stderr,chromosome_iit);
    fprintf(stderr,"Here are chromosomes in the SNPs IIT file: ");
    IIT_dump_divstrings(stderr,snps_iit);
    exit(9);
  } else if (coord_values_8p == true) {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(UINT8));
    positions4 = (UINT4 *) NULL;
    positions8 = (UINT8 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT8));
    if (positions8 == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(UINT4));
    positions8 = (UINT8 *) NULL;
    positions4 = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions4 == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  }

  /* Prepare counter */
  countermeta = Indexdb_bitpack_counter(&counterstrm,snponly_offsetsmeta,snponly_offsetsstrm,index1part);

  for (divno = 1; divno < snps_iit->ndivs; divno++) {
    divstring = IIT_divstring(snps_iit,divno);
    fprintf(stderr,"Processing positions for chromosome %s...",divstring);
    if ((chrnum = Univ_IIT_find_one(chromosome_iit,divstring)) <= 0) {
      fprintf(stderr,"not found in chromosome iit\n");
    } else {
      fprintf(stderr,"has %d snps...",IIT_nintervals(snps_iit,divno));
      chroffset = Univ_IIT_interval_low(chromosome_iit,chrnum);
      chrlength = Univ_IIT_interval_length(chromosome_iit,chrnum);
      nintervals = IIT_nintervals(snps_iit,divno);

      if (Univ_IIT_interval_type(chromosome_iit,chrnum) == circular_typeint) {
	fprintf(stderr,"and is circular...");
	nintervals_alias = 2*nintervals;
	intervals = (Labeled_interval_T *) CALLOC(nintervals_alias,sizeof(Labeled_interval_T));
	for (i = 0; i < nintervals; i++) {
	  origindex = IIT_index(snps_iit,divno,i);
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/&(snps_iit->intervals[divno][i]));
	}
	for (j = 0; i < nintervals_alias; i++, j++) {
	  origindex = IIT_index(snps_iit,divno,j);
	  interval = intervals[j]->interval;
	  copy = Interval_new(/*low*/Interval_low(interval) + chrlength,
			      /*high*/Interval_high(interval) + chrlength,
			      Interval_type(interval));
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/copy);
	}

      } else {
	nintervals_alias = nintervals;
	intervals = (Labeled_interval_T *) CALLOC(nintervals,sizeof(Labeled_interval_T));
	for (i = 0; i < nintervals; i++) {
	  origindex = IIT_index(snps_iit,divno,i);
	  intervals[i] = Labeled_interval_new(origindex,/*interval*/&(snps_iit->intervals[divno][i]));
	}
      }

      qsort(intervals,nintervals,sizeof(Labeled_interval_T),Labeled_interval_cmp);

      i = 0;
      while (i < nintervals_alias) {
	j = i + 1;
	while (j < nintervals_alias && Interval_low(intervals[j]->interval) < Interval_low(intervals[j-1]->interval) + index1part) {
	  j++;
	}
	process_snp_block(&nwarnings,/*packsizes*/NULL,/*bitpacks*/NULL,
			  snponly_offsetsmeta,snponly_offsetsstrm,countermeta,counterstrm,
			  positions4,positions8,&(intervals[i]),/*nintervals*/j-i,
			  chroffset,genome,/*snp_blocks*/NULL,
			  divstring,snps_iit,index1part,coord_values_8p);
	i = j;
      }

      for (i = nintervals; i < nintervals_alias; i++) {
	Labeled_interval_free(&(intervals[i]),/*free_interval_p*/true);
      }
      for (i = 0; i < nintervals; i++) {
	Labeled_interval_free(&(intervals[i]),/*free_interval_p*/false);
      }

      FREE(intervals);
      fprintf(stderr,"done\n");
    }
  }

  FREE(counterstrm);
  FREE(countermeta);

  /* Sort positions in each block */
  fprintf(stderr,"Sorting snp positions");
  if (coord_values_8p == true) {
    for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
      Bitpack64_block_offsets(offsets,oligoi,snponly_offsetsmeta,snponly_offsetsstrm);
      for (i = 0; i < 64; i++) {
	if ((oligoi + i) % MONITOR_INTERVAL == 0) {
	  fprintf(stderr,".");
	}
	if ((npositions = offsets[i+1] - offsets[i]) > 1) {
	  qsort(&(positions8[offsets[i]]),npositions,sizeof(UINT8),UINT8_compare);
	}
      }
    }
    fprintf(stderr,"done\n");
    return positions8;

  } else {
    for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
      Bitpack64_block_offsets(offsets,oligoi,snponly_offsetsmeta,snponly_offsetsstrm);
      for (i = 0; i < 64; i++) {
	if ((oligoi + i) % MONITOR_INTERVAL == 0) {
	  fprintf(stderr,".");
	}
	if ((npositions = offsets[i+1] - offsets[i]) > 1) {
	  qsort(&(positions4[offsets[i]]),npositions,sizeof(UINT4),UINT4_compare);
	}
      }
    }
    fprintf(stderr,"done\n");
    return positions4;
  }
}


static void
merge_positions8 (FILE *positions_high_fp, FILE *positions_low_fp,
		  UINT8 *start1, UINT8 *end1, UINT8 *start2, UINT8 *end2, 
		  Oligospace_T oligo, int index1part) {
  unsigned char position8_high;
  UINT4 position8_low;
  UINT8 *ptr1 = start1, *ptr2 = start2;
  char *nt;
#ifdef WORDS_BIGENDIAN
  UINT8 position1, position2;
#endif

  while (ptr1 < end1 && ptr2 < end2) {
#ifdef WORDS_BIGENDIAN
    position1 = Bigendian_convert_uint8(*ptr1);
    position2 = Bigendian_convert_uint8(*ptr2);
    if (position1 < position2) {
      position8_high = position1 >> POSITIONS8_HIGH_SHIFT;
      position8_low = position1 & POSITIONS8_LOW_MASK;
      FWRITE_CHAR(position8_high,positions_high_fp);
      FWRITE_UINT(position8_low,positions_low_fp);
      ptr1++;
    } else if (position2 < position1) {
      position8_high = position2 >> POSITIONS8_HIGH_SHIFT;
      position8_low = position2 & POSITIONS8_LOW_MASK;
      FWRITE_CHAR(position8_high,positions_high_fp);
      FWRITE_UINT(position8_low,positions_low_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,index1part);
      fprintf(stderr,"Problem: saw duplicate positions %u in oligo %s\n",position1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT8(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }

#else
    if (*ptr1 < *ptr2) {
      position8_high = *ptr1 >> POSITIONS8_HIGH_SHIFT;
      position8_low = *ptr1 & POSITIONS8_LOW_MASK;
      FWRITE_CHAR(position8_high,positions_high_fp);
      FWRITE_UINT(position8_low,positions_low_fp);
      ptr1++;
    } else if (*ptr2 < *ptr1) {
      position8_high = *ptr2 >> POSITIONS8_HIGH_SHIFT;
      position8_low = *ptr2 & POSITIONS8_LOW_MASK;
      FWRITE_CHAR(position8_high,positions_high_fp);
      FWRITE_UINT(position8_low,positions_low_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,index1part);
      fprintf(stderr,"Problem: saw duplicate positions %llu in oligo %s\n",(unsigned long long) *ptr1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT8(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }
#endif
  }

#ifdef WORDS_BIGENDIAN
  while (ptr1 < end1) {
    position1 = Bigendian_convert_uint8(*ptr1);
    position8_high = position1 >> POSITIONS8_HIGH_SHIFT;
    position8_low = position1 & POSITIONS8_LOW_MASK;

    FWRITE_CHAR(position8_high,positions_high_fp);
    FWRITE_UINT(position8_low,positions_low_fp);
    ptr1++;
  }

#else
  while (ptr1 < end1) {
    position8_high = *ptr1 >> POSITIONS8_HIGH_SHIFT;
    position8_low = *ptr1 & POSITIONS8_LOW_MASK;
    FWRITE_CHAR(position8_high,positions_high_fp);
    FWRITE_UINT(position8_low,positions_low_fp);
    ptr1++;
  }
#endif

#ifdef WORDS_BIGENDIAN
  while (ptr2 < end2) {
    position2 = Bigendian_convert_uint8(*ptr2);
    position8_high = position2 >> POSITIONS8_HIGH_SHIFT;
    position8_low = position2 & POSITIONS8_LOW_MASK;

    FWRITE_CHAR(position8_high,positions_high_fp);
    FWRITE_UINT(position8_low,positions_low_fp);
    ptr2++;
  }
#else
  while (ptr2 < end2) {
    position8_high = *ptr2 >> POSITIONS8_HIGH_SHIFT;
    position8_low = *ptr2 & POSITIONS8_LOW_MASK;

    FWRITE_CHAR(position8_high,positions_high_fp);
    FWRITE_UINT(position8_low,positions_low_fp);
    ptr2++;
  }
#endif

  return;
}


static void
merge_positions4 (FILE *positions_fp, UINT4 *start1, UINT4 *end1,
		  UINT4 *start2, UINT4 *end2, Oligospace_T oligo, int index1part) {
  UINT4 *ptr1 = start1, *ptr2 = start2;
  char *nt;
#ifdef WORDS_BIGENDIAN
  UINT4 position1, position2;
#endif

  while (ptr1 < end1 && ptr2 < end2) {
#ifdef WORDS_BIGENDIAN
    position1 = Bigendian_convert_uint(*ptr1);
    position2 = Bigendian_convert_uint(*ptr2);
    if (position1 < position2) {
      FWRITE_UINT(position1,positions_fp);
      ptr1++;
    } else if (position2 < position1) {
      FWRITE_UINT(position2,positions_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,index1part);
      fprintf(stderr,"Problem: saw duplicate positions %u in oligo %s\n",*ptr1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }

#else

    if (*ptr1 < *ptr2) {
      FWRITE_UINT(*ptr1,positions_fp);
      ptr1++;
    } else if (*ptr2 < *ptr1) {
      FWRITE_UINT(*ptr2,positions_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,index1part);
      fprintf(stderr,"Problem: saw duplicate positions %u in oligo %s\n",*ptr1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }
#endif
  }

#ifdef WORDS_BIGENDIAN
  while (ptr1 < end1) {
    position1 = Bigendian_convert_uint(*ptr1);
    FWRITE_UINT(position1,positions_fp);
    ptr1++;
  }
#else
  if (ptr1 < end1) {
    FWRITE_UINTS(ptr1,end1-ptr1,positions_fp);
  }
#endif

#ifdef WORDS_BIGENDIAN
  while (ptr2 < end2) {
    position2 = Bigendian_convert_uint(*ptr2);
    FWRITE_UINT(position2,positions_fp);
    ptr2++;
  }
#else
  if (ptr2 < end2) {
    FWRITE_UINTS(ptr2,end2-ptr2,positions_fp);
  }
#endif

  return;
}


/* Usage: snpindex -d <genome> -V <destdir> -v <snps_root> */


/* Note: Depends on having gmapindex sampled on mod 3 bounds */
int
main (int argc, char *argv[]) {
  char *sourcedir = NULL, *destdir = NULL, *mapdir = NULL;
  int compression_type;
  Univ_IIT_T chromosome_iit;
  IIT_T snps_iit;
  Genome_T genome;
  UINT4 *ref_offsetsmeta, *ref_offsetsstrm, *snponly_offsetsmeta, *snponly_offsetsstrm;

  size_t totalcounts, i;
#ifdef EXTRA_ALLOCATION
  Positionsptr_T npositions;
#endif

  unsigned char *ref_positions8_high;
  UINT4 *ref_positions8_low;
  UINT8 *snp_positions8, *ref_positions8;
  UINT4 *snp_positions4, *ref_positions4;
  Univcoord_T nblocks;
  Genomecomp_T *snp_blocks;
  Oligospace_T oligospace, oligoi;
  int ii;
#ifndef HAVE_MMAP
  double seconds;
#endif

  bool coord_values_8p;
  Filenames_T filenames;
  char *filename, *snp_offsetsmetafile, *snp_offsetsstrmfile, *snponly_offsetsmetafile, *snponly_offsetsstrmfile;
  FILE *genome_fp, *genomebits_fp, *positions_high_fp, *positions_low_fp;
  int ref_offsetsmeta_fd, ref_offsetsstrm_fd, ref_positions_high_fd, ref_positions_low_fd,
    snponly_offsetsmeta_fd, snponly_offsetsstrm_fd;
  size_t ref_offsetsmeta_len, ref_offsetsstrm_len, ref_positions_high_len, ref_positions_low_len,
    snponly_offsetsmeta_len, snponly_offsetsstrm_len;
  UINT4 ref_offsets[BLOCKSIZE+1], snponly_offsets[BLOCKSIZE+1];

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:k:q:V:v:w:",
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
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'snpindex --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_sourcedir = optarg; break;
    case 'd': dbroot = optarg; break;
    case 'k': required_index1part = atoi(optarg); break;
    case 'q': required_interval = atoi(optarg); break;
    case 'V': user_destdir = optarg; break;
    case 'v': snps_root = optarg; break;
    case 'w': max_warnings = atoi(optarg); break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (dbroot == NULL) {
    fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
    print_program_usage();
    exit(9);
  } else if (snps_root == NULL) {
    fprintf(stderr,"Missing name of SNP database.  Must specify with -v flag.\n");
    print_program_usage();
    exit(9);
  } else {
    sourcedir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_sourcedir,dbroot);
    fprintf(stderr,"Reading source files from %s\n",sourcedir);
  }

  if (argc > 1) {
    /* IIT file provided as an argument */
    if (Access_file_exists_p(argv[1]) == false) {
      fprintf(stderr,"SNP IIT file %s not found\n",argv[1]);
      exit(9);
    } else {
      fprintf(stderr,"Reading SNPs IIT file %s...",argv[1]);
      if ((snps_iit = IIT_read(argv[1],/*name*/NULL,/*readonlyp*/true,
			       /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false)) == NULL) {
	fprintf(stderr,"SNP IIT file %s is not valid\n",argv[1]);
	exit(9);
      }
      fprintf(stderr,"done\n");
    } 

  } else {
    mapdir = Datadir_find_mapdir(user_sourcedir,sourcedir,fileroot);
    filename = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(snps_root)+strlen(".iit")+1,sizeof(char));
    sprintf(filename,"%s/%s.iit",mapdir,snps_root);
    if (Access_file_exists_p(filename) == false) {
      fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",snps_root,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s.iit or specify IIT file as an argument\n",snps_root);
      exit(9);
    } else {
      if ((snps_iit = IIT_read(filename,/*name*/NULL,/*readonlyp*/true,
			       /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false)) == NULL) {
	fprintf(stderr,"SNP IIT file %s is not valid\n",filename);
	exit(9);
      }
      fprintf(stderr,"done\n");
    }
    FREE(filename);
    FREE(mapdir);
  }


  /* Chromosome IIT file */
  filename = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(filename,"%s/%s.chromosome.iit",sourcedir,fileroot);
  if ((chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"IIT file %s is not valid\n",filename);
    exit(9);
  } else {
    circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  }
  FREE(filename);
  coord_values_8p = Univ_IIT_coord_values_8p(chromosome_iit);

  fprintf(stderr,"Chromosomes in the genome: ");
  Univ_IIT_dump_labels(stderr,chromosome_iit);
  fprintf(stderr,"Chromosomes in the SNPs IIT file: ");
  IIT_dump_divstrings(stderr,snps_iit);

  genome = Genome_new(sourcedir,fileroot,/*snps_root*/NULL,/*genometype*/GENOME_OLIGOS,
		      /*uncompressedp*/false,/*access*/USE_MMAP_ONLY,/*sharedp*/false);

  /* Copy reference genome as basis for snp genome */
  nblocks = Genome_totallength(genome)/32U;
  snp_blocks = (Genomecomp_T *) CALLOC(nblocks*3,sizeof(Genomecomp_T));
  fprintf(stderr,"Allocating %llu*3*%d bytes for compressed genome\n",
	  (unsigned long long) nblocks,(int) sizeof(Genomecomp_T));
  memcpy(snp_blocks,Genome_blocks(genome),nblocks*3*sizeof(Genomecomp_T));

  /* Prepare for write */
  if (user_destdir == NULL) {
    destdir = sourcedir;
  } else {
    destdir = user_destdir;
  }
  fprintf(stderr,"Writing snpindex files to %s\n",destdir);

  if (max_warnings == 0) {
    show_warnings_p = false;
  }

  oligospace = power(4,index1part);
  filenames = Indexdb_get_filenames(&compression_type,&index1part,&index1interval,
				    sourcedir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
				    required_index1part,required_interval,
				    /*offsets_only_p*/false);

  /* Write offsets and compute snp_blocks */
#ifdef HAVE_MMAP
  ref_offsetsmeta = (UINT4 *) Access_mmap(&ref_offsetsmeta_fd,&ref_offsetsmeta_len,filenames->pointers_filename,/*randomp*/false);
  ref_offsetsstrm = (UINT4 *) Access_mmap(&ref_offsetsstrm_fd,&ref_offsetsstrm_len,filenames->offsets_filename,/*randomp*/false);
#else
  ref_offsetsmeta = (UINT4 *) Access_allocate_private(&ref_offsetsmeta_len,&seconds,filenames->pointers_filename,sizeof(UINT4));
  ref_offsetsstrm = (UINT4 *) Access_allocate_private(&ref_offsetsstrm_len,&seconds,filenames->offsets_filename,sizeof(UINT4));
#endif

  /* Write snponly and SNP offsets */
  snp_offsetsmetafile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(filenames->pointers_basename_ptr)+
					strlen(".")+strlen(snps_root)+1,sizeof(char));
  sprintf(snp_offsetsmetafile,"%s/%s.%s",destdir,filenames->pointers_basename_ptr,snps_root);
  snp_offsetsstrmfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(filenames->offsets_basename_ptr)+
					strlen(".")+strlen(snps_root)+1,sizeof(char));
  sprintf(snp_offsetsstrmfile,"%s/%s.%s",destdir,filenames->offsets_basename_ptr,snps_root);

  snponly_offsetsmetafile = (char *) CALLOC(strlen(snp_offsetsmetafile) + strlen(".only") + 1,sizeof(char));
  sprintf(snponly_offsetsmetafile,"%s.only",snp_offsetsmetafile);
  snponly_offsetsstrmfile = (char *) CALLOC(strlen(snp_offsetsstrmfile) + strlen(".only") + 1,sizeof(char));
  sprintf(snponly_offsetsstrmfile,"%s.only",snp_offsetsstrmfile);

  write_offsets(snponly_offsetsmetafile,snponly_offsetsstrmfile,snp_offsetsmetafile,snp_offsetsstrmfile,
		snps_iit,ref_offsetsmeta,ref_offsetsstrm,chromosome_iit,genome,snp_blocks,oligospace,index1part);

  FREE(snp_offsetsmetafile);
  FREE(snp_offsetsstrmfile);


  /* Write genomecomp */
  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".genomecomp.")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.genomecomp.%s",destdir,fileroot,snps_root);
  if ((genome_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing genome\n",filename);
    exit(9);
  }
  fprintf(stderr,"Writing genome file %s...",filename);
  FWRITE_UINTS(snp_blocks,nblocks*3,genome_fp);
  fclose(genome_fp);
  FREE(snp_blocks);
  fprintf(stderr,"done\n");


  /* Write genomebits */
  if ((genome_fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for genome\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".genomebits128.")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.genomebits128.%s",destdir,fileroot,snps_root);
  if ((genomebits_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing genome\n",filename);
    exit(9);
  }
  fprintf(stderr,"Writing genome file %s...",filename);
  Compress_unshuffle_bits128(/*out*/genomebits_fp,/*in*/genome_fp);
  fclose(genomebits_fp);
  fclose(genome_fp);
  FREE(filename);
  fprintf(stderr,"done\n");


  /* Compute snponly positions and keep in memory */
  show_warnings_p = false;	/* Already shown in write_offsets */

#ifdef HAVE_MMAP
  snponly_offsetsmeta = (UINT4 *) Access_mmap(&snponly_offsetsmeta_fd,&snponly_offsetsmeta_len,snponly_offsetsmetafile,/*randomp*/false);
  snponly_offsetsstrm = (UINT4 *) Access_mmap(&snponly_offsetsstrm_fd,&snponly_offsetsstrm_len,snponly_offsetsstrmfile,/*randomp*/false);
#else
  snponly_offsetsmeta = (UINT4 *) Access_allocate_private(&snponly_offsetsmeta_len,&seconds,snponly_offsetsmetafile,sizeof(UINT4));
  snponly_offsetsstrm = (UINT4 *) Access_allocate_private(&snponly_offsetsstrm_len,&seconds,snponly_offsetsstrmfile,sizeof(UINT4));
#endif
  if (coord_values_8p == true) {
    snp_positions8 = compute_positions(snponly_offsetsmeta,snponly_offsetsstrm,snps_iit,chromosome_iit,genome,oligospace,/*coord_values_8p*/true);
  } else {
    snp_positions4 = compute_positions(snponly_offsetsmeta,snponly_offsetsstrm,snps_iit,chromosome_iit,genome,oligospace,/*coord_values_8p*/false);
  }


  /* Read reference positions and merge */
  if (coord_values_8p == true) {
#ifdef HAVE_MMAP
    ref_positions8_high = (unsigned char *) Access_mmap(&ref_positions_high_fd,&ref_positions_high_len,
							filenames->positions_high_filename,/*randomp*/false);
    ref_positions8_low = (UINT4 *) Access_mmap(&ref_positions_low_fd,&ref_positions_low_len,
					       filenames->positions_low_filename,/*randomp*/false);
#else
    ref_positions8_high = (unsigned char *) Access_allocate_private(&ref_positions_high_len,&seconds,
								    filenames->positions_high_filename,sizeof(unsigned char));
    ref_positions8_low = (UINT4 *) Access_allocate_private(&ref_positions_low_len,&seconds,
							   filenames->positions_low_filename,sizeof(UINT4));
#endif
    /* Unpack 8-byte positions into memory */
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
    FREE(ref_positions8_high);	/* Not using shared */
    FREE(ref_positions8_low);	/* Not using shared */
    /* Access_deallocate(ref_positions8_high,ref_positions_high_shmid,ref_positions_high_key); */
    /* Access_deallocate(ref_positions8_low,ref_positions_low_shmid,ref_positions_low_key); */
#endif

  } else {
#ifdef HAVE_MMAP
    ref_positions4 = (UINT4 *) Access_mmap(&ref_positions_low_fd,&ref_positions_low_len,
					   filenames->positions_low_filename,/*randomp*/false);
#else
    ref_positions4 = (UINT4 *) Access_allocate_private(&ref_positions_low_len,&seconds,
						       filenames->positions_low_filename,sizeof(UINT4));
#endif
  }


  if (coord_values_8p == true) {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(filenames->positions_high_basename_ptr)+
			       strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s",destdir,filenames->positions_high_basename_ptr,snps_root);
    if ((positions_high_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
      fprintf(stderr,"Can't open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);
  }

  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(filenames->positions_low_basename_ptr)+
			     strlen(".")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s",destdir,filenames->positions_low_basename_ptr,snps_root);
  if ((positions_low_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);


  fprintf(stderr,"Merging snp positions with reference positions");
  if (coord_values_8p == true) {
    for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
      Bitpack64_block_offsets(snponly_offsets,oligoi,snponly_offsetsmeta,snponly_offsetsstrm);
      Bitpack64_block_offsets(ref_offsets,oligoi,ref_offsetsmeta,ref_offsetsstrm);
      for (ii = 0; ii < BLOCKSIZE; ii++) {
	if ((oligoi + ii) % MONITOR_INTERVAL == 0) {
	  fprintf(stderr,".");
	}
	merge_positions8(positions_high_fp,positions_low_fp,
			 &(snp_positions8[snponly_offsets[ii]]),&(snp_positions8[snponly_offsets[ii+1]]),
			 &(ref_positions8[ref_offsets[ii]]),&(ref_positions8[ref_offsets[ii+1]]),oligoi + ii,index1part);
      }
    }
  } else {
    for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
      Bitpack64_block_offsets(snponly_offsets,oligoi,snponly_offsetsmeta,snponly_offsetsstrm);
      Bitpack64_block_offsets(ref_offsets,oligoi,ref_offsetsmeta,ref_offsetsstrm);
      for (ii = 0; ii < BLOCKSIZE; ii++) {
	if ((oligoi + ii) % MONITOR_INTERVAL == 0) {
	  fprintf(stderr,".");
	}
	merge_positions4(positions_low_fp,&(snp_positions4[snponly_offsets[ii]]),&(snp_positions4[snponly_offsets[ii+1]]),
			 &(ref_positions4[ref_offsets[ii]]),&(ref_positions4[ref_offsets[ii+1]]),oligoi + ii,index1part);
      }
    }
  }
  fprintf(stderr,"done\n");

  fclose(positions_low_fp);
  if (coord_values_8p == true) {
    fclose(positions_high_fp);
  }


  /* Clean up */
  if (coord_values_8p == true) {
    /* For both mmap and allocated, since we have already combined positions_high and positions_low */
    FREE(ref_positions8);
  } else {
#ifdef HAVE_MMAP
    munmap((void *) ref_positions4,ref_positions_low_len);
    close(ref_positions_low_fd);
#else
    FREE(ref_positions4);	/* Not using shared */
    /* Access_deallocate(ref_positions4,ref_positions_low_shmid,ref_positions_low_key); */
#endif
  }

#ifdef HAVE_MMAP
  munmap((void *) ref_offsetsmeta,ref_offsetsmeta_len);
  munmap((void *) ref_offsetsstrm,ref_offsetsstrm_len);
  close(ref_offsetsmeta_fd);
  close(ref_offsetsstrm_fd);
#else
  FREE(ref_offsetsmeta);
  FREE(ref_offsetsstrm);
#endif

  if (coord_values_8p == true) {
    FREE(snp_positions8);
  } else {
    FREE(snp_positions4);
  }

#ifdef HAVE_MMAP
  munmap((void *) snponly_offsetsmeta,snponly_offsetsmeta_len);
  munmap((void *) snponly_offsetsstrm,snponly_offsetsstrm_len);
  close(snponly_offsetsmeta_fd);
  close(snponly_offsetsstrm_fd);
#else
  FREE(snponly_offsetsmeta);
  FREE(snponly_offsetsstrm);
#endif
  remove(snponly_offsetsmetafile);
  remove(snponly_offsetsstrmfile);
  FREE(snponly_offsetsmetafile);
  FREE(snponly_offsetsstrmfile);


  Genome_free(&genome);
  Univ_IIT_free(&chromosome_iit);
  IIT_free(&snps_iit);

  Filenames_free(&filenames);


  if (argc <= 1) {
    /* Program called using -v flag only.  No need to install. */
    /* fprintf(stderr,"IIT file already present in .maps directory\n"); */
#if 0
    /* Old code used for copying IIT file to .maps directory */
    /* To use this code, cannot free mapdir earlier */
    fprintf(stderr,"Now copying IIT file from %s/%s.iit to %s...",mapdir,snps_root,filename);
    filename1 = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(snps_root)+strlen(".iit")+1,sizeof(char));
    sprintf(filename1,"%s/%s.iit",mapdir,snps_root);
    Access_file_copy(/*dest*/filename,/*source*/filename1);
    fprintf(stderr,"done\n");
    FREE(filename1);
    FREE(mapdir);
#endif

  } else {
    /* Install IIT file */
    if (!strcmp(destdir,sourcedir)) {
      filename = (char *) CALLOC(strlen(destdir)+strlen("/")+ strlen(fileroot) + strlen(".maps/") +
				 strlen(snps_root)+strlen(".iit")+1,sizeof(char));
      sprintf(filename,"%s/%s.maps/%s.iit",destdir,fileroot,snps_root);
    } else {
      filename = (char *) CALLOC(strlen(destdir)+strlen("/")+
				 strlen(snps_root)+strlen(".iit")+1,sizeof(char));
      sprintf(filename,"%s/%s.iit",destdir,snps_root);
    }

    fprintf(stderr,"SNP genome indices created.\n");
    if (Access_file_exists_p(filename) == true) {
      if (Access_file_equal(filename,argv[1]) == false) {
	fprintf(stderr,"IIT file already present as %s, but it is different from the given file %s\n",
		filename,argv[1]);
	fprintf(stderr,"Please copy file %s as %s\n",argv[1],filename);
      } else {
	fprintf(stderr,"IIT file already present as %s, and it is the same as given file %s\n",
		filename,argv[1]);
      }

    } else {
      fprintf(stderr,"Now installing IIT file %s as %s...",
	      argv[1],filename);
      Access_file_copy(/*dest*/filename,/*source*/argv[1]);
      fprintf(stderr,"done\n");
    }

    FREE(filename);
  }

  FREE(dbversion);
  FREE(fileroot);
  FREE(sourcedir);

  return 0;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: snpindex [OPTIONS...] -d <genome> -v <snpsdb> [<iitfile>]\n\
\n\
If iitfile is provided as a non-flag argument, then use that iitfile and create SNP database\n\
as named by -v flag.  Otherwise, try to find iit file named <snpsdb>.iit in GMAP index files\n\
for <genome>.\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Options (must include -d)\n");
  fprintf(stdout,"\
  -D, --sourcedir=directory      Directory where to read genome index files (default is\n\
                                   GMAP genome directory specified at compile time)\n\
  -d, --db=STRING                Genome database\n\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less).\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  -q, --sampling=INT             Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected k-mer size\n\
  -V, --destdir=directory        Directory where to write SNP index files (default is\n\
                                   GMAP genome directory specified at compile time)\n\
  -v, --snpsdb=STRING            Name of SNP database\n\
  -w, --max-warnings=INT         Maximum number of warnings to print to stderr about\n\
                                   inconsistencies relative to the reference genome.\n\
                                   A value of 0 turns off all warnings.  A negative value\n\
                                   prints all warnings.  (default -1, meaning no limit)\n\
\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

