static char rcsid[] = "$Id: sarray-read.c 207324 2017-06-14 19:41:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "sarray-read.h"

#ifdef WORDS_BIGENDIAN
#define CONVERT(x) Bigendian_convert_uint(x)
#include "bigendian.h"
#else
#define CONVERT(x) (x)
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */

#include "mem.h"
#include "bool.h"
#include "assert.h"
#include "access.h"
#include "types.h"
#include "genomicpos.h"
#include "genome128_hr.h"
#include "bytecoding.h"
#include "bitpack64-read.h"
#include "bitpack64-readtwo.h"
#include "bitpack64-access.h"

#include "sarray-read.h"


#ifdef USE_CSA

#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
#include <emmintrin.h>
#endif

#endif


#define MAX_DEBUG1_HITS 100

/* Details of suffix array search */
#ifdef DEBUG1
#include "genome.h"
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Search through saindex */
#ifdef DEBUG1A
#define debug1a(x) x
#else
#define debug1a(x)
#endif

/* get_child */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Compressed suffix array */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Compressed suffix array: comparison with sarray */
#ifdef DEBUG3A
#define debug3a(x) x
#else
#define debug3a(x)
#endif

/* Compressed suffix array: comparison with csa phi */
#ifdef DEBUG3B
#define debug3b(x) x
#else
#define debug3b(x)
#endif

/* Compare separate buckets with a single one */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif



#define T Sarray_T
struct T {
  Univcoord_T n;
  Univcoord_T n_plus_one;

  /* Old format */
  int array_shmid;
  key_t array_key;
  Univcoord_T *array;

  int lcpchilddc_shmid;
  key_t lcpchilddc_key;
  unsigned char *lcpchilddc;

  int lcp_guide_shmid;
  key_t lcp_guide_key;
  int lcp_exceptions_shmid;
  key_t lcp_exceptions_key;
  UINT4 *lcp_guide;
  UINT4 *lcp_exceptions;
  int n_lcp_exceptions;		/* Won't be necessary if we change lcpchilddc to use guide array */
  /* int lcp_guide_interval; -- Always use 1024 */
  
  int child_guide_shmid;
  key_t child_guide_key;
  int child_exceptions_shmid;
  key_t child_exceptions_key;
  UINT4 *child_guide;
  UINT4 *child_exceptions;
  /* int n_child_exceptions; */
  int child_guide_interval; /* Always use 1024 */

#if 0
  Sarrayptr_T initindexi[4];	/* For A, C, G, T */
  Sarrayptr_T initindexj[4];	/* For A, C, G, T */
#endif

  int indexsize;
  UINT4 indexspace;		/* 4^indexsize.  Used by sarray_read to detect when we have a poly-T oligo shorter than indexsize */
#ifdef DEBUG15
  UINT4 *indexi_ptrs, *indexi_comp, *indexj_ptrs, *indexj_comp; /* bucket array: oligomer lookup into suffix array */
  UINT4 *indexij_ptrs, *indexij_comp;
#elif defined(USE_SEPARATE_BUCKETS)
  UINT4 *indexi_ptrs, *indexi_comp, *indexj_ptrs, *indexj_comp; /* bucket array: oligomer lookup into suffix array */
#else
  int indexij_ptrs_shmid;
  key_t indexij_ptrs_key;
  int indexij_comp_shmid;
  key_t indexij_comp_key;
  UINT4 *indexij_ptrs, *indexij_comp;
#endif

  Access_T array_access; int array_fd; size_t array_len;

#ifdef DEBUG15
  int indexi_ptrs_fd; size_t indexi_ptrs_len; int indexi_comp_fd; size_t indexi_comp_len;
  int indexj_ptrs_fd; size_t indexj_ptrs_len; int indexj_comp_fd; size_t indexj_comp_len;
  int indexij_ptrs_fd; size_t indexij_ptrs_len; int indexij_comp_fd; size_t indexij_comp_len;
#elif defined(USE_SEPARATE_BUCKETS)
  int indexi_ptrs_fd; size_t indexi_ptrs_len; int indexi_comp_fd; size_t indexi_comp_len;
  int indexj_ptrs_fd; size_t indexj_ptrs_len; int indexj_comp_fd; size_t indexj_comp_len;
#else
  Access_T indexij_ptrs_access; int indexij_ptrs_fd; size_t indexij_ptrs_len;
  Access_T indexij_comp_access; int indexij_comp_fd; size_t indexij_comp_len;
#endif

  Access_T lcpchilddc_access; int lcpchilddc_fd; size_t lcpchilddc_len;

  Access_T lcp_guide_access; int lcp_guide_fd; size_t lcp_guide_len;
  Access_T lcp_exceptions_access; int lcp_exceptions_fd; size_t lcp_exceptions_len;

  Access_T child_guide_access; int child_guide_fd; size_t child_guide_len;
  Access_T child_exceptions_access; int child_exceptions_fd; size_t child_exceptions_len;

};


/* For benchmarking */
Univcoord_T
Sarray_size (Sarray_T this) {
  return this->n_plus_one;
}




#if 0
/* Simplified from sarray_search_simple in sarray-write.c */
static void
sarray_search_char (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, char desired_char,
		    UINT4 *SA, UINT4 n, char *chartable) {
  Sarrayptr_T low, high, mid;
  Univcoord_T pos;
  char c;

  low = 1;
  high = n + 1;

  while (low < high) {
#if 0
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif
#ifdef WORDS_BIGENDIAN
    pos = Bigendian_convert_uint(SA[mid]);
#else
    pos = SA[mid];
#endif
    c = Genome_get_char_lex(genome,pos,n,chartable);
    if (desired_char > c) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  *initptr = low;

  low--;
  high = n;
  while (low < high) {
#if 1
    /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 || high % 2 == 1) {
      mid += 1;
    }
#else
    /* This does not work for ceiling */
    mid = low + ((high - low) / 2);
#endif
#ifdef WORDS_BIGENDIAN
    pos = Bigendian_convert_uint(SA[mid]);
#else
    pos = SA[mid];
#endif
    c = Genome_get_char_lex(genome,pos,n,chartable);
    if (desired_char >= c) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  *finalptr = high;
  return;
}
#endif



static int
log4 (int result) {
  int exponent = 0;

  while (result > 1) {
    result /= 4;
    exponent++;
  }

  return exponent;
}

static UINT4
power (int base, int exponent) {
  UINT4 result = 1;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }

  return result;
}


#if 0
void
Sarray_shmem_remove (char *dir, char *fileroot, char *snps_root, Mode_T mode, bool fwdp) {
  char *mode_prefix;
  char *sarrayfile;
  char *lcpchilddcfile;
  char *lcp_guidefile, *lcp_exceptionsfile;
  char *child_guidefile, *child_exceptionsfile;
  char *indexij_ptrsfile, *indexij_compfile;

  if (mode == STANDARD) {
    mode_prefix = ".";
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".metct.";
    } else {
      mode_prefix = ".metga.";
    }
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2iag.";
    } else {
      mode_prefix = ".a2itc.";
    }
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2itc.";
    } else {
      mode_prefix = ".a2iag.";
    }
  }

  sarrayfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sarray")+1,sizeof(char));
  sprintf(sarrayfile,"%s/%s%ssarray",dir,fileroot,mode_prefix);

  lcpchilddcfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpchilddc")+1,sizeof(char));
  sprintf(lcpchilddcfile,"%s/%s%ssalcpchilddc",dir,fileroot,mode_prefix);

  lcp_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpguide1024")+1,sizeof(char));
  sprintf(lcp_guidefile,"%s/%s%ssalcpguide1024",dir,fileroot,mode_prefix);
  lcp_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpexc")+1,sizeof(char));
  sprintf(lcp_exceptionsfile,"%s/%s%ssalcpexc",dir,fileroot,mode_prefix);

  child_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildguide1024")+1,sizeof(char));
  sprintf(child_guidefile,"%s/%s%ssachildguide1024",dir,fileroot,mode_prefix);
  child_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildexc")+1,sizeof(char));
  sprintf(child_exceptionsfile,"%s/%s%ssachildexc",dir,fileroot,mode_prefix);

  indexij_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64meta")+1,sizeof(char));
  sprintf(indexij_ptrsfile,"%s/%s%ssaindex64meta",dir,fileroot,mode_prefix);
  indexij_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64strm")+1,sizeof(char));
  sprintf(indexij_compfile,"%s/%s%ssaindex64strm",dir,fileroot,mode_prefix);

  Access_shmem_remove(indexij_ptrsfile);
  Access_shmem_remove(indexij_compfile);

  Access_shmem_remove(sarrayfile);
  Access_shmem_remove(lcpchilddcfile);
  Access_shmem_remove(lcp_guidefile);
  Access_shmem_remove(lcp_exceptionsfile);

  Access_shmem_remove(child_guidefile);
  Access_shmem_remove(child_exceptionsfile);

  FREE(child_exceptionsfile);
  FREE(child_guidefile);

  FREE(lcp_exceptionsfile);
  FREE(lcp_guidefile);

  FREE(lcpchilddcfile);

  FREE(sarrayfile);

  return;
}
#endif


Univcoord_T *
Sarray_array (T this) {
  return this->array;
}


#ifdef USE_CSA

Univcoord_T
Sarray_position (T sarray, Sarrayptr_T i) {
  Univcoord_T nhops = 0, expected_sa_i;
  Sarrayptr_T expected_i;
  __m128i converted, cmp;
  int matchbits;

  debug3(printf("Entered Sarray_position for %u:",i));
#ifdef DEBUG3A
  expected_sa_i = sarray->array[i];
#endif

  if (
#ifdef DEBUG3A
      0 && 
#endif
      sarray->array != NULL) {
    debug3(printf("Returning %u\n",sarray->array[i]));
    return sarray->array[i];
  } else {
    while ((i % sarray->sa_sampling) != 0) {
      debug3(printf(",%u",i));
#ifdef DEBUG3B
      expected_i = sarray->csa[i];
#endif

#if defined(HAVE_SSE2) && !defined(WORDS_BIGENDIAN)
      converted = _mm_sub_epi32(_mm_set1_epi32(i),epi32_convert);
      cmp = _mm_cmpgt_epi32(converted,sarray->indices0); /* To use cmpgt, sarray->indices0 is shifted down by 1 */
      matchbits = _mm_movemask_ps(_mm_castsi128_ps(cmp));
      /* assert(matchbits == 0 || matchbits == 1 || matchbits == 3 || matchbits == 7 || matchbits == 15); */
      debug3(printf("(%d)",matchbits));
      i = Bitpack64_read_one(i - sarray->index0[matchbits],sarray->csa0ptrs[matchbits],sarray->csa0comp[matchbits]);
#else
      if (i >= sarray->indexX) {
	assert(matchbits == 15);
	printf("X");
	i = Bitpack64_read_one(i-sarray->indexX,sarray->csaXptrs,sarray->csaXcomp);
      } else if (i >= sarray->indexT) {
	assert(matchbits == 7);
	printf("T");
	i = Bitpack64_read_one(i-sarray->indexT,sarray->csaTptrs,sarray->csaTcomp);
      } else if (i >= sarray->indexG) {
	assert(matchbits == 3);
	printf("G");
	i = Bitpack64_read_one(i-sarray->indexG,sarray->csaGptrs,sarray->csaGcomp);
      } else if (i >= sarray->indexC) {
	assert(matchbits == 1);
	printf("C");
	i = Bitpack64_read_one(i-sarray->indexC,sarray->csaCptrs,sarray->csaCcomp);
      } else {
	assert(matchbits == 0);
	printf("A");
	i = Bitpack64_read_one(i-sarray->indexA,sarray->csaAptrs,sarray->csaAcomp);
      }
#endif

      debug3b(assert(i == expected_i));
      nhops += 1;
    }

    debug3(printf("\n"));
    debug3(printf("Returning %u = %u - nhops %u\n",
		   sarray->array_samples[i/sarray->sa_sampling] - nhops,
		   sarray->array_samples[i/sarray->sa_sampling],nhops));
    
    debug3a(assert(sarray->array_samples[i/sarray->sa_sampling] - nhops == expected_sa_i));

    return sarray->array_samples[i/sarray->sa_sampling] - nhops;
  }
}

#elif defined(WORDS_BIGENDIAN)

Univcoord_T
Sarray_position (T sarray, Sarrayptr_T i) {
  return Bigendian_convert_uint(sarray->array[i]);
}

#else

Univcoord_T
Sarray_position (T sarray, Sarrayptr_T i) {
  return sarray->array[i];
}

#endif


T
Sarray_new (char *dir, char *fileroot, Access_mode_T sarray_access, Access_mode_T lcp_access,
	    Access_mode_T guideexc_access, Access_mode_T indexij_access, bool sharedp, Mode_T mode, bool fwdp) {
  T new;
  char *comma1;
  double seconds;
  int npages;

  bool old_format_p;
  char *sarrayfile;		/* Old format */

  char *lcpchilddcfile;
  char *lcp_guidefile, *lcp_exceptionsfile;
  char *child_guidefile, *child_exceptionsfile;
#ifdef DEBUG15
  char *indexi_ptrsfile, *indexi_compfile;
  char *indexj_ptrsfile, *indexj_compfile;
  char *indexij_ptrsfile, *indexij_compfile;
#elif defined(USE_SEPARATE_BUCKETS)
  char *indexi_ptrsfile, *indexi_compfile;
  char *indexj_ptrsfile, *indexj_compfile;
#else
  char *indexij_ptrsfile, *indexij_compfile;
#endif

  char *mode_prefix;

  if (mode == STANDARD) {
    mode_prefix = ".";
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".metct.";
    } else {
      mode_prefix = ".metga.";
    }
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2iag.";
    } else {
      mode_prefix = ".a2itc.";
    }
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2itc.";
    } else {
      mode_prefix = ".a2iag.";
    }
  }

  /* Old format */
  sarrayfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sarray")+1,sizeof(char));
  sprintf(sarrayfile,"%s/%s%ssarray",dir,fileroot,mode_prefix);

  lcpchilddcfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpchilddc")+1,sizeof(char));
  sprintf(lcpchilddcfile,"%s/%s%ssalcpchilddc",dir,fileroot,mode_prefix);

  lcp_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpguide1024")+1,sizeof(char));
  sprintf(lcp_guidefile,"%s/%s%ssalcpguide1024",dir,fileroot,mode_prefix);
  lcp_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpexc")+1,sizeof(char));
  sprintf(lcp_exceptionsfile,"%s/%s%ssalcpexc",dir,fileroot,mode_prefix);

  child_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildguide1024")+1,sizeof(char));
  sprintf(child_guidefile,"%s/%s%ssachildguide1024",dir,fileroot,mode_prefix);
  child_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildexc")+1,sizeof(char));
  sprintf(child_exceptionsfile,"%s/%s%ssachildexc",dir,fileroot,mode_prefix);

#ifdef DEBUG15
  indexi_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64meta")+1,sizeof(char));
  sprintf(indexi_ptrsfile,"%s/%s%ssaindexi64meta",dir,fileroot,mode_prefix);
  indexi_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64strm")+1,sizeof(char));
  sprintf(indexi_compfile,"%s/%s%ssaindexi64strm",dir,fileroot,mode_prefix);
  indexj_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64meta")+1,sizeof(char));
  sprintf(indexj_ptrsfile,"%s/%s%ssaindexj64meta",dir,fileroot,mode_prefix);
  indexj_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64strm")+1,sizeof(char));
  sprintf(indexj_compfile,"%s/%s%ssaindexj64strm",dir,fileroot,mode_prefix);
  indexij_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindex64meta")+1,sizeof(char));
  sprintf(indexij_ptrsfile,"%s/%s%ssaindex64meta",dir,fileroot,mode_prefix);
  indexij_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindex64strm")+1,sizeof(char));
  sprintf(indexij_compfile,"%s/%s%ssaindex64strm",dir,fileroot,mode_prefix);
#elif defined(USE_SEPARATE_BUCKETS)
  indexi_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64meta")+1,sizeof(char));
  sprintf(indexi_ptrsfile,"%s/%s%ssaindexi64meta",dir,fileroot,mode_prefix);
  indexi_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64strm")+1,sizeof(char));
  sprintf(indexi_compfile,"%s/%s%ssaindexi64strm",dir,fileroot,mode_prefix);
  indexj_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64meta")+1,sizeof(char));
  sprintf(indexj_ptrsfile,"%s/%s%ssaindexj64meta",dir,fileroot,mode_prefix);
  indexj_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64strm")+1,sizeof(char));
  sprintf(indexj_compfile,"%s/%s%ssaindexj64strm",dir,fileroot,mode_prefix);
#else
  indexij_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64meta")+1,sizeof(char));
  sprintf(indexij_ptrsfile,"%s/%s%ssaindex64meta",dir,fileroot,mode_prefix);
  indexij_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64strm")+1,sizeof(char));
  sprintf(indexij_compfile,"%s/%s%ssaindex64strm",dir,fileroot,mode_prefix);
#endif

  if (Access_file_exists_p(sarrayfile) == false) {
    fprintf(stderr,"No suffix array for genome\n");
    new = (T) NULL;

  } else if (Access_file_exists_p(lcpchilddcfile) == false) {
    fprintf(stderr,"Enhanced suffix array file %s does not exist.  The genome was built using an obsolete version\n",
	    lcpchilddcfile);
    new = (T) NULL;
    exit(9);

  } else {
    new = (T) MALLOC_KEEP(sizeof(*new));
    old_format_p = true;

    if (sarray_access == USE_MMAP_PRELOAD) {
      if (old_format_p == true) {
	fprintf(stderr,"Pre-loading suffix array...");
	new->array = (Univcoord_T *) Access_mmap_and_preload(&new->array_fd,&new->array_len,&npages,&seconds,sarrayfile,
						       sizeof(UINT4));
	new->n_plus_one = new->array_len/sizeof(UINT4); /* Should be genomiclength + 1*/
	new->n = new->n_plus_one - 1;

	comma1 = Genomicpos_commafmt(new->array_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      }
      new->array_access = MMAPPED;

    } else if (sarray_access == USE_MMAP_ONLY) {
      if (old_format_p == true) {
	new->array = (Univcoord_T *) Access_mmap(&new->array_fd,&new->array_len,sarrayfile,/*randomp*/true);
	new->n_plus_one = new->array_len/sizeof(UINT4); /* Should be genomiclength + 1*/
	new->n = new->n_plus_one - 1;
      }
      new->array_access = MMAPPED;

    } else if (sarray_access == USE_ALLOCATE) {
      if (old_format_p == true) {
	fprintf(stderr,"Allocating memory for suffix array...");
	if (sharedp == true) {
	  new->array = (Univcoord_T *) Access_allocate_shared(&new->array_access,&new->array_shmid,&new->array_key,
							&new->array_fd,&new->array_len,&seconds,sarrayfile,sizeof(UINT4));
	} else {
	  new->array = (Univcoord_T *) Access_allocate_private(&new->array_access,&new->array_len,&seconds,sarrayfile,sizeof(UINT4));
	}
	new->n_plus_one = new->array_len/sizeof(UINT4); /* Should be genomiclength + 1*/
	new->n = new->n_plus_one - 1;
	comma1 = Genomicpos_commafmt(new->array_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      }
    }

#ifdef DEBUG15
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    new->indexi_ptrs = (UINT4 *) Access_allocate_private(&new->indexi_ptrs_len,&seconds,indexi_ptrsfile,sizeof(UINT4));
    new->indexi_comp = (UINT4 *) Access_allocate_private(&new->indexi_comp_len,&seconds,indexi_compfile,sizeof(UINT4));
    new->indexj_ptrs = (UINT4 *) Access_allocate_private(&new->indexj_ptrs_len,&seconds,indexj_ptrsfile,sizeof(UINT4));
    new->indexj_comp = (UINT4 *) Access_allocate_private(&new->indexj_comp_len,&seconds,indexj_compfile,sizeof(UINT4));
    new->indexij_ptrs = (UINT4 *) Access_allocate_private(&new->indexij_ptrs_len,&seconds,indexij_ptrsfile,sizeof(UINT4));
    new->indexij_comp = (UINT4 *) Access_allocate_private(&new->indexij_comp_len,&seconds,indexij_compfile,sizeof(UINT4));
    new->indexsize = 3 + log4(((new->indexij_ptrs_len - 8)/sizeof(UINT4)/2)/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#elif defined(USE_SEPARATE_BUCKETS)
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    new->indexi_ptrs = (UINT4 *) Access_allocate_private(&new->indexi_ptrs_access,&new->indexi_ptrs_len,&seconds,indexi_ptrsfile,sizeof(UINT4));
    new->indexi_comp = (UINT4 *) Access_allocate_private(&new->indexi_comp_access,&new->indexi_comp_len,&seconds,indexi_compfile,sizeof(UINT4));
    new->indexj_ptrs = (UINT4 *) Access_allocate_private(&new->indexj_ptrs_access,&new->indexj_ptrs_len,&seconds,indexj_ptrsfile,sizeof(UINT4));
    new->indexj_comp = (UINT4 *) Access_allocate_private(&new->indexj_comp_access,&new->indexj_comp_len,&seconds,indexj_compfile,sizeof(UINT4));
    new->indexsize = 3 + log4(((new->indexi_ptrs_len - 8)/sizeof(UINT4))/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#else
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    if (indexij_access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading indexij ptrs...");
      new->indexij_ptrs = (UINT4 *) Access_mmap_and_preload(&new->indexij_ptrs_fd,&new->indexij_ptrs_len,&npages,&seconds,indexij_ptrsfile,
							    sizeof(UINT4));
      comma1 = Genomicpos_commafmt(new->indexij_ptrs_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
      FREE(comma1);

      fprintf(stderr,"Pre-loading indexij comp...");
      new->indexij_comp = (UINT4 *) Access_mmap_and_preload(&new->indexij_comp_fd,&new->indexij_comp_len,&npages,&seconds,indexij_compfile,
							    sizeof(UINT4));
      comma1 = Genomicpos_commafmt(new->indexij_comp_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
      FREE(comma1);

      new->indexij_ptrs_access = MMAPPED;
      new->indexij_comp_access = MMAPPED;

    } else if (indexij_access == USE_MMAP_ONLY) {
      new->indexij_ptrs = (UINT4 *) Access_mmap(&new->indexij_ptrs_fd,&new->indexij_ptrs_len,indexij_ptrsfile,/*randomp*/true);
      new->indexij_comp = (UINT4 *) Access_mmap(&new->indexij_comp_fd,&new->indexij_comp_len,indexij_compfile,/*randomp*/true);

      new->indexij_ptrs_access = MMAPPED;
      new->indexij_comp_access = MMAPPED;

    } else if (indexij_access == USE_ALLOCATE) {
      if (sharedp == true) {
	fprintf(stderr,"Allocating memory for indexij ptrs...");
	new->indexij_ptrs = (UINT4 *) Access_allocate_shared(&new->indexij_ptrs_access,&new->indexij_ptrs_shmid,&new->indexij_ptrs_key,
							     &new->indexij_ptrs_fd,&new->indexij_ptrs_len,&seconds,indexij_ptrsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->indexij_ptrs_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);

	fprintf(stderr,"Allocating memory for indexij comp...");
	new->indexij_comp = (UINT4 *) Access_allocate_shared(&new->indexij_comp_access,&new->indexij_comp_shmid,&new->indexij_comp_key,
							     &new->indexij_comp_fd,&new->indexij_comp_len,&seconds,indexij_compfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->indexij_comp_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      } else {
	fprintf(stderr,"Allocating memory for indexij ptrs...");
	new->indexij_ptrs = (UINT4 *) Access_allocate_private(&new->indexij_ptrs_access,&new->indexij_ptrs_len,&seconds,indexij_ptrsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->indexij_ptrs_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);

	fprintf(stderr,"Allocating memory for indexij comp...");
	new->indexij_comp = (UINT4 *) Access_allocate_private(&new->indexij_comp_access,&new->indexij_comp_len,&seconds,indexij_compfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->indexij_comp_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      }

    }

    new->indexsize = 3 + log4(((new->indexij_ptrs_len - 8)/sizeof(UINT4)/2)/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#endif
    new->indexspace = power(4,new->indexsize);

    if (lcp_access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading LCP/child/DC arrays...");
      new->lcpchilddc = (unsigned char *) Access_mmap_and_preload(&new->lcpchilddc_fd,&new->lcpchilddc_len,&npages,&seconds,
								  lcpchilddcfile,sizeof(unsigned char));
      new->lcpchilddc_access = MMAPPED;
      comma1 = Genomicpos_commafmt(new->lcpchilddc_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
      FREE(comma1);
    } else if (lcp_access == USE_MMAP_ONLY) {
      new->lcpchilddc = (unsigned char *) Access_mmap(&new->lcpchilddc_fd,&new->lcpchilddc_len,lcpchilddcfile,/*randomp*/true);
      new->lcpchilddc_access = MMAPPED;
    } else if (lcp_access == USE_ALLOCATE) {
      fprintf(stderr,"Allocating memory for lcpchildc...");
      if (sharedp == true) {
	new->lcpchilddc = (unsigned char *) Access_allocate_shared(&new->lcpchilddc_access,&new->lcpchilddc_shmid,&new->lcpchilddc_key,
								   &new->lcpchilddc_fd,&new->lcpchilddc_len,&seconds,lcpchilddcfile,sizeof(unsigned char));
      } else {
	new->lcpchilddc = (unsigned char *) Access_allocate_private(&new->lcpchilddc_access,&new->lcpchilddc_len,&seconds,lcpchilddcfile,sizeof(unsigned char));
      }
      comma1 = Genomicpos_commafmt(new->lcpchilddc_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
      FREE(comma1);
    }

    if (guideexc_access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading guide/exceptions...");
      new->lcp_guide = (UINT4 *) Access_mmap_and_preload(&new->lcp_guide_fd,&new->lcp_guide_len,&npages,&seconds,
							 lcp_guidefile,sizeof(UINT4));
      new->lcp_exceptions = (UINT4 *) Access_mmap_and_preload(&new->lcp_exceptions_fd,&new->lcp_exceptions_len,&npages,&seconds,
							 lcp_exceptionsfile,sizeof(UINT4));
      new->child_guide = (UINT4 *) Access_mmap_and_preload(&new->child_guide_fd,&new->child_guide_len,&npages,&seconds,
							 child_guidefile,sizeof(UINT4));
      new->child_exceptions = (UINT4 *) Access_mmap_and_preload(&new->child_exceptions_fd,&new->child_exceptions_len,&npages,&seconds,
							 child_exceptionsfile,sizeof(UINT4));
      new->lcp_guide_access = MMAPPED;
      new->lcp_exceptions_access = MMAPPED;
      new->child_guide_access = MMAPPED;
      new->child_exceptions_access = MMAPPED;
      fprintf(stderr,"done\n");

    } else if (guideexc_access == USE_MMAP_ONLY) {
      new->lcp_guide = (UINT4 *) Access_mmap(&new->lcp_guide_fd,&new->lcp_guide_len,
					     lcp_guidefile,/*randomp*/true);
      new->lcp_exceptions = (UINT4 *) Access_mmap(&new->lcp_exceptions_fd,&new->lcp_exceptions_len,
						  lcp_exceptionsfile,/*randomp*/true);
      new->child_guide = (UINT4 *) Access_mmap(&new->child_guide_fd,&new->child_guide_len,
					       child_guidefile,/*randomp*/true);
      new->child_exceptions = (UINT4 *) Access_mmap(&new->child_exceptions_fd,&new->child_exceptions_len,
							 child_exceptionsfile,/*randomp*/true);
      new->lcp_guide_access = MMAPPED;
      new->lcp_exceptions_access = MMAPPED;
      new->child_guide_access = MMAPPED;
      new->child_exceptions_access = MMAPPED;

    } else if (guideexc_access == USE_ALLOCATE) {
      fprintf(stderr,"Allocating memory for lcp guide...");
      if (sharedp == true) {
	new->lcp_guide = (UINT4 *) Access_allocate_shared(&new->lcp_guide_access,&new->lcp_guide_shmid,&new->lcp_guide_key,
							  &new->lcp_guide_fd,&new->lcp_guide_len,&seconds,lcp_guidefile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->lcp_guide_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      
	fprintf(stderr,"Allocating memory for lcp exceptions...");
	new->lcp_exceptions = (UINT4 *) Access_allocate_shared(&new->lcp_exceptions_access,&new->lcp_exceptions_shmid,&new->lcp_exceptions_key,
							       &new->lcp_exceptions_fd,&new->lcp_exceptions_len,&seconds,lcp_exceptionsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->lcp_exceptions_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
	
	fprintf(stderr,"Allocating memory for child guide...");
	new->child_guide = (UINT4 *) Access_allocate_shared(&new->child_guide_access,&new->child_guide_shmid,&new->child_guide_key,
							    &new->child_guide_fd,&new->child_guide_len,&seconds,child_guidefile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->child_guide_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);

	fprintf(stderr,"Allocating memory for child exceptions...");
	new->child_exceptions = (UINT4 *) Access_allocate_shared(&new->child_exceptions_access,&new->child_exceptions_shmid,&new->child_exceptions_key,
								 &new->child_exceptions_fd,&new->child_exceptions_len,&seconds,child_exceptionsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->child_exceptions_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      
      } else {
	new->lcp_guide = (UINT4 *) Access_allocate_private(&new->lcp_guide_access,&new->lcp_guide_len,&seconds,lcp_guidefile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->lcp_guide_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      
	fprintf(stderr,"Allocating memory for lcp exceptions...");
	new->lcp_exceptions = (UINT4 *) Access_allocate_private(&new->lcp_exceptions_access,&new->lcp_exceptions_len,&seconds,
								lcp_exceptionsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->lcp_exceptions_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
	
	fprintf(stderr,"Allocating memory for child guide...");
	new->child_guide = (UINT4 *) Access_allocate_private(&new->child_guide_access,&new->child_guide_len,&seconds,child_guidefile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->child_guide_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);

	fprintf(stderr,"Allocating memory for child exceptions...");
	new->child_exceptions = (UINT4 *) Access_allocate_private(&new->child_exceptions_access,&new->child_exceptions_len,&seconds,
								  child_exceptionsfile,sizeof(UINT4));
	comma1 = Genomicpos_commafmt(new->child_exceptions_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma1,seconds);
	FREE(comma1);
      }
    }

    new->n_lcp_exceptions = new->lcp_exceptions_len/(sizeof(UINT4) + sizeof(UINT4));
    new->child_guide_interval = 1024;
  }


  FREE(child_exceptionsfile);
  FREE(child_guidefile);

  FREE(lcp_exceptionsfile);
  FREE(lcp_guidefile);

  FREE(lcpchilddcfile);

#ifdef DEBUG15
  FREE(indexi_compfile);
  FREE(indexi_ptrsfile);
  FREE(indexj_compfile);
  FREE(indexj_ptrsfile);
  FREE(indexij_compfile);
  FREE(indexij_ptrsfile);
#elif defined(USE_SEPARATE_BUCKETS)
  FREE(indexi_compfile);
  FREE(indexi_ptrsfile);
  FREE(indexj_compfile);
  FREE(indexj_ptrsfile);
#else
  FREE(indexij_compfile);
  FREE(indexij_ptrsfile);
#endif

  FREE(sarrayfile);

  return new;
}


void
Sarray_free (T *old) {
  if (*old) {
#ifdef DEBUG15
    FREE((*old)->indexi_ptrs);
    FREE((*old)->indexi_comp);
    FREE((*old)->indexj_ptrs);
    FREE((*old)->indexj_comp);
    FREE((*old)->indexij_ptrs);
    FREE((*old)->indexij_comp);
#elif defined(USE_SEPARATE_BUCKETS)
    FREE((*old)->indexi_ptrs);
    FREE((*old)->indexi_comp);
    FREE((*old)->indexj_ptrs);
    FREE((*old)->indexj_comp);
#else
    if ((*old)->indexij_ptrs_access == MMAPPED) {
      munmap((void *) (*old)->indexij_ptrs,(*old)->indexij_ptrs_len);
      close((*old)->indexij_ptrs_fd);
    } else if ((*old)->indexij_ptrs_access == ALLOCATED_PRIVATE) {
      FREE((*old)->indexij_ptrs);
    } else if ((*old)->indexij_ptrs_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->indexij_ptrs,(*old)->indexij_ptrs_shmid,(*old)->indexij_ptrs_key);
    }
    if ((*old)->indexij_comp_access == MMAPPED) {
      munmap((void *) (*old)->indexij_comp,(*old)->indexij_comp_len);
      close((*old)->indexij_comp_fd);
    } else if ((*old)->indexij_comp_access == ALLOCATED_PRIVATE) {
      FREE((*old)->indexij_comp);
    } else if ((*old)->indexij_comp_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->indexij_comp,(*old)->indexij_comp_shmid,(*old)->indexij_comp_key);
    }
#endif

    if ((*old)->lcp_guide_access == MMAPPED) {
      munmap((void *) (*old)->lcp_guide,(*old)->lcp_guide_len);
      close((*old)->lcp_guide_fd);
    } else if ((*old)->lcp_guide_access == ALLOCATED_PRIVATE) {
      FREE((*old)->lcp_guide);
    } else if ((*old)->lcp_guide_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->lcp_guide,(*old)->lcp_guide_shmid,(*old)->lcp_guide_key);
    }

    if ((*old)->lcp_exceptions_access == MMAPPED) {
      munmap((void *) (*old)->lcp_exceptions,(*old)->lcp_exceptions_len);
      close((*old)->lcp_exceptions_fd);
    } else if ((*old)->lcp_exceptions_access == ALLOCATED_PRIVATE) {
      FREE((*old)->lcp_exceptions);
    } else if ((*old)->lcp_exceptions_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->lcp_exceptions,(*old)->lcp_exceptions_shmid,(*old)->lcp_exceptions_key);
    }

    if ((*old)->child_guide_access == MMAPPED) {
      munmap((void *) (*old)->child_guide,(*old)->child_guide_len);
      close((*old)->child_guide_fd);
    } else if ((*old)->child_guide_access == ALLOCATED_PRIVATE) {
      FREE((*old)->child_guide);
    } else if ((*old)->child_guide_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->child_guide,(*old)->child_guide_shmid,(*old)->child_guide_key);
    }

    if ((*old)->child_exceptions_access == MMAPPED) {
      munmap((void *) (*old)->child_exceptions,(*old)->child_exceptions_len);
      close((*old)->child_exceptions_fd);
    } else if ((*old)->child_exceptions_access == ALLOCATED_PRIVATE) {
      FREE((*old)->child_exceptions);
    } else if ((*old)->child_exceptions_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->child_exceptions,(*old)->child_exceptions_shmid,(*old)->child_exceptions_key);
    }

    if ((*old)->lcpchilddc_access == MMAPPED) {
      munmap((void *) (*old)->lcpchilddc,(*old)->lcpchilddc_len);
      close((*old)->lcpchilddc_fd);
    } else if ((*old)->lcpchilddc_access == ALLOCATED_PRIVATE) {
      FREE((*old)->lcpchilddc);
    } else if ((*old)->lcpchilddc_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->lcpchilddc,(*old)->lcpchilddc_shmid,(*old)->lcpchilddc_key);
    }

    if ((*old)->array_access == MMAPPED) {
      munmap((void *) (*old)->array,(*old)->array_len);
      close((*old)->array_fd);
    } else if ((*old)->array_access == ALLOCATED_PRIVATE) {
      FREE((*old)->array);
    } else if ((*old)->array_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->array,(*old)->array_shmid,(*old)->array_key);
    }

    FREE_KEEP(*old);
  }

  return;
}



#if 0
/* Old search method.  O(m*(log n)), where m is the querylength and n
   is the size of the suffix array searched */
static Sarrayptr_T
sarray_search_init (char *query, int querylength, int queryoffset, Compress_T query_compress, bool plusp,
		    Sarrayptr_T low, Sarrayptr_T high, Univcoord_T nmatches_low, Univcoord_T nmatches_high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  Univcoord_T nmatches_mid, fasti;
  char c;
  UINT4 sa_low, sa_mid;
  UINT4 lcp_low, lcp_mid;

  assert(querylength > 0);

  debug1(printf("sarray_search_init on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
#if 0
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif

    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

#ifdef WORDS_BIGENDIAN
    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/Bigendian_convert_uint(sarray->array[mid])-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand);
    pos = Bigendian_convert_uint(sarray->array[mid]) + fasti;
#else
    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand);
    pos = sarray->array[mid] + fasti;
#endif
    c = Genome_get_char_lex(genome,pos,sarray->n,chartable);

    if (fasti == (Univcoord_T) querylength || c > query[fasti]) {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
#ifdef WORDS_BIGENDIAN
      sa_mid = Bigendian_convert_uint(sarray->array[mid]);
#else
      sa_mid = sarray->array[mid];
#endif
      lcp_mid = Bitpack64_read_one(sa_mid,sarray->plcp_ptrs,sarray->plcp_comp) - sa_mid;
#ifdef USE_LCP
      if (lcp_mid != sarray->lcp[mid]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
    } else {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
#ifdef WORDS_BIGENDIAN
      sa_low = Bigendian_convert_uint(sarray->array[low]);
#else
      sa_low = sarray->array[low];
#endif
      lcp_low = Bitpack64_read_one(sa_low,sarray->plcp_ptrs,sarray->plcp_comp) - sa_low;
#ifdef USE_LCP
      if (lcp_low != sarray->lcp[low]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
    }

    debug1(printf("sarray_search_init with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_init ended.  Returning low %u+1\n\n",low));
  return low + 1;
}
#endif


#if 0
/* Old search method.  O(m*(log n)), where m is the querylength and n
   is the size of the suffix array searched */
static Sarrayptr_T
sarray_search_final (char *query, int querylength, int queryoffset, Compress_T query_compress, bool plusp,
		     Sarrayptr_T low, Sarrayptr_T high, Univcoord_T nmatches_low, Univcoord_T nmatches_high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  Univcoord_T nmatches_mid, fasti;
  UINT4 sa_low, sa_mid;
  UINT4 lcp_low, lcp_mid;
  char c;

  assert(querylength > 0);

  debug1(printf("sarray_search_final on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
#if 0
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif
    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

#ifdef WORDS_BIGENDIAN
    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/Bigendian_convert_uint(sarray->array[mid])-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand);
    pos = Bigendian_convert_uint(sarray->array[mid]) + fasti;
#else
    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand);
    pos = sarray->array[mid] + fasti;
#endif
    c = Genome_get_char_lex(genome,pos,sarray->n,chartable);

    if (fasti == (Univcoord_T) querylength || c < query[fasti]) {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
#ifdef WORDS_BIGENDIAN
      sa_low = Bigendian_convert_uint(sarray->array[low]);
#else
      sa_low = sarray->array[low];
#endif
      lcp_low = Bitpack64_read_one(sa_low,sarray->plcp_ptrs,sarray->plcp_comp) - sa_low;
#ifdef USE_LCP
      if (lcp_low != sarray->lcp[low]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
    } else {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
#ifdef WORDS_BIGENDIAN
      sa_mid = Bigendian_convert_uint(sarray->array[mid]);
#else
      sa_mid = sarray->array[mid];
#endif
      lcp_mid = Bitpack64_read_one(sa_mid,sarray->plcp_ptrs,sarray->plcp_comp) - sa_mid;
#ifdef USE_LCP
      if (lcp_mid != sarray->lcp[mid]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
    }

    debug1(printf("sarray_search_final with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_final ended.  Returning high %u-1\n\n",high-1));
  return high - 1;
}
#endif


int
nt_querylength (char *query, int querylength) {
  int i;
  char c;

  i = 0;
  while (i < querylength && ((c = query[i]) == 'A' || c == 'C' || c == 'G' || c == 'T')) {
    i++;
  }

  return i;
}


Oligospace_T
nt_oligo (char *query, int indexsize) {
  Oligospace_T oligo = 0U;
  int i;

  for (i = 0; i < indexsize; i++) {
    oligo *= 4;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default:
      fprintf(stderr,"Saw N in nt_oligo\n");
      abort();
    }
  }

  return oligo;
}

Oligospace_T
nt_oligo_truncate (char *query, int truncsize, int indexsize, int subst_value) {
  Oligospace_T oligo = 0U;
  int i;

  for (i = 0; i < truncsize; i++) {
    oligo *= 4;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default:
      fprintf(stderr,"Saw N in nt_oligo\n");
      abort();
    }
  }

  for ( ; i < indexsize; i++) {
    oligo *= 4;
    oligo += subst_value;
  }

  return oligo;
}



/* For child[index+1].up, just calling child[index] */
#define decode_up(index,child_bytes,child_guide,child_exceptions,child_guide_interval) index - Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval)
#define decode_down(index,child_bytes,child_guide,child_exceptions,child_guide_interval) Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval) + index + 1
#define decode_next(index,child_bytes,child_guide,child_exceptions,child_guide_interval) Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval) + index + 1

#if 0
/*                                      0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F */
static char discrim_char_before[16] = {'?','$','$','$','$','$','A','A','A','A','C','C','C','G','G','T'};
static char discrim_char_after[16]  = {'?','A','C','G','T','X','C','G','T','X','G','T','X','T','X','X'};
#endif

static bool
get_child_given_first (Sarrayptr_T *l, Sarrayptr_T *r, Sarrayptr_T i, Sarrayptr_T j, char desired_char,
		       T sarray, unsigned char *lcpchilddc, UINT4 lcp_whole, UINT4 nextl) {
  char c1, c2;
  UINT4 child_next;

  debug2(printf("Getting children for l-interval from %u to %u, char %c\n",i,j,desired_char));

#if 0
  /* First child already given */
  debug1(printf("lcp-interval %u..%u\n",i,j));
  up = decode_up(j,sarray->child_bytes,sarray->child_guide,sarray->child_exceptions,sarray->child_guide_interval);
  if (i < up && up <= j) {
    nextl = up;
    debug2(printf("nextl is up: %u\n",nextl));
  } else {
    nextl = decode_down(i,sarray->child_bytes,sarray->child_guide,sarray->child_exceptions,sarray->child_guide_interval); /* down */
    debug2(printf("nextl is down: %u\n",nextl));
  }
#endif

  /* Test first child: Use discrim_chars, rather than looking up S[SA[i] + lcp_whole] */
  c2 = Bytecoding_lcpchilddc_dc(&c1,nextl,lcpchilddc);
  debug2(printf("First child: %u to %u, discrim chars %c and %c\n",i,nextl-1,c1,c2));

  if (desired_char < c1) {
    debug2(printf("1.  Returning false, because desired %c < c1 %c\n",desired_char,c1));
    return false;
  } else if (desired_char == c1) {
    *l = i;
    *r = nextl - 1;
    debug2(printf("Returning true\n\n"));
    return true;
  } else if (desired_char < c2) {
    debug2(printf("1.  Returning false, because desired %c < c2 %c\n",desired_char,c2));
    return false;
  } else {
    /* Advance to middle children or final child */
    debug2(printf("1.  Advancing\n"));
  }

  /* Test for child[i] being down: lcp[child[i]] > lcp[i] */
  /* Test for child[i] being next_lindex: lcp[child[i]] == lcp[i] */
  /* Test middle children */
  while (nextl < j && Bytecoding_lcpchilddc_lcp_next(&child_next,nextl,/*bytes*/lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						     sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions) == lcp_whole) {
    /* Already tested for desired_char < c2 */
    if (desired_char == c2) {
      *l = nextl;
#if 0
      *r = Bytecoding_lcpchilddc_child_next(nextl,lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					    sarray->child_guide_interval) - 1; /* child[nextl] - 1 */
#else
      *r = child_next - 1;
#endif
      debug2(printf("Child: %u to %u, c2 %c\n",nextl,*r,c2));
      debug2(printf("Returning true\n\n"));
      return true;
    } else {
      debug2(printf("Child: %u",nextl));
#if 0
      nextl = Bytecoding_lcpchilddc_child_next(nextl,lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					       sarray->child_guide_interval); /* child[nextl] */
#else
      nextl = child_next;
#endif
      c2 = Bytecoding_lcpchilddc_dc(&c1,nextl,lcpchilddc);
      debug2(printf(" to %u, discrim chars %c and %c\n",nextl-1,c1,c2));

      if (desired_char < c2) {
	debug2(printf("M.  Returning false, because desired %c < c2 %c\n",desired_char,c2));
	return false;
      } else {
	debug2(printf("M.  Advancing\n"));
      }
    }
  }

  /* Test last child */
  /* Already tested for desired_char < c2 */
  debug2(printf("Final child: %u to %u, c2 %c\n",nextl,j,c2));
  if (desired_char == c2) {
    *l = nextl;
    *r = j;
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    debug2(printf("3.  Returning false, because desired %c != c2 %c\n",desired_char,c2));
    return false;
  }
}


static UINT4
find_longest_match (UINT4 nmatches, Sarrayptr_T *initptr, Sarrayptr_T *finalptr,
		    Sarrayptr_T i, Sarrayptr_T j, char *query, UINT4 querylength,
		    int queryoffset, Compress_T query_compress, T sarray, bool plusp,
		    int genestrand, char conversion[]) {
  UINT4 lcp_whole, nextl, up;
  UINT4 minlength;
  UINT4 l, r;
  Univcoord_T SA_i;

  while (nmatches < querylength) {
    if (i == j) {
      /* Singleton interval */
      debug1(printf("Singleton interval %u..%u\n",i,j));
      SA_i = Sarray_position(sarray,i);
      nmatches +=
	Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
					     /*pos5*/queryoffset+nmatches,/*pos3*/queryoffset+querylength,
					     plusp,genestrand);
      *initptr = i;
      *finalptr = j;
      return nmatches;

    } else {
      /* First child */
      debug1(printf("lcp-interval %u..%u\n",i,j));
      up = Bytecoding_lcpchilddc_child_up(j,sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					  sarray->child_guide_interval);
      if (i < up && up <= j) {
	nextl = up;
	debug2(printf("nextl is up: %u\n",nextl));
      } else {
	nextl = Bytecoding_lcpchilddc_child_next(i,sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval); /* really down */
	debug2(printf("nextl is down: %u\n",nextl));
      }

      lcp_whole = Bytecoding_lcpchilddc_lcp(nextl,sarray->lcpchilddc,sarray->lcp_exceptions,
					    sarray->n_lcp_exceptions); /* lcp(i,j) */
      debug1(printf("lcp_whole for %u..%u is %d, compared with nmatches %d\n",i,j,lcp_whole,nmatches));

      if (lcp_whole > nmatches) {
	/* Check only up to minlength, so we validate the entire interval */
	minlength = (lcp_whole < querylength) ? lcp_whole : querylength;
	debug1(printf("Looking up genome for query from %d .. %d - 1\n",nmatches,minlength));
	SA_i = Sarray_position(sarray,i);
	nmatches +=
	  Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
					       /*pos5*/queryoffset+nmatches,/*pos3*/queryoffset+minlength,
					       plusp,genestrand);
	if (nmatches < minlength) {
	  *initptr = i;
	  *finalptr = j;
	  return nmatches;

	} else if (nmatches >= querylength) {
	  debug1(printf("nmatches is now %d >= querylength %d => success\n",nmatches,querylength));
	  *initptr = i;
	  *finalptr = j;
	  return nmatches;
	}
      }
	
      debug1(printf("nmatches is now %d => desired_char is %c => %c\n",
		    nmatches,query[nmatches],conversion[query[nmatches]]));
      if (get_child_given_first(&l,&r,i,j,/*desired_char*/conversion[(int) query[nmatches]],
				sarray,sarray->lcpchilddc,lcp_whole,nextl) == false) {
	*initptr = i;
	*finalptr = j;
	return nmatches;
      } else {
	nmatches += 1;
	i = l;
	j = r;
      }
    }
  }

  *initptr = i;
  *finalptr = j;
  return nmatches;
}



/* Searches using LCP and child arrays.  Should be O(m * |Sigma|),
   where m wis the querylength and |Sigma| is the size of the alphabet
   (4 for DNA) */
/* query is a substring of the original, starting with queryoffset */
void
Sarray_read (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, bool *successp,
	     UINT4 *nmatches, char *query, UINT4 querylength, int queryoffset,
	     Compress_T query_compress, T sarray, bool plusp, int genestrand,
	     char conversion[]) {
  int effective_querylength;	/* length to first N */
  Oligospace_T oligo;
  UINT4 l, r;

#ifdef DEBUG1
  Univcoord_T SA_i, hit, child_next;
  int k = 0;
  UINT4 recount, lcp_prev, lcp_next, lcp_i, max_lcp;
  char Buffer[1000+1], c1, c2;
  bool failp;
#endif

  debug1(printf("Sarray_read on %.*s, querylength %d, plusp %d\n",querylength,query,querylength,plusp));

  /* Find initial lcp-interval */
  effective_querylength = nt_querylength(query,querylength);
  debug1(printf("sarray_search on %.*s, querylength %d, effective querylength %d, plusp %d\n",
		querylength,query,querylength,effective_querylength,plusp));


  *nmatches = 0;
  if (effective_querylength == 0) {
    *initptr = *finalptr = 0;
    *successp = false;
    return;

  } else if (effective_querylength < sarray->indexsize) {
    debug1(printf("string %.*s with effective querylength %d is shorter than indexsize",
		  querylength,query,effective_querylength));
    l = 1;
    r = sarray->n;

  } else {
    oligo = nt_oligo(query,sarray->indexsize);
#ifdef DEBUG15
    if ((l = Bitpack64_read_two(&r,oligo*2,sarray->indexij_ptrs,sarray->indexij_comp)) !=
	Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp)) {
      abort();
    } else if (r - 1 != Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp)) {
      printf("For oligo %u, separate buckets give %u and %u, while single bucket gives %u and %u\n",
	     oligo,
	     Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp),
	     Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp),
	     l,r);
      abort();
    }
    r--;			/* Because interleaved writes r+1 to maintain monotonicity */
#elif defined(USE_SEPARATE_BUCKETS)
    l = Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp);
    r = Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp);
#else
    l = Bitpack64_read_two(&r,oligo*2,sarray->indexij_ptrs,sarray->indexij_comp);
    r--;			/* Because interleaved writes r+1 to maintain monotonicity */
#endif
    debug1(printf("string %.*s is equal/longer than indexsize %d => oligo %u => interval %u..%u",
		  querylength,query,sarray->indexsize,oligo,l,r));
    if (l <= r) {
      debug1(printf(" (good)\n"));
      *nmatches = sarray->indexsize;
      /* i = l; */
      /* j = r; */
    } else {
      /* The entire lcp-interval [1,sarray->n] should also work without initindex */
      l = 1;
      r = sarray->n;
      debug1(printf(" (bad) => entire lcp-interval: %u..%u\n",l,r));
    }
  }

  if (l > r) {
    /* Did not find a match using saindex or one letter */
    *initptr = l;
    *finalptr = r;
  } else {
    *nmatches = find_longest_match(*nmatches,&(*initptr),&(*finalptr),/*i*/l,/*j*/r,
				   query,querylength,queryoffset,query_compress,sarray,
				   plusp,genestrand,conversion);
  }

  /* Search through suffix tree */
  debug1(printf("initptr gets %u, finalptr gets %u\n",*initptr,*finalptr));

  if (*nmatches < querylength) {
    *successp = false;
    debug1(printf("%s fail at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  } else {
    *successp = true;
    debug1(printf("%s success at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  }

#ifdef DEBUG1
  failp = false;

  /* Before */
  if (*nmatches > 0 && *initptr > 0U) {
    SA_i = Sarray_position(sarray,(*initptr)-1);
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand);
    printf("%d\t%u\t%u\t",recount,(*initptr)-1,SA_i/*+ 1U*/);
    c2 = Bytecoding_lcpchilddc_dc(&c1,(*initptr)-1,sarray->lcpchilddc);
    printf("%c%c\t",c1,c2);
    lcp_i = Bytecoding_lcpchilddc_lcp((*initptr)-1,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    printf("%u\t",lcp_i);
    lcp_next = Bytecoding_lcpchilddc_lcp((*initptr),/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    printf("%u\t",Bytecoding_lcpchilddc_lcp_next(&child_next,(*initptr)-1,/*bytes*/sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(SA_i,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(SA_i,recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(SA_i,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(SA_i,recount+1,Buffer);
      }
    }
    printf("%s\n",Buffer);
    if (recount >= *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false negative: recount %d at %u before init does equal expected nmatches %d\n",
	     recount,SA_i,*nmatches);
      failp = true;
    }
  }
  printf("\n");


  /* Hits */
  lcp_prev = lcp_i;
  for (k = 0; k < (int) (*finalptr - *initptr + 1) && k < MAX_DEBUG1_HITS; k++) {
    SA_i = Sarray_position(sarray,(*initptr)+k);
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand);
    printf("%d\t%u\t%u\t",recount,(*initptr)+k,SA_i/*+ 1U*/);
    c2 = Bytecoding_lcpchilddc_dc(&c1,(*initptr)+k,sarray->lcpchilddc);
    printf("%c%c\t",c1,c2);
    lcp_i = Bytecoding_lcpchilddc_lcp((*initptr)+k,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    lcp_next = Bytecoding_lcpchilddc_lcp((*initptr)+k+1,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    printf("%u\t",lcp_i);
    printf("%u\t",Bytecoding_lcpchilddc_lcp_next(&child_next,(*initptr)+k,/*bytes*/sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    max_lcp = lcp_i;
    if (lcp_prev > max_lcp) {
      max_lcp = lcp_prev;
    }
    if (lcp_next > max_lcp) {
      max_lcp = lcp_next;
    }
    if (max_lcp > 1000) {
      max_lcp = 1000;
    }

    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(SA_i,max_lcp+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(SA_i,max_lcp+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(SA_i,max_lcp+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(SA_i,max_lcp+1,Buffer);
      }
    }
    printf("%s\n",Buffer);
    if (recount != *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false positive: recount %d at %u does not equal expected nmatches %d\n",
	     recount,Sarray_position(sarray,(*initptr)),*nmatches);
      failp = true;
    }

    lcp_prev = lcp_i;
  }

  if (k < (int) (*finalptr - *initptr + 1)) {
    /* Overflow */
    printf("...\n");
    k = (int) (*finalptr - *initptr);
    hit = Sarray_position(sarray,(*initptr)+k);
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/hit-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand);
    printf("%d\t%u\t%u\t",recount,(*initptr)+k,hit /*+ 1U*/);
    c2 = Bytecoding_lcpchilddc_dc(&c1,(*initptr)+k,sarray->lcpchilddc);
    printf("%c%c\t",c1,c2);
    lcp_i = Bytecoding_lcpchilddc_lcp((*initptr)+k,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    lcp_next = Bytecoding_lcpchilddc_lcp((*initptr)+k+1,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions);
    printf("%u\t",lcp_i);
    printf("%u\t",Bytecoding_lcpchilddc_lcp_next(&child_next,(*initptr)+k,/*bytes*/sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(hit,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(hit,recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(hit,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(hit,recount+1,Buffer);
      }
    }
    printf("%s\n",Buffer);
    if (recount != *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false positive: recount %d at %u does not equal expected nmatches %d\n",
	     recount,Sarray_position(sarray,*initptr),*nmatches);
      failp = true;
    }
    /* hits[k] = sarray->array[(*initptr)++]; */
  }


  /* After */
  if (*nmatches > 0 && (SA_i = Sarray_position(sarray,(*finalptr)+1)) > 0U) {
    printf("\n");
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/SA_i-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand);
    printf("%d\t%u\t%u\t",recount,(*finalptr)+1,SA_i/*+ 1U*/);
    c2 = Bytecoding_lcpchilddc_dc(&c1,(*finalptr)+1,sarray->lcpchilddc);
    printf("%c%c\t",c1,c2);
    printf("%u\t",Bytecoding_lcpchilddc_lcp((*finalptr)+1,/*bytes*/sarray->lcpchilddc,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    printf("%u\t",Bytecoding_lcpchilddc_lcp_next(&child_next,(*finalptr)+1,/*bytes*/sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions));
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(SA_i,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(SA_i,recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(SA_i,recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(SA_i,recount+1,Buffer);
      }
    }
    printf("%s\n",Buffer);
    if (recount >= *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false negative: recount %d at %u after (*finalptr) does equal expected nmatches %d\n",
	     recount,SA_i,*nmatches);
      failp = true;
    }
  }

  if (failp == true) {
    /* Can happen because $ ranks below 0 */
    /* Can also happen with CMET or ATOI, since genome128_hr procedures find genome-to-query mismatches */
    /* abort(); */
  }
#endif

  return;
}


Univcoord_T *
Sarray_lookup (int *nhits, T sarray, char *query, UINT4 querylength, int queryoffset,
	       Compress_T query_compress, bool plusp, int genestrand,
	       char conversion[]) {
  Univcoord_T *hits;
  Sarrayptr_T initptr, finalptr, ptr;
  bool successp;
  UINT4 nmatches;
  int k;

  Sarray_read(&initptr,&finalptr,&successp,&nmatches,query,querylength,queryoffset,
	      query_compress,sarray,plusp,genestrand,conversion);

  
  if (successp == false) {
    *nhits = 0;
    return (Univcoord_T *) NULL;
  } else {
    hits = (Univcoord_T *) MALLOC((finalptr - initptr + 1)*sizeof(UINT4));
    k = 0;
    for (ptr = initptr; ptr <= finalptr; ptr++) {
      hits[k++] = Sarray_position(sarray,ptr);
    }
    *nhits = k;
    return hits;
  }
}



