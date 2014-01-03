static char rcsid[] = "$Id: sarray-write.c 116707 2013-11-27 19:41:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sarray-write.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */
#include "bool.h"
#include "mem.h"
#include "genomicpos.h"
#include "assert.h"
#include "compress.h"
#include "bitpack64-write.h"
#include "bitpack64-access.h"
#include "fopen.h"
#include "saca-k.h"
#include "genome_hr.h"
#include "popcount.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


/* #define WRITE_LCP 1 */

/* make_index */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* correctness of using Genome_consecutive_matches_pair */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif


/* For computing LCP.  Comment out because mmap is faster than fread */
/* #define READ_SA_FROM_FILE 1 */

#define MONITOR_INTERVAL 100000000 /* 100 million nt */


static void
compute_lcp (UINT4 *lcp,
#ifdef DEBUG14
	     unsigned char *s,
#endif
	     UINT4 *SA, UINT4 n) {
  UINT4 *rank, h;
  UINT4 i, j;
  char *comma;
#ifdef DEBUG14
  UINT4 horig;
#endif

  rank = (UINT4 *) CALLOC(n+1,sizeof(UINT4));
  for (i = 0; i <= n; i++) {
    rank[SA[i]] = i;
  }

  lcp[0] = 0;
  h = 0;
  for (i = 0; i <= n; i++) {
    if (rank[i] > 0) {
      j = SA[rank[i] - 1];
#ifdef DEBUG14
      horig = h;
      while (i + h < n && j + h < n && s[i+h] == s[j+h]) {
	h++;
      }
      if ((h - horig) != Genome_consecutive_matches_pair(i+horig,j+horig,/*genomelength*/n)) {
	abort();
      }
#else
      h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
#endif
      lcp[rank[i]] = h;
      if (h > 0) {
	h--;
      }
    }
    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing lcp index %s\n",comma);
      FREE(comma);
    }
  }

  FREE(rank);

  return;
}


/* Computes permuted lcp (Karkkainen, CPM 2009) */
static void
compute_plcp (UINT4 *plcp,
#ifdef READ_SA_FROM_FILE
	      FILE *fp,
#else
	      UINT4 *SA,
#endif
	      UINT4 n) {
  UINT4 *phi, h;
  UINT4 i, j;
  char *comma;
#ifdef READ_SA_FROM_FILE
  UINT4 sa, prev_sa;
#endif

  phi = plcp;			/* Use space allocated for plcp */

#ifdef READ_SA_FROM_FILE
  fprintf(stderr,"Inverting suffix array (via file)...");
  FREAD_UINT(&prev_sa,fp);
  for (i = 1; i <= n; i++) {
    FREAD_UINT(&sa,fp);
    phi[sa] = prev_sa;
    prev_sa = sa;
  }
#else
  fprintf(stderr,"Inverting suffix array (via mmap)...");
  for (i = 1; i <= n; i++) {
    phi[SA[i]] = SA[i-1];
  }
#endif
  fprintf(stderr,"done\n");
  /* Note that phi[n] is not assigned, because SA[i] == n for i == 0, and we don't look up i == 0 */


  h = 0;
  for (i = 0; i < n; i++) {
    j = phi[i];			/* To be overwritten by plcp[i] */
    h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
    plcp[i] = h;				 /* overwrites phi[i] */
    if (h > 0) {
      h--;
    }

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing permuted lcp index %s\n",comma);
      FREE(comma);
    }
  }

  /* This makes lcp[0] = 0, because lcp[0] = plcp[SA[0]] = plcp[n] */
  plcp[n] = 0;

  return;
}



static UINT4
power (int base, int exponent) {
  UINT4 result = 1U;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


#define LOW_TWO_BITS 0x3
#define RIGHT_A 0
#define RIGHT_C 1
#define RIGHT_G 2
#define RIGHT_T 3


static void
oligo_nt (char *nt, UINT4 oligo, int oligosize) {
  int i, j;
  UINT4 lowbits;

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

  return;
}


#if 0
/* May not be correct */
static Sarrayptr_T
sarray_search_nolcp_init (char *query, int querylength, Genome_T genomecomp, UINT4 *SA, Sarrayptr_T low, Sarrayptr_T high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  char c;
  int i;

  assert(querylength > 0);

  debug1(printf("sarray_search_nolcp_init on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));

    pos = SA[mid];

    i = 0;
    while (i < querylength && (c = Genome_get_char(genomecomp,pos)) == query[i]) {
      i++;
      pos++;
    }
    if (c == 'N') {
      c = 'X';
    }

    if (i == (Univcoord_T) querylength) {
      return low + 1;
    } else if (c > query[i]) {
      high = mid;
    } else {
      low = mid;
    }

    debug1(printf("sarray_search_nolcp_init with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_nolcp_init ended.  Returning low %u+1\n\n",low));
  return low + 1;
}
#endif


#if 0
/* May not be correct */
static Sarrayptr_T
sarray_search_nolcp_final (char *query, int querylength, Genome_T genomecomp, UINT4 *SA, Sarrayptr_T low, Sarrayptr_T high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  char c;
  int i;

  assert(querylength > 0);

  debug1(printf("sarray_search_nolcp_final on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));

    pos = SA[mid];

    i = 0;
    while (i < querylength && (c = Genome_get_char(genomecomp,pos)) == query[i]) {
      i++;
      pos++;
    }
    if (c == 'N') {
      c = 'X';
    }

    if (i == (Univcoord_T) querylength) {
      return high - 1;
    } else if (c < query[i]) {
      low = mid;
    } else {
      high = mid;
    }

    debug1(printf("sarray_search_nolcp_final with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_nolcp_final ended.  Returning high %u-1\n\n",high-1));
  return high - 1;
}
#endif


static UINT4
sarray_search_init (char *query, int querylength, Genome_T genomecomp, UINT4 *SA,
		    UINT4 low, UINT4 high, UINT4 nmatches_low, UINT4 nmatches_high) {
  UINT4 mid, pos;
  UINT4 nmatches_mid;
  int fasti;
  char c;
  UINT4 lcp_low, lcp_mid;

  debug1(printf("sarray_search_init on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
    /* Compute mid for UINT4s */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid;
    pos = SA[mid] + nmatches_mid;
    while (fasti < querylength && (c = Genome_get_char(genomecomp,pos)) == query[fasti]) {
      debug1(printf("Comparing query %d with pos %u\n",fasti,pos));
      fasti++;
      pos++;
    }
    if (c == 'N') {
      c = 'X';
    }

    if (fasti == querylength || c > query[fasti]) {
      high = mid;
      lcp_mid = Bitpack64_access(mid);
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
    } else {
      low = mid;
      lcp_low = Bitpack64_access(low);
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
    }

    debug1(printf("sarray_search_init with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_init ended.  Returning low %u+1\n\n",low));
  return low + 1;
}


static UINT4
sarray_find_index (bool *successp, char *query, int querylength, Genome_T genomecomp,
		   UINT4 *SA, UINT4 n) {
  UINT4 low, high, mid, pos;
  UINT4 nmatches_low, nmatches_high, nmatches_mid;
  UINT4 lcp_low, lcp_mid;

  UINT4 prevlow, prevhigh;
  int nmatches_prevlow, nmatches_prevhigh, nmatches_best = 0;

  int fasti;
  char c;

  *successp = false;
  low = prevlow = 0;
  high = prevhigh = n;
  nmatches_low = nmatches_high = 0;

  debug1(printf("sarray_search on %s, querylength %d, with low %u, high %u\n",
	       query,querylength,low,high));
  while (low + 1 < high && *successp == false) {
    /* Compute mid for UINT4s */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    nmatches_mid = (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid;
    pos = SA[mid] + nmatches_mid;
    while (fasti < querylength && (c = Genome_get_char(genomecomp,pos)) == query[fasti]) {
      debug1(printf("Comparing query %d with pos %u\n",fasti,pos));
      fasti++;
      pos++;
    }
    if (c == 'N') {
      c = 'X';
    }

    if (fasti > nmatches_best) {
      debug1(printf("fasti %d > nmatches_best %d.  Saving prevlow %u and prevhigh %u.\n",
		    fasti,nmatches_best,low,high));
      prevlow = low;
      prevhigh = high;
      nmatches_prevlow = nmatches_low;
      nmatches_prevhigh = nmatches_high;
      nmatches_best = fasti;
    }

    if (fasti == querylength) {
      *successp = true;

    } else if (c < query[fasti]) {
      low = mid;
      lcp_low = Bitpack64_access(low);
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
      debug1(printf("genome %c < query (%c) => low gets %u @ %u\n",c,query[fasti],low,SA[low]));

    } else if (c > query[fasti]) {
      high = mid;
      lcp_mid = Bitpack64_access(mid);
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
      debug1(printf("genome %c > query (%c) => high gets %u @ %u\n",c,query[fasti],high,SA[high]));

    } else {
      debug1(printf("genome %c == query (%c) => should not happen after Genome_consecutive_matches\n",
		   c,query[fasti]));
      abort();
    }

    debug1(printf("sarray_search with low %u @ %u, high %u @ %u\n",low,SA[low],high,SA[high]));
  }
  debug1(printf("\n"));

  if (*successp == false) {
    return -1U;
  } else {
    return sarray_search_init(query,querylength,genomecomp,SA,
			      low,mid,nmatches_low,nmatches_mid);
  }
}



static UINT4 *
make_index (UINT4 *oligospace, int querylength, Genome_T genomecomp, UINT4 *SA,
	    UINT4 n) {
  UINT4 *saindex;
  char *queryuc_ptr;
  UINT4 i;
  bool successp;

  *oligospace = power(4,querylength);
  saindex = (UINT4 *) CALLOC((*oligospace)+1,sizeof(UINT4));

  queryuc_ptr = (char *) CALLOC(querylength+1,sizeof(char));

  saindex[0] = 1U;
  for (i = 1; i < *oligospace; i++) {
    oligo_nt(queryuc_ptr,i,querylength);
#if 0
    if ((ptr = sarray_find_index(&successp,queryuc_ptr,querylength,genomecomp,SA,n)) == -1) {
      saindex[i] = saindex[i-1];
    } else {
      saindex[i] = ptr;
    }
#else
    saindex[i] = sarray_find_index(&successp,queryuc_ptr,querylength,genomecomp,SA,n);
#endif
    /* printf("%s\t%u\t%u\n",queryuc_ptr,i,saindex[i]); */
  }
  FREE(queryuc_ptr);

  saindex[*oligospace] = n;

  return saindex;
}


static int
compute_packsize (UINT4 *values) {
  UINT4 packsize;
  UINT4 maxvalue = 0, top;
  int i;
  int firstbit, msb;

  for (i = 0; i < 64; i++) {
    maxvalue |= values[i];
  }

#ifdef HAVE_BUILTIN_CLZ
  firstbit = __builtin_clz(maxvalue);
  packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
  asm("bsr %1,%0" : "=r"(msb) : "r"(maxvalue));
  packsize = msb + 1;
#else
  firstbit = ((top = maxvalue >> 16) ? clz_table[top] : 16 + clz_table[maxvalue]);
  packsize = 32 - firstbit;
#endif

#ifdef ALLOW_ODD_PACKSIZES
  return packsize;
#else
  return (packsize + 1) & ~1;	/* Converts packsizes to the next multiple of 2 */
#endif
}

#define LCP_BUFFER_SIZE 1000000
#define METAINFO_SIZE 1
#define BLOCKSIZE 64

#if 0
/* Uses minimal memory, but too slow.  Also, there is an error in the last block. */
void
Sarray_write_lcpptrs_naive (char *lcpptrsfile, char *lcpcompfile, UINT4 *SA,
			    UINT4 *origlcp, UINT4 /*n*/genomelength) {
  UINT4 lcp[65];
  FILE *lcpptrs_fp, *lcpcomp_fp;
  UINT4 *lcpptrs;
  int lcpptri, i;
  UINT4 positioni;

  UINT4 *lcp_buffer;
  int lcp_buffer_size = LCP_BUFFER_SIZE;
  int lcp_buffer_i;

  UINT4 nwritten;
  int packsize;


  Bitpack64_write_setup();

  /* 1 metavalue: nwritten (pointer).  Packsize can be
     computed from difference between successive pointers, if only
     even packsizes are allowed */
  lcpptrs = (UINT4 *) CALLOC((((genomelength + 1) + BLOCKSIZE - 1)/BLOCKSIZE + 1) * METAINFO_SIZE,sizeof(UINT4));
  lcpptri = 0;

  if ((lcpcomp_fp = FOPEN_WRITE_BINARY(lcpcompfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpcompfile);
    exit(9);
  }
  lcp_buffer = (UINT4 *) CALLOC(lcp_buffer_size,sizeof(UINT4));
  lcp_buffer_i = 0;

  nwritten = 0U;

  /* Last value of lcp is lcp[genomelength], with lcp[0] = 0. */
  lcp[64] = 0;
  for (positioni = 0; positioni + BLOCKSIZE <= genomelength; positioni += BLOCKSIZE) {
    /* Pointer */
    lcpptrs[lcpptri++] = nwritten;

    /* Pack block of 64 diffs */
    lcp[0] = lcp[64];
    for (i = 1; i <= BLOCKSIZE; i++) {
      lcp[i] = Genome_consecutive_matches_pair(SA[positioni+i-1],SA[positioni+i],genomelength);
      if (lcp[i] != origlcp[positioni+i]) {
	fprintf(stderr,"At index %u, positions %u and %u, origlcp = %d, new lcp = %d\n",
		positioni+i,SA[positioni+i-1],SA[positioni+i],origlcp[positioni+i],lcp[i]);
	exit(9);
      }
    }
    packsize = compute_packsize(lcp);
    lcp_buffer_i = Bitpack64_write_vert(lcpcomp_fp,lcp_buffer,lcp_buffer_size,
					lcp_buffer_i,lcp,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  if (positioni <= genomelength) {
    /* Finish last block of 64 */
    lcpptrs[lcpptri++] = nwritten;
    
    /* Note: The code below is off by 1 position.  What is assigned to lcp[1] should go into lcp[0]. */
    lcp[0] = lcp[64];

    i = 1;
    while (positioni < genomelength) {
      lcp[i] = Genome_consecutive_matches_pair(SA[positioni-1],SA[positioni],genomelength);
      if (lcp[i] != origlcp[positioni]) {
	fprintf(stderr,"At index %u, positions %u and %u, origlcp = %d, new lcp = %d\n",
		positioni,SA[positioni-1],SA[positioni],origlcp[positioni],lcp[i]);
	exit(9);
      }
      positioni++;
      i++;
    }

    while (i < BLOCKSIZE) {
      lcp[i++] = 0;
    }

    packsize = compute_packsize(lcp);
    lcp_buffer_i = Bitpack64_write_vert(lcpcomp_fp,lcp_buffer,lcp_buffer_size,
					lcp_buffer_i,lcp,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  /* Write the final pointer, which will point after the end of the
     file */
  lcpptrs[lcpptri++] = nwritten;

  if ((lcpptrs_fp = FOPEN_WRITE_BINARY(lcpptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(lcpptrs,lcpptri,lcpptrs_fp);
    FREE(lcpptrs);
    fclose(lcpptrs_fp);
  }
    
  /* Empty lcp buffer */
  if (lcp_buffer_i > 0) {
    FWRITE_UINTS(lcp_buffer,lcp_buffer_i,lcpcomp_fp);	
    lcp_buffer_i = 0;
  }
  FREE(lcp_buffer);
  fclose(lcpcomp_fp);

  return;
}
#endif



static void
Sarray_write_lcpptrs (char *lcpptrsfile, char *lcpcompfile, UINT4 *lcp,
		      UINT4 /*n*/genomelength) {
  FILE *lcpptrs_fp, *lcpcomp_fp;
  UINT4 *lcpptrs;
  int lcpptri, i;
  UINT4 positioni;

  UINT4 *lcp_buffer;
  int lcp_buffer_size = LCP_BUFFER_SIZE;
  int lcp_buffer_i;

  UINT4 values[64];

  UINT4 nwritten;
  int packsize;


  Bitpack64_write_setup();

  /* 1 metavalue: nwritten (pointer).  Packsize can be
     computed from difference between successive pointers, if only
     even packsizes are allowed */
  lcpptrs = (UINT4 *) CALLOC((((genomelength + 1) + BLOCKSIZE - 1)/BLOCKSIZE + 1) * METAINFO_SIZE,sizeof(UINT4));
  lcpptri = 0;

  if ((lcpcomp_fp = FOPEN_WRITE_BINARY(lcpcompfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpcompfile);
    exit(9);
  }
  lcp_buffer = (UINT4 *) CALLOC(lcp_buffer_size,sizeof(UINT4));
  lcp_buffer_i = 0;

  nwritten = 0U;

  /* Last value of lcp is lcp[genomelength], with lcp[0] = 0. */
  for (positioni = 0; positioni + BLOCKSIZE <= genomelength; positioni += BLOCKSIZE) {
    /* Pointer */
    lcpptrs[lcpptri++] = nwritten;

    /* Pack block of 64 diffs */
    packsize = compute_packsize(&(lcp[positioni]));
    lcp_buffer_i = Bitpack64_write_vert(lcpcomp_fp,lcp_buffer,lcp_buffer_size,
					lcp_buffer_i,&(lcp[positioni]),packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  if (positioni <= genomelength) {
    /* Finish last block of 64 */
    lcpptrs[lcpptri++] = nwritten;
    
    i = 0;
    while (positioni <= genomelength) {
      values[i++] = lcp[positioni++];
    }
    while (i < BLOCKSIZE) {
      values[i++] = 0;
    }

    packsize = compute_packsize(values);
    lcp_buffer_i = Bitpack64_write_vert(lcpcomp_fp,lcp_buffer,lcp_buffer_size,
					lcp_buffer_i,values,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  /* Write the final pointer, which will point after the end of the
     file */
  lcpptrs[lcpptri++] = nwritten;

  if ((lcpptrs_fp = FOPEN_WRITE_BINARY(lcpptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(lcpptrs,lcpptri,lcpptrs_fp);
    FREE(lcpptrs);
    fclose(lcpptrs_fp);
  }
    
  /* Empty lcp buffer */
  if (lcp_buffer_i > 0) {
    FWRITE_UINTS(lcp_buffer,lcp_buffer_i,lcpcomp_fp);	
    lcp_buffer_i = 0;
  }
  FREE(lcp_buffer);
  fclose(lcpcomp_fp);

  return;
}


static void
Sarray_write_lcpptrs_using_plcp (char *lcpptrsfile, char *lcpcompfile, UINT4 *plcp,
#ifdef READ_SA_FROM_FILE
				 FILE *fp,
#else
				 UINT4 *SA,
#endif
				 UINT4 /*n*/genomelength) {
  FILE *lcpptrs_fp, *lcpcomp_fp;
  UINT4 *lcpptrs;
  int lcpptri, i;
  UINT4 positioni;

  UINT4 *lcp_buffer;
  int lcp_buffer_size = LCP_BUFFER_SIZE;
  int lcp_buffer_i;

  UINT4 values[BLOCKSIZE];

  UINT4 nwritten;
  int packsize;
#ifdef READ_SA_FROM_FILE
  UINT4 sa;
#endif


  Bitpack64_write_setup();

  /* 1 metavalue: nwritten (pointer).  Packsize can be
     computed from difference between successive pointers, if only
     even packsizes are allowed */
  lcpptrs = (UINT4 *) CALLOC((((genomelength + 1) + BLOCKSIZE - 1)/BLOCKSIZE + 1) * METAINFO_SIZE,sizeof(UINT4));
  lcpptri = 0;

  if ((lcpcomp_fp = FOPEN_WRITE_BINARY(lcpcompfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpcompfile);
    exit(9);
  }
  lcp_buffer = (UINT4 *) CALLOC(lcp_buffer_size,sizeof(UINT4));
  lcp_buffer_i = 0;

  nwritten = 0U;

  /* Last value of lcp is lcp[genomelength], with lcp[0] = 0. */
  for (positioni = 0; positioni + BLOCKSIZE <= genomelength; positioni += BLOCKSIZE) {
    /* Pointer */
    lcpptrs[lcpptri++] = nwritten;

#ifdef READ_SA_FROM_FILE
    for (i = 0; i < BLOCKSIZE; i++) {
      FREAD_UINT(&sa,fp);
      values[i] = plcp[sa];
    }
#else
    for (i = 0; i < BLOCKSIZE; i++) {
      values[i] = plcp[SA[positioni+i]];
    }
#endif

    /* Pack block of 64 diffs */
    packsize = compute_packsize(values);
    lcp_buffer_i = Bitpack64_write_vert(lcpcomp_fp,lcp_buffer,lcp_buffer_size,
					lcp_buffer_i,values,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  if (positioni <= genomelength) {
    /* Finish last block of 64 */
    lcpptrs[lcpptri++] = nwritten;
    
    i = 0;
#ifdef READ_SA_FROM_FILE
    while (positioni <= genomelength) {
      FREAD_UINT(&sa,fp);
      values[i++] = plcp[sa];
      positioni++;
    }
#else
    while (positioni <= genomelength) {
      values[i++] = plcp[SA[positioni++]];
    }
#endif
    while (i < BLOCKSIZE) {
      values[i++] = 0;
    }

    packsize = compute_packsize(values);
    lcp_buffer_i = Bitpack64_write_vert(lcpcomp_fp,lcp_buffer,lcp_buffer_size,
					lcp_buffer_i,values,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  /* Write the final pointer, which will point after the end of the
     file */
  lcpptrs[lcpptri++] = nwritten;

  if ((lcpptrs_fp = FOPEN_WRITE_BINARY(lcpptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(lcpptrs,lcpptri,lcpptrs_fp);
    FREE(lcpptrs);
    fclose(lcpptrs_fp);
  }
    
  /* Empty lcp buffer */
  if (lcp_buffer_i > 0) {
    FWRITE_UINTS(lcp_buffer,lcp_buffer_i,lcpcomp_fp);	
    lcp_buffer_i = 0;
  }
  FREE(lcp_buffer);
  fclose(lcpcomp_fp);

  return;
}


void
Sarray_write_array (char *sarrayfile, Genome_T genomecomp, UINT4 genomelength) {
  UINT4 *SA;
  UINT4 n = genomelength;
  unsigned char *gbuffer;
  FILE *fp;


  SA = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_int_string(genomecomp,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */
  SACA_K(gbuffer,SA,n+/*virtual sentinel*/1,/*K, alphabet_size*/5,/*m*/n+1,/*level*/0);

  if ((fp = FOPEN_WRITE_BINARY(sarrayfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",sarrayfile);
    exit(9);
  } else {
    FWRITE_UINTS(SA,n+1,fp);
    fclose(fp);
  }

  FREE(gbuffer);
  FREE(SA);

  return;
}


void
Sarray_write_lcp_old (char *lcpptrsfile, char *lcpcompfile, char *sarrayfile,
		      Genome_T genomecomp, UINT4 genomelength) {
  UINT4 *SA, *lcp;
  UINT4 n = genomelength;
  int fd;
  size_t len;
#ifdef DEBUG14
  unsigned char *gbuffer;
#endif

  SA = (UINT4 *) Access_mmap(&fd,&len,sarrayfile,sizeof(UINT4),/*randomp*/true);

#ifdef DEBUG14
  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_int_string(genomecomp,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */
#endif

  fprintf(stderr,"Computing lcp...");
  lcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
#ifdef DEBUG14
  compute_lcp(lcp,gbuffer,SA,n);
#else
  compute_lcp(lcp,SA,n);
#endif
  fprintf(stderr,"done\n");
#ifdef DEBUG14
  FREE(gbuffer);
#endif

#ifdef WRITE_LCP
  if ((fp = FOPEN_WRITE_BINARY(lcpfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpfile);
    exit(9);
  } else {
    FWRITE_UINTS(lcp,n+1,fp);
    fclose(fp);
  }
#endif

  Sarray_write_lcpptrs(lcpptrsfile,lcpcompfile,lcp,n);

  FREE(lcp);
  munmap((void *) SA,len);
  close(fd);
  
  return;
}


/* Uses permuted lcp for speed and reduced memory usage */
void
Sarray_write_lcp (char *lcpptrsfile, char *lcpcompfile, char *sarrayfile,
		  Genome_T genomecomp, UINT4 genomelength) {
  UINT4 *plcp;
  UINT4 n = genomelength;
#ifdef READ_SA_FROM_FILE
  FILE *fp;
#else
  UINT4 *SA;
  int fd;
  size_t len;
#endif

  plcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));

#ifdef READ_SA_FROM_FILE
  if ((fp = FOPEN_READ_BINARY(sarrayfile)) == NULL) {
    fprintf(stderr,"Can't read from file %s\n",sarrayfile);
    exit(9);
  }
  compute_plcp(plcp,fp,n);
#else
  SA = (UINT4 *) Access_mmap(&fd,&len,sarrayfile,sizeof(UINT4),/*randomp*/false);
  compute_plcp(plcp,SA,n);
#endif


  fprintf(stderr,"Writing lcp file...");
#ifdef READ_SA_FROM_FILE
  fclose(fp);
  fp = FOPEN_READ_BINARY(sarrayfile);
  Sarray_write_lcpptrs_using_plcp(lcpptrsfile,lcpcompfile,plcp,fp,n);
#else
  Sarray_write_lcpptrs_using_plcp(lcpptrsfile,lcpcompfile,plcp,SA,n);
#endif
  fprintf(stderr,"done\n");

  FREE(plcp);

#ifdef READ_SA_FROM_FILE
  fclose(fp);
#else
  munmap((void *) SA,len);
  close(fd);
#endif
  
  return;
}


void
Sarray_write_index (char *saindexfile, char *sarrayfile, char *lcpptrsfile, char *lcpcompfile,
		    Genome_T genomecomp, UINT4 genomelength) {
  UINT4 n = genomelength, oligospace;
  UINT4 *saindex, *SA, *lcpptrs, *lcpcomp;
  FILE *fp;
  int sa_fd, lcpcomp_fd;
  size_t sa_len, lcpptrs_len, lcpcomp_len;
  double seconds;


  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);

  lcpptrs = (UINT4 *) Access_allocated(&lcpptrs_len,&seconds,lcpptrsfile,sizeof(UINT4));
  lcpcomp = (UINT4 *) Access_mmap(&lcpcomp_fd,&lcpcomp_len,lcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  lcpcomp = (UINT4 *) Access_mmap(&lcpcomp_fd,&lcpcomp_len,lcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  Bitpack64_access_setup(lcpptrs,lcpcomp);

  fprintf(stderr,"Computing saindex...");
  saindex = make_index(&oligospace,/*querylength*/12,genomecomp,SA,n);
  fprintf(stderr,"done\n");

  if ((fp = FOPEN_WRITE_BINARY(saindexfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",saindexfile);
    exit(9);
  } else {
    FWRITE_UINTS(saindex,oligospace+1,fp);
    fclose(fp);
  }

  FREE(saindex);
  FREE(lcpptrs);
  munmap((void *) lcpcomp,lcpcomp_len);
  close(lcpcomp_fd);
  munmap((void *) SA,sa_len);
  close(sa_fd);

  return;
}


void
Sarray_lcp_uncompress (char *lcpptrsfile, char *lcpcompfile, char *sarrayfile, UINT4 genomelength) {
  UINT4 n = genomelength, pos;
  UINT4 *SA, *lcpptrs, *lcpcomp;
  int sa_fd, lcpcomp_fd;
  size_t sa_len, lcpptrs_len, lcpcomp_len;
  double seconds;

  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);
  lcpptrs = (UINT4 *) Access_allocated(&lcpptrs_len,&seconds,lcpptrsfile,sizeof(UINT4));
  lcpcomp = (UINT4 *) Access_mmap(&lcpcomp_fd,&lcpcomp_len,lcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  lcpcomp = (UINT4 *) Access_mmap(&lcpcomp_fd,&lcpcomp_len,lcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  Bitpack64_access_setup(lcpptrs,lcpcomp);

  printf("i\tSA\tLCP\n");
  for (pos = 0; pos <= n; pos++) {
    printf("%u\t%u\t%u\n",pos,SA[pos],Bitpack64_access(pos));
  }

  FREE(lcpptrs);
  munmap((void *) lcpcomp,lcpcomp_len);
  close(lcpcomp_fd);
  munmap((void *) SA,sa_len);
  close(sa_fd);

  return;
}
