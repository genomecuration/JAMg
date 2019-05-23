static char rcsid[] = "$Id: extension-search.c 218736 2019-03-25 16:53:44Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "extension-search.h"

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
#include "univcoord.h"
#include "oligo.h"

#include "list.h"
#include "maxent_hr.h"

#include "genome128_hr.h"
#include "genome128_consec.h"

#include "univdiagdef.h"
#include "univdiag.h"
#include "path-solve.h"

#if defined(HAVE_SSE2)
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#ifdef HAVE_AVX512
#include <immintrin.h>
#endif


/* conservative behaves better for high-quality sequences and is faster */
/* #define LIBERAL 1 */
#define LIBERAL_KMER 8

/* Some limit is needed to prevent GSNAP from running very slowly */
#define MAX_HITS_FOR_BEST_ELT 1000


/* #define CHECK_OLIGOS 1 */

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Sorting of diagonals */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* filter_elts */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif


static Mode_T mode;

static Genome_T genomebits;
static Genome_T genomebits_alt;

static Univ_IIT_T chromosome_iit;
static int circular_typeint;
static bool *circularp;

static Genome_T genomebits;
static Genome_T genomebits_alt;
static Indexdb_T indexdb;
static Indexdb_T indexdb2;


static char conversion_fwd[128];
static char conversion_rev[128];

static Chrpos_T overall_max_distance_genome;
static Chrpos_T overall_end_distance_genome;
static Chrpos_T shortsplicedist;
static Chrpos_T min_intronlength;
static Chrpos_T max_deletionlen;
static Chrpos_T max_insertionlen_default;
static int max_end_deletions;
static int max_middle_insertions_default;

static int index1part;
static int local1part;
static int leftreadshift;
static Oligospace_T oligobase_mask;


#ifdef LARGE_GENOMES
#define GETPOS(high,low) (((Univcoord_T) high << 32) + low)
#endif


void
Extension_search_setup (Mode_T mode_in, Univ_IIT_T chromosome_iit_in, int circular_typeint_in, bool *circularp_in,
			Genome_T genomebits_in, Genome_T genomebits_alt_in, Indexdb_T indexdb_in, Indexdb_T indexdb2_in,
			Chrpos_T shortsplicedist_in, Chrpos_T shortsplicedist_novelend,
			int min_intronlength_in, int max_deletionlength, int max_end_deletions_in,
			int max_middle_insertions_in, int max_end_insertions,
			int index1part_in, int local1part_in) {
  int i;

  mode = mode_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  indexdb = indexdb_in;
  indexdb2 = indexdb2_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  circularp = circularp_in;

  for (i = 0; i < 128; i++) {
    conversion_fwd[i] = i;
    conversion_rev[i] = i;
  }
  if (mode == STANDARD) {
    /* Don't change conversion */
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    conversion_fwd['C'] = 'T';	/* CT */
    conversion_rev['G'] = 'A';	/* GA */
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    conversion_fwd['A'] = 'G';	/* AG */
    conversion_rev['T'] = 'C';	/* TC */
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    conversion_fwd['T'] = 'C';	/* TC */
    conversion_rev['A'] = 'G';	/* AG */
  }

  chromosome_iit = chromosome_iit_in;
  circular_typeint = circular_typeint_in;
  shortsplicedist = shortsplicedist_in;

  min_intronlength = min_intronlength_in;
  max_deletionlen = max_deletionlength;
  max_end_deletions = max_end_deletions_in;
  max_middle_insertions_default = max_middle_insertions_in;
  if (max_middle_insertions_in > max_end_insertions) {
    max_insertionlen_default = max_middle_insertions_in;
  } else {
    max_insertionlen_default = max_end_insertions;
  }

  if (shortsplicedist > max_deletionlen) {
    overall_max_distance_genome = shortsplicedist;
  } else {
    overall_max_distance_genome = max_deletionlen;
  }

  if (shortsplicedist_novelend > max_deletionlen) {
    overall_end_distance_genome = shortsplicedist_novelend;
  } else {
    overall_end_distance_genome = max_deletionlen;
  }

  index1part = index1part_in;
  local1part = local1part_in;

#ifdef HAVE_64_BIT
  leftreadshift = 64 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#else
  leftreadshift = 32 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#endif

  return;
}


/* Simplified version of Spanningelt_T */
#define T Elt_T


#if 0
/* Works only for queryfwd.  Use Elt_read functions instead */
static T
Elt_new (int querypos, int nmatches, Univcoord_T *all_diagonals, int n_all_diagonals) {
  T new = (T) MALLOC(sizeof(*new));
#ifdef DEBUG
  int i;
#endif

  new->qstart = querypos;
  new->qend = querypos + nmatches - 1;
  new->nmatches = nmatches;
  debug(printf("Making an Elt with querystart %d, nmatches %d => queryend %d\n",
	       new->qstart,new->nmatches,new->qend));

  new->all_diagonals = all_diagonals;
  new->n_all_diagonals = n_all_diagonals;

  new->diagonals = &(new->all_diagonals[0]);
  new->ndiagonals = n_all_diagonals;

#ifdef DEBUG
  printf("Diagonals:");
  for (i = 0; i < n_all_diagonals; i++) {
    printf(" %llu",all_diagonals[i]);
  }
  printf("\n");
#endif

  new->lowi = 0;
  new->highi = n_all_diagonals;

  return new;
}
#endif


static void
Elt_free (T *old) {

  FREE((*old)->all_diagonals);
  FREE(*old);
  return;
}

void
Elt_gc (List_T *set) {
  List_T p;
  T elt;

  for (p = *set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    Elt_free(&elt);
  }
  /* List_free(&(*set)); -- allocated by Listpool_push */
  return;
}



#if 0
static int
nt_querylength (char *query, int querylength) {
  int i;
  char c;

  i = 0;
  while (i < querylength && ((c = query[i]) == 'A' || c == 'C' || c == 'G' || c == 'T')) {
    i++;
  }

  return i;
}
#endif

#ifdef CHECK_OLIGOS
static Oligospace_T
nt_oligo (char *query, int indexsize) {
  Oligospace_T oligo = 0U;
  int i;

  for (i = 0; i < indexsize; i++) {
    oligo *= 4;
    
    printf("%c",query[i]);
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
  printf("\n");

  return oligo;
}

#define LOW_TWO_BITS 0x3


static char *
oligo_nt (UINT4 oligo, int oligosize) {
  char *nt = MALLOC((oligosize+1)*sizeof(char));
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
  nt[oligosize] = '\0';

  return nt;
}

#endif


#if !defined(HAVE_SSE4_2)
#define count_leading_zeroes_32(diff) ((diff >> 16) ? clz_table[diff >> 16] : 16 + clz_table[diff])
#elif defined(HAVE_LZCNT)
#define count_leading_zeroes_32(diff) _lzcnt_u32(diff)
#elif defined(HAVE_BUILTIN_CLZ)
#define count_leading_zeroes_32(diff) __builtin_clz(diff)
#else
#define count_leading_zeroes_32(diff) ((diff >> 16) ? clz_table[diff >> 16] : 16 + clz_table[diff])
#endif

#if !defined(HAVE_SSE4_2)
#define count_trailing_zeroes_32(diff) mod_37_bit_position[(-diff & diff) % 37]
#elif defined(HAVE_TZCNT)
#define count_trailing_zeroes_32(diff) _tzcnt_u32(diff)
#elif defined(HAVE_BUILTIN_CTZ)
#define count_trailing_zeroes_32(diff) __builtin_ctz(diff)
#else
/* lowbit = -diff & diff */
#define count_trailing_zeroes_32(diff) mod_37_bit_position[(-diff & diff) % 37]
#endif


/* query is a substring of the original, starting with queryoffset */
/* Performs the same function as Sarray_lookup, which returns the diagonals */
static T
Elt_read_queryfwd (
#ifdef LARGE_GENOMES
		   unsigned char *positions_high,
#endif
		   UINT4 *positions, int n, int diagterm, int querylength, int querystart,
		   Compress_T query_compress, bool plusp, int genestrand) {
  T new;
  int max_nmatches;
  Univcoord_T *best_diagonals, *out, diagonal;
  int *nmatches;
  int i;
#ifdef LIBERAL
  int queryend;
  int nmismatches3;
#endif
  

  debug(printf("Got %d positions\n",n));
  if (n == 0) {
    return (T) NULL;
  } else if (querystart >= querylength - index1part) {
    return (T) NULL;
  } else {
    max_nmatches = 0;
    nmatches = (int *) MALLOC(n*sizeof(int));
    for (i = 0; i < n; i++) {
#ifdef LARGE_GENOMES
      assert(GETPOS(positions_high[i],positions[i]) >= (Univcoord_T) -diagterm);
      diagonal = GETPOS(positions_high[i],positions[i]) + diagterm;
      debug(printf("plusp %d, diagonal is %llu, querystart is %d\n",plusp,diagonal,querystart));
#else
      assert(positions[i] >= (Univcoord_T) -diagterm);
      diagonal = positions[i] + diagterm;
      debug(printf("plusp %d, diagonal is %u, querystart is %d\n",plusp,diagonal,querystart));
#endif

#ifdef LIBERAL	
      /* Uses Genome_first_kmer_right, which is liberal */
      queryend = Genome_first_kmer_right(&nmismatches3,genomebits,query_compress,/*left*/diagonal,
					 /*pos5*/querystart,/*pos3*/querylength,
					 plusp,genestrand,/*query_unk_mismatch_local_p*/false,
					 /*kmer*/LIBERAL_KMER);
      if ((nmatches[i] = queryend - querystart) > max_nmatches) {
	max_nmatches = nmatches[i];
      }
#else
      /* Uses Genome_consecutive_matches_rightward, which is conservative */
      nmatches[i] = index1part +
	Genome_consecutive_matches_rightward(genomebits,query_compress,/*left*/diagonal,
					     /*pos5*/querystart+index1part,/*pos3*/querylength,
					     plusp,genestrand);
      if (nmatches[i] > max_nmatches) {
	max_nmatches = nmatches[i];
      }
#endif
    }

    out = best_diagonals = MALLOC(n*sizeof(Univcoord_T));
    for (i = 0; i < n; i++) {
      if (nmatches[i] >= max_nmatches - /*index1interval*/ 3 + 1) {
#ifdef LARGE_GENOMES
	*out++ = GETPOS(positions_high[i],positions[i]) + diagterm;
#else
	*out++ = positions[i] + diagterm;
#endif
      }
    }
    FREE(nmatches);
   

    new = (T) MALLOC(sizeof(*new));
    new->qstart = querystart;
    new->qend = querystart + max_nmatches;
    new->nmatches = max_nmatches;

    assert(new->qend <= (int) querylength);
    debug(printf("Making an Elt with querystart %d, nmatches %d => queryend %d\n",
		 new->qstart,new->nmatches,new->qend));

    new->all_diagonals = new->diagonals = best_diagonals;
    new->n_all_diagonals = new->ndiagonals = out - best_diagonals;

#ifdef DEBUG
    printf("Diagonals:");
    for (i = 0; i < new->ndiagonals; i++) {
      printf(" %llu",new->diagonals[i]);
    }
    printf("\n");
#endif

    new->lowi = 0;
    new->highi = new->n_all_diagonals;

    return new;
  }
}


static T
Elt_read_queryrev (
#ifdef LARGE_GENOMES
		   unsigned char *positions_high,
#endif
		   UINT4 *positions, int n, int diagterm, int queryend,
		   Compress_T query_compress, bool plusp, int genestrand) {
  T new;
  int max_nmatches;
  Univcoord_T *best_diagonals, *out, diagonal;
  int *nmatches;
  int i;
#ifdef LIBERAL
  int querystart;
  int nmismatches5;
#endif
  

  debug(printf("Got %d positions\n",n));
  if (n == 0) {
    return (T) NULL;
  } else if (queryend < index1part) {
    return (T) NULL;
  } else {
    max_nmatches = 0;
    nmatches = (int *) MALLOC(n*sizeof(int));
    for (i = 0; i < n; i++) {
#ifdef LARGE_GENOMES
      assert(GETPOS(positions_high[i],positions[i]) >= (Univcoord_T) -diagterm);
      diagonal = GETPOS(positions_high[i],positions[i]) + diagterm;
      debug(printf("plusp %d, diagonal is %llu, queryend is %d\n",plusp,diagonal,queryend));
#else
      assert(positions[i] >= (Univcoord_T) -diagterm);
      diagonal = positions[i] + diagterm;
      debug(printf("plusp %d, diagonal is %u, queryend is %d\n",plusp,diagonal,queryend));
#endif

#ifdef LIBERAL	
      /* Uses Genome_first_kmer_left, which is liberal */
      querystart = Genome_first_kmer_left(&nmismatches5,genomebits,query_compress,/*left*/diagonal,
					  /*pos5*/0,/*pos3*/queryend,
					  plusp,genestrand,/*query_unk_mismatch_local_p*/false,
					  /*kmer*/LIBERAL_KMER);
      if ((nmatches[i] = queryend - querystart) > max_nmatches) {
	max_nmatches = nmatches[i];
      }
#else
      /* Uses Genome_consecutive_matches_rightward, which is conservative */
      nmatches[i] = index1part +
	Genome_consecutive_matches_leftward(genomebits,query_compress,/*left*/diagonal,
					    /*pos5*/0,/*pos3*/queryend - index1part,plusp,genestrand);
      if (nmatches[i] > max_nmatches) {
	max_nmatches = nmatches[i];
      }
#endif
    }

    out = best_diagonals = MALLOC(n*sizeof(Univcoord_T));
    for (i = 0; i < n; i++) {
      if (nmatches[i] >= max_nmatches - /*index1interval*/ 3 + 1) {
#ifdef LARGE_GENOMES
	*out++ = GETPOS(positions_high[i],positions[i]) + diagterm;
#else
	*out++ = positions[i] + diagterm;
#endif
      }
    }
    FREE(nmatches);
   

    new = (T) MALLOC(sizeof(*new));
    new->qend = queryend;
    new->qstart = queryend - max_nmatches;
    new->nmatches = max_nmatches;

    assert(new->qstart >= 0);
    debug(printf("Making an Elt with queryend %d, nmatches %d => querystart %d\n",
		 new->qend,new->nmatches,new->qstart));

    new->all_diagonals = new->diagonals = best_diagonals;
    new->n_all_diagonals = new->ndiagonals = out - best_diagonals;

#ifdef DEBUG
    printf("Diagonals:");
    for (i = 0; i < new->ndiagonals; i++) {
      printf(" %llu",new->diagonals[i]);
    }
    printf("\n");
#endif

    new->lowi = 0;
    new->highi = new->n_all_diagonals;

    return new;
  }
}


void
Elt_set_queryfwd_extend_left (List_T list, Compress_T query_compress, bool plusp, int genestrand) {
  List_T p;
  T this;
  int max_nmatches;
  Univcoord_T *best_diagonals, *out, diagonal;
  int *nmatches;
  int n, i;
#ifdef LIBERAL
  int querystart;
  int nmismatches5;
#endif
  

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    n = this->n_all_diagonals;
    debug(printf("Starting with %d positions\n",n));
    assert(n > 0);

    max_nmatches = 0;
    nmatches = (int *) MALLOC(n*sizeof(int));
    for (i = 0; i < n; i++) {
      diagonal = this->all_diagonals[i];
      
#ifdef LIBERAL	
      /* Uses Genome_first_kmer_left, which is liberal */
      querystart = Genome_first_kmer_left(&nmismatches5,genomebits,query_compress,/*left*/diagonal,
					  /*pos5*/0,/*pos3*/this->qend,
					  plusp,genestrand,/*query_unk_mismatch_local_p*/false,
					  /*kmer*/LIBERAL_KMER);
      if ((nmatches[i] = this->queryend - querystart) > max_nmatches) {
	max_nmatches = nmatches[i];
      }
#else
      /* Uses Genome_consecutive_matches_rightward, which is conservative */
      nmatches[i] = index1part +
	Genome_consecutive_matches_leftward(genomebits,query_compress,/*left*/diagonal,
					    /*pos5*/0,/*pos3*/this->qend - index1part,plusp,genestrand);
      if (nmatches[i] > max_nmatches) {
	max_nmatches = nmatches[i];
      }
#endif
    }
    
    out = best_diagonals = MALLOC(n*sizeof(Univcoord_T));
    for (i = 0; i < n; i++) {
      if (nmatches[i] >= max_nmatches - /*index1interval*/ 3 + 1) {
	*out++ = this->all_diagonals[i];
      }
    }
    FREE(nmatches);
    
    this->qstart = this->qend - max_nmatches;
    this->nmatches = max_nmatches;
    
    debug(printf("Revising an Elt with queryend %d, nmatches %d => querystart %d\n",
		 this->qend,this->nmatches,this->qstart));
    FREE(this->all_diagonals);
    this->all_diagonals = this->diagonals = best_diagonals;
    this->n_all_diagonals = this->ndiagonals = out - best_diagonals;
    
#ifdef DEBUG
    printf("Diagonals:");
    for (i = 0; i < this->ndiagonals; i++) {
      printf(" %llu",this->diagonals[i]);
    }
    printf("\n");
#endif
    
    this->lowi = 0;
    this->highi = this->n_all_diagonals;
  }

  return;
}


void
Elt_set_queryrev_extend_right (List_T list, int querylength,
			       Compress_T query_compress, bool plusp, int genestrand) {
  List_T p;
  T this;
  int max_nmatches;
  Univcoord_T *best_diagonals, *out, diagonal;
  int *nmatches;
  int n, i;
#ifdef LIBERAL
  int queryend;
  int nmismatches3;
#endif
  

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    n = this->n_all_diagonals;
    debug(printf("Starting with %d positions\n",n));
    assert(n > 0);

    max_nmatches = 0;
    nmatches = (int *) MALLOC(n*sizeof(int));
    for (i = 0; i < n; i++) {
      diagonal = this->all_diagonals[i];
      
#ifdef LIBERAL	
      /* Uses Genome_first_kmer_right, which is liberal */
      queryend = Genome_first_kmer_right(&nmismatches3,genomebits,query_compress,/*left*/diagonal,
					 /*pos5*/this->querystart,/*pos3*/querylength,
					 plusp,genestrand,/*query_unk_mismatch_local_p*/false,
					 /*kmer*/LIBERAL_KMER);
      if ((nmatches[i] = queryend - this->querystart) > max_nmatches) {
	max_nmatches = nmatches[i];
      }
#else
      /* Uses Genome_consecutive_matches_rightward, which is conservative */
      nmatches[i] = index1part +
	Genome_consecutive_matches_rightward(genomebits,query_compress,/*left*/diagonal,
					     /*pos5*/this->qstart + index1part,/*pos3*/querylength,
					     plusp,genestrand);
      if (nmatches[i] > max_nmatches) {
	max_nmatches = nmatches[i];
      }
#endif
    }
    
    out = best_diagonals = MALLOC(n*sizeof(Univcoord_T));
    for (i = 0; i < n; i++) {
      if (nmatches[i] >= max_nmatches - /*index1interval*/ 3 + 1) {
	*out++ = this->all_diagonals[i];
      }
    }
    FREE(nmatches);
    
    
    this->qend = this->qstart + max_nmatches;
    this->nmatches = max_nmatches;
    
    debug(printf("Revising an Elt with querystart %d, nmatches %d => queryend %d\n",
		 this->qstart,this->nmatches,this->qend));
    FREE(this->all_diagonals);
    this->all_diagonals = this->diagonals = best_diagonals;
    this->n_all_diagonals = this->ndiagonals = out - best_diagonals;
    
#ifdef DEBUG
    printf("Diagonals:");
    for (i = 0; i < this->ndiagonals; i++) {
      printf(" %llu",this->diagonals[i]);
    }
    printf("\n");
#endif
    
    this->lowi = 0;
    this->highi = this->n_all_diagonals;
  }

  return;
}


static int
binary_search_univcoord (int lowi, int highi, Univcoord_T *diagonals, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,diagonals[lowi],middlei,diagonals[middlei],
		   highi-1,diagonals[highi-1],goal));
    if (goal < diagonals[middlei]) {
      highi = middlei;
    } else if (goal > diagonals[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


static void
Elt_filter_diagonals (T this, Univcoord_T low, Univcoord_T high) {
  int lowi, highi;

  debug(printf("Entered Elt_filter_diagonals on %d..%d with %d diagonals with low %u and high %u, nmatches %d\n",
	       this->qstart,this->qend,this->n_all_diagonals,low,high,this->nmatches));

  /* low_adj and high_adj are inclusive */
  lowi = binary_search_univcoord(/*lowi*/0,/*highi*/this->n_all_diagonals,this->all_diagonals,/*goal*/low);
  highi = binary_search_univcoord(lowi,/*highi*/this->n_all_diagonals,this->all_diagonals,/*goal*/high + 1) - 1;
  if ((this->ndiagonals = highi - lowi + 1) == 0) {
    this->diagonals = (Univcoord_T *) NULL;

  } else {
    this->diagonals = &(this->all_diagonals[lowi]);
  }

  debug(printf("Setting lowi %d and highi %d\n",lowi,highi));

  return;
}


#ifdef DEBUG
static void
Elt_dump (T elt) {
  int k;

  printf("Elt %d..%d with %d diagonals:\n",elt->qstart,elt->qend,elt->ndiagonals);
  for (k = 0; k < elt->ndiagonals; k++) {
    printf("  %u\n",elt->diagonals[k]);
  }
  printf("\n");

  return;
}

static void
Elt_dump_set (List_T set) {
  List_T p;

  for (p = set; p != NULL; p = List_next(p)) {
    Elt_dump((T) List_head(p));
  }

  return;
}
#endif


#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


static bool
elt_startp (T elt, int middle_qstart, int middle_qend) {
  if (elt->qstart >= middle_qstart && elt->qend <= middle_qend) {
    debug13(printf("Not allowing left elt that is subsumed by middle elt: q %d..%d\n",
		   elt->qstart,elt->qend));
    return false;
  } else if (elt->qend >= middle_qend) {
    debug13(printf("Not allowing left elt that extends right of middle elt: q %d..%d\n",
		   elt->qstart,elt->qend));
    return false;
  } else if ((elt->qend - middle_qstart) > (middle_qend - middle_qstart) / 2) {
    debug13(printf("Not allowing left elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->qstart,elt->qend));
    return false;
  } else if ((elt->qend - middle_qstart) > (elt->qend - elt->qstart) / 2) {
    debug13(printf("Not allowing left elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->qstart,elt->qend));
    return false;
  } else {
    return true;
  }
}


static bool
elt_endp (T elt, int middle_qstart, int middle_qend) {
  if (elt->qstart >= middle_qstart && elt->qend <= middle_qend) {
    debug13(printf("Not allowing right elt that is subsumed by middle elt: qpos %d..%d\n",
		   elt->qstart,elt->qend));
    return false;
  } else if (elt->qstart <= middle_qstart) {
    debug13(printf("Not allowing right elt that extends left of middle elt: qpos %d..%d\n",
		   elt->qstart,elt->qend));
    return false;
  } else if ((middle_qend - elt->qstart) > (middle_qend - middle_qstart) / 2) {
    debug13(printf("Not allowing right elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->qstart,elt->qend));
    return false;
  } else if ((middle_qend - elt->qstart) > (elt->qend - elt->qstart) / 2) {
    debug13(printf("Not allowing right elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->qstart,elt->qend));
    return false;
  } else {
    return true;
  }
}


#if 0
static void
filter_elts (T *elt_array, int best_i, int nelts) {
  T elt, middle;
  int i;

  middle = elt_array[best_i];

  debug13(printf("Best elt is %d..%d with %d diagonals\n",middle->qstart,middle->qend,middle->ndiagonals));

  /* Filter elts and diagonals for right side */
  debug13(printf("Filtering elts on right side\n"));
  for (i = best_i + 1; i < nelts; i++) {
    elt = elt_array[i];
    debug13(printf("Handling right elt %d: %p\n",i,elt));
    if (elt->qstart >= middle->qstart && elt->qend <= middle->qend) {
      debug13(printf("Not allowing right elt that is subsumed by middle elt: qpos %d..%d\n",
		     elt->qstart,elt->qend));
      /* Elt_free(&elt); */
      elt_array[i] = (Elt_T) NULL;
    } else if (elt->qstart <= middle->qstart) {
      debug13(printf("Not allowing right elt that extends left of middle elt: qpos %d..%d\n",
		     elt->qstart,elt->qend));
      /* Elt_free(&elt); */
      elt_array[i] = (Elt_T) NULL;
    } else if ((middle->qend - elt->qstart) > (middle->qend - middle->qstart) / 2) {
      debug13(printf("Not allowing right elt that mainly overlaps middle elt: q %d..%d\n",
		     elt->qstart,elt->qend));
      /* Elt_free(&elt); */
      elt_array[i] = (Elt_T) NULL;
    } else if ((middle->qend - elt->qstart) > (elt->qend - elt->qstart) / 2) {
      debug13(printf("Not allowing right elt that mainly overlaps middle elt: q %d..%d\n",
		     elt->qstart,elt->qend));
      /* Elt_free(&elt); */
      elt_array[i] = (Elt_T) NULL;
    }
  }

  /* Filter elts and diagonals for left side */
  debug13(printf("Filtering elts on left side\n"));
  for (i = best_i - 1; i >= 0; --i) {
    elt = elt_array[i];
    debug13(printf("Handling left elt %d: %p\n",i,elt));
    if (elt->qstart >= middle->qstart && elt->qend <= middle->qend) {
      debug13(printf("Not allowing left elt that is subsumed by middle elt: q %d..%d\n",
		     elt->qstart,elt->qend));
      /* Elt_free(&elt); */
      elt_array[i] = (Elt_T) NULL;
    } else if (elt->qend >= middle->qend) {
      debug13(printf("Not allowing left elt that extends right of middle elt: q %d..%d\n",
		     elt->qstart,elt->qend));
      /* Elt_free(&elt); */
      elt_array[i] = (Elt_T) NULL;
    } else if ((elt->qend - middle->qstart) > (middle->qend - middle->qstart) / 2) {
      debug13(printf("Not allowing left elt that mainly overlaps middle elt: q %d..%d\n",
		     elt->qstart,elt->qend));
      /* Elt_free(&elt); */
      elt_array[i] = (Elt_T) NULL;
    } else if ((elt->qend - middle->qstart) > (elt->qend - elt->qstart) / 2) {
      debug13(printf("Not allowing left elt that mainly overlaps middle elt: q %d..%d\n",
		     elt->qstart,elt->qend));
      /* Elt_free(&elt); */
      elt_array[i] = (Elt_T) NULL;
    }
  }

  return;
}
#endif


#if 0
/* Currently returning a score, but this is not used by the caller */
static int
get_diagonals (Univdiag_T *middle_diagonal, List_T *best_right_diagonals, List_T *best_left_diagonals, 
	       List_T *all_right_diagonals, List_T *all_left_diagonals, int querylength,
	       Univcoord_T chroffset, Univcoord_T chrhigh, Univcoord_T goal,
	       Univdiagpool_T univdiagpool, T *elt_array, int best_i, int nelts) {
  int best_score_right, best_score_left, best_score, score;
  T elt;
  Univcoord_T low, high, left_position, right_position;

  int min_qstart, max_qend; /* Checked using localdb */
  Univcoord_T *diagonals;
  int ndiagonals;

  int i, j;
  List_T p;

  Univdiag_T *diagonal_array, diagonal, prev_diagonal;
  int qpos;
  List_T left_diagonals, right_diagonals;

  int max_insertionlen;

  if (max_middle_insertions_default >= 0) {
    max_insertionlen = max_insertionlen_default;
  } else {
    max_insertionlen = querylength;
  }

  debug13(printf("\n***Entered get_diagonals with goal %u, chroffset %u, chrhigh %u\n",goal,chroffset,chrhigh));

  /* Create diagonals.  We give a bonus of +1 for being on the same
     diagonal.  This means that we should count consecutive regions
     within each diagonal as 2 points.  Then an indel or gap will
     give only 1 point, or a relative penalty. */

  elt = elt_array[best_i];
  left_position = goal + elt->qstart;
  right_position = goal + elt->qend;
  min_qstart = elt->qstart;
  max_qend = elt->qend;

  *middle_diagonal = Univdiag_new(univdiagpool,elt->qstart,elt->qend,/*univdiagonal*/goal);
  (*middle_diagonal)->intscore = 2*(elt->qend - elt->qstart + 1);
  debug13(printf("Creating middle diagonal: query %d..%d, diagonal %u = goal %u - chroffset %u\n",
		 elt->qstart,elt->qend,goal - chroffset,goal,chroffset));


  /* Filter diagonals for right side */
  right_diagonals = (List_T) NULL;

  debug13(printf("Filtering elts and diagonals for right side.  Subtracting %d and adding %d from goal %u\n",
		 max_insertionlen,overall_end_distance_genome,goal));
  low = subtract_bounded(goal,/*minusterm*/max_insertionlen,chroffset);
  high = add_bounded(goal,/*plusterm*/overall_end_distance_genome,chrhigh);

  for (i = best_i + 1; i < nelts; i++) {
    if ((elt = elt_array[i]) != NULL) {
      debug13(printf("Handling right elt %d: %p\n",i,elt));
      Elt_filter_diagonals(elt,low,high);
      if (elt->ndiagonals > 0) {
	diagonals = elt->diagonals;
	ndiagonals = elt->ndiagonals;
      } else {
	debug13(printf("Elt %d..%d has no diagonals\n",elt->qstart,elt->qend));
	/* ? TODO: Restrict so we search only short elts */
	ndiagonals = 0;
      }

      /* successp = false; */
      for (j = 0; j < ndiagonals; j++) {
	debug13(printf("position %d is %u\n",j,diagonals[j]));
	if (diagonals[j] + elt->qend <= right_position) {
	  debug13(printf("Not creating right diagonal that does not advance past right_position %u: query %d..%d, diagonal %u\n",
			 right_position - chroffset,elt->qstart,elt->qend,diagonals[j] - chroffset));
	} else {
	  debug13(printf("Creating right diagonal: query %d..%d, diagonal %u\n",
			 elt->qstart,elt->qend,diagonals[j] - chroffset));
	  right_diagonals = Univdiagpool_push(right_diagonals,univdiagpool,elt->qstart,elt->qend,
					      /*univdiagonal*/diagonals[j]);
	  /* In resetting right_position, we are assuming that this diagonal will be used */
	  /* right_position = diagonals[j] + elt->qend; */
	  if (elt->qend > max_qend) {
	    /* In case we succeeded with the original elt diagonals */
	    max_qend = elt->qend;
	  }
	  /* successp = true; */
	}
      }

#if 0
      /* In resetting low and high, we are assuming that this diagonal will be used */
      if (successp == true) {
	/* Success: Update low and high for next search */
	low = subtract_bounded(diagonals[0],/*minusterm*/max_insertionlen,chroffset);
	high = add_bounded(diagonals[ndiagonals-1],/*plusterm*/overall_end_distance_genome,chrhigh);
	debug13(printf("Updated low and high based on %u and %u\n",diagonals[0],diagonals[ndiagonals-1]));
      }
#endif

    }
  }

  debug13(printf("max_qend = %d\n",max_qend));
  right_diagonals = List_reverse(right_diagonals);



  /* Filter diagonals for left side */
  left_diagonals = (List_T) NULL;

  debug13(printf("Filtering elts and diagonals for left side.  Subtracting %d and adding %d from goal %u\n",
		 overall_end_distance_genome,max_insertionlen,goal));
  low = subtract_bounded(goal,/*minusterm*/overall_end_distance_genome,chroffset);
  high = add_bounded(goal,/*plusterm*/max_insertionlen,chrhigh);

  for (i = best_i - 1; i >= 0; --i) {
    if ((elt = elt_array[i]) != NULL) {
      debug13(printf("Handling left elt %d: %p\n",i,elt));
      Elt_filter_diagonals(elt,low,high);
      if (elt->ndiagonals > 0) {
	diagonals = elt->diagonals;
	ndiagonals = elt->ndiagonals;
      } else {
	debug13(printf("Elt %d..%d has no diagonals\n",elt->qstart,elt->qend));
	/* ? TODO: Restrict so we search only short elts */
	ndiagonals = 0;
      }

      /* successp = false; */
      for (j = ndiagonals - 1; j >= 0; --j) {
	debug13(printf("position %d is %u\n",j,diagonals[j]));
	if (diagonals[j] + elt->qstart >= left_position) {
	  debug13(printf("Not creating left diagonal that does not advance past left_position %u: query %d..%d, diagonal %u\n",
			 left_position - chroffset,elt->qstart,elt->qend,diagonals[j] - chroffset));
	} else {
	  debug13(printf("Creating left diagonal: query %d..%d, diagonal %u\n",
			 elt->qstart,elt->qend,diagonals[j] - chroffset));
	  left_diagonals = Univdiagpool_push(left_diagonals,univdiagpool,elt->qstart,elt->qend,
					     /*univdiagonal*/diagonals[j]);
	  /* In resetting left_position, we are assuming that this diagonal will be used */
	  /* left_position = diagonals[j] + elt->qstart; */
	  if (elt->qstart < min_qstart) {
	    /* In case we succeeded with the original elt diagonals */
	    min_qstart = elt->qstart;
	  }
	  /* successp = true; */
	}
      }

#if 0
      /* In resetting right_position, we are assuming that this diagonal will be used */
      if (successp == true) {
	/* Success: Update low and high for next search */
	low = subtract_bounded(diagonals[0],/*minusterm*/overall_end_distance_genome,chroffset);
	high = add_bounded(diagonals[ndiagonals-1],/*plusterm*/max_insertionlen,chrhigh);
	debug13(printf("Updated low and high based on %u and %u\n",diagonals[0],diagonals[ndiagonals-1]));
      }
#endif

    }
  }

  debug13(printf("min_qstart = %d\n",min_qstart));
  left_diagonals = List_reverse(left_diagonals);



  /* A.  Compute right diagonals */
  /* A1.  Scoring for dynamic programming */
  diagonal_array = (Univdiag_T *) List_to_array_n(&ndiagonals,right_diagonals);
  /* List_free(&right_diagonals); -- allocated by Univdiagpool_push */
#ifdef DEBUG12
  printf("Right side before consolidating\n");
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    printf("%d..%d at %u\n",diagonal->qstart,diagonal->qend,diagonal->univdiagonal);
  }
#endif

  *all_right_diagonals = (List_T) NULL;
  qsort(diagonal_array,ndiagonals,sizeof(Univdiag_T),Univdiag_diagonal_cmp);
  i = 0;
  while (i < ndiagonals) {
    j = i;
    while (j < ndiagonals && diagonal_array[j]->univdiagonal == diagonal_array[i]->univdiagonal) {
      j++;
    }
    if (j == i) {
      *all_right_diagonals = Univdiagpool_push_existing(*all_right_diagonals,univdiagpool,diagonal_array[i]);
    } else {
      *all_right_diagonals = Univdiagpool_push(*all_right_diagonals,univdiagpool,diagonal_array[i]->qstart,
					       diagonal_array[j-1]->qend,diagonal_array[i]->univdiagonal);
#if 0
      for (k = i; k < j; k++) {
	/* Univdiag_free(&(diagonal_array[k])); -- allocated by Univdiagpool_push */
      }
#endif
    }
    i = j;
  }
  FREE(diagonal_array);

  /* TODO: May be able to skip this sorting step */
  diagonal_array = (Univdiag_T *) List_to_array_n(&ndiagonals,*all_right_diagonals);
  qsort(diagonal_array,ndiagonals,sizeof(Univdiag_T),Univdiag_ascending_cmp);
#ifdef DEBUG12
  printf("Right side after consolidating and sorting\n");
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    printf("%d..%d at %u\n",diagonal->qstart,diagonal->qend,diagonal->univdiagonal);
  }
#endif


  debug13(printf("%d diagonals:\n",ndiagonals));
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    debug13(printf("%d: %d..%d at %u\n",i,diagonal->qstart,diagonal->qend,diagonal->univdiagonal));

    low = subtract_bounded(diagonal->univdiagonal,overall_max_distance_genome,chroffset);
    high = add_bounded(diagonal->univdiagonal,max_insertionlen,chrhigh);
    qpos = diagonal->qstart;
    best_score = 0;

    for (j = i - 1; j >= 0; --j) {
      prev_diagonal = diagonal_array[j];
      debug13(printf("  %d: %d..%d at %u  ",j,prev_diagonal->qstart,prev_diagonal->qend,prev_diagonal->univdiagonal));

      if (prev_diagonal->qend >= qpos) {
	debug13(printf("Skipping because qend %d >= qpos %d\n",prev_diagonal->qend,qpos));
      } else if (prev_diagonal->univdiagonal < low) {
	debug13(printf("Skipping because diagonal %u < low_chrpos %u\n",prev_diagonal->univdiagonal,low));
      } else if (prev_diagonal->univdiagonal > high) {
	debug13(printf("Skipping because diagonal %u > high_chrpos %u\n",prev_diagonal->univdiagonal,high));
      } else {
	score = prev_diagonal->intscore;
	if (prev_diagonal->univdiagonal == diagonal->univdiagonal) {
	  score += 1;
	}
	if (score <= best_score) {
	  debug13(printf("Skipping because score %d <= best_score %d\n",score,best_score));
	} else {
	  best_score = score;
	  /* diagonal->prev = prev_diagonal; -- Not used */
	  debug13(printf("Updating best score to be %d.  Prev diagonal is %d..%d at %llu\n",
			 best_score,prev_diagonal->qstart,prev_diagonal->qend,prev_diagonal->univdiagonal));
	}
      }
    }

    /* Handle links to middle diagonal */
    prev_diagonal = *middle_diagonal;
    debug13(printf("  Middle: %d..%d at %u  ",prev_diagonal->qstart,prev_diagonal->qend,prev_diagonal->univdiagonal));
    if (prev_diagonal->qend >= qpos) {
      debug13(printf("Skipping because qend %d >= qpos %d\n",prev_diagonal->qend,qpos));
    } else if (prev_diagonal->univdiagonal < low) {
      debug13(printf("Skipping because diagonal %u < low_chrpos %u\n",prev_diagonal->univdiagonal,low));
    } else if (prev_diagonal->univdiagonal > high) {
      debug13(printf("Skipping because diagonal %u > high_chrpos %u\n",prev_diagonal->univdiagonal,high));
    } else {
      score = prev_diagonal->intscore;
      if (prev_diagonal->univdiagonal == diagonal->univdiagonal) {
	score += 1;		/* This bonus means we should double count contiguous region within each segment */
      }
      if (score <= best_score) {
	debug13(printf("Skipping because score %d <= best_score %d\n",score,best_score));
      } else {
	best_score = score;
	/* diagonal->prev = (Univdiag_T) NULL; */
	debug13(printf("Updating best score (for link to middle diagonal) to be %d\n",best_score));
      }
    }

#if 0
    diagonal->intscore = best_score + 2*diagonal->nconsecutive;  /* Rewards longer diagonals, so not desirable */
#else
    /* Exposed a bug in the algorithm (NM_080646_112_) with a NULL alt_path, but now resolved */
    diagonal->intscore = best_score;
#endif
    debug13(printf("Right diagonal %d..%d at %u (gap of %u) gets score %d\n",
		   diagonal->qstart,diagonal->qend,diagonal->univdiagonal,
		   diagonal->univdiagonal - prev_diagonal->univdiagonal,diagonal->intscore));
  }
  FREE(diagonal_array);


  /* A2.  Optimizing for dynamic programming */
  best_score_right = 0;
  *best_right_diagonals = (List_T) NULL;
  for (p = *all_right_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    debug13(printf("Considering right diagonal %d..%d at %u (vs middle diagonal %u)\n",
		   diagonal->qstart,diagonal->qend,diagonal->univdiagonal,(*middle_diagonal)->univdiagonal));
    if (diagonal->intscore > best_score_right) {
      best_score_right = diagonal->intscore;
      /* List_free(&(*best_right_diagonals)); -- allocated by Univdiagpool_push */
      *best_right_diagonals = Univdiagpool_push_existing(NULL,univdiagpool,diagonal);
    } else if (diagonal->intscore == best_score_right) {
      *best_right_diagonals = Univdiagpool_push_existing(*best_right_diagonals,univdiagpool,diagonal);
    }
  }


  /* C.  Compute left diagonals */
  /* C1.  Scoring for dynamic programming */
  diagonal_array = (Univdiag_T *) List_to_array_n(&ndiagonals,left_diagonals);
  /* List_free(&left_diagonals); -- allocated by Univdiagpool_push */
#ifdef DEBUG12
  printf("Left side before consolidating\n");
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    printf("%d..%d at %u\n",diagonal->qstart,diagonal->qend,diagonal->univdiagonal);
  }
#endif

  *all_left_diagonals = (List_T) NULL;
  qsort(diagonal_array,ndiagonals,sizeof(Univdiag_T),Univdiag_diagonal_cmp);
  i = 0;
  while (i < ndiagonals) {
    j = i;
    while (j < ndiagonals && diagonal_array[j]->univdiagonal == diagonal_array[i]->univdiagonal) {
      j++;
    }
    if (j == i) {
      *all_left_diagonals = Univdiagpool_push_existing(*all_left_diagonals,univdiagpool,diagonal_array[i]);
    } else {
      *all_left_diagonals = Univdiagpool_push(*all_left_diagonals,univdiagpool,diagonal_array[i]->qstart,
					      diagonal_array[j-1]->qend,diagonal_array[i]->univdiagonal);
#if 0
      for (k = i; k < j; k++) {
	/* Univdiag_free(&(diagonal_array[k])); -- allocated by Univdiagpool_push */
      }
#endif
    }
    i = j;
  }
  FREE(diagonal_array);

  /* TODO: May be able to skip this sorting step */
  diagonal_array = (Univdiag_T *) List_to_array_n(&ndiagonals,*all_left_diagonals);
  qsort(diagonal_array,ndiagonals,sizeof(Univdiag_T),Univdiag_descending_cmp);
#ifdef DEBUG12
  printf("Left side after consolidating and sorting\n");
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    printf("%d..%d at %u\n",diagonal->qstart,diagonal->qend,diagonal->univdiagonal);
  }
#endif


  debug13(printf("%d diagonals:\n",ndiagonals));
  for (i = 0; i < ndiagonals; i++) {
    diagonal = diagonal_array[i];
    debug13(printf("%d: %d..%d at %u\n",i,diagonal->qstart,diagonal->qend,diagonal->univdiagonal));

    low = subtract_bounded(diagonal->univdiagonal,max_insertionlen,chroffset);
    high = add_bounded(diagonal->univdiagonal,overall_max_distance_genome,chrhigh);
    qpos = diagonal->qend;
    best_score = 0;

    for (j = i - 1; j >= 0; --j) {
      prev_diagonal = diagonal_array[j];
      debug13(printf("  %d: %d..%d at %u  ",j,prev_diagonal->qstart,prev_diagonal->qend,prev_diagonal->univdiagonal));

      if (prev_diagonal->qstart <= qpos) {
	debug13(printf("Skipping because qstart %d <= qpos %d\n",prev_diagonal->qstart,qpos));
      } else if (prev_diagonal->univdiagonal < low) {
	debug13(printf("Skipping because diagonal %u < low %u\n",prev_diagonal->univdiagonal,low));
      } else if (prev_diagonal->univdiagonal > high) {
	debug13(printf("Skipping because diagonal %u > high %u\n",prev_diagonal->univdiagonal,high));
      } else {
	score = prev_diagonal->intscore;
	if (prev_diagonal->univdiagonal == diagonal->univdiagonal) {
	  score += 1;
	}
	if (score <= best_score) {
	  debug13(printf("Skipping because score %d <= best_score %d\n",score,best_score));
	} else {
	  best_score = score;
	  /* diagonal->prev = prev_diagonal; -- Not used */
	  debug13(printf("Updating best score to be %d.  Prev diagonal is %d..%d at %llu\n",
			 best_score,prev_diagonal->qstart,prev_diagonal->qend,prev_diagonal->univdiagonal));
	}
      }
    }

    /* Handle links to middle diagonal */
    prev_diagonal = *middle_diagonal;
    debug13(printf("  Middle: %d..%d at %u  ",prev_diagonal->qstart,prev_diagonal->qend,prev_diagonal->univdiagonal));
    if (prev_diagonal->qstart <= qpos) {
      debug13(printf("Skipping because qstart %d <= qpos %d\n",prev_diagonal->qstart,qpos));
    } else if (prev_diagonal->univdiagonal < low) {
      debug13(printf("Skipping because diagonal %u < low_chrpos %u\n",prev_diagonal->univdiagonal,low));
    } else if (prev_diagonal->univdiagonal > high) {
      debug13(printf("Skipping because diagonal %u > high_chrpos %u\n",prev_diagonal->univdiagonal,high));
    } else {
      score = prev_diagonal->intscore;
      if (prev_diagonal->univdiagonal == diagonal->univdiagonal) {
	score += 1;		/* This bonus means we should double count contiguous region within each segment */
      }
      if (score <= best_score) {
	debug13(printf("Skipping because score %d <= best_score %d\n",prev_diagonal->intscore,best_score));
      } else {
	best_score = score;
	/* diagonal->prev = (Univdiag_T) NULL; */
	debug13(printf("Updating best score (for link to middle diagonal) to be %d\n",best_score));
      }
    }

#if 0
    diagonal->intscore = best_score + 2*diagonal->nconsecutive;  /* Rewards longer diagonals, so not desirable */
#else
    /* Exposes a bug in the algorithm (NM_080646_112_) */
    diagonal->intscore = best_score;
#endif
    debug13(printf("Left diagonal %d..%d at %u (gap of %u) gets score %d\n",
		   diagonal->qstart,diagonal->qend,diagonal->univdiagonal,
		   prev_diagonal->univdiagonal - diagonal->univdiagonal,diagonal->intscore));
  }
  FREE(diagonal_array);


  /* C2.  Optimizing for dynamic programming */
  best_score_left = 0;
  *best_left_diagonals = (List_T) NULL;
  for (p = *all_left_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    if (diagonal->intscore > best_score_left) {
      best_score_left = diagonal->intscore;
      /* List_free(&(*best_left_diagonals)); -- allocated by Univdiagpool_push */
      *best_left_diagonals = Univdiagpool_push_existing(NULL,univdiagpool,diagonal);
    } else if (diagonal->intscore == best_score_left) {
      *best_left_diagonals = Univdiagpool_push_existing(*best_left_diagonals,univdiagpool,diagonal);
    }
  }

#if 0
  printf("Best on the left\n");
  for (p = *best_left_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    printf("Score %d: %d..%d at %u\n",diagonal->intscore,diagonal->qstart,diagonal->qend,diagonal->diagonal);
  }
#endif


  if (best_score_left == 0 && best_score_right == 0) {
    return (*middle_diagonal)->intscore;
  } else if (best_score_left == 0) {
    return best_score_right;
  } else if (best_score_right == 0) {
    return best_score_left;
  } else {
    /* middle_diagonal score is double counted */
    return best_score_left + best_score_right - (*middle_diagonal)->intscore;
  }
}
#endif



#if 0
static void
diagonals_dump (Univcoord_T *diagonals, int ndiagonals) {
  int i;

  for (i = 0; i < ndiagonals; i++) {
    printf(" %u",diagonals[i]);
  }
  printf("\n");

  return;
}
#endif


/* Caller needs to call Elt_gc(&plus_set) and Elt_gc(&minus_set) */
static void
get_elt_sets_queryfwd (List_T *plus_set, List_T *minus_set, T *best_plus_elt, T *best_minus_elt,
		       Stage1_T stage1, int querylength,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, 
		       int nmismatches_allowed, int genestrand, Listpool_T listpool) {
  T elt;
  int best_plus_nmatches, best_minus_nmatches;
  int max_plus_qpos, max_minus_qpos, plus_qpos[3], minus_qpos[3];

  int mod;
  int niter_plus, niter_minus;

#ifdef LARGE_GENOMES
  unsigned char *positions_high;
#endif
  UINT4 *positions;
  int npositions;

  int query_lastpos = querylength - index1part;
  int queryoffset;


  debug(printf("\nStarting get_elt_sets_queryfwd with querylength %d and nmismatches_allowed %d, genestrand %d\n",
	       querylength,nmismatches_allowed,genestrand));

#if 0
  /* Allow calls to Extension_search_queryfwd in addition to other methods */
  *hits_gplus = *hits_gminus = (List_T) NULL;
#endif

  if (nmismatches_allowed < 0) {
    nmismatches_allowed = 0;
#if 0
  } else {
    /* It is possible that this makes GSNAP too slow */
    nmismatches_allowed = querylength;
#endif
  }

  /* I.  Race from plus and minus start to end */
  best_plus_nmatches = best_minus_nmatches = 0;
  *best_plus_elt = *best_minus_elt = (T) NULL;
  max_plus_qpos = max_minus_qpos = 0;
  plus_qpos[0] = minus_qpos[0] = 0;
  plus_qpos[1] = minus_qpos[1] = 1;
  plus_qpos[2] = minus_qpos[2] = 2;
  niter_plus = niter_minus = 0;

  while (niter_plus <= nmismatches_allowed && niter_minus <= nmismatches_allowed &&
	 max_plus_qpos < query_lastpos && max_minus_qpos < query_lastpos) {
    for (mod = 0; mod < 3; mod++) {
      if ((queryoffset = plus_qpos[mod]) >= query_lastpos) {
	/* Skip */
	plus_qpos[mod] += 3;
      } else if (stage1->plus_validp[queryoffset] == false) {
	/* Skip */
	plus_qpos[mod] += 3;
      } else {
	if (stage1->plus_retrievedp[queryoffset] == true) {
#ifdef LARGE_GENOMES
	  positions_high = stage1->plus_positions_high[queryoffset];
#endif
	  positions = stage1->plus_positions[queryoffset];
	  npositions = stage1->plus_npositions[queryoffset];
	  debug(printf("Already have %d plus positions at %d\n",npositions,queryoffset));
	} else {
	  /* These should be lefts */
	  assert(stage1->plus_positions[queryoffset] == NULL);
	  stage1->plus_retrievedp[queryoffset] = true;
#ifdef LARGE_GENOMES
	  npositions = stage1->plus_npositions[queryoffset] =
	    Indexdb_largeptr_with_diagterm(&stage1->plus_positions_high[queryoffset],&stage1->plus_positions[queryoffset],
					   indexdb,stage1->plus_oligos[queryoffset],/*diagterm*/-queryoffset);
	  positions_high = stage1->plus_positions_high[queryoffset];
#else
	  npositions = stage1->plus_npositions[queryoffset] =
	    Indexdb_ptr_with_diagterm(&stage1->plus_positions[queryoffset],indexdb,
				      stage1->plus_oligos[queryoffset],/*diagterm*/-queryoffset);
#endif
	  positions = stage1->plus_positions[queryoffset];
	}

	if ((elt = Elt_read_queryfwd(
#ifdef LARGE_GENOMES
				     positions_high,
#endif
				     positions,npositions,/*diagterm*/-queryoffset,
				     querylength,/*querystart*/queryoffset,
				     query_compress_fwd,/*plusp*/true,genestrand)) == NULL) {
	  plus_qpos[mod] += 3;
	} else {
	  if (elt->nmatches > best_plus_nmatches && elt->n_all_diagonals <= MAX_HITS_FOR_BEST_ELT) {
	    *best_plus_elt = elt;
	    best_plus_nmatches = elt->nmatches;
	  }
	  *plus_set = Listpool_push(*plus_set,listpool,(void *) elt);
	  plus_qpos[mod] += elt->nmatches;
	  niter_plus++;
	}
      }

      if (plus_qpos[mod] > max_plus_qpos) {
	max_plus_qpos = plus_qpos[mod];
      }


      /* querypos_rc uses standard Stage1_T convention, but we have switched to sarray convention */
      /* querypos_rc = query_lastpos - queryoffset; */
      if ((queryoffset = minus_qpos[mod]) >= query_lastpos) {
	/* Skip */
	minus_qpos[mod] += 3;
      } else if (stage1->minus_validp[queryoffset] == false) {
	/* Skip */
	minus_qpos[mod] += 3;
      } else {
	if (stage1->minus_retrievedp[queryoffset] == true) {
#ifdef LARGE_GENOMES
	  positions_high = stage1->minus_positions_high[queryoffset];
#endif
	  positions = stage1->minus_positions[queryoffset];
	  npositions = stage1->minus_npositions[queryoffset];
	  debug(printf("Already have %d minus positions at %d\n",npositions,queryoffset));
	} else {
	  /* These should be lefts */
	  assert(stage1->minus_positions[queryoffset] == NULL);
	  stage1->minus_retrievedp[queryoffset] = true;
	  /* In standard stage1 convention, diagterm would be
	     -querylength + querypos + index1part = -(query_lastpos - querypos) = -queryoffset */
#ifdef LARGE_GENOMES
	  npositions = stage1->minus_npositions[queryoffset] =
	    Indexdb_largeptr_with_diagterm(&stage1->minus_positions_high[queryoffset],&stage1->minus_positions[queryoffset],
					   indexdb2,stage1->minus_oligos[queryoffset],/*diagterm*/-queryoffset);
	  positions_high = stage1->minus_positions_high[queryoffset];
#else
	  npositions = stage1->minus_npositions[queryoffset] =
	    Indexdb_ptr_with_diagterm(&stage1->minus_positions[queryoffset],indexdb2,
				      stage1->minus_oligos[queryoffset],/*diagterm*/-queryoffset);
#endif
	  positions = stage1->minus_positions[queryoffset];
	}
	
	if ((elt = Elt_read_queryfwd(
#ifdef LARGE_GENOMES
				     positions_high,
#endif
				     positions,npositions,/*diagterm*/-queryoffset,
				     querylength,/*querystart*/queryoffset,
				     query_compress_rev,/*plusp*/false,genestrand)) == NULL) {
	  minus_qpos[mod] += 3;
	} else {
	  if (elt->nmatches > best_minus_nmatches && elt->n_all_diagonals <= MAX_HITS_FOR_BEST_ELT) {
	    *best_minus_elt = elt;
	    best_minus_nmatches = elt->nmatches;
	  }
	  *minus_set = Listpool_push(*minus_set,listpool,(void *) elt);
	  minus_qpos[mod] += elt->nmatches;
	  niter_minus++;
	}
      }
	
      if (minus_qpos[mod] > max_minus_qpos) {
	max_minus_qpos = minus_qpos[mod];
      }
    }

#ifdef DEBUG
    printf("\n");
    for (mod = 0; mod < 3; mod++) {
      printf("mod %d, plus_qpos %d\n",mod,plus_qpos[mod]);
      printf("mod %d, minus_qpos %d\n",mod,minus_qpos[mod]);
    }
    printf("max_plus_qpos %d\n",max_plus_qpos);
    printf("max_maxus_qpos %d\n",max_minus_qpos);
    printf("\n");
#endif

#ifndef LIBERAL
    /* Skip the presumed mismatch */
    max_plus_qpos += 1;
    max_minus_qpos += 1;
#endif

    plus_qpos[0] = max_plus_qpos;
    plus_qpos[1] = max_plus_qpos + 1;
    plus_qpos[2] = max_plus_qpos + 2;

    minus_qpos[0] = max_minus_qpos;
    minus_qpos[1] = max_minus_qpos + 1;
    minus_qpos[2] = max_minus_qpos + 2;
  }

  *plus_set = List_reverse(*plus_set);
  *minus_set = List_reverse(*minus_set);

#ifdef DEBUG
  printf("queryfwd plus set:\n");
  Elt_dump_set(*plus_set);
  printf("\n");
  printf("queryfwd minus set:\n");
  Elt_dump_set(*minus_set);
  printf("\n");
#endif

  if (max_minus_qpos >= query_lastpos) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYFWD PLUS: minus side won, so skip plus side\n"));
    *best_plus_elt = (Elt_T) NULL;

  } else if (max_plus_qpos < query_lastpos || *best_plus_elt == NULL) {
    debug(printf("QUERYFWD PLUS: Still could not find large pieces: plus_qpos %d < query_lastpos %d, best_plus_elt %p\n",
		 max_plus_qpos,query_lastpos,*best_plus_elt));
    *best_plus_elt = (Elt_T) NULL;

  } else if ((*best_plus_elt)->ndiagonals == 0) {
      /* Could happen if there are too many diagonals */
    debug(printf("QUERYFWD PLUS: Best elt has no diagonals\n"));
    *best_plus_elt = (Elt_T) NULL;

  } else {
    debug(printf("QUERYFWD PLUS HAS BEST ELT: "));
    debug(Elt_dump(*best_plus_elt));
    debug(printf("\n"));
  }

  if (max_plus_qpos >= query_lastpos) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYFWD MINUS: plus side won, so skip minus side\n"));
    *best_minus_elt = (Elt_T) NULL;

  } else if (max_minus_qpos < query_lastpos || *best_minus_elt == NULL) {
    debug(printf("QUERYFWD MINUS: Still could not find large pieces: minus_qpos %d < query_lastpos %d, best_minus_elt %p\n",
		 max_minus_qpos,query_lastpos,*best_minus_elt));
    *best_minus_elt = (Elt_T) NULL;

  } else if ((*best_minus_elt)->ndiagonals == 0) {
    /* Could happen if there are too many diagonals */
    debug(printf("QUERYFWD MINUS: Best elt has no diagonals\n"));
    *best_minus_elt = (Elt_T) NULL;

  } else {
    debug(printf("QUERYFWD MINUS HAS BEST ELT: "));
    debug(Elt_dump(*best_minus_elt));
    debug(printf("\n"));
  }

  return;
}


/* Note: A second try, starting at queryrev, may solve only a few more
   cases presented to this method (i.e., cases that cannot be solved
   by kmer-end search).  It therefore may be best to push these cases
   to the next method. */

/* Caller needs to call Elt_gc(&plus_set) and Elt_gc(&minus_set) */
static void
get_elt_sets_queryrev (List_T *plus_set, List_T *minus_set, T *best_plus_elt, T *best_minus_elt,
		       Stage1_T stage1, int querylength,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, 
		       int nmismatches_allowed, int genestrand, Listpool_T listpool) {
  T elt;
  int best_plus_nmatches, best_minus_nmatches;
  int min_plus_qpos, min_minus_qpos, plus_qpos[3], minus_qpos[3];

  int mod;
  int niter_plus, niter_minus;

#ifdef LARGE_GENOMES
  unsigned char *positions_high;
#endif
  UINT4 *positions;
  int npositions;

  int query_lastpos = querylength - index1part;
  int queryoffset;

  debug(printf("\nStarting get_elt_sets_queryrev with querylength %d and nmismatches_allowed %d, genestrand %d\n",
	       querylength,nmismatches_allowed,genestrand));

#if 0
  /* Allow calls to Extension_search_queryrev in addition to other methods */
  *hits_gplus = *hits_gminus = (List_T) NULL;
#endif

  if (nmismatches_allowed < 0) {
    nmismatches_allowed = 0;
#if 0
  } else {
    /* It is possible that this makes GSNAP too slow */
    nmismatches_allowed = querylength;
#endif
  }

  /* I.  Race from plus and minus end to start */
  best_plus_nmatches = best_minus_nmatches = 0;
  *best_plus_elt = *best_minus_elt = (T) NULL;
  min_plus_qpos = min_minus_qpos = query_lastpos;
  plus_qpos[0] = minus_qpos[0] = query_lastpos;
  plus_qpos[1] = minus_qpos[1] = query_lastpos-1;
  plus_qpos[2] = minus_qpos[2] = query_lastpos-2;
  niter_plus = niter_minus = 0;

  while (niter_plus <= nmismatches_allowed && niter_minus <= nmismatches_allowed &&
	 min_plus_qpos > 0 && min_minus_qpos > 0) {
    for (mod = 0; mod < 3; mod++) {
      if ((queryoffset = plus_qpos[mod]) < 0) {
	/* Skip */
	plus_qpos[mod] -= 3;
      } else if (stage1->plus_validp[queryoffset] == false) {
	/* Skip */
	plus_qpos[mod] -= 3;
      } else {
	if (stage1->plus_retrievedp[queryoffset] == true) {
#ifdef LARGE_GENOMES
	  positions_high = stage1->plus_positions_high[queryoffset];
#endif
	  positions = stage1->plus_positions[queryoffset];
	  npositions = stage1->plus_npositions[queryoffset];
	} else {
	  /* These should be lefts */
	  assert(stage1->plus_positions[queryoffset] == NULL);
	  stage1->plus_retrievedp[queryoffset] = true;
#ifdef LARGE_GENOMES
	  npositions = stage1->plus_npositions[queryoffset] =
	    Indexdb_largeptr_with_diagterm(&stage1->plus_positions_high[queryoffset],&stage1->plus_positions[queryoffset],
					   indexdb,stage1->plus_oligos[queryoffset],/*diagterm*/-queryoffset);
	  positions_high = stage1->plus_positions_high[queryoffset];
#else
	  npositions = stage1->plus_npositions[queryoffset] =
	    Indexdb_ptr_with_diagterm(&stage1->plus_positions[queryoffset],indexdb,
				      stage1->plus_oligos[queryoffset],/*diagterm*/-queryoffset);
#endif
	  positions = stage1->plus_positions[queryoffset];
	}
	
	if ((elt = Elt_read_queryrev(
#ifdef LARGE_GENOMES
				     positions_high,
#endif
				     positions,npositions,/*diagterm*/-queryoffset,
				     /*queryend*/queryoffset + index1part,
				     query_compress_fwd,/*plusp*/true,genestrand)) == NULL) {
	  plus_qpos[mod] -= 3;
	} else {
	  if (elt->nmatches > best_plus_nmatches && elt->n_all_diagonals <= MAX_HITS_FOR_BEST_ELT) {
	    *best_plus_elt = elt;
	    best_plus_nmatches = elt->nmatches;
	  }
	  *plus_set = Listpool_push(*plus_set,listpool,(void *) elt);
	  plus_qpos[mod] -= elt->nmatches;
	  niter_plus++;
	}
      }

      if (plus_qpos[mod] < min_plus_qpos) {
	min_plus_qpos = plus_qpos[mod];
      }


      /* querypos_rc uses standard Stage1_T convention, but we have switched to sarray convention */
      /* querypos_rc = query_lastpos - queryoffset; */
      if ((queryoffset = minus_qpos[mod]) < 0) {
	/* Skip */
	minus_qpos[mod] -= 3;
      } else if (stage1->minus_validp[queryoffset] == false) {
	/* Skip */
	minus_qpos[mod] -= 3;
      } else {
	if (stage1->minus_retrievedp[queryoffset] == true) {
#ifdef LARGE_GENOMES
	  positions_high = stage1->minus_positions_high[queryoffset];
#endif
	  positions = stage1->minus_positions[queryoffset];
	  npositions = stage1->minus_npositions[queryoffset];
	} else {
	  /* These should be lefts */
	  assert(stage1->minus_positions[queryoffset] == NULL);
	  stage1->minus_retrievedp[queryoffset] = true;
	  /* In standard stage1 convention, diagterm would be
	     -querylength + querypos + index1part = -(query_lastpos - querypos) = -queryoffset */
#ifdef LARGE_GENOMES
	  npositions = stage1->minus_npositions[queryoffset] =
	    Indexdb_largeptr_with_diagterm(&stage1->minus_positions_high[queryoffset],&stage1->minus_positions[queryoffset],
					   indexdb2,stage1->minus_oligos[queryoffset],/*diagterm*/-queryoffset);
	  positions_high = stage1->minus_positions_high[queryoffset];
#else
	  npositions = stage1->minus_npositions[queryoffset] =
	    Indexdb_ptr_with_diagterm(&stage1->minus_positions[queryoffset],indexdb2,
				      stage1->minus_oligos[queryoffset],/*diagterm*/-queryoffset);
#endif
	  positions = stage1->minus_positions[queryoffset];
	}
	
	if ((elt = Elt_read_queryrev(
#ifdef LARGE_GENOMES
				     positions_high,
#endif
				     positions,npositions,/*diagterm*/-queryoffset,
				     /*queryend*/queryoffset + index1part,
				     query_compress_rev,/*plusp*/false,genestrand)) == NULL) {
	  minus_qpos[mod] -= 3;
	} else {
	  if (elt->nmatches > best_minus_nmatches && elt->n_all_diagonals <= MAX_HITS_FOR_BEST_ELT) {
	    *best_minus_elt = elt;
	    best_minus_nmatches = elt->nmatches;
	  }
	  *minus_set = Listpool_push(*minus_set,listpool,(void *) elt);
	  minus_qpos[mod] -= elt->nmatches;
	  niter_minus++;
	}
      }

      if (minus_qpos[mod] < min_minus_qpos) {
	min_minus_qpos = minus_qpos[mod];
      }
    }

#ifdef DEBUG
    printf("\n");
    for (mod = 0; mod < 3; mod++) {
      printf("mod %d, plus_qpos %d\n",mod,plus_qpos[mod]);
      printf("mod %d, minus_qpos %d\n",mod,minus_qpos[mod]);
    }
    printf("min_plus_qpos %d\n",min_plus_qpos);
    printf("min_minus_qpos %d\n",min_minus_qpos);
    printf("\n");
#endif

#ifndef LIBERAL
    /* Skip the presumed mismatch */
    min_plus_qpos -= 1;
    min_minus_qpos -= 1;
#endif

    plus_qpos[0] = min_plus_qpos;
    plus_qpos[1] = min_plus_qpos - 1;
    plus_qpos[2] = min_plus_qpos - 2;

    minus_qpos[0] = min_minus_qpos;
    minus_qpos[1] = min_minus_qpos - 1;
    minus_qpos[2] = min_minus_qpos - 2;
  }

#if 0
  /* Not needed for queryrev */
  *plus_set = List_reverse(*plus_set);
  *minus_set = List_reverse(*minus_set);
#endif

#ifdef DEBUG
  printf("queryrev plus set:\n");
  Elt_dump_set(*plus_set);
  printf("\n");
  printf("queryrev minus set:\n");
  Elt_dump_set(*minus_set);
  printf("\n");
#endif

  if (min_minus_qpos <= 0) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYREV PLUS: minus side won, so skip plus side\n"));
    *best_plus_elt = (Elt_T) NULL;

  } else if (min_plus_qpos > 0 || *best_plus_elt == NULL) {
    debug(printf("QUERYREV PLUS: Still could not find large pieces: plus_qpos %d > 0, best_plus_elt %p\n",
		 min_plus_qpos,*best_plus_elt));
    *best_plus_elt = (Elt_T) NULL;

  } else if ((*best_plus_elt)->ndiagonals == 0) {
    /* Could happen if there are too many diagonals */
    debug(printf("QUERYREV PLUS: Best elt has no diagonals\n"));
    *best_plus_elt = (Elt_T) NULL;

  } else {
    debug(printf("QUERYREV PLUS HAS BEST ELT: "));
    debug(Elt_dump(*best_plus_elt));
    debug(printf("\n"));
  }

  if (min_plus_qpos <= 0) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYREV MINUS: plus side won, so skip minus side\n"));
    *best_minus_elt = (Elt_T) NULL;

  } else if (min_minus_qpos > 0 || *best_minus_elt == NULL) {
    debug(printf("QUERYREV MINUS: Still could not find large pieces: minus_qpos %d > 0, best_minus_elt %p\n",
		 min_minus_qpos,*best_minus_elt));
    *best_minus_elt = (Elt_T) NULL;

  } else if ((*best_minus_elt)->ndiagonals == 0) {
    /* Could happen if there are too many diagonals */
    debug(printf("QUERYREV MINUS: Best elt has no diagonals\n"));
    *best_minus_elt = (Elt_T) NULL;
    
  } else {
    debug(printf("QUERYREV MINUS HAS BEST ELT: "));
    debug(Elt_dump(*best_minus_elt));
    debug(printf("\n"));
  }

  return;
}


static List_T
process_seed (int *found_score_overall, int *found_score_within_trims,
	      Chrnum_T *chrnum, Univcoord_T *chroffset, Univcoord_T *chrhigh,
	      Chrpos_T *chrlength, List_T hits,
	      Univcoord_T middle_univdiagonal, int middle_qstart, int middle_qend,
	      List_T queryfwd_set, List_T queryrev_set, T queryfwd_best_elt, T queryrev_best_elt,
	      Stage1_T stage1, char *queryptr, int querylength, Compress_T query_compress,
	      int nmismatches_allowed, bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
	      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	      Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
	      int level) {

  List_T left_diagonals = NULL, right_diagonals = NULL, p;
  Univdiag_T *left_diagonal_array, *right_diagonal_array;
  Univcoord_T univdiagonal, start_low, start_high, end_low, end_high;
  int qstart, qend;
  Elt_T elt;
  int nleft, nright, i, j, k;
  bool foundp;


  debug(printf("PROCESS SEED at %u, qstart %d, qend %d\n",middle_univdiagonal,middle_qstart,middle_qend));
  *chrhigh = Univ_IIT_update_chrnum(&(*chrnum),&(*chroffset),*chrhigh,&(*chrlength),chromosome_iit,
				    middle_univdiagonal,querylength,circular_typeint);
  debug(printf("Got chrnum %d\n",*chrnum));
  if (*chroffset > middle_univdiagonal + middle_qstart) {
    middle_qstart = (int) ((*chroffset) - middle_univdiagonal);
    debug(printf("Updated middle_qstart to be %d\n",middle_qstart));
  }
  if (middle_univdiagonal + middle_qend >= *chrhigh) {
    middle_qend = (int) ((*chrhigh) - middle_univdiagonal);
    debug(printf("Updated middle_qend to be %d\n",middle_qend));
  }

  start_low = subtract_bounded(middle_univdiagonal,/*minusterm*/overall_end_distance_genome,*chroffset);
  start_high = add_bounded(middle_univdiagonal,/*plusterm*/max_insertionlen_default,*chrhigh);

  end_low = subtract_bounded(middle_univdiagonal,/*minusterm*/max_insertionlen_default,*chroffset);
  end_high = add_bounded(middle_univdiagonal,/*plusterm*/overall_end_distance_genome,*chrhigh);

  debug(printf("Computing start %u..%u and end %u..%u.  Chromosome bounds: %u..%u\n",
	       start_low,start_high,end_low,end_high,*chroffset,*chrhigh));


  for (p = queryfwd_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
    if (elt != queryfwd_best_elt && elt != queryrev_best_elt) {
      if (elt_startp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_diagonals(elt,start_low,start_high);
	for (k = 0; k < elt->ndiagonals; k++) {
	  left_diagonals = Univdiagpool_push(left_diagonals,univdiagpool,elt->qstart,elt->qend,elt->diagonals[k]);
	}
      }
      
      if (elt_endp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_diagonals(elt,end_low,end_high);
	for (k = 0; k < elt->ndiagonals; k++) {
	  right_diagonals = Univdiagpool_push(right_diagonals,univdiagpool,elt->qstart,elt->qend,elt->diagonals[k]);
	}
      }
    }
  }
  
  for (p = queryrev_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
    if (elt != queryfwd_best_elt && elt != queryrev_best_elt) {
      if (elt_startp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_diagonals(elt,start_low,start_high);
	for (k = 0; k < elt->ndiagonals; k++) {
	  left_diagonals = Univdiagpool_push(left_diagonals,univdiagpool,elt->qstart,elt->qend,elt->diagonals[k]);
	}
      }
      
      if (elt_endp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_diagonals(elt,end_low,end_high);
	for (k = 0; k < elt->ndiagonals; k++) {
	  right_diagonals = Univdiagpool_push(right_diagonals,univdiagpool,elt->qstart,elt->qend,elt->diagonals[k]);
	}
      }
    }
  }

  if (left_diagonals != NULL) {
    left_diagonal_array = (Univdiag_T *) List_to_array_n(&nleft,left_diagonals);
    qsort(left_diagonal_array,nleft,sizeof(Univdiag_T),Univdiag_diagonal_cmp);
    /* List_free(&left_diagonals); -- allocated by Univdiagpool_push */

    left_diagonals = (List_T) NULL;
    i = 0;
    while (i < nleft) {
      univdiagonal = left_diagonal_array[i]->univdiagonal;
      qstart = left_diagonal_array[i]->qstart;
      qend = left_diagonal_array[i]->qend;
      
      j = i+1;
      while (j < nleft && left_diagonal_array[j]->univdiagonal == univdiagonal) {
	if (left_diagonal_array[j]->qstart < qstart) {
	  qstart = left_diagonal_array[j]->qstart;
	}
	if (left_diagonal_array[j]->qend > qend) {
	  qend = left_diagonal_array[j]->qend;
	}
	j++;
      }
      
      left_diagonal_array[i]->qstart = qstart;
      left_diagonal_array[i]->qend = qend;
      left_diagonals = Univdiagpool_push_existing(left_diagonals,univdiagpool,left_diagonal_array[i]);
      i = j;
    }
    FREE(left_diagonal_array);
  }

  if (right_diagonals != NULL) {
    /* Sort the right diagonals and unique them */
    right_diagonal_array = (Univdiag_T *) List_to_array_n(&nright,right_diagonals);
    qsort(right_diagonal_array,nright,sizeof(Univdiag_T),Univdiag_diagonal_cmp);
    /* List_free(&right_diagonals); -- allocated by Univdiagpool_push */
    
    right_diagonals = (List_T) NULL;
    i = 0;
    while (i < nright) {
      univdiagonal = right_diagonal_array[i]->univdiagonal;
      qstart = right_diagonal_array[i]->qstart;
      qend = right_diagonal_array[i]->qend;
      
      j = i+1;
      while (j < nright && right_diagonal_array[j]->univdiagonal == univdiagonal) {
	if (right_diagonal_array[j]->qstart < qstart) {
	  qstart = right_diagonal_array[j]->qstart;
	}
	if (right_diagonal_array[j]->qend > qend) {
	  qend = right_diagonal_array[j]->qend;
	}
	j++;
      }
      
      right_diagonal_array[i]->qstart = qstart;
      right_diagonal_array[i]->qend = qend;
      right_diagonals = Univdiagpool_push_existing(right_diagonals,univdiagpool,right_diagonal_array[i]);
      i = j;
    }
    FREE(right_diagonal_array);
  }


  hits = Path_solve_from_diagonals(&foundp,&(*found_score_overall),&(*found_score_within_trims),hits,
				   middle_univdiagonal,middle_qstart,middle_qend,
				   right_diagonals,left_diagonals,queryptr,querylength,
				   stage1->mismatch_positions_alloc,stage1->spliceinfo,
				   stage1->stream_alloc,stage1->streamsize_alloc,query_compress,
				   *chrnum,*chroffset,*chrhigh,*chrlength,plusp,genestrand,
				   nmismatches_allowed,paired_end_p,first_read_p,
				   intlistpool,univcoordlistpool,listpool,univdiagpool,
				   hitlistpool,/*method*/EXT,level);

  /* List_free(&right_diagonals); -- allocated by Univdiagpool_push */
  /* List_free(&left_diagonals); -- allocated by Univdiagpool_push */

  return hits;
}


#define min(a,b) (a < b) ? a : b
#define max(a,b) (a > b) ? a : b


static List_T
extend_seeds (int *found_score_overall, int *found_score_within_trims,
	      List_T hits, T queryfwd_best_elt, T queryrev_best_elt,
	      List_T queryfwd_set, List_T queryrev_set,

	      Stage1_T stage1, char *queryptr, int querylength, Compress_T query_compress,
	      int nmismatches_allowed, bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
	      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	      Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
	      int level) {

  Univcoord_T queryfwd_left, queryrev_left;
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh, univdiagonal;
  Chrpos_T chrlength;
  int qstart, qend;
  int queryfwd_nseeds, queryrev_nseeds, i, j;


  chrnum = 1;
  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,/*chrnum*/1,circular_typeint);

  if (queryfwd_best_elt == NULL) {
    queryfwd_nseeds = 0;
  } else {
    queryfwd_nseeds = queryfwd_best_elt->n_all_diagonals;
  }

  if (queryrev_best_elt == NULL) {
    queryrev_nseeds = 0;
  } else {
    queryrev_nseeds = queryrev_best_elt->n_all_diagonals;
  }
    
  /* Go through the union of the seeds */
  i = j = 0;
  while (i < queryfwd_nseeds && j < queryrev_nseeds) {
    queryfwd_left = queryfwd_best_elt->all_diagonals[i];
    queryrev_left = queryrev_best_elt->all_diagonals[j];

    if (queryfwd_left == queryrev_left) {
      /* Combine the seeds */
      univdiagonal = queryfwd_left;
      qstart = min(queryfwd_best_elt->qstart,queryrev_best_elt->qstart);
      qend = max(queryfwd_best_elt->qend,queryrev_best_elt->qend);
      i++; j++;

    } else if (queryfwd_left < queryrev_left) {
      univdiagonal = queryfwd_left;
      qstart = queryfwd_best_elt->qstart;
      qend = queryfwd_best_elt->qend;
      i++;

    } else {
      univdiagonal = queryrev_left;
      qstart = queryrev_best_elt->qstart;
      qend = queryrev_best_elt->qend;
      j++;
    }

    /* Disallow diagonals from queryfwd_best_elt or queryrev_best_elt */
    hits = process_seed(&(*found_score_overall),&(*found_score_within_trims),
			&chrnum,&chroffset,&chrhigh,&chrlength,hits,univdiagonal,qstart,qend,
			queryfwd_set,queryrev_set,queryfwd_best_elt,queryrev_best_elt,
			stage1,queryptr,querylength,
			query_compress,nmismatches_allowed,plusp,genestrand,paired_end_p,first_read_p,
			intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,level);
  }

  if (i < queryfwd_nseeds) {
    qstart = queryfwd_best_elt->qstart;
    qend = queryfwd_best_elt->qend;
    while (i < queryfwd_nseeds) {
      univdiagonal = queryfwd_best_elt->all_diagonals[i++];
      
      /* Allow diagonals from queryrev_best_elt */
      hits = process_seed(&(*found_score_overall),&(*found_score_within_trims),
			  &chrnum,&chroffset,&chrhigh,&chrlength,hits,univdiagonal,qstart,qend,
			  queryfwd_set,queryrev_set,queryfwd_best_elt,/*queryrev_best_elt*/NULL,
			  stage1,queryptr,querylength,
			  query_compress,nmismatches_allowed,plusp,genestrand,paired_end_p,first_read_p,
			  intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,level);
    }
  }

  if (j < queryrev_nseeds) {
    qstart = queryrev_best_elt->qstart;
    qend = queryrev_best_elt->qend;
    while (j < queryrev_nseeds) {
      univdiagonal = queryrev_best_elt->all_diagonals[j++];
      
      /* Allow diagonals from queryfwd_best_elt */
      hits = process_seed(&(*found_score_overall),&(*found_score_within_trims),
			  &chrnum,&chroffset,&chrhigh,&chrlength,hits,univdiagonal,qstart,qend,
			  queryfwd_set,queryrev_set,/*queryfwd_best_elt*/NULL,queryrev_best_elt,
			  stage1,queryptr,querylength,
			  query_compress,nmismatches_allowed,plusp,genestrand,paired_end_p,first_read_p,
			  intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,level);
    }
  }

  return hits;
}


void
Extension_search (int *found_score_overall, int *found_score_within_trims,
		  List_T *hits_gplus, List_T *hits_gminus, Stage1_T stage1,
		  char *queryuc_ptr, char *queryrc, int querylength,
		  Compress_T query_compress_fwd, Compress_T query_compress_rev, 
#if 0
		  Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		  Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		  Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
		  int nmismatches_allowed, int genestrand, bool paired_end_p, bool first_read_p,
		  int level) {

  T queryfwd_best_plus_elt, queryrev_best_plus_elt, queryfwd_best_minus_elt, queryrev_best_minus_elt;


  get_elt_sets_queryfwd(&stage1->queryfwd_plus_set,&stage1->queryfwd_minus_set,
			&queryfwd_best_plus_elt,&queryfwd_best_minus_elt,
			stage1,querylength,query_compress_fwd,query_compress_rev,
			nmismatches_allowed,genestrand,listpool);
  
  get_elt_sets_queryrev(&stage1->queryrev_plus_set,&stage1->queryrev_minus_set,
			&queryrev_best_plus_elt,&queryrev_best_minus_elt,
			stage1,querylength,query_compress_fwd,query_compress_rev,
			nmismatches_allowed,genestrand,listpool);

  *hits_gplus = extend_seeds(&(*found_score_overall),&(*found_score_within_trims),*hits_gplus,
			     queryfwd_best_plus_elt,queryrev_best_plus_elt,
			     stage1->queryfwd_plus_set,stage1->queryrev_plus_set,
			     stage1,/*queryptr*/queryuc_ptr,querylength,
			     /*query_compress*/query_compress_fwd,nmismatches_allowed,
			     /*plusp*/true,genestrand,paired_end_p,first_read_p,
			     intlistpool,univcoordlistpool,listpool,univdiagpool,
			     hitlistpool,level);

  *hits_gminus = extend_seeds(&(*found_score_overall),&(*found_score_within_trims),*hits_gminus,
			      queryfwd_best_minus_elt,queryrev_best_minus_elt,
			      stage1->queryfwd_minus_set,stage1->queryrev_minus_set,
			      stage1,/*queryptr*/queryrc,querylength,
			      /*query_compress*/query_compress_rev,nmismatches_allowed,
			      /*plusp*/false,genestrand,paired_end_p,first_read_p,
			      intlistpool,univcoordlistpool,listpool,univdiagpool,
			      hitlistpool,level);

  return;
}


