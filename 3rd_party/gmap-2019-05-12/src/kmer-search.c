static char rcsid[] = "$Id: kmer-search.c 218524 2019-03-04 22:47:09Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "kmer-search.h"

#include "assert.h"
#include "mem.h"
#include "bool.h"
#include "types.h"
#include "chrnum.h"
#include "reader.h"
#include "oligo.h"

#ifndef LARGE_GENOMES
#include "merge-diagonals-simd-uint4.h"
#elif !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
#include "merge-diagonals-heap.h" /* For Merge_diagonals_large */
#include "merge-diagonals-simd-uint4.h"
#else
#include "merge-diagonals-simd-uint8.h" /* For Merge_diagonals_large */
#include "merge-diagonals-simd-uint4.h"
#endif

#ifdef LARGE_GENOMES
#include "intersect-large.h"
#endif
#include "intersect.h"
#include "transcript.h"
#include "popcount.h"
#include "genome128_hr.h"
#include "substring.h"
#include "junction.h"

#include "univdiag.h"
#include "univdiagdef.h"
#include "genome128_consec.h"
#include "splice.h"
#include "indel.h"
#include "intron.h"
#include "maxent_hr.h"
#include "knownsplicing.h"
#include "sedgesort.h"


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


/* General flow */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Processing merged diagonals */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Exon and path overlap */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Transcriptome and one end */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Transcriptome and one end details*/
#ifdef DEBUG2A
#define debug2a(x) x
#else
#define debug2a(x)
#endif

/* One end */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Both_ends */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Both end details */
#ifdef DEBUG4A
#define debug4a(x) x
#else
#define debug4a(x)
#endif


/* Trimming */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Oligos and remapping */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Binary search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif


/* Known splicing for genome */
#ifdef DEBUG4S
#define debug4s(x) x
#else
#define debug4s(x)
#endif

/* Find best path genome */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif


/* Merging faster than count table */
#define USE_MERGE 1

#define MAX_NEIGHBORS 3		/* Cannot be 0 */
#define SUBOPT 3

#define LONG_END 6
#define ALLOWED_END_MISMATCHES 2 /* For long ends */
#define ALLOWED_TRANSCRIPTOME_TRIM 3


static Mode_T mode;

static Univ_IIT_T chromosome_iit;
static int circular_typeint;
static bool *circularp;

static Genome_T genomebits;
static Genome_T genomebits_alt;

static Indexdb_T indexdb;
static Indexdb_T indexdb2;
static Indexdb_T indexdb_tr;

static bool novelsplicingp;
static Univcoord_T *splicesites;
static Splicetype_T *splicetypes;
static Chrpos_T *splicedists;
static int nsplicesites = 0;


static Chrpos_T min_intronlength;
static Chrpos_T max_deletionlen;
static Chrpos_T shortsplicedist;
static Chrpos_T overall_max_distance;

static int index1part_tr;
static int index1part;
static int index1interval;
static int local1part;

static int leftreadshift;
static Oligospace_T oligobase_mask;


/* Indexdb: oligo + diagterm -> streams of diagonals */
/* Merge: an array of diagonals with duplicates */
/* Path: genomic endpoints + gaps + trnums */

/* All calls to Substring_new are for transcriptome.  May need to make call to Univ_IIT_update_chrnum */

/* Transcriptome ends */
/* Ultrafast: check ends only.  Require extension to ends and nmismatches <= nmismatches_allowed */

/* Transcriptome complete */
/* Most prevalent: merged diagonals -> most prevalent diagonals (called loci), unique */

/* Algorithm: Check for indels.  Require extension to ends - 3 bp.
   Create lefts,querystarts,queryends,adjustments.  Convert to genomic
   paths */

/* single_hits: paths -> hits */

/* Genome */
/* Ultrafast: check ends only */
/* find_local_sets: merged->diagonals -> middle_diagonal, left_diagonals, right_diagonals */
/* Algorithm 1a: find_best_path_genome -> complete_path (list of Univdiag_T) */
/* Algorithm 1b: Kmer_search_genome (solve_via_segments_genome): complete_path -> hits */
/* Algorithm 2: Stage3end_complete_path_run_gmap: complete_path -> hits */


#ifdef DEBUG8
#ifdef HAVE_SSE4_1
static void
print_vector_dec (__m128i x) {
  printf("%d %d %d %d\n",
	 _mm_extract_epi32(x,0),_mm_extract_epi32(x,1),_mm_extract_epi32(x,2),_mm_extract_epi32(x,3));
  return;
}
#endif
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



#if 0
/* poly_p needs to be defined from [-index1part] through [querylength] */
static int
trim_ends (int *trim5, int *trim3, int *nmismatches5, int *nmismatches3,
	   bool *poly_p, Compress_T query_compress, Univcoord_T left, int querylength,
	   Genome_T omebits, bool plusp, int index1part) {
  int pos, i;
  int nmismatches;
  /* int queryend = querylength; */
  int *mismatch_positions_alloc, *mismatch_positions, *ptr, *end, cmp;

  __m128i _mismatchp, _nmatches, _values, _prev, _next, _threshold, _zero;

  /* n = nmismatches.  L = querylength (at mismatch_positions[n]) */

  /*                  X X X -1  mm_0                ... mm_{n-1} L X X X */

  /*  For trim 5 */
  /*  Values                    mm_0 mm_1 mm_2 mm_3 */
  /*  Prev                        -1 mm_0 mm_2 mm_3 */
  /*  So ptr = &(nmismatch_positions[0] */

  /*  Values                                                     L X X X */
  /*  Prev                                                mm_{n-1} L X X */
  /*  So ptr <= &(nmismatch_positions[n]) */

  /*  For trim 3 */
  /*  Values                                mm_{n-4} mm_{n-3} mm_{n-2} mm_{n-1} */
  /*  Next                                  mm_{n-3} mm_{n-2} mm_{n-1}        L */
  /*  So ptr = &(nmismatch_positions[n-4]) */

  /*  Values          X X  X   -1 */
  /*  Next            X X -1 mm_0 */
  /*  So ptr >= &(nmismatch_positions[-4]) */

  mismatch_positions_alloc = (int *) MALLOC((querylength+8)*sizeof(int));
  mismatch_positions = &(mismatch_positions_alloc[4]);
  debug8(printf("Calling Genome_mismatches_left_trim with plusp %d, left %u, pos5 %d - %d, pos3 %d - %d\n",
		plusp,left,querylength,querylength,querylength,0));
  nmismatches = Genome_mismatches_left_trim(mismatch_positions,/*max_mismatches*/querylength,
					    /*ome*/omebits,/*ome_alt*/NULL,query_compress,
					    left,/*pos5*/0,/*pos3*/querylength,
					    plusp,genestrand);

  mismatch_positions[-4] = -999;
  mismatch_positions[-3] = -999;
  mismatch_positions[-2] = -999;

  mismatch_positions[-1] = -1;
  mismatch_positions[nmismatches] = querylength;

  mismatch_positions[nmismatches+1] = -999;
  mismatch_positions[nmismatches+2] = -999;
  mismatch_positions[nmismatches+3] = -999;

#ifdef DEBUG8
  printf("%d mismatches:",nmismatches);
  for (i = -1; i <= nmismatches; i++) {
    printf(" %d",mismatch_positions[i]);
  }
  printf("\n");

  printf("poly_p at");
  for (i = 0; i < querylength; i++) {
    if (poly_p[i] == true) {
      printf(" %d",i);
    }
  }
  printf("\n");
#endif


  _zero = _mm_setzero_si128();
  _threshold = _mm_set1_epi32(index1part+1); /* Because we measure nmatches+1 */

  /* TODO: With poly oligo, we do need to check for ptr >= end */
  /* Don't have to check for ptr >= end, because we are guaranteed to get a hit */
  ptr = &(mismatch_positions[0]);
  end = &(mismatch_positions[nmismatches]);
  /* last_value = -1; */
  cmp = 0x0000;
  while (ptr <= end && cmp == 0x0000) {
    while (ptr <= end && cmp == 0x0000) {
      _values = _mm_loadu_si128((__m128i *) ptr);
      _prev = _mm_loadu_si128((__m128i *) &(ptr[-1]));
      _nmatches = _mm_sub_epi32(_values,_prev); /* really nmatches + 1 */
      _mismatchp = _mm_cmplt_epi32(_nmatches,_threshold);
      cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mismatchp,_zero));
#ifdef DEBUG8
      printf("Values  : "); print_vector_dec(_values);
      printf("Prev    : "); print_vector_dec(_prev);
      printf("nmatches: "); print_vector_dec(_nmatches);
      printf("mismatchp: "); print_vector_dec(_mismatchp);
      printf("cmp %08X\n",cmp);
      printf("\n");
#endif
      /* last_value = ptr[3]; */
      ptr += 4;
    }

    debug8(printf("cmp %08X, trailing zeroes: %d\n",cmp,count_trailing_zeroes_32(cmp)));
    i = count_trailing_zeroes_32(cmp)/4;
    pos = ptr[-5 + i] + 1;
    if (pos < 0) {
      debug8(printf("REACHED UNEXPECTED POS %d\n",pos));
    } else if (plusp == true /*&& pos < querylength*/ && poly_p[pos] == true) {
      debug8(printf("1 POLY OLIGO AT %d DOES NOT COUNT\n",pos));
      cmp = 0x0000;
    } else if (plusp == false /*&& querylength - index1part - pos >= 0*/ &&
	       poly_p[querylength - index1part - pos] == true) {
      debug8(printf("2 POLY OLIGO AT %d DOES NOT COUNT\n",pos));
      cmp = 0x0000;
    } else {
      *trim5 = pos;
      *nmismatches5 = &(ptr[-4]) - &(mismatch_positions[0]) + i;
    }
  }

  if (ptr > end) {
    debug8(printf("REACHED END FOR TRIM5\n"));
    /* Re-do without looking at poly_p */
    ptr = &(mismatch_positions[0]);
    cmp = 0x0000;
    while (cmp == 0x0000) {
      _values = _mm_loadu_si128((__m128i *) ptr);
      _prev = _mm_loadu_si128((__m128i *) &(ptr[-1]));
      _nmatches = _mm_sub_epi32(_values,_prev); /* really nmatches + 1 */
      _mismatchp = _mm_cmplt_epi32(_nmatches,_threshold);
      cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mismatchp,_zero));
#ifdef DEBUG8
      printf("Values  : "); print_vector_dec(_values);
      printf("Prev    : "); print_vector_dec(_prev);
      printf("nmatches: "); print_vector_dec(_nmatches);
      printf("mismatchp: "); print_vector_dec(_mismatchp);
      printf("cmp %08X\n",cmp);
      printf("\n");
#endif
      /* last_value = ptr[3]; */
      ptr += 4;
    }

    debug8(printf("cmp %08X, trailing zeroes: %d\n",cmp,count_trailing_zeroes_32(cmp)));
    i = count_trailing_zeroes_32(cmp)/4;
    *trim5 = ptr[-5 + i] + 1;
    *nmismatches5 = &(ptr[-4]) - &(mismatch_positions[0]) + i;
  }

  debug8(printf("TRIM5 %d, NMISMATCHES5 %d\n\n",*trim5,*nmismatches5));


  /* Don't have to check ptr, because we are guaranteed to get a hit */
  ptr = &(mismatch_positions[nmismatches]);
  end = &(mismatch_positions[0]);
  /* last_value = querylength; */
  cmp = 0x0000;
  while (ptr >= end && cmp == 0x0000) {
    while (ptr >= end && cmp == 0x0000) {
      ptr -= 4;
      _next = _mm_loadu_si128((__m128i *) &(ptr[+1]));
      _values = _mm_loadu_si128((__m128i *) ptr);
      _nmatches = _mm_sub_epi32(_next,_values); /* really nmatches + 1 */
      _mismatchp = _mm_cmplt_epi32(_nmatches,_threshold);
      cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mismatchp,_zero));
#ifdef DEBUG8
      printf("Next    : "); print_vector_dec(_next);
      printf("Values  : "); print_vector_dec(_values);
      printf("nmatches: "); print_vector_dec(_nmatches);
      printf("mismatchp: "); print_vector_dec(_mismatchp);
      printf("cmp %08X\n",cmp);
      printf("\n");
#endif
      /* last_value = ptr[0]; */
    }

    debug8(printf("cmp %08X, leading zeroes - 16: %d\n",cmp,count_leading_zeroes_32(cmp) - 16));
    i = (count_leading_zeroes_32(cmp) - 16)/4;
    pos = ptr[4-i];
    if (pos < 0) {
      debug8(printf("REACHED UNEXPECTED POS %d\n",pos));
    } else if (plusp == true /*&& pos - index1part >= 0*/ && poly_p[pos - index1part] == true) {
      debug8(printf("3 POLY OLIGO AT %d DOES NOT COUNT\n",pos));
      cmp = 0x0000;
    } else if (plusp == false /*&& pos > 0*/ && poly_p[querylength - pos] == true) {
      debug8(printf("4 POLY OLIGO AT %d DOES NOT COUNT\n",pos));
      cmp = 0x0000;
    } else {
      *trim3 = querylength - pos;
      *nmismatches3 = &(mismatch_positions[nmismatches]) - &(ptr[4]) + i;
    }
  }

  if (ptr < end) {
    debug8(printf("REACHED END FOR TRIM3\n"));
    ptr = &(mismatch_positions[nmismatches]);
    cmp = 0x0000;
    while (cmp == 0x0000) {
      ptr -= 4;
      _next = _mm_loadu_si128((__m128i *) &(ptr[+1]));
      _values = _mm_loadu_si128((__m128i *) ptr);
      _nmatches = _mm_sub_epi32(_next,_values); /* really nmatches + 1 */
      _mismatchp = _mm_cmplt_epi32(_nmatches,_threshold);
      cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mismatchp,_zero));
#ifdef DEBUG8
      printf("Next    : "); print_vector_dec(_next);
      printf("Values  : "); print_vector_dec(_values);
      printf("nmatches: "); print_vector_dec(_nmatches);
      printf("mismatchp: "); print_vector_dec(_mismatchp);
      printf("cmp %08X\n",cmp);
      printf("\n");
#endif
      /* last_value = ptr[0]; */
    }

    debug8(printf("cmp %08X, leading zeroes - 16: %d\n",cmp,count_leading_zeroes_32(cmp) - 16));
    i = (count_leading_zeroes_32(cmp) - 16)/4;
    pos = ptr[4-i];
    *trim3 = querylength - pos;
    *nmismatches3 = &(mismatch_positions[nmismatches]) - &(ptr[4]) + i;
  }

  debug8(printf("TRIM3 %d, NMISMATCHES3 %d\n\n",*trim3,*nmismatches3));

  FREE(mismatch_positions_alloc);

  return nmismatches;
}
#endif


/* Overwrite values array */
/* Repetitive oligos can calse false indel diagonals */
static int
most_prevalent_uint (int *nloci, UINT4 *values, int nvalues) {
  int max_count, count;
  UINT4 *out, *ptr, *end, *first;

  assert(nvalues > 0);

  ptr = out = &(values[0]);	/* Reset */
  end = &(values[nvalues]);

  max_count = 1;
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u",*first));
    if (ptr + max_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[max_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[max_count - 1] != ptr[max_count - 2]) {
	/* Can jump immediately */
	ptr += max_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      ptr += max_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
      debug0(printf(" => Count %d, index %d => printing\n",ptr - first,first - values));
      if ((count = ptr - first) > max_count) {
	out = &(values[0]);	/* Reset */
	max_count = count;
      }
      *out++ = *first;
    }
  }

  *nloci = out - &(values[0]);
  return max_count;
}


#if 0
static int
most_prevalent_count (UINT4 *values, int nvalues) {
  int max_count, count;
  UINT4 *ptr, *end, *first;

  assert(nvalues > 0);

  ptr = values;
  end = &(values[nvalues]);

  max_count = 1;
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u",*first));
    if (ptr + max_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[max_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[max_count - 1] != ptr[max_count - 2]) {
	/* Can jump immediately */
	ptr += max_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      ptr += max_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
      debug0(printf(" => Count %d\n",ptr - first));
      if ((count = ptr - first) > max_count) {
	max_count = count;
      }
    }
  }

  return max_count;
}
#endif



#if 0
static int
most_prevalent_uint_simd (int *nloci, UINT4 *values, int nvalues) {
  int max_count, count;
  UINT4 *ptr, *end, last_value;
  int i;

  assert(nvalues > 0);

  ptr = &(values[0]);
  end = &(values[n]);

  max_count = 1;
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u",*first));

    if (ptr + max_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[max_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      _first = _mm_set1_epi32(*first);
      _values = _mm_loadu_si128((__m128i *) ptr);

      match = _mm_cmpeq_epi32(_values,_first);
      matchbits = _mm_movemask_ps(_mm_castsi128_ps(matc));

#if !defined(HAVE_SSE4_2)
      nmatches = __builtin_popcount(matchbits);
#elif defined(HAVE_POPCNT)
      nmatches = _popcnt32(matchbits);
#elif defined HAVE_MM_POPCNT
      nmatches = _mm_popcnt_u32(matchbits);
#else
      nmatches = __builtin_popcount(matchbits);
#endif
      count = nmatches;

      while (nmatches == 4) {
	ptr += 4;
	_values = _mm_loadu_si128((__m128i *) ptr);

	match = _mm_cmpeq_epi32(_values,_first);
	matchbits = _mm_movemask_ps(_mm_castsi128_ps(matc));

#if !defined(HAVE_SSE4_2)
	nmatches = __builtin_popcount(matchbits);
#elif defined(HAVE_POPCNT)
	nmatches = _popcnt32(matchbits);
#elif defined HAVE_MM_POPCNT
	nmatches = _mm_popcnt_u32(matchbits);
#else
	nmatches = __builtin_popcount(matchbits);
#endif
	count += nmatches;
      }
    }
  }

  *out++ = count;
  ptr += nmatches;
}

#endif


struct Path_T {
  int total_nmatches;

  Transcript_T initial_transcript;

  List_T transcripts;		/* Set during disambiguate_paths */
  List_T transcripts_other;	/* Set during disambiguate_paths */

  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;

  int nexons;
  Uintlist_T endpoints;		/* Chrpos_T */
  Intlist_T gaps;

  int trim_low;
  int trim_high;

  bool concordantp;
};

static void
Path_free (Path_T *old) {
  List_T p;
  Transcript_T transcript;

  /* Don't need to free initial_transcript, since it has been pushed onto transcripts or transcripts_other of a representative path */
  assert((*old)->initial_transcript == NULL);

  for (p = (*old)->transcripts; p != NULL; p = List_next(p)) {
    transcript = (Transcript_T) List_head(p);
    Transcript_free(&transcript);
  }
  List_free(&(*old)->transcripts);

  for (p = (*old)->transcripts_other; p != NULL; p = List_next(p)) {
    transcript = (Transcript_T) List_head(p);
    Transcript_free(&transcript);
  }
  List_free(&(*old)->transcripts_other);

  Uintlist_free(&(*old)->endpoints);
  Intlist_free(&(*old)->gaps);
  FREE(*old);

  return;
}


static Path_T
Path_new (int total_nmatches, int trnum, int trstart, int trend, int trim_low, int trim_high,
	  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
	  Uintlist_T endpoints, Intlist_T gaps, bool concordantp) {
  Path_T new = (Path_T) MALLOC(sizeof(*new));
  
  new->total_nmatches = total_nmatches;

  new->initial_transcript = Transcript_new(trnum,trstart,trend);

  new->transcripts = (List_T) NULL;
  new->transcripts_other = (List_T) NULL;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;

  new->nexons = Uintlist_length(endpoints)/2;
  new->endpoints = endpoints;
  new->gaps = gaps;

  new->trim_low = trim_low;
  new->trim_high = trim_high;

  new->concordantp = concordantp;

  return new;
}


static int
Path_cmp (const void *x, const void *y) {
  Path_T a = * (Path_T *) x;
  Path_T b = * (Path_T *) y;
  Chrpos_T a0, b0;
  
  if (a->total_nmatches < b->total_nmatches) {
    return -1;
  } else if (b->total_nmatches < a->total_nmatches) {
    return +1;

  } else if (a->chrnum < b->chrnum) {
    return -1;
  } else if (b->chrnum < a->chrnum) {
    return +1;

  } else {
    a0 = Uintlist_head(a->endpoints);
    b0 = Uintlist_head(b->endpoints);

    if (a0 < b0) {
      return -1;
    } else if (a0 > b0) {
      return +1;

    } else {
      a0 = Uintlist_last_value(a->endpoints);
      b0 = Uintlist_last_value(b->endpoints);
      if (a0 < b0) {
	return -1;
      } else if (a0 > b0) {
	return +1;
      } else {
	return 0;
      }
    }
  }
}


static bool
exon_overlap_fwd_p (Uintlist_T p, Uintlist_T q) {
  Chrpos_T low, high;
  int overlap, length1, length2;

  if (Uintlist_head(Uintlist_next(p)) < Uintlist_head(q)) {
    return false;
  } else if (Uintlist_head(Uintlist_next(q)) < Uintlist_head(p)) {
    return false;
  } else if (Uintlist_head(p) == Uintlist_head(q)) {
    return true;
  } else if (Uintlist_head(Uintlist_next(p)) == Uintlist_head(Uintlist_next(q))) {
    return true;
  } else {
    low = Uintlist_head(p);
    if (Uintlist_head(q) > low) {
      low = Uintlist_head(q);
    }

    high = Uintlist_head(Uintlist_next(p));
    if (Uintlist_head(Uintlist_next(q)) < high) {
      high = Uintlist_head(Uintlist_next(q));
    }
    debug1(printf("overlap: %u..%u\n",low,high));
    overlap = high - low - 1;

    low = Uintlist_head(p);
    high = Uintlist_head(Uintlist_next(p));
    length1 = high - low - 1;
    p = Uintlist_next(Uintlist_next(p));

    low = Uintlist_head(q);
    high = Uintlist_head(Uintlist_next(q));
    length2 = high - low - 1;
    q = Uintlist_next(Uintlist_next(q));

    if (length1 > length2) {
      if (overlap < 0.8 * (double) length1) {
	return false;
      }
    } else {
      if (overlap < 0.8 * (double) length2) {
	return false;
      }
    }

    return true;
  }
}

static bool
exon_overlap_rev_p (Uintlist_T p, Uintlist_T q) {
  Chrpos_T low, high;
  int overlap, length1, length2;

  if (Uintlist_head(p) < Uintlist_head(Uintlist_next(q))) {
    return false;
  } else if (Uintlist_head(q) < Uintlist_head(Uintlist_next(p))) {
    return false;
  } else if (Uintlist_head(Uintlist_next(p)) == Uintlist_head(Uintlist_next(q))) {
    return true;
  } else if (Uintlist_head(p) == Uintlist_head(q)) {
    return true;
  } else {
    low = Uintlist_head(Uintlist_next(p));
    if (Uintlist_head(Uintlist_next(q)) > low) {
      low = Uintlist_head(Uintlist_next(q));
    }

    high = Uintlist_head(p);
    if (Uintlist_head(q) < high) {
      high = Uintlist_head(q);
    }
    debug1(printf("overlap: %u..%u\n",low,high));
    overlap = high - low - 1;

    low = Uintlist_head(Uintlist_next(p));
    high = Uintlist_head(p);
    length1 = high - low - 1;
    p = Uintlist_next(Uintlist_next(p));

    low = Uintlist_head(Uintlist_next(q));
    high = Uintlist_head(q);
    length2 = high - low - 1;
    q = Uintlist_next(Uintlist_next(q));

    if (length1 > length2) {
      if (overlap < 0.8 * (double) length1) {
	return false;
      }
    } else {
      if (overlap < 0.8 * (double) length2) {
	return false;
      }
    }

    return true;
  }
}


static bool
Path_equivalent_p (Path_T a, Path_T b) {
  if (a->total_nmatches == b->total_nmatches) {
    return true;
  } else {
    return false;
  }
}


static bool
Path_overlap_p (Path_T a, Path_T b) {
  Uintlist_T p, q;
  Chrpos_T low, high;
  int overlap, length1, length2;

  if (a->chrnum != b->chrnum) {
    return false;

  } else if (Uintlist_last_value(a->endpoints) < Uintlist_head(b->endpoints)) {
    /* No overlap */
    return false;

  } else if (Uintlist_last_value(b->endpoints) < Uintlist_head(a->endpoints)) {
    /* No overlap */
    return false;

  } else {
    /* Compute amount of overlap */
    overlap = length1 = length2 = 0;
    p = a->endpoints;
    q = b->endpoints;

    while (p != NULL && q != NULL) {
      debug1(printf("%u..%u vs %u..%u\n",
		    Uintlist_head(p),Uintlist_head(Uintlist_next(p)),
		    Uintlist_head(q),Uintlist_head(Uintlist_next(q))));

      if (Uintlist_head(Uintlist_next(p)) < Uintlist_head(q)) {
	debug1(printf("p is entirely before q\n"));
	low = Uintlist_head(p);
	high = Uintlist_head(Uintlist_next(p));
	length1 += high - low - 1;
	p = Uintlist_next(Uintlist_next(p));

      } else if (Uintlist_head(Uintlist_next(q)) < Uintlist_head(p)) {
	debug1(printf("q is entirely before p\n"));
	low = Uintlist_head(q);
	high = Uintlist_head(Uintlist_next(q));
	length2 += high - low - 1;
	q = Uintlist_next(Uintlist_next(q));

      } else {
	low = Uintlist_head(p);
	if (Uintlist_head(q) > low) {
	  low = Uintlist_head(q);
	}

	high = Uintlist_head(Uintlist_next(p));
	if (Uintlist_head(Uintlist_next(q)) < high) {
	  high = Uintlist_head(Uintlist_next(q));
	}
	debug1(printf("overlap: %u..%u\n",low,high));
	overlap += high - low - 1;

	if (Uintlist_head(p) == Uintlist_head(q)) {
	  low = Uintlist_head(p);
	  high = Uintlist_head(Uintlist_next(p));
	  length1 += high - low - 1;
	  p = Uintlist_next(Uintlist_next(p));

	  low = Uintlist_head(q);
	  high = Uintlist_head(Uintlist_next(q));
	  length2 += high - low - 1;
	  q = Uintlist_next(Uintlist_next(q));

	} else if (Uintlist_head(p) < Uintlist_head(q)) {
	  low = Uintlist_head(p);
	  high = Uintlist_head(Uintlist_next(p));
	  length1 += high - low - 1;
	  p = Uintlist_next(Uintlist_next(p));

	} else {
	  low = Uintlist_head(q);
	  high = Uintlist_head(Uintlist_next(q));
	  length2 += high - low - 1;
	  q = Uintlist_next(Uintlist_next(q));
	}
      }
    }

    while (p != NULL) {
      low = Uintlist_head(p);
      high = Uintlist_head(Uintlist_next(p));
      length1 += high - low - 1;
      p = Uintlist_next(Uintlist_next(p));
    }
      
    while (q != NULL) {
      low = Uintlist_head(q);
      high = Uintlist_head(Uintlist_next(q));
      length2 += high - low - 1;
      q = Uintlist_next(Uintlist_next(q));
    }

    debug1(printf("Overlap is %d\n",overlap));
    if (length1 > length2) {
      if (overlap < 0.8 * (double) length1) {
	return false;
      }
    } else {
      if (overlap < 0.8 * (double) length2) {
	return false;
      }
    }

    /* Look for a common low point */
    p = a->endpoints;
    q = b->endpoints;
    while (p != NULL && q != NULL) {
      if (Uintlist_head(p) == Uintlist_head(q)) {
	return true;
      } else if (Uintlist_head(p) < Uintlist_head(q)) {
	p = Uintlist_next(Uintlist_next(p));
      } else {
	q = Uintlist_next(Uintlist_next(q));
      }
    }

    /* Look for a common high point */
    p = a->endpoints;
    q = b->endpoints;
    while (p != NULL && q != NULL) {
      if (Uintlist_head(Uintlist_next(p)) == Uintlist_head(Uintlist_next(q))) {
	return true;
      } else if (Uintlist_head(Uintlist_next(p)) < Uintlist_head(Uintlist_next(q))) {
	p = Uintlist_next(Uintlist_next(p));
      } else {
	q = Uintlist_next(Uintlist_next(q));
      }
    }
      
    return false;
  }
}


/* IDEA: Instead of trying to predict overlap, try absorbing paths
   until they can no longer be combined.  Implement a combine_paths
   procedure */

/* Each path is in ascending order: low to high coordinates.  First entry is trnum, then trstart, then chrnum. */
static List_T
disambiguate_paths (Path_T *paths, int npaths) {
  List_T all_paths = NULL;
  int i, j, k, l;
  Uintlist_T final_endpoints, pstart, pend, p, *ptrs;
  Intlist_T final_gaps, gstart, g;
  Path_T path, representative_path, temp;
  Univcoord_T overall_start, overall_end, initialexon_high, finalexon_low;
  int trim_low, trim_high;
  bool foundp;
  int nconcordant;

  qsort(paths,npaths,sizeof(Uintlist_T),Path_cmp);
  for (i = 0; i < npaths; i++) {
    debug2(printf("Path %d (%d matches): %s\n",i,paths[i]->total_nmatches,Uintlist_to_string(paths[i]->endpoints)));
  }

  ptrs = (Uintlist_T *) MALLOC(npaths*sizeof(Uintlist_T));

  i = 0;
  while (i < npaths) {
    debug2(printf("Path %d: %s\n",i,Uintlist_to_string(paths[i]->endpoints)));
    /* last_value = Uintlist_last_value(paths[i]); */
    j = i + 1;
    while (j < npaths && Path_equivalent_p(paths[j],paths[i]) == true &&
	   Path_overlap_p(paths[j],paths[i]) == true) {
      j++;
    }
    debug2(printf("Cluster ends at %d\n",j));
    
    if (j - i == 1) {
      /* Cluster size is 1 */
      representative_path = paths[i];
      if (representative_path->concordantp == true) {
	representative_path->transcripts = List_push(NULL,(void *) representative_path->initial_transcript);
	representative_path->initial_transcript = (Transcript_T) NULL;
	all_paths = List_push(all_paths,(void *) representative_path);
	debug2(printf("Putting initial transcript for %d onto %d\n",i,i));
      } else {
	/* Shortcut for freeing initial transcript */
	Transcript_free(&representative_path->initial_transcript);
	Path_free(&representative_path);
      }

    } else {
      nconcordant = 0;
      for (k = i; k < j; k++) {
	if (paths[k]->concordantp == true) {
	  nconcordant += 1;
	}
      }

      if (nconcordant == 0) {
	for (k = i; k < j; k++) {
	  /* Shortcut for freeing initial transcript */
	  Transcript_free(&paths[k]->initial_transcript);
	  Path_free(&(paths[k]));
	}

      } else {
	k = i;
	l = i + nconcordant;
	while (k < i + nconcordant && l < j) {
	  while (k < i + nconcordant && paths[k]->concordantp == true) {
	    k++;
	  }
	  while (l < j && paths[l]->concordantp == false) {
	    l++;
	  }
	  if (k < i + nconcordant && l < j) {
	    /* Swap */
	    temp = paths[k];
	    paths[k] = paths[l];
	    paths[l] = temp;
	    k++; l++;
	  }
	}
	
	representative_path = paths[i];
	representative_path->transcripts = List_push(NULL,(void *) representative_path->initial_transcript);
	representative_path->initial_transcript = (Transcript_T) NULL;
	debug2(printf("Putting initial transcript for %d onto %d\n",i,i));

	for (k = i + nconcordant; k < j; k++) {
	  path = paths[k];
	  representative_path->transcripts_other =
	    List_push(representative_path->transcripts_other,(void *) path->initial_transcript);
	  path->initial_transcript = (Transcript_T) NULL;
	  Path_free(&path);
	}

	if (nconcordant == 1) {
	  all_paths = List_push(all_paths,(void *) representative_path);

	} else {
	  /* Need to find first common exon across all clusters */
	  debug2(printf("Finding initialexon_high\n"));
	  foundp = false;
	  ptrs[i] = paths[i]->endpoints;
	  while (ptrs[i] != NULL && foundp == false) {
	    foundp = true;
	    for (k = i + 1; k < i + nconcordant; k++) {
	      ptrs[k] = paths[k]->endpoints;
	      while (ptrs[k] != NULL && exon_overlap_fwd_p(ptrs[i],ptrs[k]) == false) {
		ptrs[k] = Uintlist_next(Uintlist_next(ptrs[k]));
	      }
	      if (ptrs[k] == NULL) {
		foundp = false;
	      }
	    }
	    if (foundp == false) {
	      ptrs[i] = Uintlist_next(Uintlist_next(ptrs[i]));
	    }
	  }
      
	  debug2(printf("foundp is %d\n",foundp));
	  if (foundp == true) {
	    pstart = ptrs[i];
	    overall_start = Uintlist_head(ptrs[i]);
	    initialexon_high = Uintlist_head(Uintlist_next(ptrs[i]));
	    for (k = i + 1; k < i + nconcordant; k++) {
	      if (Uintlist_head(ptrs[k]) > overall_start) {
		overall_start = Uintlist_head(ptrs[k]);
	      }
	      if (Uintlist_head(Uintlist_next(ptrs[k])) < initialexon_high) {
		initialexon_high = Uintlist_head(Uintlist_next(ptrs[k]));
	      }
	    }
	    debug2(printf("overall_start is %u\n",overall_start));
	    debug2(printf("initialexon_high is %u\n",initialexon_high));
	  }

	  debug2(printf("Finding finalexon_low\n"));
	  /* Reverse all endpoints */
	  for (k = i; k < i + nconcordant; k++) {
	    paths[k]->endpoints = Uintlist_reverse(paths[k]->endpoints);
	  }

	  foundp = false;
	  ptrs[i] = paths[i]->endpoints;
	  while (ptrs[i] != NULL && foundp == false) {
	    foundp = true;
	    for (k = i + 1; k < i + nconcordant; k++) {
	      ptrs[k] = paths[k]->endpoints;
	      while (ptrs[k] != NULL && exon_overlap_rev_p(ptrs[i],ptrs[k]) == false) {
		ptrs[k] = Uintlist_next(Uintlist_next(ptrs[k]));
	      }
	      if (ptrs[k] == NULL) {
		foundp = false;
	      }
	    }
	    if (foundp == false) {
	      ptrs[i] = Uintlist_next(Uintlist_next(ptrs[i]));
	    }
	  }
      
	  debug2(printf("foundp is %d\n",foundp));
	  if (foundp == true) {
	    pend = Uintlist_next(ptrs[i]);
	    overall_end = Uintlist_head(ptrs[i]);
	    finalexon_low = Uintlist_head(Uintlist_next(ptrs[i]));
	    for (k = i + 1; k < i + nconcordant; k++) {
	      if (Uintlist_head(ptrs[k]) < overall_end) {
		overall_end = Uintlist_head(ptrs[k]);
	      }
	      if (Uintlist_head(Uintlist_next(ptrs[k])) > finalexon_low) {
		finalexon_low = Uintlist_head(Uintlist_next(ptrs[k]));
	      }
	    }
	    debug2(printf("finalexon_low is %u\n",finalexon_low));
	  }

	  debug2(printf("initialexon high %u, finalexon_low %u\n",initialexon_high,finalexon_low));
	  /* Restore all endpoints */
	  for (k = i; k < i + nconcordant; k++) {
	    paths[k]->endpoints = Uintlist_reverse(paths[k]->endpoints);
	  }


	  /* Handle cases where one transcript may extend into intron and
	     another splices.  Want largest elt in each paths[k] that is
	     less than initialexon_high, and smallest elt that is greater
	     than finalexon_low. */
#if 0
	  overall_start = 0;
	  overall_end = (Univcoord_T) -1;
	  for (k = i; k < i + nconcordant; k++) {
	    p = paths[k]->endpoints = Uintlist_reverse(paths[k]->endpoints);
	    /* Skip entire exons */
	    while (Uintlist_head(Uintlist_next(p)) > initialexon_high) {
	      p = Uintlist_next(Uintlist_next(p));
	    }
	    p = Uintlist_next(p);
	    if (Uintlist_head(p) > overall_start) {
	      overall_start = Uintlist_head(p);
	    }

	    p = paths[k]->endpoints = Uintlist_reverse(paths[k]->endpoints);
	    /* Skip entire exons */
	    while (Uintlist_head(Uintlist_next(p)) < finalexon_low) {
	      p = Uintlist_next(Uintlist_next(p));
	    }
	    p = Uintlist_next(p);
	    if (Uintlist_head(p) < overall_end) {
	      overall_end = Uintlist_head(p);
	    }
	  }
#endif

	  debug2(printf("OVERALL START %u, END %u\n",overall_start,overall_end));

	  /* Put overall_start and overall_end on the right exons */
	  trim_low = representative_path->trim_low;
	  trim_high = representative_path->trim_high;
	  debug2(printf("INITIAL TRIMS: low %d, high %d\n",trim_low,trim_high));

	  debug2(printf("\n"));

	  p = representative_path->endpoints;
	  g = representative_path->gaps;

	  /* Skip entire exons */
	  while (p != pstart) {
	    debug2(printf("Skipping exon %u..%u which is before initial exon\n",
			  Uintlist_head(p),Uintlist_head(Uintlist_next(p))));
	    trim_low += Uintlist_head(Uintlist_next(p)) - Uintlist_head(p);
	    debug2(printf("Adding to trim_low (whole exon): %d = %u - %u\n",
			  Uintlist_head(Uintlist_next(p)) - Uintlist_head(p),
			  Uintlist_head(Uintlist_next(p)),Uintlist_head(p)));
	    p = Uintlist_next(Uintlist_next(p));
	    g = Intlist_next(g);
	  }
	  gstart = g;

	  trim_low += overall_start - Uintlist_head(pstart);
	  debug2(printf("Adding to trim_low (partial exon): %d = %u - %u\n",
			overall_start - Uintlist_head(pstart),
			overall_start,Uintlist_head(pstart)));
	  Uintlist_head_set(pstart,overall_start);

	  while (p != pend) {
	    p = Uintlist_next(Uintlist_next(p));
	  }

	  trim_high += Uintlist_head(Uintlist_next(p)) - overall_end;
	  debug2(printf("Adding to trim_high (partial exon): %d = %u - %u\n",
			Uintlist_head(Uintlist_next(p)) - overall_end,
			Uintlist_head(Uintlist_next(p)),overall_end));
	  Uintlist_head_set(Uintlist_next(p),overall_end);

	  p = Uintlist_next(Uintlist_next(p));
	  while (p != NULL) {
	    trim_high += Uintlist_head(Uintlist_next(p)) - Uintlist_head(p);
	    debug2(printf("Adding to trim_high (whole exon): %d = %u - %u\n",
			  Uintlist_head(Uintlist_next(p)) - Uintlist_head(p),
			  Uintlist_head(Uintlist_next(p)),Uintlist_head(p)));
	    p = Uintlist_next(Uintlist_next(p));
	  }

	  representative_path->trim_low = trim_low;
	  representative_path->trim_high = trim_high;
	  debug2(printf("FINAL TRIMS: low %d, high %d\n",trim_low,trim_high));


	  for (k = i+1; k < i + nconcordant; k++) {
	    path = paths[k];
	    debug2(printf("Putting initial transcript for %d onto %d\n",k,i));
	    representative_path->transcripts =
	      List_push(representative_path->transcripts,(void *) paths[k]->initial_transcript);
	    paths[k]->initial_transcript = (Transcript_T) NULL;
	    Path_free(&path);
	  }

	  final_endpoints = (Uintlist_T) NULL;
	  final_gaps = (Intlist_T) NULL;
	  for (p = pstart, g = gstart; p != pend; p = Uintlist_next(Uintlist_next(p))) {
	    final_endpoints = Uintlist_push(final_endpoints,Uintlist_head(p));
	    final_endpoints = Uintlist_push(final_endpoints,Uintlist_head(Uintlist_next(p)));
	    if (Uintlist_next(Uintlist_next(p)) != NULL) {
	      final_gaps = Intlist_push(final_gaps,Intlist_head(g));
	      g = Intlist_next(g);
	    }
	  }
	  /* Handle pend */
	  final_endpoints = Uintlist_push(final_endpoints,Uintlist_head(p));
	  final_endpoints = Uintlist_push(final_endpoints,Uintlist_head(Uintlist_next(p)));
#if 0
	  if (Uintlist_next(Uintlist_next(p)) != NULL) {
	    final_gaps = Intlist_push(final_gaps,Intlist_head(g));
	    g = Intlist_next(g);
	  }
#endif

	  assert(Intlist_length(final_gaps) == Uintlist_length(final_endpoints)/2 - 1);

	  final_endpoints = Uintlist_reverse(final_endpoints);
	  final_gaps = Intlist_reverse(final_gaps);
	  debug2(printf("FINAL ENDPOINTS: %s.  FINAL GAPS: %s.  TRIM LOW: %d, TRIM HIGH: %d\n\n",
			Uintlist_to_string(final_endpoints),Intlist_to_string(final_gaps),trim_low,trim_high));
	  Uintlist_free(&(representative_path->endpoints));
	  Intlist_free(&(representative_path->gaps));
	  representative_path->endpoints = final_endpoints;
	  representative_path->gaps = final_gaps;

	  all_paths = List_push(all_paths,(void *) representative_path);
	}
      }
    }

    i = j;
  }

  FREE(ptrs);

#if 0
  printf("PERFORMING CHECK OF INITIAL TRANSCRIPTS\n");
  for (a = all_paths; a != NULL; a = List_next(a)) {
    path = (Path_T) List_head(a);
    if (path->initial_transcript != NULL) {
      printf("Still have an initial_transcript that was not put into transcripts or transcripts_other in a representative path\n");
      abort();
    }
    for (b = path->transcripts; b != NULL; b = List_next(b)) {
      printf("Transcript %p\n",List_head(b));
    }
    for (b = path->transcripts_other; b != NULL; b = List_next(b)) {
      printf("Transcript other %p\n",List_head(b));
    }
  }
#endif

  return all_paths;
}


static List_T
single_hits_gplus (int *found_score_overall, int *found_score_within_trims,
		   List_T hits, Path_T *gplus_paths, int ngplus,
		   int querylength, Compress_T query_compress_fwd,
		   int nmismatches_allowed, int sensedir,
		   Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  Stage3end_T hit;
  List_T all_paths, substrings, junctions, q, r;
  Uintlist_T p;
  Intlist_T g;
  Path_T path, representative_path;
  int querystart, queryend;
  Univcoord_T alignstart, alignend, prev_alignend, left;
  Chrpos_T splice_distance;
  int nmismatches_whole, nmismatches, nindels, ninserts;
  Substring_T substring;
  Junction_T junction;
  bool abortp;

  debug(printf("Entered single_hits_gplus with %d gplus paths\n",ngplus));

  if (ngplus == 0) {
    all_paths = (List_T) NULL;
  } else if (ngplus == 1) {
    representative_path = gplus_paths[0];
    debug2(printf("Putting initial transcript for 0 onto 0\n"));
    representative_path->transcripts = List_push(NULL,representative_path->initial_transcript);
    representative_path->initial_transcript = (Transcript_T) NULL;
    all_paths = List_push(NULL,(void *) representative_path);
  } else {
    all_paths = disambiguate_paths(gplus_paths,ngplus);
  }

  for (q = all_paths; q != NULL; q = List_next(q)) {
    abortp = false;
    nmismatches_whole = 0;
    substrings = junctions = (List_T) NULL;

    path = (Path_T) List_head(q);
    debug2(printf("path %p has %d exons.  Endpoints %s\n",
		  path,path->nexons,Uintlist_to_string(path->endpoints)));
    /* gplus: Do not reverse endpoints or gaps */
    querystart = path->trim_low; /* gplus */
    prev_alignend = 0;
    for (p = path->endpoints, g = path->gaps; p != NULL; p = Uintlist_next(Uintlist_next(p))) {
      alignstart = Uintlist_head(p) + path->chroffset;
      alignend = Uintlist_head(Uintlist_next(p)) + path->chroffset;
      ninserts = 0;
      if (prev_alignend > 0) {
	if ((nindels = Intlist_head(g)) > 0) {
	  debug2(printf("Setting value of deletionpos to be %u\n",alignstart - nindels));
	  junction = Junction_new_deletion(nindels,/*deletionpos*/alignstart - nindels); /* gplus */
	} else if (nindels < 0) {
	  junction = Junction_new_insertion(-nindels);
	  ninserts = -nindels;
	  debug2(printf("Revising ninserts to be %d\n",ninserts));
	} else {
	  splice_distance = alignstart - prev_alignend; /* gplus */
	  debug2(printf("splice_distance %u = %u - %u\n",splice_distance,alignstart,prev_alignend));
	  junction = Junction_new_splice(splice_distance,sensedir,/*donor_prob*/2.0,/*acceptor_prob*/2.0);
	}
	junctions = Listpool_push(junctions,listpool,(void *) junction);
	g = Intlist_next(g);
      }
      querystart += ninserts;
      queryend = querystart + (alignend - alignstart); /* gplus */

      /* gplus */
      left = alignstart - querystart;
      debug2(printf("tr gplus query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
		    querystart,queryend,alignstart,alignend,alignstart - path->chroffset,alignend - path->chroffset,
		    left - path->chroffset,alignend - path->chroffset,querylength,queryend));
      nmismatches_whole += nmismatches =
	Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_fwd,left,
					  /*pos5*/querystart,/*pos3*/queryend,/*plusp*/true,/*genestrand*/0);
      debug2(printf("nmismatches fwd from %d to %d: %d\n",querystart,queryend,nmismatches));

      if ((substring = Substring_new(nmismatches,left,querystart,queryend,querylength,
				     /*plusp*/true,/*genestrand*/0,query_compress_fwd,
				     path->chrnum,path->chroffset,path->chrhigh,path->chrlength,
				     /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				     /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				     sensedir)) == NULL) {
	debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
	abortp = true;
      } else {
	debug2(printf(" => Created substring for %d..%d\n",querystart,queryend));
	substrings = Listpool_push(substrings,listpool,(void *) substring);
      }
	
      querystart = queryend;
      prev_alignend = alignend;
    }

    if (abortp == true || nmismatches_whole > nmismatches_allowed) {
      debug2(printf("ABORTING STAGE3END:  abortp %d.  nmismatches_whole %d vs nmismatches_allowed %d\n\n",
		    abortp,nmismatches_whole,nmismatches_allowed));
      for (r = substrings; r != NULL; r = List_next(r)) {
	substring = (Substring_T) List_head(r);
	Substring_free(&substring);
      }
      /* List_free(&substrings); -- allocated by Listpool_push */

      for (r = junctions; r != NULL; r = List_next(r)) {
	junction = (Junction_T) List_head(r);
	Junction_free(&junction);
      }
      /* List_free(&junctions); -- allocated by Listpool_push */

    } else {
      substrings = List_reverse(substrings);
      junctions = List_reverse(junctions);
      hit = Stage3end_new_precomputed(&(*found_score_overall),&(*found_score_within_trims),
				      /*nmismatches_bothdiff*/nmismatches_whole,
				      substrings,junctions,path->transcripts,path->transcripts_other,
				      querylength,path->chrnum,path->chroffset,path->chrhigh,path->chrlength,
				      /*gplusp*/true,/*genestrand*/0,/*sensedir*/SENSE_NULL,
				      listpool,/*method*/TR,level);
      debug(printf("Created new transcript hit %p with %d transcripts\n",hit,List_length(path->transcripts)));
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      assert(path->initial_transcript == NULL);
      path->transcripts = (List_T) NULL; /* hit owns it now */
      path->transcripts_other = (List_T) NULL; /* hit owns it now */
    }
  }

  for (q = all_paths; q != NULL; q = List_next(q)) {
    path = (Path_T) List_head(q);
    Path_free(&path);
  }
  List_free(&all_paths);

  return hits;
}


static List_T
single_hits_gminus (int *found_score_overall, int *found_score_within_trims,
		    List_T hits, Path_T *gminus_paths, int ngminus,
		    int querylength, Compress_T query_compress_rev,
		    int nmismatches_allowed, int sensedir,
		    Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  Stage3end_T hit;
  List_T all_paths, substrings, junctions, q, r;
  Uintlist_T p;
  Intlist_T g;
  Path_T path, representative_path;
  int querystart, queryend;
  Univcoord_T alignstart, alignend, prev_alignend, left;
  Chrpos_T splice_distance;
  int nmismatches_whole, nmismatches, nindels, ninserts;
  Substring_T substring;
  Junction_T junction;
  bool abortp;

  debug(printf("Entered single_hits_gminus with %d gminus paths\n",ngminus));

  if (ngminus == 0) {
    all_paths = (List_T) NULL;
  } else if (ngminus == 1) {
    representative_path = gminus_paths[0];
    debug2(printf("Putting initial transcript for 0 onto 0\n"));
    representative_path->transcripts = List_push(NULL,representative_path->initial_transcript);
    representative_path->initial_transcript = (Transcript_T) NULL;
    all_paths = List_push(NULL,(void *) representative_path);
  } else {
    all_paths = disambiguate_paths(gminus_paths,ngminus);
  }

  for (q = all_paths; q != NULL; q = List_next(q)) {
    abortp = false;
    nmismatches_whole = 0;
    substrings = junctions = (List_T) NULL;

    path = (Path_T) List_head(q);
    debug2(printf("path %p has %d exons.  Endpoints %s\n",
		  path,path->nexons,Uintlist_to_string(path->endpoints)));
    /* gminus: Reverse endpoints and gaps */
    path->endpoints = Uintlist_reverse(path->endpoints);
    path->gaps = Intlist_reverse(path->gaps);
    querystart = path->trim_high; /* gminus */
    prev_alignend = 0;
    for (p = path->endpoints, g = path->gaps; p != NULL; p = Uintlist_next(Uintlist_next(p))) {
      alignstart = Uintlist_head(p) + path->chroffset;
      alignend = Uintlist_head(Uintlist_next(p)) + path->chroffset;
      ninserts = 0;
      if (prev_alignend > 0) {
	if ((nindels = Intlist_head(g)) > 0) {
	  debug2(printf("Setting value of deletionpos to be %u\n",alignstart));
	  junction = Junction_new_deletion(nindels,/*deletionpos*/alignstart);  /* gminus */
	} else if (nindels < 0) {
	  junction = Junction_new_insertion(-nindels);
	  ninserts = -nindels;
	  debug2(printf("Revising ninserts to be %d\n",ninserts));
	} else {
	  splice_distance = prev_alignend - alignstart; /* gminus */
	  debug2(printf("splice_distance %u = %u - %u\n",splice_distance,prev_alignend,alignstart));
	  junction = Junction_new_splice(splice_distance,sensedir,/*donor_prob*/2.0,/*acceptor_prob*/2.0);
	}
	junctions = Listpool_push(junctions,listpool,(void *) junction);
	g = Intlist_next(g);
      }
      querystart += ninserts;
      queryend = querystart + (alignstart - alignend); /* gminus */

      /* gminus */
      left = alignend - (querylength - queryend);
      debug2(printf("tr gminus query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
		    querylength - queryend,querylength - querystart,alignstart,alignend,alignstart - path->chroffset,alignend - path->chroffset,
		    left - path->chroffset,alignend - path->chroffset,querylength,queryend));
      nmismatches_whole += nmismatches =
	Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_rev,left,
					  /*pos5*/querylength - queryend,/*pos3*/querylength - querystart,
					  /*plusp*/false,/*genestrand*/0);
      debug2(printf("mismatches rev from %d to %d: %d\n",querylength - queryend,querylength - querystart,nmismatches));

      if ((substring = Substring_new(nmismatches,left,querystart,queryend,querylength,
				     /*plusp*/false,/*genestrand*/0,query_compress_rev,
				     path->chrnum,path->chroffset,path->chrhigh,path->chrlength,
				     /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				     /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				     sensedir)) == NULL) {
	debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
	abortp = true;
      } else {
	debug2(printf(" => Created substring for %d..%d\n",querystart,queryend));
	substrings = Listpool_push(substrings,listpool,(void *) substring);
      }

      querystart = queryend;
      prev_alignend = alignend;
    }

    if (abortp == true || nmismatches_whole > nmismatches_allowed) {
      debug2(printf("ABORTING STAGE3END:  abortp %d.  nmismatches_whole %d vs nmismatches_allowed %d\n\n",
		    abortp,nmismatches_whole,nmismatches_allowed));
      for (r = substrings; r != NULL; r = List_next(r)) {
	substring = (Substring_T) List_head(r);
	Substring_free(&substring);
      }
      /* List_free(&substrings); -- allocated by Listpool_push */

      for (r = junctions; r != NULL; r = List_next(r)) {
	junction = (Junction_T) List_head(r);
	Junction_free(&junction);
      }
      /* List_free(&junctions); -- allocated by Listpool_push */

    } else {
      substrings = List_reverse(substrings);
      junctions = List_reverse(junctions);
      hit = Stage3end_new_precomputed(&(*found_score_overall),&(*found_score_within_trims),
				      /*nmismatches_bothdiff*/nmismatches_whole,
				      substrings,junctions,path->transcripts,path->transcripts_other,
				      querylength,path->chrnum,path->chroffset,path->chrhigh,path->chrlength,
				      /*gplusp*/false,/*genestrand*/0,/*sensedir*/SENSE_NULL,
				      listpool,/*method*/TR,level);
      debug(printf("Created new transcript hit %p with %d transcripts\n",hit,List_length(path->transcripts)));
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      path->transcripts = (List_T) NULL; /* hit owns it now */
      path->transcripts_other = (List_T) NULL; /* hit owns it now */
    }
  }

  for (q = all_paths; q != NULL; q = List_next(q)) {
    path = (Path_T) List_head(q);
    Path_free(&path);
  }
  List_free(&all_paths);

  return hits;
}



#if 0
static int
Path_trnum_cmp (const void *x, const void *y) {
  Path_T a = * (Path_T *) x;
  Path_T b = * (Path_T *) y;
  Transcript_T transcript_a, transcript_b;
  int trnum_a, trnum_b;

  transcript_a = (Transcript_T) a->initial_transcript;
  transcript_b = (Transcript_T) b->initial_transcript;
  
  trnum_a = Transcript_num(transcript_a);
  trnum_b = Transcript_num(transcript_b);

  if (trnum_a < trnum_b) {
    return -1;
  } else if (trnum_b < trnum_a) {
    return +1;
  } else {
    return 0;
  }
}
#endif


#if 0
static void
filter_concordant_paths (Path_T *paths5, int npaths5, Path_T *paths3, int npaths3) {
  int i, j, k, l, a, b;
  Transcript_T transcript5, transcript3;
  int trnum;
  
  debug2(printf("Entered Kmer_filter_concordant_paths with %d paths5 and %d paths3\n",
		npaths5,npaths3));

  qsort(paths5,npaths5,sizeof(Path_T),Path_trnum_cmp);
  qsort(paths3,npaths3,sizeof(Path_T),Path_trnum_cmp);
  i = j = 0;
  while (i < npaths5 && j < npaths3) {
    transcript5 = paths5[i]->initial_transcript;
    transcript3 = paths3[j]->initial_transcript;

    if ((trnum = Transcript_num(transcript5)) < Transcript_num(transcript3)) {
      i++;
    } else if (Transcript_num(transcript3) < Transcript_num(transcript5)) {
      j++;
    } else {
      k = i + 1;
      while (k < npaths5 && Transcript_num(paths5[k]->initial_transcript) == trnum) {
	k++;
      }

      l = j + 1;
      while (l < npaths3 && Transcript_num(paths3[l]->initial_transcript) == trnum) {
	l++;
      }

      for (a = i; a < k; a++) {
	transcript5 = paths5[a]->initial_transcript;
	for (b = j; b < l; b++) {
	  transcript3 = paths3[b]->initial_transcript;
	  if (Transcript_concordant_p(transcript5,transcript3) == true) {
	    debug2(printf(" => concordant\n"));
	    paths5[a]->concordantp = true;
	    paths3[b]->concordantp = true;
	  }
	}
      }
      i = k;
      j = l;
    }
  }

  return;
}
#endif


#ifdef HAVE_AVX2
#define all_zero_p(diff) _mm256_testz_si256(diff,diff)
#elif defined(HAVE_SSE4_1)
#define all_zero_p(diff) _mm_testz_si128(diff,diff)
#elif defined(HAVE_SSE2)
#define all_zero_p(diff) _mm_movemask_ps(_mm_castsi128_ps(diff)) == 0
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


/* Expecting transcriptome to not be a large genome */
static void
search_transcriptome_complete (Path_T **tplus_gplus_paths, int *n_tplus_gplus,
			       Path_T **tplus_gminus_paths, int *n_tplus_gminus,
			       Path_T **tminus_gminus_paths, int *n_tminus_gminus,
			       Path_T **tminus_gplus_paths, int *n_tminus_gplus,
			       
			       Trcoord_T *tplus_positions_5, int n_tplus_positions_5, int tplus_diagterm_5,
			       Trcoord_T *tminus_positions_5, int n_tminus_positions_5, int tminus_diagterm_5,
			       Trcoord_T *tplus_positions_3, int n_tplus_positions_3, int tplus_diagterm_3, 
			       Trcoord_T *tminus_positions_3, int n_tminus_positions_3, int tminus_diagterm_3,
			       
			       char *queryuc_ptr, int querylength,
			       Trcoord_T **tplus_stream_array, int *tplus_streamsize_array, int *tplus_diagterm_array,
			       Trcoord_T **tminus_stream_array, int *tminus_streamsize_array, int *tminus_diagterm_array,
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       Univ_IIT_T transcript_iit, Transcriptome_T transcriptome, Genome_T transcriptomebits, 
			       int nmismatches_allowed, bool concordantp) {
  int max_tplus_count, max_tminus_count;

  Univcoord_T left, left5, left3;
  Trcoord_T *tplus_diagonals, *tminus_diagonals;
  int n_tplus_loci, n_tminus_loci, n_tplus_diagonals, n_tminus_diagonals,
    total_ndiagonals_tplus = 0, total_ndiagonals_tminus = 0, i;
  
  Trcoord_T *diagonals;
  int ndiagonals;
  int tplus_streami = 0, tminus_streami = 0;

  Reader_T reader;
  int querypos;
  Oligostate_T last_state = INIT;
  Oligospace_T forward = 0, revcomp = 0, forward_oligo, revcomp_oligo;

#ifdef HAVE_AVX2
  __m256i _diff, _exonbounds, _trstart;
  int matchbits;
  int *ptr;
#elif defined(HAVE_SSE2)
  __m128i _diff, _exonbounds, _trstart;
  int matchbits;
  int *ptr;
#endif

  int trnum;
  int transcript_genestrand;
  Univcoord_T troffset, trhigh;
  Chrpos_T trlength, chrnum;
  int nexons, exoni;
  int overall_adj, adj, adj5, adj3, nindels;

  Uintlist_T lefts, q;
  Intlist_T querystarts, queryends, adjustments, r, s, t;
  int trim5, trim3;
  int nmismatches5, nmismatches3;

  bool want_lowest_coordinate_p;
  int indel_pos_5, indel_pos_3;
  int total_nmismatches, total_nmismatches_5, total_nmismatches_3;
  int total_nmatches, total_nmatches_5, total_nmatches_3;
  
  int *exonbounds, transcript_diag, overall_trstart, overall_trend, trstart;
#ifdef DEBUG2
  int trend;
#endif
  Chrpos_T *exonstarts, genomicpos, chrlength;
  Univcoord_T alignstart, alignend, chroffset, chrhigh;
  int querystart, queryend, exonpos, align_residual, exon_residual, len;
  bool abortp;

  Uintlist_T endpoints;
  Intlist_T gaps;
  Path_T *gplus_paths, *gminus_paths;
  int ngplus, ngminus;


  debug2(printf("\n\n***Complete transcriptome\n"));


  if (n_tplus_positions_5 > 0) {
    tplus_stream_array[tplus_streami] = tplus_positions_5;
    tplus_streamsize_array[tplus_streami] = n_tplus_positions_5;
    tplus_diagterm_array[tplus_streami] = tplus_diagterm_5;
    tplus_streami++;
  }

  if (n_tplus_positions_3 > 0) {
    tplus_stream_array[tplus_streami] = tplus_positions_3;
    tplus_streamsize_array[tplus_streami] = n_tplus_positions_3;
    tplus_diagterm_array[tplus_streami] = tplus_diagterm_3;
    tplus_streami++;
  }

  if (n_tminus_positions_5 > 0) {
    tminus_stream_array[tminus_streami] = tminus_positions_5;
    tminus_streamsize_array[tminus_streami] = n_tminus_positions_5;
    tminus_diagterm_array[tminus_streami] = tminus_diagterm_5;
    tminus_streami++;
  }

  if (n_tminus_positions_3 > 0) {
    tminus_stream_array[tminus_streami] = tminus_positions_3;
    tminus_streamsize_array[tminus_streami] = n_tminus_positions_3;
    tminus_diagterm_array[tminus_streami] = tminus_diagterm_3;
    tminus_streami++;
  }


  reader = Reader_new(queryuc_ptr,/*querystart*/1,/*queryend*/querylength - 1);
  last_state = INIT;
  forward = revcomp = 0;

  while ((last_state = Oligo_next_5(last_state,&querypos,&forward,&revcomp,reader,/*genestrand*/0)) != DONE) {
    forward_oligo = forward & oligobase_mask;
    /* These are lefts */
    ndiagonals = Indexdb_ptr_with_diagterm(&diagonals,/*plus_indexdb*/indexdb_tr,forward_oligo,
					   /*diagterm for left*/-querypos);
    if (ndiagonals > 0) {
      tplus_stream_array[tplus_streami] = diagonals;
      tplus_streamsize_array[tplus_streami] = ndiagonals;
      tplus_diagterm_array[tplus_streami] = -querypos;
      total_ndiagonals_tplus += ndiagonals;
      tplus_streami++;
    }

    revcomp_oligo = (revcomp >> leftreadshift) & oligobase_mask;
    /* These are lefts */
    ndiagonals = Indexdb_ptr_with_diagterm(&diagonals,/*minus_indexdb*/indexdb_tr,revcomp_oligo,
					   /*diagterm for left*/-querylength + querypos + index1part_tr);
    if (ndiagonals > 0) {
      tminus_stream_array[tminus_streami] = diagonals;
      tminus_streamsize_array[tminus_streami] = ndiagonals;
      tminus_diagterm_array[tminus_streami] = -querylength + querypos + index1part_tr;
      total_ndiagonals_tminus += ndiagonals;
      tminus_streami++;
    }
  }
  Reader_free(&reader);

  tplus_diagonals = Merge_diagonals(&n_tplus_diagonals,tplus_stream_array,
				    tplus_streamsize_array,tplus_diagterm_array,/*nstreams*/tplus_streami);

  if (n_tplus_diagonals == 0) {
    max_tplus_count = 0;
    n_tplus_loci = 0;
  } else {
    max_tplus_count = most_prevalent_uint(&n_tplus_loci,tplus_diagonals,n_tplus_diagonals);
  }

#ifdef DEBUG2
  printf("max_tplus_count: %d\n",max_tplus_count);
  for (i = 0; i < n_tplus_loci; i++) {
    printf("PLUS_LOCUS %u\n",tplus_diagonals[i]);
  }
#endif

  tminus_diagonals = Merge_diagonals(&n_tminus_diagonals,tminus_stream_array,
				     tminus_streamsize_array,tminus_diagterm_array,/*nstreams*/tminus_streami);

  if (n_tminus_diagonals == 0) {
    max_tminus_count = 0;
    n_tminus_loci = 0;
  } else {
    max_tminus_count = most_prevalent_uint(&n_tminus_loci,tminus_diagonals,n_tminus_diagonals);
  }


#ifdef DEBUG2
  printf("max_tminus_count: %d\n",max_tminus_count);
  for (i = 0; i < n_tminus_loci; i++) {
    printf("MINUS_LOCUS %u\n",tminus_diagonals[i]);
  }
#endif

  ngplus = ngminus = 0;	/* ngplus and ngminus partition n_tplus_loci */
  if (max_tplus_count < (querylength - index1part_tr)/5 &&
      max_tplus_count < max_tminus_count - SUBOPT) {
    /* Not enough kmers match or less than minus */
    debug2(printf("Not enough tplus kmers match (%d)\n",max_tplus_count));
    gplus_paths = gminus_paths = (Path_T *) NULL;

  } else if (n_tplus_loci == 0) {
    gplus_paths = gminus_paths = (Path_T *) NULL;

  } else {
    gplus_paths = (Path_T *) MALLOC(n_tplus_loci*sizeof(Path_T));
    gminus_paths = (Path_T *) MALLOC(n_tplus_loci*sizeof(Path_T));
    for (i = 0; i < n_tplus_loci; i++) {

      left = (Univcoord_T) tplus_diagonals[i] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
#if 0
      total_nmismatches = trim_ends(&trim5,&trim3,&nmismatches5,&nmismatches3,poly_p,
				    query_compress_fwd,left,querylength,transcriptomebits,/*plusp*/true,index1part_tr);
#else
      trim5 = Genome_first_kmer_left(&nmismatches5,transcriptomebits,query_compress_fwd,left,/*pos5*/0,/*pos3*/querylength,
				     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
      queryend = Genome_first_kmer_right(&nmismatches3,transcriptomebits,query_compress_fwd,left,/*pos5*/0,/*pos3*/querylength,
					 /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
      trim3 = querylength - queryend;
      total_nmismatches = Genome_count_mismatches_substring(transcriptomebits,NULL,query_compress_fwd,
							    left,/*pos5*/trim5,/*pos3*/queryend,
							    /*plusp*/true,/*genestrand*/0);
#endif

      total_nmatches = querylength - total_nmismatches;
      debug2(printf("tplus diagonal %u => trimmed %d..%d, nmismatches: %d total, %d at 5', and %d at 3'.  total_nmatches: %d\n",
		    tplus_diagonals[i],trim5,querylength - trim3,total_nmismatches,nmismatches5,nmismatches3,
		    total_nmatches));

      /* Look for conditions where we need not look for an indel */
      if (total_nmismatches < 3) {
	trim5 = trim3 = 0;
      } else {
	if (trim5 < LONG_END) {
	  if (nmismatches5 == 1) {
	    trim5 = 0;
	  }
	} else {
	  if (nmismatches5 <= ALLOWED_END_MISMATCHES) {
	    trim5 = 0;
  }
	}

	if (trim3 < LONG_END) {
	  if (nmismatches3 == 1) {
	    trim3 = 0;
	  }
	} else {
	  if (nmismatches3 <= ALLOWED_END_MISMATCHES) {
	    trim3 = 0;
	  }
	}
      }

      /* Should replace with a rank/select structure */
      trnum = Univ_IIT_get_one(transcript_iit,left,left);
      chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      if (transcript_genestrand > 0) {
	want_lowest_coordinate_p = true;
      } else {
	want_lowest_coordinate_p = false;
      }
      debug2(printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand));

      overall_adj = adj5 = adj3 = 0;
      total_nmismatches_5 = total_nmismatches_3 = querylength;
      debug2(printf("trim5 %d, trim3 %d\n",trim5,trim3));
      if (trim5 > 0) {
	indel_pos_5 = Indel_solve_end_low(&adj5,&total_nmismatches_5,left,/*firstbound*/trim5,querylength,
					  query_compress_fwd,transcriptomebits,NULL,
					  /*max_end_insertions*/3,/*max_end_deletions*/3,
					  /*max_mismatches*/total_nmismatches-1,/*plusp*/true,/*genestrand*/0,
					  want_lowest_coordinate_p);
	total_nmatches_5 = querylength - total_nmismatches_5;
	/* Compensate for nmismatches from indel to be fair to a straight alignment */
	if (adj5 > 0) {
	  total_nmismatches_5 += adj5;
	} else {
	  total_nmismatches_5 -= adj5;
	}
	debug2(printf("Indel_solve_end_low returns indel_pos %d, adj5 %d, with total_nmismatches %d\n",
		      indel_pos_5,adj5,total_nmismatches_5));
      }
      if (trim3 > 0) {
	indel_pos_3 = Indel_solve_end_high(&adj3,&total_nmismatches_3,left,/*lasttbound*/querylength - trim3,querylength,
					   query_compress_fwd,transcriptomebits,NULL,
					   /*max_end_insertions*/3,/*max_end_deletions*/3,
					   /*max_mismatches*/total_nmismatches-1,/*plusp*/true,/*genestrand*/0,
					   want_lowest_coordinate_p);
	total_nmatches_3 = querylength - total_nmismatches_3;
	/* Compensate for nmismatches from indel to be fair to a straight alignment */
	if (adj3 > 0) {
	  total_nmismatches_3 += adj3;
	} else {
	  total_nmismatches_3 -= adj3;
	}
	debug2(printf("Indel_solve_end_high returns indel_pos %d, adj3 %d, with total_nmismatches %d\n",
		      indel_pos_3,adj3,total_nmismatches_3));
      }

      if (total_nmismatches_5 < total_nmismatches_3) {
	if (total_nmismatches_5 >= total_nmismatches) {
	  overall_adj = 0;
	} else {
	  total_nmatches = total_nmatches_5;
	  overall_adj = adj5;
	}
      } else {
	if (total_nmismatches_3 >= total_nmismatches) {
	  overall_adj = 0;
	} else {
	  total_nmatches = total_nmatches_3;
	  overall_adj = adj3;
	}
      }

      debug2(printf("overall_adj %d, adj5 %d, adj3 %d\n",overall_adj,adj5,adj3));
      abortp = false;
      if (overall_adj == 0) {
	trim5 = trim3 = 0;      /* Reset trims, since we cannot find an indel */
	lefts = Uintlist_push(NULL,left);
	querystarts = Intlist_push(NULL,trim5);
	queryends = Intlist_push(NULL,querylength - trim3);
	adjustments = (Intlist_T) NULL;
	debug2(printf("1 querystarts: %s\n",Intlist_to_string(querystarts)));
	debug2(printf("1 queryends: %s\n",Intlist_to_string(queryends)));

      } else if (overall_adj > 0) {
	/* Deletion */
	if (total_nmismatches_5 < total_nmismatches_3) {
	  if (total_nmismatches_5 > nmismatches_allowed) {
	    debug2(printf("ABORTING INDEL BECAUSE CANNOT EXTEND TO END OF TRANSCRIPT: %d mismatches\n",total_nmismatches_5));
	    abortp = true;
	  } else {
	    overall_adj = adj5;
	    adj3 = 0;
	    left5 = left-adj5;

	    trim5 = trim3 = 0;      /* Reset trims, since we have found the (single) indel */
	    lefts = Uintlist_push(Uintlist_push(NULL,left),left5);
	    querystarts = Intlist_push(Intlist_push(NULL,indel_pos_5),trim5);
	    queryends = Intlist_push(Intlist_push(NULL,querylength - trim3),indel_pos_5);
	    adjustments = Intlist_push(NULL,adj5);
	    debug2(printf("2 querystarts: %s\n",Intlist_to_string(querystarts)));
	    debug2(printf("2 queryends: %s\n",Intlist_to_string(queryends)));
	  }
	} else {
	  if (total_nmismatches_3 > nmismatches_allowed) {
	    debug2(printf("ABORTING INDEL BECAUSE CANNOT EXTEND TO END OF TRANSCRIPT: %d mismatches\n",total_nmismatches_3));
	    abortp = true;
	  } else {
	    overall_adj = adj3;
	    adj5 = 0;
	    left3 = left+adj3;

	    trim5 = trim3 = 0;      /* Reset trims, since we have found the (single) indel */
	    lefts = Uintlist_push(Uintlist_push(NULL,left3),left);
	    querystarts = Intlist_push(Intlist_push(NULL,indel_pos_3),trim5);
	    queryends = Intlist_push(Intlist_push(NULL,querylength - trim3),indel_pos_3);
	    adjustments = Intlist_push(NULL,adj3);
	    debug2(printf("3 querystarts: %s\n",Intlist_to_string(querystarts)));
	    debug2(printf("3 queryends: %s\n",Intlist_to_string(queryends)));
	  }
	}

      } else {
	/* Insertion */
	if (total_nmismatches_5 < total_nmismatches_3) {
	  if (total_nmismatches_5 > nmismatches_allowed) {
	    debug2(printf("ABORTING INDEL BECAUSE CANNOT EXTEND TO END OF TRANSCRIPT: %d mismatches\n",total_nmismatches_5));
	    abortp = true;
	  } else {
	    overall_adj = adj5;
	    adj3 = 0;
	    left5 = left-adj5;

	    trim5 = trim3 = 0;      /* Reset trims, since we have found the (single) indel */
	    lefts = Uintlist_push(Uintlist_push(NULL,left),left5);
	    querystarts = Intlist_push(Intlist_push(NULL,indel_pos_5 - adj5),trim5);
	    queryends = Intlist_push(Intlist_push(NULL,querylength - trim3),indel_pos_5);
	    adjustments = Intlist_push(NULL,adj5);
	    debug2(printf("4 querystarts: %s\n",Intlist_to_string(querystarts)));
	    debug2(printf("4 queryends: %s\n",Intlist_to_string(queryends)));
	  }
	} else {
	  if (total_nmismatches_3 > nmismatches_allowed) {
	    debug2(printf("ABORTING INDEL BECAUSE CANNOT EXTEND TO END OF TRANSCRIPT: %d mismatches\n",total_nmismatches_3));
	    abortp = true;
	  } else {
	    overall_adj = adj3;
	    adj5 = 0;
	    left3 = left+adj3;

	    trim5 = trim3 = 0;      /* Reset trims, since we have found the (single) indel */
	    lefts = Uintlist_push(Uintlist_push(NULL,left3),left);
	    querystarts = Intlist_push(Intlist_push(NULL,indel_pos_3 - adj3),trim5);
	    queryends = Intlist_push(Intlist_push(NULL,querylength - trim3),indel_pos_3);
	    adjustments = Intlist_push(NULL,adj3);
	    debug2(printf("5 querystarts: %s\n",Intlist_to_string(querystarts)));
	    debug2(printf("5 queryends: %s\n",Intlist_to_string(queryends)));
	  }
	}
      }

      if (abortp == false) {
	nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);
	Univ_IIT_interval_bounds_linear(&troffset,&trhigh,&trlength,transcript_iit,trnum);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

	if (transcript_genestrand > 0) {
	  /* Case 1: tplus, gplus.  transcript_plusp == true && transcript_genestrand > 0 */
	  debug2(printf("Case 1: tplus, gplus\n"));
	  /* sensedir = SENSE_FORWARD; */
	  /* plusp = true; */

	  endpoints = (Uintlist_T) NULL;
	  gaps = (Intlist_T) NULL;

	  for (q = lefts, r = querystarts, s = queryends, t = adjustments; q != NULL;
	       q = Uintlist_next(q), r = Intlist_next(r), s = Intlist_next(s)) {
	    /* Find initial exonpos */
	    querystart = Intlist_head(r); /* tplus: compute from querystart to queryend */
	    queryend = Intlist_head(s);
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = Uintlist_head(q) - troffset;
	    trstart = transcript_diag + querystart;	/* tplus */
#ifdef DEBUG2
	    trend = transcript_diag + queryend; /* tplus */
	    printf("transcript_diag %d, trstart %d, trend %d, left %u, overall_adj %d\n",
		   transcript_diag,trstart,trend,Uintlist_head(q),overall_adj);
#endif

	    if (q == lefts) {
	      exoni = 0;
#ifdef HAVE_AVX2
	      _trstart = _mm256_set1_epi32(trstart);
	      ptr = exonbounds;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	      while (all_zero_p(_diff)) {
		exoni += 8;
		ptr += 8;
		_exonbounds = _mm256_loadu_si256((__m256i *) ptr);
		_diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	      }
	      matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	      exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	      _trstart = _mm_set1_epi32(trstart);
	      ptr = exonbounds;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	      while (all_zero_p(_diff)) {
		exoni += 4;
		ptr += 4;
		_exonbounds = _mm_loadu_si128((__m128i *) ptr);
		_diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	      }
	      matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	      exoni += count_trailing_zeroes_32(matchbits);
#else
	      while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
		exoni++;
	      }
#endif
	    } else {
	      while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
		exoni++;
	      }
	    }

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] + exonpos - 1; /* gplus */
	      debug2(printf("genomicpos %u = %u + %u - 1\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		queryend = querystart + len; /* tplus */
	
		alignstart = chroffset + genomicpos;
		alignend = alignstart + len; /* gplus */
		debug2(left = alignstart - querystart); /* tplus == gminus */
		debug2(printf("tr complete case 1 query %d..%d, align %u..%u [%u..%u], left %u = %u - %d\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignstart - chroffset,querystart));
	
		/* tplus: push alignstart first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("2 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  querystart += len; /* tplus */
		  genomicpos = exonstarts[exoni] - 1; /* gplus */
		  debug2(printf("genomicpos %u = %u - 1\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    
	      if (Uintlist_next(q) != NULL) {
		nindels = adj = Intlist_head(t);
		debug2(printf("INDEL: last piece %d | nindels %d | next piece %d\n",len,nindels,exon_residual));
		if (nindels < 0) {
		  nindels = -nindels;
		}
		debug2(printf("exon residual %d vs nindels %d\n",exon_residual,nindels));
		if (len <= nindels || exon_residual <= nindels) {
		  debug2(printf("ABORTING BECAUSE INDEL OVERLAPS SPLICE SITE\n"));
		  abortp = true;
		} else {
		  gaps = Intlist_push(gaps,(int) adj);
		}
		t = Intlist_next(t);
	      }
	    }
	  }
	  
	  if (abortp == true) {
	    debug2(printf("ABORTING PATH\n"));
	    Uintlist_free(&endpoints);
	    Intlist_free(&gaps);
	  } else {
	    endpoints = Uintlist_reverse(endpoints); /* gplus */
	    gaps = Intlist_reverse(gaps);	     /* gplus */
	    debug2(printf("PLUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			  ngplus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	    overall_trstart = transcript_diag + Intlist_head(querystarts); /* tplus */
	    overall_trend = transcript_diag + queryend; /* tplus: queryend is the last one processed */
	    gplus_paths[ngplus++] = Path_new(total_nmatches,trnum,overall_trstart,overall_trend,
					     /*trim_low*/trim5,/*trim_high*/trim3,
					     chrnum,chroffset,chrhigh,chrlength,endpoints,gaps,concordantp);
	  }

	} else {
	  /* Case 2: tplus, gminus.  transcript_plusp == true && transcript_genestrand < 0 */
	  debug2(printf("Case 2: tplus, gminus\n"));
	  /* sensedir = SENSE_FORWARD; */
	  /* plusp = false; */
    
	  /* abortp = false; */
	  endpoints = (Uintlist_T) NULL;
	  gaps = (Intlist_T) NULL;
	  
	  for (q = lefts, r = querystarts, s = queryends, t = adjustments; q != NULL;
	       q = Uintlist_next(q), r = Intlist_next(r), s = Intlist_next(s)) {
	    /* Find initial exonpos */
	    querystart = Intlist_head(r); /* tplus: compute from querystart to queryend */
	    queryend = Intlist_head(s);
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = Uintlist_head(q) - troffset;
	    trstart = transcript_diag + querystart;	/* tplus */
#ifdef DEBUG2
	    trend = transcript_diag + queryend; /* tplus */
	    printf("transcript_diag %d, trstart %d, trend %d, left %u, overall_adj %d\n",
		   transcript_diag,trstart,trend,Uintlist_head(q),overall_adj);
#endif

	    if (q == lefts) {
	      exoni = 0;
#ifdef HAVE_AVX2
	      _trstart = _mm256_set1_epi32(trstart);
	      ptr = exonbounds;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	      while (all_zero_p(_diff)) {
		exoni += 8;
		ptr += 8;
		_exonbounds = _mm256_loadu_si256((__m256i *) ptr);
		_diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	      }
	      matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	      exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	      _trstart = _mm_set1_epi32(trstart);
	      ptr = exonbounds;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	      while (all_zero_p(_diff)) {
		exoni += 4;
		ptr += 4;
		_exonbounds = _mm_loadu_si128((__m128i *) ptr);
		_diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	      }
	      matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	      exoni += count_trailing_zeroes_32(matchbits);
#else
	      while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
		exoni++;
	      }
#endif
	    } else {
	      while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
		exoni++;
	      }
	    }

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] - exonpos; /* gminus */
	      debug2(printf("genomicpos %u = %u - %u\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		queryend = querystart + len; /* tplus */

		alignstart = chroffset + genomicpos;
		alignend = alignstart - len; /* gminus */
		debug2(left = alignend - (querylength - queryend)); /* tplus != gminus */
		debug2(printf("tr complete case 2 query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignend - chroffset,querylength,queryend));

		/* tplus: push alignstart first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("4 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  querystart += len; /* tplus */
		  genomicpos = exonstarts[exoni]; /* gminus */
		  debug2(printf("genomicpos %u = %u\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }

	      if (Uintlist_next(q) != NULL) {
		nindels = adj = Intlist_head(t);
		debug2(printf("INDEL: last piece %d | nindels %d | next piece %d\n",len,nindels,exon_residual));
		if (nindels < 0) {
		  nindels = -nindels;
		}
		debug2(printf("exon residual %d vs nindels %d\n",exon_residual,nindels));
		if (len <= nindels || exon_residual <= nindels) {
		  debug2(printf("ABORTING BECAUSE INDEL OVERLAPS SPLICE SITE\n"));
		  abortp = true;
		} else {
		  gaps = Intlist_push(gaps,(int) adj);
		}
		t = Intlist_next(t);
	      }
	    }
	  }
      
	  if (abortp == true) {
	    debug2(printf("ABORTING PATH\n"));
	    Uintlist_free(&endpoints);
	    Intlist_free(&gaps);
	  } else {
	    /* gminus: Do not reverse endpoints or gaps */
	    debug2(printf("MINUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			  ngminus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	    overall_trstart = transcript_diag + Intlist_head(querystarts); /* tplus */
	    overall_trend = transcript_diag + queryend; /* tplus: queryend is the last one processed */
	    gminus_paths[ngminus++] = Path_new(total_nmatches,trnum,overall_trstart,overall_trend,
					       /*trim_low*/trim3,/*trim_high*/trim5,
					       chrnum,chroffset,chrhigh,chrlength,endpoints,gaps,concordantp);
	  }

	}

	Intlist_free(&adjustments);
	Intlist_free(&queryends);
	Intlist_free(&querystarts);
	Uintlist_free(&lefts);
      }
    }

    debug2(printf("tplus loci %d = ngplus %d + ngminus %d\n",n_tplus_loci,ngplus,ngminus));
  }
#if !defined(LARGE_GENOMES) || defined(HAVE_AVX512) || defined(HAVE_AVX2)
  FREE_ALIGN(tplus_diagonals);
#else
  FREE(tplus_diagonals);
#endif
  *tplus_gplus_paths = gplus_paths;
  *n_tplus_gplus = ngplus;
  *tplus_gminus_paths = gminus_paths;
  *n_tplus_gminus = ngminus;


  ngplus = ngminus = 0;	/* ngplus and ngminus partition n_tminus_loci */
  if (max_tminus_count < (querylength - index1part_tr)/5 &&
      max_tminus_count < max_tplus_count - SUBOPT) {
    /* Not enough kmers match or less than plus */
    debug2(printf("Not enough tminus kmers match (%d)\n",max_tminus_count));
    gplus_paths = gminus_paths = (Path_T *) NULL;

  } else if (n_tminus_loci == 0) {
    gplus_paths = gminus_paths = (Path_T *) NULL;

  } else {
    gplus_paths = (Path_T *) MALLOC(n_tminus_loci*sizeof(Path_T));
    gminus_paths = (Path_T *) MALLOC(n_tminus_loci*sizeof(Path_T));
    for (i = 0; i < n_tminus_loci; i++) {

      left = (Univcoord_T) tminus_diagonals[i] /*- querylength*/; /* NEW FORMULA for queryend of querylength */
#if 0
      total_nmismatches = trim_ends(&trim5,&trim3,&nmismatches5,&nmismatches3,poly_p,
				    query_compress_rev,left,querylength,transcriptomebits,/*plusp*/false,index1part_tr);
#else
      trim5 = Genome_first_kmer_left(&nmismatches5,transcriptomebits,query_compress_rev,left,/*pos5*/0,/*pos3*/querylength,
				     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
      queryend = Genome_first_kmer_right(&nmismatches3,transcriptomebits,query_compress_rev,left,/*pos5*/0,/*pos3*/querylength,
					 /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);

      trim3 = querylength - queryend;
      total_nmismatches = Genome_count_mismatches_substring(transcriptomebits,NULL,query_compress_rev,
							    left,/*pos5*/trim5,/*pos3*/queryend,
							    /*plusp*/false,/*genestrand*/0);
#endif

      total_nmatches = querylength - total_nmismatches;
      debug2(printf("tminus diagonal %u => trimmed %d..%d, nmismatches: %d total, %d at 5', and %d at 3', total_nmatches: %d\n",
		    tminus_diagonals[i],trim5,querylength - trim3,total_nmismatches,nmismatches5,nmismatches3,
		    total_nmatches));

      /* Look for conditions where we need not look for an indel */
      if (total_nmismatches < 3) {
	trim5 = trim3 = 0;
      } else {
	if (trim5 < LONG_END) {
	  if (nmismatches5 == 1) {
	    trim5 = 0;
	  }
	} else {
	  if (nmismatches5 <= ALLOWED_END_MISMATCHES) {
	    trim5 = 0;
	  }
	}

	if (trim3 < LONG_END) {
	  if (nmismatches3 == 1) {
	    trim3 = 0;
	  }
	} else {
	  if (nmismatches3 <= ALLOWED_END_MISMATCHES) {
	    trim3 = 0;
	  }
	}
      }

      /* Should replace with a rank/select structure */
      trnum = Univ_IIT_get_one(transcript_iit,left,left);
      chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      if (transcript_genestrand > 0) {
	want_lowest_coordinate_p = true;
      } else {
	want_lowest_coordinate_p = false;
      }
      debug2(printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand));

      overall_adj = adj5 = adj3 = 0;
      total_nmismatches_5 = total_nmismatches_3 = querylength;
      debug2(printf("trim5 %d, trim3 %d\n",trim5,trim3));
      if (trim5 > 0) {
	indel_pos_5 = Indel_solve_end_low(&adj5,&total_nmismatches_5,left,
					  /*firstbound*/trim5,querylength,
					  query_compress_rev,transcriptomebits,NULL,
					  /*max_end_insertions*/3,/*max_end_deletions*/3,
					  /*max_mismatches*/total_nmismatches-1,/*plusp*/false,/*genestrand*/0,
					  want_lowest_coordinate_p);
	total_nmatches_5 = querylength - total_nmismatches_5;
	/* Compensate for nmismatches from indel to be fair to a straight alignment */
	if (adj5 > 0) {
	  total_nmismatches_5 += adj5;
	} else {
	  total_nmismatches_5 -= adj5;
	}
	debug2(printf("Indel_solve_end_low returns indel_pos %d, adj5 %d, with total_nmismatches %d\n",
		      indel_pos_5,adj5,total_nmismatches_5));
      }
      if (trim3 > 0) {
	indel_pos_3 = Indel_solve_end_high(&adj3,&total_nmismatches_3,left,
					   /*lasttbound*/querylength - trim3,querylength,
					   query_compress_rev,transcriptomebits,NULL,
					   /*max_end_insertions*/3,/*max_end_deletions*/3,
					   /*max_mismatches*/total_nmismatches-1,/*plusp*/false,/*genestrand*/0,
					   want_lowest_coordinate_p);
	total_nmatches_3 = querylength - total_nmismatches_3;
	/* Compensate for nmismatches from indel to be fair to a straight alignment */
	if (adj3 > 0) {
	  total_nmismatches_3 += adj3;
	} else {
	  total_nmismatches_3 -= adj3;
	}
	debug2(printf("Indel_solve_end_high returns indel_pos %d, adj3 %d, with total_nmismatches %d\n",
		      indel_pos_3,adj3,total_nmismatches_3));
      }

      if (total_nmismatches_5 < total_nmismatches_3) {
	if (total_nmismatches_5 >= total_nmismatches) {
	  overall_adj = 0;
	} else {
	  total_nmatches = total_nmatches_5;
	  overall_adj = adj5;
	}
      } else {
	if (total_nmismatches_3 >= total_nmismatches) {
	  overall_adj = 0;
	} else {
	  total_nmatches = total_nmatches_3;
	  overall_adj = adj3;
	}
      }

      debug2(printf("overall_adj %d, adj5 %d, adj3 %d\n",overall_adj,adj5,adj3));
      abortp = false;
      if (overall_adj == 0) {
	trim5 = trim3 = 0;      /* Reset trims, since we cannot find an indel */
	lefts = Uintlist_push(NULL,left);
	querystarts = Intlist_push(NULL,trim3);
	queryends = Intlist_push(NULL,querylength - trim5);
	adjustments = (Intlist_T) NULL;
	debug2(printf("6 querystarts: %s\n",Intlist_to_string(querystarts)));
	debug2(printf("6 queryends: %s\n",Intlist_to_string(queryends)));

      } else if (overall_adj > 0) {
	/* Deletion */
	if (total_nmismatches_5 < total_nmismatches_3) {
	  if (total_nmismatches_5 > nmismatches_allowed) {
	    debug2(printf("ABORTING INDEL BECAUSE CANNOT EXTEND TO END OF TRANSCRIPT: %d mismatches\n",total_nmismatches_5));
	    abortp = true;
	  } else {
	    overall_adj = adj5;
	    adj3 = 0;
	    left5 = left-adj5;

	    trim5 = trim3 = 0;      /* Reset trims, since we have found the (single) indel */
	    lefts = Uintlist_push(Uintlist_push(NULL,left),left5);
	    querystarts = Intlist_push(Intlist_push(NULL,trim3),querylength - indel_pos_5);
	    queryends = Intlist_push(Intlist_push(NULL,querylength - indel_pos_5),querylength - trim5);
	    adjustments = Intlist_push(NULL,adj5);
	    debug2(printf("7 querystarts: %s\n",Intlist_to_string(querystarts)));
	    debug2(printf("7 queryends: %s\n",Intlist_to_string(queryends)));
	  }
	} else {
	  if (total_nmismatches_3 > nmismatches_allowed) {
	    debug2(printf("ABORTING INDEL BECAUSE CANNOT EXTEND TO END OF TRANSCRIPT: %d mismatches\n",total_nmismatches_3));
	    abortp = true;
	  } else {
	    overall_adj = adj3;
	    adj5 = 0;
	    left3 = left+adj3;

	    trim5 = trim3 = 0;      /* Reset trims, since we have found the (single) indel */
	    lefts = Uintlist_push(Uintlist_push(NULL,left3),left);
	    querystarts = Intlist_push(Intlist_push(NULL,trim3),querylength - indel_pos_3);
	    queryends = Intlist_push(Intlist_push(NULL,querylength - indel_pos_3),querylength - trim5);
	    adjustments = Intlist_push(NULL,adj3);
	    debug2(printf("8 querystarts: %s\n",Intlist_to_string(querystarts)));
	    debug2(printf("8 queryends: %s\n",Intlist_to_string(queryends)));
	  }
	}
	
      } else {
	/* Insertion */
	if (total_nmismatches_5 < total_nmismatches_3) {
	  if (total_nmismatches_5 > nmismatches_allowed) {
	    debug2(printf("ABORTING INDEL BECAUSE CANNOT EXTEND TO END OF TRANSCRIPT: %d mismatches\n",total_nmismatches_5));
	    abortp = true;
	  } else {
	    overall_adj = adj5;
	    adj3 = 0;
	    left5 = left-adj5;

	    trim5 = trim3 = 0;      /* Reset trims, since we have found the (single) indel */
	    lefts = Uintlist_push(Uintlist_push(NULL,left),left5);
	    querystarts = Intlist_push(Intlist_push(NULL,trim3),querylength - indel_pos_5);
	    queryends = Intlist_push(Intlist_push(NULL,querylength - indel_pos_5 + adj5),querylength - trim5);
	    adjustments = Intlist_push(NULL,adj5);
	    debug2(printf("9 querystarts: %s\n",Intlist_to_string(querystarts)));
	    debug2(printf("9 queryends: %s\n",Intlist_to_string(queryends)));
	  }
	} else {
	  if (total_nmismatches_3 > nmismatches_allowed) {
	    debug2(printf("ABORTING INDEL BECAUSE CANNOT EXTEND TO END OF TRANSCRIPT: %d mismatches\n",total_nmismatches_3));
	    abortp = true;
	  } else {
	    overall_adj = adj3;
	    adj5 = 0;
	    left3 = left+adj3;

	    trim5 = trim3 = 0;      /* Reset trims, since we have found the (single) indel */
	    lefts = Uintlist_push(Uintlist_push(NULL,left3),left);
	    querystarts = Intlist_push(Intlist_push(NULL,trim3),querylength - indel_pos_3);
	    queryends = Intlist_push(Intlist_push(NULL,querylength - indel_pos_3 + adj3),querylength - trim5);
	    adjustments = Intlist_push(NULL,adj3);
	    debug2(printf("10 querystarts: %s\n",Intlist_to_string(querystarts)));
	    debug2(printf("10 queryends: %s\n",Intlist_to_string(queryends)));
	  }
	}
      }

      if (abortp == false) {
	nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);
	Univ_IIT_interval_bounds_linear(&troffset,&trhigh,&trlength,transcript_iit,trnum);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	debug2(printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand));

	if (transcript_genestrand > 0) {
	  /* Case 3: tminus, gplus.  transcript_plusp == false && transcript_genestrand > 0 */
	  debug2(printf("Case 3: tminus, gplus\n"));
	  /* sensedir = SENSE_ANTI; */
	  /* plusp = false; */

	  endpoints = (Uintlist_T) NULL;
	  gaps = (Intlist_T) NULL;

	  for (q = lefts, r = querystarts, s = queryends, t = adjustments; q != NULL;
	       q = Uintlist_next(q), r = Intlist_next(r), s = Intlist_next(s)) {
	    /* Find initial exonpos */
	    queryend = Intlist_head(s); /* tminus: compute from queryend to querystart */
	    querystart = Intlist_head(r);
	    transcript_diag = Uintlist_head(q) - troffset;
	    trstart = transcript_diag + (querylength - queryend); /* tminus */
#ifdef DEBUG2
	    trend = transcript_diag + (querylength - querystart); /* tminus */
	    printf("transcript_diag %d, trstart %d, trend %d, left %u, overall_adj %d\n",
		   transcript_diag,trstart,trend,Uintlist_head(q),overall_adj);
#endif	  

	    if (q == lefts) {
	      exoni = 0;
#ifdef HAVE_AVX2
	      _trstart = _mm256_set1_epi32(trstart);
	      ptr = exonbounds;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	      while (all_zero_p(_diff)) {
		exoni += 8;
		ptr += 8;
		_exonbounds = _mm256_loadu_si256((__m256i *) ptr);
		_diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	      }
	      matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	      exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	      _trstart = _mm_set1_epi32(trstart);
	      ptr = exonbounds;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	      while (all_zero_p(_diff)) {
		exoni += 4;
		ptr += 4;
		_exonbounds = _mm_loadu_si128((__m128i *) ptr);
		_diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	      }
	      matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	      exoni += count_trailing_zeroes_32(matchbits);
#else
	      while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
		exoni++;
	      }
#endif
	    } else {
	      while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
		exoni++;
	      }
	    }

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] + exonpos - 1; /* gplus */
	      debug2(printf("genomicpos %u = %u + %u - 1\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d,align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		querystart = queryend - len; /* tminus */
      
		alignend = chroffset + genomicpos;
		alignstart = alignend + len; /* gplus */
		debug2(left = alignend - (querylength - queryend)); /* tminus != gplus */
		debug2(printf("tr complete case 3 query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignend - chroffset,querylength,queryend));

		/* tminus: plus alignend first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("6 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  queryend -= len; /* tminus */
		  genomicpos = exonstarts[exoni] - 1; /* gplus */
		  debug2(printf("genomicpos %u = %u - 1\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }

	      if (Uintlist_next(q) != NULL) {
		nindels = adj = Intlist_head(t);
		debug2(printf("INDEL: last piece %d | nindels %d | next piece %d\n",len,nindels,exon_residual));
		if (nindels < 0) {
		  nindels = -nindels;
		}
		debug2(printf("exon residual %d vs nindels %d\n",exon_residual,nindels));
		if (len <= nindels || exon_residual <= nindels) {
		  debug2(printf("ABORTING BECAUSE INDEL OVERLAPS SPLICE SITE\n"));
		  abortp = true;
		} else {
		  gaps = Intlist_push(gaps,(int) adj);
		}
		t = Intlist_next(t);
	      }
	    }
	  }

	  if (abortp == true) {
	    debug2(printf("ABORTING PATH\n"));
	    Uintlist_free(&endpoints);
	    Intlist_free(&gaps);
	  } else {
	    endpoints = Uintlist_reverse(endpoints); /* gplus */
	    gaps = Intlist_reverse(gaps);	     /* gplus */
	    debug2(printf("MINUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			  ngminus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	    overall_trstart = transcript_diag + (querylength - Intlist_head(queryends)); /* tminus */
	    overall_trend = transcript_diag + (querylength - querystart); /* tminus: querystart is the last one processed */
	    gminus_paths[ngminus++] = Path_new(total_nmatches,trnum,overall_trend,overall_trstart,
					       /*trim_low*/trim3,/*trim_high*/trim5,chrnum,chroffset,chrhigh,chrlength,
					       endpoints,gaps,concordantp);
	  }

	} else {
	  /* Case 4: tminus, gminus.  transcript_plusp == false && transcript_genestrand < 0 */
	  debug2(printf("Case 4: tminus, gminus\n"));
	  /* sensedir = SENSE_ANTI; */
	  /* plusp = true; */

	  /* abortp = false; */
	  endpoints = (Uintlist_T) NULL;
	  gaps = (Intlist_T) NULL;

	  for (q = lefts, r = querystarts, s = queryends, t = adjustments; q != NULL;
	       q = Uintlist_next(q), r = Intlist_next(r), s = Intlist_next(s)) {
	    /* Find initial exonpos */
	    queryend = Intlist_head(s); /* tminus: compute from queryend to querystart */
	    querystart = Intlist_head(r);
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = Uintlist_head(q) - troffset;
	    trstart = transcript_diag + (querylength - queryend); /* tminus */
#ifdef DEBUG2
	    trend = transcript_diag + (querylength - querystart); /* tminus */
	    printf("transcript_diag %d, trstart %d, trend %d, left %u, overall_adj %d\n",
		   transcript_diag,trstart,trend,Uintlist_head(q),overall_adj);
#endif

	    if (q == lefts) {
	      exoni = 0;
#ifdef HAVE_AVX2
	      _trstart = _mm256_set1_epi32(trstart);
	      ptr = exonbounds;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	      while (all_zero_p(_diff)) {
		exoni += 8;
		ptr += 8;
		_exonbounds = _mm256_loadu_si256((__m256i *) ptr);
		_diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	      }
	      matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	      exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	      _trstart = _mm_set1_epi32(trstart);
	      ptr = exonbounds;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	      while (all_zero_p(_diff)) {
		exoni += 4;
		ptr += 4;
		_exonbounds = _mm_loadu_si128((__m128i *) ptr);
		_diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	      }
	      matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	      exoni += count_trailing_zeroes_32(matchbits);
#else
	      while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
		exoni++;
	      }
#endif
	    } else {
	      while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
		exoni++;
	      }
	    }

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] - exonpos; /* gminus */
	      debug2(printf("genomicpos %u = %u - %u\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		querystart = queryend - len; /* tminus */
	
		alignend = chroffset + genomicpos;
		alignstart = alignend - len; /* gminus */
		debug2(left = alignstart - querystart); /* tminus == gminus */
		debug2(printf("tr complete case 4 query %d..%d, align %u..%u [%u..%u], left %u = %u - %d\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignstart - chroffset,querystart));

		/* tminus: push alignend first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("8 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  queryend -= len; /* tminus */
		  genomicpos = exonstarts[exoni]; /* gminus */
		  debug2(printf("genomicpos %u = %u\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }

	      if (Uintlist_next(q) != NULL) {
		nindels = adj = Intlist_head(t);
		debug2(printf("INDEL: last piece %d | nindels %d | next piece %d\n",len,nindels,exon_residual));
		if (nindels < 0) {
		  nindels = -nindels;
		}
		debug2(printf("exon residual %d vs nindels %d\n",exon_residual,nindels));
		if (len <= nindels || exon_residual <= nindels) {
		  debug2(printf("ABORTING BECAUSE INDEL OVERLAPS SPLICE SITE\n"));
		  abortp = true;
		} else {
		  gaps = Intlist_push(gaps,(int) adj);
		}
		t = Intlist_next(t);
	      }
	    }
	  }

	  if (abortp == true) {
	    debug2(printf("ABORTING PATH\n"));
	    Uintlist_free(&endpoints);
	    Intlist_free(&gaps);
	  } else {
	    /* gminus: Do not reverse endpoints or gaps */
	    debug2(printf("PLUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			  ngplus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	    overall_trstart = transcript_diag + (querylength - Intlist_head(queryends)); /* tminus */
	    overall_trend = transcript_diag + (querylength - querystart); /* tminus: querystart is the last one processed */
	    gplus_paths[ngplus++] = Path_new(total_nmatches,trnum,overall_trend,overall_trstart,
					     /*trim_low*/trim3,/*trim_high*/trim5,chrnum,chroffset,chrhigh,chrlength,
					     endpoints,gaps,concordantp);
	  }

	}

	Intlist_free(&adjustments);
	Intlist_free(&queryends);
	Intlist_free(&querystarts);
	Uintlist_free(&lefts);
      }
    }

    debug2(printf("tminus loci %d = ngplus %d + ngminus %d\n",n_tminus_loci,ngplus,ngminus));
  }

#if !defined(LARGE_GENOMES) || defined(HAVE_AVX512) || defined(HAVE_AVX2)
  FREE_ALIGN(tminus_diagonals);
#else
  FREE(tminus_diagonals);
#endif
  *tminus_gminus_paths = gminus_paths;
  *n_tminus_gminus = ngminus;
  *tminus_gplus_paths = gplus_paths;
  *n_tminus_gplus = ngplus;

  return;
}



static int
search_transcriptome_ends (Path_T **tplus_gplus_paths, int *n_tplus_gplus,
			   Path_T **tplus_gminus_paths, int *n_tplus_gminus,
			   Path_T **tminus_gminus_paths, int *n_tminus_gminus,
			   Path_T **tminus_gplus_paths, int *n_tminus_gplus,
				
			   Trcoord_T **tplus_positions_5, int *n_tplus_positions_5, int *tplus_diagterm_5,
			   Trcoord_T **tminus_positions_5, int *n_tminus_positions_5, int *tminus_diagterm_5,
			   Trcoord_T **tplus_positions_3, int *n_tplus_positions_3, int *tplus_diagterm_3,
			   Trcoord_T **tminus_positions_3, int *n_tminus_positions_3, int *tminus_diagterm_3,
			   
			   char *queryuc_ptr, int querylength,
			   Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Univ_IIT_T transcript_iit, Transcriptome_T transcriptome, Genome_T transcriptomebits, 
			   int nmismatches_allowed, bool concordantp) {
  Path_T *gplus_paths, *gminus_paths;
  int ngplus, ngminus;
				
  Univcoord_T left, lefta, leftb, left0, left1;
  Trcoord_T *tplus_diagpairs, *tminus_diagpairs;
  int n_tplus_diagpairs, n_tminus_diagpairs, i, k;
  bool tplus_exactp, tminus_exactp;

  int adj, nindels;
  int trim5 = 0, trim3 = 0;	/* Assuming that alignment goes to the ends (found that way) */
  int trim_a5, trim_a3, trim_b5, trim_b3;
  int nmismatches, total_nmatches, nmismatches_a5, nmismatches_a3, nmismatches_b5, nmismatches_b3;
  int best_nmismatches_i, best_nmismatches_j;
  bool want_lowest_coordinate_p;

  Reader_T reader;
  int querypos, query_lastpos, indel_pos, ninserts;
  Oligospace_T forward, revcomp, forward_oligo, revcomp_oligo;

#ifdef HAVE_AVX2
  __m256i _diff, _exonbounds, _trstart;
  int matchbits;
  int *ptr;
#elif defined(HAVE_SSE2)
  __m128i _diff, _exonbounds, _trstart;
  int matchbits;
  int *ptr;
#endif

  int trnum;
  int transcript_genestrand;
  Univcoord_T troffset, trhigh;
  Chrpos_T trlength, chrnum;
  int nexons, exoni;

  int *exonbounds, transcript_diag, overall_trstart, overall_trend, trstart;
#ifdef DEBUG2
  int trend;
#endif
  Chrpos_T *exonstarts, genomicpos, chrlength;
  Univcoord_T alignstart, alignend, chroffset, chrhigh;
  int querystart, queryend, exonpos, align_residual, exon_residual, len;
  bool abortp;

  Uintlist_T endpoints;
  Intlist_T gaps;

  *tplus_positions_5 = *tminus_positions_5 = *tplus_positions_3 = *tminus_positions_3 = (Trcoord_T *) NULL;
  *n_tplus_positions_5 = *n_tminus_positions_5 = *n_tplus_positions_3 = *n_tminus_positions_3 = 0;
			

  debug3(printf("%s\n",queryuc_ptr));
  reader = Reader_new(queryuc_ptr,/*querystart*/0,/*queryend*/querylength);
  forward = revcomp = 0;
  if (Oligo_next_5(/*last_state*/INIT,&querypos,&forward,&revcomp,reader,/*genestrand*/0) != DONE) {
    forward_oligo = forward & oligobase_mask;
    *tplus_diagterm_5 = -querypos;
    *n_tplus_positions_5 = Indexdb_ptr_with_diagterm(&(*tplus_positions_5),/*plus_indexdb*/indexdb_tr,forward_oligo,
						   /*diagterm for left*/-querypos);

    revcomp_oligo = (revcomp >> leftreadshift) & oligobase_mask;
    *tminus_diagterm_5 = -querylength + querypos + index1part_tr;
    *n_tminus_positions_5 = Indexdb_ptr_with_diagterm(&(*tminus_positions_5),/*minus_indexdb*/indexdb_tr,revcomp_oligo,
						    /*diagterm for left*/-querylength + querypos + index1part_tr);
    debug3(printf("5' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  *n_tplus_positions_5,*n_tminus_positions_5));
  }
  Reader_free(&reader);

  query_lastpos = querylength - index1part_tr;
  reader = Reader_new(queryuc_ptr,/*querystart*/query_lastpos - /*index1interval*/1 + 1,/*queryend*/querylength);
  forward = revcomp = 0;
  if (Oligo_next_5(/*last_state*/INIT,&querypos,&forward,&revcomp,reader,/*genestrand*/0) != DONE) {
    forward_oligo = forward & oligobase_mask;
    *tplus_diagterm_3 = -querypos;
    *n_tplus_positions_3 = Indexdb_ptr_with_diagterm(&(*tplus_positions_3),/*plus_indexdb*/indexdb_tr,forward_oligo,
						     *tplus_diagterm_3);

    revcomp_oligo = (revcomp >> leftreadshift) & oligobase_mask;
    *tminus_diagterm_3 = -querylength + querypos + index1part_tr;
    *n_tminus_positions_3 = Indexdb_ptr_with_diagterm(&(*tminus_positions_3),/*minus_indexdb*/indexdb_tr,revcomp_oligo,
						      *tminus_diagterm_3);
    debug3(printf("3' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  *n_tplus_positions_5,*n_tminus_positions_5));
  }
  Reader_free(&reader);
	 
  tplus_diagpairs = Intersect_approx(&tplus_exactp,&n_tplus_diagpairs,
				     *tplus_positions_5,*n_tplus_positions_5,*tplus_diagterm_5,
				     *tplus_positions_3,*n_tplus_positions_3,*tplus_diagterm_3,
				     /*maxdistance*/3);
  tminus_diagpairs = Intersect_approx(&tminus_exactp,&n_tminus_diagpairs,
				      *tminus_positions_5,*n_tminus_positions_5,*tminus_diagterm_5,
				      *tminus_positions_3,*n_tminus_positions_3,*tminus_diagterm_3,
				      /*maxdistance*/3);
  debug2(printf("***Ultrafast transcriptome: exactp %d and %d.  %d plus and %d minus diagpairs\n",
		tplus_exactp,tminus_exactp,n_tplus_diagpairs,n_tminus_diagpairs));


  if (n_tplus_diagpairs == 0) {
    FREE(tplus_diagpairs);	/* Occupies memory even if n_tplus_diagpairs == 0 */
    *tplus_gplus_paths = *tplus_gminus_paths = (Path_T *) NULL;
    *n_tplus_gplus = *n_tplus_gminus = 0;

  } else if (tplus_exactp == true) {
    ngplus = ngminus = 0;
    gplus_paths = (Path_T *) CALLOC(n_tplus_diagpairs,sizeof(Path_T));
    gminus_paths = (Path_T *) CALLOC(n_tplus_diagpairs,sizeof(Path_T));

    for (i = 0, k = 0; i < n_tplus_diagpairs; i++, k += 2) {
      debug2(printf("tplus diagpairs: %u and %u\n",tplus_diagpairs[k],tplus_diagpairs[k+1]));
      if (tplus_diagpairs[k] == tplus_diagpairs[k+1]) {
	left = (Univcoord_T) tplus_diagpairs[k] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
	debug2a(printf("ULTRAFAST PLUS DIAGONAL %u\n",tplus_diagpairs[k]));
	if ((nmismatches = Genome_count_mismatches_substring(transcriptomebits,NULL,query_compress_fwd,left,
							     /*pos5*/0,/*pos3*/querylength,
							     /*plusp*/true,/*genestrand*/0)) <= nmismatches_allowed) {
	  debug2a(printf(" => nmismatches %d\n",nmismatches));
	  total_nmatches = querylength - nmismatches;

	  /* Should replace with a rank/select structure */
	  trnum = Univ_IIT_get_one(transcript_iit,left,left);
	  chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);

	  nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);
	  Univ_IIT_interval_bounds_linear(&troffset,&trhigh,&trlength,transcript_iit,trnum);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  debug2a(printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand));

	  if (transcript_genestrand > 0) {
	    /* Case 1 ultrafast nogap: tplus, gplus.  transcript_plusp == true && transcript_genestrand > 0 */
	    debug2(printf("Case 1 ultrafast nogap: tplus, gplus\n"));
	    /* sensedir = SENSE_FORWARD; */
	    /* plusp = true; */
	
	    abortp = false;
	    endpoints = (Uintlist_T) NULL;
	    gaps = (Intlist_T) NULL;

	    /* Find initial exonpos */
	    querystart = 0;	/* tplus: compute from querystart to queryend */
	    queryend = querylength;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    trstart = transcript_diag /*+ querystart (0)*/; /* tplus */
#ifdef DEBUG2
	    trend = transcript_diag + queryend; /* tplus */
	    printf("transcript_diag %d = left %u - troffset %u, trstart %d, trend %d\n",
		   transcript_diag,left,troffset,trstart,trend);
#endif

	    exoni = 0;
#ifdef HAVE_AVX2
	    _trstart = _mm256_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	    _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 8;
	      ptr += 8;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	    _trstart = _mm_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	    _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 4;
	      ptr += 4;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#else
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
#endif

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));
	  
	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] + exonpos - 1; /* gplus */
	      debug2(printf("genomicpos %u = %u + %u - 1\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		queryend = querystart + len; /* tplus */
	
		alignstart = chroffset + genomicpos;
		alignend = alignstart + len; /* gplus */
		debug2(left = alignstart /*- querystart (0)*/); /* gplus == gminus */
		debug2(printf("tr ultrafast nogap case 1 query %d..%d, align %u..%u [%u..%u], left %u = %u - %d\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignstart - chroffset,querystart));
	
		/* tplus: push alignstart first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("2 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  querystart += len; /* tplus */
		  genomicpos = exonstarts[exoni] - 1; /* gplus */
		  debug2(printf("genomicpos %u = %u - 1\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }
	  
	    if (abortp == true) {
	      debug2(printf("ABORTING PATH\n"));
	      Uintlist_free(&endpoints);
	      Intlist_free(&gaps);
	    } else {
	      endpoints = Uintlist_reverse(endpoints); /* gplus */
	      gaps = Intlist_reverse(gaps);	     /* gplus */
	      debug2(printf("PLUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			    ngplus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	      overall_trstart = transcript_diag	 /*+ querystart (0)*/; /* tplus */
	      overall_trend = transcript_diag + queryend; /* tplus */
	      gplus_paths[ngplus++] = Path_new(total_nmatches,trnum,overall_trstart,overall_trend,
					       /*trim_low*/0,/*trim_high*/0,chrnum,chroffset,chrhigh,chrlength,
					       endpoints,gaps,concordantp);
	    }

	  } else {
	    /* Case 2 ultrafast nogap: tplus, gminus.  transcript_plusp == true && transcript_genestrand < 0 */
	    debug2(printf("Case 2 ultrafast nogap: tplus, gminus\n"));
	    /* sensedir = SENSE_FORWARD; */
	    /* plusp = false; */
    
	    abortp = false;
	    endpoints = (Uintlist_T) NULL;
	    gaps = (Intlist_T) NULL;

	    /* Find initial exonpos */
	    querystart = 0;	/* tplus: compute from querystart to queryend */
	    queryend = querylength;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    trstart = transcript_diag /*+ querystart (0) */; /* tplus */
#ifdef DEBUG2
	    trend = transcript_diag + queryend; /* tplus */
	    printf("transcript_diag %d = left %u - troffset %u, trstart %d, trend %d\n",
		   transcript_diag,left,troffset,trstart,trend);
#endif

	    exoni = 0;
#ifdef HAVE_AVX2
	    _trstart = _mm256_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	    _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 8;
	      ptr += 8;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	    _trstart = _mm_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	    _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 4;
	      ptr += 4;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#else
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
#endif

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] - exonpos; /* gminus */
	      debug2(printf("genomicpos %u = %u - %u\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		queryend = querystart + len; /* tplus */

		alignstart = chroffset + genomicpos;
		alignend = alignstart - len; /* gminus */
		debug2(left = alignend /*- (querylength - queryend) (0)*/); /* tplus != gminus */
		debug2(printf("tr ultrafast nogap case 2 query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignend - chroffset,querylength,queryend));

		/* tplus: push alignstart first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("4 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  querystart += len; /* tplus */
		  genomicpos = exonstarts[exoni]; /* gminus */
		  debug2(printf("genomicpos %u = %u\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }
      
	    if (abortp == true) {
	      debug2(printf("ABORTING PATH\n"));
	      Uintlist_free(&endpoints);
	      Intlist_free(&gaps);
	    } else {
	      /* gminus: Do not reverse endpoints or gaps */
	      debug2(printf("MINUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			    ngminus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	      overall_trstart = transcript_diag	 /*+ querystart (0)*/; /* tplus */
	      overall_trend = transcript_diag + queryend; /* tplus */
	      gminus_paths[ngminus++] = Path_new(total_nmatches,trnum,overall_trstart,overall_trend,
						 /*trim_low*/0,/*trim_high*/0,chrnum,chroffset,chrhigh,chrlength,
						 endpoints,gaps,concordantp);
	    }
	  }
	}
      }
    }

    FREE(tplus_diagpairs);	/* Occupies memory even if n_tplus_diagpairs == 0 */
    *tplus_gplus_paths = gplus_paths;
    *n_tplus_gplus = ngplus;
    *tplus_gminus_paths = gminus_paths;
    *n_tplus_gminus = ngminus;

  } else if (tminus_exactp == true) {
    FREE(tplus_diagpairs);	/* Occupies memory even if n_tplus_diagpairs == 0 */
    *tplus_gplus_paths = *tplus_gminus_paths = (Path_T *) NULL;
    *n_tplus_gplus = *n_tplus_gminus = 0;

  } else {
    /* Handle approximate diagpairs below */
  }
    

  if (n_tminus_diagpairs == 0) {
    FREE(tminus_diagpairs);	/* Occupies memory even if n_tminus_diagpairs == 0 */
    *tminus_gplus_paths = *tminus_gminus_paths = (Path_T *) NULL;
    *n_tminus_gplus = *n_tminus_gminus = 0;

  } else if (tminus_exactp == true) {
    ngplus = ngminus = 0;
    gplus_paths = (Path_T *) CALLOC(n_tminus_diagpairs,sizeof(Path_T));
    gminus_paths = (Path_T *) CALLOC(n_tminus_diagpairs,sizeof(Path_T));

    for (i = 0, k = 0; i < n_tminus_diagpairs; i++, k += 2) {
      debug2(printf("tminus diagpairs: %u and %u\n",tminus_diagpairs[k],tminus_diagpairs[k+1]));
      if (tminus_diagpairs[k] == tminus_diagpairs[k+1]) {
	left = (Univcoord_T) tminus_diagpairs[k] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
	debug2a(printf("ULTRAFAST MINUS DIAGONAL %u\n",tminus_diagpairs[k]));
	if ((nmismatches = Genome_count_mismatches_substring(transcriptomebits,NULL,query_compress_rev,left,
							     /*pos5*/0,/*pos3*/querylength,
							     /*plusp*/false,/*genestrand*/0)) <= nmismatches_allowed) {
	  debug2a(printf(" => nmismatches %d\n",nmismatches));
	  total_nmatches = querylength - nmismatches;

	  /* Should replace with a rank/select structure */
	  trnum = Univ_IIT_get_one(transcript_iit,left,left);
	  Univ_IIT_interval_bounds_linear(&troffset,&trhigh,&trlength,transcript_iit,trnum);

	  /* Need chroffset to distinguish paralogs */
	  nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);
	  chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  debug2a(printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand));

	  if (transcript_genestrand > 0) {
	    /* Case 3 ultrafast nogap: tminus, gplus.  transcript_plusp == false && transcript_genestrand > 0 */
	    debug2(printf("Case 3 ultrafast nogap: tminus, gplus\n"));
	    /* sensedir = SENSE_ANTI; */
	    /* plusp = false; */

	    abortp = false;
	    endpoints = (Uintlist_T) NULL;
	    gaps = (Intlist_T) NULL;

	    /* Find initial exonpos */
	    queryend = querylength; /* tminus: compute from queryend to querystart */
	    querystart = 0;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    trstart = transcript_diag /*+ (querylength - queryend) (0)*/; /* tminus */
#ifdef DEBUG2
	    trend = transcript_diag + (querylength /*- querystart (0)*/); /* tminus */
	    printf("transcript_diag %d = left %u - troffset %u, trstart %d, trend %d\n",
		   transcript_diag,left,troffset,trstart,trend);
#endif	  

	    exoni = 0;
#ifdef HAVE_AVX2
	    _trstart = _mm256_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	    _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 8;
	      ptr += 8;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	    _trstart = _mm_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	    _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 4;
	      ptr += 4;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#else
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
#endif

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] + exonpos - 1; /* gplus */
	      debug2(printf("genomicpos %u = %u + %u - 1\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d,align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		querystart = queryend - len; /* tminus */
      
		alignend = chroffset + genomicpos;
		alignstart = alignend + len; /* gplus */
		debug2(left = alignend /*- (querylength - queryend) (0)*/);  /* tminus != gplus */
		debug2(printf("tr ultrafast nogap case 3 query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignend - chroffset,querylength,queryend));

		/* tminus: plus alignend first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("6 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  queryend -= len; /* tminus */
		  genomicpos = exonstarts[exoni] - 1; /* gplus */
		  debug2(printf("genomicpos %u = %u - 1\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    if (abortp == true) {
	      debug2(printf("ABORTING PATH\n"));
	      Uintlist_free(&endpoints);
	      Intlist_free(&gaps);
	    } else {
	      endpoints = Uintlist_reverse(endpoints); /* gplus */
	      gaps = Intlist_reverse(gaps);	     /* gplus */
	      debug2(printf("MINUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			    ngminus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	      overall_trstart = transcript_diag /*+ (querylength - queryend) (0)*/; /* tminus */
	      overall_trend = transcript_diag + (querylength /*- querystart (0)*/); /* tminus */
	      gminus_paths[ngminus++] = Path_new(total_nmatches,trnum,overall_trend,overall_trstart,
						 /*trim_low*/0,/*trim_high*/0,chrnum,chroffset,chrhigh,chrlength,
						 endpoints,gaps,concordantp);
	    }

	  } else {
	    /* Case 4 ultrafast nogap: tminus, gminus.  transcript_plusp == false && transcript_genestrand < 0 */
	    debug2(printf("Case 4 ultrafast nogap: tminus, gminus\n"));
	    /* sensedir = SENSE_ANTI; */
	    /* plusp = true; */

	    abortp = false;
	    endpoints = (Uintlist_T) NULL;
	    gaps = (Intlist_T) NULL;

	    /* Find initial exonpos */
	    queryend = querylength; /* tminus: compute from queryend to querystart */
	    querystart = 0;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    trstart = transcript_diag /*+ (querylength - queryend) (0)*/; /* tminus */
#ifdef DEBUG2
	    trend = transcript_diag + (querylength /*- querystart (0)*/); /* tminus */
	    printf("transcript_diag %d = left %u - troffset %u, trstart %d, trend %d\n",
		   transcript_diag,left,troffset,trstart,trend);
#endif	  

	    exoni = 0;
#ifdef HAVE_AVX2
	    _trstart = _mm256_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	    _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 8;
	      ptr += 8;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	    _trstart = _mm_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	    _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 4;
	      ptr += 4;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#else
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
#endif

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] - exonpos; /* gminus */
	      debug2(printf("genomicpos %u = %u - %u\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		querystart = queryend - len; /* tminus */
	
		alignend = chroffset + genomicpos;
		alignstart = alignend - len; /* gminus */
		debug2(left = alignstart /*- querystart (0)*/); /* tminus == gminus */
		debug2(printf("tr ultrafast nogap case 4 query %d..%d, align %u..%u [%u..%u], left %u = %u - %d\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignstart - chroffset,querystart));

		/* tminus: push alignend first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("8 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  queryend -= len; /* tminus */
		  genomicpos = exonstarts[exoni]; /* gminus */
		  debug2(printf("genomicpos %u = %u\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    if (abortp == true) {
	      debug2(printf("ABORTING PATH\n"));
	      Uintlist_free(&endpoints);
	      Intlist_free(&gaps);
	    } else {
	      /* gminus: Do not reverse endpoints or gaps */
	      debug2(printf("PLUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			    ngplus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	      overall_trstart = transcript_diag /*+ (querylength - queryend) (0)*/; /* tminus */
	      overall_trend = transcript_diag + (querylength /*- querystart (0)*/); /* tminus */
	      gplus_paths[ngplus++] = Path_new(total_nmatches,trnum,overall_trend,overall_trstart,
					       /*trim_low*/0,/*trim_high*/0,chrnum,chroffset,chrhigh,chrlength,
					       endpoints,gaps,concordantp);
	    }
	  }
	}
      }
    }
    FREE(tminus_diagpairs);	/* Occupies memory even if n_tminus_diagpairs == 0 */
    *tminus_gplus_paths = gplus_paths;
    *n_tminus_gplus = ngplus;
    *tminus_gminus_paths = gminus_paths;
    *n_tminus_gminus = ngminus;

  } else if (tplus_exactp == true) {
    FREE(tminus_diagpairs);	/* Occupies memory even if n_tminus_diagpairs == 0 */
    *tminus_gplus_paths = *tminus_gminus_paths = (Path_T *) NULL;
    *n_tminus_gplus = *n_tminus_gminus = 0;

  } else {
    /* Handle approximate diagpairs below */
  }

  /* Handle approximate diagpairs */
  if (tplus_exactp == false && tminus_exactp == false) {

    /* Case for n_tplus_diagpairs == 0 already handled above */
    if (n_tplus_diagpairs > 0) {
      ngplus = ngminus = 0;
      gplus_paths = (Path_T *) CALLOC(n_tplus_diagpairs,sizeof(Path_T));
      gminus_paths = (Path_T *) CALLOC(n_tplus_diagpairs,sizeof(Path_T));

      for (i = 0, k = 0; i < n_tplus_diagpairs; i++, k += 2) {
	debug2(printf("tplus diagpairs: %u and %u\n",tplus_diagpairs[k],tplus_diagpairs[k+1]));
	lefta = (Univcoord_T) tplus_diagpairs[k] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
	leftb = (Univcoord_T) tplus_diagpairs[k+1] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */

#if 0
	trim_ends(&trim_a5,&trim_a3,&nmismatches_a5,&nmismatches_a3,poly_p,
		  query_compress_fwd,lefta,querylength,transcriptomebits,/*plusp*/true,index1part_tr);
#else
	trim_a5 = Genome_first_kmer_left(&nmismatches_a5,transcriptomebits,query_compress_fwd,lefta,/*pos5*/0,/*pos3*/querylength,
					 /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
	trim_a3 = querylength - Genome_first_kmer_right(&nmismatches_a3,transcriptomebits,query_compress_fwd,lefta,/*pos5*/0,/*pos3*/querylength,
							/*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
#endif

#if 0
	trim_ends(&trim_b5,&trim_b3,&nmismatches_b5,&nmismatches_b3,poly_p,
		  query_compress_fwd,leftb,querylength,transcriptomebits,/*plusp*/true,15);
#else
	trim_b5 = Genome_first_kmer_left(&nmismatches_b5,transcriptomebits,query_compress_fwd,leftb,/*pos5*/0,/*pos3*/querylength,
					 /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
	trim_b3 = querylength - Genome_first_kmer_right(&nmismatches_b3,transcriptomebits,query_compress_fwd,leftb,/*pos5*/0,/*pos3*/querylength,
							/*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
#endif

	debug2a(printf("trimmed a %d..%d, trimmed b %d..%d\n",
		       trim_a5,querylength - trim_a3,trim_b5,querylength - trim_b3));
	if (trim_a5 == 0 || trim_b3 == querylength) {
	  left0 = lefta;
	  left1 = leftb;
	  adj = (int) (leftb - lefta);
	  debug2a(printf("Setting a first, b second, adj %d\n",adj));
	} else if (trim_a3 == querylength || trim_b5 == 0) {
	  left0 = leftb;
	  left1 = lefta;
	  adj = (int) (lefta - leftb);
	  debug2a(printf("Setting b first, a second, adj %d\n",adj));
	} else {
	  adj = 0;
	}

	trnum = Univ_IIT_get_one(transcript_iit,left0,left0);
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if (transcript_genestrand > 0) {
	  want_lowest_coordinate_p = true;
	} else {
	  want_lowest_coordinate_p = false;
	}
	debug2(printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand));

	if (adj > 0 &&
	    (indel_pos = Indel_resolve_middle_deletion(&best_nmismatches_i,&best_nmismatches_j,
						       /*left*/left0,/*indels*/-adj,
						       /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						       /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						       /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_fwd,
						       /*querystart*/0,/*queryend*/querylength,querylength,
						       nmismatches_allowed,/*plusp*/true,/*genestrand*/0,
						       want_lowest_coordinate_p)) > 0) {
	  ninserts = 0;
	  nindels = adj;
	  total_nmatches = querylength - best_nmismatches_i - best_nmismatches_j;
	} else if (adj < 0 &&
		   (indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
							       /*left*/left0,/*indels*/-adj,
							       /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							       /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							       /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_fwd,
							       /*querystart*/0,/*queryend*/querylength,querylength,
							       nmismatches_allowed,/*plusp*/true,/*genestrand*/0,
							       want_lowest_coordinate_p)) > 0) {
	  ninserts = nindels = -adj;
	  total_nmatches = querylength - best_nmismatches_i - best_nmismatches_j - ninserts;
	} else {
	  adj = 0;
	}

	if (adj != 0) {
	  nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);
	  Univ_IIT_interval_bounds_linear(&troffset,&trhigh,&trlength,transcript_iit,trnum);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

	  if (transcript_genestrand > 0) {
	    /* Case 1 ultrafast gap: tplus, gplus.  transcript_plusp == true && transcript_genestrand > 0 */
	    debug2(printf("Case 1 ultrafast gap: tplus, gplus\n"));
	    /* sensedir = SENSE_FORWARD; */
	    /* plusp = true; */

	    abortp = false;
	    endpoints = (Uintlist_T) NULL;
	    gaps = (Intlist_T) NULL;
	    
	    /* first part */
	    left = left0;
	    querystart = 0;
	    queryend = indel_pos;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    overall_trstart = trstart = transcript_diag + querystart;	/* tplus */
#ifdef DEBUG2
	    trend = transcript_diag + queryend; /* tplus */
	    printf("transcript_diag %d, trstart %d, trend %d, left %u, adj %d\n",
		   transcript_diag,trstart,trend, left,adj);
#endif

	    debug2(printf("Searching for exoni for trstart %d\n",trstart));
	    exoni = 0;
#ifdef HAVE_AVX2
	    _trstart = _mm256_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	    _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 8;
	      ptr += 8;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	    _trstart = _mm_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	    _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 4;
	      ptr += 4;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#else
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
#endif

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] + exonpos - 1; /* gplus */
	      debug2(printf("genomicpos %u = %u + %u - 1\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		queryend = querystart + len; /* tplus */
	
		alignstart = chroffset + genomicpos;
		alignend = alignstart + len; /* gplus */
		debug2(left = alignstart - querystart); /* tplus == gminus */
		debug2(printf("tr ultrafast gap case 1 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - %d\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignstart - chroffset,querystart));
	
		/* tplus: push alignstart first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("2 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  querystart += len; /* tplus */
		  genomicpos = exonstarts[exoni] - 1; /* gplus */
		  debug2(printf("genomicpos %u = %u - 1\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }
	    
	    /* gap */
	    debug2(printf("INDEL: last piece %d | nindels %d | next piece %d\n",len,nindels,exon_residual));
	    debug2(printf("exon residual %d vs nindels %d\n",exon_residual,nindels));
	    if (len <= nindels || exon_residual <= nindels) {
	      debug2(printf("ABORTING BECAUSE INDEL OVERLAPS SPLICE SITE\n"));
	      abortp = true;
	    } else {
	      gaps = Intlist_push(gaps,(int) adj);
	    }

	    /* second part */
	    left = left1;
	    querystart = indel_pos + ninserts;
	    queryend = querylength;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    trstart = transcript_diag + querystart;	/* tplus */
	    overall_trend = transcript_diag + queryend; /* tplus */
	    debug2(printf("transcript_diag %d, trstart %d, trend %d, left %u\n",
			  transcript_diag,trstart,overall_trend,left));

	    debug2(printf("Searching for exoni for trstart %d\n",trstart));
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] + exonpos - 1; /* gplus */
	      debug2(printf("genomicpos %u = %u + %u - 1\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		queryend = querystart + len; /* tplus */
	
		alignstart = chroffset + genomicpos;
		alignend = alignstart + len; /* gplus */
		debug2(left = alignstart - querystart); /* tplus == gminus */
		debug2(printf("tr ultrafast gap case 1 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - %d\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignstart - chroffset,querystart));
	
		/* tplus: push alignstart first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("2 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  querystart += len; /* tplus */
		  genomicpos = exonstarts[exoni] - 1; /* gplus */
		  debug2(printf("genomicpos %u = %u - 1\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    if (abortp == true) {
	      debug2(printf("ABORTING PATH\n"));
	      Uintlist_free(&endpoints);
	      Intlist_free(&gaps);
	    } else {
	      endpoints = Uintlist_reverse(endpoints); /* gplus */
	      gaps = Intlist_reverse(gaps);	     /* gplus */
	      debug2(printf("PLUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			    ngplus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	      gplus_paths[ngplus++] = Path_new(total_nmatches,trnum,overall_trstart,overall_trend,
					       /*trim_low*/trim5,/*trim_high*/trim3,
					       chrnum,chroffset,chrhigh,chrlength,endpoints,gaps,concordantp);
	    }

	  } else {
	    /* Case 2 ultrafast gap: tplus, gminus.  transcript_plusp == true && transcript_genestrand < 0 */
	    debug2(printf("Case 2 ultrafast gap: tplus, gminus\n"));
	    /* sensedir = SENSE_FORWARD; */
	    /* plusp = false; */
    
	    abortp = false;
	    endpoints = (Uintlist_T) NULL;
	    gaps = (Intlist_T) NULL;
	  
	    /* first part */
	    left = left0;
	    querystart = 0; /* tplus: compute from querystart to queryend */
	    queryend = indel_pos;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    overall_trstart = trstart = transcript_diag + querystart;	/* tplus */
#ifdef DEBUG2
            trend = transcript_diag + queryend; /* tplus */
	    printf("transcript_diag %d, trstart %d, trend %d, left %u, adj %d\n",
		   transcript_diag,trstart,trend,left,adj);
#endif

	    exoni = 0;
#ifdef HAVE_AVX2
	    _trstart = _mm256_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	    _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 8;
	      ptr += 8;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	    _trstart = _mm_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	    _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 4;
	      ptr += 4;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#else
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
#endif

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] - exonpos; /* gminus */
	      debug2(printf("genomicpos %u = %u - %u\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		queryend = querystart + len; /* tplus */

		alignstart = chroffset + genomicpos;
		alignend = alignstart - len; /* gminus */
		debug2(left = alignend - (querylength - queryend)); /* tplus != gminus */
		debug2(printf("tr ultrafast gap case 2 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignend - chroffset,querylength,queryend));

		/* tplus: push alignstart first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("4 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  querystart += len; /* tplus */
		  genomicpos = exonstarts[exoni]; /* gminus */
		  debug2(printf("genomicpos %u = %u\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    /* gap */
	    debug2(printf("INDEL: last piece %d | nindels %d | next piece %d\n",len,nindels,exon_residual));
	    debug2(printf("exon residual %d vs nindels %d\n",exon_residual,nindels));
	    if (len <= nindels || exon_residual <= nindels) {
	      debug2(printf("ABORTING BECAUSE INDEL OVERLAPS SPLICE SITE\n"));
	      abortp = true;
	    } else {
	      gaps = Intlist_push(gaps,(int) adj);
	    }

	    /* second part */
	    left = left1;
	    querystart = indel_pos + ninserts;
	    queryend = querylength;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    trstart = transcript_diag + querystart;	/* tplus */
	    overall_trend = transcript_diag + queryend; /* tplus */
	    debug2(printf("transcript_diag %d, trstart %d, trend %d, left %u\n",
			  transcript_diag,trstart,overall_trend,left));
	      
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] - exonpos; /* gminus */
	      debug2(printf("genomicpos %u = %u - %u\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		queryend = querystart + len; /* tplus */

		alignstart = chroffset + genomicpos;
		alignend = alignstart - len; /* gminus */
		debug2(left = alignend - (querylength - queryend)); /* tplus != gminus */
		debug2(printf("tr ultrfast gap case 2 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignend - chroffset,querylength,queryend));

		/* tplus: push alignstart first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("4 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  querystart += len; /* tplus */
		  genomicpos = exonstarts[exoni]; /* gminus */
		  debug2(printf("genomicpos %u = %u\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    if (abortp == true) {
	      debug2(printf("ABORTING PATH\n"));
	      Uintlist_free(&endpoints);
	      Intlist_free(&gaps);
	    } else {
	      /* gminus: Do not reverse endpoints or gaps */
	      debug2(printf("MINUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			    ngminus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	      gminus_paths[ngminus++] = Path_new(total_nmatches,trnum,overall_trstart,overall_trend,
						 /*trim_low*/trim3,/*trim_high*/trim5,
						 chrnum,chroffset,chrhigh,chrlength,endpoints,gaps,concordantp);
	    }
	  }
	}
      }

      FREE(tplus_diagpairs);
      *tplus_gplus_paths = gplus_paths;
      *n_tplus_gplus = ngplus;
      *tplus_gminus_paths = gminus_paths;
      *n_tplus_gminus = ngminus;
    }

    /* Case for n_tminus_diagpairs == 0 already handled above */
    if (n_tminus_diagpairs > 0) {
      ngplus = ngminus = 0;
      gplus_paths = (Path_T *) CALLOC(n_tminus_diagpairs,sizeof(Path_T));
      gminus_paths = (Path_T *) CALLOC(n_tminus_diagpairs,sizeof(Path_T));

      for (i = 0, k = 0; i < n_tminus_diagpairs; i++, k += 2) {
	debug2(printf("tminus diagpairs: %u and %u\n",tminus_diagpairs[k],tminus_diagpairs[k+1]));
	lefta = (Univcoord_T) tminus_diagpairs[k] /*- querylength*/; /* NEW FORMULA for queryend of querylength */
	leftb = (Univcoord_T) tminus_diagpairs[k+1] /*- querylength*/; /* NEW FORMULA for queryend of querylength */

#if 0
	trim_ends(&trim_a5,&trim_a3,&nmismatches_a5,&nmismatches_a3,poly_p,
		  query_compress_rev,lefta,querylength,transcriptomebits,/*plusp*/false,index1part_tr);
#else
	trim_a5 = Genome_first_kmer_left(&nmismatches_a5,transcriptomebits,query_compress_rev,lefta,/*pos5*/0,/*pos3*/querylength,
					 /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
	trim_a3 = querylength - Genome_first_kmer_right(&nmismatches_a3,transcriptomebits,query_compress_rev,lefta,/*pos5*/0,/*pos3*/querylength,
							/*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
#endif

#if 0
	trim_ends(&trim_b5,&trim_b3,&nmismatches_b5,&nmismatches_b3,poly_p,
		  query_compress_rev,leftb,querylength,transcriptomebits,/*plusp*/false,15);
#else
	trim_b5 = Genome_first_kmer_left(&nmismatches_b5,transcriptomebits,query_compress_rev,leftb,/*pos5*/0,/*pos3*/querylength,
					 /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
	trim_b3 = querylength - Genome_first_kmer_right(&nmismatches_b3,transcriptomebits,query_compress_rev,leftb,/*pos5*/0,/*pos3*/querylength,
							/*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
#endif

	debug2a(printf("trimmed a %d..%d, trimmed b %d..%d\n",
		       trim_a5,querylength - trim_a3,trim_b5,querylength - trim_b3));
	if (trim_a5 == 0 || trim_b3 == querylength) {
	  left0 = lefta;
	  left1 = leftb;
	  adj = (int) (leftb - lefta);
	  debug2a(printf("Setting a first, b second, adj %d\n",adj));
	} else if (trim_a3 == querylength || trim_b5 == 0) {
	  left0 = leftb;
	  left1 = lefta;
	  adj = (int) (lefta - leftb);
	  debug2a(printf("Setting b first, a second, adj %d\n",adj));
	} else {
	  adj = 0;
	}

	trnum = Univ_IIT_get_one(transcript_iit,left0,left0);
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if (transcript_genestrand > 0) {
	  want_lowest_coordinate_p = true;
	} else {
	  want_lowest_coordinate_p = false;
	}
	debug2(printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand));

	if (adj > 0 &&
	    (indel_pos = Indel_resolve_middle_deletion(&best_nmismatches_i,&best_nmismatches_j,
						       /*left*/left0,/*indels*/-adj,
						       /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						       /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						       /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_rev,
						       /*querystart*/0,/*queryend*/querylength,querylength,
						       nmismatches_allowed,/*plusp*/false,/*genestrand*/0,
						       want_lowest_coordinate_p)) > 0) {
	  ninserts = 0;
	  nindels = adj;
	  total_nmatches = querylength - best_nmismatches_i - best_nmismatches_j;
	} else if (adj < 0 &&
		   (indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
							       /*left*/left0,/*indels*/-adj,
							       /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							       /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							       /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_rev,
							       /*querystart*/0,/*queryend*/querylength,querylength,
							       nmismatches_allowed,/*plusp*/false,/*genestrand*/0,
							       want_lowest_coordinate_p)) > 0) {
	  ninserts = nindels = -adj;
	  total_nmatches = querylength - best_nmismatches_i - best_nmismatches_j - ninserts;
	} else {
	  adj = 0;
	}

	if (adj != 0) {
	  nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);
	  Univ_IIT_interval_bounds_linear(&troffset,&trhigh,&trlength,transcript_iit,trnum);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

	  if (transcript_genestrand > 0) {
	    /* Case 3 ultrafast gap: tminus, gplus.  transcript_plusp == false && transcript_genestrand > 0 */
	    debug2(printf("Case 3 ultrafast gap: tminus, gplus\n"));
	    /* sensedir = SENSE_ANTI; */
	    /* plusp = false; */

	    abortp = false;
	    endpoints = (Uintlist_T) NULL;
	    gaps = (Intlist_T) NULL;

	    /* first part */
	    left = left0;
	    queryend = querylength; /* tminus: compute from queryend to querystart */
	    querystart = querylength - indel_pos;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    overall_trstart = trstart = transcript_diag + (querylength - queryend); /* tminus */
#ifdef DEBUG2
	    trend = transcript_diag + (querylength - querystart); /* tminus */
	    printf("transcript_diag %d, trstart %d, trend %d, left %u, adj %d\n",
		   transcript_diag,trstart,trend,left,adj);
#endif
	  
	    exoni = 0;
#ifdef HAVE_AVX2
	    _trstart = _mm256_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	    _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 8;
	      ptr += 8;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	    _trstart = _mm_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	    _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 4;
	      ptr += 4;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#else
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
#endif

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] + exonpos - 1; /* gplus */
	      debug2(printf("genomicpos %u = %u + %u - 1\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d,align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		querystart = queryend - len; /* tminus */
      
		alignend = chroffset + genomicpos;
		alignstart = alignend + len; /* gplus */
		debug2(left = alignend - (querylength - queryend)); /* tminus != gplus */
		debug2(printf("tr ultrafast gap case 3 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignend - chroffset,querylength,queryend));

		/* tminus: plus alignend first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("6 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  queryend -= len; /* tminus */
		  genomicpos = exonstarts[exoni] - 1; /* gplus */
		  debug2(printf("genomicpos %u = %u - 1\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    /* gap */
	    debug2(printf("INDEL: last piece %d | nindels %d | next piece %d\n",len,nindels,exon_residual));
	    debug2(printf("exon residual %d vs nindels %d\n",exon_residual,nindels));
	    if (len <= nindels || exon_residual <= nindels) {
	      debug2(printf("ABORTING BECAUSE INDEL OVERLAPS SPLICE SITE\n"));
	      abortp = true;
	    } else {
	      gaps = Intlist_push(gaps,(int) adj);
	    }

	    /* second part */
	    left = left1;
	    queryend = querylength - indel_pos - ninserts; /* tminus: compute from queryend to querystart */
	    querystart = 0;
	    transcript_diag = left - troffset;
	    trstart = transcript_diag + (querylength - queryend); /* tminus */
	    overall_trend = transcript_diag + (querylength - querystart); /* tminus */
	    debug2(printf("transcript_diag %d, trstart %d, trend %d, left %u\n",
			  transcript_diag,trstart,overall_trend,left));
	  
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
	  
	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] + exonpos - 1; /* gplus */
	      debug2(printf("genomicpos %u = %u + %u - 1\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d,align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		querystart = queryend - len; /* tminus */
      
		alignend = chroffset + genomicpos;
		alignstart = alignend + len; /* gplus */
		debug2(left = alignend - (querylength - queryend)); /* tminus != gplus */
		debug2(printf("tr ultrafast gap case 3 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignend - chroffset,querylength,queryend));

		/* tminus: plus alignend first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("6 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  queryend -= len; /* tminus */
		  genomicpos = exonstarts[exoni] - 1; /* gplus */
		  debug2(printf("genomicpos %u = %u - 1\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    if (abortp == true) {
	      debug2(printf("ABORTING PATH\n"));
	      Uintlist_free(&endpoints);
	      Intlist_free(&gaps);
	    } else {
	      endpoints = Uintlist_reverse(endpoints); /* gplus */
	      gaps = Intlist_reverse(gaps);	     /* gplus */
	      debug2(printf("MINUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			    ngminus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	      gminus_paths[ngminus++] = Path_new(total_nmatches,trnum,overall_trend,overall_trstart,
						 /*trim_low*/trim3,/*trim_high*/trim5,chrnum,chroffset,chrhigh,chrlength,
						 endpoints,gaps,concordantp);
	    }

	  } else {
	    /* Case 4 ultrafast gap: tminus, gminus.  transcript_plusp == false && transcript_genestrand < 0 */
	    debug2(printf("Case 4 ultrafast gap: tminus, gminus\n"));
	    /* sensedir = SENSE_ANTI; */
	    /* plusp = true; */

	    abortp = false;
	    endpoints = (Uintlist_T) NULL;
	    gaps = (Intlist_T) NULL;

	    /* first part */
	    left = left0;
	    queryend = querylength; /* tminus: compute from queryend to querystart */
	    querystart = querylength - indel_pos;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    overall_trstart = trstart = transcript_diag + (querylength - queryend); /* tminus */
#ifdef DEBUG2
	    trend = transcript_diag + (querylength - querystart); /* tminus */
	    printf("transcript_diag %d, trstart %d, trend %d, left %u, adj %d\n",
			  transcript_diag,trstart,trend,left,adj);
#endif

	    exoni = 0;
#ifdef HAVE_AVX2
	    _trstart = _mm256_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	    _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 8;
	      ptr += 8;
	      _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	      _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	    _trstart = _mm_set1_epi32(trstart);
	    ptr = exonbounds;
	    _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	    _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    while (all_zero_p(_diff)) {
	      exoni += 4;
	      ptr += 4;
	      _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	      _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	    }
	    matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	    exoni += count_trailing_zeroes_32(matchbits);
#else
	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }
#endif

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] - exonpos; /* gminus */
	      debug2(printf("genomicpos %u = %u - %u\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		querystart = queryend - len; /* tminus */
	
		alignend = chroffset + genomicpos;
		alignstart = alignend - len; /* gminus */
		debug2(left = alignstart - querystart); /* tminus == gminus */
		debug2(printf("tr ultrafast gap case 4 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - %d\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignstart - chroffset,querystart));

		/* tminus: push alignend first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("8 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  queryend -= len; /* tminus */
		  genomicpos = exonstarts[exoni]; /* gminus */
		  debug2(printf("genomicpos %u = %u\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    /* gap */
	    debug2(printf("INDEL: last piece %d | nindels %d | next piece %d\n",len,nindels,exon_residual));
	    debug2(printf("exon residual %d vs nindels %d\n",exon_residual,nindels));
	    if (len <= nindels || exon_residual <= nindels) {
	      debug2(printf("ABORTING BECAUSE INDEL OVERLAPS SPLICE SITE\n"));
	      abortp = true;
	    } else {
	      gaps = Intlist_push(gaps,(int) adj);
	    }

	    /* second part */
	    left = left1;
	    queryend = querylength - indel_pos - ninserts; /* tminus: compute from queryend to querystart */
	    querystart = 0;
	    debug2(printf("\nAnalyzing query %d..%d\n",querystart,queryend));
	    transcript_diag = left - troffset;
	    trstart = transcript_diag + (querylength - queryend); /* tminus */
	    overall_trend = transcript_diag + (querylength - querystart); /* tminus */
	    debug2(printf("transcript_diag %d, trstart %d, trend %d, left %u\n",
			  transcript_diag,trstart,overall_trend,left));

	    while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	      exoni++;
	    }

	    if (exoni >= nexons) {
	      debug2(printf("ABORTING BECAUSE OF EXONI %d >= NEXONS %d\n",exoni,nexons));
	      abortp = true;
	    } else {
	      debug2(printf("exoni %d out of nexons %d\n",exoni,nexons));
	      if (exoni == 0) {
		exonpos = trstart;
	      } else {
		exonpos = trstart - exonbounds[exoni - 1];
	      }
	      exon_residual = exonbounds[exoni] - trstart;
	      debug2(printf("exoni %d => exon_residual %d, exonpos %d\n",exoni,exon_residual,exonpos));

	      /* Generate one substring per exon */
	      genomicpos = exonstarts[exoni] - exonpos; /* gminus */
	      debug2(printf("genomicpos %u = %u - %u\n",genomicpos,exonstarts[exoni],exonpos));
	      align_residual = queryend - querystart;
	      while (align_residual > 0) {
		debug2(printf("abortp %d, align_residual %d, exon_residual %d\n",abortp,align_residual,exon_residual));
		len = (align_residual <= exon_residual) ? align_residual : exon_residual;
		querystart = queryend - len; /* tminus */
	
		alignend = chroffset + genomicpos;
		alignstart = alignend - len; /* gminus */
		debug2(left = alignstart - querystart); /* tminus == gminus */
		debug2(printf("tr ultrafast gap case 4 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - %d\n",
			      querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
			      left - chroffset,alignstart - chroffset,querystart));

		/* tminus: push alignend first */
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignend - chroffset));
		endpoints = Uintlist_push(endpoints,(Chrpos_T) (alignstart - chroffset));

		if (align_residual <= exon_residual) {
		  exon_residual -= align_residual; /* Need to know this if we have a following deletion */
		  align_residual = 0;
		} else if (++exoni >= nexons) {
		  /* Aligning past the end of the transcript */
		  debug2(printf("8 ABORTING BECAUSE ALIGNING PAST END OF TRANSCRIPT\n"));
		  align_residual = 0;
		  abortp = true;
		} else {
		  /* nintrons += 1; */
		  align_residual -= len;
		  queryend -= len; /* tminus */
		  genomicpos = exonstarts[exoni]; /* gminus */
		  debug2(printf("genomicpos %u = %u\n",genomicpos,exonstarts[exoni]));
		  exon_residual = exonbounds[exoni] - exonbounds[exoni - 1];
		  debug2(printf("INTRON\n"));
		  gaps = Intlist_push(gaps,0);
		}
	      }
	    }

	    if (abortp == true) {
	      debug2(printf("ABORTING PATH\n"));
	      Uintlist_free(&endpoints);
	      Intlist_free(&gaps);
	    } else {
	      /* gminus: Do not reverse endpoints or gaps */
	      debug2(printf("PLUS PATH %d.  ENDPOINTS: %s.  GAPS: %s\n\n",
			    ngplus,Uintlist_to_string(endpoints),Intlist_to_string(gaps)));
	      gplus_paths[ngplus++] = Path_new(total_nmatches,trnum,overall_trend,overall_trstart,
					       /*trim_low*/trim3,/*trim_high*/trim5,chrnum,chroffset,chrhigh,chrlength,
					       endpoints,gaps,concordantp);
	    }
	  }
	}
      }

      FREE(tminus_diagpairs);
      *tminus_gplus_paths = gplus_paths;
      *n_tminus_gplus = ngplus;
      *tminus_gminus_paths = gminus_paths;
      *n_tminus_gminus = ngminus;
    }
  }

  debug2(printf("Counts: %d %d %d %d\n",
		*n_tplus_gplus,*n_tplus_gminus,*n_tminus_gplus,*n_tminus_gminus));

  return (*n_tplus_gplus) + (*n_tplus_gminus) + (*n_tminus_gplus) + (*n_tminus_gminus);
}


void
Kmer_search_transcriptome_single (int *found_score_overall, int *found_score_within_trims,
				  List_T *hits_gplus, List_T *hits_gminus,

				  char *queryuc_ptr, int querylength,
				  Trcoord_T **tplus_stream_array, int *tplus_streamsize_array, int *tplus_diagterm_array,
				  Trcoord_T **tminus_stream_array, int *tminus_streamsize_array, int *tminus_diagterm_array,

				  Compress_T query_compress_fwd, Compress_T query_compress_rev,
				  Univ_IIT_T transcript_iit, Transcriptome_T transcriptome, Genome_T transcriptomebits, 
				  int nmismatches_allowed,
				  Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  
  Path_T *tplus_gplus_paths, *tplus_gminus_paths, *tminus_gminus_paths, *tminus_gplus_paths;
  int n_tplus_gplus, n_tplus_gminus, n_tminus_gminus, n_tminus_gplus;

  Trcoord_T *tplus_positions_5, *tminus_positions_5, *tplus_positions_3, *tminus_positions_3;
  int n_tplus_positions_5, n_tminus_positions_5, n_tplus_positions_3, n_tminus_positions_3;
  int tplus_diagterm_5, tminus_diagterm_5, tplus_diagterm_3, tminus_diagterm_3;

  /* Set concordantp to be true, meaning that no further concordance is necessary */
  if (search_transcriptome_ends(&tplus_gplus_paths,&n_tplus_gplus,
				&tplus_gminus_paths,&n_tplus_gminus,
				&tminus_gminus_paths,&n_tminus_gminus,
				&tminus_gplus_paths,&n_tminus_gplus,
				
				&tplus_positions_5,&n_tplus_positions_5,&tplus_diagterm_5,
				&tminus_positions_5,&n_tminus_positions_5,&tminus_diagterm_5,
				&tplus_positions_3,&n_tplus_positions_3,&tplus_diagterm_3,
				&tminus_positions_3,&n_tminus_positions_3,&tminus_diagterm_3,

				queryuc_ptr,querylength,
				query_compress_fwd,query_compress_rev,
				transcript_iit,transcriptome,transcriptomebits, 
				nmismatches_allowed,/*concordantp*/true) > 0) {
  } else {
    /* Will reassign paths from search_transcriptome_complete */
    FREE(tplus_gplus_paths);
    FREE(tplus_gminus_paths);
    FREE(tminus_gplus_paths);
    FREE(tminus_gminus_paths);

    search_transcriptome_complete(&tplus_gplus_paths,&n_tplus_gplus,
				  &tplus_gminus_paths,&n_tplus_gminus,
				  &tminus_gminus_paths,&n_tminus_gminus,
				  &tminus_gplus_paths,&n_tminus_gplus,
				
				  tplus_positions_5,n_tplus_positions_5,tplus_diagterm_5,
				  tminus_positions_5,n_tminus_positions_5,tminus_diagterm_5,
				  tplus_positions_3,n_tplus_positions_3,tplus_diagterm_3,
				  tminus_positions_3,n_tminus_positions_3,tminus_diagterm_3,
				 
				  queryuc_ptr,querylength,
				  tplus_stream_array,tplus_streamsize_array,tplus_diagterm_array,
				  tminus_stream_array,tminus_streamsize_array,tminus_diagterm_array,

				  query_compress_fwd,query_compress_rev,
				  transcript_iit,transcriptome,transcriptomebits, 
				  nmismatches_allowed,/*concordantp*/true);
  }


  debug(printf("Converting tplus paths to hits\n"));
  *hits_gplus = single_hits_gplus(&(*found_score_overall),&(*found_score_within_trims),
				  *hits_gplus,tplus_gplus_paths,n_tplus_gplus,
				  querylength,query_compress_fwd,nmismatches_allowed,
				  /*sensedir*/SENSE_FORWARD,listpool,hitlistpool,level);
  *hits_gminus = single_hits_gminus(&(*found_score_overall),&(*found_score_within_trims),
				    *hits_gminus,tplus_gminus_paths,n_tplus_gminus,
				    querylength,query_compress_rev,nmismatches_allowed,
				    /*sensedir*/SENSE_FORWARD,listpool,hitlistpool,level);

  debug(printf("Converting tminus paths to hits\n"));
  *hits_gplus = single_hits_gplus(&(*found_score_overall),&(*found_score_within_trims),
				  *hits_gplus,tminus_gplus_paths,n_tminus_gplus,
				  querylength,query_compress_fwd,nmismatches_allowed,
				  /*sensedir*/SENSE_ANTI,listpool,hitlistpool,level);
  *hits_gminus = single_hits_gminus(&(*found_score_overall),&(*found_score_within_trims),
				    *hits_gminus,tminus_gminus_paths,n_tminus_gminus,
				    querylength,query_compress_rev,nmismatches_allowed,
				    /*sensedir*/SENSE_ANTI,listpool,hitlistpool,level);
  FREE(tplus_gplus_paths);
  FREE(tplus_gminus_paths);
  FREE(tminus_gminus_paths);
  FREE(tminus_gplus_paths);
  
  return;
}


#if 0
void
Kmer_search_transcriptome_paired (int *found_score_5, int *found_score_3,
				  List_T *hits5_gplus, List_T *hits5_gminus,
				  List_T *hits3_gplus, List_T *hits3_gminus,
				  
				  char *queryuc_ptr_5, int querylength5,
				  char *queryuc_ptr_3, int querylength3,
				  int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,

				  Trcoord_T **tplus_stream_array_5, int *tplus_streamsize_array_5, int *tplus_diagterm_array_5,
				  Trcoord_T **tminus_stream_array_5, int *tminus_streamsize_array_5, int *tminus_diagterm_array_5,
				  Trcoord_T **tplus_stream_array_3, int *tplus_streamsize_array_3, int *tplus_diagterm_array_3,
				  Trcoord_T **tminus_stream_array_3, int *tminus_streamsize_array_3, int *tminus_diagterm_array_3,

				  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
				  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

				  Univ_IIT_T transcript_iit, Transcriptome_T transcriptome, Genome_T transcriptomebits, 
				  int nmismatches_allowed_5, int nmismatches_allowed_3, 
				  int kmer_search_sizelimit, int level) {
  
  Path_T *tplus_gplus_paths_5, *tplus_gminus_paths_5, *tminus_gminus_paths_5, *tminus_gplus_paths_5,
    *tplus_gplus_paths_3, *tplus_gminus_paths_3, *tminus_gminus_paths_3, *tminus_gplus_paths_3;
  int n_tplus_gplus_5, n_tplus_gminus_5, n_tminus_gminus_5, n_tminus_gplus_5,
    n_tplus_gplus_3, n_tplus_gminus_3, n_tminus_gminus_3, n_tminus_gplus_3;

  /* 5 and 3 here mean the 5' and 3' ends of each read, not the 5' and 3' reads */
  Trcoord_T *tplus_positions_5, *tminus_positions_5, *tplus_positions_3, *tminus_positions_3;
  int n_tplus_positions_5, n_tminus_positions_5, n_tplus_positions_3, n_tminus_positions_3;
  int tplus_diagterm_5, tminus_diagterm_5, tplus_diagterm_3, tminus_diagterm_3;


  /* Solve for 5' read */
  /* Set concordantp to be false, meaning that we want to find concordance later */
  if (search_transcriptome_ends(&tplus_gplus_paths_5,&n_tplus_gplus_5,
				&tplus_gminus_paths_5,&n_tplus_gminus_5,
				&tminus_gminus_paths_5,&n_tminus_gminus_5,
				&tminus_gplus_paths_5,&n_tminus_gplus_5,
				
				&tplus_positions_5,&n_tplus_positions_5,&tplus_diagterm_5,
				&tminus_positions_5,&n_tminus_positions_5,&tminus_diagterm_5,
				&tplus_positions_3,&n_tplus_positions_3,&tplus_diagterm_3,
				&tminus_positions_3,&n_tminus_positions_3,&tminus_diagterm_3,
				
				queryuc_ptr_5,querylength5,
				query5_compress_fwd,query5_compress_rev,
				transcript_iit,transcriptome,transcriptomebits, 
				nmismatches_allowed_5,/*concordantp*/false) > 0) {
  } else {
    /* Will reassign paths from search_transcriptome_complete */
    FREE(tplus_gplus_paths_5);
    FREE(tplus_gminus_paths_5);
    FREE(tminus_gplus_paths_5);
    FREE(tminus_gminus_paths_5);

    search_transcriptome_complete(&tplus_gplus_paths_5,&n_tplus_gplus_5,
				  &tplus_gminus_paths_5,&n_tplus_gminus_5,
				  &tminus_gminus_paths_5,&n_tminus_gminus_5,
				  &tminus_gplus_paths_5,&n_tminus_gplus_5,
				
				  tplus_positions_5,n_tplus_positions_5,tplus_diagterm_5,
				  tminus_positions_5,n_tminus_positions_5,tminus_diagterm_5,
				  tplus_positions_3,n_tplus_positions_3,tplus_diagterm_3,
				  tminus_positions_3,n_tminus_positions_3,tminus_diagterm_3,
				 
				  queryuc_ptr_5,querylength5,
				  tplus_stream_array_5,tplus_streamsize_array_5,tplus_diagterm_array_5,
				  tminus_stream_array_5,tminus_streamsize_array_5,tminus_diagterm_array,
				  query5_compress_fwd,query5_compress_rev,
				  transcript_iit,transcriptome,transcriptomebits, 
				  nmismatches_allowed_5,/*concordantp*/false);
  }


  /* Solve for 3' read */
  /* Set concordantp to be false, meaning that we want to find concordance later */
  if (search_transcriptome_ends(&tplus_gplus_paths_3,&n_tplus_gplus_3,
				&tplus_gminus_paths_3,&n_tplus_gminus_3,
				&tminus_gminus_paths_3,&n_tminus_gminus_3,
				&tminus_gplus_paths_3,&n_tminus_gplus_3,
				
				&tplus_positions_5,&n_tplus_positions_5,&tplus_diagterm_5,
				&tminus_positions_5,&n_tminus_positions_5,&tminus_diagterm_5,
				&tplus_positions_3,&n_tplus_positions_3,&tplus_diagterm_3,
				&tminus_positions_3,&n_tminus_positions_3,&tminus_diagterm_3,
				
				queryuc_ptr_3,querylength3,
				query3_compress_fwd,query3_compress_rev,
				transcript_iit,transcriptome,transcriptomebits, 
				nmismatches_allowed_3,/*concordantp*/false) > 0) {
  } else {
    /* Will reassign paths from search_transcriptome_complete */
    FREE(tplus_gplus_paths_3);
    FREE(tplus_gminus_paths_3);
    FREE(tminus_gplus_paths_3);
    FREE(tminus_gminus_paths_3);

    search_transcriptome_complete(&tplus_gplus_paths_3,&n_tplus_gplus_3,
				  &tplus_gminus_paths_3,&n_tplus_gminus_3,
				  &tminus_gminus_paths_3,&n_tminus_gminus_3,
				  &tminus_gplus_paths_3,&n_tminus_gplus_3,
				
				  tplus_positions_5,n_tplus_positions_5,tplus_diagterm_5,
				  tminus_positions_5,n_tminus_positions_5,tminus_diagterm_5,
				  tplus_positions_3,n_tplus_positions_3,tplus_diagterm_3,
				  tminus_positions_3,n_tminus_positions_3,tminus_diagterm_3,
				 
				  queryuc_ptr_3,querylength3,
				  tplus_stream_array_3,tplus_streamsize_array_3,tplus_diagterm_array_3,
				  tminus_stream_array_3,tminus_streamsize_array_3,tminus_diagterm_array_3,
				  query3_compress_fwd,query3_compress_rev,
				  transcript_iit,transcriptome,transcriptomebits, 
				  nmismatches_allowed_3,/*concordantp*/false);
  }


  debug(printf("Filtering tplus and tminus transcriptome hits for concordant pairs\n"));
  filter_concordant_paths(tplus_gplus_paths_5,n_tplus_gplus_5,
			  tplus_gplus_paths_3,n_tplus_gplus_3);
  filter_concordant_paths(tplus_gminus_paths_5,n_tplus_gminus_5,
			  tplus_gminus_paths_3,n_tplus_gminus_3);
  filter_concordant_paths(tminus_gminus_paths_5,n_tminus_gminus_5,
			  tminus_gminus_paths_3,n_tminus_gminus_3);
  filter_concordant_paths(tminus_gplus_paths_5,n_tminus_gplus_5,
			  tminus_gplus_paths_3,n_tminus_gplus_3);


  debug(printf("Converting tplus paths for queryseq5 to hits\n"));
  *hits5_gplus = single_hits_gplus(&(*found_score_5),*hits5_gplus,tplus_gplus_paths_5,n_tplus_gplus_5,
				   querylength5,query5_compress_fwd,nmismatches_allowed_5,
				   /*sensedir*/SENSE_FORWARD,listpool,hitlistpool,level);
  *hits5_gminus = single_hits_gminus(&(*found_score_5),*hits5_gminus,tplus_gminus_paths_5,n_tplus_gminus_5,
				     querylength5,query5_compress_rev,nmismatches_allowed_5,
				     /*sensedir*/SENSE_FORWARD,listpool,hitlistpool,level);
  debug(printf("Converting tplus paths for queryseq3 to hits\n"));
  *hits3_gplus = single_hits_gplus(&(*found_score_3),*hits3_gplus,tplus_gplus_paths_3,n_tplus_gplus_3,
				   querylength3,query3_compress_fwd,nmismatches_allowed_3,
				   /*sensedir*/SENSE_FORWARD,listpool,hitlistpool,level);
  *hits3_gminus = single_hits_gminus(&(*found_score_3),*hits3_gminus,tplus_gminus_paths_3,n_tplus_gminus_3,
				     querylength3,query3_compress_rev,nmismatches_allowed_3,
				     /*sensedir*/SENSE_FORWARD,listpool,hitlistpool,level);
  FREE(tplus_gplus_paths_5);
  FREE(tplus_gminus_paths_5);
  FREE(tplus_gplus_paths_3);
  FREE(tplus_gminus_paths_3);
  
  debug(printf("Converting tminus paths for queryseq5 to hits\n"));
  *hits5_gplus = single_hits_gplus(&(*found_score_5),*hits5_gplus,tminus_gplus_paths_5,n_tminus_gplus_5,
				   querylength5,query5_compress_fwd,nmismatches_allowed_5,
				   /*sensedir*/SENSE_ANTI,listpool,hitlistpool,level);
  *hits5_gminus = single_hits_gminus(&(*found_score_5),*hits5_gminus,tminus_gminus_paths_5,n_tminus_gminus_5,
				     querylength5,query5_compress_rev,nmismatches_allowed_5,
				     /*sensedir*/SENSE_ANTI,listpool,hitlistpool,level);
  debug(printf("Converting tminus paths for queryseq3 to hits\n"));
  *hits3_gplus = single_hits_gplus(&(*found_score_3),*hits3_gplus,tminus_gplus_paths_3,n_tminus_gplus_3,
				   querylength3,query3_compress_fwd,nmismatches_allowed_3,
				   /*sensedir*/SENSE_ANTI,listpool,hitlistpool,level);
  *hits3_gminus = single_hits_gminus(&(*found_score_3),*hits3_gminus,tminus_gminus_paths_3,n_tminus_gminus_3,
				     querylength3,query3_compress_rev,nmismatches_allowed_3,
				     /*sensedir*/SENSE_ANTI,listpool,hitlistpool,level);
  
  FREE(tminus_gplus_paths_5);
  FREE(tminus_gminus_paths_5);
  FREE(tminus_gplus_paths_3);
  FREE(tminus_gminus_paths_3);

  return;
}
#endif



/* Simplified from search_transcriptome_ends */
List_T
Kmer_remap_transcriptome (char *remap_sequence, int remap_seqlength, 
			  Chrnum_T chrnum, Chrpos_T lowbound, Chrpos_T highbound,
			  Univ_IIT_T transcript_iit, 
			  Genome_T transcriptomebits, Transcriptome_T transcriptome) {
  List_T transcripts = NULL;
  Compress_T remap_compress_fwd, remap_compress_rev;

  UINT4 *tplus_positions_5 = NULL, *tminus_positions_5 = NULL,
    *tplus_positions_3 = NULL, *tminus_positions_3 = NULL;
  int n_tplus_positions_5 = 0, n_tminus_positions_5 = 0, n_tplus_positions_3 = 0, n_tminus_positions_3 = 0;
  int tplus_diagterm_5 = 0, tminus_diagterm_5 = 0, tplus_diagterm_3 = 0, tminus_diagterm_3 = 0;
				
  Univcoord_T left;
  UINT4 *tplus_diagonals, *tminus_diagonals;
  int n_tplus_diagonals, n_tminus_diagonals, i;

  Reader_T reader;
  int querypos, query_lastpos;
  Oligospace_T forward, revcomp, forward_oligo, revcomp_oligo;

#ifdef HAVE_AVX2
  __m256i _diff, _exonbounds, _trstart;
  int matchbits;
  int *ptr;
#elif defined(HAVE_SSE2)
  __m128i _diff, _exonbounds, _trstart;
  int matchbits;
  int *ptr;
#endif

  int trnum;
  int trstart, trend;
  int transcript_genestrand;
  Univcoord_T troffset, trhigh;
  Chrpos_T trlength, chrstart, chrend;
  int nmismatches;
  int *exonbounds;
  Chrpos_T *exonstarts;
  int nexons, exoni;
  int exonpos;


  debug9(printf("%s\n",remap_sequence));
  remap_compress_fwd = Compress_new_fwd(remap_sequence,remap_seqlength);
  remap_compress_rev = Compress_new_rev(remap_sequence,remap_seqlength);

  reader = Reader_new(remap_sequence,/*querystart*/0,/*queryend*/remap_seqlength);
  forward = revcomp = 0;
  if (Oligo_next_5(/*last_state*/INIT,&querypos,&forward,&revcomp,reader,/*genestrand*/0) != DONE) {
    forward_oligo = forward & oligobase_mask;
    debug9(printf("%s\n",Oligo_one_nt(forward_oligo,index1part_tr)));

    tplus_diagterm_5 = -querypos;
    n_tplus_positions_5 = Indexdb_ptr_with_diagterm(&tplus_positions_5,/*plus_indexdb*/indexdb_tr,forward_oligo,
						    tplus_diagterm_5);

    revcomp_oligo = (revcomp >> leftreadshift) & oligobase_mask;
    tminus_diagterm_5 = -remap_seqlength + querypos + index1part_tr;
    n_tminus_positions_5 = Indexdb_ptr_with_diagterm(&tminus_positions_5,/*minus_indexdb*/indexdb_tr,revcomp_oligo,
						     tminus_diagterm_5);

  }
  Reader_free(&reader);

  query_lastpos = remap_seqlength - index1part_tr;
  reader = Reader_new(remap_sequence,/*querystart*/query_lastpos - /*index1interval*/1 + 1,/*queryend*/remap_seqlength);
  forward = revcomp = 0;
  if (Oligo_next_5(/*last_state*/INIT,&querypos,&forward,&revcomp,reader,/*genestrand*/0) != DONE) {
    forward_oligo = forward & oligobase_mask;
    debug9(printf("%s\n",Oligo_one_nt(forward_oligo,index1part_tr)));
    tplus_diagterm_3 = -querypos;
    n_tplus_positions_3 = Indexdb_ptr_with_diagterm(&tplus_positions_3,/*plus_indexdb*/indexdb_tr,forward_oligo,
						    tplus_diagterm_3);

    revcomp_oligo = (revcomp >> leftreadshift) & oligobase_mask;
    tminus_diagterm_3 = -remap_seqlength + querypos + index1part_tr;
    n_tminus_positions_3 = Indexdb_ptr_with_diagterm(&tminus_positions_3,/*minus_indexdb*/indexdb_tr,revcomp_oligo,
						     tminus_diagterm_3);
  }
  Reader_free(&reader);
	 
  tplus_diagonals = Intersect_exact(&n_tplus_diagonals,
				    tplus_positions_5,n_tplus_positions_5,tplus_diagterm_5,
				    tplus_positions_3,n_tplus_positions_3,tplus_diagterm_3);
  tminus_diagonals = Intersect_exact(&n_tminus_diagonals,
				     tminus_positions_5,n_tminus_positions_5,tminus_diagterm_5,
				     tminus_positions_3,n_tminus_positions_3,tminus_diagterm_3);

  debug9(printf("***Remap transcriptome: %d plus and %d minus diagonals\n",n_tplus_diagonals,n_tminus_diagonals));

  for (i = 0; i < n_tplus_diagonals; i++) {
    left = (Univcoord_T) tplus_diagonals[i] /*- remap_seqlength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
    debug9(printf("REMAP PLUS DIAGONAL %u\n",tplus_diagonals[i]));
    if ((nmismatches = Genome_count_mismatches_substring(transcriptomebits,NULL,remap_compress_fwd,left,
							 /*pos5*/0,/*pos3*/remap_seqlength,
							 /*plusp*/true,/*genestrand*/0)) > 0) {
      debug9(printf("nmismatches %d\n",nmismatches));
    } else {
      /* Should replace with a rank/select structure */
      trnum = Univ_IIT_get_one(transcript_iit,left,left);
#ifdef DEBUG9
      printf("chrnum = %d\n",Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum));
      printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand);
#endif
      if (Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum) == chrnum) {
	Univ_IIT_interval_bounds_linear(&troffset,&trhigh,&trlength,transcript_iit,trnum);
	nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);
      
	trstart = left - troffset;
	trend = trstart + remap_seqlength;
	debug9(printf("trstart %d, trend %d\n",trstart,trend));

	exoni = 0;
#ifdef HAVE_AVX2
	_trstart = _mm256_set1_epi32(trstart);
	ptr = exonbounds;
	_exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	_diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	while (all_zero_p(_diff)) {
	  exoni += 8;
	  ptr += 8;
	  _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	  _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	}
	matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	_trstart = _mm_set1_epi32(trstart);
	ptr = exonbounds;
	_exonbounds = _mm_loadu_si128((__m128i *) ptr);
	_diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	while (all_zero_p(_diff)) {
	  exoni += 4;
	  ptr += 4;
	  _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	  _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	}
	matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	exoni += count_trailing_zeroes_32(matchbits);
#else
	while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	  exoni++;
	}
#endif

	if (exoni < nexons) {
	  if (exoni == 0) {
	    exonpos = trstart;
	  } else {
	    exonpos = trstart - exonbounds[exoni - 1];
	  }
	  if (transcript_genestrand > 0) {
	    chrstart = exonstarts[exoni] + exonpos - 1; /* gplus */
	  } else {
	    chrstart = exonstarts[exoni] - exonpos; /* gminus */
	  }
	    
	  while (/*exoni < nexons && */exonbounds[exoni] <= trend) {
	    exoni++;
	  }

	  if (exoni < nexons) {
	    if (exoni == 0) {
	      exonpos = trend;
	    } else {
	      exonpos = trend - exonbounds[exoni - 1];
	    }
	    if (transcript_genestrand > 0) {
	      chrend = exonstarts[exoni] + exonpos - 1; /* gplus */
	      assert(chrstart < chrend);
	      debug9(printf("Checking if interval %u..%u overlaps with desired %u..%u",
			    chrstart,chrend,lowbound,highbound));
	      if (chrstart <= highbound && chrend >= lowbound) {
		debug9(printf(" => yes"));
		transcripts = List_push(transcripts,(void *) Transcript_new(trnum,trstart,trend));
	      }
	      debug9(printf("\n"));

	    } else {
	      chrend = exonstarts[exoni] - exonpos; /* gminus */
	      assert(chrend < chrstart);
	      debug9(printf("Checking if interval %u..%u overlaps with desired %u..%u",
			    chrend,chrstart,lowbound,highbound));
	      if (chrend <= highbound && chrstart >= lowbound) {
		debug9(printf(" => yes"));
		transcripts = List_push(transcripts,(void *) Transcript_new(trnum,trstart,trend));
	      }
	      debug9(printf("\n"));
	    }
	  }
	}
      }
    }
  }
  FREE(tplus_diagonals);	/* Occupies memory even if n_tminus_diagonals == 0 */

  for (i = 0; i < n_tminus_diagonals; i++) {
    left = (Univcoord_T) tminus_diagonals[i] /*- remap_seqlength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
    debug9(printf("REMAP MINUS DIAGONAL %u\n",tminus_diagonals[i]));
    if ((nmismatches = Genome_count_mismatches_substring(transcriptomebits,NULL,remap_compress_rev,left,
							 /*pos5*/0,/*pos3*/remap_seqlength,
							 /*plusp*/false,/*genestrand*/0)) > 0) {
      debug9(printf("nmismatches %d\n",nmismatches));
    } else {
      /* Should replace with a rank/select structure */
      trnum = Univ_IIT_get_one(transcript_iit,left,left);
#ifdef DEBUG9
      printf("chrnum = %d\n",Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum));
      printf("trnum %d, transcript_genestrand %d\n",trnum,transcript_genestrand);
#endif
      if (Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum) == chrnum) {
	Univ_IIT_interval_bounds_linear(&troffset,&trhigh,&trlength,transcript_iit,trnum);
	nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);

	trend = left - troffset;
	trstart = trend + remap_seqlength;
	debug9(printf("trstart %d, trend %d\n",trstart,trend));

	exoni = 0;
#ifdef HAVE_AVX2
	_trstart = _mm256_set1_epi32(trend);
	ptr = exonbounds;
	_exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	_diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	while (all_zero_p(_diff)) {
	  exoni += 8;
	  ptr += 8;
	  _exonbounds = _mm256_loadu_si256((__m256i *) ptr);
	  _diff = _mm256_cmpgt_epi32(_exonbounds,_trstart);
	}
	matchbits = _mm256_movemask_ps(_mm256_castsi256_ps(_diff));
	exoni += count_trailing_zeroes_32(matchbits);
#elif defined(HAVE_SSE2)
	_trstart = _mm_set1_epi32(trend);
	ptr = exonbounds;
	_exonbounds = _mm_loadu_si128((__m128i *) ptr);
	_diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	while (all_zero_p(_diff)) {
	  exoni += 4;
	  ptr += 4;
	  _exonbounds = _mm_loadu_si128((__m128i *) ptr);
	  _diff = _mm_cmpgt_epi32(_exonbounds,_trstart);
	}
	matchbits = _mm_movemask_ps(_mm_castsi128_ps(_diff));
	exoni += count_trailing_zeroes_32(matchbits);
#else
	while (/*exoni < nexons && */exonbounds[exoni] <= trend) {
	  exoni++;
	}
#endif

	if (exoni < nexons) {
	  if (exoni == 0) {
	    exonpos = trend;
	  } else {
	    exonpos = trend - exonbounds[exoni - 1];
	  }
	  if (transcript_genestrand > 0) {
	    chrend = exonstarts[exoni] + exonpos - 1; /* gplus */
	  } else {
	    chrend = exonstarts[exoni] - exonpos; /* gminus */
	  }
	  while (/*exoni < nexons && */exonbounds[exoni] <= trstart) {
	    exoni++;
	  }

	  if (exoni < nexons) {
	    if (exoni == 0) {
	      exonpos = trstart;
	    } else {
	      exonpos = trstart - exonbounds[exoni - 1];
	    }
	    if (transcript_genestrand > 0) {
	      chrstart = exonstarts[exoni] + exonpos - 1; /* gplus */
	      assert(chrend < chrstart);
	      debug9(printf("Checking if interval %u..%u overlaps with desired %u..%u",
			    chrend,chrstart,lowbound,highbound));
	      if (chrend <= highbound && chrstart >= lowbound) {
		debug9(printf(" => yes"));
		transcripts = List_push(transcripts,(void *) Transcript_new(trnum,trstart,trend));
	      }
	      debug9(printf("\n"));
	    } else {
	      chrstart = exonstarts[exoni] - exonpos; /* gminus */
	      assert(chrstart < chrend);
	      debug9(printf("Checking if interval %u..%u overlaps with desired %u..%u",
			    chrstart,chrend,lowbound,highbound));
	      if (chrstart <= highbound && chrend >= lowbound) {
		debug9(printf(" => yes"));
		transcripts = List_push(transcripts,(void *) Transcript_new(trnum,trstart,trend));
	      }
	      debug9(printf("\n"));
	    }
	  }
	}
      }
    }
  }
  FREE(tminus_diagonals);	/* Occupies memory even if n_tminus_diagonals == 0 */

  Compress_free(&remap_compress_rev);
  Compress_free(&remap_compress_fwd);

  return transcripts;
}



#if 0
/* Divides hits into transcriptome (if remapped) and genome (otherwise) */
void
Kmer_divide_remapped_hits (List_T *transcriptome_hits, List_T *genome_hits, List_T hits,
			   Genome_T genomecomp,
			   Univ_IIT_T transcript_iit, Genome_T transcriptomebits, 
			   Transcriptome_T transcriptome, int genestrand) {
  Stage3end_T hit;
  List_T transcripts, p;

  *transcriptome_hits = *genome_hits = (List_T) NULL;

  for (p = hits; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_transcripts(hit) != NULL) {
      *transcriptome_hits = Hitlist_push(*transcriptome_hits,hitlistpool,(void *) hit);

    } else if ((transcripts = Kmer_remap_transcriptome(hit,transcript_iit,transcriptomebits,
						       transcriptome,genomecomp)) != NULL) {
      Stage3end_set_transcripts(hit,transcripts);
      *transcriptome_hits = Hitlist_push(*transcriptome_hits,hitlistpool,(void *) hit);

    } else {
      *genome_hits = Hitlist_push(*genome_hits,hitlistpool,(void *) hit);
    }
  }
  Hitlist_free(&hits);

  return;
}
#endif


static int
binary_search_univcoord (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi-1,positions[highi-1],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


static List_T
combine_ends_plus (int *found_score_overall, int *found_score_within_trims,
		   int *nhits, List_T hits, int querylength,
		   Univcoord_T left0, Univcoord_T left1,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		   int *mismatch_positions_alloc, Spliceinfo_T spliceinfo, Compress_T query_compress_fwd,
		   int genestrand, int nmismatches_allowed,
		   Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  Stage3end_T hit;

  Univcoord_T left, prev_left;
  Substring_T substring1, substring2;
  Univcoord_T deletionpos;
  int splice_distance;
  int querystart, queryend;
  bool abortp;
#ifdef DEBUG2
  Univcoord_T alignstart, alignend;
#endif

  int nmismatches1, nmismatches2, nmismatches_bothdiff;
  int introntype, sensedir;

  int best_knowni_i, best_knowni_j, best_nmismatches_i, best_nmismatches_j, best_nmismatches_indel;
  double best_prob_i, best_prob_j, donor_prob, acceptor_prob;
  int j;

  List_T junctions = NULL, substrings, junctions1, junctions2;
  int ninserts = 0;
  int splice_pos_1 = -1, splice_pos_2 = -1;
  int adj, nindels, indel_pos;

  debug2(printf("Entered combine_ends_plus with left %u and left %u\n",left0,left1));

  *nhits = 0;
  if ((adj = (int) (left1 - left0)) == 0) {
    left = left0;
    chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,left,
				 querylength,circular_typeint);
      
    if ((hit = Stage3end_new_substitution(&(*found_score_overall),&(*found_score_within_trims),
					  left,/*genomiclength*/querylength,
					  querylength,mismatch_positions_alloc,query_compress_fwd,
					  /*plusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
					  nmismatches_allowed,chrnum,chroffset,chrhigh,chrlength,
					  listpool,/*method*/KMER_APPROX,level)) != NULL) {
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      *nhits += 1;
    }

    if (novelsplicingp == true) {
      if ((hit = Stage3end_new_substitution(&(*found_score_overall),&(*found_score_within_trims),
					    left,/*genomiclength*/querylength,
					    querylength,mismatch_positions_alloc,query_compress_fwd,
					    /*plusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
					    nmismatches_allowed,chrnum,chroffset,chrhigh,chrlength,
					    listpool,/*method*/KMER_APPROX,level)) != NULL) {
	hits = Hitlist_push(hits,hitlistpool,(void *) hit);
	*nhits += 1;
      }
    }

  } else if (adj < 0) {
    nindels = -adj;
    if (nindels > 3) {
      adj = 0;
    } else if ((indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
							   /*left*/left0,/*indels*/-adj,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_fwd,
							   /*querystart*/0,/*queryend*/querylength,querylength,
							   nmismatches_allowed,/*plusp*/true,genestrand,
							   /*want_lowest_coordinate_p*/true)) <= 0) {
      adj = 0;
    } else {
      ninserts = nindels;
      /* total_nmatches = querylength - best_nmismatches_i - best_nmismatches_j - ninserts; */
      junctions = Listpool_push(NULL,listpool,(void *) Junction_new_insertion(nindels));
    }

  } else if (left1 > left0 + max_deletionlen) {
    if (novelsplicingp == false && splicesites == NULL) {
      /* Cannot be a splice */
      adj = 0;

    } else if (circularp[chrnum] == true) {
      /* Cannot be a splice */
      adj = 0;

    } else {
      /* Splice */
      debug2(printf("Appears to be a splice (plus)\n"));
    
      prev_left = left0;
      left = left1;
      assert(prev_left < left);

      spliceinfo->segmenti_donor_nknown = spliceinfo->segmenti_antiacceptor_nknown = 0;
      if (nsplicesites > 0 &&
	  Knownsplicing_splicesite_p(prev_left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search_univcoord(0,nsplicesites,splicesites,prev_left);
	while (j < nsplicesites && splicesites[j] < prev_left + querylength) {
	  if (splicetypes[j] == DONOR) {
	    debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	    spliceinfo->segmenti_donor_knownpos[spliceinfo->segmenti_donor_nknown] = splicesites[j] - prev_left;
	    spliceinfo->segmenti_donor_knowni[spliceinfo->segmenti_donor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIACCEPTOR) {
	    debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	    spliceinfo->segmenti_antiacceptor_knownpos[spliceinfo->segmenti_antiacceptor_nknown] = splicesites[j] - prev_left;
	    spliceinfo->segmenti_antiacceptor_knowni[spliceinfo->segmenti_antiacceptor_nknown++] = j;
	  }
	  j++;
	}
      }
      spliceinfo->segmenti_donor_knownpos[spliceinfo->segmenti_donor_nknown] = querylength + 100;
      spliceinfo->segmenti_antiacceptor_knownpos[spliceinfo->segmenti_antiacceptor_nknown] = querylength + 100;
      
      spliceinfo->segmentj_acceptor_nknown = spliceinfo->segmentj_antidonor_nknown = 0;
      if (nsplicesites > 0 &&
	  Knownsplicing_splicesite_p(left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search_univcoord(0,nsplicesites,splicesites,left);
	while (j < nsplicesites && splicesites[j] < left + querylength) {
	  if (splicetypes[j] == ACCEPTOR) {
	    debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	    spliceinfo->segmentj_acceptor_knownpos[spliceinfo->segmentj_acceptor_nknown] = splicesites[j] - left;
	    spliceinfo->segmentj_acceptor_knowni[spliceinfo->segmentj_acceptor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIDONOR) {
	    debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	    spliceinfo->segmentj_antidonor_knownpos[spliceinfo->segmentj_antidonor_nknown] = splicesites[j] - left;
	    spliceinfo->segmentj_antidonor_knowni[spliceinfo->segmentj_antidonor_nknown++] = j;
	  }
	  j++;
	}
      }
      spliceinfo->segmentj_acceptor_knownpos[spliceinfo->segmentj_acceptor_nknown] = querylength + 100;
      spliceinfo->segmentj_antidonor_knownpos[spliceinfo->segmentj_antidonor_nknown] = querylength + 100;
      
      splice_distance = (int) (left - prev_left);
      if ((splice_pos_1 = Splice_resolve_sense(&nindels,&indel_pos,&best_knowni_i,&best_knowni_j,
					       &best_nmismatches_i,&best_nmismatches_j,&best_nmismatches_indel,
					       &best_prob_i,&best_prob_j,
					       /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
					       /*querystart*/0,/*queryend*/querylength,querylength,query_compress_fwd,
					       spliceinfo,nmismatches_allowed,/*plusp*/true,genestrand,
					       /*max_deletionlen*/0,/*max_insertionlen*/0,
					       /*allow_indel_p*/false)) > 0) {
	junctions1 = Listpool_push(NULL,listpool,
				   (void *) Junction_new_splice(splice_distance,/*sensedir*/SENSE_FORWARD,
								/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
      }
      
      if ((splice_pos_2 = Splice_resolve_antisense(&nindels,&indel_pos,&best_knowni_i,&best_knowni_j,
						   &best_nmismatches_i,&best_nmismatches_j,&best_nmismatches_indel,
						   &best_prob_i,&best_prob_j,
						   /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
						   /*querystart*/0,/*queryend*/querylength,querylength,query_compress_fwd,
						   spliceinfo,nmismatches_allowed,/*plusp*/true,genestrand,
						   /*max_deletionlen*/0,/*max_insertionlen*/0,
						   /*allow_indel_p*/false)) > 0) {
	junctions2 = Listpool_push(NULL,listpool,
				   (void *) Junction_new_splice(splice_distance,/*sensedir*/SENSE_ANTI,
								/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
      }

      adj = 0;
    }
      
    
  } else if ((indel_pos = Indel_resolve_middle_deletion_or_splice(&introntype,&best_nmismatches_i,&best_nmismatches_j,
								  /*left*/left0,/*indels*/-adj,
								  /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_fwd,
								  /*querystart*/0,/*queryend*/querylength,querylength,
								  nmismatches_allowed,/*plusp*/true,genestrand,
								  min_intronlength,/*want_lowest_coordinate_p*/true)) <= 0) {
    adj = 0;

  } else if ((Chrpos_T) (nindels = adj) < min_intronlength) {
    /* Cannot be an intron, so must be a deletion */
    deletionpos = left0 + indel_pos; /* gplus */
    junctions = Listpool_push(NULL,listpool,(void *) Junction_new_deletion(nindels,deletionpos));
    
  } else if ((sensedir = Intron_canonical_sensedir(introntype)) == SENSE_NULL) {
    /* No intron dinucleotids found, so must be a deletion */
    deletionpos = left0 + indel_pos; /* gplus */
    junctions = Listpool_push(NULL,listpool,(void *) Junction_new_deletion(nindels,deletionpos));
    
  } else if (novelsplicingp == false && splicesites == NULL) {
    /* Cannot be a splice */
    adj = 0;

  } else if (circularp[chrnum] == true) {
    /* Cannot be a splice */
    adj = 0;

  } else if (sensedir == SENSE_FORWARD) {
    /* Rare condition, so not tested */
    /* splice_distance = (int) (left1 - left0); */
    deletionpos = left0 + indel_pos; /* gplus */
    donor_prob = Maxent_hr_donor_prob(deletionpos,chroffset);
    acceptor_prob = Maxent_hr_acceptor_prob(deletionpos+nindels,chroffset);
    junctions1 = Listpool_push(NULL,listpool,(void *) Junction_new_splice(/*splice_distance*/nindels,SENSE_FORWARD,donor_prob,acceptor_prob));
    splice_pos_1 = indel_pos;
    adj = 0;
    
  } else {
    /* Rare condition, so not tested */
    /* splice_distance = (int) (left1 - left0); */
    deletionpos = left0 + indel_pos; /* gplus */
    donor_prob = Maxent_hr_antidonor_prob(deletionpos+nindels,chroffset);
    acceptor_prob = Maxent_hr_antiacceptor_prob(deletionpos,chroffset);
    junctions2 = Listpool_push(NULL,listpool,(void *) Junction_new_splice(/*splice_distance*/nindels,SENSE_ANTI,donor_prob,acceptor_prob));
    splice_pos_2 = indel_pos;
    adj = 0;
  }
  
  debug2(printf("splice_pos_1 %d, splice_pos_2 %d\n",splice_pos_1,splice_pos_2));
  if (splice_pos_1 > 0) {
    abortp = false;
    substrings = (List_T) NULL;
    
    /* second part */
    left = left1;
    querystart = splice_pos_1;
    queryend = querylength;
#ifdef DEBUG2
    alignstart = left + querystart;
    alignend = left + queryend;
    printf("(1) gplus 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches2 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_fwd,left,
					/*pos5*/querystart,/*pos3*/queryend,/*plusp*/true,genestrand);
    debug2(printf("nmismatches fwd from %d to %d: %d\n",querystart,queryend,nmismatches2));
    
    if ((substring2 = Substring_new(nmismatches2,left,querystart,queryend,querylength,
				    /*plusp*/true,genestrand,query_compress_fwd,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_FORWARD)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querystart,queryend));
      substrings = Listpool_push(NULL,listpool,(void *) substring2);
    }
    
    /* first part */
    left = left0;
    querystart = 0;
    queryend = splice_pos_1;
#ifdef DEBUG2
    alignstart = left /*+ querystart (0)*/;
    alignend = left + queryend;
    printf("(1) gplus 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches1 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_fwd,left,
					/*pos5*/querystart,/*pos3*/queryend,/*plusp*/true,genestrand);
    debug2(printf("nmismatches fwd from %d to %d: %d\n",querystart,queryend,nmismatches1));
    
    if ((substring1 = Substring_new(nmismatches1,left,querystart,queryend,querylength,
				    /*plusp*/true,genestrand,query_compress_fwd,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_FORWARD)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querystart,queryend));
      substrings = Listpool_push(substrings,listpool,(void *) substring1);
    }
    
    /* nmismatches_whole = */ nmismatches_bothdiff = nmismatches1 + nmismatches2;
    if (abortp == true ||
	(hit = Stage3end_new_precomputed(&(*found_score_overall),&(*found_score_within_trims),nmismatches_bothdiff,
					 substrings,junctions1,/*transcripts*/NULL,/*transcripts_other*/NULL,
					 querylength,chrnum,chroffset,chrhigh,chrlength,
					 /*gplusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
					 listpool,/*method*/KMER_APPROX,level)) == NULL) {
      Junction_gc(&junctions1);
      Substring_free(&substring2);
      Substring_free(&substring1);
      /* List_free(&substrings); -- allocated by Listpool_push */
    } else {
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      *nhits += 1;
    }
  }
  
  if (splice_pos_2 > 0) {
    abortp = false;
    substrings = (List_T) NULL;
    
    /* second part */
    left = left1;
    querystart = splice_pos_2;
    queryend = querylength;
#ifdef DEBUG2
    alignstart = left + querystart;
    alignend = left + queryend;
    printf("(2) gplus 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches2 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_fwd,left,
					/*pos5*/querystart,/*pos3*/queryend,/*plusp*/true,genestrand);
    debug2(printf("nmismatches fwd from %d to %d: %d\n",querystart,queryend,nmismatches2));
    
    if ((substring2 = Substring_new(nmismatches2,left,querystart,queryend,querylength,
				    /*plusp*/true,genestrand,query_compress_fwd,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_ANTI)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querystart,queryend));
      substrings = Listpool_push(NULL,listpool,(void *) substring2);
    }
    
    /* first part */
    left = left0;
    querystart = 0;
    queryend = splice_pos_2;
#ifdef DEBUG2
    alignstart = left /*+ querystart (0)*/;
    alignend = left + queryend;
    printf("(2) gplus 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches1 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_fwd,left,
					/*pos5*/querystart,/*pos3*/queryend,/*plusp*/true,genestrand);
    debug2(printf("nmismatches fwd from %d to %d: %d\n",querystart,queryend,nmismatches1));
    
    if ((substring1 = Substring_new(nmismatches1,left,querystart,queryend,querylength,
				    /*plusp*/true,genestrand,query_compress_fwd,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_ANTI)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querystart,queryend));
      substrings = Listpool_push(substrings,listpool,(void *) substring1);
    }
    
    /* nmismatches_whole = */ nmismatches_bothdiff = nmismatches1 + nmismatches2;
    if (abortp == true ||
	(hit = Stage3end_new_precomputed(&(*found_score_overall),&(*found_score_within_trims),nmismatches_bothdiff,
					 substrings,junctions2,/*transcripts*/NULL,/*transcripts_other*/NULL,
					 querylength,chrnum,chroffset,chrhigh,chrlength,
					 /*gplusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
					 listpool,/*method*/KMER_APPROX,level)) == NULL) {
      Junction_gc(&junctions2);
      Substring_free(&substring2);
      Substring_free(&substring1);
      /* List_free(&substrings); -- allocated by Listpool_push */
    } else {
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      *nhits += 1;
    }
  }
  
  if (adj != 0) {
    /* Indel */
    abortp = false;
    substrings = (List_T) NULL;
    
    /* second part */
    left = left1;
    querystart = indel_pos + ninserts;
    queryend = querylength;
#ifdef DEBUG2
    alignstart = left + querystart;
    alignend = left + queryend;
    printf("(3) gplus 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches2 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_fwd,left,
					/*pos5*/querystart,/*pos3*/queryend,/*plusp*/true,genestrand);
    debug2(printf("nmismatches fwd from %d to %d: %d\n",querystart,queryend,nmismatches2));
    
    if ((substring2 = Substring_new(nmismatches2,left,querystart,queryend,querylength,
				    /*plusp*/true,genestrand,query_compress_fwd,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_NULL)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querystart,queryend));
      substrings = Listpool_push(NULL,listpool,(void *) substring2);
    }
    
    /* first part */
    left = left0;
    querystart = 0;
    queryend = indel_pos;
#ifdef DEBUG2
    alignstart = left /*+ querystart (0)*/;
    alignend = left + queryend;
    printf("(3) gplus 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches1 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_fwd,left,
					/*pos5*/querystart,/*pos3*/queryend,/*plusp*/true,genestrand);
    debug2(printf("nmismatches fwd from %d to %d: %d\n",querystart,queryend,nmismatches1));
    
    if ((substring1 = Substring_new(nmismatches1,left,querystart,queryend,querylength,
				    /*plusp*/true,genestrand,query_compress_fwd,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_NULL)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querystart,queryend));
      substrings = Listpool_push(substrings,listpool,(void *) substring1);
    }
    
    /*nmismatches_whole = */ nmismatches_bothdiff = nmismatches1 + nmismatches2;
    if (abortp == true ||
	(hit = Stage3end_new_precomputed(&(*found_score_overall),&(*found_score_within_trims),nmismatches_bothdiff,
					 substrings,junctions,/*transcripts*/NULL,/*transcripts_other*/NULL,
					 querylength,chrnum,chroffset,chrhigh,chrlength,
					 /*gplusp*/true,genestrand,/*sensedir*/SENSE_NULL,
					 listpool,/*method*/KMER_APPROX,level)) == NULL) {
      Junction_gc(&junctions);
      Substring_free(&substring2);
      Substring_free(&substring1);
      /* List_free(&substrings); -- allocated by Listpool_push */
    } else {
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      *nhits += 1;
    }
  }

  return hits;
}


static List_T
combine_ends_minus (int *found_score_overall, int *found_score_within_trims,
		    int *nhits, List_T hits, int querylength,
		    Univcoord_T left0, Univcoord_T left1,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		    int *mismatch_positions_alloc, Spliceinfo_T spliceinfo, Compress_T query_compress_rev,
		    int genestrand, int nmismatches_allowed,
		    Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  Stage3end_T hit;

  Univcoord_T left, prev_left;
  Substring_T substring1, substring2;
  Univcoord_T deletionpos;
  int splice_distance;
  int querystart, queryend;
  bool abortp;
#ifdef DEBUG2
  Univcoord_T alignstart, alignend;
#endif

  int nmismatches1, nmismatches2, nmismatches_bothdiff;
  int introntype, sensedir;

  int best_knowni_i, best_knowni_j, best_nmismatches_i, best_nmismatches_j, best_nmismatches_indel;
  double best_prob_i, best_prob_j, donor_prob, acceptor_prob;
  int j;

  List_T junctions = NULL, substrings, junctions1, junctions2;
  int ninserts = 0;
  int splice_pos_1 = -1, splice_pos_2 = -1;
  int adj, nindels, indel_pos;

  *nhits = 0;

  debug2(printf("left0 %u, left1 %u\n",left0,left1));
  if ((adj = (int) left1 - left0) == 0) {
    left = left1;
    chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,left,
				 querylength,circular_typeint);
      
    if ((hit = Stage3end_new_substitution(&(*found_score_overall),&(*found_score_within_trims),
					  left,/*genomiclength*/querylength,
					  querylength,mismatch_positions_alloc,query_compress_rev,
					  /*plusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
					  nmismatches_allowed,chrnum,chroffset,chrhigh,chrlength,
					  listpool,/*method*/KMER_APPROX,level)) != NULL) {
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      *nhits += 1;
    }

    if (novelsplicingp == true) {
      if ((hit = Stage3end_new_substitution(&(*found_score_overall),&(*found_score_within_trims),
					    left,/*genomiclength*/querylength,
					    querylength,mismatch_positions_alloc,query_compress_rev,
					    /*plusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
					    nmismatches_allowed,chrnum,chroffset,chrhigh,chrlength,
					    listpool,/*method*/KMER_APPROX,level)) != NULL) {
	hits = Hitlist_push(hits,hitlistpool,(void *) hit);
	*nhits += 1;
      }
    }
      
  } else if (adj < 0) {
    /* TODO: Check if should be left1 */
    nindels = -adj;
    if (nindels > 3) {
      adj = 0;
    } else if ((indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
							   /*left*/left0,/*indels*/-adj,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_rev,
							   /*querystart*/0,/*queryend*/querylength,querylength,
							   nmismatches_allowed,/*plusp*/false,genestrand,
							   /*want_lowest_coordinate_p*/true)) <= 0) {
      adj = 0;
    } else {
      ninserts = nindels;
      /* total_nmatches = querylength - best_nmismatches_i - best_nmismatches_j - ninserts; */
      junctions = Listpool_push(NULL,listpool,(void *) Junction_new_insertion(nindels));
    }
    
  } else if (left1 > left0 + max_deletionlen) {
    if (novelsplicingp == false && splicesites == NULL) {
      /* Cannot be a splice */
      adj = 0;

    } else if (circularp[chrnum] == true) {
      /* Cannot be a splice */
      adj = 0;

    } else {
      /* Splice */
      debug2(printf("Appears to be a splice (minus)\n"));
    
      prev_left = left0;
      left = left1;
    
      spliceinfo->segmenti_donor_nknown = spliceinfo->segmenti_antiacceptor_nknown = 0;
      if (nsplicesites > 0 &&
	  Knownsplicing_splicesite_p(prev_left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search_univcoord(0,nsplicesites,splicesites,prev_left);
	while (j < nsplicesites && splicesites[j] < prev_left + querylength) {
	  if (splicetypes[j] == DONOR) {
	    debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	    spliceinfo->segmenti_donor_knownpos[spliceinfo->segmenti_donor_nknown] = splicesites[j] - prev_left;
	    spliceinfo->segmenti_donor_knowni[spliceinfo->segmenti_donor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIACCEPTOR) {
	    debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	    spliceinfo->segmenti_antiacceptor_knownpos[spliceinfo->segmenti_antiacceptor_nknown] = splicesites[j] - prev_left;
	    spliceinfo->segmenti_antiacceptor_knowni[spliceinfo->segmenti_antiacceptor_nknown++] = j;
	  }
	  j++;
	}
      }
      spliceinfo->segmenti_donor_knownpos[spliceinfo->segmenti_donor_nknown] = querylength + 100;
      spliceinfo->segmenti_antiacceptor_knownpos[spliceinfo->segmenti_antiacceptor_nknown] = querylength + 100;
    
      spliceinfo->segmentj_acceptor_nknown = spliceinfo->segmentj_antidonor_nknown = 0;
      if (nsplicesites > 0 &&
	  Knownsplicing_splicesite_p(left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search_univcoord(0,nsplicesites,splicesites,left);
	while (j < nsplicesites && splicesites[j] < left + querylength) {
	  if (splicetypes[j] == ACCEPTOR) {
	    debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	    spliceinfo->segmentj_acceptor_knownpos[spliceinfo->segmentj_acceptor_nknown] = splicesites[j] - left;
	    spliceinfo->segmentj_acceptor_knowni[spliceinfo->segmentj_acceptor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIDONOR) {
	    debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	    spliceinfo->segmentj_antidonor_knownpos[spliceinfo->segmentj_antidonor_nknown] = splicesites[j] - left;
	    spliceinfo->segmentj_antidonor_knowni[spliceinfo->segmentj_antidonor_nknown++] = j;
	  }
	  j++;
	}
      }
      spliceinfo->segmentj_acceptor_knownpos[spliceinfo->segmentj_acceptor_nknown] = querylength + 100;
      spliceinfo->segmentj_antidonor_knownpos[spliceinfo->segmentj_antidonor_nknown] = querylength + 100;
    
      splice_distance = (int) (left - prev_left);
      if ((splice_pos_1 = Splice_resolve_sense(&nindels,&indel_pos,&best_knowni_i,&best_knowni_j,
					       &best_nmismatches_i,&best_nmismatches_j,&best_nmismatches_indel,
					       &best_prob_i,&best_prob_j,
					       /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
					       /*querystart*/0,/*queryend*/querylength,querylength,query_compress_rev,
					       spliceinfo,nmismatches_allowed,/*plusp*/false,genestrand,
					       /*max_deletionlen*/0,/*max_insertionlen*/0,
					       /*allow_indel_p*/false)) > 0) {
	junctions1 = Listpool_push(NULL,listpool,
				   (void *) Junction_new_splice(splice_distance,/*sensedir*/SENSE_FORWARD,
								/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
      }
    
      if ((splice_pos_2 = Splice_resolve_antisense(&nindels,&indel_pos,&best_knowni_i,&best_knowni_j,
						   &best_nmismatches_i,&best_nmismatches_j,&best_nmismatches_indel,
						   &best_prob_i,&best_prob_j,
						   /*segmenti_left*/prev_left,/*segmentj_left*/left,chroffset,chroffset,
						   /*querystart*/0,/*queryend*/querylength,querylength,query_compress_rev,
						   spliceinfo,nmismatches_allowed,/*plusp*/false,genestrand,
						   /*max_deletionlen*/0,/*max_insertionlen*/0,
						   /*allow_indel_p*/false)) > 0) {
	junctions2 = Listpool_push(NULL,listpool,
				   (void *) Junction_new_splice(splice_distance,/*sensedir*/SENSE_ANTI,
								/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
      }
    
      adj = 0;
    }
    
  } else if ((indel_pos = Indel_resolve_middle_deletion_or_splice(&introntype,&best_nmismatches_i,&best_nmismatches_j,
								  /*left*/left0,/*indels*/-adj,
								  /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress_rev,
								  /*querystart*/0,/*queryend*/querylength,querylength,
								  nmismatches_allowed,/*plusp*/false,genestrand,
								  min_intronlength,/*want_lowest_coordinate_p*/true)) <= 0) {
    adj = 0;
  } else if ((Chrpos_T) (nindels = adj) < min_intronlength) {
    /* Cannot be an intron, so must be a deletion */
    deletionpos = left1 + indel_pos - nindels; /* gminus */
    junctions = Listpool_push(NULL,listpool,(void *) Junction_new_deletion(nindels,deletionpos));
    
  } else if ((sensedir = Intron_canonical_sensedir(introntype)) == SENSE_NULL) {
    /* No intron dinucleotids found, so must be a deletion */
    deletionpos = left1 + indel_pos - nindels; /* gminus */
    junctions = Listpool_push(NULL,listpool,(void *) Junction_new_deletion(nindels,deletionpos));
    
  } else if (novelsplicingp == false && splicesites == NULL) {
    /* Cannot be a splice */
    adj = 0;

  } else if (circularp[chrnum] == true) {
    /* Cannot be a splice */
    adj = 0;

  } else if (sensedir == SENSE_FORWARD) {
    /* Rare condition, so not tested */
    /* splice_distance = (int) (left1 - left0); */
    deletionpos = left1 + indel_pos - nindels; /* gminus */
    donor_prob = Maxent_hr_antidonor_prob(deletionpos+nindels,chroffset);
    acceptor_prob = Maxent_hr_antiacceptor_prob(deletionpos,chroffset);
    junctions1 = Listpool_push(NULL,listpool,
			       (void *) Junction_new_splice(/*splice_distance*/nindels,SENSE_FORWARD,donor_prob,acceptor_prob));
    splice_pos_1 = indel_pos;
    adj = 0;
    
  } else {
    /* Rare condition, so not tested */
    /* splice_distance = (int) (left1 - left0); */
    deletionpos = left1 + indel_pos - nindels; /* gminus */
    donor_prob = Maxent_hr_donor_prob(deletionpos,chroffset);
    acceptor_prob = Maxent_hr_acceptor_prob(deletionpos+nindels,chroffset);
    junctions2 = Listpool_push(NULL,listpool,
			       (void *) Junction_new_splice(/*splice_distance*/nindels,SENSE_ANTI,donor_prob,acceptor_prob));
    splice_pos_2 = indel_pos;
    adj = 0;
  }
  
  debug2(printf("splice_pos_1 %d, splice_pos_2 %d\n",splice_pos_1,splice_pos_2));
  if (splice_pos_1 > 0) {
    abortp = false;
    substrings = (List_T) NULL;
    
    /* first part */
    left = left0;
    queryend = splice_pos_1;
    querystart = 0;
#ifdef DEBUG2
    alignstart = left + queryend;
    alignend = left /*+ querystart (0)*/;
    printf("(1) gminus 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches1 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_rev,left,
					/*pos5*/querystart,/*pos3*/queryend,
					/*plusp*/false,genestrand);
    debug2(printf("mismatches rev from %d to %d: %d\n",querylength - queryend,querylength - querystart,nmismatches1));

    if ((substring1 = Substring_new(nmismatches1,left,querylength - queryend,querylength - querystart,querylength,
				    /*plusp*/false,genestrand,query_compress_rev,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_FORWARD)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querylength - queryend,querylength - querystart));
      substrings = Listpool_push(NULL,listpool,(void *) substring1);
    }
    
    /* second part */
    left = left1;
    queryend = querylength;
    querystart = splice_pos_1;
#ifdef DEBUG2
    alignstart = left + queryend;
    alignend = left + querystart;
    printf("(1) gminus 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches2 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_rev,left,
					/*pos5*/querystart,/*pos3*/queryend,
					/*plusp*/false,genestrand);
    debug2(printf("mismatches rev from %d to %d: %d\n",querystart,queryend,nmismatches2));
    if ((substring2 = Substring_new(nmismatches2,left,querylength - queryend,querylength - querystart,querylength,
				    /*plusp*/false,genestrand,query_compress_rev,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_FORWARD)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querylength - queryend,querylength - querystart));
      substrings = Listpool_push(substrings,listpool,(void *) substring2);
    }
    
    /* nmismatches_whole = */ nmismatches_bothdiff = nmismatches1 + nmismatches2;
    if (abortp == true ||
	(hit = Stage3end_new_precomputed(&(*found_score_overall),&(*found_score_within_trims),nmismatches_bothdiff,
					 substrings,junctions1,/*transcripts*/NULL,/*transcripts_other*/NULL,
					 querylength,chrnum,chroffset,chrhigh,chrlength,
					 /*gplusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
					 listpool,/*method*/KMER_APPROX,level)) == NULL) {
      Junction_gc(&junctions1);
      Substring_free(&substring2);
      Substring_free(&substring1);
      /* List_free(&substrings); -- allocated by Listpool_push */
    } else {
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      *nhits += 1;
    }
  }
  
  if (splice_pos_2 > 0) {
    abortp = false;
    substrings = (List_T) NULL;
    
    /* first part */
    left = left0;
    queryend = splice_pos_2;
    querystart = 0;
#ifdef DEBUG2
    alignstart = left + queryend;
    alignend = left /*+ querystart (0)*/;
    printf("(2) gminus 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches1 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_rev,left,
					/*pos5*/querystart,/*pos3*/queryend,
					/*plusp*/false,genestrand);
    debug2(printf("mismatches rev from %d to %d: %d\n",querylength - queryend,querylength - querystart,nmismatches1));
    if ((substring1 = Substring_new(nmismatches1,left,querylength - queryend,querylength - querystart,querylength,
				    /*plusp*/false,genestrand,query_compress_rev,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_ANTI)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querylength - queryend,querylength - querystart));
      substrings = Listpool_push(NULL,listpool,(void *) substring1);
    }
    
    /* second part */
    left = left1;
    queryend = querylength;
    querystart = splice_pos_2;
#ifdef DEBUG2
    alignstart = left + queryend;
    alignend = left + querystart;
    printf("(2) gminus 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches2 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_rev,left,
					/*pos5*/querystart,/*pos3*/queryend,
					/*plusp*/false,genestrand);
    debug2(printf("mismatches rev from %d to %d: %d\n",querystart,queryend,nmismatches2));
    if ((substring2 = Substring_new(nmismatches2,left,querylength - queryend,querylength - querystart,querylength,
				    /*plusp*/false,genestrand,query_compress_rev,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_ANTI)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querylength - queryend,querylength - querystart));
      substrings = Listpool_push(substrings,listpool,(void *) substring2);
    }
    
    /* nmismatches_whole = */ nmismatches_bothdiff = nmismatches1 + nmismatches2;
    if (abortp == true ||
	(hit = Stage3end_new_precomputed(&(*found_score_overall),&(*found_score_within_trims),nmismatches_bothdiff,
					 substrings,junctions2,/*transcripts*/NULL,/*transcripts_other*/NULL,
					 querylength,chrnum,chroffset,chrhigh,chrlength,
					 /*gplusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
					 listpool,/*method*/KMER_APPROX,level)) == NULL) {
      Junction_gc(&junctions2);
      Substring_free(&substring2);
      Substring_free(&substring1);
      /* List_free(&substrings); -- allocated by Listpool_push */
    } else {
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      *nhits += 1;
    }
  }
  
  if (adj != 0) {
    /* Indel */
    abortp = false;
    substrings = (List_T) NULL;
    
    /* first part */
    left = left0;
    queryend = indel_pos;
    querystart = 0;
#ifdef DEBUG2
    alignstart = left + queryend;
    alignend = left /*+ querystart (0)*/;
    printf("(3) gminus 1st part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches1 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_rev,left,
					/*pos5*/querystart,/*pos3*/queryend,
					/*plusp*/false,genestrand);
    debug2(printf("mismatches rev from %d to %d: %d\n",querylength - queryend,querylength - querystart,nmismatches1));
    if ((substring1 = Substring_new(nmismatches1,left,querylength - queryend,querylength - querystart,querylength,
				    /*plusp*/false,genestrand,query_compress_rev,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_NULL)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querylength - queryend,querylength - querystart));
      substrings = Listpool_push(NULL,listpool,(void *) substring1);
    }
    
    /* second part */
    left = left1;
    queryend = querylength;
    querystart = indel_pos + ninserts;
#ifdef DEBUG2
    alignstart = left + queryend;
    alignend = left + querystart;
    printf("(3) gminus 2nd part query %d..%d, align %u..%u [%u..%u], left %u = %u - (%d - %d)\n",
	   querystart,queryend,alignstart,alignend,alignstart - chroffset,alignend - chroffset,
	   left - chroffset,alignend - chroffset,querylength,queryend);
#endif
    nmismatches2 =
      Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress_rev,left,
					/*pos5*/querystart,/*pos3*/queryend,
					/*plusp*/false,genestrand);
    debug2(printf("mismatches rev from %d to %d: %d\n",querystart,queryend,nmismatches2));
    if ((substring2 = Substring_new(nmismatches2,left,querylength - queryend,querylength - querystart,querylength,
				    /*plusp*/false,genestrand,query_compress_rev,
				    chrnum,chroffset,chrhigh,chrlength,
				    /*splice_querystart_p*/false,/*splicetype_querystart*/NO_SPLICE,/*ambig_prob_querystart*/0.0,
				    /*splice_queryend_p*/false,/*splicetype_queryend*/NO_SPLICE,/*ambig_prob_queryend*/0.0,
				    /*sensedir*/SENSE_NULL)) == NULL) {
      debug2(printf("ABORTING BECAUSE OF MISMATCHES\n"));
      abortp = true;
    } else {
      debug2(printf(" => Created substring for %d..%d\n",querylength - queryend,querylength - querystart));
      substrings = Listpool_push(substrings,listpool,(void *) substring2);
    }
    
    /* nmismatches_whole = */ nmismatches_bothdiff = nmismatches1 + nmismatches2;
    if (abortp == true ||
	(hit = Stage3end_new_precomputed(&(*found_score_overall),&(*found_score_within_trims),nmismatches_bothdiff,
					 substrings,junctions,/*transcripts*/NULL,/*transcripts_other*/NULL,
					 querylength,chrnum,chroffset,chrhigh,chrlength,
					 /*gplusp*/false,genestrand,/*sensedir*/SENSE_NULL,
					 listpool,/*method*/KMER_APPROX,level)) == NULL) {
      Junction_gc(&junctions);
      Substring_free(&substring2);
      Substring_free(&substring1);
      /* List_free(&substrings); -- allocated by Listpool_push */
    } else {
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
      *nhits += 1;
    }
  }

  return hits;
}


/* Searches the ends, just like ultrafast transcriptome */
/* Need max_hits, because repetitive reads can give many exact matches in a genome */
void
Kmer_search_genome_ends_exact (bool *abort_exact_p, int *found_score_overall, int *found_score_within_trims,
			       List_T *hits_gplus, List_T *hits_gminus,
			       Stage1_T stage1, int querylength, int *mismatch_positions_alloc,
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       int genestrand, int nmismatches_allowed,
			       Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  int nhits = 0;
  Stage3end_T hit;
  int mod5, mod3, adj;

  Univcoord_T left, *diagonals;
  int ndiagonals, i;

  Chrnum_T chrnum;
  Chrpos_T chrlength;
  Univcoord_T chroffset, chrhigh;


  /* max_hits = 1000000; */

  *abort_exact_p = false;

  /* Plus */
  for (mod5 = 0; mod5 < index1interval; mod5++) {
    adj = (querylength - index1part) % index1interval;
    mod3 = (index1interval + adj - mod5) % index1interval;

#ifdef LARGE_GENOMES
    diagonals = Intersect_exact_large(&ndiagonals,
				      stage1->plus_rawpositions_high_5[mod5],stage1->plus_rawpositions_5[mod5],
				      stage1->plus_nrawpositions_5[mod5],stage1->plus_diagterms_5[mod5],
				      stage1->plus_rawpositions_high_3[mod3],stage1->plus_rawpositions_3[mod3],
				      stage1->plus_nrawpositions_3[mod3],stage1->plus_diagterms_3[mod3]);
#else
    diagonals = Intersect_exact(&ndiagonals,
				stage1->plus_rawpositions_5[mod5],stage1->plus_nrawpositions_5[mod5],stage1->plus_diagterms_5[mod5],
				stage1->plus_rawpositions_3[mod3],stage1->plus_nrawpositions_3[mod3],stage1->plus_diagterms_3[mod3]);
#endif
    debug4(printf("plus mod5 %d, mod3 %d: %d diagonals\n",mod5,mod3,ndiagonals));

    i = 0;
    while (/*nhits <= max_hits && */ i < ndiagonals) {
      left = diagonals[i] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
      debug4a(printf("ULTRAFAST PLUS DIAGONAL %u\n",left));
      debug4a(printf(" => nmismatches %d\n",nmismatches));
      chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,left,
				   querylength,circular_typeint);
	    
      if ((hit = Stage3end_new_substitution(&(*found_score_overall),&(*found_score_within_trims),
					    left,/*genomiclength*/querylength,
					    querylength,mismatch_positions_alloc,query_compress_fwd,
					    /*plusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
					    nmismatches_allowed,chrnum,chroffset,chrhigh,chrlength,
					    listpool,/*method*/KMER_EXACT,level)) != NULL) {
	*hits_gplus = Hitlist_push(*hits_gplus,hitlistpool,(void *) hit);
	nhits++;
      }

      if (novelsplicingp == true) {
	if ((hit = Stage3end_new_substitution(&(*found_score_overall),&(*found_score_within_trims),
					      left,/*genomiclength*/querylength,
					      querylength,mismatch_positions_alloc,query_compress_fwd,
					      /*plusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
					      nmismatches_allowed,chrnum,chroffset,chrhigh,chrlength,
					      listpool,/*method*/KMER_EXACT,level)) != NULL) {
	  *hits_gplus = Hitlist_push(*hits_gplus,hitlistpool,(void *) hit);
	  nhits++;
	}
      }
      i++;
    }
    FREE(diagonals);


    /* Minus */
#ifdef LARGE_GENOMES
    diagonals = Intersect_exact_large(&ndiagonals,
				      stage1->minus_rawpositions_high_5[mod5],stage1->minus_rawpositions_5[mod5],
				      stage1->minus_nrawpositions_5[mod5],stage1->minus_diagterms_5[mod5],
				      stage1->minus_rawpositions_high_3[mod3],stage1->minus_rawpositions_3[mod3],
				      stage1->minus_nrawpositions_3[mod3],stage1->minus_diagterms_3[mod3]);
#else
    diagonals = Intersect_exact(&ndiagonals,
				stage1->minus_rawpositions_5[mod5],stage1->minus_nrawpositions_5[mod5],stage1->minus_diagterms_5[mod5],
				stage1->minus_rawpositions_3[mod3],stage1->minus_nrawpositions_3[mod3],stage1->minus_diagterms_3[mod3]);
#endif
    debug4(printf("minus mod5 %d, mod3 %d: %d diagonals\n",mod5,mod3,ndiagonals));

    i = 0;
    while (/*nhits <= max_hits && */ i < ndiagonals) {
      left = diagonals[i] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
      debug4a(printf("ULTRAFAST MINUS DIAGONAL %u\n",left));
      debug4a(printf(" => nmismatches %d\n",nmismatches));
      chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,left,
				   querylength,circular_typeint);
	    
      if ((hit = Stage3end_new_substitution(&(*found_score_overall),&(*found_score_within_trims),
					    left,/*genomiclength*/querylength,
					    querylength,mismatch_positions_alloc,query_compress_rev,
					    /*plusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
					    nmismatches_allowed,chrnum,chroffset,chrhigh,chrlength,
					    listpool,/*method*/KMER_EXACT,level)) != NULL) {
	*hits_gminus = Hitlist_push(*hits_gminus,hitlistpool,(void *) hit);
	nhits++;
      }

      if (novelsplicingp == true) {
	if ((hit = Stage3end_new_substitution(&(*found_score_overall),&(*found_score_within_trims),
					      left,/*genomiclength*/querylength,
					      querylength,mismatch_positions_alloc,query_compress_rev,
					      /*plusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
					      nmismatches_allowed,chrnum,chroffset,chrhigh,chrlength,
					      listpool,/*method*/KMER_EXACT,level)) != NULL) {
	  *hits_gminus = Hitlist_push(*hits_gminus,hitlistpool,(void *) hit);
	  nhits++;
	}
      }
      i++;
    }
    FREE(diagonals);

  }

#if 0
  if (nhits > max_hits) {
    debug4(printf("Kmer_search_ends_exact aborting because nhits %d > max_hits %d\n",nhits,max_hits));
    Stage3end_gc(*hits_gplus);
    Hitlist_free(&(*hits_gplus));
    Stage3end_gc(*hits_gminus);
    Hitlist_free(&(*hits_gminus));
    *abort_exact_p = true;	/* Indicates that we should not try to run approx algorithm */
  }
#endif

  debug4(printf("Kmer_search_genome_ends_exact returning %d plus and %d minus hits\n",
		List_length(*hits_gplus),List_length(*hits_gminus)));

  return;
}


/* Performs a merge.  About 10% faster than the alternative of trying all combinations of mod5 and mod3. */
/* Need max_hits, because repetitive reads can give many exact matches in a genome */
void
Kmer_search_genome_ends_approx (int *found_score_overall, int *found_score_within_trims,
				List_T *hits_gplus, List_T *hits_gminus, Stage1_T stage1,
				Compress_T query_compress_fwd, Compress_T query_compress_rev,
				int querylength, int genestrand, int nmismatches_allowed, int sizelimit,
				Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  int nhits = 0, nhits0;
  int gplus_streami_5 = 0, gminus_streami_5 = 0, gplus_streami_3 = 0, gminus_streami_3 = 0;


  int total_npositions_gplus_5 = 0, total_npositions_gminus_5 = 0,
    total_npositions_gplus_3 = 0, total_npositions_gminus_3 = 0;

  Univcoord_T *gplus_diagonals_5, *gminus_diagonals_5, *gplus_diagonals_3, *gminus_diagonals_3;
  int n_gplus_diagonals_5, n_gminus_diagonals_5, n_gplus_diagonals_3, n_gminus_diagonals_3;

  bool gplus_exactp, gminus_exactp;
  Univcoord_T *gplus_diagpairs, *gminus_diagpairs;
  int n_gplus_diagpairs, n_gminus_diagpairs;

  int i, k;
  Univcoord_T lefta, leftb, left0, left1;

  int adj;
  int trim_a5, trim_a3, trim_b5, trim_b3;
  int nmismatches_a5, nmismatches_a3, nmismatches_b5, nmismatches_b3;

  Chrnum_T chrnum;
  Chrpos_T chrlength;
  Univcoord_T chroffset, chrhigh;


  /* max_hits = 1000000; */

  for (i = 0; i < index1interval; i++) {
    if (stage1->plus_nrawpositions_5[i] > 0) {
#ifdef LARGE_GENOMES
      stage1->gplus_stream_high_array_5[gplus_streami_5] = stage1->plus_rawpositions_high_5[i];
      stage1->gplus_stream_low_array_5[gplus_streami_5] = stage1->plus_rawpositions_5[i];
#else
      stage1->gplus_stream_array_5[gplus_streami_5] = stage1->plus_rawpositions_5[i];
#endif
      stage1->gplus_streamsize_array_5[gplus_streami_5] = stage1->plus_nrawpositions_5[i];
      stage1->gplus_diagterm_array_5[gplus_streami_5] = stage1->plus_diagterms_5[i];
      total_npositions_gplus_5 += stage1->plus_nrawpositions_5[i];
      gplus_streami_5++;
    }

    if (stage1->minus_nrawpositions_5[i] > 0) {
#ifdef LARGE_GENOMES
      stage1->gminus_stream_high_array_5[gminus_streami_5] = stage1->minus_rawpositions_high_5[i];
      stage1->gminus_stream_low_array_5[gminus_streami_5] = stage1->minus_rawpositions_5[i];
#else
      stage1->gminus_stream_array_5[gminus_streami_5] = stage1->minus_rawpositions_5[i];
#endif
      stage1->gminus_streamsize_array_5[gminus_streami_5] = stage1->minus_nrawpositions_5[i];
      stage1->gminus_diagterm_array_5[gminus_streami_5] = stage1->minus_diagterms_5[i];
      total_npositions_gminus_5 += stage1->minus_nrawpositions_5[i];
      gminus_streami_5++;
    }

    if (stage1->plus_nrawpositions_3[i] > 0) {
#ifdef LARGE_GENOMES
      stage1->gplus_stream_high_array_3[gplus_streami_3] = stage1->plus_rawpositions_high_3[i];
      stage1->gplus_stream_low_array_3[gplus_streami_3] = stage1->plus_rawpositions_3[i];
#else
      stage1->gplus_stream_array_3[gplus_streami_3] = stage1->plus_rawpositions_3[i];
#endif
      stage1->gplus_streamsize_array_3[gplus_streami_3] = stage1->plus_nrawpositions_3[i];
      stage1->gplus_diagterm_array_3[gplus_streami_3] = stage1->plus_diagterms_3[i];
      total_npositions_gplus_3 += stage1->plus_nrawpositions_3[i];
      gplus_streami_3++;
    }

    if (stage1->minus_nrawpositions_3[i] > 0) {
#ifdef LARGE_GENOMES
      stage1->gminus_stream_high_array_3[gminus_streami_3] = stage1->minus_rawpositions_high_3[i];
      stage1->gminus_stream_low_array_3[gminus_streami_3] = stage1->minus_rawpositions_3[i];
#else
      stage1->gminus_stream_array_3[gminus_streami_3] = stage1->minus_rawpositions_3[i];
#endif
      stage1->gminus_streamsize_array_3[gminus_streami_3] = stage1->minus_nrawpositions_3[i];
      stage1->gminus_diagterm_array_3[gminus_streami_3] = stage1->minus_diagterms_3[i];
      total_npositions_gminus_3 += stage1->minus_nrawpositions_3[i];
      gminus_streami_3++;
    }
  }

  debug4(printf("Comparing total_npositions_gplus_5 %d against sizelmit %d\n",total_npositions_gplus_5,sizelimit));
  if (0 && total_npositions_gplus_5 > sizelimit) {
    gplus_diagonals_5 = (Univcoord_T *) NULL;
    n_gplus_diagonals_5 = 0;
  } else {
#ifdef LARGE_GENOMES
    gplus_diagonals_5 = Merge_diagonals_large(&n_gplus_diagonals_5,
					      stage1->gplus_stream_high_array_5,stage1->gplus_stream_low_array_5,
					      stage1->gplus_streamsize_array_5,stage1->gplus_diagterm_array_5,
					      /*nstreams*/gplus_streami_5);
#else
    gplus_diagonals_5 = Merge_diagonals(&n_gplus_diagonals_5,stage1->gplus_stream_array_5,
					stage1->gplus_streamsize_array_5,stage1->gplus_diagterm_array_5,
					/*nstreams*/gplus_streami_5);
#endif
  }


  debug4(printf("Comparing total_npositions_gminus_5 %d against sizelmit %d\n",total_npositions_gminus_5,sizelimit));
  if (0 && total_npositions_gminus_5 > sizelimit) {
    gminus_diagonals_5 = (Univcoord_T *) NULL;
    n_gminus_diagonals_5 = 0;
  } else {
#ifdef LARGE_GENOMES
    gminus_diagonals_5 = Merge_diagonals_large(&n_gminus_diagonals_5,
					       stage1->gminus_stream_high_array_5,stage1->gminus_stream_low_array_5,
					       stage1->gminus_streamsize_array_5,stage1->gminus_diagterm_array_5,
					       /*nstreams*/gminus_streami_5);
#else
    gminus_diagonals_5 = Merge_diagonals(&n_gminus_diagonals_5,stage1->gminus_stream_array_5,
					 stage1->gminus_streamsize_array_5,stage1->gminus_diagterm_array_5,
					 /*nstreams*/gminus_streami_5);
#endif
  }


  debug4(printf("Comparing total_npositions_gplus_3 %d against sizelmit %d\n",total_npositions_gplus_3,sizelimit));
  if (0 && total_npositions_gplus_3 > sizelimit) {
    gplus_diagonals_3 = (Univcoord_T *) NULL;
    n_gplus_diagonals_3 = 0;
  } else {
#ifdef LARGE_GENOMES
    gplus_diagonals_3 = Merge_diagonals_large(&n_gplus_diagonals_3,
					      stage1->gplus_stream_high_array_3,stage1->gplus_stream_low_array_3,
					      stage1->gplus_streamsize_array_3,stage1->gplus_diagterm_array_3,
					      /*nstreams*/gplus_streami_3);
#else
    gplus_diagonals_3 = Merge_diagonals(&n_gplus_diagonals_3,stage1->gplus_stream_array_3,
					stage1->gplus_streamsize_array_3,stage1->gplus_diagterm_array_3,
					/*nstreams*/gplus_streami_3);
#endif

  }


  debug4(printf("Comparing total_npositions_gminus_3 %d against sizelmit %d\n",total_npositions_gminus_3,sizelimit));
  if (0 && total_npositions_gminus_3 > sizelimit) {
    gminus_diagonals_3 = (Univcoord_T *) NULL;
    n_gminus_diagonals_3 = 0;
  } else {
#ifdef LARGE_GENOMES
    gminus_diagonals_3 = Merge_diagonals_large(&n_gminus_diagonals_3,
					       stage1->gminus_stream_high_array_3,stage1->gminus_stream_low_array_3,
					       stage1->gminus_streamsize_array_3,stage1->gminus_diagterm_array_3,
					       /*nstreams*/gminus_streami_3);
#else
    gminus_diagonals_3 = Merge_diagonals(&n_gminus_diagonals_3,stage1->gminus_stream_array_3,
					 stage1->gminus_streamsize_array_3,stage1->gminus_diagterm_array_3,
					 /*nstreams*/gminus_streami_3);
#endif
  }

  gplus_diagpairs = Intersect_approx_simple(&gplus_exactp,&n_gplus_diagpairs,
					    gplus_diagonals_5,n_gplus_diagonals_5,
					    gplus_diagonals_3,n_gplus_diagonals_3,
					    /*maxdistance*/overall_max_distance);
  gminus_diagpairs = Intersect_approx_simple(&gminus_exactp,&n_gminus_diagpairs,
					     gminus_diagonals_5,n_gminus_diagonals_5,
					     gminus_diagonals_3,n_gminus_diagonals_3,
					     /*maxdistance*/overall_max_distance);
  debug4(printf("***Intersect ends approx: exactp %d and %d.  %d plus and %d minus diagpairs\n",
		gplus_exactp,gminus_exactp,n_gplus_diagpairs,n_gminus_diagpairs));

#if !defined(LARGE_GENOMES) || defined(HAVE_AVX512) || defined(HAVE_AVX2)
  FREE_ALIGN(gplus_diagonals_5);
  FREE_ALIGN(gminus_diagonals_5);
  FREE_ALIGN(gplus_diagonals_3);
  FREE_ALIGN(gminus_diagonals_3);
#else
  FREE(gplus_diagonals_5);
  FREE(gminus_diagonals_5);
  FREE(gplus_diagonals_3);
  FREE(gminus_diagonals_3);
#endif

  i = k = 0;
  while (/*nhits <= max_hits && */ i < n_gplus_diagpairs) {
    lefta = (Univcoord_T) gplus_diagpairs[k] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
    leftb = (Univcoord_T) gplus_diagpairs[k+1] /*- querylength*/ /*- 0*/; /* NEW FORMULA for querystart of 0 */
    /* printf("plus lefta %llu, leftb %llu\n",lefta,leftb); */
	  
#if 0
    trim_ends(&trim_a5,&trim_a3,&nmismatches_a5,&nmismatches_a3,poly_p,
	      query_compress_fwd,lefta,querylength,genomebits,/*plusp*/true,index1part);
#else
    trim_a5 = Genome_first_kmer_left(&nmismatches_a5,genomebits,query_compress_fwd,lefta,/*pos5*/0,/*pos3*/querylength,
				     /*plusp*/true,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part);
    trim_a3 = querylength - Genome_first_kmer_right(&nmismatches_a3,genomebits,query_compress_fwd,lefta,/*pos5*/0,/*pos3*/querylength,
						    /*plusp*/true,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part);
#endif
	  
#if 0
    trim_ends(&trim_b5,&trim_b3,&nmismatches_b5,&nmismatches_b3,poly_p,
	      query_compress_fwd,leftb,querylength,genomebits,/*plusp*/true,index1part);
#else
    trim_b5 = Genome_first_kmer_left(&nmismatches_b5,genomebits,query_compress_fwd,leftb,/*pos5*/0,/*pos3*/querylength,
				     /*plusp*/true,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part);
    trim_b3 = querylength - Genome_first_kmer_right(&nmismatches_b3,genomebits,query_compress_fwd,leftb,/*pos5*/0,/*pos3*/querylength,
						    /*plusp*/true,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part);
#endif
	  
    debug4a(printf("trimmed a %d..%d, trimmed b %d..%d\n",
		   trim_a5,querylength - trim_a3,trim_b5,querylength - trim_b3));
    if (trim_a5 == 0 || trim_b3 == querylength) {
      left0 = lefta;
      left1 = leftb;
      adj = (int) (leftb - lefta);
      debug4a(printf("Setting a first, b second, adj %d\n",adj));
	    
    } else if (trim_a3 == querylength || trim_b5 == 0) {
      left0 = leftb;
      left1 = lefta;
      adj = (int) (lefta - leftb);
      debug4a(printf("Setting b first, a second, adj %d\n",adj));
	    
    } else {
      adj = 0;
    }
	  
    /* printf("plus adj %d, left0 %llu, left1 %llu\n",adj,left0,left1); */
    if (adj != 0) {
      if (left1 > left0) {
	chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,left0,
				     querylength,circular_typeint);
	      
	if (left1 + querylength <= chrhigh) {
	  *hits_gplus = combine_ends_plus(&(*found_score_overall),&(*found_score_within_trims),
					  &nhits0,*hits_gplus,querylength,left0,left1,
					  chrnum,chroffset,chrhigh,chrlength,stage1->mismatch_positions_alloc,
					  stage1->spliceinfo,query_compress_fwd,genestrand,nmismatches_allowed,
					  listpool,hitlistpool,level);
	  nhits += nhits0;
	}
	      
      } else {
	chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,left1,
				     querylength,circular_typeint);
	      
	if (left0 + querylength <= chrhigh) {
	  *hits_gplus = combine_ends_plus(&(*found_score_overall),&(*found_score_within_trims),
					  &nhits0,*hits_gplus,querylength,left0,left1,
					  chrnum,chroffset,chrhigh,chrlength,stage1->mismatch_positions_alloc,
					  stage1->spliceinfo,query_compress_fwd,genestrand,nmismatches_allowed,
					  listpool,hitlistpool,level);
	  nhits += nhits0;
	}
      }
    }
    i++;
    k += 2;
  }
      

  i = k = 0;
  while (/*nhits <= max_hits && */ i < n_gminus_diagpairs) {
    lefta = (Univcoord_T) gminus_diagpairs[k] /*- querylength*/; /* NEW FORMULA for queryend of querylength */
    leftb = (Univcoord_T) gminus_diagpairs[k+1] /*- querylength*/; /* NEW FORMULA for queryend of querylength */
    /* printf("minus lefta %llu, leftb %llu\n",lefta,leftb); */
	  
#if 0
    trim_ends(&trim_a5,&trim_a3,&nmismatches_a5,&nmismatches_a3,poly_p,
	      query_compress_rev,lefta,querylength,genomebits,/*plusp*/false,index1part);
#else
    trim_a5 = Genome_first_kmer_left(&nmismatches_a5,genomebits,query_compress_rev,lefta,/*pos5*/0,/*pos3*/querylength,
				     /*plusp*/false,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part);
    trim_a3 = querylength - Genome_first_kmer_right(&nmismatches_a3,genomebits,query_compress_rev,lefta,/*pos5*/0,/*pos3*/querylength,
						    /*plusp*/false,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part);
#endif
	  
#if 0
    trim_ends(&trim_b5,&trim_b3,&nmismatches_b5,&nmismatches_b3,poly_p,
	      query_compress_rev,leftb,querylength,genomebits,/*plusp*/false,index1part);
#else
    trim_b5 = Genome_first_kmer_left(&nmismatches_b5,genomebits,query_compress_rev,leftb,/*pos5*/0,/*pos3*/querylength,
				     /*plusp*/false,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part);
    trim_b3 = querylength - Genome_first_kmer_right(&nmismatches_b3,genomebits,query_compress_rev,leftb,/*pos5*/0,/*pos3*/querylength,
						    /*plusp*/false,genestrand,/*query_unk_mismatch_local_p*/false,/*kmer*/index1part);
#endif
	  
    debug4a(printf("trimmed a %d..%d, trimmed b %d..%d\n",
		   trim_a5,querylength - trim_a3,trim_b5,querylength - trim_b3));
    if (trim_a5 == 0 || trim_b3 == querylength) {
      left0 = lefta;
      left1 = leftb;
      adj = (int) (leftb - lefta);
      debug4a(printf("Setting a first, b second, adj %d\n",adj));
	    
    } else if (trim_a3 == querylength || trim_b5 == 0) {
      left0 = leftb;
      left1 = lefta;
      adj = (int) (lefta - leftb);
      debug4a(printf("Setting b first, a second, adj %d\n",adj));
	    
    } else {
      adj = 0;
    }
	  
    /* printf("minus adj %d, left0 %llu, left1 %llu\n",adj,left0,left1); */
    if (adj != 0) {
      if (left1 > left0) {
	chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,left0,
				     querylength,circular_typeint);
	      
	if (left1 + querylength <= chrhigh) {
	  *hits_gminus = combine_ends_minus(&(*found_score_overall),&(*found_score_within_trims),
					    &nhits0,*hits_gminus,querylength,left0,left1,
					    chrnum,chroffset,chrhigh,chrlength,stage1->mismatch_positions_alloc,
					    stage1->spliceinfo,query_compress_rev,genestrand,nmismatches_allowed,
					    listpool,hitlistpool,level);
	  nhits += nhits0;
	}
	      
      } else {
	chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,left1,
				     querylength,circular_typeint);
	      
	if (left0 + querylength <= chrhigh) {
	  *hits_gminus = combine_ends_minus(&(*found_score_overall),&(*found_score_within_trims),
					    &nhits0,*hits_gminus,querylength,left0,left1,
					    chrnum,chroffset,chrhigh,chrlength,stage1->mismatch_positions_alloc,
					    stage1->spliceinfo,query_compress_rev,genestrand,nmismatches_allowed,
					    listpool,hitlistpool,level);
	  nhits += nhits0;
	}
      }
    }

    i++;
    k += 2;
  }

  FREE(gminus_diagpairs);
  FREE(gplus_diagpairs);

#if 0
  if (nhits > max_hits) {
    Stage3end_gc(*hits_gplus);
    Hitlist_free(&(*hits_gplus));
    Stage3end_gc(*hits_gminus);
    Hitlist_free(&(*hits_gminus));
  }
#endif

  debug4(printf("Done with Kmer_search_genome_ends_approx\n"));
  return;
}


#if 0
static Oligospace_T
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
#endif


void
Kmer_search_setup (Mode_T mode_in,
		   int index1part_tr_in, int index1part_in, int index1interval_in, int local1part_in,
		   int min_intronlength_in, int max_deletionlength,
		   Univ_IIT_T chromosome_iit_in, int circular_typeint_in, bool *circularp_in, 
		   Genome_T genomebits_in, Genome_T genomebits_alt_in,
		   Indexdb_T indexdb_in, Indexdb_T indexdb2_in, Indexdb_T indexdb_tr_in, Chrpos_T shortsplicedist_in,
		   bool novelsplicingp_in, Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		   Chrpos_T *splicedists_in, int nsplicesites_in) {

  mode = mode_in;

  chromosome_iit = chromosome_iit_in;
  circular_typeint = circular_typeint_in;
  circularp = circularp_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

  index1part_tr = index1part_tr_in;
  index1part = index1part_in;
  index1interval = index1interval_in;
  local1part = local1part_in;

#ifdef HAVE_64_BIT
  leftreadshift = 64 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#else
  leftreadshift = 32 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#endif

  indexdb = indexdb_in;
  indexdb2 = indexdb2_in;
  indexdb_tr = indexdb_tr_in;

  min_intronlength = min_intronlength_in;
  max_deletionlen = max_deletionlength;

  shortsplicedist = shortsplicedist_in;

  overall_max_distance = shortsplicedist;
  if ((Chrpos_T) max_deletionlen > overall_max_distance) {
    overall_max_distance = max_deletionlen;
  }

  novelsplicingp = novelsplicingp_in;
  splicesites = splicesites_in;
  splicetypes = splicetypes_in;
  splicedists = splicedists_in;
  nsplicesites = nsplicesites_in;

  return;
}


