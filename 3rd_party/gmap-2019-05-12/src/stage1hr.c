static char rcsid[] = "$Id: stage1hr.c 219219 2019-05-12 22:27:20Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "stage1hr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset() */
#include <math.h>
#include <ctype.h>		/* for tolower() */
#include "assert.h"
#include "mem.h"
#include "types.h"		/* Needed for HAVE_64_BIT */
#include "univcoord.h"

#include "reader.h"
#include "oligo.h"
#include "indexdb.h"

#include "list.h"
#include "intlist.h"
#include "splice.h"
#include "indel.h"
#include "stage3hr.h"
#include "concordance.h"
#include "substring.h"
#include "complement.h"
#include "compress.h"
#include "genome128_hr.h"
#include "genome_sites.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "iitdef.h"
#include "univinterval.h"
#ifdef LARGE_GENOMES
#include "uint8list.h"
#else
#include "uintlist.h"
#endif
#include "record.h"

#ifdef HAVE_64_BIT
/* Oligospace_T is 64 bits */
#include "uint8table_rh.h"
#else
/* Oligospace_T is 32 bits */
#include "uinttable_rh.h"
#endif


#include "orderstat.h"
#include "path-solve.h"
#include "kmer-search.h"
#include "extension-search.h"
#include "segment-search.h"
#include "terminal.h"
#include "distant-rna.h"
#include "distant-dna.h"

#include "univdiag.h"
#include "univdiagdef.h"

#include "comp.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#endif


#ifdef HAVE_64_BIT
#ifdef LARGE_GENOMES
#else
#define DIAGONAL_ADD_QUERYPOS 1
#endif
#endif

/* Three methods for performing a multiway merge.  Need to define one below. */

#define SPEED 1


/* Note FORMULA: formulas for querypos <-> diagonal (diagterm in call to Indexdb_read) are:

plus: diagonal = position + querylength - querypos
minus: diagonal = position + querypos + index1part

For minus, the index1part is needed in call to Indexdb_read because
position is stored at beginning of plus oligomer, which corresponds to
end of minus oligomer.  As a result, we have the following formulas:

high genomic position = diagonal (corresponds to querypos =
querylength for plus, and querypos = 0 for minus)

low genomic position = diagonal - querylength (corresponds to querypos
= 0 for plus, and querypos = querylength for minus)

Holds when we use Reader_T to read from 5' end of forward query and 3'
end of revcomp query simultaneously.  If we create a queryrc sequence,
then we can use just the plus formula, and convert the query
coordinates later.

*/


#define MAX_ITER 3		/* For looping through Segment_search */

#define NO_EXTENSIONS_BEFORE_ZERO 1

#define ALLOW_MIDDLE_ALIGNMENTS 1

/* #define EXTRACT_GENOMICSEG 1 */
#ifdef EXTRACT_GENOMICSEG
#define MAX_INDEXSIZE 8
#endif


/* MAX_NALIGNMENTS of 2 vs 1 gets 1600 improvements in 275,000 reads */
/* MAX_NALIGNMENTS of 3 vs 2 gets 96 improvements in 275,000 reads */
#define MAX_NALIGNMENTS 3

#define MAX_ALLOCATION 200

#define PAIRMAX_ADDITIONAL 10000 /* Allows for finding of unpaired GMAP alignments beyond pairmax */

#define MIN_SIZELIMIT 100

/* static int kmer_search_sizelimit = 100; */
/* static int stage1hr_sizelimit = 3000; */
/* static int extension_search_sizelimit = 3000; */

#define POLY_A 0x00000000
#define POLY_T 0xFFFFFFFF
static Oligospace_T poly_a;
static Oligospace_T poly_t;

static bool use_only_transcriptome_p;
static bool remap_transcriptome_p = false;

static Univ_IIT_T transcript_iit;
static Transcriptome_T transcriptome;
static Genome_T transcriptomebits;

static Indexdb_T indexdb_fwd;
static Indexdb_T indexdb_rev;
static Indexdb_T indexdb_tr;
static Localdb_T localdb;

static int index1part;
static int index1part_tr;
static int index1interval;

static double user_maxlevel_float;
static double user_mincoverage_float;

static Univ_IIT_T chromosome_iit;
static int circular_typeint;

static int nchromosomes;
static Genome_T genomecomp;
static Genome_T genomebits;
static Genome_T genomebits_alt;

static int leftreadshift;
static Oligospace_T oligobase_mask; /* same as kmer_mask */


/* Mode */
static Mode_T mode;
static bool snpp;
static int maxpaths_search;
static int maxpaths_report;

/* For spliceable (really "joinable", if we consider indels) */
static Chrpos_T overall_max_distance;

/* Other distances */
bool novelsplicingp;
static Chrpos_T shortsplicedist;


/* Penalties */
static int subopt_levels;

static bool find_dna_chimeras_p;
static bool distances_observed_p;

static Chrpos_T min_intronlength;
static Chrpos_T expected_pairlength;
static Chrpos_T pairlength_deviation;

static int distantsplicing_penalty;
static int min_distantsplicing_end_matches;
static int min_distantsplicing_identity;


#define A_CHAR 0x0
#define C_CHAR 0x1
#define G_CHAR 0x2
#define T_CHAR 0x3


/* Originally allowed only 1, to print only unique translocations.
   But need to allow enough to avoid missing some translocations. */
/* For transcript splicing, need to increase MAXCHIMERAPATHS */
/* #define MAXCHIMERAPATHS 100 */
#define MAXCHIMERAPATHS 10000

#define NREQUIRED_FAST 2	/* For candidate generation using
				   multimiss.  A value of 2 implies 
				   specificity of a 24-mer, which
				   should be low for a human-sized
				   genome */

#define MAX_INDEX1INTERVAL 3



/* Overall flow */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Stage1_init */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Filling oligos */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* consolidate_paired_results and choose_among_paired */ 
#ifdef DEBUG16
#define debug16(x) x
#else
#define debug16(x)
#endif


#define T Stage1_T
static T
Stage1_new (int querylength) {
  T new = (T) MALLOC(sizeof(*new));
  int overhang = index1interval-1;
  int mod;

  new->plus_validp = (bool *) CALLOC(querylength+overhang,sizeof(bool));
  new->minus_validp = (bool *) CALLOC(querylength+overhang,sizeof(bool));

  new->plus_oligos = (Oligospace_T *) MALLOC((querylength+overhang)*sizeof(Oligospace_T));
  new->minus_oligos = (Oligospace_T *) MALLOC((querylength+overhang)*sizeof(Oligospace_T));

  new->retrievedp_allocated = (bool *) CALLOC(2 * (querylength+overhang),sizeof(bool));
  new->plus_retrievedp = &(new->retrievedp_allocated[overhang]);
  new->minus_retrievedp = &(new->retrievedp_allocated[(querylength+overhang)+overhang]);

#ifdef LARGE_GENOMES
  new->positions_high_allocated = (unsigned char **) CALLOC(2 * (querylength+overhang),sizeof(unsigned char *));
  new->plus_positions_high = &(new->positions_high_allocated[overhang]);
  new->minus_positions_high = &(new->positions_high_allocated[(querylength+overhang)+overhang]);
#endif

  new->positions_allocated = (UINT4 **) CALLOC(2 * (querylength+overhang),sizeof(UINT4 *));
  new->plus_positions = &(new->positions_allocated[overhang]);
  new->minus_positions = &(new->positions_allocated[(querylength+overhang)+overhang]);

  new->npositions_allocated = (int *) CALLOC(2 * (querylength+overhang),sizeof(int));
  new->plus_npositions = &(new->npositions_allocated[overhang]);
  new->minus_npositions = &(new->npositions_allocated[(querylength+overhang)+overhang]);

  /* Need to allocate (max_mismatches + 2), where max_mismatches is provided to Genome_mismatches_left or Genome_mismatches_right */
  new->mismatch_positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  new->positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  new->spliceinfo = Spliceinfo_new(querylength);

  /* Memory allocated for Segment_identify in segment-search.c, and
     Merge_diagonals in kmer-search.c (which needs four sets of
     arrays) */
#ifdef LARGE_GENOMES
  new->stream_high_alloc = (unsigned char **) MALLOC(4*querylength*sizeof(unsigned char *));
  new->gplus_stream_high_array_5 = &(new->stream_high_alloc[0]);
  new->gminus_stream_high_array_5 = &(new->stream_high_alloc[querylength]);
  new->gplus_stream_high_array_3 = &(new->stream_high_alloc[2*querylength]);
  new->gminus_stream_high_array_3 = &(new->stream_high_alloc[3*querylength]);

  new->stream_low_alloc = (UINT4 **) MALLOC(4*querylength*sizeof(UINT4 *));
  new->gplus_stream_low_array_5 = &(new->stream_low_alloc[0]);
  new->gminus_stream_low_array_5 = &(new->stream_low_alloc[querylength]);
  new->gplus_stream_low_array_3 = &(new->stream_low_alloc[2*querylength]);
  new->gminus_stream_low_array_3 = &(new->stream_low_alloc[3*querylength]);
#endif

  new->stream_alloc = (Univcoord_T **) MALLOC(4*querylength*sizeof(Univcoord_T *));
  new->gplus_stream_array_5 = &(new->stream_alloc[0]);
  new->gminus_stream_array_5 = &(new->stream_alloc[querylength]);
  new->gplus_stream_array_3 = &(new->stream_alloc[2*querylength]);
  new->gminus_stream_array_3 = &(new->stream_alloc[3*querylength]);

#ifdef LARGE_GENOMES
  new->tplus_stream_array = new->gplus_stream_low_array_5;
  new->tminus_stream_array = new->gminus_stream_low_array_5;
#else
  new->tplus_stream_array = new->gplus_stream_array_5;
  new->tminus_stream_array = new->gminus_stream_array_5;
#endif

  new->streamsize_alloc = (int *) MALLOC(4*querylength*sizeof(int));
  new->tplus_streamsize_array = new->gplus_streamsize_array_5 = &(new->streamsize_alloc[0]);
  new->tminus_streamsize_array = new->gminus_streamsize_array_5 = &(new->streamsize_alloc[querylength]);
  new->gplus_streamsize_array_3 = &(new->streamsize_alloc[2*querylength]);
  new->gminus_streamsize_array_3 = &(new->streamsize_alloc[3*querylength]);

  new->querypos_diagterm_alloc = (int *) MALLOC(4*querylength*sizeof(int));
  new->tplus_diagterm_array = new->gplus_diagterm_array_5 = &(new->querypos_diagterm_alloc[0]);
  new->tminus_diagterm_array = new->gminus_diagterm_array_5 = &(new->querypos_diagterm_alloc[querylength]);
  new->gplus_diagterm_array_3 = &(new->querypos_diagterm_alloc[2*querylength]);
  new->gminus_diagterm_array_3 = &(new->querypos_diagterm_alloc[3*querylength]);


  for (mod = 0; mod < 2*index1interval; mod++) {
#ifdef LARGE_GENOMES
    new->plus_rawpositions_high_5[mod] = (unsigned char *) NULL;
    new->minus_rawpositions_high_5[mod] = (unsigned char *) NULL;
    new->plus_rawpositions_high_3[mod] = (unsigned char *) NULL;
    new->minus_rawpositions_high_3[mod] = (unsigned char *) NULL;
#endif
    new->plus_rawpositions_5[mod] = (UINT4 *) NULL;
    new->minus_rawpositions_5[mod] = (UINT4 *) NULL;
    new->plus_rawpositions_3[mod] = (UINT4 *) NULL;
    new->minus_rawpositions_3[mod] = (UINT4 *) NULL;
  }

  /* Uses Listpool_T procedures */
  new->queryfwd_plus_set = (List_T) NULL;
  new->queryfwd_minus_set = (List_T) NULL;
  new->queryrev_plus_set = (List_T) NULL;
  new->queryrev_minus_set = (List_T) NULL;

  return new;
}

#ifdef DEBUG
static void
Stage1_dump (T this, int querylength) {
  int query_lastpos = querylength - index1part, querypos;
  int i;

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->plus_retrievedp[querypos] == true) {
      printf("plus %d (%d):",querypos,this->plus_npositions[querypos]);
      for (i = 0; i < this->plus_npositions[querypos]; i++) {
	printf(" %u",this->plus_positions[querypos][i]);
      }
      printf("\n");
    }
    if (this->minus_retrievedp[querypos] == true) {
      printf("minus %d (%d):",querypos,this->minus_npositions[querypos]);
      for (i = 0; i < this->minus_npositions[querypos]; i++) {
	printf(" %u",this->minus_positions[querypos][i]);
      }
      printf("\n");
    }
  }
  printf("\n");
  return;
}
#endif


static void
Stage1_init (T this, char *queryuc_ptr, int querylength, int genestrand) {
  Reader_T reader;
  Oligostate_T last_state;
  Oligospace_T forward, revcomp, forward_oligo, revcomp_oligo;
  int querypos, query_lastpos;
  int mod;


  debug1(printf("%s\n",queryuc_ptr));
  reader = Reader_new(queryuc_ptr,/*querystart*/0,/*queryend*/querylength);
  last_state = INIT;
  forward = revcomp = 0;
  mod = 0;

  while (mod < index1interval &&
	 (last_state = Oligo_next_5(last_state,&querypos,&forward,&revcomp,reader,genestrand)) != DONE) {
    while (mod < index1interval && mod < querypos) {
      debug1(printf("Skipping rawpositions_5 %d, because querypos is not the expected one\n",mod));
#ifdef LARGE_GENOMES
      this->plus_rawpositions_high_5[mod] = (unsigned char *) NULL;
      this->minus_rawpositions_high_5[mod] = (unsigned char *) NULL;
#endif
      this->plus_rawpositions_5[mod] = (UINT4 *) NULL;
      this->minus_rawpositions_5[mod] = (UINT4 *) NULL;
      this->plus_nrawpositions_5[mod] = 0;
      this->minus_nrawpositions_5[mod] = 0;
      mod++;
    }

    if (mod < index1interval) {
      forward_oligo = forward & oligobase_mask;
      this->plus_diagterms_5[mod] = -querypos;
#ifdef LARGE_GENOMES      
      this->plus_nrawpositions_5[mod] = 
	Indexdb_largeptr_with_diagterm(&(this->plus_rawpositions_high_5[mod]),&(this->plus_rawpositions_5[mod]),
				       /*plus_indexdb*/indexdb_fwd,forward_oligo,this->plus_diagterms_5[mod]);
#else
      this->plus_nrawpositions_5[mod] =
	Indexdb_ptr_with_diagterm(&(this->plus_rawpositions_5[mod]),/*plus_indexdb*/indexdb_fwd,forward_oligo,
				  this->plus_diagterms_5[mod]);
#endif
      debug1(printf("(1) plus_nrawpositions_5[%d] = %d, oligo %016lX\n",mod,this->plus_nrawpositions_5[mod],forward_oligo));
      
      revcomp_oligo = (revcomp >> leftreadshift) & oligobase_mask;
      this->minus_diagterms_5[mod] = -querylength + querypos + index1part;
#ifdef LARGE_GENOMES
      this->minus_nrawpositions_5[mod] =
	Indexdb_largeptr_with_diagterm(&(this->minus_rawpositions_high_5[mod]),&(this->minus_rawpositions_5[mod]),
				       /*minus_indexdb*/indexdb_rev,revcomp_oligo,this->minus_diagterms_5[mod]);
#else
      this->minus_nrawpositions_5[mod] =
	Indexdb_ptr_with_diagterm(&(this->minus_rawpositions_5[mod]),/*minus_indexdb*/indexdb_rev,revcomp_oligo,
				  this->minus_diagterms_5[mod]);
#endif
      debug1(printf("(2) minus_nrawpositions_5[%d] = %d, oligo %016lX\n",mod,this->minus_nrawpositions_5[mod],revcomp_oligo));
      
      debug1(printf("5' end: %s %s: %d plus positions, %d minus positions, genestrand %d\n",
		    Oligo_one_nt(forward_oligo,index1part),Oligo_one_nt(revcomp_oligo,index1part),
		    this->plus_nrawpositions_5[mod],this->minus_nrawpositions_5[mod],genestrand));
    }

    mod++;
  }

  while (mod < index1interval) {
    debug1(printf("Skipping rawpositions_5 %d, because last_state was DONE\n",mod));
#ifdef LARGE_GENOMES
    this->plus_rawpositions_high_5[mod] = (unsigned char *) NULL;
    this->minus_rawpositions_high_5[mod] = (unsigned char *) NULL;
#endif
    this->plus_rawpositions_5[mod] = (UINT4 *) NULL;
    this->minus_rawpositions_5[mod] = (UINT4 *) NULL;
    this->plus_nrawpositions_5[mod] = 0;
    this->minus_nrawpositions_5[mod] = 0;
    mod++;
  }
  Reader_free(&reader);


  query_lastpos = querylength - index1part;
  reader = Reader_new(queryuc_ptr,/*querystart*/0,/*queryend*/querylength);
  last_state = INIT;
  forward = revcomp = 0;
  mod = 0;
  while (mod < index1interval &&
	 (last_state = Oligo_next_3(last_state,&querypos,&forward,&revcomp,reader,genestrand)) != DONE) {
    while (mod < index1interval && mod < query_lastpos - querypos) {
      debug1(printf("Skipping rawpositions_3 %d, because querypos is not the expected one\n",mod));
#ifdef LARGE_GENOMES
      this->plus_rawpositions_high_3[mod] = (unsigned char *) NULL;
      this->minus_rawpositions_high_3[mod] = (unsigned char *) NULL;
#endif
      this->plus_rawpositions_3[mod] = (UINT4 *) NULL;
      this->minus_rawpositions_3[mod] = (UINT4 *) NULL;
      this->plus_nrawpositions_3[mod] = 0;
      this->minus_nrawpositions_3[mod] = 0;
      mod++;
    }

    if (mod < index1interval) {
      forward_oligo = (forward >> leftreadshift) & oligobase_mask;
      this->plus_diagterms_3[mod] = -querypos;
#ifdef LARGE_GENOMES
      this->plus_nrawpositions_3[mod] =
	Indexdb_largeptr_with_diagterm(&(this->plus_rawpositions_high_3[mod]),&(this->plus_rawpositions_3[mod]),
				       /*plus_indexdb*/indexdb_fwd,forward_oligo,this->plus_diagterms_3[mod]);
#else
      this->plus_nrawpositions_3[mod] =
	Indexdb_ptr_with_diagterm(&(this->plus_rawpositions_3[mod]),/*plus_indexdb*/indexdb_fwd,forward_oligo,
				  this->plus_diagterms_3[mod]);
#endif
      debug1(printf("(3) plus_nrawpositions_3[%d] = %d, oligo %016lX\n",mod,this->plus_nrawpositions_3[mod],forward_oligo));
      
      revcomp_oligo = revcomp & oligobase_mask;
      this->minus_diagterms_3[mod] = -querylength + querypos + index1part;
#ifdef LARGE_GENOMES
      this->minus_nrawpositions_3[mod] =
	Indexdb_largeptr_with_diagterm(&(this->minus_rawpositions_high_3[mod]),&(this->minus_rawpositions_3[mod]),
				       /*minus_indexdb*/indexdb_rev,revcomp_oligo,this->minus_diagterms_3[mod]);
#else
      this->minus_nrawpositions_3[mod] =
	Indexdb_ptr_with_diagterm(&(this->minus_rawpositions_3[mod]),/*minus_indexdb*/indexdb_rev,revcomp_oligo,
				  this->minus_diagterms_3[mod]);
#endif
      debug1(printf("(4) minus_nrawpositions_3[%d] = %d, oligo %016lX\n",mod,this->minus_nrawpositions_3[mod],revcomp_oligo));

      debug1(printf("3' end: %s %s: %d plus positions, %d minus positions, genestrand %d\n",
		    Oligo_one_nt(forward_oligo,index1part),Oligo_one_nt(revcomp_oligo,index1part),
		    this->plus_nrawpositions_3[mod],this->minus_nrawpositions_3[mod],genestrand));
    }

    mod++;
  }

  while (mod < index1interval) {
    /* printf("Skipping rawpositions_3 %d, because last_state was DONE\n",mod); */
#ifdef LARGE_GENOMES
    this->plus_rawpositions_high_3[mod] = (unsigned char *) NULL;
    this->minus_rawpositions_high_3[mod] = (unsigned char *) NULL;
#endif
    this->plus_rawpositions_3[mod] = (UINT4 *) NULL;
    this->minus_rawpositions_3[mod] = (UINT4 *) NULL;
    this->plus_nrawpositions_3[mod] = 0;
    this->minus_nrawpositions_3[mod] = 0;
    mod++;
  }
  Reader_free(&reader);

  return;
}




/* onep indicates whether Kmer_search_genome_one_end was run, which
   fills beyond the first index1interval */
static void
Stage1_rearrange (T this, int querylength, bool onep) {
  int query_lastpos = querylength - index1part;
  int mod;

  for (mod = 0; mod < index1interval; mod++) {
    this->plus_retrievedp[mod] = true;
#ifdef LARGE_GENOMES
    this->plus_positions_high[mod] = this->plus_rawpositions_high_5[mod];
#endif
    this->plus_positions[mod] = this->plus_rawpositions_5[mod];
    this->plus_npositions[mod] = this->plus_nrawpositions_5[mod];

    this->plus_retrievedp[query_lastpos-mod] = true;
#ifdef LARGE_GENOMES
    this->plus_positions_high[query_lastpos-mod] = this->plus_rawpositions_high_3[mod];
#endif
    this->plus_positions[query_lastpos-mod] = this->plus_rawpositions_3[mod];
    this->plus_npositions[query_lastpos-mod] = this->plus_nrawpositions_3[mod];

    /* Using new sarray and segment-based conventions */
    this->minus_retrievedp[query_lastpos-mod] = true;
#ifdef LARGE_GENOMES
    this->minus_positions_high[query_lastpos-mod] = this->minus_rawpositions_high_5[mod];
#endif
    this->minus_positions[query_lastpos-mod] = this->minus_rawpositions_5[mod];
    this->minus_npositions[query_lastpos-mod] = this->minus_nrawpositions_5[mod];

    this->minus_retrievedp[mod] = true;
#ifdef LARGE_GENOMES
    this->minus_positions_high[mod] = this->minus_rawpositions_high_3[mod];
#endif
    this->minus_positions[mod] = this->minus_rawpositions_3[mod];
    this->minus_npositions[mod] = this->minus_nrawpositions_3[mod];

#if 0
    printf("Initializing plus_positions[%d] to be %p, with count of %d\n",
	   mod,this->plus_positions[mod],this->plus_npositions[mod]);
    printf("Initializing plus_positions[%d] to be %p, with count of %d\n",
	 query_lastpos-mod,this->plus_positions[query_lastpos-mod],this->plus_npositions[query_lastpos-mod]);

    printf("Initializing minus_positions[%d] to be %p, with count of %d\n",
	   mod,this->minus_positions[mod],this->minus_npositions[mod]);
    printf("Initializing minus_positions[%d] to be %p, with count of %d\n",
	 query_lastpos-mod,this->minus_positions[query_lastpos-mod],this->minus_npositions[query_lastpos-mod]);
#endif
  }

  if (onep == true) {
    /* Not currently doing one_kmer */
    for (mod = 0; mod < index1interval; mod++) {
      this->plus_retrievedp[index1part+mod] = true;
#ifdef LARGE_GENOMES
      this->plus_positions_high[index1part+mod] = this->plus_rawpositions_high_5[index1interval+mod];
#endif
      this->plus_positions[index1part+mod] = this->plus_rawpositions_5[index1interval+mod];
      this->plus_npositions[index1part+mod] = this->plus_nrawpositions_5[index1interval+mod];

      this->plus_retrievedp[query_lastpos-index1part-mod] = true;
#ifdef LARGE_GENOMES
      this->plus_positions_high[query_lastpos-index1part-mod] = this->plus_rawpositions_high_3[index1interval+mod];
#endif
      this->plus_positions[query_lastpos-index1part-mod] = this->plus_rawpositions_3[index1interval+mod];
      this->plus_npositions[query_lastpos-index1part-mod] = this->plus_nrawpositions_3[index1interval+mod];

      /* Using new sarray and segment-based conventions */
      this->minus_retrievedp[query_lastpos-index1part-mod] = true;
#ifdef LARGE_GENOMES
      this->minus_positions_high[query_lastpos-index1part-mod] = this->minus_rawpositions_high_5[index1interval+mod];
#endif
      this->minus_positions[query_lastpos-index1part-mod] = this->minus_rawpositions_5[index1interval+mod];
      this->minus_npositions[query_lastpos-index1part-mod] = this->minus_nrawpositions_5[index1interval+mod];

      this->minus_retrievedp[index1part+mod] = true;
#ifdef LARGE_GENOMES
      this->minus_positions_high[index1part+mod] = this->minus_rawpositions_high_3[index1interval+mod];
#endif
      this->minus_positions[index1part+mod] = this->minus_rawpositions_3[index1interval+mod];
      this->minus_npositions[index1part+mod] = this->minus_nrawpositions_3[index1interval+mod];

#if 0
      printf("Initializing plus_positions[%d] to be %p, with count of %d\n",
	     index1part+mod,this->plus_positions[index1part+mod],this->plus_npositions[index1part+mod]);
      printf("Initializing plus_positions[%d] to be %p, with count of %d\n",
	     query_lastpos-index1part-mod,this->plus_positions[query_lastpos-index1part-mod],
	     this->plus_npositions[query_lastpos-index1part-mod]);

      printf("Initializing minus_positions[%d] to be %p, with count of %d\n",
	     index1part+mod,this->minus_positions[index1part+mod],this->minus_npositions[index1part+mod]);
      printf("Initializing minus_positions[%d] to be %p, with count of %d\n",
	     query_lastpos-index1part-mod,this->minus_positions[query_lastpos-index1part-mod],
	     this->minus_npositions[query_lastpos-index1part-mod]);
#endif
    }
  }

  /* Stage1_dump(this,querylength); */

  return;
}


void
Stage1_fill_all_oligos (T this, char *queryuc_ptr, int querylength, int genestrand) {
  Reader_T reader;
  int query_lastpos, querypos, querypos_rc;
  Oligostate_T last_state = INIT;
  Oligospace_T forward = 0, revcomp = 0, oligo;
#ifdef HAVE_64_BIT
  Uint8table_T plus_seenp, minus_seenp;
#else
  Uinttable_T plus_seenp, minus_seenp;
#endif

#ifdef HAVE_64_BIT
  plus_seenp = Uint8table_new(/*hint*/querylength,/*save_contents_p*/false);
  minus_seenp = Uint8table_new(/*hint*/querylength,/*save_contents_p*/false);
#else
  plus_seenp = Uinttable_new(/*hint*/querylength,/*save_contents_p*/false);
  minus_seenp = Uinttable_new(/*hint*/querylength,/*save_contents_p*/false);
#endif

  query_lastpos = querylength - index1part;
  reader = Reader_new(queryuc_ptr,/*querystart*/0,/*queryend*/querylength);

  /* Note: leftshifting is done here, rather than in Oligo_lookup */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  /* debug(printf("oligobase_mask: %010llX\n",oligobase_mask)); */
  querypos = 0;
  while ((last_state = Oligo_next_5(last_state,&querypos,&forward,&revcomp,reader,genestrand)) != DONE) {
    if (last_state != VALID) {
      /* querypos is not defined when last_state != VALID */
      debug8(printf("oligo at plus %d, minus %d is not valid\n",querypos,querypos_rc));
    } else {
      querypos_rc = query_lastpos - querypos;
      oligo = this->plus_oligos[querypos] = forward & oligobase_mask;
      debug8(printf("Putting oligo %016lX at plus %d\n",oligo,querypos));
#ifdef HAVE_64_BIT
      if (Uint8table_get(plus_seenp,oligo) != NULL) {
	debug8(printf("oligo at plus %d already seen, so marking as invalid\n",querypos));
	this->plus_validp[querypos] = false;
      } else {
	this->plus_validp[querypos] = true;
	Uint8table_put(plus_seenp,oligo,(void *) true);
      }
#else
      if (Uinttable_get(plus_seenp,oligo) != NULL) {
	debug8(printf("oligo at plus %d already seen, so marking as invalid\n",querypos));
	this->plus_validp[querypos] = false;
      } else {
	this->plus_validp[querypos] = true;
	Uinttable_put(plus_seenp,oligo,(void *) true);
      }
#endif

      oligo = this->minus_oligos[querypos_rc] = (revcomp >> leftreadshift) & oligobase_mask;
      debug8(printf("Putting oligo %016lX at minus %d\n",oligo,querypos_rc));
#ifdef HAVE_64_BIT
      if (Uint8table_get(minus_seenp,oligo) != NULL) {
	debug8(printf("oligo at minus %d already seen, so marking as invalid\n",querypos));
	this->minus_validp[querypos_rc] = false;
      } else {
	this->minus_validp[querypos_rc] = true;
	Uint8table_put(minus_seenp,oligo,(void *) true);
      }
#else
      if (Uinttable_get(minus_seenp,oligo) != NULL) {
	debug8(printf("oligo at minus %d already seen, so marking as invalid\n",querypos));
	this->minus_validp[querypos_rc] = false;
      } else {
	this->minus_validp[querypos_rc] = true;
	Uinttable_put(minus_seenp,oligo,(void *) true);
      }
#endif
    }
  }

  Reader_free(&reader);

#ifdef HAVE_64_BIT
  Uint8table_free(&minus_seenp);
  Uint8table_free(&plus_seenp);
#else
  Uinttable_free(&minus_seenp);
  Uinttable_free(&plus_seenp);
#endif

  return;
}


void
Stage1_fill_position_ptrs_plus (T this, int querystart, int queryend, int genestrand) {
  int query_lastpos, querypos;
  Indexdb_T plus_indexdb;

  if (genestrand == +2) {
    plus_indexdb = indexdb_rev;
  } else {
    plus_indexdb = indexdb_fwd;
  }

  query_lastpos = queryend - index1part;

  /* Assumes that plus_oligos have been filled in (by Stage1_fill_all_oligos */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  /* debug(printf("oligobase_mask: %010llX\n",oligobase_mask)); */
  for (querypos = querystart; querypos <= query_lastpos; querypos++) {
    if (this->plus_retrievedp[querypos] == true) {
      /* No need to do anything */
    } else if (this->plus_validp[querypos] == false) {
      /* Oligo not valid */
    } else {
#ifdef LARGE_GENOMES
      this->plus_npositions[querypos] = 
	Indexdb_largeptr_with_diagterm(&this->plus_positions_high[querypos],&this->plus_positions[querypos],
				       plus_indexdb,this->plus_oligos[querypos],/*diagterm*/-querypos);
#else
      this->plus_npositions[querypos] = 
	Indexdb_ptr_with_diagterm(&this->plus_positions[querypos],plus_indexdb,
				  this->plus_oligos[querypos],/*diagterm*/-querypos);
#endif
      this->plus_retrievedp[querypos] = true;
    }
  }

  return;
}

void
Stage1_fill_position_ptrs_minus (T this, int querystart, int queryend, int genestrand) {
  int query_lastpos, querypos;
  Indexdb_T minus_indexdb;

  if (genestrand == +2) {
    minus_indexdb = indexdb_fwd;
  } else {
    minus_indexdb = indexdb_rev;
  }

  query_lastpos = queryend - index1part;

  /* Assumes that minus_oligos have been filled in (by Stage1_fill_all_oligos */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  /* debug(printf("oligobase_mask: %010llX\n",oligobase_mask)); */
  for (querypos = querystart; querypos <= query_lastpos; querypos++) {
    if (this->minus_retrievedp[querypos] == true) {
      /* No need to do anything */
    } else if (this->minus_validp[querypos] == false) {
      /* Oligo not valid */
    } else {
#ifdef LARGE_GENOMES
      this->minus_npositions[querypos] =
	Indexdb_largeptr_with_diagterm(&this->minus_positions_high[querypos],&this->minus_positions[querypos],
				       minus_indexdb,this->minus_oligos[querypos],/*diagterm*/-querypos);
#else
      this->minus_npositions[querypos] =
	Indexdb_ptr_with_diagterm(&this->minus_positions[querypos],minus_indexdb,
				  this->minus_oligos[querypos],/*diagterm*/-querypos);
#endif
      this->minus_retrievedp[querypos] = true;
    }
  }

  return;
}


#if 0
static void
Stage1_ptr_all_positions (T this, int querylength, int genestrand) {
  int query_lastpos, querypos, querypos_rc;
  Indexdb_T plus_indexdb, minus_indexdb;

  if (genestrand == +2) {
    plus_indexdb = indexdb_rev;
    minus_indexdb = indexdb_fwd;
  } else {
    plus_indexdb = indexdb_fwd;
    minus_indexdb = indexdb_rev;
  }

  query_lastpos = querylength - index1part;

  /* Note: leftshifting is done here, rather than in Oligo_lookup */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  /* debug(printf("oligobase_mask: %010llX\n",oligobase_mask)); */
  for (querypos = 0, querypos_rc = query_lastpos; querypos <= query_lastpos; querypos++, --querypos_rc) {
    if (this->plus_retrievedp[querypos] == true) {
      /* No need to do anything */
    } else if (this->plus_validp[querypos] == false) {
      /* Oligo not valid */
    } else {
#ifdef LARGE_GENOMES
      this->plus_npositions[querypos] = 
	Indexdb_largeptr_with_diagterm(&this->plus_positions_high[querypos],&this->plus_positions[querypos],
				       plus_indexdb,this->plus_oligos[querypos],/*diagterm*/-querypos);
#else
      this->plus_npositions[querypos] = 
	Indexdb_ptr_with_diagterm(&this->plus_positions[querypos],plus_indexdb,
				  this->plus_oligos[querypos],/*diagterm*/-querypos);
#endif
      this->plus_retrievedp[querypos] = true;
    }

    if (this->minus_retrievedp[querypos_rc] == true) {
      /* No need to do anything */
    } else if (this->minus_validp[querypos_rc] == false) {
      /* Oligo not valid */
    } else {
#ifdef LARGE_GENOMES
      this->minus_npositions[querypos_rc] =
	Indexdb_largeptr_with_diagterm(&this->minus_positions_high[querypos_rc],&this->minus_positions[querypos_rc],
				       minus_indexdb,this->minus_oligos[querypos_rc],/*diagterm*/-querypos_rc);
#else
      this->minus_npositions[querypos_rc] =
	Indexdb_ptr_with_diagterm(&this->minus_positions[querypos_rc],minus_indexdb,
				  this->minus_oligos[querypos_rc],/*diagterm*/-querypos_rc);
#endif
      this->minus_retrievedp[querypos_rc] = true;
    }
  }

  return;
}
#endif


static int
determine_sizelimit (T this, int querylength) {
  int cutoff, *set, count;
  int n;
  int query_lastpos, querypos, querypos_rc;

  assert(querylength >= index1part);

  query_lastpos = querylength - index1part;
  set = (int *) MALLOC(2*(query_lastpos+1)*sizeof(int));
  n = 0;
  for (querypos = 0, querypos_rc = query_lastpos; querypos <= query_lastpos; querypos++, --querypos_rc) {
    if (this->plus_validp[querypos] == true) {
      set[n++] = count = this->plus_npositions[querypos];
    }

    if (this->minus_validp[querypos_rc] == true) {
      set[n++] = count = this->minus_npositions[querypos_rc];
    }
  }

  if (n < 5) {
    cutoff = MIN_SIZELIMIT;
  } else if ((cutoff = Orderstat_int_pct_inplace(set,n,/*pct*/0.60)) < MIN_SIZELIMIT) {
    cutoff = MIN_SIZELIMIT;
  }
  FREE(set);

  return cutoff;
}





void
Stage1_free (T *old) {

  /* Stage1hr_check(*old); */

  if (*old) {
#if 0
    /* Now pointing to data structure, and not copying values */
    for (i = 0; i <= query_lastpos; i++) {
      if ((*old)->plus_retrievedp[i] == true) {
#ifdef LARGE_GENOMES
	FREE((*old)->plus_positions_high[i]);
#endif
	FREE((*old)->plus_positions[i]);
      }

      if ((*old)->minus_retrievedp[i] == true) {
#ifdef LARGE_GENOMES
	FREE((*old)->minus_positions_high[i]);
#endif
	FREE((*old)->minus_positions[i]);
      }
    }
#endif

    FREE((*old)->mismatch_positions_alloc);
    FREE((*old)->positions_alloc);
    Spliceinfo_free(&(*old)->spliceinfo);

#ifdef LARGE_GENOMES
    FREE((*old)->stream_high_alloc);
    FREE((*old)->stream_low_alloc);
#endif
    FREE((*old)->stream_alloc);
    FREE((*old)->streamsize_alloc);
    FREE((*old)->querypos_diagterm_alloc);


#ifdef LARGE_GENOMES
    FREE((*old)->positions_high_allocated);
#endif
    FREE((*old)->positions_allocated);
    FREE((*old)->npositions_allocated);

    FREE((*old)->retrievedp_allocated);

    FREE((*old)->minus_oligos);
    FREE((*old)->plus_oligos);

    FREE((*old)->minus_validp);
    FREE((*old)->plus_validp);

    Elt_gc(&(*old)->queryfwd_plus_set);
    Elt_gc(&(*old)->queryfwd_minus_set);
    Elt_gc(&(*old)->queryrev_plus_set);
    Elt_gc(&(*old)->queryrev_minus_set);

    FREE(*old);
  }

  return;
}


/************************************************************************/

static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}


#if 0
static void
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}
#endif

/************************************************************************/



static void
stage3list_gc (List_T *old) {
  List_T p;
  Stage3end_T hit;

  for (p = *old; p != NULL; p = p->rest) {
    hit = (Stage3end_T) p->first;
    Stage3end_free(&hit);
  }
  Hitlist_free(&(*old));

  return;
}


static void
remap_to_transcriptome (List_T hits) {
  List_T p;
  Stage3end_T hit;
  List_T transcripts;

  char *remap_sequence;
  int remap_seqlength;

  for (p = hits; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (Stage3end_transcripts(hit) != NULL) {
      /* Already mapped to transcriptome */
    } else {
      remap_sequence = Stage3end_substrings_genomic_sequence(&remap_seqlength,hit,genomecomp);
      if ((transcripts = Kmer_remap_transcriptome(remap_sequence,remap_seqlength,
						  Stage3end_chrnum(hit),Stage3end_chrpos_low(hit),Stage3end_chrpos_high(hit),
						  transcript_iit,transcriptomebits,transcriptome)) != NULL) {
	/* Following procedure frees the list hit->transcripts */
	Stage3end_set_transcripts(hit,transcripts);
      }
      FREE(remap_sequence);
    }
  }

  return;
}



static int
single_read (int *found_score_overall, int *found_score_within_trims,
	     List_T *hits_gplus, List_T *hits_gminus,
	     T this, char *queryuc_ptr, char *queryrc, int querylength,
	     Compress_T query_compress_fwd, Compress_T query_compress_rev, int genestrand,
#if 0
	     Oligoindex_array_T oligoindices_minor,
	     Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
	     Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
#endif
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
	     bool paired_end_p, bool first_read_p) {
  int sizelimit;

  int done_level, nmismatches_allowed;
  /* int max_terminal_mismatches; */
  int max_splice_mismatches;
  int min_trim;

  bool abort_exact_p = false, elt_set_extended_p = false;
  Stage3end_T hit;

  List_T p;
  bool onep = false;


  /* For Segment_search */
  struct Record_T *plus_records = NULL, *minus_records = NULL;
  int plus_nrecords = 0, minus_nrecords = 0;

  List_T startfrags_plus = NULL, endfrags_plus = NULL,
    startfrags_minus = NULL, endfrags_minus = NULL;


  debug(printf("Entered single_read with queryuc_ptr %s\n",queryuc_ptr));
  *hits_gplus = *hits_gminus = (List_T) NULL;

#if 0
  /* Use user_maxlevel_float for final filtering, not for search */
  if (user_maxlevel_float < 0.0) {
    nmismatches_allowed = querylength/index1part;
  } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
    nmismatches_allowed = (int) rint(user_maxlevel_float * (double) querylength);
  } else {
    nmismatches_allowed = (int) user_maxlevel_float;
  }
#else
  nmismatches_allowed = querylength/index1part;
#endif

  *found_score_overall = *found_score_within_trims = querylength;

  if (querylength < index1part + index1interval - 1) {
    return /*level*/0;
  }

  if (indexdb_tr != NULL) {
    debug(printf("Running Kmer_search_transcriptome_single\n"));
    Kmer_search_transcriptome_single(&(*found_score_overall),&(*found_score_within_trims),
				     &(*hits_gplus),&(*hits_gminus),queryuc_ptr,querylength,
				     this->tplus_stream_array,this->tplus_streamsize_array,this->tplus_diagterm_array,
				     this->tminus_stream_array,this->tminus_streamsize_array,this->tminus_diagterm_array,
				     query_compress_fwd,query_compress_rev,
				     transcript_iit,transcriptome,transcriptomebits,
				     nmismatches_allowed,listpool,hitlistpool,/*level*/0);

  }

  if (use_only_transcriptome_p == true || *found_score_within_trims <= nmismatches_allowed) {
    return TR;
  } else {
    debug(printf("LEVEL %d: RUNNING KMER_SEARCH_GENOME_ENDS\n",1));
    Stage1_init(this,queryuc_ptr,querylength,genestrand);

    Kmer_search_genome_ends_exact(&abort_exact_p,&(*found_score_overall),&(*found_score_within_trims),
				  &(*hits_gplus),&(*hits_gminus),this,querylength,this->mismatch_positions_alloc,
				  query_compress_fwd,query_compress_rev,genestrand,nmismatches_allowed,
				  listpool,hitlistpool,/*level*/1);
    
    debug(printf("After Kmer_search_genome_ends, we have found_score %d (vs allowed %d), %d plus and %d minus hits\n",
		 *found_score_within_trims,nmismatches_allowed,List_length(*hits_gplus),List_length(*hits_gminus)));
  }
  
  if (*found_score_within_trims <= nmismatches_allowed) {
    if (indexdb_tr != NULL && remap_transcriptome_p == true) {
      remap_to_transcriptome(*hits_gplus);
      remap_to_transcriptome(*hits_gminus);
    }
    return KMER_EXACT;

  } else {
    Stage1_rearrange(this,querylength,onep);
    Stage1_fill_all_oligos(this,queryuc_ptr,querylength,genestrand);
    debug(printf("LEVEL %d: RUNNING EXTENSION_SEARCH\n",2));
    Extension_search(&(*found_score_overall),&(*found_score_within_trims),&(*hits_gplus),&(*hits_gminus),
		     this,queryuc_ptr,queryrc,querylength,query_compress_fwd,query_compress_rev,
		     /*pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,*/
		     intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
		     nmismatches_allowed,genestrand,paired_end_p,first_read_p,/*level*/2);
    debug(printf("After Extension search, we have found_score_within_trims %d (vs allowed %d), %d plus and %d minus hits\n",
		 *found_score_within_trims,nmismatches_allowed,List_length(*hits_gplus),List_length(*hits_gminus)));
  }

  /* Need to use found_score_overall here and afterwards in order to find fusions */
  debug(printf("Comparing found_score_overall %d with nmismatches_allowed %d\n",
	       *found_score_overall,nmismatches_allowed));
  if (abort_exact_p == true || *found_score_overall <= nmismatches_allowed) {
    if (indexdb_tr != NULL && remap_transcriptome_p == true) {
      remap_to_transcriptome(*hits_gplus);
      remap_to_transcriptome(*hits_gminus);
    }
    return EXT;

  } else {
    debug(printf("LEVEL %d: RUNNING KMER_SEARCH_GENOME_ENDS_APPROX\n",3));
    Kmer_search_genome_ends_approx(&(*found_score_overall),&(*found_score_within_trims),
				   &(*hits_gplus),&(*hits_gminus),this,
				   query_compress_fwd,query_compress_rev,querylength,
				   genestrand,nmismatches_allowed,
				   /*sizelimit*/3000,listpool,hitlistpool,/*level*/3);
  }

  debug(printf("found_score_overall %d\n",*found_score_overall));
  if (*found_score_overall <= nmismatches_allowed) {
    if (indexdb_tr != NULL && remap_transcriptome_p == true) {
      debug(printf("Remapping to transcriptome\n"));
      remap_to_transcriptome(*hits_gplus);
      remap_to_transcriptome(*hits_gminus);
      debug(printf("Done with remapping to transcriptome\n"));
    }
    return KMER_APPROX;

  } else {
    if (querylength >= index1part) {
      Stage1_fill_position_ptrs_plus(this,/*querystart*/0,/*queryend*/querylength,genestrand);
      Stage1_fill_position_ptrs_minus(this,/*querystart*/0,/*queryend*/querylength,genestrand);
      /* Need sizelimit to constrain segment search */
      sizelimit = determine_sizelimit(this,querylength);
      
      debug(printf("LEVEL %d: RUNNING SEGMENT_SEARCH\n",4));
      debug(printf("Starting Segment_identify on plus strand\n"));
      plus_records = Segment_identify(&plus_nrecords,
#ifdef LARGE_GENOMES
				      this->plus_positions_high,
#endif
				      this->plus_positions,this->plus_npositions,this->plus_validp,
#ifdef LARGE_GENOMES
				      this->stream_high_alloc,this->stream_low_alloc,
#else
				      this->stream_alloc,
#endif
				      this->streamsize_alloc,this->querypos_diagterm_alloc,
				      querylength,sizelimit);
      debug(printf("Done\n"));

      debug(printf("Starting Segment_identify on minus strand\n"));
      minus_records = Segment_identify(&minus_nrecords,
#ifdef LARGE_GENOMES
				       this->minus_positions_high,
#endif
				       this->minus_positions,this->minus_npositions,this->minus_validp,
#ifdef LARGE_GENOMES
				       this->stream_high_alloc,this->stream_low_alloc,
#else
				       this->stream_alloc,
#endif
				       this->streamsize_alloc,this->querypos_diagterm_alloc,
				       querylength,sizelimit);
      debug(printf("Done\n"));
      
      debug(printf("Starting Segment_search_all\n"));
      Segment_search_all(&(*found_score_overall),&(*found_score_within_trims),&(*hits_gplus),&(*hits_gminus),
			 plus_records,plus_nrecords,minus_records,minus_nrecords,
			 
			 queryuc_ptr,queryrc,querylength,this->mismatch_positions_alloc,
			 this->spliceinfo,this->stream_alloc,this->streamsize_alloc,
			 query_compress_fwd,query_compress_rev,genestrand,paired_end_p,first_read_p,
			 intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
			 /*method*/SEGMENT1,/*level*/4);
      debug(printf("Done\n"));

      FREE(minus_records);
      FREE(plus_records);
    }
  }


  debug(printf("Comparing found_score_overall %d with nmismatches_allowed %d\n",
	       *found_score_overall,nmismatches_allowed));
  if (*found_score_overall <= nmismatches_allowed) {
    if (indexdb_tr != NULL && remap_transcriptome_p == true) {
      remap_to_transcriptome(*hits_gplus);
      remap_to_transcriptome(*hits_gminus);
    }
    return SEGMENT1;
    
  } else {
    /* Distant RNA splicing.  No longer remap to transcriptome */
    if ((done_level = (*found_score_overall) + subopt_levels) > nmismatches_allowed) {
      done_level = nmismatches_allowed;
    }
    debug(printf("done_level %d = found_score %d + subopt_levels %d\n",done_level,*found_score_within_trims,subopt_levels));
    
    min_trim = querylength;	/* min of max(trim5,trim3) over each hit */
    for (p = *hits_gplus; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      if (Stage3end_max_trim(hit) < min_trim) {
	min_trim = Stage3end_max_trim(hit);
      }
    }
    for (p = *hits_gminus; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      if (Stage3end_max_trim(hit) < min_trim) {
	min_trim = Stage3end_max_trim(hit);
      }
    }
    
    debug(printf("Comparing min_trim %d with min_distantsplicing_end_matches %d\n",
		 min_trim,min_distantsplicing_end_matches));
    if (min_trim < min_distantsplicing_end_matches) {
      /* Don't find distant splicing */
      debug(printf("Skipping distant splicing because min_trim %d < min_distantsplicing_end_matches %d\n",
		   min_trim,min_distantsplicing_end_matches));
      
    } else if ((max_splice_mismatches = done_level - distantsplicing_penalty) < 0) {
      debug(printf("Skipping distant splicing because done_level %d - distantsplicing_penalty %d < 0\n",
		   done_level,distantsplicing_penalty));

    } else if (novelsplicingp == false) {
      /* TODO: Implement distant DNA fusions */
      debug(printf("Skipping distant splicing novelsplicingp is false\n"));

    } else {
      debug(printf("Candidate for distant splicing because min_trim %d >= min_distantsplicing_end_matches %d\n",
		   min_trim,min_distantsplicing_end_matches));
      debug(printf("Have sets: %d %d %d %d\n",
		   List_length(this->queryfwd_plus_set),List_length(this->queryfwd_minus_set),
		   List_length(this->queryrev_plus_set),List_length(this->queryrev_minus_set)));
      debug(printf("LEVEL %d: RUNNING DISTANT_RNA_SOLVE\n",5));

#if 0
      Elt_set_queryfwd_extend_left(this->queryfwd_plus_set,/*query_compress*/query_compress_fwd,
				   /*plusp*/true,genestrand);
      Elt_set_queryfwd_extend_left(this->queryfwd_minus_set,/*query_compress*/query_compress_rev,
				   /*plusp*/false,genestrand);
      Elt_set_queryrev_extend_right(this->queryrev_plus_set,querylength,/*query_compress*/query_compress_fwd,
				    /*plusp*/true,genestrand);
      Elt_set_queryrev_extend_right(this->queryrev_minus_set,querylength,/*query_compress*/query_compress_rev,
				    /*plusp*/false,genestrand);
      elt_set_extended_p = true;
#endif

      Distant_rna_solve(&(*found_score_overall),&(*found_score_within_trims),&(*hits_gplus),&(*hits_gminus),
			&startfrags_plus,&endfrags_plus,&startfrags_minus,&endfrags_minus,

			this->queryfwd_plus_set,this->queryfwd_minus_set,
			this->queryrev_plus_set,this->queryrev_minus_set,
			
			this->mismatch_positions_alloc,this->positions_alloc,
			query_compress_fwd,query_compress_rev,queryuc_ptr,queryrc,querylength,
			max_splice_mismatches,genestrand,/*first_read_p*/true,listpool,hitlistpool,
			/*level*/5);
    }
  }
  

  if (*found_score_overall <= nmismatches_allowed) {
    Substring_list_gc(&startfrags_plus);
    Substring_list_gc(&endfrags_plus);
    Substring_list_gc(&startfrags_minus);
    Substring_list_gc(&endfrags_minus);
    return DISTANT_RNA;

  } else if (find_dna_chimeras_p == true) {
    debug(printf("LEVEL %d: RUNNING DISTANT_DNA_SOLVE\n",6));
    startfrags_plus = Substring_sort_siteN_halves(startfrags_plus,listpool,/*ascendingp*/true);
    endfrags_plus = Substring_sort_siteN_halves(endfrags_plus,listpool,/*ascendingp*/true);
    startfrags_minus = Substring_sort_siteN_halves(startfrags_minus,listpool,/*ascendingp*/true);
    endfrags_minus = Substring_sort_siteN_halves(endfrags_minus,listpool,/*ascendingp*/true);

    Distant_dna_solve(&(*found_score_overall),&(*found_score_within_trims),&(*hits_gplus),&(*hits_gminus),
		      startfrags_plus,endfrags_plus,startfrags_minus,endfrags_minus,
		      this->mismatch_positions_alloc,
		      query_compress_fwd,query_compress_rev,querylength,
		      first_read_p,listpool,hitlistpool,/*level*/6);

    if (*found_score_overall <= nmismatches_allowed) {
      Substring_list_gc(&startfrags_plus);
      Substring_list_gc(&endfrags_plus);
      Substring_list_gc(&startfrags_minus);
      Substring_list_gc(&endfrags_minus);
      return DISTANT_DNA;

    } else {
      debug(printf("Consider terminals.  Have %d hits plus and %d hits minus\n",List_length(*hits_gplus),List_length(*hits_gminus)));
      /* Terminals.  Use hits == NULL as criterion, rather than
	 found_score, since Terminal_solve cannot improve found_score */
      if (*hits_gplus == NULL && *hits_gminus == NULL) {
	/* No need to recompute min_trim, if hits_gplus and hits_gminus are NULL */

#if 0
	if (elt_set_extended_p == false) {
	  Elt_set_queryfwd_extend_left(this->queryfwd_plus_set,/*query_compress*/query_compress_fwd,
				       /*plusp*/true,genestrand);
	  Elt_set_queryfwd_extend_left(this->queryfwd_minus_set,/*query_compress*/query_compress_rev,
				       /*plusp*/false,genestrand);
	  Elt_set_queryrev_extend_right(this->queryrev_plus_set,querylength,/*query_compress*/query_compress_fwd,
					/*plusp*/true,genestrand);
	  Elt_set_queryrev_extend_right(this->queryrev_minus_set,querylength,/*query_compress*/query_compress_rev,
					/*plusp*/false,genestrand);
	  elt_set_extended_p = true;
	}
#endif

	/* max_terminal_mismatches = done_level; */
	debug(printf("LEVEL %d: RUNNING TERMINAL_SOLVE_PLUS\n",7));
	if ((*hits_gplus = Terminal_solve_plus(&(*found_score_overall),&(*found_score_within_trims),
					       this->queryfwd_plus_set,this->queryrev_plus_set,
					       query_compress_fwd,querylength,
					       genestrand,listpool,hitlistpool,/*level*/7)) == NULL) {
	  debug(printf("Got no terminals\n"));
	} else {
	  debug(printf("Got %d terminals\n",List_length(*hits_gplus)));
	}
	
	debug(printf("LEVEL %d: RUNNING TERMINAL_SOLVE_MINUS\n",7));
	if ((*hits_gminus = Terminal_solve_minus(&(*found_score_overall),&(*found_score_within_trims),
						 this->queryfwd_minus_set,this->queryrev_minus_set,
						 query_compress_rev,querylength,
						 genestrand,listpool,hitlistpool,/*level*/7)) == NULL) {
	  debug(printf("Got no terminals\n"));
	} else {
	  debug(printf("Got %d terminals\n",List_length(*hits_gminus)));
	}
      }
    }
  }

  Substring_list_gc(&startfrags_plus);
  Substring_list_gc(&endfrags_plus);
  Substring_list_gc(&startfrags_minus);
  Substring_list_gc(&endfrags_minus);
  return TERMINAL;
}


Stage3end_T *
Stage1_single_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    Shortread_T queryseq,
#if 0
		    Oligoindex_array_T oligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
#endif
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool) {
  Stage3end_T *stage3array;
  T this;
  List_T hits, geneplus_hits, geneminus_hits, hits_gplus, hits_gminus;

  int found_score_overall, found_score_within_trims;
  int querylength, max_mismatches, min_coverage;
  char *queryuc_ptr, *queryrc, *quality_string;
  Compress_T query_compress_fwd, query_compress_rev;


  querylength = Shortread_fulllength(queryseq);
  queryuc_ptr = Shortread_fullpointer_uc(queryseq);

  queryrc = (char *) MALLOC((querylength+1)*sizeof(char));
  make_complement_buffered(queryrc,queryuc_ptr,querylength);

  query_compress_fwd = Compress_new_fwd(queryuc_ptr,querylength);
  query_compress_rev = Compress_new_rev(queryuc_ptr,querylength);

  if (mode == STANDARD || mode == CMET_STRANDED || mode == ATOI_STRANDED || mode == TTOC_STRANDED) {
    this = Stage1_new(querylength);
    single_read(&found_score_overall,&found_score_within_trims,&hits_gplus,&hits_gminus,this,
		queryuc_ptr,queryrc,querylength,query_compress_fwd,query_compress_rev,/*genestrand*/0,
		/*oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,*/
		intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
		/*paired_end_p*/false,/*first_read_p*/true);
    Stage1_free(&this);
    hits = List_append(hits_gplus,hits_gminus);


  } else if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
    this = Stage1_new(querylength);
    single_read(&found_score_overall,&found_score_within_trims,&hits_gplus,&hits_gminus,this,
		queryuc_ptr,queryrc,querylength,query_compress_fwd,query_compress_rev,/*genestrand*/+1,
		/*oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,*/
		intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
		/*paired_end_p*/false,/*first_read_p*/true);
    Stage1_free(&this);
    geneplus_hits = List_append(hits_gplus,hits_gminus);

    this = Stage1_new(querylength);
    single_read(&found_score_overall,&found_score_within_trims,&hits_gplus,&hits_gminus,this,
		queryuc_ptr,queryrc,querylength,query_compress_fwd,query_compress_rev,/*genestrand*/+2,
		/*oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,*/
		intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
		/*paired_end_p*/false,/*first_read_p*/true);
    Stage1_free(&this);
    geneminus_hits = List_append(hits_gplus,hits_gminus);

    hits = List_append(geneplus_hits,geneminus_hits);

  } else {
    fprintf(stderr,"Do not recognize mode %d\n",mode);
    abort();
  }

  /* Previously used user_maxlevel_float in searching, now just for filtering */
  if (user_maxlevel_float < 0.0) {
    max_mismatches = querylength;
  } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
    max_mismatches = (int) rint(user_maxlevel_float * (double) querylength);
  } else {
    max_mismatches = (int) user_maxlevel_float;
  }

  if (user_mincoverage_float < 0.0) {
    min_coverage = 0;
  } else if (user_mincoverage_float > 0.0 && user_mincoverage_float < 1.0) {
    min_coverage = (int) rint(user_mincoverage_float * (double) querylength);
  } else {
    min_coverage = (int) user_mincoverage_float;
  }

  debug(printf("Have %d hits before filter_coverage\n",List_length(hits)));
  hits = Stage3end_filter(hits,hitlistpool,max_mismatches,min_coverage);
  debug(printf("Have %d hits after filter_coverage\n",List_length(hits)));

  if (hits == NULL) {
    *npaths_primary = *npaths_altloc = 0;
    stage3array = (Stage3end_T *) NULL;
    
  } else {
    hits = Stage3end_remove_circular_alias(hits,hitlistpool); /* Contains a call to unalias_circular */
    hits = Stage3end_remove_duplicates(hits,hitlistpool); /* Aliases can cause duplicates */
    debug(printf("Have %d hits after remove_duplicates\n",List_length(hits)));
    
    hits = Stage3end_optimal_score(hits,hitlistpool,querylength,/*finalp*/false);
    debug(printf("Have %d hits after optimal_score\n",List_length(hits)));
    hits = Stage3end_remove_overlaps(hits,hitlistpool,/*finalp*/true);
    debug(printf("Have %d hits after remove_overlaps\n",List_length(hits)));
    hits = Stage3end_optimal_score(hits,hitlistpool,querylength,/*finalp*/true);
    debug(printf("Have %d hits after final optimal_score\n",List_length(hits)));

    quality_string = Shortread_quality_string(queryseq);
    Stage3end_count_hits(&(*npaths_primary),&(*npaths_altloc),hits);
    stage3array = (Stage3end_T *) List_to_array_out(hits,NULL); Hitlist_free(&hits); /* Return value */
    stage3array = Stage3end_eval_and_sort(/*npaths*/(*npaths_primary) + (*npaths_altloc),
					  &(*first_absmq),&(*second_absmq),
					  stage3array,queryuc_ptr,quality_string,/*displayp*/true);
  }
  
  FREE(queryrc);

  Compress_free(&query_compress_fwd);
  Compress_free(&query_compress_rev);

  return stage3array;
}




static Pairtype_T
choose_among_paired (int *best_nmatches_paired, int *best_nmatches_5, int *best_nmatches_3,
		     List_T hitpairs, List_T samechr, List_T conc_transloc) {
  Pairtype_T final_pairtype = UNPAIRED;
  List_T p;
  Stage3pair_T hitpair;
  int nmatches, nmatches5, nmatches3;

  debug16(printf("choose: %d hitpairs, %d conc_transloc, %d samechr\n",
		 List_length(hitpairs),List_length(conc_transloc),List_length(samechr)));

  *best_nmatches_paired = 0;
  *best_nmatches_5 = 0;
  *best_nmatches_3 = 0;
  for (p = hitpairs; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if ((nmatches = Stage3pair_nmatches_to_trims(&nmatches5,&nmatches3,hitpair)) > *best_nmatches_paired) {
      final_pairtype = CONCORDANT;
      *best_nmatches_paired = nmatches;
      *best_nmatches_5 = nmatches5;
      *best_nmatches_3 = nmatches3;
    }
  }

  *best_nmatches_paired += 1; /* penalty for choosing translocation over others */

  for (p = conc_transloc; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if ((nmatches = Stage3pair_nmatches_to_trims(&nmatches5,&nmatches3,hitpair)) > *best_nmatches_paired) {
      final_pairtype = CONCORDANT_TRANSLOCATIONS;
      *best_nmatches_paired = nmatches;
      *best_nmatches_5 = nmatches5;
      *best_nmatches_3 = nmatches3;
    }
  }

  for (p = samechr; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if ((nmatches = Stage3pair_nmatches_to_trims(&nmatches5,&nmatches3,hitpair)) > *best_nmatches_paired) {
      final_pairtype = PAIRED_UNSPECIFIED;
      *best_nmatches_paired = nmatches;
      *best_nmatches_5 = nmatches5;
      *best_nmatches_3 = nmatches3;
    }
  }

  debug16(printf("best_nmatches_paired among paired = %d = %d + %d\n",
		 *best_nmatches_paired,*best_nmatches_5,*best_nmatches_3));
  debug16(printf("final pairtype: %s\n",Pairtype_string(final_pairtype)));

  return final_pairtype;
}



static List_T
paired_read_segment_search (bool *abort_pairing_p, List_T *samechr, List_T *conc_transloc,
			    int last_level_5, int last_level_3, 
			    int found_score_overall_5, int found_score_within_trims_5,
			    int found_score_overall_3, int found_score_within_trims_3, T this5, T this3,
			    Ladder_T ladder5_plus, Ladder_T ladder5_minus, Ladder_T ladder3_plus, Ladder_T ladder3_minus,
			    char *queryuc_ptr_5, char *queryrc5, int querylength5, char *queryuc_ptr_3, char *queryrc3, int querylength3,
			    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			    int genestrand, Chrpos_T pairmax_linear,
#if 0
			    Oligoindex_array_T oligoindices_minor,
			    Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
#endif
			    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool) {
  List_T hitpairs;
  int last_level;
  
  /* int sensedir; -- Used by Pair_solve_via_gmap */

  /* int plus_specifici_5[3], minus_specifici_5[3], plus_specifici_3[3], minus_specifici_3[3]; */
  Univcoord_T *gplus5_diagonals, *gminus5_diagonals, *gplus3_diagonals, *gminus3_diagonals;
  int gplus5_ndiagonals, gminus5_ndiagonals, gplus3_ndiagonals, gminus3_ndiagonals;

  /* concordant_score is for all hitpairs; adjacent score is for those hitpairs within adjacent_pairlength */
  int adjacent_score, concordant_score;
  int done_level_5, done_level_3, nmismatches_allowed, nmismatches_allowed_5, nmismatches_allowed_3;
  /* int max_terminal_mismatches; */
  int max_splice_mismatches;
  int min_trim_5, min_trim_3;

  int maxpairedpaths;

  List_T hits5_gplus = NULL, hits5_gminus = NULL, hits3_gplus = NULL, hits3_gminus = NULL;

#if 0
  /* If we prioritize transcriptome hits over genomic hits */
  List_T transcriptome_hits5 = NULL, transcriptome_hits3 = NULL;
  List_T genome_hits5 = NULL, genome_hits5 = NULL;
#endif


  /* For Segment_search */
  struct Record_T *plus_records5 = NULL, *minus_records5 = NULL, *plus_records3 = NULL, *minus_records3 = NULL;
  int plus_nrecords5 = 0, minus_nrecords5 = 0, plus_nrecords3 = 0, minus_nrecords3 = 0;

  List_T startfrags5_plus = NULL, endfrags5_plus = NULL, startfrags5_minus = NULL, endfrags5_minus = NULL,
    startfrags3_plus = NULL, endfrags3_plus = NULL, startfrags3_minus = NULL, endfrags3_minus = NULL;

  int level;


  debug(printf("Entered paired_read_segment_search with queryuc_ptr_5 %s\n",queryuc_ptr_5));
  debug(printf("Entered paired_read_segment_search with queryuc_ptr_3 %s\n",queryuc_ptr_3));

  *abort_pairing_p = false;
  hitpairs = (List_T) NULL;
  *samechr = (List_T) NULL;
  /* *conc_transloc = (List_T) NULL; -- May contain results from concordant alignment of single ends */

  if (last_level_5 < last_level_3) {
    level = last_level = last_level_5;
  } else {
    level = last_level = last_level_3;
  }
  debug(printf("last_level is %s <- %s and %s\n",
	       Method_string(last_level),Method_string(last_level_5),Method_string(last_level_3)));


  /* Take the larger of maxpaths_search and 10*maxpaths_report */
  maxpairedpaths = maxpaths_search;
  if (maxpairedpaths < 10*maxpaths_report) {
    maxpairedpaths = 10*maxpaths_report;
  }

  adjacent_score = concordant_score = querylength5 + querylength3;

#if 0
  /* Use user_maxlevel_float for final filtering, not for search */
  if (user_maxlevel_float < 0.0) {
    nmismatches_allowed_5 = querylength5/index1part_tr;
    nmismatches_allowed_3 = querylength3/index1part_tr;
  } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
    nmismatches_allowed_5 = (int) rint(user_maxlevel_float * (double) querylength5);
    nmismatches_allowed_3 = (int) rint(user_maxlevel_float * (double) querylength3);
  } else {
    nmismatches_allowed_5 = nmismatches_allowed_3 = (int) user_maxlevel_float;
  }
#else
  nmismatches_allowed_5 = querylength5/index1part_tr;
  nmismatches_allowed_3 = querylength3/index1part_tr;
#endif

  nmismatches_allowed = nmismatches_allowed_5 + nmismatches_allowed_3;


  /* found_score_5 = querylength5; */
  /* found_score_3 = querylength3; */
  /* found_score = querylength5 + querylength3; */

	
  if (last_level < KMER_EXACT) {
    if (last_level_5 < KMER_EXACT) {
      debug(printf("LEVEL %d: RUNNING KMER_SEARCH_GENOME_ENDS ON 5' READ\n",level));
      Stage1_init(this5,queryuc_ptr_5,querylength5,genestrand);
    }
    
    if (last_level_3 < KMER_EXACT) {
      debug(printf("LEVEL %d: RUNNING KMER_SEARCH_GENOME_ENDS ON 3' READ\n",level));
      Stage1_init(this3,queryuc_ptr_3,querylength3,genestrand);
    }
  }

  if (last_level < EXT) {
    if (last_level_5 < EXT && querylength5 >= index1part) {
      Stage1_rearrange(this5,querylength5,/*one5p*/false);
      Stage1_fill_all_oligos(this5,queryuc_ptr_5,querylength5,genestrand);
    }

    if (last_level_3 < EXT && querylength3 >= index1part) {
      Stage1_rearrange(this3,querylength3,/*one3p*/false);
      Stage1_fill_all_oligos(this3,queryuc_ptr_3,querylength3,genestrand);
    }
  }

  /* SEGMENT2 algorithm is different from SEGMENT1, so we should do this in all cases */
  if (querylength5 >= index1part) {
    debug(printf("LEVEL %d: RUNNING SEGMENT_SEARCH_GENOME (ANCHORED) ON 5' READ\n",level));
    if (last_level_5 < SEGMENT1) {
      Stage1_fill_position_ptrs_plus(this5,/*querystart*/0,/*queryend*/querylength5,genestrand);
      Stage1_fill_position_ptrs_minus(this5,/*querystart*/0,/*queryend*/querylength5,genestrand);
    }
    
    gplus3_diagonals = Ladder_genomicstarts(&gplus3_ndiagonals,ladder3_plus);
    gminus3_diagonals = Ladder_genomicends(&gminus3_ndiagonals,ladder3_minus);

    plus_records5 = Segment_identify_lower(&plus_nrecords5,
#ifdef LARGE_GENOMES
					   this5->plus_positions_high,
#endif
					   this5->plus_positions,this5->plus_npositions,this5->plus_validp,
					   this5->stream_alloc,this5->streamsize_alloc,this5->querypos_diagterm_alloc,
					   /*max_pairlength*/pairmax_linear,querylength5,
					   gplus3_diagonals,gplus3_ndiagonals);
    
    minus_records5 = Segment_identify_higher(&minus_nrecords5,
#ifdef LARGE_GENOMES
					     this5->minus_positions_high,
#endif
					     this5->minus_positions,this5->minus_npositions,this5->minus_validp,
					     this5->stream_alloc,this5->streamsize_alloc,this5->querypos_diagterm_alloc,
					     /*max_pairlength*/pairmax_linear,querylength5,
					     gminus3_diagonals,gminus3_ndiagonals);
    FREE(gminus3_diagonals);
    FREE(gplus3_diagonals);
    
    Segment_search_all(&found_score_overall_5,&found_score_within_trims_5,&hits5_gplus,&hits5_gminus,
		       plus_records5,plus_nrecords5,minus_records5,minus_nrecords5,
		       
		       queryuc_ptr_5,queryrc5,querylength5,this5->mismatch_positions_alloc,
		       this5->spliceinfo,this5->stream_alloc,this5->streamsize_alloc,
		       query5_compress_fwd,query5_compress_rev,genestrand,/*paired_end_p*/true,/*first_read_p*/true,
		       intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
		       /*method*/SEGMENT2,level);
  }
  
  if (querylength3 >= index1part) {
    debug(printf("LEVEL %d: RUNNING SEGMENT_SEARCH_GENOME (ANCHORED) ON 3' READ\n",level));
    if (last_level_3 < SEGMENT1) {
      Stage1_fill_position_ptrs_plus(this3,/*querystart*/0,/*queryend*/querylength3,genestrand);
      Stage1_fill_position_ptrs_minus(this3,/*querystart*/0,/*queryend*/querylength3,genestrand);
    }
    
    gplus5_diagonals = Ladder_genomicstarts(&gplus5_ndiagonals,ladder5_plus);
    gminus5_diagonals = Ladder_genomicends(&gminus5_ndiagonals,ladder5_minus);
    
    plus_records3 = Segment_identify_higher(&plus_nrecords3,
#ifdef LARGE_GENOMES
					    this3->plus_positions_high,
#endif
					    this3->plus_positions,this3->plus_npositions,this3->plus_validp,
					    this3->stream_alloc,this3->streamsize_alloc,this3->querypos_diagterm_alloc,
					    /*max_pairlength*/pairmax_linear,querylength3,
					    gplus5_diagonals,gplus5_ndiagonals);
    
    minus_records3 = Segment_identify_lower(&minus_nrecords3,
#ifdef LARGE_GENOMES
					    this3->minus_positions_high,
#endif
					    this3->minus_positions,this3->minus_npositions,this3->minus_validp,
					    this3->stream_alloc,this3->streamsize_alloc,this3->querypos_diagterm_alloc,
					    /*max_pairlength*/pairmax_linear,querylength3,
					    gminus5_diagonals,gminus5_ndiagonals);
    FREE(gminus5_diagonals);
    FREE(gplus5_diagonals);
    
    Segment_search_all(&found_score_overall_3,&found_score_within_trims_3,&hits3_gplus,&hits3_gminus,
		       plus_records3,plus_nrecords3,minus_records3,minus_nrecords3,
		       
		       queryuc_ptr_3,queryrc3,querylength3,this3->mismatch_positions_alloc,
		       this3->spliceinfo,this3->stream_alloc,this3->streamsize_alloc,
		       query3_compress_fwd,query3_compress_rev,genestrand,/*paired_end_p*/true,/*first_read_p*/false,
		       intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
		       /*method*/SEGMENT2,level);
  }
  
  debug(printf("found scores: %d and %d (vs %d and %d allowed)\n",
	       found_score_within_trims_5,found_score_within_trims_3,nmismatches_allowed_5,nmismatches_allowed_3));
  debug(printf("hits5 plus: %d new.  hits5 minus: %d new.\n",List_length(hits5_gplus),List_length(hits5_gminus)));
  debug(printf("hits3 plus: %d new.  hits3 minus: %d new.\n",List_length(hits3_gplus),List_length(hits3_gminus)));
  hitpairs = Concordance_pair_up_genome(&(*abort_pairing_p),&adjacent_score,&concordant_score,&(*conc_transloc),
					hitpairs,hits5_gplus,hits5_gminus,hits3_gplus,hits3_gminus,
					ladder5_plus,ladder5_minus,ladder3_plus,ladder3_minus,
					querylength5,querylength3,
#if 0
					queryuc_ptr_5,queryuc_ptr_3,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					listpool,hitlistpool,maxpairedpaths,genestrand); level++;
  Hitlist_free(&hits5_gplus); Hitlist_free(&hits5_gminus); Hitlist_free(&hits3_gplus); Hitlist_free(&hits3_gminus);
  debug(printf("(4) After level %d, we have adjacent_score %d, concordant_score %d (vs allowed %d), %d pairs\n",
	       level,adjacent_score,concordant_score,nmismatches_allowed,List_length(hitpairs)));
    
  FREE(plus_records3); FREE(minus_records3); FREE(plus_records5); FREE(minus_records5);
  

  if (*abort_pairing_p == true || concordant_score <= nmismatches_allowed) {
    return hitpairs;

  } else if (last_level < DISTANT_RNA) {
    /* Find distant splicing.  Use concordant_score, rather than
       adjacent_score as criterion, since distant splicing cannot find
       adjacent pairs */
      
    /* distant splicing on 5' read */
    if (last_level_5 < DISTANT_RNA && found_score_within_trims_5 >= found_score_within_trims_3) {
      if ((done_level_5 = found_score_within_trims_5 + subopt_levels) > nmismatches_allowed_5) {
	done_level_5 = nmismatches_allowed_5;
      }
      debug(printf("done_level_5 %d = found_score_within_trims_5 %d + subopt_levels %d\n",
		   done_level_5,found_score_within_trims_5,subopt_levels));
      
      /* min of max(trim5,trim3) over each hit */      
      min_trim_5 = Ladder_minimax_trim(ladder5_plus,ladder5_minus,querylength5);
      
      debug(printf("For 5' end, comparing min_trim %d with min_distantsplicing_end_matches %d\n",
		   min_trim_5,min_distantsplicing_end_matches));
      if (min_trim_5 < min_distantsplicing_end_matches) {
	/* Don't find distant splicing */
	debug(printf("For 5' end, skipping distant splicing because min_trim %d < min_distantsplicing_end_matches %d\n",
		     min_trim_5,min_distantsplicing_end_matches));
	
      } else if ((max_splice_mismatches = done_level_5 - distantsplicing_penalty) >= 0) {
	if (novelsplicingp == true) {
	  debug(printf("For 5' end, candidate for distant splicing because min_trim %d >= min_distantsplicing_end_matches %d\n",
		       min_trim_5,min_distantsplicing_end_matches));
	  debug(printf("Have sets: %d %d %d %d\n",
		       List_length(this5->queryfwd_plus_set),List_length(this5->queryfwd_minus_set),
		       List_length(this5->queryrev_plus_set),List_length(this5->queryrev_minus_set)));
	  Distant_rna_solve(&found_score_overall_5,&found_score_within_trims_5,&hits5_gplus,&hits5_gminus,
			    &startfrags5_plus,&endfrags5_plus,&startfrags5_minus,&endfrags5_minus,

			    this5->queryfwd_plus_set,this5->queryfwd_minus_set,
			    this5->queryrev_plus_set,this5->queryrev_minus_set,
			    
			    this5->mismatch_positions_alloc,this5->positions_alloc,
			    query5_compress_fwd,query5_compress_rev,queryuc_ptr_5,queryrc5,querylength5,
			    max_splice_mismatches,genestrand,/*first_read_p*/true,listpool,hitlistpool,level);
	} else {
	  /* TODO: Implement distant DNA fusions */
	}
      }
    }
    
    /* Distant splicing on 3' read */
    if (last_level_3 < DISTANT_RNA && found_score_within_trims_3 >= found_score_within_trims_5) {
      if ((done_level_3 = found_score_within_trims_3 + subopt_levels) > nmismatches_allowed_3) {
	done_level_3 = nmismatches_allowed_3;
      }
      debug(printf("done_level_3 %d = found_score_within_trims_3 %d + subopt_levels %d\n",
		   done_level_3,found_score_within_trims_3,subopt_levels));
      
      /* min of max(trim5,trim3) over each hit */      
      min_trim_3 = Ladder_minimax_trim(ladder3_plus,ladder3_minus,querylength3);
      
      debug(printf("For 3' end, comparing min_trim %d with min_distantsplicing_end_matches %d\n",
		   min_trim_3,min_distantsplicing_end_matches));
      if (min_trim_3 < min_distantsplicing_end_matches) {
	/* Don't find distant splicing */
	debug(printf("For 3'end, skipping distant splicing because min_trim %d < min_distantsplicing_end_matches %d\n",
		     min_trim_3,min_distantsplicing_end_matches));
	
      } else if ((max_splice_mismatches = done_level_3 - distantsplicing_penalty) >= 0) {
	if (novelsplicingp == true) {
	  debug(printf("For 3' end, candidate for distant splicing because min_trim %d >= min_distantsplicing_end_matches %d\n",
		       min_trim_3,min_distantsplicing_end_matches));
	  debug(printf("Have sets: %d %d %d %d\n",
		       List_length(this3->queryfwd_plus_set),List_length(this3->queryfwd_minus_set),
		       List_length(this3->queryrev_plus_set),List_length(this3->queryrev_minus_set)));
	  Distant_rna_solve(&found_score_overall_3,&found_score_within_trims_3,&hits3_gplus,&hits3_gminus,
			    &startfrags3_plus,&endfrags3_plus,&startfrags3_minus,&endfrags3_minus,

			    this3->queryfwd_plus_set,this3->queryfwd_minus_set,
			    this3->queryrev_plus_set,this3->queryrev_minus_set,
			    
			    this3->mismatch_positions_alloc,this3->positions_alloc,
			    query3_compress_fwd,query3_compress_rev,queryuc_ptr_3,queryrc3,querylength3,
			    max_splice_mismatches,genestrand,/*first_read_p*/false,listpool,hitlistpool,level);
	} else {
	  /* TODO: Implement distant DNA */
	}
      }
    }
    
    debug(printf("found scores: %d and %d (vs %d and %d allowed)\n",
		 found_score_within_trims_5,found_score_within_trims_3,nmismatches_allowed_5,nmismatches_allowed_3));
    debug(printf("hits5 plus: %d new.  hits5 minus: %d new.\n",List_length(hits5_gplus),List_length(hits5_gminus)));
    debug(printf("hits3 plus: %d new.  hits3 minus: %d new.\n",List_length(hits3_gplus),List_length(hits3_gminus)));
    hitpairs = Concordance_pair_up_distant(&(*abort_pairing_p),&concordant_score,&(*samechr),&(*conc_transloc),hitpairs,
					   hits5_gplus,hits5_gminus,hits3_gplus,hits3_gminus,
					   ladder5_plus,ladder5_minus,ladder3_plus,ladder3_minus,
					   querylength5,querylength3,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,
					   query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,hitlistpool,maxpairedpaths,genestrand); level++;
    Hitlist_free(&hits5_gplus); Hitlist_free(&hits5_gminus); Hitlist_free(&hits3_gplus); Hitlist_free(&hits3_gminus);
    debug(printf("(6) After level %d, we have adjacent_score %d, concordant_score %d (vs allowed %d), %d pairs\n",
		 level,adjacent_score,concordant_score,nmismatches_allowed,List_length(hitpairs)));
  }
  
  
  if (last_level < DISTANT_RNA && (*abort_pairing_p == true || concordant_score <= nmismatches_allowed)) {
    Substring_list_gc(&startfrags3_plus);
    Substring_list_gc(&endfrags3_plus);
    Substring_list_gc(&startfrags3_minus);
    Substring_list_gc(&endfrags3_minus);
    Substring_list_gc(&startfrags5_plus);
    Substring_list_gc(&endfrags5_plus);
    Substring_list_gc(&startfrags5_minus);
    Substring_list_gc(&endfrags5_minus);
    return hitpairs;

  } else if (last_level < TERMINAL) {
    /* Terminals */
    if (last_level_5 < TERMINAL) {
      /* Need to recompute, because of new distant hits or because distant splicing was not run */
      min_trim_5 = Ladder_minimax_trim(ladder5_plus,ladder5_minus,querylength5);

      if (min_trim_5 < min_distantsplicing_end_matches) {
	debug(printf("For 5' end, skipping terminals because min_trim %d < min_distantsplicing_end_matches %d\n",
		     min_trim_5,min_distantsplicing_end_matches));
      } else {
	/* max_terminal_mismatches = done_level_5; */
	if ((hits5_gplus = Terminal_solve_plus(&found_score_overall_5,&found_score_within_trims_5,
					       this5->queryfwd_plus_set,this5->queryrev_plus_set,
					       query5_compress_fwd,querylength5,
					       genestrand,listpool,hitlistpool,level)) == NULL) {
	  debug(printf("Got no terminals\n"));
	} else {
	  debug(printf("Got %d terminals\n",List_length(hits5_gplus)));
	}
	
	if ((hits5_gminus = Terminal_solve_minus(&found_score_overall_5,&found_score_within_trims_5,
						 this5->queryfwd_minus_set,this5->queryrev_minus_set,
						 query5_compress_rev,querylength5,
						 genestrand,listpool,hitlistpool,level)) == NULL) {
	  debug(printf("Got no terminals\n"));
	} else {
	  debug(printf("Got %d terminals\n",List_length(hits5_gminus)));
	}
      }
    }
    
    if (last_level_3 < DISTANT_RNA) {
      /* Need to recompute, because of new distant hits or because distant splicing was not run */
      min_trim_3 = Ladder_minimax_trim(ladder3_plus,ladder3_minus,querylength3);
      
      if (min_trim_3 < min_distantsplicing_end_matches) {
	debug(printf("For 3'end, skipping terminals because min_trim %d < min_distantsplicing_end_matches %d\n",
		     min_trim_3,min_distantsplicing_end_matches));
      } else {
	/* max_terminal_mismatches = done_level_3; */
	if ((hits3_gplus = Terminal_solve_plus(&found_score_overall_3,&found_score_within_trims_3,
					       this3->queryfwd_plus_set,this3->queryrev_plus_set,
					       query3_compress_fwd,querylength3,
					       genestrand,listpool,hitlistpool,level)) == NULL) {
	  debug(printf("Got no terminals\n"));
	} else {
	  debug(printf("Got %d terminals\n",List_length(hits3_gplus)));
	}
	
	if ((hits3_gminus = Terminal_solve_minus(&found_score_overall_3,&found_score_within_trims_3,
						 this3->queryfwd_minus_set,this3->queryrev_minus_set,
						 query3_compress_rev,querylength3,
						 genestrand,listpool,hitlistpool,level)) == NULL) {
	  debug(printf("Got no terminals\n"));
	} else {
	  debug(printf("Got %d terminals\n",List_length(hits3_gminus)));
	}
      }
    }
    
    hitpairs = Concordance_pair_up_genome(&(*abort_pairing_p),&adjacent_score,&concordant_score,&(*conc_transloc),
					  hitpairs,hits5_gplus,hits5_gminus,hits3_gplus,hits3_gminus,
					  ladder5_plus,ladder5_minus,ladder3_plus,ladder3_minus,
					  querylength5,querylength3,
#if 0
					  queryuc_ptr_5,queryuc_ptr_3,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					  listpool,hitlistpool,maxpairedpaths,genestrand); level++;
    Hitlist_free(&hits5_gplus); Hitlist_free(&hits5_gminus); Hitlist_free(&hits3_gplus); Hitlist_free(&hits3_gminus);
    debug(printf("(7) After level %d, we have adjacent_score %d, concordant_score %d (vs allowed %d), %d pairs\n",
		 level,adjacent_score,concordant_score,nmismatches_allowed,List_length(hitpairs)));
  }

  Substring_list_gc(&startfrags3_plus);
  Substring_list_gc(&endfrags3_plus);
  Substring_list_gc(&startfrags3_minus);
  Substring_list_gc(&endfrags3_minus);
  Substring_list_gc(&startfrags5_plus);
  Substring_list_gc(&endfrags5_plus);
  Substring_list_gc(&startfrags5_minus);
  Substring_list_gc(&endfrags5_minus);

  debug(printf("Exiting paired_read_segment_search\n"));
  return hitpairs;
}



static List_T
paired_read (bool *abort_pairing_p, List_T *hits5, List_T *hits3, List_T *samechr, List_T *conc_transloc,
	     char *queryuc_ptr_5, char *queryrc5, int querylength5, char *queryuc_ptr_3, char *queryrc3, int querylength3,
	     Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
	     Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
	     int genestrand, Chrpos_T pairmax_linear,
#if 0
	     Oligoindex_array_T oligoindices_minor,
	     Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
	     Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
#endif
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool) {

  List_T hits5_gplus, hits5_gminus, hits3_gplus, hits3_gminus;
  int adjacent_score, concordant_score;
  List_T hitpairs;
  int maxpairedpaths;
  int level5, level3;
  int found_score_overall_5, found_score_within_trims_5, found_score_overall_3, found_score_within_trims_3;

  T this5, this3;
  Ladder_T ladder5_plus, ladder5_minus, ladder3_plus, ladder3_minus;


  *abort_pairing_p = false;
  hitpairs = (List_T) NULL;
  *samechr = (List_T) NULL;
  *conc_transloc = (List_T) NULL;

  /* Take the larger of maxpaths_search and 10*maxpaths_report */
  maxpairedpaths = maxpaths_search;
  if (maxpairedpaths < 10*maxpaths_report) {
    maxpairedpaths = 10*maxpaths_report;
  }

  adjacent_score = concordant_score = querylength5 + querylength3;


  /* TODO: Return level from single_read.  Then if we don't have concordant hitpairs, take the min(level5,level3) and start paired_readfrom t
here */

  this5 = Stage1_new(querylength5);
  level5 = single_read(&found_score_overall_5,&found_score_within_trims_5,&hits5_gplus,&hits5_gminus,this5,
		       queryuc_ptr_5,queryrc5,querylength5,query5_compress_fwd,query5_compress_rev,genestrand,
		       /*oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,*/
		       intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
		       /*paired_end_p*/true,/*first_read_p*/true);
  debug(printf("Got %d plus and %d minus hits for 5' end\n",List_length(hits5_gplus),List_length(hits5_gminus)));

  this3 = Stage1_new(querylength3);
  level3 = single_read(&found_score_overall_3,&found_score_within_trims_3,&hits3_gplus,&hits3_gminus,this3,
		       queryuc_ptr_3,queryrc3,querylength3,query3_compress_fwd,query3_compress_rev,genestrand,
		       /*oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,*/
		       intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool,
		       /*paired_end_p*/true,/*first_read_p*/false);
  debug(printf("Got %d plus and %d minus hits for 3' end\n",List_length(hits3_gplus),List_length(hits3_gminus)));

  ladder5_plus = Ladder_new(NULL,hitlistpool,/*end5p*/true);
  ladder5_minus = Ladder_new(NULL,hitlistpool,/*end5p*/true);
  ladder3_plus = Ladder_new(NULL,hitlistpool,/*end5p*/false);
  ladder3_minus = Ladder_new(NULL,hitlistpool,/*end5p*/false);

  hitpairs = Concordance_pair_up_genome(&(*abort_pairing_p),&adjacent_score,&concordant_score,&(*conc_transloc),
					hitpairs,hits5_gplus,hits5_gminus,hits3_gplus,hits3_gminus,
					ladder5_plus,ladder5_minus,ladder3_plus,ladder3_minus,
					querylength5,querylength3,
#if 0
					queryuc_ptr_5,queryuc_ptr_3,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					listpool,hitlistpool,maxpairedpaths,genestrand);
  Hitlist_free(&hits5_gplus); Hitlist_free(&hits5_gminus); Hitlist_free(&hits3_gplus); Hitlist_free(&hits3_gminus);
  debug(printf("After initial concordance, have %d concordant hitpairs and %d concordant transloc\n\n",
	       List_length(hitpairs),List_length(*conc_transloc)));

  /* DEBUGGING */
  if (hitpairs == NULL) {
    hitpairs = paired_read_segment_search(&(*abort_pairing_p),&(*samechr),&(*conc_transloc),
					  level5,level3,found_score_overall_5,found_score_within_trims_5,
					  found_score_overall_3,found_score_within_trims_3,this5,this3,
					  ladder5_plus,ladder5_minus,ladder3_plus,ladder3_minus,
					  queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  genestrand,pairmax_linear,
					  intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool);
  }

  Stage1_free(&this3); Stage1_free(&this5);
  Ladder_to_hits(&(*hits5),&(*hits3),&ladder5_plus,&ladder5_minus,&ladder3_plus,&ladder3_minus,
		 hitlistpool);

  debug(printf("Exiting paired_read\n"));
  return hitpairs;
}



/* Have three lists: hitpairs, samechr, and conc_transloc => result */
static Stage3pair_T *
consolidate_paired_results (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
			    Stage3end_T **stage3array5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
			    Stage3end_T **stage3array3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
			    List_T hitpairs, List_T samechr, List_T conc_transloc, List_T hits5, List_T hits3,

			    char *queryuc_ptr_5, int querylength5, char *queryuc_ptr_3, int querylength3,
			    char *quality_string_5, char *quality_string_3,
#if 0
			    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
#endif
			    int max_mismatches_5, int max_mismatches_3, int min_coverage_5, int min_coverage_3,
			    Listpool_T listpool, Hitlistpool_T hitlistpool) {
  Stage3pair_T *stage3pairarray, stage3pair, newpair;
  Stage3end_T hit5, hit3;
  List_T result, singlehits5, singlehits3, p;
  int genestrand;
  Pairtype_T pairtype;
  int best_nmatches_paired, best_nmatches_paired_5, best_nmatches_paired_3;

  
  *final_pairtype = choose_among_paired(&best_nmatches_paired,&best_nmatches_paired_5,&best_nmatches_paired_3,
					hitpairs,samechr,conc_transloc);
  debug16(printf("Entered consolidate_paired_results with final_pairtype %d\n",*final_pairtype));

  /* DEBUGGING */
  if (*final_pairtype == CONCORDANT) {
    /* Have concordant results */
    debug16(printf("Have %d concordant results, %d samechr, and %d conc_transloc\n",
		   List_length(hitpairs),List_length(samechr),List_length(conc_transloc)));
    for (p = samechr; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    Hitlist_free(&samechr);

    for (p = conc_transloc; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    Hitlist_free(&conc_transloc);
	  
    if (1) {
      debug16(printf("Via transcriptome: Before removing overlaps, %d results\n",List_length(hitpairs)));
      result = Stage3pair_optimal_score(hitpairs,hitlistpool,querylength5,querylength3,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,hitlistpool,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/true);
      
#ifdef RESOLVE_INSIDE_GENERAL
      result = Stage3pair_resolve_insides(result,queryuc_ptr_5,queryuc_ptr_3,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  pairpool,dynprogL,dynprogM,dynprogR,
					  oligoindices_minor,diagpool,cellpool,listpool,
					  hitlistpool);
#endif

      result = Stage3pair_resolve_multimapping(result,hitlistpool);
      /* result = Stage3pair_sort_distance(result); */
      debug16(printf("After removing overlaps, %d results\n",List_length(result)));

    } else if (true /*gmap_improvement_p == false*/) {
      debug16(printf("No GMAP improvement: Before removing overlaps, %d results\n",List_length(hitpairs)));
      result = Stage3pair_optimal_score(hitpairs,hitlistpool,querylength5,querylength3,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,hitlistpool,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/true);

#ifdef RESOLVE_INSIDE_GENERAL
      result = Stage3pair_resolve_insides(result,queryuc_ptr_5,queryuc_ptr_3,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  pairpool,dynprogL,dynprogM,dynprogR,
					  oligoindices_minor,diagpool,cellpool,listpool,
					  hitlistpool);
#endif

      result = Stage3pair_resolve_multimapping(result,hitlistpool);
      /* result = Stage3pair_sort_distance(result); */
      debug16(printf("After removing overlaps, %d results\n",List_length(result)));

    } else {
      debug16(printf("GMAP improvement: Before removing overlaps, %d results\n",List_length(hitpairs)));
      result = Stage3pair_optimal_score(hitpairs,hitlistpool,querylength5,querylength3,/*finalp*/false);
      result = Stage3pair_remove_overlaps(result,hitlistpool,/*translocp*/false,/*finalp*/false);
      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/false);
      result = Stage3pair_resolve_multimapping(result,hitlistpool);
      /* result = Stage3pair_sort_distance(result); */
      debug16(printf("After removing overlaps, %d results\n",List_length(result)));

#ifdef TODO
      result = align_pair_with_gmap(&(*final_pairtype),result,
				    queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
				    oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				    cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/true);
      if (Stage3pair_sense_consistent_p(result) == false) {
	debug16(printf("sense is not consistent\n"));
	result = align_pair_with_gmap(&(*final_pairtype),result,
				      queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
				      oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				      cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/true);
      }
#endif

      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,hitlistpool,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/true);

#ifdef RESOLVE_INSIDE_GENERAL
      result = Stage3pair_resolve_insides(result,queryuc_ptr_5,queryuc_ptr_3,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  pairpool,dynprogL,dynprogM,dynprogR,
					  oligoindices_minor,diagpool,cellpool,listpool,
					  hitlistpool);
#endif

      result = Stage3pair_resolve_multimapping(result,hitlistpool);
    }

  } else if (*final_pairtype == PAIRED_UNSPECIFIED) {
    /* Have paired results */
    debug16(printf("Have paired unspecified\n"));
    for (p = hitpairs; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    Hitlist_free(&hitpairs);

    for (p = conc_transloc; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    Hitlist_free(&conc_transloc);

    if (true /*gmap_improvement_p == false*/) {
      debug16(printf("No GMAP improvement: Before removing overlaps, %d results\n",List_length(samechr)));
      result = Stage3pair_optimal_score(samechr,hitlistpool,querylength5,querylength3,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,hitlistpool,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/true);
      result = Stage3pair_resolve_multimapping(result,hitlistpool);
    } else {
      debug16(printf("GMAP improvement: Before removing overlaps, %d results\n",List_length(samechr)));
      result = Stage3pair_optimal_score(samechr,hitlistpool,querylength5,querylength3,/*finalp*/false);
      result = Stage3pair_remove_overlaps(result,hitlistpool,/*translocp*/false,/*finalp*/false);
      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/false);
      result = Stage3pair_resolve_multimapping(result,hitlistpool);

#ifdef TODO
      result = align_pair_with_gmap(&(*final_pairtype),result,
				    queryuc_ptr_5,queryuc_ptr_3,querylength5,querylength3,
				    oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,
				    cutoff_level_5,cutoff_level_3,/*expect_concordant_p*/false);
      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/true);
      result = Stage3pair_remove_overlaps(result,hitlistpool,/*translocp*/false,/*finalp*/true);
      result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/true);
#endif

      if (Stage3pair_concordantp(result) == true) {
	debug16(printf("Found remaining concordant solution, so removing non-concordant ones\n"));
	*final_pairtype = CONCORDANT;
	result = Stage3pair_filter_nonconcordant(result,hitlistpool);

#ifdef RESOLVE_INSIDE_GENERAL
	result = Stage3pair_resolve_insides(result,queryuc_ptr_5,queryuc_ptr_3,
					    query5_compress_fwd,query5_compress_rev,
					    query3_compress_fwd,query3_compress_rev,
					    pairpool,dynprogL,dynprogM,dynprogR,
					    oligoindices_minor,diagpool,cellpool,listpool,
					    hitlistpool);
#endif

	debug16(printf("Concordant results: %d\n",List_length(result)));
      } else {
	*final_pairtype = PAIRED_UNSPECIFIED;
      }

      result = Stage3pair_resolve_multimapping(result,hitlistpool);
    }

    *final_pairtype = CONCORDANT; /* Because paired terminals are each concordant */


  } else if (*final_pairtype == CONCORDANT_TRANSLOCATIONS) {
    debug16(printf("Have %d concordant translocation results\n",List_length(conc_transloc)));
    for (p = hitpairs; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    Hitlist_free(&hitpairs);

    for (p = samechr; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    Hitlist_free(&samechr);

    result = Stage3pair_optimal_score(conc_transloc,hitlistpool,querylength5,querylength3,/*finalp*/true);
    result = Stage3pair_remove_overlaps(result,hitlistpool,/*translocp*/true,/*finalp*/true);
    result = Stage3pair_optimal_score(result,hitlistpool,querylength5,querylength3,/*finalp*/true);

#ifdef RESOLVE_INSIDE_GENERAL
    result = Stage3pair_resolve_insides(result,queryuc_ptr_5,queryuc_ptr_3,
					query5_compress_fwd,query5_compress_rev,
					query3_compress_fwd,query3_compress_rev,
					pairpool,dynprogL,dynprogM,dynprogR,
					oligoindices_minor,diagpool,cellpool,listpool,
					hitlistpool);
#endif

    result = Stage3pair_resolve_multimapping(result,hitlistpool);
    debug16(printf("Finally, have %d concordant translocation results\n",List_length(result)));

  } else {
    debug16(printf("Have unpaired results\n"));
    /* Need to free conc_transloc, since we can get here with multiple results */
    for (p = conc_transloc; p != NULL; p = List_next(p)) {
      stage3pair = (Stage3pair_T) List_head(p);
      Stage3pair_free(&stage3pair);
    }
    Hitlist_free(&conc_transloc);

    result = (List_T) NULL;
  }

  debug16(printf("After eval: %d in result, %d hitpairs, %d conc_transloc, %d samechr\n",
		 List_length(result),List_length(hitpairs),List_length(conc_transloc),List_length(samechr)));

  if (result == NULL) {
    debug16(printf("Entering Stage3end_optimal_score for 5' end\n"));
    singlehits5 = Stage3end_optimal_score(hits5,hitlistpool,querylength5,/*finalp*/true);
    /* singlehits5 = Stage3end_reject_trimlengths(singlehits5); */
    singlehits5 = Stage3end_linearize_5(singlehits5);
    singlehits5 = Stage3end_remove_overlaps(singlehits5,hitlistpool,/*finalp*/true);
    singlehits5 = Stage3end_optimal_score(singlehits5,hitlistpool,querylength5,/*finalp*/true);
    singlehits5 = Stage3end_resolve_multimapping(singlehits5,hitlistpool);

    debug16(printf("Entering Stage3end_optimal_score for 3' end\n"));
    singlehits3 = Stage3end_optimal_score(hits3,hitlistpool,querylength3,/*finalp*/true);
    /* singlehits3 = Stage3end_reject_trimlengths(singlehits3); */
    singlehits3 = Stage3end_linearize_3(singlehits3);
    singlehits3 = Stage3end_remove_overlaps(singlehits3,hitlistpool,/*finalp*/true);
    singlehits3 = Stage3end_optimal_score(singlehits3,hitlistpool,querylength3,/*finalp*/true);
    singlehits3 = Stage3end_resolve_multimapping(singlehits3,hitlistpool);

    debug16(printf("5' end has %d hits and 3' end has %d hits\n",
		   List_length(singlehits5),List_length(singlehits3)));

    if (List_length(singlehits5) == 1 && List_length(singlehits3) == 1) {
      hit5 = (Stage3end_T) List_head(singlehits5);
      hit3 = (Stage3end_T) List_head(singlehits3);
      if ((genestrand = Stage3end_genestrand(hit5)) == Stage3end_genestrand(hit3) &&
	  (pairtype = Stage3_determine_pairtype(hit5,hit3,/*stage3pair*/NULL)) != UNPAIRED) {
	/* Convert unpaired uniq to a paired uniq */
	debug16(printf("Converting unpaired uniq to paired uniq, with initial pairtype %s\n",Pairtype_string(pairtype)));
	if ((newpair = Stage3pair_new(hit5,hit3,genestrand,pairtype,
#if 0
				      queryuc_ptr_5,queryuc_ptr_3,
				      query5_compress_fwd,query5_compress_rev,
				      query3_compress_fwd,query3_compress_rev,
				      pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,
				      diagpool,cellpool,
#endif
				      listpool,/*expect_concordant_p*/pairtype == CONCORDANT ? true : false,
				      /*transcriptome_guided_p*/false)) != NULL) {
	  result = Hitlist_push(NULL,hitlistpool,(void *) newpair);
	  
	  *nhits5_primary = *nhits5_altloc = 0;
	  *nhits3_primary = *nhits3_altloc = 0;
	  *stage3array5 = *stage3array3 = (Stage3end_T *) NULL;
	  
	  if (Stage3pair_altlocp(newpair) == true) {
	    *npaths_primary = 0;
	    *npaths_altloc = 1;
	  } else {
	    *npaths_primary = 1;
	    *npaths_altloc = 0;
	  }
	  if (pairtype == CONCORDANT) {
	    debug16(printf("final pairtype is CONCORDANT\n"));
	    *final_pairtype = CONCORDANT;
	    
#ifdef RESOLVE_INSIDE_GENERAL
	    result = Stage3pair_resolve_insides(result,queryuc_ptr_5,queryuc_ptr_3,
						query5_compress_fwd,query5_compress_rev,
						query3_compress_fwd,query3_compress_rev,
						pairpool,dynprogL,dynprogM,dynprogR,
						oligoindices_minor,diagpool,cellpool,listpool,
						hitlistpool);
#endif
	    
	  } else {
	    debug16(printf("final pairtype is PAIRED_UNSPECIFIED\n"));
	    *final_pairtype = PAIRED_UNSPECIFIED;
	  }

	  stage3pairarray = (Stage3pair_T *) CALLOC_OUT(1,sizeof(Stage3pair_T));
	  stage3pairarray[0] = (Stage3pair_T) List_head(result);
	  Hitlist_free(&result);
	  
#if 0
	  Stage3pair_privatize(stage3pairarray,/*npairs*/1);
#endif
	  Stage3pair_eval_and_sort(/*npaths*/(*npaths_primary) + (*npaths_altloc),
				   &(*first_absmq),&(*second_absmq),stage3pairarray,
				   queryuc_ptr_5,queryuc_ptr_3,quality_string_5,quality_string_3);
	  
	  stage3list_gc(&singlehits3);
	  stage3list_gc(&singlehits5);
	  debug16(printf("1 Exiting consolidate_paired_results with final_pairtype %d\n",*final_pairtype));
	  return stage3pairarray;
	}
      }
    }

    /* Fall through: halfmapping or unpaired */
    *npaths_primary = *npaths_altloc = 0;
    *final_pairtype = UNPAIRED;
	  
    singlehits5 = Stage3end_filter(singlehits5,hitlistpool,max_mismatches_5,min_coverage_5);
    if (singlehits5 == NULL) {
      *nhits5_primary = *nhits5_altloc = 0;
      *stage3array5 = (Stage3end_T *) NULL;
    } else {
#if 0
      singlehits5 = Stage3end_unalias_circular(singlehits5);
#else
      singlehits5 = Stage3end_remove_circular_alias(singlehits5,hitlistpool); /* Contains a call to unalias_circular */
      singlehits5 = Stage3end_remove_duplicates(singlehits5,hitlistpool); /* Aliases can cause duplicates */
      Stage3end_count_hits(&(*nhits5_primary),&(*nhits5_altloc),singlehits5);
#endif
      *stage3array5 = (Stage3end_T *) List_to_array_out(singlehits5,NULL); Hitlist_free(&singlehits5); /* Return value */
    }

    singlehits3 = Stage3end_filter(singlehits3,hitlistpool,max_mismatches_3,min_coverage_3);
    if (singlehits3 == NULL) {
      *nhits3_primary = *nhits3_altloc = 0;
      *stage3array3 = (Stage3end_T *) NULL;
    } else {
#if 0
      singlehits3 = Stage3end_unalias_circular(singlehits3);
#else
      singlehits3 = Stage3end_remove_circular_alias(singlehits3,hitlistpool); /* Contains a call to unalias_circular */
      singlehits3 = Stage3end_remove_duplicates(singlehits3,hitlistpool); /* Aliases can cause duplicates */
      Stage3end_count_hits(&(*nhits3_primary),&(*nhits3_altloc),singlehits3);
#endif
      *stage3array3 = (Stage3end_T *) List_to_array_out(singlehits3,NULL); Hitlist_free(&singlehits3); /* Return value */
    }

    if ((*nhits5_primary) + (*nhits5_altloc) > 0) {
      if ((*nhits3_primary) + (*nhits3_altloc) == 1) {
      /* Use single 3' hit to guide sorting of multiple 5' hits */
        *stage3array5 = Stage3end_eval_and_sort_guided((*nhits5_primary) + (*nhits5_altloc),
   	                                               &(*first_absmq5),&(*second_absmq5),/*guide*/(*stage3array3)[0],
						       *stage3array5,queryuc_ptr_5,quality_string_5,/*displayp*/true);
      } else {
        *stage3array5 = Stage3end_eval_and_sort((*nhits5_primary) + (*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
 						*stage3array5,queryuc_ptr_5,quality_string_5,/*displayp*/true);
      }
    }

    if ((*nhits3_primary) + (*nhits3_altloc) > 0) {
      if ((*nhits5_primary) + (*nhits5_altloc) == 1) {
	/* Use single 5' hit to guide sorting of multiple 3' hits */
        *stage3array3 = Stage3end_eval_and_sort_guided((*nhits3_primary) + (*nhits3_altloc),
                                                       &(*first_absmq3),&(*second_absmq3),/*guide*/(*stage3array5)[0],
						       *stage3array3,queryuc_ptr_3,quality_string_3,/*displayp*/true);
      } else {
        *stage3array3 = Stage3end_eval_and_sort((*nhits3_primary) + (*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
						*stage3array3,queryuc_ptr_3,quality_string_3,/*displayp*/true);
      }
    }
    debug16(printf("Result is NULL, and we have %d hits on 5' end and %d hits on 3' end\n",
		   (*nhits5_primary) + (*nhits5_altloc),(*nhits3_primary) + (*nhits3_altloc)));
    debug16(printf("2 Exiting consolidate_paired_results with final_pairtype %d\n",*final_pairtype));
    return (Stage3pair_T *) NULL;

  } else {
    debug16(printf("final pairtype is %d\n",*final_pairtype));
    debug16(printf("Result is not NULL (%d paths), and we fall through to concordant, paired, or transloc pairs\n",
		   List_length(result)));

    result = Stage3pair_filter(result,hitlistpool,max_mismatches_5,max_mismatches_3,min_coverage_5,min_coverage_3);
    if (result == NULL) {
      *npaths_primary = *npaths_altloc = 0;
      stage3list_gc(&hits3);
      stage3list_gc(&hits5);
      *nhits5_primary = *nhits5_altloc = 0;
      *nhits3_primary = *nhits3_altloc = 0;
      *stage3array5 = *stage3array3 = (Stage3end_T *) NULL;
      debug16(printf("3 Exiting consolidate_paired_results with final_pairtype %d\n",*final_pairtype));
      return (Stage3pair_T *) NULL;

    } else {
      /* result != NULL */
      /* Concordant, paired, or transloc pairs found.  Remove single hits. */
      Stage3pair_count_hits(&(*npaths_primary),&(*npaths_altloc),result);
      stage3pairarray = (Stage3pair_T *) List_to_array_out(result,NULL); Hitlist_free(&result); /* Return value */
#if 0
      Stage3pair_privatize(stage3pairarray,/*npaths*/(*npaths_primary) + (*npaths_altloc));
#endif
      Stage3pair_eval_and_sort(/*npaths*/(*npaths_primary) + (*npaths_altloc),
			       &(*first_absmq),&(*second_absmq),
			       stage3pairarray,queryuc_ptr_5,queryuc_ptr_3,quality_string_5,quality_string_3);
      stage3list_gc(&hits3);
      stage3list_gc(&hits5);

      *nhits5_primary = *nhits5_altloc = 0;
      *nhits3_primary = *nhits3_altloc = 0;
      *stage3array5 = *stage3array3 = (Stage3end_T *) NULL;
      debug16(printf("4 Exiting consolidate_paired_results with final_pairtype %d\n",*final_pairtype));

      return stage3pairarray;
    }
  }
}


Stage3pair_T *
Stage1_paired_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
		    Stage3end_T **stage3array5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
		    Stage3end_T **stage3array3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
		    Shortread_T queryseq5, Shortread_T queryseq3, Chrpos_T pairmax_linear,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool) {
  Stage3pair_T *stage3pairarray;
  bool abort_pairing_p, geneplus_abort_pairing_p, geneminus_abort_pairing_p;
  List_T hitpairs, geneplus_hitpairs, geneminus_hitpairs;
  List_T hits5, hits3, geneplus_hits5, geneplus_hits3, geneminus_hits5, geneminus_hits3;
  List_T samechr, geneplus_samechr, geneminus_samechr;
  List_T conc_transloc, geneplus_conc_transloc, geneminus_conc_transloc;

  int querylength5, querylength3;
  int max_mismatches_5, max_mismatches_3;
  int min_coverage_5, min_coverage_3;
  char *queryuc_ptr_5, *queryuc_ptr_3, *queryrc5, *queryrc3, *quality_string_5, *quality_string_3;
  Compress_T query5_compress_fwd, query5_compress_rev, query3_compress_fwd, query3_compress_rev;


  querylength5 = Shortread_fulllength(queryseq5);
  querylength3 = Shortread_fulllength(queryseq3);

  /* Previously used user_maxlevel_float in searching, now just for filtering */
  if (user_maxlevel_float < 0.0) {
    max_mismatches_5 = querylength5;
    max_mismatches_3 = querylength3;
  } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
    max_mismatches_5 = (int) rint(user_maxlevel_float * (double) querylength5);
    max_mismatches_3 = (int) rint(user_maxlevel_float * (double) querylength3);
  } else {
    max_mismatches_5 = max_mismatches_3 = (int) user_maxlevel_float;
  }

  if (user_mincoverage_float < 0.0) {
    min_coverage_5 = min_coverage_3 = 0;
  } else if (user_mincoverage_float > 0.0 && user_mincoverage_float < 1.0) {
    min_coverage_5 = (int) rint(user_mincoverage_float * (double) querylength5);
    min_coverage_3 = (int) rint(user_mincoverage_float * (double) querylength3);
  } else {
    min_coverage_5 = min_coverage_3 = (int) user_mincoverage_float;
  }

  queryuc_ptr_5 = Shortread_fullpointer_uc(queryseq5);
  queryuc_ptr_3 = Shortread_fullpointer_uc(queryseq3);
  quality_string_5 = Shortread_quality_string(queryseq5);
  quality_string_3 = Shortread_quality_string(queryseq3);

  queryrc5 = (char *) MALLOC((querylength5+1)*sizeof(char));
  queryrc3 = (char *) MALLOC((querylength3+1)*sizeof(char));
  make_complement_buffered(queryrc5,queryuc_ptr_5,querylength5);
  make_complement_buffered(queryrc3,queryuc_ptr_3,querylength3);

  query5_compress_fwd = Compress_new_fwd(queryuc_ptr_5,querylength5);
  query5_compress_rev = Compress_new_rev(queryuc_ptr_5,querylength5);
  query3_compress_fwd = Compress_new_fwd(queryuc_ptr_3,querylength3);
  query3_compress_rev = Compress_new_rev(queryuc_ptr_3,querylength3);


  if (mode == STANDARD || mode == CMET_STRANDED || mode == ATOI_STRANDED || mode == TTOC_STRANDED) {
    hitpairs = paired_read(&abort_pairing_p,&hits5,&hits3,&samechr,&conc_transloc,
			   queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
			   query5_compress_fwd,query5_compress_rev,
			   query3_compress_fwd,query3_compress_rev,
			   /*genestrand*/0,pairmax_linear,
			   intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool);
    if (abort_pairing_p == true) {
      /* Too many concordant.  Could put into a separate split output */
      debug(printf("abort_pairing_p is true\n"));
      *final_pairtype = CONCORDANT;
    } else {
      *final_pairtype = CONCORDANT;
    }

    stage3pairarray =
      consolidate_paired_results(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),&(*final_pairtype),
				 &(*stage3array5),&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
				 &(*stage3array3),&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
				 hitpairs,samechr,conc_transloc,hits5,hits3,
				 
				 queryuc_ptr_5,querylength5,queryuc_ptr_3,querylength3,
				 quality_string_5,quality_string_3,
				 /*query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,*/
				 
				 max_mismatches_5,max_mismatches_3,min_coverage_5,min_coverage_3,
				 /*oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,*/
				 listpool,hitlistpool);

  } else if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
    geneplus_hitpairs = paired_read(&geneplus_abort_pairing_p,&geneplus_hits5,&geneplus_hits3,
				    &geneplus_samechr,&geneplus_conc_transloc,
				    queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
				    query5_compress_fwd,query5_compress_rev,
				    query3_compress_fwd,query3_compress_rev,
				    /*genestrand*/+1,pairmax_linear,
				    intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool);

    geneminus_hitpairs = paired_read(&geneminus_abort_pairing_p,&geneminus_hits5,&geneminus_hits3,
				     &geneminus_samechr,&geneminus_conc_transloc,
				     queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
				     query5_compress_fwd,query5_compress_rev,
				     query3_compress_fwd,query3_compress_rev,
				     /*genestrand*/+2,pairmax_linear,
				     intlistpool,univcoordlistpool,listpool,univdiagpool,hitlistpool);

    if (geneplus_abort_pairing_p == true || geneminus_abort_pairing_p == true) {
      /* Too many concordant.  Could put into a separate split output */
      debug(printf("abort_pairing_p is true\n"));
      *final_pairtype = CONCORDANT;
    } else {
      *final_pairtype = CONCORDANT;
    }

    hits5 = List_append(geneplus_hits5,geneminus_hits5);
    hits3 = List_append(geneplus_hits3,geneminus_hits3);
    hitpairs = List_append(geneplus_hitpairs,geneminus_hitpairs);
    samechr = List_append(geneplus_samechr,geneminus_samechr);
    conc_transloc = List_append(geneplus_conc_transloc,geneminus_conc_transloc);
    
    stage3pairarray =
      consolidate_paired_results(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),&(*final_pairtype),
				 &(*stage3array5),&(*nhits5_primary),&(*nhits5_altloc),&(*first_absmq5),&(*second_absmq5),
				 &(*stage3array3),&(*nhits3_primary),&(*nhits3_altloc),&(*first_absmq3),&(*second_absmq3),
				 hitpairs,samechr,conc_transloc,hits5,hits3,
				 
				 queryuc_ptr_5,querylength5,queryuc_ptr_3,querylength3,
				 quality_string_5,quality_string_3,
				 /*query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,*/
				 
				 max_mismatches_5,max_mismatches_3,min_coverage_5,min_coverage_3,
				 /*oligoindices_minor,pairpool,diagpool,cellpool,dynprogL,dynprogM,dynprogR,*/
				 listpool,hitlistpool);

  } else {
    fprintf(stderr,"Do not recognize mode %d\n",mode);
    abort();
  }

  Compress_free(&query5_compress_fwd);
  Compress_free(&query5_compress_rev);
  Compress_free(&query3_compress_fwd);
  Compress_free(&query3_compress_rev);

  FREE(queryrc3);
  FREE(queryrc5);

  debug(printf("Returning with final_pairtype %d\n",*final_pairtype));
  return stage3pairarray;
}


void
Stage1hr_setup (Univ_IIT_T transcript_iit_in, Transcriptome_T transcriptome_in, Genome_T transcriptomebits_in,
		bool use_only_transcriptome_p_in,

		Indexdb_T indexdb_fwd_in, Indexdb_T indexdb_rev_in, Indexdb_T indexdb_tr_in, Localdb_T localdb_in,

		int index1part_tr_in, int index1part_in, int index1interval_in, 
		double user_maxlevel_float_in, double user_mincoverage_float_in,

		Univ_IIT_T chromosome_iit_in, int nchromosomes_in,
		Genome_T genomecomp_in, Genome_T genomebits_in, Genome_T genomebits_alt_in,
		Mode_T mode_in, int maxpaths_search_in, int maxpaths_report_in,

		bool find_dna_chimeras_p_in, bool distances_observed_p_in, int subopt_levels_in,
		int max_middle_insertions, int max_middle_deletions,

		bool novelsplicingp_in, Chrpos_T shortsplicedist_in,
		Chrpos_T min_intronlength_in, Chrpos_T expected_pairlength_in, Chrpos_T pairlength_deviation_in,

		int distantsplicing_penalty_in,
		int min_distantsplicing_end_matches_in, int min_distantsplicing_identity_in) {

  use_only_transcriptome_p = use_only_transcriptome_p_in;

  transcript_iit = transcript_iit_in;
  transcriptome = transcriptome_in;
  transcriptomebits = transcriptomebits_in;

  indexdb_fwd = indexdb_fwd_in;
  indexdb_rev = indexdb_rev_in;
  indexdb_tr = indexdb_tr_in;
  localdb = localdb_in;

  index1part_tr = index1part_tr_in;
  index1part = index1part_in;
  index1interval = index1interval_in;

  user_maxlevel_float = user_maxlevel_float_in;
  user_mincoverage_float = user_mincoverage_float_in;

  chromosome_iit = chromosome_iit_in;
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  nchromosomes = nchromosomes_in;

  genomecomp = genomecomp_in;
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

#ifdef HAVE_64_BIT
  leftreadshift = 64 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#else
  leftreadshift = 32 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#endif

  poly_a = POLY_A & oligobase_mask;
  poly_t = POLY_T & oligobase_mask;

  mode = mode_in;
  maxpaths_search = maxpaths_search_in;
  maxpaths_report = maxpaths_report_in;

  find_dna_chimeras_p = find_dna_chimeras_p_in;
  distances_observed_p = distances_observed_p_in;

  subopt_levels = subopt_levels_in;

  novelsplicingp = novelsplicingp_in;
  shortsplicedist = shortsplicedist_in;

  overall_max_distance = shortsplicedist;
  if (max_middle_deletions > (int) overall_max_distance) {
    overall_max_distance = max_middle_deletions;
  }
  if (max_middle_insertions > (int) overall_max_distance) {
    overall_max_distance = max_middle_insertions;
  }

  min_intronlength = min_intronlength_in;
  expected_pairlength = expected_pairlength_in;
  pairlength_deviation = pairlength_deviation_in;

  distantsplicing_penalty = distantsplicing_penalty_in;
  min_distantsplicing_end_matches = min_distantsplicing_end_matches_in;
  min_distantsplicing_identity = min_distantsplicing_identity_in;

  if (genomebits_alt_in != NULL) {
    snpp = true;
  } else {
    snpp = false;
  }

  return;
}
