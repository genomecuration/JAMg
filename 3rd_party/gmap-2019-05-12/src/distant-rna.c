static char rcsid[] = "$Id: distant-rna.c 218675 2019-03-16 01:25:48Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "distant-rna.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For rint */
#include "mem.h"
#include "assert.h"
#include "types.h"

#include "genomicpos.h"
#include "univdiagdef.h"
#include "substring.h"

#include "genome_sites.h"
#include "genome128_hr.h"
#include "maxent.h"
#include "maxent_hr.h"

#include "extension-search.h"	/* For handling Elt_T objects */


/* Originally allowed only 1, to print only unique translocations.
   But need to allow enough to avoid missing some translocations. */
/* For transcript splicing, need to increase MAXCHIMERAPATHS */
/* #define MAXCHIMERAPATHS 100 */
#define MAXCHIMERAPATHS 10000
#define MIN_SPLICE_EXON_LENGTH 0
#define MIN_FRAGLENGTH 25

static Univ_IIT_T chromosome_iit;
static int circular_typeint;

static Genome_T genomebits;
static Genome_T genomebits_alt;

static int index1part = 15;
static int index1interval = 3;

static int subopt_levels;

static Chrpos_T shortsplicedist;
static int min_distantsplicing_end_matches;
static int min_distantsplicing_identity;

static int localsplicing_penalty;
static int distantsplicing_penalty;


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG4E
#define debug4e(x) x
#else
#define debug4e(x)
#endif

#ifdef DEBUG4L
#define debug4l(x) x
#else
#define debug4l(x)
#endif

#ifdef DEBUG4LD
#define debug4ld(x) x
#else
#define debug4ld(x)
#endif




/* Do not compare against true or false */
/* Moderate criterion */
static int
sufficient_splice_prob_distant (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < min_distantsplicing_end_matches) {
    return 0;
  } else if (support < 30) {
    return (spliceprob > 0.95);
  } else if (support < 35) {
    return (spliceprob > 0.90);
  } else if (support < 40) {
    return (spliceprob > 0.85);
  } else {
    return (spliceprob > 0.70);
  }
}

#if 0
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
#endif


/* Streamlined version */
/* Produces lists of distant_donors and distant_acceptors that are substrings */
/* Note: Call to Genome_donor_positions and similar functions can
   yield positions at (splice_pos_start - 1) through (splice_pos_end +
   1).  Therefore, need to have splice_pos_start > 0 and
   splice_pos_end < querylength. */
/* TODO: Change to lists of Stage3end_T objects, including GMAP.
   Change definition of a chimera to be two Stage3end_T objects, instead
   of two substrings. */
/* mismatch_positions_alloc and positions_alloc are MALLOC((querylength+1)*sizeof(int)) */
static void
find_spliceends_rna_plus (List_T **distant_donors, List_T **distant_antidonors,
			  List_T **distant_acceptors, List_T **distant_antiacceptors,
			  List_T *distant_startfrags, List_T *distant_endfrags,

			  List_T elt_set,
#ifdef DEBUG4E
			  char *queryptr,
#endif
			  int querylength, int query_lastpos,
			  int *mismatch_positions_alloc, int *positions_alloc,
			  Compress_T query_compress, Listpool_T listpool,
			  int max_mismatches_allowed, int genestrand) {
#ifdef DEBUG4E
  char *gbuffer;
#endif

  Elt_T elt;
  List_T p;
  int k;

  Substring_T substring;
  Univcoord_T segment_left;
  int nmismatches, nmismatches_end_trim, i;
  int splice_pos_qstart, splice_pos_qend, splice_pos;
  bool splice_querystart_p, splice_queryend_p;
  Splicetype_T splicetype_querystart, splicetype_queryend;
  double ambig_prob_querystart, ambig_prob_queryend, prob;

  int total_nmismatches_splice, nmismatches_frag;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *mismatch_positions_splice;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;


  debug4e(printf("Entering find_spliceends_rna_plus with %d elts\n",List_length(elt_set)));

  for (p = elt_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);

    debug4e(printf("Comparing elt->qend %d < query_lastpos %d\n",elt->qend,query_lastpos));
    if (elt->qend < query_lastpos /*&& (elt->qstart < index1part || segment->spliceable_low_p == true)*/) {
      /* Find splices on genomic right */
      for (k = 0; k < elt->n_all_diagonals; k++) {
	segment_left = elt->all_diagonals[k];
	chrnum = Univ_IIT_get_one(chromosome_iit,segment_left,segment_left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

	debug4e(printf("find_spliceends_rna_plus: Checking up to %d mismatches at diagonal %llu (elt qstart..qend %d..%d)\n",
		       max_mismatches_allowed,(unsigned long long) segment_left,elt->qstart,elt->qend));
	debug4e(
		gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		printf("genome 0..: %s\n",gbuffer);
		printf("query  0..: %s\n",queryptr);
		FREE(gbuffer);
		);

	/* Trim at end (away from splice, or qstart) */
	debug4e(printf("Trimming at end: %d..%d\n",0,elt->qstart));

	/* SENSE FORWARD */
	splice_querystart_p = Substring_qstart_trim(&splice_pos_qstart,&splicetype_querystart,&ambig_prob_querystart,
						    segment_left,/*qend: was elt->qstart*/elt->qend,
						    /*plusp*/true,genestrand,mismatch_positions_alloc,query_compress,chroffset,
						    /*sensedir*/SENSE_FORWARD);
	debug4e(printf("sense: splice_pos_qstart %d with splice_querystart_p %d\n",splice_pos_qstart,splice_querystart_p));
	if (splice_pos_qstart >= 0 && (splice_pos_qstart < 6 || splice_querystart_p == true)) {
	  /* genomic left anchor or middle anchor */
	  debug4e(printf("Searching genomic right (calling Genome_mismatches_left): plus genomic left anchor or middle anchor: %d..%d\n",
			 elt->qstart,querylength));
	  mismatch_positions_splice = mismatch_positions_alloc;
	  total_nmismatches_splice = Genome_mismatches_left(mismatch_positions_alloc,max_mismatches_allowed,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    /*left*/segment_left,/*pos5*/splice_pos_qstart,/*pos3*/querylength,
							    /*plusp*/true,genestrand);
	  debug4e(
		  printf("%d mismatches at splice (%d allowed) at:",total_nmismatches_splice,max_mismatches_allowed);
		  for (i = 0; i <= total_nmismatches_splice; i++) {
		    printf(" %d",mismatch_positions_splice[i]);
		  }
		  printf("\n");
		  );

	  if (total_nmismatches_splice > max_mismatches_allowed) {
	    splice_pos_qend = mismatch_positions_splice[max_mismatches_allowed];
	    nmismatches_frag = max_mismatches_allowed; /* should be total_nmismatches_splice - 1 */
	  } else {
	    splice_pos_qend = mismatch_positions_splice[total_nmismatches_splice];
	    nmismatches_frag = total_nmismatches_splice;
	  }

	  if (splice_pos_qend - splice_pos_qstart >= MIN_FRAGLENGTH) {
	    
	    if ((substring = Substring_new_startfrag(nmismatches_frag,
						     /*querystart*/splice_pos_qstart,/*queryend*/splice_pos_qend,
						     /*left*/segment_left,query_compress,
						     querylength,/*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	      debug4e(printf("=> plus startfrag: at %d, %d..%d\n",Substring_siteN_pos(substring),
			     Substring_querystart(substring),Substring_queryend(substring)));
	      debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	      *distant_startfrags = Listpool_push(*distant_startfrags,listpool,(void *) substring);
	    }

	    if (segment_left + splice_pos_qstart >= DONOR_MODEL_LEFT_MARGIN) {
	      /* (1) Originally on plus strand.  No complement */
 	      debug4e(printf("Search for donor splice sites from %d up to %d\n",splice_pos_qstart,splice_pos_qend));
	      donori_nsites = Genome_donor_positions_novel(positions_alloc,segment_left,splice_pos_qstart,splice_pos_qend);
	      donori_positions = positions_alloc;
	      debug4e(
		      printf("Donor dinucleotides:");
		      for (i = 0; i < donori_nsites; i++) {
			printf(" %d",donori_positions[i]);
		      }
		      printf("\n");
		      );
	      
	      i = 0;
	      nmismatches = 0;
	      while (i < donori_nsites && nmismatches <= max_mismatches_allowed) {
		splice_pos = donori_positions[i];
		while (nmismatches <= total_nmismatches_splice &&
		       mismatch_positions_splice[nmismatches] < splice_pos) { /* Changed from <= to < */
		  debug4e(printf("  mismatch at %d\n",mismatch_positions_splice[nmismatches]));
		  nmismatches++;
		}
		debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
		
		if (splice_pos - splice_pos_qstart >= MIN_SPLICE_EXON_LENGTH) {
		  assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,segment_left,
									  /*pos5*/splice_pos_qstart,/*pos3*/splice_pos,/*plusp*/true,genestrand));
		  if (segment_left + querylength < chrhigh) {
		    prob = Maxent_hr_donor_prob(segment_left + splice_pos,chroffset);
		    debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
				   splice_pos,nmismatches,prob,sufficient_splice_prob_distant(splice_pos,nmismatches,prob)));
		    if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		      debug4e(printf("Novel donor for segment at %llu, splice_pos %d (%d mismatches)\n",
				     (unsigned long long) segment_left,splice_pos,nmismatches));

		      if ((substring = Substring_new_donor(nmismatches,/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/-1,
							   /*querystart*/splice_pos_qstart,/*queryend*/splice_pos,/*sitepos*/splice_pos,prob,
							   splice_querystart_p,splicetype_querystart,ambig_prob_querystart,
							   /*splice_queryend_p*/true,/*splicetype_querystart*/DONOR,/*ambig_prob_queryend*/prob,
							   /*left*/segment_left,query_compress,querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
							   chrnum,chroffset,chrhigh,chrlength)) != NULL) {
			debug4e(printf("=> plus donor: %f at %d (%d mismatches) %d..%d\n",
				       prob,Substring_siteD_pos(substring),nmismatches,
				       Substring_querystart(substring),Substring_queryend(substring)));
			debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
			(*distant_donors)[nmismatches] = Listpool_push((*distant_donors)[nmismatches],
								       listpool,(void *) substring);
		      }
		    }
		  }
		}

		i++;
	      }
	    }
	  }
	}

	/* SENSE ANTI */
	splice_querystart_p = Substring_qstart_trim(&splice_pos_qstart,&splicetype_querystart,&ambig_prob_querystart,
						    segment_left,/*qend: was elt->qstart*/elt->qend,
						    /*plusp*/true,genestrand,mismatch_positions_alloc,query_compress,chroffset,
						    /*sensedir*/SENSE_ANTI);
	debug4e(printf("antisense: splice_pos_qstart %d with splice_querystart_p %d\n",splice_pos_qstart,splice_querystart_p));
	if (splice_pos_qstart >= 0 && (splice_pos_qstart < 6 || splice_querystart_p == true)) {
	  /* genomic left anchor or middle anchor */
	  debug4e(printf("Searching genomic right (calling Genome_mismatches_left): plus genomic left anchor or middle anchor: %d..%d\n",
			 elt->qstart,querylength));
	  mismatch_positions_splice = mismatch_positions_alloc;
	  total_nmismatches_splice = Genome_mismatches_left(mismatch_positions_alloc,max_mismatches_allowed,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    /*left*/segment_left,/*pos5*/splice_pos_qstart,/*pos3*/querylength,
							    /*plusp*/true,genestrand);
	  debug4e(
		  printf("%d mismatches at splice (%d allowed) at:",total_nmismatches_splice,max_mismatches_allowed);
		  for (i = 0; i <= total_nmismatches_splice; i++) {
		    printf(" %d",mismatch_positions_splice[i]);
		  }
		  printf("\n");
		  );
	  
	  if (total_nmismatches_splice > max_mismatches_allowed) {
	    splice_pos_qend = mismatch_positions_splice[max_mismatches_allowed];
	    nmismatches_frag = max_mismatches_allowed; /* should be total_nmismatches_splice - 1 */
	  } else {
	    splice_pos_qend = mismatch_positions_splice[total_nmismatches_splice];
	    nmismatches_frag = total_nmismatches_splice;
	  }

	  if (splice_pos_qend - splice_pos_qstart >= MIN_FRAGLENGTH) {
	    if ((substring = Substring_new_startfrag(nmismatches_frag,
						     /*querystart*/splice_pos_qstart,/*queryend*/splice_pos_qend,
						     /*left*/segment_left,query_compress,
						     querylength,/*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	      debug4e(printf("=> plus startfrag: at %d, %d..%d\n",Substring_siteN_pos(substring),
			     Substring_querystart(substring),Substring_queryend(substring)));
	      debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	      *distant_startfrags = Listpool_push(*distant_startfrags,listpool,(void *) substring);
	    }

	    if (segment_left + splice_pos_qstart >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	      /* (2) Splicing originally on minus strand.  Complement */
	      debug4e(printf("Search for antiacceptor splice sites from %d up to %d\n",splice_pos_qstart,splice_pos_qend));
	      antiacceptori_nsites = Genome_antiacceptor_positions_novel(positions_alloc,segment_left,splice_pos_qstart,splice_pos_qend);
	      antiacceptori_positions = positions_alloc;
	      debug4e(
		      printf("Antiacceptor dinucleotides:");
		      for (i = 0; i < antiacceptori_nsites; i++) {
			printf(" %d",antiacceptori_positions[i]);
		      }
		      printf("\n");
		      );
	      
	      i = 0;
	      nmismatches = 0;
	      while (i < antiacceptori_nsites && nmismatches <= max_mismatches_allowed) {
		splice_pos = antiacceptori_positions[i];
		while (nmismatches <= total_nmismatches_splice &&
		       mismatch_positions_splice[nmismatches] < splice_pos) { /* Changed from <= to < */
		  debug4e(printf("  mismatch at %d\n",mismatch_positions_splice[nmismatches]));
		  nmismatches++;
		}
		debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
		
		if (splice_pos - splice_pos_qstart >= MIN_SPLICE_EXON_LENGTH) {
		  assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,segment_left,
									  /*pos5*/splice_pos_qstart,/*pos3*/splice_pos,/*plusp*/true,genestrand));
		  if (segment_left + querylength < chrhigh) {
		    prob = Maxent_hr_antiacceptor_prob(segment_left + splice_pos,chroffset);
		    debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
				   splice_pos,nmismatches,prob,sufficient_splice_prob_distant(splice_pos,nmismatches,prob)));
		    if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		      debug4e(printf("Novel antiacceptor for segment at %llu, splice_pos %d (%d mismatches)\n",
				     (unsigned long long) segment_left,splice_pos,nmismatches));
		      if ((substring = Substring_new_acceptor(nmismatches,/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/-1,
							      /*querystart*/splice_pos_qstart,/*queryend*/splice_pos,/*sitepos*/splice_pos,prob,
							      splice_querystart_p,splicetype_querystart,ambig_prob_querystart,
							      /*splice_queryend_p*/true,/*splicetype_querystart*/ANTIACCEPTOR,/*ambig_prob_queryend*/prob,
							      /*left*/segment_left,query_compress,querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
							      chrnum,chroffset,chrhigh,chrlength)) != NULL) {
			debug4e(printf("=> plus antiacceptor : %f at %d (%d mismatches) %d..%d\n",
				       prob,Substring_siteA_pos(substring),nmismatches,
				       Substring_querystart(substring),Substring_queryend(substring)));
			debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
			(*distant_antiacceptors)[nmismatches] = Listpool_push((*distant_antiacceptors)[nmismatches],
									      listpool,(void *) substring);
		      }
		    }
		  }
		}
		
		i++;
	      }
	    }
	  }
	}

      }
    }


    debug4e(printf("Comparing elt->qstart %d > index1part %d\n",elt->qstart,index1part));
    if (elt->qstart > index1part /*&& (elt->qend > query_lastpos || segment->spliceable_high_p == true)*/) {
      /* Find splices on genomic left */
      for (k = 0; k < elt->n_all_diagonals; k++) {
	segment_left = elt->all_diagonals[k];
	chrnum = Univ_IIT_get_one(chromosome_iit,segment_left,segment_left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	
	debug4e(printf("find_spliceends_rna_plus: Checking up to %d mismatches at diagonal %llu (elt qstart..qend %d..%d)\n",
		       max_mismatches_allowed,(unsigned long long) segment_left,elt->qstart,elt->qend));
	debug4e(
		gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		printf("genome 0..: %s\n",gbuffer);
		printf("query  0..: %s\n",queryptr);
		FREE(gbuffer);
		);

	/* Trim at end (away from splice, or qend) */
	debug4e(printf("Trimming at end: %d..%d\n",elt->qend,querylength));

	/* SENSE FORWARD */
	splice_queryend_p = Substring_qend_trim(&splice_pos_qend,&splicetype_queryend,&ambig_prob_queryend,
						segment_left,/*qstart:was elt->qend*/elt->qstart,querylength,
						/*plusp*/true,genestrand,mismatch_positions_alloc,query_compress,chroffset,
						/*sensedir*/SENSE_FORWARD);
	debug4e(printf("sense: splice_pos_qend %d with splice_queryend_p %d\n",splice_pos_qend,splice_queryend_p));
	if (splice_pos_qend >= 0 && (splice_pos_qend > querylength - 6 || splice_queryend_p == true)) {
	  /* genomic right anchor or middle anchor */
	  debug4e(printf("Searching genomic left (calling Genome_mismatches_right): plus genomic right anchor or middle anchor: %d..%d\n",
			 0,splice_pos_qend));
	  mismatch_positions_splice = mismatch_positions_alloc;
	  total_nmismatches_splice = Genome_mismatches_right(mismatch_positions_alloc,max_mismatches_allowed,
							     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							     /*left*/segment_left,/*pos5*/0,/*pos3*/splice_pos_qend,
							     /*plusp*/true,genestrand);
	  debug4e(
		  printf("%d mismatches at splice (%d allowed) at:",total_nmismatches_splice,max_mismatches_allowed);
		  for (i = 0; i <= total_nmismatches_splice; i++) {
		    printf(" %d",mismatch_positions_splice[i]);
		  }
		  printf("\n");
		  );
	  
	  if (total_nmismatches_splice > max_mismatches_allowed) {
	    splice_pos_qstart = mismatch_positions_splice[max_mismatches_allowed] + 1;
	    nmismatches_frag = max_mismatches_allowed; /* should be total_nmismatches_splice - 1 */
	  } else {
	    splice_pos_qstart = mismatch_positions_splice[total_nmismatches_splice] + 1;
	    nmismatches_frag = total_nmismatches_splice;
	  }
	  
	  if (splice_pos_qend - splice_pos_qstart >= MIN_FRAGLENGTH) {
	    if ((substring = Substring_new_endfrag(nmismatches_frag,
						   /*querystart*/splice_pos_qstart,/*queryend*/splice_pos_qend,
						   /*left*/segment_left,query_compress,
						   querylength,/*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	      debug4e(printf("=> plus endfrag: at %d, %d..%d\n",Substring_siteN_pos(substring),
			     Substring_querystart(substring),Substring_queryend(substring)));
	      debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	      *distant_endfrags = Listpool_push(*distant_endfrags,listpool,(void *) substring);
	    }

	    if (segment_left + splice_pos_qend >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	      /* (1) Splicing originally on plus strand.  No complement. */
	      debug4e(printf("Search for acceptor splice sites from %d down to %d\n",splice_pos_qend,splice_pos_qstart));
	      acceptorj_nsites = Genome_acceptor_positions_novel(positions_alloc,segment_left,splice_pos_qstart,splice_pos_qend);
	      acceptorj_positions = positions_alloc;
	      debug4e(
		      printf("Acceptor dinucleotides:");
		      for (i = 0; i < acceptorj_nsites; i++) {
			printf(" %d",acceptorj_positions[i]);
		      }
		      printf("\n");
		      );
	    
	      i = acceptorj_nsites - 1;
	      nmismatches = 0;
	      while (i >= 0 && nmismatches <= max_mismatches_allowed) {
		splice_pos = acceptorj_positions[i];
		while (nmismatches <= total_nmismatches_splice &&
		       mismatch_positions_splice[nmismatches] >= splice_pos) { /* Must be >= */
		  debug4e(printf("  mismatch at %d\n",mismatch_positions_splice[nmismatches]));
		  nmismatches++;
		}
		debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
	      
		if (splice_pos_qend - splice_pos >= MIN_SPLICE_EXON_LENGTH) {
		  assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,segment_left,
									  /*pos5*/splice_pos,/*pos3*/splice_pos_qend,/*plusp*/true,genestrand));
		  if (segment_left + querylength < chrhigh) {
		    prob = Maxent_hr_acceptor_prob(segment_left + splice_pos,chroffset);
		    debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
				   splice_pos,nmismatches,prob,sufficient_splice_prob_distant(querylength - splice_pos,nmismatches,prob)));
		    if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		      debug4e(printf("Novel acceptor for segment at %llu, splice_pos %d (%d mismatches)\n",
				     (unsigned long long) segment_left,splice_pos,nmismatches));
		      if ((substring = Substring_new_acceptor(nmismatches,/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/-1,
							      /*querystart*/splice_pos,/*queryend*/splice_pos_qend,/*sitepos*/splice_pos,prob,
							      /*splice_querystart_p*/true,/*splicetype_querystart*/ACCEPTOR,/*ambig_prob_querystart*/prob,
							      splice_queryend_p,splicetype_queryend,ambig_prob_queryend,
							      /*left*/segment_left,query_compress,querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_FORWARD,
							      chrnum,chroffset,chrhigh,chrlength)) != NULL) {
			debug4e(printf("=> plus acceptor: %f at %d (%d mismatches) %d..%d\n",
				       prob,Substring_siteA_pos(substring),nmismatches,
				       Substring_querystart(substring),Substring_queryend(substring)));
			debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
			(*distant_acceptors)[nmismatches] = Listpool_push((*distant_acceptors)[nmismatches],
									  listpool,(void *) substring);
		      }
		    }
		  }
		}
		
		i--;
	      }
	    }
	  }
	}

	  
	/* SENSE ANTI */
	splice_queryend_p = Substring_qend_trim(&splice_pos_qend,&splicetype_queryend,&ambig_prob_queryend,
						segment_left,/*qstart:was elt->qend*/elt->qstart,querylength,
						/*plusp*/true,genestrand,mismatch_positions_alloc,query_compress,chroffset,
						/*sensedir*/SENSE_ANTI);
	debug4e(printf("antisense: splice_pos_qend %d with splice_queryend_p %d\n",splice_pos_qend,splice_queryend_p));
	if (splice_pos_qend >= 0 && (splice_pos_qend > querylength - 6 || splice_queryend_p == true)) {
	  /* genomic right anchor or middle anchor */
	  debug4e(printf("Searching genomic left (calling Genome_mismatches_right): plus genomic right anchor or middle anchor: %d..%d\n",
			 0,splice_pos_qend));
	  mismatch_positions_splice = mismatch_positions_alloc;
	  total_nmismatches_splice = Genome_mismatches_right(mismatch_positions_alloc,max_mismatches_allowed,
							     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							     /*left*/segment_left,/*pos5*/0,/*pos3*/splice_pos_qend,
							     /*plusp*/true,genestrand);
	  debug4e(
		  printf("%d mismatches at splice (%d allowed) at:",total_nmismatches_splice,max_mismatches_allowed);
		  for (i = 0; i <= total_nmismatches_splice; i++) {
		    printf(" %d",mismatch_positions_splice[i]);
		  }
		  printf("\n");
		  );
	  
	  if (total_nmismatches_splice > max_mismatches_allowed) {
	    splice_pos_qstart = mismatch_positions_splice[max_mismatches_allowed] + 1;
	    nmismatches_frag = max_mismatches_allowed; /* should be total_nmismatches_splice - 1 */
	  } else {
	    splice_pos_qstart = mismatch_positions_splice[total_nmismatches_splice] + 1;
	    nmismatches_frag = total_nmismatches_splice;
	  }
	  
	  if (splice_pos_qend - splice_pos_qstart >= MIN_FRAGLENGTH) {
	    if ((substring = Substring_new_endfrag(nmismatches_frag,
						   /*querystart*/splice_pos_qstart,/*queryend*/splice_pos_qend,
						   /*left*/segment_left,query_compress,
						   querylength,/*plusp*/true,genestrand,chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	      debug4e(printf("=> plus endfrag: at %d, %d..%d\n",Substring_siteN_pos(substring),
			     Substring_querystart(substring),Substring_queryend(substring)));
	      debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	      *distant_endfrags = Listpool_push(*distant_endfrags,listpool,(void *) substring);
	    }

	    if (segment_left + splice_pos_qend >= DONOR_MODEL_RIGHT_MARGIN) {
	      /* (2) Splicing originally on minus strand.  Complement.  */
	      debug4e(printf("Search for antidonor splice sites from %d down to %d\n",splice_pos_qend,splice_pos_qstart));
	      antidonorj_nsites = Genome_antidonor_positions_novel(positions_alloc,segment_left,splice_pos_qstart,splice_pos_qend);
	      antidonorj_positions = positions_alloc;
	      debug4e(
		      printf("Antidonor dinucleotides:");
		      for (i = 0; i < antidonorj_nsites; i++) {
			printf(" %d",antidonorj_positions[i]);
		      }
		      printf("\n");
		      );
	    
	      i = antidonorj_nsites - 1;
	      nmismatches = 0;
	      while (i >= 0 && nmismatches <= max_mismatches_allowed) {
		splice_pos = antidonorj_positions[i];
		while (nmismatches <= total_nmismatches_splice &&
		       mismatch_positions_splice[nmismatches] >= splice_pos) { /* Must be >= */
		  debug4e(printf("  mismatch at %d\n",mismatch_positions_splice[nmismatches]));
		  nmismatches++;
		}
		debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
	      
		if (splice_pos_qend - splice_pos >= MIN_SPLICE_EXON_LENGTH) {
		  assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,segment_left,
									  /*pos5*/splice_pos,/*pos3*/splice_pos_qend,/*plusp*/true,genestrand));
		  if (segment_left + querylength < chrhigh) {
		    prob = Maxent_hr_antidonor_prob(segment_left + splice_pos,chroffset);
		    debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
				   splice_pos,nmismatches,prob,sufficient_splice_prob_distant(querylength - splice_pos,nmismatches,prob)));
		    if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		      debug4e(printf("Novel antidonor for segmenti at %llu, splice_pos %d (%d mismatches)\n",
				     (unsigned long long) segment_left,splice_pos,nmismatches));
		      if ((substring = Substring_new_donor(nmismatches,/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/-1,
							   /*querystart*/splice_pos,/*queryend*/splice_pos_qend,/*sitepos*/splice_pos,prob,
							   /*splice_querystart_p*/true,/*splicetype_querystart*/ANTIDONOR,/*ambig_prob_querystart*/prob,
							   splice_queryend_p,splicetype_queryend,ambig_prob_queryend,
							   /*left*/segment_left,query_compress,querylength,/*plusp*/true,genestrand,/*sensedir*/SENSE_ANTI,
							   chrnum,chroffset,chrhigh,chrlength)) != NULL) {
			debug4e(printf("=> plus antidonor: %f at %d (%d mismatches) %d..%d\n",
				       prob,Substring_siteD_pos(substring),nmismatches,
				       Substring_querystart(substring),Substring_queryend(substring)));
			debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
			(*distant_antidonors)[nmismatches] = Listpool_push((*distant_antidonors)[nmismatches],
									   listpool,(void *) substring);
		      }
		    }
		  }
		}
	      
		i--;
	      }
	    }
	  }
	}

      }
    }
  }

  return;
}


static void
find_spliceends_rna_minus (List_T **distant_donors, List_T **distant_antidonors,
			   List_T **distant_acceptors, List_T **distant_antiacceptors,
			   List_T *distant_startfrags, List_T *distant_endfrags,

			   List_T elt_set,
#ifdef DEBUG4E
			   char *queryptr,
#endif
			   int querylength, int query_lastpos,
			   int *mismatch_positions_alloc, int *positions_alloc,
			   Compress_T query_compress, Listpool_T listpool,
			   int max_mismatches_allowed, int genestrand) {
#ifdef DEBUG4E
  char *gbuffer;
#endif

  Elt_T elt;
  List_T p;
  int k;

  Substring_T substring;
  Univcoord_T segment_left;
  int nmismatches, nmismatches_end_trim, i;
  int splice_pos_qstart, splice_pos_qend, splice_pos;
  bool splice_querystart_p, splice_queryend_p;
  Splicetype_T splicetype_querystart, splicetype_queryend;
  double ambig_prob_querystart, ambig_prob_queryend, prob;

  int total_nmismatches_splice, nmismatches_frag;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *mismatch_positions_splice;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;


  debug4e(printf("Entering find_spliceends_rna_minus with %d elts\n",List_length(elt_set)));

  for (p = elt_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);

    debug4e(printf("Comparing elt->qend %d < query_lastpos %d\n",elt->qend,query_lastpos));
    if (elt->qend < query_lastpos /*&& (elt->qstart < index1part || segment->spliceable_low_p == true)*/) {
      /* Find splices on genomic right */
      for (k = 0; k < elt->n_all_diagonals; k++) {
	segment_left = elt->all_diagonals[k];
	chrnum = Univ_IIT_get_one(chromosome_iit,segment_left,segment_left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

	debug4e(printf("find_spliceends_rna_minus: Checking up to %d mismatches at diagonal %llu (elt qstart..qend %d..%d)\n",
		       max_mismatches_allowed,(unsigned long long) segment_left,elt->qstart,elt->qend));
	debug4e(
		gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		printf("genome 0..: %s\n",gbuffer);
		printf("query  0..: %s\n",queryptr);
		FREE(gbuffer);
		);

	/* Trim at end (away from splice, or qstart) */
	debug4e(printf("Trimming at end: %d..%d\n",0,elt->qstart));

	/* SENSE ANTI */
	splice_queryend_p = Substring_qstart_trim(&splice_pos_qstart,&splicetype_queryend,&ambig_prob_queryend,
						  segment_left,/*qend: was elt->qstart*/elt->qend,
						  /*plusp*/false,genestrand,mismatch_positions_alloc,query_compress,chroffset,
						  /*sensedir*/SENSE_ANTI);
	debug4e(printf("antisense: splice_pos_qstart %d with splice_queryend_p %d\n",splice_pos_qstart,splice_queryend_p));
	if (splice_pos_qstart >= 0 && (splice_pos_qstart < 6 || splice_queryend_p == true)) {
	  /* genomic left anchor or middle anchor */
	  debug4e(printf("Searching genomic right (calling Genome_mismatches_left): minus genomic left anchor or middle anchor: %d..%d\n",
			 elt->qstart,querylength));
	  mismatch_positions_splice = mismatch_positions_alloc;
	  total_nmismatches_splice = Genome_mismatches_left(mismatch_positions_alloc,max_mismatches_allowed,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    /*left*/segment_left,/*pos5*/splice_pos_qstart,/*pos3*/querylength,
							    /*plusp*/false,genestrand);
	  debug4e(
		  printf("%d mismatches at splice (%d allowed) at:",total_nmismatches_splice,max_mismatches_allowed);
		  for (i = 0; i <= total_nmismatches_splice; i++) {
		    printf(" %d",mismatch_positions_splice[i]);
		  }
		  printf("\n");
		  );
	  
	  if (total_nmismatches_splice > max_mismatches_allowed) {
	    splice_pos_qend = mismatch_positions_splice[max_mismatches_allowed];
	    nmismatches_frag = max_mismatches_allowed; /* should be total_nmismatches_splice - 1 */
	  } else {
	    splice_pos_qend = mismatch_positions_splice[total_nmismatches_splice];
	    nmismatches_frag = total_nmismatches_splice;
	  }

	  if (splice_pos_qend - splice_pos_qstart >= MIN_FRAGLENGTH) {
	    if ((substring = Substring_new_endfrag(nmismatches_frag,
						   /*querystart*/querylength - splice_pos_qend,/*queryend*/querylength - splice_pos_qstart,
						   /*left*/segment_left,query_compress,
						   querylength,/*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	      debug4e(printf("=> minus endfrag: at %d, %d..%d\n",Substring_siteN_pos(substring),
			     Substring_querystart(substring),Substring_queryend(substring)));
	      debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	      *distant_endfrags = Listpool_push(*distant_endfrags,listpool,(void *) substring);
	    }
	  
	    if (segment_left + splice_pos_qstart >= DONOR_MODEL_LEFT_MARGIN) {
	      /* (1) Originally on plus strand.  No complement */
	      debug4e(printf("Search for donor splice sites from %d up to %d\n",splice_pos_qstart,splice_pos_qend));
	      donori_nsites = Genome_donor_positions_novel(positions_alloc,segment_left,splice_pos_qstart,splice_pos_qend);
	      donori_positions = positions_alloc;
	      debug4e(
		      printf("Donor dinucleotides:");
		      for (i = 0; i < donori_nsites; i++) {
			printf(" %d",donori_positions[i]);
		      }
		      printf("\n");
		      );
	    
	      i = 0;
	      nmismatches = 0;
	      while (i < donori_nsites && nmismatches <= max_mismatches_allowed) {
		splice_pos = donori_positions[i];
		while (nmismatches <= total_nmismatches_splice &&
		       mismatch_positions_splice[nmismatches] < splice_pos) { /* Changed from <= to < */
		  debug4e(printf("  mismatch at %d\n",mismatch_positions_splice[nmismatches]));
		  nmismatches++;
		}
		debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
	      
		if (splice_pos - splice_pos_qstart >= MIN_SPLICE_EXON_LENGTH) {
		  assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,segment_left,
									  /*pos5*/splice_pos_qstart,/*pos3*/splice_pos,/*plusp*/false,genestrand));
		  if (segment_left + querylength < chrhigh) {
		    prob = Maxent_hr_donor_prob(segment_left + splice_pos,chroffset);
		    debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
				   splice_pos,nmismatches,prob,sufficient_splice_prob_distant(splice_pos,nmismatches,prob)));
		    if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		      debug4e(printf("Novel donor for segment at %llu, splice_pos %d (%d mismatches)\n",
				     (unsigned long long) segment_left,splice_pos,nmismatches));
		      if ((substring = Substring_new_donor(nmismatches,/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/-1,
							   /*querystart*/querylength - splice_pos,/*queryend*/querylength - splice_pos_qstart,
							   /*sitepos*/querylength - splice_pos,prob,
							   /*splice_querystart_p*/true,/*splicetype_querystart*/ANTIDONOR,/*ambig_prob_querystart*/prob,
							   splice_queryend_p,splicetype_queryend,ambig_prob_queryend,
							   /*left*/segment_left,query_compress,querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
							   chrnum,chroffset,chrhigh,chrlength)) != NULL) {
			debug4e(printf("=> minus donor: %f at %d (%d mismatches) %d..%d\n",
				       prob,Substring_siteD_pos(substring),nmismatches,
				       Substring_querystart(substring),Substring_queryend(substring)));
			debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
			(*distant_donors)[nmismatches] = Listpool_push((*distant_donors)[nmismatches],
								       listpool,(void *) substring);
		      }
		    }
		  }
		}

		i++;
	      }
	    }
	  }
	}


	/* SENSE FORWARD */
	splice_queryend_p = Substring_qstart_trim(&splice_pos_qstart,&splicetype_queryend,&ambig_prob_queryend,
						  segment_left,/*qend:was elt->qstart*/elt->qend,
						  /*plusp*/false,genestrand,mismatch_positions_alloc,query_compress,chroffset,
						  /*sensedir*/SENSE_FORWARD);
	debug4e(printf("sense: splice_pos_qstart %d with splice_queryend_p %d\n",splice_pos_qstart,splice_queryend_p));
	if (splice_pos_qstart >= 0 && (splice_pos_qstart < 6 || splice_queryend_p == true)) {
	  /* genomic left anchor or middle anchor */
	  debug4e(printf("Searching genomic right (calling Genome_mismatches_left): minus genomic left anchor or middle anchor: %d..%d\n",
			 elt->qstart,querylength));
	  mismatch_positions_splice = mismatch_positions_alloc;
	  total_nmismatches_splice = Genome_mismatches_left(mismatch_positions_alloc,max_mismatches_allowed,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							    /*left*/segment_left,/*pos5*/splice_pos_qstart,/*pos3*/querylength,
							    /*plusp*/false,genestrand);
	  debug4e(
		  printf("%d mismatches at splice (%d allowed) at:",total_nmismatches_splice,max_mismatches_allowed);
		  for (i = 0; i <= total_nmismatches_splice; i++) {
		    printf(" %d",mismatch_positions_splice[i]);
		  }
		  printf("\n");
		  );
	  
	  if (total_nmismatches_splice > max_mismatches_allowed) {
	    splice_pos_qend = mismatch_positions_splice[max_mismatches_allowed];
	    nmismatches_frag = max_mismatches_allowed; /* should be total_nmismatches_splice - 1 */
	  } else {
	    splice_pos_qend = mismatch_positions_splice[total_nmismatches_splice];
	    nmismatches_frag = total_nmismatches_splice;
	  }

	  if (splice_pos_qend - splice_pos_qstart >= MIN_FRAGLENGTH) {
	    if ((substring = Substring_new_endfrag(nmismatches_frag,
						   /*querystart*/querylength - splice_pos_qend,/*queryend*/querylength - splice_pos_qstart,
						   /*left*/segment_left,query_compress,
						   querylength,/*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	      debug4e(printf("=> minus endfrag: at %d, %d..%d\n",Substring_siteN_pos(substring),
			     Substring_querystart(substring),Substring_queryend(substring)));
	      debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	      *distant_endfrags = Listpool_push(*distant_endfrags,listpool,(void *) substring);
	    }

	    if (segment_left + splice_pos_qstart >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	      /* (2) Splicing originally on minus strand.  Complement */
	      debug4e(printf("Search for antiacceptor splice sites from %d up to %d\n",splice_pos_qstart,splice_pos_qend));
	      antiacceptori_nsites = Genome_antiacceptor_positions_novel(positions_alloc,segment_left,splice_pos_qstart,splice_pos_qend);
	      antiacceptori_positions = positions_alloc;
	      debug4e(
		      printf("Antiacceptor dinucleotides:");
		      for (i = 0; i < antiacceptori_nsites; i++) {
			printf(" %d",antiacceptori_positions[i]);
		      }
		      printf("\n");
		      );
	    
	      i = 0;
	      nmismatches = 0;
	      while (i < antiacceptori_nsites && nmismatches <= max_mismatches_allowed) {
		splice_pos = antiacceptori_positions[i];
		while (nmismatches <= total_nmismatches_splice &&
		       mismatch_positions_splice[nmismatches] < splice_pos) { /* Changed from <= to < */
		  debug4e(printf("  mismatch at %d\n",mismatch_positions_splice[nmismatches]));
		  nmismatches++;
		}
		debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
	      
		if (splice_pos - splice_pos_qstart >= MIN_SPLICE_EXON_LENGTH) {
		  assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,segment_left,
									  /*pos5*/splice_pos_qstart,/*pos3*/splice_pos,/*plusp*/false,genestrand));
		  if (segment_left + querylength < chrhigh) {
		    prob = Maxent_hr_antiacceptor_prob(segment_left + splice_pos,chroffset);
		    debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
				   splice_pos,nmismatches,prob,sufficient_splice_prob_distant(splice_pos,nmismatches,prob)));
		    if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		      debug4e(printf("Novel antiacceptor for segment at %llu, splice_pos %d (%d mismatches)\n",
				     (unsigned long long) segment_left,splice_pos,nmismatches));
		      if ((substring = Substring_new_acceptor(nmismatches,/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/-1,
							      /*querystart*/querylength - splice_pos,/*queryend*/querylength - splice_pos_qstart,
							      /*sitepos*/querylength - splice_pos,prob,
							      /*splice_querystart_p*/true,/*splicetype_querystart*/ACCEPTOR,/*ambig_prob_querystart*/prob,
							      splice_queryend_p,splicetype_queryend,ambig_prob_queryend,
							      /*left*/segment_left,query_compress,querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
							      chrnum,chroffset,chrhigh,chrlength)) != NULL) {
			debug4e(printf("=> minus antiacceptor : %f at %d (%d mismatches) %d..%d\n",
				       prob,Substring_siteA_pos(substring),nmismatches,
				       Substring_querystart(substring),Substring_queryend(substring)));
			debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
			(*distant_antiacceptors)[nmismatches] = Listpool_push((*distant_antiacceptors)[nmismatches],
									      listpool,(void *) substring);
		      }
		    }
		  }
		}
	      
		i++;
	      }
	    }
	  }
	}

      }
    }


    debug4e(printf("Comparing elt->qstart %d > index1part %d\n",elt->qstart,index1part));
    if (elt->qstart > index1part /*&& (elt->qend > query_lastpos || segment->spliceable_high_p == true)*/) {
      /* Find splices on genomic left */
      for (k = 0; k < elt->n_all_diagonals; k++) {
	segment_left = elt->all_diagonals[k];
	chrnum = Univ_IIT_get_one(chromosome_iit,segment_left,segment_left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);

	debug4e(printf("find_spliceends_rna_minus: Checking up to %d mismatches at diagonal %llu (elt qstart..qend %d..%d)\n",
		       max_mismatches_allowed,(unsigned long long) segment_left,elt->qstart,elt->qend));
	debug4e(
		gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		printf("genome 0..: %s\n",gbuffer);
		printf("query  0..: %s\n",queryptr);
		FREE(gbuffer);
		);

	/* Trim at end (away from splice, or qend) */
	debug4e(printf("Trimming at end: %d..%d\n",elt->qend,querylength));

	/* SENSE ANTI */
	splice_querystart_p = Substring_qend_trim(&splice_pos_qend,&splicetype_querystart,&ambig_prob_querystart,
						  segment_left,/*qstart:was elt->qend*/elt->qstart,querylength,
						  /*plusp*/false,genestrand,mismatch_positions_alloc,query_compress,chroffset,
						  /*sensedir*/SENSE_ANTI);
	debug4e(printf("antisense: splice_pos_qend %d with splice_querystart_p %d\n",splice_pos_qend,splice_querystart_p));
	if (splice_pos_qend >= 0 && (splice_pos_qend > querylength - 6 || splice_querystart_p == true)) {
	  /* genomic right anchor or middle anchor */
	  debug4e(printf("Searching genomic left (calling Genome_mismatches_right): minus genomic right anchor or middle anchor: %d..%d\n",
			 0,splice_pos_qend));
	  mismatch_positions_splice = mismatch_positions_alloc;
	  total_nmismatches_splice = Genome_mismatches_right(mismatch_positions_alloc,max_mismatches_allowed,
							     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							     /*left*/segment_left,/*pos5*/0,/*pos3*/splice_pos_qend,
							     /*plusp*/false,genestrand);
	  debug4e(
		  printf("%d mismatches at splice (%d allowed) at:",total_nmismatches_splice,max_mismatches_allowed);
		  for (i = 0; i <= total_nmismatches_splice; i++) {
		    printf(" %d",mismatch_positions_splice[i]);
		  }
		  printf("\n");
		  );
	  
	  if (total_nmismatches_splice > max_mismatches_allowed) {
	    splice_pos_qstart = mismatch_positions_splice[max_mismatches_allowed] + 1;
	    nmismatches_frag = max_mismatches_allowed; /* should be total_nmismatches_splice - 1 */
	  } else {
	    splice_pos_qstart = mismatch_positions_splice[total_nmismatches_splice] + 1;
	    nmismatches_frag = total_nmismatches_splice;
	  }

	  if (splice_pos_qend - splice_pos_qstart >= MIN_FRAGLENGTH) {
	    if ((substring = Substring_new_startfrag(nmismatches_frag,
						     /*querystart*/querylength - splice_pos_qend,/*queryend*/querylength - splice_pos_qstart,
						     /*left*/segment_left,query_compress,
						     querylength,/*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	      debug4e(printf("=> minus startfrag: at %d, %d..%d\n",Substring_siteN_pos(substring),
			     Substring_querystart(substring),Substring_queryend(substring)));
	      debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	      *distant_startfrags = Listpool_push(*distant_startfrags,listpool,(void *) substring);
	    }

	    if (segment_left + splice_pos_qend >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	      /* (1) Splicing originally on plus strand.  No complement. */
	      debug4e(printf("Search for acceptor splice sites from %d down to %d\n",splice_pos_qend,splice_pos_qstart));
	      acceptorj_nsites = Genome_acceptor_positions_novel(positions_alloc,segment_left,splice_pos_qstart,splice_pos_qend);
	      acceptorj_positions = positions_alloc;
	      debug4e(
		      printf("Acceptor dinucleotides:");
		      for (i = 0; i < acceptorj_nsites; i++) {
			printf(" %d",acceptorj_positions[i]);
		      }
		      printf("\n");
		      );
	    
	      i = acceptorj_nsites - 1;
	      nmismatches = 0;
	      while (i >= 0 && nmismatches <= max_mismatches_allowed) {
		splice_pos = acceptorj_positions[i];
		while (nmismatches <= total_nmismatches_splice &&
		       mismatch_positions_splice[nmismatches] >= splice_pos) { /* Must be >= */
		  debug4e(printf("  mismatch at %d\n",mismatch_positions_splice[nmismatches]));
		  nmismatches++;
		}
		debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
	      
		if (splice_pos_qend - splice_pos >= MIN_SPLICE_EXON_LENGTH) {
		  assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,segment_left,
									  /*pos5*/splice_pos,/*pos3*/splice_pos_qend,/*plusp*/false,genestrand));
		  if (segment_left + querylength < chrhigh) {
		    prob = Maxent_hr_acceptor_prob(segment_left + splice_pos,chroffset);
		    debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
				   splice_pos,nmismatches,prob,sufficient_splice_prob_distant(querylength - splice_pos,nmismatches,prob)));
		    if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		      debug4e(printf("Novel acceptor for segment at %llu, splice_pos %d (%d mismatches)\n",
				     (unsigned long long) segment_left,splice_pos,nmismatches));
		      if ((substring = Substring_new_acceptor(nmismatches,/*acceptor_coord*/segment_left + splice_pos,/*acceptor_knowni*/-1,
							      /*querystart*/querylength - splice_pos_qend,/*queryend*/querylength - splice_pos,
							      /*sitepos*/querylength - splice_pos,prob,
							      splice_querystart_p,splicetype_querystart,ambig_prob_querystart,
							      /*splice_queryend_p*/true,/*splicetype_queryend*/ANTIACCEPTOR,/*ambig_prob_queryend*/prob,
							      /*left*/segment_left,query_compress,querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_ANTI,
							      chrnum,chroffset,chrhigh,chrlength)) != NULL) {
			debug4e(printf("=> minus acceptor: %f at %d (%d mismatches) %d..%d\n",
				       prob,Substring_siteA_pos(substring),nmismatches,
				       Substring_querystart(substring),Substring_queryend(substring)));
			debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
			(*distant_acceptors)[nmismatches] = Listpool_push((*distant_acceptors)[nmismatches],
									  listpool,(void *) substring);
		      }
		    }
		  }
		}
	      
		i--;
	      }
	    }
	  }
	}
	  

	/* SENSE FORWARD */
	splice_querystart_p = Substring_qend_trim(&splice_pos_qend,&splicetype_querystart,&ambig_prob_querystart,
						  segment_left,/*qstart:was elt->qend*/elt->qstart,querylength,
						  /*plusp*/false,genestrand,mismatch_positions_alloc,query_compress,chroffset,
						  /*sensedir*/SENSE_FORWARD);
	debug4e(printf("sense: splice_pos_qend %d with splice_querystart_p %d\n",splice_pos_qend,splice_querystart_p));
	if (splice_pos_qend >= 0 && (splice_pos_qend > querylength - 6 || splice_querystart_p == true)) {
	  /* genomic right anchor or middle anchor */
	  debug4e(printf("Searching genomic left (calling Genome_mismatches_right): minus genomic right anchor or middle anchor: %d..%d\n",
			 0,splice_pos_qend));
	  mismatch_positions_splice = mismatch_positions_alloc;
	  total_nmismatches_splice = Genome_mismatches_right(mismatch_positions_alloc,max_mismatches_allowed,
							     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							     /*left*/segment_left,/*pos5*/0,/*pos3*/splice_pos_qend,
							     /*plusp*/false,genestrand);
	  debug4e(
		  printf("%d mismatches at splice (%d allowed) at:",total_nmismatches_splice,max_mismatches_allowed);
		  for (i = 0; i <= total_nmismatches_splice; i++) {
		    printf(" %d",mismatch_positions_splice[i]);
		  }
		  printf("\n");
		  );
	  
	  if (total_nmismatches_splice > max_mismatches_allowed) {
	    splice_pos_qstart = mismatch_positions_splice[max_mismatches_allowed] + 1;
	    nmismatches_frag = max_mismatches_allowed; /* should be total_nmismatches_splice - 1 */
	  } else {
	    splice_pos_qstart = mismatch_positions_splice[total_nmismatches_splice] + 1;
	    nmismatches_frag = total_nmismatches_splice;
	  }

	  if (splice_pos_qend - splice_pos_qstart >= MIN_FRAGLENGTH) {
	    if ((substring = Substring_new_startfrag(nmismatches_frag,
						     /*querystart*/querylength - splice_pos_qend,/*queryend*/querylength - splice_pos_qstart,
						     /*left*/segment_left,query_compress,
						     querylength,/*plusp*/false,genestrand,chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	      debug4e(printf("=> minus startfrag: at %d, %d..%d\n",Substring_siteN_pos(substring),
			     Substring_querystart(substring),Substring_queryend(substring)));
	      debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	      *distant_startfrags = Listpool_push(*distant_startfrags,listpool,(void *) substring);
	    }

	    if (segment_left + splice_pos_qend >= DONOR_MODEL_RIGHT_MARGIN) {
	      /* (2) Splicing originally on minus strand.  Complement.  */
	      debug4e(printf("Search for antidonor splice sites from %d down to %d\n",splice_pos_qend,splice_pos_qstart));
	      antidonorj_nsites = Genome_antidonor_positions_novel(positions_alloc,segment_left,splice_pos_qstart,splice_pos_qend);
	      antidonorj_positions = positions_alloc;
	      debug4e(
		      printf("Antidonor dinucleotides:");
		      for (i = 0; i < antidonorj_nsites; i++) {
			printf(" %d",antidonorj_positions[i]);
		      }
		      printf("\n");
		      );
	    
	      i = antidonorj_nsites - 1;
	      nmismatches = 0;
	      while (i >= 0 && nmismatches <= max_mismatches_allowed) {
		splice_pos = antidonorj_positions[i];
		while (nmismatches <= total_nmismatches_splice &&
		       mismatch_positions_splice[nmismatches] >= splice_pos) { /* Must be >= */
		  debug4e(printf("  mismatch at %d\n",mismatch_positions_splice[nmismatches]));
		  nmismatches++;
		}
		debug4e(printf(" splice pos %d, nmismatches %d\n",splice_pos,nmismatches));
	      
		if (splice_pos_qend - splice_pos >= MIN_SPLICE_EXON_LENGTH) {
		  assert(nmismatches == Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,segment_left,
									  /*pos5*/splice_pos,/*pos3*/splice_pos_qend,/*plusp*/false,genestrand));
		  if (segment_left + querylength < chrhigh) {
		    prob = Maxent_hr_antidonor_prob(segment_left + splice_pos,chroffset);
		    debug4e(printf("splice pos %d, nmismatches %d, prob %f, sufficient %d\n",
				   splice_pos,nmismatches,prob,sufficient_splice_prob_distant(querylength - splice_pos,nmismatches,prob)));
		    if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		      debug4e(printf("Novel antidonor for segmenti at %llu, splice_pos %d (%d mismatches)\n",
				     (unsigned long long) segment_left,splice_pos,nmismatches));
		      if ((substring = Substring_new_donor(nmismatches,/*donor_coord*/segment_left + splice_pos,/*donor_knowni*/-1,
							   /*querystart*/querylength - splice_pos_qend,/*queryend*/querylength - splice_pos,
							   /*sitepos*/querylength - splice_pos,prob,
							   splice_querystart_p,splicetype_querystart,ambig_prob_querystart,
							   /*splice_queryend_p*/true,/*splicetype_queryend*/DONOR,/*ambig_prob_queryend*/prob,
							   /*left*/segment_left,query_compress,querylength,/*plusp*/false,genestrand,/*sensedir*/SENSE_FORWARD,
							   chrnum,chroffset,chrhigh,chrlength)) != NULL) {
			debug4e(printf("=> minus antidonor: %f at %d (%d mismatches) %d..%d\n",
				       prob,Substring_siteD_pos(substring),nmismatches,
				       Substring_querystart(substring),Substring_queryend(substring)));
			debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
			(*distant_antidonors)[nmismatches] = Listpool_push((*distant_antidonors)[nmismatches],
									   listpool,(void *) substring);
		      }
		    }
		  }
		}
	      
		i--;
	      }
	    }
	  }
	}

      }
    }
  }

  return;
}


#if 0
/* Requires known splice sites */
static bool
intragenic_splice_p (Chrpos_T splicedistance, Substring_T donor, Substring_T acceptor) {
  int knowni;

  if ((knowni = Substring_splicesitesD_knowni(donor)) >= 0) {
    if (splicedists[knowni] >= splicedistance) {
      return true;
    }
  }

  if ((knowni = Substring_splicesitesA_knowni(acceptor)) >= 0) {
    if (splicedists[knowni] >= splicedistance) {
      return true;
    }
  }

  return false;
}
#endif


static void
find_splicepairs_rna (int *found_score_overall, int *found_score_within_trims,
		      List_T *hits_plus, List_T *hits_minus,

		      List_T *donors_plus, List_T *antidonors_plus,
		      List_T *acceptors_plus, List_T *antiacceptors_plus,
		      List_T *donors_minus, List_T *antidonors_minus,
		      List_T *acceptors_minus, List_T *antiacceptors_minus,

		      int querylength, int nmismatches_allowed, bool first_read_p,
		      Listpool_T listpool, Hitlistpool_T hitlistpool, int level) {
  List_T localhits_plus = NULL, localhits_minus = NULL, distanthits_plus = NULL, distanthits_minus = NULL,
    p, q, qsave;
  int ndistantsplicepairs = 0;
  Substring_T donor, acceptor;
  int min_endlength_1, min_endlength_2, nmismatches1, nmismatches2, pos;
  Chrpos_T distance;
  Univcoord_T donor_genomicstart, acceptor_genomicstart;
  bool shortdistancep;
  double nonidentity = 1.0 - min_distantsplicing_identity;
  Chrnum_T chrnum;
  Stage3end_T splice;

  debug(printf("Starting find_splicepairs_rna with nonidentity %f\n",nonidentity));
  debug4l(printf("Starting find_splicepairs_rna with nonidentity %f\n",nonidentity));

  if (nonidentity == 0.0) {
    nmismatches_allowed = 0;
  }

  for (nmismatches1 = 0; nmismatches1 <= nmismatches_allowed; nmismatches1++) {
    nmismatches2 = nmismatches_allowed - nmismatches1;

    if (nonidentity == 0.0) {
      min_endlength_1 = min_endlength_2 = min_distantsplicing_end_matches;
    } else {
      min_endlength_1 = rint((double) nmismatches1/nonidentity);
      if (min_endlength_1 < min_distantsplicing_end_matches) {
	min_endlength_1 = min_distantsplicing_end_matches;
      }
      min_endlength_2 = rint((double) nmismatches2/nonidentity);
      if (min_endlength_2 < min_distantsplicing_end_matches) {
	min_endlength_2 = min_distantsplicing_end_matches;
      }
    }

    debug4l(printf("  nmismatches1 = %d, nmismatches2 = %d, min_endlength_1 = %d, min_endlength_2 = %d\n",
		   nmismatches1,nmismatches2,min_endlength_1,min_endlength_2));

    /************************************************************************
     *   Same strands
     ************************************************************************/

    /* 1.  End 1 to End 2.  Same strands. */
    p = donors_plus[nmismatches1];
    q = acceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_rna (%d+%d mismatches): donors+ (%d) to acceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end1-end2: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of donor < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	/* Generate all pairs at this splice_pos */
	qsave = q;
	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(donor)) != Substring_chrnum(acceptor)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((acceptor_genomicstart = Substring_genomicstart(acceptor)) > (donor_genomicstart = Substring_genomicstart(donor))) {
		distance = acceptor_genomicstart - donor_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
#if 0
		} else if (distances_observed_p == true && intragenic_splice_p(distance,donor,acceptor) == true) {
		  shortdistancep = true;
#endif
		} else {
		  shortdistancep = false;
		}
	      } else {
		distance = donor_genomicstart - acceptor_genomicstart;
		shortdistancep = false; /* scramble */
	      }
	      debug4ld(printf("1-2. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  shortdistancep = %d\n",
			      pos,min_endlength_1,querylength-min_endlength_2,
			      (unsigned long long) Substring_genomicstart(donor),
			      (unsigned long long) Substring_genomicstart(acceptor),shortdistancep));

	      if (shortdistancep) {
		if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
						   donor,acceptor,distance,/*shortdistancep*/true,querylength,
						   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
						   /*sensedir*/SENSE_FORWARD,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
		  /* local splice */
		  localhits_plus = Hitlist_push(localhits_plus,hitlistpool,(void *) splice);
		}
	      } else if (ndistantsplicepairs <= MAXCHIMERAPATHS) {
		if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
						   donor,acceptor,distance,/*shortdistancep*/false,querylength,
						   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
						   /*sensedir*/SENSE_FORWARD,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
		  /* distant splice */
		  distanthits_plus = Hitlist_push(distanthits_plus,hitlistpool,(void *) splice);
		  ndistantsplicepairs++;
		}
	      }

	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }

    /* 4. End 3 to End 4.  Same strands. */
    p = donors_minus[nmismatches1];
    q = acceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_rna (%d+%d mismatches): donors- (%d) to acceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end3-end4: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of donor < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(donor)) != Substring_chrnum(acceptor)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((acceptor_genomicstart = Substring_genomicstart(acceptor)) > (donor_genomicstart = Substring_genomicstart(donor))) {
		distance = acceptor_genomicstart - donor_genomicstart;
		shortdistancep = false; /* scramble */
	      } else {
		distance = donor_genomicstart - acceptor_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
#if 0
		} else if (distances_observed_p == true && intragenic_splice_p(distance,donor,acceptor) == true) {
		  shortdistancep = true;
#endif
		} else {
		  shortdistancep = false;
		}
	      }
	      debug4ld(printf("3-4. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  shortdistancep = %d.\n",
			      pos,min_endlength_1,querylength-min_endlength_2,
			      (unsigned long long) Substring_genomicstart(donor),
			      (unsigned long long) Substring_genomicstart(acceptor),shortdistancep));
	      if (shortdistancep) {
		if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
						   donor,acceptor,distance,/*shortdistancep*/true,querylength,
						   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
						   /*sensedir*/SENSE_FORWARD,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
		  localhits_minus = Hitlist_push(localhits_minus,hitlistpool,(void *) splice);
		}
	      } else if (ndistantsplicepairs <= MAXCHIMERAPATHS) {
		if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
						   donor,acceptor,distance,/*shortdistancep*/false,querylength,
						   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
						   /*sensedir*/SENSE_FORWARD,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
		  distanthits_minus = Hitlist_push(distanthits_minus,hitlistpool,(void *) splice);
		  ndistantsplicepairs++;
		}
	      }
	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }

    /* 5. End 5 to End 6.  Same strands. */
    p = antidonors_plus[nmismatches1];
    q = antiacceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_rna (%d+%d mismatches): antidonors+ (%d) to antiacceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end5-end6: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_2) {
	debug4ld(printf("chimera_pos of donor < min_endlength_2\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_1\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really an continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(donor)) != Substring_chrnum(acceptor)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((acceptor_genomicstart = Substring_genomicstart(acceptor)) > (donor_genomicstart = Substring_genomicstart(donor))) {
		distance = acceptor_genomicstart - donor_genomicstart;
		shortdistancep = false; /* scramble */
	      } else {
		distance = donor_genomicstart - acceptor_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
#if 0
		} else if (distances_observed_p == true && intragenic_splice_p(distance,donor,acceptor) == true) {
		  shortdistancep = true;
#endif
		} else {
		  shortdistancep = false;
		}
	      }

	      debug4ld(printf("5-6. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  shortdistancep = %d\n",
			      pos,min_endlength_2,querylength-min_endlength_1,
			      (unsigned long long) Substring_genomicstart(donor),
			      (unsigned long long) Substring_genomicstart(acceptor),shortdistancep));
	      if (shortdistancep) {
		if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
						   donor,acceptor,distance,/*shortdistancep*/true,querylength,
						   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
						   /*sensedir*/SENSE_ANTI,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
		  localhits_plus = Hitlist_push(localhits_plus,hitlistpool,(void *) splice);
		}
	      } else if (ndistantsplicepairs <= MAXCHIMERAPATHS) {
		if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
						   donor,acceptor,distance,/*shortdistancep*/false,querylength,
						   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
						   /*sensedir*/SENSE_ANTI,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
		  distanthits_plus = Hitlist_push(distanthits_plus,hitlistpool,(void *) splice);
		  ndistantsplicepairs++;
		}
	      }
	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }

    /* 8. End 7 to End 8.  Same strands. */
    p = antidonors_minus[nmismatches1];
    q = antiacceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_rna (%d+%d mismatches): antidonors- (%d) to antiacceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end7-end8: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_2) {
	debug4ld(printf("chimera_pos of donor < min_endlength_2\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_1\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;

	while (p != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL /* && *nsplicepairs <= MAXCHIMERAPATHS */ && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
	      if ((chrnum = Substring_chrnum(donor)) != Substring_chrnum(acceptor)) {
		distance = 0U;
		shortdistancep = false;
	      } else if ((acceptor_genomicstart = Substring_genomicstart(acceptor)) > (donor_genomicstart = Substring_genomicstart(donor))) {
		distance = acceptor_genomicstart - donor_genomicstart;
		if (distance <= shortsplicedist) {
		  shortdistancep = true;
#if 0
		} else if (distances_observed_p == true && intragenic_splice_p(distance,donor,acceptor) == true) {
		  shortdistancep = true;
#endif
		} else {
		  shortdistancep = false;
		}
	      } else {
		distance = donor_genomicstart - acceptor_genomicstart;
		shortdistancep = false; /* scramble */
	      }
	      debug4ld(printf("7-8. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  shortdistancep = %d.\n",
			      pos,min_endlength_2,querylength-min_endlength_1,
			      (unsigned long long) Substring_genomicstart(donor),
			      (unsigned long long) Substring_genomicstart(acceptor),shortdistancep));
	      if (shortdistancep) {
		if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
						   donor,acceptor,distance,/*shortdistancep*/true,querylength,
						   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
						   /*sensedir*/SENSE_ANTI,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
		  localhits_minus = Hitlist_push(localhits_minus,hitlistpool,(void *) splice);
		}
	      } else if (ndistantsplicepairs <= MAXCHIMERAPATHS) {
		if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
						   donor,acceptor,distance,/*shortdistancep*/false,querylength,
						   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
						   /*sensedir*/SENSE_ANTI,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
		  distanthits_minus = Hitlist_push(distanthits_minus,hitlistpool,(void *) splice);
		  ndistantsplicepairs++;
		}
	      }
	    }
	    q = q->rest;

	  }
	  p = p->rest;
	}
      }
    }
  }


  if (localhits_plus != NULL || localhits_minus != NULL) {
    /* A local splice takes precedence over distant splices */
    Stage3end_gc(distanthits_plus);
    Stage3end_gc(distanthits_minus);

    *hits_plus = List_append(*hits_plus,localhits_plus);
    *hits_minus = List_append(*hits_minus,localhits_minus);
    return;

  } else {
    *hits_plus = List_append(*hits_plus,distanthits_plus);
    *hits_minus = List_append(*hits_minus,distanthits_minus);
  }


  for (nmismatches1 = 0; nmismatches1 <= nmismatches_allowed; nmismatches1++) {
    nmismatches2 = nmismatches_allowed - nmismatches1;

    if (nonidentity == 0.0) {
      min_endlength_1 = min_endlength_2 = min_distantsplicing_end_matches;
    } else {
      min_endlength_1 = rint((double) nmismatches1/nonidentity);
      if (min_endlength_1 < min_distantsplicing_end_matches) {
	min_endlength_1 = min_distantsplicing_end_matches;
      }
      min_endlength_2 = rint((double) nmismatches2/nonidentity);
      if (min_endlength_2 < min_distantsplicing_end_matches) {
	min_endlength_2 = min_distantsplicing_end_matches;
      }
    }

    debug4l(printf("  nmismatches1 = %d, nmismatches2 = %d, min_endlength_1 = %d, min_endlength_2 = %d\n",
		   nmismatches1,nmismatches2,min_endlength_1,min_endlength_2));



    /************************************************************************
     *   Different strands
     ************************************************************************/

    /* 2. End 1 to End 4.  Different strands. */
    p = donors_plus[nmismatches1];
    q = acceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_rna (%d+%d mismatches): donors+ (%d) to acceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end1-end4: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of donor < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) - pos) > (Substring_genomicstart(donor) + pos)) {
	      distance = (Substring_genomicstart(acceptor) - pos) - (Substring_genomicstart(donor) + pos);
	    } else {
	      distance = (Substring_genomicstart(donor) + pos) - (Substring_genomicstart(acceptor) - pos);
	    }
	    debug4ld(printf("1-4. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  Different strands, so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    (unsigned long long) Substring_genomicstart(donor),
			    (unsigned long long) Substring_genomicstart(acceptor)));
	    if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
					       donor,acceptor,distance,/*shortdistancep*/false,querylength,
					       /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
					       /*sensedir*/SENSE_FORWARD,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
	      if (Stage3end_plusp(splice) == true) {
		/* Determined by substring_for_concordance */
		*hits_plus = Hitlist_push(*hits_plus,hitlistpool,(void *) splice);
	      } else {
		*hits_minus = Hitlist_push(*hits_minus,hitlistpool,(void *) splice);
	      }
	      ndistantsplicepairs++;
	    }
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }

    /* 3. End 3 to End 2.  Different strands. */
    p = donors_minus[nmismatches1];
    q = acceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_rna (%d+%d mismatches): donors- (%d) to acceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end3-end2: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_1) {
	debug4ld(printf("chimera_pos of donor < min_endlength_1\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_2\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if (Substring_genomicstart(acceptor) > Substring_genomicstart(donor)) {
	      distance = (Substring_genomicstart(acceptor) + pos) - (Substring_genomicstart(donor) - pos);
	    } else {
	      distance = (Substring_genomicstart(donor) - pos) - (Substring_genomicstart(acceptor) + pos);
	    }
	    debug4ld(printf("3-2. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  Different strands so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    (unsigned long long) Substring_genomicstart(donor),
			    (unsigned long long) Substring_genomicstart(acceptor)));
	    if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
					       donor,acceptor,distance,/*shortdistancep*/false,querylength,
					       /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
					       /*sensedir*/SENSE_FORWARD,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
	      if (Stage3end_plusp(splice) == true) {
		/* Determined by substring_for_concordance */
		*hits_plus = Hitlist_push(*hits_plus,hitlistpool,(void *) splice);
	      } else {
		*hits_minus = Hitlist_push(*hits_minus,hitlistpool,(void *) splice);
	      }
	      ndistantsplicepairs++;
	    }
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }


    /* 6. End 5 to End 8.  Different strands. */
    p = antidonors_plus[nmismatches1];
    q = antiacceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_rna (%d+%d mismatches): antidonors+ (%d) to antiacceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end5-end8: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_2) {
	debug4ld(printf("chimera_pos of donor < min_endlength_2\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_1\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) - pos) > (Substring_genomicstart(donor) + pos)) {
	      distance = (Substring_genomicstart(acceptor) - pos) - (Substring_genomicstart(donor) + pos);
	    } else {
	      distance = (Substring_genomicstart(donor) + pos) - (Substring_genomicstart(acceptor) - pos);
	    }
	    debug4ld(printf("5-8. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  Different strands so not shortdistance.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    (unsigned long long) Substring_genomicstart(donor),
			    (unsigned long long) Substring_genomicstart(acceptor)));
	    if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
					       donor,acceptor,distance,/*shortdistancep*/false,querylength,
					       /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
					       /*sensedir*/SENSE_ANTI,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
	      if (Stage3end_plusp(splice) == true) {
		/* Determined by substring_for_concordance */
		*hits_plus = Hitlist_push(*hits_plus,hitlistpool,(void *) splice);
	      } else {
		*hits_minus = Hitlist_push(*hits_minus,hitlistpool,(void *) splice);
	      }
	      ndistantsplicepairs++;
	    }
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }

    /* 7. End 7 to End 6.  Different strands. */
    p = antidonors_minus[nmismatches1];
    q = antiacceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_rna (%d+%d mismatches): antidonors- (%d) to antiacceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    while (p != NULL && q != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("end7-end6: donor at %llu %d..%d and acceptor at %llu %d..%d\n",
		      (unsigned long long) Substring_genomicstart(donor),
		      Substring_querystart(donor),Substring_queryend(donor),
		      (unsigned long long) Substring_genomicstart(acceptor),
		      Substring_querystart(acceptor),Substring_queryend(acceptor)));

      if ((pos = Substring_siteD_pos(donor)) < min_endlength_2) {
	debug4ld(printf("chimera_pos of donor < min_endlength_2\n"));
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	debug4ld(printf("chimera_pos of donor > querylength - min_endlength_1\n"));
	p = p->rest;
      } else if (pos < Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d < chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	p = p->rest;
      } else if (pos > Substring_siteA_pos(acceptor)) {
	debug4ld(printf("chimera_pos of donor %d > chimera_pos of acceptor %d\n",pos,Substring_siteA_pos(acceptor)));
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteD_pos(((Substring_T) p->first)) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && ndistantsplicepairs <= MAXCHIMERAPATHS && Substring_siteA_pos(((Substring_T) q->first)) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %llu, pos %d\n",(unsigned long long) Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) + pos) > (Substring_genomicstart(donor) - pos)) {
	      distance = (Substring_genomicstart(acceptor) + pos) - (Substring_genomicstart(donor) - pos);
	    } else {
	      distance = (Substring_genomicstart(donor) - pos) - (Substring_genomicstart(acceptor) + pos);
	    }
	    debug4ld(printf("7-6. Pushing a candidate at splice_pos %d (%d..%d), donor %llu to acceptor %llu.  Different strands so not shortdistance.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    (unsigned long long) Substring_genomicstart(donor),
			    (unsigned long long) Substring_genomicstart(acceptor)));
	    if ((splice = Stage3end_new_splice(&(*found_score_overall),&(*found_score_within_trims),
					       donor,acceptor,distance,/*shortdistancep*/false,querylength,
					       /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
					       /*sensedir*/SENSE_ANTI,listpool,/*method*/DISTANT_RNA,level)) != NULL) {
	      if (Stage3end_plusp(splice) == true) {
		/* Determined by substring_for_concordance */
		*hits_plus = Hitlist_push(*hits_plus,hitlistpool,(void *) splice);
	      } else {
		*hits_minus = Hitlist_push(*hits_minus,hitlistpool,(void *) splice);
	      }
	      ndistantsplicepairs++;
	    }
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  }

#if 0
  debug4l(printf("ndistantsplicepairs %d, maxchimerapaths %d\n",ndistantsplicepairs,MAXCHIMERAPATHS));
  if (ndistantsplicepairs > MAXCHIMERAPATHS) {
    /* Can afford to ignore these if MAXCHIMERAPATHS is set high enough */
    stage3list_gc(&distantsplicing);
    return distantsplicing_orig;
  } else {
    return List_append(distantsplicing_orig,distantsplicing);
  }
#endif

  return;
}


/* done_level should probably be renamed final_level.  opt_level
   should probably be renamed found_level or opt_level. */
void
Distant_rna_solve (int *found_score_overall, int *found_score_within_trims,
		   List_T *hits_plus, List_T *hits_minus,

		   List_T *startfrags_plus, List_T *endfrags_plus,
		   List_T *startfrags_minus, List_T *endfrags_minus,
			   
		   List_T queryfwd_plus_set, List_T queryfwd_minus_set,
		   List_T queryrev_plus_set, List_T queryrev_minus_set,

		   int *mismatch_positions_alloc, int *positions_alloc,
		   Compress_T query_compress_fwd, Compress_T query_compress_rev,
		   char *queryuc_ptr, char *queryrc,
		   int querylength, int max_splice_mismatches,
		   int genestrand, bool first_read_p, Listpool_T listpool,
		   Hitlistpool_T hitlistpool, int level) {
  int query_lastpos = querylength - index1part;
  int done_level, nmismatches;
  int i;
  int nsplicepairs = 0;
  List_T *donors_plus, *antidonors_plus, *acceptors_plus, *antiacceptors_plus,
    *donors_minus, *antidonors_minus, *acceptors_minus, *antiacceptors_minus;
  /* bool ambiguousp; */
#ifdef DEBUG13
  int missing_hit, missing_gmap;
#endif


  /* 9 (RNA).  Find distant splicing for RNA iteratively using both known and novel splice sites */
  debug(printf("*** Stage 9 (RNA).  Distant_rna_solve, allowing %d mismatches ***\n",max_splice_mismatches));

  donors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
  antidonors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
  acceptors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
  antiacceptors_plus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
  donors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
  antidonors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
  acceptors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));
  antiacceptors_minus = (List_T *) CALLOCA(max_splice_mismatches+1,sizeof(List_T));

  *startfrags_plus = *endfrags_plus = (List_T) NULL;
  *startfrags_minus = *endfrags_minus = (List_T) NULL;

  debug(printf("Starting find_spliceends_rna (plus)\n"));
  find_spliceends_rna_plus(&donors_plus,&antidonors_plus,&acceptors_plus,&antiacceptors_plus,
			   &(*startfrags_plus),&(*endfrags_plus),queryfwd_plus_set,

#ifdef DEBUG4E
			   /*queryptr*/queryuc_ptr,
#endif
			   querylength,query_lastpos,
			   mismatch_positions_alloc,positions_alloc,/*query_compress*/query_compress_fwd,
			   listpool,max_splice_mismatches,genestrand);

  find_spliceends_rna_plus(&donors_plus,&antidonors_plus,&acceptors_plus,&antiacceptors_plus,
			   &(*startfrags_plus),&(*endfrags_plus),queryrev_plus_set,
#ifdef DEBUG4E
			   /*queryptr*/queryuc_ptr,
#endif
			   querylength,query_lastpos,
			   mismatch_positions_alloc,positions_alloc,/*query_compress*/query_compress_fwd,
			   listpool,max_splice_mismatches,genestrand);
  debug(printf("Finished find_spliceends_rna (plus)\n"));


  debug(printf("Starting find_spliceends_rna (minus)\n"));
  find_spliceends_rna_minus(&antidonors_minus,&donors_minus,&antiacceptors_minus,&acceptors_minus,
			    &(*startfrags_minus),&(*endfrags_minus),
			    queryfwd_minus_set,
#ifdef DEBUG4E
			    /*queryptr*/queryrc,
#endif
			    querylength,query_lastpos,
			    mismatch_positions_alloc,positions_alloc,/*query_compress*/query_compress_rev,
			    listpool,max_splice_mismatches,genestrand);

  find_spliceends_rna_minus(&antidonors_minus,&donors_minus,&antiacceptors_minus,&acceptors_minus,
			    &(*startfrags_minus),&(*endfrags_minus),
			    queryrev_minus_set,
#ifdef DEBUG4E
			    /*queryptr*/queryrc,
#endif
			    querylength,query_lastpos,
			    mismatch_positions_alloc,positions_alloc,/*query_compress*/query_compress_rev,
			    listpool,max_splice_mismatches,genestrand);
  debug(printf("Finished find_spliceends_rna (minus)\n"));


#if 0
  opt_level = ((*found_score_overall) < opt_level) ? (*found_score_overall) : opt_level;
#endif
  done_level = (*found_score_overall) + subopt_levels;

  nmismatches = 0;
  /* ambiguousp = false; */
  debug(printf("Comparing nmismatches level %d against done_level %d - distantsplicing_penalty %d\n",
	       nmismatches,done_level,distantsplicing_penalty));
  while (nmismatches <= done_level - distantsplicing_penalty && nmismatches <= max_splice_mismatches &&
	 nsplicepairs < MAXCHIMERAPATHS /*&& ambiguousp == false*/) {
    debug(printf("*** Stage 9 (RNA).  Distant splicing, allowing %d mismatches ***\n",nmismatches));

    debug4e(printf("Sorting splice ends\n"));
    donors_plus[nmismatches] = Substring_sort_siteD_halves(donors_plus[nmismatches],listpool,/*ascendingp*/true);
    acceptors_plus[nmismatches] = Substring_sort_siteA_halves(acceptors_plus[nmismatches],listpool,/*ascendingp*/true);
    
    antidonors_plus[nmismatches] = Substring_sort_siteD_halves(antidonors_plus[nmismatches],listpool,/*ascendingp*/false);
    antiacceptors_plus[nmismatches] = Substring_sort_siteA_halves(antiacceptors_plus[nmismatches],listpool,/*ascendingp*/false);
    
    donors_minus[nmismatches] = Substring_sort_siteD_halves(donors_minus[nmismatches],listpool,/*ascendingp*/false);
    acceptors_minus[nmismatches] = Substring_sort_siteA_halves(acceptors_minus[nmismatches],listpool,/*ascendingp*/false);
    
    antidonors_minus[nmismatches] = Substring_sort_siteD_halves(antidonors_minus[nmismatches],listpool,/*ascendingp*/true);
    antiacceptors_minus[nmismatches] = Substring_sort_siteA_halves(antiacceptors_minus[nmismatches],listpool,/*ascendingp*/true);
    
    debug4e(printf("Splice ends at %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		   nmismatches,
		   List_length(donors_plus[nmismatches]),List_length(acceptors_plus[nmismatches]),
		   List_length(antidonors_plus[nmismatches]),List_length(antiacceptors_plus[nmismatches]),
		   List_length(donors_minus[nmismatches]),List_length(acceptors_minus[nmismatches]),
		   List_length(antidonors_minus[nmismatches]),List_length(antiacceptors_minus[nmismatches])));

    find_splicepairs_rna(&(*found_score_overall),&(*found_score_within_trims),&(*hits_plus),&(*hits_minus),
			 donors_plus,antidonors_plus,acceptors_plus,antiacceptors_plus,
			 donors_minus,antidonors_minus,acceptors_minus,antiacceptors_minus,
			 querylength,nmismatches,first_read_p,listpool,hitlistpool,level);

    nmismatches++;
  }

  for (i = 0; i <= max_splice_mismatches; i++) {
    Substring_list_gc(&(donors_plus[i]));
    Substring_list_gc(&(antidonors_plus[i]));
    Substring_list_gc(&(acceptors_plus[i]));
    Substring_list_gc(&(antiacceptors_plus[i]));
    Substring_list_gc(&(donors_minus[i]));
    Substring_list_gc(&(antidonors_minus[i]));
    Substring_list_gc(&(acceptors_minus[i]));
    Substring_list_gc(&(antiacceptors_minus[i]));
  }
  FREEA(donors_plus);
  FREEA(antidonors_plus);
  FREEA(acceptors_plus);
  FREEA(antiacceptors_plus);
  FREEA(donors_minus);
  FREEA(antidonors_minus);
  FREEA(acceptors_minus);
  FREEA(antiacceptors_minus);

  return;
}


void
Distant_rna_setup (Univ_IIT_T chromosome_iit_in, int circular_typeint_in,
		   Genome_T genomebits_in, Genome_T genomebits_alt_in,
		   int index1part_in, int index1interval_in,
		   int subopt_levels_in, Chrpos_T shortsplicedist_in,
		   int min_distantsplicing_end_matches_in, int min_distantsplicing_identity_in,
		   int localsplicing_penalty_in, int distantsplicing_penalty_in) {

  chromosome_iit = chromosome_iit_in;
  circular_typeint = circular_typeint_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

  index1part = index1part_in;
  index1interval = index1interval_in;

  subopt_levels = subopt_levels_in;
  shortsplicedist = shortsplicedist_in;
  min_distantsplicing_end_matches = min_distantsplicing_end_matches_in;
  min_distantsplicing_identity = min_distantsplicing_identity_in;
  
  localsplicing_penalty = localsplicing_penalty_in;
  distantsplicing_penalty = distantsplicing_penalty_in;

  return;
}


