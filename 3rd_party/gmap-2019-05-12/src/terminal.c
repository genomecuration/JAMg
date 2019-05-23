static char rcsid[] = "$Id: terminal.c 218675 2019-03-16 01:25:48Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "terminal.h"

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "assert.h"
#include "types.h"

#include "genomicpos.h"
#include "substring.h"

#include "genome128_hr.h"
#include "extension-search.h"	/* For handling Elt_T objects */


#define MAXTERMINALS 1000

static Univ_IIT_T chromosome_iit;
static int circular_typeint;

static Genome_T genomebits;
static Genome_T genomebits_alt;

static int index1part = 15;
static int index1interval = 3;

static int subopt_levels;


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


/* Modified from find_spliceends_rna */
static List_T
find_terminals (List_T terminals, List_T elt_set,
#ifdef DEBUG4E
		char *queryptr,
#endif
		int querylength, int query_lastpos,
		Compress_T query_compress, Hitlistpool_T hitlistpool,
		bool plusp, int genestrand) {
#ifdef DEBUG4E
  char *gbuffer;
#endif

  Elt_T elt;
  List_T p;
  int k;

  Substring_T hit;
  Univcoord_T segment_left;
  int first_querypos, last_querypos;

  /* int nmismatches_left, nmismatches_right; */
#ifdef DEBUG4E
  int i;
#endif

  int *mismatch_positions;
  int *positions_alloc;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    mismatch_positions = (int *) ALLOCA((querylength+1)*sizeof(int));
    positions_alloc = (int *) ALLOCA((querylength+1)*sizeof(int));
  } else {
    mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
    positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  }
#else
  mismatch_positions = (int *) MALLOC((querylength+1)*sizeof(int));
  positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
#endif

  debug4e(printf("Entering find_spliceends_rna (plusp %d) with %d elts\n",plusp,List_length(elt_set)));

  for (p = elt_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
#if 0
    /* Not sure if this is still necessary */
    if ((first_querypos = segment->querystart - (index1interval - 1)) < 0) {
      first_querypos = 0;
    }
    if ((last_querypos = segment->queryend + (index1interval - 1)) > querylength) {
      last_querypos = querylength;
    }
#else
    first_querypos = elt->qstart;
    last_querypos = elt->qend;
#endif

    if (last_querypos < query_lastpos /*&& (first_querypos < index1part || segment->spliceable_low_p == true)*/) {
      /* Find splices on genomic right */
      for (k = 0; k < elt->n_all_diagonals; k++) {
	segment_left = elt->all_diagonals[k];

	debug4e(printf("find_terminals: Checking up to %d mismatches at diagonal %llu (querypos %d..%d), plusp %d\n",
		       max_mismatches_allowed,(unsigned long long) segment_left,
		       first_querypos,last_querypos,plusp));
	debug4e(
		gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		printf("genome 0..: %s\n",gbuffer);
		printf("query  0..: %s\n",queryptr);
		FREE(gbuffer);
		);

#if 0
	if (plusp) {
	  /* genomic left anchor or middle anchor */
	  debug4e(printf("Searching genomic right: plus genomic left anchor or middle anchor: %d..%d\n",
			 first_querypos,querylength));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						    /*left*/segment_left,/*pos5*/first_querypos,/*pos3*/querylength,
						    plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );
	  
	  terminal_start = first_querypos + 1;
	  if (nmismatches_left <= max_mismatches_allowed) {
	    terminal_end = querylength - 1 - 1;
	  } else if ((terminal_end = mismatch_positions[nmismatches_left-1]) > querylength - 1 - 1) {
	    terminal_end = querylength - 1 - 1;
	  }

	} else {
	  /* Minus */
	  /* genomic left anchor or middle anchor */
	  debug4e(printf("Searching genomic right: minus genomic left anchor or middle anchor: %d..%d\n",
			 querylength - last_querypos,querylength));
	  nmismatches_left = Genome_mismatches_left(mismatch_positions,max_mismatches_allowed,
						    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						    /*left*/segment_left,/*pos5*/querylength - last_querypos,
						    /*pos3*/querylength,plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on left (%d allowed) at:",
			 nmismatches_left,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_left; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );

	  terminal_start = querylength - last_querypos + 1;
	  if (nmismatches_left <= max_mismatches_allowed) {
	    terminal_end = querylength - 1 - 1;
	  } else if ((terminal_end = mismatch_positions[nmismatches_left-1]) > querylength - 1 - 1) {
	    terminal_end = querylength - 1 - 1;
	  }
	}

	/* if (terminal_start <= terminal_end) */
#endif

	chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,segment_left,
				     querylength,circular_typeint);
	  
	if ((hit = Substring_new(/*nmismatches*/-1,segment_left,/*querystart*/0,/*queryend*/querylength,querylength,
				 plusp,genestrand,query_compress,chrnum,chroffset,chrhigh,chrlength,
				 /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				 /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				 /*orig_sensedir*/SENSE_NULL)) != NULL) {
	  debug4e(printf("=> %s terminal: (%d mismatches) %d..%d\n",
			 plusp == true ? "plus" : "minus",Substring_nmismatches_bothdiff(hit),
			 Substring_querystart(hit),Substring_queryend(hit)));
	  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	  terminals = Hitlist_push(terminals,hitlistpool,(void *) hit);
	}
      }
    }

    if (first_querypos > index1part /*&& (last_querypos > query_lastpos || segment->spliceable_high_p == true)*/) {
      /* Find splices on genomic left */
      for (k = 0; k < elt->n_all_diagonals; k++) {
	segment_left = elt->all_diagonals[k];

	debug4e(printf("find_terminals: Checking up to %d mismatches at diagonal %llu (querypos %d..%d), plusp %d\n",
		       max_mismatches_allowed,(unsigned long long) segment_left,
		       first_querypos,last_querypos,plusp));
	debug4e(
		gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
		Genome_fill_buffer_blocks(segment_left,querylength,gbuffer);
		printf("genome 0..: %s\n",gbuffer);
		printf("query  0..: %s\n",queryptr);
		FREE(gbuffer);
		);

#if 0
	if (plusp) {
	  /* genomic right anchor or middle anchor */
	  debug4e(printf("Searching genomic left: plus genomic right anchor or middle anchor: %d..%d\n",
			 0,last_querypos));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						      /*left*/segment_left,/*pos5*/0,/*pos3*/last_querypos,
						      plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );
	
	  terminal_end = last_querypos - 1 - 1;
	  if (nmismatches_right <= max_mismatches_allowed) {
	    terminal_start = 1;
	  } else if ((terminal_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    terminal_start = 1;
	  }

	} else {
	  /* Minus */
	  /* genomic right anchor or middle anchor*/
	  debug4e(printf("Searching genomic left: minus genomic right anchor or middle anchor: %d..%d\n",
			 0,querylength - first_querypos));
	  nmismatches_right = Genome_mismatches_right(mismatch_positions,max_mismatches_allowed,
						      /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						      /*left*/segment_left,/*pos5*/0,/*pos3*/querylength - first_querypos,
						      plusp,genestrand);
	  debug4e(
		  printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		  for (i = 0; i <= nmismatches_right; i++) {
		    printf(" %d",mismatch_positions[i]);
		  }
		  printf("\n");
		  );
	  
	  terminal_end = querylength - first_querypos - 1 - 1;
	  if (nmismatches_right <= max_mismatches_allowed) {
	    terminal_start = 1;
	  } else if ((terminal_start = mismatch_positions[nmismatches_right-1]) < 1) {
	    terminal_start = 1;
	  }
	}

	/* if (terminal_start <= terminal_end) */
#endif

	chrnum = Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,segment_left,
				     querylength,circular_typeint);
	if ((hit = Substring_new(/*nmismatches*/-1,segment_left,/*querystart*/0,/*queryend*/querylength,querylength,
				 plusp,genestrand,query_compress,chrnum,chroffset,chrhigh,chrlength,
				 /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				 /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				 /*orig_sensedir*/SENSE_NULL)) != NULL) {
	  debug4e(printf("=> %s terminal: (%d mismatches) %d..%d\n",
			 plusp == true ? "plus" : "minus",Substring_nmismatches_bothdiff(hit),
			 Substring_querystart(hit),Substring_queryend(hit)));
	  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
	  terminals = Hitlist_push(terminals,hitlistpool,(void *) hit);
	}
      }
    }
  }


#ifdef HAVE_ALLOCA
  if (querylength <= MAX_STACK_READLENGTH) {
    FREEA(mismatch_positions);
    FREEA(positions_alloc);
  } else {
    FREE(mismatch_positions);
    FREE(positions_alloc);
  }
#else
  FREE(mismatch_positions);
  FREE(positions_alloc);
#endif

  return terminals;
}


/* done_level should probably be renamed final_level.  opt_level
   should probably be renamed found_level or opt_level. */
List_T
Terminal_solve_plus (int *found_score_overall, int *found_score_within_trims,
		     List_T queryfwd_plus_set, List_T queryrev_plus_set,

		     Compress_T query_compress_fwd, int querylength, 
		     int genestrand, Listpool_T listpool,
		     Hitlistpool_T hitlistpool, int level) {
  List_T terminals = NULL, terminals_plus, p;
  Stage3end_T hit;
  int query_lastpos = querylength - index1part;
  int nterminals;
  Substring_T substring;
#ifdef DEBUG13
  int missing_hit, missing_gmap;
#endif


  /* 9 (Term).  Find terminals */
  debug(printf("*** Stage 9 (Term).  Terminal_solve, allowing %d mismatches ***\n",max_terminal_mismatches));

  terminals_plus = (List_T) NULL;

  debug(printf("Starting find_terminals (plus)\n"));
  terminals_plus = find_terminals(terminals_plus,queryfwd_plus_set,

#ifdef DEBUG4E
				  /*queryptr*/queryuc_ptr,
#endif
				  querylength,query_lastpos,
				  /*query_compress*/query_compress_fwd,
				  hitlistpool,/*plusp*/true,genestrand);

  terminals_plus = find_terminals(terminals_plus,queryrev_plus_set,
#ifdef DEBUG4E
				  /*queryptr*/queryuc_ptr,
#endif
				  querylength,query_lastpos,
				  /*query_compress*/query_compress_fwd,
				  hitlistpool,/*plusp*/true,genestrand);
  debug(printf("Finished find_terminals (plus)\n"));

  debug4e(printf("Terminals: plus %d\n",List_length(terminals_plus)));

#if 0
  opt_level = ((*found_score_within_trims) < opt_level) ? (*found_score_within_trims) : opt_level;
#endif
  /* done_level = (*found_score) + subopt_levels; */

#if 0
  debug4e(printf("Sorting terminals\n"));
  /* TODO: Prioritize inner terminals */
  terminals_left = Substring_sort_nmatches(terminals_left);
  terminals_right = Substring_sort_nmatches(terminals_right);
#endif

  debug(printf("*** Stage 9 (Term)\n"));

  nterminals = 0;
  for (p = terminals_plus; p != NULL && nterminals < MAXTERMINALS; p = List_next(p)) {
    substring = (Substring_T) p->first;
    if ((hit = Stage3end_new_terminal(&(*found_score_overall),&(*found_score_within_trims),
				      substring,querylength,/*gplusp*/true,genestrand,/*sensedir*/SENSE_NULL,
				      listpool,/*method*/TERMINAL,level)) != NULL) {
      terminals = Hitlist_push(terminals,hitlistpool,(void *) hit);
      nterminals++;
    }
  }


  /* Excess distant splicing should be freed already in find_splicepairs_rna */
  terminals = Stage3end_remove_overlaps(terminals,hitlistpool,/*finalp*/false);
#if 0
  /* Do not filter terminals by optimal score, since criterion is concordance, not length */
  debug(printf("Entering Stage3end_optimal_score with %d hits\n",List_length(terminals)));
  terminals = Stage3end_optimal_score(terminals,query_compress_fwd,query_compress_rev,querylength,
				      /*keep_gmap_p*/true,/*finalp*/false);
  debug(printf("Exiting Stage3end_optimal_score with %d hits\n",List_length(terminals)));
#endif
  
  if (terminals) {
#if 0
    opt_level = ((*found_score_within_trims) < opt_level) ? (*found_score) : opt_level;
    if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
      done_level = user_maxlevel;
    }
#endif
    /* done_level = (*found_score) + subopt_levels; */
    debug(printf("9 (Term)> found_score = %d, opt_level %d, done_level %d\n",
		 *found_score_within_trims,opt_level,done_level));
  }


  Substring_list_gc(&terminals_plus);

  debug(printf("%d terminals\n",List_length(terminals)));

  return terminals;
}


List_T
Terminal_solve_minus (int *found_score_overall, int *found_score_within_trims,
		      List_T queryfwd_minus_set, List_T queryrev_minus_set,

		      Compress_T query_compress_rev, int querylength,
		      int genestrand, Listpool_T listpool,
		      Hitlistpool_T hitlistpool, int level) {
  List_T terminals = NULL, terminals_minus, p;
  Stage3end_T hit;
  int query_lastpos = querylength - index1part;
  int nterminals;
  Substring_T substring;
#ifdef DEBUG13
  int missing_hit, missing_gmap;
#endif


  /* 9 (Term).  Find terminals */
  debug(printf("*** Stage 9 (Term).  Terminal_solve, allowing %d mismatches ***\n",max_terminal_mismatches));

  terminals_minus = (List_T) NULL;

  debug(printf("Starting find_terminals (minus)\n"));
  terminals_minus = find_terminals(terminals_minus,queryfwd_minus_set,
#ifdef DEBUG4E
				   /*queryptr*/queryrc,
#endif
				   querylength,query_lastpos,
				   /*query_compress*/query_compress_rev,
				   hitlistpool,/*plusp*/false,genestrand);

  terminals_minus = find_terminals(terminals_minus,queryrev_minus_set,
#ifdef DEBUG4E
				   /*queryptr*/queryrc,
#endif
				   querylength,query_lastpos,
				   /*query_compress*/query_compress_rev,
				   hitlistpool,/*plusp*/false,genestrand);
  debug(printf("Finished find_terminals (minus)\n"));


  debug4e(printf("Terminals: minus %d\n",List_length(terminals_minus)));

#if 0
  opt_level = ((*found_score_within_trims) < opt_level) ? (*found_score_within_trims) : opt_level;
#endif
  /* done_level = (*found_score_within_trims) + subopt_levels; */

#if 0
  debug4e(printf("Sorting terminals\n"));
  /* TODO: Prioritize inner terminals */
  terminals_left = Substring_sort_nmatches(terminals_left);
  terminals_right = Substring_sort_nmatches(terminals_right);
#endif

  debug(printf("*** Stage 9 (Term)\n"));

  nterminals = 0;
  for (p = terminals_minus; p != NULL && nterminals < MAXTERMINALS; p = List_next(p)) {
    substring = (Substring_T) p->first;
    if ((hit = Stage3end_new_terminal(&(*found_score_overall),&(*found_score_within_trims),
				      substring,querylength,/*gplusp*/false,genestrand,/*sensedir*/SENSE_NULL,
				      listpool,/*method*/TERMINAL,level)) != NULL) {
      terminals = Hitlist_push(terminals,hitlistpool,(void *) hit);
      nterminals++;
    }
  }


  /* Excess distant splicing should be freed already in find_splicepairs_rna */
  terminals = Stage3end_remove_overlaps(terminals,hitlistpool,/*finalp*/false);
#if 0
  /* Do not filter terminals by optimal score, since criterion is concordance, not length */
  debug(printf("Entering Stage3end_optimal_score with %d hits\n",List_length(terminals)));
  terminals = Stage3end_optimal_score(terminals,hitlistpool,querylength,/*finalp*/false);
  debug(printf("Exiting Stage3end_optimal_score with %d hits\n",List_length(terminals)));
#endif
  
  if (terminals) {
#if 0
    opt_level = ((*found_score_within_trims) < opt_level) ? (*found_score_within_trims) : opt_level;
    if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
      done_level = user_maxlevel;
    }
#endif
    /* done_level = (*found_score_within_trims) + subopt_levels; */
    debug(printf("9 (Term)> found_score = %d, opt_level %d, done_level %d\n",
		 *found_score_within_trims,opt_level,done_level));
  }


  Substring_list_gc(&terminals_minus);

  debug(printf("%d terminals\n",List_length(terminals)));

  return terminals;
}


void
Terminal_setup (Univ_IIT_T chromosome_iit_in, int circular_typeint_in,
		Genome_T genomebits_in, Genome_T genomebits_alt_in,
		int index1part_in, int index1interval_in, int subopt_levels_in) {

  chromosome_iit = chromosome_iit_in;
  circular_typeint = circular_typeint_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

  index1part = index1part_in;
  index1interval = index1interval_in;

  subopt_levels = subopt_levels_in;

  return;
}


