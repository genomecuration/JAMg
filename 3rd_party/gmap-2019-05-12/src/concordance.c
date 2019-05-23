static char rcsid[] = "$Id: concordance.c 218690 2019-03-19 17:22:06Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "concordance.h"

#include <stdlib.h>		/* For qsort */
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "assert.h"
#include "mem.h"
#include "univdiag.h"
#include "univdiagdef.h"
#include "transcript.h"
#include "stage3hrdef.h"


static int subopt_levels;
static bool novelsplicingp;

static Chrpos_T pairmax_transcriptome;
static Chrpos_T pairmax_linear;	/* For two ends that both lack a splice */
static Chrpos_T pairmax_circular;

static Chrpos_T expected_pairlength;
static Chrpos_T pairlength_deviation;
static Chrpos_T adjacent_pairlength; /* For two ends, one of which has a splice */

static bool *circularp;
static bool merge_samechr_p;


#define MAX_HITS 1000

#define T Stage3end_T

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Details */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif



#if 0
/* Previously called by stage1hr.c, but now we use Transcript_concordance */
/* Since we haven't called merge procedure yet, we can just look at the Intlist_head */
List_T
Stage3hr_filter_concordant_tr (List_T *disjoint, List_T hits, List_T mates,
			       Hitlistpool_T hitlistpool) {
  List_T common = NULL, p, q;
  T hit, mate;
  int trnum;
  Chrpos_T start1, end2;
  bool concordantp;

  *disjoint = (List_T) NULL;
  for (p = hits; p != NULL; p = List_next(p)) {
    hit = (T) List_head(p);
    trnum = Intlist_head(hit->trnums);
    concordantp = false;

    for (q = mates; q != NULL; q = List_next(q)) {
      mate = (T) List_head(q);
      if (Intlist_head(mate->trnums) == trnum) {
	start1 = (Chrpos_T) Intlist_head(hit->trstarts);
	end2 = (Chrpos_T) Intlist_head(mate->trends);
	if (start1 < end2) {
	  if (end2 <= start1 + pairmax_transcriptome) {
	    concordantp = true;
	  }
	} else {
	  if (start1 <= end2 + pairmax_transcriptome) {
	    concordantp = true;
	  }
	}
      }
    }

    if (concordantp == true) {
      debug7(printf("trnum %d is concordant\n",trnum));
      common = Hitlist_push(common,hitlistpool,(void *) hit);
    } else {
      debug7(printf("trnum %d is discordant\n",trnum));
      *disjoint = Hitlist_push(*disjoint,hitlistpool,(void *) hit);
    }
  }

  debug7(printf("Returning %d common hits\n",List_length(common)));
  return common;
}
#endif


#if 0
static char *
print_sense (int sense) {
  if (sense == SENSE_NULL) {
    return "sense:null";
  } else if (sense == SENSE_ANTI) {
    return "sense:anti";
  } else if (sense == SENSE_FORWARD) {
    return "sense:fwd";
  } else {
    abort();
  }
}
#endif


static int
do_transcriptome_plus (int *concordant_score_overall, List_T *hitpairs,
		       T *hits5, int nhits5, T *hits3, int nhits3, int genestrand,
#if 0
		       char *queryuc_ptr_5, char *queryuc_ptr_3,
		       Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		       Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		       Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		       Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		       Listpool_T listpool, Hitlistpool_T hitlistpool) {
  int nconcordant = 0;
  int i, j;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  int min_insertlength;
  int score;
  Chrnum_T chrnum;

  if (nhits5 > 0 && nhits3 > 0) {
    i = j = 0;
    while (i < nhits5) {
      hit5 = hits5[i];
      chrnum = hit5->effective_chrnum;
      
#ifdef DEBUG
      printf("plus/plus: i=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
	     i,nhits5,hit5->effective_chrnum,hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
	     print_sense(hit5->sensedir_for_concordance),Method_string(hit5->method),hit5->circularalias,hit5);
      if (j >= 0 && j < nhits3) {
	printf("    j=%d/%d #%d:%u..%u %p",j,nhits3,hits3[j]->effective_chrnum,
	       hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,hits3[j]);
      }
      printf("\n");
#endif
      
      while (j >= 0 && hits3[j]->effective_chrnum == chrnum) {
	debug(printf("  backup: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]->circularalias,hits3[j]));
	j--;
      }
      j++;		/* Finish backup */
      
      /* No need for advance */
      
      while (j < nhits3 && hits3[j]->effective_chrnum == chrnum) {
	debug(printf("  samechr: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]->circularalias,hits3[j]));
	hit3 = hits3[j];
	
	if (Transcript_intersect_p(&min_insertlength,hit5->transcripts,hit3->transcripts) == false) {
	  debug(printf(" => transcripts do not intersect:\n"));
	  debug(Transcript_print_nums(hit5->transcripts));
	  debug(Transcript_print_nums(hit5->transcripts));
	  debug(printf("\n"));
	  
	} else {
	  debug(printf(" => concordant transcriptome with min_insertlength %d\n",min_insertlength));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/CONCORDANT,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/true,/*transcriptome_guided_p*/true)) != NULL) {
	    stage3pair->insertlength = min_insertlength;
	    
	    /* No need to update found_score */
	    debug(printf(" Have new pair with scores %d + %d\n",stage3pair->hit5->score,stage3pair->hit3->score));
	    
	    if ((score = (stage3pair->hit5->score_overall + stage3pair->hit3->score_overall)) < *concordant_score_overall) {
	      *concordant_score_overall = score;
	      assert(stage3pair->hit5->nmatches_plus_spliced_trims <= stage3pair->hit5->querylength);
	      assert(stage3pair->hit3->nmatches_plus_spliced_trims <= stage3pair->hit3->querylength);
	      debug(printf(" => Updating concordant_score to be %d = %d + %d\n",
			   *concordant_score_overall,stage3pair->hit5->score_overall,stage3pair->hit3->score_overall));
	    }
	    *hitpairs = Hitlist_push(*hitpairs,hitlistpool,(void *) stage3pair);
	    nconcordant++;
	  }
	}
	debug(printf("\n"));
	
	j++;
      }
      j--;		/* Finish advance */
      
      i++;
    }
  }

  return nconcordant;
}


static int
do_transcriptome_minus (int *concordant_score_overall, List_T *hitpairs,
			T *hits5, int nhits5, T *hits3, int nhits3, int genestrand,
#if 0
			char *queryuc_ptr_5, char *queryuc_ptr_3,
			Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
			Listpool_T listpool, Hitlistpool_T hitlistpool) {
  int nconcordant = 0;
  int i, j;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  int min_insertlength;
  int score;
  Chrnum_T chrnum;

  if (nhits3 > 0 && nhits5 > 0) {
    i = j = 0;
    while (i < nhits3) {
      hit3 = hits3[i];
      chrnum = hit3->effective_chrnum;
      
#ifdef DEBUG
      printf("minus/minus: i=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
	     i,nhits3,hit3->effective_chrnum,hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset,
	     print_sense(hit3->sensedir_for_concordance),Method_string(hit3->method),hit3->circularalias,hit3);
      if (j >= 0 && j < nhits5) {
	printf("    j=%d/%d #%d:%u..%u %p",j,nhits5,hits5[j]->effective_chrnum,
	       hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,hits5[j]);
      }
      printf("\n");
#endif
      
      while (j >= 0 && hits5[j]->effective_chrnum == chrnum) {
	debug(printf("  backup: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]->circularalias,hits5[j]));
	j--;
      }
      j++;			/* Finish backup */
      
      /* No need for advance */
      
      while (j < nhits5 && hits5[j]->effective_chrnum == chrnum) {
	debug(printf("  overlap: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]->circularalias,hits5[j]));
	hit5 = hits5[j];
	
	if (Transcript_intersect_p(&min_insertlength,hit5->transcripts,hit3->transcripts) == false) {
	  debug(printf(" => transcripts do not intersect:\n"));
	  debug(Transcript_print_nums(hit5->transcripts));
	  debug(Transcript_print_nums(hit5->transcripts));
	  debug(printf("\n"));
	  
	} else {
	  debug(printf(" => concordant transcriptome with min_insertlength %d\n",min_insertlength));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/CONCORDANT,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/true,/*transcriptome_guided_p*/true)) != NULL) {
	    stage3pair->insertlength = min_insertlength;
	    
	    if ((score = (stage3pair->hit5->score_overall + stage3pair->hit3->score_overall)) < *concordant_score_overall) {
	      *concordant_score_overall = score;
	      assert(stage3pair->hit5->nmatches_plus_spliced_trims <= stage3pair->hit5->querylength);
	      assert(stage3pair->hit3->nmatches_plus_spliced_trims <= stage3pair->hit3->querylength);
	      debug(printf(" => Updating concordant_score to be %d = %d + %d\n",
			   *concordant_score_overall,stage3pair->hit5->score_overall,stage3pair->hit3->score_overall));
	    }
	    *hitpairs = Hitlist_push(*hitpairs,hitlistpool,(void *) stage3pair);
	    nconcordant++;
	  }
	}
	debug(printf("\n"));
	
	j++;
      }
      j--;		/* Finish advance */
      
      i++;
    }
  }

  return nconcordant;
}



List_T
Concordance_pair_up_transcriptome (bool *abort_pairing_p, int *concordant_score_overall, List_T hitpairs,

				   Ladder_T ladder5_plus, Ladder_T ladder5_minus,
				   Ladder_T ladder3_plus, Ladder_T ladder3_minus,
#if 0
				   char *queryuc_ptr_5, char *queryuc_ptr_3,
				   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
				   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
				   Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				   Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
				   Listpool_T listpool, Hitlistpool_T hitlistpool,
				   int maxpairedpaths, int genestrand) {
  int nconcordant = 0, n;  /* was List_length(hitpairs), but this hurts minus if there are too many plus hits */

  T *hits5, *hits3;
  int nhits5, nhits3;
  int max_frontier_score, frontier_score, cutoff_level, level, score5, score3;


  debug(printf("Starting Concordance_pair_up_transcriptome\n"));

  /* Initial value, which resets to frontier_score + subopt_levels upon our first hitpair */
  cutoff_level = Ladder_cutoff(ladder5_plus) + Ladder_cutoff(ladder3_plus);
  if ((level = Ladder_cutoff(ladder5_minus) + Ladder_cutoff(ladder3_minus)) > cutoff_level) {
    cutoff_level = level;
  }
  max_frontier_score = cutoff_level;

  frontier_score = 0;
  while (*abort_pairing_p == false && frontier_score <= cutoff_level) {
    debug(printf("frontier_score = %d\n",frontier_score));
    for (score5 = 0; score5 <= frontier_score; score5++) {
      score3 = frontier_score - score5;
      debug(printf("score5 = %d, score3 = %d\n",score5,score3));

      /* plus/plus: hits5_plus against hits3_plus (really on minus) */
      if (score5 <= Ladder_cutoff(ladder5_plus) && score3 <= Ladder_cutoff(ladder3_plus)) {
	hits5 = Ladder_hits_for_score(&nhits5,ladder5_plus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,ladder3_plus,hitlistpool,score3);
	debug(printf("at score %d+%d, nhits5_plus = %d, nhits3_plus = %d\n",score5,score3,nhits5,nhits3));

	if ((n = do_transcriptome_plus(&(*concordant_score_overall),&hitpairs,
				       hits5,nhits5,hits3,nhits3,genestrand,
#if 0
				       queryuc_ptr_5,queryuc_ptr_3,
				       query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				       pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
				       listpool,hitlistpool)) > 0) {
	  if (nconcordant == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
	  if ((nconcordant += n) > maxpairedpaths) {
	    debug(printf(" -- %d concordant paths exceeds %d",nconcordant,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }

      /* minus/minus: hits3_minus (really on plus) against hits5_minus */
      if (score5 <= Ladder_cutoff(ladder5_minus) && score3 <= Ladder_cutoff(ladder3_minus)) {
	hits5 = Ladder_hits_for_score(&nhits5,ladder5_minus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,ladder3_minus,hitlistpool,score3);
	debug(printf("at score %d+%d, nhits5_minus = %d, nhits3_minus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_transcriptome_minus(&(*concordant_score_overall),&hitpairs,
					hits5,nhits5,hits3,nhits3,genestrand,
#if 0
					queryuc_ptr_5,queryuc_ptr_3,
					query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					listpool,hitlistpool)) > 0) {
	  if (nconcordant == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
	  if ((nconcordant += n) > maxpairedpaths) {
	    debug(printf(" -- %d concordant paths exceeds %d",nconcordant,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }
    }

    frontier_score++;
  }

  debug(printf("Finished with Concordance_pair_up_transcriptome: %d concordant\n",List_length(hitpairs)));

  return hitpairs;
}


static int
do_genome_plus (int *adjacent_score, int *concordant_score_overall, int *concordant_score_within_trims,
		List_T *local_hitpairs, List_T *distant_hitpairs, T *hits5, int nhits5, T *hits3, int nhits3, int genestrand,
		int querylength5, int querylength3, 
#if 0
		char *queryuc_ptr_5, char *queryuc_ptr_3,
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		Listpool_T listpool, Hitlistpool_T hitlistpool) {

  int nconcordant = 0;
  int i, j;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  Univcoord_T insert_start;
  int score;
  int pairmax;

  *concordant_score_within_trims = querylength5 + querylength3;

  if (nhits5 > 0 && nhits3 > 0) {
    i = j = 0;
    while (i < nhits5) {
      hit5 = hits5[i];
      if (circularp[hit5->effective_chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }
      
      insert_start = hit5->genomicend - querylength5;
#ifdef DEBUG
      printf("plus/plus: i=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
	     i,nhits5,hit5->effective_chrnum,hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
	     print_sense(hit5->sensedir_for_concordance),Method_string(hit5->method),hit5->circularalias,hit5);
      if (j >= 0 && j < nhits3) {
	printf("    j=%d/%d #%d:%u..%u %p",j,nhits3,hits3[j]->effective_chrnum,
	       hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,hits3[j]);
      }
      printf("\n");
#endif

      while (j >= 0 && 
	     hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	debug(printf("  backup: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]->circularalias,hits3[j]));
	j--;
      }
      j++;		/* Finish backup */

      while (j < nhits3 && 
	     hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	debug(printf("  advance: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]->circularalias,hits3[j]));
	j++;
      }
      
      while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	debug(printf("  overlap: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]->circularalias,hits3[j]));
	hit3 = hits3[j];
	
	if (hit5->effective_chrnum != hit3->effective_chrnum) {
	  debug(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
	  
	} else if (SENSE_INCONSISTENT_P(hit5->sensedir_for_concordance,hit3->sensedir_for_concordance)) {
	  /* Use sensedir_for_concordance here and not sensedir */
	  debug(printf(" => sense inconsistent: %d | %d = %d",hit5->sensedir_for_concordance,hit3->sensedir_for_concordance,hit5->sensedir_for_concordance|hit3->sensedir_for_concordance));
	  
	} else if (hit3->genomicend < hit5->genomicstart) {
	  debug(printf(" => scramble because end3 %llu < start5 %llu\n",
		       (unsigned long long) hit3->genomicend,(unsigned long long) hit5->genomicstart));
	  
	} else {
	  debug(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
		       hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/CONCORDANT,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/true,/*transcriptome_guided_p*/false)) != NULL) {
	    if (Stage3pair_distant_splice_p(stage3pair) == true) {
	      *distant_hitpairs = Hitlist_push(*distant_hitpairs,hitlistpool,(void *) stage3pair);
	    } else {
	      if ((score = stage3pair->hit5->score_overall + stage3pair->hit3->score_overall) < *concordant_score_overall) {
		debug(printf(" => Updating concordant_score_overall to be %d\n",score));
		*concordant_score_overall = score;
	      }
	      if ((score = stage3pair->hit5->score_within_trims + stage3pair->hit3->score_within_trims) < *concordant_score_within_trims) {
		*concordant_score_within_trims = score;
		assert(stage3pair->hit5->nmatches_plus_spliced_trims <= stage3pair->hit5->querylength);
		assert(stage3pair->hit3->nmatches_plus_spliced_trims <= stage3pair->hit3->querylength);
		debug(printf(" => Updating concordant_score_within_trims to be %d = %d + %d\n",
			     *concordant_score_within_trims,stage3pair->hit5->score_within_trims,
			     stage3pair->hit3->score_within_trims));
	      }
	      if (stage3pair->insertlength <= adjacent_pairlength && score < *adjacent_score) {
		*adjacent_score = score;
	      }
	      *local_hitpairs = Hitlist_push(*local_hitpairs,hitlistpool,(void *) stage3pair);
	      nconcordant++;
	    }
	  }
	}
	debug(printf("\n"));
	
	j++;
      }
      j--;		/* Finish advance */
      
      i++;
    }
  }

  return nconcordant;
}



static int
do_genome_minus (int *adjacent_score, int *concordant_score_overall, int *concordant_score_within_trims,
		 List_T *local_hitpairs, List_T *distant_hitpairs,
		 T *hits5, int nhits5, T *hits3, int nhits3, int genestrand,
		 int querylength5, int querylength3,
#if 0
		 char *queryuc_ptr_5, char *queryuc_ptr_3,
		 Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		 Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		 Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		 Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		 Listpool_T listpool, Hitlistpool_T hitlistpool) {

  *concordant_score_within_trims = querylength5 + querylength3;

  int nconcordant = 0;
  int i, j;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  Univcoord_T insert_start;
  int score;
  int pairmax;


  if (nhits3 > 0 && nhits5 > 0) {
    i = j = 0;
    while (i < nhits3) {
      hit3 = hits3[i];
      if (circularp[hit3->effective_chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }
      insert_start = hit3->genomicstart - querylength3;
#ifdef DEBUG
      printf("minus/minus: i=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
	     i,nhits3,hit3->effective_chrnum,hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset,
	     print_sense(hit3->sensedir_for_concordance),Method_string(hit3->method),hit3->circularalias,hit3);
      if (j >= 0 && j < nhits5) {
	printf("    j=%d/%d #%d:%u..%u %p",j,nhits5,hits5[j]->effective_chrnum,
	       hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,hits5[j]);
      }
      printf("\n");
#endif

      while (j >= 0 && 
	     hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	debug(printf("  backup: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]->circularalias,hits5[j]));
	j--;
      }
      j++;			/* Finish backup */
      
      while (j < nhits5 && 
	     hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	debug(printf("  advance: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]->circularalias,hits5[j]));
	j++;
      }
      
      while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	debug(printf("  overlap: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]->circularalias,hits5[j]));
	hit5 = hits5[j];

	if (hit3->effective_chrnum != hit5->effective_chrnum) {
	  debug(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));

	} else if (SENSE_INCONSISTENT_P(hit3->sensedir_for_concordance,hit5->sensedir_for_concordance)) {
	  /* Use sensedir_for_concordance here and not sensedir */
	  debug(printf(" => sense inconsistent: %d | %d = %d",hit5->sensedir_for_concordance,hit3->sensedir_for_concordance,hit5->sensedir_for_concordance|hit3->sensedir_for_concordance));
		
	} else if (hit5->genomicstart < hit3->genomicend) {
	  debug(printf(" => scramble because start5 %llu < end3 %llu\n",
		       (unsigned long long) hit5->genomicstart,(unsigned long long) hit3->genomicend));
		
	} else {
	  debug(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
		       hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/CONCORDANT,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/true,/*transcriptome_guided_p*/false)) != NULL) {
	    if (Stage3pair_distant_splice_p(stage3pair) == true) {
	      *distant_hitpairs = Hitlist_push(*distant_hitpairs,hitlistpool,(void *) stage3pair);
	    } else {
	      if ((score = stage3pair->hit5->score_overall + stage3pair->hit3->score_overall) < *concordant_score_overall) {
		debug(printf(" => Updating concordant_score_overall to be %d\n",score));
		*concordant_score_overall = score;
	      }
	      if ((score = stage3pair->hit5->score_within_trims + stage3pair->hit3->score_within_trims) < *concordant_score_within_trims) {
		*concordant_score_within_trims = score;
		assert(stage3pair->hit5->nmatches_plus_spliced_trims <= stage3pair->hit5->querylength);
		assert(stage3pair->hit3->nmatches_plus_spliced_trims <= stage3pair->hit3->querylength);
		debug(printf(" => Updating concordant_score_within_trims to be %d = %d + %d\n",
			     *concordant_score_within_trims,stage3pair->hit5->score_within_trims,stage3pair->hit3->score_within_trims));
	      }
	      if (stage3pair->insertlength <= adjacent_pairlength && score < *adjacent_score) {
		*adjacent_score = score;
	      }
	      *local_hitpairs = Hitlist_push(*local_hitpairs,hitlistpool,(void *) stage3pair);
	      nconcordant++;
	    }
	  }
	}
	debug(printf("\n"));
	      
	j++;
      }
      j--;		/* Finish advance */
	    
      i++;
    }
  }

  return nconcordant;
}



/* Finds concordant pairs if nconcordant is 0 */
List_T
Concordance_pair_up_genome (bool *abort_pairing_p, int *adjacent_score, int *concordant_score_overall,
			    List_T *distant_hitpairs, List_T hitpairs,

			    List_T hitlist5_gplus, List_T hitlist5_gminus,
			    List_T hitlist3_gplus, List_T hitlist3_gminus,
			    
			    Ladder_T ladder5_plus, Ladder_T ladder5_minus,
			    Ladder_T ladder3_plus, Ladder_T ladder3_minus,
			    int querylength5, int querylength3,
#if 0
			    char *queryuc_ptr_5, char *queryuc_ptr_3,
			    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			    Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
			    Listpool_T listpool, Hitlistpool_T hitlistpool,
			    int maxpairedpaths, int genestrand) {
  int nconcordant = 0, n; /* Also used as an indicator of when to reset cutoff_level */
  Ladder_T newladder5_plus, newladder5_minus, newladder3_plus, newladder3_minus;

  T *hits5, *hits3;
  int nhits5, nhits3;

  int frontier_score, score5, score3;
  int concordant_score_within_trims;
  int cutoff_level, level;
  

  debug(printf("Starting Concordance_pair_up_genome\n"));

#if 0
  /* Rely instead on Ladder_cutoff, in case there are a few good hits */
  if (0 && List_length(hitlist5_gplus) > MAX_HITS) {
    Stage3end_gc(hitlist5_gplus);
    newladder5_plus = Ladder_new((List_T) NULL,hitlistpool,/*end5p*/true);
  } else {
    newladder5_plus = Ladder_new(hitlist5_gplus,hitlistpool,/*end5p*/true);
    Ladder_merge(ladder5_plus,newladder5_plus);
  }

  if (0 && List_length(hitlist5_gminus) > MAX_HITS) {
    Stage3end_gc(hitlist5_gminus);
    newladder5_minus = Ladder_new((List_T) NULL,hitlistpool,/*end5p*/true);
  } else {
    newladder5_minus = Ladder_new(hitlist5_gminus,hitlistpool,/*end5p*/true);
    Ladder_merge(ladder5_minus,newladder5_minus);
  }

  if (0 && List_length(hitlist3_gplus) > MAX_HITS) {
    Stage3end_gc(hitlist3_gplus);
    newladder3_plus = Ladder_new((List_T) NULL,hitlistpool,/*end5p*/false);
  } else {
    newladder3_plus = Ladder_new(hitlist3_gplus,hitlistpool,/*end5p*/false);
    Ladder_merge(ladder3_plus,newladder3_plus);
  }

  if (0 && List_length(hitlist3_gminus) > MAX_HITS) {
    Stage3end_gc(hitlist3_gminus);
    newladder3_minus = Ladder_new((List_T) NULL,hitlistpool,/*end5p*/false);
  } else {
    newladder3_minus = Ladder_new(hitlist3_gminus,hitlistpool,/*end5p*/false);
    Ladder_merge(ladder3_minus,newladder3_minus);
  }

#else
  newladder5_plus = Ladder_new(hitlist5_gplus,hitlistpool,/*end5p*/true);
  newladder5_minus = Ladder_new(hitlist5_gminus,hitlistpool,/*end5p*/true);
  newladder3_plus = Ladder_new(hitlist3_gplus,hitlistpool,/*end5p*/false);
  newladder3_minus = Ladder_new(hitlist3_gminus,hitlistpool,/*end5p*/false);

  Ladder_merge(ladder5_plus,newladder5_plus,hitlistpool);
  Ladder_merge(ladder5_minus,newladder5_minus,hitlistpool);
  Ladder_merge(ladder3_plus,newladder3_plus,hitlistpool);
  Ladder_merge(ladder3_minus,newladder3_minus,hitlistpool);
#endif

  
  /* Initial value, which resets to frontier_score + subopt_levels upon our first hitpair */
  cutoff_level = Ladder_cutoff(newladder5_plus) + Ladder_cutoff(ladder3_plus);
  if ((level = Ladder_cutoff(newladder5_minus) + Ladder_cutoff(ladder3_minus)) > cutoff_level) {
    cutoff_level = level;
  }
  if ((level = Ladder_cutoff(ladder5_plus) + Ladder_cutoff(newladder3_plus)) > cutoff_level) {
    cutoff_level = level;
  }
  if ((level = Ladder_cutoff(ladder5_minus) + Ladder_cutoff(newladder3_minus)) > cutoff_level) {
    cutoff_level = level;
  }
  /* max_frontier_score = cutoff_level; */

  frontier_score = 0;
  /* Eliminating check on frontier_score against *found_score, because too greedy and can miss correct result */
  while (*abort_pairing_p == false && frontier_score <= cutoff_level) {
    debug(printf("frontier_score = %d, cutoff_level %d\n",frontier_score,cutoff_level));
    for (score5 = 0; score5 <= frontier_score; score5++) {
      score3 = frontier_score - score5;
      debug(printf("score5 = %d, score3 = %d\n",score5,score3));

      /* New hits 5 vs all hits 3 */
      /* plus/plus: hits5_plus against hits3_plus (really on minus) */
      if (score5 <= Ladder_cutoff(newladder5_plus) && score3 <= Ladder_cutoff(ladder3_plus)) {
	hits5 = Ladder_hits_for_score(&nhits5,newladder5_plus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,ladder3_plus,hitlistpool,score3);
	debug(printf("at score %d+%d, nnewhits5_plus = %d, nhits3_plus = %d\n",score5,score3,nhits5,nhits3));

	if ((n = do_genome_plus(&(*adjacent_score),&(*concordant_score_overall),&concordant_score_within_trims,
				&hitpairs,&(*distant_hitpairs),
				hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
				queryuc_ptr_5,queryuc_ptr_3,
				query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
				listpool,hitlistpool)) > 0) {
#if 0
	  /* Too greedy */
	  if (nconcordant == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
#else
	  if (concordant_score_within_trims + subopt_levels < cutoff_level) {
	    cutoff_level = concordant_score_within_trims + subopt_levels;
	    debug(printf("Updated cutoff level to be %d\n",cutoff_level));
	  }
#endif
	  if ((nconcordant += n) > maxpairedpaths) {
	    debug(printf(" -- %d concordant paths exceeds %d",nconcordant,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }

      /* minus/minus: hits3_minus (really on plus) against newhits5_minus */
      if (score5 <= Ladder_cutoff(newladder5_minus) && score3 <= Ladder_cutoff(ladder3_minus)) {
	hits5 = Ladder_hits_for_score(&nhits5,newladder5_minus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,ladder3_minus,hitlistpool,score3);
	debug(printf("at score %d+%d, nnewhits5_minus = %d, nhits3_minus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_genome_minus(&(*adjacent_score),&(*concordant_score_overall),&concordant_score_within_trims,
				 &hitpairs,&(*distant_hitpairs),
				 hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
				 queryuc_ptr_5,queryuc_ptr_3,
				 query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				 pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif				 
				 listpool,hitlistpool)) > 0) {
#if 0
	  /* Too greedy */
	  if (nconcordant == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
#else
	  if (concordant_score_within_trims + subopt_levels < cutoff_level) {
	    cutoff_level = concordant_score_within_trims + subopt_levels;
	    debug(printf("Updated cutoff level to be %d\n",cutoff_level));
	  }
#endif
	  if ((nconcordant += n) > maxpairedpaths) {
	    debug(printf(" -- %d concordant paths exceeds %d",nconcordant,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }

      /* All hits 5 vs new hits 3 */
      /* plus/plus: hits5_plus against newhits3_plus (really on minus) */
      if (score5 <= Ladder_cutoff(ladder5_plus) && score3 <= Ladder_cutoff(newladder3_plus)) {
	hits5 = Ladder_hits_for_score(&nhits5,ladder5_plus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,newladder3_plus,hitlistpool,score3);
	debug(printf("at score %d+%d, nhits5_plus = %d, nnewhits3_plus = %d\n",score5,score3,nhits5,nhits3));

	if ((n = do_genome_plus(&(*adjacent_score),&(*concordant_score_overall),&concordant_score_within_trims,
				&hitpairs,&(*distant_hitpairs),
				hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
				queryuc_ptr_5,queryuc_ptr_3,
				query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
				listpool,hitlistpool)) > 0) {
#if 0
	  /* Too greedy */
	  if (nconcordant == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
#else
	  if (concordant_score_within_trims + subopt_levels < cutoff_level) {
	    cutoff_level = concordant_score_within_trims + subopt_levels;
	    debug(printf("Updated cutoff level to be %d\n",cutoff_level));
	  }
#endif
	  if ((nconcordant += n) > maxpairedpaths) {
	    debug(printf(" -- %d concordant paths exceeds %d",nconcordant,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }

      /* minus/minus: newhits3_minus (really on plus) against hits5_minus */
      if (score5 <= Ladder_cutoff(ladder5_minus) && score3 <= Ladder_cutoff(newladder3_minus)) {
	hits5 = Ladder_hits_for_score(&nhits5,ladder5_minus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,newladder3_minus,hitlistpool,score3);

	debug(printf("at score %d+%d, nhits5_minus = %d, nnewhits3_minus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_genome_minus(&(*adjacent_score),&(*concordant_score_overall),&concordant_score_within_trims,
				 &hitpairs,&(*distant_hitpairs),
				 hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
				 queryuc_ptr_5,queryuc_ptr_3,
				 query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				 pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
				 listpool,hitlistpool)) > 0) {
#if 0
	  /* Too greedy */
	  if (nconcordant == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
#else
	  if (concordant_score_within_trims + subopt_levels < cutoff_level) {
	    cutoff_level = concordant_score_within_trims + subopt_levels;
	    debug(printf("Updated cutoff level to be %d\n",cutoff_level));
	  }
#endif
	  if ((nconcordant += n) > maxpairedpaths) {
	    debug(printf(" -- %d concordant paths exceeds %d",nconcordant,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }
    }

    frontier_score++;
  }


  Ladder_free(&newladder5_plus);
  Ladder_free(&newladder5_minus);
  Ladder_free(&newladder3_plus);
  Ladder_free(&newladder3_minus);

  Ladder_gc_duplicates(ladder5_plus);
  Ladder_gc_duplicates(ladder5_minus);
  Ladder_gc_duplicates(ladder3_plus);
  Ladder_gc_duplicates(ladder3_minus);

  debug(printf("Finished with Concordance_pair_up_genome: %d concordant\n",List_length(hitpairs)));

  return hitpairs;
}


static int
do_distant_plus_plus (int *concordant_score, List_T *hitpairs,
		      List_T *conc_transloc, List_T *samechr,
		      T *hits5, int nhits5, T *hits3, int nhits3, int genestrand,
		      int querylength5, int querylength3,
#if 0
		      char *queryuc_ptr_5, char *queryuc_ptr_3,
		      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		      Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		      Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		      Listpool_T listpool, Hitlistpool_T hitlistpool) {
  int nfound = 0;
  int i, j;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  int score;
  Univcoord_T insert_start;
  Chrpos_T pairmax;

  if (nhits5 > 0 && nhits3 > 0) {
    i = j = 0;
    while (i < nhits5) {
      hit5 = hits5[i];
      if (circularp[hit5->effective_chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }
      
      insert_start = hit5->genomicend - querylength5;
#ifdef DEBUG
      printf("plus/plus: i=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
	     i,nhits5,hit5->effective_chrnum,hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
	     print_sense(hit5->sensedir_for_concordance),Method_string(hit5->method),hit5->circularalias,hit5);
      if (j >= 0 && j < nhits3) {
	printf("    j=%d/%d #%d:%u..%u %p",j,nhits3,hits3[j]->effective_chrnum,
	       hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,hits3[j]);
      }
      printf("\n");
#endif
      
      while (j >= 0 && 
	     hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	debug(printf("  backup: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]->circularalias,hits3[j]));
	j--;
      }
      j++;		/* Finish backup */
      
      while (j < nhits3 && 
	     hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	debug(printf("  advance: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]->circularalias,hits3[j]));
	j++;
      }
      
      while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	debug(printf("  overlap: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]->circularalias,hits3[j]));
	hit3 = hits3[j];
	
	if (hit5->effective_chrnum != hit3->effective_chrnum) {
	  debug(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
	  
	} else if (hit5->distant_splice_p == true && hit3->distant_splice_p == true /* && hit5->other_chrnum != hit3->other_chrnum */) {
	  /* Could potentially miss an alignment if the two ends overlap */
	  debug(printf(" => double splice translocations"));
	  
	} else if (hit5->distant_splice_p == true || hit3->distant_splice_p == true) {
	  debug(printf(" => conc_transloc effchr %d (chrnum5 %d, chrnum3 %d)",
		       hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/CONCORDANT_TRANSLOCATIONS,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/true,/*transcriptome_guided_p*/false)) != NULL) {
	    *conc_transloc = Hitlist_push(*conc_transloc,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
	  
	} else if (SENSE_INCONSISTENT_P(hit5->sensedir_for_concordance,hit3->sensedir_for_concordance)) {
	  /* Use sensedir_for_concordance here and not sensedir */
	  debug(printf(" => sense inconsistent: %d | %d = %d",hit5->sensedir_for_concordance,hit3->sensedir_for_concordance,hit5->sensedir_for_concordance|hit3->sensedir_for_concordance));
	  
	} else if (hit3->genomicend < hit5->genomicstart) {
	  debug(printf(" => scramble because end3 %llu < start5 %llu\n",
		       (unsigned long long) hit3->genomicend,(unsigned long long) hit5->genomicstart));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/PAIRED_SCRAMBLE,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/false,/*transcriptome_guided_p*/false)) != NULL) {
	    *samechr = Hitlist_push(*samechr,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
	  
	} else {
	  /* Unexpected? */
	  debug(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
		       hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/CONCORDANT,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/true,/*transcriptome_guided_p*/false)) != NULL) {
	    
	    if ((score = (stage3pair->hit5->querylength - stage3pair->hit5->nmatches_plus_spliced_trims) + 
		 (stage3pair->hit3->querylength - stage3pair->hit3->nmatches_plus_spliced_trims)) < *concordant_score) {
	      *concordant_score = score;
	      assert(stage3pair->hit5->nmatches_plus_spliced_trims <= stage3pair->hit5->querylength);
	      assert(stage3pair->hit3->nmatches_plus_spliced_trims <= stage3pair->hit3->querylength);
	      debug(printf(" => Updating concordant_score to be %d = %d + %d\n",
			   *concordant_score,stage3pair->hit5->querylength - stage3pair->hit5->nmatches_plus_spliced_trims,
			   stage3pair->hit3->querylength - stage3pair->hit3->nmatches_plus_spliced_trims));
	    }
	    *hitpairs = Hitlist_push(*hitpairs,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
	}
	debug(printf("\n"));
	
	j++;
      }
      j--;		/* Finish advance */
      
      i++;
    }
  }
  
  return nfound;
}


static int
do_distant_minus_minus (int *concordant_score, List_T *hitpairs,
			List_T *conc_transloc, List_T *samechr,
			T *hits5, int nhits5, T *hits3, int nhits3, int genestrand,
			int querylength5, int querylength3,
#if 0
			char *queryuc_ptr_5, char *queryuc_ptr_3,
			Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
			Listpool_T listpool, Hitlistpool_T hitlistpool) {

  int nfound = 0;
  int i, j;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  int score;
  Univcoord_T insert_start;
  Chrpos_T pairmax;

  if (nhits3 > 0 && nhits5 > 0) {
    i = j = 0;
    while (i < nhits3) {
      hit3 = hits3[i];
      if (circularp[hit3->effective_chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }
      insert_start = hit3->genomicstart - querylength3;
#ifdef DEBUG
      printf("minus/minus: i=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
	     i,nhits3,hit3->effective_chrnum,hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset,
	     print_sense(hit3->sensedir_for_concordance),Method_string(hit3->method),hit3->circularalias,hit3);
      if (j >= 0 && j < nhits5) {
	printf("    j=%d/%d #%d:%u..%u %p",j,nhits5,hits5[j]->effective_chrnum,
	       hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,hits5[j]);
      }
      printf("\n");
#endif
      
      while (j >= 0 && 
	     hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	debug(printf("  backup: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]->circularalias,hits5[j]));
	j--;
      }
      j++;			/* Finish backup */
      
      while (j < nhits5 && 
	     hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	debug(printf("  advance: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p\n",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]->circularalias,hits5[j]));
	j++;
      }
      
      while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	debug(printf("  overlap: j=%d/%d #%d:%u..%u %s %s circularalias:%d %p",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]->circularalias,hits5[j]));
	hit5 = hits5[j];
	
	/* Do want to see pairs previously seen */
	if (hit3->effective_chrnum != hit5->effective_chrnum) {
	  debug(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
	  
	} else if (hit5->distant_splice_p == true && hit3->distant_splice_p == true /* && hit5->other_chrnum != hit3->other_chrnum */) {
	  /* Could potentially miss an alignment if the two ends overlap */
	  debug(printf(" => double splice translocations"));
	  
	} else if (hit5->distant_splice_p == true || hit3->distant_splice_p == true) {
	  debug(printf(" => conc_transloc effchr %d (chrnum5 %d, chrnum3 %d)",
		       hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/CONCORDANT_TRANSLOCATIONS,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/true,/*transcriptome_guided_p*/false)) != NULL) {
	    *conc_transloc = Hitlist_push(*conc_transloc,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
		
	} else if (SENSE_INCONSISTENT_P(hit3->sensedir_for_concordance,hit5->sensedir_for_concordance)) {
	  /* Use sensedir_for_concordance here and not sensedir */
	  debug(printf(" => sense inconsistent: %d | %d = %d",hit5->sensedir_for_concordance,hit3->sensedir_for_concordance,hit5->sensedir_for_concordance|hit3->sensedir_for_concordance));
		
	} else if (hit5->genomicstart < hit3->genomicend) {
	  debug(printf(" => scramble because start5 %llu < end3 %llu\n",
		       (unsigned long long) hit5->genomicstart,(unsigned long long) hit3->genomicend));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/PAIRED_SCRAMBLE,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/false,/*transcriptome_guided_p*/false)) != NULL) {
	    *samechr = Hitlist_push(*samechr,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
		
	} else {
	  /* Unexpected? */
	  debug(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
		       hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/CONCORDANT,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/true,/*transcriptome_guided_p*/false)) != NULL) {
		  
	    if ((score = (stage3pair->hit5->querylength - stage3pair->hit5->nmatches_plus_spliced_trims) + 
		 (stage3pair->hit3->querylength - stage3pair->hit3->nmatches_plus_spliced_trims)) < *concordant_score) {
	      *concordant_score = score;
	      assert(stage3pair->hit5->nmatches_plus_spliced_trims <= stage3pair->hit5->querylength);
	      assert(stage3pair->hit3->nmatches_plus_spliced_trims <= stage3pair->hit3->querylength);
	      debug(printf(" => Updating concordant_score to be %d = %d + %d\n",
			   *concordant_score,stage3pair->hit5->querylength - stage3pair->hit5->nmatches_plus_spliced_trims,
			   stage3pair->hit3->querylength - stage3pair->hit3->nmatches_plus_spliced_trims));
	    }
	    *hitpairs = Hitlist_push(*hitpairs,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
	}
	debug(printf("\n"));
	      
	j++;
      }
      j--;		/* Finish advance */
	    
      i++;
    }
  }
  
  return nfound;
}


static int
do_distant_plus_minus (List_T *samechr,
		       T *hits5, int nhits5, T *hits3, int nhits3, int genestrand,
		       int querylength5, int querylength3,
#if 0
		       char *queryuc_ptr_5, char *queryuc_ptr_3,
		       Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		       Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		       Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		       Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		       Listpool_T listpool, Hitlistpool_T hitlistpool) {

  int nfound = 0;
  int i, j;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  Univcoord_T insert_start;
  Chrpos_T pairmax;

  if (nhits5 > 0 && nhits3 > 0) {
    i = j = 0;
    while (i < nhits5) {
      hit5 = hits5[i];
      if (circularp[hit5->effective_chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }
      insert_start = hit5->genomicend - querylength5;
#ifdef DEBUG
      printf("plus/minus: i=%d/%d #%d:%u..%u %s %s %p",
	     i,nhits5,hit5->effective_chrnum,hit5->genomicstart - hit5->chroffset,hit5->genomicend - hit5->chroffset,
	     print_sense(hit5->sensedir_for_concordance),Method_string(hit5->method),hit5);
      if (j >= 0 && j < nhits3) {
	printf("    j=%d/%d #%d:%u..%u %p",j,nhits3,hits3[j]->effective_chrnum,
	       hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,hits3[j]);
      }
      printf("\n");
#endif
	    
      while (j >= 0 && 
	     hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	debug(printf("  backup: j=%d/%d #%d:%u..%u %s %s %p\n",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]));
	j--;
      }
      j++;		/* Finish backup */

      while (j < nhits3 && 
	     hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	debug(printf("  advance: j=%d/%d #%d:%u..%u %s %s %p\n",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]));
	j++;
      }

      while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	debug(printf("  overlap: j=%d/%d #%d:%u..%u %s %s %p",
		     j,nhits3,hits3[j]->effective_chrnum,hits3[j]->genomicstart - hits3[j]->chroffset,hits3[j]->genomicend - hits3[j]->chroffset,
		     print_sense(hits3[j]->sensedir_for_concordance),Method_string(hits3[j]->method),hits3[j]));
	hit3 = hits3[j];
		
	if (hit5->effective_chrnum != hit3->effective_chrnum) {
	  debug(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		
	} else if (hit5->distant_splice_p == true && hit3->distant_splice_p == true /* && hit5->other_chrnum != hit3->other_chrnum */) {
	  debug(printf(" => double splice translocations"));
		
	} else if (SENSE_INCONSISTENT_FOR_INVERSION_P(hit5->sensedir_for_concordance,hit3->sensedir_for_concordance)) {
	  /* Use sensedir_for_concordance here and not sensedir */
	  debug(printf(" => sense inconsistent for inversion"));
#if 0
	} else if (hits3[j]->genomicstart + querylength3 <= insert_start) {
	  debug(printf(" => scramble"));
	  if (nsamechr <= maxpairedpaths &&
	      (stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/PAIRED_SCRAMBLE,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/false,/*transcriptome_guided_p*/false)) != NULL) {
	    *samechr = Hitlist_push(*samechr,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
#endif
	} else {
	  debug(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
		       hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/PAIRED_INVERSION,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/false,/*transcriptome_guided_p*/false)) != NULL) {
	    *samechr = Hitlist_push(*samechr,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
	}
	debug(printf("\n"));

	j++;
      }
      j--;		/* Finish advance */

      i++;
    }
  }

  return nfound;
}


static int
do_distant_minus_plus (List_T *samechr,
		       T *hits5, int nhits5, T *hits3, int nhits3, int genestrand,
		       int querylength5, int querylength3,
#if 0
		       char *queryuc_ptr_5, char *queryuc_ptr_3,
		       Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		       Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		       Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		       Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
		       Listpool_T listpool, Hitlistpool_T hitlistpool) {
  int nfound = 0;
  int i, j;
  Stage3pair_T stage3pair;
  T hit5, hit3;
  Univcoord_T insert_start;
  Chrpos_T pairmax;

  if (nhits3 > 0 && nhits5 > 0) {
    i = j = 0;
    while (i < nhits3) {
      hit3 = hits3[i];
      if (circularp[hit3->effective_chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }
      insert_start = hit3->genomicstart - querylength3;
#ifdef DEBUG
      printf("minus/plus: i=%d/%d #%d:%u..%u %s %s %p",
	     i,nhits3,hit3->effective_chrnum,hit3->genomicstart - hit3->chroffset,hit3->genomicend - hit3->chroffset,
	     print_sense(hit3->sensedir_for_concordance),Method_string(hit3->method),hit3);
      if (j >= 0 && j < nhits5) {
	printf("    j=%d/%d #%d:%u..%u %p",j,nhits5,hits5[j]->effective_chrnum,
	       hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,hits5[j]);
      }
      printf("\n");
#endif

      while (j >= 0 && 
	     hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	debug(printf("  backup: j=%d/%d #%d:%u..%u %s %s %p\n",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]));
	j--;
      }
      j++;			/* Finish backup */

      while (j < nhits5 && 
	     hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	debug(printf("  advance: j=%d/%d #%d:%u..%u %s %s %p\n",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]));
	j++;
      }

      while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	debug(printf("  overlap: j=%d/%d #%d:%u..%u %s %s %p",
		     j,nhits5,hits5[j]->effective_chrnum,hits5[j]->genomicstart - hits5[j]->chroffset,hits5[j]->genomicend - hits5[j]->chroffset,
		     print_sense(hits5[j]->sensedir_for_concordance),Method_string(hits5[j]->method),hits5[j]));
	hit5 = hits5[j];

	if (hit3->effective_chrnum != hit5->effective_chrnum) {
	  debug(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		
	} else if (hit5->distant_splice_p == true && hit3->distant_splice_p == true /* && hit5->other_chrnum != hit3->other_chrnum */) {
	  debug(printf(" => double splice translocations"));
		
	} else if (SENSE_INCONSISTENT_FOR_INVERSION_P(hit3->sensedir_for_concordance,hit5->sensedir_for_concordance)) {
	  /* Use sensedir_for_concordance here and not sensedir */
	  debug(printf(" => sense inconsistent for inversion"));
#if 0
	} else if (hits5[j]->genomicend + querylength5 <= insert_start) {
	  debug(printf(" => scramble"));
	  if ((*nsamechr) <= maxpairedpaths &&
	      (stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/PAIRED_SCRAMBLE,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/false,/*transcriptome_guided_p*/false)) != NULL) {
	    *samechr = Hitlist_push(*samechr,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
#endif
	} else {
	  debug(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
		       hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
	  if ((stage3pair = Stage3pair_new(hit5,hit3,genestrand,/*pairtype*/PAIRED_INVERSION,
#if 0
					   queryuc_ptr_5,queryuc_ptr_3,
					   query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					   pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					   listpool,/*expect_concordant_p*/false,/*transcriptome_guided_p*/false)) != NULL) {
	    *samechr = Hitlist_push(*samechr,hitlistpool,(void *) stage3pair);
	    nfound++;
	  }
	}
	debug(printf("\n"));
	      
	j++;
      }
      j--;		/* Finish advance */
	    
      i++;
    }
  }

  return nfound;
}



List_T
Concordance_pair_up_distant (bool *abort_pairing_p, int *concordant_score,
			     List_T *samechr, List_T *conc_transloc, List_T hitpairs, 

			     List_T hitlist5_gplus, List_T hitlist5_gminus,
			     List_T hitlist3_gplus, List_T hitlist3_gminus,

			     Ladder_T ladder5_plus, Ladder_T ladder5_minus,
			     Ladder_T ladder3_plus, Ladder_T ladder3_minus,
			     int querylength5, int querylength3,
#if 0
			     char *queryuc_ptr_5, char *queryuc_ptr_3,
			     Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			     Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			     Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool,
#endif
			     Listpool_T listpool, Hitlistpool_T hitlistpool,
			     int maxpairedpaths, int genestrand) {
  int nfound = 0, n;
  int max_frontier_score, frontier_score, score5, score3;
  int cutoff_level, level;
  Ladder_T newladder5_plus, newladder5_minus, newladder3_plus, newladder3_minus;

  T *hits5, *hits3;
  int nhits5, nhits3;


  debug(printf("Starting Concordance_pair_up_distant\n"));

#if 0
  /* Rely instead on Ladder_cutoff, in case there are a few good hits */
  if (0 && List_length(hitlist5_gplus) > MAX_HITS) {
    Stage3end_gc(hitlist5_gplus);
    newladder5_plus = Ladder_new((List_T) NULL,hitlistpool,/*end5p*/true);
  } else {
    newladder5_plus = Ladder_new(hitlist5_gplus,hitlistpool,/*end5p*/true);
    Ladder_merge(ladder5_plus,newladder5_plus);
  }

  if (0 && List_length(hitlist5_gminus) > MAX_HITS) {
    Stage3end_gc(hitlist5_gminus);
    newladder5_minus = Ladder_new((List_T) NULL,hitlistpool,/*end5p*/true);
  } else {
    newladder5_minus = Ladder_new(hitlist5_gminus,hitlistpool,/*end5p*/true);
    Ladder_merge(ladder5_minus,newladder5_minus);
  }

  if (0 && List_length(hitlist3_gplus) > MAX_HITS) {
    Stage3end_gc(hitlist3_gplus);
    newladder3_plus = Ladder_new((List_T) NULL,hitlistpool,/*end5p*/false);
  } else {
    newladder3_plus = Ladder_new(hitlist3_gplus,hitlistpool,/*end5p*/false);
    Ladder_merge(ladder3_plus,newladder3_plus);
  }

  if (0 && List_length(hitlist3_gminus) > MAX_HITS) {
    Stage3end_gc(hitlist3_gminus);
    newladder3_minus = Ladder_new((List_T) NULL,hitlistpool,/*end5p*/false);
  } else {
    newladder3_minus = Ladder_new(hitlist3_gminus,hitlistpool,/*end5p*/false);
    Ladder_merge(ladder3_minus,newladder3_minus);
  }

#else
  newladder5_plus = Ladder_new(hitlist5_gplus,hitlistpool,/*end5p*/true);
  newladder5_minus = Ladder_new(hitlist5_gminus,hitlistpool,/*end5p*/true);
  newladder3_plus = Ladder_new(hitlist3_gplus,hitlistpool,/*end5p*/false);
  newladder3_minus = Ladder_new(hitlist3_gminus,hitlistpool,/*end5p*/false);

  Ladder_merge(ladder5_plus,newladder5_plus,hitlistpool);
  Ladder_merge(ladder5_minus,newladder5_minus,hitlistpool);
  Ladder_merge(ladder3_plus,newladder3_plus,hitlistpool);
  Ladder_merge(ladder3_minus,newladder3_minus,hitlistpool);
#endif


  cutoff_level = Ladder_cutoff(newladder5_plus) + Ladder_cutoff(ladder3_plus);
  if ((level = Ladder_cutoff(newladder5_minus) + Ladder_cutoff(ladder3_minus)) > cutoff_level) {
    cutoff_level = level;
  }
  if ((level = Ladder_cutoff(ladder5_plus) + Ladder_cutoff(newladder3_plus)) > cutoff_level) {
    cutoff_level = level;
  }
  if ((level = Ladder_cutoff(ladder5_minus) + Ladder_cutoff(newladder3_minus)) > cutoff_level) {
    cutoff_level = level;
  }
  if ((level = Ladder_cutoff(ladder5_plus) + Ladder_cutoff(ladder3_minus)) > cutoff_level) {
    cutoff_level = level;
  }
  if ((level = Ladder_cutoff(ladder5_minus) + Ladder_cutoff(ladder3_plus)) > cutoff_level) {
    cutoff_level = level;
  }
  max_frontier_score = cutoff_level;

  frontier_score = 0;
  while (frontier_score <= cutoff_level) {
    debug(printf("frontier_score = %d\n",frontier_score));
    for (score5 = 0; score5 <= frontier_score; score5++) {
      score3 = frontier_score - score5;
      debug(printf("score5 = %d, score3 = %d\n",score5,score3));

      /* New hits 5 vs all hits 3 */
      /* plus/plus: hits5_plus against hits3_plus (really on minus) */
      if (score5 <= Ladder_cutoff(newladder5_plus) && score3 <= Ladder_cutoff(ladder3_plus)) {
	hits5 = Ladder_hits_for_score(&nhits5,newladder5_plus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,ladder3_plus,hitlistpool,score3);
	debug(printf("at score %d+%d, nnewhits5_plus = %d, nhits3_plus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_distant_plus_plus(&(*concordant_score),&hitpairs,&(*conc_transloc),&(*samechr),
				      hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
				      queryuc_ptr_5,queryuc_ptr_3,
				      query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				      pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
				      listpool,hitlistpool)) > 0) {
	  if (nfound == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
	  if ((nfound += n) > maxpairedpaths) {
	    debug(printf(" -- %d found paths exceeds %d",nfound,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }

      /* minus/minus: hits3_minus (really on plus) against hits5_minus */
      if (score5 <= Ladder_cutoff(newladder5_minus) && score3 <= Ladder_cutoff(ladder3_minus)) {
	hits5 = Ladder_hits_for_score(&nhits5,newladder5_minus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,ladder3_minus,hitlistpool,score3);
	debug(printf("at score %d+%d, nnewhits5_minus = %d, nhits3_minus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_distant_minus_minus(&(*concordant_score),&hitpairs,&(*conc_transloc),&(*samechr),
					hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
					queryuc_ptr_5,queryuc_ptr_3,
					query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					listpool,hitlistpool)) > 0) {
	  if (nfound == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
	  if ((nfound += n) > maxpairedpaths) {
	    debug(printf(" -- %d found paths exceeds %d",nfound,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }
	

      /* All hits 5 vs new hits 3 */
      /* plus/plus: hits5_plus against hits3_plus (really on minus) */
      if (score5 <= Ladder_cutoff(ladder5_plus) && score3 <= Ladder_cutoff(newladder3_plus)) {
	hits5 = Ladder_hits_for_score(&nhits5,ladder5_plus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,newladder3_plus,hitlistpool,score3);
	debug(printf("at score %d+%d, nnewhits5_plus = %d, nhits3_plus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_distant_plus_plus(&(*concordant_score),&hitpairs,&(*conc_transloc),&(*samechr),
				      hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
				      queryuc_ptr_5,queryuc_ptr_3,
				      query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				      pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
				      listpool,hitlistpool)) > 0) {

	  if (nfound == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
	  if ((nfound += n) > maxpairedpaths) {
	    debug(printf(" -- %d found paths exceeds %d",nfound,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }

      /* minus/minus: hits3_minus (really on plus) against hits5_minus */
      if (score5 <= Ladder_cutoff(ladder5_minus) && score3 <= Ladder_cutoff(newladder3_minus)) {
	hits5 = Ladder_hits_for_score(&nhits5,ladder5_minus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,newladder3_minus,hitlistpool,score3);
	debug(printf("at score %d+%d, nnewhits5_minus = %d, nhits3_minus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_distant_minus_minus(&(*concordant_score),&hitpairs,&(*conc_transloc),&(*samechr),
					hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
					queryuc_ptr_5,queryuc_ptr_3,
					query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
					pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
					listpool,hitlistpool)) > 0) {
	  if (nfound == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
	  if ((nfound += n) > maxpairedpaths) {
	    debug(printf(" -- %d found paths exceeds %d",nfound,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }


      /* All hits5 vs all hits3 */
      /* plus/minus (inversions): hits5_plus against hits3_minus */
      if (score5 <= Ladder_cutoff(ladder5_plus) && score3 <= Ladder_cutoff(ladder3_minus)) {
	hits5 = Ladder_hits_for_score(&nhits5,ladder5_plus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,ladder3_minus,hitlistpool,score3);
	debug(printf("at score %d+%d, nhits5_plus = %d, nhits3_minus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_distant_plus_minus(&(*samechr),hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
				       queryuc_ptr_5,queryuc_ptr_3,
				       query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				       pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
				       listpool,hitlistpool)) > 0) {
	  if (nfound == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
	  if ((nfound += n) > maxpairedpaths) {
	    debug(printf(" -- %d found paths exceeds %d",nfound,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }

      /* minus/plus (inversions): hits3_plus against hits5_minus */
      if (score5 <= Ladder_cutoff(ladder5_minus) && score3 <= Ladder_cutoff(ladder3_plus)) {
	hits5 = Ladder_hits_for_score(&nhits5,ladder5_minus,hitlistpool,score5);
	hits3 = Ladder_hits_for_score(&nhits3,ladder3_plus,hitlistpool,score3);
	debug(printf("at score %d+%d, nhits5_minus = %d, nhits3_plus = %d\n",score5,score3,nhits5,nhits3));
	if ((n = do_distant_minus_plus(&(*samechr),hits5,nhits5,hits3,nhits3,genestrand,querylength5,querylength3,
#if 0
				       queryuc_ptr_5,queryuc_ptr_3,
				       query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
				       pairpool,dynprogL,dynprogM,dynprogR,oligoindices_minor,diagpool,cellpool,
#endif
				       listpool,hitlistpool)) > 0) {
	  if (nfound == 0) {
	    if ((cutoff_level = frontier_score + subopt_levels) > max_frontier_score) {
	      cutoff_level = max_frontier_score;
	    }
	  }
	  if ((nfound += n) > maxpairedpaths) {
	    debug(printf(" -- %d found paths exceeds %d",nfound,maxpairedpaths));
	    *abort_pairing_p = true;
	  }
	}
      }
    }

    frontier_score++;
  }


  Ladder_free(&newladder5_plus);
  Ladder_free(&newladder5_minus);
  Ladder_free(&newladder3_plus);
  Ladder_free(&newladder3_minus);

  Ladder_gc_duplicates(ladder5_plus);
  Ladder_gc_duplicates(ladder5_minus);
  Ladder_gc_duplicates(ladder3_plus);
  Ladder_gc_duplicates(ladder3_minus);

  debug(printf("Finished with Concordance_pair_up_distant: %d concordant, %d samechr, %d conc_transloc\n",
		List_length(hitpairs),List_length(*samechr),List_length(*conc_transloc)));

  return hitpairs;
}



#if 0
static int
path_ascending_cmp (const void *a, const void *b) {
  List_T x = * (List_T *) a;
  List_T y = * (List_T *) b;
  Univdiag_T diagonalx, diagonaly;
  Univcoord_T path_start_x, path_start_y, path_end_x, path_end_y;

  diagonalx = List_head(x);
  diagonaly = List_head(y);
  path_start_x = diagonalx->univdiagonal;
  path_start_y = diagonaly->univdiagonal;

  if (path_start_x < path_start_y) {
    return -1;
  } else if (path_start_y < path_start_x) {
    return +1;
  } else {
    diagonalx = (Univdiag_T) List_last_value(x);
    diagonaly = (Univdiag_T) List_last_value(y);
    path_end_x = diagonalx->univdiagonal;
    path_end_y = diagonaly->univdiagonal;
    if (path_end_x < path_end_y) {
      return -1;
    } else if (path_end_y < path_end_x) {
      return +1;
    } else {
      return 0;
    }
  }
}
#endif


#if 0
static int
hit_ascending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->high < y->high) {
    return -1;
  } else if (y->high < x->high) {
    return +1;
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else {
    return 0;
  }
}
#endif


void
Concordance_setup (int subopt_levels_in, bool novelsplicingp_in,
		   Chrpos_T pairmax_transcriptome_in,
		   Chrpos_T pairmax_linear_in, Chrpos_T pairmax_circular_in,
		   Chrpos_T expected_pairlength_in, Chrpos_T pairlength_deviation_in, 
		   bool *circularp_in, bool merge_samechr_p_in) {
  subopt_levels = subopt_levels_in;
  novelsplicingp = novelsplicingp_in;

  pairmax_transcriptome = pairmax_transcriptome_in;
  pairmax_linear = pairmax_linear_in;
  pairmax_circular = pairmax_circular_in;

  expected_pairlength = expected_pairlength_in;
  pairlength_deviation = pairlength_deviation_in;
  adjacent_pairlength = expected_pairlength + pairlength_deviation;

  circularp = circularp_in;
  merge_samechr_p = merge_samechr_p_in;

  return;
}
