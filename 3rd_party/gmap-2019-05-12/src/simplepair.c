static char rcsid[] = "$Id: simplepair.c 218286 2019-01-23 16:46:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "simplepair.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include <math.h>		/* For rint(), abs() */
#include <ctype.h>		/* For toupper */

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "comp.h"
#include "complement.h"
#include "transcript.h"
#include "method.h"



#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* hardclip_pairarray */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif


static bool novelsplicingp;
static IIT_T splicesites_iit;
static Univ_IIT_T transcript_iit;

static bool sam_insert_0M_p = false;
static bool md_lowercase_variant_p;
static bool snps_p;

static bool cigar_extended_p;
static Cigar_action_T cigar_action;



#define T Simplepair_T


Chrpos_T
Simplepair_head_genomepos (List_T pairs) {
  return ((T) List_head(pairs))->genomepos;
}

Chrpos_T
Simplepair_last_genomepos (List_T pairs) {
  return ((T) List_last_value(pairs))->genomepos;
}


/* For GSNAP, generate for output thread only.  For GMAP, pairs needed by worker threads are made in pairpool.c */
T
Simplepair_new_out (int querypos, Chrpos_T genomepos, char cdna, char comp, char genome) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  
  new->querypos = querypos;
  new->genomepos = genomepos;
  new->cdna = cdna;
  new->comp = comp;
  new->genome = genome;
  new->genomealt = genome;

  switch (comp) {
  case FWD_CANONICAL_INTRON_COMP /*'>'*/:
  case REV_CANONICAL_INTRON_COMP /*'<'*/:
  case NONINTRON_COMP /*'='*/:
  case SHORTGAP_COMP /*'~'*/:
  case INTRONGAP_COMP /*'.'*/:
  case FWD_GCAG_INTRON_COMP /*')'*/:
  case REV_GCAG_INTRON_COMP /*'('*/:
  case FWD_ATAC_INTRON_COMP /*']'*/:
  case REV_ATAC_INTRON_COMP /*'['*/:
  case DUALBREAK_COMP /*'#'*/:
    new->gapp = true; break;
  default: new->gapp = false;
  }

  return new;
}

void
Simplepair_free_out (T *old) {
  if (*old) {
    FREE_OUT(*old);
  }
  return;
}


static void
dump_one (T this, bool zerobasedp) {

  debug1(printf("%p ",this));

  if (this->gapp == true) {
    printf("*** Gap: type: ");
    switch (this->comp) {
    case FWD_CANONICAL_INTRON_COMP: printf("> GT-AG"); break;
    case FWD_GCAG_INTRON_COMP: printf(") GC-AG"); break;
    case FWD_ATAC_INTRON_COMP: printf("] AT-AC"); break;
    case REV_ATAC_INTRON_COMP: printf("[ AT-AC"); break;
    case REV_GCAG_INTRON_COMP: printf("( GC-AG"); break;
    case REV_CANONICAL_INTRON_COMP: printf("< GT-AG"); break;
    case SHORTGAP_COMP: printf("~ shortgap"); break;
    case NONINTRON_COMP: printf("= nonintron"); break;
    default: printf("? unknown"); break;
    }

    printf(" ***");

  } else {
    printf("%d %d %c ",
	   this->querypos + !zerobasedp,this->genomepos + !zerobasedp,this->cdna);

    /* Subtract 1 because dynprogindices start at +1 and -1 */
    putchar(this->comp);
    printf(" %c",this->genome);
    if (this->genomealt != this->genome) {
      printf(" alt:%c",this->genomealt);
    }
  }

  return;
}


/* Useful for debugging */
void
Simplepair_dump_list (List_T pairs, bool zerobasedp) {
  T this, prev = NULL, old = NULL;
  List_T p;

  printf("***Start of list***\n");
  for (p = pairs; p != NULL; p = List_next(p)) {
    this = List_head(p);
    dump_one(this,zerobasedp);
    printf("\n");

    if (this->querypos != -1) {
      if (old != NULL) {
	if (old->querypos > prev->querypos) {
	  if (prev->querypos < this->querypos) {
	    fprintf(stderr,"%d %d %d\n",old->querypos,prev->querypos,this->querypos);
	    abort();
	  }
	} else if (old->querypos < prev->querypos) {
	  if (prev->querypos > this->querypos) {
	    fprintf(stderr,"%d %d %d\n",old->querypos,prev->querypos,this->querypos);
	    abort();
	  }
	}
      }

      old = prev;
      prev = this;
    }

  }
  printf("***End of list***\n");
  return;
}  


List_T
Simplepair_strip_gaps_at_head (List_T pairs) {
  T pair;

  while (pairs != NULL && (pair = pairs->first) != NULL &&
	 (pair->gapp == true || pair->cdna == ' ' || pair->genome == ' ')) {
    pairs = List_pop_out(pairs,(void **) &pair);
    Simplepair_free_out(&pair);
  }

  return pairs;
}

List_T
Simplepair_strip_gaps_at_tail (List_T pairs) {
  T pair;

  if (pairs != NULL) {
    pairs = List_reverse(pairs);

    while (pairs != NULL && (pair = pairs->first) != NULL &&
	   (pair->gapp == true || pair->cdna == ' ' || pair->genome == ' ')) {
      pairs = List_pop_out(pairs,(void **) &pair);
      Simplepair_free_out(&pair);
    }

    pairs = List_reverse(pairs);
  }

  return pairs;
}



static struct T *
hardclip_pairarray (int *clipped_npairs, int hardclip_start, int hardclip_end,
		    struct T *pairs, int npairs, int querylength) {
  struct T *clipped_pairs, *ptr;
  int i, starti;

  debug10(printf("Entered hardclip_pairarray with hardclip_start %d, hardclip_end %d, querylength %d\n",
		 hardclip_start,hardclip_end,querylength));
  debug10(Simplepair_dump_array(pairs,npairs,true));
  debug10(printf("Starting with %d pairs\n",npairs));

  i = 0;
  ptr = pairs;
  while (i < npairs && ptr->querypos < hardclip_start) {
    i++;
    ptr++;
  }
  while (i < npairs && (ptr->gapp == true || ptr->cdna == ' ' || ptr->genome == ' ')) {
    i++;
    ptr++;
  }

  if (i >= npairs) {
    /* hardclip_start passes right end of read, so invalid */
    debug10(printf("i = %d, so passed end of read\n",i));
    hardclip_start = 0;
  } else if (hardclip_start > 0) {
    hardclip_start = ptr->querypos;
  }

  starti = i;
  debug10(printf("starti is %d\n",starti));

  clipped_pairs = ptr;

  while (i < npairs && ptr->querypos < querylength - hardclip_end) {
    i++;
    ptr++;
  }

  i--;
  ptr--;
  while (i >= starti && (ptr->gapp == true || ptr->cdna == ' ' || ptr->genome == ' ')) {
    i--;
    ptr--;
  }
  
  if (i < 0) {
    /* hardclip_end passes left end of read, so invalid */
    debug10(printf("i = %d, so passed left end of read\n",i));
    hardclip_end = 0;
  } else if (hardclip_end > 0) {
    hardclip_end = querylength - 1 - ptr->querypos;
  }

  if (hardclip_start == 0 && hardclip_end == 0) {
    debug10(printf("Unable to hard clip\n"));
    *clipped_npairs = npairs;
    clipped_pairs = pairs;
  } else {
    *clipped_npairs = i - starti + 1;
  }

  debug10(printf("Ending with %d pairs\n",*clipped_npairs));
  debug10(printf("Exiting hardclip_pairs with hardclip_start %d, hardclip_end %d\n",
		 hardclip_start,hardclip_end));

  return clipped_pairs;
}


Chrpos_T
Simplepair_genomicpos_low (int hardclip_low, int hardclip_high,
			   struct T *pairarray, int npairs, int querylength,
			   bool plusp, bool hide_soft_clips_p) {
  struct T *clipped_pairs;
  int clipped_npairs;
  T lowpair, highpair;

#if 0
  if (clipdir >= 0) {
    if (plusp == true) {
      if (first_read_p == true) {
	hardclip_high = hardclip5;
	hardclip_low = 0;
      } else {
	hardclip_high = 0;
	hardclip_low = hardclip3;
      }
    } else {
      if (first_read_p == true) {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
    }
  } else {
    if (plusp == true) {
      if (first_read_p == true) {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
    } else {
      if (first_read_p == true) {
	hardclip_high = hardclip5;
	hardclip_low = 0;
      } else {
	hardclip_high = 0;
	hardclip_low = hardclip3;
      }
    }
  }
#endif

  if (plusp == true) {
    clipped_pairs = hardclip_pairarray(&clipped_npairs,hardclip_low,hardclip_high,
				       pairarray,npairs,querylength);
    lowpair = &(clipped_pairs[0]);
    highpair = &(clipped_pairs[clipped_npairs-1]);
    if (hide_soft_clips_p == true) {
      assert(lowpair->querypos == 0);
      assert(highpair->querypos == querylength - 1);
      /* *chrpos_high = highpair->genomepos + 1 - (querylength - 1 - highpair->querypos); */
      return lowpair->genomepos + 1 - lowpair->querypos;
    } else {
      /* *chrpos_high = highpair->genomepos + 1; */
      return lowpair->genomepos + 1;
    }
  } else {
    /* Swap hardclip_low and hardclip_high */
    clipped_pairs = hardclip_pairarray(&clipped_npairs,hardclip_high,hardclip_low,
				       pairarray,npairs,querylength);
    lowpair = &(clipped_pairs[clipped_npairs-1]);
    highpair = &(clipped_pairs[0]);
    if (hide_soft_clips_p == true) {
      assert(lowpair->querypos == querylength - 1);
      assert(highpair->querypos == 0);
      /* *chrpos_high = highpair->genomepos + 1 - highpair->querypos; */
      return lowpair->genomepos + 1 + (querylength - 1 - lowpair->querypos);
    } else {
      /* *chrpos_high = highpair->genomepos + 1; */
      return lowpair->genomepos + 1;
    }
  }
}


static List_T
push_token (List_T tokens, char *token) {
  char *copy;

  copy = (char *) MALLOC_OUT((strlen(token)+1) * sizeof(char));
  strcpy(copy,token);
  return List_push_out(tokens,(void *) copy);
}

static void
tokens_free (List_T *tokens) {
  List_T p;
  char *token;

  for (p = *tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FREE_OUT(token);
  }
  List_free_out(&(*tokens));

  return;
}

static void
print_tokens (Filestring_T fp, List_T tokens) {
  List_T p;
  char *token;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FPRINTF(fp,"%s",token);
    /* FREE_OUT(token); -- Now freed within Stage3end_free or Stage3_free */
  }

  return;
}

static int
cigar_length (List_T tokens) {
  int length = 0, tokenlength;
  List_T p;
  char *token;
  char type;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    type = token[strlen(token)-1];
    /* Should include 'H', but that gets added according to hardclip_low and hardclip_high */
    if (type == 'S' || type == 'I' || type == 'M' || type == 'X' || type == '=') {
      sscanf(token,"%d",&tokenlength);
      length += tokenlength;
    }
  }

  return length;
}


static List_T
clean_cigar (List_T tokens, bool watsonp) {
  List_T clean, unique = NULL, p;
  char token[11], *curr_token, *last_token;
  int length = 0;
  char type, last_type = ' ';
  bool duplicatep = false;

  for (p = tokens; p != NULL; p = List_next(p)) {
    curr_token = (char *) List_head(p);
    type = curr_token[strlen(curr_token)-1];
    if (type == last_type) {
      length += atoi(last_token);
      FREE_OUT(last_token);
      duplicatep = true;
    } else {
      if (last_type == ' ') {
	/* Skip */
      } else if (duplicatep == false) {
	unique = List_push_out(unique,(void *) last_token);
      } else {
	length += atoi(last_token);
	FREE_OUT(last_token);
	sprintf(token,"%d%c",length,last_type);
	unique = push_token(unique,token);
      }
      last_type = type;
      duplicatep = false;
      length = 0;
    }
    last_token = curr_token;
  }
  if (last_type == ' ') {
    /* Skip */
  } else if (duplicatep == false) {
    unique = List_push_out(unique,(void *) last_token);
  } else {
    length += atoi(last_token);
    FREE_OUT(last_token);
    sprintf(token,"%d%c",length,last_type);
    unique = push_token(unique,token);
  }
  List_free_out(&tokens);


  if (sam_insert_0M_p == false) {
    /* Return result */
    if (watsonp) {
      /* Put tokens in forward order */
      return unique;
    } else {
      /* Keep tokens in reverse order */
      return List_reverse(unique);
    }

  } else {
    /* Insert "0M" between adjacent I and D operations */
    last_type = ' ';
    clean = (List_T) NULL;
    for (p = unique; p != NULL; p = List_next(p)) {
      curr_token = (char *) List_head(p);
      type = curr_token[strlen(curr_token)-1];
      if (last_type == 'I' && type == 'D') {
	clean = push_token(clean,"0M");
      } else if (last_type == 'D' && type == 'I') {
	clean = push_token(clean,"0M");
      }
      clean = List_push_out(clean,(void *) curr_token);
      last_type = type;
    }
    List_free_out(&unique);

    /* Return result */
    if (watsonp) {
      /* Put tokens in forward order */
      return List_reverse(clean);
    } else {
      /* Keep tokens in reverse order */
      return clean;
    }
  }
}


static List_T
compute_cigar_standard (bool *intronp, int *hardclip_start, int *hardclip_end, struct T *pairs, int npairs, int querylength_given,
			bool watsonp,
#ifdef CONVERT_INTRONS_TO_DELETIONS
			int sensedir,
#endif
			int chimera_part) {
  List_T tokens = NULL;
  char token[11];
  int Mlength = 0, Ilength = 0, Dlength = 0;
  bool in_exon = false, deletionp;
  struct T *ptr, *prev, *this = NULL;
  int exon_queryend = -1;
  Chrpos_T exon_genomestart = 0;
  Chrpos_T exon_genomeend, genome_gap;
  int query_gap;
  int last_querypos = -1;
  Chrpos_T last_genomepos = (Chrpos_T) -1;
  int i;

  /* *chimera_hardclip_start = *chimera_hardclip_high = 0; */
  *intronp = false;

  ptr = pairs;

  if (chimera_part == +1) {
    if (ptr->querypos > *hardclip_start) {
      if (ptr->querypos > 0) {
	/* Clip to beginning */
	*hardclip_start = ptr->querypos;
	sprintf(token,"%dH",*hardclip_start);
	tokens = push_token(tokens,token);
      }
    } else {
      if (*hardclip_start > 0) {
	/* Clip to hard clip boundary */
	sprintf(token,"%dH",*hardclip_start);
	tokens = push_token(tokens,token);
      }
    }
  } else {
    if (*hardclip_start > 0) {
      sprintf(token,"%dH",*hardclip_start);
      tokens = push_token(tokens,token);
    }
    if (ptr->querypos > (*hardclip_start)) {
      sprintf(token,"%dS",ptr->querypos - (*hardclip_start));
      tokens = push_token(tokens,token);
    }
  }

  this = (T) NULL;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

#if 0
    /* Cigar_print_tokens(stdout,tokens); */
    Pair_dump_one(this,true);
    printf("\n");
#endif

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
#if 0
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
#endif
	
	if (Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(tokens,token);
	} else if (Ilength > 0) {
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(tokens,token);
	} else if (Dlength > 0) {
	  sprintf(token,"%dD",Dlength);
	  tokens = push_token(tokens,token);
	}

	Mlength = Ilength = Dlength = 0;

	in_exon = false;
      }

    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* exon_querystart = this->querypos + 1; */
	exon_genomestart = this->genomepos + 1;

	if (prev != NULL) {
	  /* Gap */
	  /* abs() gives a large value when flag -m64 is specified */
	  /* genome_gap = abs(intron_end - intron_start) + 1; */
	  if (watsonp) {
	    /* intron_end = exon_genomestart - 1; */
	    /* genome_gap = (intron_end - intron_start) + 1; */
	    genome_gap = exon_genomestart - exon_genomeend - 1;
	  } else {
	    /* intron_end = exon_genomestart + 1; */
	    /* genome_gap = (intron_start - intron_end) + 1; */
	    genome_gap = exon_genomeend - exon_genomestart - 1;
	  }

	  deletionp = false;
#ifdef CONVERT_INTRONS_TO_DELETIONS
	  if (sensedir == SENSE_FORWARD) {
	    if (prev->comp == FWD_CANONICAL_INTRON_COMP ||
		prev->comp == FWD_GCAG_INTRON_COMP ||
		prev->comp == FWD_ATAC_INTRON_COMP) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else {
	      sprintf(token,"%uD",genome_gap);
	      deletionp = true;
	    }
	  } else if (sensedir == SENSE_ANTI) {
	    if (prev->comp == REV_CANONICAL_INTRON_COMP ||
		prev->comp == REV_GCAG_INTRON_COMP ||
		prev->comp == REV_ATAC_INTRON_COMP) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else {
	      sprintf(token,"%uD",genome_gap);
	      deletionp = true;
	    }
	  } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN){
	    sprintf(token,"%uN",genome_gap);
	    *intronp = true;
	  } else {
	    sprintf(token,"%uD",genome_gap);
	    deletionp = true;
	  }
#else
	  sprintf(token,"%uN",genome_gap);
	  *intronp = true;
#endif
	  tokens = push_token(tokens,token);

	  /* Check for dual gap.  Doesn't work for hard clipping. */
	  /* assert(exon_queryend >= 0); */

	  query_gap = this->querypos - exon_queryend;
	  assert(query_gap >= 0);
	  if (query_gap > 0) {
	    if (deletionp == true && sam_insert_0M_p == true) {
	      /* Put zero matches between deletion and insertion, since some programs will complain */
	      sprintf(token,"0M");
	      tokens = push_token(tokens,token);
	    }

	    sprintf(token,"%uI",query_gap);
	    tokens = push_token(tokens,token);
	  }
	}

	in_exon = true;
      }

      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (this->genome == ' ') {
	  /* Insertion relative to genome */
	  if (Mlength > 0) {
	    sprintf(token,"%dM",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Dlength > 0) {
	    /* unlikely */
	    sprintf(token,"%dD",Dlength);
	    tokens = push_token(tokens,token);
	    Dlength = 0;
	  }
	  Ilength++;
	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (Mlength > 0) {
	    sprintf(token,"%dM",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Ilength > 0) {
	    sprintf(token,"%dI",Ilength);
	    tokens = push_token(tokens,token);
	    Ilength = 0;
	  }
	  Dlength++;
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	/* Count even if unknown base */

	if (Ilength > 0) {
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(tokens,token);
	  Ilength = 0;
	} else if (Dlength > 0) {
	  sprintf(token,"%dD",Dlength);
	  tokens = push_token(tokens,token);
	  Dlength = 0;
	}
	Mlength++;

      }
    }

    if (this != NULL) {
      if (this->cdna != ' ') {
	last_querypos = this->querypos;
      }
      if (this->genome != ' ') {
	last_genomepos = this->genomepos;
      }
    }
  }

  /* prev = this; */
  /* exon_queryend = last_querypos + 1; */
  /* exon_genomeend = last_genomepos + 1; */

  if (Mlength > 0) {
    sprintf(token,"%dM",Mlength);
    tokens = push_token(tokens,token);
  } else if (Ilength > 0) {
    sprintf(token,"%dI",Ilength);
    tokens = push_token(tokens,token);
  } else if (Dlength > 0) {
    sprintf(token,"%dD",Dlength);
    tokens = push_token(tokens,token);
  }


  /* Terminal clipping */
  if (chimera_part == -1) {
    if (last_querypos < querylength_given - 1 - (*hardclip_end)) {
      if (last_querypos < querylength_given - 1) {
	/* Clip to end */
	*hardclip_end = querylength_given - 1 - last_querypos;
	sprintf(token,"%dH",*hardclip_end);
	tokens = push_token(tokens,token);
      }
    } else {
      if (*hardclip_end > 0) {
	/* Clip to hard clip boundary */
	sprintf(token,"%dH",*hardclip_end);
	tokens = push_token(tokens,token);
      }
    }
  } else {
    if (last_querypos < querylength_given - 1 - (*hardclip_end)) {
      sprintf(token,"%dS",querylength_given - 1 - (*hardclip_end) - last_querypos);
      tokens = push_token(tokens,token);
    }
    if (*hardclip_end > 0) {
      sprintf(token,"%dH",*hardclip_end);
      tokens = push_token(tokens,token);
    }
  }

  return clean_cigar(tokens,watsonp);
}


static List_T
compute_cigar_extended (bool *intronp, int *hardclip_start, int *hardclip_end, struct T *pairs, int npairs, int querylength_given,
			bool watsonp,
#ifdef CONVERT_INTRONS_TO_DELETIONS
			int sensedir,
#endif
			int chimera_part) {
  List_T tokens = NULL;
  char token[11];
  int Elength = 0, Xlength = 0, Ilength = 0, Dlength = 0;
  bool in_exon = false, deletionp;
  struct T *ptr, *prev, *this = NULL;
  int exon_queryend = -1;
  Chrpos_T exon_genomestart = 0;
  Chrpos_T exon_genomeend, genome_gap;
  int query_gap;
  int last_querypos = -1;
  Chrpos_T last_genomepos = (Chrpos_T) -1;
  int i;

  /* *chimera_hardclip_start = *chimera_hardclip_high = 0; */
  *intronp = false;

  ptr = pairs;

  if (chimera_part == +1) {
    if (ptr->querypos > *hardclip_start) {
      if (ptr->querypos > 0) {
	/* Clip to beginning */
	*hardclip_start = ptr->querypos;
	sprintf(token,"%dH",*hardclip_start);
	tokens = push_token(tokens,token);
      }
    } else {
      if (*hardclip_start > 0) {
	/* Clip to hard clip boundary */
	sprintf(token,"%dH",*hardclip_start);
	tokens = push_token(tokens,token);
      }
    }
  } else {
    if (*hardclip_start > 0) {
      sprintf(token,"%dH",*hardclip_start);
      tokens = push_token(tokens,token);
    }
    if (ptr->querypos > (*hardclip_start)) {
      sprintf(token,"%dS",ptr->querypos - (*hardclip_start));
      tokens = push_token(tokens,token);
    }
  }

  this = (T) NULL;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

#if 0
    /* Cigar_print_tokens(stdout,tokens); */
    Pair_dump_one(this,true);
    printf("\n");
#endif

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
#if 0
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
#endif
	
	if (Elength > 0) {
	  sprintf(token,"%d=",Elength);
	  tokens = push_token(tokens,token);
	} else if (Xlength > 0) {
	  sprintf(token,"%dX",Xlength);
	  tokens = push_token(tokens,token);
	} else if (Ilength > 0) {
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(tokens,token);
	} else if (Dlength > 0) {
	  sprintf(token,"%dD",Dlength);
	  tokens = push_token(tokens,token);
	}

	Elength = Xlength = Ilength = Dlength = 0;

	in_exon = false;
      }

    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* exon_querystart = this->querypos + 1; */
	exon_genomestart = this->genomepos + 1;

	if (prev != NULL) {
	  /* Gap */
	  /* abs() gives a large value when flag -m64 is specified */
	  /* genome_gap = abs(intron_end - intron_start) + 1; */
	  if (watsonp) {
	    /* intron_end = exon_genomestart - 1; */
	    /* genome_gap = (intron_end - intron_start) + 1; */
	    genome_gap = exon_genomestart - exon_genomeend - 1;
	  } else {
	    /* intron_end = exon_genomestart + 1; */
	    /* genome_gap = (intron_start - intron_end) + 1; */
	    genome_gap = exon_genomeend - exon_genomestart - 1;
	  }

	  deletionp = false;
#ifdef CONVERT_INTRONS_TO_DELETIONS
	  if (sensedir == SENSE_FORWARD) {
	    if (prev->comp == FWD_CANONICAL_INTRON_COMP ||
		prev->comp == FWD_GCAG_INTRON_COMP ||
		prev->comp == FWD_ATAC_INTRON_COMP) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else {
	      sprintf(token,"%uD",genome_gap);
	      deletionp = true;
	    }
	  } else if (sensedir == SENSE_ANTI) {
	    if (prev->comp == REV_CANONICAL_INTRON_COMP ||
		prev->comp == REV_GCAG_INTRON_COMP ||
		prev->comp == REV_ATAC_INTRON_COMP) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else {
	      sprintf(token,"%uD",genome_gap);
	      deletionp = true;
	    }
	  } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN){
	    sprintf(token,"%uN",genome_gap);
	    *intronp = true;
	  } else {
	    sprintf(token,"%uD",genome_gap);
	    deletionp = true;
	  }
#else
	  sprintf(token,"%uN",genome_gap);
	  *intronp = true;
#endif
	  tokens = push_token(tokens,token);

	  /* Check for dual gap.  Doesn't work for hard clipping. */
	  /* assert(exon_queryend >= 0); */

	  query_gap = this->querypos - exon_queryend;
	  assert(query_gap >= 0);
	  if (query_gap > 0) {
	    if (deletionp == true && sam_insert_0M_p == true) {
	      /* Put zero matches between deletion and insertion, since some programs will complain */
	      sprintf(token,"0M");
	      tokens = push_token(tokens,token);
	    }

	    sprintf(token,"%uI",query_gap);
	    tokens = push_token(tokens,token);
	  }
	}

	in_exon = true;
      }

      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (this->genome == ' ') {
	  /* Insertion relative to genome */
	  if (Elength > 0) {
	    sprintf(token,"%d=",Elength);
	    tokens = push_token(tokens,token);
	    Elength = 0;
	  } else if (Xlength > 0) {
	    sprintf(token,"%dX",Xlength);
	    tokens = push_token(tokens,token);
	    Xlength = 0;
	  } else if (Dlength > 0) {
	    /* unlikely */
	    sprintf(token,"%dD",Dlength);
	    tokens = push_token(tokens,token);
	    Dlength = 0;
	  }
	  Ilength++;
	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (Elength > 0) {
	    sprintf(token,"%d=",Elength);
	    tokens = push_token(tokens,token);
	    Elength = 0;
	  } else if (Xlength > 0) {
	    sprintf(token,"%dX",Xlength);
	    tokens = push_token(tokens,token);
	    Xlength = 0;
	  } else if (Ilength > 0) {
	    sprintf(token,"%dI",Ilength);
	    tokens = push_token(tokens,token);
	    Ilength = 0;
	  }
	  Dlength++;
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	/* Count even if unknown base */

	if (Ilength > 0) {
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(tokens,token);
	  Ilength = 0;
	} else if (Dlength > 0) {
	  sprintf(token,"%dD",Dlength);
	  tokens = push_token(tokens,token);
	  Dlength = 0;
	}

	if (prev == NULL || prev->gapp || prev->comp == INDEL_COMP || prev->comp == SHORTGAP_COMP) {
	  if (this->cdna == this->genome) {
	    Elength++;
	  } else {
	    Xlength++;
	  }

	} else if (prev->cdna == prev->genome) {
	  if (this->cdna == this->genome) {
	    Elength++;
	  } else {
	    if (Elength > 0) {
	      sprintf(token,"%d=",Elength);
	      tokens = push_token(tokens,token);
	      Elength = 0;
	    }
	    Xlength++;
	  }

	} else {
	  if (this->cdna != this->genome) {
	    Xlength++;
	  } else {
	    if (Xlength > 0) {
	      sprintf(token,"%dX",Xlength);
	      tokens = push_token(tokens,token);
	      Xlength = 0;
	    }
	    Elength++;
	  }
	}
      }
    }

    if (this != NULL) {
      if (this->cdna != ' ') {
	last_querypos = this->querypos;
      }
      if (this->genome != ' ') {
	last_genomepos = this->genomepos;
      }
    }
  }

  /* prev = this; */
  /* exon_queryend = last_querypos + 1; */
  /* exon_genomeend = last_genomepos + 1; */

  if (Elength > 0) {
    sprintf(token,"%d=",Elength);
    tokens = push_token(tokens,token);
  } else if (Xlength > 0) {
    sprintf(token,"%dX",Xlength);
    tokens = push_token(tokens,token);
  } else if (Ilength > 0) {
    sprintf(token,"%dI",Ilength);
    tokens = push_token(tokens,token);
  } else if (Dlength > 0) {
    sprintf(token,"%dD",Dlength);
    tokens = push_token(tokens,token);
  }


  /* Terminal clipping */
  if (chimera_part == -1) {
    if (last_querypos < querylength_given - 1 - (*hardclip_end)) {
      if (last_querypos < querylength_given - 1) {
	/* Clip to end */
	*hardclip_end = querylength_given - 1 - last_querypos;
	sprintf(token,"%dH",*hardclip_end);
	tokens = push_token(tokens,token);
      }
    } else {
      if (*hardclip_end > 0) {
	/* Clip to hard clip boundary */
	sprintf(token,"%dH",*hardclip_end);
	tokens = push_token(tokens,token);
      }
    }
  } else {
    if (last_querypos < querylength_given - 1 - (*hardclip_end)) {
      sprintf(token,"%dS",querylength_given - 1 - (*hardclip_end) - last_querypos);
      tokens = push_token(tokens,token);
    }
    if (*hardclip_end > 0) {
      sprintf(token,"%dH",*hardclip_end);
      tokens = push_token(tokens,token);
    }
  }

  return clean_cigar(tokens,watsonp);
}


static List_T
compute_cigar (bool *intronp, int *hardclip_start, int *hardclip_end, struct T *pairs, int npairs, int querylength_given,
	       bool watsonp, int chimera_part) {
  if (cigar_extended_p == true) {
    return compute_cigar_extended(&(*intronp),&(*hardclip_start),&(*hardclip_end),pairs,npairs,querylength_given,
				  watsonp,chimera_part);
  } else {
    return compute_cigar_standard(&(*intronp),&(*hardclip_start),&(*hardclip_end),pairs,npairs,querylength_given,
				  watsonp,chimera_part);
  }
}


typedef enum {IN_MATCHES, IN_MISMATCHES, IN_DELETION} MD_state_T;

static char complCode[128] = COMPLEMENT_LC;

static List_T
compute_md_string (int *nmismatches_refdiff, int *nmismatches_bothdiff, int *nindels,
		   struct T *pairs, int npairs, bool watsonp, List_T cigar_tokens) {
  List_T md_tokens = NULL, p;
  char *cigar_token, token[11], *first_token, type;
  T this;
  int nmatches = 0, length;
  MD_state_T state = IN_MISMATCHES;
  int i, k = 0;

  *nmismatches_refdiff = *nmismatches_bothdiff = *nindels = 0;

  debug4(Pair_dump_array(pairs,npairs,true));
  debug4(printf("watsonp %d\n",watsonp));

  if (watsonp == true) {
    for (p = cigar_tokens; p != NULL; p = List_next(p)) {
      cigar_token = (char *) List_head(p);
      debug4(printf("token is %s\n",cigar_token));
      type = cigar_token[strlen(cigar_token)-1];
      length = atoi(cigar_token);
    
      if (type == 'H') {
	/* k += length; */

      } else if (type == 'S') {
	/* k += length; */

      } else if (type == 'M' || type == 'X' || type == '=') {
	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  debug4(printf("M %d/%d comp %c\n",i,length,this->comp));
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	    nmatches++;
	    state = IN_MATCHES;

	  } else if (this->comp == MISMATCH_COMP) {
	    if (state == IN_MATCHES) {
	      sprintf(token,"%d",nmatches);
	      md_tokens = push_token(md_tokens,token);
	      nmatches = 0;
	    } else if (state == IN_DELETION) {
	      md_tokens = push_token(md_tokens,"0");
	    }
	    state = IN_MISMATCHES;

	    *nmismatches_refdiff += 1;
	    if (md_lowercase_variant_p && this->cdna == this->genomealt) {
	      /* A mismatch against the reference only => alternate variant */
	      sprintf(token,"%c",tolower(this->genome));
	    } else {
	      /* A true mismatch against both variants */
	      *nmismatches_bothdiff += 1;
	      sprintf(token,"%c",this->genome);
	    }
	    md_tokens = push_token(md_tokens,token);

	  } else {
	    fprintf(stderr,"Unexpected comp '%c'\n",this->comp);
	    abort();
	  }
	}

      } else if (type == 'I') {
	while (k < npairs && pairs[k].comp == INDEL_COMP && pairs[k].genome == ' ') {
	  *nindels += 1;
	  k++;
	}
	state = IN_MATCHES;

      } else if (type == 'N') {
	while (k < npairs && pairs[k].gapp == true) {
	  k++;
	}

      } else if (type == 'D') {
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    md_tokens = push_token(md_tokens,token);
	    nmatches = 0;
	  }
	}

	if (state != IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}
	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  sprintf(token,"%c",this->genome);
	  md_tokens = push_token(md_tokens,token);
	  *nindels += 1;
	}

	state = IN_DELETION;

      } else {
	fprintf(stderr,"Don't recognize type %c\n",type);
	abort();
      }
    }

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      md_tokens = push_token(md_tokens,token);
    }

    md_tokens = List_reverse(md_tokens);

  } else {
    cigar_tokens = List_reverse(cigar_tokens);
    for (p = cigar_tokens; p != NULL; p = List_next(p)) {
      cigar_token = (char *) List_head(p);
      debug4(printf("token is %s\n",cigar_token));
      type = cigar_token[strlen(cigar_token)-1];
      length = atoi(cigar_token);
    
      if (type == 'H') {
	/* k += length; */

      } else if (type == 'S') {
	/* k += length; */

      } else if (type == 'M' || type == 'X' || type == '=') {
	if (state == IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}

	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  debug4(printf("M %d/%d comp %c\n",i,length,this->comp));
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	    nmatches++;
	    state = IN_MATCHES;

	  } else if (this->comp == MISMATCH_COMP) {
	    if (state == IN_MATCHES) {
	      sprintf(token,"%d",nmatches);
	      md_tokens = push_token(md_tokens,token);
	      nmatches = 0;
	    }
	    state = IN_MISMATCHES;

	    *nmismatches_refdiff += 1;

	    if (md_lowercase_variant_p && this->cdna == this->genomealt) {
	      /* A mismatch against the reference only => alternate variant */
	      sprintf(token,"%c",tolower(complCode[(int) this->genome]));
	    } else {
	      *nmismatches_bothdiff += 1;
	      sprintf(token,"%c",complCode[(int) this->genome]);
	    }
	    md_tokens = push_token(md_tokens,token);


	  } else {
	    fprintf(stderr,"Unexpected comp '%c'\n",this->comp);
	    abort();
	  }
	}

      } else if (type == 'I') {
	if (state == IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}

	while (k < npairs && pairs[k].comp == INDEL_COMP && pairs[k].genome == ' ') {
	  *nindels += 1;
	  k++;
	}
	state = IN_MATCHES;

      } else if (type == 'N') {
#if 0
	/* Ignore deletion adjacent to intron, to avoid double ^^ */
	if (state == IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}
#endif

	while (k < npairs && pairs[k].gapp == true) {
	  k++;
	}

      } else if (type == 'D') {
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    md_tokens = push_token(md_tokens,token);
	    nmatches = 0;
	  }
	} else if (state == IN_MISMATCHES) {
	  md_tokens = push_token(md_tokens,"0");
	}

	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  sprintf(token,"%c",complCode[(int) this->genome]);
	  md_tokens = push_token(md_tokens,token);
	  *nindels += 1;
	}
	state = IN_DELETION;

      } else {
	fprintf(stderr,"Don't recognize type %c\n",type);
	abort();
      }
    }

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      md_tokens = push_token(md_tokens,token);
    }

    /* Restore cigar_tokens */
    cigar_tokens = List_reverse(cigar_tokens);
  }

  assert(k == npairs);

  /* Insert initial 0 token if necessary */
  if (md_tokens != NULL) {
    first_token = (char *) List_head(md_tokens);
    if (!isdigit(first_token[0])) {
      md_tokens = push_token(md_tokens,"0");
    }
  }

  return md_tokens;
}


/* Modeled after Shortread_print_chopped */
static void
print_chopped (Filestring_T fp, char *contents, int querylength,
	       int hardclip_start, int hardclip_end) {
  int i;

  for (i = hardclip_start; i < querylength - hardclip_end; i++) {
    PUTC(contents[i],fp);
  }
  return;
}

/* Differs from Shortread version, in that hardclip_high and hardclip_low are not reversed */
static void
print_chopped_revcomp (Filestring_T fp, char *contents, int querylength,
		       int hardclip_start, int hardclip_end) {
  int i;

  for (i = querylength - 1 - hardclip_end; i >= hardclip_start; --i) {
    PUTC(complCode[(int) contents[i]],fp);
  }
  return;
}


static void
print_chopped_end (Filestring_T fp, char *contents, int querylength,
		   int hardclip_start, int hardclip_end) {
  int i;

  for (i = 0; i < hardclip_start; i++) {
    PUTC(contents[i],fp);
  }

  /* No separator */

  for (i = querylength - hardclip_end; i < querylength; i++) {
    PUTC(contents[i],fp);
  }

  return;
}

/* Differs from Shortread version, in that hardclip_high and hardclip_low are not reversed */
static void
print_chopped_end_revcomp (Filestring_T fp, char *contents, int querylength,
			   int hardclip_start, int hardclip_end) {
  int i;

  for (i = querylength - 1; i >= querylength - hardclip_end; --i) {
    PUTC(complCode[(int) contents[i]],fp);
  }

  /* No separator */

  for (i = hardclip_start - 1; i >= 0; --i) {
    PUTC(complCode[(int) contents[i]],fp);
  }

  return;
}


static void
print_chopped_end_quality (Filestring_T fp, char *quality, int querylength,
			   int hardclip_start, int hardclip_end) {
  int i;

  if (hardclip_start > 0) {
    for (i = 0; i < hardclip_start; i++) {
      PUTC(quality[i],fp);
    }
    return;

  } else {
    for (i = querylength - hardclip_end; i < querylength; i++) {
      PUTC(quality[i],fp);
    }
    return;
  }
}

/* Differs from Shortread version, in that hardclip_high and hardclip_low are not reversed */
static void
print_chopped_end_quality_reverse (Filestring_T fp, char *quality, int querylength,
				   int hardclip_start, int hardclip_end) {
  int i;

  if (hardclip_start > 0) {
    for (i = hardclip_start - 1; i >= 0; --i) {
      PUTC(quality[i],fp);
    }
    return;

  } else {
    for (i = querylength - 1; i >= querylength - hardclip_end; --i) {
      PUTC(quality[i],fp);
    }
    return;
  }
}



/* Modeled after Shortread_print_quality */
static void
print_quality (Filestring_T fp, char *quality, int querylength,
	       int hardclip_start, int hardclip_end, int shift) {
  int i;
  int c;

  if (quality == NULL) {
    PUTC('*',fp);
  } else {
    for (i = hardclip_start; i < querylength - hardclip_end; i++) {
      if ((c = quality[i] + shift) <= 32) {
	fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		shift,quality[i]);
	abort();
      } else {
	PUTC(c,fp);
      }
    }
  }
  return;
}


static void
print_quality_revcomp (Filestring_T fp, char *quality, int querylength,
		       int hardclip_start, int hardclip_end, int shift) {
  int i;
  int c;

  if (quality == NULL) {
    PUTC('*',fp);
  } else {
    for (i = querylength - 1 - hardclip_end; i >= hardclip_start; --i) {
      if ((c = quality[i] + shift) <= 32) {
	fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		shift,quality[i]);
	abort();
      } else {
	PUTC(c,fp);
      }
    }
  }

  return;
}

/* Derived from print_gff3_cdna_match */
/* Assumes pairarray has been hard clipped already */
static void
print_overlap_sam_line (Filestring_T fp, char *abbrev, char *acc1, char *acc2, char *chrstring,
			bool watsonp, int sensedir, List_T cigar_tokens, List_T md_tokens,
			int nmismatches_refdiff, int nmismatches_bothdiff, int nindels,
			bool intronp, char *queryseq_ptr, char *quality_string,
			int hardclip_start, int hardclip_end,
			int querylength, int quality_shift,
			int pathnum, int npaths_primary, int npaths_altloc, int absmq_score, int second_absmq, unsigned int flag,
			Chrpos_T chrpos, Chrpos_T chrlength,
			
			Shortread_T queryseq, int pair_mapq_score, int end_mapq_score,
			Stage3pair_T stage3pair, Stage3end_T stage3hit,	bool first_read_p,
			char *sam_read_group_id) {
  bool invertp = false;

  /* Should already be checked when Stage3_T or Stage3end_T object was created */
  if (cigar_action == CIGAR_ACTION_IGNORE) {
    /* Don't check */
  } else if (cigar_length(cigar_tokens) + hardclip_start + hardclip_end == querylength) {
    /* Okay */
  } else if (cigar_action == CIGAR_ACTION_WARNING) {
    fprintf(stderr,"Warning: for %s, CIGAR length %d plus hardclips %d and %d do not match sequence length %d\n",
	    acc1,cigar_length(cigar_tokens),hardclip_start,hardclip_end,querylength);
  } else if (cigar_action == CIGAR_ACTION_NOPRINT) {
    fprintf(stderr,"Warning: for %s, CIGAR length %d plus hardclips %d and %d do not match sequence length %d\n",
	    acc1,cigar_length(cigar_tokens),hardclip_start,hardclip_end,querylength);
    return;
  } else {
    /* CIGAR_ACTION_ABORT */
    fprintf(stderr,"Error: for %s, CIGAR length %d plus hardclips %d and %d do not match sequence length %d\n",
	    acc1,cigar_length(cigar_tokens),hardclip_start,hardclip_end,querylength);
    abort();
  }

  /* 1. QNAME or Accession */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s\t",acc1);
  } else {
    FPRINTF(fp,"%s,%s\t",acc1,acc2);
  }

  /* 2. Flags */
  FPRINTF(fp,"%u\t",flag);

  /* 3. RNAME or Chrstring */
  /* 4. POS or Chrlow */
  /* Taken from GMAP part of SAM_chromosomal_pos */
  if (chrpos > chrlength) {
    FPRINTF(fp,"%s\t%u\t",chrstring,chrpos - chrlength /*+ 1*/);
  } else {
    FPRINTF(fp,"%s\t%u\t",chrstring,chrpos /*+ 1*/);
  }

  /* 5. MAPQ or Mapping quality */
  FPRINTF(fp,"%d\t",pair_mapq_score);

  /* 6. CIGAR */
  print_tokens(fp,cigar_tokens);

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  FPRINTF(fp,"\t*\t0");		/* Because sequence and mate are merged */

  /* 9. ISIZE: Insert size */
  FPRINTF(fp,"\t0");		/* Because sequence and mate are merged */

  /* 10. SEQ: queryseq and 11. QUAL: quality_scores */
  FPRINTF(fp,"\t");
  if (watsonp == true) {
    print_chopped(fp,queryseq_ptr,querylength,hardclip_start,hardclip_end);
    FPRINTF(fp,"\t");
    print_quality(fp,quality_string,querylength,hardclip_start,hardclip_end,
		  quality_shift);
  } else {
    print_chopped_revcomp(fp,queryseq_ptr,querylength,hardclip_start,hardclip_end);
    FPRINTF(fp,"\t");
    print_quality_revcomp(fp,quality_string,querylength,hardclip_start,hardclip_end,
			  quality_shift);
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH and XI */
  if (hardclip_start > 0 || hardclip_end > 0) {
    FPRINTF(fp,"\tXH:Z:");
    if (watsonp == true) {
      print_chopped_end(fp,queryseq_ptr,querylength,hardclip_start,hardclip_end);
    } else {
      print_chopped_end_revcomp(fp,queryseq_ptr,querylength,hardclip_start,hardclip_end);
    }

    if (quality_string != NULL) {
      FPRINTF(fp,"\tXI:Z:");
      if (watsonp == true) {
	print_chopped_end_quality(fp,quality_string,querylength,hardclip_start,hardclip_end);
      } else {
	print_chopped_end_quality_reverse(fp,quality_string,querylength,hardclip_start,hardclip_end);
      }
    }
  }

  if (queryseq != NULL) {
    /* 12. TAGS: XB */
    Shortread_print_barcode(fp,queryseq);

    /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
    Shortread_print_chop(fp,queryseq,invertp);
  }

  /* 12. TAGS: MD string */
  FPRINTF(fp,"\tMD:Z:");
  print_tokens(fp,md_tokens);

  /* 12. TAGS: NH */
  FPRINTF(fp,"\tNH:i:%d",npaths_primary + npaths_altloc);
  
  /* 12. TAGS: HI */
  FPRINTF(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  FPRINTF(fp,"\tNM:i:%d",nmismatches_refdiff + nindels);

  if (snps_p) {
    /* 12. TAGS: XW and XV */
    FPRINTF(fp,"\tXW:i:%d",nmismatches_bothdiff);
    FPRINTF(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }


  /* 12. TAGS: SM */
  FPRINTF(fp,"\tSM:i:%d",end_mapq_score);

  /* 12. TAGS: XQ */
  FPRINTF(fp,"\tXQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  FPRINTF(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  if (novelsplicingp == false && splicesites_iit == NULL) {
    /* Do not print XS field */

  } else if (sensedir == SENSE_FORWARD) {
    if (watsonp == true) {
      FPRINTF(fp,"\tXS:A:+");
    } else {
      FPRINTF(fp,"\tXS:A:-");
    }

  } else if (sensedir == SENSE_ANTI) {
    if (watsonp == true) {
      FPRINTF(fp,"\tXS:A:-");
    } else {
      FPRINTF(fp,"\tXS:A:+");
    }

  } else if (intronp == false) {
    /* Skip.  No intron in this end and mate is not revealing. */

#if 0
  } else if (force_xs_direction_p == true) {
    /* Don't print XS field for SENSE_NULL */
    /* Could not determine sense, so just report arbitrarily as + */
    /* This option provided for users of Cufflinks, which cannot handle XS:A:? */
    FPRINTF(fp,"\tXS:A:+");
    
  } else {
    /* Non-canonical.  Don't report. */
    FPRINTF(fp,"\tXS:A:?");
#endif
  }


  /* 12. TAGS: XX, XY */
  if (stage3pair == NULL) {
    /* Single-end */
    if (Stage3end_transcripts(stage3hit) != NULL) {
      FPRINTF(fp,"\tXX:Z:");
      Transcript_print_info(fp,Stage3end_transcripts(stage3hit),transcript_iit,invertp);
    }

  } else {
    /* Paired-end */
    if (first_read_p == true) {
      if (Stage3pair_transcripts5(stage3pair) != NULL) {
	FPRINTF(fp,"\tXX:Z:");
	Transcript_print_info(fp,Stage3pair_transcripts5(stage3pair),transcript_iit,invertp);
      }
      Transcript_print_diff(fp,Stage3end_transcripts(stage3hit),Stage3pair_transcripts5(stage3pair),
			    transcript_iit,invertp);

    } else {
      if (Stage3pair_transcripts3(stage3pair) != NULL) {
	FPRINTF(fp,"\tXX:Z:");
	Transcript_print_info(fp,Stage3pair_transcripts3(stage3pair),transcript_iit,invertp);
      }
      Transcript_print_diff(fp,Stage3end_transcripts(stage3hit),Stage3pair_transcripts3(stage3pair),
			    transcript_iit,invertp);
    }
  }

  /* TAGS: XZ */
  if (Stage3end_transcripts_other(stage3hit) != NULL) {
    FPRINTF(fp,"\tXZ:Z:");
    Transcript_print_info(fp,Stage3end_transcripts_other(stage3hit),transcript_iit,invertp);
  }

  /* 12. TAGS: XC (circular).  Not currently handled for GMAP */

  /* 12. TAGS: XG.  Not currently handled for GMAP */
  Method_samprint(fp,Stage3end_method(stage3hit));


  FPRINTF(fp,"\n");

  return;
}


/* Modified from Pair_print_sam */
void
Simplepair_overlap_print_sam (Filestring_T fp, char *abbrev, struct T *pairarray, int npairs,
			      char *acc1, char *acc2, Chrnum_T chrnum, Univ_IIT_T chromosome_iit,
			      char *queryseq_ptr, char *quality_string,
			      int hardclip_low, int hardclip_high, int querylength_given,
			      bool watsonp, int sensedir,
			      int quality_shift, bool first_read_p, int pathnum, int npaths_primary, int npaths_altloc,
			      int absmq_score, int second_absmq, Chrpos_T chrpos, Chrpos_T chrlength,

			      Shortread_T queryseq, unsigned int flag,
			      int pair_mapq_score, int end_mapq_score,
			      Stage3pair_T stage3pair, Stage3end_T stage3hit,
			      
			      char *sam_read_group_id) {
  char *chrstring = NULL;
  int chimera_part = 0;
#if 0
  char *mate_chrstring, *mate_chrstring_alloc = NULL;
#elif defined(GSNAP)
#else
  unsigned int flag;
#endif

  List_T cigar_tokens, md_tokens = NULL;
  int nmismatches_refdiff, nmismatches_bothdiff, nindels;
  bool intronp;
  int hardclip_start, hardclip_end;
  /* int hardclip_start_zero = 0, hardclip_end_zero = 0; */
  struct T *clipped_pairarray;
  int clipped_npairs;
  bool cigar_tokens_alloc;


  if (chrnum == 0) {
    /* chrstring = Sequence_accession(usersegment); */
    fprintf(stderr,"Error in Simplepair_overlap_print_sam: chrnum == 0\n");
    exit(9);
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

#if 0
  /* Previously this was code for GSNAP.  However for GSNAP, this is called only when the sequence and mate are merged */
  if (mate_chrpos_low == 0U) {
    mate_chrstring = "*";
  } else if (mate_chrnum == 0) {
    abort();
  } else if (/* chrpos > 0U && chrnum > 0 && */ mate_chrnum == chrnum) {
    mate_chrstring = "=";
  } else {
    mate_chrstring = mate_chrstring_alloc = Chrnum_to_string(mate_chrnum,chromosome_iit);
  }
#elif defined(GSNAP)
  /* flag is given as a parameter */
#else
  flag = compute_sam_flag_nomate(npaths_primary + npaths_altloc,first_read_p,watsonp,sam_paired_p);
#endif

  debug4(printf("Entered SAM_print_pairs with watsonp %d, first_read_p %d, hardclip_low %d, and hardclip_high %d\n",
		watsonp,first_read_p,hardclip_low,hardclip_high));

  if (watsonp == true) {
    hardclip_start = hardclip_low;
    hardclip_end = hardclip_high;
  } else {
    hardclip_start = hardclip_high;
    hardclip_end = hardclip_low;
  }
  debug4(printf("hardclip_start %d, hardclip_end %d\n",hardclip_start,hardclip_end));


  /* Because merged_overlap_p is true */
  /* clipped_pairarray = pairarray; */
  /* clipped_npairs = npairs; */
  clipped_pairarray = hardclip_pairarray(&clipped_npairs,hardclip_start,hardclip_end,
					 pairarray,npairs,querylength_given);
  cigar_tokens = compute_cigar(&intronp,&hardclip_start,&hardclip_end,clipped_pairarray,clipped_npairs,querylength_given,
			       watsonp,chimera_part);
  cigar_tokens_alloc = true;


  /* Cigar updates hardclip5 and hardclip3 for chimeras */
  md_tokens = compute_md_string(&nmismatches_refdiff,&nmismatches_bothdiff,&nindels,
				clipped_pairarray,clipped_npairs,watsonp,cigar_tokens);

  print_overlap_sam_line(fp,abbrev,acc1,acc2,chrstring,
			 watsonp,sensedir,cigar_tokens,md_tokens,
			 nmismatches_refdiff,nmismatches_bothdiff,nindels,
			 intronp,queryseq_ptr,quality_string,hardclip_start,hardclip_end,
			 querylength_given,quality_shift,pathnum,npaths_primary,npaths_altloc,
			 absmq_score,second_absmq,flag,chrpos,chrlength,
			 queryseq,pair_mapq_score,end_mapq_score,stage3pair,stage3hit,first_read_p,
			 sam_read_group_id);

  /* Print procedures free the character strings */
  tokens_free(&md_tokens);
  if (cigar_tokens_alloc == true) {
    tokens_free(&cigar_tokens);
  }

  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Simplepair_setup (bool novelsplicingp_in, IIT_T splicesites_iit_in,
		  Univ_IIT_T transcript_iit_in, bool sam_insert_0M_p_in,
		  bool md_lowercase_variant_p_in, bool snps_p_in,
		  bool cigar_extended_p_in, Cigar_action_T cigar_action_in) {

  novelsplicingp = novelsplicingp_in;
  splicesites_iit = splicesites_iit_in;
  transcript_iit = transcript_iit_in;

  sam_insert_0M_p = sam_insert_0M_p_in;
  md_lowercase_variant_p = md_lowercase_variant_p_in;
  snps_p = snps_p_in;

  cigar_extended_p = cigar_extended_p_in;
  cigar_action = cigar_action_in;

  return;
}
