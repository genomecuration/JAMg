static char rcsid[] = "$Id: transcript.c 212789 2018-01-26 14:08:00Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "transcript.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For qsort */
#include "mem.h"


#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


static Chrpos_T pairmax_transcriptome;
static Chrpos_T expected_pairlength;


void
Transcript_setup (Chrpos_T pairmax_transcriptome_in, Chrpos_T expected_pairlength_in) {
  pairmax_transcriptome = pairmax_transcriptome_in;
  expected_pairlength = expected_pairlength_in;
  return;
}


#define T Transcript_T


int
Transcript_num (T this) {
  return this->num;
}


void
Transcript_free (T *old) {

  FREE(*old);

  return;
}

void
Transcript_gc (List_T *list) {
  List_T p;
  T this;

  for (p = *list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    Transcript_free(&this);
  }
  List_free(&(*list));

  return;
}


T
Transcript_new (int num, int start, int end) {
  T new = (T) MALLOC(sizeof(*new));

  new->num = num;
  new->start = start;
  new->end = end;

  return new;
}


static T
Transcript_copy (T old) {
  T new = (T) MALLOC(sizeof(*new));

  new->num = old->num;
  new->start = old->start;
  new->end = old->end;

  return new;
}

List_T
Transcript_copy_list (List_T old) {
  List_T new = NULL, p;
  T this;

  for (p = old; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    new = List_push(new,(void *) Transcript_copy(this));
  }
  return List_reverse(new);
}

static bool
Transcript_equal (T x, T y) {
  if (x->num == y->num &&
      x->start == y->start &&
      x->end == y->end) {
    return true;
  } else {
    return false;
  }
}

bool
Transcript_in_list_p (T x, List_T list) {
  List_T p;
  T y;

  for (p = list; p != NULL; p = List_next(p)) {
    y = (T) List_head(p);
    if (x->num == y->num &&
	x->start == y->start &&
	x->end == y->end) {
      return true;
    }
  }
  return false;
}


void
Transcript_print_nums (List_T list) {
  List_T p;
  T this;

  printf(" Trnums:");
  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    printf(" %d:%d",this->num,this->start);
  }
  return;
}

void
Transcript_print_list (List_T list) {
  List_T p;
  T this;

  printf(" Trnums:");
  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    printf(" %d",this->num);
  }
  printf("\n");

  printf(" Trstarts:");
  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    printf(" %d",this->start);
  }
  printf("\n");

  printf(" Trends:");
  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    printf(" %d",this->end);
  }
  printf("\n");

  return;
}


static int
Transcript_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->num < y->num) {
    return -1;
  } else if (y->num < x->num) {
    return +1;
  } else if (x->start < y->start) {
    return -1;
  } else if (y->start < x->start) {
    return +1;
  } else if (x->end < y->end) {
    return -1;
  } else if (y->end < x->end) {
    return +1;
  } else {
    return 0;
  }
}


bool
Transcript_intersect_p (int *min_insertlength, List_T transcripts5, List_T transcripts3) {
  T *array5, *array3;
  int ntranscripts5, ntranscripts3, i, j, k, l;
  int a, b;
  int trnum5, trnum3;
  Chrpos_T start5, end5, end3;
  int insertlength;

  *min_insertlength = -1;

  if (transcripts5 != NULL && transcripts3 != NULL) {
    array5 = (T *) List_to_array_n(&ntranscripts5,transcripts5);
    array3 = (T *) List_to_array_n(&ntranscripts3,transcripts3);
    qsort(array5,ntranscripts5,sizeof(T),Transcript_cmp);
    qsort(array3,ntranscripts3,sizeof(T),Transcript_cmp);

    i = k = 0;
    while (i < ntranscripts5 && k < ntranscripts3) {
      trnum5 = array5[i]->num;
      j = i+1;
      while (j < ntranscripts5 && array5[j]->num == trnum5) {
	j++;
      }

      trnum3 = array3[k]->num;
      l = k+1;
      while (l < ntranscripts3 && array3[l]->num == trnum3) {
	l++;
      }
      
      if (trnum5 < trnum3) {
	i = j;

      } else if (trnum3 < trnum5) {
	k = l;

      } else {
	for (a = i; a < j; a++) {
	  start5 = (Chrpos_T) array5[a]->start;
	  end5 = (Chrpos_T) array5[a]->end;

	  if (start5 < end5) {
	    /* Plus on transcript */
	    for (b = k; b < l; b++) {
	      end3 = (Chrpos_T) array3[b]->end;
	      debug2(printf("Checking for concordance on trnum %d: start5 %u, end3 %u",trnum5,start5,end3));
	      if (start5 < end3 && end3 <= start5 + pairmax_transcriptome) {
		insertlength = (int) (end3 - start5);
		if (*min_insertlength < 0 || insertlength < *min_insertlength) {
		  *min_insertlength = insertlength;
		}
	      }
	    }

	  } else {
	    /* Minus on transcript */
	    for (b = k; b < l; b++) {
	      end3 = (Chrpos_T) array3[b]->end;
	      debug2(printf("Checking for concordance on trnum %d: end3 %u, start5 %u",trnum5,end3,start5));
	      if (end3 < start5 && start5 <= end3 + pairmax_transcriptome) {
		insertlength = (int) (start5 - end3);
		if (*min_insertlength < 0 || insertlength < *min_insertlength) {
		  *min_insertlength = insertlength;
		}
	      }
	    }
	  }
	}

	i = j;
	k = l;
      }
    }

    FREE(array3);
    FREE(array5);
  }

  if (*min_insertlength < 0) {
    return false;
  } else {
    return true;
  }
}


bool
Transcript_concordant_p (T transcript5, T transcript3) {
  int trnum5, trnum3;
  Chrpos_T start5, end5, end3;

  trnum5 = transcript5->num;
  trnum3 = transcript3->num;
  if (trnum5 != trnum3) {
    return false;
  } else {
    start5 = (Chrpos_T) transcript5->start;
    end5 = (Chrpos_T) transcript5->end;

    if (start5 < end5) {
      /* Plus on transcript */
      end3 = (Chrpos_T) transcript3->end;
      debug2(printf("Checking for concordance on trnum %d: start5 %u, end3 %u",trnum5,start5,end3));
      if (start5 < end3 && end3 <= start5 + pairmax_transcriptome) {
	return true;
      } else {
	return false;
      }
	
    } else {
      /* Minus on transcript */
      end3 = (Chrpos_T) transcript3->end;
      debug2(printf("Checking for concordance on trnum %d: end3 %u, start5 %u",trnum5,end3,start5));
      if (end3 < start5 && start5 <= end3 + pairmax_transcriptome) {
	return true;
      } else {
	return false;
      }
    }
  }
}



void
Transcript_concordance (List_T *newtranscripts5, List_T *newtranscripts3, List_T transcripts5, List_T transcripts3) {
  T *array5, *array3;
  int ntranscripts5, ntranscripts3, i, j, k, l;
  int besta, bestb, a, b;
  int trnum5, trnum3;
  Chrpos_T start5, end5, end3;
  Chrpos_T best_absdifflength, absdifflength, insertlength;


  *newtranscripts5 = (List_T) NULL;
  *newtranscripts3 = (List_T) NULL;

  if (transcripts5 != NULL && transcripts3 != NULL) {
    array5 = (T *) List_to_array_n(&ntranscripts5,transcripts5);
    array3 = (T *) List_to_array_n(&ntranscripts3,transcripts3);
    qsort(array5,ntranscripts5,sizeof(T),Transcript_cmp);
    qsort(array3,ntranscripts3,sizeof(T),Transcript_cmp);

    i = k = 0;
    while (i < ntranscripts5 && k < ntranscripts3) {
      trnum5 = array5[i]->num;
      j = i+1;
      while (j < ntranscripts5 && array5[j]->num == trnum5) {
	j++;
      }

      trnum3 = array3[k]->num;
      l = k+1;
      while (l < ntranscripts3 && array3[l]->num == trnum3) {
	l++;
      }
      
      if (trnum5 < trnum3) {
	i = j;

      } else if (trnum3 < trnum5) {
	k = l;

      } else {
	besta = bestb = -1;
	best_absdifflength = 0;

	for (a = i; a < j; a++) {
	  start5 = (Chrpos_T) array5[a]->start;
	  end5 = (Chrpos_T) array5[a]->end;

	  if (start5 < end5) {
	    /* Plus on transcript */
	    for (b = k; b < l; b++) {
	      end3 = (Chrpos_T) array3[b]->end;
	      debug2(printf("Checking for concordance on trnum %d: start5 %u, end3 %u",trnum5,start5,end3));
	      if (start5 < end3 && end3 <= start5 + pairmax_transcriptome) {
		if ((insertlength = end3 - start5) < expected_pairlength) {
		  absdifflength = expected_pairlength - insertlength;
		} else {
		  absdifflength = insertlength - expected_pairlength;
		}
		if (besta < 0 || absdifflength < best_absdifflength) {
		  besta = a;
		  bestb = b;
		  best_absdifflength = absdifflength;
		}
	      }
	    }

	  } else {
	    /* Minus on transcript */
	    for (b = k; b < l; b++) {
	      end3 = (Chrpos_T) array3[b]->end;
	      debug2(printf("Checking for concordance on trnum %d: end3 %u, start5 %u",trnum5,end3,start5));
	      if (end3 < start5 && start5 <= end3 + pairmax_transcriptome) {
		if ((insertlength = start5 - end3) < expected_pairlength) {
		  absdifflength = expected_pairlength - insertlength;
		} else {
		  absdifflength = insertlength - expected_pairlength;
		}
		if (besta < 0 || absdifflength < best_absdifflength) {
		  besta = a;
		  bestb = b;
		  best_absdifflength = absdifflength;
		}
	      }
	    }
	  }
	}

	if (besta >= 0) {
	  *newtranscripts5 = List_push(*newtranscripts5,(void *) Transcript_copy(array5[besta]));
	  *newtranscripts3 = List_push(*newtranscripts3,(void *) Transcript_copy(array3[bestb]));
	}

	i = j;
	k = l;
      }
    }

    FREE(array3);
    FREE(array5);
  }

  return;
}



void
Transcript_print_info (Filestring_T fp, List_T transcripts, Univ_IIT_T transcript_iit, bool invertp) {
  T this, *array;
  char *label;
  bool allocp;
  int n, i;

  if (transcripts != NULL) {
    array = (T *) List_to_array_n(&n,transcripts);
    qsort(array,n,sizeof(T),Transcript_cmp);

    this = array[0];
    label = Univ_IIT_label(transcript_iit,this->num,&allocp);
    if (invertp == false) {
      if (this->start < this->end) {
	FPRINTF(fp,"+%s:%d..%d",label,this->start+1,this->end);
      } else {
	FPRINTF(fp,"-%s:%d..%d",label,this->start,this->end+1);
      }
    } else {
      if (this->start < this->end) {
	FPRINTF(fp,"-%s:%d..%d",label,this->end,this->start+1);
      } else {
	FPRINTF(fp,"+%s:%d..%d",label,this->end+1,this->start);
      }
    }
    if (allocp) {
      FREE(label);
    }

    for (i = 1; i < n; i++) {
      this = array[i];
      label = Univ_IIT_label(transcript_iit,this->num,&allocp);
      if (invertp == false) {
	if (this->start < this->end) {
	  FPRINTF(fp,",+%s:%d..%d",label,this->start+1,this->end);
	} else {
	  FPRINTF(fp,",-%s:%d..%d",label,this->start,this->end+1);
	}
      } else {
	if (this->start < this->end) {
	  FPRINTF(fp,",-%s:%d..%d",label,this->end,this->start+1);
	} else {
	  FPRINTF(fp,",+%s:%d..%d",label,this->end+1,this->start);
	}
      }
      if (allocp) {
	FREE(label);
      }
    }

    FREE(array);
  }

  return;
}


/* Only for SAM output */
void
Transcript_print_diff (Filestring_T fp, List_T transcripts, List_T common_transcripts,
		       Univ_IIT_T transcript_iit, bool invertp) {
  bool disjointp, equalp, firstp;
  List_T p, q;
  T this;
  char *label;
  bool allocp;

  p = transcripts;
  disjointp = false;
  while (disjointp == false && p != NULL) {
    this = (T) List_head(p);
    q = common_transcripts;
    equalp = false;
    while (equalp == false && q != NULL) {
      if (Transcript_equal(this,(T) List_head(q)) == true) {
	equalp = true;
      }
      q = List_next(q);
    }
    if (equalp == false) {
      disjointp = true;
    }
    p = List_next(p);
  }

  if (disjointp == true) {
    FPRINTF(fp,"\tXY:Z:");
    firstp = true;

    for (p = transcripts; p != NULL; p = List_next(p)) {
      this = (T) List_head(p);

      equalp = false;
      for (q = common_transcripts; q != NULL; q = List_next(q)) {
	if (Transcript_equal(this,(T) List_head(q)) == true) {
	  equalp = true;
	}
      }

      if (equalp == false) {
	label = Univ_IIT_label(transcript_iit,this->num,&allocp);
	if (firstp == false) {
	  FPRINTF(fp,",");
	}
	if (invertp == false) {
	  if (this->start < this->end) {
	    FPRINTF(fp,"+%s:%d..%d",label,this->start+1,this->end);
	  } else {
	    FPRINTF(fp,"-%s:%d..%d",label,this->start,this->end+1);
	  }
	} else {
	  if (this->start < this->end) {
	    FPRINTF(fp,"-%s:%d..%d",label,this->end,this->start+1);
	  } else {
	    FPRINTF(fp,"+%s:%d..%d",label,this->end+1,this->start);
	  }
	}
	if (allocp) {
	  FREE(label);
	}
	firstp = false;
      }
    }
  }

  return;
}


