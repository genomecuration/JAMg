static char rcsid[] = "$Id: ladder.c 218473 2019-02-22 23:39:06Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "ladder.h"
#include "stage3hrdef.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>		/* For memcpy */
#include "assert.h"
#include "mem.h"
#include "sedgesort.h"


#define MAX_HITS 1000

static int
genomicstart_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;

  if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
  } else {
    return 0;
  }
}

static int
genomicend_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;

  if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}

static int
stream_length (Stage3end_T *stream) {
  int length = 0;

  while (*stream++ != NULL) {
    length += 1;
  }

  return length;
}


#define T Ladder_T
struct T {
  int maxscore;
  bool end5p;

  bool *sortedp;
  Stage3end_T **byscore;	/* Each array is terminated by an element of NULL, to be compatible with streams */
  int *nhits;

  List_T *streams;		/* Each stream has type (Stage3end_T *), terminated by an element of
				   NULL.  Managed by Hitlistpool_T procedures */

  List_T duplicates;		/* Contains results from
				   Stage3end_remove_duplicates_array.
				   Hits are potentially shared between
				   newladder and ladder, so free only
				   from ladder after concordance is
				   performed.  Managed by Hitlistpool_T procedures */
};


#ifdef CHECK_ASSERTIONS
static void
Ladder_check (T this) {
  int score;
  int nhits, i;
  List_T p;
  Stage3end_T *stream;

  for (score = 0; score <= this->maxscore; score++) {
    if (this->sortedp[score] == true) {
      assert(this->streams[score] == (List_T) NULL);
      
      if (this->nhits[score] == 0) {
	assert(this->byscore[score] == (Stage3end_T *) NULL);
      } else {
	assert(this->byscore[score] != (Stage3end_T *) NULL);
	for (i = 0; i < this->nhits[score]; i++) {
	  assert(this->byscore[score][i] != (Stage3end_T) NULL);
	}
	assert(this->byscore[score][this->nhits[score]] == (Stage3end_T) NULL);
      }

    } else {
      assert(this->byscore[score] == (Stage3end_T *) NULL);

      if (this->nhits[score] == 0) {
	assert(this->streams[score] == (List_T) NULL);
      } else {
	assert(this->streams[score] != (List_T) NULL);

	nhits = 0;
	for (p = this->streams[score]; p != NULL; p = List_next(p)) {
	  stream = (Stage3end_T *) List_head(p);
	  while (*stream++ != NULL) {
	    nhits++;
	  }
	}
	assert(nhits == this->nhits[score]);
      }
    }
  }

  return;
}
#endif


int
Ladder_maxscore (T this) {
  return this->maxscore;
}


int
Ladder_cutoff (T this) {
  int score = -1;
  int total_nhits = 0;

  return this->maxscore;

  /* Commented out because we want to see at least one score */
  while (score+1 <= this->maxscore && total_nhits /*+ this->nhits[score+1]*/ <= MAX_HITS) {
    total_nhits += this->nhits[++score];
  }

  return score;
}


Stage3end_T *
Ladder_hits_for_score (int *nhits, T this, Hitlistpool_T hitlistpool, int score) {
  List_T q;
  Stage3end_T *hits, *in, *out, *stream;
  int n;

#ifdef CHECK_ASSERTIONS
  Ladder_check(this);
#endif

  if ((n = this->nhits[score]) == 0) {
    *nhits = 0;
    return (Stage3end_T *) NULL;

  } else if (this->sortedp[score] == true) {
    *nhits = n;
    return this->byscore[score];

  } else {
    out = hits = (Stage3end_T *) MALLOC((n + 1) * sizeof(Stage3end_T));
    for (q = this->streams[score]; q != NULL; q = List_next(q)) {
      in = stream = (Stage3end_T *) List_head(q);
      while (*in != (Stage3end_T) NULL) {
	*out++ = *in++;
      }
      FREE(stream);
    }
    Hitlist_free(&this->streams[score]);
    this->streams[score] = (List_T) NULL;

    /* Need to remove duplicates at this point to make concordance go
       faster and to avoid possible overflows */
    this->byscore[score] = Stage3end_remove_duplicates_array(&(*nhits),&this->duplicates,hits,n,
							     hitlistpool);
    this->nhits[score] = *nhits;
    this->byscore[score][*nhits] = (Stage3end_T) NULL; /* Terminate with a NULL */

    /* Sort for concordance */
    if (this->end5p == true) {
      qsort(this->byscore[score],*nhits,sizeof(Stage3end_T),genomicend_cmp);
    } else {
      qsort(this->byscore[score],*nhits,sizeof(Stage3end_T),genomicstart_cmp);
    }
      
    this->sortedp[score] = true;

#ifdef CHECK_ASSERTIONS
    Ladder_check(this);
#endif

    return this->byscore[score];
  }

}

int *
Ladder_nhits (T this) {
  return this->nhits;
}

static List_T
Ladder_hitlist (T this, List_T hitlist, Hitlistpool_T hitlistpool) {
  int score;
  List_T p;
  Stage3end_T *stream;

  for (score = 0; score <= this->maxscore; score++) {
    if (this->nhits[score] == 0) {
      /* Skip */
    } else if (this->sortedp[score] == true) {
      stream = this->byscore[score];
      while (*stream != NULL) {
	hitlist = Hitlist_push(hitlist,hitlistpool,(void *) *stream++);
      }

    } else {
      for (p = this->streams[score]; p != NULL; p = List_next(p)) {
	stream = (Stage3end_T *) List_head(p);
	while (*stream != NULL) {
	  hitlist = Hitlist_push(hitlist,hitlistpool,(void *) *stream++);
	}
      }
    }
  }

  return hitlist;
}


Univcoord_T *
Ladder_genomicstarts (int *ndiagonals, T this) {
  Univcoord_T *diagonals;
  int score, k = 0;
  List_T p;
  Stage3end_T *stream;

  *ndiagonals = 0;
  for (score = 0; score <= this->maxscore; score++) {
    *ndiagonals += this->nhits[score];
  }

  /* Need extra entry for Sedgesort */
  diagonals = (Univcoord_T *) MALLOC(((*ndiagonals) + 1) * sizeof(Univcoord_T));

  for (score = 0; score <= this->maxscore; score++) {
    if (this->nhits[score] == 0) {
      /* Skip */
    } else if (this->sortedp[score] == true) {
      stream = this->byscore[score];
      while (*stream != NULL) {
	diagonals[k++] = Stage3end_genomicstart(*stream++);
      }

    } else {
      for (p = this->streams[score]; p != NULL; p = List_next(p)) {
	stream = (Stage3end_T *) List_head(p);
	while (*stream != NULL) {
	  diagonals[k++] = Stage3end_genomicstart(*stream++);
	}
      }
    }
  }
  
#ifdef LARGE_GENOMES
  Sedgesort_uint8(diagonals,*ndiagonals);
#else
  Sedgesort_uint4(diagonals,*ndiagonals);
#endif

  return diagonals;
}


Univcoord_T *
Ladder_genomicends (int *ndiagonals, T this) {
  Univcoord_T *diagonals;
  int score, k = 0;
  List_T p;
  Stage3end_T *stream;

  *ndiagonals = 0;
  for (score = 0; score <= this->maxscore; score++) {
    *ndiagonals += this->nhits[score];
  }

  /* Need extra entry for Sedgesort */
  diagonals = (Univcoord_T *) MALLOC(((*ndiagonals) + 1) * sizeof(Univcoord_T));

  for (score = 0; score <= this->maxscore; score++) {
    if (this->nhits[score] == 0) {
      /* Skip */
    } else if (this->sortedp[score] == true) {
      stream = this->byscore[score];
      while (*stream != NULL) {
	diagonals[k++] = Stage3end_genomicend(*stream++);
      }

    } else {
      for (p = this->streams[score]; p != NULL; p = List_next(p)) {
	stream = (Stage3end_T *) List_head(p);
	while (*stream != NULL) {
	  diagonals[k++] = Stage3end_genomicend(*stream++);
	}
      }
    }
  }
  
#ifdef LARGE_GENOMES
  Sedgesort_uint8(diagonals,*ndiagonals);
#else
  Sedgesort_uint4(diagonals,*ndiagonals);
#endif

  return diagonals;
}




void
Ladder_gc_hits (T this) {
  int score;
  List_T p;
  Stage3end_T *stream, hit;

  for (score = 0; score <= this->maxscore; score++) {
    if (this->nhits[score] == 0) {
      /* Skip */
    } else if (this->sortedp[score] == true) {
      stream = this->byscore[score];
      while (*stream != NULL) {
	hit = *stream++;
	Stage3end_free(&hit);
      }

    } else {
      for (p = this->streams[score]; p != NULL; p = List_next(p)) {
	stream = (Stage3end_T *) List_head(p);
	while (*stream != NULL) {
	  hit = *stream++;
	  Stage3end_free(&hit);
	}
      }
    }
  }

  return;
}


int
Ladder_minimax_trim (T ladder_plus, T ladder_minus, int querylength) {
  int min_trim = querylength, trim;
  T this;
  int score;
  List_T p;
  Stage3end_T *stream, hit;

  this = ladder_plus;
  for (score = 0; score <= this->maxscore; score++) {
    if (this->nhits[score] == 0) {
      /* Skip */
    } else if (this->sortedp[score] == true) {
      stream = this->byscore[score];
      while ((hit = *stream++) != NULL) {
	if ((trim = Stage3end_max_trim(hit)) < min_trim) {
	  min_trim = trim;
	}
      }

    } else {
      for (p = this->streams[score]; p != NULL; p = List_next(p)) {
	stream = (Stage3end_T *) List_head(p);
	while ((hit = *stream++) != NULL) {
	  if ((trim = Stage3end_max_trim(hit)) < min_trim) {
	    min_trim = trim;
	  }
	}
      }
    }
  }


  this = ladder_minus;
  for (score = 0; score <= this->maxscore; score++) {
    if (this->nhits[score] == 0) {
      /* Skip */
    } else if (this->sortedp[score] == true) {
      stream = this->byscore[score];
      while ((hit = *stream++) != NULL) {
	if ((trim = Stage3end_max_trim(hit)) < min_trim) {
	  min_trim = trim;
	}
      }

    } else {
      for (p = this->streams[score]; p != NULL; p = List_next(p)) {
	stream = (Stage3end_T *) List_head(p);
	while ((hit = *stream++) != NULL) {
	  if ((trim = Stage3end_max_trim(hit)) < min_trim) {
	    min_trim = trim;
	  }
	}
      }
    }
  }

  return min_trim;
}



void
Ladder_free (T *old) {
  int score;
  List_T p;
  Stage3end_T *stream;

  for (score = 0; score <= (*old)->maxscore; score++) {
    FREE((*old)->byscore[score]);
    for (p = (*old)->streams[score]; p != NULL; p = List_next(p)) {
      stream = (Stage3end_T *) List_head(p);
      FREE(stream);
    }
    Hitlist_free(&(*old)->streams[score]);
  }
  FREE((*old)->sortedp);
  FREE((*old)->byscore);
  FREE((*old)->nhits);
  FREE((*old)->streams);
  Hitlist_free(&(*old)->duplicates);

  FREE(*old);

  return;
}


void
Ladder_gc_duplicates (T this) {
  List_T p;
  Stage3end_T hit;

  for (p = this->duplicates; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    Stage3end_free(&hit);
  }
  Hitlist_free(&this->duplicates);
  assert(this->duplicates == NULL);
  
  return;
}



void
Ladder_to_hits (List_T *hits5, List_T *hits3,
		T *ladder5_plus, T *ladder5_minus,
		T *ladder3_plus, T *ladder3_minus,
		Hitlistpool_T hitlistpool) {
  *hits5 = Ladder_hitlist(*ladder5_plus,NULL,hitlistpool);
  *hits5 = Ladder_hitlist(*ladder5_minus,*hits5,hitlistpool);

  *hits3 = Ladder_hitlist(*ladder3_plus,NULL,hitlistpool);
  *hits3 = Ladder_hitlist(*ladder3_minus,*hits3,hitlistpool);

  Ladder_gc_duplicates(*ladder5_plus);
  Ladder_gc_duplicates(*ladder5_minus);
  Ladder_gc_duplicates(*ladder3_plus);
  Ladder_gc_duplicates(*ladder3_minus);

  Ladder_free(&(*ladder5_plus));
  Ladder_free(&(*ladder5_minus));
  Ladder_free(&(*ladder3_plus));
  Ladder_free(&(*ladder3_minus));
  
  return;
}



T
Ladder_new (List_T hitlist, Hitlistpool_T hitlistpool, bool end5p) {
  T new = (T) MALLOC(sizeof(*new));
  List_T p;
  Stage3end_T *stream, hit;
  int maxscore, score;

  new->end5p = end5p;

  maxscore = 0;
  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    if (hit->score_within_trims > maxscore) {
      maxscore = hit->score_within_trims;
    }
  }

  new->maxscore = maxscore;
  new->sortedp = (bool *) CALLOC(maxscore+1,sizeof(bool));
  new->byscore = (Stage3end_T **) CALLOC(maxscore+1,sizeof(Stage3end_T *));
  new->nhits = (int *) CALLOC(maxscore+1,sizeof(int));
  new->streams = (List_T *) CALLOC(maxscore+1,sizeof(List_T));
  
  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    score = hit->score_within_trims;
    new->nhits[score]++;
  }

  for (score = 0; score <= maxscore; score++) {
    if (new->nhits[score] > 0) {
      stream = (Stage3end_T *) MALLOC((new->nhits[score]+1) * sizeof(Stage3end_T));
      stream[new->nhits[score]] = (Stage3end_T) NULL; /* Terminate with a NULL */
      new->streams[score] = Hitlist_push(NULL,hitlistpool,(void *) stream); /* Stream consists of a single array */
      new->nhits[score] = 0;	/* Reset for use below */
    }
  }

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (Stage3end_T) List_head(p);
    score = hit->score_within_trims;
    assert(score >= 0);
    stream = (Stage3end_T *) List_head(new->streams[score]);
    stream[new->nhits[score]++] = hit;
  }

  new->duplicates = (List_T) NULL;

  return new;
}


void
Ladder_merge (T dest, T source, Hitlistpool_T hitlistpool) {
  int score;
  bool *sortedp;
  Stage3end_T **byscore, *stream, *oldstream;
  int *nhits, length;
  List_T *streams, p;


  /* Revise byscore and nhits */
  if (source->maxscore <= dest->maxscore) {
    for (score = 0; score <= source->maxscore; score++) {
      if (source->nhits[score] > 0) {
	if (dest->sortedp[score] == true) {
	  assert(dest->streams[score] == NULL);
	  dest->streams[score] = Hitlist_push(NULL,hitlistpool,(void *) dest->byscore[score]);
	  dest->byscore[score] = (Stage3end_T *) NULL;
	  dest->sortedp[score] = false;
	}

	if (source->sortedp[score] == true) {
	  stream = (Stage3end_T *) MALLOC((source->nhits[score]+1)*sizeof(Stage3end_T));
	  memcpy(stream,source->byscore[score],(source->nhits[score]+1)*sizeof(Stage3end_T));
	  dest->streams[score] = Hitlist_push(dest->streams[score],hitlistpool,(void *) stream);
	} else {
	  for (p = source->streams[score]; p != NULL; p = List_next(p)) {
	    oldstream = (Stage3end_T *) List_head(p);
	    length = stream_length(oldstream);
	    stream = (Stage3end_T *) MALLOC((length+1)*sizeof(Stage3end_T));
	    memcpy(stream,oldstream,(length+1)*sizeof(Stage3end_T));
	    dest->streams[score] = Hitlist_push(dest->streams[score],hitlistpool,(void *) stream);
	  }
	}
	dest->nhits[score] += source->nhits[score];
      }
    }

  } else {
    sortedp = (bool *) CALLOC(source->maxscore+1,sizeof(bool));
    byscore = (Stage3end_T **) CALLOC(source->maxscore+1,sizeof(Stage3end_T *));
    nhits = (int *) CALLOC(source->maxscore+1,sizeof(int));
    streams = (List_T *) CALLOC(source->maxscore+1,sizeof(List_T));

    for (score = 0; score <= dest->maxscore; score++) {
      if (source->nhits[score] == 0) {
	/* Copy information from dest */
	sortedp[score] = dest->sortedp[score];
	byscore[score] = dest->byscore[score];
	nhits[score] = dest->nhits[score];
	streams[score] = dest->streams[score];

      } else {
	nhits[score] = dest->nhits[score] + source->nhits[score];
	if (dest->sortedp[score] == true) {
	  assert(dest->streams[score] == NULL);
	  streams[score] = Hitlist_push(NULL,hitlistpool,(void *) dest->byscore[score]);
	  byscore[score] = (Stage3end_T *) NULL;
	  sortedp[score] = false;
	} else {
	  streams[score] = dest->streams[score];
	}

	if (source->sortedp[score] == true) {
	  streams[score] = Hitlist_push(streams[score],hitlistpool,(void *) source->byscore[score]);
	} else {
	  for (p = source->streams[score]; p != NULL; p = List_next(p)) {
	    oldstream = (Stage3end_T *) List_head(p);
	    length = stream_length(oldstream);
	    stream = (Stage3end_T *) MALLOC((length+1)*sizeof(Stage3end_T));
	    memcpy(stream,oldstream,(length+1)*sizeof(Stage3end_T));
	    streams[score] = Hitlist_push(streams[score],hitlistpool,(void *) stream);
	  }
	}
      }
    }

    for ( ; score <= source->maxscore; score++) {
      if ((nhits[score] = source->nhits[score]) > 0) {
	if ((sortedp[score] = source->sortedp[score]) == true) {
	  byscore[score] = (Stage3end_T *) MALLOC((nhits[score]+1)*sizeof(Stage3end_T));
	  memcpy(byscore[score],source->byscore[score],(nhits[score]+1)*sizeof(Stage3end_T));
	} else {
	  for (p = source->streams[score]; p != NULL; p = List_next(p)) {
	    oldstream = (Stage3end_T *) List_head(p);
	    length = stream_length(oldstream);
	    stream = (Stage3end_T *) MALLOC((length+1)*sizeof(Stage3end_T));
	    memcpy(stream,oldstream,(length+1)*sizeof(Stage3end_T));
	    streams[score] = Hitlist_push(streams[score],hitlistpool,(void *) stream);
	  }
	}
      }
    }

    FREE(dest->sortedp);
    FREE(dest->byscore);
    FREE(dest->nhits);
    FREE(dest->streams);

    dest->maxscore = source->maxscore;
    dest->sortedp = sortedp;
    dest->byscore = byscore;
    dest->nhits = nhits;
    dest->streams = streams;
  }

  /* Because source was just produced by Ladder_new */
  assert(source->duplicates == (List_T) NULL);
  /* dest->duplicates = List_append(dest->duplicates,source->duplicates); */

#ifdef CHECK_ASSERTIONS
  Ladder_check(source);
  Ladder_check(dest);
#endif

  return;
}



#if 0
int
Ladder_cutoff_level (T this) {
#if 0
  /* Reset cutoff_level to check up to MAX_HITS */
  debug5(printf("Computing cutoff level up to max_score %d and MAX_HITS %d\n",max_score,MAX_HITS));
  score = 0;
  total_nhits = nhits[score];
  while (score+1 <= max_score && total_nhits + nhits[score+1] < MAX_HITS) {
    total_nhits += nhits[score+1];
    score++;
    debug5(printf("Allowing score to go to %d, because nhits = %d < max_hits %d\n",score,total_nhits,MAX_HITS));
  }
  debug5(printf("Setting cutoff_level to be %d\n\n",score));
  cutoff_level = score;
#else
  cutoff_level = max_score;
#endif

  return cutoff_level;
}
#endif


