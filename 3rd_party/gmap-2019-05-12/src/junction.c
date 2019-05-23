static char rcsid[] = "$Id: junction.c 218255 2019-01-22 17:19:58Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "junction.h"
#include "mem.h"
#include "assert.h"
#include "complement.h"
#include "maxent_hr.h"
#include "sense.h"


#define MIN_INTRONLEN 30


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Debugging procedures */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif



#define T Junction_T
struct T {
  Junctiontype_T type;
  int nindels;			/* Should be positive */
  Univcoord_T deletionpos;

  Chrpos_T splice_distance;
  int sensedir;
  double donor_prob;
  double acceptor_prob;
};


#ifdef DEBUG1
void
Junction_print (T this) {
  if (this == NULL) {
    printf("No junction\n");
  } else if (this->type == INS_JUNCTION) {
    printf("Insertion of %d\n",this->nindels);
  } else if (this->type == DEL_JUNCTION) {
    printf("Deletion of %d at %llu\n",this->nindels,(unsigned long long) this->deletionpos);
  } else if (this->type == SPLICE_JUNCTION) {
    if (this->splice_distance == 0) {
      printf("Splice ambiguous with sense %d, prob %f and %f\n",
	     this->sensedir,this->donor_prob,this->acceptor_prob);
    } else {
      printf("Splice with sense %d of %u, prob %f and %f\n",
	     this->sensedir,this->splice_distance,this->donor_prob,this->acceptor_prob);
    }
  }
  return;
}
#endif

#ifdef DEBUG1
void
Junction_print_list (List_T list) {
  T this;
  List_T p;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    if (this == NULL) {
      printf("None,");
    } else if (this->type == INS_JUNCTION) {
      printf("Ins:%d,",this->nindels);
    } else if (this->type == DEL_JUNCTION) {
      printf("Del:%d,",this->nindels);
    } else if (this->type == SPLICE_JUNCTION) {
      if (this->splice_distance == 0) {
	printf("Amb:%f-%f,",this->donor_prob,this->acceptor_prob);
      } else {
	printf("Splice:%u,",this->splice_distance);
      }
    }
  }

  return;
}
#endif

void
Junction_free (T *old) {
  FREE(*old);
  return;
}

void
Junction_gc (List_T *list) {
  List_T p;
  T old;

  for (p = *list; p != NULL; p = List_next(p)) {
    old = (T) List_head(p);
    Junction_free(&old);
  }
  /* List_free(&(*list)); -- Allocated by Listpool_push */
  return;
}

T
Junction_new_insertion (int nindels) {
  T new = (T) MALLOC(sizeof(*new));

  assert(nindels > 0);

  new->type = INS_JUNCTION;
  new->nindels = nindels;
  new->deletionpos = 0;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = 0.0;
  new->acceptor_prob = 0.0;

  return new;
}

T
Junction_new_deletion (int nindels, Univcoord_T deletionpos) {
  T new = (T) MALLOC(sizeof(*new));

  assert(nindels > 0);

  new->type = DEL_JUNCTION;
  new->nindels = nindels;
  new->deletionpos = deletionpos;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = 0.0;
  new->acceptor_prob = 0.0;

  return new;
}

T
Junction_new_splice (Chrpos_T splice_distance, int sensedir, double donor_prob, double acceptor_prob) {
  T new = (T) MALLOC(sizeof(*new));

  /* A zero splice distance is created for ambiguous splices */
  assert((int) splice_distance >= 0);

  new->type = SPLICE_JUNCTION;
  new->nindels = 0;
  new->deletionpos = 0;

  new->splice_distance = splice_distance;
  new->sensedir = sensedir;
  new->donor_prob = donor_prob;
  new->acceptor_prob = acceptor_prob;

  return new;
}


T
Junction_new_chimera (int sensedir, double donor_prob, double acceptor_prob) {
  T new = (T) MALLOC(sizeof(*new));

  new->type = CHIMERA_JUNCTION;
  new->nindels = 0;
  new->deletionpos = 0;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = donor_prob;
  new->acceptor_prob = acceptor_prob;

  return new;
}


T
Junction_new_generic (Univcoord_T left1, Univcoord_T left2, int querypos1, int querypos2,
		      Univcoord_T chroffset, bool plusp, int sensedir) {
  int nindels;
  double donor_prob, acceptor_prob;
  Chrpos_T splice_distance;

  debug(printf("Entered Junction_new_generic with left %u, querypos %d to left %u, querypos %d, plusp %d, sensedir %d\n",
	       left1,querypos1,left2,querypos2,plusp,sensedir));

  if (left1 + querypos1 > left2 + querypos2) {
    nindels = (left1 + querypos1) - (left2 + querypos2);
    debug(printf("Returning an insertion with %d indels\n",nindels));
    return Junction_new_insertion(nindels);

  } else if (left1 + querypos1 == left2 + querypos2) {
    debug(printf("Returning NULL\n"));
    return (Junction_T) NULL;
    
  } else if ((nindels = (left2 + querypos2) - (left1 + querypos1)) < MIN_INTRONLEN) {
    debug(printf("Returning a deletion with %d indels\n",nindels));
    return Junction_new_deletion(nindels,/*deletionpos*/left1+querypos1);

  } else {
    splice_distance = left2 - left1;
    if (sensedir == SENSE_FORWARD) {
      if (plusp == true) {
	donor_prob = Maxent_hr_donor_prob(left1 + querypos1,chroffset);
	acceptor_prob = Maxent_hr_acceptor_prob(left2 + querypos2,chroffset);
	debug(printf("Returning a splice with prob %f + %f\n",donor_prob,acceptor_prob));
	return Junction_new_splice(splice_distance,SENSE_FORWARD,donor_prob,acceptor_prob);
      } else {
	donor_prob = Maxent_hr_antidonor_prob(left2 + querypos2,chroffset);
	acceptor_prob = Maxent_hr_antiacceptor_prob(left1 + querypos1,chroffset);
	debug(printf("Returning a splice with prob %f + %f\n",donor_prob,acceptor_prob));
	return Junction_new_splice(splice_distance,SENSE_FORWARD,donor_prob,acceptor_prob);
      }
      
    } else {
      if (plusp == true) {
	donor_prob = Maxent_hr_antidonor_prob(left2 + querypos2,chroffset);
	acceptor_prob = Maxent_hr_antiacceptor_prob(left1 + querypos1,chroffset);
	debug(printf("Returning a splice with prob %f + %f\n",donor_prob,acceptor_prob));
	return Junction_new_splice(splice_distance,SENSE_ANTI,donor_prob,acceptor_prob);
      } else {
	donor_prob = Maxent_hr_donor_prob(left1 + querypos1,chroffset);
	acceptor_prob = Maxent_hr_acceptor_prob(left2 + querypos2,chroffset);
	debug(printf("Returning a splice with prob %f + %f\n",donor_prob,acceptor_prob));
	return Junction_new_splice(splice_distance,SENSE_ANTI,donor_prob,acceptor_prob);
      }
    }
  }
}


T
Junction_copy (T old) {
  T new = (T) MALLOC(sizeof(*new));

  new->type = old->type;
  new->nindels = old->nindels;
  new->deletionpos = old->deletionpos;

  new->splice_distance = old->splice_distance;
  new->sensedir = old->sensedir;
  new->donor_prob = old->donor_prob;
  new->acceptor_prob = old->acceptor_prob;

  return new;
}


List_T
Junction_copy_list (List_T old, Listpool_T listpool) {
  List_T new = NULL, p;

  for (p = old; p != NULL; p = List_next(p)) {
    new = Listpool_push(new,listpool,(void *) Junction_copy((T) List_head(p)));
  }
  return List_reverse(new);
}


Junctiontype_T
Junction_type (T this) {
  return this->type;
}

char *
Junction_typestring (T this) {
  switch (this->type) {
  case NO_JUNCTION: return "None";
  case INS_JUNCTION: return "Insertion";
  case DEL_JUNCTION: return "Deletion";
  case SPLICE_JUNCTION: return "Splice";
  case CHIMERA_JUNCTION: return "Chimera";
  case AMB_JUNCTION: return "Amb";
  case END_JUNCTION: return "End";
  }
  return (char *) NULL;
}

double
Junction_prob (T this) {
  return this->donor_prob + this->acceptor_prob;
}

int
Junction_sensedir (T this) {
  return this->sensedir;
}

double
Junction_donor_prob (T this) {
  return this->donor_prob;
}

double
Junction_acceptor_prob (T this) {
  return this->acceptor_prob;
}

double
Junction_splice_score (T this) {
  return this->donor_prob + this->acceptor_prob;
}

int
Junction_nindels (T this) {
  return this->nindels;
}

int
Junction_adj (T this) {
  if (this->type == DEL_JUNCTION) {
    return +this->nindels;
  } else if (this->type == INS_JUNCTION) {
    return -this->nindels;
  } else {
    return 0;
  }
}

int
Junction_ninserts (T this) {
  if (this->type == INS_JUNCTION) {
    return this->nindels;
  } else {
    return 0;
  }
}

int
Junction_total_ninserts (List_T list) {
  int ninserts = 0;
  T this;
  List_T p;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    if (this->type == INS_JUNCTION) {
      ninserts += this->nindels;
    }
  }

  return ninserts;
}



static char complCode[128] = COMPLEMENT_LC;

static char *
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

  return sequence;
}


Univcoord_T
Junction_deletionpos (T this) {
  return this->deletionpos;
}

void
Junction_set_deletionpos (T this, Univcoord_T deletionpos) {
  this->deletionpos = deletionpos;
  return;
}

char *
Junction_deletion_string (T this, Genome_T genome, bool plusp) {
  char *deletion_string;
  
  /* printf("Entered Junction_deletion_string with plusp %d\n",plusp); */
  /* printf("deletionpos = %u\n",this->deletionpos); */

  deletion_string = (char *) MALLOC((this->nindels+1)*sizeof(char));
  Genome_fill_buffer_simple(genome,this->deletionpos,this->nindels,deletion_string);
  if (plusp == false) {
    make_complement_inplace(deletion_string,this->nindels);
  }

  /* printf("string = %s\n",deletion_string); */
  return deletion_string;
}


Chrpos_T
Junction_splice_distance (T this) {
  return this->splice_distance;
}

void
Junction_set_unambiguous (T this, Chrpos_T distance, double donor_prob, double acceptor_prob) {
  this->splice_distance = distance;
  this->donor_prob = donor_prob;
  this->acceptor_prob = acceptor_prob;

  return;
}

void
Junction_set_ambiguous (T this) {
  this->splice_distance = 0;

  return;
}


