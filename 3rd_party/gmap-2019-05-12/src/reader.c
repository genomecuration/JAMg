static char rcsid[] = "$Id: reader.c 218153 2019-01-17 05:38:29Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "reader.h"
#include "mem.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


#define T Reader_T
struct T {
  int querystart;
  int queryend;
  int querystart_save;
  int queryend_save;

  char *startinit;
  char *startptr;
  char *endptr;

  char *startbound;		/* Saved for reset */
  char *endbound;
};

int
Reader_querystart (T this) {
  return this->querystart;
}

int
Reader_queryend (T this) {
  return this->queryend;
}

/* Same as current querypos + oligosize */
int
Reader_startpos (T this) {
  return (this->startptr - this->startinit);
}

int
Reader_endpos (T this) {
  return (this->endptr - this->startinit);
}


void
Reader_reset_start (T this, int querypos) {
  char *sequence;

  sequence = this->startinit;
  this->startptr = &(sequence[querypos]);
  return;
}

void
Reader_reset_end (T this, int querypos) {
  char *sequence;

  sequence = this->startinit;
  this->endptr = &(sequence[querypos]);
  return;
}


void
Reader_reset_ends (T this) {
  /*
  char *sequence;
  sequence = this->startinit;
  this->startptr = &(sequence[this->querystart]);
  this->endptr = &(sequence[this->queryend-1]);
  */

  this->querystart = this->querystart_save;
  this->queryend = this->queryend_save;

  this->startptr = this->startbound;
  this->endptr = this->endbound;

  return;
}



T
Reader_new (char *sequence, int querystart, int queryend) {
  T new = (T) MALLOC(sizeof(*new));

  new->querystart_save = new->querystart = querystart;
  new->queryend_save = new->queryend = queryend;

  new->startinit = sequence;
  new->startbound = new->startptr = &(sequence[querystart]);
  new->endbound = new->endptr = &(sequence[queryend-1]);
  return new;
}

void
Reader_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

#ifndef GSNAP
char
Reader_getc (T this, cDNAEnd_T cdnaend, int blocksize) {
  debug(fprintf(stderr,"Read_getc has startptr %d and endptr %d\n",
		this->startptr-this->startinit,this->endptr-this->startinit));
  if (this->startptr - this->endptr >= blocksize) {
    return '\0';
  } else if (cdnaend == FIVE) {
    return *(this->startptr++);
  } else { 
    return *(this->endptr--);
  }
}
#endif

char
Reader_getc_5 (T this) {
  debug(printf("Read_getc has startptr %d and endptr %d\n",
	       this->startptr - this->startinit,this->endptr - this->startinit));
  if (this->startptr > this->endbound) {
    return '\0';
  } else {
    return *(this->startptr++);
  }
}

char
Reader_getc_3 (T this) {
  debug(printf("Read_getc has startptr %d and endptr %d\n",
	       this->startptr - this->startinit,this->endptr - this->startinit));
  if (this->endptr < this->startbound) {
    return '\0';
  } else {
    return *(this->endptr--);
  }
}


#if 0
void
Reader_set_5 (T this, int querypos, int oligosize) {
  this->startptr = this->startinit + querypos;
  debug(printf("Reader_set now has startpos %d\n",this->startptr - this->startinit));
  return;
}
#endif

#if 0
void
Reader_set_3 (T this, int querypos, int oligosize) {
  this->endptr = this->startinit + querypos + oligosize - 1;
  debug(printf("Reader_set now has endpos %d\n",this->endptr - this->startinit));
  return;
}
#endif

#if 0
void
Reader_skip_5 (T this, int nskip) {
  this->startptr += nskip;
  return;
}
#endif

#if 0
void
Reader_skip_3 (T this, int nskip) {
  this->endptr -= nskip;
  return;
}
#endif


#if 0
/* For debugging */
static Oligospace_T
nt_oligo (char *query, int indexsize) {
  Oligospace_T oligo = 0U;
  int i;

  debug9(printf("Reader_check: "));
  for (i = 0; i < indexsize; i++) {
    oligo *= 4;
    
    debug9(printf("%c",query[i]));
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
  debug9(printf(" => oligo %016lX\n",oligo));

  return oligo;
}

Oligospace_T
Reader_check (T this, int querypos, int indexsize) {

  debug9(printf("Read: %s\n",this->startinit));
  return nt_oligo(&(this->startinit[querypos]),indexsize);
}
#endif


/* For testing */
/*
static void
process_input (FILE *input) {
  bool ssfile;
  int queryno = 0, seqlength, skiplength;
  char *query, initc;

  if ((initc = Reader_input_init(input)) == '>') {
    ssfile = false;
    query = Reader_input_header(input);
  } else {
    ssfile = true;
    query = (char *) CALLOC(strlen("NO_HEADER")+1,sizeof(char));
    strcpy(query,"NO_HEADER");
  }
  fprintf(stderr,"Read sequence %s\n",query);
  Reader_input_sequence(&seqlength,&skiplength,input,initc);
  print_contents(&(Sequence[0]),SEQUENCELEN);
  fprintf(stderr,"\nSeqlength = %d\n",seqlength);
  FREE(query);
  queryno++;

  while ((query = Reader_input_header(input)) != NULL) {
    fprintf(stderr,"Read sequence %s\n",query);
    Reader_input_sequence(&seqlength,&skiplength,input,initc);
    print_contents(&(Sequence[0]),SEQUENCELEN);
    fprintf(stderr,"\nSeqlength = %d\n",seqlength);
    FREE(query);
    queryno++;
  }

  return;
}

int
main (int argc, char *argv[]) {
  process_input(stdin);
  return 0;
}
*/
