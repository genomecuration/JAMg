static char rcsid[] = "$Id: uinttable_rh.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* #define STANDALONE 1 */

#include "uinttable_rh.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>		/* For memset */
#include <stdlib.h>		/* For qsort */

#ifdef STANDALONE
#include <stdlib.h>
#define CALLOC calloc
#define MALLOC malloc
#define FREE free
#define HAVE_64_BIT 1
#define UINT8 unsigned long long
#else
#include "mem.h"
#include "assert.h"
#include "list.h"
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Implementation of a hash table using open addressing and Robin Hood probing */
/* Key may not be -1U, because that value is used to indicate empty */

#define EMPTY (unsigned int) -1	/* Needed so we can use 0U as a key */

/* Choose powers of 2 to have a fast hash function */
static unsigned int powers_of_2[] = { 2U, 4U, 8U, 16U, 32U, 64U, 128U, 256U, 512U, 1024U,
				      2048U, 4096U, 8192U, 16384U, 32768U, 65536U, 131072U,
				      262144U, 524288U, 1048576U, 2097152U, 4194304U, 8388608U,
				      16777216U, 33554432U, 67108864U, 134217728U,
				      268435456U, 536870912U, 1073741824U, 2147483648U};
#define LOAD_FACTOR 3

#define T Uinttable_T
struct T {
  int size;
  int longest_dist;

  int nbuckets;
  unsigned int mask;			/* hash by taking the lower order bits */

  unsigned int *keys;
  void **contents;
  void **saved_contents;	/* Stored in order, to make Uinttable_values faster */
};


#if 0
/* For arbitrary nbuckets */
static inline int
reduce (unsigned int x, int nbuckets) {
  return x % nbuckets;
}
#else

/* For powers of 2 */
static inline int
reduce (unsigned int x, unsigned int mask) {
  return x & mask;
}
#endif


/* n is the number of total elements, although some may be repeated */
T 
Uinttable_new (int hint, bool save_contents_p) {
  T new = (T) MALLOC(sizeof(*new));
  int level;

  for (level = 0; powers_of_2[level] < (unsigned int) (hint * LOAD_FACTOR); level++) {
  }
  new->nbuckets = powers_of_2[level];
  new->mask = (unsigned int) (new->nbuckets - 1);
  debug(printf("Creating table with %d buckets\n",new->nbuckets));

  new->size = 0;
  new->longest_dist = 0;
  new->keys = (unsigned int *) MALLOC(new->nbuckets * sizeof(unsigned int));
  memset(new->keys,0xff,new->nbuckets * sizeof(unsigned int));

  new->contents = (void **) MALLOC(new->nbuckets * sizeof(void *));
  if (save_contents_p == false) {
    new->saved_contents = (void **) NULL;
  } else {
    new->saved_contents = (void **) MALLOC(new->nbuckets * sizeof(void *));
  }

  return new;
}
  
int
Uinttable_length (T this) {
  return this->size;
}


void *
Uinttable_get (T this, const unsigned int key) {
  int probe_dist = 0;
  int bucketi;

  bucketi = reduce(key,this->mask);
  while (probe_dist <= this->longest_dist && this->keys[bucketi] != EMPTY) {
    if (this->keys[bucketi] == key) {
      return this->contents[bucketi];
    } else {
      probe_dist++;
      if (++bucketi >= this->nbuckets) {
	bucketi = 0;
      }
    }
  }

  /* Unsuccessful search */
  return (void *) NULL;
}


static void
insert_aux (T this, unsigned int probe_key, void *probe_contents) {
  int probe_dist = 0, occupant_dist;
  int bucketi, occupant_bucketi;
  unsigned int occupant_key;
  void *occupant_contents;

  debug(printf("Entered insert_aux with key %u, and contents %p\n",probe_key,probe_contents));
  bucketi = reduce(probe_key,this->mask);
  while ((occupant_key = this->keys[bucketi]) != EMPTY) {
    occupant_bucketi = reduce(occupant_key,this->mask);
    if ((occupant_dist = bucketi - occupant_bucketi) < 0) {
      occupant_dist += this->nbuckets;
    }

    if (occupant_dist < probe_dist) {
      /* Put probe here and make the occupant move */
      debug(printf("Putting key %u and contents %p into bucketi\n",
		   probe_key,probe_contents,bucketi));

      occupant_contents = this->contents[bucketi];
      this->contents[bucketi] = probe_contents;
      probe_contents = occupant_contents;

      /* occupant_key = this->keys[bucketi]; -- obtained above */
      this->keys[bucketi] = probe_key;
      probe_key = occupant_key;

      if (probe_dist > this->longest_dist) {
	this->longest_dist = probe_dist;
      }
      probe_dist = occupant_dist;
    }

    probe_dist += 1;
    if (++bucketi >= this->nbuckets) {
      bucketi = 0;
    }
  }

  debug(printf("Putting key %u and contents %p into bucketi\n",
	       probe_key,probe_contents,bucketi));
  this->contents[bucketi] = probe_contents;
  this->keys[bucketi] = probe_key;

  if (probe_dist > this->longest_dist) {
    this->longest_dist = probe_dist;
  }

  return;
}


static void
grow (T this) {
  int old_nbuckets = this->nbuckets, i;
  unsigned int *old_keys = this->keys;
  void **old_contents = this->contents, **old_saved_contents;

  this->nbuckets *= 2;
  this->mask = (unsigned int) (this->nbuckets - 1);
  debug(printf("Growing hash table to %d buckets\n",this->nbuckets));

  this->keys = (unsigned int *) MALLOC(this->nbuckets * sizeof(unsigned int));
  memset(this->keys,0xff,this->nbuckets * sizeof(unsigned int));

  this->contents = (void **) MALLOC(this->nbuckets * sizeof(void *));
  if (this->saved_contents != NULL) {
    old_saved_contents = this->saved_contents;
    this->saved_contents = (void **) MALLOC(this->nbuckets * sizeof(void *));
    memcpy(this->saved_contents,old_saved_contents,this->size * sizeof(void *));
    FREE(old_saved_contents);
  }
  
  this->longest_dist = 0;
  for (i = 0; i < old_nbuckets; i++) {
    if (old_keys[i] != EMPTY) {
      insert_aux(this,old_keys[i],old_contents[i]);
    }
  }

  FREE(old_contents);
  FREE(old_keys);
  return;
}


void
Uinttable_put (T this, unsigned int key, void *contents) {
  assert(key != EMPTY);

  if (this->size == this->nbuckets) {
    grow(this);
  }
  insert_aux(this,key,contents);
  /* this->saved_contents[this->size] = contents; */
  this->size += 1;

  return;
}

void
Uinttable_put_and_save (T this, unsigned int key, void *contents) {
  assert(key != EMPTY);

  if (this->size == this->nbuckets) {
    grow(this);
  }
  insert_aux(this,key,contents);
  this->saved_contents[this->size] = contents;
  this->size += 1;

  return;
}


/* Valid only when caller specified save_contents_p and performs Uinttable_put_and_save */
void **
Uinttable_saved_values (int *nvalues, T this) {

#if 0
  void **valuearray;
  int i, k = 0;

  valuearray = (void **) MALLOC(this->size * sizeof(void *));
  for (i = 0; i < this->nbuckets; i++) {
    if (this->keys[i] != EMPTY) {
      valuearray[k++] = this->contents[i];
    }
  }
  *nvalues = this->size;
  return valuearray;

#else
  /* Faster, since we don't need to scan all buckets */
  *nvalues = this->size;
  return this->saved_contents;
#endif
}

static int
uint_compare (const void *a, const void *b) {
  unsigned int x = * (unsigned int *) a;
  unsigned int y = * (unsigned int *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}


unsigned int *
Uinttable_keys (T this, bool sortp) {
  unsigned int *keyarray, key;
  int i, k = 0;

  keyarray = (unsigned int *) CALLOC(this->size+1,sizeof(unsigned int));
  for (i = 0; i < this->nbuckets; i++) {
    if ((key = this->keys[i]) != EMPTY) {
      keyarray[k++] = key;
    }
  }

  if (sortp == true) {
    qsort(keyarray,this->size,sizeof(unsigned int),uint_compare);
  }

  return keyarray;
}

void 
Uinttable_free (T *old) {
  debug(printf("Freeing table\n"));
  FREE((*old)->saved_contents);
  FREE((*old)->contents);
  FREE((*old)->keys);
  FREE(*old);
  return;
}


#ifdef STANDALONE
int
main (int argc, char *argv[]) {
  T table = Uinttable_new(10);
  void *elt;
  int i;

  printf("First key: %u\n",table->keys[0]);

  for (i = 0; i < 50; i++) {
    Uinttable_put(table,i,(void *) i);
  }

  for (i = 0; i < 60; i++) {
    elt = Uinttable_get(table,i);
    printf("%d => %d\n",i,(int) elt);
  }

  Uinttable_free(&table);

  return 0;
}
#endif

