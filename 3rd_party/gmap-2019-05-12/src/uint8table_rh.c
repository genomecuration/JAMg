static char rcsid[] = "$Id: uint8table_rh.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "uint8table_rh.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>		/* For memset */
#include <stdlib.h>		/* For qsort */

#include "mem.h"
#include "assert.h"
#include "list.h"
#include "types.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#ifdef HAVE_64_BIT		/* Defined in types.h, included in uint8table.h */

/* Implementation of a hash table using open addressing and Robin Hood probing */
/* Key may not be -1U, because that value is used to indicate empty */

#define EMPTY (UINT8) -1	/* Needed so we can use 0 as a key */

/* Choose powers of 2 to have a fast hash function */
static UINT8 powers_of_2[] = { 2UL, 4UL, 8UL, 16UL, 32UL, 64UL, 128UL, 256UL, 512UL, 1024UL,
			       2048UL, 4096UL, 8192UL, 16384UL, 32768UL, 65536UL, 131072UL,
			       262144UL, 524288UL, 1048576UL, 2097152UL, 4194304UL, 8388608UL,
			       16777216UL, 33554432UL, 67108864UL, 134217728UL,
			       268435456UL, 536870912UL, 1073741824UL, 2147483648UL};
#define LOAD_FACTOR 3

#define T Uint8table_T
struct T {
  int size;
  int longest_dist;

  int nbuckets;
  UINT8 mask;			/* hash by taking the lower order bits */

  UINT8 *keys;
  void **contents;
  void **saved_contents;	/* Stored in order, to make Uint8table_values faster */
};


#if 0
/* For arbitrary nbuckets */
static inline int
reduce (UINT8 x, int nbuckets) {
  return x % nbuckets;
}
#else

/* For powers of 2 */
static inline int
reduce (UINT8 x, UINT8 mask) {
  return x & mask;
}
#endif


/* n is the number of total elements, although some may be repeated */
T 
Uint8table_new (int hint, bool save_contents_p) {
  T new = (T) MALLOC(sizeof(*new));
  int level;

  for (level = 0; powers_of_2[level] < (UINT8) (hint * LOAD_FACTOR); level++) {
  }
  new->nbuckets = powers_of_2[level];
  new->mask = (UINT8) (new->nbuckets - 1);
  debug(printf("Creating table with %d buckets\n",new->nbuckets));

  new->size = 0;
  new->longest_dist = 0;
  new->keys = (UINT8 *) MALLOC(new->nbuckets * sizeof(UINT8));
  memset(new->keys,0xff,new->nbuckets * sizeof(UINT8));

  new->contents = (void **) MALLOC(new->nbuckets * sizeof(void *));
  if (save_contents_p == false) {
    new->saved_contents = (void **) NULL;
  } else {
    new->saved_contents = (void **) MALLOC(new->nbuckets * sizeof(void *));
  }

  return new;
}
  
int
Uint8table_length (T this) {
  return this->size;
}


void *
Uint8table_get (T this, const UINT8 key) {
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
insert_aux (T this, UINT8 probe_key, void *probe_contents) {
  int probe_dist = 0, occupant_dist;
  int bucketi, occupant_bucketi;
  UINT8 occupant_key;
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
  UINT8 *old_keys = this->keys;
  void **old_contents = this->contents, **old_saved_contents;

  this->nbuckets *= 2;
  this->mask = (UINT8) (this->nbuckets - 1);
  debug(printf("Growing hash table to %d buckets\n",this->nbuckets));

  this->keys = (UINT8 *) MALLOC(this->nbuckets * sizeof(UINT8));
  memset(this->keys,0xff,this->nbuckets * sizeof(UINT8));

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
Uint8table_put (T this, UINT8 key, void *contents) {
  assert(key != EMPTY);

  if (this->size == this->nbuckets) {
    grow(this);
  }
  insert_aux(this,key,contents);
  /* this->saved_contents[this->size] = contents; */
  this->size += 1;

  return;
}

/* Valid only when caller and performed Uinttable_put_and_save */
void
Uint8table_put_and_save (T this, UINT8 key, void *contents) {
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
Uint8table_saved_values (int *nvalues, T this) {

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
uint8_compare (const void *a, const void *b) {
  UINT8 x = * (UINT8 *) a;
  UINT8 y = * (UINT8 *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}


UINT8 *
Uint8table_keys (T this, bool sortp) {
  UINT8 *keyarray, key;
  int i, k = 0;

  keyarray = (UINT8 *) CALLOC(this->size+1,sizeof(UINT8));
  for (i = 0; i < this->nbuckets; i++) {
    if ((key = this->keys[i]) != EMPTY) {
      keyarray[k++] = key;
    }
  }

  if (sortp == true) {
    qsort(keyarray,this->size,sizeof(UINT8),uint8_compare);
  }

  return keyarray;
}

void 
Uint8table_free (T *old) {
  debug(printf("Freeing table\n"));
  FREE((*old)->saved_contents);
  FREE((*old)->contents);
  FREE((*old)->keys);
  FREE(*old);
  return;
}

#endif

