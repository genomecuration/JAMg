static char rcsid[] = "$Id: uinttableuint.c 210071 2017-09-23 00:40:40Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "uinttableuint.h"
#include <stdio.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>		/* For qsort */
#include <string.h>		/* For strcmp */
#include "mem.h"
#include "assert.h"

#define T Uinttableuint_T
struct T {
  int size;
  int length;
  unsigned int timestamp;
  struct binding {
    struct binding *link;
    unsigned int key;
    unsigned int value;
    unsigned int timeindex;
  } **buckets;
};



T 
Uinttableuint_new (int hint) {
  T table;
  int i;
  static int primes[] = { 509, 509, 1021, 2053, 4093,
			  8191, 16381, 32771, 65521, INT_MAX };

  assert(hint >= 0);
  for (i = 1; primes[i] < hint; i++) {
  }
  table = (T) MALLOC(sizeof(*table) +
		     primes[i-1]*sizeof(table->buckets[0]));
  table->size = primes[i-1];
  table->buckets = (struct binding **)(table + 1);
  for (i = 0; i < table->size; i++) {
    table->buckets[i] = NULL;
  }
  table->length = 0;
  table->timestamp = 0;
  return table;
}

unsigned int
Uinttableuint_get (T table, const unsigned int key) {
  int i;
  struct binding *p;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  i = key % table->size;
  /* printf("Doing Uinttableuint_get on %s at bucket %d\n",(char *) key, i); */
  for (p = table->buckets[i]; p; p = p->link) {
    /* printf("  Comparing %s with %s at %p, key = %p\n",(char *) key, (char *) p->key, p, p->key); */
    if (key == p->key) {
      break;
    }
  }
  return p ? p->value : 0;
}

unsigned int
Uinttableuint_put (T table, const unsigned int key, unsigned int value) {
  int i;
  struct binding *p;
  unsigned int prev;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  i = key % table->size;
  for (p = table->buckets[i]; p; p = p->link) {
    if (key == p->key) {
      break;
    }
  }
  if (p == NULL) {
    NEW(p);
    p->key = key;
    /* printf("Doing Uinttable_put at %p, key = %p\n",p,p->key); */
    p->link = table->buckets[i];
    table->buckets[i] = p;
    table->length++;
    prev = 0;
  } else {
    prev = p->value;
  }
  p->value = value;
  p->timeindex = table->timestamp;
  table->timestamp++;
  return prev;
}

int 
Uinttableuint_length (T table) {
  assert(table);
  return table->length;
}

void 
Uinttableuint_map (T table,
	       void (*apply)(const unsigned int key, unsigned int *value, void *cl),
	       void *cl) {
  int i;
  struct binding *p;

  assert(table);
  assert(apply);
  for (i = 0; i < table->size; i++)
    for (p = table->buckets[i]; p; p = p->link) {
      apply(p->key, &p->value, cl);
    }
}

unsigned int
Uinttableuint_remove (T table, const unsigned int key) {
  int i;
  struct binding **pp;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  table->timestamp++;
  i = key % table->size;
  for (pp = &table->buckets[i]; *pp; pp = &(*pp)->link) {
    if (key == (*pp)->key) {
      struct binding *p = *pp;
      unsigned int value = p->value;
      *pp = p->link;
      FREE(p);
      table->length--;
      return value;
    }
  }
  return 0;
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
Uinttableuint_keys (T table, bool sortp) {
  unsigned int *keyarray;
  int i, j = 0;
  struct binding *p;

  assert(table);
  keyarray = (unsigned int *) CALLOC(table->length+1,sizeof(unsigned int));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      keyarray[j++] = p->key;
    }
  }

  if (sortp == true) {
    qsort(keyarray,table->length,sizeof(unsigned int),uint_compare);
  }

  return keyarray;
}


static int
timeindex_cmp (const void *x, const void *y) {
  struct binding *a = * (struct binding **) x;
  struct binding *b = * (struct binding **) y;

  if (a->timeindex < b->timeindex) {
    return -1;
  } else if (a->timeindex > b->timeindex) {
    return +1;
  } else {
    return 0;
  }
}


unsigned int *
Uinttableuint_keys_by_timeindex (T table) {
  unsigned int *keyarray;
  int i, j = 0;
  struct binding **buckets, *p;

  assert(table);
  buckets = (struct binding **) CALLOC(table->length+1,sizeof(struct binding *));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      buckets[j++] = p;
    }
  }
  qsort(buckets,table->length,sizeof(struct binding *),timeindex_cmp);

  keyarray = (unsigned int *) CALLOC(table->length,sizeof(unsigned int));
  for (j = 0; j < table->length; j++) {
    p = buckets[j];
    keyarray[j] = p->key;
  }
  FREE(buckets);

  return keyarray;
}


unsigned int *
Uinttableuint_values (T table) {
  unsigned int *valuearray;
  int i, j = 0;
  struct binding *p;

  assert(table);
  valuearray = (unsigned int *) CALLOC(table->length,sizeof(unsigned int));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      valuearray[j++] = p->value;
    }
  }
  return valuearray;
}

void 
Uinttableuint_free (T *table) {
  assert(table && *table);
  if ((*table)->length > 0) {
    int i;
    struct binding *p, *q;
    for (i = 0; i < (*table)->size; i++) {
      for (p = (*table)->buckets[i]; p; p = q) {
	q = p->link;
	FREE(p);
      }
    }
  }
  FREE(*table);
}
