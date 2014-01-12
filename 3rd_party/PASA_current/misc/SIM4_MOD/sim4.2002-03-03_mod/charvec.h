/* genvec charvec char ; 1999-10-05 22:59:09 */
#ifndef HAS_GEN_charvec_H
#define HAS_GEN_charvec_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct charvec {
    char *a;
    unsigned int len;
    unsigned int max;
    void *((*alloc)(void*, size_t));
    void (*free)(void*);
} charvec_t;

#define charvec_INIT(a,f)  {0, 0, 0, a, f}

charvec_t* charvec_new(void* ((*)(void*, size_t)), void ((*)(void*)));
charvec_t* charvec_free(charvec_t *t);
int charvec_init(charvec_t *t, void* ((*)(void*,size_t)), void((*)(void*)));
int charvec_fini(charvec_t *t);
int charvec_need(charvec_t *t, unsigned int n);
int charvec_more(charvec_t *t, unsigned int n);
int charvec_append(charvec_t *t, char e);
int charvec_fit(charvec_t *t);

#ifndef GENVEC_INBOUNDS
#define GENVEC_INBOUNDS(t,n) ((0<=(n))&&((n)<(t)->len))
#endif

#ifndef GENVEC_GET
#define GENVEC_GET(t,n) (assert(GENVEC_INBOUNDS(t,n)) , (t)->a[n])
#endif

#endif

