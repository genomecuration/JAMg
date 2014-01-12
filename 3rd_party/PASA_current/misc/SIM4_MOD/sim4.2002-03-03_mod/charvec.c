/* genvec charvec char ; 1999-10-05 22:59:09 */
#include "charvec.h"

charvec_t* charvec_new(void* ((*ra)(void*,size_t)), void (*fr)(void*))
{
    charvec_t *vec = ra(0, sizeof(*vec));
    if (vec) {
        if (charvec_init(vec, ra, fr))
            return vec;
        fr(vec);
    }
    return 0;
}

charvec_t* charvec_free(charvec_t *t)
{
    charvec_fini(t);
    t->free(t);
    return 0;
}

int charvec_init(charvec_t *t, void* ((*a)(void*,size_t)), void (*f)(void*))
{
    assert(t);
    t->a = 0;
    t->len = 0;
    t->max = 0;
    t->alloc = a;
    t->free = f;
    return charvec_need(t, 0);
}

int charvec_fini(charvec_t *t)
{
    assert(t);

    if (t->a && t->free) { t->free(t->a); t->a = 0; t->max = 0; }
    t->len = 0;
    return 1;
}

#ifndef BASE_ALLOC
#define BASE_ALLOC 30
#endif
enum { BASE = BASE_ALLOC };

int charvec_need(charvec_t *t, unsigned int n)
{
    assert(t);
    if (t->a == 0) {
        assert(t->alloc);
        t->len = 0; 
        t->max = n;
        t->a = t->alloc(0, n * sizeof(char));
        return t->a != 0;
    }

    if (n > t->max) { 
        unsigned int i = BASE + n + (n >> 3); 
        void *p = t->alloc(t->a, i * sizeof(char));
        if (!p)
            return 0;
        t->max = i;
        t->a = p;
    }
    return 1;
}

int charvec_more(charvec_t *t, unsigned int n)
{
    assert(t);
    return charvec_need(t, n + t->len);
}

int charvec_append(charvec_t *t, char e)
{
    assert(t);
    if (!charvec_more(t, 1))
        return 0;
    t->a[t->len++] = e;
    return 1;
}

int charvec_fit(charvec_t *t)
{
    assert(t);
    assert(t->alloc);

    {
    unsigned int i = t->len;
    void *p = t->alloc(t->a, i * sizeof(char));
    if (!p)
        return 0;
    t->max = i;
    t->a = p;
    return 1;
    }
}

