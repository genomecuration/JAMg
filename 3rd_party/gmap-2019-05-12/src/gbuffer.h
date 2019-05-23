/* $Id: gbuffer.h 218286 2019-01-23 16:46:55Z twu $ */
#ifndef GBUFFER_INCLUDED
#define GBUFFER_INCLUDED

#define T Gbuffer_T
typedef struct T *T;


extern int
Gbuffer_gbufferlen (T this);
extern char *
Gbuffer_chars1 (T this);
extern char *
Gbuffer_chars2 (T this);
extern char *
Gbuffer_chars3 (T this);

extern void
Gbuffer_free_contents (T this);
extern void
Gbuffer_alloc_contents (T this, int gbufferlen);

extern void
Gbuffer_free (T *old);
extern T
Gbuffer_new ();


#undef T
#endif

