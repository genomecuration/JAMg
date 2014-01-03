// Copyright (c) 2003  by  Mihaela Pertea.

#ifndef _MISC_H
#define _MISC_H

#define TRUE 1
#define FALSE 0


/* Function prototypes */
extern void * check_calloc (size_t, size_t);
extern FILE * check_fopen (const char *, const char *);
extern void * check_malloc (size_t);
extern void * check_realloc (void *, size_t);

#endif
