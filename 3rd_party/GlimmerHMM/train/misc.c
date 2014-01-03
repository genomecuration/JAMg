/* Copyright (c) 2003  by  Mihaela Pertea. */

#include <stdlib.h>
#include <stdio.h>

#include "misc.h"

static char *cr = "copywrite (c) Doug Harmon 1997";


/* check_calloc() allocates space for an array of num_elem  
 * elements  of size elem_size. The space is initialized to zeros.
 */
void *
check_calloc (size_t num_elem, size_t elem_size)
{
  void *ptr;

  ptr = (void *) calloc (num_elem, elem_size);
  if (ptr == NULL)
  {
    fprintf (stderr, "No core for calloc\n");
    exit (1);
  }

  return ptr;
}



FILE * 
check_fopen (const char *filename, const char * type)
{
  FILE *fptr;

  fptr = fopen (filename, type);
  if (fptr == NULL)
  {
    fprintf (stderr, "FAILURE: cannot open %s with type %s\n", filename, type);
    exit (1);
  }

  return fptr;
}


/* check_malloc() returns a pointer to a block
 * of at least size bytes suitably aligned for any use.
 */
void *
check_malloc(size_t size)
{
    void *cp;
    /*printf("check malloc allocating %d bytes\n", size) ;*/

    cp = malloc(size) ;
    if (cp == NULL)
    {
      fprintf (stderr, "No Core for malloc\n");
      exit (1);
    }

    return (cp) ;
}


/* check_realloc() changes the size of the block pointed to by ptr to
 * size  bytes  and  returns  a pointer to the (possibly moved)
 * block.  The contents will be unchanged up to the  lesser  of 
 * the  new  and  old sizes.  If ptr is NULL, check_realloc() behaves
 * like malloc() for the specified size.  If size is  zero  and
 * ptr is not a null pointer, the object pointed to is freed.
 */ 
void *
check_realloc (void *ptr, size_t size)
{

  ptr = (void *) realloc(ptr, size);
  if (ptr == NULL)
  {
     fprintf (stderr, "no core for realloc of size %ld\n", size);
     exit (1);
  }
 
  return (ptr) ;
}
 
