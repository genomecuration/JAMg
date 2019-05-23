static char rcsid[] = "$Id: getline.c 215545 2018-05-31 00:32:58Z twu $";

/* #define STANDALONE 1 */

#include "getline.h"

#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include <strings.h>		/* For rindex */

#ifdef STANDALONE
#define MALLOC malloc
#define FREE free
#else
#include "mem.h"
#endif


#define BUFSIZE 1024

/* Reads arbitrarily long lines from fp.  Strips the '\n' character from the end */
char *
Getline (FILE *fp) {
  size_t size, length;
  char *buffer, *ptr;

  size = BUFSIZE;
  buffer = (char *) MALLOC(BUFSIZE*sizeof(char));
  if (fgets(buffer,BUFSIZE,fp) == NULL) {
    FREE(buffer);
    return (char *) NULL;

  } else {
    length = strlen(buffer);
    while (!feof(fp) && buffer[length-1] != '\n') {
      size += BUFSIZE;

      ptr = buffer;
      buffer = (char *) MALLOC(size*sizeof(char));
      memcpy(buffer,ptr,length*sizeof(char));
      FREE(ptr);

      ptr = fgets(&(buffer[length]),BUFSIZE,fp);
      length += strlen(ptr);
    }

    if (buffer[length-1] == '\n') {
      buffer[length-1] = '\0';
      length -= 1;
    }

    return buffer;
  }
}

/* length does not include the "\n" character */
char *
Getline_wlength (int *string_length, FILE *fp) {
  size_t size, length;
  char *buffer, *ptr;

  size = BUFSIZE;
  buffer = (char *) MALLOC(BUFSIZE*sizeof(char));
  if (fgets(buffer,BUFSIZE,fp) == NULL) {
    FREE(buffer);
    *string_length = 0;
    return (char *) NULL;

  } else {
    length = strlen(buffer);
    while (!feof(fp) && buffer[length-1] != '\n') {
      size += BUFSIZE;

      ptr = buffer;
      buffer = (char *) MALLOC(size*sizeof(char));
      memcpy(buffer,ptr,length*sizeof(char));
      FREE(ptr);

      ptr = fgets(&(buffer[length]),BUFSIZE,fp);
      length += strlen(ptr);
    }

    if (buffer[length-1] == '\n') {
      buffer[length-1] = '\0';
      length -= 1;
    }

    *string_length = length;
    return buffer;
  }
}


/* Reads arbitrarily long lines from fp.  Does not strip the '\n'
   character from the end.  Can be used if we want to reproduce the
   input file exactly.  line_length includes the '\n' */
char *
Getline_wlinefeed (int *string_length, FILE *fp) {
  size_t size, length;
  char *buffer, *ptr;

  size = BUFSIZE;
  buffer = (char *) MALLOC(BUFSIZE*sizeof(char));
  if (fgets(buffer,BUFSIZE,fp) == NULL) {
    FREE(buffer);
    return (char *) NULL;

  } else {
    length = strlen(buffer);
    while (!feof(fp) && buffer[length-1] != '\n') {
      size += BUFSIZE;

      ptr = buffer;
      buffer = (char *) MALLOC(size*sizeof(char));
      memcpy(buffer,ptr,length*sizeof(char));
      FREE(ptr);

      ptr = fgets(&(buffer[length]),BUFSIZE,fp);
      length += strlen(ptr);
    }

#if 0
    if (buffer[length-1] == '\n') {
      buffer[length-1] = '\0';
      length -= 1;
    }
#endif

    *string_length = length;
    return buffer;
  }
}


#ifdef STANDALONE
int
main (int argc, char *argv[]) {
  FILE *fp;
  char *line;
  int length;

  fp = fopen(argv[1],"r");
  while ((line = Getline_wlength(&length,fp)) != NULL) {
    printf("%d: %s\n",length,line);
    FREE(line);
  }
  fclose(fp);

  return 0;
}
#endif
