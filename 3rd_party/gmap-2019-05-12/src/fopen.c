static char rcsid[] = "$Id: fopen.c 214794 2018-04-21 00:54:37Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For USE_FOPEN_BINARY, USE_FOPEN_TEXT */
#endif

#include "fopen.h"
#include <stdlib.h>
#include <string.h>
#include "mem.h"


FILE *
Fopen_read_text (char *read_files_command, char *filename) {
  FILE *fp;
  char *command;

  if (read_files_command != NULL) {
    command = (char *) MALLOC((strlen(read_files_command) + strlen(" ") + strlen(filename) + 1) * sizeof(char));
    sprintf(command,"%s %s",read_files_command,filename);

    /* Note: popen does not take "rt" as fopen might */
    if ((fp = popen(command,"r")) == NULL) {
      fprintf(stderr,"Cannot open file %s with command %s\n",filename,read_files_command);
    }
    FREE(command);

  } else {
#if USE_FOPEN_TEXT
    if ((fp = fopen(filename,"rt")) == NULL) {
      fprintf(stderr,"Cannot open file %s\n",filename);
    }
#else
    if ((fp = fopen(filename,"r")) == NULL) {
      fprintf(stderr,"Cannot open file %s\n",filename);
    }
#endif
  }

  return fp;
}

