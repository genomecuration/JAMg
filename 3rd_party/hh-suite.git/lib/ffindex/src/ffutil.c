/*
 * Ffindex
 * written by Andy Hauser <hauser@genzentrum.lmu.de>.
 * Please add your name here if you distribute modified versions.
 * 
 * Ffindex is provided under the Create Commons license "Attribution-ShareAlike
 * 3.0", which basically captures the spirit of the Gnu Public License (GPL).
 * 
 * See:
 * http://creativecommons.org/licenses/by-sa/3.0/
 * 
 * Ffindex is a very simple database for small files. The files are stored
 * concatenated in one big data file, seperated by '\0'. A second file
 * contains a plain text index, giving name, offset and length of of the small
 * files.
 */

#include "ffutil.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

int fferror_print(char *sourcecode_filename, int line, const char *function_name, const char *message)
{
  int myerrno = errno;
  char* errstr = strerror(myerrno);
  fprintf(stderr, "%s:%d %s: %s: %s\n", sourcecode_filename , line, function_name, message, errstr);
  return myerrno;
}


/* remove \n, assumes UNIX line endings! */
char* ffnchomp(char *s, size_t len)
{
  // prevent underflow because of the size_t
  // we want to chomp off the last element
  len = len > 1 ? len - 1 : 0;
  if(s[len] == '\n')
    s[len] = '\0';

  return s;
}

size_t ffcount_lines(const char *filename)
{
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    return 0;
  }

  size_t lines = 0;
  while (!feof(fp)) {
    char ch = (char) fgetc(fp);
    if (ch == '\n') {
      lines++;
    }
  }

  fclose(fp);

  return lines;
}

void ffmerge_splits(const char* data_filename, const char* index_filename, int splits, int remove_temporary)
{
  if (!data_filename)
    return;

  if (!index_filename)
    return;

  char merge_command[FILENAME_MAX * 5];
  char tmp_filename[FILENAME_MAX];

  for (int i = 1; i < splits; i++)
  {
    snprintf(merge_command, FILENAME_MAX,
             "ffindex_build -as -d %s.%d -i %s.%d %s %s",
             data_filename, i, index_filename, i, data_filename, index_filename);

    int ret = system(merge_command);
    if (ret == 0 && remove_temporary)
    {
      snprintf(tmp_filename, FILENAME_MAX, "%s.%d", index_filename, i);
      remove(tmp_filename);

      snprintf(tmp_filename, FILENAME_MAX, "%s.%d", data_filename, i);
      remove(tmp_filename);
    }
  }
}

/* vim: ts=2 sw=2 et
*/
