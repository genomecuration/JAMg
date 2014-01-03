/*
 * Ffindex
 * written by Andreas Hauser <andy@splashground.de>.
 * Please add your name here if you distribute modified versions.
 * 
 * Ffindex is provided under the Create Commons license "Attribution-ShareAlike
 * 3.0", which basically captures the spirit of the Gnu Public License (GPL).
 * 
 * See:
 * http://creativecommons.org/licenses/by-sa/3.0/
*/

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include "ffindex.h"
#include "ffutil.h"

#define MAX_FILENAME_LIST_FILES 4096


void usage(char *program_name)
{
    fprintf(stderr, "USAGE: %s OUT_DATA_FILE OUT_INDEX_FILE -d IN_DATA_FILE -i IN_INDEX_FILE\n"
                    "\t-d FFDATA_FILE\tinput ffindex data file to be optimized\n"
                    "\t-i FFINDEX_FILE\tinput ffindex index file to be optimized\n"
                    "\t-v\t\tprint version and other info then exit\n"
                    FFINDEX_COPYRIGHT,
                    program_name, MAX_FILENAME_LIST_FILES, FFINDEX_MAX_ENTRY_NAME_LENTH, FFINDEX_MAX_INDEX_ENTRIES_DEFAULT);
}

int main(int argn, char **argv)
{
  int version = 0;
  int opt, err = EXIT_SUCCESS;
  char* list_ffindex_data[MAX_FILENAME_LIST_FILES];
  char* list_ffindex_index[MAX_FILENAME_LIST_FILES];
  size_t list_ffindex_data_index = 0;
  size_t list_ffindex_index_index = 0;
  size_t list_filenames_index = 0;
  while ((opt = getopt(argn, argv, "vd:i:")) != -1)
  {
    switch (opt)
    {
      case 'd':
        list_ffindex_data[list_ffindex_data_index++] = optarg;
        break;
      case 'i':
        list_ffindex_index[list_ffindex_index_index++] = optarg;
        break;
      case 'v':
        version = 1;
        break;
      default:
        usage(argv[0]);
        return EXIT_FAILURE;
    }
  }

  if(version == 1)
  {
    /* Don't you dare running it on a platform where byte != 8 bits */
    printf("%s version %.2f, off_t = %zd bits\n", argv[0], FFINDEX_VERSION, sizeof(off_t) * 8);
    return EXIT_SUCCESS;
  }

  if(argn - optind < 2)
  {
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  if(list_ffindex_data_index != list_ffindex_index_index)
  {
    fprintf(stderr, "ERROR: -d and -i must be specified pairwise\n");
    return EXIT_FAILURE;
  }

  char *data_filename  = argv[optind++];
  char *index_filename = argv[optind++];
  FILE *data_file = NULL, *index_file = NULL;

  size_t offset = 0;

  /* open index and data file, seek to end if needed */
  char *mode = "w";
  if(append)
    mode = "a";
  err = ffindex_index_open(data_filename, index_filename, mode, &data_file, &index_file, &offset);
  if(err != EXIT_SUCCESS)
    return EXIT_FAILURE;

  /* For each list_file insert */
  if(list_filenames_index > 0)
    for(int i = 0; i < list_filenames_index; i++)
    {
      FILE *list_file = fopen(list_filenames[i], "r");
      if( list_file == NULL) { perror(list_filenames[i]); return EXIT_FAILURE; }
      if(ffindex_insert_list_file(data_file, index_file, &offset, list_file) < 0)
      {
        perror(list_filenames[i]);
        err = -1;
      }
    }

  /* Append other ffindexes */
  if(list_ffindex_data_index > 0)
  {
    for(int i = 0; i < list_ffindex_data_index; i++)
    {
      FILE* data_file_to_add  = fopen(list_ffindex_data[i], "r");  if(  data_file_to_add == NULL) { perror(list_ffindex_data[i]); return EXIT_FAILURE; }
      FILE* index_file_to_add = fopen(list_ffindex_index[i], "r"); if( index_file_to_add == NULL) { perror(list_ffindex_index[i]); return EXIT_FAILURE; }
      size_t data_size;
      char *data_to_add = ffindex_mmap_data(data_file_to_add, &data_size);
      ffindex_index_t* index_to_add = ffindex_index_parse(index_file_to_add, 0);
      ffindex_insert_ffindex(data_file, index_file, &offset, data_to_add, index_to_add);
    }
  }


  /* Insert files and directories into the index */
  for(int i = optind; i < argn; i++)
  {
    char *path = argv[i];
    struct stat sb;
    if(stat(path, &sb) == -1)
    {
      fferror_print(__FILE__, __LINE__, __func__, path);
      continue;
    }

    if(S_ISDIR(sb.st_mode))
    {
      err = ffindex_insert_dir(data_file, index_file, &offset, path);
      if(err < 0)fferror_print(__FILE__, __LINE__, __func__, path);
    }
    else if(S_ISREG(sb.st_mode))
    {
      ffindex_insert_file(data_file, index_file, &offset, path, path);
    }
  }
  fclose(data_file);

  /* Sort the index entries and write back */
  if(sort)
  {
    rewind(index_file);
    ffindex_index_t* index = ffindex_index_parse(index_file, 0);
    if(index == NULL)
    {
      fferror_print(__FILE__, __LINE__, __func__, index_filename);
      exit(EXIT_FAILURE);
    }
    fclose(index_file);
    ffindex_sort_index_file(index);
    index_file = fopen(index_filename, "w");
    if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }
    err += ffindex_write(index, index_file);
  }

  return err;
}

/* vim: ts=2 sw=2 et
 */
