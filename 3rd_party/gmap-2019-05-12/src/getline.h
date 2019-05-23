#ifndef GETLINE_INCLUDED
#define GETLINE_INCLUDED

#include <stdio.h>

extern char *
Getline (FILE *fp);
extern char *
Getline_wlength (int *string_length, FILE *fp);
extern char *
Getline_wlinefeed (int *string_length, FILE *fp);

#endif


