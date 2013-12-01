/*****************************************************************************\
 lexica.c

After "Hello World", I believe that lexical analysis is the next best program
to write. It often shows off the strengths and weakness of a language as well
as the canonical loops for dynamic lists and maps.

Copyright (C) Ian Korf 2002-2013.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The MIT License (MIT) - opensource.org/licenses/MIT

\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zoe.h"

/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char **argv) {
	FILE    *file = NULL;
	char     line[8192], *token = NULL;
	int      i, len, string_start, string_found;
	int     *iptr = NULL;
	void    *val  = NULL;
	zoeHash hash = NULL;
	zoeTVec keys = NULL;
	
	zoeSetProgramName(argv[0]);
	zoeSetOption("-table", 0);
	zoeSetOption("-stats", 0);
	zoeParseOptions(&argc, argv);
	
	/* commandline */
	if (argc != 2) {
		zoeS(stderr, "usage: %s [-table -stats] <file>\n", argv[0]);
		exit(1);
	};
	
	/* create data structures */
	hash = zoeNewHash();
	
	/* open file */
	if ((file = fopen(argv[1], "r")) == NULL) zoeExit("file error (%s)", argv[1]);
	
	/* parse file into tokens, store unique tokens */
	while (fgets(line, sizeof(line), file) != NULL) {
		len = strlen(line);
		string_start = 0;
		string_found = 0;
		for (i = 0; i < len; i++) {
			if (isalnum((int)line[i])) {
				if (string_found) continue;
				string_found = 1;
				string_start = i;
			} else {
				if (!string_found) continue;
				string_found = 0;
				line[i] = '\0';
				token = line + string_start;
				
				val = zoeGetHash(hash, token);
				if (val == NULL) {
					iptr = malloc(sizeof(int));
					*iptr = 1;
					zoeSetHash(hash, token, iptr);
				} else {
					iptr = val;
					(*iptr)++;   /* uglier: (*(int*)val)++; */
				}
				
			}
		}
	}
	(void)fclose(file);
	
	
	/* print out word usage */
	if (zoeOption("-table")) {
		keys = zoeKeysOfHash(hash);
		qsort(keys->elem, keys->size, sizeof(void*), zoeTcmp);
		for (i = 0; i < keys->size; i++) {
			zoeS(stdout, "%s\t%d\n",
				(char*)keys->elem[i],
				*(int*)zoeGetHash(hash, (char*)keys->elem[i])
			);
		}
		zoeDeleteTVec(keys);
	}

	if (zoeOption("-stats"))  zoeStatHash(hash);
	
	/* clean house */
	zoeDeleteHash(hash);
	
	return 0;
}
