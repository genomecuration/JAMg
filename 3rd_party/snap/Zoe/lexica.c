/*****************************************************************************\
 lexica.c

After "Hello World", I believe that lexical analysis is the next best program
to write. It often shows off the strengths and weakness of a language as well
as the canonical loops for dynamic lists and maps.

Copyright (C) 2002-2013 Ian Korf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

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
