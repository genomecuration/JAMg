/*****************************************************************************\
 seqstats.c

Produces statistics on sequence files

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
#include "zoe.h"


int ASCII_SET = 256;

int main (int argc, char *argv[]) {
	zoeFile file;
	zoeFastaFile seq = NULL;
	double entropy = 0;
	float f;
	double d;
	int file_count = 0;
	long total_length = 0;
	int this_arg;
	int i;
	long count[256];
	long GC = 0, AT = 0;
	long longest = 0, shortest = INT_MAX;
	
	for(i=0;i<ASCII_SET;i++) count[i] = 0;
	
	/* set the program name */
	zoeSetProgramName(argv[0]);
	
	if (argc == 1) {
		zoeS(stderr, "usage: %s <fasta files>\n", argv[0]);
		exit(1);
	}
	
	for (this_arg = 1; this_arg < argc; this_arg++) {
	
		/* open the file and get the sequences */
		file = zoeOpenFile(argv[this_arg]);
		while ((seq = zoeReadFastaFile(file.stream)) != NULL) {
			file_count++;
			
			/* count the letters */
			for (i = 0; i < seq->length; i++) {
				count[(int) seq->seq[i]]++;
				switch (seq->seq[i]) {
					case 'a': case 'A': case 't': case 'T': AT++; break;
					case 'g': case 'G': case 'c': case 'C': GC++; break;
				}
			}
			total_length += seq->length;
			if (seq->length > longest) longest = seq->length;
			if (seq->length < shortest) shortest = seq->length;
			
			zoeDeleteFastaFile(seq);
		}
		zoeCloseFile(file);
	}
	
	/* output some stats */
	zoeS(stdout, "%d files\n", file_count);
	zoeS(stdout, "%g total letters\n", (float)total_length);
	zoeS(stdout, "range %g to %g\n", (float)shortest, (float)longest);
	f = total_length/file_count;
	zoeS(stdout, "%g average\n", f);
	zoeS(stdout, "%g GC\n", (double)GC / (double)(GC+AT));
	
	for (i=0;i<ASCII_SET;i++) {
		if (count[i] != 0) {
			f = (double)count[i]/total_length;
			d = - ((double)count[i]/total_length)
				* zoeLog2((double)count[i]/total_length);
			zoeS(stdout, "%c\t%g\t%f\n", i, (float)count[i], f);
			entropy += d;
		}
	}
	zoeS(stdout, "entropy = %g bits\n", entropy);
	return 0;
}
