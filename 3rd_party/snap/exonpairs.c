/*****************************************************************************\
 exonpairs.c

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

char * ZOE = NULL; /* environment variable */

static char usage[]  = "\n\
usage: exonpairs [options] <HMM file> <FASTA file> [options]\n\
options:\n\
  -min-intron <int>     minimum Intron length [30]\n\
  -max-intron <int>     maximum Intron length [10000]\n\
  -einit-length <int>   minimum Einit length in bp [10]\n\
  -eterm-length <int>   minimum Eterm length in bp [10]\n\
  -exon-length <int>    minimum Exon length in bp [30]\n\
  -intron-score <float> minimum Intron score in bits [-5]\n\
  -einit-score <float>  minimum Einit score in bits [-5]\n\
  -eterm-score <float>  minimum Eterm score in bits [-5]\n\
  -exon-score <float>   minimum Exon score in bits [-5]\n\
  -pair-score <float>   minimum pair score (exon-intron-exon) [10]\n\
  -flank-length <int>   length of flanking sequence per exon [20]\n\
  -lcmask               treat lowercase as N\n\
";

int   PADDING             = 48;
int   MIN_INTRON_LENGTH   = 40;
int   MAX_INTRON_LENGTH   = 10000;
int   MIN_EINIT_LENGTH    = 10;
int   MIN_ETERM_LENGTH    = 10;
int   MIN_EXON_LENGTH     = 30;
float MIN_INTRON_SCORE    = -5;
float MIN_EINIT_SCORE     = -5;
float MIN_ETERM_SCORE     = -5;
float MIN_EXON_SCORE      = -5;
float MIN_PAIR_SCORE      = 10;
int   FLANK_LENGTH        = 20;

static zoeFeatureVec get_exons (const zoeTrellis t) {
	int               i, j, skip, length;

	zoeFeatureVec     temp;
	zoeFeatureVec     exons = zoeNewFeatureVec();
	zoeFeature        exon;
	zoeFeatureFactory factory = t->factory[Exon];
	
	for (i = PADDING; i < t->dna->length - PADDING; i++) {
		temp = factory->create(factory, i);
		if (temp) {
			for (j = 0; j < temp->size; j++) {
				exon = temp->elem[j];				
				length = exon->end - exon->start +1;	
				exon->score = zoeScoreExon(t, exon, PADDING, 0);
				
				skip = 0;
				switch (exon->label) {
					case Esngl: skip = 1; break;
					case Einit:
						if (exon->score < MIN_EINIT_SCORE)  skip = 1; break;
						if (length      < MIN_EINIT_LENGTH) skip = 1; break;
					case Eterm:
						if (exon->score < MIN_ETERM_SCORE)  skip = 1; break;
						if (length      < MIN_ETERM_LENGTH) skip = 1; break;
					case Exon:
						if (exon->score < MIN_EXON_SCORE)   skip = 1; break;
						if (length      < MIN_EXON_LENGTH)  skip = 1; break;
					default: zoeExit("no, not possible");
				}
				if (skip) continue;
				
				zoePushFeatureVec(exons, exon);
			}
			zoeDeleteFeatureVec(temp);
		}
	}
	
	return exons;
}

static void report_pair (const zoeDNA dna, const zoeFeature e1, const zoeFeature e2,
		const zoeFeature intron, strand_t strand, zoeHash reported) {
	int start;
	char seq1[FLANK_LENGTH +1], seq2[FLANK_LENGTH +1];
	char location[64];
	int c1, c2;
		
	start = e1->end - FLANK_LENGTH +1;
	strncpy(seq1, dna->seq + start, FLANK_LENGTH);
	seq1[FLANK_LENGTH] = '\0';
	
	start = e2->start;
	strncpy(seq2, dna->seq + start, FLANK_LENGTH);
	seq2[FLANK_LENGTH] = '\0';
	
	/* check if this junction has already been reported */
	sprintf(location, "%d-%d", e1->end, e2->start);
	if (zoeGetHash(reported, location)) return; /* already reported */
	zoeSetHash(reported, location, (void*)1);
	
	/* output */
	if (strand == '+') {
		c1 = e1->end   - PADDING +1;
		c2 = e2->start - PADDING +1;
	} else {
		c2 = dna->length - PADDING - e1->end;
		c1 = dna->length - PADDING - e2->start;
	}
	printf(">%c%d..%d (%d)\n%s\n%s\n", strand, c1, c2, c2-c1+1, seq1, seq2);

}

static void exon_pairs (const zoeDNA dna, const zoeHMM hmm, strand_t strand) {
	zoeTrellis      trellis;
	zoeFeatureVec   exons;
	zoeFeature      e1, e2, intron;
	int             i, j, d;
	score_t         pair_score;
	score_t        *intron_score;
	zoeHash         reported;
	
	/* init */
	zoeSetTrellisMeter(0); /* quiet */
	trellis = zoeNewTrellis(dna, hmm, NULL);
	exons = get_exons(trellis);
	intron = zoeNewFeature(Intron, 0, 0, '+', 0, 0, 0, 0, 0);
	reported = zoeNewHash();
	
	/* precompute intron scores */
	intron_score = malloc(trellis->dna->length * sizeof(score_t));
	for (i = 0; i < PADDING; i++) intron_score[i] = 0;
	for (i = PADDING; i < trellis->dna->length; i++) {
		intron_score[i] = intron_score[i-1]
			+ trellis->scanner[Int0]->score(trellis->scanner[Int0], i)
			- trellis->exp_score;
	}

	/* main loop */
	for (i = 0; i < exons->size; i++) {
		e1 = exons->elem[i];

		for (j = i+1; j < exons->size; j++) {
			e2 = exons->elem[j];
			
			/* exon type */
			if ((e1->label == Einit && e2->label == Einit) ||
				(e1->label == Eterm) ||
				(e1->label == Exon && e2->label == Einit)) continue;
			
			/* exon phase */
			if ((e1->inc3 == 0 && e2->inc5 != 0) ||
				(e1->inc3 == 1 && e2->inc5 != 2) ||
				(e1->inc3 == 2 && e2->inc5 != 1)) continue;
					
			/* intron length */
			intron->start = e1->end +1;
			intron->end   = e2->start -1;
			d = intron->end - intron->start +1;
			if (d < MIN_INTRON_LENGTH) continue;
			if (d > MAX_INTRON_LENGTH) break;
			
			/* intron score */
			intron->score = intron_score[intron->end] - intron_score[intron->start]
				+ zoeScoreDuration(trellis->hmm->dmap[Int0], intron->end - intron->start +1);
			if (intron->score < MIN_INTRON_SCORE) continue;
			
			/* pair score */
			pair_score = e1->score + e2->score + intron->score;
			if (pair_score < MIN_PAIR_SCORE) continue;
			
			/* output */
			report_pair(trellis->dna, e1, e2, intron, strand, reported);
		}
		
	}
	
	fprintf(stderr, "%d exons found on %c strand\n", reported->keys->size, strand);
	
	/* clean up */
	zoeDeleteFeature(intron);
	zoeDeleteHash(reported);
	zoeDeleteTrellis(trellis);
	zoeDeleteFeatureVec(exons);
	
}
/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char *argv[]) {
	zoeFile         dna_file;
	zoeFastaFile    fasta;
	zoeDNA          dna, rc;
	zoeHMM          hmm;
	
	/* command line */
	zoeSetProgramName(argv[0]);
	zoeSetOption("-lcmask", 0);
	zoeSetOption("-min-intron", 1);
	zoeSetOption("-max-intron", 1);
	zoeSetOption("-intron-score", 1);
	zoeSetOption("-einit-length", 1);
	zoeSetOption("-eterm-length", 1);
	zoeSetOption("-exon-length", 1);
	zoeSetOption("-einit-score", 1);
	zoeSetOption("-eterm-score", 1);
	zoeSetOption("-exon-score", 1);
	zoeSetOption("-pair-score", 1);
	zoeSetOption("-flank-length", 1);
	zoeParseOptions(&argc, argv);
	if (argc != 3) {
		zoeE("%s", usage);
		exit(1);
	}
	
	/* options */
	if (zoeOption("-min-intron"))   MIN_INTRON_LENGTH = atoi(zoeOption("-min-intron"));
	if (zoeOption("-max-intron"))   MAX_INTRON_LENGTH = atoi(zoeOption("-max-intron"));
	if (zoeOption("-einit-length")) MIN_EINIT_LENGTH  = atoi(zoeOption("-einit-length"));
	if (zoeOption("-eterm-length")) MIN_ETERM_LENGTH  = atoi(zoeOption("-eterm-length"));
	if (zoeOption("-exon-length"))  MIN_EXON_LENGTH   = atoi(zoeOption("-exon-length"));
	if (zoeOption("-intron-score")) MIN_INTRON_SCORE  = atof(zoeOption("-intron-score"));
	if (zoeOption("-einit-score"))  MIN_EINIT_SCORE   = atof(zoeOption("-einit-score"));
	if (zoeOption("-eterm-score"))  MIN_ETERM_SCORE   = atof(zoeOption("-eterm-score"));
	if (zoeOption("-exon-score"))   MIN_EXON_SCORE    = atof(zoeOption("-exon-score"));
	if (zoeOption("-pair-score"))   MIN_PAIR_SCORE    = atof(zoeOption("-pair-score"));
	if (zoeOption("-flank-length")) FLANK_LENGTH      = atoi(zoeOption("-flank-length"));
	
	/* HMM and FASTA */
	ZOE = getenv("ZOE");
	hmm = zoeGetHMM(argv[1]);
	dna_file = zoeOpenFile(argv[2]);
	fasta = zoeReadFastaFile(dna_file.stream);
	dna = zoeNewDNA(fasta->def, fasta->seq);
	zoeDeleteFastaFile(fasta);
	if (zoeOption("-lcmask")) zoeLCsmooth(dna, 10, 10,100);
	rc = zoeAntiDNA(">anti", dna);
	
	/* do it */
	exon_pairs(dna, hmm, '+');
	exon_pairs(rc,  hmm, '-');
			
	return 0;
}


