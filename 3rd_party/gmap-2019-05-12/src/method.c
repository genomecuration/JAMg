static char rcsid[] = "$Id: method.c 218675 2019-03-16 01:25:48Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "method.h"

#include <stdlib.h>

char *
Method_string (Method_T method) {
  /* May want to turn on NO_COMPARE in substring.c also */
  switch (method) {
  case TR: return "tr";
  case KMER_EXACT: return "exact";
  case KMER_ONE_END: return "one";
  case KMER_APPROX: return "approx";
  case EXT: return "ext";
  case EXT_GMAP: return "ext-gmap";
  case SEGMENT1: return "seg1";
  case SEGMENT2: return "seg2";
  case DISTANT_RNA: return "dist-rna";
  case DISTANT_DNA: return "dist-dna";
  case TERMINAL: return "term";

  case KMER_MIDDLE: return "mid";
  case SEGMENT_GMAP: return "seg-gmap";
  case KMER: return "kmer";
  case KMER_GMAP: return "kmer-gmap";
  case REGION_GMAP: return "gmap";
  default: abort();
  }
}


void
Method_samprint (Filestring_T fp, Method_T method) {
  switch (method) {
  case TR: FPRINTF(fp,"\tXG:Z:tr"); break;
  case KMER_EXACT: FPRINTF(fp,"\tXG:Z:exact"); break;
  case KMER_ONE_END: FPRINTF(fp,"\tXG:Z:one"); break;
  case KMER_APPROX: FPRINTF(fp,"\tXG:Z:approx"); break;
  case EXT: FPRINTF(fp,"\tXG:Z:ext"); break;
  case EXT_GMAP: FPRINTF(fp,"\tXG:Z:ext-gmap"); break;
  case SEGMENT1: FPRINTF(fp,"\tXG:Z:seg1"); break;
  case SEGMENT2: FPRINTF(fp,"\tXG:Z:seg2"); break;
  case DISTANT_RNA: FPRINTF(fp,"\tXG:Z:dist-rna"); break;
  case DISTANT_DNA: FPRINTF(fp,"\tXG:Z:dist-dna"); break;
  case TERMINAL: FPRINTF(fp,"\tXG:Z:term"); break;

  case KMER_MIDDLE: FPRINTF(fp,"\tXG:Z:mid"); break;
  case SEGMENT_GMAP: FPRINTF(fp,"\tXG:Z:seg-gmap"); break;
  case KMER: FPRINTF(fp,"\tXG:Z:kmer"); break;
  case KMER_GMAP: FPRINTF(fp,"\tXG:Z:kmer-gmap"); break;
  case REGION_GMAP: FPRINTF(fp,"\tXG:Z:gmap"); break;
  default: abort();
  }

  return;
}

void
Method_print (Filestring_T fp, Method_T method) {
  switch (method) {
  case TR: FPRINTF(fp,"\tmethod:tr"); break;
  case KMER_EXACT: FPRINTF(fp,"\tmethod:exact"); break;
  case KMER_ONE_END: FPRINTF(fp,"\tmethod:one"); break;
  case KMER_APPROX: FPRINTF(fp,"\tmethod:approx"); break;
  case EXT: FPRINTF(fp,"\tmethod:ext"); break;
  case EXT_GMAP: FPRINTF(fp,"\tmethod:ext-gmap"); break;
  case SEGMENT1: FPRINTF(fp,"\tmethod:seg1"); break;
  case SEGMENT2: FPRINTF(fp,"\tmethod:seg2"); break;
  case DISTANT_RNA: FPRINTF(fp,"\tmethod:dist-rna"); break;
  case DISTANT_DNA: FPRINTF(fp,"\tmethod:dist-dna"); break;
  case TERMINAL: FPRINTF(fp,"\tmethod:term"); break;

  case KMER_MIDDLE: FPRINTF(fp,"\tmethod:mid"); break;
  case SEGMENT_GMAP: FPRINTF(fp,"\tmethod:seg-gmap"); break;
  case KMER: FPRINTF(fp,"\tmethod:kmer"); break;
  case KMER_GMAP: FPRINTF(fp,"\tmethod:kmer-gmap"); break;
  case REGION_GMAP: FPRINTF(fp,"\tmethod:gmap"); break;
  default: abort();
  }

  return;
}

