static char rcsid[] = "$Id: gmapindex.c 218147 2019-01-16 21:28:41Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include <ctype.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <sys/mman.h>		/* For munmap */
#include "types.h"
#include "bool.h"
#include "assert.h"
#include "mem.h"
#include "fopen.h"
#include "getline.h"
#include "filesuffix.h"

#include "table.h"
#ifdef HAVE_64_BIT
#include "tableuint8.h"
typedef Tableuint8_T Table_chrpos_T;
#else
typedef Tableuint_T Table_chrpos_T;
#endif
#include "tableuint.h"
#include "compress.h"
#include "chrom.h"
#include "segmentpos.h"
#include "univinterval.h"
#include "iit-write-univ.h"
#include "iit-read-univ.h"
#include "genome.h"
#include "genome128_hr.h"
#include "genome-write.h"
#include "indexdb-write.h"
#include "localdb-write.h"
#include "compress-write.h"
#include "intlist.h"
#include "indexdb.h"		/* For Indexdb_filenames_T */
#include "indexdbdef.h"		/* For compression types */
#include "bitpack64-write.h"
#include "parserange.h"

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#if !defined(HAVE_SSE4_2)
/* Skip popcnt */
#elif defined(HAVE_POPCNT)
#include <immintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip mm_popcnt */
#elif defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip lzcnt and tzcnt */
#elif defined(HAVE_LZCNT) || defined(HAVE_TZCNT)
#include <immintrin.h>
#endif


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#define BUFFERSIZE 8192

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Program variables */
typedef enum {NONE, AUXFILES, GENOME, UNSHUFFLE, COUNT, OFFSETS, POSITIONS, 
	      LOCALDB} Action_T;
static Action_T action = NONE;
static char *sourcedir = ".";
static char *destdir = ".";
static char *fileroot = NULL;
static int compression_types = BITPACK64_COMPRESSION;
static int compression_type;
static int index1part = 15;
static int index1interval = 3;	/* Interval for storing 15-mers */
static int local1part = 8;
static int local1interval = 1;
static bool genome_lc_p = false;
static bool rawp = false;
static bool writefilep = false;
/* static bool sortchrp = true;	? Sorting now based on order in .coords file */
static int wraplength = 0;
static bool mask_lowercase_p = false;
static int nmessages = 50;

static char *mitochondrial_string = NULL;
static Sorttype_T divsort = CHROM_SORT;
static char *sortfilename = NULL;
static bool huge_offsets_p = false;


#if 0
/* Using non-SIMD code, so no need to check compiler assumptions */
static void
check_compiler_assumptions () {
  unsigned int x = rand(), y = rand();
#ifdef HAVE_SSE2
  int z;
  __m128i a;
#ifdef HAVE_SSE4_1
  char negx, negy;
#endif
#endif


  fprintf(stderr,"Checking compiler assumptions for popcnt: ");
  fprintf(stderr,"%08X ",x);
#ifdef HAVE_LZCNT
  fprintf(stderr,"_lzcnt_u32=%d ",_lzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stderr,"__builtin_clz=%d ",__builtin_clz(x));
#endif
#ifdef HAVE_TZCNT
  fprintf(stderr,"_tzcnt_u32=%d ",_tzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stderr,"__builtin_ctz=%d ",__builtin_ctz(x));
#endif

#ifdef HAVE_POPCNT
  fprintf(stderr,"_popcnt32=%d ",_popcnt32(x));
#endif
#if defined(HAVE_MM_POPCNT)
  fprintf(stderr,"_mm_popcnt_u32=%d ",_mm_popcnt_u32(x));
#endif
#if defined(HAVE_BUILTIN_POPCOUNT)
  fprintf(stderr,"__builtin_popcount=%d ",__builtin_popcount(x));
#endif

  fprintf(stderr,"\n");

#ifdef HAVE_SSE2
  fprintf(stderr,"Checking compiler assumptions for SSE2: ");
  fprintf(stderr,"%08X %08X",x,y);
  a = _mm_xor_si128(_mm_set1_epi32(x),_mm_set1_epi32(y));
  z = _mm_cvtsi128_si32(a);
  fprintf(stderr," xor=%08X\n",z);

#ifdef HAVE_SSE4_1
  if ((negx = (char) x) > 0) {
    negx = -negx;
  }
  if ((negy = (char) y) > 0) {
    negy = -negy;
  }

  fprintf(stderr,"Checking compiler assumptions for SSE4.1: ");
  fprintf(stderr,"%d %d",negx,negy);
  a = _mm_max_epi8(_mm_set1_epi8(negx),_mm_set1_epi8(negy));
  z = _mm_extract_epi8(a,0);
  fprintf(stderr," max=%d => ",z);
  if (negx > negy) {
    if (z == (int) negx) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  } else {
    if (z == (int) negy) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  }

#endif

#endif

  fprintf(stderr,"Finished checking compiler assumptions\n");

  return;
}
#endif


#if 0
/************************************************************************
 *   Reading strain from file
 ************************************************************************/

static char *
read_strain_from_strainfile (char *strainfile) {
  FILE *fp;
  char *refstrain = NULL, *line, *strain, *straintype;
  int line_length;

  if (strainfile != NULL) {
    fp = fopen(strainfile,"r");
    if (fp == NULL) {
      fprintf(stderr,"Cannot open strain file %s\n",strainfile);
    } else {
      while ((line = Getline_wlength(&line_length,fp)) != NULL) {
	if (line[0] == '#') {
	  /* Skip */
	} else {
	  strain = (char *) MALLOC((line_length+1)*sizeof(char));
	  straintype = (char *) MALLOC((line_length+1)*sizeof(char));
	  if (sscanf(line,"%s %s",strain,straintype) == 2) {
	    if (!strcmp(straintype,"reference") || !strcmp(straintype,"Reference") || 
		!strcmp(straintype,"REFERENCE")) {
	      if (refstrain != NULL) {
		fprintf(stderr,"More than one reference strain seen in %s\n",strainfile);
		exit(9);
	      }
	      refstrain = (char *) CALLOC(strlen(strain)+1,sizeof(char));
	      strcpy(refstrain,strain);
	    }
	  }
	  FREE(straintype);
	  FREE(strain);
	}
	FREE(line);
      }

      fclose(fp);
    }
  }

  if (refstrain != NULL) {
    return refstrain;
  } else {
    refstrain = (char *) CALLOC(strlen("reference")+1,sizeof(char));
    strcpy(refstrain,"reference");
    return refstrain;
  }
}

static char *
read_strain_from_coordsfile (char *coordsfile) {
  FILE *fp;
  char *refstrain = NULL, *line, *strain, *ptr;
  int line_length;

  if (coordsfile != NULL) {
    fp = fopen(coordsfile,"r");
    if (fp == NULL) {
      fprintf(stderr,"Cannot open coords file %s\n",coordsfile);
    } else {
      while ((line = Getline_wlength(&line_length,fp)) != NULL) {
	if (line[0] == '#') {
	  if ((ptr = strstr(line,"Reference strain:")) != NULL) {
	    strain = (char *) MALLOC((line_length+1)*sizeof(char));
	    if (sscanf(ptr,"Reference strain: %s",strain) == 1) {
	      if (refstrain != NULL) {
		fprintf(stderr,"More than one reference strain seen in %s\n",coordsfile);
		exit(9);
	      }
	      refstrain = (char *) CALLOC(strlen(strain)+1,sizeof(char));
	      strcpy(refstrain,strain);
	    }
	    FREE(strain);
	  }
	}	    
	FREE(line);
      }

      fclose(fp);
    }
  }

  if (refstrain != NULL) {
    return refstrain;
  } else {
    refstrain = (char *) CALLOC(strlen("reference")+1,sizeof(char));
    strcpy(refstrain,"reference");
    return refstrain;
  }
}
#endif




/************************************************************************
 *   Creating aux file
 ************************************************************************/

/* accsegmentpos_table: char *accession -> Segmentpos_T segmentpos
   chrlength_table:     Chrom_T chrom -> Chrpos_T chrlength
*/

static void
chrlength_update (Table_chrpos_T chrlength_table, Chrom_T chrom, Univcoord_T segend) {
  Univcoord_T oldsegend;

#ifdef HAVE_64_BIT
  if ((oldsegend = (Univcoord_T) Tableuint8_get(chrlength_table,chrom)) == 0) {
    /* Initial entry for this chromosome */
    Tableuint8_put(chrlength_table,chrom,segend);

  } else if (segend > oldsegend) {
    /* Revise */
    Tableuint8_put(chrlength_table,chrom,segend);
  }
#else
  if ((oldsegend = (Univcoord_T) Tableuint_get(chrlength_table,chrom)) == 0) {
    /* Initial entry for this chromosome */
    Tableuint_put(chrlength_table,chrom,segend);

  } else if (segend > oldsegend) {
    /* Revise */
    Tableuint_put(chrlength_table,chrom,segend);
  }
#endif

  return;
}

static void
store_accession (Table_T accsegmentpos_table, Table_chrpos_T chrlength_table, Tableuint_T chrorder_table,
		 char *accession, char *chr_string, Chrpos_T chrpos1, 
		 Chrpos_T chrpos2, bool revcompp, Chrpos_T seglength, 
		 int contigtype, unsigned int universal_coord, bool circularp,
		 Chrpos_T alt_scaffold_start, Chrpos_T alt_scaffold_end) {
  Chrom_T chrom;
  Segmentpos_T segmentpos;
  unsigned int order;

  if (chrorder_table != NULL) {
    order = Tableuint_get(chrorder_table,chr_string);
    chrom = Chrom_from_string(chr_string,mitochondrial_string,order,circularp,
			      alt_scaffold_start,alt_scaffold_end);
  } else {
    chrom = Chrom_from_string(chr_string,mitochondrial_string,/*order*/universal_coord,circularp,
			      alt_scaffold_start,alt_scaffold_end);
  }

  segmentpos = Segmentpos_new(chrom,chrpos1,chrpos2,revcompp,seglength,contigtype);
  Table_put(accsegmentpos_table,(void *) accession,(void *) segmentpos);

  /* Update chrlength */
  if (chrpos2 > chrpos1 + seglength) {
    chrlength_update(chrlength_table,chrom,chrpos2);
  } else {
    chrlength_update(chrlength_table,chrom,chrpos1+seglength);
  }

  return;
}


/* We assume that header has already been read.  We need to check each
   new line for a new header */
static Chrpos_T
count_sequence () {
  Chrpos_T seglength = 0U;
  int c;
  char Buffer[BUFFERSIZE], *p;
  bool newline = true;

  while (1) {
    /* Start of new line */
    if (newline == true) {
      if ((c = getc(stdin)) == EOF || c == '>') {
	return seglength;
      } else {
	seglength += 1;
      }
    }

    if (fgets(Buffer,BUFFERSIZE,stdin) == NULL) {
      return seglength;
    } else {
      if ((p = rindex(Buffer,'\n')) != NULL) {
	*p = '\0';
	newline = true;
      } else {
	newline = false;
      }
      seglength += (Chrpos_T) strlen(Buffer);
    }
  }
}

static void
skip_sequence (Chrpos_T seglength) {
  int c;
  char Buffer[BUFFERSIZE];

  while (seglength > BUFFERSIZE) {
    if (fread(Buffer,sizeof(char),BUFFERSIZE,stdin) < BUFFERSIZE) {
      fprintf(stderr,"End of file reached.  Expecting %u more characters\n",seglength);
      exit(9);
    }
    seglength -= BUFFERSIZE;
  }

  if (seglength > 0U) {
    if (fread(Buffer,sizeof(char),seglength,stdin) < seglength) {
      fprintf(stderr,"End of file reached.  Expecting %u more characters\n",seglength);
      exit(9);
    }
  }

  if ((c = getchar()) != EOF && c != '\n') {
    fprintf(stderr,"Expecting linefeed at end of sequence.  Saw %d (%c) instead\n",c,c);
    exit(9);
  }

  if ((c = getchar()) != EOF && c != '>') {
    fprintf(stderr,"Expecting new FASTA line.  Saw %d (%c) instead\n",c,c);
    exit(9);
  }

  return;
}

static bool
process_sequence_aux (Chrpos_T *seglength, Table_T accsegmentpos_table, Table_chrpos_T chrlength_table,
		      Tableuint_T chrorder_table, Table_T alt_scaffold_info_table,
		      char *fileroot, int ncontigs) {
  char *line, *accession_p, *accession, 
    *chrpos_string, *primary_string, *chr_string, *coords, *ptr, *p;
  int line_length;
  Chrpos_T chrpos1, chrpos2, lower, upper, primary_chrstart, primary_chrend;
  Univcoord_T universal_coord = 0U;
  bool revcompp, circularp;
  Chrpos_T alt_scaffold_start, alt_scaffold_end;
  char *alt_scaffold_chr, *altscaffold_info, *primary_chr;
  int nitems;

  /* Store sequence info */
  if ((line = Getline_wlength(&line_length,stdin)) == NULL) {
    return false;
  } else {
    accession_p = (char *) MALLOC((line_length+1)*sizeof(char));
    chrpos_string = (char *) MALLOC((line_length+1)*sizeof(char));
  }

  nitems = sscanf(line,"%s %s %llu",accession_p,chrpos_string,&universal_coord);

  if (nitems < 2) {
    fprintf(stderr,"Can't parse line %s\n",line);
    exit(1);
  } else {
    if (ncontigs < nmessages) {
      fprintf(stderr,"Logging contig %s at %s in genome %s\n",accession_p,chrpos_string,fileroot);
    } else if (ncontigs == nmessages) {
      fprintf(stderr,"More than %d contigs.  Will stop printing messages\n",nmessages);
    }

    if ((coords = rindex(chrpos_string,':')) == NULL) {
      fprintf(stderr,"Can't parse chromosomal coordinates %s\n",chrpos_string);
      exit(1);
    } else {
      chr_string = (char *) MALLOC((coords-chrpos_string+1)*sizeof(char));
      strncpy(chr_string,chrpos_string,(coords-chrpos_string)*sizeof(char));
      chr_string[coords-chrpos_string] = '\0';

      coords = &(coords[1]);	/* To skip ':' */
      if (sscanf(coords,"%u..%u",&chrpos1,&chrpos2) == 2) {
	/* 1:3..5, one-based, inclusive => (2,5), zero-based, boundaries */
	if (chrpos1 <= chrpos2) {
	  chrpos1--;
	  revcompp = false;
	  lower = chrpos1;
	  upper = chrpos2;
	} else {
	  chrpos2--;
	  revcompp = true;
	  lower = chrpos2;
	  upper = chrpos1;
	}
      } else if (sscanf(coords,"%u",&chrpos1) == 1) {
	/* 1:3, one-based, inclusive => (3,3), zero-based, boundaries */
	revcompp = false;
	lower = upper = chrpos1;
      } else {
	fprintf(stderr,"Can't parse chromosomal coordinates %s\n",coords);
	exit(1);
      }
    }

#if 0
    /* No longer supporting strains/types */
    p = line;
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to first space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past first space */
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to second space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past second space */

    if (*p == '\0') {
      contigtype = 0;		/* Empty type string */
    } else {
#if 0
      if ((ptr = rindex(p,'\n')) != NULL) {
	/* No longer necessary with call to Getline */
	while (isspace((int) *ptr)) { ptr--; } /* Erase empty space */
	ptr++;
	*ptr = '\0';
      }
#endif

      /* Erase empty space at end of line */
      ptr = &(line[line_length-1]);
      while (isspace((int) *ptr)) { ptr--; }
      ptr++;
      *ptr = '\0';

      if ((contigtype = Tableuint_get(contigtype_table,(void *) p)) == 0) {
	debug(printf("Storing type %s.\n",p));
	/* Store types as 1-based */
	contigtype = Tableuint_length(contigtype_table) + 1;
	typestring = (char *) CALLOC(strlen(p)+1,sizeof(char));
	strcpy(typestring,p);
	Tableuint_put(contigtype_table,(void *) typestring,contigtype);
	*contigtypelist = List_push(*contigtypelist,typestring);
      }
    }
#endif

    /* Look for circular.  Code modeled after parsing for strain above. */
    p = line;
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to first space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past first space */
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to second space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past second space */
    while (*p != '\0' && !isspace((int) *p)) { p++; } /* Skip to third space */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Skip past third space */

    circularp = false;
    alt_scaffold_start = alt_scaffold_end = 0;
    if (*p == '\0') {
      /* circularp = false; */
      if (ncontigs < nmessages) {
	fprintf(stderr," => primary (linear) chromosome");
      }
    } else {
#if 0
      if ((ptr = rindex(p,'\n')) != NULL) {
	/* No longer necessary with call to Getline */
	while (isspace((int) *ptr)) { ptr--; } /* Erase empty space */
	ptr++;
	*ptr = '\0';
      }
#endif

      /* Erase empty space at end of line */
      ptr = &(line[line_length-1]);
      while (isspace((int) *ptr)) { ptr--; }
      ptr++;
      *ptr = '\0';

      if (!strcmp(p,"circular")) {
	if (ncontigs < nmessages) {
	  fprintf(stderr," => circular chromosome");
	}
	circularp = true;
      } else if (!strcmp(p,"linear")) {
	/* Ignore.  Not a altloc chromosome. */
	if (ncontigs < nmessages) {
	  fprintf(stderr," => primary (linear) chromosome");
	}

      } else {
	/* Alternate scaffold.  Store primary information as a new type. */
	if (ncontigs < nmessages) {
	  fprintf(stderr," => alternate scaffold for segment %s",p);
	}

	if (sscanf(p,"%u..%u",&alt_scaffold_start,&alt_scaffold_end) == 2) {
	  while (*p != '\0' && *p != ',') {
	    p++;
	  }
	  if (*p == '\0') {
	    fprintf(stderr,"Expected to find a comma in the altscaffold tag\n");
	  } else {
	    p++;
	    if (Parserange_simple(&primary_chr,&revcompp,&primary_chrstart,&primary_chrend,p) == false) {
	      fprintf(stderr,"Cannot parse %s\n",p);
	      abort();
	    } else {
	      /* Maximum length of an unsigned char is 10 digits */
	      primary_string = (char *) MALLOC((strlen(primary_chr)+10+2+10+1+1+10+2+10+1+1)*sizeof(char));
	      sprintf(primary_string,"%u..%u,%s:%u..%u\n",
		      alt_scaffold_start,alt_scaffold_end,primary_chr,primary_chrstart,primary_chrend);
	    }
	    FREE(primary_chr);

	    altscaffold_info = (char *) CALLOC(strlen(primary_string)+1,sizeof(char));
	    strcpy(altscaffold_info,primary_string);
	    
	    alt_scaffold_chr = (char *) CALLOC(strlen(chr_string)+1,sizeof(char));
	    strcpy(alt_scaffold_chr,chr_string);

	    Table_put(alt_scaffold_info_table,(void *) alt_scaffold_chr,(void *) altscaffold_info);
	    FREE(primary_string);
	  }
	}
      }
    }

    /* The '>' character was already stripped off by the last call to count_sequence() */
    accession = (char *) CALLOC(strlen(accession_p)+1,sizeof(char));
    strcpy(accession,accession_p);

    if (ncontigs < nmessages) {
      fprintf(stderr,"\n");
    }

    FREE(chrpos_string);
    FREE(accession_p);
    FREE(line);
  }

  if (rawp == true) {
    *seglength = upper - lower;
    fprintf(stderr,"Skipping %u characters\n",*seglength);
    skip_sequence(*seglength);
  } else {
    *seglength = count_sequence();
    if (*seglength != upper - lower) {
      fprintf(stderr,"%s has expected sequence length %u-%u=%u but actual length %u\n",
	      accession,upper,lower,upper-lower,*seglength);
    }
  }

  if (nitems < 3) {
    universal_coord = 0U;
  }

  store_accession(accsegmentpos_table,chrlength_table,chrorder_table,
		  accession,chr_string,lower,upper,revcompp,
		  *seglength,/*contigtype*/0,universal_coord,circularp,
		  alt_scaffold_start,alt_scaffold_end);
  FREE(chr_string);
    
  return true;
}


/************************************************************************
 *   Creating genome and related files
 ************************************************************************/

/* Modifies chrlength_table to store offsets, rather than chrlengths */
static void
write_chromosome_file (char *genomesubdir, char *fileroot, Table_chrpos_T chrlength_table,
		       bool coord_values_8p) {
  FILE *textfp, *chrsubsetfp;
  char *divstring, *textfile, *chrsubsetfile, *iitfile, *chr_string, emptystring[1];
  int n, i;
  int typeint;
  Chrom_T *chroms;
  Univcoord_T chroffset = 0;
  Chrpos_T chrlength, alt_scaffold_start, alt_scaffold_end;
  List_T divlist = NULL;
  List_T intervallist = NULL, chrtypelist = NULL, labellist = NULL, annotlist = NULL, p;
  Table_T intervaltable, labeltable, annottable;
  Univinterval_T interval;

  emptystring[0] = '\0';

  if (divsort == NO_SORT) {
    fprintf(stderr,"divsort == NO_SORT\n");
#ifdef HAVE_64_BIT
    chroms = (Chrom_T *) Tableuint8_keys_by_timeindex(chrlength_table,0U);
    n = Tableuint8_length(chrlength_table);
#else
    chroms = (Chrom_T *) Tableuint_keys_by_timeindex(chrlength_table,0U);
    n = Tableuint_length(chrlength_table);
#endif

  } else {
    /* Get chromosomes in order */
#ifdef HAVE_64_BIT
    chroms = (Chrom_T *) Tableuint8_keys(chrlength_table,0U);
    n = Tableuint8_length(chrlength_table);
#else
    chroms = (Chrom_T *) Tableuint_keys(chrlength_table,0U);
    n = Tableuint_length(chrlength_table);
#endif
    switch (divsort) {
    case ALPHA_SORT: qsort(chroms,n,sizeof(Chrom_T),Chrom_compare_alpha); break;
    case NUMERIC_ALPHA_SORT: qsort(chroms,n,sizeof(Chrom_T),Chrom_compare_numeric_alpha); break;
    case CHROM_SORT: qsort(chroms,n,sizeof(Chrom_T),Chrom_compare_chrom); break;
    case FILENAME_SORT: qsort(chroms,n,sizeof(Chrom_T),Chrom_compare_order); break;
    default: abort();
    }
  }
  fprintf(stderr,"Have a total of %d chromosomes\n",n);
    

  /* Write chromosome text file and chrsubset file */
  textfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".chromosome")+1,sizeof(char));
  sprintf(textfile,"%s/%s.chromosome",genomesubdir,fileroot);
  fprintf(stderr,"Writing chromosome file %s\n",textfile);

  /* Use binary, not text, so files are Unix-compatible */
  if ((textfp = FOPEN_WRITE_BINARY(textfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",textfile);
    exit(9);
  }
  FREE(textfile);

  chrsubsetfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				  strlen(fileroot)+strlen(".chrsubset")+1,sizeof(char));
  sprintf(chrsubsetfile,"%s/%s.chrsubset",genomesubdir,fileroot);
  /* Use binary, not text, so files are Unix-compatible */
  if ((chrsubsetfp = FOPEN_WRITE_BINARY(chrsubsetfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",chrsubsetfile);
    exit(9);
  }
  FREE(chrsubsetfile);
  fprintf(chrsubsetfp,">all\n");
  fprintf(chrsubsetfp,"\n");

  chrtypelist = List_push(chrtypelist,"circular"); /* typeint 1 */
  chrtypelist = List_push(chrtypelist,"");	   /* typeint 0 */


  for (i = 0; i < n; i++) {
#ifdef HAVE_64_BIT
    chrlength = (Chrpos_T) Tableuint8_get(chrlength_table,chroms[i]);
#else
    chrlength = (Chrpos_T) Tableuint_get(chrlength_table,chroms[i]);
#endif
    assert(chroffset <= chroffset + (Univcoord_T) chrlength - 1);
    chr_string = Chrom_string(chroms[i]);
    if (i < nmessages) {
      fprintf(stderr,"Chromosome %s has universal coordinates %llu..%llu\n",
	      chr_string,(unsigned long long) chroffset+1,
	      (unsigned long long) chroffset + 1 + (Univcoord_T) chrlength - 1);
    } else if (i == nmessages) {
      fprintf(stderr,"More than %d contigs.  Will stop printing messages\n",nmessages);
    }
      
    if (n <= 100) {
      fprintf(chrsubsetfp,">%s\n",chr_string);
      fprintf(chrsubsetfp,"+%s\n",chr_string);
    }

    fprintf(textfp,"%s\t%llu..%llu\t%u",
	    chr_string,(unsigned long long) chroffset+1,
	    (unsigned long long) chroffset + (Univcoord_T) chrlength,chrlength);
    if (Chrom_circularp(chroms[i]) == true) {
      fprintf(textfp,"\tcircular");
      typeint = 1;
    } else {
      typeint = 0;
    }
    fprintf(textfp,"\n");

    intervallist = List_push(intervallist,(void *) Univinterval_new(chroffset,chroffset + (Univcoord_T) chrlength - 1,typeint));
    labellist = List_push(labellist,(void *) chr_string);
    annotlist = List_push(annotlist,(void *) emptystring); /* No annotations */

    /* Now chrlength_table holds chroffsets, not chrlengths */
#ifdef HAVE_64_BIT
    Tableuint8_put(chrlength_table,chroms[i],chroffset);
#else
    Tableuint_put(chrlength_table,chroms[i],chroffset);
#endif
    if (Chrom_circularp(chroms[i]) == true) {
      chroffset += (Univcoord_T) chrlength;
      chroffset += (Univcoord_T) chrlength;
    } else if (Chrom_altlocp(&alt_scaffold_start,&alt_scaffold_end,chroms[i]) == true) {
      /* Handle same as primary, non-circular */
      chroffset += chrlength;
    } else {
      chroffset += (Univcoord_T) chrlength;
    }
  }
  FREE(chroms);
  intervallist = List_reverse(intervallist);
  labellist = List_reverse(labellist);

  fclose(chrsubsetfp);
  fclose(textfp);

  /* Write chromosome IIT file */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(NULL,(void *) divstring);

  intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
  annottable = Table_new(65522,Table_string_compare,Table_string_hash);

  Table_put(intervaltable,(void *) divstring,intervallist);
  Table_put(labeltable,(void *) divstring,labellist);
  Table_put(annottable,(void *) divstring,annotlist);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  fprintf(stderr,"Writing chromosome IIT file %s\n",iitfile);
  IIT_write_univ(iitfile,divlist,chrtypelist,
		 intervaltable,labeltable,annottable,
		 coord_values_8p,/*label_pointers_8p*/false,/*annot_pointers_8p*/false);
  FREE(iitfile);

  List_free(&divlist);
  FREE(divstring);

  Table_free(&annottable);
  Table_free(&labeltable);
  Table_free(&intervaltable);

  List_free(&annotlist);

  /* Do not free strings in labellist, since they are not allocated */
  List_free(&labellist);

  /* chrtypelist has no dynamically allocated strings */
  List_free(&chrtypelist);

  for (p = intervallist; p != NULL; p = List_next(p)) {
    interval = (Univinterval_T) List_head(p);
    Univinterval_free(&interval);
  }
  List_free(&intervallist);

  return;
}

static void
write_alt_scaffold_file (char *genomesubdir, char *fileroot, Table_T alt_scaffold_info_table,
			 Univ_IIT_T chromosome_iit, bool coord_values_8p) {
  char *iitfile, *string;
  char **alt_chrs, *altloc_chr, *altscaffold_info, *primary_chr, *p, *q;
  char *divstring, *annot;
  Univinterval_T interval;
  Univcoord_T chroffset;
  Chrpos_T alt_scaffold_start, alt_scaffold_end, primary_start, primary_end;
  int n, i;
  int index;
  List_T intervallist = NULL, labellist = NULL, annotlist = NULL, chrtypelist = NULL, divlist = NULL;
  Table_T intervaltable, labeltable, annottable;
  /* bool revcompp; */


  if ((n = Table_length(alt_scaffold_info_table)) == 0) {
    printf("No alternate scaffolds observed\n");

  } else {
    alt_chrs = (char **) Table_keys(alt_scaffold_info_table,0U);

    for (i = 0; i < n; i++) {
      altloc_chr = alt_chrs[i];
      if ((index = Univ_IIT_find_linear(chromosome_iit,altloc_chr)) < 0) {
	fprintf(stderr,"Error: Could not find alternate chromosome %s in IIT file\n",altloc_chr);
	abort();
      } else {
	interval = Univ_IIT_interval(chromosome_iit,index);
	chroffset = Univinterval_low(interval);
      }

      altscaffold_info = (char *) Table_get(alt_scaffold_info_table,(void *) altloc_chr);
      primary_chr = (char *) MALLOC((strlen(altscaffold_info)+1)*sizeof(char));
#if 0
      fprintf(textfp,"\t%s",altscaffold_info); /* contains a tab character */
#endif
      if (sscanf(altscaffold_info,"%u..%u",&alt_scaffold_start,&alt_scaffold_end) < 2) {
	fprintf(stderr,"Could not parse altscaffold_info %s\n",altscaffold_info);
	abort();
      } else {
	p = altscaffold_info;
	while (*p != '\0' && *p != ',') {
	  p++;
	}
	if (*p == '\0') {
	  fprintf(stderr,"Expected to find a comma in the altscaffold tag\n");
	} else {
	  p++;
	  q = p;
	  while (*q != '\0' && *q != ':') {
	    q++;
	  }
	  strncpy(primary_chr,p,q-p);
	  primary_chr[q-p] = '\0';
	  q++;
	  p = q;

	  if (sscanf(p,"%u..%u\n",&primary_start,&primary_end) < 2) {
	    fprintf(stderr,"Could not parse altscaffold_info %s\n",p);
	    abort();
	  } else if ((index = Univ_IIT_find_linear(chromosome_iit,primary_chr)) < 0) {
	    fprintf(stderr,"Warning: Could not find primary chromosome %s in IIT file, so ignoring alt scaffold %s:%u..%u\n",
		    primary_chr,altloc_chr,alt_scaffold_start,alt_scaffold_end);
	  } else {
	    intervallist = List_push(intervallist, 
				     (void *) Univinterval_new(chroffset + alt_scaffold_start-1, /* zero-based */
							       chroffset + alt_scaffold_end,
							       /*type*/0));
	    labellist = List_push(labellist,(void *) altloc_chr);
	  
	    string = (char *) MALLOC((strlen(primary_chr)+1+10+2+10+1)*sizeof(char));
	    sprintf(string,"%s:%u..%u",primary_chr,primary_start,primary_end);
	    annot = (char *) CALLOC(strlen(string)+1,sizeof(char));
	    strcpy(annot,string);
	    annotlist = List_push(annotlist,(void *) annot);
	    FREE(string);
	  }
	}
      }

      FREE(primary_chr);
    }

    /* Write alt_scaffold IIT file */
    divstring = (char *) CALLOC(1,sizeof(char));
    divstring[0] = '\0';
    divlist = List_push(NULL,divstring);

    intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
    labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
    annottable = Table_new(65522,Table_string_compare,Table_string_hash);

    Table_put(intervaltable,(void *) divstring,intervallist);
    Table_put(labeltable,(void *) divstring,labellist);
    Table_put(annottable,(void *) divstring,annotlist);

    chrtypelist = List_push(chrtypelist,"");	   /* typeint 0 */

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".altscaffold.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.altscaffold.iit",genomesubdir,fileroot);
    fprintf(stderr,"Writing altscaffold IIT file %s\n",iitfile);
    IIT_write_univ(iitfile,divlist,chrtypelist,
		   intervaltable,labeltable,annottable,
		   coord_values_8p,/*label_pointers_8p*/false,/*annot_pointers_8p*/false);
    FREE(iitfile);

    List_free(&chrtypelist); /* chrtypelist has a single static string */

    List_free(&divlist);
    FREE(divstring);

    Table_free(&annottable);
    Table_free(&labeltable);
    Table_free(&intervaltable);

    List_free(&annotlist);
  }

  return;
}


static Table_T current_accsegmentpos_table;

#if 0
/* x is really a char ** */
static int
bysegmentpos_compare_alpha (const void *x, const void *y) {
  char *acc1 = * (char **) x;
  char *acc2 = * (char **) y;
  Segmentpos_T a = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc1);
  Segmentpos_T b = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc2);

  return Segmentpos_compare_alpha(&a,&b);
}
#endif


static int
bysegmentpos_compare (const void *x, const void *y) {
  char *acc1 = * (char **) x;
  char *acc2 = * (char **) y;
  Segmentpos_T a = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc1);
  Segmentpos_T b = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc2);

  if (divsort == ALPHA_SORT) {
    return Segmentpos_compare_alpha(&a,&b);
  } else if (divsort == NUMERIC_ALPHA_SORT) {
    return Segmentpos_compare_numeric_alpha(&a,&b);
  } else if (divsort == CHROM_SORT) {
    return Segmentpos_compare_chrom(&a,&b);
  } else if (divsort == FILENAME_SORT) {
    return Segmentpos_compare_order(&a,&b);
  } else {
    abort();
  }
}

static void
write_contig_file (char *genomesubdir, char *fileroot, Table_T accsegmentpos_table, 
		   Table_chrpos_T chrlength_table, List_T contigtypelist, bool coord_values_8p) {
  FILE *textfp;
  char *textfile, *iitfile, *annot;
  int naccessions, i;
  char **accessions, *divstring, seglength[12]; /* 2^32 = 4*10^9 */
  Segmentpos_T segmentpos;
  Chrom_T chrom;
  Univcoord_T chroffset, universalpos1, universalpos2;
  List_T divlist = NULL, intervallist = NULL, labellist = NULL, annotlist = NULL, p;
  Table_T intervaltable, labeltable, annottable;
  Univinterval_T interval;
#if 0
  void **keys;
  int *values, ntypes;
#endif
  
  if (divsort == NO_SORT) {
    accessions = (char **) Table_keys_by_timeindex(accsegmentpos_table,NULL);
    naccessions = Table_length(accsegmentpos_table);
  } else {
    /* Get accessions in order */
    accessions = (char **) Table_keys(accsegmentpos_table,NULL);
    naccessions = Table_length(accsegmentpos_table);
    current_accsegmentpos_table = accsegmentpos_table;
    qsort(accessions,naccessions,sizeof(char *),bysegmentpos_compare);
  }

#if 0
  /* Get types in order */
  keys = Tableuint_keys(contigtype_table,NULL);
  values = Tableuint_values(contigtype_table,0);
  ntypes = Tableuint_length(contigtype_table);
  contigtypes = (char **) CALLOC(ntypes+1,sizeof(char *)); /* Add 1 for type 0 */
  contigtypes[0] = "";
  for (j = 0; j < ntypes; j++) {
    contigtypes[values[j]] = keys[j];
  }
  FREE(values);
  FREE(keys);
#endif

  /* Write contig text file */
  textfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".contig")+1,sizeof(char));
  sprintf(textfile,"%s/%s.contig",genomesubdir,fileroot);
  /* Use binary, not text, so files are Unix-compatible */
  if ((textfp = FOPEN_WRITE_BINARY(textfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",textfile);
    exit(9);
  }
  FREE(textfile);

  for (i = 0; i < naccessions; i++) {
    segmentpos = (Segmentpos_T) Table_get(accsegmentpos_table,(void *) accessions[i]);
    chrom = Segmentpos_chrom(segmentpos);
#ifdef HAVE_64_BIT
    chroffset = (Univcoord_T) Tableuint8_get(chrlength_table,chrom);
#else
    chroffset = (Univcoord_T) Tableuint_get(chrlength_table,chrom);
#endif
    universalpos1 = chroffset + (Univcoord_T) Segmentpos_chrpos1(segmentpos);
    universalpos2 = chroffset + (Univcoord_T) Segmentpos_chrpos2(segmentpos);

    /* Print as 1-based, inclusive [a,b] */
    if (Segmentpos_revcompp(segmentpos) == true) {
      fprintf(textfp,"%s\t%llu..%llu\t%s:%u..%u\t%u",
	      accessions[i],(unsigned long long) universalpos2+1,(unsigned long long) universalpos1,
	      Chrom_string(chrom),Segmentpos_chrpos2(segmentpos)+1,Segmentpos_chrpos1(segmentpos),
	      Segmentpos_length(segmentpos));
    } else {
      fprintf(textfp,"%s\t%llu..%llu\t%s:%u..%u\t%u",
	      accessions[i],(unsigned long long) universalpos1+1,(unsigned long long) universalpos2,
	      Chrom_string(chrom),Segmentpos_chrpos1(segmentpos)+1,Segmentpos_chrpos2(segmentpos),
	      Segmentpos_length(segmentpos));
    }


#if 0
    if (Segmentpos_type(segmentpos) > 0) {
      fprintf(textfp,"\t%s",contigtypes[Segmentpos_type(segmentpos)]);
    }
#endif

    fprintf(textfp,"\n");

    /* Store as 0-based, inclusive [a,b] */
    labellist = List_push(labellist,(void *) accessions[i]);
    if (Segmentpos_revcompp(segmentpos) == true) {
      intervallist = List_push(intervallist, 
			       (void *) Univinterval_new(universalpos2-1,universalpos1,
							 Segmentpos_type(segmentpos)));

      /* IIT version 1.  Indicate revcomp with '-' as first character, which indicates that the
	 contig was reverse complement */
      sprintf(seglength,"-%u",Segmentpos_length(segmentpos));

    } else {
      intervallist = List_push(intervallist, 
			       (void *) Univinterval_new(universalpos1,universalpos2-1,
							 Segmentpos_type(segmentpos)));

      /* IIT version 1.  Indicate revcomp with '-' as first character */
      sprintf(seglength,"%u",Segmentpos_length(segmentpos));
    }

    annot = (char *) CALLOC(strlen(seglength)+1,sizeof(char));
    strcpy(annot,seglength);
    annotlist = List_push(annotlist,(void *) annot);

  }
  fclose(textfp);

#if 0
  FREE(contigtypes);
#endif
  FREE(accessions);
  intervallist = List_reverse(intervallist);
  /* contigtypelist = List_reverse(contigtypelist); -- Done by caller */ 
  labellist = List_reverse(labellist);
  annotlist = List_reverse(annotlist);

  /* Write contig IIT file */
  divstring = (char *) CALLOC(1,sizeof(char));
  divstring[0] = '\0';
  divlist = List_push(NULL,divstring);

  intervaltable = Table_new(65522,Table_string_compare,Table_string_hash);
  labeltable = Table_new(65522,Table_string_compare,Table_string_hash);
  annottable = Table_new(65522,Table_string_compare,Table_string_hash);

  Table_put(intervaltable,(void *) divstring,intervallist);
  Table_put(labeltable,(void *) divstring,labellist);
  Table_put(annottable,(void *) divstring,annotlist);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);

#if 0
  debug(
	for (p = contigtypelist; p != NULL; p = List_next(p)) {
	  printf("Type %s\n",(char *) List_head(p));
	}
	);
#endif
  IIT_write_univ(iitfile,divlist,contigtypelist,
		 intervaltable,labeltable,annottable,
		 coord_values_8p,/*label_pointers_8p*/false,/*annot_pointers_8p*/false);
  FREE(iitfile);

  List_free(&divlist);
  FREE(divstring);

  Table_free(&annottable);
  Table_free(&labeltable);
  Table_free(&intervaltable);

  for (p = annotlist; p != NULL; p = List_next(p)) {
    annot = (char *) List_head(p);
    FREE(annot);
  }
  List_free(&annotlist);

  /* Labels (accessions) are freed by accsegmentpos_table_gc */
  List_free(&labellist);

  for (p = intervallist; p != NULL; p = List_next(p)) {
    interval = (Univinterval_T) List_head(p);
    Univinterval_free(&interval);
  }
  List_free(&intervallist);

  return;
}

/************************************************************************/


#if 0
static void
stringlist_gc (List_T *list) {
  List_T p;
  char *string;

  for (p = *list; p != NULL; p = List_next(p)) {
    string = (char *) List_head(p);
    FREE(string);
  }
  List_free(&(*list));
  return;
}
#endif


static void
chrlength_table_gc (Table_chrpos_T *chrlength_table) {
  /* Don't free chrom entries in table, because they are freed by Segmentpos_free */
#ifdef HAVE_64_BIT
  Tableuint8_free(&(*chrlength_table));
#else
  Tableuint_free(&(*chrlength_table));
#endif

  return;
}

static void
accsegmentpos_table_gc (Table_T *accsegmentpos_table) {
  int n, i = 0;
  char *accession;
  Segmentpos_T segmentpos;
  void **keys, **values;

  /* For some reason, this was failing on some computers, perhaps
     because we weren't checking for n > 0 */
  if ((n = Table_length(*accsegmentpos_table)) > 0) {
    keys = Table_keys(*accsegmentpos_table,NULL);
    values = Table_values(*accsegmentpos_table,NULL);
    for (i = 0; i < n; i++) {
      accession = (char *) keys[i];
      FREE(accession);
    }
    for (i = 0; i < n; i++) {
      segmentpos = (Segmentpos_T) values[i];
      Segmentpos_free(&segmentpos);
    }
    FREE(values);
    FREE(keys);
  }

  Table_free(&(*accsegmentpos_table));
  return;
}

static void
alt_scaffold_info_table_gc (Table_T *alt_scaffold_info_table) {
  int n, i = 0;
  char *altloc_chr, *altscaffold_info;
  void **keys, **values;

  if ((n = Table_length(*alt_scaffold_info_table)) > 0) {
    keys = Table_keys(*alt_scaffold_info_table,NULL);
    values = Table_values(*alt_scaffold_info_table,NULL);
    for (i = 0; i < n; i++) {
      altloc_chr = (char *) keys[i];
      FREE(altloc_chr);
    }
    for (i = 0; i < n; i++) {
      altscaffold_info = (char *) values[i];
      FREE(altscaffold_info);
    }
    FREE(values);
    FREE(keys);
  }

  Table_free(&(*alt_scaffold_info_table));
  return;
}


#if 0
static char *
remove_slashes (char *buffer) {
  char *copy, *p;

  copy = (char *) CALLOC(strlen(buffer)+1,sizeof(char));
  strcpy(copy,buffer);
  p = copy;
  while (*p != '\0') {
    if (*p == '/') {
      *p = '_';
    }
    p++;
  }
  return copy;
}
#endif


static int
add_compression_type (char *string) {
  if (!strcmp(string,"none")) {
    compression_types = NO_COMPRESSION;
    return 0;
  } else if (!strcmp(string,"all")) {
    compression_types = BITPACK64_COMPRESSION;
    return 1;
  } else {
    if (!strcmp(string,"bitpack64")) {
      compression_types |= BITPACK64_COMPRESSION;
    } else {
      fprintf(stderr,"Don't recognize compression type %s\n",string);
      fprintf(stderr,"Allowed values are: none, all, bitpack64\n");
      exit(9);
    }
    return 1;
  }
}



#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int 
main (int argc, char *argv[]) {
  int ncontigs;
  Table_T accsegmentpos_table, alt_scaffold_info_table;

  FILE *fp;
  char *key, **keys, *chrname, *chrname_alt, *line;
  int line_length;
  Tableuint_T chrorder_table = NULL;
  unsigned int order;
  int i;

  Table_chrpos_T chrlength_table;
  List_T contigtypelist = NULL, p;
  Univ_IIT_T chromosome_iit, contig_iit;
  char *typestring;
  Univcoord_T genomelength, totalnts;
  char *chromosomefile, *iitfile, *positionsfile_high, *positionsfile_low, interval_char, local_interval_char;

  Indexdb_filenames_T indexdb_filenames;
  Chrpos_T seglength;
  bool coord_values_8p;




#ifdef HAVE_64_BIT
  UINT8 noffsets;
#endif

  int c;
  extern int optind;
  extern char *optarg;
  char *string;

  while ((c = getopt(argc,argv,"F:D:d:z:k:q:j:ArlGUNHOPQWw:e:Ss:n:m9")) != -1) {
    switch (c) {
    case 'F': sourcedir = optarg; break;
    case 'D': destdir = optarg; break;
    case 'd': fileroot = optarg; break;

    case 'z':
      compression_types = NO_COMPRESSION; /* Initialize */
      string = strtok(optarg,",");
      if (add_compression_type(string) != 0) {
	while ((string = strtok(NULL,",")) != NULL && add_compression_type(string) != 0) {
	}
      }
      break;

    case 'k': index1part = atoi(optarg);
      if (index1part > MAXIMUM_KMER) {
	fprintf(stderr,"The choice of k-mer size must be %d or less\n",MAXIMUM_KMER);
	exit(9);
      }
      break;
    case 'q': index1interval = atoi(optarg); break;

    case 'j': local1part = atoi(optarg);
#if 0
      if (local1part > MAXIMUM_LOCAL1PART) {
	fprintf(stderr,"The choice of local size must be %d or less\n",MAXIMUM_LOCAL1PART);
	exit(9);
      }
#endif
      break;

    case 'A': action = AUXFILES; break;
    case 'r': rawp = true; break;
    case 'l': genome_lc_p = true; break;
    case 'G': action = GENOME; break;
    case 'U': action = UNSHUFFLE; break;
    case 'N': action = COUNT; break;
    case 'H': huge_offsets_p = true; break;

    case 'O': action = OFFSETS; break;
    case 'P': action = POSITIONS; break;

    case 'Q': action = LOCALDB; break;

    case 'W': writefilep = true; break;
    case 'w': wraplength = atoi(optarg); break;
    case 'e': nmessages = atoi(optarg); break;

    case 's': 
      if (!strcmp(optarg,"none")) {
	divsort = NO_SORT;
      } else if (!strcmp(optarg,"alpha")) {
	divsort = ALPHA_SORT;
      } else if (!strcmp(optarg,"numeric-alpha")) {
	divsort = NUMERIC_ALPHA_SORT;
      } else if (!strcmp(optarg,"chrom")) {
	divsort = CHROM_SORT;
      } else if (!strcmp(optarg,"names")) {
	divsort = FILENAME_SORT;
      } else {
	fprintf(stderr,"Don't recognize sort type %s.  Allowed values are none, alpha, or chrom.",optarg);
	exit(9);
      }
      break;

    case 'n': sortfilename = optarg; break;

    case 'm': mask_lowercase_p = true; break;

    case '9': /* check_compiler_assumptions(); */ return 0; break;

    default: fprintf(stderr,"Unknown flag %c\n",c); exit(9);

    }
  }
  argc -= optind;
  argv += optind;

  if (index1interval == 3) {
    interval_char = '3';
  } else if (index1interval == 2) {
    interval_char = '2';
  } else if (index1interval == 1) {
    interval_char = '1';
  } else {
    fprintf(stderr,"Selected indexing interval %d is not allowed.  Only values allowed are 3, 2, or 1\n",index1interval);
    exit(9);
  }

  if (local1interval == 3) {
    local_interval_char = '3';
  } else if (local1interval == 2) {
    local_interval_char = '2';
  } else if (local1interval == 1) {
    local_interval_char = '1';
  } else {
    fprintf(stderr,"Selected indexing interval %d is not allowed.  Only values allowed are 3, 2, or 1\n",local1interval);
    exit(9);
  }

  if (action != UNSHUFFLE) {
    if (fileroot == NULL) {
      fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
      exit(9);
    }
  }

  if (action == AUXFILES) {
    /* Usage: cat <fastafile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -A
       Requires <fastafile> in appropriate format
       Writes <destdir>/<dbname>.chromosome and <destdir>/<dbname>.contig files 
       and corresponding .iit files */

    if (divsort == FILENAME_SORT) {
      if (sortfilename == NULL) {
	fprintf(stderr,"For sorting by names file, need to provide file to -n flag");
	exit(9);
      } else if ((fp = fopen(sortfilename,"r")) == NULL) {
	fprintf(stderr,"Unable to open file %s provided to -n flag",sortfilename);
	exit(9);
      } else {
	chrorder_table = Tableuint_new(65522,Table_string_compare,Table_string_hash);
	order = 1;
	while ((line = Getline_wlength(&line_length,fp)) != NULL) {
	  chrname = (char *) MALLOC((line_length+1)*sizeof(char));
	  chrname_alt = (char *) MALLOC((line_length+1)*sizeof(char));
   	  if (sscanf(line,"%s %s",chrname,chrname_alt) == 2) {	  
	    key = (char *) CALLOC(strlen(chrname_alt)+1,sizeof(char));
	    strcpy(key,chrname_alt);
	    Tableuint_put(chrorder_table,(void *) key,order);
          } else if (sscanf(line,"%s",chrname) == 1) {
	    key = (char *) CALLOC(strlen(chrname)+1,sizeof(char));
	    strcpy(key,chrname);
	    Tableuint_put(chrorder_table,(void *) key,order);
          } else {
	    fprintf(stderr,"Unable to parse line %s\n",line);
	  }
	  FREE(chrname_alt);
	  FREE(chrname);
	  FREE(line);

          order += 1;
	}
	fclose(fp);
      }
    }

    if (getc(stdin) != '>') {
      fprintf(stderr,"Expected file to start with '>'\n");
      exit(9);
    }

    /* Holds contigs.  keys are strings; values are structs. */
    accsegmentpos_table = Table_new(65522,Table_string_compare,Table_string_hash);
    /* Hold chromosomes.  keys are Chrom_Ts; values are uints. */
#ifdef HAVE_64_BIT
    chrlength_table = Tableuint8_new(65522,Chrom_compare_table,Chrom_hash_table);
#else
    chrlength_table = Tableuint_new(65522,Chrom_compare_table,Chrom_hash_table);
#endif
    alt_scaffold_info_table = Table_new(65522,Table_string_compare,Table_string_hash);


#if 0
    /* No longer supporting strains */
    /* keys are strings; values are ints */
    contigtype_table = Tableuint_new(100,Table_string_compare,Table_string_hash);

    refstrain = read_strain_from_coordsfile(coordsfile);
    fprintf(stderr,"Reference strain is %s\n",refstrain);
    contigtypelist = List_push(NULL,refstrain);
#endif
    /* The zeroth type is empty */
    typestring = (char *) CALLOC(1,sizeof(char));
    typestring[0] = '\0';
    contigtypelist = List_push(NULL,typestring);

    ncontigs = 0;
    totalnts = 0U;

    while (process_sequence_aux(&seglength,accsegmentpos_table,chrlength_table,
				chrorder_table,alt_scaffold_info_table,
				fileroot,ncontigs) == true) {
      if (totalnts + seglength < totalnts) {
	/* Exceeds 32 bits */
	fprintf(stderr,"The total length of genomic sequence exceeds 2^32 = 4,294,967,296 bp, which the GMAP index format cannot handle\n");
	exit(9);
      } else {
	totalnts = totalnts + seglength;
      }
      ncontigs++;
    }
    fprintf(stderr,"Total genomic length = %llu bp\n",(unsigned long long) totalnts);

    if (ncontigs == 0) {
      fprintf(stderr,"No contig information was provided to gmapindex\n");
      exit(9);
    }

#ifdef HAVE_64_BIT
    if (totalnts > 4294967295) {
      coord_values_8p = true;
    } else {
      coord_values_8p = false;
    }
#else
    coord_values_8p = false;
#endif

    write_chromosome_file(destdir,fileroot,chrlength_table,coord_values_8p);
    write_contig_file(destdir,fileroot,accsegmentpos_table,chrlength_table,contigtypelist,coord_values_8p);

    /* Note: We are using destdir here, not sourcedir */
    chromosomefile = (char *) CALLOC(strlen(destdir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",destdir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    write_alt_scaffold_file(destdir,fileroot,alt_scaffold_info_table,chromosome_iit,coord_values_8p);
    Univ_IIT_free(&chromosome_iit);
    FREE(chromosomefile);

    for (p = contigtypelist; p != NULL; p = List_next(p)) {
      typestring = (char *) List_head(p);
      FREE(typestring);
    }
    List_free(&contigtypelist);

    chrlength_table_gc(&chrlength_table);
    accsegmentpos_table_gc(&accsegmentpos_table);

    alt_scaffold_info_table_gc(&alt_scaffold_info_table);

    if (chrorder_table != NULL) {
      keys = (char **) Tableuint_keys(chrorder_table,NULL);
      for (i = 0; i < Tableuint_length(chrorder_table); i++) {
	FREE(keys[i]);
      }
      FREE(keys);
      Tableuint_free(&chrorder_table);
    }

  } else if (action == GENOME) {
    /* Usage: cat <fastafile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -G
       Requires <fastafile> in appropriate format and <sourcedir>/<dbname>.chromosome.iit 
       and <sourcedir>/<dbname>.contig.iit files.
       Creates <destdir>/<dbname>.genomecomp (horizontal format in blocks of 32).
       Then the unshuffle command creates another version (vertical format in blocks of 128).  */

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
    FREE(chromosomefile);

    iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			      strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.contig.iit",sourcedir,fileroot);
    if ((contig_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",iitfile);
      exit(9);
    }
    FREE(iitfile);
    
    if (Univ_IIT_ntypes(contig_iit) == 1) {
      /* index1part needed only if writing an uncompressed genome using a file */
      Genome_write_comp32(destdir,fileroot,stdin,contig_iit,
			  chromosome_iit,genome_lc_p,rawp,writefilep,
			  genomelength,index1part,nmessages);
    } else if (Univ_IIT_ntypes(contig_iit) > 1) {
      fprintf(stderr,"GMAPINDEX no longer supports alternate strains\n");
      abort();
    }

    Univ_IIT_free(&chromosome_iit);
    Univ_IIT_free(&contig_iit);

  } else if (action == UNSHUFFLE) {
    /* Usage: cat <fastafile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -G
       Requires <fastafile> in appropriate format and <sourcedir>/<dbname>.chromosome.iit 
       and <sourcedir>/<dbname>.contig.iit files.
       Creates <destdir>/<dbname>.genomebits128 (in blocks of 128) */

    Compress_unshuffle_bits128(stdout,stdin);

  } else if (action == COUNT) {
    /* Usage: cat <genomefile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -C */

#ifndef HAVE_64_BIT
    printf("0\n");
#else
    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    noffsets = Indexdb_count_offsets(stdin,chromosome_iit,index1part,index1interval,
				     genome_lc_p,fileroot,mask_lowercase_p);
    printf("%llu\n",(unsigned long long) noffsets);

    Univ_IIT_free(&chromosome_iit);
#endif

  } else if (action == OFFSETS) {
    /* Usage: gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -O <genomefile>
       Creates <destdir>/<dbname>.ref153offsets64meta and ref153offsets64strm */

    if (argc == 0) {
      fp = stdin;
    } else if ((fp = fopen(argv[0],"rb")) == NULL) {
      fprintf(stderr,"Could not open file %s\n",argv[0]);
      exit(9);
    }

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    if (huge_offsets_p == false) {
      fprintf(stderr,"Offset compression types:");
      if ((compression_types & BITPACK64_COMPRESSION) != 0) {
	fprintf(stderr," bitpack64");
      }
      fprintf(stderr,"\n");
      
      Indexdb_write_offsets(destdir,interval_char,fp,chromosome_iit,
			    index1part,index1interval,
			    genome_lc_p,fileroot,mask_lowercase_p);
    } else {
      fprintf(stderr,"Offset compression types:");
      if ((compression_types & BITPACK64_COMPRESSION) != 0) {
	fprintf(stderr," bitpack64");
      }
      fprintf(stderr,"\n");
      
      Indexdb_write_offsets_huge(destdir,interval_char,fp,chromosome_iit,
				 index1part,index1interval,
				 genome_lc_p,fileroot,mask_lowercase_p);
    }

    if (argc > 0) {
      fclose(fp);
    }

    Univ_IIT_free(&chromosome_iit);

  } else if (action == POSITIONS) {
    /* Usage: gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -P <genomefile>
       Requires <sourcedir>/<dbname>.ref153offsets64meta and .ref153offsets64strm
       Creates <destdir>/<dbname>.ref153positions */

    if (argc == 0) {
      fp = stdin;
    } else if ((fp = fopen(argv[0],"rb")) == NULL) {
      fprintf(stderr,"Could not open file %s\n",argv[0]);
      exit(9);
    }

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
    FREE(chromosomefile);

    indexdb_filenames = Indexdb_get_filenames(&compression_type,&index1part,&index1interval,
					      sourcedir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
					      /*required_index1part*/index1part,
					      /*required_interval*/index1interval,/*offsets_only_p*/true);

    if (Univ_IIT_coord_values_8p(chromosome_iit) == true) {
      coord_values_8p = true;

      positionsfile_high = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen(IDX_FILESUFFIX)+
					   /*for kmer*/2+/*for interval char*/1+
					   strlen(POSITIONS_HIGH_FILESUFFIX)+1,sizeof(char));
      sprintf(positionsfile_high,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,POSITIONS_HIGH_FILESUFFIX);

      positionsfile_low = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen(IDX_FILESUFFIX)+
					   /*for kmer*/2+/*for interval char*/1+
					   strlen(POSITIONS_FILESUFFIX)+1,sizeof(char));
      sprintf(positionsfile_low,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,POSITIONS_FILESUFFIX);

    } else {
      coord_values_8p = false;
      
      positionsfile_high = (char *) NULL;

      positionsfile_low = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen(IDX_FILESUFFIX)+
					   /*for kmer*/2+/*for interval char*/1+
					   strlen(POSITIONS_FILESUFFIX)+1,sizeof(char));
      sprintf(positionsfile_low,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,POSITIONS_FILESUFFIX);
    }

    if (huge_offsets_p == false) {
      Indexdb_write_positions(positionsfile_high,positionsfile_low,indexdb_filenames->pointers_filename,
			      indexdb_filenames->offsets_filename,fp,chromosome_iit,
			      index1part,index1interval,genomelength,
			      genome_lc_p,writefilep,fileroot,mask_lowercase_p,
			      coord_values_8p);
    } else {
      Indexdb_write_positions_huge(positionsfile_high,positionsfile_low,
				   indexdb_filenames->pages_filename,indexdb_filenames->pointers_filename,
				   indexdb_filenames->offsets_filename,fp,chromosome_iit,
				   index1part,index1interval,genomelength,
				   genome_lc_p,writefilep,fileroot,mask_lowercase_p,
				   coord_values_8p);
    }

    if (argc > 0) {
      fclose(fp);
    }

    Indexdb_filenames_free(&indexdb_filenames);

    FREE(positionsfile_high);
    FREE(positionsfile_low);
    Univ_IIT_free(&chromosome_iit);

  } else if (action == LOCALDB) {
    /* Usage: gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -Q <genomefile>
       Creates <destdir>/<dbname>.ref061locoffsets64table, ref061locoffsets64meta, ref061locoffsets64strm, ref061locpositions */

    if (argc == 0) {
      fp = stdin;
    } else if ((fp = fopen(argv[0],"rb")) == NULL) {
      fprintf(stderr,"Could not open file %s\n",argv[0]);
      exit(9);
    }

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(chromosomefile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",chromosomefile);
      exit(9);
    }
    FREE(chromosomefile);

    if (huge_offsets_p == false) {
      fprintf(stderr,"Offset compression types:");
      if ((compression_types & BITPACK64_COMPRESSION) != 0) {
	fprintf(stderr," bitpack64");
      }
      fprintf(stderr,"\n");
      
      Localdb_write(destdir,local_interval_char,fp,chromosome_iit,
		    local1part,local1interval,
		    genome_lc_p,fileroot,mask_lowercase_p);
    } else {
#if 0
      /* TODO: Implement */
      fprintf(stderr,"Offset compression types:");
      if ((compression_types & BITPACK64_COMPRESSION) != 0) {
	fprintf(stderr," bitpack64");
      }
      fprintf(stderr,"\n");
      
      Localdb_write_huge(destdir,local_interval_char,fp,chromosome_iit,
			 local1part,local1interval,
			 genome_lc_p,fileroot,mask_lowercase_p);
#endif
    }

    if (argc > 0) {
      fclose(fp);
    }

    Univ_IIT_free(&chromosome_iit);
  }

  return 0;
}

