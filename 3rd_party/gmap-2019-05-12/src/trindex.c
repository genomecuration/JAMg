static char rcsid[] = "$Id: trindex.c 215170 2018-05-12 00:44:00Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>		/* For open */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include <ctype.h>		/* For toupper */
#include <sys/mman.h>		/* For munmap */
#include <math.h>		/* For qsort */
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#include "mem.h"
#include "fopen.h"
#include "access.h"
#include "types.h"
#include "iit-read.h"
#include "bitpack64-write.h"

#include "datadir.h"
#include "getopt.h"


/* Creates the following files: */
/*   <destdir>/<dbroot>.tr.divs -- text file */
/*   <destdir>/<dbroot>.tr.divints */
/*   <destdir>/<dbroot>.tr.exoninfo */
/*   <destdir>/<dbroot>.tr.offsets64meta */
/*   <destdir>/<dbroot>.tr.offsets64strm */


static struct option long_options[] = {
  /* Input options */
  {"transcriptdir", required_argument, 0, 'C'},	/* user_destdir */
  {"transcriptdb", required_argument, 0, 'c'}, /* transcriptome_dbroot */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void print_program_usage ();


static char *user_destdir = NULL;
static char *transcriptome_dbroot = NULL;
static char *transcriptome_fileroot = NULL;
static char *genomesubdir = NULL;


#define BUF_SIZE 1024

int
main (int argc, char *argv[]) {
  char *destdir = NULL, *dbversion;
  FILE *fp;
  char *pointersfile, *offsetsfile, *filename;
  char *divstring;
  int *divints, divint;
  unsigned int *offsets;
  
  char *genes_file;
  IIT_T genes_iit;
  int in, out;
  char buf[BUF_SIZE];
  int nread;

  int ntranscripts, indx;
  int transcript_genestrand;
  int nexons;
  int *exonbounds;
  Chrpos_T *exonstarts;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;


  while ((opt = getopt_long(argc,argv,"C:c:",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;

      if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
	
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'trindex --help'",long_name);
	exit(9);
      }
      break;

    case 'C': user_destdir = optarg; break;
    case 'c': transcriptome_dbroot = optarg; break;

    case '?': fprintf(stderr,"For usage, run 'trindex --help'\n"); exit(9);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;

  
  genomesubdir = Datadir_find_genomesubdir(&transcriptome_fileroot,&dbversion,user_destdir,
					   transcriptome_dbroot);
  FREE(dbversion);

  if (user_destdir == NULL) {
    destdir = genomesubdir;
  } else {
    destdir = user_destdir;
  }

  if (argc <= 0) {
    fprintf(stderr,"Need to specify a genes IIT file\n");
    exit(9);
  }

  genes_file = argv[0];
  if ((genes_iit = IIT_read(genes_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
    fprintf(stderr,"Reading genes file %s...",genes_file);
  } else {
#if 0
    mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,genome_fileroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(genes_file)+1,sizeof(char));
    sprintf(iitfile,"%s/%s",mapdir,genes_file);
    if ((genes_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
      fprintf(stderr,"Reading genes file %s...",iitfile);
      FREE(iitfile);
      FREE(mapdir);
    } else {
      fprintf(stderr,"Genes file %s.iit not found locally or in %s.  Available files:\n",genes_file,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s or specify a full directory path\n",genes_file);
      exit(9);
    }
#endif
  }


  /* Copy genes file back out */
  in = open(genes_file, O_RDONLY);

  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(transcriptome_fileroot)+
			     strlen(".tr.iit")+1,sizeof(char));
  sprintf(filename,"%s/%s.tr.iit",destdir,transcriptome_fileroot);
  out = open(filename,/*openFlags*/O_CREAT | O_WRONLY | O_TRUNC,
	     /*filePerms*/S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);

  while ((nread = read(in,buf,BUF_SIZE)) > 0) {
    if (write(out,buf,nread) != nread) {
      fprintf(stderr,"Could not write whole buffer\n");
      abort();
    }
  }

  close(out);
  close(in);
  FREE(filename);


  /* Write divs as a text file */
  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(transcriptome_fileroot)+
			     strlen(".tr.divs")+1,sizeof(char));
  sprintf(filename,"%s/%s.tr.divs",destdir,transcriptome_fileroot);
  fp = fopen(filename,"w");
  FREE(filename);

  for (divint = 1; divint < IIT_ndivs(genes_iit); divint++) {
    divstring = IIT_divstring(genes_iit,divint);
    fprintf(fp,"%s\n",divstring);
  }
  fclose(fp);


  /* Write exoninfo and compute offsets */
  ntranscripts = IIT_total_nintervals(genes_iit);
  divints = MALLOC(ntranscripts*sizeof(int));

  offsets = MALLOC((ntranscripts+1)*sizeof(unsigned int));
  offsets[0] = 0;

  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(transcriptome_fileroot)+
				 strlen(".tr.exoninfo")+1,sizeof(char));
  sprintf(filename,"%s/%s.tr.exoninfo",destdir,transcriptome_fileroot);
  fp = fopen(filename,"wb");
  FREE(filename);

  for (indx = 1; indx <= ntranscripts; indx++) {
    nexons = IIT_gene_exons_array(&transcript_genestrand,&divint,&exonbounds,&exonstarts,genes_iit,indx);
    FWRITE_INTS(exonbounds,nexons,fp);
    FWRITE_UINTS(exonstarts,nexons,fp);
    FREE(exonstarts);
    FREE(exonbounds);

    if (transcript_genestrand > 0) {
      divints[indx-1] = divint;
    } else {
      divints[indx-1] = -divint;
    }
    offsets[indx] = offsets[indx-1] + (unsigned int) nexons;
  }
  fclose(fp);


  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(transcriptome_fileroot)+
				 strlen(".tr.divints")+1,sizeof(char));
  sprintf(filename,"%s/%s.tr.divints",destdir,transcriptome_fileroot);
  fp = fopen(filename,"wb");
  FREE(filename);

  FWRITE_INTS(divints,ntranscripts,fp);
  fclose(fp);
  FREE(divints);



  /* Write offsets */
  pointersfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(transcriptome_fileroot)+
				 strlen(".tr.offsets64meta")+1,sizeof(char));
  sprintf(pointersfile,"%s/%s.tr.offsets64meta",destdir,transcriptome_fileroot);
  
  offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(transcriptome_fileroot)+
				strlen(".tr.offsets64strm")+1,sizeof(char));
  sprintf(offsetsfile,"%s/%s.tr.offsets64strm",destdir,transcriptome_fileroot);

  fprintf(stderr,"Writing %d offsets compressed via bitpack64...",ntranscripts+1);
  Bitpack64_write_differential(pointersfile,offsetsfile,offsets,ntranscripts);
  fprintf(stderr,"done\n");
  FREE(offsetsfile);
  FREE(pointersfile);

  FREE(offsets);

  IIT_free(&genes_iit);

  fprintf(stderr,"Wrote transcriptome files to %s/%s.*\n",destdir,transcriptome_fileroot);
  FREE(transcriptome_fileroot);
  FREE(genomesubdir);

  return 0;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: trindex [OPTIONS...] -c <transcriptome> <genes IIT file>\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Options (must include -c)\n");
  fprintf(stdout,"\
  -C, --transcriptdir=directory  Directory where to write transcriptome index files (default is\n\
                                   GMAP genome directory specified at compile time)\n\
  -c, --transcriptdb=STRING      Transcriptome database\n\
  --help                         Show this help message\n\
");

  fprintf(stdout,"\n");
  fprintf(stdout,"\
  Note: Before running trindex, you need a genes IIT file, which gives the transcript coordinates\n\
  for some genome index <genome>.  You can use the utility programs provided in the GMAP package:\n\
  ensembl_genes, gff3_genes, gtf_genes, or psl_genes to generate this file, or you can align\n\
  transcripts using GMAP --format=map_exons.  You then cat this file to iit_store -o <filename>\n\
  to generate the genes IIT file\n\
\n\
  You also need to create a GMAP index for the transcripts in the genes IIT file.  If you already\n\
  have the transcripts in a FASTA file, you can do\n\
\n\
      gmap_build -d <transcriptome> -q 1 <fasta file>\n\
\n\
  But if you don't have the transcripts, but only the genes IIT file, you can generate the\n\
  FASTA file by doing\n\
\n\
      get-genome  --sequence --dump -d <genome> -m <genes IIT file> > <fasta file>\n\
\n\
  and then doing\n\
\n\
      gmap_build -d <transcriptome> -q 1 <fasta file>\n\
\n\
  Once you have built the GMAP index for the transcriptome, you can run this program by doing\n\
\n\
      trindex -c <transcriptome> <genes IIT file>\n\
\n\
  Finally, you can then perform transcriptome-guided genomic alignment by doing\n\
\n\
      gsnap -c <transcriptome> -d <genome> <reads>\n\
\n\
");
  
  return;
}
