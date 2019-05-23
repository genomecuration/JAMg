static char rcsid[] = "$Id: transcriptome.c 212997 2018-02-02 06:38:09Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "transcriptome.h"

#include <stdio.h>
#include <stddef.h>
#include <string.h>		/* Needed for strlen and strcpy */

#include "mem.h"
#include "access.h"
#include "list.h"
#include "genomicpos.h"
#include "bitpack64-readtwo.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif



#define T Transcriptome_T
struct T {

  Chrnum_T *chrnums;

  Access_T offsetsmeta_access;
  int offsetsmeta_shmid;
  key_t offsetsmeta_key;
  int offsetsmeta_fd;
  size_t offsetsmeta_len;
  UINT4 *offsetsmeta;

  Access_T offsetsstrm_access;
  int offsetsstrm_shmid;
  key_t offsetsstrm_key;
  int offsetsstrm_fd;
  size_t offsetsstrm_len;
  UINT4 *offsetsstrm;

  Access_T exoninfo_access;
  int exoninfo_shmid;
  key_t exoninfo_key;
  int exoninfo_fd;
  size_t exoninfo_len;

  /* exoninfo has exonbounds (int) and exonstarts (unsigned int) interleaved */
  unsigned int *exoninfo;
};


#define BUFFERSIZE 1024000

void
Transcriptome_free (T *old) {
  if (*old) {

    if ((*old)->exoninfo_access == ALLOCATED_PRIVATE) {
      FREE((*old)->exoninfo);
    } else if ((*old)->exoninfo_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->exoninfo,(*old)->exoninfo_shmid,(*old)->exoninfo_key);
    } else {
      abort();
    }

    if ((*old)->offsetsstrm_access == ALLOCATED_PRIVATE) {
      FREE((*old)->offsetsstrm);
    } else if ((*old)->offsetsstrm_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->offsetsstrm,(*old)->offsetsstrm_shmid,(*old)->offsetsstrm_key);
    } else {
      abort();
    }
      
    if ((*old)->offsetsmeta_access == ALLOCATED_PRIVATE) {
      FREE((*old)->offsetsmeta);
    } else if ((*old)->offsetsmeta_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->offsetsmeta,(*old)->offsetsmeta_shmid,(*old)->offsetsmeta_key);
    } else {
      abort();
    }

    FREE((*old)->chrnums);
    FREE(*old);
  }
  return;
}


T
Transcriptome_new (char *genomesubdir, char *fileroot, Univ_IIT_T chromosome_iit,
		   bool sharedp) {
  T new = (T) MALLOC(sizeof(*new));
  char Buffer[BUFFERSIZE];
  int len;

  FILE *fp;
  char *filename, *divstring, **divs;
  List_T divlist = NULL;
  Chrnum_T chrnum;
  bool *baddivp;

  char *comma;
  double seconds;
  int ntranscripts, trnum, i;
  int ndivs, divint;


  if (chromosome_iit == NULL) {
    /* Read transcript chromosome numbers independently of genome */

    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+strlen(".tr.divints")+1,sizeof(char));
    sprintf(filename,"%s/%s.tr.divints",genomesubdir,fileroot);
    ntranscripts = Access_filesize(filename)/sizeof(int);

    fp = fopen(filename,"rb");
    FREE(filename);

    new->chrnums = (Chrnum_T *) MALLOC((ntranscripts+1)*sizeof(Chrnum_T));
    new->chrnums[0] = 0;
    for (trnum = 1; trnum <= ntranscripts; trnum++) {
      FREAD_INT(&divint,fp);
      new->chrnums[trnum] = divint;
    }
    fclose(fp);


  } else {
    /* Read transcript chromosome numbers in context of genome chromosome numbers */
    /* Read divstrings */
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+strlen(".tr.divs")+1,sizeof(char));
    sprintf(filename,"%s/%s.tr.divs",genomesubdir,fileroot);
    fp = fopen(filename,"r");
    FREE(filename);

    ndivs = 0;
    while (fgets(Buffer,BUFFERSIZE,fp) != NULL) {
      len = strlen(Buffer);
      if (Buffer[len-1] != '\n') {
	fprintf(stderr,"Buffer overflow on reading div %s\n",Buffer);
	exit(9);
      } else {
	Buffer[len-1] = '\0';
      }
      divstring = (char *) MALLOC(len*sizeof(char));
      strcpy(divstring,Buffer);
      divlist = List_push(divlist,(void *) divstring);
      ndivs++;
    }
    divlist = List_reverse(divlist);
    divs = (char **) List_to_array(divlist,NULL);
    List_free(&divlist);
    fclose(fp);


    /* Translate signed divints to signed chrnums */
    baddivp = (bool *) CALLOC(ndivs,sizeof(bool));

    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+strlen(".tr.divints")+1,sizeof(char));
    sprintf(filename,"%s/%s.tr.divints",genomesubdir,fileroot);
    ntranscripts = Access_filesize(filename)/sizeof(int);

    fp = fopen(filename,"rb");
    FREE(filename);

    new->chrnums = (Chrnum_T *) MALLOC((ntranscripts+1)*sizeof(Chrnum_T));
    new->chrnums[0] = 0;
    for (trnum = 1; trnum <= ntranscripts; trnum++) {
      FREAD_INT(&divint,fp);
      if (divint > 0) {
	if ((chrnum = Univ_IIT_find_one(chromosome_iit,divs[divint-1])) < 0) {
	  if (baddivp[divint-1] == false) {
	    fprintf(stderr,"Cannot find genomic chromosome for divstring %s\n",divs[divint-1]);
	    baddivp[divint-1] = true;
	  }
	  new->chrnums[trnum] = 0;
	} else {
	  new->chrnums[trnum] = chrnum;
	}
      } else {
	divint = -divint;
	if ((chrnum = Univ_IIT_find_one(chromosome_iit,divs[divint-1])) < 0) {
	  if (baddivp[divint-1] == false) {
	    fprintf(stderr,"Cannot find genomic chromosome for divstring %s\n",divs[divint-1]);
	    baddivp[divint-1] = true;
	  }
	  new->chrnums[trnum] = 0;
	} else {
	  new->chrnums[trnum] = -chrnum;
	}
      }
      /* printf("For trnum %d, %s, assigned chrnum %d\n",trnum,divs[divint-1],new->chrnums[trnum]); */
    }
    FREE(baddivp);
    fclose(fp);

    for (i = 0; i < ndivs; i++) {
      FREE(divs[i]);
    }
    FREE(divs);
  }


  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+strlen(".tr.offsets64meta")+1,sizeof(char));
  sprintf(filename,"%s/%s.tr.offsets64meta",genomesubdir,fileroot);
  if (sharedp == true) {
    new->offsetsmeta = (UINT4 *) Access_allocate_shared(&new->offsetsmeta_access,&new->offsetsmeta_shmid,&new->offsetsmeta_key,
							&new->offsetsmeta_fd,&new->offsetsmeta_len,&seconds,
							filename,sizeof(UINT4));
  } else {
    new->offsetsmeta = (UINT4 *) Access_allocate_private(&new->offsetsmeta_access,&new->offsetsmeta_len,&seconds,
							filename,sizeof(UINT4));
  }
  comma = Genomicpos_commafmt(new->offsetsmeta_len);
  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
  FREE(comma);
  FREE(filename);



  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+strlen(".tr.offsets64strm")+1,sizeof(char));
  sprintf(filename,"%s/%s.tr.offsets64strm",genomesubdir,fileroot);
  if (sharedp == true) {
    new->offsetsstrm = (UINT4 *) Access_allocate_shared(&new->offsetsstrm_access,&new->offsetsstrm_shmid,&new->offsetsstrm_key,
							&new->offsetsstrm_fd,&new->offsetsstrm_len,&seconds,
							filename,sizeof(UINT4));
  } else {
    new->offsetsstrm = (UINT4 *) Access_allocate_private(&new->offsetsstrm_access,&new->offsetsstrm_len,&seconds,
							 filename,sizeof(UINT4));
  }
  comma = Genomicpos_commafmt(new->offsetsstrm_len);
  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
  FREE(comma);
  FREE(filename);



  /* Read exoninfo file */
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+strlen(".tr.exoninfo")+1,sizeof(char));
  sprintf(filename,"%s/%s.tr.exoninfo",genomesubdir,fileroot);
  if (sharedp == true) {
    new->exoninfo = (UINT4 *) Access_allocate_shared(&new->exoninfo_access,&new->exoninfo_shmid,&new->exoninfo_key,
						     &new->exoninfo_fd,&new->exoninfo_len,&seconds,
						     filename,sizeof(UINT4));
  } else {
    new->exoninfo = (UINT4 *) Access_allocate_private(&new->exoninfo_access,&new->exoninfo_len,&seconds,
						      filename,sizeof(UINT4));
  }
  comma = Genomicpos_commafmt(new->exoninfo_len);
  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
  FREE(comma);
  FREE(filename);
    

  return new;
}


Chrnum_T
Transcriptome_chrnum (int *transcript_genestrand, T this, int trnum) {
  Chrnum_T chrnum;

  if ((chrnum = this->chrnums[trnum]) == 0) {
    return 0;
  } else if (chrnum > 0) {
    *transcript_genestrand = +1;
    return chrnum;
  } else {
    *transcript_genestrand = -1;
    return -chrnum;
  }
}



/* Assumes that we have already checked that this->chrnums != 0 */
int
Transcriptome_exons (int **exonbounds, Chrpos_T **exonstarts, T this, int trnum) {
  UINT4 offset0, offset1;
  int nexons;

  offset0 = Bitpack64_read_two(&offset1,(Oligospace_T) (trnum - 1),this->offsetsmeta,this->offsetsstrm);
  nexons = (int) (offset1 - offset0);
  *exonbounds = &(((int *) this->exoninfo)[2*offset0]);
  *exonstarts = &(this->exoninfo[2*offset0 + nexons]);

  return nexons;
}


bool
Transcriptome_genomic_bounded_p (int trnum, Chrnum_T chrbound, Chrpos_T lowbound, Chrpos_T highbound, T this) {
  UINT4 offset0, offset1;
  int nexons;
  int *exonbounds;
  Chrnum_T chrnum;
  Chrpos_T *exonstarts;
  Univcoord_T start, end;

  if ((chrnum = this->chrnums[trnum]) == 0) {
    return false;
  } else if (chrnum > 0) {
    if (chrnum != chrbound) {
      return false;
    }
  } else {
    if (-chrnum != chrbound) {
      return false;
    }
  }

  offset0 = Bitpack64_read_two(&offset1,(Oligospace_T) (trnum - 1),this->offsetsmeta,this->offsetsstrm);
  nexons = (int) (offset1 - offset0);

  exonbounds = &(((int *) this->exoninfo)[2*offset0]);
  exonstarts = &(this->exoninfo[2*offset0 + nexons]);

  start = exonstarts[0];
  end = exonstarts[nexons - 1] + exonbounds[nexons - 1];
  if (start < end) {
    if (end < lowbound) {
      return false;
    } else if (start > highbound) {
      return false;
    } else {
      return true;
    }

  } else {
    if (start < lowbound) {
      return false;
    } else if (end > highbound) {
      return false;
    } else {
      return true;
    }
  }
}

