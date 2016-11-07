static char rcsid[] = "$Id: indexdb-write.c 184465 2016-02-18 00:09:34Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "indexdb-write.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include <ctype.h>		/* For toupper */
#include <sys/mman.h>		/* For munmap */

#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For lseek and close */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
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
#include "types.h"		/* For Oligospace_T */

#include "compress-write.h"	/* For Compress_get_char */
#include "interval.h"
#include "complement.h"
#include "access.h"
#include "genomicpos.h"
#include "bool.h"
#include "indexdbdef.h"
#include "iit-read-univ.h"
#include "indexdb.h"

#include "bitpack64-read.h"
#include "bitpack64-write.h"
#include "bitpack64-access.h"
#include "bitpack64-incr.h"


#define MAX_BITPACK_BLOCKSIZE 64

/* Another MONITOR_INTERVAL is in compress.c */
#define MONITOR_INTERVAL 100000000 /* 100 million nt */
#define OFFSETS_BUFFER_SIZE 1000000
#define WRITE_CHUNK 1000000
#define POSITIONS8_HIGH_SHIFT 32
#define POSITIONS8_LOW_MASK 0xFFFFFFFF

/* #define ALLOW_ODD_PACKSIZES 1 */

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Writing of positions */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif




/************************************************************************
 *   Write procedures -- called by gmapindex/pmapindex
 ************************************************************************/

static Oligospace_T
power (int base, int exponent) {
  Oligospace_T result = 1;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


#if 0
static void
check_bitpack (char *offsetsmetafile, char *offsetsstrmfile,
	       Oligospace_T oligospace, Blocksize_T blocksize) {
  UINT4 *offsetsmeta, *offsetsstrm, *info;
  int offsetsmeta_fd, offsetsstrm_fd;
  size_t offsetsmeta_len, offsetsstrm_len;
  UINT4 *blockptr;
  Oligospace_T oligo;
  int packsize;
#ifndef HAVE_MMAP
  double seconds;
#endif

#ifdef HAVE_MMAP
  offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,offsetsmetafile,/*randomp*/false);
  offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,offsetsstrmfile,/*randomp*/false);
#else
  offsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_len,&seconds,offsetsmetafile,sizeof(UINT4));
  offsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_len,&seconds,offsetsstrmfile,sizeof(UINT4));
#endif

  for (oligo = 0; oligo < oligospace; oligo += blocksize) {
    info = &(offsetsmeta[oligo/BLOCKSIZE * METAINFO_SIZE]);
    blockptr = &(offsetsstrm[info[0]]);
    packsize = (info[2] - info[0])/2;
    if (packsize > 32) {
      fprintf(stderr,"Error:\n");
      fprintf(stderr,"oligo: %08X\n",oligo);
      fprintf(stderr,"nwritten: %u\n",info[0]);
      fprintf(stderr,"packsize: %d\n",packsize);
      fprintf(stderr,"offset: %u\n",info[1]);
      abort();
    }
  }

#ifdef HAVE_MMAP
  munmap((void *) offsetsstrm,offsetsstrm_len);
  munmap((void *) offsetsmeta,offsetsmeta_len);
#else
  FREE(offsetsstrm);
  FREE(offsetsmeta);
#endif

  return;
}
#endif


#if 0
/* Old method */

#ifndef PMAP
static void
check_offsets_from_bitpack (char *offsetsmetafile, char *offsetsstrmfile, Positionsptr_T *offsets,
			    Oligospace_T oligospace, Blocksize_T blocksize) {
  UINT4 *offsetsmeta;
  UINT4 *offsetsstrm;
  Positionsptr_T offsets_decoded[MAX_BITPACK_BLOCKSIZE+1];
  int offsetsmeta_fd, offsetsstrm_fd;
  size_t offsetsmeta_len, offsetsstrm_len;
  Oligospace_T oligoi, i;
#ifndef HAVE_MMAP
  int shmid;
  key_t key;
  double seconds;
#endif


#ifdef HAVE_MMAP
  offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,offsetsmetafile,/*randomp*/false);
  offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,offsetsstrmfile,/*randomp*/false);
#else
  offsetsmeta = (UINT4 *) Access_allocate_private(&shmid,&key,&offsetsmeta_len,&seconds,offsetsmetafile,sizeof(UINT4));
  offsetsstrm = (UINT4 *) Access_allocate_private(&shmid,&key,&offsetsstrm_len,&seconds,offsetsstrmfile,sizeof(UINT4));
#endif

  for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
    Bitpack64_block_offsets(offsets_decoded,oligoi,offsetsmeta,offsetsstrm);
    for (i = 0; i <= blocksize; i++) {
      if (offsets_decoded[i] != offsets[oligoi+i]) {
	fprintf(stderr,"Problem with bitpack at oligo %llu+%llu = %llu: %llu != %llu.  Please inform twu@gene.com\n",
		(unsigned long long) oligoi,(unsigned long long) i,
		(unsigned long long) oligoi+i,(unsigned long long) offsets_decoded[i],(unsigned long long) offsets[oligoi+i]);
	exit(9);
      }
    }
  }



#ifdef HAVE_MMAP
  munmap((void *) offsetsstrm,offsetsstrm_len);
  munmap((void *) offsetsmeta,offsetsmeta_len);
#else
  FREE(offsetsstrm);
  FREE(offsetsmeta);
#endif

  return;
}
#endif
#endif


#if 0
/* Old method */

#ifndef PMAP
static void
check_offsets_from_bitpack_huge (char *offsetspagesfile, char *offsetsmetafile, char *offsetsstrmfile,
				 Hugepositionsptr_T *offsets, Oligospace_T oligospace, Blocksize_T blocksize) {
  UINT4 *offsetspages;
  UINT4 *offsetsmeta;
  UINT4 *offsetsstrm;
  Hugepositionsptr_T offsets64[65];
  int offsetsmeta_fd, offsetsstrm_fd;
  size_t offsetspages_len, offsetsmeta_len, offsetsstrm_len;
  Oligospace_T oligoi, i;
  double seconds;


  offsetspages = (UINT4 *) Access_allocate_private(&offsetspages_len,&seconds,offsetspagesfile,sizeof(UINT4));
#ifdef HAVE_MMAP
  offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,offsetsmetafile,/*randomp*/false);
  offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,offsetsstrmfile,/*randomp*/false);
#else
  offsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_len,&seconds,offsetsmetafile,sizeof(UINT4));
  offsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_len,&seconds,offsetsstrmfile,sizeof(UINT4));
#endif

  for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
    Bitpack64_block_offsets_huge(offsets64,oligoi,offsetspages,offsetsmeta,offsetsstrm);
    for (i = 0; i <= 64; i++) {
      if (offsets64[i] != offsets[oligoi+i]) {
#ifdef OLIGOSPACE_NOT_LONG
	fprintf(stderr,"\nProblem with bitpack64 at oligo %u+%u = %u: uncompressed %u != expected %u.  Your compiler may be defective.  Please inform twu@gene.com\n",
		oligoi,i,oligoi+i,offsets64[i],offsets[oligoi+i]);
#else
	fprintf(stderr,"\nProblem with bitpack64 at oligo %llu+%llu = %llu: uncompressed %llu != expected %llu.  Your compiler may be defective.  Please inform twu@gene.com\n",
		(unsigned long long) oligoi,(unsigned long long) i,(unsigned long long) oligoi+i,
		(unsigned long long) offsets64[i],(unsigned long long) offsets[oligoi+i]);
#endif
      }
    }
  }

#ifdef HAVE_MMAP
  munmap((void *) offsetsstrm,offsetsstrm_len);
  munmap((void *) offsetsmeta,offsetsmeta_len);
#else
  FREE(offsetsstrm);
  FREE(offsetsmeta);
#endif
  FREE(offsetspages);

  return;
}
#endif
#endif





#ifdef HAVE_64_BIT
UINT8
Indexdb_count_offsets (FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
		       Width_T index1part_aa, bool watsonp,
#else
		       Width_T index1part,
#endif
		       Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
  UINT8 noffsets = 0;
  int c;
  char *comma;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Shortoligomer_T low = 0U;
#else
  int between_counter = 0, in_counter = 0;
#endif
  int circular_typeint;
  int nchrs, chrnum;


  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#endif

  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Counting positions in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Counting positions in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

/* Don't need 64-bit distinction, because we need oligomers only for PMAP */
#ifdef PMAP
    /* Need to compute low, so we can determine if we have a stop codon */
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      low = 0U; 
      break;
    default: 
      if (genome_lc_p == true) {
	low = 0U;
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
#else
    switch (uppercaseCode[c]) {
    case 'A': break;
    case 'C': break;
    case 'G': break;
    case 'T': break;
    case 'X': case 'N': in_counter = 0; break;
    default: 
      if (genome_lc_p == true) {
	in_counter = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %llu\n",c,(unsigned long long) position);
	abort();
      }
    }
#endif

#ifdef PMAP
    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));

    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	noffsets += 1;
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }

#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	noffsets += 1;
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
#ifndef PMAP
      in_counter = 0;
#else
      /* high = */ low = /* carry = */ 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#endif

      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

  return noffsets;
}
#endif


#define BLOCKSIZE 64

void
Indexdb_write_offsets (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
		       Alphabet_T alphabet, Width_T index1part_aa, bool watsonp,
#else
		       Width_T index1part,
#endif
		       Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  char *pointersfile, *offsetsfile;
  char *filesuffix;

  /* If offsets[oligospace] > 2^32, then will will want to allocate and write 8-mers for offsets file */
  char *packsizes;
  UINT4 **bitpacks;

  char *comma;
  int c, nchrs, chrnum;
  Oligospace_T oligospace, bmerspace;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
#ifdef HAVE_64_BIT
  Oligospace_T oligo = 0ULL;
#else
  Shortoligomer_T high = 0U, low = 0U, carry;
#endif
  Oligospace_T bmer;

#ifdef PMAP
  int alphabet_size;
  int frame = -1, between_counter[3], in_counter[3];
  Shortoligomer_T aaindex;
  int index1part_nt = 3*index1part_aa;
#else
  Oligospace_T masked, mask;
  int between_counter, in_counter;
#endif
#ifdef DEBUG1
  char *aa;
#endif

  int circular_typeint;

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  oligospace = power(alphabet_size,index1part_aa);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;

#else
#ifdef HAVE_64_BIT
  mask = ~(~0ULL << 2*index1part);
#else
  mask = ~(~0U << 2*index1part);
#endif
  oligospace = power(4,index1part);
  between_counter = in_counter = 0;
#endif

  bmerspace = oligospace/BLOCKSIZE;
  fprintf(stderr,"Allocating %llu*%d bytes for packsizes\n",bmerspace,(int) sizeof(char));
  packsizes = (char *) CALLOC_NO_EXCEPTION(bmerspace,sizeof(char));
  fprintf(stderr,"Allocating %llu*%d bytes for bitpacks\n",bmerspace,(int) sizeof(UINT4 *));
  bitpacks = (UINT4 **) CALLOC_NO_EXCEPTION(bmerspace,sizeof(UINT4 *));

  if (packsizes == NULL || bitpacks == NULL) {
#ifdef PMAP
    fprintf(stderr,"Unable to allocate %llu+%llu bytes of memory, needed to build offsets with %d-mers\n",
	    bmerspace*sizeof(char),bmerspace*sizeof(UINT4 *),index1part_aa);
#else
    fprintf(stderr,"Unable to allocate %llu+%llu bytes of memory, needed to build offsets with %d-mers\n",
	    bmerspace*sizeof(char),bmerspace*sizeof(UINT4 *),index1part);
#endif
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }


  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef HAVE_64_BIT
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N':
      oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %llu\n",c,(unsigned long long) position);
	abort();
      }
    }

#else
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      high = low = carry = 0U; 
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
#ifdef PMAP
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
	in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#endif

#ifdef PMAP
    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));

    if (in_counter[frame] > 0) {
      if (watsonp == true) {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_fwd(oligo) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
#endif
      } else {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_rev(oligo) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
#endif
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
#ifdef HAVE_64_BIT
	aaindex = Alphabet_get_aa_index(oligo,watsonp,index1part_nt);
#else
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
#endif
	bmer = aaindex/BLOCKSIZE;
	if (Bitpack64_access_filledp(aaindex,(int) packsizes[bmer],bitpacks[bmer]) == true) {
	  bitpacks[bmer] = Bitpack64_realloc_one((int) packsizes[bmer],bitpacks[bmer]);
	  packsizes[bmer] += 2;
	}
	Bitpack64_incr_bitpack(aaindex,(int) packsizes[bmer],bitpacks[bmer]);

	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }

#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	bmer = masked/BLOCKSIZE;
	if (Bitpack64_access_filledp(masked,(int) packsizes[bmer],bitpacks[bmer]) == true) {
	  bitpacks[bmer] = Bitpack64_realloc_one((int) packsizes[bmer],bitpacks[bmer]);
	  packsizes[bmer] += 2;
	}
	Bitpack64_incr_bitpack(masked,(int) packsizes[bmer],bitpacks[bmer]);

	between_counter = 0;
      }
      in_counter--;
    }
#endif

    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
#ifndef PMAP
      oligo = 0;
      in_counter = 0;
#elif defined(HAVE_64_BIT)
      oligo = 0ULL;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#endif

      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }


#ifdef PMAP
  if (watsonp == true) {
    filesuffix = FWD_FILESUFFIX;
  } else {
    filesuffix = REV_FILESUFFIX;
  }
#else
  filesuffix = IDX_FILESUFFIX;
#endif

  pointersfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				 strlen(".")+strlen(filesuffix)+
				 /*for kmer*/2+/*for interval char*/1+
				 strlen("offsets64meta")+1,sizeof(char));
#ifdef PMAP
  sprintf(pointersfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"offsets64meta");
#else
  sprintf(pointersfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,index1part,interval_char,"offsets64meta");
#endif
	
  offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				strlen(".")+strlen(filesuffix)+
				/*for kmer*/2+/*for interval char*/1+
				strlen("offsets64strm")+1,sizeof(char));
#ifdef PMAP
  sprintf(offsetsfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"offsets64strm");
#else
  sprintf(offsetsfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,index1part,interval_char,"offsets64strm");
#endif

  fprintf(stderr,"Writing %llu offsets compressed via bitpack64...",oligospace+1UL);
  Bitpack64_write_differential_bitpacks(pointersfile,offsetsfile,packsizes,bitpacks,oligospace);
  fprintf(stderr,"done\n");

  FREE(offsetsfile);
  FREE(pointersfile);

  for (bmer = 0; bmer < bmerspace; bmer++) {
    FREE(bitpacks[bmer]);
  }
  FREE(bitpacks);
  FREE(packsizes);

  return;
}



#ifdef HAVE_64_BIT

void
Indexdb_write_offsets_huge (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			    Alphabet_T alphabet, Width_T index1part_aa, bool watsonp,
#else
			    Width_T index1part,
#endif
			    Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  char *pagesfile, *pointersfile, *offsetsfile;
  char *filesuffix;

  /* If offsets[oligospace] > 2^32, then we allocate and write 8-mers for offsets file */
  char *packsizes;
  UINT4 **bitpacks;

  char *comma;
  int c, nchrs, chrnum;
  Oligospace_T oligospace, bmerspace;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
#ifdef HAVE_64_BIT
  Oligospace_T oligo = 0ULL;
#else
  Shortoligomer_T high = 0U, low = 0U, carry;
#endif
  Oligospace_T bmer;

#ifdef PMAP
  int alphabet_size;
  int frame = -1, between_counter[3], in_counter[3];
  Oligospace_T aaindex;
  int index1part_nt = 3*index1part_aa;
#else
  Oligospace_T masked, mask;
  int between_counter, in_counter;
#endif
#ifdef DEBUG1
  char *aa;
#endif

  int circular_typeint;

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  oligospace = power(alphabet_size,index1part_aa);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;

#else
#ifdef HAVE_64_BIT
  mask = ~(~0ULL << 2*index1part);
#else
  mask = ~(~0U << 2*index1part);
#endif
  oligospace = power(4,index1part);
  between_counter = in_counter = 0;
#endif

  bmerspace = oligospace/BLOCKSIZE;

  fprintf(stderr,"Allocating %llu*%d bytes for packsizes\n",bmerspace,(int) sizeof(char));
  packsizes = (char *) CALLOC_NO_EXCEPTION(bmerspace,sizeof(char));
  fprintf(stderr,"Allocating %llu*%d bytes for bitpacks\n",bmerspace,(int) sizeof(UINT4 *));
  bitpacks = (UINT4 **) CALLOC_NO_EXCEPTION(bmerspace,sizeof(UINT4 *));

  if (packsizes == NULL || bitpacks == NULL) {
#ifdef PMAP
    fprintf(stderr,"Unable to allocate %llu+%llu bytes of memory, needed to build offsets with %d-mers\n",
	    bmerspace*sizeof(char),bmerspace*sizeof(UINT4 *),index1part_aa);
#else
    fprintf(stderr,"Unable to allocate %llu+%llu bytes of memory, needed to build offsets with %d-mers\n",
	    bmerspace*sizeof(char),bmerspace*sizeof(UINT4 *),index1part);
#endif
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }


  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {

#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef HAVE_64_BIT
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N':
      oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %llu\n",c,(unsigned long long) position);
	abort();
      }
    }

#else
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      high = low = carry = 0U; 
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
#ifdef PMAP
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
	in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#endif

#ifdef PMAP
    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));

    if (in_counter[frame] > 0) {
      if (watsonp == true) {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_fwd(oligo) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
#endif
      } else {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_rev(oligo) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
#endif
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
#ifdef HAVE_64_BIT
	aaindex = Alphabet_get_aa_index(oligo,watsonp,index1part_nt);
#else
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
#endif
	bmer = aaindex/BLOCKSIZE;
	if (Bitpack64_access_filledp(aaindex,(int) packsizes[bmer],bitpacks[bmer]) == true) {
	  bitpacks[bmer] = Bitpack64_realloc_one((int) packsizes[bmer],bitpacks[bmer]);
	  packsizes[bmer] += 2;
	}
	Bitpack64_incr_bitpack(aaindex,(int) packsizes[bmer],bitpacks[bmer]);

	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }

#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	bmer = masked/BLOCKSIZE;
	if (Bitpack64_access_filledp(masked,(int) packsizes[bmer],bitpacks[bmer]) == true) {
	  bitpacks[bmer] = Bitpack64_realloc_one((int) packsizes[bmer],bitpacks[bmer]);
	  packsizes[bmer] += 2;
	}
	Bitpack64_incr_bitpack(masked,(int) packsizes[bmer],bitpacks[bmer]);

	between_counter = 0;
      }
      in_counter--;
    }
#endif

    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
#ifndef PMAP
      oligo = 0;
      in_counter = 0;
#elif defined(HAVE_64_BIT)
      oligo = 0ULL;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#endif

      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }


#ifdef PMAP
  if (watsonp == true) {
    filesuffix = FWD_FILESUFFIX;
  } else {
    filesuffix = REV_FILESUFFIX;
  }
#else
  filesuffix = IDX_FILESUFFIX;
#endif

  pagesfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			      strlen(".")+strlen(filesuffix)+
			      /*for kmer*/2+/*for interval char*/1+
			      strlen("offsets64pages")+1,sizeof(char));
#ifdef PMAP
  sprintf(pagesfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"offsets64pages");
#else
  sprintf(pagesfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,index1part,interval_char,"offsets64pages");
#endif

  pointersfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				 strlen(".")+strlen(filesuffix)+
				 /*for kmer*/2+/*for interval char*/1+
				 strlen("offsets64meta")+1,sizeof(char));
#ifdef PMAP
  sprintf(pointersfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"offsets64meta");
#else
  sprintf(pointersfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,index1part,interval_char,"offsets64meta");
#endif
	
  offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				strlen(".")+strlen(filesuffix)+
				/*for kmer*/2+/*for interval char*/1+
				strlen("offsets64strm")+1,sizeof(char));
#ifdef PMAP
  sprintf(offsetsfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"offsets64strm");
#else
  sprintf(offsetsfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,index1part,interval_char,"offsets64strm");
#endif

  fprintf(stderr,"Writing %llu offsets compressed via bitpack64...",oligospace+1UL);
  Bitpack64_write_differential_huge_bitpacks(pagesfile,pointersfile,offsetsfile,packsizes,bitpacks,oligospace);

  FREE(offsetsfile);
  FREE(pointersfile);
  FREE(pagesfile);
  
  for (bmer = 0; bmer < bmerspace; bmer++) {
    FREE(bitpacks[bmer]);
  }
  FREE(bitpacks);
  FREE(packsizes);

  return;
}
#endif


UINT4 *
Indexdb_bitpack_counter (UINT4 **counterstrm, UINT4 *offsetsmeta, UINT4 *offsetsstrm,
#ifdef PMAP
			 int alphabet_size, Width_T index1part_aa
#else
			 Width_T index1part
#endif
			 ) {
  UINT4 *countermeta;
  UINT4 nvectors;
  UINT4 *offsets = NULL, diff, maxdiff, k;
  Oligospace_T oligospace, oligoi;
  Blocksize_T blocksize, packsize, firstbit, i;
#ifdef HAVE_BUILTIN_CLZ
#elif defined(HAVE_ASM_BSR)
  int msb;
#endif


#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
#else
  oligospace = power(4,index1part);
#endif
  blocksize = 64;
  offsets = (UINT4 *) MALLOC((blocksize + 1) * sizeof(UINT4));


  countermeta = (UINT4 *) MALLOC((oligospace/BLOCKSIZE + 1) * sizeof(UINT4));

  fprintf(stderr,"Expanding offsetsstrm into counters...");
#if defined(LARGE_GENOMES)
  fprintf(stderr,"Cannot do expand offsets on large genomes\n");
  exit(9);
#else
  k = 0;
  nvectors = 0;
  for (oligoi = 0ULL; oligoi < oligospace; oligoi += blocksize) {
    countermeta[k++] = nvectors;

    Bitpack64_block_offsets(offsets,oligoi,offsetsmeta,offsetsstrm);

    maxdiff = 0;
    for (i = 1; i <= blocksize; i++) {
      if ((diff = offsets[i] - offsets[i-1]) > maxdiff) {
	maxdiff = diff;
      }
    }

    if (maxdiff == 0) {
      packsize = 0;
    } else {
#ifdef HAVE_BUILTIN_CLZ
      firstbit = __builtin_clz(maxdiff);
      packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
      asm("bsr %1,%0" : "=r"(msb) : "r"(maxdiff));
      packsize = msb + 1;
#else
      firstbit = ((maxdiff >> 16) ? clz_table[maxdiff >> 16] : 16 + clz_table[maxdiff]);
      packsize = 32 - firstbit;
#endif
    }

    packsize = (packsize + 1) & ~1; /* Converts packsizes to the next multiple of 2 */
    /* nvectors += (packsize * blocksize)/128; */
    nvectors += packsize / 2;	/* For blocksize of 64 => *64/128 */
  }
#endif
  countermeta[k] = nvectors;
  fprintf(stderr,"done\n");

  fprintf(stderr,"Allocating %u bytes for counterstrm\n",nvectors * 16);
  *counterstrm = (UINT4 *) CALLOC(nvectors * 4,sizeof(UINT4)); /* 4 words per 128-bit vector */

  FREE(offsets);
  return countermeta;
}


#ifdef HAVE_64_BIT
UINT4 *
Indexdb_bitpack_counter_huge (UINT4 **counterstrm, UINT4 *offsetspages, UINT4 *offsetsmeta, UINT4 *offsetsstrm,
#ifdef PMAP
			      int alphabet_size, Width_T index1part_aa
#else
			      Width_T index1part
#endif
			      ) {
  UINT4 *countermeta;
  UINT4 nvectors;
  UINT8 *offsets = NULL;
  UINT4 diff, maxdiff, k;
  Oligospace_T oligospace, oligoi;
  /* Oligospace_T bmerspace; */
  /* UINT4 *pageptr; */
  /* int npages; */
  Blocksize_T packsize, firstbit, i;
#ifdef HAVE_BUILTIN_CLZ
#elif defined(HAVE_ASM_BSR)
  int msb;
#endif


#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
#else
  oligospace = power(4,index1part);
#endif
  /* bmerspace = oligospace/BLOCKSIZE; */
  offsets = (UINT8 *) MALLOC((BLOCKSIZE + 1) * sizeof(UINT8));


  countermeta = (UINT4 *) MALLOC((oligospace/BLOCKSIZE + 1) * sizeof(UINT4));

  fprintf(stderr,"Expanding offsetsstrm into counters...");
  k = 0;
  nvectors = 0;
  for (oligoi = 0ULL; oligoi < oligospace; oligoi += BLOCKSIZE) {
    countermeta[k++] = nvectors;

    Bitpack64_block_offsets_huge(offsets,oligoi,offsetspages,offsetsmeta,offsetsstrm);

    maxdiff = 0;
    for (i = 1; i <= BLOCKSIZE; i++) {
      if ((diff = (UINT4) (offsets[i] - offsets[i-1])) > maxdiff) {
	maxdiff = diff;
      }
    }

    if (maxdiff == 0) {
      packsize = 0;
    } else {
#ifdef HAVE_BUILTIN_CLZ
      firstbit = __builtin_clz(maxdiff);
      packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
      asm("bsr %1,%0" : "=r"(msb) : "r"(maxdiff));
      packsize = msb + 1;
#else
      firstbit = ((maxdiff >> 16) ? clz_table[maxdiff >> 16] : 16 + clz_table[maxdiff]);
      packsize = 32 - firstbit;
#endif
    }

    packsize = (packsize + 1) & ~1; /* Converts packsizes to the next multiple of 2 */
    /* nvectors += (packsize * blocksize)/128; */
    nvectors += packsize / 2;	/* For blocksize of 64 => *64/128 */
  }

  countermeta[k] = nvectors;
  fprintf(stderr,"done\n");

  fprintf(stderr,"Allocating %u bytes for counterstrm\n",nvectors * 16);
  *counterstrm = (UINT4 *) CALLOC(nvectors * 4,sizeof(UINT4)); /* 4 words per 128-bit vector */

  FREE(offsets);
  return countermeta;
}
#endif



#if 0
static bool
need_to_sort_p (Univcoord_T *positions, int length) {
  Univcoord_T prevpos;
  int j;

  prevpos = positions[0];
  for (j = 1; j < length; j++) {
    if (positions[j] <= prevpos) {
      return true;
    } else {
      prevpos = positions[j];
    }
  }
  return false;
}
#endif


#if 0
static void
positions_move_absolute_1 (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*sizeof(unsigned char);

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %zu*%zu=%zu\n",
	    (size_t) ptr,sizeof(unsigned char),(size_t) offset);
    perror("Error in indexdb.c, positions_move_absolute_1");
    exit(9);
  }
  return;
}
#endif

#if 0
static void
positions_move_absolute_4 (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*sizeof(UINT4);

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %zu*%zu=%zu\n",
	    (size_t) ptr,sizeof(UINT4),(size_t) offset);
    perror("Error in indexdb.c, positions_move_absolute_4");
    exit(9);
  }
  return;
}
#endif

#if 0
/* Replaced by positions_move_absolute_1 and positions_move_absolute_4 */
static void
positions_move_absolute_8 (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*((off_t) sizeof(UINT8));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %zu*%zu=%zu\n",
	    (size_t) ptr,sizeof(UINT8),(size_t) offset);
    perror("Error in indexdb.c, positions_move_absolute_8");
    exit(9);
  }
  return;
}
#endif


#if 0
/* Old method */

/* Works directly in file, so we don't need to allocate memory */
static void
compute_positions_in_file (int positions_high_fd, int positions_low_fd, Positionsptr_T *offsets,
			   FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			   Width_T index1part_aa, bool watsonp,
#else
			   Width_T index1part,
#endif
			   Width_T index1interval, bool genome_lc_p, char *fileroot,
			   bool mask_lowercase_p, bool coord_values_8p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  UINT8 adjposition8;
  UINT4 adjposition4;
  Chrpos_T chrpos = 0U;
  char *comma;
  int c, nchrs, chrnum;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Storedoligomer_T aaindex;
  Width_T index1part_nt = 3*index1part_aa;
#else
  Oligospace_T oligo = 0U, masked, mask;
  int between_counter = 0, in_counter = 0;
#endif
  int circular_typeint;

#ifdef ADDSENTINEL
  Oligospace_T oligospace, oligoi;
  oligospace = power(4,index1part);
#endif

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  /* oligospace = power(alphabet_size,index1part_aa); */
#else
#ifdef HAVE_64_BIT
  mask = ~(~0ULL << 2*index1part);
#else
  mask = ~(~0U << 2*index1part);
#endif
  /* oligospace = power(4,index1part); */
#endif

  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1; break;
    case 'G': low = (low << 2) | 2; break;
    case 'T': low = (low << 2) | 3; break;
    case 'X': case 'N': 
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    case 'X': case 'N': oligo = 0U; in_counter = 0; break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U; in_counter = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %llu\n",c,(unsigned long long) position);
	abort();
      }
    }
#endif

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */
    
#ifdef PMAP
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
	if (coord_values_8p == true) {
	  positions_move_absolute_1(positions_high_fd,offsets[aaindex]);
	  positions_move_absolute_4(positions_low_fd,offsets[aaindex]);
	  offsets[aaindex] += 1;
	  if (watsonp == true) {
	    adjposition8 =  position - index1part_nt + 1;
	  } else {
	    adjposition8 = position;
	  }
	  WRITE_CHAR((unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT),positions_high_fd);
	  WRITE_UINT((UINT4) (adjposition8 & POSITIONS8_LOW_MASK),positions_low_fd);
	} else {
	  positions_move_absolute_4(positions_low_fd,offsets[aaindex]);
	  offsets[aaindex] += 1;
	  if (watsonp == true) {
	    adjposition4 = (UINT4) (position - index1part_nt + 1);
	  } else {
	    adjposition4 = (UINT4) position;
	  }
	  WRITE_UINT(adjposition4,positions_low_fd);
	}
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	if (coord_values_8p == true) {
	  positions_move_absolute_1(positions_high_fd,offsets[masked]);
	  positions_move_absolute_4(positions_low_fd,offsets[masked]);
	  offsets[masked] += 1;
	  adjposition8 = position - index1part + 1;
	  WRITE_CHAR((unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT),positions_high_fd);
	  WRITE_UINT((UINT4) (adjposition8 & POSITIONS8_LOW_MASK),positions_low_fd);
	} else {
	  positions_move_absolute_4(positions_low_fd,offsets[masked]);
	  offsets[masked] += 1;
	  adjposition4 = (UINT4) (position - index1part + 1);
	  WRITE_UINT(adjposition4,positions_low_fd);
	}
	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
#ifdef PMAP
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      oligo = 0U; in_counter = 0;
#endif
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

#ifdef ADDSENTINEL
  if (coord_values_8p == true) {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      positions_move_absolute_1(positions_high_fd,offsets[oligoi]);
      positions_move_absolute_4(positions_low_fd,offsets[oligoi]);
      WRITE_CHAR((unsigned char) (-1UL >> POSITIONS8_HIGH_SHIFT),positions_high_fd);
      WRITE_UINT((UINT4) (-1UL & POSITIONS8_LOW_SHIFT),positions_low_fd);
    }
  } else {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      positions_move_absolute_4(positions_low_fd,offsets[oligoi]);
      WRITE_UINT(-1U,positions_low_fd);
    }
  }
#endif

  return;
}
#endif


#if 0
/* Old method */

/* Requires sufficient memory to hold all positions */
static void
compute_positions_in_memory (UINT4 *positions4, unsigned char *positions8_high, UINT4 *positions8_low,
			     Positionsptr_T *offsets, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			     Width_T index1part_aa, bool watsonp,
#else
			     Width_T index1part,
#endif
			     Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
			     bool coord_values_8p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
  char *comma;
  int c, nchrs, chrnum;
  UINT8 adjposition8;
  Positionsptr_T offset;

#ifdef HAVE_64_BIT
  Oligospace_T oligo = 0ULL;
#else
  Shortoligomer_T high = 0U, low = 0U, carry;
#endif

#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Oligospace_T aaindex;
  Width_T index1part_nt = 3*index1part_aa;
  debug1(char *aa);
#else
  Oligospace_T masked, mask;
  int between_counter, in_counter;
  debug1(char *nt);
#endif

  int circular_typeint;

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;

#else
#ifdef HAVE_64_BIT
  mask = ~(~0ULL << 2*index1part);
#else
  mask = ~(~0U << 2*index1part);
#endif
  between_counter = in_counter = 0;
#endif


  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef HAVE_64_BIT
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N':
      oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %llu\n",c,(unsigned long long) position);
	abort();
      }
    }

#else
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      high = low = carry = 0U; 
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
#ifdef PMAP
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
	in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#endif

#ifdef PMAP
    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));
    
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_fwd(oligo) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#endif
      } else {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_rev(oligo) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#endif
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
#ifdef HAVE_64_BIT
	aaindex = Alphabet_get_aa_index(oligo,watsonp,index1part_nt);
#else
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
#endif
	offset = Bitpack64_read_one(aaindex,offsetsmeta,offsetsstrm) + Bitpack64_incr(aaindex,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  if (watsonp == true) {
	    adjposition8 = position - index1part_nt + 1;
	    positions8_high[offset] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offset] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position-index1part_nt+1UL);
	  } else {
	    positions8_high[offset] = (unsigned char) (position >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offset] = (UINT4) (position & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position);
	  }
	} else {
	  if (watsonp == true) {
	    positions4[offset] = (UINT4) (position - index1part_nt + 1);
	    debug1(adjposition = position - index1part_nt + 1);
	  } else {
	    positions4[offset] = (UINT4) position;
	    debug1(adjposition = position);
	  }
	}

	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }

#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	offset = Bitpack64_read_one(masked,offsetsmeta,offsetsstrm) + Bitpack64_incr(masked,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  adjposition8 = position - index1part + 1;
	  positions8_high[offset] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	  positions8_low[offset] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	} else {
	  positions4[offset] = (UINT4) (position - index1part + 1);
	}
	debug1(nt = shortoligo_nt(masked,index1part);
	       printf("Storing %s at %llu, chrpos %u\n",nt,(unsigned long long) (position-index1part+1U),chrpos-index1part+1U);
	       FREE(nt));

	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
      debug1(printf("Skipping because position %u is at chrbound\n",position));
#ifndef PMAP
      oligo = 0;
      in_counter = 0;
#elif defined(HAVE_64_BIT)
      oligo = 0ULL;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#endif

      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

  return;
}
#endif


/* Requires sufficient memory to hold all positions */
static void
compute_positions (UINT4 *positions4, unsigned char *positions8_high, UINT4 *positions8_low,
		   UINT4 *offsetsmeta, UINT4 *offsetsstrm, UINT4 *countermeta, UINT4 *counterstrm,
		   FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
		   Width_T index1part_aa, bool watsonp,
#else
		   Width_T index1part,
#endif
		   Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
		   bool coord_values_8p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
  char *comma;
  int c, nchrs, chrnum;
  UINT8 adjposition8;
  Positionsptr_T offset;

#ifdef HAVE_64_BIT
  Oligospace_T oligo = 0ULL;
#else
  Shortoligomer_T high = 0U, low = 0U, carry;
#endif

#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Oligospace_T aaindex;
  Width_T index1part_nt = 3*index1part_aa;
  debug1(char *aa);
#else
  Oligospace_T masked, mask;
  int between_counter = 0, in_counter = 0;
  debug1(char *nt);
#endif

  int circular_typeint;

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
#ifdef HAVE_64_BIT
  mask = ~(~0ULL << 2*index1part);
#else
  mask = ~(~0U << 2*index1part);
#endif
#endif


  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef HAVE_64_BIT
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N':
      oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %llu\n",c,(unsigned long long) position);
	abort();
      }
    }

#else
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      high = low = carry = 0U; 
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
#ifdef PMAP
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
	in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#endif

#ifdef PMAP
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_fwd(oligo) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#endif
      } else {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_rev(oligo) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#endif
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
#ifdef HAVE_64_BIT
	aaindex = Alphabet_get_aa_index(oligo,watsonp,index1part_nt);
#else
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
#endif
	offset = Bitpack64_read_one(aaindex,offsetsmeta,offsetsstrm) + Bitpack64_incr(aaindex,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  if (watsonp == true) {
	    adjposition8 = position - index1part_nt + 1;
	    positions8_high[offset] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offset] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position-index1part_nt+1UL);
	  } else {
	    positions8_high[offset] = (unsigned char) (position >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offset] = (UINT4) (position & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position);
	  }
	} else {
	  if (watsonp == true) {
	    positions4[offset] = (UINT4) (position - index1part_nt + 1);
	    debug1(adjposition = position - index1part_nt + 1);
	  } else {
	    positions4[offset] = (UINT4) position;
	    debug1(adjposition = position);
	  }
	}

	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }

#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	offset = Bitpack64_read_one(masked,offsetsmeta,offsetsstrm) + Bitpack64_incr(masked,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  adjposition8 = position - index1part + 1;
	  positions8_high[offset] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	  positions8_low[offset] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	} else {
	  positions4[offset] = (UINT4) (position - index1part + 1);
	}

	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
      debug1(printf("Skipping because position %u is at chrbound\n",position));
#ifndef PMAP
      oligo = 0;
      in_counter = 0;
#elif defined(HAVE_64_BIT)
      oligo = 0ULL;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#endif

      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

  return;
}


#ifdef HAVE_64_BIT
/* Requires sufficient memory to hold all positions */
/* Main difference is that we use Hugepositionsptr_T * for offsets.
   Also, huge genome implies large genome, so need only positions8. */
static void
compute_positions_huge (UINT4 *positions4, unsigned char *positions8_high, UINT4 *positions8_low,
			UINT4 *offsetspages, UINT4 *offsetsmeta, UINT4 *offsetsstrm,
			UINT4 *countermeta, UINT4 *counterstrm,
			FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			Width_T index1part_aa, bool watsonp,
#else
			Width_T index1part,
#endif
			Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
			bool coord_values_8p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
  char *comma;
  int c, nchrs, chrnum;
  UINT8 adjposition8;
  Hugepositionsptr_T offset;

#ifdef HAVE_64_BIT
  Oligospace_T oligo = 0ULL;
#else
  Shortoligomer_T high = 0U, low = 0U, carry;
#endif

#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Oligospace_T aaindex;
  Width_T index1part_nt = 3*index1part_aa;
  debug1(char *aa);
#else
  Oligospace_T masked, mask;
  int between_counter = 0, in_counter = 0;
  debug1(char *nt);
#endif

  int circular_typeint;


  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
#ifdef HAVE_64_BIT
  mask = ~(~0ULL << 2*index1part);
#else
  mask = ~(~0U << 2*index1part);
#endif
#endif


  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef HAVE_64_BIT
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N':
      oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U;
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %llu\n",c,(unsigned long long) position);
	abort();
      }
    }

#else
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      high = low = carry = 0U; 
#ifdef PMAP
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
#ifdef PMAP
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
	in_counter = 0;
#endif
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#endif
    
#ifdef PMAP
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_fwd(oligo) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#endif
      } else {
#ifdef HAVE_64_BIT
	if (Alphabet_get_codon_rev(oligo) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#else
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
#endif
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
#ifdef HAVE_64_BIT
	aaindex = Alphabet_get_aa_index(oligo,watsonp,index1part_nt);
#else
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
#endif
	offset = Bitpack64_read_one_huge(aaindex,offsetspages,offsetsmeta,offsetsstrm) +
	  Bitpack64_incr(aaindex,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  if (watsonp == true) {
	    adjposition8 = position - index1part_nt + 1;
	    positions8_high[offset] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offset] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position-index1part_nt+1);
	  } else {
	    positions8_high[offset] = (unsigned char) (position >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offset] = (UINT4) (position & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position);
	  }
	} else {
	  if (watsonp == true) {
	    positions4[offset] = (UINT4) (position-index1part_nt+1U);
	    debug1(adjposition = position-index1part_nt+1U);
	  } else {
	    positions4[offset] = (UINT4) position;
	    debug1(adjposition = position);
	  }
	}

	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }

#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	offset = Bitpack64_read_one_huge(masked,offsetspages,offsetsmeta,offsetsstrm) +
	  Bitpack64_incr(masked,countermeta,counterstrm);
	if (coord_values_8p == true) {
	  adjposition8 = position - index1part + 1;
	  positions8_high[offset] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	  positions8_low[offset] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	} else {
	  positions4[offset] = (UINT4) (position - index1part + 1);
	}

	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
      debug1(printf("Skipping because position %u is at chrbound\n",position));
#ifndef PMAP
      oligo = 0;
      in_counter = 0;
#elif defined(HAVE_64_BIT)
      oligo = 0ULL;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#endif

      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

  return;
}
#endif


#define WRITE_CHUNK 1000000

void
Indexdb_write_positions (char *positionsfile_high, char *positionsfile_low, char *offsetsmetafile, char *offsetsstrmfile,
			 FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			 Alphabet_T alphabet, int index1part_aa, bool watsonp,
#else
			 int index1part,
#endif
			 int index1interval, Univcoord_T genomelength, bool genome_lc_p, bool writefilep,
			 char *fileroot, bool mask_lowercase_p, bool coord_values_8p) {
  FILE *positions_high_fp, *positions_low_fp; /* For building positions in memory */
  /* int positions_high_fd, positions_low_fd; */ /* For building positions in file */
  int offsetsmeta_fd, offsetsstrm_fd;
  size_t offsetsmeta_len, offsetsstrm_len;
  UINT4 *offsetsmeta, *offsetsstrm, *countermeta, *counterstrm;
  Oligospace_T oligospace;
  Positionsptr_T totalcounts, count;
  UINT4 *positions4;
  unsigned char *positions8_high;
  UINT4 *positions8_low;
  size_t filesize;
#ifdef PMAP
  int alphabet_size;
#endif


#ifdef HAVE_MMAP
  offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,offsetsmetafile,/*randomp*/false);
  offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,offsetsstrmfile,/*randomp*/false);
#else
  offsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_len,&seconds,offsetsmetafile,sizeof(UINT4));
  offsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_len,&seconds,offsetsstrmfile,sizeof(UINT4));
#endif

#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  oligospace = power(alphabet_size,index1part_aa);
#else
  oligospace = power(4,index1part);
#endif

  totalcounts = Bitpack64_read_one(oligospace,offsetsmeta,offsetsstrm);
  if (totalcounts == 0) {
    if (
#ifdef PMAP
	genomelength > (Univcoord_T) (index1part_aa * 3)
#else
	genomelength > (Univcoord_T) index1part
#endif
	) {
      fprintf(stderr,"Something is wrong with the offsets file.  Total counts is zero.\n");
      exit(9);
    } else {
#if 0
      if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	exit(9);
      } else if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
      fclose(positions_high_fp);
      fclose(positions_low_fp);
#endif
#ifdef HAVE_MMAP
      munmap((void *) offsetsstrm,offsetsstrm_len);
      munmap((void *) offsetsmeta,offsetsmeta_len);
#else
      FREE(offsetsstrm);
      FREE(offsetsmeta);
#endif
      return;
    }
  } else {
#ifdef PMAP
    countermeta = Indexdb_bitpack_counter(&counterstrm,offsetsmeta,offsetsstrm,alphabet_size,index1part_aa);
#else
    countermeta = Indexdb_bitpack_counter(&counterstrm,offsetsmeta,offsetsstrm,index1part);
#endif
  }

  if (writefilep == true) {
    fprintf(stderr,"User requested build of positions in file.  Not supported\n");
    abort();
#if 0
    positions_high_fd = Access_fileio_rw(positionsfile_high);
    positions_low_fd = Access_fileio_rw(positionsfile_low);
#ifdef PMAP
    compute_positions_in_file(positions_high_fd,positions_low_fd,offsets,sequence_fp,chromosome_iit,
			      index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,mask_lowercase_p,
			      coord_values_8p);
#else
    compute_positions_in_file(positions_high_fd,positions_low_fd,offsets,sequence_fp,chromosome_iit,
			      index1part,index1interval,genome_lc_p,fileroot,mask_lowercase_p,
			      coord_values_8p);
#endif
    close(positions_high_fd);
    close(positions_low_fd);
#endif

  } else if (coord_values_8p == true) {
    fprintf(stderr,"Trying to allocate %u*(%d+%d) bytes of memory for positions...",totalcounts,(int) sizeof(unsigned char),(int) sizeof(UINT4));
    positions8_high = (unsigned char *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(unsigned char));
    positions8_low = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions8_high == NULL || positions8_low == NULL) {
      fprintf(stderr,"failed.  Building positions in file.  Not supported\n");
      abort();
#if 0
      positions_high_fd = Access_fileio_rw(positionsfile_high);
      positions_low_fd = Access_fileio_rw(positionsfile_low);
#ifdef PMAP
      compute_positions_in_file(positions_high_fd,positions_low_fd,offsets,sequence_fp,chromosome_iit,
				index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,
				mask_lowercase_p,/*coord_values_8p*/true);
#else
      compute_positions_in_file(positions_high_fd,positions_low_fd,offsets,sequence_fp,chromosome_iit,
				index1part,index1interval,genome_lc_p,fileroot,
				mask_lowercase_p,/*coord_values_8p*/true);
#endif
      close(positions_high_fd);
      close(positions_low_fd);

      if ((filesize = Access_filesize(positionsfile_high)) != (size_t) (totalcounts * sizeof(unsigned char))) {
	fprintf(stderr,"Error after file-based build: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
		positionsfile_high,totalcounts*sizeof(unsigned char),filesize);
	abort();
      } else if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	fprintf(stderr,"Error after file-based build: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	abort();
      }
#endif

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	exit(9);
      } else if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
#ifdef PMAP
      compute_positions(/*positions4*/NULL,positions8_high,positions8_low,
			offsetsmeta,offsetsstrm,countermeta,counterstrm,
			sequence_fp,chromosome_iit,index1part_aa,watsonp,
			index1interval,genome_lc_p,fileroot,
			mask_lowercase_p,/*coord_values_8p*/true);
#else
      compute_positions(/*positions4*/NULL,positions8_high,positions8_low,
			offsetsmeta,offsetsstrm,countermeta,counterstrm,
			sequence_fp,chromosome_iit,index1part,
			index1interval,genome_lc_p,fileroot,
			mask_lowercase_p,/*coord_values_8p*/true);
#endif
      fprintf(stderr,"Writing %u genomic positions to files %s and %s ...\n",
	      totalcounts,positionsfile_high,positionsfile_low);
      FWRITE_CHARS(positions8_high,totalcounts,positions_high_fp);
      FWRITE_UINTS(positions8_low,totalcounts,positions_low_fp);

      fclose(positions_high_fp);
      fclose(positions_low_fp);

      if ((filesize = Access_filesize(positionsfile_high)) != (size_t) (totalcounts * sizeof(unsigned char))) {
	fprintf(stderr,"Error: expected file size for %s is %zu, but observed only %zu.  Trying now to write with smaller chunks.\n",
		positionsfile_high,totalcounts*sizeof(unsigned char),filesize);
	if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_CHARS(&(positions8_high[count]),count,positions_high_fp);
	}
	if (count < totalcounts) {
	  FWRITE_CHARS(&(positions8_high[count]),totalcounts - count,positions_high_fp);
	}
	fclose(positions_high_fp);

	if ((filesize = Access_filesize(positionsfile_high)) != (size_t) (totalcounts * sizeof(unsigned char))) {
	  fprintf(stderr,"Error persists: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_high,totalcounts*sizeof(unsigned char),filesize);
	  abort();
	}
      }

      if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	fprintf(stderr,"Error: expected file size for %s is %zu, but observed only %zu.  Trying now to write with smaller chunks.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_UINTS(&(positions8_low[count]),count,positions_low_fp);
	}
	if (count < totalcounts) {
	  FWRITE_UINTS(&(positions8_low[count]),totalcounts - count,positions_low_fp);
	}
	fclose(positions_low_fp);

	if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	  fprintf(stderr,"Error persists: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	  abort();
	}
      }

      FREE(positions8_high);
      FREE(positions8_low);
    }

  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory for positions...",totalcounts,(int) sizeof(UINT4));
    positions4 = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions4 == NULL) {
      fprintf(stderr,"failed.  Building positions in file.  Not supported.\n");
      abort();
#if 0
      positions_low_fd = Access_fileio_rw(positionsfile_low);
#ifdef PMAP
      compute_positions_in_file(/*positions_high_fd*/0,positions_low_fd,offsets,sequence_fp,chromosome_iit,
				index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,mask_lowercase_p,
				/*coord_values_8*/false);
#else
      compute_positions_in_file(/*positions_high_fd*/0,positions_low_fd,offsets,sequence_fp,chromosome_iit,
				index1part,index1interval,genome_lc_p,fileroot,mask_lowercase_p,
				/*coord_values_8p*/false);
#endif
      close(positions_low_fd);

      if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	fprintf(stderr,"Error after file-based build: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	abort();
      }
#endif

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
#ifdef PMAP
      compute_positions(positions4,/*positions8_high*/NULL,/*positions8_low*/NULL,
			offsetsmeta,offsetsstrm,countermeta,counterstrm,
			sequence_fp,chromosome_iit,
			index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,
			mask_lowercase_p,/*coord_values_8p*/false);
#else
      compute_positions(positions4,/*positions8_high*/NULL,/*positions8_low*/NULL,
			offsetsmeta,offsetsstrm,countermeta,counterstrm,
			sequence_fp,chromosome_iit,
			index1part,index1interval,genome_lc_p,fileroot,
			mask_lowercase_p,/*coord_values_8p*/false);
#endif
      fprintf(stderr,"Writing %u genomic positions to file %s ...\n",
	      totalcounts,positionsfile_low);
      FWRITE_UINTS(positions4,totalcounts,positions_low_fp);
      fclose(positions_low_fp);
      fprintf(stderr,"done\n");

      if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	fprintf(stderr,"Error: expected file size for %s is %zu, but observed only %zu.  Trying now to write with smaller chunks.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_UINTS(&(positions4[count]),count,positions_low_fp);
	}
	if (count < totalcounts) {
	  FWRITE_UINTS(&(positions4[count]),totalcounts - count,positions_low_fp);
	}
	fclose(positions_low_fp);

	if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	  fprintf(stderr,"Error persists: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	  abort();
	}
      }

      FREE(positions4);
    }
  }

  FREE(counterstrm);
  FREE(countermeta);

#ifdef HAVE_MMAP
  munmap((void *) offsetsstrm,offsetsstrm_len);
  munmap((void *) offsetsmeta,offsetsmeta_len);
#else
  FREE(offsetsstrm);
  FREE(offsetsmeta);
#endif

  return;
}



#ifdef HAVE_64_BIT
void
Indexdb_write_positions_huge (char *positionsfile_high, char *positionsfile_low,
			      char *offsetspagesfile, char *offsetsmetafile, char *offsetsstrmfile,
			      FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			      Alphabet_T alphabet, int index1part_aa, bool watsonp,
#else
			      int index1part,
#endif
			      int index1interval, Univcoord_T genomelength, bool genome_lc_p, bool writefilep,
			      char *fileroot, bool mask_lowercase_p, bool coord_values_8p) {
  FILE *positions_high_fp, *positions_low_fp; /* For building positions in memory */
  int offsetspages_fd, offsetsmeta_fd, offsetsstrm_fd;
  size_t offsetspages_len, offsetsmeta_len, offsetsstrm_len;
  UINT4 *offsetspages, *offsetsmeta, *offsetsstrm, *countermeta, *counterstrm;
  Oligospace_T oligospace;
  Hugepositionsptr_T totalcounts;
  UINT4 *positions4;
  unsigned char *positions8_high;
  UINT4 *positions8_low;
  size_t filesize;
#ifdef PMAP
  int alphabet_size;
#endif
  UINT4 count;


#ifdef HAVE_MMAP
  offsetspages = (UINT4 *) Access_mmap(&offsetspages_fd,&offsetspages_len,offsetspagesfile,/*randomp*/false);
  offsetsmeta = (UINT4 *) Access_mmap(&offsetsmeta_fd,&offsetsmeta_len,offsetsmetafile,/*randomp*/false);
  offsetsstrm = (UINT4 *) Access_mmap(&offsetsstrm_fd,&offsetsstrm_len,offsetsstrmfile,/*randomp*/false);
#else
  offsetspages = (UINT4 *) Access_allocate_private(&offsetspages_len,&seconds,offsetspagesfile,sizeof(UINT4));
  offsetsmeta = (UINT4 *) Access_allocate_private(&offsetsmeta_len,&seconds,offsetsmetafile,sizeof(UINT4));
  offsetsstrm = (UINT4 *) Access_allocate_private(&offsetsstrm_len,&seconds,offsetsstrmfile,sizeof(UINT4));
#endif
  fprintf(stderr,"offsetspages has len %zu\n",offsetspages_len);
  fprintf(stderr,"offsetsmeta has len %zu\n",offsetsmeta_len);
  fprintf(stderr,"offsetsstrm has len %zu\n",offsetsstrm_len);


#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  oligospace = power(alphabet_size,index1part_aa);
#else
  oligospace = power(4,index1part);
  oligospace = power(4,index1part);
#endif

  totalcounts = Bitpack64_read_one_huge(oligospace,offsetspages,offsetsmeta,offsetsstrm);
  if (totalcounts == 0) {
    if (
#ifdef PMAP
	genomelength > (Univcoord_T) (index1part_aa * 3)
#else
	genomelength > (Univcoord_T) index1part
#endif
	) {
      fprintf(stderr,"Something is wrong with the offsets file.  Total counts is zero.\n");
      exit(9);
    } else {
#if 0
      if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	exit(9);
      } else if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
      fclose(positions_high_fp);
      fclose(positions_low_fp);
#endif
#ifdef HAVE_MMAP
      munmap((void *) offsetsstrm,offsetsstrm_len);
      munmap((void *) offsetsmeta,offsetsmeta_len);
      munmap((void *) offsetspages,offsetspages_len);
#else
      FREE(offsetsstrm);
      FREE(offsetsmeta);
      FREE(offsetspages);
#endif
      return;
    }

  } else {
#ifdef PMAP
    countermeta = Indexdb_bitpack_counter_huge(&counterstrm,offsetspages,offsetsmeta,offsetsstrm,alphabet_size,index1part_aa);
#else
    countermeta = Indexdb_bitpack_counter_huge(&counterstrm,offsetspages,offsetsmeta,offsetsstrm,index1part);
#endif
  }

  if (writefilep == true) {
    fprintf(stderr,"User requested build of positions in file.  Not supported\n");
    abort();

  } else if (coord_values_8p == true) {
    fprintf(stderr,"Trying to allocate %llu*(%d+%d) bytes of memory for positions...",
	    totalcounts,(int) sizeof(unsigned char),(int) sizeof(UINT8));
    positions8_high = (unsigned char *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(unsigned char));
    positions8_low = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions8_high == NULL || positions8_low == NULL) {
      fprintf(stderr,"failed.  Building positions in file.  Not supported\n");
      abort();

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	exit(9);
      } else if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
#ifdef PMAP
      compute_positions_huge(/*positions4*/NULL,positions8_high,positions8_low,
			     offsetspages,offsetsmeta,offsetsstrm,countermeta,counterstrm,
			     sequence_fp,chromosome_iit,index1part_aa,watsonp,
			     index1interval,genome_lc_p,fileroot,
			     mask_lowercase_p,/*coord_values_8p*/true);
#else
      compute_positions_huge(/*positions4*/NULL,positions8_high,positions8_low,
			     offsetspages,offsetsmeta,offsetsstrm,countermeta,counterstrm,
			     sequence_fp,chromosome_iit,index1part,index1interval,
			     genome_lc_p,fileroot,mask_lowercase_p,/*coord_values_8p*/true);
#endif
      fprintf(stderr,"Writing %llu genomic positions to files %s and %s...\n",
	      (unsigned long long) totalcounts,positionsfile_high,positionsfile_low);
      FWRITE_CHARS(positions8_high,totalcounts,positions_high_fp);
      FWRITE_UINTS(positions8_low,totalcounts,positions_low_fp);

      fclose(positions_high_fp);
      fclose(positions_low_fp);

      if ((filesize = Access_filesize(positionsfile_high)) != (size_t) (totalcounts * sizeof(unsigned char))) {
	fprintf(stderr,"Error: expected file size for %s is %zu, but observed only %zu.  Trying now to write with smaller chunks.\n",
	positionsfile_high,(size_t) (totalcounts*sizeof(unsigned char)),filesize);
	if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_CHARS(&(positions8_high[count]),count,positions_high_fp);
	}
	if (count < totalcounts) {
	  FWRITE_CHARS(&(positions8_high[count]),totalcounts - count,positions_high_fp);
	}
	fclose(positions_high_fp);

	if ((filesize = Access_filesize(positionsfile_high)) != (size_t) (totalcounts * sizeof(unsigned char))) {
	  fprintf(stderr,"Error persists: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
	          positionsfile_high,(size_t) (totalcounts*sizeof(unsigned char)),filesize);
	  abort();
	}
      }

      if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	fprintf(stderr,"Error: expected file size for %s is %zu, but observed only %zu.  Trying now to write with smaller chunks.\n",
		positionsfile_low,(size_t) (totalcounts*sizeof(UINT4)),filesize);
	if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_UINTS(&(positions8_low[count]),count,positions_low_fp);
	}
	if (count < totalcounts) {
	  FWRITE_UINTS(&(positions8_low[count]),totalcounts - count,positions_low_fp);
	}
	fclose(positions_low_fp);

	if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	  fprintf(stderr,"Error persists: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_low,(size_t) (totalcounts*sizeof(UINT4)),filesize);
	  abort();
	}
      }

      FREE(positions8_high);
      FREE(positions8_low);
    }

  } else {
    fprintf(stderr,"Something is wrong.  We have a huge genome (offsets need 8 bytes), but positions can fit in 4 bytes\n");
    fprintf(stderr,"Please report this bug to twu@gene.com\n");
    abort();

    fprintf(stderr,"Trying to allocate %llu*%d bytes of memory for positions...",totalcounts,(int) sizeof(UINT4));
    positions4 = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions4 == NULL) {
      fprintf(stderr,"failed.  Not able to proceed.\n");
      abort();

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
#ifdef PMAP
      compute_positions_huge(positions4,/*positions8_high*/NULL,/*positions8_low*/NULL,
			     offsetspages,offsetsmeta,offsetsstrm,countermeta,counterstrm,
			     sequence_fp,chromosome_iit,index1part_aa,watsonp,index1interval,
			     genome_lc_p,fileroot,mask_lowercase_p,/*coord_values_8p*/false);
#else
      compute_positions_huge(positions4,/*positions8_high*/NULL,/*positions8_low*/NULL,
			     offsetspages,offsetsmeta,offsetsstrm,countermeta,counterstrm,
			     sequence_fp,chromosome_iit,index1part,index1interval,
			     genome_lc_p,fileroot,mask_lowercase_p,/*coord_values_8p*/false);
#endif
      fprintf(stderr,"Writing %llu genomic positions to file %s ...\n",
	      (unsigned long long) totalcounts,positionsfile_low);
      FWRITE_UINTS(positions4,totalcounts,positions_low_fp);
      fclose(positions_low_fp);
      fprintf(stderr,"done\n");

      if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	fprintf(stderr,"Error: expected file size for %s is %zu, but observed only %zu.  Trying now to write with smaller chunks.\n",
		positionsfile_low,(size_t) (totalcounts*sizeof(UINT4)),filesize);
	if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_UINTS(&(positions4[count]),count,positions_low_fp);
	}
	if (count < totalcounts) {
	  FWRITE_UINTS(&(positions4[count]),totalcounts - count,positions_low_fp);
	}
	fclose(positions_low_fp);

	if ((filesize = Access_filesize(positionsfile_low)) != (size_t) (totalcounts * sizeof(UINT4))) {
	  fprintf(stderr,"Error persists: expected file size for %s is %zu, but observed only %zu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_low,(size_t) (totalcounts*sizeof(UINT4)),filesize);
	  abort();
	}
      }

      FREE(positions4);
    }
  }

  FREE(counterstrm);
  FREE(countermeta);

#ifdef HAVE_MMAP
  munmap((void *) offsetsstrm,offsetsstrm_len);
  munmap((void *) offsetsmeta,offsetsmeta_len);
  munmap((void *) offsetspages,offsetspages_len);
#else
  FREE(offsetsstrm);
  FREE(offsetsmeta);
  FREE(offsetspages);
#endif

  return;
}
#endif


