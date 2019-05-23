static char rcsid[] = "$Id: localdb-write.c 214305 2018-03-19 23:40:43Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "localdb-write.h"

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
#include "filesuffix.h"

#include "compress-write.h"	/* For Compress_get_char */
#include "complement.h"

#include "uintlist.h"

#include "epu16-bitpack64-write.h"
#include "epu16-bitpack64-access.h"
#include "epu16-bitpack64-incr.h"


/* Another MONITOR_INTERVAL is in compress.c */
#define MONITOR_INTERVAL 100000000 /* 100 million nt */


static Oligospace_T
power (int base, int exponent) {
  Oligospace_T result = 1;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


#define BLOCKSIZE 64

/* Want to call with part 6, interval 1*/
void
Localdb_write (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
	       Alphabet_T alphabet, Width_T local1part_aa, bool watsonp,
#else
	       Width_T local1part,
#endif
	       Width_T local1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  FILE *region_fp, *ptrs_fp, *comp_fp, *positions_fp;
  char *regionfile, *pointersfile, *offsetsfile, *positionsfile;
  char *filesuffix;

  /* If offsets[oligospace] > 2^32, then will will want to allocate and write 8-mers for offsets file */
  UINT4 ascending, total_npositions, npositions;
  char *packsizes;
  UINT2 **bitpacks;

  char *comma;
  int c, nchrs, chrnum;
  Oligospace_T oligospace, bmerspace;
  Univcoord_T position = 0, adjposition, next_chrbound;
  Chrpos_T chrpos = 0U;
#ifdef HAVE_64_BIT
  Oligospace_T oligo = 0ULL;
#else
  Shortoligomer_T high = 0U, low = 0U, carry;
#endif
  Oligospace_T bmer;

  Uintlist_T *positions, p;
  UINT2 pos;


#ifdef PMAP
  int alphabet_size;
  int frame = -1, between_counter[3], in_counter[3];
  Shortoligomer_T aaindex;
  int local1part_nt = 3*local1part_aa;
#else
  Oligospace_T masked, mask;
  int between_counter, in_counter;
#endif
#ifdef DEBUG1
  char *aa;
#endif

  int circular_typeint;


#ifdef PMAP
  if (watsonp == true) {
    filesuffix = FWD_FILESUFFIX;
  } else {
    filesuffix = REV_FILESUFFIX;
  }
#else
  filesuffix = IDX_FILESUFFIX;
#endif

  regionfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen(filesuffix)+
				  /*for kmer*/2+/*for interval char*/1+
				  strlen("loctable")+1,sizeof(char));
#ifdef PMAP
  sprintf(regionfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,local1part_aa,interval_char,"loctable");
#else
  sprintf(regionfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,local1part,interval_char,"loctable");
#endif
	
  if ((region_fp = FOPEN_WRITE_BINARY(regionfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",regionfile);
    exit(9);
  } else {
    FREE(regionfile);
  }


  pointersfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				 strlen(".")+strlen(filesuffix)+
				 /*for kmer*/2+/*for interval char*/1+
				 strlen("locoffsets64meta")+1,sizeof(char));
#ifdef PMAP
  sprintf(pointersfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,local1part_aa,interval_char,"locoffsets64meta");
#else
  sprintf(pointersfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,local1part,interval_char,"locoffsets64meta");
#endif
	
  if ((ptrs_fp = FOPEN_WRITE_BINARY(pointersfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",pointersfile);
    exit(9);
  } else {
    FREE(pointersfile);
  }


  offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				strlen(".")+strlen(filesuffix)+
				/*for kmer*/2+/*for interval char*/1+
				strlen("locoffsets64strm")+1,sizeof(char));
#ifdef PMAP
  sprintf(offsetsfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,local1part_aa,interval_char,"locoffsets64strm");
#else
  sprintf(offsetsfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,local1part,interval_char,"locoffsets64strm");
#endif

  if ((comp_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",offsetsfile);
    exit(9);
  } else {
    FREE(offsetsfile);
  }


  positionsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen(filesuffix)+
				  /*for kmer*/2+/*for interval char*/1+
				  strlen("locpositions")+1,sizeof(char));
#ifdef PMAP
  sprintf(positionsfile,"%s/%s.%s.%s%d%c%s",
	  destdir,fileroot,Alphabet_string(alphabet),filesuffix,local1part_aa,interval_char,"locpositions");
#else
  sprintf(positionsfile,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,filesuffix,local1part,interval_char,"locpositions");
#endif
	
  if ((positions_fp = FOPEN_WRITE_BINARY(positionsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",positionsfile);
    exit(9);
  } else {
    FREE(positionsfile);
  }


  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  oligospace = power(alphabet_size,local1part_aa);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;

#else
#ifdef HAVE_64_BIT
  mask = ~(~0ULL << 2*local1part);
#else
  mask = ~(~0U << 2*local1part);
#endif
  oligospace = power(4,local1part);
  between_counter = in_counter = 0;
#endif

  positions = (Uintlist_T *) MALLOC(oligospace * sizeof(Uintlist_T));
  memset(positions,0,oligospace*sizeof(Uintlist_T));


  bmerspace = oligospace/BLOCKSIZE;
  fprintf(stderr,"Allocating %llu*%d bytes for packsizes\n",bmerspace,(int) sizeof(char));
  packsizes = (char *) CALLOC_NO_EXCEPTION(bmerspace,sizeof(char));
  fprintf(stderr,"Allocating %llu*%d bytes for bitpacks\n",bmerspace,(int) sizeof(UINT2 *));
  bitpacks = (UINT2 **) CALLOC_NO_EXCEPTION(bmerspace,sizeof(UINT2 *));

  if (packsizes == NULL || bitpacks == NULL) {
#ifdef PMAP
    fprintf(stderr,"Unable to allocate %llu+%llu bytes of memory, needed to build offsets with %d-mers\n",
	    bmerspace*sizeof(char),bmerspace*sizeof(UINT2 *),local1part_aa);
#else
    fprintf(stderr,"Unable to allocate %llu+%llu bytes of memory, needed to build offsets with %d-mers\n",
	    bmerspace*sizeof(char),bmerspace*sizeof(UINT2 *),local1part);
#endif
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }


  ascending = 0;
  total_npositions = 0;

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
      fprintf(stderr,"Indexing offsets/positions of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,local1part_aa,local1interval,comma);
#else
      fprintf(stderr,"Indexing offsets/positions of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,local1part,local1interval,comma);
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
    if (in_counter[frame] == local1part_aa + 1) {
      if (between_counter[frame] >= local1interval) {
#ifdef HAVE_64_BIT
	aaindex = Alphabet_get_aa_index(oligo,watsonp,local1part_nt);
#else
	aaindex = Alphabet_get_aa_index(high,low,watsonp,local1part_nt);
#endif
	bmer = aaindex/BLOCKSIZE;
	if (Epu16_bitpack64_access_filledp(aaindex,(int) packsizes[bmer],bitpacks[bmer]) == true) {
	  bitpacks[bmer] = Epu16_bitpack64_realloc_one((int) packsizes[bmer],bitpacks[bmer]);
	  packsizes[bmer] += 2;
	}
	Epu16_bitpack64_incr(aaindex,(int) packsizes[bmer],bitpacks[bmer]);

	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }

#else
    if (in_counter == local1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= local1interval
#else
	  (chrpos-local1part+1U) % local1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	positions[masked] = Uintlist_push(positions[masked],(position - local1part + 1) % 65536);

	bmer = masked/BLOCKSIZE;

	if (Epu16_bitpack64_access_filledp(masked,(int) packsizes[bmer],bitpacks[bmer]) == true) {
	  bitpacks[bmer] = Epu16_bitpack64_realloc_one((int) packsizes[bmer],bitpacks[bmer]);
	  packsizes[bmer] += 2;
	}
	Epu16_bitpack64_incr(masked,(int) packsizes[bmer],bitpacks[bmer]);

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

    position += 1;
    if (position > (Univcoord_T) (local1part - 1) && (adjposition = position - local1part + 1) % 65536 == 0) {
      /* Write offsets */
      FWRITE_UINT(ascending,region_fp);
      ascending += Epu16_bitpack64_append_differential(&npositions,ptrs_fp,comp_fp,
						       packsizes,bitpacks,oligospace);
      memset(packsizes,0,bmerspace*sizeof(char));
      for (bmer = 0; bmer < bmerspace; bmer++) {
	FREE(bitpacks[bmer]);
      }
      memset(bitpacks,0,bmerspace*sizeof(UINT2 *));


      /* Write positions */
      FWRITE_UINT(total_npositions,region_fp);
      for (masked = 0; masked < oligospace; masked++) {
	positions[masked] = Uintlist_reverse(positions[masked]);
	for (p = positions[masked]; p != NULL; p = Uintlist_next(p)) {
	  pos = (UINT2) Uintlist_head(p);
	  FWRITE_USHORT(pos,positions_fp);
	  total_npositions++;
	}
	Uintlist_free(&positions[masked]);
      }
      memset(positions,0,oligospace*sizeof(Uintlist_T));
    }
  }

  if ((adjposition = position - local1part + 1) % 65536 != 0) {
    /* Write offsets */
    FWRITE_UINT(ascending,region_fp);
    ascending += Epu16_bitpack64_append_differential(&npositions,ptrs_fp,comp_fp,
						     packsizes,bitpacks,oligospace);
    /* memset(packsizes,0,bmerspace*sizeof(char)); */
    for (bmer = 0; bmer < bmerspace; bmer++) {
      FREE(bitpacks[bmer]);
    }
    /* memset(bitpacks,0,bmerspace*sizeof(UINT2 *)); */


    /* Write positions */
    FWRITE_UINT(total_npositions,region_fp);
    for (masked = 0; masked < oligospace; masked++) {
      positions[masked] = Uintlist_reverse(positions[masked]);
      for (p = positions[masked]; p != NULL; p = Uintlist_next(p)) {
	pos = Uintlist_head(p);
	FWRITE_USHORT(pos,positions_fp);
	total_npositions++;
      }
      Uintlist_free(&positions[masked]);
    }
  }

  /* Write final values for end of strm and end of positions (used by cmetindex and atoiindex) */
  FWRITE_UINT(ascending,region_fp);
  FWRITE_UINT(total_npositions,region_fp);

  fclose(positions_fp);
  fclose(comp_fp);
  fclose(ptrs_fp);
  fclose(region_fp);

  FREE(bitpacks);
  FREE(packsizes);

  FREE(positions);

  fprintf(stderr,"Done\n");

  return;
}


