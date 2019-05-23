static char rcsid[] = "$Id: knownsplicing.c 218286 2019-01-23 16:46:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "knownsplicing.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For qsort */

#include "assert.h"
#include "mem.h"
#include "univcoord.h"

#include "iitdef.h"
#include "interval.h"
#include "genome128_hr.h"


#define SPLICEDIST_EXTRA 10 /* Reported intron length does not include dinucleotides */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Full dump */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Splicecomp */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


#define clear_start(diff,startdiscard) (diff & (~0U << (startdiscard)))
#define clear_end(diff,enddiscard) (diff & ~(~0U << (enddiscard)))


static Genomecomp_T *splicecomp;

void
Knownsplicing_setup (Genomecomp_T *splicecomp_in) {
  splicecomp = splicecomp_in;

  return;
}


/* Previously was in splicetrie.c */
/* Modified from count_mismatches_limit */
bool
Knownsplicing_splicesite_p (Univcoord_T left, int pos5, int pos3) {
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki;
  Genomecomp_T *endblock, *ptr;
  Genomecomp_T splicesitep;

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;
  startdiscard = (left+pos5) % 32;
  enddiscard = (left+pos3) % 32;
  
  if (endblocki == startblocki) {
    debug4(printf("** Single block **\n"));
    splicesitep = splicecomp[startblocki];
    splicesitep = clear_start(splicesitep,startdiscard);
    splicesitep = clear_end(splicesitep,enddiscard);

    return (splicesitep ? true : false);

  } else if (endblocki == startblocki + 1) {
    /* Only two blocks to check */

    if (32 - startdiscard >= enddiscard) {
      /* Two blocks to check and more bits counted in startblock */
      debug4(printf("* Two blocks, start block first **\n"));

      /* 1/2: Startblock */
      splicesitep = splicecomp[startblocki];
      splicesitep = clear_start(splicesitep,startdiscard);
      if (splicesitep) {
	return true;
      }
      
      /* 2/2: Endblock */
      splicesitep = splicecomp[endblocki];
      splicesitep = clear_end(splicesitep,enddiscard);

      return (splicesitep ? true : false);

    } else {
      /* Two blocks to check and more bits counted in endblock */
      debug4(printf("** Two blocks, end block first **\n"));

      /* 1/2: Endblock */
      splicesitep = splicecomp[endblocki];
      splicesitep = clear_end(splicesitep,enddiscard);
      if (splicesitep) {
	return true;
      }

      /* 2/2: Startblock */
      splicesitep = splicecomp[startblocki];
      splicesitep = clear_start(splicesitep,startdiscard);
      return (splicesitep ? true : false);
    }

  } else {
    /* More than 2 blocks to check */
    debug4(printf("** More than two blocks **\n"));

    ptr = &(splicecomp[startblocki+1]);
    endblock = &(splicecomp[endblocki]);

    while (ptr < endblock) {
      if (*ptr) {
	return true;
      }
      ptr++;
    }

    if (enddiscard >= 32 - startdiscard) {
      /* More bits in end block */
      debug4(printf("** Final block, end block first **\n"));

      /* n/n: Go first to end block */
      splicesitep = *ptr;
      splicesitep = clear_end(splicesitep,enddiscard);
      if (splicesitep) {
	return true;
      }
      
      /* 1/n: Go second to start block */
      splicesitep = splicecomp[startblocki];
      splicesitep = clear_start(splicesitep,startdiscard);
      /* debug4(printf("adding start mask %08X\n",clear_start_mask(startdiscard))); */
      
      return (splicesitep ? true : false);

    } else {
      /* 1/n: Go first to start block */
      debug4(printf("** Final block, start block first **\n"));

      splicesitep = splicecomp[startblocki];
      splicesitep = clear_start(splicesitep,startdiscard);
      if (splicesitep) {
	return true;
      }

      /* n/n: Go second to end block */
      splicesitep = splicecomp[endblocki];
      splicesitep = clear_end(splicesitep,enddiscard);

      return (splicesitep ? true : false);
    }
  }
}


/*               87654321 */
#define UINT4_LEFT_A 0x00000000
#define UINT4_LEFT_C 0x40000000
#define UINT4_LEFT_G 0x80000000
#define UINT4_LEFT_T 0xC0000000
#define UINT4_HIGH2  0xC0000000


/* Puts leftmost character into lowest bits */
/* For right splicestrings, we want the leftmost character in the highest bits */

#if 0
static Genomecomp_T
compress16 (bool *saw_n_p, char *buffer) {
  Genomecomp_T low = 0U;
  int c;
  int i;

  /* *saw_n_p = false; -- Want to check both ref and alt, so rely on caller to set */
  for (i = 0; i < 16; i++) {
    c = buffer[i];
    low >>= 2;
    switch (c) {
    case 'A': break;
    case 'C': low |= UINT4_LEFT_C; break;
    case 'G': low |= UINT4_LEFT_G; break;
    case 'T': low |= UINT4_LEFT_T; break;
    default: *saw_n_p = true; break;
    }
  }

  return low;
}
#endif


static bool
look_for_n (bool saw_n_p, char *buffer) {
  int c;
  int i;

  /* saw_n_p = false; -- Want to check both ref and alt, so rely on caller to set */
  for (i = 0; i < 16; i++) {
    c = buffer[i];
    switch (c) {
    case 'A': case 'C': case 'G': case 'T': break;
    default: saw_n_p = true; break;
    }
  }

  return saw_n_p;
}


#if 0
static Genomecomp_T
uint4_reverse (Genomecomp_T forward) {
  Genomecomp_T reverse = 0U;
  int c;
  int i;

  for (i = 0; i < 16; i++) {
    c = forward & 0x03;
    reverse <<= 2;
    reverse |= c;
    forward >>= 2;
  }

  return reverse;
}
#endif


Univcoord_T *
Knownsplicing_retrieve_via_splicesites (bool *distances_observed_p, Genomecomp_T **splicecomp,
					Splicetype_T **splicetypes, Chrpos_T **splicedists,
					int *nsplicesites, IIT_T splicing_iit, int *splicing_divint_crosstable,
					int donor_typeint, int acceptor_typeint, Univ_IIT_T chromosome_iit,
					Genome_T genome, Genome_T genomealt, Chrpos_T shortsplicedist) {
  Univcoord_T *splicesites, chroffset, chrhigh, position;
  Chrpos_T chrlength, chrpos;
  Univcoord_T last_donor, last_antidonor, last_acceptor, last_antiacceptor;
  int last_donor_k, last_antidonor_k, last_acceptor_k, last_antiacceptor_k;
  int distance;
  int *splicesites1;
  int divno, nsplicesites1, i, k;
  Chrnum_T chrnum;
  Interval_T *intervals, interval;
  struct Interval_T *interval_structs;
  char gbuffer_ref[17], gbuffer_alt[17], *chr;
  char *restofheader;
  /* char *annot; */
  bool firstp = true, saw_n_p, allocp, alloc_header_p;
  int ntoolong = 0;

#ifdef DEBUG2
  List_T p;
#endif

  int nblocks;

  k = 0;
  for (chrnum = 1; chrnum <= Univ_IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicing_divint_crosstable[chrnum]) > 0) {
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrlength = Univ_IIT_length(chromosome_iit,chrnum); */
      splicesites1 = IIT_get_with_divno(&nsplicesites1,splicing_iit,divno,
					0U,chrlength-1U,/*sortp*/false);
      if (nsplicesites1 > 0) {
	if (firstp == true) {
	  /* annot = */ IIT_annotation(&restofheader,splicing_iit,splicesites1[0],&alloc_header_p);
	  if (restofheader[0] == '\0') {
	    fprintf(stderr,"splice distances absent...");
	    *distances_observed_p = false;
	  } else if (sscanf(restofheader,"%d",&distance) < 1) {
	    fprintf(stderr,"splice distances absent...");
	    *distances_observed_p = false;
	  } else {
	    fprintf(stderr,"splice distances present...");
	    *distances_observed_p = true;
	  }
	  if (alloc_header_p == true) {
	    FREE(restofheader);
	  }
	  firstp = false;
	}

	intervals = (Interval_T *) CALLOC(nsplicesites1,sizeof(Interval_T));
	for (i = 0; i < nsplicesites1; i++) {
	  intervals[i] = &(splicing_iit->intervals[divno][i]);
	}
	qsort(intervals,nsplicesites1,sizeof(Interval_T),Interval_cmp_low);

	last_donor = last_antidonor = last_acceptor = last_antiacceptor = 0U;
	for (i = 0; i < nsplicesites1; i++) {
	  interval = intervals[i];
	  chrpos = Interval_low(interval);
	  position = chrpos + chroffset;

	  if (position >= chrhigh) {
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	    fprintf(stderr,"\nSplice site %s:%u extends beyond chromosome length %u.  Discarding...",
		    chr,chrpos,chrlength);
	    if (allocp) FREE(chr);

	  } else if (Interval_type(interval) == donor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position != last_donor) {
		last_donor = position;
		k++;
	      }
	    } else {
	      if (position != last_antidonor) {
		last_antidonor = position;
		k++;
	      }
	    }
	  } else if (Interval_type(interval) == acceptor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position != last_acceptor) {
		last_acceptor = position;
		k++;
	      }
	    } else {
	      if (position != last_antiacceptor) {
		last_antiacceptor = position;
		k++;
	      }
	    }
	  }
	}
	FREE(intervals);
	FREE(splicesites1);
      }
    }
  }

  *nsplicesites = k;
  debug(printf("total unique splicesites: %d\n",*nsplicesites));
  fprintf(stderr,"%d unique splicesites...",*nsplicesites);

  /* The above procedure determines number of unique splicesites */

  if (*nsplicesites == 0) {
    *splicecomp = (Genomecomp_T *) NULL;
    splicesites = (Univcoord_T *) NULL;
    *splicetypes = (Splicetype_T *) NULL;
    *splicedists = (Chrpos_T *) NULL;
    return splicesites;
  }

  nblocks = (Genome_totallength(genome)+31)/32U;
  *splicecomp = (Genomecomp_T *) CALLOC(nblocks,sizeof(Genomecomp_T));
  splicesites = (Univcoord_T *) CALLOC((*nsplicesites) + 1,sizeof(Univcoord_T));
  *splicetypes = (Splicetype_T *) CALLOC(*nsplicesites,sizeof(Splicetype_T));
  *splicedists = (Chrpos_T *) CALLOC(*nsplicesites,sizeof(Chrpos_T));

  /* Use interval_structs instead of intervals, because we want to
     copy information and avoid creating many Interval_T objects */

  k = 0;
  for (chrnum = 1; chrnum <= Univ_IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicing_divint_crosstable[chrnum]) > 0) {
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      /* chrlength = Univ_IIT_length(chromosome_iit,chrnum); */
      splicesites1 = IIT_get_with_divno(&nsplicesites1,splicing_iit,divno,
					0U,chrlength-1U,/*sortp*/false);
      if (nsplicesites1 > 0) {
	interval_structs = (struct Interval_T *) CALLOC(nsplicesites1,sizeof(struct Interval_T));
	for (i = 0; i < nsplicesites1; i++) {
	  /* intervals[i] = &(splicing_iit->intervals[divno][i]); */
	  /* Copy so we can store distance information in Interval_high */
	  Interval_copy_existing(&(interval_structs[i]),&(splicing_iit->intervals[divno][i]));
	  if (*distances_observed_p == false) {
	    /* No, want to have essentially zero distance */
	    /* Interval_store_length(intervals[i],shortsplicedist); */
	  } else {
	    /* annot = */ IIT_annotation(&restofheader,splicing_iit,splicesites1[i],&alloc_header_p);
	    if (sscanf(restofheader,"%d",&distance) != 1) {
	      fprintf(stderr,"splicesites file missing distance in entry %s...exiting\n",
		      IIT_label(splicing_iit,splicesites1[i],&allocp));
	      exit(9);
	    } else if (distance < 0) {
	      fprintf(stderr,"splicesites file has a negative distance %d in entry %s...exiting\n",
		      distance,IIT_label(splicing_iit,splicesites1[i],&allocp));
	      exit(9);
	    } else if (distance > (int) shortsplicedist) {
	      ntoolong++;
	      Interval_store_length(&(interval_structs[i]),distance + SPLICEDIST_EXTRA);  /* Previously stored shortsplicedist */
	    } else {
	      Interval_store_length(&(interval_structs[i]),distance + SPLICEDIST_EXTRA);
	    }
	    if (alloc_header_p == true) {
	      FREE(restofheader);
	    }
	  }
	}

	qsort(interval_structs,nsplicesites1,sizeof(struct Interval_T),Interval_cmp_low_struct);

	last_donor = last_antidonor = last_acceptor = last_antiacceptor = 0U;
	for (i = 0; i < nsplicesites1; i++) {
	  interval = &(interval_structs[i]);
	  chrpos = Interval_low(interval);
	  position = chrpos + chroffset;

	  if (position >= chrhigh) {
#if 0
	    /* Warning given previously */
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	    fprintf(stderr,"\nSplice site %s:%u extends beyond chromosome length %u.  Discarding...",
		    chr,chrpos,chrlength);
	    if (allocp) FREE(chr);
#endif

	  } else if (Interval_type(interval) == donor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position == last_donor) {
		if (Interval_length(interval) > (*splicedists)[last_donor_k]) {
		  (*splicedists)[last_donor_k] = Interval_length(interval);
		}

	      } else {
		last_donor_k = k;
		last_donor = splicesites[k] = position;

		assert(position/32U < (unsigned int) nblocks);
		(*splicecomp)[position/32U] |= (1 << (position % 32));
		(*splicetypes)[k] = DONOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position-16,16,gbuffer_ref);
		saw_n_p = look_for_n(saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genome,genomealt,position-16,16,gbuffer_alt);
		  saw_n_p = look_for_n(saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			  chr,chrpos);
		  if (allocp) FREE(chr);

		} else {
		  k++;
		}
	      }

	    } else {
	      if (position == last_antidonor) {
		if (Interval_length(interval) > (*splicedists)[last_antidonor_k]) {
		  (*splicedists)[last_antidonor_k] = Interval_length(interval);
		}

	      } else {
		last_antidonor_k = k;
		last_antidonor = splicesites[k] = position;

		assert(position/32U < (unsigned int) nblocks);
		(*splicecomp)[position/32U] |= (1 << (position % 32));
		(*splicetypes)[k] = ANTIDONOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position,16,gbuffer_ref);
		saw_n_p = look_for_n(saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genome,genomealt,position,16,gbuffer_alt);
		  saw_n_p = look_for_n(saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			  chr,chrpos);
		  if (allocp) FREE(chr);
		} else {
		  k++;
		}
	      }
	    }

	  } else if (Interval_type(interval) == acceptor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position == last_acceptor) {
		if (Interval_length(interval) > (*splicedists)[last_acceptor_k]) {
		  (*splicedists)[last_acceptor_k] = Interval_length(interval);
		}

	      } else {
		last_acceptor_k = k;
		last_acceptor = splicesites[k] = position;

		assert(position/32U < (unsigned int) nblocks);
		(*splicecomp)[position/32U] |= (1 << (position % 32));
		(*splicetypes)[k] = ACCEPTOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position,16,gbuffer_ref);
		saw_n_p = look_for_n(saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genome,genomealt,position,16,gbuffer_alt);
		  saw_n_p = look_for_n(saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			  chr,chrpos);
		  if (allocp) FREE(chr);
		} else {
		  k++;
		}
	      }

	    } else {
	      if (position == last_antiacceptor) {
		if (Interval_length(interval) > (*splicedists)[last_antiacceptor_k]) {
		  (*splicedists)[last_antiacceptor_k] = Interval_length(interval);
		}

	      } else {
		last_antiacceptor_k = k;
		last_antiacceptor = splicesites[k] = position;

		assert(position/32U < (unsigned int) nblocks);
		(*splicecomp)[position/32U] |= (1 << (position % 32));
		(*splicetypes)[k] = ANTIACCEPTOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position-16,16,gbuffer_ref);
		saw_n_p = look_for_n(saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genome,genomealt,position-16,16,gbuffer_alt);
		  saw_n_p = look_for_n(saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...",
			  chr,chrpos);
		  if (allocp) FREE(chr);
		} else {
		  k++;
		}
	      }

	    }
	  }
	}

	FREE(interval_structs);
	FREE(splicesites1);
      }
    }
  }

  *nsplicesites = k;
  splicesites[*nsplicesites] = (Univcoord_T) -1; /* Marker for comparison in identify_all_segments */
  fprintf(stderr,"\n%d splicesites are valid...",*nsplicesites);

#ifdef DEBUG2
  for (k = 0; k < *nsplicesites; k++) {
    printf("%d: %u %s\n",k,splicesites[k],Splicetype_string((*splicetypes)[k]));
  }
#endif

  if (ntoolong > 0) {
    fprintf(stderr,"%d entries with distance > %d specified for local splice distance...",ntoolong,shortsplicedist);
  }

  return splicesites;
}

