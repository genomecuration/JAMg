static char rcsid[] = "$Id: path-solve.c 218691 2019-03-19 17:38:22Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-solve.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */

#include "mem.h"
#include "assert.h"
#include "access.h"

#include "comp.h"
#include "splice.h"
#include "indel.h"
#include "intron.h"
#include "maxent_hr.h"
#include "genome128_hr.h"
#include "knownsplicing.h"

#include "univdiagdef.h"
#include "substring.h"
#include "junction.h"

#include "stage3hr.h"
#include "sedgesort.h"
#ifdef LARGE_GENOMES
#include "uint8table_rh.h"
#else
#include "uinttable_rh.h"
#endif


#define MIN_INTRONLEN 9
#define USE_LOCALDB 1
#define MAX_GMAP_PATHS 10
#define MAX_DEPTH_MIDDLE 3
#define MAX_DEPTH_LOCAL 2


/* known splicing */
#ifdef DEBUG4S
#define debug4s(x) x
#else
#define debug4s(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* best_path_genome */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif

/* via gmap */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif


static Genome_T genomebits;
static Genome_T genomebits_alt;

static bool *circularp;

static int index1part;
static int index1interval;
static int local1part;
static Localdb_T localdb;
static Localdb_T localdb2;

static int endtrim_allowed = 8;

static Chrpos_T min_intronlength;

/* Splicing */
static bool novelsplicingp;
static Univcoord_T *splicesites;
static Splicetype_T *splicetypes;
static Chrpos_T *splicedists;
static int nsplicesites;

/* For localdb */
static int max_insertionlen;
static Chrpos_T max_deletionlen;
static Chrpos_T overall_end_distance_genome;

/* For splice plus indel */
static int max_splice_deletionlen = 3;
static int max_splice_insertionlen = 3;


static int
binary_search_univcoord (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;
  
  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));
  
  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi-1,positions[highi-1],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }
  
  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


#define T Path_T
typedef struct T *T;
struct T {
  Intlist_T endpoints;
  Univcoordlist_T lefts;
  Intlist_T nmismatches;
  List_T junctions;
  
  double last_medial_splice_prob;
  double last_distal_splice_prob;
  int last_distal_knowni;
  
  bool splice5p;
  bool splice3p;
  Splicetype_T splicetype5;
  Splicetype_T splicetype3;
  double ambig_prob_5;
  double ambig_prob_3;

  Substring_T alts_substring;
};


#ifdef DEBUG13
static void
Path_print (T path, Univcoord_T chroffset) {
  printf("%s  %s  %s  alts:%p  splice5p:%d  splice3p:%d  jcns:",
	 Uintlist_to_string_offset(path->lefts,chroffset),Intlist_to_string(path->endpoints),
	 Intlist_to_string(path->nmismatches),path->alts_substring,path->splice5p,path->splice3p);
  Junction_print_list(path->junctions);
  printf("\n");
  return;
}
#endif


static T
Path_copy_5 (T old, bool splice5p, Splicetype_T splicetype5, double ambig_prob_5,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->endpoints = Intlistpool_copy(old->endpoints,intlistpool);
  new->lefts = Univcoordlistpool_copy(old->lefts,univcoordlistpool);
  new->nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
  new->junctions = Junction_copy_list(old->junctions,listpool);
  
  new->last_medial_splice_prob = old->last_medial_splice_prob;
  new->last_distal_splice_prob = old->last_distal_splice_prob;
  new->last_distal_knowni = old->last_distal_knowni;
  
  new->splice5p = splice5p;
  new->splicetype5 = splicetype5;
  new->ambig_prob_5 = ambig_prob_5;

  new->splice3p = old->splice3p;
  new->splicetype3 = old->splicetype3;
  new->ambig_prob_3 = old->ambig_prob_3;

  new->alts_substring = old->alts_substring;
  
  return new;
}

static T
Path_copy_3 (T old,  bool splice3p, Splicetype_T splicetype3, double ambig_prob_3,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->endpoints = Intlistpool_copy(old->endpoints,intlistpool);
  new->lefts = Univcoordlistpool_copy(old->lefts,univcoordlistpool);
  new->nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
  new->junctions = Junction_copy_list(old->junctions,listpool);
  
  new->last_medial_splice_prob = old->last_medial_splice_prob;
  new->last_distal_splice_prob = old->last_distal_splice_prob;
  new->last_distal_knowni = old->last_distal_knowni;
  
  new->splice5p = old->splice5p;
  new->splicetype5 = old->splicetype5;
  new->ambig_prob_5 = old->ambig_prob_5;

  new->splice3p = splice3p;
  new->splicetype3 = splicetype3;
  new->ambig_prob_3 = ambig_prob_3;

  new->alts_substring = old->alts_substring;
  
  return new;
}

static void
Path_free (T *old) {
  /* Intlist_free(&(*old)->endpoints); -- allocated by Intlistpool_push */
  /* Uintlist_free(&(*old)->lefts); -- allocated by Uintlistpool_push */
  /* Intlist_free(&(*old)->nmismatches); -- allocated by Intlistpool_push */
  Junction_gc(&(*old)->junctions);
  
  if ((*old)->alts_substring != NULL) {
    Substring_free(&(*old)->alts_substring);
  }
  
  FREE(*old);
  return;
}

static void
Path_gc (List_T *list) {
  List_T p;
  T old;
  
  for (p = *list; p != NULL; p = List_next(p)) {
    old = (T) List_head(p);
    Path_free(&old);
  }
  /* List_free(&(*list)); -- allocated by Listpool_push */
  
  return;
}


static int
Path_nsegments (T this) {
  int nsegments;

  nsegments = Univcoordlist_length(this->lefts);
  if (this->alts_substring != NULL) {
    nsegments += 1;
  }
  return nsegments;
}


/* Chooses result primarily with fewest mismatches, then with shortest
   length.  Used for single-end reads or outer parts of paired-end
   reads */
static T
find_best_splice (List_T newpaths
#ifdef DEBUG13
		  ,Univcoord_T chroffset
#endif
		  ) {
  T best_path = NULL, path;
  Chrpos_T shortest_splice_dist, splice_dist;
  int qdistal;
  int best_nmismatches, nmismatches;
  List_T p;
  
  debug13(printf("Entered find_best_splice with %d paths\n",List_length(newpaths)));

  /* First, require that all paths represent splices over the same query interval */
  path = (T) List_head(newpaths);
  for (p = newpaths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    debug13(Path_print(path,chroffset));
    if (path->alts_substring != NULL) {
      return (T) NULL;
    } else if (path->junctions == NULL) {
      return (T) NULL;
    } else if (Junction_type((Junction_T) List_head(path->junctions)) != SPLICE_JUNCTION) {
      return (T) NULL;
    }
  }
  debug13(printf("\n"));

  /* Then, pick the best one */
  path = (T) List_head(newpaths);
  qdistal = Intlist_head(path->endpoints);

  assert(Intlist_head(path->nmismatches) >= 0);
  assert(Intlist_head(Intlist_next(path->nmismatches)) >= 0);

  best_path = path;
  best_nmismatches = Intlist_head(path->nmismatches) + Intlist_head(Intlist_next(path->nmismatches));
  shortest_splice_dist = Junction_splice_distance((Junction_T) List_head(path->junctions));
  debug13(Path_print(best_path,chroffset));
  debug13(printf("path has %d nmismatches and %u splice_dist\n",best_nmismatches,shortest_splice_dist));

  for (p = List_next(newpaths); p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    if (Intlist_head(path->endpoints) != qdistal) {
      debug13(printf("qdistals are different\n"));
      return (T) NULL;
    }

    debug13(Path_print(path,chroffset));
    debug13(printf("path has %d nmismatches and %u splice_dist\n",
		   Intlist_head(path->nmismatches) + Intlist_head(Intlist_next(path->nmismatches)),
		   Junction_splice_distance((Junction_T) List_head(path->junctions))));
    if ((nmismatches = Intlist_head(path->nmismatches) + Intlist_head(Intlist_next(path->nmismatches))) <
	best_nmismatches) {
      best_path = path;
      best_nmismatches = nmismatches;
      shortest_splice_dist = Junction_splice_distance((Junction_T) List_head(path->junctions));

    } else if (nmismatches == best_nmismatches &&
	       (splice_dist = Junction_splice_distance((Junction_T) List_head(path->junctions))) <
	       shortest_splice_dist) {
      best_path = path;
      /* best_nmismatches = nmismatches; */
      shortest_splice_dist = splice_dist;
    }
  }

  debug13(printf("best path has %d nmismatches and %u splice_dist\n",best_nmismatches,shortest_splice_dist));
  return best_path;
}
      

/* Allows for differences in mismatches.  Returns splice_qpos, or -1
   if no common splice_qpos exists.  Used for inner parts of
   paired_end_reads */
static int
find_common_splice_qpos (double *common_splice_prob, int *common_nmismatches, List_T newpaths) {
  int splice_qpos, qdistal, nmismatches, ncommon = 0;
  Univcoord_T last_left;
  List_T p;
  T path;
  
  path = (T) List_head(newpaths);
  if (path->alts_substring != NULL) {
    /* Already has an alt */
    return -1;
  } else if (path->junctions == NULL) {
    /* Only a single segment */
    return -1;
  } else if (Junction_type((Junction_T) List_head(path->junctions)) != SPLICE_JUNCTION) {
    /* Not a splice */
    return -1;
  } else {
    qdistal = Intlist_head(path->endpoints);
    splice_qpos = Intlist_head(Intlist_next(path->endpoints));
    nmismatches = Intlist_head(Intlist_next(path->nmismatches));
    last_left = Univcoordlist_head(path->lefts);
    ncommon += 1;
  }
  
  for (p = List_next(newpaths); p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    if (path->alts_substring != NULL) {
      /* Already has an alt */
      return -1;
    } else if (path->junctions == NULL) {
      /* Only a single segment */
      return -1;
    } else if (Junction_type((Junction_T) List_head(path->junctions)) != SPLICE_JUNCTION) {
      /* Not a splice */
      return -1;
    } else if (Intlist_head(path->endpoints) != qdistal) {
      return -1;
    } else if (Intlist_head(Intlist_next(path->endpoints)) != splice_qpos) {
      return -1;
    } else if (Intlist_head(Intlist_next(path->nmismatches)) != nmismatches) {
      return -1;
    } else if (Univcoordlist_head(path->lefts) == last_left) {
      /* Duplicate, possible because two diagonals were extended similarly */
    } else {
      ncommon += 1;
    }
  }
  
  if (ncommon == 1) {
    return -1;
  } else {
    *common_splice_prob = path->last_medial_splice_prob;
    *common_nmismatches = nmismatches;
    return splice_qpos;
  }
}


/* qstart and qend are the genome-normalized coordinates, so qstart
   marks the left coordinate and qend marks the right coordinate.  For
   a plus-strand alignment, qstart = querystart and qend = queryend.
   For a minus-strand alignment qstart = querylength - querystart and
   qend = querylength - queryend. */

/* Need to convert from qstart and qend to querystart and queryend
   when creating Substring_T objects */

static T
Path_new_for_qstart_extension (Univcoord_T univdiagonal, int qstart, int qend,
			       bool splice5p, Splicetype_T splicetype5,
			       double ambig_prob_5, Intlistpool_T intlistpool,
			       Univcoordlistpool_T univcoordlistpool) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->endpoints = Intlistpool_push(Intlistpool_push(NULL,intlistpool,qend),intlistpool,qstart);
  new->lefts = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal);
  new->nmismatches = Intlistpool_push(NULL,intlistpool,-1);
  new->junctions = (List_T) NULL;
  
  new->last_medial_splice_prob = 0.0;
  new->last_distal_splice_prob = 0.0;
  new->last_distal_knowni = -1;
  
  new->splice5p = splice5p;
  new->splicetype5 = splicetype5;
  new->ambig_prob_5 = ambig_prob_5;

  new->splice3p = false;
  new->splicetype3 = NO_SPLICE;
  new->ambig_prob_3 = 0.0;

  new->alts_substring = (Substring_T) NULL;
  
  return new;
}


/* Returns a List_T of Path_T objects, or NULL if the diagonal cannot
   be attached.  All Path_T objects are copies of the original path */
static List_T
attach_qstart_diagonal (T path, Univcoord_T low_left, int low_qstart, int low_qend,
			Chrnum_T chrnum, Univcoord_T chroffset, int querylength,
			Spliceinfo_T spliceinfo, int *mismatch_positions_alloc, Compress_T query_compress,
			bool plusp, int genestrand, int max_mismatches_allowed,
			Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, int sensedir) {
  List_T newpaths = NULL;
  T newpath;
  Univcoord_T left;
  int nspliceends, i;
  bool splice5p;
  Splicetype_T splicetype5;
  int *trim_qstarts, *ambig_qstarts, qend, result;
  double *ambig_probs_5, ambig_prob_5;
  int nindels, indel_pos, deletionpos, splice_qpos;
  Chrpos_T splice_distance;
  double donor_prob, acceptor_prob, best_prob_i, best_prob_j, zero_prob = 0.0;
  int best_nmismatches_i, best_nmismatches_j, best_nmismatches_indel, nmismatches,
    best_knowni_i, best_knowni_j, j;
  int introntype;
#ifdef DEBUG13
  int qstart, ninserts;
#endif
  
  
  left = Univcoordlist_head(path->lefts);

  qend = Intlist_head(Intlist_next(path->endpoints)) /*+ ninserts*/;
  
#ifdef DEBUG13
  if (path->junctions == NULL) {
    ninserts = 0;
  } else {
    ninserts = Junction_ninserts(List_head(path->junctions));
  }
#endif
  debug13(qstart = Intlist_head(path->endpoints) + ninserts);
  debug13(printf("Entering attach_qstart_diagonal with sensedir %d, low_left %u %d..%d and left %u %d..%d (diff %d)\n",
		 sensedir,low_left - chroffset,low_qstart,low_qend,left - chroffset,qstart,qend,left - low_left));
  splice5p = Substring_trimmed_qstarts(&result,&splicetype5,&ambig_qstarts,&ambig_probs_5,
				       low_left,/*qstart:0,*//*qend*/low_qend,
				       plusp,genestrand,mismatch_positions_alloc,query_compress,
				       chroffset,sensedir);

  if (result < 0) {
    nspliceends = 0;
  } else if (splice5p == false) {
    nspliceends = 1;
    trim_qstarts = &result;
    ambig_probs_5 = &zero_prob;
  } else {
    nspliceends = result;
    trim_qstarts = ambig_qstarts;
  }

  debug13(printf("Got %d spliceends\n",nspliceends));
  for (i = 0; i < nspliceends; i++) {
    ambig_prob_5 = ambig_probs_5[i];
    low_qstart = trim_qstarts[i];

    if (low_left == left) {
      if (low_qstart >= Intlist_head(path->endpoints)) {
	debug13(printf("Mismatch fails, since new endpoint %d >= old endpoint %d\n",low_qstart,Intlist_head(path->endpoints)));
      } else {
	/* Mismatch: Revise the endpoint */
	debug13(printf("Mismatch extends from %d to %d\n",Intlist_head(path->endpoints),low_qstart));
	newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
    
	/* Determine nmismatches */
	nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							/*pos5*/low_qstart,
							/*pos3*/Intlist_head(Intlist_next(newpath->endpoints)) /*+ ninserts*/,
							plusp,genestrand);
	debug13(printf("Counting mismatches from %d to %d => %d\n",
		       low_qstart,Intlist_head(Intlist_next(newpath->endpoints)),nmismatches));
	Intlist_head_set(newpath->nmismatches,nmismatches); /* was -1 */
	Intlist_head_set(newpath->endpoints,low_qstart);
	newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
      }
    
    } else if (low_left > left + max_insertionlen) {
      /* Impossible */
      debug13(printf("Impossible\n"));

    } else if (low_left > left) {
      /* Insertion */
      nindels = low_left - left;
      if ((indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
						      /*left*/low_left,/*indels*/+nindels,
						      /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						      /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						      /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						      low_qstart,qend,querylength,
						      max_mismatches_allowed,/*plusp:true*/true,genestrand,
						      /*want_lowest_coordinate_p*/true)) <= 0) {
	debug13(printf("Insertion fails\n"));
      
      } else {
	newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
	
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart);
	newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(newpath->nmismatches,best_nmismatches_j);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	
	newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left);
	debug13(printf("Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       low_qstart,qend,indel_pos,nindels,best_nmismatches_i,best_nmismatches_j));
	newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
      }
    
    } else if (low_left + max_deletionlen >= left) {
      /* Deletion (or short intron) */
      nindels = left - low_left;
      if ((indel_pos = Indel_resolve_middle_deletion_or_splice(&introntype,&best_nmismatches_i,&best_nmismatches_j,
							       /*left*/low_left,/*indels*/-nindels,
							       /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							       low_qstart,qend,querylength,
							       max_mismatches_allowed,/*plusp:true*/true,genestrand,
							       min_intronlength,/*want_lowest_coordinate_p*/true)) <= 0) {
	debug13(printf("Deletion or short intron fails\n"));
	
      } else {
	newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
	
	assert(nindels >= 0);
	if ((Chrpos_T) nindels < min_intronlength || (novelsplicingp == false && splicesites == NULL) || circularp[chrnum] == true) {
	  /* Cannot be an intron, so must be a deletion */
	  deletionpos = low_left + indel_pos;
	  newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(nindels,deletionpos));
	  
	} else if ((sensedir = Intron_canonical_sensedir(introntype)) == SENSE_NULL) {
	  /* No intron dinucleotides found, so must be a deletion */
	  deletionpos = low_left + indel_pos;
	  newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(nindels,deletionpos));
	  
	} else {
	  deletionpos = low_left + indel_pos;
	  if (plusp == true) {
	    if (sensedir == SENSE_FORWARD) {
	      donor_prob = Maxent_hr_donor_prob(deletionpos,chroffset);
	      acceptor_prob = Maxent_hr_acceptor_prob(deletionpos+nindels,chroffset);
	      newpath->junctions = Listpool_push(newpath->junctions,listpool,
						 Junction_new_splice(/*splice_distance*/nindels,SENSE_FORWARD,donor_prob,acceptor_prob));
	      
	    } else {
	      donor_prob = Maxent_hr_antidonor_prob(deletionpos+nindels,chroffset);
	      acceptor_prob = Maxent_hr_antiacceptor_prob(deletionpos,chroffset);
	      newpath->junctions = Listpool_push(newpath->junctions,listpool,
						 Junction_new_splice(/*splice_distance*/nindels,SENSE_ANTI,donor_prob,acceptor_prob));
	    }
	    
	  } else {
	    /* Sense is reversed on minus strand */
	    if (sensedir == SENSE_ANTI) {
	      /* check */
	      donor_prob = Maxent_hr_antidonor_prob(deletionpos+nindels,chroffset);
	      acceptor_prob = Maxent_hr_antiacceptor_prob(deletionpos,chroffset);
	      newpath->junctions = Listpool_push(newpath->junctions,listpool,
						 Junction_new_splice(/*splice_distance*/nindels,SENSE_FORWARD,donor_prob,acceptor_prob));
	      
	    } else {
	      /* check */
	      donor_prob = Maxent_hr_donor_prob(deletionpos,chroffset);
	      acceptor_prob = Maxent_hr_acceptor_prob(deletionpos+nindels,chroffset);
	      newpath->junctions = Listpool_push(newpath->junctions,listpool,
						 Junction_new_splice(/*splice_distance*/nindels,SENSE_ANTI,donor_prob,acceptor_prob));
	    }
	  }
	}
	
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart);
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(newpath->nmismatches,best_nmismatches_j);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	
	newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left);
	debug13(printf("Deletion or short splice in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       low_qstart,qend,indel_pos,nindels,best_nmismatches_i,best_nmismatches_j));
	newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
      }
    
    } else {
      /* Splice */
      spliceinfo->segmenti_donor_nknown = spliceinfo->segmenti_antiacceptor_nknown = 0;
      if (nsplicesites > 0 &&
	  Knownsplicing_splicesite_p(low_left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search_univcoord(0,nsplicesites,splicesites,low_left);
	while (j < nsplicesites && splicesites[j] < low_left + querylength) {
	  if (splicetypes[j] == DONOR) {
	    debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	    spliceinfo->segmenti_donor_knownpos[spliceinfo->segmenti_donor_nknown] = splicesites[j] - low_left;
	    spliceinfo->segmenti_donor_knowni[spliceinfo->segmenti_donor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIACCEPTOR) {
	    debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	    spliceinfo->segmenti_antiacceptor_knownpos[spliceinfo->segmenti_antiacceptor_nknown] = splicesites[j] - low_left;
	    spliceinfo->segmenti_antiacceptor_knowni[spliceinfo->segmenti_antiacceptor_nknown++] = j;
	  }
	  j++;
	}
      }
      spliceinfo->segmenti_donor_knownpos[spliceinfo->segmenti_donor_nknown] = querylength + 100;
      spliceinfo->segmenti_antiacceptor_knownpos[spliceinfo->segmenti_antiacceptor_nknown] = querylength + 100;
      
      spliceinfo->segmentj_acceptor_nknown = spliceinfo->segmentj_antidonor_nknown = 0;
      if (nsplicesites > 0 &&
	  Knownsplicing_splicesite_p(left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search_univcoord(0,nsplicesites,splicesites,left);
	while (j < nsplicesites && splicesites[j] < left + querylength) {
	  if (splicetypes[j] == ACCEPTOR) {
	    debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	    spliceinfo->segmentj_acceptor_knownpos[spliceinfo->segmentj_acceptor_nknown] = splicesites[j] - left;
	    spliceinfo->segmentj_acceptor_knowni[spliceinfo->segmentj_acceptor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIDONOR) {
	    debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	    spliceinfo->segmentj_antidonor_knownpos[spliceinfo->segmentj_antidonor_nknown] = splicesites[j] - left;
	    spliceinfo->segmentj_antidonor_knowni[spliceinfo->segmentj_antidonor_nknown++] = j;
	  }
	  j++;
	}
      }
      spliceinfo->segmentj_acceptor_knownpos[spliceinfo->segmentj_acceptor_nknown] = querylength + 100;
      spliceinfo->segmentj_antidonor_knownpos[spliceinfo->segmentj_antidonor_nknown] = querylength + 100;
      
      if (sensedir == SENSE_FORWARD) {
	if ((splice_qpos = Splice_resolve_sense(&nindels,&indel_pos,&best_knowni_i,&best_knowni_j,
						&best_nmismatches_i,&best_nmismatches_j,&best_nmismatches_indel,
						&best_prob_i,&best_prob_j,
						/*segmenti_left*/low_left,/*segmentj_left*/left,chroffset,chroffset,
						low_qstart,qend,querylength,query_compress,
						spliceinfo,max_mismatches_allowed,plusp,genestrand,
						max_splice_deletionlen,max_splice_insertionlen,
						/*allow_indel_p*/true)) <= 0) {
	  debug13(printf("Splice_resolve_sense: fails\n"));
	
	} else if (nindels != 0 && splice_qpos < indel_pos) { /* splice is distal, indel is medial */
	  /* Push indel (based on left) then splice */
	  splice_distance = left - low_left + nindels;

	  debug13(printf("Splice_resolve_sense: %d indels at %d then splice at %d\n",nindels,indel_pos,splice_qpos));
	  newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
	  newpath->last_medial_splice_prob = best_prob_j;
	  newpath->last_distal_splice_prob = best_prob_i;
	  newpath->last_distal_knowni = best_knowni_i;
	  
	  Intlist_head_set(newpath->endpoints,indel_pos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart);
	  
	  /* Indel first */
	  if (nindels > 0) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	  } else {
	    deletionpos = low_left + indel_pos;
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(-nindels,deletionpos));
	  }
	  
	  /* Splice second */
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  }
	  
	  /* For qstart, push j first, then push i */
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_indel);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,left + nindels); /* ? check */
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
	  
	} else if (nindels != 0 && indel_pos < splice_qpos) { /* indel is distal, splice is medial */
	  /* Push splice then indel (based on low_left) */
	  splice_distance = left - low_left + nindels;

	  debug13(printf("Splice_resolve_sense: splice at %d then %d indels at %d\n",splice_qpos,nindels,indel_pos));
	  newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
	  Intlist_head_set(newpath->endpoints,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,indel_pos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart);
	  
	  /* Splice first */
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  }
	  
	  /* Indel second */
	  if (nindels > 0) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	  } else {
	    deletionpos = low_left + indel_pos;
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(-nindels,deletionpos));
	  }
	  
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_j);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_indel);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left - nindels); /* ? check */
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);

	} else {
	  /* Splice only */
	  splice_distance = left - low_left;

	  newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
	  newpath->last_medial_splice_prob = best_prob_j;
	  newpath->last_distal_splice_prob = best_prob_i;
	  newpath->last_distal_knowni = best_knowni_i;
	  
	  Intlist_head_set(newpath->endpoints,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart);
	  
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  }
	  
	  /* For qstart, push j first, then push i */
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_j);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
	}
      
      } else {
	if ((splice_qpos = Splice_resolve_antisense(&nindels,&indel_pos,&best_knowni_i,&best_knowni_j,
						    &best_nmismatches_i,&best_nmismatches_j,&best_nmismatches_indel,
						    &best_prob_i,&best_prob_j,
						    /*segmenti_left*/low_left,/*segmentj_left*/left,chroffset,chroffset,
						    low_qstart,qend,querylength,query_compress,
						    spliceinfo,max_mismatches_allowed,plusp,genestrand,
						    max_splice_deletionlen,max_splice_insertionlen,
						    /*allow_indel_p*/true)) <= 0) {
	  debug13(printf("Splice_resolve_antisense: fails\n"));

	} else if (nindels != 0 && splice_qpos < indel_pos) { /* splice is distal, indel is medial */
	  /* Push indel (based on left) then splice */
	  splice_distance = left - low_left + nindels;

	  debug13(printf("Splice_resolve_antisense: %d indels at %d then splice at %d\n",nindels,indel_pos,splice_qpos));
	  newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
	  newpath->last_medial_splice_prob = best_prob_j;
	  newpath->last_distal_splice_prob = best_prob_i;
	  newpath->last_distal_knowni = best_knowni_i;
	  
	  Intlist_head_set(newpath->endpoints,indel_pos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart);
	  
	  /* Indel first */
	  if (nindels > 0) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	  } else {
	    deletionpos = low_left + indel_pos;
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(-nindels,deletionpos));
	  }
	  
	  /* Splice second */
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  }
	  
	  /* For qstart, push j first, then push i */
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_indel);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,left + nindels); /* ? check */
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);

	} else if (nindels != 0 && indel_pos < splice_qpos) { /* indel is distal, splice is medial */
	  /* Push splice then indel (based on low_left) */
	  splice_distance = left - low_left + nindels;

	  debug13(printf("Splice_resolve_antisense: splice at %d then %d indels at %d\n",splice_qpos,nindels,indel_pos));
	  newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
	  Intlist_head_set(newpath->endpoints,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,indel_pos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart);
	  
	  /* Splice first */
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  }
	  
	  /* Indel second */
	  if (nindels > 0) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	  } else {
	    deletionpos = low_left + indel_pos;
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(-nindels,deletionpos));
	  }
	  
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_j);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_indel);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left - nindels); /* ? check */
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);

	} else {
	  /* Splice only */
	  splice_distance = left - low_left;

	  newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,intlistpool,univcoordlistpool,listpool);
	  newpath->last_medial_splice_prob = best_prob_j;
	  newpath->last_distal_splice_prob = best_prob_i;
	  newpath->last_distal_knowni = best_knowni_i;
	  
	  Intlist_head_set(newpath->endpoints,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart);
	  
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   low_qstart,qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  }
	  
	  /* For qstart, push j first, then push i */
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_j);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,low_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
	}
      }
    }
  }

  if (splice5p == true) {
    FREE(ambig_qstarts);
    FREE(ambig_probs_5);
  }

  return newpaths;
}


static Substring_T
combine_into_qstart_alts (List_T newpaths, int common_splice_qpos, double common_splice_prob,
			  int querylength, bool plusp, int genestrand, Compress_T query_compress,
			  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			  int sensedir) {
  T newpath;
  List_T p;
  Univcoord_T *alts_coords;
  int *alts_nmismatches, *alts_knowni;
  double *alts_probs;
  int n, k;

  n = List_length(newpaths);
  alts_coords = (Univcoord_T *) MALLOC_OUT(n*sizeof(Univcoord_T));
  alts_knowni = (int *) MALLOC_OUT(n*sizeof(int));
  alts_nmismatches = (int *) MALLOC_OUT(n*sizeof(int));
  alts_probs = (double *) MALLOC_OUT(n*sizeof(double));

  for (p = newpaths, k = 0; p != NULL; p = List_next(p), k++) {
    newpath = List_head(p);
    alts_coords[k] = Univcoordlist_head(newpath->lefts) + common_splice_qpos;
    alts_nmismatches[k] = Intlist_head(newpath->nmismatches);
    alts_knowni[k] = newpath->last_distal_knowni;
    alts_probs[k] = newpath->last_distal_splice_prob;
    debug13(printf("Start alt %u [%u]\n",alts_coords[k],alts_coords[k] - chroffset));
  }

  /* qstart = Intlist_head(newpath->endpoints); -- Should all be the same */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      return Substring_new_alts_D(/*querystart*/0,/*queryend*/common_splice_qpos,common_splice_qpos,
				  querylength,plusp,genestrand,query_compress,
				  chrnum,chroffset,chrhigh,chrlength,
				  alts_coords,alts_knowni,alts_nmismatches,alts_probs,n,
				  common_splice_prob,/*substring1p*/true);
    } else {
      return Substring_new_alts_A(/*querystart*/querylength - common_splice_qpos,
				  /*queryend*/querylength,common_splice_qpos,
				  querylength,plusp,genestrand,query_compress,
				  chrnum,chroffset,chrhigh,chrlength,
				  alts_coords,alts_knowni,alts_nmismatches,alts_probs,n,
				  common_splice_prob,/*substring1p*/false);
    }
  } else {
    if (plusp == true) {
      return Substring_new_alts_A(/*querystart*/0,/*queryend*/common_splice_qpos,common_splice_qpos,
				  querylength,plusp,genestrand,query_compress,
				  chrnum,chroffset,chrhigh,chrlength,
				  alts_coords,alts_knowni,alts_nmismatches,alts_probs,n,
				  common_splice_prob,/*substring1p*/true);
    } else {
      return Substring_new_alts_D(/*querystart*/querylength - common_splice_qpos,
				  /*queryend*/querylength,common_splice_qpos,
				  querylength,plusp,genestrand,query_compress,
				  chrnum,chroffset,chrhigh,chrlength,
				  alts_coords,alts_knowni,alts_nmismatches,alts_probs,n,
				  common_splice_prob,/*substring1p*/false);
    }
  }
}


/* Returns a list of Path_T objects */
static List_T
compute_qstart_paths (int depth, bool *terminalp, T path, List_T diagonals,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      int querylength, Spliceinfo_T spliceinfo, int *mismatch_positions_alloc, Compress_T query_compress,
		      bool plusp, int genestrand, int max_mismatches_allowed,
		      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		      Listpool_T listpool, Univdiagpool_T univdiagpool, int sensedir, bool innerp) {
  List_T all_child_paths, child_paths, newpaths, newdiagonals, p;
  T best_child_path, child_path;
  bool all_terminal_p, child_terminal_p;
  int common_splice_qpos, qstart, ignore_int;
  Univcoord_T ignore_univcoord;
  Junction_T junction;
  double common_splice_prob;
  int common_nmismatches;
  Univdiag_T diagonal;

#ifdef DEBUG13
  printf("Entered compute_qstart_paths with innerp %d.  Current path: ",innerp);
  Path_print(path,chroffset);
  printf("Diagonals:\n");
  for (p = diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    printf("%u %d..%d\n",diagonal->univdiagonal - chroffset,diagonal->qstart,diagonal->qend);
  }
  printf("\n");
#endif

  if (diagonals == NULL) {
    debug13(printf("diagonals is NULL, so exiting.  terminalp is true\n"));
    *terminalp = true;
    return Listpool_push(NULL,listpool,(void *) path);

  } else if (depth > MAX_DEPTH_MIDDLE) {
    debug13(printf("reached maximum depth, so exiting\n"));
    *terminalp = false;
    return Listpool_push(NULL,listpool,(void *) path);

  } else {
    /* Partition diagonals into those that can be added, and the rest */
    newpaths = (List_T) NULL;
    newdiagonals = (List_T) NULL;

    qstart = Intlist_head(path->endpoints);
    for (p = diagonals; p != NULL; p = List_next(p)) {
      diagonal = (Univdiag_T) List_head(p);
      if (diagonal->qstart >= qstart) {
	/* Skip */
      } else if ((child_paths = attach_qstart_diagonal(path,/*low_left*/diagonal->univdiagonal,
						       /*low_qstart*/diagonal->qstart,
						       /*low_qend*/diagonal->qend,
						       chrnum,chroffset,querylength,
						       spliceinfo,mismatch_positions_alloc,query_compress,
						       plusp,genestrand,max_mismatches_allowed,intlistpool,
						       univcoordlistpool,listpool,sensedir)) != NULL) {
	/* *extendedp = true; */
	newpaths = List_append(child_paths,newpaths);
      } else {
	/* Save for the next iteration */
	newdiagonals = Univdiagpool_push_existing(newdiagonals,univdiagpool,diagonal);
      }
    }

    if (newpaths == NULL) {
      debug13(printf("newpaths is NULL, so exiting.  terminalp is true\n"));
      *terminalp = true;
      return Listpool_push(NULL,listpool,(void *) path);

    } else if (List_next(newpaths) == NULL) {
      debug13(printf("newpaths has one entry, so exiting\n"));
      *terminalp = false;
      /* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
      Path_free(&path);
      return newpaths;

    } else {
      all_terminal_p = true;
      all_child_paths = (List_T) NULL;
      for (p = newpaths; p != NULL; p = List_next(p)) {
	debug13(printf("Calling recursively on "));
	debug13(Path_print(List_head(p),chroffset));
	debug13(printf("\n"));
	child_paths = compute_qstart_paths(depth+1,&child_terminal_p,List_head(p),newdiagonals,
					   chrnum,chroffset,chrhigh,chrlength,querylength,
					   spliceinfo,mismatch_positions_alloc,query_compress,plusp,genestrand,
					   max_mismatches_allowed,intlistpool,univcoordlistpool,
					   listpool,univdiagpool,sensedir,innerp);
	debug13(printf("Returning recursively\n"));
	all_child_paths = List_append(all_child_paths,child_paths);
	if (child_terminal_p == false) {
	  all_terminal_p = false;
	}
      }

      if (all_terminal_p == false) {
	debug13(printf("Not all children are terminals, so exiting\n"));
	/* List_free(&newpaths); -- allocated by Listpool_push */
	/* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
	*terminalp = false;
	Path_free(&path);
	return all_child_paths;

      } else if (innerp == false &&
#ifdef DEBUG13
		 (best_child_path = find_best_splice(all_child_paths,chroffset))
#else
		 (best_child_path = find_best_splice(all_child_paths))
#endif
		 != NULL) {
	debug13(printf("Best child path:\n"));
	debug13(Path_print(best_child_path,chroffset));
	for (p = all_child_paths; p != NULL; p = List_next(p)) {
	  child_path = (T) List_head(p);
	  if (child_path != best_child_path) {
	    Path_free(&child_path);
	  }
	}
	*terminalp = false;
	Path_free(&path);
	return Listpool_push(NULL,listpool,(void *) best_child_path);

      } else if ((common_splice_qpos = find_common_splice_qpos(&common_splice_prob,&common_nmismatches,all_child_paths)) < 0) {
	debug13(printf("No common splice pos found, so exiting\n"));
	/* List_free(&newpaths); -- allocated by Listpool_push */
	/* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
	*terminalp = false;
	Path_free(&path);
	return all_child_paths;

      } else if (common_splice_prob < 0.9) {
	debug13(printf("Common splice prob %f is poor, so exiting\n",common_splice_prob));
	/* List_free(&newpaths); -- allocated by Listpool_push */
	/* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
	*terminalp = false;
	Path_free(&path);
	return all_child_paths;

#if 0
      } else if ((best_child_path = dominating_child_path(all_child_paths)) != NULL) {
	debug13(printf("Found a single best path\n"));
	for (p = all_child_paths; p != NULL; p = List_next(p)) {
	  if ((child_path = List_head(p)) != best_child_path) {
	    Path_free(&child_path);
	  }
	}
	*terminalp = false;
	Path_free(&path);
	return Listpool_push(NULL,listpool,(void *) best_child_path);
#endif

      } else {
	debug13(printf("All %d children are terminals with a common splice pos of %d, so combining into alt\n",
		       List_length(all_child_paths),common_splice_qpos));
#ifdef DEBUG13
	best_child_path = find_best_splice(all_child_paths,chroffset);
#else
	best_child_path = find_best_splice(all_child_paths);
#endif
	assert(best_child_path != NULL);
	debug13(printf("Best child path:\n"));
	debug13(Path_print(best_child_path,chroffset));

	best_child_path->alts_substring = combine_into_qstart_alts(all_child_paths,common_splice_qpos,common_splice_prob,
								   querylength,plusp,genestrand,query_compress,
								   chrnum,chroffset,chrhigh,chrlength,sensedir);

	best_child_path->endpoints = Intlistpool_pop(best_child_path->endpoints,&ignore_int);
	best_child_path->nmismatches = Intlistpool_pop(best_child_path->nmismatches,&ignore_int);
	best_child_path->lefts = Univcoordlistpool_pop(best_child_path->lefts,&ignore_univcoord);
	best_child_path->junctions = Listpool_pop(best_child_path->junctions,(void **) &junction);
	Junction_free(&junction);
	
	debug13(printf("(1) Path after alt: "));
	debug13(Path_print(best_child_path,chroffset));

	for (p = all_child_paths; p != NULL; p = List_next(p)) {
	  child_path = (T) List_head(p);
	  if (child_path != best_child_path) {
	    Path_free(&child_path);
	  }
	}
	*terminalp = false;
	Path_free(&path);

	/* List_free(&newpaths); -- allocated by Listpool_push */
	/* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
	return Listpool_push(NULL,listpool,(void *) best_child_path);
      }
    }
  }
}


static void
check_for_ascending_values (Univcoord_T *diagonals, int ndiagonals) {
  int i;

  for (i = 0; i < ndiagonals - 1; i++) {
    if (diagonals[i+1] < diagonals[i]) {
      abort();
    }
  }
  return;
}


static List_T
add_qstart_local (bool *addedp, T path, Univcoord_T *local_diagonals, int ndiagonals, int pos3,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		  int querylength, Spliceinfo_T spliceinfo, int *mismatch_positions_alloc,
		  Compress_T query_compress, bool plusp, int genestrand, int max_mismatches_allowed,
		  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, int sensedir, bool innerp) {
  List_T newpaths = NULL, paths, p;
  T best_child_path, child_path;
  int common_splice_qpos;
  int common_nmismatches;
  double common_splice_prob;
  int nhits, i;

  void *ignore;
  int ignore_int;
  Univcoord_T ignore_univcoord;
  Junction_T junction;


  debug13(printf("Entered add_qstart_local.  Current path: "));
  debug13(Path_print(path,chroffset));
  debug13(printf("\n"));
#ifdef CHECK_ASSERTIONS
  check_for_ascending_values(local_diagonals,ndiagonals);
#endif

  *addedp = false;

  i = ndiagonals - 1;
  while (i >= 0 /*&& nhits < 2*/) {
    debug13(printf("attaching qstart local %u %d..%d\n",local_diagonals[i] - chroffset,0,pos3));
    if ((paths = attach_qstart_diagonal(path,/*low_left*/local_diagonals[i],
					/*low_qstart*/0,/*low_qend*/pos3,
					chrnum,chroffset,querylength,
					spliceinfo,mismatch_positions_alloc,query_compress,
					plusp,genestrand,max_mismatches_allowed,intlistpool,
					univcoordlistpool,listpool,sensedir)) != NULL) {
      newpaths = List_append(paths,newpaths);
      *addedp = true;
    }
    i--;
  }

  if ((nhits = List_length(newpaths)) == 0) {
    debug13(printf("newpaths is NULL, so exiting\n"));
    return Listpool_push(NULL,listpool,(void *) path);

  } else if (nhits == 1) {
    debug13(printf("newpaths has one entry, so exiting\n"));
    Path_free(&path);
    return newpaths;

  } else if (innerp == false &&
#ifdef DEBUG13
	     (best_child_path = find_best_splice(newpaths,chroffset))
#else
	     (best_child_path = find_best_splice(newpaths))
#endif
	     != NULL) {
    for (p = newpaths; p != NULL; p = List_next(p)) {
      child_path = (T) List_head(p);
      if (child_path != best_child_path) {
	Path_free(&child_path);
      }
    }
    Path_free(&path);
    return Listpool_push(NULL,listpool,(void *) best_child_path);

  } else if ((common_splice_qpos = find_common_splice_qpos(&common_splice_prob,&common_nmismatches,newpaths)) < 0) {
    debug13(printf("No common splice pos found, so exiting\n"));
    Path_free(&path);
    return newpaths;

  } else if (common_splice_prob < 0.9) {
    debug13(printf("Common splice prob %f is poor, so exiting\n",common_splice_prob));
    Path_free(&path);
    return newpaths;

  } else {
    debug13(printf("Creating an alt substring\n"));
#if 0
    Intlist_head_set(path->endpoints,common_splice_qpos);
    Intlist_head_set(path->nmismatches,common_nmismatches);
#else
    Path_free(&path);
    path = (Path_T) List_head(newpaths);
    path->endpoints = Intlistpool_pop(path->endpoints,&ignore_int);
    path->nmismatches = Intlistpool_pop(path->nmismatches,&ignore_int);
    path->lefts = Univcoordlistpool_pop(path->lefts,&ignore_univcoord);
    if (path->junctions != NULL) {
      path->junctions = Listpool_pop(path->junctions,(void **) &junction);
      Junction_free(&junction);
    }
#endif

    path->alts_substring = combine_into_qstart_alts(newpaths,common_splice_qpos,common_splice_prob,
						    querylength,plusp,genestrand,query_compress,
						    chrnum,chroffset,chrhigh,chrlength,sensedir);
    debug13(printf("(2) Resulting path: "));
    debug13(Path_print(path,chroffset));
    newpaths = Listpool_pop(newpaths,(void **) &ignore);
    Path_gc(&newpaths);
    return Listpool_push(NULL,listpool,(void *) path);
  }
}


static T
Path_new_for_qend_extension (Univcoord_T univdiagonal, int qstart, int qend, bool splice3p, Splicetype_T splicetype3,
			     double ambig_prob_3, Intlistpool_T intlistpool,
			     Univcoordlistpool_T univcoordlistpool) {
  T new = (T) MALLOC(sizeof(*new));

  new->endpoints = Intlistpool_push(Intlistpool_push(NULL,intlistpool,qstart),intlistpool,qend);
  new->lefts = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal);
  new->nmismatches = Intlistpool_push(NULL,intlistpool,-1);
  new->junctions = (List_T) NULL;

  new->last_medial_splice_prob = 0.0;
  new->last_distal_splice_prob = 0.0;
  new->last_distal_knowni = -1;

  new->splice5p = false;
  new->splicetype5 = NO_SPLICE;
  new->ambig_prob_5 = 0.0;

  new->splice3p = splice3p;
  new->splicetype3 = splicetype3;
  new->ambig_prob_3 = ambig_prob_3;

  new->alts_substring = (Substring_T) NULL;

  return new;
}


/* Returns a List_T of Path_T objects, or NULL if the diagonal cannot
   be attached.  All Path_T objects are copies of the original path */
static List_T
attach_qend_diagonal (T path, Univcoord_T high_left, int high_qstart, int high_qend,
		      Chrnum_T chrnum, Univcoord_T chroffset, int querylength,
		      Spliceinfo_T spliceinfo, int *mismatch_positions_alloc, Compress_T query_compress,
		      bool plusp, int genestrand, int max_mismatches_allowed,
		      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		      Listpool_T listpool, int sensedir) {
  List_T newpaths = NULL;
  T newpath;
  Univcoord_T left;
  int nspliceends, i;
  bool splice3p;
  Splicetype_T splicetype3;
  int *trim_qends, *ambig_qends, qstart, ninserts, result;
  double *ambig_probs_3, ambig_prob_3;
  int nindels, indel_pos, deletionpos, splice_qpos;
  Chrpos_T splice_distance;
  double donor_prob, acceptor_prob, best_prob_i, best_prob_j, zero_prob = 0.0;
  int best_nmismatches_i, best_nmismatches_j, best_nmismatches_indel, nmismatches,
    best_knowni_i, best_knowni_j, j;
  int introntype;
#ifdef DEBUG13
  int qend;
#endif


  left = Univcoordlist_head(path->lefts);

#if 0
  ninserts = Junction_total_ninserts(path->junctions);
#else
  if (path->junctions == NULL) {
    ninserts = 0;
  } else {
    ninserts = Junction_ninserts(List_head(path->junctions));
  }
#endif

  
  qstart = Intlist_head(Intlist_next(path->endpoints)) + ninserts;

  debug13(qend = Intlist_head(path->endpoints) /*+ ninserts*/);
  debug13(printf("Entering attach_qend_diagonal with sensedir %d, left %u %d..%d and high_left %u %d..%d (diff %d)\n",
		 sensedir,left - chroffset,qstart,qend,high_left - chroffset,high_qstart,high_qend,high_left - left));

  splice3p = Substring_trimmed_qends(&result,&splicetype3,&ambig_qends,&ambig_probs_3,
				     high_left,/*qstart*/high_qstart,
				     /*qend:querylength,*/querylength,plusp,genestrand,
				     mismatch_positions_alloc,query_compress,
				     chroffset,sensedir);

  if (result < 0) {
    nspliceends = 0;
  } else if (splice3p == false) {
    nspliceends = 1;
    trim_qends = &result;
    ambig_probs_3 = &zero_prob;
  } else {
    nspliceends = result;
    trim_qends = ambig_qends;
  }

  debug13(printf("Got %d spliceends\n",nspliceends));
  for (i = 0; i < nspliceends; i++) {
    ambig_prob_3 = ambig_probs_3[i];
    high_qend = trim_qends[i];

    if (high_left == left) {
      if (high_qend <= Intlist_head(path->endpoints)) {
	debug13(printf("Mismatch fails, since new endpoint %d <= old endpoint %d\n",high_qend,Intlist_head(path->endpoints)));
      } else {
	/* Mismatch: Revise the endpoint */
	debug13(printf("Mismatch extends from %d to %d\n",Intlist_head(path->endpoints),high_qend));
	newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);

	/* Determine nmismatches */
	if (path->junctions == NULL) {
	  ninserts = 0;
	} else {
	  ninserts = Junction_ninserts(List_head(path->junctions));
	}
	nmismatches = Genome_count_mismatches_substring(genomebits,genomebits_alt,query_compress,left,
							/*pos5*/Intlist_head(Intlist_next(newpath->endpoints)) + ninserts,
							/*pos3*/high_qend,plusp,genestrand);
	debug13(printf("Counting mismatches from %d to %d => %d\n",
		       Intlist_head(Intlist_next(newpath->endpoints)),high_qend,nmismatches));
	Intlist_head_set(newpath->nmismatches,nmismatches); /* was -1 */
	Intlist_head_set(newpath->endpoints,high_qend);
	newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
      }

    } else if (high_left + max_insertionlen < left) {
      /* Impossible */
      debug13(printf("Impossible\n"));

    } else if (high_left < left) {
      /* Insertion */
      nindels = left - high_left;
      if ((indel_pos = Indel_resolve_middle_insertion(&best_nmismatches_i,&best_nmismatches_j,
						      /*left*/left,/*indels*/+nindels,
						      /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						      /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						      /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						      qstart,high_qend,querylength,
						      max_mismatches_allowed,/*plusp:true*/true,genestrand,
						      /*want_lowest_coordinate_p*/true)) <= 0) {
	debug13(printf("Insertion fails\n"));

      } else {
	newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);
	
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend);
	newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(newpath->nmismatches,best_nmismatches_i);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	
	newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left);
	debug13(printf("Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       qstart,high_qend,indel_pos,nindels,best_nmismatches_i,best_nmismatches_j));
	newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
      }

    } else if (high_left <= left + max_deletionlen) {
      /* Deletion (or short intron) */
      nindels = high_left - left;
      if ((indel_pos = Indel_resolve_middle_deletion_or_splice(&introntype,&best_nmismatches_i,&best_nmismatches_j,
							       /*left*/left,/*indels*/-nindels,
							       /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							       qstart,high_qend,querylength,
							       max_mismatches_allowed,/*plusp:true*/true,genestrand,
							       min_intronlength,/*want_lowest_coordinate_p*/true)) <= 0) {
	debug13(printf("Deletion or short intron fails\n"));
      
      } else {
	newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);
	
	assert(nindels >= 0);
	if ((Chrpos_T) nindels < min_intronlength || (novelsplicingp == false && splicesites == NULL) || circularp[chrnum] == true) {
	  /* Cannot be an intron, so must be a deletion */
	  deletionpos = left + indel_pos;
	  newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(nindels,deletionpos));
	  
	} else if ((sensedir = Intron_canonical_sensedir(introntype)) == SENSE_NULL) {
	  /* No intron dinucleotides found, so must be a deletion */
	  deletionpos = left + indel_pos;
	  newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(nindels,deletionpos));
	  
	} else {
	  deletionpos = left + indel_pos;
	  if (plusp == true) {
	    if (sensedir == SENSE_FORWARD) {
	      donor_prob = Maxent_hr_donor_prob(deletionpos,chroffset);
	      acceptor_prob = Maxent_hr_acceptor_prob(deletionpos+nindels,chroffset);
	      newpath->junctions = Listpool_push(newpath->junctions,listpool,
						 Junction_new_splice(/*splice_distance*/nindels,SENSE_FORWARD,donor_prob,acceptor_prob));
	      
	    } else {
	      donor_prob = Maxent_hr_antidonor_prob(deletionpos+nindels,chroffset);
	      acceptor_prob = Maxent_hr_antiacceptor_prob(deletionpos,chroffset);
	      newpath->junctions = Listpool_push(newpath->junctions,listpool,
						 Junction_new_splice(/*splice_distance*/nindels,SENSE_ANTI,donor_prob,acceptor_prob));
	    }
	    
	  } else {
	    /* Sense is reversed on minus strand */
	    if (sensedir == SENSE_ANTI) {
	      /* check */
	      donor_prob = Maxent_hr_antidonor_prob(deletionpos+nindels,chroffset);
	      acceptor_prob = Maxent_hr_antiacceptor_prob(deletionpos,chroffset);
	      newpath->junctions = Listpool_push(newpath->junctions,listpool,
						 Junction_new_splice(/*splice_distance*/nindels,SENSE_FORWARD,donor_prob,acceptor_prob));
	      
	    } else {
	      /* check */
	      donor_prob = Maxent_hr_donor_prob(deletionpos,chroffset);
	      acceptor_prob = Maxent_hr_acceptor_prob(deletionpos+nindels,chroffset);
	      newpath->junctions = Listpool_push(newpath->junctions,listpool,
						 Junction_new_splice(/*splice_distance*/nindels,SENSE_ANTI,donor_prob,acceptor_prob));
	    }
	  }
	}
	
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend);
	
	/* For qend, push i first, then push j */
	Intlist_head_set(newpath->nmismatches,best_nmismatches_i);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	
	newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left);
	debug13(printf("Deletion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       qstart,high_qend,indel_pos,nindels,best_nmismatches_i,best_nmismatches_j));
	newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
      }
    
    } else {
      /* Splice */
      spliceinfo->segmenti_donor_nknown = spliceinfo->segmenti_antiacceptor_nknown = 0;
      if (nsplicesites > 0 &&
	  Knownsplicing_splicesite_p(left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search_univcoord(0,nsplicesites,splicesites,left);
	while (j < nsplicesites && splicesites[j] < left + querylength) {
	  if (splicetypes[j] == DONOR) {
	    debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	    spliceinfo->segmenti_donor_knownpos[spliceinfo->segmenti_donor_nknown] = splicesites[j] - left;
	    spliceinfo->segmenti_donor_knowni[spliceinfo->segmenti_donor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIACCEPTOR) {
	    debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	    spliceinfo->segmenti_antiacceptor_knownpos[spliceinfo->segmenti_antiacceptor_nknown] = splicesites[j] - left;
	    spliceinfo->segmenti_antiacceptor_knowni[spliceinfo->segmenti_antiacceptor_nknown++] = j;
	  }
	  j++;
	}
      }
      spliceinfo->segmenti_donor_knownpos[spliceinfo->segmenti_donor_nknown] = querylength + 100;
      spliceinfo->segmenti_antiacceptor_knownpos[spliceinfo->segmenti_antiacceptor_nknown] = querylength + 100;
      
      spliceinfo->segmentj_acceptor_nknown = spliceinfo->segmentj_antidonor_nknown = 0;
      if (nsplicesites > 0 &&
	  Knownsplicing_splicesite_p(high_left,/*pos5*/1,/*pos3*/querylength) == true) {
	j = binary_search_univcoord(0,nsplicesites,splicesites,high_left);
	while (j < nsplicesites && splicesites[j] < high_left + querylength) {
	  if (splicetypes[j] == ACCEPTOR) {
	    debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	    spliceinfo->segmentj_acceptor_knownpos[spliceinfo->segmentj_acceptor_nknown] = splicesites[j] - high_left;
	    spliceinfo->segmentj_acceptor_knowni[spliceinfo->segmentj_acceptor_nknown++] = j;
	  } else if (splicetypes[j] == ANTIDONOR) {
	    debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	    spliceinfo->segmentj_antidonor_knownpos[spliceinfo->segmentj_antidonor_nknown] = splicesites[j] - high_left;
	    spliceinfo->segmentj_antidonor_knowni[spliceinfo->segmentj_antidonor_nknown++] = j;
	  }
	  j++;
	}
      }
      spliceinfo->segmentj_acceptor_knownpos[spliceinfo->segmentj_acceptor_nknown] = querylength + 100;
      spliceinfo->segmentj_antidonor_knownpos[spliceinfo->segmentj_antidonor_nknown] = querylength + 100;
      
      if (sensedir == SENSE_FORWARD) {
	/* Sense */
	if ((splice_qpos = Splice_resolve_sense(&nindels,&indel_pos,&best_knowni_i,&best_knowni_j,
						&best_nmismatches_i,&best_nmismatches_j,&best_nmismatches_indel,
						&best_prob_i,&best_prob_j,
						/*segmenti_left*/left,/*segmentj_left*/high_left,chroffset,chroffset,
						qstart,high_qend,querylength,query_compress,
						spliceinfo,max_mismatches_allowed,plusp,genestrand,
						max_splice_deletionlen,max_splice_insertionlen,
						/*allow_indel_p*/true)) <= 0) {
	  debug13(printf("Splice_resolve_sense: fails\n"));

	} else if (nindels != 0 && indel_pos < splice_qpos) { /* indel is medial, splice is distal */
	  /* Push indel (based on left) then splice (verified) */
	  splice_distance = high_left - left + nindels;
      
	  debug13(printf("Splice_resolve_sense: %d indels at %d then splice at %d\n",nindels,indel_pos,splice_qpos));
	  newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);
	  newpath->last_medial_splice_prob = best_prob_i;
	  newpath->last_distal_splice_prob = best_prob_j;
	  newpath->last_distal_knowni = best_knowni_j;
	  
	  Intlist_head_set(newpath->endpoints,indel_pos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend);
	  
	  /* Indel first */
	  if (nindels > 0) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	  } else {
	    deletionpos = left + indel_pos;
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(-nindels,deletionpos));
	  }
	  
	  /* Splice second */
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  }
	  
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_indel);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,left - nindels); /* ? check ! */
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);

	} else if (nindels != 0 && splice_qpos < indel_pos) {
	  /* Push splice then indel (based on high_left) (verified) */
	  splice_distance = high_left - left + nindels;

	  debug13(printf("Splice_resolve_sense: splice at %d then %d indels at %d\n",splice_qpos,nindels,indel_pos));
	  newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);
	  Intlist_head_set(newpath->endpoints,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,indel_pos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend);
	  
	  /* Splice first */
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  }
	  
	  /* Indel second */
	  if (nindels > 0) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	  } else {
	    deletionpos = left + indel_pos;
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(-nindels,deletionpos));
	  }
	  
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_i);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_indel);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left + nindels); /* ? check ! */
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);

	} else {
	  /* Splice only */
	  splice_distance = high_left - left;

	  newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);
	  newpath->last_medial_splice_prob = best_prob_i;
	  newpath->last_distal_splice_prob = best_prob_j;
	  newpath->last_distal_knowni = best_knowni_j;
	  
	  Intlist_head_set(newpath->endpoints,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend);
	  
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_FORWARD,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_sense: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  }
	  
	  /* For qend, push i first, then push j */
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_i);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
	}
      
      } else {
	/* Antisense */
	if ((splice_qpos = Splice_resolve_antisense(&nindels,&indel_pos,&best_knowni_i,&best_knowni_j,
						    &best_nmismatches_i,&best_nmismatches_j,&best_nmismatches_indel,
						    &best_prob_i,&best_prob_j,
						    /*segmenti_left*/left,/*segmentj_left*/high_left,chroffset,chroffset,
						    qstart,high_qend,querylength,query_compress,
						    spliceinfo,max_mismatches_allowed,plusp,genestrand,
						    max_splice_deletionlen,max_splice_insertionlen,
						    /*allow_indel_p*/true)) <= 0) {
	  debug13(printf("Splice_resolve_antisense: fails\n"));
	
	} else if (nindels != 0 && indel_pos < splice_qpos) { /* indel is medial, splice is distal */
	  /* Push indel (based on left) then splice (verified) */
	  splice_distance = high_left - left + nindels;

	  debug13(printf("Splice_resolve_antisense: %d indels at %d then splice at %d\n",nindels,indel_pos,splice_qpos));
	  newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);
	  newpath->last_medial_splice_prob = best_prob_i;
	  newpath->last_distal_splice_prob = best_prob_j;
	  newpath->last_distal_knowni = best_knowni_j;
	  
	  Intlist_head_set(newpath->endpoints,indel_pos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend);
	  
	  /* Indel first */
	  if (nindels > 0) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	  } else {
	    deletionpos = left + indel_pos;
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(-nindels,deletionpos));
	  }
	  
	  /* Splice second */
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  }
	  
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_indel);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_i);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,left - nindels); /* ? check ! */
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);

	} else if (nindels != 0 && splice_qpos < indel_pos) {
	  /* Push splice then indel (based on high_left) */
	  splice_distance = high_left - left + nindels;  /* verified */

	  debug13(printf("Splice_resolve_antisense: splice at %d then %d indels at %d\n",splice_qpos,nindels,indel_pos));
	  newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);
	  Intlist_head_set(newpath->endpoints,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,indel_pos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend);
	  
	  /* Splice first */
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  }
	  
	  /* Indel second */
	  if (nindels > 0) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_insertion(nindels));
	  } else {
	    deletionpos = left + indel_pos;
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,Junction_new_deletion(-nindels,deletionpos));
	  }
	  
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_i);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_indel);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left + nindels); /* ? check ! */
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);

	} else {
	  /* Splice only */
	  splice_distance = high_left - left;

	  newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,intlistpool,univcoordlistpool,listpool);
	  newpath->last_medial_splice_prob = best_prob_i;
	  newpath->last_distal_splice_prob = best_prob_j;
	  newpath->last_distal_knowni = best_knowni_j;
	  
	  Intlist_head_set(newpath->endpoints,splice_qpos);
	  newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend);
	  
	  if (plusp == true) {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_j,/*acceptor_prob*/best_prob_i));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_j,best_prob_i));
	  } else {
	    newpath->junctions = Listpool_push(newpath->junctions,listpool,
					       Junction_new_splice(splice_distance,SENSE_ANTI,/*donor_prob*/best_prob_i,/*acceptor_prob*/best_prob_j));
	    debug13(printf("Splice_resolve_antisense: splice_qpos in range %d..%d is at %d with distance %u, nmismatches %d+%d, and probs %f and %f\n",
			   qstart,high_qend,splice_qpos,splice_distance,best_nmismatches_i,best_nmismatches_j,best_prob_i,best_prob_j));
	  }
	  
	  /* For qend, push i first, then push j */
	  Intlist_head_set(newpath->nmismatches,best_nmismatches_i);
	  newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,best_nmismatches_j);
	  
	  newpath->lefts = Univcoordlistpool_push(newpath->lefts,univcoordlistpool,high_left);
	  newpaths = Listpool_push(newpaths,listpool,(void *) newpath);
	}
      }
    }
  }

  if (splice3p == true) {
    FREE(ambig_qends);
    FREE(ambig_probs_3);
  }

  return newpaths;
}


static Substring_T
combine_into_qend_alts (List_T newpaths, int common_splice_qpos, double common_splice_prob,
			int querylength, bool plusp, int genestrand, Compress_T query_compress,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			int sensedir) {
  T newpath;
  List_T p;
  Univcoord_T *alts_coords;
  int *alts_nmismatches, *alts_knowni;
  double *alts_probs;
  int n, k;

  n = List_length(newpaths);
  alts_coords = (Univcoord_T *) MALLOC_OUT(n*sizeof(Univcoord_T));
  alts_knowni = (int *) MALLOC_OUT(n*sizeof(int));
  alts_nmismatches = (int *) MALLOC_OUT(n*sizeof(int));
  alts_probs = (double *) MALLOC_OUT(n*sizeof(double));

  for (p = newpaths, k = 0; p != NULL; p = List_next(p), k++) {
    newpath = List_head(p);
    alts_coords[k] = Univcoordlist_head(newpath->lefts) + common_splice_qpos;
    alts_nmismatches[k] = Intlist_head(newpath->nmismatches);
    alts_knowni[k] = newpath->last_distal_knowni;
    alts_probs[k] = newpath->last_distal_splice_prob;
    debug13(printf("End alt %u [%u\\n",alts_coords[k],alts_coords[k] - chroffset));
  }

  /* qend = Intlist_head(newpath->endpoints); */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      return Substring_new_alts_A(/*querystart*/common_splice_qpos,/*queryend*/querylength,
				  common_splice_qpos,querylength,plusp,genestrand,query_compress,
				  chrnum,chroffset,chrhigh,chrlength,
				  alts_coords,alts_knowni,alts_nmismatches,alts_probs,n,
				  common_splice_prob,/*substring1p*/false);
    } else {
      return Substring_new_alts_D(/*querystart*/0,/*queryend*/querylength - common_splice_qpos,
				  common_splice_qpos,querylength,plusp,genestrand,query_compress,
				  chrnum,chroffset,chrhigh,chrlength,
				  alts_coords,alts_knowni,alts_nmismatches,alts_probs,n,
				  common_splice_prob,/*substring1p*/true);
    }

  } else {
    if (plusp == true) {
      return Substring_new_alts_D(/*querystart*/common_splice_qpos,/*queryend*/querylength,
				  common_splice_qpos,querylength,plusp,genestrand,query_compress,
				  chrnum,chroffset,chrhigh,chrlength,
				  alts_coords,alts_knowni,alts_nmismatches,alts_probs,n,
				  common_splice_prob,/*substring1p*/false);
    } else {
      return Substring_new_alts_A(/*querystart*/0,/*queryend*/querylength - common_splice_qpos,
				  common_splice_qpos,querylength,plusp,genestrand,query_compress,
				  chrnum,chroffset,chrhigh,chrlength,
				  alts_coords,alts_knowni,alts_nmismatches,alts_probs,n,
				  common_splice_prob,/*substring1p*/true);
    }
  }
}


/* Returns a list of Path_T objects */
static List_T
compute_qend_paths (int depth, bool *terminalp, T path, List_T diagonals,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		    int querylength, Spliceinfo_T spliceinfo, int *mismatch_positions_alloc, Compress_T query_compress,
		    bool plusp, int genestrand, int max_mismatches_allowed,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Univdiagpool_T univdiagpool, int sensedir, bool innerp) {
  List_T all_child_paths, child_paths, newpaths, newdiagonals, p;
  T best_child_path, child_path;
  bool all_terminal_p, child_terminal_p;
  int common_splice_qpos, qend, ignore_int;
  Univcoord_T ignore_univcoord;
  Junction_T junction;
  double common_splice_prob;
  int common_nmismatches;
  Univdiag_T diagonal;

#ifdef DEBUG13
  printf("Entered compute_qend_paths with innerp %d.  Current path: ",innerp);
  Path_print(path,chroffset);
  printf("Diagonals:\n");
  for (p = diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(p);
    printf("%u %d..%d\n",diagonal->univdiagonal - chroffset,diagonal->qstart,diagonal->qend);
  }
  printf("\n");
#endif

  if (diagonals == NULL) {
    debug13(printf("diagonals is NULL, so exiting.  terminalp is true\n"));
    *terminalp = true;
    return Listpool_push(NULL,listpool,(void *) path);

  } else if (depth > MAX_DEPTH_MIDDLE) {
    debug13(printf("reached maximum depth, so exiting\n"));
    *terminalp = false;
    return Listpool_push(NULL,listpool,(void *) path);

  } else {
    /* Partition diagonals into those that can be added, and the rest */
    newpaths = (List_T) NULL;
    newdiagonals = (List_T) NULL;

    qend = Intlist_head(path->endpoints);
    for (p = diagonals; p != NULL; p = List_next(p)) {
      diagonal = (Univdiag_T) List_head(p);
      if (diagonal->qend <= qend) {
	/* Skip */
      } else if ((child_paths = attach_qend_diagonal(path,/*high_left*/diagonal->univdiagonal,
						     /*high_qstart*/diagonal->qstart,
						     /*high_qend*/diagonal->qend,
						     chrnum,chroffset,querylength,
						     spliceinfo,mismatch_positions_alloc,query_compress,
						     plusp,genestrand,max_mismatches_allowed,intlistpool,
						     univcoordlistpool,listpool,sensedir)) != NULL) {
	/* *extendedp = true; */
	newpaths = List_append(child_paths,newpaths);
      } else {
	/* Save for the next iteration */
	newdiagonals = Univdiagpool_push_existing(newdiagonals,univdiagpool,diagonal);
      }
    }

    if (newpaths == NULL) {
      debug13(printf("newpaths is NULL, so exiting.  terminalp is true\n"));
      *terminalp = true;
      return Listpool_push(NULL,listpool,(void *) path);

    } else if (List_next(newpaths) == NULL) {
      debug13(printf("newpaths has one entry, so exiting.  terminalp is false\n"));
      *terminalp = false;
      /* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
      Path_free(&path);
      return newpaths;

    } else {
      all_terminal_p = true;
      all_child_paths = (List_T) NULL;
      for (p = newpaths; p != NULL; p = List_next(p)) {
	debug13(printf("Calling recursively on "));
	debug13(Path_print(List_head(p),chroffset));
	debug13(printf("\n"));
	child_paths = compute_qend_paths(depth+1,&child_terminal_p,List_head(p),newdiagonals,
					 chrnum,chroffset,chrhigh,chrlength,querylength,
					 spliceinfo,mismatch_positions_alloc,query_compress,plusp,genestrand,
					 max_mismatches_allowed,intlistpool,univcoordlistpool,
					 listpool,univdiagpool,sensedir,innerp);
	debug13(printf("Returning recursively\n"));
	all_child_paths = List_append(all_child_paths,child_paths);
	if (child_terminal_p == false) {
	  all_terminal_p = false;
	}
      }

      if (all_terminal_p == false) {
	debug13(printf("Not all children are terminals, so exiting\n"));
	/* List_free(&newpaths); -- allocated by Listpool_push */
	/* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
	*terminalp = false;
	Path_free(&path);
	return all_child_paths;

      } else if (innerp == false &&
#ifdef DEBUG13
		 (best_child_path = find_best_splice(all_child_paths,chroffset))
#else
		 (best_child_path = find_best_splice(all_child_paths))
#endif
		 != NULL) {
	debug13(printf("Best child path:\n"));
	debug13(Path_print(best_child_path,chroffset));
	for (p = all_child_paths; p != NULL; p = List_next(p)) {
	  child_path = (T) List_head(p);
	  if (child_path != best_child_path) {
	    Path_free(&child_path);
	  }
	}
	*terminalp = false;
	Path_free(&path);
	return Listpool_push(NULL,listpool,(void *) best_child_path);

      } else if ((common_splice_qpos = find_common_splice_qpos(&common_splice_prob,&common_nmismatches,all_child_paths)) < 0) {
	debug13(printf("No common splice pos found, so exiting\n"));
	/* List_free(&newpaths); -- allocated by Listpool_push */
	/* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
	*terminalp = false;
	Path_free(&path);
	return all_child_paths;

      } else if (common_splice_prob < 0.9) {
	debug13(printf("Common splice prob %f is poor, so exiting\n",common_splice_prob));
	/* List_free(&newpaths); -- allocated by Listpool_push */
	/* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
	*terminalp = false;
	Path_free(&path);
	return all_child_paths;

#if 0
      } else if ((best_child_path = dominating_child_path(all_child_paths)) != NULL) {
	debug13(printf("Found a single best path\n"));
	for (p = all_child_paths; p != NULL; p = List_next(p)) {
	  if ((child_path = List_head(p)) != best_child_path) {
	    Path_free(&child_path);
	  }
	}
	*terminalp = false;
	Path_free(&path);
	return Listpool_push(NULL,listpool,(void *) best_child_path);
#endif

      } else {
	debug13(printf("All %d children are terminals with a common splice pos of %d, so combining into alt\n",
		       List_length(all_child_paths),common_splice_qpos));
#ifdef DEBUG13
	best_child_path = find_best_splice(all_child_paths,chroffset);
#else
	best_child_path = find_best_splice(all_child_paths);
#endif
	assert(best_child_path != NULL);
	debug13(printf("Best child path:\n"));
	debug13(Path_print(best_child_path,chroffset));

	best_child_path->alts_substring = combine_into_qend_alts(all_child_paths,common_splice_qpos,common_splice_prob,
								 querylength,plusp,genestrand,query_compress,
								 chrnum,chroffset,chrhigh,chrlength,sensedir);

	best_child_path->endpoints = Intlistpool_pop(best_child_path->endpoints,&ignore_int);
	best_child_path->nmismatches = Intlistpool_pop(best_child_path->nmismatches,&ignore_int);
	best_child_path->lefts = Univcoordlistpool_pop(best_child_path->lefts,&ignore_univcoord);
	best_child_path->junctions = Listpool_pop(best_child_path->junctions,(void **) &junction);
	Junction_free(&junction);

	debug13(printf("(3) Path after alt: "));
	debug13(Path_print(best_child_path,chroffset));

	for (p = all_child_paths; p != NULL; p = List_next(p)) {
	  child_path = (T) List_head(p);
	  if (child_path != best_child_path) {
	    Path_free(&child_path);
	  }
	}

	*terminalp = false;
	Path_free(&path);

	/* List_free(&newpaths); -- allocated by Listpool_push */
	/* List_free(&newdiagonals); -- allocated by Univdiagpool_push */
	return Listpool_push(NULL,listpool,(void *) best_child_path);
      }
    }
  }
}


static List_T
add_qend_local (bool *addedp, T path, Univcoord_T *local_diagonals, int ndiagonals, int pos5,
		Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		int querylength, Spliceinfo_T spliceinfo, int *mismatch_positions_alloc,
		Compress_T query_compress, bool plusp, int genestrand, int max_mismatches_allowed,
		Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		Listpool_T listpool, int sensedir, bool innerp) {
  List_T newpaths = NULL, paths, p;
  T best_child_path, child_path;
  int common_splice_qpos;
  int common_nmismatches;
  double common_splice_prob;
  int nhits, i;

  void *ignore;
  int ignore_int;
  Univcoord_T ignore_univcoord;
  Junction_T junction;


  debug13(printf("Entered add_qend_local.  Current path: "));
  debug13(Path_print(path,chroffset));
  debug13(printf("\n"));
#ifdef CHECK_ASSERTIONS
  check_for_ascending_values(local_diagonals,ndiagonals);
#endif

  *addedp = false;

  i = 0;
  while (i < ndiagonals /*&& nhits < 2*/) {
    debug13(printf("attaching qend local %u %d..%d\n",local_diagonals[i] - chroffset,pos5,querylength));
    if ((paths = attach_qend_diagonal(path,/*high_left*/local_diagonals[i],
				      /*high_qstart*/pos5,/*high_qend*/querylength,
				      chrnum,chroffset,querylength,
				      spliceinfo,mismatch_positions_alloc,query_compress,
				      plusp,genestrand,max_mismatches_allowed,intlistpool,
				      univcoordlistpool,listpool,sensedir)) != NULL) {
      newpaths = List_append(paths,newpaths);
      *addedp = true;
    }
    i++;
  }

  if ((nhits = List_length(newpaths)) == 0) {
    debug13(printf("newpaths is NULL, so exiting\n"));
    return Listpool_push(NULL,listpool,(void *) path);

  } else if (nhits == 1) {
    debug13(printf("newpaths has one entry, so exiting\n"));
    Path_free(&path);
    return newpaths;

  } else if (innerp == false &&
#ifdef DEBUG13
	     (best_child_path = find_best_splice(newpaths,chroffset))
#else
	     (best_child_path = find_best_splice(newpaths))
#endif
	     != NULL) {
    debug13(printf("Found a best child path\n"));
    for (p = newpaths; p != NULL; p = List_next(p)) {
      child_path = (T) List_head(p);
      if (child_path != best_child_path) {
	Path_free(&child_path);
      }
    }
    Path_free(&path);
    return Listpool_push(NULL,listpool,(void *) best_child_path);

  } else if ((common_splice_qpos = find_common_splice_qpos(&common_splice_prob,&common_nmismatches,newpaths)) < 0) {
    debug13(printf("No common splice pos found, so exiting\n"));
    Path_free(&path);
    return newpaths;

  } else if (common_splice_prob < 0.9) {
    debug13(printf("Common splice prob %f is poor, so exiting\n",common_splice_prob));
    Path_free(&path);
    return newpaths;

  } else {
    debug13(printf("Creating an alt substring\n"));
#if 0
    Intlist_head_set(path->endpoints,common_splice_qpos);
    Intlist_head_set(path->nmismatches,common_nmismatches);
#else
    Path_free(&path);
    path = (Path_T) List_head(newpaths);
    path->endpoints = Intlistpool_pop(path->endpoints,&ignore_int);
    path->nmismatches = Intlistpool_pop(path->nmismatches,&ignore_int);
    path->lefts = Univcoordlistpool_pop(path->lefts,&ignore_univcoord);
    if (path->junctions != NULL) {
      path->junctions = Listpool_pop(path->junctions,(void **) &junction);
      Junction_free(&junction);
    }
#endif

    path->alts_substring = combine_into_qend_alts(newpaths,common_splice_qpos,common_splice_prob,
						  querylength,plusp,genestrand,query_compress,
						  chrnum,chroffset,chrhigh,chrlength,sensedir);
    debug13(printf("(4) Resulting path: "));
    debug13(Path_print(path,chroffset));
    newpaths = Listpool_pop(newpaths,(void **) &ignore);
    Path_gc(&newpaths);
    return Listpool_push(NULL,listpool,(void *) path);
  }
}


/* Sometimes merging of left and right paths can result in anomalies */
static bool
endpoints_acceptable_p (Intlist_T endpoints, List_T junctions) {
  Intlist_T p;
  List_T q;
  Junction_T junction;
  int last_endpoint;

  debug13(printf("Evaluating endpoints for acceptability: %s\n",Intlist_to_string(endpoints)));

  last_endpoint = 0;
  /* Skip first endpoint (which we force to be 0) and last endpoint (which we force to be querylength) */
  for (p = Intlist_next(endpoints), q = junctions; Intlist_next(p) != NULL; p = Intlist_next(p), q = List_next(q)) {
    junction = (Junction_T) List_head(q);
    if (last_endpoint + Junction_ninserts(junction) >= Intlist_head(p)) {
      debug13(printf("Endpoint %d + %d >= %d, so unacceptable\n",
		     last_endpoint,Junction_ninserts(junction),Intlist_head(p)));
      return false;
    } else {
      debug13(printf("Endpoint %d + %d < %d, so acceptable\n",
		     last_endpoint,Junction_ninserts(junction),Intlist_head(p)));
    }
    last_endpoint = Intlist_head(p);
  }

  return true;
}


/* Always solves against plus strand of genome.  Just provide either
   queryuc/query_compress_fwd (coords measured from beginning of
   sequence) or queryrc/query_compress_rev (coords measured from end
   of sequence).  All coordinates measured from low end.
   Sense/antisense is with respect to the plus strand.  But to
   interface with Stage3end_new_substring command, need to flip
   coordinates for case where queryrc aligns to plus strand. */

/* chrnum is fixed from middle_diagonal */
static List_T
solve_via_segments_genome (bool *foundp, int *found_score_overall, int *found_score_within_trims,
			   List_T hits, List_T qstart_paths, List_T qend_paths, int sensedir,
			   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			   int querylength, Compress_T query_compress, bool plusp, int genestrand,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Hitlistpool_T hitlistpool, Method_T method, int level) {
  Stage3end_T hit;
  int best_nsegments, nsegments;
  int best_nmatches, nmatches, middle_nmismatches;
  List_T a;

  Intlist_T endpoints, q;
  Univcoordlist_T lefts, r;
  Intlist_T nmismatches, s;
  Junction_T junction;
  List_T junctions, t;

  T best_qstart_path, best_qend_path, qstart_path, qend_path;
  int middle_qstart_r, middle_qend_r, middle_qstart_l, middle_qend_l;
  bool left_chop_p, right_chop_p;
  Splicetype_T left_splicetype, right_splicetype;
  double left_ambig_prob, right_ambig_prob;
  int int_ignore;
  Univcoord_T univcoord_ignore;
  void *void_ignore;


#ifdef DEBUG13
  printf("\n");
  printf("*** Entered solve_via_segments_genome: sensedir %d\n",sensedir);
#endif


  /* Test left paths separately */
  best_qstart_path = (T) NULL;
  best_nsegments = 0;
  best_nmatches = 0;
  for (a = qstart_paths; a != NULL; a = List_next(a)) {
    qstart_path = (T) List_head(a);
    debug13(printf("\n"));
    debug13(printf("+ Left path %p: ",qstart_path));
    debug13(Path_print(qstart_path,chroffset));

    /* No need to copy */
    q = qstart_path->endpoints;
    r = qstart_path->lefts;
    s = qstart_path->nmismatches;
    t = qstart_path->junctions;

    if (endpoints_acceptable_p(q,t) == false) {
      debug13(printf("=> Unacceptable\n"));

    } else if ((nsegments = Path_nsegments(qstart_path)) > best_nsegments) {
      best_qstart_path = qstart_path;
      best_nsegments = nsegments;
      best_nmatches = Stage3end_nmatches_substrings(/*endpoints*/q,/*lefts*/r,/*nmismatches*/s,/*junctions*/t,
						    querylength,query_compress,
						    /*qend_alts*/NULL,/*qstart_alts*/qstart_path->alts_substring,
						    plusp,genestrand,chrnum,chroffset,chrhigh,chrlength,
						    qstart_path->splice5p,qstart_path->splice3p,
						    listpool);
      debug13(printf("(A) Left path without extension: %d segments, %d matches\n",nsegments,best_nmatches));
    } else if (nsegments == best_nsegments &&
	       (nmatches = Stage3end_nmatches_substrings(/*endpoints*/q,/*lefts*/r,/*nmismatches*/s,/*junctions*/t,
							 querylength,query_compress,
							 /*qend_alts*/NULL,/*qstart_alts*/qstart_path->alts_substring,
							 plusp,genestrand,chrnum,chroffset,chrhigh,chrlength,
							 qstart_path->splice5p,qstart_path->splice3p,
							 listpool)) > best_nmatches) {
      best_qstart_path = qstart_path;
      best_nmatches = nmatches;
      debug13(printf("(A) Left path without extension: %d segments, %d matches\n",nsegments,nmatches));
    } else {
      debug13(printf("(A) Left path without extension: %d segments, %d matches\n",nsegments,nmatches));
    }
  }

  if (best_qstart_path == NULL) {
    /* Skip */
  } else if (best_qstart_path->alts_substring != NULL) {
    left_chop_p = false;
  } else if (best_qstart_path->junctions != NULL &&
      Junction_type(junction = (Junction_T) List_head(best_qstart_path->junctions)) == SPLICE_JUNCTION &&
      Junction_splice_score(junction) < 1.8) {
    left_chop_p = true;
    if (plusp == true) {
      if (sensedir == SENSE_FORWARD) {
	left_splicetype = ACCEPTOR;
	left_ambig_prob = Junction_acceptor_prob(junction);
      } else if (sensedir == SENSE_ANTI) {
	left_splicetype = ANTIDONOR;
	left_ambig_prob = Junction_donor_prob(junction);
      }
    } else {
      if (sensedir == SENSE_FORWARD) {
	left_splicetype = DONOR;
	left_ambig_prob = Junction_donor_prob(junction);
      } else if (sensedir == SENSE_ANTI) {
	left_splicetype = ANTIACCEPTOR;
	left_ambig_prob = Junction_acceptor_prob(junction);
      }
    }

  } else {
    left_chop_p = false;
  }


  /* Test right paths separately */
  best_qend_path = (T) NULL;
  best_nsegments = 0;
  best_nmatches = 0;
  for (a = qend_paths; a != NULL; a = List_next(a)) {
    qend_path = (T) List_head(a);
    debug13(printf("\n"));
    debug13(printf("+ Right path %p: ",qend_path));
    debug13(Path_print(qend_path,chroffset));

    q = Intlist_reverse(Intlistpool_copy(qend_path->endpoints,intlistpool));
    r = Univcoordlist_reverse(Univcoordlistpool_copy(qend_path->lefts,univcoordlistpool));
    s = Intlist_reverse(Intlistpool_copy(qend_path->nmismatches,intlistpool));
    t = List_reverse(Listpool_copy(qend_path->junctions,listpool));

    if (endpoints_acceptable_p(q,t) == false) {
      debug13(printf("=> Unacceptable\n"));
    } else if ((nsegments = Path_nsegments(qend_path)) > best_nsegments) {
      best_qend_path = qend_path;
      best_nsegments = nsegments;
      best_nmatches = Stage3end_nmatches_substrings(/*endpoints*/q,/*lefts*/r,/*nmismatches*/s,/*junctions*/t,
						    querylength,query_compress,
						    /*qend_alts*/qend_path->alts_substring,/*qstart_alts*/NULL,
						    plusp,genestrand,chrnum,chroffset,chrhigh,chrlength,
						    qend_path->splice5p,qend_path->splice3p,
						    listpool);
      debug13(printf("(A) Right path without extension: %d segments, %d matches\n",nsegments,best_nmatches));
    } else if (nsegments == best_nsegments &&
	       (nmatches = Stage3end_nmatches_substrings(/*endpoints*/q,/*lefts*/r,/*nmismatches*/s,/*junctions*/t,
							 querylength,query_compress,
							 /*qend_alts*/qend_path->alts_substring,/*qstart_alts*/NULL,
							 plusp,genestrand,chrnum,chroffset,chrhigh,chrlength,
							 qend_path->splice5p,qend_path->splice3p,
							 listpool)) > best_nmatches) {
      best_qend_path = qend_path;
      best_nmatches = nmatches;
      debug13(printf("(A) Right path without extension: %d segments, %d matches\n",nsegments,nmatches));
    } else {
      debug13(printf("(A) Right path without extension: %d segments, %d matches\n",nsegments,nmatches));
    }

    /* Intlist_free(&q); -- allocated by Intlistpool_push */
    /* Univcoordlist_free(&r); -- allocated by Uintlistpool_push */
    /* Intlist_free(&s); -- allocated by Intlistpool_push */
    /* List_free(&t); -- Allocated by Listpool_push */
  }

  if (best_qend_path == NULL) {
    /* Skip */
  } else if (best_qend_path->alts_substring != NULL) {
    right_chop_p = false;
  } else if (best_qend_path->junctions != NULL &&
	     Junction_type(junction = (Junction_T) List_head(best_qend_path->junctions)) == SPLICE_JUNCTION &&
      Junction_splice_score(junction) < 1.8) {
    right_chop_p = true;
    if (plusp == true) {
      if (sensedir == SENSE_FORWARD) {
	right_splicetype = DONOR;
	right_ambig_prob = Junction_donor_prob(junction);
      } else if (sensedir == SENSE_ANTI) {
	right_splicetype = ANTIACCEPTOR;
	right_ambig_prob = Junction_acceptor_prob(junction);
      }
    } else {
      if (sensedir == SENSE_FORWARD) {
	right_splicetype = ACCEPTOR;
	right_ambig_prob = Junction_acceptor_prob(junction);
      } else if (sensedir == SENSE_ANTI) {
	right_splicetype = ANTIDONOR;
	right_ambig_prob = Junction_donor_prob(junction);
      }
    }
  } else {
    right_chop_p = false;
  }


  /* Combine best_qstart_path and best_qend_path */
  if (best_qstart_path == NULL || best_qend_path == NULL) {
    /* All paths must have gone over chromosome boundary */
    return hits;
  } else {
    debug13(printf("++ Best left path %p: ",best_qstart_path));
    debug13(Path_print(best_qstart_path,chroffset));
    debug13(printf("\n"));
    debug13(printf("++ Best right path %p: ",best_qend_path));
    debug13(Path_print(best_qend_path,chroffset));
    debug13(printf("\n\n"));

    middle_qend_l = Intlist_last_value(best_qstart_path->endpoints);
    middle_qstart_l = Intlist_penultimate_value(best_qstart_path->endpoints);
    /* This is the junction at the start of the middle segment */
    if (best_qstart_path->junctions != NULL &&
	(junction = (Junction_T) List_last_value(best_qstart_path->junctions)) != NULL) {
      middle_qstart_l += Junction_ninserts(junction);
    }
  
    middle_qstart_r = Intlist_last_value(best_qend_path->endpoints);
    middle_qend_r = Intlist_penultimate_value(best_qend_path->endpoints);
  }

  if (middle_qstart_l >= middle_qend_r) {
    /* Middle segment got eliminated, so abort */
    return hits;
  } else if (middle_qstart_l == middle_qstart_r) {
    /* Take nmismatches from qend_path */
    middle_nmismatches = Intlist_last_value(best_qend_path->nmismatches);
  } else if (middle_qend_l == middle_qend_r) {
    /* Take nmismatches from qstart_path */
    middle_nmismatches = Intlist_last_value(best_qstart_path->nmismatches);
  } else {
    /* Use -1 for the middle diagonal */
    middle_nmismatches = -1;
  }

  endpoints = Intlist_reverse(Intlistpool_copy_but_last(best_qend_path->endpoints,intlistpool));
  endpoints = Intlist_append(Intlistpool_copy_but_last(best_qstart_path->endpoints,intlistpool),endpoints);

  junctions = List_reverse(Listpool_copy(best_qend_path->junctions,listpool));
  junctions = List_append(Listpool_copy(best_qstart_path->junctions,listpool),junctions);
  
  if (endpoints_acceptable_p(endpoints,junctions) == false) {
    /* Intlist_free(&endpoints); -- allocated by Intlistpool_push */
    /* List_free(&junctions) -- allocated by Listpool_push */
    return hits;
  } else {
    nmismatches = Intlist_reverse(Intlistpool_copy_but_last(best_qend_path->nmismatches,intlistpool));
    nmismatches = Intlistpool_push(nmismatches,intlistpool,middle_nmismatches);
    nmismatches = Intlist_append(Intlistpool_copy_but_last(best_qstart_path->nmismatches,intlistpool),nmismatches);
  
    lefts = Univcoordlist_reverse(Univcoordlistpool_copy(best_qend_path->lefts,univcoordlistpool));
    lefts = Univcoordlist_append(Univcoordlistpool_copy_but_last(best_qstart_path->lefts,univcoordlistpool),lefts);
  }
  
  
#ifdef DEBUG13
  printf("\n");
  printf("Combined qstart_alts %p, qend_alts %p\n",best_qstart_path->alts_substring,best_qend_path->alts_substring);
  printf("Combined splice5p %d, splice3p %d\n",best_qstart_path->splice5p,best_qend_path->splice3p);
  printf("Combined ambig_prob_5 %f, ambig_prob_3 %f\n",best_qstart_path->ambig_prob_5,qend_path->ambig_prob_3);
  printf("Combined endpoints: %s\n",Intlist_to_string(endpoints));
  printf("Combined lefts: %s\n",Uintlist_to_string_offset(lefts,chroffset));
  printf("Combined mismatches: %s\n",Intlist_to_string(nmismatches));
  printf("Combined junctions: ");
  Junction_print_list(junctions);
  printf("\n");
#endif

  if ((hit = Stage3end_new_substrings(&(*found_score_overall),&(*found_score_within_trims),
				      endpoints,lefts,nmismatches,junctions,querylength,query_compress,
				      /*qend_alts*/best_qend_path->alts_substring,/*qstart_alts*/best_qstart_path->alts_substring,
				      plusp,genestrand,sensedir,chrnum,chroffset,chrhigh,chrlength,
				      best_qstart_path->splice5p,best_qstart_path->splicetype5,best_qstart_path->ambig_prob_5,
				      best_qend_path->splice3p,best_qend_path->splicetype3,best_qend_path->ambig_prob_3,
				      listpool,method,level)) != NULL) {
    *foundp = true;
    hits = Hitlist_push(hits,hitlistpool,(void *) hit);
  }

  /* Need to chop, because a bad splice can result in an overreach */
  if (left_chop_p == true && right_chop_p == true) {
    endpoints = Intlistpool_pop(Intlist_reverse(Intlistpool_pop(Intlist_reverse(endpoints),&int_ignore)),&int_ignore);
    lefts = Univcoordlistpool_pop(Univcoordlist_reverse(Univcoordlistpool_pop(Univcoordlist_reverse(lefts),&univcoord_ignore)),&univcoord_ignore);
    nmismatches = Intlistpool_pop(Intlist_reverse(Intlistpool_pop(Intlist_reverse(nmismatches),&int_ignore)),&int_ignore);
    junctions = Listpool_pop(List_reverse(Listpool_pop(List_reverse(junctions),&void_ignore)),&void_ignore);
#ifdef DEBUG13
    printf("\n");
    printf("Chopped probs: %f and %f\n",left_ambig_prob,right_ambig_prob);
    printf("Chopped endpoints: %s\n",Intlist_to_string(endpoints));
    printf("Chopped lefts: %s\n",Uintlist_to_string_offset(lefts,chroffset));
    printf("Chopped mismatches: %s\n",Intlist_to_string(nmismatches));
    printf("Chopped junctions: ");
    Junction_print_list(junctions);
    printf("\n");
#endif

    if ((hit = Stage3end_new_substrings(&(*found_score_overall),&(*found_score_within_trims),
					endpoints,lefts,nmismatches,junctions,
					querylength,query_compress,/*qend_alts*/NULL,/*qstart_alts*/NULL,
					plusp,genestrand,sensedir,chrnum,chroffset,chrhigh,chrlength,
					/*splice5p*/true,left_splicetype,left_ambig_prob,
					/*splice3p*/true,right_splicetype,right_ambig_prob,
					listpool,method,level)) != NULL) {
      *foundp = true;
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
    }

  } else if (left_chop_p == true) {
    endpoints = Intlistpool_pop(endpoints,&int_ignore);
    lefts = Univcoordlistpool_pop(lefts,&univcoord_ignore);
    nmismatches = Intlistpool_pop(nmismatches,&int_ignore);
    junctions = Listpool_pop(junctions,&void_ignore);
#ifdef DEBUG13
    printf("\n");
    printf("Chopped probs: %f and %f\n",left_ambig_prob,best_qend_path->ambig_prob_3);
    printf("Chopped endpoints: %s\n",Intlist_to_string(endpoints));
    printf("Chopped lefts: %s\n",Uintlist_to_string_offset(lefts,chroffset));
    printf("Chopped mismatches: %s\n",Intlist_to_string(nmismatches));
    printf("Chopped junctions: ");
    Junction_print_list(junctions);
    printf("\n");
#endif

    if ((hit = Stage3end_new_substrings(&(*found_score_overall),&(*found_score_within_trims),
					endpoints,lefts,nmismatches,junctions,
					querylength,query_compress,/*qend_alts*/best_qend_path->alts_substring,/*qstart_alts*/NULL,
					plusp,genestrand,sensedir,chrnum,chroffset,chrhigh,chrlength,
					/*splice5p*/true,left_splicetype,left_ambig_prob,
					best_qend_path->splice3p,best_qend_path->splicetype3,best_qend_path->ambig_prob_3,
					listpool,method,level)) != NULL) {
      *foundp = true;
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
    }

  } else if (right_chop_p == true) {
    endpoints = Intlist_reverse(Intlistpool_pop(Intlist_reverse(endpoints),&int_ignore));
    lefts = Univcoordlist_reverse(Univcoordlistpool_pop(Univcoordlist_reverse(lefts),&univcoord_ignore));
    nmismatches = Intlist_reverse(Intlistpool_pop(Intlist_reverse(nmismatches),&int_ignore));
    junctions = List_reverse(Listpool_pop(List_reverse(junctions),&void_ignore));
#ifdef DEBUG13
    printf("\n");
    printf("Chopped probs: %f and %f\n",best_qstart_path->ambig_prob_5,right_ambig_prob);
    printf("Chopped endpoints: %s\n",Intlist_to_string(endpoints));
    printf("Chopped lefts: %s\n",Uintlist_to_string_offset(lefts,chroffset));
    printf("Chopped mismatches: %s\n",Intlist_to_string(nmismatches));
    printf("Chopped junctions: ");
    Junction_print_list(junctions);
    printf("\n");
#endif

    if ((hit = Stage3end_new_substrings(&(*found_score_overall),&(*found_score_within_trims),
					endpoints,lefts,nmismatches,junctions,
					querylength,query_compress,/*qend_alts*/NULL,/*qstart_alts*/best_qstart_path->alts_substring,
					plusp,genestrand,sensedir,chrnum,chroffset,chrhigh,chrlength,
					best_qstart_path->splice5p,best_qstart_path->splicetype5,best_qstart_path->ambig_prob_5,
					/*splice3p*/true,right_splicetype,right_ambig_prob,
					listpool,method,level)) != NULL) {
      *foundp = true;
      hits = Hitlist_push(hits,hitlistpool,(void *) hit);
    }
  }

  /* Intlist_free(&endpoints); -- allocated by Intlistpool_push */
  /* Univcoordlist_free(&lefts); -- allocated by Uintlistpool_push */
  /* Intlist_free(&nmismatches); -- allocated by Intlistpool_push */
  /* List_free(&junctions) -- allocated by Listpool_push */

  return hits;
}


static void
check_for_descending_qend (int prev_qend, Univcoord_T prev_left, List_T qstart_diagonals) {
  List_T p;
  Univdiag_T diagonal;

  for (p = qstart_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(qstart_diagonals);
    if (diagonal->qend > prev_qend) {
      abort();
    } else if (diagonal->qend == prev_qend && diagonal->univdiagonal > prev_left) {
      abort();
    }
    prev_qend = diagonal->qend;
    prev_left = diagonal->univdiagonal;
  }
  return;
}

static void
check_for_ascending_qstart (int prev_qstart, Univcoord_T prev_left, List_T qend_diagonals) {
  List_T p;
  Univdiag_T diagonal;

  for (p = qend_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(qend_diagonals);
    if (diagonal->qstart < prev_qstart) {
      abort();
    } else if (diagonal->qstart == prev_qstart && diagonal->univdiagonal < prev_left) {
      abort();
    }
    prev_qstart = diagonal->qstart;
    prev_left = diagonal->univdiagonal;
  }
  return;
}
  

typedef struct Diagonals_store_T *Diagonals_store_T;
struct Diagonals_store_T {
  int ndiagonals;
  Univcoord_T *diagonals;
};

static Diagonals_store_T
Diagonals_store_new (Univcoord_T *diagonals, int ndiagonals) {
  Diagonals_store_T new = (Diagonals_store_T) MALLOC(sizeof(*new));

  new->ndiagonals = ndiagonals;
  new->diagonals = diagonals;
  return new;
}

static void
Diagonals_store_free (Diagonals_store_T *old) {
#ifndef LARGE_GENOMES
  if ((*old)->diagonals) {
    FREE_ALIGN((*old)->diagonals);	/* Aligned from Localdb_get_diagonals */
  }
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
  if ((*old)->diagonals) {
    FREE_ALIGN((*old)->diagonals);	/* Aligned from Localdb_get_diagonals */
  }
#else
  FREE((*old)->diagonals);
#endif
  FREE(*old);

  return;
}


static List_T
compute_qstart_local (int depth, Path_T path, char *queryptr, int querylength,
		      Localdb_T localdb1, Univcoordtable_T *local_diagonals_cache,
		      List_T *all_tables, int *mismatch_positions_alloc, Spliceinfo_T spliceinfo, 
		      Univcoord_T **stream_alloc, int *streamsize_alloc, Compress_T query_compress, 
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      int nmismatches_allowed, bool plusp, int genestrand, int sensedir, bool innerp) {
  List_T results, newpaths, p;
  bool addedp;
  Univcoordtable_T table;
  Univcoord_T *local_diagonals, left, low, high;
  Diagonals_store_T diagonals_store;
  int ndiagonals;
  int qstart;

  if (localdb1 == NULL) {
    /* Localdb not possible */
    return Listpool_push(NULL,listpool,(void *) path);

  } else if (path->alts_substring != NULL) {
    /* Not possible to add localdb because of the ambiguity */
    return Listpool_push(NULL,listpool,(void *) path);

  } else if ((qstart = Intlist_head(path->endpoints)) < endtrim_allowed) {
    /* Not worth it to add localdb because already near the end */
    return Listpool_push(NULL,listpool,(void *) path);

  } else if (depth > MAX_DEPTH_LOCAL) {
    /* Too much depth */
    return Listpool_push(NULL,listpool,(void *) path);

  } else {
    left = Univcoordlist_head(path->lefts);
    debug13(printf("(1) Calling localdb_get_diagonals at %d..%d at left %u\n",0,qstart,left - chroffset));

    if ((table = local_diagonals_cache[qstart]) == NULL) {
      table = local_diagonals_cache[qstart] = Univcoordtable_new(/*hint*/10,/*save_contents_p*/true);
      *all_tables = Listpool_push(*all_tables,listpool,(void *) table);

      if (circularp[chrnum] == true) {
	low = subtract_bounded(left + /*pos5*/0,max_deletionlen,chroffset);
      } else {
	low = subtract_bounded(left + /*pos5*/0,overall_end_distance_genome,chroffset);
      }
      high = add_bounded(left + /*pos3*/qstart,max_insertionlen,chrhigh);
      high = subtract_bounded(high,/*local1part*/8,chroffset);
      debug13(printf("  overall_end_distance %u, low %u, high %u\n",overall_end_distance_genome,low - chroffset,high - chroffset));
      local_diagonals = Localdb_get_diagonals(&ndiagonals,localdb1,queryptr,
					      /*pos5*/0,/*pos3*/qstart,low,high,
					      plusp,genestrand,stream_alloc,streamsize_alloc,
					      /*remove_repetitive_p*/false);
      Univcoordtable_put_and_save(table,left,(void *) Diagonals_store_new(local_diagonals,ndiagonals));
      
    } else if ((diagonals_store = Univcoordtable_get(table,left)) == NULL) {
      if (circularp[chrnum] == true) {
	low = subtract_bounded(left + /*pos5*/0,max_deletionlen,chroffset);
      } else {
	low = subtract_bounded(left + /*pos5*/0,overall_end_distance_genome,chroffset);
      }
      high = add_bounded(left + /*pos3*/qstart,max_insertionlen,chrhigh);
      high = subtract_bounded(high,/*local1part*/8,chroffset);
      debug13(printf("  overall_end_distance %u, low %u, high %u\n",overall_end_distance_genome,low - chroffset,high - chroffset));
      local_diagonals = Localdb_get_diagonals(&ndiagonals,localdb1,queryptr,
					      /*pos5*/0,/*pos3*/qstart,low,high,
					      plusp,genestrand,stream_alloc,streamsize_alloc,
					      /*remove_repetitive_p*/false);
      Univcoordtable_put_and_save(table,left,(void *) Diagonals_store_new(local_diagonals,ndiagonals));
      
    } else {
      local_diagonals = diagonals_store->diagonals;
      ndiagonals = diagonals_store->ndiagonals;
    } 
    debug13(printf("Got %d localdb diagonals\n",ndiagonals));
    
    if (ndiagonals == 0) {
      /* No extensions found */
      return Listpool_push(NULL,listpool,(void *) path);

    } else {
      newpaths = add_qstart_local(&addedp,path,local_diagonals,ndiagonals,/*pos3*/qstart,
				  chrnum,chroffset,chrhigh,chrlength,querylength,
				  spliceinfo,mismatch_positions_alloc,query_compress,
				  plusp,genestrand,nmismatches_allowed,intlistpool,univcoordlistpool,
				  listpool,sensedir,innerp);
      if (addedp == false) {
	/* No extensions could be added */
	return newpaths;
      } else {
	results = (List_T) NULL;
	for (p = newpaths; p != NULL; p = List_next(p)) {
	  results = List_append(results,
				compute_qstart_local(depth+1,(T) List_head(p),queryptr,querylength,
						     localdb1,local_diagonals_cache,
						     &(*all_tables),mismatch_positions_alloc,spliceinfo, 
						     stream_alloc,streamsize_alloc,query_compress, 
						     chrnum,chroffset,chrhigh,chrlength,intlistpool,
						     univcoordlistpool,listpool,nmismatches_allowed,
						     plusp,genestrand,sensedir,innerp));
	}
	return results;
      }
    }
  }
}


static List_T
compute_qend_local (int depth, Path_T path, char *queryptr, int querylength,
		    Localdb_T localdb1, Univcoordtable_T *local_diagonals_cache,
		    List_T *all_tables, int *mismatch_positions_alloc, Spliceinfo_T spliceinfo, 
		    Univcoord_T **stream_alloc, int *streamsize_alloc, Compress_T query_compress, 
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		    int nmismatches_allowed, bool plusp, int genestrand, int sensedir, bool innerp) {
  List_T results, newpaths, p;
  bool addedp;
  Univcoordtable_T table;
  Univcoord_T *local_diagonals, left, low, high;
  Diagonals_store_T diagonals_store;
  int ndiagonals;
  int qend;

  if (localdb1 == NULL) {
    /* Localdb not possible */
    return Listpool_push(NULL,listpool,(void *) path);
	
  } else if (path->alts_substring != NULL) {
    /* Not possible to add localdb because of the ambiguity */
    return Listpool_push(NULL,listpool,(void *) path);
	
  } else if (querylength - (qend = Intlist_head(path->endpoints)) < endtrim_allowed) {
    /* Not worth it to add localdb because already near the end */
    return Listpool_push(NULL,listpool,(void *) path);
	
  } else if (depth > MAX_DEPTH_LOCAL) {
    /* Too much depth */
    return Listpool_push(NULL,listpool,(void *) path);

  } else {
    left = Univcoordlist_head(path->lefts);
    debug13(printf("(2) Calling localdb_get_diagonals at %d..%d at left %u\n",qend,querylength,left - chroffset));
	
    if ((table = local_diagonals_cache[qend]) == NULL) {
      table = local_diagonals_cache[qend] = Univcoordtable_new(/*hint*/10,/*save_contents_p*/true);
      *all_tables = Listpool_push(*all_tables,listpool,(void *) table);
      
      low = subtract_bounded(left + /*pos5*/qend,max_insertionlen,chroffset);
      if (circularp[chrnum] == true) {
	high = add_bounded(left + /*pos3*/querylength,max_deletionlen,chrhigh);
      } else {
	high = add_bounded(left + /*pos3*/querylength,overall_end_distance_genome,chrhigh);
      }
      high = subtract_bounded(high,/*local1part*/8,chroffset);
      debug13(printf("  overall_end_distance %u, low %u, high %u\n",overall_end_distance_genome,low - chroffset,high - chroffset));
      local_diagonals = Localdb_get_diagonals(&ndiagonals,localdb1,queryptr,
					      /*pos5*/qend,/*pos3*/querylength,low,high,
					      plusp,genestrand,stream_alloc,streamsize_alloc,
					      /*remove_repetitive_p*/false);
      Univcoordtable_put_and_save(table,left,(void *) Diagonals_store_new(local_diagonals,ndiagonals));
      
    } else if ((diagonals_store = Univcoordtable_get(table,left)) == NULL) {
      low = subtract_bounded(left + /*pos5*/qend,max_insertionlen,chroffset);
      if (circularp[chrnum] == true) {
	high = add_bounded(left + /*pos3*/querylength,max_deletionlen,chrhigh);
      } else {
	high = add_bounded(left + /*pos3*/querylength,overall_end_distance_genome,chrhigh);
      }
      high = subtract_bounded(high,/*local1part*/8,chroffset);
      debug13(printf("  overall_end_distance %u, low %u, high %u\n",overall_end_distance_genome,low - chroffset,high - chroffset));
      local_diagonals = Localdb_get_diagonals(&ndiagonals,localdb1,queryptr,
					      /*pos5*/qend,/*pos3*/querylength,low,high,
					      plusp,genestrand,stream_alloc,streamsize_alloc,
					      /*remove_repetitive_p*/false);
      Univcoordtable_put_and_save(table,left,(void *) Diagonals_store_new(local_diagonals,ndiagonals));
      
    } else {
      local_diagonals = diagonals_store->diagonals;
      ndiagonals = diagonals_store->ndiagonals;
    } 
    debug13(printf("Got %d localdb diagonals\n",ndiagonals));
    
    if (ndiagonals == 0) {
      /* No extensions found */
      return Listpool_push(NULL,listpool,(void *) path);
      
    } else {
      newpaths = add_qend_local(&addedp,path,local_diagonals,ndiagonals,/*pos5*/qend,
				chrnum,chroffset,chrhigh,chrlength,querylength,
				spliceinfo,mismatch_positions_alloc,query_compress,
				plusp,genestrand,nmismatches_allowed,intlistpool,univcoordlistpool,
				listpool,sensedir,innerp);
      if (addedp == false) {
	/* No extensions could be added */
	return newpaths;
      } else {
	results = (List_T) NULL;
	for (p = newpaths; p != NULL; p = List_next(p)) {
	  results = List_append(results,
				compute_qend_local(depth+1,(T) List_head(p),queryptr,querylength,
						   localdb1,local_diagonals_cache,
						   &(*all_tables),mismatch_positions_alloc,spliceinfo, 
						   stream_alloc,streamsize_alloc,query_compress, 
						   chrnum,chroffset,chrhigh,chrlength,intlistpool,
						   univcoordlistpool,listpool,nmismatches_allowed,
						   plusp,genestrand,sensedir,innerp));
	}
	return results;
      }
    }
  }
}



List_T
Path_solve_from_diagonals (bool *foundp, int *found_score_overall, int *found_score_within_trims, List_T hits,
			   Univcoord_T middle_diagonal_univdiagonal, int middle_diagonal_qstart, int middle_diagonal_qend,
			   List_T qend_diagonals, List_T qstart_diagonals, char *queryptr, int querylength,
			   int *mismatch_positions_alloc, Spliceinfo_T spliceinfo, 
			   Univcoord_T **stream_alloc, int *streamsize_alloc, Compress_T query_compress, 
			   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			   bool plusp, int genestrand, int nmismatches_allowed, bool paired_end_p, bool first_read_p,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Univdiagpool_T univdiagpool, Hitlistpool_T hitlistpool,
			   Method_T method, int level) {
  List_T p;
  int gap, max_qend, min_qstart, qstart, qend, pos5, pos3;
  int *trim_qstarts, *trim_qends;
  Localdb_T localdb1;
  Univdiag_T diagonal;

  Univcoord_T min_low, max_high, low, high;
  int ndiagonals;

  List_T sense_qstart_paths, antisense_qstart_paths, sense_qend_paths, antisense_qend_paths, paths;
  Path_T path;
  bool innerp, terminalp;

  Univcoord_T *local_diagonals;
  
  Univcoordtable_T *local_diagonals_cache, table; /* local_diagonals are often re-used from sense to antisense */
  Diagonals_store_T *diagonals_stores;
  List_T all_tables;
  int nvalues, i;

  int result;
  int nspliceends;
  bool splice5p, splice3p;
  Splicetype_T splicetype5, splicetype3;
  int *ambig_qstarts, *ambig_qends;
  double *ambig_probs_5, *ambig_probs_3, zero_prob = 0.0;



  debug13(printf("Entered Path_solve_from_diagonals, first_read_p %d, with middle_diagonal #%d %u %d..%d, %d qstart diagonals, and %d qend diagonals\n",
		 first_read_p,chrnum,middle_diagonal_univdiagonal - chroffset,middle_diagonal_qstart,middle_diagonal_qend,
		 List_length(qstart_diagonals),List_length(qend_diagonals)));
#ifdef CHECK_ASSERTIONS
  check_for_descending_qend(middle_diagonal_qend,middle_diagonal_univdiagonal,qstart_diagonals);
  check_for_ascending_qstart(middle_diagonal_qstart,middle_diagonal_univdiagonal,qend_diagonals);
#endif
  localdb1 = (plusp == true) ? localdb : localdb2;


  if (qstart_diagonals != NULL && localdb1 != NULL) {
    /* Analyze qstart diagonals for gaps in the query */
    /* If NULL, then covered by the call to localdb starting at 0 */
    max_qend = 0;
    min_low = (Univcoord_T) -1;
    for (p = qstart_diagonals; p != NULL; p = List_next(p)) {
      diagonal = (Univdiag_T) List_head(p);
      if (diagonal->qend > max_qend) {
	max_qend = diagonal->qend;
      }
      if (diagonal->univdiagonal < min_low) {
	min_low = diagonal->univdiagonal;
      }
    }
    if ((gap = middle_diagonal_qstart - max_qend) > 8) {
      debug13(printf("Gap for qstart is %d = %d - %d\n",gap,middle_diagonal_qstart,max_qend));
      pos3 = middle_diagonal_qstart;
      low = min_low + /*pos5*/max_qend;
      high = middle_diagonal_univdiagonal + pos3 - /*local1part*/8;
      if (low < high) {
	local_diagonals = Localdb_get_diagonals(&ndiagonals,localdb1,queryptr,/*pos5*/max_qend,pos3,
						low,high,plusp,genestrand,stream_alloc,streamsize_alloc,
						/*remove_repetitive_p*/true);
	for (i = 0; i < ndiagonals; i++) {
	  qstart_diagonals = Univdiagpool_push(qstart_diagonals,univdiagpool,/*qstart*/max_qend,
					       /*qend*/pos3,local_diagonals[i]);
	}
#ifndef LARGE_GENOMES
	FREE_ALIGN(local_diagonals);
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
	FREE_ALIGN(local_diagonals);
#else
	FREE(local_diagonals);
#endif
      }
    }
  }


  if (qend_diagonals != NULL && localdb1 != NULL) {
    /* Analyze qend diagonals for gaps in the query */
    /* If NULL, then covered by the call to localdb ending at querylength */
    min_qstart = querylength;
    max_high = (Univcoord_T) 0;
    for (p = qend_diagonals; p != NULL; p = List_next(p)) {
      diagonal = (Univdiag_T) List_head(p);
      if (diagonal->qstart < min_qstart) {
	min_qstart = diagonal->qstart;
      }
      if (diagonal->univdiagonal > max_high) {
	max_high = diagonal->univdiagonal;
      }
    }
    if ((gap = min_qstart - middle_diagonal_qend) > 8) {
      debug13(printf("Gap for qend is %d = %d - %d\n",gap,min_qstart,middle_diagonal_qend));
      pos5 = middle_diagonal_qend;
      low = middle_diagonal_univdiagonal + pos5;
      high = max_high + /*pos3*/min_qstart - /*local1part*/8;
      if (low < high) {
	local_diagonals = Localdb_get_diagonals(&ndiagonals,localdb1,queryptr,pos5,/*pos3*/min_qstart,
						low,high,plusp,genestrand,stream_alloc,streamsize_alloc,
						/*remove_repetitive_p*/true);
	for (i = 0; i < ndiagonals; i++) {
	  qend_diagonals = Univdiagpool_push(qend_diagonals,univdiagpool,/*qstart*/pos5,
					     /*qend*/min_qstart,local_diagonals[i]);
	}
#ifndef LARGE_GENOMES
	FREE_ALIGN(local_diagonals);
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
	FREE_ALIGN(local_diagonals);
#else
	FREE(local_diagonals);
#endif
      }
    }
  }


  /* Qstart */
  local_diagonals_cache = (Univcoordtable_T *) CALLOC(querylength,sizeof(Univcoordtable_T));
  all_tables = (List_T) NULL;
  if (paired_end_p == false) {
    innerp = false;
  } else if (first_read_p == plusp) {
    innerp = false;
  } else {
    innerp = true;
  }
  

  /* Qstart sense.  Also handles the non-splicing case. */
  debug13(printf("Calling Substring_trimmed_qstarts with %d..%d\n",middle_diagonal_qstart,middle_diagonal_qend));
  splice5p = Substring_trimmed_qstarts(&result,&splicetype5,&ambig_qstarts,&ambig_probs_5,
				       /*left*/middle_diagonal_univdiagonal,
				       /*qstart:middle_diagonal_qstart,*/
				       /*qend*/middle_diagonal_qend,plusp,genestrand,
				       mismatch_positions_alloc,query_compress,
				       chroffset,/*sensedir*/SENSE_FORWARD);
  debug13(printf("splice5p %d, result %d\n",splice5p,result));
  
  if (result < 0) {
    nspliceends = 0;
    /* trim_qstarts = (int *) NULL; */
  } else if (splice5p == false) {
    nspliceends = 1;
    trim_qstarts = &result;
    ambig_probs_5 = &zero_prob;
  } else {
    nspliceends = result;
    trim_qstarts = ambig_qstarts;
  }
#ifdef CHECK_ASSERTIONS
  if (nspliceends > 1) {
    assert(trim_qstarts[0] < trim_qstarts[1]);
  }
#endif
  
  sense_qstart_paths = (List_T) NULL;
  for (i = 0; i < nspliceends; i++) {
    qstart = trim_qstarts[i];
    path = Path_new_for_qstart_extension(middle_diagonal_univdiagonal,qstart/*instead of middle_diagonal_qstart*/,middle_diagonal_qend,
					 splice5p,splicetype5,ambig_probs_5[i],intlistpool,univcoordlistpool);
    debug13(printf("For SENSE_FORWARD, plusp %d, revised qstart from %d to %d: ",plusp,middle_diagonal_qstart,qstart));
    debug13(Path_print(path,chroffset));
    
    paths = compute_qstart_paths(/*depth*/0,&terminalp,path,qstart_diagonals,
				 chrnum,chroffset,chrhigh,chrlength,querylength,
				 spliceinfo,mismatch_positions_alloc,query_compress,
				 plusp,genestrand,nmismatches_allowed,intlistpool,univcoordlistpool,
				 listpool,univdiagpool,/*sensedir*/SENSE_FORWARD,innerp);
    for (p = paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      debug13(printf("Path for qstart_local (sense): "));
      debug13(Path_print(path,chroffset));
      sense_qstart_paths = List_append(sense_qstart_paths,
				       compute_qstart_local(/*depth*/0,path,queryptr,querylength,localdb1,local_diagonals_cache,
							    &all_tables,mismatch_positions_alloc,spliceinfo,
							    stream_alloc,streamsize_alloc,query_compress,
							    chrnum,chroffset,chrhigh,chrlength,intlistpool,
							    univcoordlistpool,listpool,nmismatches_allowed,
							    plusp,genestrand,/*sensedir*/SENSE_FORWARD,innerp));
    }
    /* List_free(&paths); -- allocated by Listpool_push */
  }
  
#if 0
  if (any_extended_p == false) {
    if (trim_qstarts == NULL) {
      path = Path_new_for_qstart_extension(middle_diagonal_univdiagonal,middle_diagonal_qstart,middle_diagonal_qend,
					   /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,intlistpool,univcoordlistpool);

    } else {
      path = Path_new_for_qstart_extension(middle_diagonal_univdiagonal,trim_qstarts[0],middle_diagonal_qend,
					   splice5p,splicetype5,ambig_probs_5[0],intlistpool,univcoordlistpool);
    }
    sense_qstart_paths = Listpool_push(NULL,listpool,(void *) path);
  }
#endif

  if (splice5p == true) {
    FREE(ambig_qstarts);
    FREE(ambig_probs_5);
  }

#ifdef DEBUG13
  for (p = sense_qstart_paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    printf("** Sense left path: ");
    Path_print(path,chroffset);
  }
#endif


  /* Qstart antisense */
  if (novelsplicingp == true) {
    debug13(printf("Calling Substring_trimmed_qstarts with %d..%d\n",middle_diagonal_qstart,middle_diagonal_qend));
    splice5p = Substring_trimmed_qstarts(&result,&splicetype5,&ambig_qstarts,&ambig_probs_5,
					 /*left*/middle_diagonal_univdiagonal,
					 /*qstart:middle_diagonal_qstart,*/
					 /*qend*/middle_diagonal_qend,plusp,genestrand,
					 mismatch_positions_alloc,query_compress,
					 chroffset,/*sensedir*/SENSE_ANTI);
    debug13(printf("splice5p %d, result %d\n",splice5p,result));
  
    if (result < 0) {
      nspliceends = 0;
      /* trim_qstarts = (int *) NULL; */
    } else if (splice5p == false) {
      nspliceends = 1;
      trim_qstarts = &result;
      ambig_probs_5 = &zero_prob;
    } else {
      nspliceends = result;
      trim_qstarts = ambig_qstarts;
    }
#ifdef CHECK_ASSERTIONS
    if (nspliceends > 1) {
      assert(trim_qstarts[0] < trim_qstarts[1]);
    }
#endif

    antisense_qstart_paths = (List_T) NULL;
    for (i = 0; i < nspliceends; i++) {
      qstart = trim_qstarts[i];
      path = Path_new_for_qstart_extension(middle_diagonal_univdiagonal,qstart/*instead of middle_diagonal_qstart*/,middle_diagonal_qend,
					   splice5p,splicetype5,ambig_probs_5[i],intlistpool,univcoordlistpool);
      debug13(printf("For SENSE_ANTI, plusp %d, revised qstart from %d to %d: ",plusp,middle_diagonal_qstart,qstart));
      debug13(Path_print(path,chroffset));

      paths = compute_qstart_paths(/*depth*/0,&terminalp,path,qstart_diagonals,
				   chrnum,chroffset,chrhigh,chrlength,querylength,
				   spliceinfo,mismatch_positions_alloc,query_compress,
				   plusp,genestrand,nmismatches_allowed,intlistpool,univcoordlistpool,
				   listpool,univdiagpool,/*sensedir*/SENSE_ANTI,innerp);
      for (p = paths; p != NULL; p = List_next(p)) {
	path = (T) List_head(p);
	debug13(printf("Path for qstart_local (antisense): "));
	debug13(Path_print(path,chroffset));
	antisense_qstart_paths = List_append(antisense_qstart_paths,
					     compute_qstart_local(/*depth*/0,path,queryptr,querylength,
								  localdb1,local_diagonals_cache,
								  &all_tables,mismatch_positions_alloc,spliceinfo,
								  stream_alloc,streamsize_alloc,query_compress,
								  chrnum,chroffset,chrhigh,chrlength,intlistpool,
								  univcoordlistpool,listpool,nmismatches_allowed,
								  plusp,genestrand,/*sensedir*/SENSE_ANTI,innerp));
      }
      /* List_free(&paths); -- allocated by Listpool_push */
    }

#if 0
    if (any_extended_p == false) {
      if (trim_qstarts == NULL) {
	path = Path_new_for_qstart_extension(middle_diagonal_univdiagonal,middle_diagonal_qstart,middle_diagonal_qend,
					     /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,intlistpool,univcoordlistpool);
      } else {
	path = Path_new_for_qstart_extension(middle_diagonal_univdiagonal,trim_qstarts[0],middle_diagonal_qend,
					     splice5p,splicetype5,ambig_probs_5[0],intlistpool,univcoordlistpool);
      }
      antisense_qstart_paths = Listpool_push(NULL,listpool,(void *) path);
    }
#endif
    
    if (splice5p == true) {
      FREE(ambig_qstarts);
      FREE(ambig_probs_5);
    }

#ifdef DEBUG13
    for (p = antisense_qstart_paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      printf("** Antisense left path: ");
      Path_print(path,chroffset);
    }
#endif
  }

  for (p = all_tables; p != NULL; p = List_next(p)) {
    table = (Univcoordtable_T) List_head(p);
    diagonals_stores = (Diagonals_store_T *) Univcoordtable_saved_values(&nvalues,table);
    for (i = 0; i < nvalues; i++) {
      Diagonals_store_free(&(diagonals_stores[i]));
    }
    /* FREE(diagonals_stores); -- Not allocated by Univcoordtable_saved_values */
    Univcoordtable_free(&table);
  }
  /* List_free(&all_tables); -- allocated by Listpool_push */
  FREE(local_diagonals_cache);


  /* Qend */
  local_diagonals_cache = (Univcoordtable_T *) CALLOC(querylength,sizeof(Univcoordtable_T));
  all_tables = (List_T) NULL;
  if (paired_end_p == false) {
    innerp = false;
  } else if (first_read_p == plusp) {
    innerp = true;
  } else {
    innerp = false;
  }

  /* Qend sense.  Also handles the non-splicing case */
  debug13(printf("Calling Substring_trimmed_qends with %d..%d\n",middle_diagonal_qstart,middle_diagonal_qend));
  splice3p = Substring_trimmed_qends(&result,&splicetype3,&ambig_qends,&ambig_probs_3,
				     /*left*/middle_diagonal_univdiagonal,
				     /*qstart*/middle_diagonal_qstart,
				     /*qend:middle_diagonal_qend,*/
				     querylength,plusp,genestrand,
				     mismatch_positions_alloc,query_compress,
				     chroffset,/*sensedir*/SENSE_FORWARD);
  debug13(printf("splice3p %d, result %d\n",splice3p,result));
  
  if (result < 0) {
    nspliceends = 0;
    /* trim_qends = (int *) NULL; */
  } else if (splice3p == false) {
    nspliceends = 1;
    trim_qends = &result;
    ambig_probs_3 = &zero_prob;
  } else {
    nspliceends = result;
    trim_qends = ambig_qends;
  }
#ifdef CHECK_ASSERTIONS
  if (nspliceends > 1) {
    assert(trim_qends[0] > trim_qends[1]);
  }
#endif

  sense_qend_paths = (List_T) NULL;
  for (i = 0; i < nspliceends; i++) {
    qend = trim_qends[i];
    path = Path_new_for_qend_extension(middle_diagonal_univdiagonal,middle_diagonal_qstart,qend/*instead of middle_diagonal_qend*/,
				       splice3p,splicetype3,ambig_probs_3[i],intlistpool,univcoordlistpool);
    debug13(printf("For SENSE_FORWARD, plusp %d, revised qend from %d to %d: ",plusp,middle_diagonal_qend,qend));
    debug13(Path_print(path,chroffset));
    
    paths = compute_qend_paths(/*depth*/0,&terminalp,path,qend_diagonals,
			       chrnum,chroffset,chrhigh,chrlength,querylength,
			       spliceinfo,mismatch_positions_alloc,query_compress,
			       plusp,genestrand,nmismatches_allowed,intlistpool,univcoordlistpool,
			       listpool,univdiagpool,/*sensedir*/SENSE_FORWARD,innerp);
    
    for (p = paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      debug13(printf("Path for qend_local (sense): "));
      debug13(Path_print(path,chroffset));
      sense_qend_paths = List_append(sense_qend_paths,
				     compute_qend_local(/*depth*/0,path,queryptr,querylength,
							localdb1,local_diagonals_cache,
							&all_tables,mismatch_positions_alloc,spliceinfo,
							stream_alloc,streamsize_alloc,query_compress,
							chrnum,chroffset,chrhigh,chrlength,intlistpool,
							univcoordlistpool,listpool,nmismatches_allowed,
							plusp,genestrand,/*sensedir*/SENSE_FORWARD,innerp));
    }
    /* List_free(&paths); -- allocated by Listpool_push */
  }
  
#if 0
  if (any_extended_p == false) {
    if (trim_qends == NULL) {
      path = Path_new_for_qend_extension(middle_diagonal_univdiagonal,middle_diagonal_qstart,middle_diagonal_qend,
					 /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,intlistpool,univcoordlistpool);
    } else {
      path = Path_new_for_qend_extension(middle_diagonal_univdiagonal,middle_diagonal_qstart,trim_qends[0],
					 splice3p,splicetype3,ambig_probs_3[0],intlistpool,univcoordlistpool);
    }
    sense_qend_paths = Listpool_push(NULL,listpool,(void *) path);
  }
#endif

  if (splice3p == true) {
    FREE(ambig_qends);
    FREE(ambig_probs_3);
  }

#ifdef DEBUG13
  for (p = sense_qend_paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    printf("** Sense right path: ");
    Path_print(path,chroffset);
  }
#endif


  /* Qend antisense */
  if (novelsplicingp == true) {
    debug13(printf("Calling Substring_trimmed_qends with %d..%d\n",middle_diagonal_qstart,middle_diagonal_qend));
    splice3p = Substring_trimmed_qends(&result,&splicetype3,&ambig_qends,&ambig_probs_3,
				       /*left*/middle_diagonal_univdiagonal,
				       /*qstart*/middle_diagonal_qstart,
				       /*qend:middle_diagonal_qend,*/
				       querylength,plusp,genestrand,
				       mismatch_positions_alloc,query_compress,
				       chroffset,/*sensedir*/SENSE_ANTI);
    debug13(printf("splice3p %d, result %d\n",splice3p,result));

    if (result < 0) {
      nspliceends = 0;
      /* trim_qends = (int *) NULL; */
    } else if (splice3p == false) {
      nspliceends = 1;
      trim_qends = &result;
      ambig_probs_3 = &zero_prob;
    } else {
      nspliceends = result;
      trim_qends = ambig_qends;
    }
#ifdef CHECK_ASSERTIONS
    if (nspliceends > 1) {
      assert(trim_qends[0] > trim_qends[1]);
    }
#endif

    antisense_qend_paths = (List_T) NULL;
    for (i = 0; i < nspliceends; i++) {
      qend = trim_qends[i];
      path = Path_new_for_qend_extension(middle_diagonal_univdiagonal,middle_diagonal_qstart,qend/*instead of middle_diagonal_qend*/,
					 splice3p,splicetype3,ambig_probs_3[i],intlistpool,univcoordlistpool);
      debug13(printf("For SENSE_ANTI, plusp %d, revised qend from %d to %d: ",plusp,middle_diagonal_qend,qend));
      debug13(Path_print(path,chroffset));

      paths = compute_qend_paths(/*depth*/0,&terminalp,path,qend_diagonals,
				 chrnum,chroffset,chrhigh,chrlength,querylength,
				 spliceinfo,mismatch_positions_alloc,query_compress,
				 plusp,genestrand,nmismatches_allowed,intlistpool,univcoordlistpool,
				 listpool,univdiagpool,/*sensedir*/SENSE_ANTI,innerp);
      for (p = paths; p != NULL; p = List_next(p)) {
	path = (T) List_head(p);
	debug13(printf("Path for qend_local (antisense): "));
	debug13(Path_print(path,chroffset));
	antisense_qend_paths = List_append(antisense_qend_paths,
					   compute_qend_local(/*depth*/0,path,queryptr,querylength,
							      localdb1,local_diagonals_cache,
							      &all_tables,mismatch_positions_alloc,spliceinfo,
							      stream_alloc,streamsize_alloc,query_compress,
							      chrnum,chroffset,chrhigh,chrlength,intlistpool,
							      univcoordlistpool,listpool,nmismatches_allowed,
							      plusp,genestrand,/*sensedir*/SENSE_ANTI,innerp));
      }
      /* List_free(&paths); -- allocated by Listpool_push */
    }
    
#if 0
    if (any_extended_p == false) {
      if (trim_qends == NULL) {
	path = Path_new_for_qend_extension(middle_diagonal_univdiagonal,middle_diagonal_qstart,middle_diagonal_qend,
					   /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,intlistpool,univcoordlistpool);
      } else {
	path = Path_new_for_qend_extension(middle_diagonal_univdiagonal,middle_diagonal_qstart,trim_qends[0],
					   splice3p,splicetype3,ambig_probs_3[0],intlistpool,univcoordlistpool);
      }
      antisense_qend_paths = Listpool_push(NULL,listpool,(void *) path);
    }
#endif
    
    if (splice3p == true) {
      FREE(ambig_qends);
      FREE(ambig_probs_3);
    }

#ifdef DEBUG13
    for (p = antisense_qend_paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      printf("** Antisense right path: ");
      Path_print(path,chroffset);
      printf("\n");
    }
#endif
  }

  for (p = all_tables; p != NULL; p = List_next(p)) {
    table = (Univcoordtable_T) List_head(p);
    diagonals_stores = (Diagonals_store_T *) Univcoordtable_saved_values(&nvalues,table);
    for (i = 0; i < nvalues; i++) {
      Diagonals_store_free(&(diagonals_stores[i]));
    }
    /* FREE(diagonals_stores); -- Not allocated by Univcoordtable_saved_values */
    Univcoordtable_free(&table);
  }
  /* List_free(&all_tables); -- allocated by Listpool_push */
  FREE(local_diagonals_cache);


  *foundp = false;
  hits = solve_via_segments_genome(&(*foundp),&(*found_score_overall),&(*found_score_within_trims),hits,
				   sense_qstart_paths,sense_qend_paths,/*sensedir*/SENSE_FORWARD,
				   chrnum,chroffset,chrhigh,chrlength,querylength,query_compress,
				   plusp,genestrand,intlistpool,univcoordlistpool,
				   listpool,hitlistpool,method,level);
  Path_gc(&sense_qstart_paths);
  Path_gc(&sense_qend_paths);

  if (novelsplicingp == true) {
    hits = solve_via_segments_genome(&(*foundp),&(*found_score_overall),&(*found_score_within_trims),hits,
				     antisense_qstart_paths,antisense_qend_paths,/*sensedir*/SENSE_ANTI,
				     chrnum,chroffset,chrhigh,chrlength,querylength,query_compress,
				     plusp,genestrand,intlistpool,univcoordlistpool,
				     listpool,hitlistpool,method,level);
    Path_gc(&antisense_qstart_paths);
    Path_gc(&antisense_qend_paths);
  }

  return hits;
}


void
Path_solve_setup (Genome_T genomebits_in, Genome_T genomebits_alt_in, bool *circularp_in,
		  Localdb_T localdb_in, Localdb_T localdb2_in,
		  Chrpos_T shortsplicedist_novelend, int min_intronlength_in, int max_deletionlength,
		  int max_insertionlen_in, 
		  bool novelsplicingp_in, Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
		  Chrpos_T *splicedists_in, int nsplicesites_in,
		  int index1part_in, int index1interval_in, int local1part_in) {

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  circularp = circularp_in;

  localdb = localdb_in;
  localdb2 = localdb2_in;

  min_intronlength = min_intronlength_in;
  max_deletionlen = max_deletionlength;

  novelsplicingp = novelsplicingp_in;
  splicesites = splicesites_in;
  splicetypes = splicetypes_in;
  splicedists = splicedists_in;
  nsplicesites = nsplicesites_in;

  /* For localdb */
  max_insertionlen = max_insertionlen_in;
  if (shortsplicedist_novelend > max_deletionlen) {
    overall_end_distance_genome = shortsplicedist_novelend;
  } else {
    overall_end_distance_genome = max_deletionlen;
  }

  index1part = index1part_in;
  index1interval = index1interval_in;
  local1part = local1part_in;

  return;
}


