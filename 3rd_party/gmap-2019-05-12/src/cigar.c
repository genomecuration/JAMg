static char rcsid[] = "$Id: cigar.c 218383 2019-02-15 23:21:46Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>		/* For isupper */

#include "cigar.h"
#include "mem.h"
#include "complement.h"
#include "stage3hr.h"



#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* print_extended_cigar */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Cigar_compute_main */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Cigar_compute_supplemental */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif



static bool hide_soft_clips_p;
static bool cigar_extended_p = false;
static bool merge_samechr_p;
static bool md_lowercase_variant_p;
static bool sam_hardclip_use_S_p;
static bool sam_insert_0M_p;


void
Cigar_setup (bool cigar_extended_p_in, bool hide_soft_clips_p_in,
	     bool merge_samechr_p_in, bool md_lowercase_variant_p_in,
	     bool sam_hardclip_use_S_p_in, bool sam_insert_0M_p_in) {
  cigar_extended_p = cigar_extended_p_in;
  hide_soft_clips_p = hide_soft_clips_p_in;
  merge_samechr_p = merge_samechr_p_in;
  md_lowercase_variant_p = md_lowercase_variant_p_in;
  sam_hardclip_use_S_p = sam_hardclip_use_S_p_in;
  sam_insert_0M_p = sam_insert_0M_p_in;
  return;
}


/* Modified from print_cigar and assuming type == M */
static int
length_cigar_M (int stringlength, int querypos, int querylength,
		int hardclip_low, int hardclip_high, bool plusp, bool lastp) {
  int length = 0;
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;

  if (plusp == true) {
    debug1(printf("\nEntering length_cigar_M with stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Ignoring initial %dH\n",cliplength));
	/* FPRINTF(fp,"%dH",cliplength); */
      }
      matchlength = endpos - startpos;
      if (matchlength <= 0) {
	/* Skip */
      } else {
	/* type != 'E' */
	debug1(printf("  Adding length of  %d\n",matchlength));
	length += matchlength;	/* FPRINTF(fp,"%d%c",matchlength,type); */
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Ignoring final %dH\n",cliplength));
	/* FPRINTF(fp,"%dH",cliplength); */
      }
    }

  } else {
    debug1(printf("\nEntering length_cigar_M with stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Ignoring initial %dH\n",cliplength));
	/* FPRINTF(fp,"%dH",cliplength); */
      }
      matchlength = startpos - endpos;
      if (matchlength <= 0) {
	/* Skip */
      } else {
	/* type != 'E' */
	debug1(printf("  Adding length of %d\n",matchlength));
	length += matchlength;  /* FPRINTF(fp,"%d%c",matchlength,type); */
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Ignoring final %dH\n",cliplength));
	/* FPRINTF(fp,"%dH",cliplength); */
      }
    }
  }

  return length;
}



/* Returns true if prints the requested type */
static bool
print_cigar (Filestring_T fp, char type, int stringlength, int querypos, int querylength,
	     int hardclip_low, int hardclip_high, bool plusp, bool lastp, int trimlength) {
  bool printp = false;
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  
  if (plusp == true) {
    debug1(printf("\nEntering print_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	if (sam_hardclip_use_S_p) {
	  FPRINTF(fp,"%dS",cliplength);
	} else {
	  FPRINTF(fp,"%dH",cliplength);
	}
      }
      matchlength = endpos - startpos;
      if (matchlength <= 0 && (type != 'M' || sam_insert_0M_p == false)) {
	/* Skip, except for 0M if desired */
      } else if (type != 'E') {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	FPRINTF(fp,"%d%c",matchlength,type);
	printp = true;
      } else if (matchlength == trimlength) {
	debug1(printf("  Pushing %dS\n",matchlength));
	FPRINTF(fp,"%dS",matchlength);
      } else {
	debug1(printf("  Pushing %dH because matchlength %d != trimlength %d\n",
		      matchlength,matchlength,trimlength));
	if (sam_hardclip_use_S_p) {
	  FPRINTF(fp,"%dS",matchlength);
	} else {
	  FPRINTF(fp,"%dH",matchlength);
	}
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	if (sam_hardclip_use_S_p) {
	  FPRINTF(fp,"%dS",cliplength);
	} else {
	  FPRINTF(fp,"%dH",cliplength);
	}
      }
    }

  } else {
    debug1(printf("\nEntering print_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	if (sam_hardclip_use_S_p) {
	  FPRINTF(fp,"%dS",cliplength);
	} else {
	  FPRINTF(fp,"%dH",cliplength);
	}
      }
      matchlength = startpos - endpos;
      if (matchlength <= 0 && (type != 'M' || sam_insert_0M_p == false)) {
	/* Skip, except for 0M if desired */
      } else if (type != 'E') {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	FPRINTF(fp,"%d%c",matchlength,type);
	printp = true;
      } else if (matchlength == trimlength) {
	debug1(printf("  Pushing %dS\n",matchlength));
	FPRINTF(fp,"%dS",matchlength);
      } else {
	debug1(printf("  Pushing %dH because matchlength %d != trimlength %d\n",
		      matchlength,matchlength,trimlength));
	if (sam_hardclip_use_S_p) {
	  FPRINTF(fp,"%dS",matchlength);
	} else {
	  FPRINTF(fp,"%dH",matchlength);
	}
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	if (sam_hardclip_use_S_p) {
	  FPRINTF(fp,"%dS",cliplength);
	} else {
	  FPRINTF(fp,"%dH",cliplength);
	}
      }
    }
  }

  return printp;
}



/* Returns true if we printed an = or X entry */
/* Based on print_md_string */
static bool
print_extended_cigar (Filestring_T fp, char *genomicfwd_refdiff,
		      int stringlength, int querypos, int querylength,
		      int hardclip_low, int hardclip_high, bool plusp, bool lastp) {
  bool printp = false;
  int nmatches = 0, nmismatches = 0;
  int starti, endi, i;
  bool hardclip_end_p = false;
  int cliplength, endpos;

  if (plusp == true) {
    debug2(printf("\nEntering print_extended_cigar with querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus: %s ref\n",
		  querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff));
    if (hardclip_low == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_low > querypos) {
      /* startpos = hardclip_low; */
      starti = hardclip_low - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_low %d - querypos %d\n",
		    starti,hardclip_low,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      endi = (querylength - hardclip_high) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_high %d) - querypos %d\n",
		    endi,querylength,hardclip_high,querypos));
    } else {
      endpos = querypos + stringlength;
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd_refdiff == NULL) {
      if (endi > starti) {
	nmatches += (endi - starti);
      }

    } else if (md_lowercase_variant_p == false) {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  if (nmismatches > 0 /*|| hardclip_end_p == true*/) {
	    FPRINTF(fp,"%dX",nmismatches);
	    nmismatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmatches++;

	} else {
	  /* A true mismatch against both variants */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmismatches++;
	}
      }

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  if (nmismatches > 0 /*|| hardclip_end_p == true*/) {
	    FPRINTF(fp,"%dX",nmismatches);
	    nmismatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmatches++;

#if 0
	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmismatches++;
#endif

	} else {
	  /* A true mismatch against both variants */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmismatches++;
	}
      }
    }

    if (nmatches > 0) {
      FPRINTF(fp,"%d=",nmatches);
      printp = true;
    } else if (nmismatches > 0) {
      FPRINTF(fp,"%dX",nmismatches);
      printp = true;
    }

    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	if (sam_hardclip_use_S_p) {
	  FPRINTF(fp,"%dS",cliplength);
	} else {
	  FPRINTF(fp,"%dH",cliplength);
	}
      }
    }

  } else {
    debug2(printf("\nEntering print_extended_cigar with querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus: %s ref\n",
		  querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff));
    querypos = querylength - querypos - stringlength;
    debug2(printf("  Revising querypos to be %d\n",querypos));

    if (hardclip_low == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_low > querypos) {
      /* startpos = hardclip_low; */
      starti = hardclip_low - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_low %d - querypos %d\n",
		    starti,hardclip_low,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      endi = (querylength - hardclip_high) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_high %d) - querypos %d\n",
		    endi,querylength,hardclip_high,querypos));
    } else {
      endpos = querypos + stringlength;
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd_refdiff == NULL) {
      if (endi > starti) {
	nmatches += (endi - starti);
      }

    } else if (md_lowercase_variant_p == false) {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  if (nmismatches > 0 /*|| hardclip_end_p == true*/) {
	    FPRINTF(fp,"%dX",nmismatches);
	    nmismatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmatches++;

	} else {
	  /* A true mismatch against both variants */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmismatches++;
	}
      }

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  if (nmismatches > 0 /*|| hardclip_end_p == true*/) {
	    FPRINTF(fp,"%dX",nmismatches);
	    nmismatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmatches++;

#if 0
	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmismatches++;
#endif

	} else {
	  /* A true mismatch against both variants */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	    printp = true;
	  }
	  nmismatches++;
	}
      }
    }

    if (nmatches > 0) {
      FPRINTF(fp,"%d=",nmatches);
      printp = true;
    } else if (nmismatches > 0) {
      FPRINTF(fp,"%dX",nmismatches);
      printp = true;
    }

    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	if (sam_hardclip_use_S_p) {
	  FPRINTF(fp,"%dS",cliplength);
	} else {
	  FPRINTF(fp,"%dH",cliplength);
	}
      }
    }
  }

  return printp;
}


static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}


/* Returns true if we printed an M, =, or X entry */
static bool
print_cigar_M (Filestring_T fp, Substring_T substring, int substring_length, int substring_start,
	       int stringlength, int querypos, int querylength,
	       int hardclip_low, int hardclip_high, bool plusp, bool lastp, int trimlength) {
  bool printp;
  char *genomicfwd_refdiff, *genomicdir_refdiff;
  
  if (cigar_extended_p == false) {
    return print_cigar(fp,/*type*/'M',stringlength,querypos,querylength,
		       hardclip_low,hardclip_high,plusp,lastp,trimlength);

  } else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == NULL) {
    return print_extended_cigar(fp,/*genomicfwd_refdiff*/NULL,/*stringlength*/substring_length,
				/*querypos*/substring_start,querylength,
				hardclip_low,hardclip_high,plusp,lastp);
  } else if (plusp == true) {
    return print_extended_cigar(fp,&(genomicdir_refdiff[substring_start]),/*stringlength*/substring_length,
				/*querypos*/substring_start,querylength,
				hardclip_low,hardclip_high,plusp,lastp);
  } else {
    genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
    printp = print_extended_cigar(fp,genomicfwd_refdiff,/*stringlength*/substring_length,
				  /*querypos*/substring_start,querylength,
				  hardclip_low,hardclip_high,plusp,lastp);
    FREEA(genomicfwd_refdiff);
    return printp;
  }
}


#if 0
/* Copy also in pair.c for GMAP */
static bool
check_cigar_types (Intlist_T cigar_types) {
  Intlist_T p;
  int type;
  bool M_present_p = false;

  for (p = cigar_types; p != NULL; p = Intlist_next(p)) {
    type = Intlist_head(p);
    if (type == 'M') {
      M_present_p = true;
#if 0
    } else if (type == 'H' && last_type == 'S') {
      debug1(printf("check_cigar_types detects adjacent S and H, so returning false\n"));
      return false;
    } else if (type == 'S' && last_type == 'H') {
      debug1(printf("check_cigar_types detects adjacent S and H, so returning false\n"));
      return false;
#endif
    }
  }

  return M_present_p;
}
#endif



/* Modified from Cigar_print_substrings */
int
Cigar_length_substrings (Stage3end_T stage3end, int querylength, int hardclip_low, int hardclip_high) {
  int length = 0;
  List_T startp, startq, endp;
  /* List_T prevp, nextp, finalp; */

  Substring_T substring, substringL, substringH;
  /* Junction_T post_junction; */
  /* int type; */

  List_T substrings_LtoH, junctions_LtoH;
  List_T p, q;
  /* int substring_start, substring_length; */

  bool plusp;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif

  
  if ((substrings_LtoH = Stage3end_substrings_LtoH(stage3end)) == NULL) {
    return 0;
  } else {
    plusp = Stage3end_plusp(stage3end);
    substrings_LtoH = Stage3end_substrings_LtoH(stage3end);
    junctions_LtoH = Stage3end_junctions_LtoH(stage3end);
    substringL = (Substring_T) List_head(substrings_LtoH);
    substringH = (Substring_T) List_last_value(substrings_LtoH);
  }

  if (Substring_has_alts_p(substringL) == true) {
    /* prevp = substrings_LtoH; */
    startp = List_next(substrings_LtoH);
    startq = List_next(junctions_LtoH);
  } else {
    /* prevp = (List_T) NULL; */
    startp = substrings_LtoH;
    startq = junctions_LtoH;
  }
  if (Substring_has_alts_p(substringH) == true) {
    endp = List_last_item(substrings_LtoH);
  } else {
    endp = (List_T) NULL;
  }

  debug(printf("End has %d substrings\n",List_length(substrings_LtoH)));

  p = startp;
  q = startq;

  if (plusp == true) {
    /* Plus */
    while (p != endp && Substring_queryend((Substring_T) List_head(p)) < hardclip_low) {
      /* Skip, because substring entirely in hard-clipped region */
      debug(printf("Skipping %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));
      /* prevp = p; */
      p = List_next(p);
      q = List_next(q);
    }

    if (p == NULL) {
      debug(printf("Empty substring\n"));

    } else {
      substring = (Substring_T) List_head(p);
      if (List_next(p) == endp || Substring_queryend(substring) >= querylength - hardclip_high) {
	/* Single substring */
	debug(printf("Single substring %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	
	if (hide_soft_clips_p == true) {
	  /* substring_start = Substring_querystart_pretrim(substring); */
	  /* substring_length = Substring_match_length_pretrim(substring); */
	  length += length_cigar_M(Substring_querystart(substring) + Substring_match_length(substring) +
				   (querylength - Substring_queryend(substring)),/*querypos*/0,querylength,
				   hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
	} else {
	  /* substring_start = Substring_querystart(substring); */
	  /* substring_length = Substring_match_length(substring); */
#if 0
	  print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		      /*querypos*/0,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/false,/*trimlength*/0);
#endif
	  length += length_cigar_M(Substring_match_length(substring),
				   /*querypos*/Substring_querystart(substring),querylength,
				   hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
#if 0
	  print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		      /*querypos*/Substring_queryend(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
#endif
	}
	/* finalp = p; */
	/* nextp = List_next(p); */
	
      } else {
	/* First substring, plus */
	debug(printf("First substring, plus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	
	/* post_junction = (Junction_T) List_head(q); */
	
	if (hide_soft_clips_p == true) {
	  /* substring_start = Substring_querystart_pretrim(substring); */
	  /* substring_length = Substring_match_length_pretrim(substring); */
	  length += length_cigar_M(Substring_querystart(substring) +
				   Substring_match_length(substring),
				   /*querypos*/0,querylength,hardclip_low,hardclip_high,
				   /*plusp*/true,/*lastp*/false);
	} else {
	  /* substring_start = Substring_querystart(substring); */
	  /* substring_length = Substring_match_length(substring); */
#if 0
	  print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		      /*querypos*/0,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/false,/*trimlength*/0);
#endif
	  length += length_cigar_M(Substring_match_length(substring),
				   /*querypos*/Substring_querystart(substring),querylength,
				   hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	}
	p = List_next(p);
      
	while (p != endp && Substring_queryend((Substring_T) List_head(p)) < querylength - hardclip_high) {
#if 0
	  if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    debug1(printf("1. Pushing %dD\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	    nindels += Junction_nindels(post_junction);
	  } else if (type == INS_JUNCTION) {
	    debug1(printf("1. Pushing %dI\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dI",Junction_nindels(post_junction));
	    nindels += Junction_nindels(post_junction);
	  } else if (type == SPLICE_JUNCTION) {
	    debug1(printf("1. Pushing %dN\n",Junction_splice_distance(post_junction)));
	    FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	  }
#endif
	  q = List_next(q);
#if 0
	  if (q == NULL) {
	  } else {
	    post_junction = (Junction_T) List_head(q);
	  }
#endif
	  
	  substring = (Substring_T) List_head(p);
	  if (List_next(p) == endp) {
	    /* Last substring, plus, not hard-clipped */
	    debug(printf("Last substring, plus, not hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
			 Substring_queryend((Substring_T) List_head(p))));
	    
	    if (hide_soft_clips_p == true) {
	      /* substring_start = Substring_querystart_pretrim(substring); */
	      /* substring_length = Substring_match_length_pretrim(substring); */
	      length += length_cigar_M(Substring_match_length(substring) +
				       (querylength - Substring_queryend(substring)),
				       /*querypos*/Substring_querystart(substring),querylength,
				       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
	    } else {
	      /* substring_start = Substring_querystart(substring); */
	      /* substring_length = Substring_match_length(substring); */
	      length += length_cigar_M(Substring_match_length(substring),
				       /*querypos*/Substring_querystart(substring),querylength,
				       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
#if 0
	      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
			  /*querypos*/Substring_queryend(substring),querylength,
			  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
#endif
	    }
	    /* finalp = p; */
	    /* nextp = List_next(p); */
	    
	  } else {
	    /* Middle substring, plus */
	    debug(printf("Middle substring, plus %d..%d\n",Substring_querystart((Substring_T) List_head(p)), 
			 Substring_queryend((Substring_T) List_head(p))));
	    /* substring_start = Substring_querystart(substring); */
	    /* substring_length = Substring_match_length(substring); */
	    
	    length += length_cigar_M(Substring_match_length(substring),
				     /*querypos*/Substring_querystart(substring),querylength,
				     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	  }
	  p = List_next(p);
	}
      
	if (p != endp) {
#if 0
	  if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    debug1(printf("2. Pushing %dD\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	    nindels += Junction_nindels(post_junction);
	  } else if (type == INS_JUNCTION) {
	    debug1(printf("2. Pushing %dI\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dI",Junction_nindels(post_junction));
	    nindels += Junction_nindels(post_junction);
	  } else if (type == SPLICE_JUNCTION) {
	    debug1(printf("2. Pushing %dN\n",Junction_splice_distance(post_junction)));
	    FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	  }
#endif
	  
	  /* Last substring, plus, hard-clipped */
	  substring = (Substring_T) List_head(p);
	  debug(printf("Last substring, plus, hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		       Substring_queryend((Substring_T) List_head(p))));
	  if (hide_soft_clips_p == true) {
	    /* substring_start = Substring_querystart_pretrim(substring); */
	    /* substring_length = Substring_match_length_pretrim(substring); */
	    length += length_cigar_M(Substring_match_length(substring) +
				     (querylength - Substring_queryend(substring)),
				     /*querypos*/Substring_querystart(substring),querylength,
				     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
	  } else {
	    /* substring_start = Substring_querystart(substring); */
	    /* substring_length = Substring_match_length(substring); */
	    length += length_cigar_M(Substring_match_length(substring),
				     /*querypos*/Substring_querystart(substring),querylength,
				     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
#if 0
	    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
#endif
	  }
	  /* finalp = p; */
	  /* nextp = List_next(p); */
	  
	}
      }
    }

  } else {
    /* Minus */
    while (p != endp && Substring_querystart((Substring_T) List_head(p)) >= querylength - hardclip_low) {
      /* Skip, because substring entirely in hard-clipped region */
      debug(printf("Skipping %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));
      /* prevp = p; */
      p = List_next(p);
      q = List_next(q);
    }

    if (p == NULL) {
      debug(printf("Empty substring\n"));

    } else {
      substring = (Substring_T) List_head(p);
      if (List_next(p) == endp || Substring_querystart(substring) < hardclip_high) {
	/* Single substring */
	debug(printf("Single substring %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	
	if (hide_soft_clips_p == true) {
	  /* substring_start = Substring_querystart_pretrim(substring); */
	  /* substring_length = Substring_match_length_pretrim(substring); */
	  length += length_cigar_M((querylength - Substring_queryend(substring)) + 
				   Substring_match_length(substring) + Substring_querystart(substring),
				   /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
				   /*plusp*/false,/*lastp*/true);
	} else {
	  /* substring_start = Substring_querystart(substring); */
	  /* substring_length = Substring_match_length(substring); */
#if 0
	  print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		      /*querypos*/querylength,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
#endif
	  length += length_cigar_M(Substring_match_length(substring),
				   /*querypos*/Substring_queryend(substring),querylength,
				   hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
#if 0
	  print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		      /*querypos*/Substring_querystart(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
#endif
	}
	/* finalp = p; */
	/* nextp = List_next(p); */
	
      } else {
	/* First substring, minus */
	debug(printf("First substring, minus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	
	/* post_junction = (Junction_T) List_head(q); */
	
	if (hide_soft_clips_p == true) {
	  /* substring_start = Substring_querystart_pretrim(substring); */
	  /* substring_length = Substring_match_length_pretrim(substring); */
	  length += length_cigar_M((querylength - Substring_queryend(substring)) +
				   Substring_match_length(substring),
				   /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
				   /*plusp*/false,/*lastp*/false);
	} else {
	  /* substring_start = Substring_querystart(substring); */
	  /* substring_length = Substring_match_length(substring); */
#if 0
	  print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		      /*querypos*/querylength,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
#endif
	  length += length_cigar_M(Substring_match_length(substring),
				   /*querypos*/Substring_queryend(substring),querylength,
				   hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	}
	p = List_next(p);
	
	while (p != endp && Substring_querystart((Substring_T) List_head(p)) >= hardclip_high) {
#if 0
	  if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    debug1(printf("3. Pushing %dD\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	    nindels += Junction_nindels(post_junction);
	  } else if (type == INS_JUNCTION) {
	    debug1(printf("3. Pushing %dI\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dI",Junction_nindels(post_junction));
	    nindels += Junction_nindels(post_junction);
	  } else if (type == SPLICE_JUNCTION) {
	    debug1(printf("3. Pushing %dN\n",Junction_splice_distance(post_junction)));
	    FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	  }
#endif
	  q = List_next(q);
#if 0
	  if (q == NULL) {
	  } else {
	    post_junction = (Junction_T) List_head(q);
	  }
#endif
	  
	  substring = (Substring_T) List_head(p);
	  if (List_next(p) == endp) {
	    /* Last substring, minus, not hard-clipped */
	    debug(printf("Last substring, minus, not hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
			 Substring_queryend((Substring_T) List_head(p))));
	    
	    if (hide_soft_clips_p == true) {
	      /* substring_start = Substring_querystart_pretrim(substring); */
	      /* substring_length = Substring_match_length_pretrim(substring); */
	      length += length_cigar_M(Substring_match_length(substring) +
				       Substring_querystart(substring),
				       /*querypos*/Substring_queryend(substring),querylength,
				       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	    } else {
	      /* substring_start = Substring_querystart(substring); */
	      /* substring_length = Substring_match_length(substring); */
	      length += length_cigar_M(Substring_match_length(substring),
				       /*querypos*/Substring_queryend(substring),querylength,
				       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
#if 0
	      print_cigar(fp,/*type*/'S',Substring_querystart(substring),
			  /*querypos*/Substring_querystart(substring),querylength,hardclip_low,hardclip_high,
			  /*plusp*/false,/*lastp*/true,/*trimlength*/0);
#endif
	    }
	    /* finalp = p; */
	    /* nextp = List_next(p); */
	    
	  } else {
	    /* Middle substring, minus */
	    debug(printf("Middle substring, minus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
			 Substring_queryend((Substring_T) List_head(p))));
	    /* substring_start = Substring_querystart(substring); */
	    /* substring_length = Substring_match_length(substring); */
	    
	    length += length_cigar_M(Substring_match_length(substring),
				     /*querypos*/Substring_queryend(substring),querylength,
				     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	  }
	  p = List_next(p);
	}
	
	if (p != endp) {
#if 0
	  if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    debug1(printf("4. Pushing %dD\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	    nindels += Junction_nindels(post_junction);
	  } else if (type == INS_JUNCTION) {
	    debug1(printf("4. Pushing %dI\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dI",Junction_nindels(post_junction));
	    nindels += Junction_nindels(post_junction);
	  } else if (type == SPLICE_JUNCTION) {
	    debug1(printf("4. Pushing %dN\n",Junction_splice_distance(post_junction)));
	    FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	  }
#endif
	  
	  /* Last substring, minus, hard-clipped */
	  substring = (Substring_T) List_head(p);
	  debug(printf("Last substring, minus, hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		       Substring_queryend((Substring_T) List_head(p))));
	  
	  if (hide_soft_clips_p == true) {
	    /* substring_start = Substring_querystart_pretrim(substring); */
	    /* substring_length = Substring_match_length_pretrim(substring); */
	    length += length_cigar_M(Substring_querystart(substring),
				     /*querypos*/Substring_queryend(substring),querylength,
				     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	  } else {
	    /* substring_start = Substring_querystart(substring); */
	    /* substring_length = Substring_match_length(substring); */
	    length += length_cigar_M(Substring_match_length(substring),
				     /*querypos*/Substring_queryend(substring),querylength,
				     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
#if 0
	    print_cigar(fp,/*type*/'S',Substring_querystart(substring),
			/*querypos*/Substring_querystart(substring),querylength,hardclip_low,hardclip_high,
			/*plusp*/false,/*lastp*/true,/*trimlength*/0);
#endif
	  }
	  /* finalp = p; */
	  /* nextp = List_next(p); */
	  
	}
      }
    }
  }

  return length;
}


static void
Cigar_print_substrings (int *nindels, List_T *startp, List_T *startq, List_T *prevp, List_T *nextp, List_T *finalp, List_T *endp,
			Filestring_T fp, Stage3end_T stage3end, int querylength,
			int hardclip_low, int hardclip_high) {
  Substring_T substring, substringL, substringH;
  Junction_T post_junction;
  int type;

  List_T substrings_LtoH, junctions_LtoH;
  List_T p, q;
  int substring_start, substring_length;

  bool M_printedp = false;
  bool plusp;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif

  
  debug(printf("\n\n***Entered Cigar_print_substrings with hardclip_low %d and hardclip_high %d\n",hardclip_low,hardclip_high));

  *nindels = 0;
  if ((substrings_LtoH = Stage3end_substrings_LtoH(stage3end)) == NULL) {
    FPRINTF(fp,"*");
    return;
  } else {
    plusp = Stage3end_plusp(stage3end);
    substrings_LtoH = Stage3end_substrings_LtoH(stage3end);
    junctions_LtoH = Stage3end_junctions_LtoH(stage3end);
    substringL = (Substring_T) List_head(substrings_LtoH);
    substringH = (Substring_T) List_last_value(substrings_LtoH);
  }


  if (Substring_has_alts_p(substringL) == true) {
    *prevp = substrings_LtoH;
    *startp = List_next(substrings_LtoH);
    *startq = List_next(junctions_LtoH);
  } else {
    *prevp = (List_T) NULL;
    *startp = substrings_LtoH;
    *startq = junctions_LtoH;
  }
  if (Substring_has_alts_p(substringH) == true) {
    *endp = List_last_item(substrings_LtoH);
  } else {
    *endp = (List_T) NULL;
  }

  debug(printf("End has %d substrings\n",List_length(substrings_LtoH)));

  p = *startp;
  q = *startq;
  if (plusp == true) {
    /* Plus */
    while (p != *endp && Substring_queryend((Substring_T) List_head(p)) < hardclip_low) {
      /* Skip, because substring entirely in hard-clipped region */
      debug(printf("Skipping %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));
      *prevp = p;
      p = List_next(p);
      q = List_next(q);
    }

    if (p == NULL) {
      debug(printf("Empty substring\n"));
      FPRINTF(fp,"*");

    } else {
      substring = (Substring_T) List_head(p);
      if (List_next(p) == *endp || Substring_queryend(substring) >= querylength - hardclip_high) {
	/* Single substring */
	debug(printf("Single substring %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	
	if (hide_soft_clips_p == true) {
	  substring_start = Substring_querystart_pretrim(substring);
	  substring_length = Substring_match_length_pretrim(substring);
	  M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
				      Substring_querystart(substring) + Substring_match_length(substring) +
				      (querylength - Substring_queryend(substring)),/*querypos*/0,querylength,
				      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	} else {
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);
	  print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		      /*querypos*/0,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	  M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
				      Substring_match_length(substring),
				      /*querypos*/Substring_querystart(substring),querylength,
				      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	  print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		      /*querypos*/Substring_queryend(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	}
	*finalp = p;
	*nextp = List_next(p);
	
      } else {
	/* First substring, plus */
	debug(printf("First substring, plus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	
	post_junction = (Junction_T) List_head(q);
	
	if (hide_soft_clips_p == true) {
	  substring_start = Substring_querystart_pretrim(substring);
	  substring_length = Substring_match_length_pretrim(substring);
	  M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
				      Substring_querystart(substring) +
				      Substring_match_length(substring),
				      /*querypos*/0,querylength,hardclip_low,hardclip_high,
				      /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	} else {
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);
	  print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		      /*querypos*/0,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	  M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
				      Substring_match_length(substring),
				      /*querypos*/Substring_querystart(substring),querylength,
				      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	}
	p = List_next(p);
	
	while (p != *endp && Substring_queryend((Substring_T) List_head(p)) < querylength - hardclip_high) {
	  if (M_printedp == false) {
	    debug1(printf("Skipping initial indel or splice\n"));

	  } else if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    debug1(printf("1. Pushing %dD\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	    *nindels += Junction_nindels(post_junction);

	  } else if (type == INS_JUNCTION) {
	    debug1(printf("1. Pushing %dI\n",Junction_nindels(post_junction)));
#if 0
	    FPRINTF(fp,"%dI",Junction_nindels(post_junction));
#else
	    print_cigar(fp,/*type*/'I',/*stringlength*/Junction_nindels(post_junction),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
#endif
	    *nindels += Junction_nindels(post_junction);

	  } else if (type == SPLICE_JUNCTION) {
	    debug1(printf("1. Pushing %dN\n",Junction_splice_distance(post_junction)));
	    FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	  }
	  q = List_next(q);
	  if (q == NULL) {
	  } else {
	    post_junction = (Junction_T) List_head(q);
	  }
	  
	  substring = (Substring_T) List_head(p);
	  if (List_next(p) == *endp) {
	    /* Last substring, plus, not hard-clipped */
	    debug(printf("Last substring, plus, not hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
			 Substring_queryend((Substring_T) List_head(p))));
	    
	    if (hide_soft_clips_p == true) {
	      substring_start = Substring_querystart_pretrim(substring);
	      substring_length = Substring_match_length_pretrim(substring);
	      M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					  Substring_match_length(substring) +
					  (querylength - Substring_queryend(substring)),
					  /*querypos*/Substring_querystart(substring),querylength,
					  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	    } else {
	      substring_start = Substring_querystart(substring);
	      substring_length = Substring_match_length(substring);
	      M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,Substring_match_length(substring),
					  /*querypos*/Substring_querystart(substring),querylength,
					  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
			  /*querypos*/Substring_queryend(substring),querylength,
			  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	    }
	    *finalp = p;
	    *nextp = List_next(p);
	    
	  } else {
	    /* Middle substring, plus */
	    debug(printf("Middle substring, plus %d..%d\n",Substring_querystart((Substring_T) List_head(p)), 
			 Substring_queryend((Substring_T) List_head(p))));
	    substring_start = Substring_querystart(substring);
	    substring_length = Substring_match_length(substring);
	    
	    M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					Substring_match_length(substring),
					/*querypos*/Substring_querystart(substring),querylength,
					hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	  }
	  p = List_next(p);
	}
	
	if (p != *endp) {
	  if (M_printedp == false) {
	    debug1(printf("Skipping initial indel or splice\n"));

	  } else if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    debug1(printf("2. Pushing %dD\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	    *nindels += Junction_nindels(post_junction);

	  } else if (type == INS_JUNCTION) {
	    debug1(printf("2. Pushing %dI\n",Junction_nindels(post_junction)));
#if 0
	    FPRINTF(fp,"%dI",Junction_nindels(post_junction));
#else
	    print_cigar(fp,/*type*/'I',/*stringlength*/Junction_nindels(post_junction),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
#endif
	    *nindels += Junction_nindels(post_junction);

	  } else if (type == SPLICE_JUNCTION) {
	    debug1(printf("2. Pushing %dN\n",Junction_splice_distance(post_junction)));
	    FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	  }
	  
	  /* Last substring, plus, hard-clipped */
	  substring = (Substring_T) List_head(p);
	  debug(printf("Last substring, plus, hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		       Substring_queryend((Substring_T) List_head(p))));
	  if (hide_soft_clips_p == true) {
	    substring_start = Substring_querystart_pretrim(substring);
	    substring_length = Substring_match_length_pretrim(substring);
	    M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					Substring_match_length(substring) +
					(querylength - Substring_queryend(substring)),
					/*querypos*/Substring_querystart(substring),querylength,
					hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	  } else {
	    substring_start = Substring_querystart(substring);
	    substring_length = Substring_match_length(substring);
	    M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					Substring_match_length(substring),
					/*querypos*/Substring_querystart(substring),querylength,
					hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	  }
	  *finalp = p;
	  *nextp = List_next(p);
	  
	}
      }
    }

  } else {
    /* Minus */
    while (p != *endp && Substring_querystart((Substring_T) List_head(p)) >= querylength - hardclip_low) {
      /* Skip, because substring entirely in hard-clipped region */
      debug(printf("Skipping %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));
      *prevp = p;
      p = List_next(p);
      q = List_next(q);
    }

    if (p == NULL) {
      debug(printf("Empty substring\n"));
      FPRINTF(fp,"*");

    } else {
      substring = (Substring_T) List_head(p);
      if (List_next(p) == *endp || Substring_querystart(substring) < hardclip_high) {
	/* Single substring */
	debug(printf("Single substring %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	
	if (hide_soft_clips_p == true) {
	  substring_start = Substring_querystart_pretrim(substring);
	  substring_length = Substring_match_length_pretrim(substring);
	  M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
				      (querylength - Substring_queryend(substring)) + 
				      Substring_match_length(substring) + Substring_querystart(substring),
				      /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
				      /*plusp*/false,/*lastp*/true,/*trimlength*/0);
	} else {
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);
	  print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		      /*querypos*/querylength,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	  M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
				      Substring_match_length(substring),
				      /*querypos*/Substring_queryend(substring),querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	  print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		      /*querypos*/Substring_querystart(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	}
	*finalp = p;
	*nextp = List_next(p);
	
      } else {
	/* First substring, minus */
	debug(printf("First substring, minus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		     Substring_queryend((Substring_T) List_head(p))));
	
	post_junction = (Junction_T) List_head(q);
	
	if (hide_soft_clips_p == true) {
	  substring_start = Substring_querystart_pretrim(substring);
	  substring_length = Substring_match_length_pretrim(substring);
	  M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
				      (querylength - Substring_queryend(substring)) +
				      Substring_match_length(substring),
				      /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
				      /*plusp*/false,/*lastp*/false,/*trimlength*/0);
	} else {
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);
	  print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		      /*querypos*/querylength,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	  M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
				      Substring_match_length(substring),
				      /*querypos*/Substring_queryend(substring),querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	}
	p = List_next(p);
	
	while (p != *endp && Substring_querystart((Substring_T) List_head(p)) >= hardclip_high) {
	  if (M_printedp == false) {
	    debug1(printf("Skipping initial indel or splice\n"));

	  } else if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    debug1(printf("3. Pushing %dD\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	    *nindels += Junction_nindels(post_junction);

	  } else if (type == INS_JUNCTION) {
	    debug1(printf("3. Pushing %dI\n",Junction_nindels(post_junction)));
#if 0
	    FPRINTF(fp,"%dI",Junction_nindels(post_junction));
#else
	    print_cigar(fp,/*type*/'I',/*stringlength*/Junction_nindels(post_junction),
			/*querypos*/Substring_querystart(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
#endif
	    *nindels += Junction_nindels(post_junction);

	  } else if (type == SPLICE_JUNCTION) {
	    debug1(printf("3. Pushing %dN\n",Junction_splice_distance(post_junction)));
	    FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	  }
	  q = List_next(q);
	  if (q == NULL) {
	  } else {
	    post_junction = (Junction_T) List_head(q);
	  }
	  
	  substring = (Substring_T) List_head(p);
	  if (List_next(p) == *endp) {
	    /* Last substring, minus, not hard-clipped */
	    debug(printf("Last substring, minus, not hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
			 Substring_queryend((Substring_T) List_head(p))));
	    
	    if (hide_soft_clips_p == true) {
	      substring_start = Substring_querystart_pretrim(substring);
	      substring_length = Substring_match_length_pretrim(substring);
	      M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					  Substring_match_length(substring) +
					  Substring_querystart(substring),
					  /*querypos*/Substring_queryend(substring),querylength,
					  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	    } else {
	      substring_start = Substring_querystart(substring);
	      substring_length = Substring_match_length(substring);
	      M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					  Substring_match_length(substring),
					  /*querypos*/Substring_queryend(substring),querylength,
					  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	      print_cigar(fp,/*type*/'S',Substring_querystart(substring),
			  /*querypos*/Substring_querystart(substring),querylength,hardclip_low,hardclip_high,
			  /*plusp*/false,/*lastp*/true,/*trimlength*/0);
	    }
	    *finalp = p;
	    *nextp = List_next(p);
	    
	  } else {
	    /* Middle substring, minus */
	    debug(printf("Middle substring, minus %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
			 Substring_queryend((Substring_T) List_head(p))));
	    substring_start = Substring_querystart(substring);
	    substring_length = Substring_match_length(substring);
	    
	    M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					Substring_match_length(substring),
					/*querypos*/Substring_queryend(substring),querylength,
					hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	  }
	  p = List_next(p);
	}
	
	if (p != *endp) {
	  if (M_printedp == false) {
	    debug1(printf("Skipping initial indel or splice\n"));

	  } else if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	    debug1(printf("4. Pushing %dD\n",Junction_nindels(post_junction)));
	    FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	    *nindels += Junction_nindels(post_junction);

	  } else if (type == INS_JUNCTION) {
	    debug1(printf("4. Pushing %dI\n",Junction_nindels(post_junction)));
#if 0
	    FPRINTF(fp,"%dI",Junction_nindels(post_junction));
#else
	    print_cigar(fp,/*type*/'I',/*stringlength*/Junction_nindels(post_junction),
			/*querypos*/Substring_querystart(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
#endif
	    *nindels += Junction_nindels(post_junction);

	  } else if (type == SPLICE_JUNCTION) {
	    debug1(printf("4. Pushing %dN\n",Junction_splice_distance(post_junction)));
	    FPRINTF(fp,"%uN",Junction_splice_distance(post_junction));
	  }
	  
	  /* Last substring, minus, hard-clipped */
	  substring = (Substring_T) List_head(p);
	  debug(printf("Last substring, minus, hard-clipped %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		       Substring_queryend((Substring_T) List_head(p))));
	  
	  if (hide_soft_clips_p == true) {
	    substring_start = Substring_querystart_pretrim(substring);
	    substring_length = Substring_match_length_pretrim(substring);
	    M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					Substring_match_length(substring) +
					Substring_querystart(substring),
					/*querypos*/Substring_queryend(substring),querylength,
					hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	  } else {
	    substring_start = Substring_querystart(substring);
	    substring_length = Substring_match_length(substring);
	    M_printedp |= print_cigar_M(fp,substring,substring_length,substring_start,
					Substring_match_length(substring),
					/*querypos*/Substring_queryend(substring),querylength,
					hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	    print_cigar(fp,/*type*/'S',Substring_querystart(substring),
			/*querypos*/Substring_querystart(substring),querylength,hardclip_low,hardclip_high,
			/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	  }
	  *finalp = p;
	  *nextp = List_next(p);
	  
	}
      }
    }
  }

  return;
}


static void
Cigar_print_halfdonor (int *part_hardclip_low, int *part_hardclip_high,
		       Filestring_T fp, Substring_T donor, Stage3end_T this,
		       int querylength) {
  bool sensep;
  int hardclip_low, hardclip_high;
  bool plusp;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif


  plusp = Substring_plusp(donor);

  if (Stage3end_sensedir(this) == SENSE_ANTI) {
    sensep = false;
  } else {
    sensep = true;
  }

  if (sensep == true) {
    if (plusp == true) {
      hardclip_low = 0;
      hardclip_high = querylength - Substring_queryend(donor);
    } else {
      hardclip_high = 0;
      hardclip_low = querylength - Substring_queryend(donor);
    }
    
  } else {
    if (plusp == true) {
      hardclip_high = 0;
      hardclip_low = Substring_querystart(donor);
    } else {
      hardclip_low = 0;
      hardclip_high = Substring_querystart(donor);
    }
  }

  *part_hardclip_low = hardclip_low;
  *part_hardclip_high = hardclip_high;

  if (sensep == true) {
    /* Doesn't hold for DNA-Seq chimeras */
    /* assert(Substring_siteD_pos(donor) == Substring_queryend(donor)); */
    if (plusp == true) {
      /* sensep true, plusp true */
      /* FPRINTF(fp,"donor sensep true, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    Substring_querystart(donor) + 
		    Substring_match_length(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_queryend(donor));

      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_queryend(donor));
      }

    } else {
      /* sensep true, plusp false */
      /* FPRINTF(fp,"donor sensep false, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_queryend(donor));
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(donor) +
		    Substring_querystart(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_queryend(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/0);
      }
    }

  } else {
    /* Doesn't hold for DNA-Seq chimeras */
    /* assert(Substring_siteD_pos(donor) == Substring_querystart(donor)); */
    if (plusp == true) {
      /* sensep false, plusp true */
      /* FPRINTF(fp,"donor sensep false, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_querystart(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor) + (querylength - Substring_queryend(donor)),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_querystart(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      }

    } else {
      /* sensep false, plusp false */
      /* FPRINTF(fp,"donor sensep true, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',(querylength - Substring_queryend(donor)) + Substring_match_length(donor),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/Substring_trim_querystart(donor));

      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/Substring_trim_querystart(donor));
      }
    }
  }

  return;
}


void
Cigar_print_halfacceptor (int *part_hardclip_low, int *part_hardclip_high,
			  Filestring_T fp, Substring_T acceptor, Stage3end_T this,
			  int querylength) {
  bool sensep;
  int hardclip_low, hardclip_high;
  bool plusp;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif


  plusp = Substring_plusp(acceptor);

  if (Stage3end_sensedir(this) == SENSE_ANTI) {
    sensep = false;
  } else {
    sensep = true;
  }

  if (sensep == true) {
    if (plusp == true) {
      hardclip_high = 0;
      hardclip_low = Substring_querystart(acceptor);
    } else {
      hardclip_low = 0;
      hardclip_high = Substring_querystart(acceptor);
    }
    
  } else {
    if (plusp == true) {
      hardclip_low = 0;
      hardclip_high = querylength - Substring_queryend(acceptor);
    } else {
      hardclip_high = 0;
      hardclip_low = querylength - Substring_queryend(acceptor);
    }
  }

  *part_hardclip_low = hardclip_low;
  *part_hardclip_high = hardclip_high;

  if (sensep == true) {
    /* Doesn't hold for DNA-Seq chimeras */
    /* assert(Substring_siteA_pos(acceptor) == Substring_querystart(acceptor)); */
    if (plusp == true) {
      /* sensep true, plusp true */
      /* FPRINTF(fp,"acceptor sensep true, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',Substring_querystart(acceptor) + Substring_match_length(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_queryend(acceptor));
      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_queryend(acceptor));
      }

    } else {
      /* sensep true, plusp false */
      /* FPRINTF(fp,"acceptor sensep true, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_queryend(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor) + Substring_querystart(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_queryend(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/0);
      }
    }

  } else {
    /* sensep false, plusp true */
    /* Doesn't hold for DNA-Seq chimeras */
    /* assert(Substring_siteA_pos(acceptor) == Substring_queryend(acceptor)); */
    if (plusp == true) {
      /* FPRINTF(fp,"acceptor sensep false, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_querystart(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor) + (querylength - Substring_queryend(acceptor)),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_querystart(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_queryend(acceptor));
      }

    } else {
      /* sensep false, plusp false */
      /* FPRINTF(fp,"acceptor sensep false, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',(querylength - Substring_queryend(acceptor)) + Substring_match_length(acceptor),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/Substring_trim_querystart(acceptor));
      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/Substring_trim_querystart(acceptor));
      }
    }
  }

  return;
}


Filestring_T
Cigar_compute_main (int *part_hardclip_low, int *part_hardclip_high,
		    int *nindels, List_T *startp, List_T *startq,
		    List_T *prevp, List_T *nextp, List_T *finalp, List_T *endp,
		    bool *plusp, Chrnum_T *chrnum, Chrpos_T *chrpos,
		    Stage3end_T stage3end, int querylength, bool first_read_p, Stage3end_T mate,
		    int hardclip_low, int hardclip_high, bool hide_soft_clips_p) {
  Filestring_T cigar_fp;
  Hittype_T hittype;
  int circularpos;
  int length_querystart, length_queryend;
  Substring_T substring, donor, acceptor;

  debug3(printf("Entered Cigar_compute_main with hardclip_low %d and hardclip_high %d\n",
		hardclip_low,hardclip_high));

  cigar_fp = Filestring_new(/*id*/0);

  if (stage3end == NULL) {
    debug3(printf("stage3end is NULL\n"));
    FPRINTF(cigar_fp,"*");
    *part_hardclip_low = hardclip_low;
    *part_hardclip_high = hardclip_high;
    *plusp = true;
    *chrnum = 0;
    *chrpos = 0;

  } else if ((hittype = Stage3end_hittype(stage3end)) == TRANSLOC_SPLICE ||
	     (hittype == SAMECHR_SPLICE && merge_samechr_p == false)) {
    debug3(printf("Have a translocation\n"));

    if (mate == NULL) {
      /* No notion of concordance, so pick longer side to be main */
      debug3(printf("Mate is NULL, so pick longer side\n"));

      donor = Stage3end_substring_donor(stage3end);
      acceptor = Stage3end_substring_acceptor(stage3end);
      if (Substring_match_length(donor) > Substring_match_length(acceptor)) {
	/* Donor is longer, so concordant and main */
	Cigar_print_halfdonor(&(*part_hardclip_low),&(*part_hardclip_high),
			      cigar_fp,donor,stage3end,querylength);
	*plusp = Substring_plusp(donor);
	*chrnum = Substring_chrnum(donor);
	*chrpos = Substring_compute_chrpos(donor,*part_hardclip_low,
					   /**part_hardclip_high,*/hide_soft_clips_p);

      } else {
	/* Acceptor is longer, so concordant and main */
	Cigar_print_halfacceptor(&(*part_hardclip_low),&(*part_hardclip_high),
				 cigar_fp,acceptor,stage3end,querylength);
	*plusp = Substring_plusp(acceptor);
	*chrnum = Substring_chrnum(acceptor);
	*chrpos = Substring_compute_chrpos(acceptor,*part_hardclip_low,
					   /**part_hardclip_high,*/hide_soft_clips_p);
      }

    } else if (Stage3end_donor_concordant_p(stage3end,first_read_p) == true) {
      /* donor is concordant and main */
      substring = Stage3end_substring_donor(stage3end);
      debug3(printf("Donor is concordant and main with plusp %d\n",Substring_plusp(substring)));

      Cigar_print_halfdonor(&(*part_hardclip_low),&(*part_hardclip_high),
			    cigar_fp,/*donor*/substring,stage3end,querylength);
      *plusp = Substring_plusp(substring);
      *chrnum = Substring_chrnum(substring);
      *chrpos = Substring_compute_chrpos(substring,*part_hardclip_low,
					 /**part_hardclip_high,*/hide_soft_clips_p);
    } else {
      /* acceptor is concordant and main */
      substring = Stage3end_substring_acceptor(stage3end);
      debug3(printf("Acceptor is concordant and main, with plusp %d\n",Substring_plusp(substring)));

      Cigar_print_halfacceptor(&(*part_hardclip_low),&(*part_hardclip_high),
			       cigar_fp,/*acceptor*/substring,stage3end,querylength);
      *plusp = Substring_plusp(substring);
      *chrnum = Substring_chrnum(substring);
      *chrpos = Substring_compute_chrpos(substring,*part_hardclip_low,
					 /**part_hardclip_high,*/hide_soft_clips_p);
    }

  } else if ((circularpos = Stage3end_circularpos(stage3end)) <= 0) {
    Cigar_print_substrings(&(*nindels),&(*startp),&(*startq),
			   &(*prevp),&(*nextp),&(*finalp),&(*endp),
			   cigar_fp,stage3end,querylength,hardclip_low,hardclip_high);

    /* Always want the substring with the lowest coordinate */
    if ((substring = Stage3end_substring_low(stage3end,hardclip_low)) == NULL) {
      *part_hardclip_low = *part_hardclip_high = 0;
      *plusp = true;
      *chrnum = 0;
      *chrpos = 0;
    } else {
      *part_hardclip_low = hardclip_low;
      *part_hardclip_high = hardclip_high;
      *plusp = Substring_plusp(substring);
      *chrnum = Substring_chrnum(substring);
      *chrpos = Substring_compute_chrpos(substring,hardclip_low,/*hardclip_high,*/hide_soft_clips_p);
    }

  } else {
    /* Circular alignment: Pick longer side to be the main */
    debug3(printf("Have a circular alignment with circularpos %d\n",circularpos));
    debug3(printf("For length_querystart, will enforce hardclips of %d and %d\n",
		  hardclip_low,querylength-circularpos));
    debug3(printf("For length_queryend, will enforce hardclips of %d and %d\n",
		  circularpos,hardclip_high));

    length_querystart = Cigar_length_substrings(stage3end,querylength,
						hardclip_low,/*hardclip_high*/querylength-circularpos);
    length_queryend = Cigar_length_substrings(stage3end,querylength,
					      /*hardclip_low*/circularpos,hardclip_high);
    debug3(printf("length_querystart is %d, length_queryend is %d\n",length_querystart,length_queryend));

    if (length_querystart > length_queryend) {
      debug3(printf("querystart is longer\n"));
      Cigar_print_substrings(&(*nindels),&(*startp),&(*startq),
			     &(*prevp),&(*nextp),&(*finalp),&(*endp),
			     cigar_fp,stage3end,querylength,
			     hardclip_low,/*hardclip_high*/querylength-circularpos);

      if ((substring = Stage3end_substring_low(stage3end,hardclip_low)) == NULL) {
	*part_hardclip_low = *part_hardclip_high = 0;
	*plusp = true;
	*chrnum = 0;
	*chrpos = 0;
      } else {
	*part_hardclip_low = hardclip_low;
	*part_hardclip_high = querylength - circularpos;
	*plusp = Substring_plusp(substring);
	*chrnum = Substring_chrnum(substring);
	*chrpos = Substring_compute_chrpos(substring,hardclip_low,
					   /*hardclip_high:querylength-circularpos,*/
					   hide_soft_clips_p);
	if (*chrpos > Substring_chrlength(substring)) {
	  *chrpos -= Substring_chrlength(substring);
	}
      }

    } else {
      debug3(printf("queryend is longer\n"));
      Cigar_print_substrings(&(*nindels),&(*startp),&(*startq),
			     &(*prevp),&(*nextp),&(*finalp),&(*endp),
			     cigar_fp,stage3end,querylength,
			     /*hardclip_low*/circularpos,hardclip_high);

      if ((substring = Stage3end_substring_low(stage3end,/*hardclip_low*/circularpos)) == NULL) {
	*part_hardclip_low = *part_hardclip_high = 0;
	*plusp = true;
	*chrnum = 0;
	*chrpos = 0;
      } else {
	*part_hardclip_low = circularpos;
	*part_hardclip_high = hardclip_high;
	*plusp = Substring_plusp(substring);
	*chrnum = Substring_chrnum(substring);
	*chrpos = Substring_compute_chrpos(substring,/*hardclip_low*/circularpos,
					   /*hardclip_high,*/hide_soft_clips_p);
	if (*chrpos > Substring_chrlength(substring)) {
	  *chrpos -= Substring_chrlength(substring);
	}
      }
    }
  }
  debug3(printf("\n\n"));

  if (sam_hardclip_use_S_p == true) {
    *part_hardclip_low = *part_hardclip_high = 0;
  }

  debug3(printf("Exiting Cigar_compute_main with plusp %d, chrnum %d, chrpos %u\n",*plusp,*chrnum,*chrpos));
  Filestring_stringify(cigar_fp);
  return cigar_fp;
}


Filestring_T
Cigar_compute_supplemental (int *part_hardclip_low, int *part_hardclip_high,
			    int *nindels, List_T *startp, List_T *startq,
			    List_T *prevp, List_T *nextp, List_T *finalp, List_T *endp,
			    bool *plusp, Chrnum_T *chrnum, Chrpos_T *chrpos,
			    Stage3end_T stage3end, int querylength, bool first_read_p, Stage3end_T mate,
			    int hardclip_low, int hardclip_high, bool hide_soft_clips_p) {
  Filestring_T cigar_fp;
  Hittype_T hittype;
  int circularpos;
  int length_querystart, length_queryend;
  Substring_T substring, donor, acceptor;

  debug4(printf("Entered Cigar_compute_supplemental with hardclip_low %d and hardclip_high %d\n",
		hardclip_low,hardclip_high));

  cigar_fp = Filestring_new(/*id*/0);

  if (stage3end == NULL) {
    /* Shouldn't be calling for supplemental */
    abort();

  } else if ((hittype = Stage3end_hittype(stage3end)) == TRANSLOC_SPLICE ||
	     (hittype == SAMECHR_SPLICE && merge_samechr_p == false)) {
    debug4(printf("Have a translocation\n"));

    if (mate == NULL) {
      /* No notion of concordance, so pick shorter side to be supplemental */
      debug4(printf("Mate is NULL, so pick shorter side\n"));

      donor = Stage3end_substring_donor(stage3end);
      acceptor = Stage3end_substring_acceptor(stage3end);
      if (Substring_match_length(donor) > Substring_match_length(acceptor)) {
	/* Acceptor is shorter, so supplemental */
	Cigar_print_halfacceptor(&(*part_hardclip_low),&(*part_hardclip_high),
				 cigar_fp,acceptor,stage3end,querylength);
	*plusp = Substring_plusp(acceptor);
	*chrnum = Substring_chrnum(acceptor);
	*chrpos = Substring_compute_chrpos(acceptor,*part_hardclip_low,
					   /**part_hardclip_high,*/hide_soft_clips_p);

      } else {
	/* Donor is shorter, so supplemental */
	Cigar_print_halfdonor(&(*part_hardclip_low),&(*part_hardclip_high),
			      cigar_fp,donor,stage3end,querylength);
	*plusp = Substring_plusp(donor);
	*chrnum = Substring_chrnum(donor);
	*chrpos = Substring_compute_chrpos(donor,*part_hardclip_low,
					   /**part_hardclip_high,*/hide_soft_clips_p);
      }

    } else if (Stage3end_donor_concordant_p(stage3end,first_read_p) == true) {
      /* acceptor is other and supplemental */
      substring = Stage3end_substring_acceptor(stage3end);
      debug3(printf("Acceptor is concordant and other, with plusp %d\n",Substring_plusp(substring)));

      Cigar_print_halfacceptor(&(*part_hardclip_low),&(*part_hardclip_high),
			       cigar_fp,/*acceptor*/substring,stage3end,querylength);
      *plusp = Substring_plusp(substring);
      *chrnum = Substring_chrnum(substring);
      *chrpos = Substring_compute_chrpos(substring,*part_hardclip_low,
					 /**part_hardclip_high,*/hide_soft_clips_p);
    } else {
      /* donor is other and supplemental */
      substring = Stage3end_substring_donor(stage3end);
      debug3(printf("Donor is concordant and other, with plusp %d\n",Substring_plusp(substring)));

      Cigar_print_halfdonor(&(*part_hardclip_low),&(*part_hardclip_high),
			    cigar_fp,/*donor*/substring,stage3end,querylength);
      *plusp = Substring_plusp(substring);
      *chrnum = Substring_chrnum(substring);
      *chrpos = Substring_compute_chrpos(substring,*part_hardclip_low,
					 /**part_hardclip_high,*/hide_soft_clips_p);
    }
      
  } else if ((circularpos = Stage3end_circularpos(stage3end)) <= 0) {
    /* Shouldn't be calling for supplemental */
    abort();

  } else {
    /* Circular alignment: Pick shorter side to be the supplemental */
    debug4(printf("Have a circular alignment with circularpos %d\n",circularpos));
    debug4(printf("For length_querystart, will enforce hardclips of %d and %d\n",
		  hardclip_low,querylength-circularpos));
    debug4(printf("For length_queryend, will enforce hardclips of %d and %d\n",
		  circularpos,hardclip_high));

    length_querystart = Cigar_length_substrings(stage3end,querylength,
						hardclip_low,/*hardclip_high*/querylength-circularpos);
    length_queryend = Cigar_length_substrings(stage3end,querylength,
					      /*hardclip_low*/circularpos,hardclip_high);
    debug4(printf("length_querystart is %d, length_queryend is %d\n",length_querystart,length_queryend));

    if (length_querystart > length_queryend) {
      debug4(printf("queryend is shorter\n"));
      Cigar_print_substrings(&(*nindels),&(*startp),&(*startq),
			     &(*prevp),&(*nextp),&(*finalp),&(*endp),
			     cigar_fp,stage3end,querylength,
			     /*hardclip_low*/circularpos,hardclip_high);

      if ((substring = Stage3end_substring_low(stage3end,/*hardclip_low*/circularpos)) == NULL) {
	*part_hardclip_low = *part_hardclip_high = 0;
	*plusp = true;
	*chrnum = 0;
	*chrpos = 0;
      } else {
	*part_hardclip_low = circularpos;
	*part_hardclip_high = hardclip_high;
	*plusp = Substring_plusp(substring);
	*chrnum = Substring_chrnum(substring);
	*chrpos = Substring_compute_chrpos(substring,/*hardclip_low*/circularpos,
					   /*hardclip_high,*/hide_soft_clips_p);
	if (*chrpos > Substring_chrlength(substring)) {
	  *chrpos -= Substring_chrlength(substring);
	}
      }

    } else {
      debug4(printf("querystart is shorter\n"));
      Cigar_print_substrings(&(*nindels),&(*startp),&(*startq),
			     &(*prevp),&(*nextp),&(*finalp),&(*endp),
			     cigar_fp,stage3end,querylength,
			     hardclip_low,/*hardclip_high*/querylength-circularpos);

      if ((substring = Stage3end_substring_low(stage3end,hardclip_low)) == NULL) {
	*part_hardclip_low = *part_hardclip_high = 0;
	*plusp = true;
	*chrnum = 0;
	*chrpos = 0;
      } else {
	*part_hardclip_low = hardclip_low;
	*part_hardclip_high = querylength - circularpos;
	*plusp = Substring_plusp(substring);
	*chrnum = Substring_chrnum(substring);
	*chrpos = Substring_compute_chrpos(substring,hardclip_low,
					   /*hardclip_high:querylength-circularpos,*/
					   hide_soft_clips_p);
	if (*chrpos > Substring_chrlength(substring)) {
	  *chrpos -= Substring_chrlength(substring);
	}
      }
    }
  }
  debug4(printf("\n\n"));

  if (sam_hardclip_use_S_p == true) {
    *part_hardclip_low = *part_hardclip_high = 0;
  }

  debug3(printf("Exiting Cigar_compute_supplemental with plusp %d, chrnum %d, chrpos %u\n",*plusp,*chrnum,*chrpos));
  Filestring_stringify(cigar_fp);
  return cigar_fp;
}


