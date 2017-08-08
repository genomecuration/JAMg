static char rcsid[] = "$Id: cigar.c 207314 2017-06-14 19:28:08Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>		/* For isupper */

#include "cigar.h"
#include "mem.h"
#include "complement.h"



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

#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


static bool hide_soft_clips_p;
static bool cigar_extended_p = false;
static bool merge_samechr_p;
static bool md_lowercase_variant_p;



#if 0
static void
print_tokens_stdout (List_T tokens) {
  List_T p;
  char *token;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    printf("%s",token);
  }

  return;
}
#endif



#if 0
/* Derived from print_tokens_gff3 */
static void
print_tokens_sam (Filestring_T fp, List_T tokens) {
  List_T p;
  char *token;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FPRINTF(fp,"%s",token);
    FREE(token);
  }

  return;
}
#endif

#if 0
static List_T
push_token (List_T tokens, char *token) {
  char *copy;

  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);
  return List_push(tokens,(void *) copy);
}
#endif


#if 0
/* Currently used for insertions and deletions */
static List_T
compute_cigar_old (List_T tokens, char type, int stringlength, int querypos, int querylength,
		   int hardclip_low, int hardclip_high, bool plusp, bool firstp, bool lastp) {
  char token[10];
  
  debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plusp %d\n",
		type,stringlength,querypos,querylength,hardclip_low,hardclip_high,plusp));

  if (firstp == true) {
    debug1(printf("firstp is true\n"));
    if (plusp == true) {
      if (hardclip_low > 0) {
	sprintf(token,"%dH",hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (querypos > hardclip_low) {
	sprintf(token,"%dS",querypos - hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    } else {
      if (hardclip_high > 0) {
	sprintf(token,"%dH",hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (querypos < querylength - hardclip_high) {
	sprintf(token,"%dS",querypos - hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  if (type == 'D' || type == 'N') {
    if (querypos < hardclip_low || querypos >= querylength - hardclip_high) {
      stringlength = 0;
    }

  } else if (plusp == true) {
    debug1(printf("Comparing querypos %d..%d against %d..%d\n",
		  querypos,querypos + stringlength,hardclip_low,querylength - hardclip_high));
    if (/* querypos < hardclip_low && */querypos + stringlength < hardclip_low) {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 1: stringlength 0\n"));
    } else if (querypos < hardclip_low) {
      if (querypos + stringlength < querylength - hardclip_high) {
	/* Print part after hardclip_low */
	stringlength = (querypos + stringlength) - hardclip_low;
	debug1(printf("Case 2: stringlength %d\n",stringlength));
      } else {
	/* Print part between hardclip_low and hardclip_high */
	stringlength = (querylength - hardclip_high) - hardclip_low;
	debug1(printf("Case 3: stringlength %d\n",stringlength));
      }
    } else if (querypos < querylength - hardclip_high) {
      if (querypos + stringlength >= querylength - hardclip_high) {
	/* Print up to hardclip_high */
	stringlength = (querylength - hardclip_high) - querypos;
	debug1(printf("Case 4: stringlength %d\n",stringlength));
      } else {
	/* Print full stringlength */
	debug1(printf("Case 5: stringlength %d\n",stringlength));
      }
    } else {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 6: stringlength 0\n"));
    }

  } else {
    debug1(printf("Comparing querypos %d..%d against %d..%d\n",
		  querypos,querypos - stringlength,hardclip_low,querylength - hardclip_high));
    if (/* querypos >= querylength - hardclip_high && */ querypos - stringlength >= querylength - hardclip_high) {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 1: stringlength 0\n"));
    } else if (querypos >= querylength - hardclip_high) {
      if (querypos - stringlength >= hardclip_low) {
	/* Print part after hardclip_high */
	stringlength = (querylength - hardclip_high) - (querypos - stringlength);
	debug1(printf("Case 2: stringlength %d\n",stringlength));
      } else {
	/* Print part between hardclip_low and hardclip_high */
	stringlength = (querylength - hardclip_high) - hardclip_low;
	debug1(printf("Case 3: stringlength %d\n",stringlength));
      }
    } else if (querypos >= hardclip_low) {
      if (querypos - stringlength < hardclip_low) {
	/* Print up to hardclip_low */
	stringlength = querypos - hardclip_low;
	debug1(printf("Case 4: stringlength %d\n",stringlength));
      } else {
	/* Print full stringlength */
	debug1(printf("Case 5: stringlength %d\n",stringlength));
      }
    } else {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 5: stringlength 0\n"));
    }
  }

  if (stringlength > 0) {
    sprintf(token,"%d%c",stringlength,type);
    debug1(printf("Pushing token %s\n",token));
    tokens = push_token(tokens,token);
  }

  if (lastp == true) {
    debug1(printf("lastp is true\n"));
    if (plusp == true) {
      querypos += stringlength;
      if (querypos < querylength - 1 - hardclip_high) {
	sprintf(token,"%dS",querylength - 1 - hardclip_high - querypos);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (hardclip_high > 0) {
	sprintf(token,"%dH",hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    } else {
      querypos -= stringlength;
      if (querypos > hardclip_low) {
	sprintf(token,"%dS",hardclip_low - querypos);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (hardclip_low > 0) {
	sprintf(token,"%dH",hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  return tokens;
}
#endif


#if 0
/* Currently used for insertions and deletions */
static List_T
compute_cigar (List_T tokens, char type, int stringlength, int querypos, int querylength,
	       int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  char token[10];
  
  if (plusp == true) {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
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
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	sprintf(token,"%d%c",matchlength,type);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }

  } else {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
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
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	sprintf(token,"%d%c",matchlength,type);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  return tokens;
}
#endif


#if 0
/* Modified from compute_cigar */
static Intlist_T
compute_cigar_types_only (Intlist_T types, char type, int stringlength, int querypos, int querylength,
			  int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  
  if (plusp == true) {
    debug1(printf("\nEntering compute_cigar_types_only with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
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
	types = Intlist_push(types,'H');
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	types = Intlist_push(types,type);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
    }

  } else {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
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
	types = Intlist_push(types,'H');
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	types = Intlist_push(types,type);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
    }
  }

  return types;
}
#endif


static void
print_cigar (Filestring_T fp, char type, int stringlength, int querypos, int querylength,
	     int hardclip_low, int hardclip_high, bool plusp, bool lastp, int trimlength) {
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
	FPRINTF(fp,"%dH",cliplength);
      }
      matchlength = endpos - startpos;
      if (matchlength <= 0) {
	/* Skip */
      } else if (type != 'E') {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	FPRINTF(fp,"%d%c",matchlength,type);
      } else if (matchlength == trimlength) {
	debug1(printf("  Pushing %dS\n",matchlength));
	FPRINTF(fp,"%dS",matchlength);
      } else {
	debug1(printf("  Pushing %dH because matchlength %d != trimlength %d\n",
		      matchlength,matchlength,trimlength));
	FPRINTF(fp,"%dH",matchlength);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	FPRINTF(fp,"%dH",cliplength);
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
	FPRINTF(fp,"%dH",cliplength);
      }
      matchlength = startpos - endpos;
      if (matchlength <= 0) {
	/* Skip */
      } else if (type != 'E') {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	FPRINTF(fp,"%d%c",matchlength,type);
      } else if (matchlength == trimlength) {
	debug1(printf("  Pushing %dS\n",matchlength));
	FPRINTF(fp,"%dS",matchlength);
      } else {
	debug1(printf("  Pushing %dH because matchlength %d != trimlength %d\n",
		      matchlength,matchlength,trimlength));
	FPRINTF(fp,"%dH",matchlength);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	FPRINTF(fp,"%dH",cliplength);
      }
    }
  }

  return;
}





/* Based on print_md_string */
static void
print_extended_cigar (Filestring_T fp, char *genomicfwd_refdiff,
		      int stringlength, int querypos, int querylength,
		      int hardclip_low, int hardclip_high, bool plusp, bool lastp) {
  int nmatches = 0, nmismatches = 0;
  int starti, endi, i;
  bool hardclip_end_p = false;
  int cliplength, endpos;

  if (plusp == true) {
    debug2(printf("\nEntering print_extended_cigar with querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus: %s ref, %s both\n",
		  querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff,genomicfwd_bothdiff));
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
	  if (nmismatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%dX",nmismatches);
	    nmismatches = 0;
	    hardclip_end_p = false;
	  }
	  nmatches++;

	} else {
	  /* A true mismatch against both variants */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	  }
	  nmismatches++;
	}
      }

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  if (nmismatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%dX",nmismatches);
	    nmismatches = 0;
	    hardclip_end_p = false;
	  }
	  nmatches++;

#if 0
	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	  }
	  nmismatches++;
#endif

	} else {
	  /* A true mismatch against both variants */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	  }
	  nmismatches++;
	}
      }
    }

    if (nmatches > 0) {
      FPRINTF(fp,"%d=",nmatches);
    } else if (nmismatches > 0) {
      FPRINTF(fp,"%dX",nmismatches);
    }

    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	FPRINTF(fp,"%dH",cliplength);
      }
    }

  } else {
    debug2(printf("\nEntering print_extended_cigar with querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus: %s ref, %s both\n",
		  querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff,genomicfwd_bothdiff));
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
	  if (nmismatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%dX",nmismatches);
	    nmismatches = 0;
	    hardclip_end_p = false;
	  }
	  nmatches++;

	} else {
	  /* A true mismatch against both variants */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	  }
	  nmismatches++;
	}
      }

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  if (nmismatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%dX",nmismatches);
	    nmismatches = 0;
	    hardclip_end_p = false;
	  }
	  nmatches++;

#if 0
	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	  }
	  nmismatches++;
#endif

	} else {
	  /* A true mismatch against both variants */
	  if (nmatches > 0 || hardclip_end_p == true) {
	    FPRINTF(fp,"%d=",nmatches);
	    nmatches = 0;
	    hardclip_end_p = false;
	  }
	  nmismatches++;
	}
      }
    }

    if (nmatches > 0) {
      FPRINTF(fp,"%d=",nmatches);
    } else if (nmismatches > 0) {
      FPRINTF(fp,"%dX",nmismatches);
    }

    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	FPRINTF(fp,"%dH",cliplength);
      }
    }
  }

  return;
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


static void
print_cigar_M (Filestring_T fp, Substring_T substring, int substring_length, int substring_start,
	       int stringlength, int querypos, int querylength,
	       int hardclip_low, int hardclip_high, bool plusp, bool lastp, int trimlength) {
  char *genomicfwd_refdiff, *genomicdir_refdiff;
  
  if (cigar_extended_p == false) {
    print_cigar(fp,/*type*/'M',stringlength,querypos,querylength,
		hardclip_low,hardclip_high,plusp,lastp,trimlength);
  } else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == NULL) {
    print_extended_cigar(fp,/*genomicfwd_refdiff*/NULL,/*stringlength*/substring_length,
			 /*querypos*/substring_start,querylength,
			 hardclip_low,hardclip_high,plusp,lastp);
  } else if (plusp == true) {
    print_extended_cigar(fp,&(genomicdir_refdiff[substring_start]),/*stringlength*/substring_length,
			 /*querypos*/substring_start,querylength,
			 hardclip_low,hardclip_high,plusp,lastp);
  } else {
    genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
    print_extended_cigar(fp,genomicfwd_refdiff,/*stringlength*/substring_length,
			 /*querypos*/substring_start,querylength,
			 hardclip_low,hardclip_high,plusp,lastp);
    FREEA(genomicfwd_refdiff);
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



void
Cigar_print_substrings (int *nindels, List_T *startp, List_T *startq, List_T *prevp, List_T *nextp, List_T *finalp, List_T *endp,
			Filestring_T fp, Stage3end_T stage3end,
			int querylength, int hardclip_low, int hardclip_high) {
  Substring_T substring, substringL, substringH;
  Junction_T post_junction;
  int type;

  List_T substrings_LtoH, junctions_LtoH;
  List_T p, q;
  int substring_start, substring_length;

  bool plusp;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif

  
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


  if (Substring_ambiguous_p(substringL) == true) {
    *prevp = substrings_LtoH;
    *startp = List_next(substrings_LtoH);
    *startq = List_next(junctions_LtoH);
  } else {
    *prevp = (List_T) NULL;
    *startp = substrings_LtoH;
    *startq = junctions_LtoH;
  }
  if (Substring_ambiguous_p(substringH) == true) {
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

    substring = (Substring_T) List_head(p);
    if (List_next(p) == *endp ||	Substring_queryend(substring) >= querylength - hardclip_high) {
      /* Single substring */
      debug(printf("Single substring %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));

      if (hide_soft_clips_p == true) {
	substring_start = Substring_querystart_orig(substring);
	substring_length = Substring_match_length_orig(substring);
	print_cigar_M(fp,substring,substring_length,substring_start,
		      Substring_querystart(substring) + Substring_match_length(substring) +
		      (querylength - Substring_queryend(substring)),/*querypos*/0,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
      } else {
	substring_start = Substring_querystart(substring);
	substring_length = Substring_match_length(substring);
	print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar_M(fp,substring,substring_length,substring_start,
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
	substring_start = Substring_querystart_orig(substring);
	substring_length = Substring_match_length_orig(substring);
	print_cigar_M(fp,substring,substring_length,substring_start,
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
	print_cigar_M(fp,substring,substring_length,substring_start,
		      Substring_match_length(substring),
		      /*querypos*/Substring_querystart(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
      }
      p = List_next(p);
      
      while (p != *endp && Substring_queryend((Substring_T) List_head(p)) < querylength - hardclip_high) {
	if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	  debug1(printf("1. Pushing %dD\n",Junction_nindels(post_junction)));
	  FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	  *nindels += Junction_nindels(post_junction);
	} else if (type == INS_JUNCTION) {
	  debug1(printf("1. Pushing %dI\n",Junction_nindels(post_junction)));
	  FPRINTF(fp,"%dI",Junction_nindels(post_junction));
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
	    substring_start = Substring_querystart_orig(substring);
	    substring_length = Substring_match_length_orig(substring);
	    print_cigar_M(fp,substring,substring_length,substring_start,
			  Substring_match_length(substring) +
			  (querylength - Substring_queryend(substring)),
			  /*querypos*/Substring_querystart(substring),querylength,
			  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	  } else {
	    substring_start = Substring_querystart(substring);
	    substring_length = Substring_match_length(substring);
	    print_cigar_M(fp,substring,substring_length,substring_start,Substring_match_length(substring),
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

	  print_cigar_M(fp,substring,substring_length,substring_start,
			Substring_match_length(substring),
			/*querypos*/Substring_querystart(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	}
	p = List_next(p);
      }
      
      if (p != *endp) {
	if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	  debug1(printf("2. Pushing %dD\n",Junction_nindels(post_junction)));
	  FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	  *nindels += Junction_nindels(post_junction);
	} else if (type == INS_JUNCTION) {
	  debug1(printf("2. Pushing %dI\n",Junction_nindels(post_junction)));
	  FPRINTF(fp,"%dI",Junction_nindels(post_junction));
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
	  substring_start = Substring_querystart_orig(substring);
	  substring_length = Substring_match_length_orig(substring);
	  print_cigar_M(fp,substring,substring_length,substring_start,
			Substring_match_length(substring) +
			(querylength - Substring_queryend(substring)),
			/*querypos*/Substring_querystart(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true,/*trimlength*/0);
	} else {
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);
	  print_cigar_M(fp,substring,substring_length,substring_start,
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

    substring = (Substring_T) List_head(p);
    if (List_next(p) == *endp || Substring_querystart(substring) < hardclip_high) {
      /* Single substring */
      debug(printf("Single substring %d..%d\n",Substring_querystart((Substring_T) List_head(p)),
		   Substring_queryend((Substring_T) List_head(p))));

      if (hide_soft_clips_p == true) {
	substring_start = Substring_querystart_orig(substring);
	substring_length = Substring_match_length_orig(substring);
	print_cigar_M(fp,substring,substring_length,substring_start,
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
	print_cigar_M(fp,substring,substring_length,substring_start,
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
	substring_start = Substring_querystart_orig(substring);
	substring_length = Substring_match_length_orig(substring);
	print_cigar_M(fp,substring,substring_length,substring_start,
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
	print_cigar_M(fp,substring,substring_length,substring_start,
		      Substring_match_length(substring),
		      /*querypos*/Substring_queryend(substring),querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
      }
      p = List_next(p);

      while (p != *endp && Substring_querystart((Substring_T) List_head(p)) >= hardclip_high) {
	if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	  debug1(printf("3. Pushing %dD\n",Junction_nindels(post_junction)));
	  FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	  *nindels += Junction_nindels(post_junction);
	} else if (type == INS_JUNCTION) {
	  debug1(printf("3. Pushing %dI\n",Junction_nindels(post_junction)));
	  FPRINTF(fp,"%dI",Junction_nindels(post_junction));
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
	    substring_start = Substring_querystart_orig(substring);
	    substring_length = Substring_match_length_orig(substring);
	    print_cigar_M(fp,substring,substring_length,substring_start,
			  Substring_match_length(substring) +
			  Substring_querystart(substring),
			  /*querypos*/Substring_queryend(substring),querylength,
			  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	  } else {
	    substring_start = Substring_querystart(substring);
	    substring_length = Substring_match_length(substring);
	    print_cigar_M(fp,substring,substring_length,substring_start,
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

	  print_cigar_M(fp,substring,substring_length,substring_start,
			Substring_match_length(substring),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false,/*trimlength*/0);
	}
	p = List_next(p);
      }

      if (p != *endp) {
	if ((type = Junction_type(post_junction)) == DEL_JUNCTION) {
	  debug1(printf("4. Pushing %dD\n",Junction_nindels(post_junction)));
	  FPRINTF(fp,"%dD",Junction_nindels(post_junction));
	  *nindels += Junction_nindels(post_junction);
	} else if (type == INS_JUNCTION) {
	  debug1(printf("4. Pushing %dI\n",Junction_nindels(post_junction)));
	  FPRINTF(fp,"%dI",Junction_nindels(post_junction));
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
	  substring_start = Substring_querystart_orig(substring);
	  substring_length = Substring_match_length_orig(substring);
	  print_cigar_M(fp,substring,substring_length,substring_start,
			Substring_match_length(substring) +
			Substring_querystart(substring),
			/*querypos*/Substring_queryend(substring),querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true,/*trimlength*/0);
	} else {
	  substring_start = Substring_querystart(substring);
	  substring_length = Substring_match_length(substring);
	  print_cigar_M(fp,substring,substring_length,substring_start,
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

  return;
}


void
Cigar_print_halfdonor (Filestring_T fp, Substring_T donor, Stage3end_T this,
		       int querylength, int *hardclip_low, int *hardclip_high,
		       bool use_hardclip_p) {
  bool sensep;
  int transloc_hardclip_low, transloc_hardclip_high;
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

  if (use_hardclip_p == true) {
    if (sensep == true) {
      if (plusp == true) {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = querylength - Substring_queryend(donor);
      } else {
	transloc_hardclip_high = 0;
	transloc_hardclip_low = querylength - Substring_queryend(donor);
      }

    } else {
      if (plusp == true) {
	transloc_hardclip_high = 0;
	transloc_hardclip_low = Substring_querystart(donor);
      } else {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = Substring_querystart(donor);
      }
    }

    if (transloc_hardclip_low > *hardclip_low) {
      *hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > *hardclip_high) {
      *hardclip_high = transloc_hardclip_high;
    }
  }


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
		    /*querypos*/0,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(donor));

      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		    /*querypos*/0,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(donor));
      }

    } else {
      /* sensep true, plusp false */
      /* FPRINTF(fp,"donor sensep false, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_right(donor));
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(donor) +
		    Substring_querystart(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_right(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,*hardclip_low,*hardclip_high,
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
		    /*querypos*/0,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_left(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor) + (querylength - Substring_queryend(donor)),
		    /*querypos*/Substring_querystart(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/0,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_left(donor));
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      }

    } else {
      /* sensep false, plusp false */
      /* FPRINTF(fp,"donor sensep true, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',(querylength - Substring_queryend(donor)) + Substring_match_length(donor),
		    /*querypos*/querylength,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/Substring_trim_left(donor));

      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/Substring_trim_left(donor));
      }
    }
  }

  return;
}


void
Cigar_print_halfacceptor (Filestring_T fp, Substring_T acceptor, Stage3end_T this,
			  int querylength, int *hardclip_low, int *hardclip_high,
			  bool use_hardclip_p) {
  bool sensep;
  int transloc_hardclip_low, transloc_hardclip_high;
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

  if (use_hardclip_p == true) {
    if (sensep == true) {
      if (plusp == true) {
	transloc_hardclip_high = 0;
	transloc_hardclip_low = Substring_querystart(acceptor);
      } else {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = Substring_querystart(acceptor);
      }

    } else {
      if (plusp == true) {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = querylength - Substring_queryend(acceptor);
      } else {
	transloc_hardclip_high = 0;
	transloc_hardclip_low = querylength - Substring_queryend(acceptor);
      }
    }

    if (transloc_hardclip_low > *hardclip_low) {
      *hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > *hardclip_high) {
      *hardclip_high = transloc_hardclip_high;
    }
  }


  if (sensep == true) {
    /* Doesn't hold for DNA-Seq chimeras */
    /* assert(Substring_siteA_pos(acceptor) == Substring_querystart(acceptor)); */
    if (plusp == true) {
      /* sensep true, plusp true */
      /* FPRINTF(fp,"acceptor sensep true, plusp true\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',Substring_querystart(acceptor) + Substring_match_length(acceptor),
		    /*querypos*/0,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(acceptor));
      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(acceptor));
      }

    } else {
      /* sensep true, plusp false */
      /* FPRINTF(fp,"acceptor sensep true, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_right(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor) + Substring_querystart(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/Substring_trim_right(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,*hardclip_low,*hardclip_high,
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
		    /*querypos*/0,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_left(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor) + (querylength - Substring_queryend(acceptor)),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/0);
      } else {
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/true,/*lastp*/false,/*trimlength*/Substring_trim_left(acceptor));
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/true,/*lastp*/true,
		    /*trimlength*/Substring_trim_right(acceptor));
      }

    } else {
      /* sensep false, plusp false */
      /* FPRINTF(fp,"acceptor sensep false, plusp false\n"); */
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',(querylength - Substring_queryend(acceptor)) + Substring_match_length(acceptor),
		    /*querypos*/querylength,querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/false,/*lastp*/false,/*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/Substring_trim_left(acceptor));
      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    *hardclip_low,*hardclip_high,/*plusp*/false,/*lastp*/false,
		    /*trimlength*/0);
	print_cigar(fp,/*type*/'E',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,*hardclip_low,*hardclip_high,
		    /*plusp*/false,/*lastp*/true,/*trimlength*/Substring_trim_left(acceptor));
      }
    }
  }

  return;
}




static void
print_exon_exon_cigar (Filestring_T fp, Stage3end_T this, int querylength) {
  Substring_T donor, acceptor;
  int sensedir;

  /* Shouldn't have any overlap on a distant splice */
  int hardclip_low = 0, hardclip_high = 0;

  sensedir = Stage3end_sensedir(this);

  if (sensedir == SENSE_FORWARD) {
    donor = Stage3end_substring_donor(this);
    Cigar_print_halfdonor(fp,donor,this,querylength,&hardclip_low,&hardclip_high,/*use_hardclip_p*/true);

  } else if (Stage3end_sensedir(this) == SENSE_ANTI) {
    acceptor = Stage3end_substring_acceptor(this);
    Cigar_print_halfacceptor(fp,acceptor,this,querylength,&hardclip_low,&hardclip_high,/*use_hardclip_p*/true);

  } else {
    /* SENSE_NULL (DNA distant chimera) */
    acceptor = Stage3end_substring_acceptor(this);
    Cigar_print_halfacceptor(fp,acceptor,this,querylength,&hardclip_low,&hardclip_high,/*use_hardclip_p*/true);
  }

  return;
}



void
Cigar_print_mate (Filestring_T fp, Stage3end_T mate, int mate_querylength, int mate_hardclip_low, int mate_hardclip_high) {
  Hittype_T hittype;
  int nindels;
  List_T startp, startq, prevp, nextp, finalp, endp;


  if (mate == NULL) {
    FPRINTF(fp,"*");		/* CIGAR for nomapping */

  } else if ((hittype = Stage3end_hittype(mate)) == GMAP) {
    Pair_print_tokens(fp,Stage3end_cigar_tokens(mate));

  } else if (hittype == TRANSLOC_SPLICE || (hittype == SAMECHR_SPLICE && merge_samechr_p == false)) {
    print_exon_exon_cigar(fp,mate,mate_querylength);

  } else {
    Cigar_print_substrings(&nindels,&startp,&startq,&prevp,&nextp,&finalp,&endp,
			   fp,mate,mate_querylength,mate_hardclip_low,mate_hardclip_high);
  }

  return;
}




